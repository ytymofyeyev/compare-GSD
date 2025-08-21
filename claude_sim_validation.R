# Simulation script to validate results of compare_GSD.Rmd

library(dplyr)
library(ggplot2)
library(survival)
library(gsDesign)
library(parallel)
library(tidyr)
library(purrr)

# Source the utility functions used in the original analysis
source("R/utils_compare_GSD.R")
source("R/utils_adoptr_basic.R")

# Set seed for reproducibility
set.seed(12345)

#-----------------------------------------------------
# 1. Define simulation parameters from the Rmd file
#-----------------------------------------------------
study_name <- "STUDYNAME"
time_unit <- "mo"
hr_grid <- seq(0.6, 1.0, 0.1) # Use fewer points for simulation
alpha <- 0.025
n_sims <- 1000  # Number of simulations per scenario

# Primary outcome assumptions
primary_outcome <- list(
  list(
    pref = "",
    eta = -log(1 - 0.05) / 12,  # Dropout hazard
    lambdaC = -log(0.69) / 24,  # Control arm hazard
    hr = 0.60                   # Hazard ratio assumption
  )
)

# Define the designs
beta_GSD <- 1 - 0.903 # type II error for fixed GSD

# Fixed GSD design
design_GSD <- list(
  tag = 'GSD without ERE',
  ratio = 1,
  test.type = 4,
  timing = c(0.75, 1),
  beta = 1-0.903,
  R = 24,
  gamma = 540/24,
  sfu = gsDesign::sfPower,
  sfupar = 2,
  sfl = gsDesign::sfPower,
  sflpar = 10.2
)
class(design_GSD) <- "gsd_fixed"

# Adaptive GSD based on conditional power
design_aGSD <- design_GSD
design_aGSD$tag <- 'Adaptive GSD'
design_aGSD$hrCP <- 0.60
design_aGSD$maxinflation <- 1.30
design_aGSD$cpadj <- c(0.001, 0.8)
design_aGSD$betastar <- 0.20
design_aGSD$z2 <- gsDesign::z2NC
class(design_aGSD) <- "adaptive_mp"

# Discrete version of the adaptive GSD
design_aGSD_d <- design_GSD
design_aGSD_d$tag <- 'Adaptive GSD (discrete)'
design_aGSD_d$maxinflation <- 1.30
design_aGSD_d$z1_mesh_in_f1_e1 <- c(0, 0.6, 0.68, 0.8)
design_aGSD_d$I_ratio_to_n2max <- c(1, 0.925, 0.68, 0.455)
design_aGSD_d$approx_method <- "constant"
design_aGSD_d$order <- 12
class(design_aGSD_d) <- "adaptive_pwl"

# Combine all designs
design_spec <- list(
  design_GSD,
  design_aGSD,
  design_aGSD_d
)

#-----------------------------------------------------
# 2. Process designs using the same approach as in Rmd
#-----------------------------------------------------
designs <- 
  design_spec |>
  tidyr::expand_grid(primary_outcome) |>
  tidyr::unnest_wider(primary_outcome) |>
  exec_calc(alpha, hr_grid)

#-----------------------------------------------------
# 3. Define functions for data simulation and analysis
#-----------------------------------------------------

# Function to simulate survival data for a single trial
simulate_survival_trial <- function(design, hr, n_subjects) {
  # Extract parameters
  lambdaC <- design$lambdaC  # Control arm hazard
  eta <- design$eta          # Dropout hazard
  R <- design$R              # Recruitment duration
  
  # Calculate treatment arm hazard
  lambdaT <- hr * lambdaC
  
  # Simulate data
  # Allocate subjects to treatment arms
  treatment <- rep(c(0, 1), c(n_subjects/2, n_subjects/2))
  
  # Simulate recruitment times (uniform over recruitment period)
  recruitment_times <- runif(n_subjects, 0, R)
  
  # Simulate event times based on treatment assignment
  event_times <- numeric(n_subjects)
  for (i in 1:n_subjects) {
    if (treatment[i] == 0) {
      event_times[i] <- rexp(1, rate = lambdaC)
    } else {
      event_times[i] <- rexp(1, rate = lambdaT)
    }
  }
  
  # Simulate dropout times
  dropout_times <- rexp(n_subjects, rate = eta)
  
  # Determine observed times and event indicators
  observed_times <- pmin(event_times, dropout_times)
  event_indicator <- as.numeric(event_times <= dropout_times)
  
  # Create dataset
  data <- data.frame(
    id = 1:n_subjects,
    treatment = treatment,
    recruitment = recruitment_times,
    time = observed_times,
    event = event_indicator,
    calendar_time = recruitment_times + observed_times
  )
  
  return(data)
}

# Function to analyze data at a specific number of events
analyze_at_events <- function(data, n_events) {
  # Sort by calendar time
  data <- data[order(data$calendar_time), ]
  
  # Get data up to the specified number of events
  events_so_far <- cumsum(data$event)
  cutoff_index <- which(events_so_far >= n_events)[1]
  
  # Handle case when not enough events
  if (length(cutoff_index) == 0 || is.na(cutoff_index)) {
    return(list(
      z = NA,
      p = NA,
      hr = NA,
      events = sum(data$event),
      duration = max(data$calendar_time)
    ))
  }
  
  data_analysis <- data[1:cutoff_index, ]
  cutoff_time <- max(data_analysis$calendar_time)
  
  # Calculate observed time for each subject
  data$obs_time <- pmin(data$calendar_time, cutoff_time) - data$recruitment
  data$obs_event <- data$event * (data$calendar_time <= cutoff_time)
  
  # Check if we have enough events for meaningful analysis
  if (sum(data$obs_event) < 10) {
    return(list(
      z = NA,
      p = NA,
      hr = NA,
      events = sum(data$obs_event),
      duration = cutoff_time
    ))
  }
  
  # Perform log-rank test
  tryCatch({
    fit <- survival::coxph(
      Surv(obs_time, obs_event) ~ treatment, 
      data = data
    )
    
    # Get z-score and p-value
    z <- -summary(fit)$coefficients[1, 4]  # Negative because we want z > 0 to favor treatment
    p <- summary(fit)$coefficients[1, 5] / 2  # One-sided p-value
    hr <- exp(summary(fit)$coefficients[1, 1])
    
    return(list(
      z = z,
      p = p,
      hr = hr,
      events = sum(data$obs_event),
      duration = cutoff_time
    ))
  }, error = function(e) {
    return(list(
      z = NA,
      p = NA,
      hr = NA,
      events = sum(data$obs_event),
      duration = cutoff_time
    ))
  })
}

# Function to determine number of events for stage 2 based on z1 (for adaptive designs)
determine_events_stage2 <- function(design, z1, base_events) {
  if (inherits(design, "gsd_fixed")) {
    return(base_events[2])
  } else if (inherits(design, "adaptive_mp")) {
    # Implementation of Mehta & Pocock approach
    if (is.na(z1)) {
      return(base_events[2])  # Default to original plan if z1 is NA
    }
    
    if (z1 >= design$gsd$upper$bound[1]) {
      return(base_events[2])  # Early stop for efficacy
    } else if (z1 < design$gsd$lower$bound[1]) {
      return(base_events[2])  # Continue with original plan after futility
    } else {
      # Promising zone - adjust based on conditional power
      cp <- conditional_power(z1, design$hrCP, design$gsd$n.I[1], design$gsd$n.I[2])
      if (cp < design$cpadj[1]) {
        return(base_events[2])  # Unfavorable zone
      } else if (cp > design$cpadj[2]) {
        return(base_events[2])  # Favorable zone
      } else {
        # Promising zone - calculate events needed for conditional power = 1 - betastar
        events_needed <- calculate_events_for_cp(z1, design$hrCP, 
                                                design$gsd$n.I[1], 
                                                1 - design$betastar)
        return(min(events_needed, design$maxinflation * base_events[2]))
      }
    }
  } else if (inherits(design, "adaptive_pwl")) {
    # Discrete version using piecewise linear function
    if (is.na(z1)) {
      return(base_events[2])  # Default to original plan if z1 is NA
    }
    
    if (z1 >= design$gsd$upper$bound[1]) {
      return(base_events[2])  # Early stop for efficacy
    } else if (z1 < design$gsd$lower$bound[1]) {
      return(base_events[2])  # Continue with original plan after futility
    } else {
      # Use the piecewise linear function to determine inflation factor
      inflation_factor <- approxfun(
        x = design$z1_mesh_in_f1_e1,
        y = design$I_ratio_to_n2max,
        method = design$approx_method,
        rule = 2
      )(z1)
      
      # Calculate final number of events
      # The inflation is inversely related to the I_ratio_to_n2max value
      return(base_events[2] * (1 + (design$maxinflation - 1) * (1 - inflation_factor)))
    }
  }
}

# Helper functions for conditional power
conditional_power <- function(z1, hr, n1, n2) {
  theta <- -log(hr)
  information_fraction <- n1 / n2
  drift <- theta * sqrt(n2)
  
  cp <- 1 - pnorm((qnorm(1-alpha) - drift * sqrt(1 - information_fraction) - 
                     z1 * sqrt(information_fraction)) / 
                    sqrt(1 - information_fraction))
  return(cp)
}

calculate_events_for_cp <- function(z1, hr, n1, target_cp) {
  theta <- -log(hr)
  
  # Function to find n2 that gives target conditional power
  find_n2 <- function(n2) {
    information_fraction <- n1 / n2
    drift <- theta * sqrt(n2)
    
    cp <- 1 - pnorm((qnorm(1-alpha) - drift * sqrt(1 - information_fraction) - 
                       z1 * sqrt(information_fraction)) / 
                      sqrt(1 - information_fraction))
    return(cp - target_cp)
  }
  
  # Find n2 using root-finding
  n2_min <- n1  # Minimum possible n2
  n2_max <- n1 * 10  # Some large upper bound
  
  tryCatch({
    result <- uniroot(find_n2, c(n2_min, n2_max), tol = 0.0001)
    return(result$root)
  }, error = function(e) {
    # If root finding fails, return maximum inflation
    return(n1 * 10)
  })
}

# Function to run the combination test for final analysis
combination_test <- function(z1, z2, n1, n2, w1, w2) {
  if (is.na(z1) || is.na(z2)) {
    return(NA)
  }
  
  # Calculate the statistic for the incremental data
  if (n2 <= n1) {
    return(z1)  # If no additional events, use interim z-value
  }
  
  # Calculate the incremental z-statistic
  z2_star <- sqrt(n2/(n2-n1)) * z2 - sqrt(n1/(n2-n1)) * z1
  
  # Calculate the combined statistic
  z_comb <- (w1 * z1 + w2 * z2_star) / sqrt(w1^2 + w2^2)
  
  return(z_comb)
}

#-----------------------------------------------------
# 4. Run simulations for each design and HR scenario
#-----------------------------------------------------

simulate_design <- function(design, hr, n_sims) {
  cat("Simulating design:", design$tag, "with HR =", hr, "\n")
  
  # Calculate required sample size based on events
  n_subjects <- ceiling(2 * design$gsd$n.I[2] / 0.8)  # Assuming 80% event rate
  
  # Create empty results data frame
  results <- data.frame(
    sim_id = 1:n_sims,
    reject_ia = logical(n_sims),
    reject_fa = logical(n_sims),
    reject_overall = logical(n_sims),
    duration_ia = numeric(n_sims),
    duration_fa = numeric(n_sims),
    events_ia = numeric(n_sims),
    events_fa = numeric(n_sims)
  )
  
  # Base numbers of events from GSD
  base_events <- design$gsd$n.I
  
  # Run simulations
  for (i in 1:n_sims) {
    if (i %% 100 == 0) cat("  Simulation", i, "of", n_sims, "\n")
    
    # Simulate trial data
    trial_data <- simulate_survival_trial(design, hr, n_subjects)
    
    # Analyze at interim
    ia_analysis <- analyze_at_events(trial_data, base_events[1])
    results$events_ia[i] <- ia_analysis$events
    results$duration_ia[i] <- ia_analysis$duration
    
    # Check if we stop at interim for efficacy
    if (!is.na(ia_analysis$z) && ia_analysis$z >= design$gsd$upper$bound[1]) {
      results$reject_ia[i] <- TRUE
      results$reject_overall[i] <- TRUE
      results$events_fa[i] <- results$events_ia[i]
      results$duration_fa[i] <- results$duration_ia[i]
      next
    } else {
      results$reject_ia[i] <- FALSE
    }
    
    # Check if we stop at interim for futility (non-binding)
    if (!is.na(ia_analysis$z) && ia_analysis$z < design$gsd$lower$bound[1]) {
      # For simulation purposes, we'll continue but track futility
      futility_stop <- TRUE
    } else {
      futility_stop <- FALSE
    }
    
    # Determine number of events for stage 2
    events_stage2 <- determine_events_stage2(design, ia_analysis$z, base_events)
    
    # Analyze at final
    fa_analysis <- analyze_at_events(trial_data, events_stage2)
    results$events_fa[i] <- fa_analysis$events
    results$duration_fa[i] <- fa_analysis$duration
    
    # Apply combination test for final analysis
    w1 <- sqrt(base_events[1])
    w2 <- sqrt(base_events[2] - base_events[1])
    
    z_comb <- combination_test(
      ia_analysis$z, fa_analysis$z, 
      base_events[1], fa_analysis$events,
      w1, w2
    )
    
    # Check if we reject at final
    if (!futility_stop && !is.na(z_comb) && z_comb >= design$gsd$upper$bound[2]) {
      results$reject_fa[i] <- TRUE
      results$reject_overall[i] <- TRUE
    } else {
      results$reject_fa[i] <- FALSE
    }
  }
  
  # Calculate summary statistics
  summary <- list(
    design = design$tag,
    hr = hr,
    power_ia = mean(results$reject_ia, na.rm = TRUE),
    power_fa = mean(results$reject_fa, na.rm = TRUE),
    power_overall = mean(results$reject_overall, na.rm = TRUE),
    mean_duration_ia = mean(results$duration_ia, na.rm = TRUE),
    mean_duration_fa = mean(results$duration_fa, na.rm = TRUE),
    mean_events_ia = mean(results$events_ia, na.rm = TRUE),
    mean_events_fa = mean(results$events_fa, na.rm = TRUE),
    expected_duration = mean(results$duration_fa, na.rm = TRUE),
    expected_events = mean(results$events_fa, na.rm = TRUE)
  )
  
  return(list(results = results, summary = summary))
}

# Function to safely run simulations for a single design-HR combination
safe_simulate_design <- function(design_idx, hr_val, design_spec, primary_outcome, n_sims) {
  tryCatch({
    design <- design_spec[[design_idx]]
    # Add primary outcome parameters to design for simulation
    design <- c(design, primary_outcome[[1]])
    simulate_design(design, hr_val, n_sims)
  }, error = function(e) {
    cat("Error in simulation for design", design_spec[[design_idx]]$tag, 
        "with HR =", hr_val, ":", conditionMessage(e), "\n")
    # Return empty result structure
    list(
      results = data.frame(),
      summary = list(
        design = design_spec[[design_idx]]$tag,
        hr = hr_val,
        power_ia = NA,
        power_fa = NA,
        power_overall = NA,
        mean_duration_ia = NA,
        mean_duration_fa = NA,
        mean_events_ia = NA,
        mean_events_fa = NA,
        expected_duration = NA,
        expected_events = NA
      )
    )
  })
}

# Run simulations in parallel for all designs and HR values
sim_results <- expand.grid(
  design_index = 1:length(design_spec),
  hr = hr_grid,
  stringsAsFactors = FALSE
) %>%
  mutate(
    results = mcmapply(
      safe_simulate_design,
      design_index, hr,
      MoreArgs = list(
        design_spec = design_spec,
        primary_outcome = primary_outcome,
        n_sims = n_sims
      ),
      SIMPLIFY = FALSE,
      mc.cores = min(parallel::detectCores() - 1, 4)
    )
  )

# Extract summaries
sim_summaries <- sim_results %>%
  mutate(
    design_tag = sapply(design_index, function(idx) design_spec[[idx]]$tag),
    summary = lapply(results, function(x) x$summary)
  ) %>%
  select(design_tag, hr, summary) %>%
  unnest_wider(summary)

#-----------------------------------------------------
# 5. Compare simulation results with analytical calculations
#-----------------------------------------------------

# Function to extract analytical results from the designs object
extract_analytical_results <- function(designs, hr_grid) {
  analytical_results <- data.frame()
  
  for (i in 1:length(designs$design)) {
    design_tag <- designs$design[[i]]$tag
    
    # Extract power values for each HR
    for (hr in hr_grid) {
      hr_idx <- which(abs(designs$hr_grid - hr) < 1e-6)
      
      if (length(hr_idx) > 0) {
        result <- data.frame(
          design_tag = design_tag,
          hr = hr,
          power_ia = designs$power_ia[[i]][hr_idx],
          power_overall = designs$power[[i]][hr_idx],
          expected_duration = designs$dur_exp[[i]][hr_idx],
          expected_events = designs$n_exp[[i]][hr_idx]
        )
        
        analytical_results <- rbind(analytical_results, result)
      }
    }
  }
  
  return(analytical_results)
}

# Extract analytical results
analytical_results <- extract_analytical_results(designs, hr_grid)

# Combine analytical and simulation results for comparison
comparison <- full_join(
  sim_summaries %>% 
    select(design_tag, hr, power_ia, power_overall, expected_duration, expected_events),
  analytical_results,
  by = c("design_tag", "hr"),
  suffix = c("_sim", "_analytical")
)

# Calculate differences
comparison <- comparison %>%
  mutate(
    power_ia_diff = power_ia_sim - power_ia_analytical,
    power_overall_diff = power_overall_sim - power_overall_analytical,
    duration_diff = expected_duration_sim - expected_duration_analytical,
    events_diff = expected_events_sim - expected_events_analytical
  )

#-----------------------------------------------------
# 6. Create comparison plots
#-----------------------------------------------------

# Plot comparison of power at IA
plot_power_ia_comparison <- ggplot(comparison, aes(x = hr, color = design_tag)) +
  geom_line(aes(y = power_ia_analytical, linetype = "Analytical")) +
  geom_point(aes(y = power_ia_sim, shape = "Simulation"), size = 3) +
  labs(
    title = "Power at Interim Analysis: Analytical vs. Simulation",
    x = "Hazard Ratio",
    y = "Power at IA1",
    color = "Design",
    linetype = "Method",
    shape = "Method"
  ) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1))

# Plot comparison of overall power
plot_power_overall_comparison <- ggplot(comparison, aes(x = hr, color = design_tag)) +
  geom_line(aes(y = power_overall_analytical, linetype = "Analytical")) +
  geom_point(aes(y = power_overall_sim, shape = "Simulation"), size = 3) +
  labs(
    title = "Overall Power: Analytical vs. Simulation",
    x = "Hazard Ratio",
    y = "Overall Power",
    color = "Design",
    linetype = "Method",
    shape = "Method"
  ) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1))

# Plot comparison of expected duration
plot_duration_comparison <- ggplot(comparison, aes(x = hr, color = design_tag)) +
  geom_line(aes(y = expected_duration_analytical, linetype = "Analytical")) +
  geom_point(aes(y = expected_duration_sim, shape = "Simulation"), size = 3) +
  labs(
    title = "Expected Study Duration: Analytical vs. Simulation",
    x = "Hazard Ratio",
    y = paste0("Expected Duration (", time_unit, ")"),
    color = "Design",
    linetype = "Method",
    shape = "Method"
  ) +
  theme_bw()

# Plot comparison of expected events
plot_events_comparison <- ggplot(comparison, aes(x = hr, color = design_tag)) +
  geom_line(aes(y = expected_events_analytical, linetype = "Analytical")) +
  geom_point(aes(y = expected_events_sim, shape = "Simulation"), size = 3) +
  labs(
    title = "Expected Number of Events: Analytical vs. Simulation",
    x = "Hazard Ratio",
    y = "Expected Number of Events",
    color = "Design",
    linetype = "Method",
    shape = "Method"
  ) +
  theme_bw()

# Save results and plots
save(sim_results, sim_summaries, analytical_results, comparison, 
     file = "validate_GSD_simulations_results.RData")

# Create a PDF with all comparison plots
pdf("validate_GSD_simulations_plots.pdf", width = 10, height = 8)
print(plot_power_ia_comparison)
print(plot_power_overall_comparison)
print(plot_duration_comparison)
print(plot_events_comparison)
dev.off()

# Create a table with detailed comparisons
comparison_table <- comparison %>%
  select(design_tag, hr, 
         power_ia_sim, power_ia_analytical, power_ia_diff,
         power_overall_sim, power_overall_analytical, power_overall_diff,
         expected_duration_sim, expected_duration_analytical, duration_diff,
         expected_events_sim, expected_events_analytical, events_diff)

write.csv(comparison_table, "validate_GSD_simulations_comparison.csv", row.names = FALSE)

# Print a summary of the validation
cat("\n=== Validation Summary ===\n")
cat("Average absolute differences between simulation and analytical results:\n")

summary_stats <- comparison %>%
  group_by(design_tag) %>%
  summarize(
    avg_power_ia_diff = mean(abs(power_ia_diff), na.rm = TRUE),
    avg_power_overall_diff = mean(abs(power_overall_diff), na.rm = TRUE),
    avg_duration_diff = mean(abs(duration_diff), na.rm = TRUE),
    avg_events_diff = mean(abs(events_diff), na.rm = TRUE)
  )

print(summary_stats)

cat("\nValidation complete. Results saved to:\n")
cat("- validate_GSD_simulations_results.RData\n")
cat("- validate_GSD_simulations_plots.pdf\n")
cat("- validate_GSD_simulations_comparison.csv\n")