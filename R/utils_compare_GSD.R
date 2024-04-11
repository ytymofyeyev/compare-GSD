##### Extending/fixing {gsDesign} functions ####################################

# Gets the expected number of events at a particular time across the control and
# treatment arms. Essentially a wrapper around gsDesign::eEvents; see there for
# argument definitions
eEvents_total                  <- function(
    hr      = 1,
    ratio   = 1,
    lambdaC = 1,
    eta     = 0,
    gamma   = 1,
    R       = 1,
    S       = NULL,
    T,
    Tfinal  = NULL,
    minfup  = 0,
    digits  = 4,
    target  = 0) # Target number of events; used for root-solving by tEvents()
{
  if (T == 0) return(0)
  Qe  <- ratio/(1 + ratio)
  eDC <- gsDesign:::eEvents1(lambda = lambdaC,
                             eta    = eta,
                             gamma  = gamma*(1 - Qe),
                             R      = R,
                             S      = S,
                             T      = T,
                             Tfinal = Tfinal,
                             minfup = minfup)
  eDE <- gsDesign:::eEvents1(lambda = lambdaC*hr,
                             eta    = eta,
                             gamma  = gamma*Qe,
                             R      = R,
                             S      = S,
                             T      = T,
                             Tfinal = Tfinal,
                             minfup = minfup)
  sum(eDC$d + eDE$d) - target
}

# Gets the expected time when a given number of events will be reached
tEvents                        <- function(
    e,                                  # Target number of events    
    hr      = 1,
    ratio   = 1,
    lambdaC = 1,
    eta     = 0,
    gamma   = 1,
    R       = 1,
    S       = NULL,
    T,
    Tfinal  = NULL,
    minfup  = 0,
    tol     = .Machine$double.eps^0.25)
{
  # Solve eEvents_total for parameter T
  stats::uniroot(f        = eEvents_total,
                 interval = c(1e-04, 1e5),
                 hr       = hr,
                 ratio    = ratio,
                 lambdaC  = lambdaC,
                 eta      = eta,
                 gamma    = gamma,
                 R        = R,
                 S        = S,
                 Tfinal   = Tfinal,
                 minfup   = minfup,
                 target   = e,
                 tol      = tol)$root
}

# A corrected version of the function gsDesign::Power.ssrCP; see there for
# argument definitions. Computes the power of a two-stage SSR design
Power.ssrCP_corrected          <- function(
    x,
    theta = NULL,
    delta = NULL,
    r     = 18)
{
  if (!inherits(x, "ssrCP")) {
    stop("Power.ssrCP must be called with x of class ssrCP")
  }
  if (is.null(theta) & is.null(delta)) {
    theta        <- (0:80) / 40 * x$x$delta
    delta        <- (x$x$delta1 - x$x$delta0) / x$x$delta * theta + x$x$delta0
  } else if (!is.null(theta)) {
    delta        <- (x$x$delta1 - x$x$delta0) / x$x$delta * theta + x$x$delta0
  } else {
    theta        <- (delta - x$x$delta0) / (x$x$delta1 - x$x$delta0) * x$x$delta
  }
  en             <- theta
  Power          <- en
  mu             <- sqrt(x$x$n.I[1]) * theta
  Power          <- stats::pnorm(x$x$upper$bound[1] - mu, lower.tail = FALSE)
  en             <- (x$x$n.I[1] + x$overrun) *
    (Power + stats::pnorm(x$x$lower$bound[1] - mu))
  cpmin          <- gsDesign::condPower(z1    = x$x$lower$bound[1],
                                        n2    = x$x$n.I[2] - x$x$n.I[1],
                                        z2    = x$z2fn,
                                        theta = x$theta,
                                        x     = x$x,
                                        beta  = x$beta)
  cpmax          <- gsDesign::condPower(z1    = x$x$upper$bound[1],
                                        n2    = x$x$n.I[2] - x$x$n.I[1],
                                        z2    = x$z2fn,
                                        theta = x$theta,
                                        x     = x$x,
                                        beta  = x$beta)
  if (cpmax <= x$cpadj[1] || cpmin >= x$cpadj[2]) {
    en           <- en + x$x$n.I[2] * (stats::pnorm(x$x$upper$bound[1] - mu) -
                                         stats::pnorm(x$x$lower$bound[1] - mu))
    a            <- x$x$lower$bound[1]
    b            <- x$x$upper$bound[2]
    n2           <- x$x$n.I[2] - x$x$n.I[1]
    grid         <- gsDesign::normalGrid(mu     = (a + b) / 2,
                                         bounds = c(a, b),
                                         r      = r)
    for (i in 1:length(theta)) {
      Power[i]   <- Power[i] + sum(
        stats::dnorm(grid$z - sqrt(x$x$n.I[1]) * theta[i]) * grid$gridwgts *
          stats::pnorm(x$z2fn(grid$z, x = x$x, n2 = n2) - theta[i] * sqrt(n2),
                       lower.tail = FALSE)
      )
    }
    return(data.frame(theta = theta,
                      delta = delta,
                      Power = Power,
                      en    = en))
  }
  if (cpmin < x$cpadj[1]) {
    changepoint  <- stats::uniroot(gsDesign:::condPowerDiff,
                                   interval = c(x$x$lower$bound[1],
                                                x$x$upper$bound[1]),
                                   target   = x$cpadj[1],
                                   n2       = x$x$n.I[2] - x$x$n.I[1],
                                   z2       = x$z2fn,
                                   theta    = x$theta,
                                   x        = x$x)$root
    en           <- en + x$x$n.I[2] * (stats::pnorm(changepoint - mu) -
                                         stats::pnorm(x$x$lower$bound[1] - mu))
    a            <- x$x$lower$bound[1]
    b            <- changepoint
    n2           <- x$x$n.I[2] - x$x$n.I[1]
    grid         <- gsDesign::normalGrid(mu     = (a + b) / 2,
                                         bounds = c(a, b),
                                         r      = r)
    for (i in 1:length(theta)) {
      Power[i]   <- Power[i] + sum(
        stats::dnorm(grid$z - sqrt(x$x$n.I[1]) * theta[i]) * grid$gridwgts *
          stats::pnorm(x$z2fn(grid$z, x = x$x, n2 = n2) - theta[i] * sqrt(n2),
                       lower.tail = FALSE)
      )
    }
  } else {
    changepoint  <- x$x$lower$bound[1]
  }
  if (cpmax > x$cpadj[2]) {
    changepoint2 <- stats::uniroot(gsDesign:::condPowerDiff,
                                   interval = c(changepoint,
                                                x$x$upper$bound[1]),
                                   target   = x$cpadj[2],
                                   n2       = x$x$n.I[2] - x$x$n.I[1],
                                   z2       = x$z2fn,
                                   theta    = x$theta,
                                   x        = x$x)$root
    en           <- en + x$x$n.I[2] * (stats::pnorm(x$x$upper$bound[1] - mu) -
                                         stats::pnorm(changepoint2 - mu))
    a            <- changepoint2
    b            <- x$x$upper$bound[1]
    n2           <- x$x$n.I[2] - x$x$n.I[1]
    grid         <- gsDesign::normalGrid(mu     = (a + b) / 2,
                                         bounds = c(a, b),
                                         r      = r)
    for (i in 1:length(theta)) {
      Power[i]   <- Power[i] + sum(
        stats::dnorm(grid$z - sqrt(x$x$n.I[1]) * theta[i]) * grid$gridwgts *
          stats::pnorm(x$z2fn(grid$z, x = x$x, n2 = n2) - theta[i] * sqrt(n2),
                       lower.tail = FALSE)
      )
    }
  } else {
    #changepoint2 <- x$upper$bound[1]
    changepoint2 <- x$x$upper$bound[1] # YT corrected
  }
  if (gsDesign:::n2sizediff(z1     = changepoint,
                            target = x$maxinc * x$x$n.I[2],
                            beta   = x$beta,
                            z2     = x$z2fn,
                            theta  = x$theta,
                            x      = x$x) > 0) {
    changepoint3 <- stats::uniroot(gsDesign:::condPowerDiff,
                                   interval = c(changepoint, 15),
                                   target   = 1 - x$beta,
                                   x        = x$x,
                                   n2       =
                                     x$maxinc * x$x$n.I[2] - x$x$n.I[1],
                                   z2       = x$z2fn,
                                   theta    = x$theta)$root
    if (changepoint3 >= changepoint2) {
      en         <- en + x$maxinc * x$x$n.I[2] *
        (stats::pnorm(changepoint2 - mu) - stats::pnorm(changepoint - mu))
      a          <- changepoint
      b          <- changepoint2
      n2         <- x$maxinc * x$x$n.I[2] - x$x$n.I[1]
      grid       <- gsDesign::normalGrid(mu     = (a + b) / 2,
                                         bounds = c(a, b),
                                         r      = r)
      for (i in 1:length(theta)) {
        Power[i] <- Power[i] + sum(
          stats::dnorm(grid$z - sqrt(x$x$n.I[1]) * theta[i]) * grid$gridwgts *
            stats::pnorm(x$z2fn(grid$z, x = x$x, n2 = n2) - theta[i] * sqrt(n2),
                         lower.tail = FALSE)
        )
      }
      return(data.frame(theta = theta,
                        delta = delta,
                        Power = Power,
                        en    = en))
    }
    en           <- en + x$maxinc * x$x$n.I[2] *
      (stats::pnorm(changepoint3 - mu) - stats::pnorm(changepoint - mu))
    a            <- changepoint
    b            <- changepoint3
    n2           <- x$maxinc * x$x$n.I[2] - x$x$n.I[1]
    grid         <- gsDesign::normalGrid(mu     = (a + b) / 2,
                                         bounds = c(a, b),
                                         r      = r)
    for (i in 1:length(theta)) {
      Power[i]   <- Power[i] + sum(
        stats::dnorm(grid$z - sqrt(x$x$n.I[1]) * theta[i]) * grid$gridwgts *
          stats::pnorm(x$z2fn(grid$z, x = x$x, n2 = n2) - theta[i] * sqrt(n2),
                       lower.tail = FALSE)
      )
    }
  } else {
    changepoint3 <- changepoint
  }
  grid           <- gsDesign::normalGrid(
    mu     = (changepoint3 + changepoint2) / 2,
    bounds = c(changepoint3, changepoint2),
    r      = r
  )
  y              <- gsDesign::ssrCP(z1     = grid$z,
                                    theta  = x$theta,
                                    maxinc = x$maxinc * 2,
                                    beta   = x$beta,
                                    x      = x$x,
                                    cpadj  = c(0.05, 0.9999),
                                    z2     = x$z2fn)$dat
  for (i in 1:length(theta)) {
    grid$density <- stats::dnorm(y$z1 - sqrt(x$x$n.I[1]) *
                                   theta[i])
    Power[i]     <- Power[i] + sum(
      grid$density * grid$gridwgts *
        stats::pnorm(y$z2 - theta[i] * sqrt(y$n2 - x$x$n.I[1]),
                     lower.tail = FALSE)
    )
    en[i]        <- en[i] + sum(grid$density * grid$gridwgts * y$n2)
  }
  return(data.frame(theta = theta,
                    delta = delta,
                    Power = Power,
                    en    = en))
}

##### Minor utility functions ##################################################

# Assembles a verbal description of a design
abbreviate_sf                  <- function(
    b,          # Upper or lower from gsDesign object
    digits = 2) # Number of digits to print to
{
  sfName       <- b$sf(1, 1, b$param)$name |>
    switch("Lan-DeMets O'Brien-Fleming approximation" = "LDOF",
           b$sf(1, 1, b$param)$name)
  if (!is.null(b$param) & length(b$param) <= 4) {
    sfParam    <- paste(
      lapply(b$param, function(y) { if (is.numeric(y)) round(y, digits) }),
      collapse = " "
    )
    sfParamStr <- paste0(" (", sfParam, ")")
  } else {
    sfParamStr <- ""
  }
  paste0(sfName, sfParamStr)
}

#Aux function to call nSurv with the 'minfup' = NULL
call_nSurv <- function(
    ...)
{
  do.call(gsDesign::nSurv, args = list(..., minfup = NULL))
}

# Aux function to call gsSurv() or set a dummy {gsDesign} object
call_gsSurv                    <- function(
    ...)
{
  params <- list(... , minfup = NULL)
  if (params$test.type %in% 1:6) {
    do.call(gsDesign::gsSurv, args = params)
  } else {
    return(params)
  }
}

exec_calc                      <- function(
    designs,
    alpha,
    hr_grid)
{
  designs |>
    dplyr::mutate(R         = purrr::map_dbl(design_spec, \(x) x$R),
                  gamma     = purrr::map_dbl(design_spec, \(x) x$gamma),
                  test.type = purrr::map_dbl(design_spec, \(x) x$test.type),
                  beta      = purrr::map_dbl(design_spec, \(x) x$beta),
                  ratio     = purrr::map_dbl(design_spec, \(x) x$ratio),
                  timing    = purrr::map(design_spec, \(x) x$timing),
                  sfu       = purrr::map(design_spec, \(x) x$sfu),
                  sfupar    = purrr::map(design_spec, \(x) x$sfupar),
                  sfl       = purrr::map(design_spec, \(x) x$sfl),
                  sflpar    = purrr::map(design_spec, \(x) x$sflpar)
    ) |>
    dplyr::mutate(fixedd = purrr::pmap(.l = list(R       = R,
                                                 gamma   = gamma,
                                                 eta     = eta,
                                                 lambdaC = lambdaC,
                                                 hr      = hr,
                                                 alpha   = alpha,
                                                 beta    = beta,
                                                 ratio   = ratio),
                                       .f = call_nSurv),
                  k      = purrr::map_int(timing,
                                          \(x) ifelse(is.null(x), 2,
                                                      length(x))),
                  gsd    = purrr::pmap(.l = list(k         = k,
                                                 test.type = test.type,
                                                 R         = R,
                                                 gamma     = gamma,
                                                 eta       = eta,
                                                 lambdaC   = lambdaC,
                                                 hr        = hr,
                                                 alpha     = alpha,
                                                 beta      = beta,
                                                 ratio     = ratio,
                                                 timing    = timing,
                                                 sfu       = sfu,
                                                 sfupar    = sfupar,
                                                 sfl       = sfl,
                                                 sflpar    = sflpar),
                                       .f = call_gsSurv),
                  OC     = purrr::pmap(list(design_spec, gsd), get_design_OC,
                                       hr_grid = hr_grid))
}

##### Operating characteristic functions #######################################

# S3 class to generate operating characteristics
get_design_OC                  <- function(
    spec,
    gsd,
    ...)
{
  UseMethod("get_design_OC")
}

# Gets the operating characteristics of a fixed GSD
get_design_OC.gsd_fixed        <- function(
    design_spec,
    gsd,
    hr_grid)
{
  ratio      <- gsd$ratio
  Q          <- ratio/(1 + ratio)  
  gsd_OC     <- gsDesign::gsProbability(d = gsd, theta = -Q*log(hr_grid))
  # Compile EN, power and total study duration
  tibble::tibble(
    Design     = paste0(
      "'", design_spec$tag, ": '*", "list(",
      #            " upper: ", abbreviate_sf(gsd$upper),
      #            " lower: "    , abbreviate_sf(gsd$lower),
      #             ", ni's : ", paste0( round(gsd$n.I), collapse = ", ") 
      paste0("italic(d)[", 1:gsd$k, "] == ", round(gsd$n.I), collapse = ", "),
      ")"
      #           ", Nrand = ", round(sum(design_spec$R*design_spec$gamma))
    ),      
    theta      = gsd_OC$theta,
    HR         = exp(-theta/Q),
    en         = gsd_OC$en,
    Power      = apply(gsd_OC$upper$prob, 2, sum),
    Power_1IA  = gsd_OC$upper$prob[1, ], 
    Duration   = purrr::map2_dbl(en, HR, tEvents,
                                 ratio   = ratio, 
                                 lambdaC = gsd$lambdaC,
                                 eta     = gsd$etaC,
                                 gamma   = gsd$gamma,
                                 R       = gsd$R),
    Time_IA1   = purrr::map2_dbl(gsd$n.I[1], HR, tEvents,
                                 ratio   = ratio,
                                 lambdaC = gsd$lambdaC,
                                 eta     = gsd$etaC,
                                 gamma   = gsd$gamma,
                                 R       = gsd$R),
    Max_Time   = purrr::map2_dbl(gsd$n.I[gsd$k], HR, tEvents,
                                 ratio   = ratio,
                                 lambdaC = gsd$lambdaC,
                                 eta     = gsd$etaC,
                                 gamma   = gsd$gamma,
                                 R       = gsd$R)
  ) |> 
    dplyr::mutate( 
      z1          = -Q*log(hr_grid)*sqrt(gsd$n.I[1]),
      zf1         = ifelse(!is.null(gsd$lower), gsd$lower$bound[1], -Inf),
      n2aux       = dplyr::if_else(gsd$k == 2, gsd$n.I[2], NA),
      n2          = dplyr::if_else(gsd$upper$bound[1] < z1 | z1 < zf1,
                                   gsd$n.I[1], n2aux),
      eHR1        = exp(-1/Q * gsd$upper$bound[1]/sqrt(gsd$n.I[1])),  
      fHR1        = exp(-1/Q * zf1/sqrt(gsd$n.I[1])),
      eZ1         = gsd$upper$bound[1],
      fZ1         = zf1,
      zone        = dplyr::case_when(gsd$upper$bound[1] < z1 ~ "eff",
                                     z1 < zf1 ~ "fut",
                                     z1 <= gsd$upper$bound[1] &
                                       z1 >= zf1 ~ "cont"),
      design_type = "Fixed"
    )
}

# Gets the operating characteristics of an adaptive 2-stage design with event
# re-estimation, based on Mehta and Pocock's promising zone approach
get_design_OC.adaptive_mp      <- function(
    design_spec,
    gsd,
    hr_grid)
{
  ratio   <- gsd$ratio
  Q       <- ratio/(1 + ratio)  
  gsd_OC  <- gsDesign::gsProbability(d = gsd, theta = -Q*log(hr_grid)) 
  zIA     <- -Q*log(hr_grid)*sqrt(gsd$n.I[1])
  thetaCP <- if (!is.null(design_spec$hrCP)) -Q*log(design_spec$hrCP) else NULL
  x       <- gsDesign::ssrCP(x      = gsd,
                             z1     = zIA,
                             theta  = thetaCP,
                             beta   = design_spec$betastar,
                             cpadj  = design_spec$cpadj,
                             maxinc = design_spec$maxinflation,
                             z2     = design_spec$z2)
  powObj  <- Power.ssrCP_corrected(x, theta = -Q*log(hr_grid))
  tibble::tibble(Design = paste0(
    "'", design_spec$tag,  ": '*", "list(",
    #            " upper: ", abbreviate_sf(gsd$upper),
    #            " lower: "    , abbreviate_sf(gsd$lower),
    #            " CP at ", ifelse(is.null(design_spec$hrCP),"trend",
    #                              design_spec$hrCP),
    "italic(d)[1] == ", round(gsd$n.I[1]),
    ", italic(d)[2] == (",
    round(min(x$dat$n2[x$dat$n2 > gsd$n.I[1]])),
    " - ",
    round(max(x$dat$n2)), "))"
    #            ", Nrand = ", round(sum(design_spec$R*design_spec$gamma))
  ), 
  cbind(powObj, design_type = "Adaptive")
  ) |> dplyr::mutate(
    HR        = exp(-theta/Q),
    Power_1IA = gsd_OC$upper$prob[1, ], 
    # Compute expected duration corresponding to number of events
    Duration  = purrr::map2_dbl(en, HR, tEvents,
                                ratio   = ratio,
                                lambdaC = gsd$lambdaC, 
                                eta     = gsd$etaC,
                                gamma   = gsd$gamma,
                                R       = gsd$R),
    Time_IA1  = purrr::map2_dbl(gsd$n.I[1], HR, tEvents,
                                ratio   = ratio, 
                                lambdaC = gsd$lambdaC,
                                eta     = gsd$etaC,
                                gamma   = gsd$gamma,
                                R       = gsd$R),
    Max_Time  = purrr::map2_dbl(gsd$n.I[2]*design_spec$maxinflation, HR, tEvents,
                                ratio   = ratio, 
                                lambdaC = gsd$lambdaC,
                                eta     = gsd$etaC,
                                gamma   = gsd$gamma,
                                R       = gsd$R),
    n2        = x$dat$n2,
    z1        = x$dat$z1,
    eHR1      = exp(-1/Q * gsd$upper$bound[1]/sqrt(gsd$n.I[1])),
    zf1       = ifelse(!is.null(gsd$lower), gsd$lower$bound[1], -Inf),
    fHR1      = exp(-1/Q * zf1/sqrt(gsd$n.I[1])),
    eZ1       = gsd$upper$bound[1],
    fZ1       = zf1,
    zone      = dplyr::case_when(gsd$upper$bound[1] < z1 ~ "eff",
                                 z1 < zf1 ~ "fut",
                                 z1 <= gsd$upper$bound[1] & z1 >= zf1 ~ "cont")
  )
}

# Gets the operating characteristics of an adaptive 2-stage design with event
# re-estimation, based on the {adoptr} optimal design approach
get_design_OC.adaptive_optimal <- function(
    design_spec,
    gsd,
    hr_grid)
{
  alpha              <- gsd$alpha    
  beta               <- gsd$beta
  ratio              <- gsd$ratio
  Q                  <- ratio/(1 + ratio)
  theta_alt          <- -log(gsd$hr)*sqrt(ratio)/(1 + ratio)
  datadist           <- adoptr::Normal(two_armed = FALSE)
  null               <- adoptr::PointMassPrior(theta = 0,         mass = 1)
  alternative        <- adoptr::PointMassPrior(theta = theta_alt, mass = 1)
  power              <- adoptr::Power(dist = datadist, prior = alternative)
  toer               <- adoptr::Power(dist = datadist, prior = null)
  mss                <- adoptr::MaximumSampleSize()
  # Set function to minimize
  ess                <- adoptr::ExpectedSampleSize(dist  = datadist,
                                                   prior = alternative)
  cp                 <- adoptr::ConditionalPower(dist  = datadist,
                                                 prior = alternative)
  initial_design     <- adoptr::get_initial_design(theta = theta_alt,
                                                   alpha = alpha,
                                                   beta  = beta,
                                                   type  = "two-stage",
                                                   dist  = datadist,
                                                   order = design_spec$order)
  initial_design@n1  <- design_spec$n1_min
  initial_design@c1f <- design_spec$zFutil
  initial_design@c1e <- design_spec$zEffic
  initial_design     <- adoptr::make_fixed(initial_design, n1, c1f, c1e)
  # Call optimizer
  opt1               <- adoptr::minimize(
    ess,   # Objective function: expected number of events
    adoptr::subject_to(power >= 1 - beta, # power requirement
                       toer  <= alpha),   # type I control
    initial_design,
    opts = list(algorithm = "NLOPT_LN_COBYLA",
                xtol_rel  = 1e-05,
                maxeval   = 10000)
  )
  summary(opt1$design, "Power" = power, "ESS" = ess, "CP" = cp)
  tibble::tibble(
    Design   = paste0(
      "'", design_spec$tag,  ": Fixed '*", "list(",
      "italic(d)[1] == ", round(opt1$design@n1),
      ", italic(z)[IA1~efficacy] == ", round(opt1$design@c1e, 2),
      ", italic(z)[IA1~futility] == ", round(opt1$design@c1f, 2),
      ")"
    ), 
    theta    = -Q*log(hr_grid),
    HR       = exp(-theta/Q),
    z1       = theta*sqrt(opt1$design@n1),
    zone     = dplyr::case_when(design_spec$zEffic < z1 ~ "eff",
                                z1 < design_spec$zFutil ~ "fut",
                                z1 <= design_spec$zEffic &
                                  z1 >= design_spec$zFutil ~ "cont")
  ) |> dplyr::mutate(
    n2        = opt1$design@n1 +
      purrr::map_dbl(z1, \(z1) adoptr::n2(opt1$design, x = z1)),
    n2_max    = max(n2),
    eHR1      = exp(-opt1$design@c1e/sqrt(opt1$design@n1)/Q),
    fHR1      = exp(-opt1$design@c1f/sqrt(opt1$design@n1)/Q),
    eZ1       = opt1$design@c1e,
    fZ1       = opt1$design@c1f,
    en        = purrr::map_dbl(
      theta,
      function(d) {
        adoptr::evaluate(adoptr::ExpectedSampleSize(
          datadist, adoptr::PointMassPrior(d, 1)), opt1$design
        )
      }
    ),
    Duration  = purrr::map2_dbl(en, HR, tEvents,
                                ratio   = ratio,
                                lambdaC = gsd$lambdaC, 
                                eta     = gsd$eta,
                                gamma   = gsd$gamma,
                                R       = gsd$R),
    Time_IA1  = purrr::map2_dbl(opt1$design@n1, HR, tEvents,
                                ratio   = ratio, 
                                lambdaC = gsd$lambdaC,
                                eta     = gsd$eta,
                                gamma   = gsd$gamma,
                                R       = gsd$R),
    Max_Time  = purrr::map2_dbl(n2_max, HR, tEvents,
                                ratio   = ratio, 
                                lambdaC = gsd$lambdaC,
                                eta     = gsd$eta,
                                gamma   = gsd$gamma,
                                R       = gsd$R),
    Power_1IA = purrr::map_dbl(
      theta, \(x) {1 - stats::pnorm(q    = opt1$design@c1e,
                                    mean = x*sqrt(opt1$design@n1)) }
    )
  )
}

# Gets the operating characteristics of an adaptive 2-stage design with event
# re-estimation, based on the design approach where the event re-estimation rule
# is defined explicitly by a piecewise-linear function (including step function)
get_design_OC.adaptive_pwl <- function(
    design_spec,
    gsd,
    hr_grid)
{
  ratio       <- gsd$ratio
  Q           <- ratio / (1 + ratio)
  theta_grid  <- -Q*log(hr_grid)
  gsD_OC      <-
    gsDesign::gsProbability(d = gsd, theta = theta_grid)
  I1  <- gsd$n.I[1]
  f1  <- gsd$lower$bound[1]
  e1  <- gsd$upper$bound[1]
  n2max <- gsd$n.I[2] * design_spec$maxinflation - gsd$n.I[1]
  
  # add the ends of interval for 'z_mesh' and 'I_mesh"
  # also add dummy value to capture early stopping for futility
  # z- and I-mesh has to take into account  'f' and 'rule' parameters in 'approxfun()'
  z_mesh <-
    c(-.Machine$double.eps ^ 0.5, design_spec$z1_mesh_in_f1_e1, 1)
  I_mesh <-
    c(0, design_spec$I_ratio_to_n2max, 0)
  obj <- list(z1_nodes = f1 + (e1 - f1) * z_mesh,
              I_values = I1 + I_mesh * n2max,
              f        = 0,
              rule     = 2,
              method   = design_spec$approx_method)
  xw    <- gauss_legendre_xw(design_spec$order)
  I2_tilde_approx <-
    approxfun(
      x      = obj$z1_nodes,
      y      = obj$I_values,
      f      = obj$f,
      rule   = obj$rule,
      method = obj$method
    )
  z1Grid <- theta_grid * sqrt(I1)
  n2     <- I2_tilde_approx(z1Grid)
  
  en <- ei_2_Vec(
    theta         = theta_grid,
    f1            = f1,
    e1            = e1,
    I1            = I1,
    x             = xw$x,
    w             = xw$w,
    I2_constraint = "piecewise-linear",
    I2_tilde      = get_I_at_legendre_x(obj, xw$x)
  )
  Power <- power_2_Vec(
    theta            = theta_grid,
    f1               = f1,
    e1               = e1,
    I1               = I1,
    x                = xw$x,
    w                = xw$w,
    e2_constraint    = "linear",
    # use the combination test theory
    e2_tilde         = c(gsd$upper$bound[2] * sqrt(gsd$n.I[2] / (gsd$n.I[2] - gsd$n.I[1])),
                         (-1) * sqrt(gsd$n.I[1] / (gsd$n.I[2] - gsd$n.I[1]))),
    I2_constraint    = "piecewise-linear",
    I2_tilde         = get_I_at_legendre_x(obj, xw$x),
    binding_futility = FALSE
  )
  gsD_OC_dat <- tibble::tibble(
    Design = paste0(
      "'", design_spec$tag,  ": '*", "list(",
      #            " upper: ", abbreviate_sf(gsd$upper),
      #            " lower: "    , abbreviate_sf(gsd$lower),
      #            " CP at ", ifelse(is.null(design_spec$hrCP),"trend",
      #                              design_spec$hrCP),
      "italic(d)[1] == ", round(gsd$n.I[1]),
      ", italic(d)[2] == (", 
      round(min(n2[n2 > gsd$n.I[1]])),
      " - ",
      round(max(n2)), "))"
      #            ", Nrand = ", round(sum(design_spec$R*design_spec$gamma))
    ),   
    theta     = theta_grid,
    HR        = exp(-1 / Q * theta),
    en        = en,
    Power     = Power,
    Power_1IA = gsD_OC$upper$prob[1,],
    Duration  = purrr::map2_dbl( en, HR, tEvents, 
                                 ratio = ratio,
                                 lambdaC = gsd$lambdaC,
                                 eta = gsd$etaC,
                                 gamma = gsd$gamma,
                                 R = gsd$R),
    Time_IA1 = purrr::map2_dbl( I1, HR, tEvents,
                                ratio = ratio,
                                lambdaC = gsd$lambdaC,
                                eta = gsd$etaC,
                                gamma = gsd$gamma,
                                R = gsd$R),
    Max_Time = purrr::map2_dbl( max(obj$I_values), HR, tEvents,
                                ratio = ratio,
                                lambdaC = gsd$lambdaC,
                                eta = gsd$etaC,
                                gamma = gsd$gamma,
                                R = gsd$R)
  ) |>
    mutate(
      z1         = theta_grid * sqrt(I1),
      zf1        = f1,
      zone      = dplyr::case_when(gsd$upper$bound[1] < z1 ~ "eff",
                                   z1 < zf1 ~ "fut",
                                   z1 <= e1 & z1 >= zf1 ~ "cont"),
      n2         = n2,
      eHR1       = exp(-1 / Q * e1 / sqrt(I1)),
      fHR1       = exp(-1 / Q * f1 / sqrt(I1)),
      eZ1        = e1,
      fZ1        = f1,
      design_type = "Adaptive"
    )
  return(gsD_OC_dat)
}

##### Plotting functions #######################################################

plot_drop_out                  <- function(
    eta,
    time_unit)
{
  df <- tibble::tibble(time          = seq(1e-6, 24, length.out = 500),
                       drop_out_rate = 1 - exp(-time*eta))
  ggplot2::ggplot(df, ggplot2::aes(time, drop_out_rate)) +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = paste0("Time (", time_unit, ")"),
                  y = "Drop-out rate (%)")
}

plot_duration_by_n_and_e       <- function(
    designs,
    primary_outcome,
    time_unit)
{
  df                 <- tibble::tibble(
    ratio   = 1,
    n       = round(designs$spec[[1]]$R*designs$spec[[1]]$gamma*
                      seq(0.75, 1.25, length.out = 500)),
    eta     = primary_outcome$eta,
    lambdaC = primary_outcome$lambdaC,
    hr      = primary_outcome$hr
  ) |>
    tidyr::expand_grid(R = designs$spec[[1]]$R,
                       e = round(designs$gsd[[1]]$n.I)) |>
    dplyr::mutate(gamma = n/R) |>
    dplyr::mutate(time = purrr::pmap_dbl(.l = list(e, hr, ratio, lambdaC, eta,
                                                   gamma, R),
                                         .f = tEvents),
                  e    = factor(e))
  ggplot2::ggplot(df, ggplot2::aes(x = n, y = time, col = e)) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(x   = expression(paste("Number of subjects, ", italic(n))),
                  y   = paste0("Expected time to event target (", time_unit,
                               ")"),
                  col = "Events")
}

plot_duration_expected         <- function(
    designs,
    primary_outcome,
    time_unit)
{
  ggplot2::ggplot(dplyr::bind_rows(designs$OC),
                  ggplot2::aes(x = HR, y = Duration, col = Design)) +
    ggplot2::geom_line() +
    ggplot2::scale_colour_discrete(labels = scales::parse_format()) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::geom_vline(xintercept = primary_outcome$hr, lty = "dotted") +
    ggplot2::labs(x = expression(paste("Hazard ratio, ", italic(HR))),
                  y = paste0("Expected duration (", time_unit, ")"))
}

plot_duration_ia1              <- function(
    designs,
    primary_outcome,
    time_unit) {
  ggplot2::ggplot(dplyr::bind_rows(designs$OC),
                  ggplot2::aes(x = HR, y = Time_IA1, col = Design)) +
    ggplot2::geom_line() +
    ggplot2::scale_colour_discrete(labels = scales::parse_format()) +
    ggplot2::labs(x = expression(paste("Hazard ratio, ", italic(HR))),
                  y = paste0("Expected duration to IA1 (", time_unit, ")")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::geom_vline(xintercept = primary_outcome$hr, lty = "dotted")
}

plot_duration_fa               <- function(
    designs,
    primary_outcome,
    time_unit)
{
  ggplot2::ggplot(dplyr::bind_rows(designs$OC),
                  ggplot2::aes(x = HR, y = Max_Time, col = Design)) +
    ggplot2::geom_line() +
    ggplot2::scale_colour_discrete(labels = scales::parse_format()) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::geom_vline(xintercept = primary_outcome$hr, lty = "dotted") +
    ggplot2::labs(x = expression(paste("Hazard ratio, ", italic(HR))),
                  y = paste0("Expected duration to FA (", time_unit, ")"))
}

plot_events_by_obs_HR          <- function(
    designs)
{
  ggplot2::ggplot(dplyr::bind_rows(designs$OC), 
                  ggplot2::aes(x = HR, y = n2, col = Design, by = zone)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(ggplot2::aes(xintercept = eHR1, col = Design),
                        lty = "dashed") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = fHR1, col = Design),
                        lty = "dashed") +
    ggplot2::scale_colour_discrete(labels = scales::parse_format()) +
    ggplot2::labs(x = expression(paste("Hazard ratio observed at IA1, ",
                                       widehat(italic(HR))[1])),
                  y = "Stage 2 event target") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 1))
}

plot_events_by_obs_z          <- function(
    designs)
{
  ggplot2::ggplot(dplyr::bind_rows(designs$OC), 
                  ggplot2::aes(x = z1, y = n2, col = Design, by = zone)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(ggplot2::aes(xintercept = eZ1, col = Design),
                        lty = "dashed") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = fZ1, col = Design),
                        lty = "dashed") +
    ggplot2::scale_colour_discrete(labels = scales::parse_format()) +
    ggplot2::labs(x = expression(paste("Z-value observed at IA1, ",
                                       italic(z))[1]),
                  y = "Stage 2 event target") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 1))
}

plot_events_expected           <- function(
    designs,
    primary_outcome)
{
  ggplot2::ggplot(dplyr::bind_rows(designs$OC),
                  ggplot2::aes(x = HR, y = en, col = Design)) +
    ggplot2::geom_line() +
    ggplot2::scale_colour_discrete(labels = scales::parse_format()) +
    ggplot2::labs(x = expression(paste("Hazard ratio, ", italic(HR))),
                  y = "Expected events") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::geom_vline(xintercept = primary_outcome$hr, lty = "dotted") + 
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 1))
}

plot_power                     <- function(
    designs,
    primary_outcome)
{
  ggplot2::ggplot(dplyr::bind_rows(designs$OC),
                  ggplot2::aes(x = HR, y = Power, col = Design)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::scale_colour_discrete(labels = scales::parse_format()) +
    ggplot2::labs(x = expression(paste("Hazard ratio, ", italic(HR))),
                  y = "Power (%)") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 1)) + 
    ggplot2::geom_vline(xintercept = primary_outcome$hr, lty = "dotted")
}

plot_power_ia1                 <- function(
    designs,
    primary_outcome)
{
  ggplot2::ggplot(dplyr::bind_rows(designs$OC),
                  ggplot2::aes(x = HR, y = Power_1IA, col = Design)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::scale_colour_discrete(labels = scales::parse_format()) +
    ggplot2::labs(x = expression(paste("Hazard ratio, ", italic(HR))),
                  y = "Power at IA1 (%)") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 1)) + 
    ggplot2::geom_vline(xintercept = primary_outcome$hr, lty = "dotted")
}
