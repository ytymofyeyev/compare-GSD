# Modification of Michael Grayling's code 'adoptr_basic_v5.R' to
# implement piecewise linear specification of n2(z1) rule

require(tidyverse)

# Gets pivots and weights of Gauss-Legendre integration
gauss_legendre_xw  <- function(order) {
  j                                            <- seq_len(order - 1)
  b                                            <- j/(4*j^2 - 1)^0.5
  A                                            <- numeric(order^2)
  A[c((order + 1)*(j - 1) + 2, (order + 1)*j)] <- b
  dim(A)                                       <- c(order, order)
  sd                                           <- eigen(A, symmetric = TRUE)
  data.frame(x = rev(sd$values),
             w = 2*rev(as.vector(sd$vectors[1, ]))^2)
}

# Function to return the power of a given two-stage design for a particular
# effect theta
power_2            <- function(theta, f1, e1, I1, x, w, e2_constraint, e2_tilde,
                               I2_constraint, I2_tilde, binding_futility) {
  sqrt_I1                            <- sqrt(I1)
  E1                                 <- pnorm(theta*sqrt_I1 - e1)
  a                                  <- (e1 - f1)/2
  b                                  <- a + f1
  z1                                 <- a*x + b
  if (binding_futility) {
    if (e2_constraint == "linear" || e2_constraint == "piecewise-linear" ) {
      e2_tilde_int                   <- e2_tilde[1] + e2_tilde[2]*z1
    } else if (e2_constraint == "none") {
      e2_tilde_int                   <- e2_tilde
    }
    if (I2_constraint == "linear") {
      I2_tilde_int                   <- I2_tilde[1] + I2_tilde[2]*z1
    } else if (I2_constraint == "piecewise-linear") {
      # mapping  z1  = a*x + b 
      a                             <- (e1 - f1)/2
      b                             <- a + f1
      z1_pivots                     <- a*x + b
      I2_tilde_approx               <- approxfun(z1_pivots, I2_tilde,
                                                 method = "linear")
      I2_tilde_int                  <- I2_tilde_approx(z1)
    } else if (I2_constraint == "none") {
      I2_tilde_int                  <- I2_tilde
    }
    E2                               <-
      a*sum(w*dnorm(z1 - theta*sqrt_I1)*
              pnorm(theta*sqrt(I2_tilde_int) - e2_tilde_int))
  } else {
    if  (I2_constraint == "none") {
      I2_tilde_spline                <- splinefun(z1, I2_tilde,
                                                  method = "monoH.FC")
    } else if (I2_constraint == "piecewise-linear"){
      I2_tilde_approx                <- approxfun(z1, I2_tilde,
                                                  method = "linear", rule = 2)
    }
    s                                <- x/2 + 0.5
    z1                               <- e1 - (1 - s)/s
    e2_tilde_int                     <- e2_tilde[1] + e2_tilde[2]*z1
    if (I2_constraint == "linear") {
      I2_tilde_int                   <- I2_tilde[1] + I2_tilde[2]*z1
      I2_tilde_int[z1 < f1]          <- I2_tilde[1] + I2_tilde[2]*f1
    } else if(I2_constraint == "piecewise-linear") {
      # I2_tilde_f1                    <- I2_tilde_spline(f1)
      I2_tilde_f1                    <- I2_tilde_approx(f1)
      I2_tilde_int                   <-
        sapply(seq_along(z1),
               function(i) {
                 ifelse(z1[i] < f1, I2_tilde_f1, I2_tilde_approx(z1[i]))
               })
      I2_tilde_int[I2_tilde_int < 0] <- 0
    } else if (I2_constraint == "none") {
      I2_tilde_f1                    <- I2_tilde_spline(f1)
      I2_tilde_int                   <-
        sapply(seq_along(z1),
               function(i) {
                 ifelse(z1[i] < f1, I2_tilde_f1, I2_tilde_spline(z1[i]))
               })
      I2_tilde_int[I2_tilde_int < 0] <- 0
    }
    E2                               <-
      0.5*sum(w*dnorm(z1 - theta*sqrt_I1)*
                pnorm(theta*sqrt(I2_tilde_int) - e2_tilde_int)/s^2)
  }
  E1 + E2
}
power_2_Vec <- Vectorize(power_2, "theta")


# Function to return the Expected Information of a given two-stage design for a
# particular effect theta (treating f1 as binding)
ei_2               <- function(theta, f1, e1, I1, x, w, I2_constraint,
                               I2_tilde) {
  # I2_tilde = defines I at x - values, so need to provide corrsponding I(z1) 
  # i.e., need to map x to z1 accordingly 
  a              <- (e1 - f1)/2
  z1             <- a*x + (a + f1)
  if (I2_constraint == "linear") {
    I2_tilde_int <- I2_tilde[1] + I2_tilde[2]*z1
  } else if (I2_constraint == "piecewise-linear") {
    I2_tilde_int <- I2_tilde
  } else if (I2_constraint == "none") {
    I2_tilde_int <- I2_tilde
  }
  I1 + a*sum(w*stats::dnorm(z1 - theta*sqrt(I1))*I2_tilde_int)
}

ei_2_Vec <- Vectorize(ei_2, "theta")

# aux function that extract from 'obj' I to be evaluated at a given Legendre x
get_I_at_legendre_x <- function(obj, x, ...) {
  UseMethod("get_I_at_legendre_x")
}

get_I_at_legendre_x.default <- function(obj, x) {
  f1        <- head(obj$z1_nodes, 1)
  e1        <- tail(obj$z1_nodes, 1)
  n1        <- head(obj$I_values, 1)
  I2_tilde_approx <- approxfun(x      = c(obj$z1_nodes,Inf),
                               y      = c(obj$I_values,n1) - n1,
                               f      = obj$f,
                               rule   = obj$rule,
                               method = obj$method)  
  a         <- (e1 - f1)/2
  b         <- a + f1
  z1_pivots <- a*x + b
  return(I2_tilde_approx(z1_pivots))
}

get_I_at_legendre_x.gsDesign <- function(obj, x) {
  # obj is GSD with 1 IA
  if (!inherits(obj, "gsDesign") || obj$k  != 2)
    stop('obj must be gsdesign k =2 ')
  order <- length(x)
  return(rep(obj$n.I[2] - obj$n.I[1], order))
}

if(FALSE){
  # Example
  z1grid = seq(-1, 3, 0.05)
  
  obj <- list(z1_nodes = c(0.28, 2.44),
              I_values = c(122, 174))
  
  I2_tilde_approx <-
    approxfun(
      x = c(obj$z1_nodes, Inf),
      y = c(obj$I_values, obj$I_values[1]),
      f = 1,
      rule = 2,
      method = "constant"
    )
  plot(x = z1grid,
       y = I2_tilde_approx(z1grid),
       type = "b")
  
  
  thetaGrid  <- seq(-0.02, 0.35, 0.01)
  order <- 12
  xw    <- gauss_legendre_xw(order)
  
  ratio <- 1
  Q <- ratio / (1 + ratio)
  HRgrid    = exp(-1 / Q * thetaGrid)
  
  EN <- ei_2_Vec(
    theta         = thetaGrid,
    f1            = head(obj$z1_nodes, 1),
    e1            = tail(obj$z1_nodes, 1),
    I1            = head(obj$I_values, 1),
    x             = xw$x,
    w             = xw$w,
    I2_constraint = "piecewise-linear",
    I2_tilde      = get_I_at_legendre_x(obj, xw$x)
  )
  
  ggplot(data.frame(x = HRgrid, y = EN)) +
    geom_line(aes(x = x, y = y))
  
  # Example 2
  
  pp <- power_2_Vec(
    theta            = thetaGrid,
    f1            = head(obj$z1_nodes, 1),
    e1            = tail(obj$z1_nodes, 1),
    I1            = head(obj$I_values, 1),
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
  
  ggplot(data.frame(x = HRgrid, y = pp)) +
    geom_line(aes(x = x, y = y)) +
    xlim(c(0.58, 0.75)) +
    ylim(c(.50, 1)) 
}

