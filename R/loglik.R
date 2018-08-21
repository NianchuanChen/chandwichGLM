# We want to calculate the vector of contributions to the GLM log-likelihood,
# from individual observations, at a given set of parameter values, including
# the dispersion parameter (if any).

# Poisson

#' @export
logLik.poisglm <- function(object, ...) {
  if (object$family$family != "poisson") {
    stop("logLik.poisglm can only be used for objects in the Poisson family")
  }
  if (nargs() > 1L) {
    warning("extra arguments discarded")
  }
  # offset
  off <- object$offset
  if (is.null(off)) {
    off <- 0
  }
  # linear predictor
  eta <- object$linear.predictors
  # Poisson mean.
  mu <- object$family$linkinv(off + eta)
  # responses
  y <- object$y
  # weights supplied to call to glm()
  w <- object$prior.weights
  # return weighted loglikelihood contributions
  return(w * stats::dpois(x = y, lambda = mu, log = TRUE))
}

#' @export
pois_glm_loglik <- function(pars, glm_object) {
  new_object <- glm_object
  class(new_object) <- "poisglm"
  # Create the linear predictor eta from the parameter values in pars
  new_object$linear.predictors <- drop(new_object$x %*% pars)
  return(logLik(new_object))
}

# Binomial


#' @export
logLik.binomglm <- function(object, ...) {
  if (object$family$family != "binomial") {
    stop("logLik.binomglm can only be used for objects in the Binomial family")
  }
  if (nargs() > 1L) {
    warning("extra arguments discarded")
  }
  # offset
  off <- object$offset
  if (is.null(off)) {
    off <- 0
  }
  # linear predictor
  eta <- object$linear.predictors
  # Binomial mean
  mu <- object$family$linkinv(off + eta)
  # responses (number of successes)
  n_successes <- object$n_successes
  # number of trials
  size <- object$n_trials
  # weights supplied to call to glm()
  w <- object$prior.weights
  # return weighted loglikelihood contributions
  return(w * stats::dbinom(x = n_successes, size = size, prob = mu, log = TRUE))
}


#' @export
binom_glm_loglik <- function(pars, glm_object) {
  new_object <- glm_object
  class(new_object) <- "binomglm"
  # Create the linear predictor eta from the parameter values in pars
  new_object$linear.predictors <- drop(new_object$x %*% pars)
  return(logLik(new_object))
}
