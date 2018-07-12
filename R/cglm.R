# ================================== cglm  ====================================

#' Fit Generalized Models to cluster correlated data
#'
#' \code{cglm} has the same functionality as \code{\link[stats]{glm}}
#' except that cluster correlated data can be analyzed using the methodology
#' of \href{http://dx.doi.org/10.1093/biomet/asm015}{Chandler and Bate (2007)},
#' implemented by the \code{\link[chandwich]{chandwich}} package.
#'
#' @inheritParams stats::glm
#' @param cluster A vector or factor indicating from which cluster the
#'   respective loglikelihood contributions from \code{loglik} originate.
#'   Must have the same length as the vector returned by \code{loglik}.
#'   If \code{cluster} is not supplied then it is set inside
#'   \code{\link[chandwich]{adjust_loglik}} under the assumption that
#'   each observation forms its own cluster.
#' @param use_alg A logical scalar.  Whether or not to provide to
#'   \code{\link[chandwich]{adjust_loglik}} via the arguments
#'   \code{alg_deriv} and \code{alg_hess} functions to evalulate
#'   the derivatives (scores) of the loglikelihood with respect to the
#'   parameters and the Hessian of the loglikelihood respectively.
#' @details Three adjustments to the independence loglikelihood described in
#'   Chandler and Bate (2007) are available.  The vertical adjustment is
#'   described in Section 6 and two horizontal adjustments are described
#'   in Sections 3.2 to 3.4.  See the descriptions of \code{type} and, for the
#'   horizontal adjustments, the descriptions of \code{C_cholesky} and
#'   \code{C_spectral}, in the documentation of
#'   \code{\link[chandwich]{adjust_loglik}}.
#'
#'   The adjustments involve first and second derivatives of the loglikelihood
#'   with respect to the model parameters.  If \code{use_alg = FALSE} then
#'   these are estimated using \code{\link[numDeriv]{jacobian}} and
#'   \code{\link[stats:optim]{optimHess}}.
#'
#'   If the argument \code{cluster} is not provided then the adjustment to the
#'   loglikelihood is based on a robust sandwich estimator for the covariance
#'   matrix of the parameter estimators, under the assumption that each
#'   observation forms its own cluster.
#'
#'   If you only wish to adjust parameter covariance matrices, that is,
#'   you are not interested in an adjusted loglikelihood, then the
#'   \href{https://CRAN.R-project.org/package=sandwich}{sandwich} package
#'   will provide equivalent results and is more flexible than
#'   \code{chandwichGLM}.  The vignettes in the sandwich package give a very
#'   clear account of this general area.
#' @return An object of class \code{"chandwich"}, returned by
#'   \code{\link[chandwich]{adjust_loglik}}.  The function
#'   \code{\link{summary}} (i.e. \code{\link[chandwich]{summary.chandwich}})
#'   can be used to obtain or print a summary of the results.
#'   The function \code{\link[chandwich]{conf_intervals}} can be used to
#'   calculate confidence intervals for each parameter.
#'   \code{\link[chandwich]{plot.confint}} can be used to illustrate the
#'   an individual confidence interval in a plot.
#'   The function \code{\link[chandwich]{conf_region}}, and its plot method
#'   \code{\link[chandwich]{plot.confreg}}, can be used to plot confidence
#'   regions for pairs of parameters.
#' @references Chandler, R. E. and Bate, S. (2007). Inference for clustered
#'   data using the independence loglikelihood. \emph{Biometrika},
#'   \strong{94}(1), 167-183. \url{http://dx.doi.org/10.1093/biomet/asm015}
#' @seealso
#'   \code{\link[stats]{glm}},
#'   \code{\link[chandwich]{chandwich}},
#'   \code{\link[chandwich]{adjust_loglik}}.
#' @examples
#' ## Example from the help file for stats::glm()
#' ## Dobson (1990) Page 93: Randomized Controlled Trial :
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' d.AD <- data.frame(treatment, outcome, counts)

#' glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())
#' cglm1 <- cglm(counts ~ outcome + treatment, family = poisson(),
#'               cluster = 1:nrow(d.AD))
#' cglm2 <- cglm(counts ~ outcome + treatment, family = poisson(),
#'               cluster = 1:nrow(d.AD), use_alg = FALSE)
#' summary(glm.D93)
#' summary(cglm1)
#' summary(cglm2)
#'
#' # Confidence intervals for each parameter
#' cints <- chandwich::conf_intervals(cglm1)
#' cints
#' plot(cints, which_par = 2)
#'
#' # A confidence region for a pair of parameters
#' conf <- chandwich::conf_region(cglm1, which_pars = 1:2)
#' plot(conf, conf = c(50, 75, 95, 99))
#' @export
cglm <- function(formula, family = gaussian, data, weights, subset, na.action,
                 start = NULL, etastart, mustart, offset, control = list(...),
                 model = TRUE, method = "glm.fit", x = FALSE, y = TRUE,
                 contrasts = NULL, cluster = NULL, use_alg = TRUE, ...) {
  #
  # 1. Fit the model using stats::glm()
  #
  # Extract the arguments from the user's call, including any in ...
  glm_call <- match.call(expand.dots = TRUE)
  # Change the first component (the function to call) to stats::glm()
  glm_call[[1L]] <- quote(stats::glm)
  # Set cluster and use_alg to NULL: they are not arguments of stats::glm()
  glm_call$cluster <- NULL
  glm_call$use_alg <- NULL
  # Add x = TRUE so that the returned object will contain the design matrix x
  glm_call$x <- TRUE
  # Call stats::glm with the user's arguments
  res <- eval.parent(glm_call)
  glm_object <- eval.parent(glm_call)
  #
  # 2. Use chandwich::adjust_loglik() to adjust the log-likelihood
  #
  # Extract the name of the GLM family
  glm_family <- res$family$family
  # Set the independence log-likelihood to be adjusted
  loglik_for_chandwich <- switch(glm_family,
                                 poisson = pois_glm_loglik,
                                 NULL)
  if (use_alg) {
    # Set a function to evaluate the score matrix
    alg_deriv_for_chandwich <- switch(glm_family,
                                      poisson = no_dispersion_glm_alg_deriv,
                                      NULL)
    # Set a function to evaluate the Hessian matrix of the loglikelihood
    alg_hess_for_chandwich <- switch(glm_family,
                                     poisson = no_dispersion_glm_alg_hess,
                                     NULL)
  } else {
    alg_deriv_for_chandwich <- NULL
    alg_hess_for_chandwich <- NULL
  }
  # Perform the adjustment using by chandwich::adjust_loglik()
  # We provide mle and hess_at_mle in order that the MLE and SE of the
  # parameters matches those obtained by stats::glm()
  res <- chandwich::adjust_loglik(loglik = loglik_for_chandwich,
                                  glm_object = glm_object,
                                  cluster = cluster,
                                  p = length(glm_object$coefficients),
                                  par_names = names(glm_object$coefficients),
                                  alg_deriv = alg_deriv_for_chandwich,
                                  alg_hess = alg_hess_for_chandwich)
  return(res)
}

# Note: pars not used because we only need the derivatives at the MLE,
# which is stored in glm_object

no_dispersion_glm_alg_deriv <- function(pars, glm_object) {
  return(sandwich::estfun(glm_object))
}

no_dispersion_glm_alg_hess <- function(pars, glm_object) {
  return(-solve(vcov(glm_object)))
}