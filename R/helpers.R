#----------------------------------------------------------------------------
#' Simulate a Gaussian linear model
#'
#' Generate data from a (sparse) Gaussian linear model.
#' The covariates are correlated Gaussian variables. The user
#' may control the signal-to-noise and the number of nonzero coefficients.
#'
#' @param n number of observations
#' @param p number of covariates
#' @param p_sig number of true nonzero coefficients (signals)
#' @param SNR signal-to-noise ratio
#' @return a list with the following elements:
#' \itemize{
#' \item \code{y}: the response variable
#' \item \code{X}: the matrix of covariates
#' \item \code{beta_true}: the true regression coefficients (including an intercept)
#' \item \code{Ey_true}: the true expectation of \code{y} (\code{X\%*\%beta_true})
#' \item \code{sigma_true}: the true error standard deviation
#' }
#'
#' @details The true regression coefficients include an intercept (-1) and
#' otherwise the \code{p_sig} nonzero coefficients are half equal to 1 and
#' half equal to -1.
#'
#' @examples
#' # Simulate data:
#' dat = simulate_lm(n = 100, p = 10)
#' names(dat) # what is returned
#'
#' @importFrom stats arima.sim model.matrix
#' @export
simulate_lm = function(n, p,
                       p_sig = min(5, p/2),
                       SNR = 1){

  #----------------------------------------------------------------------------
  # Check:
  if(p_sig > p)
    stop("The number of significant covariates cannot exceed p")
  #----------------------------------------------------------------------------
  # Simulate a design matrix with correlated predictors:
  ar1 = 0.75
  X = cbind(1,
            t(apply(matrix(0, nrow = n, ncol = p), 1, function(x)
              arima.sim(n = p, list(ar = ar1), sd = sqrt(1-ar1^2)))))
  colnames(X) = colnames(model.matrix(~X - 1)) # naming

  # Shuffle the non-intercept columns:
  X[,-1] = X[,sample(2:(p+1))]
  #----------------------------------------------------------------------------
  # True coefficients:
  beta_true = c(-1,
                rep(1, ceiling(p_sig/2)),
                rep(-1, floor(p_sig/2)),
                rep(0, p - ceiling(p_sig/2) - floor(p_sig/2)))

  # True expectation of y:
  Ey_true = X%*%beta_true

  # True SD:
  sigma_true = sd(Ey_true)/sqrt(SNR)

  # Observed y:
  y = Ey_true + sigma_true*rnorm(n)

  return(list(
    y = y, X = X,
    beta_true = beta_true,
    Ey_true = Ey_true,
    sigma_true = sigma_true
  ))
}
#----------------------------------------------------------------------------
#' Simulate a Gaussian linear model with random intercepts
#'
#' Generate data from a (sparse) Gaussian linear model with random
#' intercepts, i.e., for repeated measurements of (longitudinal) data.
#' The covariates are correlated Gaussian variables. The user
#' may control the signal-to-noise, the number of nonzero coefficients, and
#' the intraclass correlation
#'
#' @param n number of subjects
#' @param p number of covariates
#' @param m number of observations per subject
#' @param rho intraclass correlation coefficient
#' @param p_sig number of true nonzero coefficients (signals)
#' @param SNR signal-to-noise ratio
#' @return a list with the following elements:
#' \itemize{
#' \item \code{Y}: the matrix of response variables
#' \item \code{X}: the matrix of covariates
#' \item \code{beta_true}: the true regression coefficients (including an intercept)
#' \item \code{Ey_true}: the true expectation of \code{y} (\code{X\%*\%beta_true})
#' \item \code{m_scale_true}: the true Mahalanobis scale factor, 1/(sigma_e^2/sigma_u^2 + m)
#' }
#'
#' @details The true regression coefficients include an intercept (-1) and
#' otherwise the \code{p_sig} nonzero coefficients are half equal to 1 and
#' half equal to -1.
#'
#' @examples
#' # Simulate data:
#' dat = simulate_lm_randint(n = 100, p = 10, m = 4)
#' names(dat) # what is returned
#'
#' @importFrom stats arima.sim
#' @export
simulate_lm_randint = function(n, p, m,
                               rho = 0.25,
                               p_sig = min(5, p/2),
                               SNR = 1){

  #----------------------------------------------------------------------------
  # Check:
  if(p_sig > p)
    stop("The number of significant covariates cannot exceed p")
  #----------------------------------------------------------------------------
  # Simulate a design matrix with correlated predictors:
  ar1 = 0.75
  X = cbind(1,
            t(apply(matrix(0, nrow = n, ncol = p), 1, function(x)
              arima.sim(n = p, list(ar = ar1), sd = sqrt(1-ar1^2)))))
  colnames(X) = colnames(model.matrix(~X - 1)) # naming

  # Shuffle the non-intercept columns:
  X[,-1] = X[,sample(2:(p+1))]
  #----------------------------------------------------------------------------
  # True coefficients:
  beta_true = c(-1,
                rep(1, ceiling(p_sig/2)),
                rep(-1, floor(p_sig/2)),
                rep(0, p - ceiling(p_sig/2) - floor(p_sig/2)))

  # True expectation of y:
  Ey_true = X%*%beta_true

  # Compute the variances:
  sigma_tot = sd(Ey_true)/sqrt(SNR)
  sigma_u = sqrt(sigma_tot^2*rho) # rho = sigma_u^2/(sigma_u^2 + sigma_e^2)
  sigma_e = sqrt(sigma_tot^2 - sigma_u^2) # sigma_tot^2 = sigma_u^2 + sigma_e^2

  # Random effects:
  u_i = rnorm(n = n, mean = 0, sd = sigma_u)

  # Observed Y:
  Y = sapply(1:n, function(i){
    Ey_true[i] + u_i[i] + sigma_e*rnorm(n = m)
  })

  # True Mahalanobos scale factor:
  m_scale_true = 1/(sigma_e^2/sigma_u^2 + m)

  return(list(
    Y = Y, X = X,
    beta_true = beta_true,
    Ey_true = Ey_true,
    m_scale_true = m_scale_true
  ))
}
#' Get posterior predictive draws and log-predictive density
#'
#' Given posterior samples from the conditional mean and
#' conditional standard deviation of a Gaussian regression model,
#' compute posterior predictive draws and the log-predictive
#' density (lpd) at the observed data points (draw-by-draw).
#'
#' @param post_y_hat \code{nsave x n} draws of the conditional mean
#' @param post_sigma \code{nsave} draws of the conditional standard deviation
#' @param yy optional \code{n}-dimensional vector of data points; if NULL,
#' the lpd is not computed
#' @return a list with the following elements:
#' \itemize{
#' \item \code{post_y_pred}: \code{nsave x n} posterior predictive draws
#' \item \code{post_lpd}: \code{nsave x n} evaluations of the log-predictive density
#' }
#' @examples
#' # Simulate data:
#' dat = simulate_lm(n = 100, p = 10)
#' y = dat$y; X = dat$X
#'
#' # Fit a Bayesian linear model:
#' fit = bayeslm::bayeslm(y ~ X[,-1], # intercept already included
#'               N = 1000, burnin = 500) # small sim for ex
#'
#' # Compute predictive draws:
#' temp = post_predict(post_y_hat = tcrossprod(fit$beta, X),
#'                     post_sigma = fit$sigma,
#'                     yy = y)
#' names(temp) # what is returned
#'
#' # Compare fitted values to the truth:
#' plot(dat$Ey_true,
#'      colMeans(temp$post_y_pred),
#'      xlab = 'True E(y | x)', ylab = 'Fitted')
#' abline(0,1)
#'
#' @export
post_predict = function(post_y_hat,
                        post_sigma,
                        yy = NULL){
  # Dimensions:
  S = length(post_sigma) # number of simulations
  n = ncol(post_y_hat)   # number of observations

  # Check:
  if(nrow(post_y_hat) != S)
    stop('nrow(post_y_hat) must equal length(post_sigma)')

  if(!is.null(yy) & length(yy) !=n)
    stop('length(yy) must equal ncol(post_y_hat)')

  # Storage:
  post_y_pred = array(NA, c(S, n)) # posterior predictive draws
  if(!is.null(yy)) {
    post_lpd = array(NA, c(S, n)) # pointwise log-likelihood evaluations
  } else post_lpd = NULL

  # Draw-by-draw:
  for(s in 1:S){
    post_y_pred[s,] = rnorm(n = n,
                            mean = post_y_hat[s,],
                            sd = post_sigma[s])
    if(!is.null(yy)) {
      post_lpd[s,] = dnorm(yy,
                           mean = post_y_hat[s,],
                           sd = post_sigma[s],
                           log = TRUE)
    }
  }

  # Return:
  list(post_y_pred = post_y_pred,
       post_lpd = post_lpd
  )
}
#' Compute the pseudo X and Y variables for LMM summarization
#'
#' Given output from a random intercept model, compute
#' the "X" and "Y" variables needed for the least squares
#' reparametrization.
#'
#' @param XX (\code{n x p}) matrix of covariates
#' @param post_y_pred (\code{nsave x m x n}) array of posterior predictive draws
#' @param post_sigma_e (\code{nsave}) draws from the posterior distribution
#' of the observation error SD
#' @param post_sigma_u (\code{nsave}) draws from the posterior distribution
#' of the random intercept SD
#' @param post_y_pred_sum (\code{nsave x n}) matrix of the posterior predictive
#' draws summed over the replicates within each subject (optional)
#' @return list of the covariates and the response
#' @import Matrix
#' @export
getXY_randint = function(XX, post_y_pred,
                         post_sigma_e,
                         post_sigma_u,
                         post_y_pred_sum = NULL){

  # Get dimensions:
  n = nrow(XX); p = ncol(XX); m = dim(post_y_pred)[2]
  S = nrow(post_y_pred) # number of posterior simulations

  # This can be slow if n is large:
  if(is.null(post_y_pred_sum)){
    post_y_pred_sum = apply(post_y_pred, c(1,3), sum)
  }

  # Estimated Mahalanobis weight matrix (block diagonal), ignoring sigma_e^2:
  Omega_hat = diag(1, m) -
    mean(1/(post_sigma_e^2/post_sigma_u^2 + m))

  # Blocked, then matrix square-root:
  #Omega_block_hat = bdiag(lapply(1:n, function(i){Omega_hat}))
  #sqrt_Omega_block_hat = Matrix::chol(Omega_block_hat)

  # Matrix square-root, then blocked:
  sqrt_Omega_hat = chol(Omega_hat)
  sqrt_Omega_block_hat = bdiag(lapply(1:n, function(i){sqrt_Omega_hat}))

  # Stacked X-matrix:
  X_stack = apply(XX, 2, function(x)
    matrix(rep(x, each = m), nrow = m))

  # Design matrix:
  X_Omega_hat = as.matrix(sqrt_Omega_block_hat%*%X_stack)

  # Response matrix, ignoring sigma_e^2::
  Omega_y_hat = colMeans(post_y_pred) -
    rep(colMeans(1/(post_sigma_e^2/post_sigma_u^2 + m)*post_y_pred_sum),
        each = m)
  # y_Omega_hat = as.matrix(crossprod(Matrix::solve(sqrt_Omega_block_hat),
  #                                   matrix(Omega_y_hat)))
  y_Omega_hat = matrix(crossprod(Matrix::solve(sqrt_Omega_hat),
                                 Omega_y_hat))

  return(list(
    X_star = X_Omega_hat,
    y_star = y_Omega_hat
  ))
}
#' Compute the pseudo X and Y variables for LMM summarization
#'
#' Given output from a random intercept model, compute
#' the "X" and "Y" variables needed for the least squares
#' reparametrization.
#'
#' @param YY \code{m x n} matrix of response variables
#' @param y_hat \code{n x 1} vector of fitted values (common across the \code{m} replicates)
#' @param m_scale the Mahalanobis scale factor 1/(sigma_e^2/sigma_u^2 + m)
#' @return The Mahalanobis loss (scalar)
loss_maha = function(YY, y_hat, m_scale){
  m = nrow(YY); n = ncol(YY);
  1/(n*m)*as.numeric(
    sum(YY^2) - 2*crossprod(y_hat, colSums(YY)) + m*sum((y_hat)^2) -
      m_scale*(sum((colSums(YY) - m*y_hat)^2))
  )
}
#----------------------------------------------------------------------------
#' Sampler for horseshoe prior parameters
#'
#' Compute one draw of horseshoe prior parameters (local precision, global
#' precision, and local and global parameter expansion terms).
#'
#' @param omega \code{n x p} matrix of errors
#' @param params list of parameters to update
#' @param sigma_e the observation error standard deviation; for (optional) scaling purposes
#' @return List of relevant components in \code{params}: \code{sigma_wt}, the \code{n x p} matrix of standard deviations,
#' and the local and global precisions and parameter-expansion terms.
#'
#' @note To avoid scaling by the observation standard deviation \code{sigma_e},
#' simply use \code{sigma_e = 1} in the functional call.
#'
sampleHS = function(omega, params,  sigma_e = 1){

  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)

  # For numerical reasons, keep from getting too small
  hsOffset = tcrossprod(rep(1,n), apply(omega, 2, function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))
  hsInput2 = omega^2 + hsOffset

  # Local scale params:
  params$tauLambdaj = matrix(rgamma(n = n*p, shape = 1, rate = params$xiLambdaj + hsInput2/2), nrow =n)
  params$xiLambdaj = matrix(rgamma(n = n*p, shape = 1, rate = params$tauLambdaj + tcrossprod(rep(1,n), params$tauLambda)), nrow =n)

  # Global scale params:
  params$tauLambda = rgamma(n = p, shape = 0.5 + n/2, colSums(params$xiLambdaj) + params$xiLambda)
  params$xiLambda = rgamma(n = p, shape = 1, rate = params$tauLambda + 1/sigma_e^2)

  params$sigma_wt = 1/sqrt(params$tauLambdaj)

  return(params)
}
#----------------------------------------------------------------------------
#' Initialize the horseshoe prior parameters
#'
#' Compute the standard deviations, local and global precisions, and
#' parameter expansion terms for the horseshoe prior initialization.
#'
#' @param omega \code{n x p} matrix of evolution errors
#' @return List of relevant components: \code{sigma_wt}, the \code{n x p} matrix of standard deviations,
#' and the local and global precisions and parameter-expansion terms.
#'
initHS = function(omega){

  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)

  # Local precision:
  tauLambdaj = 1/omega^2;
  xiLambdaj = 1/(2*tauLambdaj); # px term

  # Global precision
  tauLambda = 1/(2*colMeans(xiLambdaj));
  xiLambda = 1/(tauLambda + 1) # px term

  # Parameters to store/return:
  return(list(sigma_wt = 1/sqrt(tauLambdaj),
              tauLambdaj = tauLambdaj,
              xiLambdaj = xiLambdaj,
              tauLambda = tauLambda,
              xiLambda = xiLambda))
}
