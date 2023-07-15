#' Bayesian linear mixed models (LMMs)
#'
#' Efficient blocked Gibbs sampler for a Bayesian linear
#' regression model with random intercepts. The fixed effects
#' (regression coefficients) and the random effects (intercepts)
#' are sampled *jointly* for Monte Carlo efficiency, while the variance
#' components are sampled in a separate block. The model uses a
#' horseshoe prior for the fixed effects (regression coefficients).
#'
#' @param Y (\code{m x n}) matrix of response variables
#' @param X (\code{n x p}) matrix of covariates
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burn-in)
#' @return a list containing the following elements:
#' \itemize{
#' \item \code{coefficients} the estimated regression coefficients (posterior means)
#' \item \code{fitted.values} the fitted values (posterior means)
#' \item \code{post_y_pred} posterior predictive draws
#' \item \code{post_y_pred_sum} posterior predictive totals for each subject
#' \item \code{post_beta} posterior draws of the regression coefficients
#' \item \code{post_u} posterior draws of the random intercepts
#' \item \code{post_sigma_u} posterior draws of the random intercept standard deviation
#' \item \code{post_sigma_e} posterior draws of the observation error standard deviation
#' \item \code{post_lpd} posterior draws of the log-likelihood (i.e., the log-likelihood
#' evaluated at each posterior draw of the model parameters)
#' }
#' @examples
#' # Simulate some data:
#' dat = simulate_lm_randint(n = 100, # subjects
#'                           p = 15,  # covariates
#'                           m = 4)   # replicates per subject
#' Y = dat$Y; X = dat$X
#'
#' # Dimensions:
#' dim(Y) # m x n
#' dim(X) # n x p
#'
#' # Fit the model:
#' fit = bayeslmm(Y = Y, X = X) # should take a few seconds
#' names(fit) # what is returned
#'
#' # Estimated coefficients:
#' coef(fit)
#'
#' # Compare to ground truth:
#' plot(coef(fit), dat$beta_true,
#'      main = 'True and estimated coefficients',
#'      xlab = 'Estimated', ylab = 'True')
#' abline(0,1)
#'
#' # 90% credible intervals:
#' ci_beta = t(apply(fit$post_beta, 2,
#'                   quantile, c(0.05, 0.95)))
#'
#' # Fitted values (m x n):
#' dim(fitted(fit))
#'
#' # MCMC diagnostics:
#' plot(as.ts(fit$post_beta[,1:6]))
#'
#' @importFrom truncdist rtrunc
#' @importFrom stats rnorm mad rgamma dnorm
#' @export
bayeslmm = function(Y, X,
                    nsave = 1000, nburn = 1000){

  # Storage:
  m = nrow(Y) # Number of replicates per subject
  n = ncol(Y) # Number of subjects
  p = ncol(X) # Number of covariates

  # Check for missingness:
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna))
  if(any.missing){
    # Indices of missing values:
    na.ind = which(is.na(Y), arr.ind = TRUE);
    Y[na.ind] = mean(Y, na.rm=TRUE)
  }

  # Recurring terms:
  XtX = crossprod(X)
  XtY = crossprod(X, t(Y))
  #---------------------------------------------------------------------------
  # Initialize:
  #---------------------------------------------------------------------------
  # Initialize the regression coefficients using the subject means:
  y_bar = colMeans(Y)
  beta = coef(lm(y_bar ~ X - 1))

  # In case there are NAs...likely due to p > n
  if(any(is.na(beta))){
    #beta[is.na(beta)] = 0
    beta = chol2inv(chol(XtX + diag(mean(XtX^2), p)))%*%rowMeans(XtY)
  }

  # Form in the "Y" dimension:
  XB = rep(X%*%beta, each = m) # does not need to be a matrix

  # Initialize the random effects based on the residuals:
  u_i = colMeans(Y - XB)
  sigma_u = sd(u_i)

  # And the observation SD:
  sigma_e = sd(Y - XB - rep(u_i, each = m))

  # Parameters for the beta shrinkage priors:
  hsparams = initHS(beta)
  sigma_beta = as.numeric(hsparams$sigma_wt)

  # Store the MCMC output in separate arrays (better computation times)
  post_y_pred = array(NA, c(nsave, m, n))
  post_y_pred_sum = array(NA, c(nsave, n))
  post_beta = array(NA, c(nsave, p))
  post_u = array(NA, c(nsave, n))
  post_sigma_u = post_sigma_e = numeric(nsave)
  post_lpd = array(NA, c(nsave, n))

  # Run the MCMC:
  for(nsi in 1:(nburn + nsave)){

    # Step 0: impute if needed
    if(any.missing){
      Y[na.ind] = matrix(XB +  rep(u_i, each = m),
                         nrow = m)[na.ind] +
        sigma_e*rnorm(n = nrow(na.ind))
      XtY = crossprod(X, t(Y))
    }

    # Step 1: sample the beta coefficients
    # after integrating out the random effects

    # First, compute the matrix inverse:
    Omega = sigma_e^-2*(diag(m) -
                          1/(sigma_e^2/sigma_u^2 + m))

    # m=2
    #detSig = (sigma_u^2 + sigma_e^2)^2 - sigma_u^4
    #Omega = 1/(detSig)*matrix(c(sigma_u^2 + sigma_e^2, -sigma_u^2,
    #                            -sigma_u^2, sigma_u^2 + sigma_e^2), nrow = 2)

    # Terms needed:
    chQbeta = chol(XtX*sum(Omega) + diag(sigma_beta^-2))
    ell_beta = XtY%*%colSums(Omega) #(crossprod(X, Y[1,])*sum(Omega[,1]) + crossprod(X, Y[2,])*sum(Omega[,2]))

    # And sample:
    beta = backsolve(chQbeta,
                     forwardsolve(t(chQbeta), ell_beta) +
                       rnorm(p))
    # Fitted part:
    XB = rep(X%*%beta, each = m)

    # Step 2: sample the random effects
    postSD = 1/sqrt(m*sigma_e^-2 + sigma_u^-2)
    postMean = sigma_e^-2*colSums(Y - XB)*postSD^2
    u_i = rnorm(n = n, mean = postMean, sd = postSD)

    # Step 3: sample the random effects SD
    u_offset = any(u_i^2 < 10^-16)*max(10^-8, mad(u_i)/10^6) # offset for numerical reasons
    sse_u = sum(u_i^2) + u_offset

    sigma_u = 1/sqrt(truncdist::rtrunc(n = 1, "gamma",
                                       a = (1/100)^2, b = Inf,
                                       shape = (n+1)/2,
                                       rate = 1/2*sse_u))

    # Step 4: sample the observation SD
    Yhat = XB + rep(u_i, each = m)
    sigma_e = 1/sqrt(rgamma(n = 1, shape = n*m/2,
                            rate = sum((Y - Yhat)^2)/2))


    # Step 5: sample the beta variance parameters
    hsparams = sampleHS(beta, params = hsparams)
    sigma_beta = as.numeric(hsparams$sigma_wt)

    # Store the MCMC output:
    if(nsi > nburn){

      # Save the MCMC samples:
      post_y_pred[nsi - nburn,,] = rnorm(n = m*n,
                                         mean = Yhat,
                                         sd = sigma_e)
      post_y_pred_sum[nsi - nburn,] = colSums(post_y_pred[nsi - nburn,,])
      post_beta[nsi - nburn,] = beta
      post_u[nsi - nburn,] = u_i
      post_sigma_u[nsi - nburn] = sigma_u
      post_sigma_e[nsi - nburn] = sigma_e
      post_lpd[nsi - nburn,] = colSums(dnorm(Yna,
                                             mean = XB + rep(u_i, each = m),
                                             sd = sigma_e, log = TRUE), na.rm=TRUE)
    }
  }

  list(
    coefficients = colMeans(post_beta),
    fitted.values = colMeans(post_y_pred),
    post_y_pred = post_y_pred,
    post_y_pred_sum = post_y_pred_sum,
    post_beta = post_beta,
    post_u = post_u,
    post_sigma_u = post_sigma_u,
    post_sigma_e = post_sigma_e,
    post_lpd = post_lpd
  )
}
