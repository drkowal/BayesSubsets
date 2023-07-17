#' Branch-and-bound algorithm for linear subset search
#'
#' Search for the "best" (according to residual sum of squares)
#' linear subsets of each size. The algorithm may collect
#' the \code{n_best} "best" subsets of each size, include or
#' exclude certain variables automatically, and apply
#' forward, backward, or exhaustive search.
#'
#' @param yy vector of response variables
#' @param XX matrix of covariates
#' @param wts vector of observation weights (for weighted least squares)
#' @param n_best number of "best" subsets for each model size
#' @param to_include indices of covariates to include in *all* subsets
#' @param to_exclude indices of covariates to exclude from *all* subsets
#' @param searchtype use exhaustive search, forward selection, backward selection or sequential replacement to search
#' @return \code{inclusion_index}: the matrix of inclusion indicators (columns) for
#' each subset returned (rows)
#'
#' @examples
#' # Simulate data:
#' dat = simulate_lm(n = 100, p = 10)
#'
#' # Run branch-and-bound:
#' indicators = branch_and_bound(yy = dat$y, XX = dat$X)
#'
#' # Inspect:
#' head(indicators)
#'
#' # Dimensions:
#' dim(indicators)
#'
#' # Model sizes:
#' rowSums(indicators)
#'
#' @importFrom leaps regsubsets
#' @export
branch_and_bound = function(yy,
                          XX,
                          wts = NULL,
                          n_best = 15,
                          to_include = 1,
                          to_exclude = NULL,
                          searchtype = 'exhaustive'
){

  # Get dimensions of the covariate matrix:
  n = nrow(XX); p = ncol(XX)

  # Useful step:
  xnames = colnames(XX);
  colnames(XX) = 1:p

  # Some basic checks:
  if(length(yy) !=n)
    stop('length of yy must equal the number of rows of XX')

  if(!is.null(to_include) && !is.null(to_exclude)){
    # Quick check:
    if(any(!is.na(match(to_include, to_exclude))))
      stop('Cannot include and exclude the same variables!')

    # Reindex to account for the excluded terms:
    to_include = match(to_include, (1:p)[-to_exclude])
  }

  # And a warning if the number of predictors is too large
  if(p - length(to_exclude) > 40 && searchtype == 'exhaustive'){
    warning("Inadvisable number of predictors for exhaustive search: slow computing! Try pre-screening to exclude some variables.")
  }

  if(!is.null(wts)){
    if(length(wts) != n)
      stop('wts must have length n')
  } else wts = rep(1,n)

  # Delete the excluded columns, if any:
  if(!is.null(to_exclude))
    XX = XX[,-to_exclude]

  # Branch-and-bound search:
  fit_all = leaps::regsubsets(x = XX, y = yy,
                              weights = wts,
                              nbest = n_best, nvmax = p,
                              method = searchtype,
                              intercept = FALSE,
                              really.big = TRUE,
                              force.in = to_include)

  # Indicator matrix of variable inclusions for each subset:
  temp_inclusion_index = summary(fit_all)$which

  # Inclusion index: need to adjust for excluded variables
  inclusion_index = matrix(FALSE,   # FALSE means excluded
                       nrow = nrow(temp_inclusion_index),
                       ncol = p)

  # Update the non-excluded entries:
  inclusion_index[,match(colnames(temp_inclusion_index),
                     1:p)] = temp_inclusion_index

  # And add the variable names:
  colnames(inclusion_index) = xnames

  # If we have any auto-include, then add those at the top:
  if(!is.null(to_include)){
    force_in = rep(FALSE, p)
    force_in[to_include] = TRUE
    inclusion_index = rbind(force_in, inclusion_index)
  }

  return(inclusion_index)
}
#' Compute the predictive and empirical cross-validated squared error loss
#'
#' Use posterior predictive draws and a sampling-importance resampling (SIR)
#' algorithm to approximate the cross-validated predictive squared error loss.
#' The empirical squared error loss (i.e., the usual quantity in cross-validation)
#' is also returned. The values are computed relative to the "best"
#' subset according to minimum empirical squared error loss.
#' Specifically, these quantities are computed for a collection of
#' linear models that are fit to the Bayesian model output, where
#' each linear model features a different subset of predictors.
#'
#' @param post_y_pred \code{S x n} matrix of posterior predictive
#' at the given \code{XX} covariate values
#' @param post_lpd \code{S} evaluations of the log-likelihood computed
#' at each posterior draw of the parameters
#' @param XX \code{n x p} matrix of covariates at which to evaluate
#' @param yy \code{n}-dimensional vector of response variables
#' @param indicators \code{L x p} matrix of inclusion indicators (booleans)
#' where each row denotes a candidate subset
#' @param post_y_hat \code{S x n} matrix of posterior fitted values
#' at the given \code{XX} covariate values
#' @param K number of cross-validation folds
#' @param sir_frac fraction of the posterior samples to use for SIR
#' @return a list with two elements: \code{pred_loss} and \code{emp_loss}
#' for the predictive and empirical loss, respectively, for each subset.
#' @details The quantity \code{post_y_hat} is the conditional expectation of the
#' response for each covariate value (columns) and using the parameters sampled
#' from the posterior (rows). For Bayesian linear regression, this term is
#' \code{X \%*\% beta}. If unspecified, the algorithm will instead use \code{post_y_pred},
#' which is still correct but has lower Monte Carlo efficiency.
#'
#' @export
pp_loss = function(post_y_pred,
                   post_lpd,
                   XX,
                   yy,
                   indicators,
                   post_y_hat = NULL,
                   K = 10,
                   sir_frac = 0.5
){

  # Get dimensions of the covariate matrix:
  n = nrow(XX); p = ncol(XX)

  # And some other dimensions
  S = nrow(post_y_pred) # number of posterior simulations
  L = nrow(indicators)  # number of subsets to consider

  # Some basic checks:
  if(length(yy) !=n)
    stop('length of yy must equal the number of rows of XX')

  if(ncol(post_y_pred) != n)
    stop('number of columns of post_y_pred must equal the number of rows of XX')

  if(nrow(post_lpd) != S)
    stop('post_lpd must have the same number of rows (simulations) as post_y_pred')

  if(is.null(post_y_hat)){
    post_y_hat = post_y_pred # just use the predictive draws
  } else {
    # Check:
    if(nrow(post_y_hat) != S)
      stop('post_y_hat must have the same number of rows (simulations) as post_y_pred')
    if(ncol(post_y_hat) != n)
      stop('number of columns of post_y_hat must equal the number of rows of XX')
  }

  if(ncol(indicators) != p)
    stop('indicators must have the same number of columns as XX has rows')

  # Samples for sampling-importance resampling (SIR):
  S_sir = ceiling(S*sir_frac) # should be a small fraction of S

  # Indices of the kth *holdout* set, k=1:K
  I_k = split(sample(1:n),
              1:K)

  # For each holdout set, obtain samples from the predictive distribution
  # conditional only on the training data; then compare with test data
  # and posterior predictive distribution of the test data
  pred_loss = array(0, c(S_sir, L))
  emp_loss = numeric(L)
  for(k in 1:K){

    # Log importance weights:
    logw_k = -rowSums(as.matrix(post_lpd[,I_k[[k]]]))
    logw_k = logw_k - min(logw_k) # useful numerical adjustment

    # SIR indices:
    samp_inds = sample(1:S,
                       S_sir,
                       prob=exp(logw_k),
                       replace = FALSE)

    # Testing (out) and training (in) predictive distributions:
    post_Y_out = post_y_pred[samp_inds, I_k[[k]]]
    post_Y_in = post_y_pred[samp_inds, -I_k[[k]]]

    # Testing quantities of interest:
    n_out = length(I_k[[k]]); # number of out-of-sample (testing) points
    X_out = matrix(XX[I_k[[k]],], nrow = n_out); # corresponding x-variables
    XtX_out = crossprod(X_out) # X'X for these testing points
    X_in = XX[-I_k[[k]],] # in-sample x-variables
    postY2_out = rowMeans(post_Y_out^2) # squared predictive variables summed across n_out

    # Predictive mean (training points):
    y_hat_in = colMeans(post_y_hat[samp_inds,
                                   -I_k[[k]]])
    # Response variables on the testing data:
    yy_out = yy[I_k[[k]]]

    # Regression for each model:
    beta_path_in = matrix(0, nrow = p, ncol = L)
    for(ell in 1:L){
      beta_path_in[indicators[ell,], ell] = coef(
        lm(y_hat_in ~ X_in[,indicators[ell,]] - 1)
      )
    }

    # Out-of-sample predictive loss:
    pred_loss = pred_loss + 1/K*(
      apply(beta_path_in, 2, function(beta_ell){
        postY2_out +
          1/n_out*as.numeric(crossprod(beta_ell, XtX_out)%*%beta_ell) -
          2/n_out*tcrossprod(post_Y_out, t(X_out%*%beta_ell))
      })
    )

    # Out-of-sample empirical loss:
    emp_loss = emp_loss + 1/K*(
      apply(beta_path_in, 2, function(beta_ell){
        mean(yy_out^2) +
          1/n_out*as.numeric(crossprod(beta_ell, XtX_out)%*%beta_ell) -
          2/n_out*crossprod(yy_out, X_out%*%beta_ell)
      })
    )
  }

  # Best subset by out-of-sample empirical loss:
  ell_ref = which.min(emp_loss)

  # Adjust the empirical loss relative to this term:
  emp_loss = 100*(emp_loss - emp_loss[ell_ref])/emp_loss[ell_ref]

  # Percent difference in predictive loss relative to best model:
  pred_loss = apply(pred_loss, 2, function(ploss)
    100*(ploss - pred_loss[,ell_ref])/pred_loss[,ell_ref])

  return(
    list(pred_loss = pred_loss,
       emp_loss = emp_loss)
  )
}
#' Compute the predictive squared error loss on *new* testing points
#'
#' Use posterior predictive draws at new \code{XX} points, compute
#' the predictive squared error loss. The values are
#' computed relative to the largest subset provided, which is typically
#' the full set of covariates (and also the minimizer of the expected
#' predictive loss). These quantities are computed for a collection of
#' linear models that are fit to the Bayesian model output, where
#' each linear model features a different subset of predictors.
#'
#' @param post_y_pred \code{S x n} matrix of posterior predictive
#' at the given \code{XX} covariate values
#' @param XX \code{n x p} matrix of covariates at which to evaluate
#' @param indicators \code{L x p} matrix of inclusion indicators (booleans)
#' where each row denotes a candidate subset
#' @param post_y_hat \code{S x n} matrix of posterior fitted values
#' at the given \code{XX} covariate values
#' @return \code{pred_loss}: the predictive loss for each subset.
#' @details The quantity \code{post_y_hat} is the conditional expectation of the
#' response for each covariate value (columns) and using the parameters sampled
#' from the posterior (rows). For Bayesian linear regression, this term is
#'  \code{X \%*\% beta}. If unspecified, the algorithm will instead use \code{post_y_pred},
#' which is still correct but has lower Monte Carlo efficiency.
#'
#' @export
pp_loss_out = function(post_y_pred,
                       XX,
                       indicators,
                       post_y_hat = NULL
){

  # Get dimensions of the covariate matrix:
  n = nrow(XX); p = ncol(XX)

  # And some other dimensions
  S = nrow(post_y_pred) # number of posterior simulations
  L = nrow(indicators)  # number of subsets to consider

  if(ncol(post_y_pred) != n)
    stop('number of columns of post_y_pred must equal the number of rows of XX')

  if(is.null(post_y_hat)){
    post_y_hat = post_y_pred # just use the predictive draws
  } else {
    # Check:
    if(nrow(post_y_hat) != S)
      stop('post_y_hat must have the same number of rows (simulations) as post_y_pred')
    if(ncol(post_y_hat) != n)
      stop('number of columns of post_y_hat must equal the number of rows of XX')
  }

  if(ncol(indicators) != p)
    stop('indicators must have the same number of columns as XX has rows')

  # Fitted values:
  y_hat = colMeans(post_y_hat)

  # Useful terms:
  XtX = crossprod(XX) # X'X for these testing points
  post_y2_pred = rowMeans(post_y_pred^2) # squared predictive variables summed across n_out

  # Storage:
  pred_loss = array(0, c(S, L))
  for(ell in 1:L){
    # Coefficients for this model:
    beta_ell = rep(0, p)
    beta_ell[indicators[ell,]] = coef(lm(y_hat ~ XX[,indicators[ell,]] - 1))

    # And predictive loss:
    pred_loss[,ell] = post_y2_pred +
      1/n*as.numeric(crossprod(beta_ell, XtX)%*%beta_ell) -
      2/n*tcrossprod(post_y_pred, t(XX%*%beta_ell))
  }

  # Reference: largest subset
  ell_ref = L

  # Percent difference in predictive loss relative to best model:
  pred_loss = apply(pred_loss, 2, function(ploss)
    100*(ploss - pred_loss[,ell_ref])/pred_loss[,ell_ref])

  return(
    pred_loss = pred_loss
  )
}
#' Compute the predictive and empirical cross-validated loss for binary data.
#'
#' Use posterior predictive draws and a sampling-importance resampling (SIR)
#' algorithm to approximate the cross-validated predictive loss.
#' The empirical loss (i.e., the usual quantity in cross-validation)
#' is also returned. The values are computed relative to the "best"
#' subset according to minimum empirical loss.
#' Specifically, these quantities are computed for a collection of
#' linear models that are fit to the Bayesian model output, where
#' each linear model features a different subset of predictors.
#' The loss function may be chosen as cross-entropy or misclassification rate
#'
#' @param post_y_pred \code{S x n} matrix of posterior predictive
#' at the given \code{XX} covariate values
#' @param post_lpd \code{S} evaluations of the log-likelihood computed
#' at each posterior draw of the parameters
#' @param XX \code{n x p} matrix of covariates at which to evaluate
#' @param yy \code{n}-dimensional vector of response variables
#' @param indicators \code{L x p} matrix of inclusion indicators (booleans)
#' where each row denotes a candidate subset
#' @param loss_type loss function to be used:
#' "cross-ent" (cross-entropy) or "misclass" (misclassication rate)
#' @param post_y_hat \code{S x n} matrix of posterior fitted values
#' at the given \code{XX} covariate values
#' @param K number of cross-validation folds
#' @param sir_frac fraction of the posterior samples to use for SIR
#' @return a list with two elements: \code{pred_loss} and \code{emp_loss}
#' for the predictive and empirical loss, respectively, for each subset.
#' @details The quantity \code{post_y_hat} is the conditional expectation of the
#' response for each covariate value (columns) and using the parameters sampled
#' from the posterior (rows). For binary data, this is also the estimated
#' probability of "success".
#' If unspecified, the algorithm will instead use \code{post_y_pred},
#' which is still correct but has lower Monte Carlo efficiency.
#'
#' @importFrom stats glm
#' @export
pp_loss_binary = function(post_y_pred,
                          post_lpd,
                          XX,
                          yy,
                          indicators,
                          loss_type = 'cross-ent',
                          post_y_hat = NULL,
                          K = 10,
                          sir_frac = 0.5
){

  # Get dimensions of the covariate matrix:
  n = nrow(XX); p = ncol(XX)

  # And some other dimensions
  S = nrow(post_y_pred) # number of posterior simulations
  L = nrow(indicators)  # number of subsets to consider

  # Some basic checks:
  if(length(yy) !=n)
    stop('length of yy must equal the number of rows of XX')

  if(ncol(post_y_pred) != n)
    stop('number of columns of post_y_pred must equal the number of rows of XX')

  if(nrow(post_lpd) != S)
    stop('post_lpd must have the same number of rows (simulations) as post_y_pred')

  if(is.null(post_y_hat)){
    post_y_hat = post_y_pred # just use the predictive draws
  } else {
    # Check:
    if(nrow(post_y_hat) != S)
      stop('post_y_hat must have the same number of rows (simulations) as post_y_pred')
    if(ncol(post_y_hat) != n)
      stop('number of columns of post_y_hat must equal the number of rows of XX')
  }

  if(ncol(indicators) != p)
    stop('indicators must have the same number of columns as XX has rows')

  if(is.na(match(loss_type, c("cross-ent", "misclass"))))
    stop('loss_type must be one of "cross-ent" or "misclass"')

  # Define the loss function:
    # y_obs: observed data in {0,1}
    # g_hat: prediction on link scale (R^1), e,g X%*%beta
  if(loss_type == "cross-ent"){ # cross-entropy loss
    loss_fun = function(y_obs, g_hat)
      log(1 + exp(g_hat)) - g_hat*y_obs
  }
  if(loss_type == "misclass"){ # misclassification rate
    loss_fun = function(y_obs, g_hat)
      1.0*(y_obs != (g_hat > 0))
  }

  # Samples for sampling-importance resampling (SIR):
  S_sir = ceiling(S*sir_frac) # should be a small fraction of S

  # Indices of the kth *holdout* set, k=1:K
  I_k = split(sample(1:n),
              1:K)

  # For each holdout set, obtain samples from the predictive distribution
  # conditional only on the training data; then compare with test data
  # and posterior predictive distribution of the test data
  pred_loss = array(0, c(S_sir, L))
  emp_loss = numeric(L)
  for(k in 1:K){

    # Log importance weights:
    logw_k = -rowSums(as.matrix(post_lpd[,I_k[[k]]]))
    logw_k = logw_k - min(logw_k) # useful numerical adjustment

    # SIR indices:
    samp_inds = sample(1:S,
                       S_sir,
                       prob=exp(logw_k),
                       replace = FALSE)

    # Testing (out) and training (in) predictive distributions:
    post_Y_out = post_y_pred[samp_inds, I_k[[k]]]
    post_Y_in = post_y_pred[samp_inds, -I_k[[k]]]

    # Testing quantities of interest:
    n_out = length(I_k[[k]]); # number of out-of-sample (testing) points
    X_out = matrix(XX[I_k[[k]],], nrow = n_out); # corresponding x-variables
    X_in = XX[-I_k[[k]],] # in-sample x-variables

    # Predictive mean (training points):
    y_hat_in = colMeans(post_y_hat[samp_inds,
                                   -I_k[[k]]])
    # Response variables on the testing data:
    yy_out = yy[I_k[[k]]]

    # Regression for each model:
    beta_path_in = matrix(0, nrow = p, ncol = L)
    for(ell in 1:L){
      beta_path_in[indicators[ell,], ell] =
        suppressWarnings(
          coef(
            glm(y_hat_in ~ X_in[,indicators[ell,]] - 1,
                family = 'binomial')
          )
        )
    }

    # Out-of-sample predictive loss:
    pred_loss = pred_loss + 1/K*(
      apply(beta_path_in, 2, function(beta_ell){
        rowMeans(loss_fun(post_Y_out,
                          rep(X_out%*%beta_ell, each = S_sir)))
      })
    )

    # Out-of-sample empirical loss:
    emp_loss = emp_loss + 1/K*(
      apply(beta_path_in, 2, function(beta_ell){
        mean(loss_fun(yy_out, X_out%*%beta_ell))
      })
    )
  }

  # Best subset by out-of-sample empirical loss:
  ell_ref = which.min(emp_loss)

  # Adjust the empirical loss relative to this term:
  emp_loss = 100*(emp_loss - emp_loss[ell_ref])/emp_loss[ell_ref]

  # Percent difference in predictive loss relative to best model:
  pred_loss = apply(pred_loss, 2, function(ploss)
    100*(ploss - pred_loss[,ell_ref])/pred_loss[,ell_ref])

  return(
    list(pred_loss = pred_loss,
       emp_loss = emp_loss)
  )
}
#' Compute the predictive and empirical cross-validated Mahalanobis loss
#' under the random intercept model
#'
#' Use posterior predictive draws and a sampling-importance resampling (SIR)
#' algorithm to approximate the cross-validated predictive Mahalanobis loss.
#' The empirical Mahalanobis loss is also returned. The values are computed relative to the "best"
#' subset according to minimum empirical Mahalanobis loss.
#' Specifically, these quantities are computed for a collection of
#' linear models that are fit to the Bayesian model output, where
#' each linear model features a different subset of predictors.
#'
#' @param post_y_pred \code{S x m x n} matrix of posterior predictive
#' at the given \code{XX} covariate values for \code{m} replicates per subject
#' @param post_lpd \code{S} evaluations of the log-likelihood computed
#' at each posterior draw of the parameters
#' @param post_sigma_e (\code{nsave}) draws from the posterior distribution
#' of the observation error SD
#' @param post_sigma_u (\code{nsave}) draws from the posterior distribution
#' of the random intercept SD
#' @param XX \code{n x p} matrix of covariates at which to evaluate
#' @param YY \code{m x n} matrix of response variables (optional)
#' @param indicators \code{L x p} matrix of inclusion indicators (booleans)
#' where each row denotes a candidate subset
#' @param post_y_pred_sum (\code{nsave x n}) matrix of the posterior predictive
#' draws summed over the replicates within each subject (optional)
#' @param K number of cross-validation folds
#' @param sir_frac fraction of the posterior samples to use for SIR
#' @return a list with two elements: \code{pred_loss} and \code{emp_loss}
#' for the predictive and empirical loss, respectively, for each subset.
#'
#' @export
pp_loss_randint = function(post_y_pred,
                           post_lpd,
                           post_sigma_e,
                           post_sigma_u,
                           XX,
                           YY,
                           indicators,
                           post_y_pred_sum = NULL,
                           K = 10,
                           sir_frac = 0.5
){

  # Get dimensions of the covariate matrix:
  n = nrow(XX); p = ncol(XX)

  # And some other dimensions
  S = nrow(post_y_pred) # number of posterior simulations
  m = dim(post_y_pred)[2] # number of replicates per subject
  L = nrow(indicators)  # number of subsets to consider

  # Some basic checks:
  if(!is.null(YY) && ncol(YY) !=n)
    stop('the number of columns of YY must equal the number of rows of XX')

  if(dim(post_y_pred)[3] != n)
    stop('incorrect dimensions for post_y_pred')

  if(!is.null(post_lpd) && nrow(post_lpd) != S)
    stop('post_lpd must have the same number of rows (simulations) as post_y_pred')

  if(ncol(indicators) != p)
    stop('indicators must have the same number of columns as XX has rows')

  # This can be slow if n is large:
  if(is.null(post_y_pred_sum)){
    post_y_pred_sum = apply(post_y_pred, c(1,3), sum)
  }

  # Samples for sampling-importance resampling (SIR):
  S_sir = ceiling(S*sir_frac) # should be a small fraction of S

  # Indices of the kth *holdout* set, k=1:K
  I_k = split(sample(1:n),
              1:K)

  # For each holdout set, obtain samples from the predictive distribution
  # conditional only on the training data; then compare with test data
  # and posterior predictive distribution of the test data
  pred_loss = array(0, c(S_sir, L))
  emp_loss = numeric(L)
  for(k in 1:K){

    # Log importance weights:
    logw_k = -rowSums(as.matrix(post_lpd[,I_k[[k]]]))
    logw_k = logw_k - min(logw_k) # useful numerical adjustment

    # SIR indices:
    samp_inds = sample(1:S,
                       S_sir,
                       prob=exp(logw_k),
                       replace = FALSE)

    # Compute the objects needed for in-sample regression:
    objXY = getXY_randint(XX = XX[-I_k[[k]], ],
                          post_y_pred = post_y_pred[samp_inds, ,-I_k[[k]]],
                          post_sigma_e = post_sigma_e[samp_inds],
                          post_sigma_u = post_sigma_u[samp_inds],
                          post_y_pred_sum = post_y_pred_sum[samp_inds, -I_k[[k]]])
    X_star_in = objXY$X_star; y_star_in = objXY$y_star; rm(objXY)

    # Objects needed for out-of-sample regression:
    n_out = length(I_k[[k]])
    X_out = matrix(XX[I_k[[k]],], nrow = n_out); # corresponding x-variables
    Y_out = matrix(YY[,I_k[[k]]], nrow = m); # corresponding y-variables
    post_y_pred_out = post_y_pred[samp_inds, ,I_k[[k]]] # out-of-sample predictive points

    # Recurring terms:
    Y2_out = sum(Y_out^2, na.rm=TRUE); Y_sum_out = colSums(Y_out, na.rm=TRUE)
    post_y_pred_sum_out = post_y_pred_sum[samp_inds, I_k[[k]]] # out-of-sample predictive sum
    post_y2_pred_out = rowSums(post_y_pred_out^2)
    post_m_scale = 1/(post_sigma_e[samp_inds]^2/post_sigma_u[samp_inds]^2 + m)
    m_scale_hat = mean(post_m_scale)

    # For each subsets compute (i) the coefficients and (ii) the predictive and empirical losses
    for(ell in 1:L){
      # Estimate the in-sample coefficients:
      beta_ell = rep(0, p)
      beta_ell[indicators[ell,]] = coef(
        lm(y_star_in ~ X_star_in[,indicators[ell,]] - 1)
      )
      # Fitted values:
      XB_ell = X_out%*%beta_ell

      # Compute the *empirical* Mahalanobis loss:
      emp_loss[ell] = emp_loss[ell] + 1/K*1/(n_out*m)*(
        Y2_out - 2*crossprod(XB_ell, Y_sum_out) + m*sum((XB_ell)^2) -
          m_scale_hat*(sum((Y_sum_out - m*XB_ell)^2))
      )
      # Compare to the full version:
      # res = matrix(Y_out - rep(X_out%*%beta_ell, each = m))
      # Omega_hat = diag(1, m) - m_scale_hat
      # Omega_block_hat = bdiag(lapply(1:n_out, function(i){Omega_hat}))
      # t(res)%*%Omega_block_hat%*%res

      # Compute the *predictive* Mahalanobis loss:
      pred_loss[,ell] = pred_loss[,ell] + 1/K*1/(n_out*m)*(
        post_y2_pred_out - 2*tcrossprod(post_y_pred_sum_out, t(XB_ell)) + m*sum((XB_ell)^2) -
          post_m_scale*rowSums((post_y_pred_sum_out - rep(m*XB_ell, each = S_sir))^2)
      )
    }
  }

  # Under missingness, use predictive expectation:
  if(any(is.na(YY))) emp_loss = colMeans(pred_loss)

  # Best subset by out-of-sample empirical loss:
  ell_ref = which.min(emp_loss)

  # Adjust the empirical loss relative to this term:
  emp_loss = 100*(emp_loss - emp_loss[ell_ref])/emp_loss[ell_ref]

  # Percent difference in predictive loss relative to best model:
  pred_loss = apply(pred_loss, 2, function(ploss)
    100*(ploss - pred_loss[,ell_ref])/pred_loss[,ell_ref])

  return(
    list(pred_loss = pred_loss,
         emp_loss = emp_loss)
  )
}
#' Compute the acceptable family of linear subsets
#'
#' Given output from a Bayesian model and a candidate of
#' subsets, compute the *acceptable family* of subsets that
#' match or nearly match the predictive accuracy of the "best" subset.
#' The acceptable family may be computed for any set of covariate values
#' \code{XX}; if \code{XX = X} are the in-sample points, then
#' cross-validation is used to assess out-of-sample predictive performance.
#'
#' @param post_y_pred \code{S x n} matrix of posterior predictive
#' at the given \code{XX} covariate values
#' @param post_lpd \code{S} evaluations of the log-likelihood computed
#' at each posterior draw of the parameters (optional)
#' @param XX \code{n x p} matrix of covariates at which to evaluate
#' @param yy \code{n}-dimensional vector of response variables (optional)
#' @param indicators \code{L x p} matrix of inclusion indicators (booleans)
#' where each row denotes a candidate subset
#' @param eps_level probability required to match the predictive
#' performance of the "best" model (up to \code{eta_level})
#' @param eta_level allowable margin (%) between each acceptable model
#' and the "best" model
#' @param post_y_hat \code{S x n} matrix of posterior fitted values
#' at the given \code{XX} covariate values (optional)
#' @param K number of cross-validation folds (optional)
#' @param sir_frac fraction of the posterior samples to use for SIR (optional)
#' @param plot logical; if TRUE, include a plot to summarize the predictive
#' performance across candidate subsets
#' @return a list containing the following elements:
#' \itemize{
#' \item \code{all_accept}: indices (i.e., rows of \code{indicators})
#' that correspond to the acceptable subsets
#' \item \code{beta_hat_small} linear coefficients for the
#' smallest acceptable model
#' \item \code{beta_hat_min} linear coefficients for the
#' "best" acceptable model
#' \item \code{ell_small}: index (i.e., row of \code{indicators}) of the
#' smallest acceptable model
#' \item \code{ell_min}: index (i.e., row of \code{indicators}) of the
#' "best" acceptable model
#' }
#' @details When \code{XX = X} is the observed covariate values,
#' then \code{post_lpd} and \code{yy} must be provided. These
#' are used to compute the cross-validated predictive and empirical
#' squared errors; the predictive version relies on a sampling importance-resampling
#' procedure.
#'
#' When \code{XX} corresponds to a new set of covariate values, then set \code{post_lpd = NULL}
#' and \code{yy = NULL} (these are the default values).
#'
#' Additional details on the predictive and empirical comparisons are
#' in \code{pp_loss} and \code{pp_loss_out}.
#'
#' @importFrom graphics abline arrows lines
#' @importFrom stats approxfun binomial coef lm quantile runif
#' @export
accept_family = function(post_y_pred,
                   post_lpd = NULL,
                   XX,
                   yy = NULL,
                   indicators,
                   eps_level = 0.05,
                   eta_level = 0.00,
                   post_y_hat = NULL,
                   K = 10,
                   sir_frac = 0.5,
                   plot = TRUE
){

  # Get dimensions of the covariate matrix:
  n = nrow(XX); p = ncol(XX)

  # And some other dimensions
  S = nrow(post_y_pred) # number of posterior simulations
  L = nrow(indicators)  # number of subsets to consider

  # Some basic checks:
  if(!is.null(yy) && length(yy) !=n)
    stop('length of yy must equal the number of rows of XX')

  if(ncol(post_y_pred) != n)
    stop('number of columns of post_y_pred must equal the number of rows of XX')

  if(!is.null(post_lpd) && nrow(post_lpd) != S)
    stop('post_lpd must have the same number of rows (simulations) as post_y_pred')

  if(is.null(post_y_hat)){
    post_y_hat = post_y_pred # just use the predictive draws
  } else {
    # Check:
    if(nrow(post_y_hat) != S)
      stop('post_y_hat must have the same number of rows (simulations) as post_y_pred')
    if(ncol(post_y_hat) != n)
      stop('number of columns of post_y_hat must equal the number of rows of XX')
  }

  if(ncol(indicators) != p)
    stop('indicators must have the same number of columns as XX has rows')

  # Compute the empirical and predictive losses:
  if(is.null(post_lpd) && is.null(yy)){
    # Totally out-of-sample:
    pred_loss = pp_loss_out(post_y_pred = post_y_pred,
                            XX = XX,
                            indicators = indicators,
                            post_y_hat = post_y_hat)
    emp_loss = NULL
  } else {
    # In-sample, so use cross-validation:
    pred_emp_loss = pp_loss(post_y_pred = post_y_pred,
                            post_lpd = post_lpd,
                            XX = XX,
                            yy = yy,
                            indicators = indicators,
                            post_y_hat = post_y_hat,
                            K = K,
                            sir_frac = sir_frac)
    # Store these better:
    emp_loss = pred_emp_loss$emp_loss
    pred_loss = pred_emp_loss$pred_loss
    rm(pred_emp_loss)
  }

  # Indices of acceptable subsets:
  all_accept = which(colMeans(pred_loss <= eta_level)
                     >= eps_level)

  # Index of "best" subset, if possible:
  if(!is.null(emp_loss)){
    ell_min = which.min(emp_loss)
  } else ell_min = NULL

  # Subset sizes:
  subset_size = rowSums(indicators)

  # Minimize size of acceptable subset:
  min_size_accept = min(subset_size[all_accept])

  # Index of acceptable subsets with this size:
  ind_min_size_accept = which(subset_size[all_accept] == min_size_accept)

  # If more than one, select the minimum empirical loss (or expected predictive loss)
  if(length(ind_min_size_accept) > 1){
    if(!is.null(emp_loss)){
      # use empirical loss, if available
      ell_small = all_accept[ind_min_size_accept][which.min(emp_loss[all_accept[ind_min_size_accept]])]
    } else {
      # otherwise, use expected predictive loss
      ell_small = all_accept[ind_min_size_accept][which.min(colMeans(pred_loss)[all_accept[ind_min_size_accept]])]
    }
  } else ell_small = all_accept[ind_min_size_accept]
  if(is.infinite(ell_small)) ell_small = ell_min # in case there is none?

  # Compute the coefficients for each index:

  # fitted values (pseudo-response)
  if(!is.null(post_y_hat)){
    y_hat = colMeans(post_y_hat)
  } else y_hat = colMeans(post_y_pred)

  # coefficients for "best" subset
  if(!is.null(ell_min)){
    beta_hat_min = numeric(p)
    beta_hat_min[indicators[ell_min,]] =
      coef(lm(y_hat ~ XX[, indicators[ell_min,]] - 1))
  } else beta_hat_min = NULL

  # coefficients for smallest acceptable subset:
  beta_hat_small = numeric(p)
  beta_hat_small[indicators[ell_small,]] =
    coef(lm(y_hat ~ XX[, indicators[ell_small,]] - 1))

  # Now add the plot, if desired:
  if(plot){
    # 100(1 - 2*eps_level)% prediction interval for the loss
    pi_loss = t(apply(pred_loss, 2,
                       quantile, c(eps_level, 1 - eps_level)))

    jitter = runif(L, min = -0.25, max = 0.25)
    plot(subset_size + jitter,
         colMeans(pred_loss), type='p', ylim = range(pi_loss, 0), lwd=4,
         xlab = 'Subset Size',
         ylab = 'Difference in loss (%)',
         main = paste('Difference in ',K,'-fold loss (%)', sep=''))
    abline(v = 1:p, col='gray')
    abline(h = eta_level, lwd=2, lty=6)
    abline(v = (subset_size + jitter)[ell_small], lwd=5, col='darkgray')
    if(!is.null(ell_min)) abline(v = (subset_size + jitter)[ell_min], lwd=3, col='lightgray', lty=3)
    arrows(jitter + subset_size, pi_loss[,1],
           jitter + subset_size, pi_loss[,2],
           length=0.05, angle=90, code=3, lwd=4)
    if(!is.null(emp_loss)) lines(subset_size + jitter, emp_loss, type='p', lwd=4, col='gray', pch=4, cex = 1.5)
  }

  return(
    list(
      all_accept = all_accept,
      beta_hat_small = beta_hat_small,
      beta_hat_min = beta_hat_min,
      ell_small = ell_small,
      ell_min = ell_min
    )
  )
}
#' Compute the acceptable family for binary data
#'
#' Given output from a Bayesian model and a candidate of
#' subsets, compute the *acceptable family* of subsets that
#' match or nearly match the predictive accuracy of the "best" subset.
#' This function applies for binary data, such as logistic regression.
#'
#' @param post_y_pred \code{S x n} matrix of posterior predictive
#' at the given \code{XX} covariate values
#' @param post_lpd \code{S} evaluations of the log-likelihood computed
#' at each posterior draw of the parameters
#' @param XX \code{n x p} matrix of covariates at which to evaluate
#' @param indicators \code{L x p} matrix of inclusion indicators (booleans)
#' where each row denotes a candidate subset
#' @param eps_level probability required to match the predictive
#' performance of the "best" model (up to \code{eta_level})
#' @param eta_level allowable margin (%) between each acceptable model
#' and the "best" model
#' @param loss_type loss function to be used:
#' "cross-ent" (cross-entropy) or "misclass" (misclassication rate)
#' @param yy \code{n}-dimensional vector of response variables
#' @param post_y_hat \code{S x n} matrix of posterior fitted values
#' at the given \code{XX} covariate values (optional)
#' @param K number of cross-validation folds
#' @param sir_frac fraction of the posterior samples to use for SIR
#' @param plot logical; if TRUE, include a plot to summarize the predictive
#' performance across candidate subsets
#' @return a list containing the following elements:
#' \itemize{
#' \item \code{all_accept}: indices (i.e., rows of \code{indicators})
#' that correspond to the acceptable subsets
#' \item \code{beta_hat_small} linear coefficients for the
#' smallest acceptable model
#' \item \code{beta_hat_min} linear coefficients for the
#' "best" acceptable model
#' \item \code{ell_small}: index (i.e., row of \code{indicators}) of the
#' smallest acceptable model
#' \item \code{ell_min}: index (i.e., row of \code{indicators}) of the
#' "best" acceptable model
#' }
#' @details see \code{pp_loss_binary} for additional details
#' about the predictive and empirical comparisons.
#'
#' @export
accept_family_binary = function(post_y_pred,
                         post_lpd,
                         XX,
                         indicators,
                         eps_level = 0.05,
                         eta_level = 0.00,
                         loss_type = "cross-ent",
                         yy = NULL,
                         post_y_hat = NULL,
                         K = 10,
                         sir_frac = 0.5,
                         plot = TRUE
){

  # Get dimensions of the covariate matrix:
  n = nrow(XX); p = ncol(XX)

  # And some other dimensions
  S = nrow(post_y_pred) # number of posterior simulations
  L = nrow(indicators)  # number of subsets to consider

  # Some basic checks:
  if(length(yy) !=n)
    stop('length of yy must equal the number of rows of XX')

  if(ncol(post_y_pred) != n)
    stop('number of columns of post_y_pred must equal the number of rows of XX')

  if(nrow(post_lpd) != S)
    stop('post_lpd must have the same number of rows (simulations) as post_y_pred')

  if(is.null(post_y_hat)){
    post_y_hat = post_y_pred # just use the predictive draws
  } else {
    # Check:
    if(nrow(post_y_hat) != S)
      stop('post_y_hat must have the same number of rows (simulations) as post_y_pred')
    if(ncol(post_y_hat) != n)
      stop('number of columns of post_y_hat must equal the number of rows of XX')
  }

  if(ncol(indicators) != p)
    stop('indicators must have the same number of columns as XX has rows')

  # Check the loss types:
  if(is.na(match(loss_type, c("cross-ent", "misclass"))))
    stop('loss_type must be one of "cross-ent" or "misclass"')

  # Compute the empirical and predictive losses:
  pred_emp_loss = pp_loss_binary(post_y_pred = post_y_pred,
                                 post_lpd = post_lpd,
                                 XX = XX,
                                 yy = yy,
                                 indicators = indicators,
                                 loss_type = loss_type,
                                 post_y_hat = post_y_hat,
                                 K = K,
                                 sir_frac = sir_frac)

  # Store these better:
  emp_loss = pred_emp_loss$emp_loss
  pred_loss = pred_emp_loss$pred_loss
  rm(pred_emp_loss)

  # Indices of acceptable subsets:
  all_accept = which(colMeans(pred_loss <= eta_level)
                     >= eps_level)

  # Index of "best" subset, if possible:
  if(!is.null(emp_loss)){
    ell_min = which.min(emp_loss)
  } else ell_min = NULL

  # Subset sizes:
  subset_size = rowSums(indicators)

  # Minimize size of acceptable subset:
  min_size_accept = min(subset_size[all_accept])

  # Index of acceptable subsets with this size:
  ind_min_size_accept = which(subset_size[all_accept] == min_size_accept)

  # If more than one, select the minimum empirical loss (or expected predictive loss)
  if(length(ind_min_size_accept) > 1){
    if(!is.null(emp_loss)){
      # use empirical loss, if available
      ell_small = all_accept[ind_min_size_accept][which.min(emp_loss[all_accept[ind_min_size_accept]])]
    } else {
      # otherwise, use expected predictive loss
      ell_small = all_accept[ind_min_size_accept][which.min(colMeans(pred_loss)[all_accept[ind_min_size_accept]])]
    }
  } else ell_small = all_accept[ind_min_size_accept]
  if(is.infinite(ell_small)) ell_small = ell_min # in case there is none?

  # Compute the coefficients for each index:

  # fitted values (pseudo-response)
  if(!is.null(post_y_hat)){
    y_hat = colMeans(post_y_hat)
  } else y_hat = colMeans(post_y_pred)

  # coefficients for "best" subset
  if(!is.null(ell_min)){
    beta_hat_min = numeric(p)
    beta_hat_min[indicators[ell_min,]] = suppressWarnings(
      coef(glm(y_hat ~ XX[,indicators[ell_min,]] - 1,
               family = binomial())))
  } else beta_hat_min = NULL

  # coefficients for smallest acceptable subset:
  beta_hat_small = numeric(p)
  beta_hat_small[indicators[ell_small,]] = suppressWarnings(
    coef(glm(y_hat ~ XX[,indicators[ell_small,]] - 1,
             family = binomial())))

  # Now add the plot, if desired:
  if(plot){
    # 100(1 - 2*eps_level)% prediction interval for the loss
    pi_loss = t(apply(pred_loss, 2,
                      quantile, c(eps_level, 1 - eps_level)))

    jitter = runif(L, min = -0.25, max = 0.25)
    plot(subset_size + jitter,
         colMeans(pred_loss), type='p', ylim = range(pi_loss, 0), lwd=4,
         xlab = 'Subset Size',
         ylab = 'Difference in loss (%)',
         main = paste('Difference in ',K,'-fold loss (%)', sep=''))
    abline(h = eta_level, lwd=2, lty=6)
    abline(v = (subset_size + jitter)[ell_small], lwd=5, col='darkgray')
    if(!is.null(ell_min)) abline(v = (subset_size + jitter)[ell_min], lwd=3, col='lightgray', lty=3)
    arrows(jitter + subset_size, pi_loss[,1],
           jitter + subset_size, pi_loss[,2],
           length=0.05, angle=90, code=3, lwd=4)
    if(!is.null(emp_loss)) lines(subset_size + jitter, emp_loss, type='p', lwd=4, col='gray', pch=4, cex = 1.5)
  }

  return(
    list(
      all_accept = all_accept,
      beta_hat_small = beta_hat_small,
      beta_hat_min = beta_hat_min,
      ell_small = ell_small,
      ell_min = ell_min
    )
  )
}
#' Compute the acceptable family of linear subsets for the
#' random intercept model
#'
#' Given output from a Bayesian random intercept model and a candidate of
#' subsets, compute the *acceptable family* of subsets that
#' match or nearly match the predictive accuracy of the "best" subset.
#' The acceptable family may be computed for any set of covariate values
#' \code{XX}; if \code{XX = X} are the in-sample points, then
#' cross-validation is used to assess out-of-sample predictive performance.
#'
#' @param post_y_pred \code{S x m x n} matrix of posterior predictive
#' at the given \code{XX} covariate values for \code{m} replicates per subject
#' @param post_lpd \code{S} evaluations of the log-likelihood computed
#' at each posterior draw of the parameters
#' @param post_sigma_e (\code{nsave}) draws from the posterior distribution
#' of the observation error SD
#' @param post_sigma_u (\code{nsave}) draws from the posterior distribution
#' of the random intercept SD
#' @param XX \code{n x p} matrix of covariates at which to evaluate
#' @param YY \code{m x n} matrix of response variables (optional)
#' @param indicators \code{L x p} matrix of inclusion indicators (booleans)
#' where each row denotes a candidate subset
#' @param post_y_pred_sum (\code{nsave x n}) matrix of the posterior predictive
#' draws summed over the replicates within each subject (optional)
#' @param eps_level probability required to match the predictive
#' performance of the "best" model (up to \code{eta_level})
#' @param eta_level allowable margin (%) between each acceptable model
#' and the "best" model
#' @param K number of cross-validation folds (optional)
#' @param sir_frac fraction of the posterior samples to use for SIR (optional)
#' @param plot logical; if TRUE, include a plot to summarize the predictive
#' performance across candidate subsets
#' @return a list containing the following elements:
#' \itemize{
#' \item \code{all_accept}: indices (i.e., rows of \code{indicators})
#' that correspond to the acceptable subsets
#' \item \code{beta_hat_small} linear coefficients for the
#' smallest acceptable model
#' \item \code{beta_hat_min} linear coefficients for the
#' "best" acceptable model
#' \item \code{ell_small}: index (i.e., row of \code{indicators}) of the
#' smallest acceptable model
#' \item \code{ell_min}: index (i.e., row of \code{indicators}) of the
#' "best" acceptable model
#' }
#' @details When \code{XX = X} is the observed covariate values,
#' then \code{post_lpd} and \code{yy} must be provided. These
#' are used to compute the cross-validated predictive and empirical
#' squared errors; the predictive version relies on a sampling importance-resampling
#' procedure.
#'
#' When \code{XX} corresponds to a new set of covariate values, then set \code{post_lpd = NULL}
#' and \code{yy = NULL} (these are the default values).
#'
#' Additional details on the predictive and empirical comparisons are
#' in \code{pp_loss_randint}.
#'
#' @importFrom graphics abline arrows lines
#' @importFrom stats approxfun binomial coef lm quantile runif
#' @export
accept_family_randint = function(post_y_pred,
                                 post_lpd,
                                 post_sigma_e,
                                 post_sigma_u,
                                 XX,
                                 YY,
                                 indicators,
                                 post_y_pred_sum = NULL,
                                 eps_level = 0.05,
                                 eta_level = 0.00,
                                 K = 10,
                                 sir_frac = 0.5,
                                 plot = TRUE
){

  # Get dimensions of the covariate matrix:
  n = nrow(XX); p = ncol(XX)

  # And some other dimensions
  S = nrow(post_y_pred) # number of posterior simulations
  m = dim(post_y_pred)[2] # number of replicates per subject
  L = nrow(indicators)  # number of subsets to consider

  # Some basic checks:
  if(!is.null(YY) && ncol(YY) !=n)
    stop('the number of columns of YY must equal the number of rows of XX')

  if(dim(post_y_pred)[3] != n)
    stop('incorrect dimensions for post_y_pred')

  if(!is.null(post_lpd) && nrow(post_lpd) != S)
    stop('post_lpd must have the same number of rows (simulations) as post_y_pred')

  if(ncol(indicators) != p)
    stop('indicators must have the same number of columns as XX has rows')

  # This can be slow if n is large:
  if(is.null(post_y_pred_sum)){
    post_y_pred_sum = apply(post_y_pred, c(1,3), sum)
  }

  # Compute the empirical and predictive losses:
  pred_emp_loss = pp_loss_randint(post_y_pred = post_y_pred,
                                  post_lpd = post_lpd,
                                  post_sigma_e = post_sigma_e,
                                  post_sigma_u = post_sigma_u,
                                  XX = XX,
                                  YY = YY,
                                  indicators = indicators,
                                  post_y_pred_sum = post_y_pred_sum,
                                  K = K,
                                  sir_frac = sir_frac)

  # Store these better:
  emp_loss = pred_emp_loss$emp_loss
  pred_loss = pred_emp_loss$pred_loss
  rm(pred_emp_loss)

  # Indices of acceptable subsets:
  all_accept = which(colMeans(pred_loss <= eta_level)
                     >= eps_level)

  # Index of "best" subset, if possible:
  if(!is.null(emp_loss)){
    ell_min = which.min(emp_loss)
  } else ell_min = NULL

  # Subset sizes:
  subset_size = rowSums(indicators)

  # Minimize size of acceptable subset:
  min_size_accept = min(subset_size[all_accept])

  # Index of acceptable subsets with this size:
  ind_min_size_accept = which(subset_size[all_accept] == min_size_accept)

  # If more than one, select the minimum empirical loss (or expected predictive loss)
  if(length(ind_min_size_accept) > 1){
    if(!is.null(emp_loss)){
      # use empirical loss, if available
      ell_small = all_accept[ind_min_size_accept][which.min(emp_loss[all_accept[ind_min_size_accept]])]
    } else {
      # otherwise, use expected predictive loss
      ell_small = all_accept[ind_min_size_accept][which.min(colMeans(pred_loss)[all_accept[ind_min_size_accept]])]
    }
  } else ell_small = all_accept[ind_min_size_accept]
  if(is.infinite(ell_small)) ell_small = ell_min # in case there is none?

  # Compute the coefficients for each index:
  objXY = getXY_randint(XX = XX,
                        post_y_pred = post_y_pred,
                        post_sigma_e = post_sigma_e,
                        post_sigma_u = post_sigma_u,
                        post_y_pred_sum = post_y_pred_sum)
  X_star = objXY$X_star; y_star = objXY$y_star; rm(objXY)

  # coefficients for "best" subset
  if(!is.null(ell_min)){
    beta_hat_min = numeric(p)
    beta_hat_min[indicators[ell_min,]] =
      coef(lm(y_star ~ X_star[, indicators[ell_min,]] - 1))
  } else beta_hat_min = NULL

  # coefficients for smallest acceptable subset:
  beta_hat_small = numeric(p)
  beta_hat_small[indicators[ell_small,]] =
    coef(lm(y_star ~ X_star[, indicators[ell_small,]] - 1))

  # Now add the plot, if desired:
  if(plot){
    # 100(1 - 2*eps_level)% prediction interval for the loss
    pi_loss = t(apply(pred_loss, 2,
                      quantile, c(eps_level, 1 - eps_level)))

    jitter = runif(L, min = -0.25, max = 0.25)
    plot(subset_size + jitter,
         colMeans(pred_loss), type='p', ylim = range(pi_loss, 0), lwd=4,
         xlab = 'Subset Size',
         ylab = 'Difference in loss (%)',
         main = paste('Difference in ',K,'-fold Mahalanobis loss (%)', sep=''), #)
         cex.main = 2, cex.lab = 2, cex.axis = 2)
    abline(v = 1:p, col='gray')
    abline(h = eta_level, lwd=2, lty=6)
    abline(v = (subset_size + jitter)[ell_small], lwd=8, col='darkgray')
    if(!is.null(ell_min)) abline(v = (subset_size + jitter)[ell_min], lwd=8, col='lightgray', lty=3)
    arrows(jitter + subset_size, pi_loss[,1],
           jitter + subset_size, pi_loss[,2],
           length=0.15, angle=90, code=3, lwd=8)
    if(!is.null(emp_loss)) lines(subset_size + jitter, emp_loss, type='p', lwd=8, col='gray', pch=4, cex = 1.5)
  }

  return(
    list(
      all_accept = all_accept,
      beta_hat_small = beta_hat_small,
      beta_hat_min = beta_hat_min,
      ell_small = ell_small,
      ell_min = ell_min
    )
  )
}
#' Variable importance for the acceptable family
#'
#' Given the candidate subsets and the indicators of the acceptable family
#' of subsets, compute for each variable the proportion of acceptable subsets in
#' which that variable appears. If specified, variable co-appearances
#' can be computed as reported as well.
#' @param indicators \code{L x p} matrix of inclusion indicators (booleans)
#' where each row denotes a candidate subset
#' @param all_accept indices (i.e., rows of \code{indicators})
#' that correspond to the acceptable subsets
#' @param co logical; if TRUE, compute and return the co-variable importances
#' @param xnames the names of the x-variables
#' @return a list with the variable importances \code{vi_inc} and the co-variable
#' importances \code{vi_co}
#' @export
var_imp = function(indicators, all_accept, co = TRUE, xnames = NULL){

  # Number of covariates:
  p = ncol(indicators)

  if(length(all_accept) > nrow(indicators))
    stop('Cannot have more acceptable subsets than candidate subsets!')

  if(!is.null(xnames)){
    if(length(xnames) != p)
      stop('xnames must have length p')
  } else xnames = paste(1:p)

  # For each variable, how often is it included in an acceptable subset?
  vi_inc = colMeans(indicators[all_accept,]); names(vi_inc) = xnames

  if(co){
    vi_co = diag(vi_inc, p)
    colnames(vi_co) = rownames(vi_co) = xnames
    for(j1 in 1:(p-1)){
      for(j2 in (j1 + 1):p){
        vi_co[j1,j2] = vi_co[j2,j1] =
          mean(indicators[all_accept, j1] & indicators[all_accept, j2])
      }
    }

  } else vi_co = NULL

  return(
    list(
      vi_inc = vi_inc,
      vi_co = vi_co
    )
  )
}
#' (Adaptive) lasso for Bayesian variable selection
#'
#' Given output from a Bayesian model, compute an (adaptive)
#' lasso search path to provide linear variable selection.
#' This function only indicates the variables selected along the path.
#'
#' @param yy_hat vector of fitted values
#' @param XX matrix of covariates
#' @param wts vector of observation weights (for weighted least squares)
#' @param ww_hat vector of weights for the adaptive lasso
#' @return \code{inclusion_index}: the matrix of inclusion indicators (columns) for
#' each subset returned (rows)
#' @details If an intercept is in the \code{XX} matrix,
#' it is assumed that (i) the intercept is the first column and
#' (ii) the intercept should be included in each subset.
#' @import glmnet
#' @export
lasso_path = function(yy_hat,
                      XX,
                      wts = NULL,
                      ww_hat = NULL
){

  # Remove the intercept, if there is one:
  if(all(XX[,1] == 1)) XX = XX[,-1]

  # Get dimensions of the covariate matrix:
  n = nrow(XX); p = ncol(XX)

  # Some basic checks:
  if(length(yy_hat) !=n)
    stop('length of yy_hat must equal the number of rows of XX')

  if(!is.null(ww_hat)){
    if(length(ww_hat) != p)
      stop('length of ww_hat should be p - 1 (excludes the intercept)')
  } else ww_hat = rep(1, p)

  if(!is.null(wts)){
    if(length(wts) != n)
      stop('wts must have length n')
  } else wts = rep(1,n)

  # Fit the lasso:
  fit_lasso = glmnet(x = XX,
                     y = yy_hat,
                     weights = wts,
                     penalty.factor = ww_hat)

  # Store the estimates:
  beta_path = rbind(fit_lasso$a0,  # intercept
                    fit_lasso$beta) # coefs

  # And the lambda values:
  lambda_path = c(fit_lasso$lambda,
                  0) # Add zero for no penalty

  # Subset sizes:
  subset_size = c(apply(beta_path, 2, function(a) sum(a != 0)),
                  p+1) # Add full subset

  # Index of best subset for each size:
  size_ind = sapply(unique(subset_size), function(mi) max(which((subset_size == mi))))

  # Update the path to approximate (almost) all possible subset sizes:
  lambda_path = unique(approxfun(subset_size[size_ind], lambda_path[size_ind])(1:(p+1)))

  # Add some lambda points between zero and one:
  lambda_path = sort(unique(c(lambda_path, seq(1,0, length.out = 100))),decreasing = TRUE)

  # And re-fit:
  fit_lasso = glmnet(x = XX,
                     y = yy_hat,
                     weights = wts,
                     penalty.factor = ww_hat,
                     lambda = lambda_path)

  # Update relevant parameters:
  beta_path = rbind(fit_lasso$a0,  # intercept
                    fit_lasso$beta) # coefs

  lambda_path = fit_lasso$lambda
  subset_size = apply(beta_path, 2, function(a) sum(a != 0))

  # Index of best subset for each size:
  size_ind = sapply(unique(subset_size), function(mi) max(which((subset_size == mi))))

  # Subset to subsets of these unique sizes:
  beta_path = beta_path[,size_ind];
  lambda_path = lambda_path[size_ind];
  subset_size = subset_size[size_ind]

  inclusion_index = t(as.matrix(beta_path != 0))

  return(inclusion_index)
}
#' Projected predictive distribution for regression coefficients
#'
#' Given draws from the predictive distribution, project these
#' draws onto (a  subset of) the covariates. This produces
#' many predictive draws for the regression coefficients,
#' which provides uncertainty quantification.
#' @param post_y_pred \code{S x n} matrix of posterior predictive
#' at the given \code{XX} covariate values
#' @param XX \code{n x p} matrix of covariates
#' @param sub_x vector of inclusion indicators for the \code{p} covariates;
#' the remaining coefficients will be fixed at zero
#' @param use_ols logical; if TRUE, use ordinary least squares regression (default);
#' otherwise use logistic regression
#' @return \code{post_beta}: the \code{S x p} matrix of
#' draws from the projected predictive distribution of  the regression coefficients.
#' @export
proj_posterior = function(post_y_pred, XX, sub_x = 1:ncol(XX), use_ols = TRUE){

  # Get dimensions of the covariate matrix:
  n = nrow(XX); p = ncol(XX)
  S = nrow(post_y_pred) # number of posterior simulations

  if(ncol(post_y_pred) != n)
    stop('number of columns of post_y_pred must equal the number of rows of XX')

  if(length(sub_x) > p)
    stop('length of sub_x must be less than or equal to p')

  # Subset the columns of XX:
  XX_nz = XX[, sub_x]

  # Storage: include zeros as needed!
  post_beta = array(0, c(S, p))
  if(use_ols){
    post_beta[,sub_x] = tcrossprod(post_y_pred, t(XX_nz)) %*%
      chol2inv(chol(crossprod(XX_nz)))
  } else{
    post_beta[,sub_x] = t(sapply(1:S, function(s){
      suppressWarnings(
        coef(
          glm(post_y_pred[s,] ~ XX_nz - 1, family = binomial())
        )
      )
    }))
  }

  colnames(post_beta) = colnames(XX)

  return(post_beta)
}
#' Projected predictive distribution for regression coefficients
#' in the random intercept model
#'
#' Given draws from the predictive distribution of the random intercept model, project these
#' draws onto (a  subset of) the covariates using Mahalanobis loss. This produces
#' many predictive draws for the regression coefficients,
#' which provides uncertainty quantification.
#' @param post_y_pred \code{S x m x n} matrix of posterior predictive
#' at the given \code{XX} covariate values
#' @param XX \code{n x p} matrix of covariates
#' @param sub_x vector of inclusion indicators for the \code{p} covariates;
#' the remaining coefficients will be fixed at zero
#' @param post_sigma_e (\code{nsave}) draws from the posterior distribution
#' of the observation error SD
#' @param post_sigma_u (\code{nsave}) draws from the posterior distribution
#' of the random intercept SD
#' @param post_y_pred_sum (\code{nsave x n}) matrix of the posterior predictive
#' draws summed over the replicates within each subject (optional)
#' @return \code{post_beta}: the \code{S x p} matrix of
#' draws from the projected predictive distribution of  the regression coefficients.
#' @import Matrix
#' @export
proj_posterior_randint = function(post_y_pred, XX, sub_x = 1:ncol(XX),
                                  post_sigma_e, post_sigma_u,
                                  post_y_pred_sum = NULL){

  # Get dimensions of the covariate matrix:
  n = nrow(XX); p = ncol(XX)

  # And some other dimensions
  S = nrow(post_y_pred) # number of posterior simulations
  m = dim(post_y_pred)[2] # number of replicates per subject

  if(dim(post_y_pred)[3] != n)
    stop('incorrect dimensions for post_y_pred')

  if(length(sub_x) > p)
    stop('length of sub_x must be less than or equal to p')

  # This can be slow if n is large:
  if(is.null(post_y_pred_sum)){
    post_y_pred_sum = apply(post_y_pred, c(1,3), sum)
  }

  # Recurring term:
  m_scale_hat = mean(1/(post_sigma_e^2/post_sigma_u^2 + m))

  # Subset the columns of XX:
  XX_nz = XX[, sub_x]

  # Storage: include zeros as needed!
  post_beta = array(0, c(S, p))

  # Estimated Mahalanobis weight matrix (block diagonal), ignoring sigma_e^2:
  Omega_hat = diag(1, m) - m_scale_hat

  # Matrix square-root, then blocked:
  sqrt_Omega_hat = chol(Omega_hat)
  sqrt_Omega_block_hat = bdiag(lapply(1:n, function(i){sqrt_Omega_hat}))

  # Stacked X-matrix:
  X_stack = apply(XX_nz, 2, function(x)
    matrix(rep(x, each = m), nrow = m))

  # Design matrix:
  X_star = as.matrix(sqrt_Omega_block_hat%*%X_stack)

  post_beta[,sub_x] = t(sapply(1:S, function(s){
    # Response matrix, ignoring sigma_e^2::
    Omega_y_s = post_y_pred[s,,] - rep(m_scale_hat*post_y_pred_sum[s,], each = m)
    y_star_s = matrix(crossprod(Matrix::solve(sqrt_Omega_hat),
                                   Omega_y_s))
    # And the coefficients:
    coef(lm(y_star_s ~ X_star-1))
  }))

  return(post_beta)
}
#' Compute the optimal linear coefficients for any covariates
#'
#' Given the fitted values from a Bayesian model and (a subset of)
#' covariates of interest, compute the optimal linear
#' coefficients. This applies for continuous (\code{use_ols = TRUE})
#' or binary (\code{use_ols = FALSE}) outcomes.
#'
#' @param y_hat \code{n} vector of fitted values
#' at the given \code{XX} covariate values
#' @param XX \code{n x p} matrix of covariates (e.g., restricted to the
#' subset of \code{p} variables of interest)
#' @param use_ols logical; if TRUE, use ordinary least squares regression (default);
#' otherwise use logistic regression
#' @return the \code{p} optimal regression coefficients
#'
#' @export
get_coefs = function(y_hat, XX, use_ols = TRUE){

  # Get dimensions of the covariate matrix:
  n = nrow(XX); p = ncol(XX)

  if(length(y_hat) != n)
    stop('length of y_hat must equal the number of rows of XX')

  # Coefficients:
  if(use_ols){
    beta_hat = coef(lm(y_hat ~ XX - 1))
  } else {
    beta_hat = suppressWarnings(
      coef(glm(y_hat ~ XX - 1, family = binomial()))
    )
  }

  # To clarify the output:
  names(beta_hat) = colnames(XX)

  return(beta_hat)
}
#' Marginal pre-screening algorithm
#'
#' Given \code{S} draws of the \code{p} regression coefficients,
#' this function (i) computes (SD-standarized) posterior means,
#' (ii) orders the absolute values of this statistic, and
#' (iii) returns the indices of the largest \code{num_to_keep} values.
#'
#' @param post_beta (\code{S x p}) matrix of posterior simulations
#' of the \code{p} regression coefficients
#' @param num_to_keep number of variables to return in the screening process
#' @return \code{ind_keep}: the indices of the variables to keep
#' @details If the Bayesian model is nonlinear, \code{post_beta} may be
#' obtained instead by projecting each vector of posterior predictive draws
#' onto the covariate matrix.
#' @importFrom stats sd
#' @export
prescreen = function(post_beta, num_to_keep){

  # Posterior mean of coefficients:
  beta_hat = colMeans(post_beta)

  # Posterior SD of coefficients:
  beta_sd = apply(post_beta, 2, sd)

  # Standardized coefficients:
  beta_hat_std = beta_hat/beta_sd

  # Order the coefficients
  ind_beta_rank = order(abs(beta_hat_std), decreasing = TRUE)

  # And keep the top num_to_keep:
  ind_keep = ind_beta_rank[1:num_to_keep]

  return(ind_keep)
}
#' Marginal pre-screening algorithm given lasso output
#'
#' Given the estimated lasso coefficients, return the
#' indices of the first \code{num_to_keep} variables that enter the model
#' @param beta_hat_lasso (\code{p x L}) matrix of estimated coefficients
#' from the lasso; sparsity decreases in the column index \code{L}
#' @param num_to_keep number of variables to return in the screening process
#' @return \code{ind_keep}: the indices of the variables to keep
#' @import glmnet
#' @export
prescreen_lasso = function(beta_hat_lasso, num_to_keep){

  # Number of nonzero predictors for each row:
  num_nz = colSums(beta_hat_lasso!=0)

  # Column index of the largest model of size num_to_keep:
  col_ind = max(which((num_nz <= num_to_keep)))

  # Indices of nonzero coefficients for this model:
  ind_keep = which(beta_hat_lasso[,col_ind] != 0)

  return(ind_keep)
}
