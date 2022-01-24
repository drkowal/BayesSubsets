# BayesSubsets
Bayesian Subset Selection and Variable Importance

Given any Bayesian model for prediction or classification,
    these tools will (i) search for the "best" subsets of linear predictors;
    (ii) collect the *acceptable family* of near-optimal subsets of linear predictors
    that match or nearly-match the predictive performance of the "best" model (according to 
    predictive cross-validation); (iii) summarize the acceptable family using the "best" model,
    the smallest acceptable model, and customized variable importance metrics; and (iv) provide
    predictive uncertainty quantification for the linear coefficients. In effect, these methods
    provide interpretable (linear) summaries of the Bayesian model. The strategy of collecting 
    the acceptable family of subsets substantially reduces the inherent sensitivity in selecting 
    a single "best" model and recognizes that many different subsets may perform quite well. 
    For any subset of covariates, the estimated linear coefficients are optimal in a Bayesian 
    decision analysis sense. 
