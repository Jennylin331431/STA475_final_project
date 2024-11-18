### Helper functions ###


#Predict Survival Probabilities from a Weibull Model
#'
#' This function calculates survival probabilities for given covariate values and time points
#' based on a fitted Weibull model.
#'
#' @param object A fitted Weibull model object from survreg()
#' @param newdata A data frame containing covariate values for prediction
#'           (must have the same column names as the data used to fit the model)
#'           If the model was fitted with a formula, the data frame should contain all variables in the formula.
#'           If the model was fitted with no predictors, the data frame should contain a dummy variable taking the values NA
#' @param t A numeric vector of time points at which to predict survival probabilities
#'
#' @return A numeric vector of predicted survival probabilities corresponding to each time point in t, based on the fitted model in `object`
#'
#' @details
#' The function uses the fitted Weibull model to predict survival probabilities.
#' It calculates the linear predictor, then the cumulative hazard, and finally the survival probability.
#'
#' @examples
#' # Assuming 'fit' is a fitted Weibull model and 'new_data' contains new covariate values:
#' predict_survival_weibull(fit, new_data, t = c(1, 5, 10))
predict_survival_weibull <- function(object, # fitted Weibull model from survreg
                                     newdata, # covariate values for which to predict survival
                                     t){ # time points to predict survival at
  # mu_hat = beta_hat * z
  mu_hat <- predict(object, newdata=newdata, type="link")
  # H(t) = -logS(t) = exp((logt - beta_hat*z)/tau_hat) = (t/exp(mu_hat))^(1/tau_hat)
  cum_hazard <- (t / exp(mu_hat))^(1/object$scale)
  # S(t) = exp(-H(t))
  surv <- exp(-cum_hazard)
  return(surv)
}


#Predict Survival Probabilities from an exponential model
#'
#' This function calculates survival probabilities for given covariate values and time points
#' based on a fitted exponential model.
#'
#' @param object A fitted exponential model object from survreg()
#' @param newdata A data frame containing covariate values for prediction
#'           (must have the same column names as the data used to fit the model)
#'           If the model was fitted with a formula, the data frame should contain all variables in the formula.
#'           If the model was fitted with no predictors, the data frame should contain a dummy variable taking the values NA
#' @param t A numeric vector of time points at which to predict survival probabilities
#'
#' @return A numeric vector of predicted survival probabilities corresponding to each time point in t, based on the fitted model in `object`
#'
#' @details
#' The function uses the fitted exponential model to predict survival probabilities.
#' It calculates the linear predictor, then the cumulative hazard, and finally the survival probability.
#'
#' @examples
#' # Assuming 'fit' is a fitted exponential model and 'new_data' contains new covariate values:
#' predict_survival_exponential(fit, new_data, t = c(1, 5, 10))
predict_survival_exponential <- function(object, # fitted exponential model from survreg
                                         newdata, # covariate values for which to predict survival
                                         t){ # time points to predict survival at
  # mu_hat = beta_hat * z
  mu_hat <- predict(object, newdata=newdata, type="link")
  # H(t) = -log S(t) = exp(log(t) - beta_hat*z) = t * exp(-beta_hat*z) = t / exp(beta_hat*z)
  cum_hazard <- (t / exp(mu_hat))
  # S(t) = exp(-H(t))
  surv <- exp(-cum_hazard)
  return(surv)
}

predict_survival_loglogistic <- function(object, # fitted loglogistic model from survreg
                                         newdata, # covariate values for which to predict survival
                                         t){ # time points to predict survival at
  # mu_hat = beta_hat * z
  mu_hat <- predict(object, newdata=newdata, type="link")
  tau_hat <- object$scale
  # S(t) = 1 / (1 + exp((log(t) - mu_hat)/tau_hat))
  surv <- 1 / (1 + exp((log(t) - mu_hat)/tau_hat))
  return(surv)
}


#' Predict Survival Probabilities from a Log-Logistic Model
#'
#' This function calculates survival probabilities for given covariate values and time points
#' based on a fitted log-logistic model.
#'
#' @param object A fitted log-logistic model object from survreg()
#' @param newdata A data frame containing covariate values for prediction
#'           (must have the same column names as the data used to fit the model)
#'           If the model was fitted with a formula, the data frame should contain all variables in the formula.
#'           If the model was fitted with no predictors, the data frame should contain a dummy variable taking the values NA
#' @param t A numeric vector of time points at which to predict survival probabilities
#'
#' @return A numeric vector of predicted survival probabilities corresponding to each time point in t
#'
#' @details
#' The function uses the fitted log-logistic model to predict survival probabilities.
#' It calculates the linear predictor (mu_hat) and uses the scale parameter (tau_hat)
#' to compute the survival probability directly.
#'
#' @examples
#' # Assuming 'fit' is a fitted log-logistic model and 'new_data' contains new covariate values:
#' predict_survival_loglogistic(fit, new_data, t = c(1, 5, 10))
predict_survival_lognormal <- function(object, # fitted lognormal model from survreg
                                       newdata, # covariate values for which to predict survival
                                       t){ # time points to predict survival at
  # mu_hat = beta_hat * z
  mu_hat <- predict(object, newdata=newdata, type="link")
  tau_hat <- object$scale
  # S(t) = 1 - Phi((log(t) - mu_hat)/tau_hat)
  surv <- 1 - pnorm((log(t) - mu_hat)/tau_hat)
  return(surv)
}
