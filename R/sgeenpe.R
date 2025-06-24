#' Title Estimate the parameters using the estimated correlation matrix
#'
#' @param formula a formula expression, of the form \code{response ~ predictors}. The \code{response} is a \code{Surv} object with right censoring.
#' @param data a data frame in which to interpret the variables named in the \code{formula}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}
#' @param itermax specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{itermax} iterations and
#' the estimates will be based on the last iteration. The default \code{itermax = 100}.
#'
#' @return A dataframe containing information related to estimated parameters
#' @export
#' @importFrom survival Surv
sgeenpe = function(formula,data,eps = 1e-06,itermax =100){
  message("Parameter estimation of the independent working matrix")
  dataframe_indepence = sgeecorr(formula = formula ,data = data,corrstr = "independence",corr_matix=NULL,eps = eps,itermax =itermax)
  message("Use the estimated working matrix")
  matrix_est = geemat(formula = formula ,data = data,beta_Lambda=dataframe_indepence)
  print(matrix_est)
  message("Parameter estimation using the estimated working matrix")
  dataframe_QU = sgeecorr(formula = formula,data = data,corrstr = "def",corr_matix = matrix_est,eps = eps,itermax =itermax)
  return(dataframe_QU)
}

