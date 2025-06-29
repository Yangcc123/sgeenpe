% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sgeenpe.R
\name{sgeenpe}
\alias{sgeenpe}
\title{Title Estimate the parameters using the estimated correlation matrix}
\usage{
sgeenpe(formula, data, eps = 1e-06, itermax = 100)
}
\arguments{
\item{formula}{a formula expression, of the form \code{response ~ predictors}. The \code{response} is a \code{Surv} object with right censoring.}

\item{data}{a data frame in which to interpret the variables named in the \code{formula}.}

\item{eps}{tolerance for convergence. The default is \code{eps = 1e-6}}

\item{itermax}{specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after \code{itermax} iterations and
the estimates will be based on the last iteration. The default \code{itermax = 100}.}
}
\value{
A dataframe containing information related to estimated parameters
}
\description{
Title Estimate the parameters using the estimated correlation matrix
}
