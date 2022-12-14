% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FractionalRidgeMain.R
\name{frac_ridge}
\alias{frac_ridge}
\title{A fast and interpretable reparameterization of ridge regression}
\usage{
frac_ridge(X, y, fracs = NULL, intercept = FALSE, scale_x = FALSE)
}
\arguments{
\item{X}{A size n*p matrix (design matrix for regression, with n number of observations and p number of variables)}

\item{y}{A size n*t matrix (data matrix, with n number of observations and t number of targets)}

\item{fracs}{A float or a size f vector (user supplied fractions of L2-norm of parameters, relative to OLS solution, ranging from 0 to 1, and it cannot be set = 0)}

\item{intercept}{A logical argument indicates whether the user needs an intercept or not}

\item{scale_x}{A logical argument indicates whether the user wants to scale the data matrix X}
}
\value{
A list containing coef (estimated parameters for every desired fraction, a 3D array with shape (p, f, t)) and alphas (the alpha values associated with each solution, shape (f, p))
}
\description{
A fast and interpretable reparameterization of ridge regression
}
\examples{
x = matrix(rnorm(100*10), 100, 10)
y = rnorm(100)
fracs = c(0.1, 0.2, 0.3)
res = frac_ridge(x, y, fracs, intercept = TRUE, scale_x = FALSE)
}
