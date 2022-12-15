\name{fracridge-package}
\alias{fracridge-package}
\alias{fracridge}
\docType{package}
\title{
  A Fast and Interpretable Reparameterization of Ridge Regression Reimplemented in R
}
\description{
  A reimplementation of fast and interpretable Reparameterization of ridge regression in R.
}
\details{
  We reimplement a fast and interpretable Reparameterization of ridge regression (Fractional Ridge Regression) in R. The package includes two most important functions frac_ridge() and frac_ridge_pred(). The former function is to obtain the estimated coefficients for ridge regression given user's input fraction, design matrix X, and data y, and the latter function is to perform prediction given the estimated coefficients and the design matrix X.
}
\author{
Zhiqi Zhou, <zhiqizh@umich.edu>
Yeming Lin, <umymlin@umich.edu>
Jiale Tan, <jialetan@umich.edu>

Maintainer: Zhiqi Zhou <zhiqizh@umich.edu>
}
\references{
Rokem, A., & Kay, K. (2020). Fractional ridge regression: a fast, interpretable reparameterization of ridge regression (https://academic.oup.com/gigascience/article/9/12/giaa133/6011381).
}
\keyword{ package }
\seealso{
  https://github.com/zhiqizhouu/fracridge_r
}
\examples{
  set.seed(615)
  X = matrix(rnorm(100 * 2), nrow = 100, ncol = 2)
  beta = matrix(rnorm(2 * 1), nrow = 2, ncol = 1)
  y = X \%*\% beta
  sd_y = apply(y, 2, sd)
  noise = matrix(rnorm(100 * 1, mean = 0, sd = sd_y), nrow = 100, ncol = 1)
  y = y + noise   
  fracs = c(0.1, 0.2, 0.3)
  beta_ridge = frac_ridge(X, y, fracs)$coefs
  y_pred = frac_ridge_pred(X, beta_ridge)
}
