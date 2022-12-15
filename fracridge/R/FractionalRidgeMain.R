#' A singular value decomposition of matrix 
#' @param X A size n*p matrix (design matrix for regression, with n number of observations and p number of variables)
#' @param y A size n*t matrix (data matrix, with n number of observations and t number of targets)
#' @return A list containing singular values of X (lambda_val), unitary matrix vt from SVD (vt), OLS coefficients in the rotated space (beta_ols), and y in the rotated space (y_trsf)
#' @examples
#' x = matrix(rnorm(100*10), 100, 10)
#' y = rnorm(100)
#' res = svd_x(x, y)
#' @export 
svd_x = function(X, y){
  if (is.null(ncol(y))){
    y = as.matrix(y)
  }
  if (nrow(X) > ncol(X)){
    res = svd(t(X) %*% X)
    uu = res$u
    lambda_sq = res$d
    vt = t(res$v)
    lambdas = sqrt(lambda_sq)
    if (ncol(y) >= nrow(X)){
      y_trsf = (diag(1/lambdas) %*% vt %*% t(X)) %*% y
    }else{
      y_trsf = diag(1/lambdas) %*% vt %*% (t(X) %*% y)
    }
  }else{
    res = svd(X)
    uu = res$u
    lambdas = res$d
    vt = t(res$v)
    y_trsf = t(uu) %*% y
  }
  beta_ols = y_trsf/lambdas
  return(list(lambda_val = lambdas, vt = vt, beta_ols = beta_ols, y_trsf = y_trsf))
} 

#' A fast and interpretable reparameterization of ridge regression
#' @param X A size n*p matrix (design matrix for regression, with n number of observations and p number of variables)
#' @param y A size n*t matrix (data matrix, with n number of observations and t number of targets)
#' @param fracs A float or a size f vector (user supplied fractions of L2-norm of parameters, relative to OLS solution, ranging from 0 to 1, and it cannot be set = 0)
#' @param intercept A logical argument indicates whether the user needs an intercept or not
#' @param scale_x A logical argument indicates whether the user wants to scale the data matrix X
#' @return A list containing coef (estimated parameters for every desired fraction, a 3D array with shape (p, f, t)) and alphas (the alpha values associated with each solution, shape (f, p))
#' @examples 
#' x = matrix(rnorm(100*10), 100, 10)
#' y = rnorm(100)
#' fracs = c(0.1, 0.2, 0.3)
#' res = frac_ridge(x, y, fracs, intercept = TRUE, scale_x = FALSE)
#' @importFrom stats approxfun
#' @importFrom Rcpp evalCpp
#' @export
frac_ridge = function(X, y, fracs = NULL, intercept = FALSE, scale_x = FALSE){
  if (is.null(ncol(y))){
    y = as.matrix(y)
  }
  
  # Check if the user needs to scale a matrix
  if (scale_x == TRUE){
    X = scaleMat(X)
  }
  
  # Check if the user needs an intercept: if intercept == TRUE, append a column of 1s at the first column of X
  if (intercept == TRUE){
    X = cbind(as.vector(rep(1, nrow(X))), X)
  }
  
  # Check if the vector of fracs is null
  # if it is null, we will define a vector of fractions for the user
  if (is.null(fracs)){
    fracs = seq(0.01, 1, 0.05)
  }
  
  # do singular value decomposition of x
  res_svd = svd_x(X, y)
  # Obtain singular values and transformed ols coefs
  lambda = res_svd$lambda_val
  ols_coef_trsf = res_svd$beta_ols
  # Get the minimum and maximum singular value (which helps us to calculate the range for candidate alphas)
  lambda_p = min(lambda)
  lambda_1 = max(lambda)
  
  # Define a range of candidate alphas
  upper = 10^3*lambda_1^2
  lower = 10^(-3)*lambda_p^2
  alpha = 10^(seq(floor(log10(lower)), ceiling(log10(upper)), by = 0.2))
  alpha = c(0, alpha)
  l = length(alpha)

  # A function to calculate the scaling factor applied to coefficients in the rotated space
  # which is lambda^2/(lambda^2+alpha), where lambda are the singular values
  SF <- function(alpha,lambda){
    p = length(lambda)
    l = length(alpha)
    alpha_grid = matrix(rep(alpha, p), p, l, byrow = TRUE)
    sf_mat = sweep(alpha_grid, 1, lambda^2, "+")  # per row
    sf_mat = lambda^2/sf_mat
    return(sf_mat)
  }
  
  # Get a matrix of scaling factor based on the defined alpha and singular values
  sf_mat = SF(alpha, lambda)
  
  # Prelocate the solution
  bb = dim(y)[2]
  ff = length(fracs)
  nn = dim(X)[1]
  pp = dim(X)[2]
  if (nn >= pp){
    first_dim = pp
  }else{
    first_dim = nn
  }
  
  # Initialize coefs and alphas array
  coef = array(0, dim = c(first_dim, ff, bb))
  alphas = array(0, dim = c(ff, bb))
  # The main loop is over targets
  for (ii in 1:ncol(y)){
    # Apply the scaling factors per alpha
    beta_rr_norm = sqrt(t(sf_mat)^2 %*% ols_coef_trsf[,ii]^2)
    
    # Normalize to the length of the ols solution
    # beta_rr_norm[1] is equal to ols norm
    gammas = beta_rr_norm/beta_rr_norm[1]
    
    # Do interpolation using a range of alpha and its corresponding gammas
    # The interpolation is performed in a log transformed space to avoid numerical issues
    interp_f = approxfun(rev(gammas), rev(log(1 + alpha)))
    
    # Interpolated alpha given desired fractions
    target_alphas = interp_f(fracs)
    # Undo the log transform from the previous step
    new_alphas = exp(target_alphas) - 1
    # Place the alphas for this target
    alphas[,ii] = new_alphas
    
    # Calculate the new scaling factor using the interpolated alphas
    sf_mat_final = SF(new_alphas, lambda)
    
    # Use new scaling factors to calculate the coefficients of beta_ridge in the rotated space
    beta_ridge_rt = sf_mat_final * ols_coef_trsf[,ii]
    coef[,,ii] = beta_ridge_rt
  }
  
  # After iterating over all targets, we first reshape the coef array into a matrix in order to make the unrotation work well
  dim(coef) = c(first_dim, ff * bb)
  # Unrotate the coef using the unitary v matrix and reshape to the desired 3D array output
  beta_ridge = t(res_svd$vt) %*% coef

  dim(beta_ridge) = c(pp, ff, bb)
  
  return(list(coefs = beta_ridge, alphas = alphas))
}

#' A function to do the prediction using X (design matrix) and beta_ridge (estimated coefficients)
#' @param X A size n*p matrix (design matrix for regression, with n number of observations and p number of variables)
#' @param beta_ridge A size (p, f, t) array (f - the number of fractions in previous frac_ridge input) obtained from frac_ridge(...)$coefs
#' @param intercept A logical argument indicates whether the user needs an intercept or not
#' @return A size (p, t, f) array with each array[,,i] representing the estimated y for each user-defined fraction in frac_ridge(...)
#' @examples
#' x = matrix(rnorm(100*10), 100, 10)
#' y = rnorm(100)
#' beta_ridge = frac_ridge(x, y, 0.1, intercept = TRUE, scale_x = FALSE)$coefs
#' y_pred = frac_ridge_pred(x, beta_ridge, TRUE)
#' @export
frac_ridge_pred <- function(X, beta_ridge, intercept = FALSE){
  if (intercept == TRUE){
    X = cbind(as.vector(rep(1, nrow(X))), X)
  }
  num_frac = dim(beta_ridge)[2]
  num_target = dim(beta_ridge)[3]
  pred = array(0, dim = c(nrow(X), num_frac, num_target))
  for(i in 1:num_target){
    pred[,,i] = X %*% beta_ridge[,,i]
  } 
  # reshape the array of y_pred 
  pred = aperm(pred, c(1,3,2))
  return(pred)
}

#' A function to calculate Rsquared score using frac_ridge output
#' @param y A size n*t matrix (data matrix, with n number of observations and t number of targets)
#' @param X A size n*p matrix (design matrix for regression, with n number of observations and p number of variables)
#' @param beta_ridge A size (p, f, t) array (f - the number of fractions in previous frac_ridge input) obtained from frac_ridge(...)$coefs
#' @param intercept A logical argument indicates whether the user needs an intercept or not
#' @return r2_score A size f vector consisting of average Rsquared score for each fraction
#' @examples
#' x = matrix(rnorm(100*10), 100, 10)
#' y = rnorm(100)
#' beta_ridge = frac_ridge(x, y, 0.1, intercept = TRUE, scale_x = FALSE)$coefs
#' r2 = r2_score(y, x, beta_ridge, TRUE)
#' @export
r2_score <- function(y, X, beta_ridge, intercept = FALSE){
  if (is.null(ncol(y))){
    y = as.matrix(y)
  }
  y_pred = frac_ridge_pred(X, beta_ridge, intercept)
  num_frac = dim(y_pred)[3]
  mean_r2 = rep(0, num_frac)
  for (i in 1:num_frac){
    numer = colSums((y - y_pred[,,i])^2)
    yybar = sweep(y, 2, colMeans(y))
    denom = colSums(yybar^2)
    r2 = 1 - numer/denom
    mean_r2[i] = mean(r2)
  }
  return(mean_r2)
}
