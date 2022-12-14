#' A function to generate simulated data with n rows, p predictors, and t targets
#' @param n - rows; p - number of features; t - number of targets (number of columns of y)
gen_large_sample_data = function(n, p = 5, t = 1){
  set.seed(615)
  X = matrix(rnorm(n * p), nrow = n, ncol = p)
  beta = matrix(rnorm(p * t), nrow = p, ncol = t)
  y = X %*% beta
  sd_y = apply(y, 2, sd)
  noise = matrix(rnorm(n * t, mean = 0, sd = sd_y), nrow = n, ncol = t)
  final_y = y + noise
  return(list(X = X, y = final_y))
}

#' A naive implementation of standard ridge regression using closed form solution in which the inverse of matrix is computed for every alpha value
standard_ridge_reg = function(X, y, alpha){
  std_br = array(0, dim = c(ncol(X), ncol(y), length(alpha)))
  # if the user's input alpha is a vector, we have to iterate over each alpha
  for (i in 1:length(alpha)){
    std_br[,,i] = (solve(t(X) %*% X + diag(alpha[i], ncol(X), ncol(X)))) %*% t(X) %*% y 
  }
  if (ncol(y) == 1){
    std_br = matrix(std_br, nrow = ncol(X), ncol = length(alpha))
  }
  return(std_br)
}

#' An improved implementation of standard ridge regression using singular value decomposition of the design matrix X
#' We can avoid taking the inverse of a matrix using this method
svd_ridge_reg = function(X, y, alpha){
  svd_res = svd_x(X, y)
  bb = dim(y)[2]
  ff = length(alpha)
  nn = dim(X)[1]
  pp = dim(X)[2]
  if (nn >= pp){
    first_dim = pp
  }else{
    first_dim = nn
  }
  
  # initialize alphas and coefs matrix
  coef = array(0, dim = c(first_dim, ff, bb))
  beta_rr_sf <- function(alpha, lambda){
    p = length(lambda)
    l = length(alpha)
    alpha_grid = matrix(rep(alpha, p), p, l, byrow = TRUE)
    sf_mat = sweep(alpha_grid, 1, lambda^2, "+")  # per row
    sf_mat = lambda/sf_mat
    return(sf_mat)
  }
  lambda = svd_res$lambda_val
  vt = svd_res$vt
  ytr = svd_res$y_trsf
  sf = beta_rr_sf(alpha, lambda)
  
  for (ii in 1:ncol(y)){
    beta_rr_trsf = sf * as.vector(ytr[,ii])
    coef[,,ii] = beta_rr_trsf
  }
  dim(coef) = c(first_dim, ff * bb)
  beta_ridge = t(vt) %*% coef
  dim(beta_ridge) = c(pp, ff, bb)
  beta_ridge = aperm(beta_ridge, c(1,3,2))
  if (bb == 1){
    beta_ridge = matrix(beta_ridge, nrow = pp, ncol = ff)
  }
  return(beta_ridge)
}

#' Ridge regression for multiple targets using library(glmnet) 
#' We add a for loop to incorporate the situation when we have multiple targets
# install.packages("glmnet")
library(glmnet)
ridge_glmnet = function(X, y, lambda, intercept = FALSE, scale = FALSE){
  for (i in 1:ncol(y)){
    res = glmnet(X, y[,i], family = 'gaussian', intercept = intercept, alpha = 0, lambda = lambda, scale = scale)
  }
  return(coef(res))
}


############## SIMULATED DATA STUDY ###################

########## BASE CASE: n = 5000, p = 100, t = 1, f = 10
# install.packages("microbenchmark")
# microbenchmark is used to measure and compare the computational time
library(microbenchmark)
# fix n, t, f, and change p = 10, 20, 50, 100, 200, 300, 500, 750, 1000, 1500, 2000, 3000
n= 5000
t = 1
f = 10
p = c(10, 20, 50, 100, 200, 300, 500, 750, 1000, 1500, 2000, 3000)
p_param = matrix(NA, nrow = 4, ncol = length(p)) # we have 4 methods to compare
colnames(p_param) = p
alphas = seq(1, 100, length = f)
fracs = seq(0.1, 1, length = f)
for (i in 1:length(p)){
  Xmat = gen_large_sample_data(n, p[i], t)$X
  y = gen_large_sample_data(n, p[i], t)$y
  res_p = summary(microbenchmark(standard_ridge_reg(Xmat, y, alphas), svd_ridge_reg(Xmat, y, alphas), frac_ridge(Xmat, y, fracs), glmnet(Xmat, y, lambda = alphas), unit = 's', times = 10L))
  p_param[,i] = res_p$mean
  # print(i) ## tell which iteration you are at
}


# fix n, p, t, and change f = 1, 5, 10, 20, 50, 75, 100
Xmat = gen_large_sample_data(5000, 100, 1)$X
y = gen_large_sample_data(5000, 100, 1)$y
f = c(1, 5, 10, 20, 50, 75, 100)
f_param = matrix(NA, nrow = 4, ncol = length(f))
colnames(f_param) = f
for (i in 1:length(f)){
  alphas = seq(1, 100, length = f[i])
  fracs = seq(0.1, 1, length = f[i])
  res_f = summary(microbenchmark(standard_ridge_reg(Xmat, y, alphas), svd_ridge_reg(Xmat, y, alphas), frac_ridge(Xmat, y, fracs), glmnet(Xmat, y, lambda = alphas), unit = 's', times = 10L))
  f_param[,i] = res_f$mean
  # print(i)
}


# fix p, t, f, and change n = 25, 50, 100, 500, 1000, 5000, 10000, 15000, 20000
n = c(25, 50, 100, 500, 1000, 5000, 10000, 15000, 20000)
f = 10
t = 1
p = 100
n_param = matrix(NA, nrow = 4, ncol = length(n))
alphas = seq(1, 100, length = f)
fracs = seq(0.1, 1, length = f)
for (i in 1:length(n)){
  Xmat = gen_large_sample_data(n[i], p, t)$X
  y = gen_large_sample_data(n[i], p, t)$y
  res_n = summary(microbenchmark(standard_ridge_reg(Xmat, y, alphas), svd_ridge_reg(Xmat, y, alphas), frac_ridge(Xmat, y, fracs), glmnet(Xmat, y, lambda = alphas), unit = 's', times = 10L))
  n_param[,i] = res_n$mean
  # print(i)
}


# fix n, p, f, and change t = 1, 5, 25, 50, 100, 500, 1000, 2500, 5000
t = c(1, 5, 25, 50, 100, 500, 1000, 2500, 5000)
f = 10
n = 5000
p = 100
t_param = matrix(NA, nrow = 4, ncol = length(t))
colnames(t_param) = t
alphas = seq(1, 100, length = f)
fracs = seq(0.1, 1, length = f)
for (i in 1:length(t)){
  Xmat = gen_large_sample_data(n, p, t[i])$X
  y = gen_large_sample_data(n, p, t[i])$y
  res_t = summary(microbenchmark(standard_ridge_reg(Xmat, y, alphas), svd_ridge_reg(Xmat, y, alphas), frac_ridge(Xmat, y, fracs), ridge_glmnet(Xmat, y, lambda = alphas), unit = 's', times = 10L))
  t_param[,i] = res_t$mean
  # print(i)
}




