# Standard ridge regression
standard_ridge_reg = function(X, y, alpha){
  std_br = array(0, dim = c(ncol(X), ncol(y), length(alpha)))
  ##### incorporate situation with multitargets
  for (i in 1:length(alpha)){
    std_br[,,i] = (solve(t(X) %*% X + diag(alpha[i], ncol(X), ncol(X)))) %*% t(X) %*% y 
  }
  if (ncol(y) == 1){
    std_br = matrix(std_br, nrow = ncol(X), ncol = length(alpha))
  }
  return(std_br)
}

# SVD decomposition ridge regression
svd_ridge_reg = function(X, y, alpha){
  
  svd_res = svd_x(X, y)
  lambda = svd_res$lambda_val
  vt = svd_res$vt
  ytr = svd_res$y_trsf
  
  svd_br = array(NA, dim = c(ncol(X), ncol(y), length(alpha)))
  for (i in 1:length(alpha)){
    svd_br[,,i] = t(vt) %*% solve(diag(lambda^2)  + diag(alpha[i], ncol(X), ncol(X))) %*% diag(lambda) %*% ytr
  }
  svd_br <- aperm(svd_br, c(1,3,2))
  return(svd_br)
}

# ridge regression from glmnet
install.packages("glmnet")
library(glmnet)
ridge_glmnet = function(X, y, lambda, intercept = FALSE, scale = FALSE){
  for (i in 1:ncol(y)){
    res = glmnet(X, y[,i], family = 'gaussian', intercept = intercept, alpha = 0, lambda = lambda, scale = scale)
  }
  return(coef(res))
}


############  Real-world data study
# Data resource: Medical Cost Personal Datasets from Kaggle.
# https://www.kaggle.com/datasets/mirichoi0218/insurance
# read data
data <- read.csv("insurance.csv")

# converge to factor
data$sex.f[data$sex %in% "male"] <- 0
data$sex.f[data$sex %in% "female"] <- 1

data$smoker.f[data$smoker %in% "yes"] <- 0
data$smoker.f[data$smoker %in% "no"] <- 1

data$region.f[data$region %in% "northeast"] <- 0
data$region.f[data$region %in% "northwest"] <- 1
data$region.f[data$region %in% "southeast"] <- 2
data$region.f[data$region %in% "southwest"] <- 3

data = subset(data, select = -c(sex,smoker,region) )

# Data spliting to Training and Testing
sample <- sample(c(TRUE, FALSE), nrow(data), replace=TRUE, prob=c(0.7,0.3))
train  <- data[sample, ]
test   <- data[!sample, ]

data_train_X <- as.matrix(subset(train, select=-c(charges)))
data_train_Y <- as.matrix(subset(train, select=c(charges)))

data_test_X <- as.matrix(subset(test, select=-c(charges)))
data_test_Y <- as.matrix(subset(test, select=c(charges)))

# Calculate R^2 and L2-norm of beta weights for fractional ridge regression
fraction <- seq(0.01,1,0.01)
a = frac_ridge(data_train_X, data_train_Y, fraction)
beta_ridge <- a$br
alpha <- a$alphas

y_pred <- frac_ridge_pred(data_test_X, beta_ridge)
r2_fraction_ridge <- r2_score(data_test_Y,y_pred)
y_pred_train <- frac_ridge_pred(data_train_X, beta_ridge)
r2_fraction_ridge_train <- r2_score(data_train_Y,y_pred_train)

beta_ridge_new <- matrix(beta_ridge, nrow = dim(data_test_X)[2], ncol = length(fraction))
L2_norm_fraction <- (colSums(beta_ridge_new^2))^0.5

# Calculate R^2 and L2-norm of beta weights for standard ridge regression
beta_ridge_2 = standard_ridge_reg(data_train_X, data_train_Y,alpha)
y_pred <- data_test_X %*% beta_ridge_2
y_pred_train <- data_train_X %*% beta_ridge_2
r2_standard_ridge <- c()
r2_standard_ridge_train <- c()
for (i in 1:dim(y_pred)[2]){
  numer <- sum((data_test_Y-y_pred[,i])^2)
  denom <- sum((data_test_Y-mean(data_test_Y))^2)
  r2_standard_ridge <- append(r2_standard_ridge,(1 - numer/denom))
}

for (i in 1:dim(y_pred_train)[2]){
  numer <- sum((data_train_Y-y_pred_train[,i])^2)
  denom <- sum((data_train_Y-mean(data_train_Y))^2)
  r2_standard_ridge_train <- append(r2_standard_ridge_train,(1 - numer/denom))
}

L2_norm_standard <- (colSums(beta_ridge_2^2))^0.5

# Calculate R^2 for SVD ridge 
beta_ridge_3 = svd_ridge_reg(data_train_X, data_train_Y,alpha)
y_pred <- frac_ridge_pred(data_test_X, beta_ridge_3)
r2_svd_ridge <- r2_score(data_test_Y,y_pred)

# Calculate R^2 and L2-norm beta weights for glmnet
beta_ridge_4 = ridge_glmnet(data_train_X, data_train_Y,alpha, intercept = FALSE, scale = FALSE)
beta_ridge_4 <- beta_ridge_4[-1,]
y_pred <- data_test_X %*% beta_ridge_4
y_pred_train <- data_train_X %*% beta_ridge_4
r2_glmnet_ridge <- c()
r2_glmnet_ridge_train <- c()
for (i in 1:dim(y_pred)[2]){
  numer <- sum((data_test_Y-y_pred[,i])^2)
  denom <- sum((data_test_Y-mean(data_test_Y))^2)
  r2_glmnet_ridge <- append(r2_glmnet_ridge,(1 - numer/denom))
}
for (i in 1:dim(y_pred_train)[2]){
  numer <- sum((data_train_Y-y_pred_train[,i])^2)
  denom <- sum((data_train_Y-mean(data_train_Y))^2)
  r2_glmnet_ridge_train <- append(r2_glmnet_ridge_train,(1 - numer/denom))
}
L2_norm_glmnet <- (colSums(beta_ridge_4^2))^0.5

# ggplot R^2 and L2 norm
install.packages("ggplot2")
library("ggplot2")

data_fun <- data.frame(x = log10(alpha),      
                       value= c(r2_standard_ridge,r2_standard_ridge_train),
                       Performance = rep(c("Test" , "Training" ),each = length(r2_standard_ridge)))
p <- ggplot(data_fun,aes(x , value, col = Performance)) + geom_line() + scale_x_reverse() +
  xlab("Log10(alpha)") + ylab(~~italic(R)^2) + labs(title="Standard Ridge Regression") +
  theme(plot.title=element_text(hjust=0.5)) + geom_vline(xintercept = 0.5, linetype="dotted")
p1 <- p+theme_bw() + theme(plot.title=element_text(hjust=0.5))+ theme(legend.position="none")


data_fun <- data.frame(x = fraction,        
                       value= c(r2_fraction_ridge,r2_fraction_ridge_train),
                       Performance = rep(c("Test" , "Training" ),each = length(r2_fraction_ridge)))
p <- ggplot(data_fun,aes(x , value, col = Performance)) + geom_line() +
  xlab("Fraction")+ylab("")+ labs(title="Fractional Ridge Regression") +
  theme(plot.title=element_text(hjust=0.5)) + geom_vline(xintercept = 0.96, linetype="dotted")
p2 <- p+theme_bw() + theme(plot.title=element_text(hjust=0.5))+ theme(legend.position="none")


data_fun <- data.frame(x = log10(alpha),         
                       value= c(r2_glmnet_ridge,r2_glmnet_ridge_train),
                       Performance = rep(c("Test" , "Training" ),each = length(r2_glmnet_ridge)))
p <- ggplot(data_fun,aes(x , value, col = Performance)) + geom_line() + scale_x_reverse() +
  xlab("Log10(alpha)") +ylab("") + labs(title="Ridge Regression from Glmnet") +
  theme(plot.title=element_text(hjust=0.5)) + geom_vline(xintercept = 0.5, linetype="dotted")
p3 <- p+theme_bw() + theme(plot.title=element_text(hjust=0.5))

data_fun <- data.frame(x = log10(alpha),       
                       value= c(L2_norm_standard))
p <- ggplot(data_fun,aes(x , value)) + geom_line() + scale_x_reverse() +
  xlab("Log10(alpha)") + ylab("L2-norm of beta weights") + labs(title="Standard Ridge Regression") +
  theme(plot.title=element_text(hjust=0.5))
p4 <- p+theme_bw() + theme(plot.title=element_text(hjust=0.5))

data_fun <- data.frame(x = fraction,       
                       value= c(L2_norm_fraction))
p <- ggplot(data_fun,aes(x , value)) + geom_line()  +
  xlab("Fraction") + ylab("") + labs(title="Fractional Ridge Regression") +
  theme(plot.title=element_text(hjust=0.5))
p5 <- p+theme_bw() + theme(plot.title=element_text(hjust=0.5))


data_fun <- data.frame(x = log10(alpha),          
                       value= c(L2_norm_glmnet))
p <- ggplot(data_fun,aes(x , value)) + geom_line() + scale_x_reverse() +
  xlab("Log10(alpha)") + ylab("") + labs(title="Ridge Regression from Glmnet") +
  theme(plot.title=element_text(hjust=0.5))
p6 <- p+theme_bw() + theme(plot.title=element_text(hjust=0.5))

# merging p1,p2,p3 into one figure
install.packages("ggpubr")
library(ggpubr)
ggarrange(p1, p2, p3 + rremove("x.text"), 
          ncol = 3, nrow = 1,common.legend = TRUE, legend="top")
# merging p4,p5,p6 into one figure
ggarrange(p4, p5, p6 + rremove("x.text"), 
          ncol = 3, nrow = 1)



