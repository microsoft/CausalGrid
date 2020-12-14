
exp_data <- function(n_4=25, dim_D=1, err_sd=0.01){
  #n_4 is n/4. We get this to make sure we have the same in each chunk
  n = n_4*4
  stopifnot(dim_D %in% c(0,1,2))
  #dim_D in {0,1,2}
  X1 = cbind(runif(n_4, 0, .5), runif(n_4, 0, .5))
  X2 = cbind(runif(n_4, 0, .5), runif(n_4, .5, 1))
  X3 = cbind(runif(n_4, .5, 1), runif(n_4, 0, .5))
  X4 = cbind(runif(n_4, .5, 1), runif(n_4, .5, 1))
  X = rbind(X1, X2, X3, X4)
  
  alpha = ifelse(X[,1]>.5, ifelse(X[,2]>.5,.5,.8), ifelse(X[,2]>.5, 2, -2))
  #alpha=0
  y = alpha + rnorm(n,0,err_sd)
  if(dim_D) {
    if(dim_D==1) {
      beta = ifelse(X[,1]>.5, ifelse(X[,2]>.5,-1,2), ifelse(X[,2]>.5, 4, 6))
      #beta = ifelse(X[,1]>.5, -1,1)
      d = matrix(rnorm(n), n, 1)
      y = y + beta*d
      colnames(d) <- "d"
    }
    else {
      beta1 = ifelse(X[,1]>.5,-1, 4)
      beta2 = ifelse(X[,2]>.5, 2, 6)
      d = matrix(rnorm(2*n), n, 2)
      y = y + beta1*d[,1] + beta2*d[,2]
      colnames(d) = c("d1", "d2")
    }
    
  }
  else {
    d = NULL
  }
  y = as.matrix(y, nrow=n, ncol=1)
  colnames(y) = "y"
  colnames(X) = c("X1", "X2")
  return(list(y=y, X=X, d=d))
}

mix_data_y <- function(n=200) {
  X = data.frame(X1=c(rep(0, n/4), rep(1, n/4), rep(0, n/4), rep(1, n/4)), 
                 X2=factor(c(rep("A", n/2), rep("B", n/2))))
  alpha = c(rep(0, n/4), rep(1, n/4), rep(2, n/4), rep(3, n/4))
  y = alpha
  return(list(y=y, X=X))
}

mix_data_d <- function(n=200) {
  X = data.frame(X1=c(rep(0, n/4), rep(1, n/4), rep(0, n/4), rep(1, n/4)), 
                 X2=factor(c(rep("A", n/2), rep("B", n/4), rep("C", n/4))))
  tau = c(rep(0, n/4), rep(1, n/4), rep(2, n/4), rep(3, n/4))
  d = rep(0:1, n/2)
  y = d*tau
  return(list(y=y, X=X, d=d))
}


two_groups_data <- function(){
  
  X = matrix(factor(c(rep("M", 100), rep("F", 100))),nrow = 200 ,ncol = 1) 
  y = c(rep(5, 100), rep(50, 100))
  
  return(list(y=y, X=X))
}

two_groups_data_int <- function(){
  
  X = matrix(c(rep(1, 100), rep(2, 100), rep(0, 200)) ,nrow = 200 ,ncol = 2) 
  y = c(rep(5, 100), rep(50, 100))
  
  return(list(y=y, X=X))
}

AI_sim <- function(n=500, design=1) {
  w = rbinom(n, 1, 0.5)
  K = c(2, 10, 20)[design]
  X = matrix(rnorm(n*K), nrow=n, ncol=K)
  X_I = X>0
  if(design==1) {
    eta = X %*% matrix(c(0.5, 1), ncol=1)
    kappa = X %*% matrix(c(0.5, 0), ncol=1)
  }
  if(design==2) {
    eta = X %*% matrix(c(rep(0.5, 2), rep(1, 4), rep(0, 4)), ncol=1)
    kappa = (X*X_I) %*% matrix(c(rep(1,2), rep(0,8)), ncol=1)
  }
  if(design==3) {
    eta = X %*% matrix(c(rep(0.5, 4), rep(1, 4), rep(0, 12)), ncol=1)
    kappa = (X*X_I) %*% matrix(c(rep(1,4), rep(0,16)), ncol=1)
  }
  epsilon = rnorm(n, 0, 0.01)
  Y = eta + 0.5*(2*w-1)*kappa + epsilon
  return(list(Y, X, w, kappa))
}
