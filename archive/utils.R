
#' Formula-based t-test
#'
#' @param mean1 
#' @param mean2 
#' @param n1 
#' @param n2 
#' @param s2_1 
#' @param s2_2 
#'
#' @return
#' @export
#'
#' @examples
t_test_form <- function(mean1, mean2, n1, n2, s2_1, s2_2, var.equal=FALSE){
  if (var.equal) {
    if (n1 == n2) {
      sp = sqrt((s2_1 + s2_2)/2)
      t = (mean1 - mean2)/(sp*sqrt(2/n1))
    }
    else{
      sp = sqrt(((n1 - 1)*s2_1 + (n2 - 1)*s2_2)/(n1 + n2 - 2))
      t = (mean1 - mean2)/(sp*sqrt(1/n1 + 1/n2))
    }
    df = n1 + n2 - 2
  }
  else {#Welch's t-test
    s2_n_1 = s2_1/n1
    s2_n_2 = s2_2/n2
    s2_delta = s2_n_1 + s2_n_2
    t = (mean1 - mean2)/(sqrt(s2_delta))
    df = s2_delta^2/((s2_n_1^2/(n1 - 1)) + (s2_n_2^2/(n2 - 1)))
    
  }
  p = 2*pt(abs(t), df, lower.tail = FALSE)
  return(list("statistic" = t, "parameter" = df, "p.value" = p))
}


#update OLS algorith: https://stats.stackexchange.com/questions/23481/
get_beta_fit <- function(X, y, t, new_split, prev=fit_work){
  
}

#returns a KxN matrix where (k,n) tells you the row in X of the n smallest value of X_k
order_idx <- function(X){
  K = ncol(X)
  N = nrow(X)
  ret = matrix(NA, K, N)
  for(k in 1:K) {
    ret[k,] = sort.list(X[,k])
  }
  return(ret)
}