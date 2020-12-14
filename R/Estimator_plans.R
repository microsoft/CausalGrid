
# Estimator Fns -----------
cont_te_estimator <- function(y, d, ...) {
  if(is_vec(d)) {
    # Straight formulas is much faster than OLS
    #formula reference: http://cameron.econ.ucdavis.edu/e240a/reviewbivariate.pdf
    y_avg = mean(y)
    d_avg = mean(d)
    d_demean = d-d_avg
    sum_d_dev = sum(d_demean^2)
    param_est = sum((y-y_avg)*d_demean)/sum_d_dev
  }
  else {
    ols_fit = lm(y~d)
    param_est = coef(ols_fit)[-1]
  }
  return(list(param_est=param_est))
}

cont_te_var_estimator <- function(y, d, ...) {
  if(is_vec(d)) {
    # Straight formulas is much faster than OLS
    #formula reference: http://cameron.econ.ucdavis.edu/e240a/reviewbivariate.pdf
    y_avg = mean(y)
    d_avg = mean(d)
    d_demean = d-d_avg
    sum_d_dev = sum(d_demean^2)
    param_est = sum((y-y_avg)*d_demean)/sum_d_dev
    b0 = y_avg - param_est*d_avg
    y_hat = b0+param_est*d
    err = y - y_hat
    var_est = (sum(err^2)/(length(y)-2))/sum_d_dev
  }
  else {
    if(length(y)==0) {
      print("Ahh")
    }
    ols_fit = lm(y~d)
    param_est = coef(ols_fit)[-1]
    var_est = diag(vcov(lm(y~d)))[-1]
  }
  return(list(param_est=param_est, var_est=var_est))
}

#Handles removing factors with only 1 level
robust_lm_d <- function(y, d, X, ctrl_names) {
  ctrl_str = if(length(ctrl_names)>0) paste0("+", paste(ctrl_names, collapse="+")) else ""
  tryCatch(ols_fit <- lm(formula(paste0("y~d", ctrl_str)), data=as.data.frame(X)),
           error=function(e) {
             ctrl_names2 <- ctrl_names[sapply(ctrl_names, function(ctrl_name){length(unique(X[, ctrl_name]))}) > 1]
             ctrl_str2 <- if(length(ctrl_names2)>0) paste0("+", paste(ctrl_names2, collapse="+")) else ""
             ols_fit <<- lm(formula(paste0("y~d", ctrl_str2)), data=as.data.frame(X))
           })
  return(ols_fit)
}

cont_te_X_estimator <- function(y, d, X, ctrl_names) {
  d_ncols = if(is_vec(d)) 1 else ncol(d)
  ols_fit = robust_lm_d(y, d, X, ctrl_names)
  param_est=coef(ols_fit)[2:(1+d_ncols)]
  return(list(param_est=param_est))
}



cont_te_var_X_estimator <- function(y, d, X, ctrl_names) {
  d_ncols = if(is_vec(d)) 1 else ncol(d)
  ols_fit = robust_lm_d(y, d, X, ctrl_names)
  param_est=coef(ols_fit)[2:(1+d_ncols)]
  var_est=diag(vcov(ols_fit))[2:(1+d_ncols)]
  return(list(param_est=param_est, var_est=var_est))
}

lcl_colMeans <- function(y) {
  if(is.list(y)) #list of dataframe
    return(sapply(y, mean))
  if(is_vec(y)) #vector
    return(mean(y))
  #matrix
  return(colMeans(y))
}

lcl_colVars_est <- function(y) {
  if(is.list(y)) #list of dataframe
    return(sapply(y, function(c) var(c)/(length(c)-1)))
  if(is_vec(y)) #vector
    return(var(y)/(length(y)-1))
  #matrix
  return(apply(y, 2, function(c) var(c)/(length(c)-1)))
}

mean_var_estimator <- function(y, ...) {
  #int_str = "(Intercept)" #"const"
  #ols_fit <- lm(y~1)
  #param_est=coef(ols_fit)[int_str]
  #var_est=vcov(ols_fit)[int_str, int_str]
  # The below is much faster
  
  return(list(param_est=lcl_colMeans(y), var_est=lcl_colVars_est(y)))
}

mean_estimator <- function(y, ...) {
  return(list(param_est=lcl_colMeans(y)))
}


# Generics ---------------

#Aside from these generics, subclasses must have $dof scalar

#' Fit_InitTr
#'
#' @param obj Object
#' @param X_tr X
#' @param y_tr y
#' @param d_tr d_tr
#' @param cv_folds CV folds
#' @param verbosity verbosity
#' @param dim_cat vector of dimensions that are categorical
#'
#' @return Updated Object
#' @export
Fit_InitTr <- function(obj, X_tr, y_tr, d_tr=NULL, cv_folds, verbosity=0, dim_cat=c()) { UseMethod("Fit_InitTr", obj)}


#' Do_Residualize
#'
#' @param obj Object
#' @param y y
#' @param X X
#' @param d d
#' @param d d (Default=NULL)
#' @param sample one of 'tr' or 'est'
#'
#' @return list(y=) or list(y=, d=)
#' @export
Do_Residualize <- function(obj, y, X, d, sample) { UseMethod("Do_Residualize", obj)}

#' Param_Est
#'
#' @param obj Object
#' @param y y A N-vector
#' @param d d A N-vector or Nxm matrix (so that they can be estimated jointly)
#' @param X X A NxK matrix or data.frame
#' @param sample Sample: "trtr", "trcv", "est" 
#' @param ret_var Return Variance in the return list
#'
#' @return list(param_est=...)
#' @export
Param_Est <- function(obj, y, d=NULL, X, sample="est", ret_var=FALSE) { UseMethod("Param_Est", obj)}

# lm_X_est ---------------

lm_X_est <- function(lasso=FALSE, control_est=TRUE) {
  return(structure(list(lasso=lasso, control_est=control_est), class = c("Estimator_plan","lm_X_est")))
}

#' is.lm_X_est
#'
#' @param x Object
#'
#' @return Boolean
#' @export
is.lm_X_est <- function(x) {inherits(x, "lm_X_est")}

dummyVar_common <- function(X, dim_cat) {
  X_new = NULL
  groups = c()
  #Since regularizing, be careful about the reference class (so can't use dummyVars or one_hot easily)
  for(k in 1:ncol(X)) {
    n_cols=1
    X_k = X[[k]]
    if(k %in% dim_cat) {
      n_l = nlevels(X_k)
      level_common = dimnames(sort(table(factor(X_k)), decreasing=T))[[1]][1]
      level_common_int = match(level_common, levels(X_k))
      X_k = model.matrix(~1+C(X_k, contr.treatment(n_l, base=level_common_int)))[, 2:n_l] #col=1 is intercept
      n_cols = ncol(X_k)
    }
    if(is.null(X_new)) X_new = X_k
    else X_new = cbind(X_new, X_k)
    groups = c(groups, rep(k, n_cols))
  }
  return(list(X_new, groups))
}

lasso_select <- function(obj, X_tr, y_tr, cv_folds, verbosity, dim_cat) {
  if(length(dim_cat)<1) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("Package \"glmnet\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    if(is.data.frame(X_tr)) X_tr = as.matrix(X_tr)
    if(length(cv_folds)==1)
      lasso_fit = glmnet::cv.glmnet(X_tr, y_tr, nfolds=cv_folds)
    else
      lasso_fit = glmnet::cv.glmnet(X_tr, y_tr, foldid=cv_folds)
    c = coef(lasso_fit, s = "lambda.min")
    sel = c[2:length(c), ]!=0
  }
  else {
    list[X_new, groups] = dummyVar_common(X_tr, dim_cat)
    if (!requireNamespace("gglasso", quietly = TRUE)) {
      stop("Package \"gglasso\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    if(length(cv_folds)==1) {
      gg_fit = gglasso::cv.gglasso(X_new, y_tr, nfolds=cv_folds, loss="ls", groups)
    }
    else {
      gg_fit = gglasso::cv.gglasso(X_new, y_tr, foldid=cv_folds, loss="ls", groups)
    }
    c = coef(gg_fit, s="lambda.min")
    sel = sort(unique(groups[c[2:length(c)]!=0]))
  }
  if(verbosity>0) print(c)
  return(colnames(X_tr)[sel])
}

#' Fit_InitTr.lm_X_est
#'
#' @param obj lm_X_est object
#' @param X_tr X_tr
#' @param y_tr y_tr
#' @param d_tr d_tr
#' @param cv_folds cv_folds
#' @param verbosity verbosity
#' @param dim_cat dim_cat
#'
#' @return Updated object
#' @export
#' @method Fit_InitTr lm_X_est
Fit_InitTr.lm_X_est <- function(obj, X_tr, y_tr, d_tr=NULL, cv_folds, verbosity=0, dim_cat=c()) {
  assert_that(!is.null(d_tr))
  if(obj$lasso & length(dim_cat)<1) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("Package \"ranger\" needed for this function to work. Please install it.", call. = FALSE)
    }
  }
  
  list[M, m_mode, N, K] = get_sample_type(y_tr, X_tr, d_tr, checks=TRUE)

  if(m_mode==0 || m_mode==2) {
    if(obj$lasso)
      obj$ctrl_names = lasso_select(obj, X_tr, y_tr, cv_folds, verbosity, dim_cat)
    else
      obj$ctrl_names = colnames(X_tr)
    obj$dof = 2+length(obj$ctrl_names)
  }
  else {
    if(obj$lasso) {
      if(m_mode==1)
        obj$ctrl_names = mapply(function(X_s, y_s, cv_folds_s) lasso_select(obj, X_s, y_s, cv_folds_s, verbosity, dim_cat), X_tr, y_tr, cv_folds)
      if(m_mode==3)
        obj$ctrl_names = apply(y_tr, 2, function(y_col) lasso_select(obj, X_tr, y_col, cv_folds, verbosity, dim_cat))
    }
    else {
      obj$ctrl_names = rep(list(colnames(X_tr)), M)
    }
    obj$dof = 2+sapply(obj$ctrl_names, length)
  }
  if(verbosity>0) cat(paste("LassoCV-picked control variables: ", paste(obj$ctrl_names, collapse=" "), "\n"))

  return(obj)
}

#' Do_Residualize.lm_X_est
#'
#' @param obj obj
#' @param y y
#' @param X X
#' @param d d
#' @param sample one of 'tr' or 'est'
#'
#' @return list(y=...) or list(y=..., d=...)
#' @export
#' @method Do_Residualize lm_X_est
Do_Residualize.lm_X_est <- function(obj, y, X, d, sample) {return(list(y=y, d=d))}

#' Param_Est.lm_X_est
#'
#' @param obj obj
#' @param y y
#' @param d d
#' @param X X
#' @param sample Sample: "trtr", "trcv", "est" 
#' @param ret_var Return variance in return list
#'
#' @return list(param_est=...) or list(param_est=..., var_est=...)
#' @export
#' @method Param_Est lm_X_est
Param_Est.lm_X_est <- function(obj, y, d=NULL, X, sample="est", ret_var=FALSE) {
  assert_that(!is.null(d))
  
  if(sample=="trtr" || (sample=="est" && obj$control_est)) {
    if(ret_var) return(cont_te_var_X_estimator(y, d, X, obj$ctrl_names))
    return(cont_te_X_estimator(y, d, X, obj$ctrl_names))
  }
  if(ret_var) return(cont_te_var_estimator(y, d, X))
  return(cont_te_estimator(y, d, X))
}


# simple_est ---------------

simple_est <- function(te_fn, te_var_fn, dof=2) {
  return(structure(list(te_fn=te_fn, te_var_fn=te_var_fn, dof=dof), class = c("Estimator_plan", "simple_est")))
} 

#' is.simple_est
#'
#' @param x Object
#'
#' @return Boolean
#' @export
is.simple_est <- function(x) {inherits(x, "simple_est")}

gen_simple_est_plan <- function(has_d=TRUE) {
  if(has_d) {
    return(simple_est(cont_te_estimator, cont_te_var_estimator))
  }
  return(simple_est(mean_estimator, mean_var_estimator, dof=1))
}

#' Fit_InitTr.simple_est
#'
#' @param obj obj
#' @param X_tr X_tr
#' @param y_tr y_tr
#' @param d_tr d_tr
#' @param cv_folds cv_folds
#' @param verbosity verbosity
#' @param dim_cat dim_cat
#'
#' @return Updated object
#' @export
#' @method Fit_InitTr simple_est
Fit_InitTr.simple_est <- function(obj, X_tr, y_tr, d_tr=NULL, cv_folds, verbosity=0, dim_cat=c()) {return(obj)}

#' Do_Residualize.simple_est
#'
#' @param obj obj
#' @param y y
#' @param X X
#' @param d d
#' @param sample one of 'tr' or 'est'
#'
#' @return list(y=...) and list(y=..., d=...)
#' @export
#' @method Do_Residualize simple_est
Do_Residualize.simple_est <- function(obj, y, X, d, sample) {return(list(y=y, d=d))}

#' Param_Est.simple_est
#'
#' @param obj obj
#' @param y y
#' @param d d
#' @param X X
#' @param sample Sample: "trtr", "trcv", "est" 
#' @param ret_var Return variance in return list
#'
#' @return list(param_est=...)
#' @export
#' @method Param_Est simple_est
Param_Est.simple_est <- function(obj, y, d=NULL, X, sample="est", ret_var=FALSE) {
  if(ret_var) return(obj$te_var_fn(y, d, X))
  return(obj$te_fn(y, d, X))
}

# grid_rf ---------------

#' grid_rf
#'
#' @param num.trees number of trees in the random forest
#' @param num.threads num.threads
#' @param dof degrees-of-freedom
#' @param resid_est Residualize the Estimation sample (using fit from training)
#'
#' @return grid_rf object
#' @export
grid_rf <- function(num.trees=500, num.threads=NULL, dof=2, resid_est=TRUE) {
  return(structure(list(num.trees=num.trees, num.threads=num.threads, dof=dof, resid_est=resid_est), 
                   class = c("Estimator_plan","grid_rf")))
} 

#' is.grid_rf
#'
#' @param x Object
#'
#' @return Boolean
#' @export
is.grid_rf <- function(x) {inherits(x, "grid_rf")}

rf_fit_data <- function(obj, target, X) {
  if(is_vec(target))
    return(ranger::ranger(y=target, x=X, 
                  num.trees = obj$num.trees, num.threads = obj$num.threads))
  fits = list()
  for(m in 1:ncol(target)){
    fits[[m]] = ranger::ranger(y=target[,m], x=X, 
           num.trees = obj$num.trees, num.threads = obj$num.threads)
  }
  return(fits)
}

#' Fit_InitTr.grid_rf
#' Note that for large data, the rf_y_fit and potentially rf_d_fit objects may be large.
#' They can be null'ed out after fitting
#'
#' @param obj Object
#' @param X_tr X
#' @param y_tr y
#' @param d_tr d_tr
#' @param cv_folds CV folds
#' @param verbosity verbosity
#' @param dim_cat vector of dimensions that are categorical
#'
#' @return Updated Object
#' @export
#' @method Fit_InitTr grid_rf
Fit_InitTr.grid_rf <- function(obj, X_tr, y_tr, d_tr=NULL, cv_folds, verbosity=0, dim_cat=c()) {
  assert_that(!is.null(d_tr)) #Only residualize when having treatment
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package \"ranger\" needed for this function to work. Please install it.", call. = FALSE)
  }
  obj$rf_y_fit = rf_fit_data(obj, y_tr, X_tr)
  
  if(!is.null(d_tr)) {
    obj$rf_d_fit = rf_fit_data(obj, d_tr, X_tr)
  }
  return(obj)
}

rf_predict_data <- function(fit, target, X) {
  if(is_vec(target))
    return(predict(fit, X, type="response")$predictions)
  preds = matrix(NA, nrow=nrow(X), ncol=ncol(X))
  for(m in 1:ncol(target)){
    preds[,m] = predict(fit[[m]], X, type="response")$predictions
  }
  return(preds)
}

#' Do_Residualize.grid_rf
#'
#' @param obj Object
#' @param y y
#' @param X X
#' @param d d (Default=NULL)
#' @param sample one of 'tr' or 'est'
#'
#' @return list(y=) or list(y=, d=)
#' @export
#' @method Do_Residualize grid_rf
Do_Residualize.grid_rf <- function(obj, y, X, d, sample) { 
  if(sample=="est" && !obj$resid_est) return(list(y=y, d=d))
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package \"ranger\" needed for this function to work. Please install it.", call. = FALSE)
  }
  y_res = y - rf_predict_data(obj$rf_y_fit, y, X)
  d_res = if(is.null(d)) NULL else d - rf_predict_data(obj$rf_d_fit, d, X)
  return(list(y=y_res, d=d_res))
}

#' Param_Est.grid_rf
#'
#' @param obj Object
#' @param y y
#' @param d d
#' @param X X
#' @param sample Sample: "trtr", "trcv", "est" 
#' @param ret_var Return Variance in the return list
#'
#' @return list(param_est=...)
#' @export
#' @method Param_Est grid_rf
Param_Est.grid_rf <- function(obj, y, d=NULL, X, sample="est", ret_var=FALSE) {
  assert_that(is.flag(ret_var), sample %in% c("est", "trtr", "trcv"), !is.null(d))
  
  if(ret_var) return(cont_te_var_estimator(y, d, X))
  return(cont_te_estimator(y, d, X))
}