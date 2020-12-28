
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

#can handle MULTI_D
cont_te_X_estimator <- function(y, d, X, ctrl_names) {
  d_ncols = if(is_vec(d)) 1 else ncol(d)
  ols_fit = robust_lm_d(y, d, X, ctrl_names)
  param_est=coef(ols_fit)[2:(1+d_ncols)]
  return(list(param_est=param_est))
}


#can handle MULTI_D
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

#' Estimator Plan
#' 
#' An estimator plan is an S3 objet that defines functions for how to estimate parameters in the various stages of subgroup analysis.
#' 
#' With an outcome-only model, then we just call \code{\link{est_params}}. 
#' When there is a treatment, then we initially also call \code{\link{fit_on_train}}, \code{\link{residualize}}.
#'
#' Each plan also has a field \code{dof} for a degrees of freedom-correction to standard errors (e.g. 2 (intercept and treatment effect) or more if controls are included) 
#' Provided estimator plans include \code{\link{simple_est}}, \code{\link{lm_est}}, and \code{\link{grid_rf}}. User provided plans just need to implement the relevant functions.
#'
#' @name EstimatorPlan
NULL
#> NULL

#Aside from these generics, subclasses must have $dof scalar

#' Fit estimator_plan on training sample
#'
#' @param obj an EstimatorPlan object
#' @param X_tr N_trxK matrix
#' @param y_tr N_tr-vector
#' @param d_tr NULL, N_tr-vector, or N_trxM matrix (if multiple treatments)
#' @param cv_folds CV folds
#' @param verbosity verbosity
#' @param dim_cat vector of dimensions that are categorical
#'
#' @return Updated EstimatorPlan object
#' @export
fit_on_train <- function(obj, X_tr, y_tr, d_tr, cv_folds, verbosity=0, dim_cat=c()) { UseMethod("fit_on_train", obj)}


#' Residualize sample
#' 
#' Residualize y and d (only called when doing treatment effect estimation)
#'
#' @param obj an EstimatorPlan object
#' @param y N-vector
#' @param X NxK matrix
#' @param d NULL, N-vector, or NxM matrix (if multiple treatments)
#' @param sample one of \code{'tr'} or \code{'est'}
#'
#' @return \code{list(y=, d=)}
#' @export
residualize <- function(obj, y, X, d, sample) { UseMethod("residualize", obj)}

#' Estimate parameters
#'
#' @param obj an EstimatorPlan object
#' @param y A N-vector
#' @param d A N-vector or NxM matrix (so that they can be estimated jointly)
#' @param X A NxK matrix or data.frame
#' @param sample One of: "trtr", "trcv", "est" 
#' @param ret_var Return Variance in the return list
#'
#' @return \code{list(param_est=...)} or \code{list(param_est=...)} if \code{ret_var} 
#' @export
est_params <- function(obj, y, d=NULL, X, sample="est", ret_var=FALSE) { UseMethod("est_params", obj)}

# lm_est ---------------

#' create an lm_est object
#'
#' @param lasso  Use Lasso on the train sample to choose control variables
#' @param control_est Use controls when estimating effects on estimation sample
#'
#' @return Estimation plan
#' @export
lm_est <- function(lasso=FALSE, control_est=TRUE) {
  return(structure(list(lasso=lasso, control_est=control_est), class = c("estimator_plan","lm_est")))
}

#' Is lm_est
#'
#' @param x an R object
#'
#' @return Boolean
#' @export
#' @describeIn lm_est is lm_est 
is.lm_est <- function(x) {inherits(x, "lm_est")}

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

fit_on_train.lm_est <- function(obj, X_tr, y_tr, d_tr, cv_folds, verbosity=0, dim_cat=c()) {
  if(obj$lasso & length(dim_cat)<1) {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      stop("Package \"ranger\" needed for this function to work. Please install it.", call. = FALSE)
    }
  }

  if(obj$lasso)
    obj$ctrl_names = lasso_select(obj, X_tr, y_tr, cv_folds, verbosity, dim_cat)
  else
    obj$ctrl_names = colnames(X_tr)
  obj$dof = 2+length(obj$ctrl_names)
  if(verbosity>0) cat(paste("LassoCV-picked control variables: ", paste(obj$ctrl_names, collapse=" "), "\n"))

  return(obj)
}

residualize.lm_est <- function(obj, y, X, d, sample) {return(list(y=y, d=d))}

#' @describeIn est_params lm_est
#' @export
est_params.lm_est <- function(obj, y, d=NULL, X, sample="est", ret_var=FALSE) {
  assert_that(!is.null(d))
  
  if(sample=="trtr" || (sample=="est" && obj$control_est)) {
    if(ret_var) return(cont_te_var_X_estimator(y, d, X, obj$ctrl_names))
    return(cont_te_X_estimator(y, d, X, obj$ctrl_names))
  }
  if(ret_var) return(cont_te_var_estimator(y, d, X))
  return(cont_te_estimator(y, d, X))
}


# simple_est ---------------

#' Create a simple_est object
#'
#' @param te_fn Bare Treatment Estimation function 
#' @param te_var_fn Treatment estimatino function when variance needed
#' @param dof 
#'
#' @return Estimation Plan
#' @export
simple_est <- function(te_fn, te_var_fn, dof=2) {
  return(structure(list(te_fn=te_fn, te_var_fn=te_var_fn, dof=dof), class = c("estimator_plan", "simple_est")))
} 

#' Is simple_est
#'
#' @param x an R object
#'
#' @return Boolean
#' @export
#' @describeIn simple_est is simple_est
is.simple_est <- function(x) {inherits(x, "simple_est")}

gen_simple_est_plan <- function(has_d=TRUE) {
  if(has_d) {
    return(simple_est(cont_te_estimator, cont_te_var_estimator))
  }
  return(simple_est(mean_estimator, mean_var_estimator, dof=1))
}

fit_on_train.simple_est <- function(obj, X_tr, y_tr, d_tr, cv_folds, verbosity=0, dim_cat=c()) {
  obj$dof = 1 + if(is_vec(d_tr)) 1 else ncol(d_tr)
  return(obj)
}

residualize.simple_est <- function(obj, y, X, d, sample) {return(list(y=y, d=d))}

#' @describeIn est_params simple_est
#' @export
est_params.simple_est <- function(obj, y, d=NULL, X, sample="est", ret_var=FALSE) {
  if(ret_var) return(obj$te_var_fn(y, d, X))
  return(obj$te_fn(y, d, X))
}

# grid_rf ---------------

#' Create a grid_rf object
#' 
#' Residualizes the train sample using cross-fitting
#' Residualizes the estimation samples using a train fit
#'
#' @param num.trees number of trees in the random forest
#' @param num.threads num.threads
#' @param dof degrees-of-freedom
#' @param resid_est Residualize the Estimation sample (using fit from training)
#'
#' @return grid_rf object
#' @export
grid_rf <- function(num.trees=500, num.threads=NULL, resid_est=TRUE) {
  return(structure(list(num.trees=num.trees, num.threads=num.threads, resid_est=resid_est), 
                   class = c("estimator_plan","grid_rf")))
} 

#' Is grid_rf
#'
#' @param x an R object
#'
#' @return Boolean
#' @export
#' @describeIn grid_rf is grid_rf
is.grid_rf <- function(x) {inherits(x, "grid_rf")}

ranger_cross_fit <- function(target, X, cv_folds, num.trees, num.threads) {
  n_folds = max(cv_folds)
  fits = list()
  for(f in 1:n_folds) {
    mask = (cv_folds!=f)
    fits[[f]] = ranger::ranger(y=target[mask], x=X[mask,], 
                               num.trees = num.trees, num.threads = num.threads)
  }
  return(structure(fits, class = c("ranger_cross_fit")))
}

predict.ranger_cross_fit <- function(obj, X, cv_folds) {
  preds = rep(NA, length(cv_folds))
  n_folds = max(cv_folds)
  for(f in 1:n_folds) {
    mask =(cv_folds==f)
    preds[mask] = predict(obj[[f]], X[mask,], type="response")$predictions
  }
  return(preds)
}

#target can be a vector or matrix (in which case we return a list of fits)
rf_fit_data <- function(obj, target, X, cv_folds=NULL) {
  if(!is.null(cv_folds)) {
    if(is_vec(target))
      return(ranger_cross_fit(target, X, cv_folds, obj$num.trees, obj$num.threads))
    fits = list()
    for(m in 1:ncol(target)){
      fits[[m]] = ranger_cross_fit(target[,m], X, cv_folds, obj$num.trees, obj$num.threads)
    }
  }
  else {
    if(is_vec(target))
      return(ranger::ranger(y=target, x=X, 
                            num.trees = obj$num.trees, num.threads = obj$num.threads))
    fits = list()
    for(m in 1:ncol(target)){
      fits[[m]] = ranger::ranger(y=target[,m], x=X, num.trees = obj$num.trees, num.threads = obj$num.threads)
    }
    
  }
  return(fits)
}

fit_on_train.grid_rf <- function(obj, X_tr, y_tr, d_tr, cv_folds, verbosity=0, dim_cat=c()) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package \"ranger\" needed for this function to work. Please install it.", call. = FALSE)
  }
  obj$cv_folds = cv_folds
  obj$rf_y_xfit = rf_fit_data(obj, y_tr, X_tr, cv_folds)
  obj$rf_y_fit = rf_fit_data(obj, y_tr, X_tr)
  
  if(!is.null(d_tr)) {
    obj$rf_d_xfit = rf_fit_data(obj, d_tr, X_tr, cv_folds)
    obj$rf_d_fit = rf_fit_data(obj, d_tr, X_tr)
  }
  obj$dof=1+if(is_vec(d_tr)) 1 else ncol(d_tr)
  return(obj)
}

rf_predict_data <- function(fit, target, X, cv_folds=NULL) {
  if(!is.null(cv_folds)) {
    if(is_vec(target))
      return(predict.ranger_cross_fit(fit, X, cv_folds))
    preds = matrix(NA, nrow=nrow(X), ncol=ncol(target))
    for(m in 1:ncol(target)){
      preds[,m] = predict.ranger_cross_fit(fit[[m]], X, cv_folds)
    }
  }
  else {
    if(is_vec(target))
      return(predict(fit, X, type="response")$predictions)
    preds = matrix(NA, nrow=nrow(X), ncol=ncol(target))
    for(m in 1:ncol(target)){
      preds[,m] = predict(fit[[m]], X, type="response")$predictions
    }
  }
  return(preds)
}

residualize.grid_rf <- function(obj, y, X, d, sample) { 
  if(sample=="est" && !obj$resid_est) return(list(y=y, d=d))
  if(is.data.frame(X)) {
    X = as.matrix(data.frame(lapply(X, as.numeric))) #predict.ranger can't handle factor variables
  }
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package \"ranger\" needed for this function to work. Please install it.", call. = FALSE)
  }
  if(sample=="est") {
    y_res = y - rf_predict_data(obj$rf_y_fit, y, X)
    d_res = d - rf_predict_data(obj$rf_d_fit, d, X)
  }
  else {
    y_res = y - rf_predict_data(obj$rf_y_xfit, y, X, cv_folds=obj$cv_folds)
    d_res = d - rf_predict_data(obj$rf_d_xfit, d, X, cv_folds=obj$cv_folds)
  }
  return(list(y=y_res, d=d_res))
}

#' @describeIn est_params grid_rf
#' @export
est_params.grid_rf <- function(obj, y, d=NULL, X, sample="est", ret_var=FALSE) {
  assert_that(is.flag(ret_var), sample %in% c("est", "trtr", "trcv"), !is.null(d))
  
  if(ret_var) return(cont_te_var_estimator(y, d, X))
  return(cont_te_estimator(y, d, X))
}