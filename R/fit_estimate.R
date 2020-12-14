#Notes:
# - In order to work with both X as matrix and data.frame I used X[,k], but this is messed up 
#   with incoming Tibbles so convert those.

# If empty and err cells will be removed from the calculation, but counts of these returned
# estimator_var: takes sample for cell and returns coefficient estimate and estimated variance of estimate
emse_hat_obj <-function(y, X , d, N_est, partition=NULL, cell_factor=NULL, estimator_var=NULL, debug=FALSE, 
                    warn_on_error=FALSE, alpha=NULL, est_plan=NULL, sample="trtr") {
  if(is.null(est_plan)) {
    if(is.null(estimator_var)) est_plan = gen_simple_est_plan(has_d=!is.null(d))
    else est_plan = simple_est(estimator_var, estimator_var)
  }
  if(is.null(cell_factor)) {
    cell_factor = get_factor_from_partition(partition, X)
  }
  if(!is.null(alpha)) {
    stopifnot(alpha>=0 & alpha<=1)
  }
  list[M, m_mode, N, K] = get_sample_type(y, X, d, checks=FALSE)
  list[lvls, n_cells] = lcl_levels(cell_factor)
  cell_contribs1 = rep(0, n_cells)
  cell_contribs2 = rep(0, n_cells)
  N_eff = 0
  N_cell_empty = 0
  N_cell_err = 0
  for(cell_i in 1:n_cells) {
    list[y_cell, d_cell, X_cell, N_l] <- get_cell(y, X, d, cell_factor, cell_i, lvls)
    if(any(N_l==0)) {
      N_cell_empty = N_cell_empty+1
      next
    }
    
    list[param_est, var_est] = Param_Est_m(est_plan, y_cell, d_cell, X_cell, sample=sample, ret_var=TRUE, m_mode) 
    if(!all(is.finite(param_est)) || !all(is.finite(var_est))) {
      N_cell_err = N_cell_err+1
      msg = paste("Failed estimation: (N_l=", N_l, ", param_est=", param_est, ", var_est=", var_est,
                  ifelse(!is.null(d), paste(", var_d=", var(d_cell)) , ""),
                  ")\n")
      if(warn_on_error) warning(msg)
      next
    }
    t1_mse = sum(N_l*param_est^2)
    t2_var = sum(N_l*var_est)
    N_eff = N_eff + sum(N_l)
    cell_contribs1[cell_i] = t1_mse
    cell_contribs2[cell_i] = t2_var
    if(debug) print(paste("N=", N, "; cell_i=", cell_i, "; N_l=", N_l, "; param_est=", param_est, 
                          "; var_est=", var_est))
  }
  t1_mse = -1/N_eff*sum(cell_contribs1)
  t2_var = (1/N_eff + 1/N_est)*sum(cell_contribs2) 
  val = if(is.null(alpha)) t1_mse + t2_var else alpha*t1_mse + (1-alpha)*t2_var
  if(debug) print(paste("cell sums", val))
  if(!is.finite(val)) stop("Non-finite val")
  return(c(val, N_cell_empty, N_cell_err))
}

get_cell <- function(y, X, d, cell_factor, cell_i, lvls) {
  list[M, m_mode, N, K] = get_sample_type(y, X, d)
  if(is_sep_sample(X)) {
    y_cell = d_cell = X_cell = list()
    N_l = rep(0, M)
    for(m in 1:M) {
      cell_ind = cell_factor[[m]]==lvls[[m]][cell_i]
      y_cell[[m]] = y[[m]][cell_ind]
      d_cell[[m]] = d[[m]][cell_ind]
      X_cell[[m]] = X[[m]][cell_ind, , drop=FALSE]
      N_l[m] = sum(cell_ind)
    }
  }
  else {
    cell_ind = cell_factor==lvls[cell_i]
    N_l = if(M==1) sum(cell_ind) else rep(sum(cell_ind), M)
    y_cell = if(is_vec(y)) y[cell_ind] else y[cell_ind, , drop=FALSE]
    d_cell = if(is_vec(d)) d[cell_ind] else d[cell_ind, , drop=FALSE]
    X_cell = X[cell_ind, , drop=FALSE]
  }
  return(list(y_cell, d_cell, X_cell, N_l))
}

lcl_levels <- function(cell_factor) {
  if(!is.list(cell_factor)) {
    lvls = levels(cell_factor)
    return(list(lvls, length(lvls)))
  }
  lvls = lapply(cell_factor, levels)
  return(list(lvls, length(lvls[[1]])))
}

Param_Est_m <- function(est_plan, y_cell, d_cell, X_cell, sample=sample, ret_var=FALSE, m_mode) {
  if(!is_sep_estimators(m_mode)) { #single estimation
    return(Param_Est(est_plan, y_cell, d_cell, X_cell, sample=sample, ret_var))
  }
  if(m_mode==1){
    if(length(est_plan)!=3) {
      print("ahh")
    }
    if(ret_var) {
      rets = mapply(function(est_plan_s, y_cell_s, d_cell_s, X_cell_s) 
        unlist(Param_Est(est_plan_s, y_cell_s, d_cell_s, X_cell_s, sample=sample, ret_var)), 
        est_plan, y_cell, d_cell, X_cell, SIMPLIFY = TRUE)
      return(list(param_ests=rets[1,], var_ests=rets[2,]))
    }
    rets = mapply(function(est_plan_s, y_cell_s, d_cell_s, X_cell_s) 
      Param_Est(est_plan_s, y_cell_s, d_cell_s, X_cell_s, sample=sample, ret_var)[[1]], 
      est_plan, y_cell, d_cell, X_cell, SIMPLIFY = TRUE)
    return(list(param_ests = rets))
  }
  
  M = ncol(y_cell)
  if(ret_var) {
    rets = sapply(1:M, function(m) unlist(Param_Est(est_plan[[m]], y_cell[,m], d_cell, X_cell, sample=sample, ret_var)))
    return(list(param_ests=rets[1,], var_ests=rets[2,]))
  }
  rets = sapply(1:M, function(m) Param_Est(est_plan[[m]], y_cell[,m], d_cell, X_cell, sample=sample, ret_var)[[1]])
  return(list(param_ests = rets))
}

mse_hat_obj <-function(y, X, d, partition=NULL, cell_factor=NULL, estimator=NULL, debug=FALSE, 
                   warn_on_error=FALSE, est_plan=NULL, sample="trtr", ...) {
  if(is.null(est_plan)) {
    if(is.null(estimator)) est_plan = gen_simple_est_plan(has_d=!is.null(d))
    else est_plan = simple_est(estimator, estimator)
  }
  if(is.null(cell_factor)) {
    cell_factor = get_factor_from_partition(partition, X)
  }
  list[M, m_mode, N, K] = get_sample_type(y, X, d, checks=FALSE)
  list[lvls, n_cells] = lcl_levels(cell_factor)
  cell_contribs = rep(0, n_cells)
  N_eff = 0
  N_cell_empty = 0
  N_cell_error = 0
  for(cell_i in 1:n_cells) {
    list[y_cell, d_cell, X_cell, N_l] <- get_cell(y, X, d, cell_factor, cell_i, lvls)
    if(any(N_l==0)) {
      N_cell_empty = N_cell_empty+1
      next
    }
    list[param_est] = Param_Est_m(est_plan, y_cell, d_cell, X_cell, sample=sample, ret_var=FALSE, m_mode)
    if(!all(is.finite(param_est))) {
      N_cell_error = N_cell_error+1
      msg = paste("Failed estimation: (N_l=", N_l, 
                  ifelse(!is.null(d), paste(", var_d=", var(d_cell)), ""), 
                  ")\n")
      if(warn_on_error) warning(msg)
      next
    }
    t1_mse = sum(N_l*param_est^2)
    N_eff = N_eff + sum(N_l)
    cell_contribs[cell_i] = t1_mse
    if(debug) print(paste("cell_i=", cell_i, "; N_l=", N_l, "; param_est=", param_est))
  }
  val = -1/N_eff*sum(cell_contribs) # Use N_eff to remove from average given errors
  if(debug) print(paste("cell sums", val))
  return(c(val, N_cell_empty, N_cell_error))
}


#' estimate_cell_stats
#'
#' @param y Nx1 matrix of outcome (label/target) data
#' @param X NxK matrix of features (covariates)
#' @param d (Optional) NxP matrix (with colnames) of treatment data. If all equally important they should 
#'          be normalized to have the same variance.
#' @param partition (Optional, need this or cell_factor) partitioning returned from fit_estimate_partition
#' @param cell_factor (Optional, need this or partition) 
#' @param estimator_var (Optional) a function with signature list(param_est, var_est) = function(y, d) 
#'                      (where if no d then can pass in null). If NULL then will choose between built-in 
#'                      mean-estimator and scalar_te_estimator
#' @param est_plan Estimation plan
#' @param alpha Alpha
#' 
#' @return list
#' \item{cell_factor}{Factor with levels for each cell for X. Length N.}
#' \item{stats}{data.frame(cell_i, N_est, param_ests, var_ests, tstats, pval, ci_u, ci_l, p_fwer, p_fdr)}
#' @export
estimate_cell_stats <- function(y, X, d=NULL, partition=NULL, cell_factor=NULL, estimator_var=NULL, 
                                est_plan=NULL, alpha=0.05) {
  list[M, m_mode, N, K] = get_sample_type(y, X, d, checks=TRUE)
  X = ensure_good_X(X)
  
  if(is.null(est_plan)) {
    if(is.null(estimator_var)) est_plan = gen_simple_est_plan(has_d=!is.null(d))
    else est_plan = simple_est(NULL, estimator_var)
  }
  if(is.null(cell_factor)) {
    cell_factor = get_factor_from_partition(partition, X)
  }
  list[lvls, n_cells] = lcl_levels(cell_factor)
  param_ests = matrix(NA, nrow=n_cells, ncol=M)
  var_ests = matrix(NA, nrow=n_cells, ncol=M)
  cell_sizes = matrix(NA, nrow=n_cells, ncol=M)
  for(cell_i in 1:n_cells) {
    list[y_cell, d_cell, X_cell, N_l] <- get_cell(y, X, d, cell_factor, cell_i, lvls)
    cell_sizes[cell_i,] = N_l
    list[param_ests_c, var_ests_c] = Param_Est_m(est_plan, y_cell, d_cell, X_cell, sample="est", ret_var=TRUE, m_mode=m_mode)
    param_ests[cell_i,] = param_ests_c
    var_ests[cell_i,] = var_ests_c
  }
  dofs = t(t(cell_sizes) - get_dofs(est_plan, M, m_mode)) #subtract dofs from each row
  colnames(cell_sizes) = if(M==1) "N_est" else paste("N_est", 1:M, sep="")
  colnames(param_ests) = if(M==1) "param_ests" else paste("param_ests", 1:M, sep="")
  colnames(var_ests) = if(M==1) "var_ests" else paste("var_ests", 1:M, sep="")
  base_df = cbind(data.frame(cell_i=1:n_cells), cell_sizes, param_ests, var_ests)
  list[stat_df, pval] = exp_stats(base_df, param_ests, var_ests, dofs, alpha=alpha, M)
  p_fwer = matrix(p.adjust(pval, "hommel"), ncol=ncol(pval)) #slightly more powerful than "hochberg". Given indep these are better than "bonferroni" and "holm"
  p_fdr = matrix(p.adjust(pval, "BH"), ncol=ncol(pval)) #ours are independent so don't need "BY"
  colnames(p_fwer) = if(M==1) "p_fwer" else paste("p_fwer", 1:M, sep="")
  colnames(p_fdr) = if(M==1) "p_fdr" else paste("p_fdr", 1:M, sep="")
  stat_df = cbind(stat_df, p_fwer, p_fdr)
  return(list(cell_factor=cell_factor, stats=stat_df))
}

get_dofs <- function(est_plan, M, m_mode) {
  if(m_mode==0) 
    return(est_plan$dof)
  
  if(is_sep_estimators(m_mode))
    return(sapply(est_plan, function(plan) plan$dof))
  
  #Only 1 estimator but multiple d, so make sure right length
  dof = est_plan$dof
  if(length(dof)==1) dof=rep(dof, M)
  return(dof)
}

#' est_full_stats
#'
#' @param y y
#' @param d d
#' @param X X
#' @param est_plan  est_plan
#' @param y_es y_es
#' @param d_es d_es
#' @param X_es X_es
#' @param index_tr index_tr
#' @param alpha alpha
#'
#' @return Stats df
#' @export
est_full_stats <- function(y, d, X, est_plan, y_es=NULL, d_es=NULL, X_es=NULL, index_tr=NULL, alpha=0.05) {
  list[M, m_mode, N, K] = get_sample_type(y, X, d, checks=TRUE)
  X = ensure_good_X(X)
  
  if(is.null(y_es)) {
    list[y_tr, y_es, X_tr, X_es, d_tr, d_es, N_est] = split_sample_m(y, X, d, index_tr)
  }
  N_es = nrow_m(X_es, M)
  full_Ns = rbind(N, N_es)
  colnames(full_Ns) = if(M==1) "N_est" else paste("N_est", 1:M, sep="")
  list[full_param_ests_all, full_var_ests_all] = Param_Est_m(est_plan, y, d, X, sample="est", ret_var=TRUE, m_mode=m_mode)
  list[full_param_ests_es, full_var_ests_es] = Param_Est_m(est_plan, y_es, d_es, X_es, sample="est", ret_var=TRUE, m_mode=m_mode)
  M = length(full_param_ests_all)
  full_param_ests = rbind(full_param_ests_all, full_param_ests_es)
  colnames(full_param_ests) = if(M==1) "param_ests" else paste("param_ests", 1:M, sep="")
  full_var_ests = rbind(full_var_ests_all, full_var_ests_es)
  colnames(full_var_ests) = if(M==1) "var_ests" else paste("var_ests", 1:M, sep="")
  base_df = cbind(data.frame(sample=c("all", "est")), full_Ns, full_param_ests, full_var_ests)
  dofs = t(t(full_Ns) - get_dofs(est_plan, M, m_mode)) #subtract dofs from each row
  list[full_stat_df, pval] = exp_stats(base_df, full_param_ests, full_var_ests, dofs, alpha=alpha, M)
  return(full_stat_df)
}

exp_stats <- function(stat_df, param_ests, var_ests, dofs, alpha=0.05, M) {
  tstats = param_ests/sqrt(var_ests)
  colnames(tstats) = if(M==1) "tstats" else paste("tstats", 1:M, sep="")
  t.half.alpha = qt(1-alpha/2, df=dofs)*sqrt(var_ests)
  ci_u = param_ests + t.half.alpha
  colnames(ci_u) = if(M==1) "ci_u" else paste("ci_u", 1:M, sep="")
  ci_l = param_ests - t.half.alpha
  colnames(ci_l) = if(M==1) "ci_l" else paste("ci_l", 1:M, sep="")
  pval = 2*pt(abs(tstats), df=dofs, lower.tail=FALSE)
  colnames(pval) = if(M==1) "pval" else paste("pval", 1:M, sep="")
  stat_df = cbind(stat_df, tstats, ci_u, ci_l, pval)
  #pval_right= pt(tstats, df=dofs, lower.tail=FALSE) #right-tailed. Checking for just a positive effect (H_a is "greater")
  #pval_left = pt(tstats, df=dofs, lower.tail=TRUE) #left-tailed. Checking for just a negative effect (H_a is "less")
  return(list(stat_df, pval))
}

#' predict_te.estimated_partition
#' 
#' Predicted unit-level treatment effect
#'
#' @param obj estimated_partition object
#' @param new_X new X
#'
#' @return predicted treatment effect
#' @export
predict_te.estimated_partition <- function(obj, new_X) {
  #TODO: for mode 1 &2 maybe return a matrix rather than list
  new_X = ensure_good_X(new_X)
  new_X_range = get_X_range(new_X)
  
  cell_factor = get_factor_from_partition(obj$partition, new_X, new_X_range)
  if(obj$M==1) {
    N=nrow(new_X)
    cell_factor_df = data.frame(id=1:N, cell_i = as.integer(cell_factor))
    m_df = merge(cell_factor_df, obj$cell_stats$stats)
    m_df = m_df[order(m_df[["id"]]), ]
    return(m_df[["param_ests"]])
  }
  N = nrow_m(X, M)
  rets = list()
  for(m in 1:M) {
    cell_factor_df = data.frame(id=1:N[m], cell_i = as.integer(cell_factor[[m]]))
    m_df = merge(cell_factor_df, obj$cell_stats$stats)
    m_df = m_df[order(m_df[["id"]]), ]
    rets[[m]] = m_df[["param_ests"]]
  }
  return(rets)
}

get_importance_weights_full_k <- function(k_i, to_compute, X_d, y, d, X_tr, y_tr, d_tr, y_es, X_es, d_es, X_range, pot_break_points, verbosity, ...) {
  if(verbosity>0) cat(paste("Feature weight > ", k_i, "of", length(to_compute),"\n"))
  k = to_compute[k_i]
  X_k = drop_col_k_m(X_d, k)
  X_tr_k = drop_col_k_m(X_tr, k)
  X_es_k = drop_col_k_m(X_es, k)
  main_ret = fit_estimate_partition_int(X_k, y, d, X_tr_k, y_tr, d_tr, y_es, X_es_k, d_es, X_range[-k], pot_break_points=pot_break_points[-k], verbosity=verbosity, ...)
  nk_val = mse_hat_obj(y_es, X_es_k, d=d_es, partition=main_ret$partition, est_plan=main_ret$est_plan, sample="est")[1] #use oos version instead of main_ret$is_obj_val_seq[partition_i]
  return(list(nk_val, main_ret$partition$nsplits_by_dim))
}

# Just use mse_hat as we're working not on the Tr sample, but the est sample
# The ... params are passed to get_importance_weights_full_k -> fit_estimate_partition_int
# There's an undocumented "fast" version. Not very great as assings 0 to any feature not split on
get_importance_weights <- function(X, y, d, X_tr, y_tr, d_tr, y_es, X_es, d_es, X_range, pot_break_points, partition, est_plan, type, verbosity, pr_cl, ...) {
  if(verbosity>0) cat("Feature weights: Started.\n")
  K = length(X_range)
  if(sum(partition$nsplits_by_dim)==0) return(rep(0, K))
  full_val = mse_hat_obj(y_es, X_es, d=d_es, partition = partition, est_plan=est_plan, sample="est")[1]
  
  if(K==1) {
    null_val = mse_hat_obj(y_es, X_es, d=d_es, partition = grid_partition(partition$X_range, partition$varnames), est_plan=est_plan, sample="est")[1]
    if(verbosity>0) cat("Feature weights: Finished.\n")
    return(null_val - full_val)
  }
  
  if("fast"==type) {
    new_vals = rep(0, K)
    factors_by_dim = get_factors_from_partition(partition, X_es)
    for(k in 1:K) {
      if(partition$nsplits_by_dim[k]>0) {
        cell_factor_nk = gen_holdout_interaction_m(factors_by_dim, k, is_sep_sample(X_tr))
        new_vals[k] = mse_hat_obj(y_es, X_es, d=d_es, cell_factor = cell_factor_nk, est_plan=est_plan, sample="est")[1]
      }
    }
    if(verbosity>0) cat("Feature weights: Finished.\n")
    return(new_vals - full_val)
  }
  
  #if("full"==type)
  new_vals = rep(full_val, K)
  to_compute = which(partition$nsplits_by_dim>0)

  params = c(list(to_compute, X_d=X, y=y, d=d, X_tr=X_tr, y_tr=y_tr, d_tr=d_tr, y_es=y_es, X_es=X_es, d_es=d_es, X_range=X_range, pot_break_points=pot_break_points, est_plan=est_plan, verbosity=verbosity-1), 
             list(...))

  rets = my_apply(1:length(to_compute), get_importance_weights_full_k, verbosity==1 || !is.null(pr_cl), pr_cl, params)
  for(k_i in 1:length(to_compute)) {
    k = to_compute[k_i]
    new_vals[k] = rets[[k_i]][[1]]
  }
  if(verbosity>0) cat("Feature weights: Finished.\n")
  return(new_vals - full_val)
}

get_feature_interactions_k12 <- function(ks_i, to_compute, X_d, y, d, X_tr, y_tr, d_tr, y_es, X_es, d_es, X_range, pot_break_points, verbosity, ...) {
  if(verbosity>0) cat(paste("Feature interaction weight > ", ks_i, "of", length(to_compute),"\n"))
  ks = to_compute[[ks_i]]
  X_k = drop_col_k_m(X_d, ks)
  X_tr_k = drop_col_k_m(X_tr, ks)
  X_es_k = drop_col_k_m(X_es, ks)
  main_ret = fit_estimate_partition_int(X_k, y, d, X_tr_k, y_tr, d_tr, y_es, X_es_k, d_es, X_range[-ks], pot_break_points=pot_break_points[-ks], verbosity=verbosity, ...)
  nk_val = mse_hat_obj(y_es, X_es_k, d=d_es, partition=main_ret$partition, est_plan=main_ret$est_plan, sample="est")[1] #use oos version instead of main_ret$is_obj_val_seq[partition_i]
  return(nk_val)
  
}

get_feature_interactions <- function(X, y, d, X_tr, y_tr, d_tr, y_es, X_es, d_es, X_range, pot_break_points, partition, est_plan, verbosity, pr_cl, ...) {
  
  if(verbosity>0) cat("Feature weights: Started.\n")
  K = length(X_range)
  delta_k12 = matrix(as.integer(diag(rep(NA, K))), ncol=K) #dummy for K<3 cases
  if(sum(partition$nsplits_by_dim)==0){ 
    if(verbosity>0) cat("Feature weights: Finished.\nFeature interaction weights: Started.\nFeature interaction interactions: Finished.\n")
    return(list(delta_k=rep(0, K), delta_k12=delta_k12))
  }
  full_val = mse_hat_obj(y_es, X_es, d=d_es, partition = partition, est_plan=est_plan, sample="est")[1]
  
  if(K==1) {
    null_val = mse_hat_obj(y_es, X_es, d=d_es, partition = grid_partition(partition$X_range, partition$varnames), est_plan=est_plan, sample="est")[1]
    if(verbosity>0) cat("Feature weights: Finished.\nFeature interaction weights: Started.\nFeature interaction interactions: Finished.\n")
    return(list(delta_k=null_val - full_val, delta_k12=delta_k12))
  }
  
  #compute the single-removed values (and keep around the nsplits from each new partition)
  new_val_k = rep(full_val, K)
  to_compute_k = which(partition$nsplits_by_dim>0)
  params = c(list(to_compute=to_compute_k, X_d=X, y=y, d=d, X_tr=X_tr, y_tr=y_tr, d_tr=d_tr, y_es=y_es, X_es=X_es, d_es=d_es, X_range=X_range, pot_break_points=pot_break_points, est_plan=est_plan, verbosity=verbosity-1), 
             list(...))
  rets_k = my_apply(1:length(to_compute_k), get_importance_weights_full_k, verbosity==1 || !is.null(pr_cl), pr_cl, params)
  for(k_i in 1:length(to_compute_k)) {
    k = to_compute_k[k_i]
    new_val_k[k] = rets_k[[k_i]][[1]]
  }
  delta_k = new_val_k - full_val
  if(K==2) {
    null_val = mse_hat_obj(y_es, X_es, d=d_es, partition = grid_partition(partition$X_range, partition$varnames), est_plan=est_plan, sample="est")[1]
    delta_k12 = matrix(null_val - full_val, ncol=2) + diag(rep(NA, K))
    if(verbosity>0) cat("Feature weights: Finished.\nFeature interaction weights: Started.\nFeature interaction interactions: Finished.\n")
    return(list(delta_k=delta_k, delta_k12=delta_k12))
  } 
  if(verbosity>0) cat("Feature weights: Finished.\n")
  
  
  #Compute the pair-removed values
  if(verbosity>0) cat("Feature interaction weights: Started.\n")
  new_val_k12 = matrix(full_val, ncol=K, nrow=K)
  to_compute = list()
  for(k1 in 1:(K-1)) {
    if(partition$nsplits_by_dim[k1]==0) {
      new_val_k12[k1,] = new_val_k
      new_val_k12[,k1] = new_val_k
    }
    else {
      k1_i = which(to_compute_k==k1)
      nsplits_by_dim_k1= rets_k[[k1_i]][[2]]
      for(k2 in (k1+1):K) {
        if(nsplits_by_dim_k1[k2-1]==0) { #nsplits_by_dim_k1 is missing k1 so drop k2 back one
          new_val_k12[k1,k2] = new_val_k[k1]
          new_val_k12[k2,k1] = new_val_k[k1]
        }
        else {
          to_compute = c(list(c(k1, k2)), to_compute)
        }
      }
    }
  }
  params = c(list(to_compute=to_compute, X_d=X, y=y, d=d, X_tr=X_tr, y_tr=y_tr, d_tr=d_tr, y_es=y_es, X_es=X_es, d_es=d_es, X_range=X_range, pot_break_points=pot_break_points, est_plan=est_plan, verbosity=verbosity-1), 
             list(...))
  rets_k12 = my_apply(1:length(to_compute), get_feature_interactions_k12, verbosity==1 || !is.null(pr_cl), pr_cl, params)
  for(ks_i in 1:length(to_compute)) {
    k1 = to_compute[[ks_i]][1]
    k2 = to_compute[[ks_i]][2]
    new_val = rets_k12[[ks_i]]
    new_val_k12[k1, k2] = new_val
    new_val_k12[k2, k1] = new_val
  }
  delta_k12 = t(t((new_val_k12 - full_val) - delta_k) - delta_k) + diag(rep(NA, K))
  if(verbosity>0) cat("Feature interaction interactions: Finished.\n")
  return(list(delta_k=delta_k, delta_k12=delta_k12))
}


fit_and_residualize <- function(est_plan, X_tr, y_tr, d_tr, cv_folds, y_es, X_es, d_es, verbosity, dim_cat) {
  est_plan = Fit_InitTr(est_plan, X_tr, y_tr, d_tr, cv_folds, verbosity=verbosity, dim_cat=dim_cat)
  list[y_tr, d_tr] = Do_Residualize(est_plan, y_tr, X_tr, d_tr, sample="tr")
  list[y_es, d_es] = Do_Residualize(est_plan, y_es, X_es, d_es, sample="tr")
  return(list(est_plan, y_tr, d_tr, y_es, d_es))
}

# ... params are passed to fit_partition()
fit_estimate_partition_int <- function(X, y, d, X_tr, y_tr, d_tr, y_es, X_es, d_es, dim_cat, X_range, est_plan, honest, cv_folds, verbosity, M, m_mode, 
                                       alpha, partition_i, ...) {
  K = length(X_range)
  obj_fn = if(honest) emse_hat_obj else mse_hat_obj
  
  list[est_plan, y_tr, d_tr, y_es, d_es] = fit_and_residualize_m(est_plan, X_tr, y_tr, d_tr, cv_folds, y_es, X_es, d_es, m_mode, M, verbosity, dim_cat)
  
  if(verbosity>0) cat("Training partition on training set\n")
  fit_ret = fit_partition(y=y_tr, X=X_tr, d=d_tr, X_aux=X_es, d_aux=d_es, cv_folds=cv_folds, verbosity=verbosity, 
                          X_range=X_range, obj_fn=obj_fn, est_plan=est_plan, valid_fn=NULL, ...)
  list[partition, is_obj_val_seq, complexity_seq, partition_i, partition_seq, split_seq, lambda, cv_foldid] = fit_ret
  
  if(verbosity>0) cat("Estimating cell statistics on estimation set\n")
  cell_stats = estimate_cell_stats(y_es, X_es, d_es, partition, est_plan=est_plan, alpha=alpha)
  
  full_stat_df = est_full_stats(y, d, X, est_plan, y_es=y_es, d_es=d_es, X_es=X_es)
  
  return(list(partition=partition, is_obj_val_seq=is_obj_val_seq, complexity_seq=complexity_seq, partition_i=partition_i, partition_seq=partition_seq, split_seq=split_seq, lambda=lambda, cv_foldid=cv_foldid, cell_stats=cell_stats, full_stat_df=full_stat_df, est_plan=est_plan))
}


#' fit_estimate_partition
#' 
#' Split the data, one one side train/fit the partition and then on the other estimate subgroup effects.
#' With multiple treatment effects (M) there are 3 options (the first two have the same sample across treatment effects).
#'  1) Multiple pairs of (Y_{m},W_{m}). y,X,d are then lists of length M. Each element then has the typical size
#'     The N_m may differ across m. The number of columns of X will be the same across m.
#'  2) Multiple treatments and a single outcome. d is then a NxM matrix.
#'  3) A single treatment and multiple outcomes. y is then a NXM matrix.
#'
#' @param y N vector of outcome (label/target) data
#' @param X NxK matrix of features (covariates). Must be numerical (unordered categorical variables must be 
#'          1-hot encoded.)
#' @param d (Optional) N vector of treatment data.
#' @param max_splits Maximum number of splits even if splits continue to improve OOS fit
#' @param max_cells Maximum number of cells
#' @param min_size Minimum size of cells
#' @param cv_folds Number of CV Folds or foldids. If Multiple effect #3 and using vector, then pass in list of vectors. 
#' @param potential_lambdas potential lambdas to search through in CV
#' @param lambda.1se Use the 1se rule to pick the best lambda
#' @param partition_i Default is NA. Use this to avoid CV automated selection of the partition
#' @param tr_split - can be ratio or vector of indexes. If Multiple effect #3 and using vector then pass in list of vectors. 
#' @param verbosity If >0 prints out progress bar for each split
#' @param pot_break_points NULL or a k-dim list of vectors giving potential split points for non-categorical 
#'                         variables (can put c(0) for categorical). Similar to 'discrete splitting' in 
#'                         CausalTree though their they do separate split-points for treated and controls.
#' @param bucket_min_n Minimum number of observations needed between different split checks for continuous features
#' @param bucket_min_d_var Ensure positive variance of d for the observations between different split checks 
#'                         for continuous features
#' @param honest Whether to use the emse_hat or mse_hat. Use emse for outcome mean. For treatment effect, 
#'               use if want predictive accuracy, but if only want to identify true dimensions of heterogeneity 
#'               then don't use.
#' @param ctrl_method Method for determining additional control variables. Empty ("") for nothing, "all" or "lasso"
#' @param pr_cl Parallel Cluster (If NULL, default, then will be single-processor)
#' @param alpha Default=0.05
#' @param bump_B Number of bump bootstraps
#' @param bump_ratio For bootstraps the ratio of sample size to sample (between 0 and 1, default 1)
#' @param importance_type Options:
#'                        single - (smart) redo full fitting removing each possible dimension
#'                        interaction - (smart) redo full fitting removing each pair of dimensions
#'                         "" - Nothing
#'
#' @return An object with class \code{"estimated_partition"}.
#' \item{partition}{Parition obj defining cuts}
#' \item{cell_stats}{list(cell_factor=cell_factor, stats=stat_df) from estimate_cell_stats() using est sample}
#' \item{importance_weights}{importance_weights}
#' \item{interaction_weights}{interaction_weights}
#' \item{has_d}{has_d}
#' \item{lambda}{lambda used}
#' \item{is_obj_val_seq}{In-sample objective function values for sequence of partitions}
#' \item{complexity_seq}{Complexity #s (# cells-1) for sequence of partitions}
#' \item{partition_i}{Index of Partition selected in sequence}
#' \item{split_seq}{Sequence of splits. Note that split i corresponds to partition i+1}
#' \item{index_tr}{Index of training sample (Size of N)}
#' \item{cv_foldid}{CV foldids for the training sample (Size of N_tr)}
#' \item{varnames}{varnames (or c("X1", "X2",...) if X doesn't have colnames)}
#' \item{honest}{honest target}
#' \item{est_plan}{Estimation plan}
#' \item{full_stat_df}{full_stat_df}
#' @export
fit_estimate_partition <- function(y, X, d=NULL, max_splits=Inf, max_cells=Inf, min_size=3, 
                                   cv_folds=2, potential_lambdas=NULL, lambda.1se=FALSE, partition_i=NA,
                                   tr_split = 0.5, verbosity=0, pot_break_points=NULL, 
                                   bucket_min_n=NA, bucket_min_d_var=FALSE, honest=FALSE,
                                   ctrl_method="", pr_cl=NULL, alpha=0.05, bump_B=0, bump_ratio=1,
                                   importance_type="") {
  
  list[M, m_mode, N, K] = get_sample_type(y, X, d, checks=TRUE)
  if(is_sep_sample(X) && length(tr_split)>1) {
    assert_that(is.list(tr_split) && length(tr_split)==M)
  }
  X = ensure_good_X(X)
  X = update_names_m(X)
  
  dim_cat = which(get_dim_cat_m(X))
  
  #Split the sample
  if(length(tr_split)==1)
    index_tr = gen_split_m(N, tr_split, m_mode==1)
  else 
    index_tr = tr_split
  
  list[y_tr, y_es, X_tr, X_es, d_tr, d_es, N_est] = split_sample_m(y, X, d, index_tr)
  
  
  #Setup est_plan
  if(is.character(ctrl_method)) {
    if(ctrl_method=="") {
      est_plan = gen_simple_est_plan(has_d=!is.null(d))
    } 
    else {
      assert_that(ctrl_method %in% c("all", "LassoCV", "rf"), !is.null(d))
      est_plan = if(ctrl_method=="rf") grid_rf() else lm_X_est(lasso = ctrl_method=="LassoCV")
    }
  }
  else {
    assert_that(inherits(ctrl_method, "Estimator_plan"))
    est_plan = ctrl_method
  } 
  if(is_sep_estimators(m_mode)) {
    est_plan1 = est_plan
    est_plan = list()
    for(m in 1:M) est_plan[[m]] = est_plan1
  }
  
  X_range = get_X_range(X)
  t0 = Sys.time()
  main_ret = fit_estimate_partition_int(X=X, y=y, d=d, X_tr=X_tr, y_tr=y_tr, d_tr=d_tr, y_es=y_es, X_es=X_es, d_es=d_es, 
                                        X_range=X_range, max_splits=max_splits, max_cells=max_cells, min_size=min_size, 
                                        cv_folds=cv_folds, potential_lambdas=potential_lambdas, lambda.1se=lambda.1se, 
                                        partition_i=partition_i, verbosity=verbosity, pot_break_points=pot_break_points, 
                                        bucket_min_n=bucket_min_n, bucket_min_d_var=bucket_min_d_var, honest=honest, 
                                        pr_cl=pr_cl, alpha=alpha, bump_B=bump_B, bump_ratio=bump_ratio, M=M, m_mode=m_mode, 
                                        dim_cat=dim_cat, est_plan=est_plan, N_est=N_est)
  list[partition, is_obj_val_seq, complexity_seq, partition_i, partition_seq, split_seq, lambda, cv_foldid, cell_stats, full_stat_df, est_plan] = main_ret
  
  importance_weights <- interaction_weights <- NULL
  if(importance_type=="interaction") {
    import_ret = get_feature_interactions(X=X, y=y, d=d, X_tr=X_tr, y_tr=y_tr, d_tr=d_tr, y_es=y_es, X_es=X_es, d_es=d_es, X_range=X_range,
                                          pot_break_points=pot_break_points, partition=partition, est_plan=est_plan, verbosity=verbosity, pr_cl=pr_cl, 
                                                  max_splits=max_splits, max_cells=max_cells, min_size=min_size, cv_folds=cv_folds, potential_lambdas=potential_lambdas, lambda.1se=lambda.1se, partition_i=partition_i, 
                                                  bucket_min_n=bucket_min_n, bucket_min_d_var=bucket_min_d_var, honest=honest, alpha=alpha, bump_B=bump_B, 
                                                  bump_ratio=bump_ratio, M=M, m_mode=m_mode, dim_cat=dim_cat, N_est=N_est)
    importance_weights = import_ret$delta_k
    interaction_weights = import_ret$delta_k12
  }
  else if(importance_type %in% c("single", "fast")) {
    importance_weights = get_importance_weights(X=X, y=y, d=d, X_tr=X_tr, y_tr=y_tr, d_tr=d_tr, y_es=y_es, X_es=X_es, d_es=d_es, X_range=X_range,
                                                pot_break_points=pot_break_points, partition=partition, est_plan=est_plan, type=importance_type, verbosity=verbosity, pr_cl=pr_cl, 
                                                max_splits=max_splits, max_cells=max_cells, min_size=min_size, cv_folds=cv_folds, potential_lambdas=potential_lambdas, lambda.1se=lambda.1se, partition_i=partition_i, 
                                                bucket_min_n=bucket_min_n, bucket_min_d_var=bucket_min_d_var, honest=honest, alpha=alpha, bump_B=bump_B, 
                                                bump_ratio=bump_ratio, M=M, m_mode=m_mode, dim_cat=dim_cat, N_est=N_est)
  }
  
  tn = Sys.time()
  td = tn-t0
  if(verbosity>0) cat(paste("Entire Fit-Estimation Duration: ", format(as.numeric(td)), " ", attr(td, "units"), "\n"))

  return(structure(list(partition=partition, 
                        cell_stats=cell_stats, 
                        importance_weights=importance_weights, 
                        interaction_weights=interaction_weights, 
                        has_d=!is.null(d),
                        lambda=lambda, 
                        is_obj_val_seq=is_obj_val_seq,
                        complexity_seq=complexity_seq,
                        partition_i=partition_i,
                        split_seq=split_seq,
                        index_tr=index_tr,
                        cv_foldid=cv_foldid,
                        varnames=names(X),
                        honest=honest,
                        est_plan=est_plan,
                        full_stat_df=full_stat_df,
                        m_mode=m_mode,
                        M=M),
                   class = c("estimated_partition")))
}

#' is.estimated_partition
#'
#' @param x Object
#'
#' @return True if x is an estimated_partition
#' @export
is.estimated_partition <- function(x) {
  inherits(x, "estimated_partition")
} 

#' num_cells.estimated_partition
#'
#' @param obj Estimated Partition
#'
#' @return Number of cells
#' @export
#' @method num_cells estimated_partition
num_cells.estimated_partition <- function(obj) {
  return(num_cells(obj$partition))
}

#' change_complexity
#' 
#' Doesn't update the importance weights
#'
#' @param fit estimated_partition 
#' @param y Nx1 matrix of outcome (label/target) data
#' @param X NxK matrix of features (covariates). Must be numerical (unordered categorical 
#'          variables must be 1-hot encoded.)
#' @param d (Optional) NxP matrix (with colnames) or vector of treatment data. If all equally 
#'          important they should be normalized to have the same variance.
#' @param partition_i partition_i - 1 is the last include in split_seq included in new partition
#'
#' @return updated estimated_partition
#' @export 
change_complexity <- function(fit, y, X, d=NULL, partition_i) {
  #TODO: Refactor checks from fit_estimation_partition and put them here
  list[y_tr, y_es, X_tr, X_es, d_tr, d_es, N_est] = split_sample_m(y, X, d, fit$index_tr)
  
  fit$partition = partition_from_split_seq(fit$split_seq, fit$partition$X_range, 
                                           varnames=fit$partition$varnames, max_include=partition_i-1)
  fit$cell_stats = estimate_cell_stats(y_es, X_es, d_es, fit$partition, est_plan=fit$est_plan)
  
  return(fit)
}


#' get_desc_df.estimated_partition
#'
#' @param obj estimated_partition object 
#' @param do_str If True, use a string like "(a, b]", otherwise have two separate columns with a and b
#' @param drop_unsplit If True, drop columns for variables overwhich the partition did not split
#' @param digits digits Option (default is NULL)
#' @param import_order should we use importance ordering or input ordering (default)
#'
#' @return data.frame
#' @export 
get_desc_df.estimated_partition <- function(obj, do_str=TRUE, drop_unsplit=TRUE, digits=NULL, import_order=FALSE) {
  M = obj$M
  stats = obj$cell_stats$stats[c(F, rep(T,M), rep(T,M), rep(F,M),rep(F,M), rep(F,M), rep(F,M), rep(T,M), rep(F,M), rep(F,M))]
  part_df = get_desc_df.grid_partition(obj$partition, do_str=do_str, drop_unsplit=drop_unsplit, digits=digits)
  
  imp_weights = obj$importance_weights
  if(drop_unsplit) {
    imp_weights = imp_weights[obj$partition$nsplits_by_dim>0]
  }
  if(import_order) part_df = part_df[, order(-1* imp_weights)]
  
  return(cbind(part_df, stats))
}

#' print.estimated_partition
#'
#' @param x estimated_partition object 
#' @param do_str If True, use a string like "(a, b]", otherwise have two separate columns with a and b
#' @param drop_unsplit If True, drop columns for variables overwhich the partition did not split
#' @param digits digits options
#' @param import_order should we use importance ordering or input ordering (default)
#' @param ... Additional arguments. These won't be passed to print.data.frame
#'
#' @return string (and displayed)
#' @export 
#' @method print estimated_partition
print.estimated_partition <- function(x, do_str=TRUE, drop_unsplit=TRUE, digits=NULL, import_order=FALSE, ...) {
  return(print(get_desc_df.estimated_partition(x, do_str, drop_unsplit, digits, import_order=import_order), 
               digits=digits, ...))
}

#predict.estimated_partition <- function(object, X, d=NULL, type="response") {
# TDDO: Have to store y_hat as well as tau_hat 
#}

#libs required and suggested. Use if sourcing directly. 
#If you don't want to use the Rcpp versio of const_vect (`const_vect = const_vectr`) then you can skip Rcpp
#lapply(lib_list, require, character.only = TRUE)
CausalGrid_libs <- function(required=TRUE, suggested=TRUE, load_Rcpp=FALSE) {
  lib_list = c()
  if(required) lib_list = c(lib_list, "caret", "gsubfn", "assertthat")
  if(load_Rcpp) lib_list = c(lib_list, "Rcpp")
  if(suggested) lib_list = c(lib_list, "ggplot2", "glmnet", "gglasso", "parallel", "pbapply", "ranger")
  #Build=Rcpp. Full dev=testthat, knitr, rmarkdown, renv, rprojroot
  return(lib_list)
}

#' any_sign_effect
#' fdr - conservative
#' sim_mom_ineq - Need samples sizes to sufficiently large so that the effects are normally distributed
#'
#' @param obj obj
#' @param check_negative If true, check for a negative. If false, check for positive. 
#' @param method one of c("fdr", "sim_mom_ineq")
#' @param alpha alpha
#' @param n_sim n_sim
#'
#' @return list(are_any= boolean of whether effect is negative)
#' @export
any_sign_effect <- function(obj, check_negative=T, method="fdr", alpha=0.05, n_sim=500) {
  #TODO: could also
  assert_that(method %in% c("fdr", "sim_mom_ineq"))
  if(method=="fdr") {
    assert_that(obj$has_d, alpha>0, alpha<1)
    dofs = obj$cell_stats$stats[["N_est"]] - obj$est_plan$dof
    pval_right= pt(obj$cell_stats$stats$tstats, df=dofs, lower.tail=FALSE) #right-tailed. Checking for just a positive effect (H_a is "greater")
    pval_left = pt(obj$cell_stats$stats$tstats, df=dofs, lower.tail=TRUE) #left-tailed. Checking for just a negative effect (H_a is "less")
    pval1s = if(check_negative) pval_left else pval_right
    pval1s_fdr = p.adjust(pval1s, "BH")
    are_any = sum(pval1s_fdr<alpha) > 0
    return(list(are_any=are_any, pval1s=pval1s, pval1s_fdr=pval1s_fdr))
  }
  else {
    N_cell = nrow(obj$cell_stats$stats)
    te_se = sqrt(obj$cell_stats$stats[["var_ests"]])
    tstat_ext = if(check_negative) min(obj$cell_stats$stats[["tstats"]]) else max(obj$cell_stats$stats[["tstats"]])
    sim_tstat_exts = rep(NA, n_sim)
    for(s in 1:n_sim) {
      sim_te = rnorm(N_cell, mean=0, sd=te_se)
      sim_tstat_exts[s] =  if(check_negative) min(sim_te/te_se) else max(sim_te/te_se)
    }
    if(check_negative) {
      are_any = sum(sim_tstat_exts < quantile(sim_tstat_exts, alpha)) > 0
    }
    else {
      are_any = sum(sim_tstat_exts > quantile(sim_tstat_exts, 1-alpha)) > 0
    }
  }
  return(list(are_any=are_any))
}
