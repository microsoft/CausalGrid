# Utils

# Generics ------------------


#' Return number of cells for the object
#' 
#' Returns the number of cells for the an object
#'
#' @param obj Object
#'
#' @return Number of cells in partition (at least 1)
#' @export
num_cells <- function(obj) {
  UseMethod("num_cells", obj)
} 

# General Utils ----------------

#handles vectors and 2D structures
row_sample <- function(data, sample) {
  if(is.null(data)) return(NULL)
  if(is.null(ncol(data)))
    return(data[sample])  
  return(data[sample, , drop=FALSE])
}

is_vec <- function(X) {
  return(is.null(ncol(X)) || ncol(X)==1)
}

get_dim_cat <- function(X) {
  if(is.data.frame(X)) {
    return(sapply(X, is.factor) & !sapply(X, is.ordered))
  }
  return(rep(F, ncol(X)))
}

update_names <- function(X) {
  if(is.null(colnames(X))){
    colnames(X) = paste("X", 1:ncol(X), sep="")
  } 
  return(X)
}

#Note that the verbosity passed in here could be different than a member of the params
my_apply <- function(X, fn_k, apply_verbosity, pr_cl, params) {
  K = length(X)
  if(requireNamespace("pbapply", quietly = TRUE) & (apply_verbosity>0) & (is.null(pr_cl) || length(pr_cl)<K)) {
    rets = do.call(pbapply::pblapply, c(list(X, fn_k), params, list(cl=pr_cl)))
  }
  else if(!is.null(pr_cl)) {
    rets = do.call(parallel::parLapply, c(list(pr_cl, X, fn_k), params))
  }
  else {
    rets = do.call(lapply, c(list(X, fn_k), params))
  }
  return(rets)
}


replace_k_factor <- function(base_facts, k, new_fact) {
  base_facts[[k]] = new_fact
  return(base_facts)
}

is_factor_dim_k <- function(X, k) {
  return(is.factor(X[, k]))
}


#Standard way to check if vector is constant is const_vectr(), but is O(n).
#Checking element-by-element would often be faster, but this is inefficient in R
#and faster in C. const_vect1() and const_vect2() were two versions (first using
#'inline', second just Rcpp), but couldn't get to work in building a package.
#The Rcpp version is now in a separate file.

const_vectr <- function(x) {
  if(length(x)==0) return(TRUE)
  r = range(x)
  return(r[1]==r[2])
}

# Fold utils --------------------------

gen_folds <- function(y, nfolds) {
  assert_that(nfolds>1, msg="Need nfolds>1")
  idxOut_new  = caret::createFolds(y, k=nfolds, list=TRUE)
  
  idx_new = list()
  for(f in 1:nfolds) {
    idx_new[[f]] = sort(unlist(idxOut_new[-f], use.names=FALSE))
  }
  
  foldids = foldlists_to_foldids(idxOut_new)
  
  return(list(foldids=foldids, index=idx_new, indexOut=idxOut_new))
  
}

factor_from_idxs <-function(N, nfolds, indexOut) {
  folds = vector("numeric", N)
  for(f in 1:nfolds) {
    folds[indexOut[[f]]] = f
  }
  folds_f = as.factor(folds)
  return(folds_f)
}


foldlists_to_foldids <- function(indexOut) {
  nfolds = length(indexOut)
  N = sum(sapply(indexOut, length))
  foldids = rep(0, N)
  for(f in 1:nfolds) {
    foldids[indexOut[[f]]] = f
  }
  return(foldids)
}

foldids_to_foldlists <- function(foldids, nfolds) {
  index = list()
  indexOut = list()
  for(f in 1:nfolds){
    index[[f]] = which(foldids!=f)
    indexOut[[f]] = which(foldids==f)
  }
  return(list(index=index, indexOut=indexOut))
}


DS.SINGLE = 0
DS.MULTI_SAMPLE = 1
DS.MULTI_D = 2
DS.MULTI_Y = 3

expand_fold_info <- function(y, cv_folds, m_mode=DS.SINGLE) {
  if(length(cv_folds)==1) {
    nfolds = cv_folds
    folds_ret = gen_folds_m(y, nfolds, m_mode)
    foldids = foldlists_to_foldids_m(folds_ret, m_mode)
  }
  else {
    foldids = cv_folds
    nfolds = if(m_mode!=DS.MULTI_SAMPLE) max(foldids) else max(foldids[[1]])
    folds_ret = foldids_to_foldlists_m(foldids, nfolds, m_mode)
  }
  
  return(list(nfolds, folds_ret, foldids))
}

expand_bump_samples <- function(bump_samples, bump_ratio, N, m_mode) {
  if(length(bump_samples)==1) {
    bump_B = bump_samples
    bump_samples <- lapply(seq_len(bump_B), function(b){sample_m(bump_ratio, N, m_mode==DS.MULTI_SAMPLE)})  
  }
  return(bump_samples)
}

# Multi-sample utils ----------------------
is_sep_sample <- function(X) {
  return(is.list(X) & !is.data.frame(X))
}

is_sep_estimators <- function(m_mode) {
  return(m_mode==DS.MULTI_SAMPLE || m_mode==DS.MULTI_Y)
}

sum_m <- function(data, M_mult) {
  if(!M_mult) return(sum(data))
  return(sapply(data, sum))
}

ensure_good_X <- function(X) {
  if(is_sep_sample(X)) {
    return(lapply(X, ensure_good_X))
  }
  
  if (is.matrix(X)) {
    are_equal(mode(X), "numeric")
  }
  else {
    assert_that(is.data.frame(X), msg="X is not a matrix or data.frame")
    if (inherits(X, "tbl")) X <- as.data.frame(X) # tibble's return tibble (rather than vector) for X[,k], making is.factor(X[,k]) and others fail. Could switch to doing X[[k]] for df-like objects
    for (k in seq_len(ncol(X))) are_equal(mode(X[[k]]), "numeric")
  }
  assert_that(ncol(X) >= 1, msg="X has no columns.")
  return(X)
}


get_sample_type <- function(y, X, d=NULL, checks=FALSE) {
  if(is_sep_sample(X)) { #Different samples
    m_mode=DS.MULTI_SAMPLE
    M = length(X)
    N = sapply(X, nrow)
    K = ncol(X[[1]])
    
    if(checks) {
      check_list_dims <- function(new_type) {
        assert_that(is.list(new_type), length(new_type)==M, msg="Separate samples, but aux param isn't list of same length")
        for(m in 1:M) assert_that(length(new_type[[m]])==N[[m]], msg="Separate samples, but aux param's list elements aren't the right length.")
      }
      check_list_dims(y)
      if(!is.null(d)) check_list_dims(d)
      
      for(m in 1:M) {
        assert_that(ncol(X[[m]])==K, msg="Separate samples, but X's don't all have the same number of columns.")
      } 
    }
    
  }
  else { #Same sample
    N = nrow(X)
    K = ncol(X)
    
    if(!is.null(d) && is.matrix(d) && ncol(d)>1) {
      m_mode= DS.MULTI_D
      M = ncol(d)
      if(checks){
        assert_that(!inherits(d, "tbl"), msg="d not allowed to be a tibble") #TODO: Could silently conver
        assert_that(nrow(d)==N, length(y)==N, msg="d and N don't have the right number of rows.")
      }
    }
    else if(!is.null(d) && is.matrix(y) && ncol(y)>1) {
      m_mode= DS.MULTI_Y
      M = ncol(y)
      N = nrow(X)
      if(checks){
        assert_that(!inherits(y, "tbl"), msg="d not allowed to be a tibble") #TODO: Could silently conver
        assert_that(is.null(d) || length(d)==N, nrow(y)==N, msg="d and N don't have the right number of rows.") 
      }
    }
    else {
      m_mode= DS.SINGLE
      M=1    
      if(checks)
        assert_that(is.null(d) || length(d)==N, length(y)==N, msg="d and N don't have the right number of rows.")
    }
    
    if(M>1) N= rep(N, M)
  }
  return(list(M, m_mode, N, K))
}

check_M_K <- function(M, m_mode, K, X_aux, d_aux) {
  if(m_mode==DS.MULTI_SAMPLE) {
    assert_that(length(X_aux)==M, is.null(d_aux) || length(d_aux)==M, msg="Separate samples, but X_aux or d_aux don't have the right structure.")
    for(m in 1:M) assert_that(ncol(X_aux[[m]])==K, msg="Separate samples, but an element of X_aux doesn't have the right number of columns.")
  }
  else {
    assert_that(ncol(X_aux)==K, msg="X_aux doesn't have the right number of columns.")
    if(m_mode==DS.MULTI_D) assert_that(ncol(d_aux)==M, msg="MULTI_D case but not the right number of cols in d.")
  }
} 

# Multi-sample wrappers --------------------


valid_partition_m <- function(M_mult, valid_fn, cell_factor, d, cell_factor_aux=NULL, d_aux=NULL, min_size=0) {
  if(M_mult) {
    for(m in 1:length(cell_factor)) {
      v_ret = valid_fn(cell_factor[[m]], d[[m]], cell_factor_aux[[m]], d_aux[[m]], min_size)
      if(v_ret$fail) return(v_ret)
    }
    return(list(fail=FALSE))
  }
  return(valid_fn(cell_factor, d, cell_factor_aux, d_aux, min_size))
}

#length of N is M (so even if same sample)
nrow_m <- function(X, M) {
  if(is_sep_sample(X)) {
    return(sapply(X, nrow))
  }
  N = nrow(X)
  if(M==1) return(N)
  
  return(rep(N, M))
}

# Return M-list if mode_m==1 else sample
sample_m <- function(ratio, N, M_mult) {
  if(!M_mult) {
    if(length(N)>1) N=N[1] #for modes 2 & 3
    return(sample(N, N*ratio, replace=TRUE))
  }
  return(lapply(N, function(N_s) sample(N_s, N_s*ratio, replace=TRUE)))
}

#assumes separate samples if m_mode==DS.MULTI_SAMPLE
subsample_m <- function(y, X, d, sample) {
  M_mult = is_sep_sample(X)
  if(!M_mult) {
    return(list(row_sample(y,sample), X[sample,,drop=FALSE], row_sample(d,sample)))
  }
  return(list(mapply(function(y_s, sample_s) y_s[sample_s], y, sample, SIMPLIFY=FALSE),
              mapply(function(X_s, sample_s) X_s[sample_s,,drop=FALSE], X, sample, SIMPLIFY=FALSE),
              mapply(function(d_s, sample_s) d_s[sample_s], d, sample, SIMPLIFY=FALSE)))
}


gen_split_m <- function(N, tr_split, M_mult) {
  if(!M_mult) {
    if(length(N>1)) N=N[1] # mode 2 & 3
    return(base::sample(N, tr_split*N))
  }
  return(lapply(N, function(n) base::sample(n, tr_split*n)))
}

split_sample <- function(y, X, d, index_tr) {
  list[y_tr, y_es] = list(row_sample(y, index_tr), row_sample(y, -index_tr))
  list[d_tr, d_es] = list(row_sample(d, index_tr), row_sample(d, -index_tr))
  X_tr = X[index_tr, , drop=FALSE]
  X_es = X[-index_tr, , drop=FALSE]
  N_est = nrow(X_es)

  return(list(y_tr, y_es, X_tr, X_es, d_tr, d_es, N_est))
}

split_sample_m <- function(y, X, d, index_tr) {
  if(!is_sep_sample(X)) {
    return(split_sample(y, X, d, index_tr))
  }
  else {
    y_tr = y_es = X_tr = X_es = d_tr = d_es = list()
    N_est = rep(0, length(X))
    for(m in 1:length(X))
      list[y_tr[[m]], y_es[[m]], X_tr[[m]], X_es[[m]], d_tr[[m]], d_es[[m]], N_est[m]] = split_sample(y[[m]], X[[m]], d[[m]], index_tr[[m]])
    N_est = sapply(X_es, nrow)
  }
  return(list(y_tr, y_es, X_tr, X_es, d_tr, d_es, N_est))
}

gen_folds_m <-function(y, nfolds, m_mode) {
  if(m_mode!=DS.MULTI_SAMPLE) {
    if(is.list(y)) y = y[[1]]
    if(is.matrix(y)) y = y[,1]
    
    return(gen_folds(y, nfolds))
  }
  M = length(y)
  return(lapply(1:M, function(m) gen_folds(y[[m]], nfolds)))
}

foldlists_to_foldids_m <- function(folds_ret, m_mode) {
  if(m_mode!=DS.MULTI_SAMPLE) return(foldlists_to_foldids(folds_ret$indexOut))
  
  return(lapply(folds_ret, function(f_ret) foldlists_to_foldids(f_ret$indexOut)))
}

foldids_to_foldlists_m <- function(foldids, nfolds, m_mode) {
  if(m_mode!=DS.MULTI_SAMPLE) return(foldids_to_foldlists(foldids, nfolds))
  
  return(lapply(foldids, function(f_ids) foldids_to_foldlists(f_ids, nfolds)))
}

split_sample_folds_m <- function(y, X, d, folds_ret, f) {
  if(!is_sep_sample(X)) {
    list[y_f_tr, y_f_cv] = list(row_sample(y, folds_ret$index[[f]]), row_sample(y, folds_ret$indexOut[[f]]))
    list[d_f_tr, d_f_cv] = list(row_sample(d, folds_ret$index[[f]]), row_sample(d, folds_ret$indexOut[[f]]))
    X_f_tr = X[folds_ret$index[[f]], , drop=FALSE]
    X_f_cv = X[folds_ret$indexOut[[f]], , drop=FALSE]
  }
  else {
    y_f_tr = y_f_cv = X_f_tr = X_f_cv = d_f_tr = d_f_cv = list()
    for(m in 1:length(X))
      list[y_f_tr[[m]], y_f_cv[[m]], X_f_tr[[m]], X_f_cv[[m]], d_f_tr[[m]], d_f_cv[[m]]] = split_sample_folds_m(y[[m]], X[[m]], d[[m]], folds_ret[[m]], f)
  }
  return(list(y_f_tr, y_f_cv, X_f_tr, X_f_cv, d_f_tr, d_f_cv))
}

fit_and_residualize_m <- function(est_plan, X_tr, y_tr, d_tr, cv_folds, y_es, X_es, d_es, m_mode, M, verbosity, dim_cat) {
  if(!is_sep_estimators(m_mode))
    return(fit_and_residualize(est_plan, X_tr, y_tr, d_tr, cv_folds, y_es, X_es, d_es, verbosity, dim_cat))
  
  if(m_mode==DS.MULTI_SAMPLE) {
    for(m in 1:M)
      list[est_plan[[m]], y_tr[[m]], d_tr[[m]], y_es[[m]], d_es[[m]]] = fit_and_residualize(est_plan[[m]], X_tr[[m]], y_tr[[m]], d_tr[[m]], cv_folds[[m]], y_es[[m]], X_es[[m]], d_es[[m]], verbosity, dim_cat)
    return(list(est_plan, y_tr, d_tr, y_es, d_es))
  }
  
  #m_mode==DS.MULTI_Y
  #We overwrite the d's
  for(m in 1:M)
    list[est_plan[[m]], y_tr[,m], d_tr, y_es[,m], d_es] = fit_and_residualize(est_plan[[m]], X_tr, y_tr[,m], d_tr, cv_folds, y_es[,m], X_es, d_es, verbosity, dim_cat)
  return(list(est_plan, y_tr, d_tr, y_es, d_es))
}

Param_Est_m <- function(est_plan, y_cell, d_cell, X_cell, sample=sample, ret_var=FALSE, m_mode) {
  if(!is_sep_estimators(m_mode)) { #single estimation
    return(est_params(est_plan, y_cell, d_cell, X_cell, sample=sample, ret_var))
  }
  if(m_mode==DS.MULTI_SAMPLE){
    if(ret_var) {
      rets = mapply(function(est_plan_s, y_cell_s, d_cell_s, X_cell_s) 
        unlist(est_params(est_plan_s, y_cell_s, d_cell_s, X_cell_s, sample=sample, ret_var)), 
        est_plan, y_cell, d_cell, X_cell, SIMPLIFY = TRUE)
      return(list(param_ests=rets[1,], var_ests=rets[2,]))
    }
    rets = mapply(function(est_plan_s, y_cell_s, d_cell_s, X_cell_s) 
      est_params(est_plan_s, y_cell_s, d_cell_s, X_cell_s, sample=sample, ret_var)[[1]], 
      est_plan, y_cell, d_cell, X_cell, SIMPLIFY = TRUE)
    return(list(param_ests = rets))
  }
  
  #m_mode==DS.MULTI_Y
  M = ncol(y_cell)
  if(ret_var) {
    rets = sapply(1:M, function(m) unlist(est_params(est_plan[[m]], y_cell[,m], d_cell, X_cell, sample=sample, ret_var)))
    return(list(param_ests=rets[1,], var_ests=rets[2,]))
  }
  rets = sapply(1:M, function(m) est_params(est_plan[[m]], y_cell[,m], d_cell, X_cell, sample=sample, ret_var)[[1]])
  return(list(param_ests = rets))
}


interaction_m <- function(facts, M_mult=FALSE, drop=FALSE) {
  if(!M_mult) {
    return(interaction(facts, drop=drop))
  }
  return(lapply(facts, function(f) interaction(f, drop=drop)))
}

interaction2_m <- function(f1, f2, M_mult=FALSE, drop=FALSE) {
  if(!M_mult) {
    return(interaction(f1, f2, drop=drop))
  }
  return(mapply(function(f1_s, f2_s) interaction(f1_s, f2_s, drop=drop), f1, f2, SIMPLIFY=FALSE))
}


gen_holdout_interaction_m <- function(factors_by_dim, k, M_mult) {
  if(!M_mult) 
    return(gen_holdout_interaction(factors_by_dim, k))
  
  return(lapply(factors_by_dim, function(f_by_dim) gen_holdout_interaction(f_by_dim, k)))
}

is_factor_dim_k_m <- function(X, k, M_mult) {
  if(!M_mult)
    return(is_factor_dim_k(X, k))
  return(is_factor_dim_k(X[[1]], k))
}

droplevels_m <- function(factor, M_mult) {
  if(!M_mult) return(droplevels(factor))
  return(lapply(factor, droplevels))
}

apply_mask_m <- function(data, mask, M_mult) {
  if(is.null(data)) return(NULL)
  if(!M_mult) return(row_sample(data, mask))
  return(mapply(function(data_s, mask_s) row_sample(data_s, mask_s), data, mask, SIMPLIFY=FALSE))
}

any_const_m <- function(d_shifted, shifted, shifted_cell_factor_nk, m_mode) {
  if(m_mode==DS.SINGLE || m_mode==DS.MULTI_Y)
    return(any(by(d_shifted, shifted_cell_factor_nk, FUN=const_vect)))
  if(m_mode==DS.MULTI_SAMPLE)
    return( any(mapply(function(d_shifted_s, shifted_cell_factor_nk_s)
      any(by(d_shifted_s, shifted_cell_factor_nk_s, FUN=const_vect))
      , d_shifted, shifted_cell_factor_nk ))  )
  #m_mode==DS.MULTI_D
  return( any(apply(d_shifted, 2, function(d_shifted_s)  any(by(d_shifted_s, shifted_cell_factor_nk, FUN=const_vect)) ))  )
}
gen_cat_window_mask_m <- function(X, k, window) {
  if(is.null(X)) return(NULL)
  M_mult = is_sep_sample(X)
  if(!M_mult) return(X[, k] %in% window)
  return(lapply(X, function(X_s) X_s[, k] %in% window))
} 
gen_cat_win_split_cond_m <- function(X, win_mask, k, win_split_val) {
  M_mult = is_sep_sample(X)
  if(!M_mult)
    return(factor(X[win_mask, k] %in% win_split_val, levels=c(FALSE, TRUE)))
  return(mapply(function(X_s, win_mask_s) factor(X_s[win_mask_s, k] %in% win_split_val, levels=c(FALSE, TRUE)), X, win_mask, SIMPLIFY=FALSE))
} 

gen_cont_window_mask_m <- function(X, k, win_LB, win_UB) {
  if(is.null(X)) return(NULL)
  M_mult = is_sep_sample(X)
  if(!M_mult) return(win_LB<X[, k] & X[, k]<=win_UB)
  return(lapply(X, function(X_s) win_LB<X_s[, k] & X_s[, k]<=win_UB))
}

gen_cont_win_split_cond_m <- function(X, win_mask, k, X_k_cut) {
  M_mult = is_sep_sample(X)
  if(!M_mult)
    return(factor(X[win_mask, k] <= X_k_cut, levels=c(FALSE, TRUE)))
  return(mapply(function(X_s, win_mask_s) factor(X_s[win_mask_s, k] <= X_k_cut, levels=c(FALSE, TRUE)), X, win_mask, SIMPLIFY=FALSE))
} 

replace_k_factor_m <- function(base_facts, k, new_fact, M_mult) {
  if(!M_mult) {
    return(replace_k_factor(base_facts, k, new_fact))
  }
  for(m in 1:length(base_facts)) {
    base_facts[[m]][[k]] = new_fact[[m]]
  }
  return(base_facts)
}

update_names_m <- function(X) {
  if(!is_sep_sample(X))
    return(update_names(X))
  
  for(m in 1:length(X)) {
    X[[m]] = update_names(X[[m]])
  }
  return(X)
}

get_dim_cat_m <- function(X) {
  if(!is_sep_sample(X))
    return(get_dim_cat(X))
  
  return(get_dim_cat(X[[1]]))
}

drop_col_k_m <- function(X, k) {
  if(!is_sep_sample(X))
    return(X[,-k, drop=FALSE])
  
  for(m in 1:length(X)) {
    X[[m]] = X[[m]][,-k, drop=FALSE]
  }
  return(X)
}


get_col_k_m <- function(X, k) {
  if(!is_sep_sample(X))
    return(X[,k, drop=FALSE])
  
  for(m in 1:length(X)) {
    X[[m]] = X[[m]][,k, drop=FALSE]
  }
  return(X)
}


