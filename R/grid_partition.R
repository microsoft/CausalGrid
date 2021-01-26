# agnostic to objective function or data splits


# grid_partition -----------------

#' Grid Partition
#' 
#' A \code{\link{grid_partition}} is defines a grid over a feature-space. It can be built by composing \code{\link{partition_split}}s. 
#' 
#' The partition is typically built by a search algorithm as \code{\link{fit_partition}}.
#' @name GridPartition
NULL
#> NULL

#' Create a null grid_partition
#' 
#' Create a empty partition. Splits can be added using \code{\link{add_partition_split}}.
#' Information about a split can be retrieved using \code{\link{num_cells}}, \code{\link{get_desc_df.grid_partition}} and \code{\link{print}}
#' With data, one can determine the cell for each observation using \code{\link{predict}}
#'
#' @param X_range Such as from \code{\link{get_X_range}}
#' @param varnames Names of the X-variables
#'
#' @return Grid Partition
#' @export
grid_partition <- function(X_range, varnames=NULL) {
  K = length(X_range)
  s_by_dim = vector("list", length=K) #splits_by_dim(s_seq) #stores Xk_val's
  dim_cat = c()
  for (k in 1:K) {
    if(mode(X_range[[k]])=="character") { 
      dim_cat = c(dim_cat, k) 
      s_by_dim[[k]] = list()
    }
    else { 
      s_by_dim[[k]] = vector("numeric")
    }
  }
  nsplits_by_dim = rep(0, K)
  
  return(structure(list(s_by_dim = s_by_dim, nsplits_by_dim = nsplits_by_dim, varnames=varnames, dim_cat=dim_cat, 
                        X_range=X_range), class = c("grid_partition")))  
}

#' Is grid_partition
#' 
#' Test whether an object is an \code{grid_function}
#'
#' @param x an R object
#'
#' @return True if x is a grid_partition
#' @export
#' @describeIn grid_partition is grid_partition
is_grid_partition <- function(x) {
  inherits(x, "grid_partition")
} 


#' Get X_range
#' 
#' Gets the "range" of each variable in X. For numeric variables this is (min, max).
#' For factors this means vector of levels.  
#'
#' @param X data
#'
#' @return list of length K with each element being the "range" along that dimension
#' @export
get_X_range <- function(X) {
  if(is_sep_sample(X))
    X = do.call("rbind", X)
  if(is.matrix(X)) {
    are_equal(mode(X), "numeric")
  }
  else {
    assert_that(is.data.frame(X), msg="X is not a matrix or data.frame")
    if(inherits(X, "tbl")) X = as.data.frame(X) #tibble's return tibble (rather than vector) for X[,k], making is.factor(X[,k]) and others fail. Could switch to doing X[[k]] for df-like objects
    for(k in seq_len(ncol(X))) are_equal(mode(X[[k]]), "numeric")
  }
  assert_that(ncol(X)>=1, msg="X has no columns")
  
  X_range = list()
  K = ncol(X)
  for(k in 1:K) {
    X_k = X[, k]
    X_range[[k]] = if(is.factor(X_k)) levels(X_k) else range(X_k) #c(min, max)
  }
  return(X_range)
}


#' Get factor describing cell number fo each observation
#' 
#' Note that currently if X has values more extreme (e.g., for numeric or factor levels ) than was used to generate the partition
#' then we will return NA unless you provide and updated X_range.
#'
#' @param obj partition
#' @param X X data or list of X
#' @param X_range (Optional) overrides the partition$X_range
#'
#' @return Factor
#' @export
predict.grid_partition <- function(obj, X, X_range=NULL) {
  facts = get_factors_from_partition(obj, X, X_range=X_range)
  return(interaction_m(facts, is_sep_sample(X)))
}


#' @describeIn num_cells grid_partition
#' @export
num_cells.grid_partition <- function(obj) {
  return(prod(obj$nsplits_by_dim+1))
}

#' Print grid_partition
#' 
#' Prints a data.frame with options
#'
#' @param x partition object
#' @param do_str If True, use a string like "(a, b]", otherwise have two separate columns with a and b
#' @param drop_unsplit If True, drop columns for variables overwhich the partition did not split
#' @param digits digits Option
#' @param ... Additional arguments. Passed to data.frame
#'
#' @return string (and displayed)
#' @export
print.grid_partition <- function(x, do_str=TRUE, drop_unsplit=TRUE, digits=NULL, ...) {
  #To check: digits
  assert_that(is.flag(do_str), is.flag(drop_unsplit), msg="One of do_str or drop_unsplit are not flags")
  return(print(get_desc_df.grid_partition(x, do_str=do_str, drop_unsplit=drop_unsplit, digits=digits), 
               digits=digits, ...))
}


#' Get descriptive data.frame for grid_partition
#' 
#' A dataset with rows for each cell and columns defining partitioning
#'
#' @param partition Partition
#' @param cont_bounds_inf If True, will put continuous bounds as -Inf/Inf. Otherwise will use X_range bounds
#' @param do_str If True, use a string like "(a, b]", otherwise have two separate columns with a and b
#' @param drop_unsplit If True, drop columns for variables overwhich the partition did not split
#' @param digits digits option
#' @param unsplit_cat_star if we don't split on a categorical var, should we show as "*" (otherwise list all levels)
#'
#' @return data.frame
#' @export
get_desc_df.grid_partition <- function(partition, cont_bounds_inf=TRUE, do_str=FALSE, drop_unsplit=FALSE, 
                                       digits=NULL, unsplit_cat_star=TRUE) {
  #To check: digits
  assert_that(is.flag(cont_bounds_inf), is.flag(do_str), is.flag(drop_unsplit), is.flag(unsplit_cat_star), msg="One (cont_bounds_inf, do_str, drop_unsplit, unsplit_cat_star)of are not flags.")
  # A split at x_k means that we split to those <= and >
  
  n_segs = partition$nsplits_by_dim+1
  n_cells = prod(n_segs)
  
  if(n_cells==1 & drop_unsplit) return(as.data.frame(matrix(NA, nrow=1, ncol=0)))
  
  #Old code
  #library(tidyverse)
  #desc_df = data.frame(labels=levels(grid_fit$cell_stats$cell_factor), 
  #                     stringsAsFactors = FALSE) %>% separate(labels, names(X), "(?<=]).(?=[(])", PERL=TRUE)
  
  K = length(partition$nsplits_by_dim)
  X_range = partition$X_range
  if(cont_bounds_inf) {
    for(k in 1:K) {
      if(!k %in% partition$dim_cat) X_range[[k]] = c(-Inf, Inf)
    }
  }
  colnames=partition$varnames
  if(is.null(colnames)) colnames = paste("X", 1:K, sep="")
  
  list_of_windows = list()
  for(k in 1:K) {
    list_of_windows[[k]] = if(k %in% partition$dim_cat) get_windows_cat(partition$s_by_dim[[k]], X_range[[k]]) else get_window_cont(partition$s_by_dim[[k]], X_range[[k]])
  }
  
  format_cell_cat <- function(win, unsplit_cat_star, n_tot_dim, sep=", ") {
    if(unsplit_cat_star && n_tot_dim==1) return("*")
    return(paste(win, collapse=sep))
  }
  format_cell_cont <- function(win) {
    if(is.infinite(win[1]) && is.infinite(win[2])) return("*")
    if(is.infinite(win[1])) return(paste0("<=", format(win[2], digits=digits)))
    if(is.infinite(win[2])) return(paste0(">", format(win[1], digits=digits)))
    return(paste0("(", format(win[1], digits=digits), ", ", format(win[2], digits=digits), "]"))
  }
  
  raw_data = data.frame(row.names=1:n_cells)
  str_data = data.frame(row.names=1:n_cells)
  for(k in 1:K) {
    raw_data_k = list()
    str_data_k = c()
    for(cell_i in 1:n_cells) {
      segment_indexes = segment_indexes_from_cell_i(cell_i, n_segs)
      win = list_of_windows[[k]][[segment_indexes[k]]]
      raw_data_k[[cell_i]] = win
      str_data_k[cell_i] = if(k %in% partition$dim_cat) format_cell_cat(win, unsplit_cat_star, length(list_of_windows[[k]])) else format_cell_cont(win)
    }
    raw_data[[colnames[k]]] = cbind(raw_data_k) #make a list-column: https://stackoverflow.com/a/51308306
    str_data[[colnames[k]]] = factor(str_data_k, levels=unique(str_data_k)) #will be in low-high order
  }
  desc_df = if(do_str) str_data else raw_data
  if(drop_unsplit) desc_df = desc_df[n_segs>1]
  
  
  return(desc_df)
}

#' Adds partition_split to grid_partition
#' 
#' Update the partition with an additional split.
#'
#' @param obj Grid Partition object
#' @param s Partition Split object
#'
#' @return updated Grid Partition
#' @export
add_partition_split <- function(obj, s) {
  k = s[[1]]
  X_k_cut = s[[2]]
  
  if(k %in% obj$dim_cat) obj$s_by_dim[[k]][[obj$nsplits_by_dim[k]+1]] = X_k_cut
  else obj$s_by_dim[[k]] = sort(c(X_k_cut, obj$s_by_dim[[k]]))
  obj$nsplits_by_dim[k] = obj$nsplits_by_dim[k]+1
  
  return(obj)
}

get_factors_from_splits_dim <- function(X_k, X_k_range, s_by_dim_k) {
  if(mode(X_k_range)=="character") {
    windows = get_windows_cat(s_by_dim_k, X_k_range)
    fac = X_k
    new_name_map = levels(fac)
    new_names = c()
    for(window in windows) {
      new_name = if(length(window)>1) paste0("{", paste(window, collapse=","), "}") else window[1]
      new_name_map[levels(fac) %in% window] = new_name
      new_names = c(new_names, new_name)
    }
    levels(fac) <- new_name_map
    fac = factor(fac, new_names)
  }
  else {
    bottom_break = X_k_range[1]
    top_break = X_k_range[2]
    #if(nsplits_by_dim_k>0) {
    #bottom_split = s_by_dim_k[1]
    #if(bottom_split==bottom_break)
    bottom_break = bottom_break-1 #not needed
    #}
    top_break = top_break+1
    breaks = c(bottom_break, s_by_dim_k, top_break)
    fac = cut(X_k, breaks, labels=NULL, include.lower=TRUE) #right=FALSE makes [a,b) segments. labels=FALSE makes just numeric vector
  }
  return(fac)
}

get_factors_from_splits_dim_m <- function(X, X_k_range, s_by_dim_k, k) {
  M_mult = is_sep_sample(X)
  if(!M_mult)
    return(get_factors_from_splits_dim(X[,k], X_k_range, s_by_dim_k))
  return(lapply(X, function(X_s) get_factors_from_splits_dim(X_s[,k], X_k_range, s_by_dim_k)))
}

# for a continuous variables, splits are just values
# for a factor variable, a split is a vector of levels (strings)

dummy_X_range <- function(K) {
  X_range = list()
  for(k in 1:K) {
    X_range[[k]] = c(-Inf, Inf)
  }
  return(X_range)
}

#First element is most insignificant (fastest changing), rather than lexicographic
#cell_i and return value are 1-indexed
segment_indexes_from_cell_i <- function(cell_i, n_segments) {
  K = length(n_segments)
  size = cumprod(n_segments)
  if(cell_i > size[K])
    print("Error: too big")
  index = rep(0, K)
  cell_i_rem = cell_i-1 #convert to 0-indexing
  for(k in 1:K) {
    index[k] = cell_i_rem %% n_segments[k]
    cell_i_rem = cell_i_rem %/% n_segments[k]
  }
  index = index+1 #convert from 0-indexing
  return(index)
}

partition_from_split_seq <- function(split_seq, X_range, varnames=NULL, max_include=Inf) {
  part = grid_partition(X_range, varnames)
  for(i in seq_len(min(length(split_seq), max_include))) part = add_partition_split(part, split_seq[[i]])
  return(part)
}

get_factors_from_partition <- function(partition, X, X_range=NULL) {
  X_range = if(is.null(X_range)) partition$X_range else X_range
  factors_by_dim = list()
  if(is_sep_sample(X)) {
    K = ncol(X[[1]])
    for(m in 1:length(X)) {
      factors_by_dim_m = list()
      for(k in 1:K) {
        factors_by_dim_m[[k]] = get_factors_from_splits_dim(X[[m]][, k], X_range[[k]], partition$s_by_dim[[k]])
      }
      factors_by_dim[[m]] = factors_by_dim_m
    }
  }
  else {
    K = ncol(X)
    for(k in 1:K) {
      factors_by_dim[[k]] = get_factors_from_splits_dim(X[, k], X_range[[k]], partition$s_by_dim[[k]])
    }
  }
  return(factors_by_dim)
}

# partition_split ---------------------

#' Create partition_split
#' 
#' Describes a single partition split. Used with \code{\link{add_partition_split}}.
#'
#' @param k dimension
#' @param X_k_cut cut value
#'
#' @return Partition Split
#' @export
partition_split <- function(k, X_k_cut) {
  return(structure(list(k=k, X_k_cut=X_k_cut), class=c("partition_split")))
} 

#' Is grid_partition_split
#' 
#' Tests whether or not an object is a \code{partition_split}.
#'
#' @param x an R object
#'
#' @return Boolean
#' @export
#' @describeIn grid_partition_split is grid_partition_split
is_grid_partition_split <- function(x){ 
  inherits(x, "partition_split") 
}

#' Print partition_split
#' 
#' Prints information for a \code{partition_split}
#'
#' @param x Object
#' @param ... Additional arguments. Unused.
#'
#' @return None
#' @export
print.partition_split <- function(x, ...) {
  cat(paste0(x[[1]], ": ", x[[2]], "\n"))
}

# Search algo --------------------


#' Fit grid_partition
#' 
#' Fit partition on some data, optionally finding best lambda using CV and then re-fiting on full data.
#' 
#' Returns the partition and information about the fitting process
#'  
#' @section Multiple estimates:
#' With multiple core estimates (M) there are 3 options (the first two have the same sample across treatment effects).\enumerate{
#'  \item DS.MULTI_SAMPLE: Multiple pairs of (Y_{m},W_{m}). y,X,d are then lists of length M. Each element then has the typical size
#'     The N_m may differ across m. The number of columns of X will be the same across m.
#'  \item DS.MULTI_D: Multiple treatments and a single outcome. d is then a NxM matrix.
#'  \item DS.MULTI_Y: A single treatment and multiple outcomes. y is then a NXM matrix.
#' }
#'
#' @param y Nx1 matrix of outcome (label/target) data. With multiple core estimates see Details below.
#' @param X NxK matrix of features (covariates). With multiple core estimates see Details below.
#' @param d (Optional) NxP matrix (with colnames) of treatment data. If all equally important they 
#'          should be normalized to have the same variance. With multiple core estimates see Details below.
#' 
#' @param X_aux aux X sample to compute statistics on (OOS data)
#' @param d_aux aux d sample to compute statistics on (OOS data)
#' @param max_splits Maximum number of splits even if splits continue to improve OOS fit
#' @param max_cells Maximum number of cells even if more splits continue to improve OOS fit
#' @param min_size Minimum cell size when building full grid, cv_tr will use (F-1)/F*min_size, cv_te doesn't use any.
#' @param cv_folds Number of CV Folds or a vector of foldids. 
#'                If m_mode==DS.MULTI_SAMPLE, then a list with foldids per Dataset.
#' @param verbosity 0 print no message. 
#'                  1 prints progress bar for high-level loops. 
#'                  2 prints detailed output for high-level loops. 
#'                  Nested operations decrease verbosity by 1.
#' @param breaks_per_dim NULL (for all possible breaks); 
#'                       K-length vector with # of break (chosen by quantiles); or 
#'                       K-dim list of vectors giving potential split points for non-categorical 
#'                         variables (can put c(0) for categorical). 
#'                      Similar to 'discrete splitting' in CausalTree though their they do separate split-points 
#'                      for treated and controls.
#' @param potential_lambdas potential lambdas to search through in CV
#' @param X_range list of min/max for each dimension (e.g., from \code{\link{get_X_range}})
#' @param bucket_min_n Minimum number of observations needed between different split checks
#' @param bucket_min_d_var Ensure positive variance of d for the observations between different split checks
#' @param obj_fn Default is \code{\link{eval_mse_hat}}. User-provided must allow same signature.
#' @param est_plan \link{EstimatorPlan}.
#' @param partition_i Default NA. Use this to avoid CV
#' @param pr_cl Default NULL. Parallel cluster. Used for:\enumerate{
#'                \item CVing the optimal lambda, 
#'                \item fitting full tree (at each split going across dimensions), 
#'                \item fitting trees over the bumped samples
#'              }
#' @param bump_samples Number of bump bootstraps (default 0), or list of such length where each items is a bootstrap sample.
#'                     If m_mode==DS.MULTI_SAMPLE then each item is a sublist with such bootstrap samples over each dataset.
#' @param bump_ratio For bootstraps the ratio of sample size to sample (between 0 and 1, default 1)
#'
#' @return An object.
#'         \item{partition}{Grid Partition (type=\code{\link{grid_partition}})}
#'         \item{is_obj_val_seq}{Full sequence of in-sample objective function values}
#'         \item{complexity_seq}{Full sequence of partition complexities (num_cells - 1)}
#'         \item{partition_i}{Index of partition chosen}
#'         \item{partition_seq}{Full sequence of Grid Partitions}
#'         \item{split_seq}{Full sequence of splits (type=\code{\link{partition_split}})}
#'         \item{lambda}{lambda chosen}
#'         \item{folds_index_out}{List of the held-out observations for each fold (e.g., we might have generated them)}
#' @export
fit_partition <- function(y, X, d=NULL, X_aux=NULL, d_aux=NULL, max_splits=Inf, max_cells=Inf, 
                          min_size=3, cv_folds=2, verbosity=0, breaks_per_dim=NULL, potential_lambdas=NULL, 
                          X_range=NULL, bucket_min_n=NA, bucket_min_d_var=FALSE, obj_fn, 
                          est_plan, partition_i=NA, pr_cl=NULL, bump_samples=0, bump_ratio=1, ...) {
  #Hidden params:
  # - @param lambda_1se Use the 1se rule to pick the best lambda
  # - @param valid_fn Function to quickly check if partition could be valid. User can override.
  # - @param split_check_fn Alternative split-check function
  # - @param N_est N of samples in the Estimation dataset
  # - @param nsplits_k_warn_limit
  # - @param bump_complexity, method 1 is c(FALSE, FALSE), method 2 is c(FALSE, TRUE), and method 3 is c(TRUE)
  extra_params = list(...)
  valid_fn = split_check_fn = NULL
  lambda_1se=FALSE
  N_est=NA
  nsplits_k_warn_limit=200
  bump_complexity=list(doCV=FALSE, incl_comp_in_pick=FALSE)
  if(length(extra_params)>0) {
    if("valid_fn" %in% names(extra_params)) valid_fn = extra_params[['valid_fn']]
    if("split_check_fn" %in% names(extra_params)) split_check_fn = extra_params[['split_check_fn']]
    if("lambda_1se" %in% names(extra_params)) lambda_1se = extra_params[['lambda_1se']]
    if("N_est" %in% names(extra_params)) N_est = extra_params[['N_est']]
    if("nsplits_k_warn_limit" %in% names(extra_params)) nsplits_k_warn_limit = extra_params[['nsplits_k_warn_limit']]
    if("bump_complexity" %in% names(extra_params)) bump_complexity = extra_params[['bump_complexity']]
    good_args = c("valid_fn", "split_check_fn", "lambda_1se", "N_est","nsplits_k_warn_limit", "bump_complexity")
    bad_names = names(extra_params)[!(names(extra_params) %in% good_args)]
    assert_that(length(bad_names)==0, msg=paste(c(list("Illegal arguments:"), bad_names), collapse = " "))
  }
  
  #To check: y, X, d, N_est, X_aux, d_aux, breaks_per_dim, potential_lambdas, X_range, bucket_min_n
  assert_that(max_splits>0, max_cells>0, min_size>0, msg="max_splits, max_cells, min_size need to be positive")
  assert_that(is.flag(lambda_1se), is.flag(bucket_min_d_var), msg="One of (lambda_1se, bucket_min_d_var) are not flags.") 
  assert_that(inherits(est_plan, "estimator_plan") || (is.list(est_plan) && inherits(est_plan[[1]], "estimator_plan")), msg="estimator_plan argument (or it's first element) doesn't inherit from estimator_plan class") 
  #verbosity can be negative if decrementd from a fit_estimate call
  list[M, m_mode, N, K] = get_sample_type(y, X, d, checks=TRUE)
  if(is_sep_sample(X) && length(cv_folds)>1) {
    assert_that(is.list(cv_folds) && length(cv_folds)==M, msg="When separate samples and length(cv_folds)>1, need is.list(cv_folds) && length(cv_folds)==M.")
  }
  check_M_K(M, m_mode, K, X_aux, d_aux)
  do_cv = is.na(partition_i) && (is.null(potential_lambdas) || length(potential_lambdas)>0)
  do_bump = length(bump_samples)>1 || bump_samples > 0
  if(!do_cv) assert_that(bump_complexity$doCV==FALSE, msg="When not doing CV, can't including bumping in CV.")
  if(do_bump && bump_complexity$doCV) {
    if(length(bump_samples==1)) bump_samples = list(bump_samples, bump_samples)
    cv_bump_samples = bump_samples[[1]]
    bump_samples = bump_samples[[2]]
  }
  else cv_bump_samples=0
  
  if(is.null(X_range)) X_range = get_X_range(X)
  if(!is.list(breaks_per_dim) && length(breaks_per_dim)==1) breaks_per_dim = get_quantile_breaks(X, X_range, g=breaks_per_dim)
  if(is.null(valid_fn)) valid_fn = valid_partition
  
  if(is.null(split_check_fn) && (!is.na(bucket_min_n) | bucket_min_d_var)) {
    split_check_fn = purrr::partial(rolling_split_check, bucket_min_n=bucket_min_n, bucket_min_d_var=bucket_min_d_var)
  }
  else{
    split_check_fn = NULL
  }
  
  if(verbosity>0) cat("Grid: Started.\n")
  
  
  if(verbosity>0) cat("Grid: Fitting grid structure on full set\n")
  fit_ret = fit_partition_full(y, X, d, X_aux, d_aux, X_range=X_range, max_splits=max_splits, 
                               max_cells=max_cells, min_size=min_size,  verbosity=verbosity-1, 
                               breaks_per_dim=breaks_per_dim, N_est, split_check_fn=split_check_fn, 
                               obj_fn=obj_fn, allow_empty_aux=FALSE, 
                               allow_est_errors_aux=FALSE, min_size_aux=1, est_plan=est_plan, 
                               pr_cl=pr_cl, valid_fn=valid_fn, nsplits_k_warn_limit=nsplits_k_warn_limit)
  list[partition_seq, is_obj_val_seq, split_seq] = fit_ret
  complexity_seq = sapply(partition_seq, num_cells) - 1
  
  foldids = NA
  if(!is.na(partition_i)) {
    lambda = NA
    max_splits = partition_i-1
    if(length(partition_seq)< partition_i) {
      cat("Note: Couldn't build grid to desired granularity. Using most granular")
      partition_i = length(partition_seq)
    }
    assert_that(bump_complexity$incl_comp_in_pick==FALSE, msg="When no complexity penalization used, can't include complexity cost in bumping calculation.")
  }
  else {
    if(do_cv) {
      list[nfolds, folds_ret, foldids] = expand_fold_info(y, cv_folds, m_mode)
      list[lambda,lambda_oos, n_cell_table] = cv_pick_lambda(y=y, X=X, d=d, folds_ret=folds_ret, nfolds=nfolds, potential_lambdas=potential_lambdas, N_est=N_est, max_splits=max_splits, max_cells=max_cells, 
                              min_size=min_size, verbosity=verbosity, breaks_per_dim=breaks_per_dim, X_range=X_range, lambda_1se=lambda_1se, 
                              split_check_fn=split_check_fn, obj_fn=obj_fn,
                              est_plan=est_plan, pr_cl=pr_cl, valid_fn=valid_fn, cv_bump_samples=cv_bump_samples, bump_ratio=bump_ratio)
    }
    else {
      lambda = potential_lambdas[1]
    }
    partition_i = which.min(is_obj_val_seq + lambda*complexity_seq)
  }

  
  if(do_bump) {
    if(verbosity>0) cat("Grid > Bumping: Started.\n")
    
    if(bump_complexity$incl_comp_in_pick) { 
      best_val = is_obj_val_seq[partition_i] + lambda*complexity_seq[partition_i]
    }
    else {
      best_val = is_obj_val_seq[partition_i]
    }
    
    b_rets = gen_bumped_partitions(bump_samples, bump_ratio, N, m_mode, verbosity, pr_cl, min_size=min_size*bump_ratio, 
                                   y=y, X_d=X, d=d, X_aux=X_aux, d_aux=d_aux, X_range=X_range, max_splits=max_splits, 
                                   max_cells=max_cells,  
                                   breaks_per_dim=breaks_per_dim, N_est=N_est, split_check_fn=split_check_fn, obj_fn=obj_fn, 
                                   min_size_aux=min_size, est_plan=est_plan, 
                                   valid_fn=valid_fn)
    bump_B = length(b_rets)
    
    best_b = NA
    for(b in seq_len(bump_B)) {
      b_ret = b_rets[[b]]
      if(do_cv || bump_complexity$incl_comp_in_pick) b_complexity_seq = sapply(b_ret$partition_seq, num_cells) - 1
      if(do_cv) {
        partition_i_b = which.min(b_ret$is_obj_val_seq + lambda*b_complexity_seq)
      }
      else {
        partition_i_b = partition_i #default
        if(length(b_ret$partition_seq)<partition_i_b) {
          cat("Note: Couldn't build grid to desired granularity. Using most granular\n")
          partition_i_b = length(b_ret$partition_seq)
        }
      }
      partition_b = b_ret$partition_seq[[partition_i_b]]
      
      obj_ret = obj_fn(y, X, d, N_est=N_est, partition=partition_b, est_plan=est_plan, sample="trtr")
      if(obj_ret[2]>0 | obj_ret[3]>0) next #N_cell_empty, N_cell_error
      if(bump_complexity$incl_comp_in_pick) {
        bump_val = obj_ret[1] + lambda*b_complexity_seq[partition_i_b]
      }
      else {
        bump_val = obj_ret[1]
      }
      
      if(bump_val < best_val){
        best_val = bump_val
        best_b = b
      }
    }
    if(!is.na(best_b)) {
      if(verbosity>0) {
        cat(paste("Grid > Bumping: Finished. Picking bumped partition."))
        cat(paste(" Old (unbumped) is_obj_val_seq=[", paste(is_obj_val_seq, collapse=" "), "]."))
        cat(paste(" Old (unbumped) complexity_seq=[", paste(complexity_seq, collapse=" "), "].\n"))
      }
      list[partition_seq, is_obj_val_seq_best_b, split_seq] = b_rets[[best_b]]
      if(do_cv) {
        b_complexity_seq = sapply(b_rets[[best_b]]$partition_seq, num_cells) - 1
        partition_i = which.min(b_rets[[best_b]]$is_obj_val_seq + lambda*b_complexity_seq)
      } 
      complexity_seq = sapply(partition_seq, num_cells) - 1
      is_obj_val_seq = sapply(partition_seq, function(p){
        obj_fn(y, X, d, N_est=N_est, partition=p, est_plan=est_plan, sample="trtr")[1]
      })
    }
    else { 
      if(verbosity>0) cat(paste("Grid > Bumping: Finished. No bumped partitions better than original.\n"))
    }
  }
  
  if(verbosity>0) {
    #print(partition_seq)
    cat(paste("Grid: Finished. is_obj_val_seq=[", paste(is_obj_val_seq, collapse=" "), "]."))
    if(do_cv) {
      cat(paste(" complexity_seq=[", paste(complexity_seq, collapse=" "), "]."))
      cat(paste(" best partition=", paste(partition_i, collapse=" "), "."))
    }
    cat("\n")
  }
  partition = partition_seq[[partition_i]]
  return(list(partition=partition, is_obj_val_seq=is_obj_val_seq, complexity_seq=complexity_seq, 
              partition_i=partition_i, partition_seq=partition_seq, split_seq=split_seq, lambda=lambda, 
              foldids=foldids))
}


#' Get break-points by looking at quantiles
#' 
#' Provides a set of potential split points for data according to quantiles (if possible)
#'
#' @param X Features
#' @param X_range X-range
#' @param g # of quantiles
#' @param type Quantile type (see ?quantile and https://mathworld.wolfram.com/Quantile.html). 
#'             Types1-3 are discrete and this is good for passing to unique() when there are clumps
#'
#' @return list of potential breaks
get_quantile_breaks <- function(X, X_range, g=20, type=3) {
  if(is.null(g)) g=20 #fit_estimate has a different default that might get passed in.
  if(is_sep_sample(X)) X = X[[1]]
  X = ensure_good_X(X)
  
  breaks_per_dim = list()
  K = ncol(X)
  for(k in 1:K) {
    X_k = X[,k]
    if(is.factor(X_k)) {
      breaks_per_dim[[k]] = c(0) #Dummy
    }
    else {
      if(storage.mode(X_k)=="integer" && (X_range[[k]][2]-X_range[[k]][1])<=g) {
        vals = sort(unique(X_k))
        breaks_per_dim[[k]] = vals[-c(length(vals), 1)]
      }
      else {
        #unique(sort(X[,k])) #we will automatically skip the top point
        #if you want g cuts, then there are g+2 outer nodes
        qs = quantile(X_k, seq(0, 1, length.out=g+2), names=FALSE, type=type)
        qs = unique(qs)
        breaks_per_dim[[k]] = qs[-c(length(qs), 1)]
      }
    }
  }
  return(breaks_per_dim)
}

# if d vectors are empty doesn't return fail 
valid_partition <- function(cell_factor, d=NULL, cell_factor_aux=NULL, d_aux=NULL, min_size=0) {
  #check none of the cells are too small
  if(min_size>0) {
    if(length(cell_factor)==0) return(list(fail=TRUE, min_size=0))
    lowest_size = min(table(cell_factor))
    if(lowest_size<min_size) return(list(fail=TRUE, min_size=lowest_size))
    
    if(!is.null(cell_factor_aux)) {
      if(length(cell_factor_aux)==0) return(list(fail=TRUE, min_size_aux=0))
      lowest_size_aux = min(table(cell_factor_aux))
      if(lowest_size_aux<min_size) return(list(fail=TRUE, min_size_aux=lowest_size_aux))
    }
  }
  
  if(!is.null(d)) {
    if(!is_vec(d)) {
      for(m in 1:ncol(d)) {
        if(any(by(d[,m], cell_factor, FUN=const_vect))) {
          return(list(fail=TRUE, always_d_var=FALSE))
        }
      }
      
    }
    else {
      if(any(by(as.vector(d), cell_factor, FUN=const_vect))) {
        return(list(fail=TRUE, always_d_var=FALSE))
      }
    }
  }
  if(!is.null(d_aux)) {
    if(!is_vec(d_aux)) {
      for(m in 1:ncol(d_aux)) {
        if(any(by(d_aux[,m], cell_factor, FUN=const_vect))) {
          return(list(fail=TRUE, always_d_var=FALSE))
        }
      }
      
    }
    else {
      if(any(by(as.vector(d_aux), cell_factor_aux, FUN=const_vect))) {
        return(list(fail=TRUE, always_d_var_aux=FALSE))
      }
    }
  }
  return(list(fail=FALSE))
}


gen_bumped_partitions <- function(bump_samples, bump_ratio, N, m_mode, verbosity, pr_cl, allow_empty_aux=FALSE, allow_est_errors_aux=FALSE, ...) {
  assert_that(bump_ratio>0, bump_ratio<=1, msg="bump_ration needs to be in (0,1]")
  bump_samples = expand_bump_samples(bump_samples, bump_ratio, N, m_mode)
  bump_B = length(bump_samples)
  
  params = c(list(samples=bump_samples, verbosity=verbosity-1, allow_empty_aux=FALSE, allow_est_errors_aux=FALSE, pr_cl=NULL, m_mode=m_mode),
             list(...))
  
  b_rets = my_apply(1:bump_B, fit_partition_bump_b, verbosity==1 || !is.null(pr_cl), pr_cl, params)
  return(b_rets)
}

#if not mid-point then the all but the last are the splits
get_usable_break_points <- function(breaks_per_dim, X, X_range, dim_cat, mid_point=TRUE) {
  if(is_sep_sample(X)) X = X[[1]]
  K = ncol(X)
  #old code
  if(is.null(breaks_per_dim)) {
    breaks_per_dim = list()
    for(k in 1:K) {
      if(!k %in% dim_cat) {
        u = unique(sort(X[, k]))
        if(mid_point) {
          breaks_per_dim[[k]] = u[-length(u)] + diff(u) / 2
        }
        else {
          breaks_per_dim[[k]] = u[-length(u)] #skip last point
        }
      }
      else {
        breaks_per_dim[[k]] = c(0) #Dummy just for place=holder
      }
    }
  }
  else { #make sure they didn't include the lowest point
    for(k in 1:K) {
      if(!k %in% dim_cat) {
        n_k = length(breaks_per_dim[[k]])
        if(breaks_per_dim[[k]][n_k]==X_range[[k]][2]) {
          breaks_per_dim[[k]] = breaks_per_dim[[k]][-n_k]
        }
      }
      breaks_per_dim[[k]] = unname(breaks_per_dim[[k]]) #names messed up the get_desc_df() (though not in debugSource)
    }
  }
  return(breaks_per_dim)
}


#Typically is_obj_val_seq trends negative. If first element is min, then return c()
get_lambda_ties <- function(is_obj_val_seq, complexity_seq) {
  n_seq = length(is_obj_val_seq)
  slopes = c() #will go from strongly negative and increases and we stop before reaching 0
  hull_i = 1
  while(hull_i < n_seq) {
    i_slopes = rep(NA, n_seq)
    for(i in (hull_i+1):n_seq) {
      i_slopes[i] = (is_obj_val_seq[i] - is_obj_val_seq[hull_i])/(complexity_seq[i]- complexity_seq[hull_i])
    }
    best_slope = min(i_slopes, na.rm=TRUE)
    if(best_slope>=0) break
    slopes = c(slopes, best_slope)
    hull_i = which.min(i_slopes)
  }
  if(length(slopes)>1) {
    lambda_ties = abs(slopes) #slightly bigger will go will pick the index earlier, slightly bigger later
  }
  else {
    lambda_ties = c()
  }
  return(lambda_ties)
}

gen_cat_window_splits <- function(chr_vec) {
  n = length(chr_vec)
  splits=list()
  for(m in seq_len(floor(n/2))) {
    cs = combn(chr_vec, m, simplify=F)
    if(m==n/2) cs = cs[1:(length(cs)/2)] #or just filter by those that contain chr_vec[1]
    splits = c(splits, cs)
  }
  return(splits)
}

n_cat_window_splits <- function(window_len) {
  n_splits = 0
  for(m in seq_len(floor(window_len/2))) {
    n_choose = choose(window_len, m)
    n_splits = n_splits + if(m==window_len/2) n_choose/2 else n_choose
  }
  return(n_splits)
}

n_cat_splits <- function(s_by_dim_k, X_range_k) {
  windows = get_windows_cat(s_by_dim_k, X_range_k)
  n_splits = 0
  for(window in windows) n_splits = n_splits + n_cat_window_splits(length(window))
  return(n_splits)
}

get_windows_cat <- function(s_by_dim_k, X_k_range) {
  windows = s_by_dim_k
  windows[[length(windows)+1]] = X_k_range[!X_k_range %in% unlist(c(windows))]
  return(windows)
}

get_window_cont <- function(s_by_dim_k, X_k_range) {
  windows=list()
  n_w = length(s_by_dim_k)+1
  for(w in 1:n_w) {
    wmin = if(w==1) X_k_range[1] else s_by_dim_k[w-1]
    wmax = if(w==n_w) X_k_range[2] else s_by_dim_k[w]
    windows[[w]] = c(wmin, wmax)
  }
  return(windows)
}

gen_holdout_interaction <- function(factors_by_dim, k) {
  if(length(factors_by_dim)>1)
    return(interaction(factors_by_dim[-k]))
  return(factor(rep("|", length(factors_by_dim[[1]]))))
}

n_breaks_k <- function(breaks_per_dim, k, partition, X_range) {
  if(k %in% partition$dim_cat) return(n_cat_splits(partition$s_by_dim[[k]], X_range[[k]]))
  return(length(breaks_per_dim[[k]]))
}


rolling_split_check <- function(shifted_N, shifted_d=NULL, shifted_cell_factor_nk, m_mode, bucket_min_n=NA, bucket_min_d_var=FALSE) {
  if(!is.na(bucket_min_n) && min(shifted_N)<bucket_min_n){
    #cat("Skipped: increment not big enough\n")
    return(FALSE)
  } 
  
  if(bucket_min_d_var && !is.null(shifted_d) && any_const_m(shifted_d, shifted_cell_factor_nk, m_mode)) {
    return(FALSE)
  }
  
  return(TRUE)
}

fit_partition_full_k <- function(k, y, X_d, d, X_range, pb, debug, valid_breaks, factors_by_dim, X_aux, 
                                 factors_by_dim_aux, partition, verbosity, allow_empty_aux=TRUE, d_aux, 
                                 allow_est_errors_aux, min_size, min_size_aux, obj_fn, N_est, est_plan, 
                                 split_check_fn = NULL, breaks_per_dim, valid_fn=NULL) { #, n_cut
  assert_that(is.flag(allow_empty_aux), msg="allow_empty_aux needs to be logical flags.")
  list[M, m_mode, N, K] = get_sample_type(y, X_d, d, checks=FALSE)
  if(is.null(valid_fn)) valid_fn = valid_partition
  search_ret = list()
  best_new_val = Inf
  valid_breaks_k = valid_breaks[[k]]
  cell_factor_nk = gen_holdout_interaction_m(factors_by_dim, k, is_sep_sample(X_d))
  if(!is.null(X_aux)) {
    cell_factor_nk_aux = gen_holdout_interaction_m(factors_by_dim_aux, k, is_sep_sample(X_aux))
  }
  
  if(!is_factor_dim_k_m(X_d, k, m_mode==DS.MULTI_SAMPLE)) {
    n_pot_break_points_k = length(breaks_per_dim[[k]])
    vals = rep(NA, n_pot_break_points_k)
    prev_split_checked = X_range[[k]][1]
    win_LB = X_range[[k]][1]-1
    win_UB = if(length(partition$s_by_dim[[k]])>0) partition$s_by_dim[[k]][1] else X_range[[k]][2]
    win_mask = gen_cont_window_mask_m(X_d, k, win_LB, win_UB)
    win_mask_aux = gen_cont_window_mask_m(X_aux, k, win_LB, win_UB)
    for(X_k_cut_i in seq_len(n_pot_break_points_k)) { #cut-point is top end of segment, 
      if (verbosity>0 && !is.null(pb)) utils::setTxtProgressBar(pb, utils::getTxtProgressBar(pb)+1)
      X_k_cut = breaks_per_dim[[k]][X_k_cut_i]
      if(X_k_cut %in% partition$s_by_dim[[k]]) {
        prev_split_checked = X_k_cut
        win_LB = X_k_cut
        higher_prev_split = partition$s_by_dim[[k]][partition$s_by_dim[[k]]>X_k_cut]
        win_UB = if(length(higher_prev_split)>0) min(higher_prev_split) else X_range[[k]][2]
        win_mask = gen_cont_window_mask_m(X_d, k, win_LB, win_UB)
        win_mask_aux = gen_cont_window_mask_m(X_aux, k, win_LB, win_UB)
        next
      } 
      if(!valid_breaks_k[[1]][X_k_cut_i]) next
      new_split = partition_split(k, X_k_cut)
      tent_partition = add_partition_split(partition, new_split)
      
      tent_split_fac_k = get_factors_from_splits_dim_m(X_d, X_range[[k]], tent_partition$s_by_dim[[k]], k)
      tent_cell_factor = interaction2_m(cell_factor_nk, tent_split_fac_k, m_mode==DS.MULTI_SAMPLE)
      if(!is.null(X_aux)) {
        tent_split_fac_k_aux = get_factors_from_splits_dim_m(X_aux, X_range[[k]], tent_partition$s_by_dim[[k]], k)
        tent_cell_factor_aux = interaction2_m(cell_factor_nk_aux, tent_split_fac_k_aux, is_sep_sample(X_aux))
      }
      
      if(!is.null(split_check_fn)){
        shifted_mask = gen_cont_window_mask_m(X_d, k, prev_split_checked, X_k_cut)
        shifted_N = sum_m(shifted_mask, m_mode==DS.MULTI_SAMPLE)
        shifted_cell_factor_nk = droplevels_m(apply_mask_m(cell_factor_nk, shifted_mask, m_mode==DS.MULTI_SAMPLE), m_mode==DS.MULTI_SAMPLE)
        shifted_d = if(is.null(d)) NULL else apply_mask_m(d, shifted_mask, m_mode==DS.MULTI_SAMPLE)
        split_OK = split_check_fn(shifted_N, shifted_d, shifted_cell_factor_nk, m_mode)
        if(!split_OK) {
          valid_breaks_k[[1]][X_k_cut_i] = FALSE
          next
        }
      }
      
      # do_window_approach
      #The bucket checks don't help much. 
      #- Though I do check for non-zero var of D, that's just on the left so to check on right side too
      #- Note that though not min_size as different than bucket_min_n)
      win_split_cond = gen_cont_win_split_cond_m(X_d, win_mask, k, X_k_cut)
      win_cell_factor_nk = apply_mask_m(cell_factor_nk, win_mask, m_mode==DS.MULTI_SAMPLE)
      win_cell_factor = interaction2_m(win_cell_factor_nk, win_split_cond, m_mode==DS.MULTI_SAMPLE)
      win_d = apply_mask_m(d, win_mask, m_mode==DS.MULTI_SAMPLE)
      valid_ret = valid_partition_m(m_mode==DS.MULTI_SAMPLE, valid_fn, win_cell_factor, d=win_d, min_size=min_size)
      if(!valid_ret$fail) {
        if(!allow_empty_aux && !is.null(X_aux)) {
          win_split_cond_aux = gen_cont_win_split_cond_m(X_aux, win_mask_aux, k, X_k_cut)
          win_cell_factor_aux = interaction2_m(apply_mask_m(cell_factor_nk_aux, win_mask_aux, is_sep_sample(X_aux)), 
                                               win_split_cond_aux, is_sep_sample(X_aux), drop=allow_empty_aux)
          win_d_aux = if(!allow_est_errors_aux) apply_mask_m(d_aux, win_mask_aux, is_sep_sample(X_aux)) else NULL
          valid_ret = valid_partition_m(is_sep_sample(X_aux), valid_fn, win_cell_factor_aux, d=win_d_aux, min_size=min_size_aux)
        }
      }
      # Global approach
      # valid_ret = valid_fn(tent_cell_factor, d=d, min_size=min_size)
      # if(!valid_ret$fail) {
      #   valid_ret = valid_fn(tent_cell_factor_aux, d=d_aux, min_size=2)
      # }
      if(valid_ret$fail) {
        #cat("Invalid partition\n")
        valid_breaks_k[[1]][X_k_cut_i] = FALSE
        next
      }
      if(debug) cat(paste("k", k, ". X_k", X_k_cut, "\n"))
      obj_ret = obj_fn(y, X_d, d, N_est=N_est, cell_factor_tr = tent_cell_factor, debug=debug, est_plan=est_plan, 
                       sample="trtr")
      if(obj_ret[3]>0) { #don't need to check [2] (empty cells) as we already did that
        #cat("Estimation errors\n")
        valid_breaks_k[[1]][X_k_cut_i] = FALSE
        next
      }
      val = obj_ret[1]
      stopifnot(is.finite(val))
      prev_split_checked = X_k_cut
      
      if(val<best_new_val) {
        #if(verbosity>0) print(paste("Testing split at ", X_k_cut, ". Val=", split_res$val))
        best_new_val = val
        new_factors_by_dim = replace_k_factor_m(factors_by_dim, k, tent_split_fac_k, is_sep_sample(X_d))
        if(!is.null(X_aux)) {
          new_factors_by_dim_aux = replace_k_factor_m(factors_by_dim_aux, k, tent_split_fac_k_aux, is_sep_sample(X_aux))
        }
        else new_factors_by_dim_aux = NULL
        search_ret = list(val=val, new_split=new_split, new_factors_by_dim=new_factors_by_dim, 
                          new_factors_by_dim_aux=new_factors_by_dim_aux)
      }
    }
    
  }
  else { #categorical variable
    windows = get_windows_cat(partition$s_by_dim[[k]], X_range[[k]])
    for(window_i in seq_len(length(windows))) {
      window = windows[[window_i]]
      win_mask = gen_cat_window_mask_m(X_d, k, window)
      win_mask_aux = gen_cat_window_mask_m(X_aux, k, window)
      pot_splits = gen_cat_window_splits(window)
      for(win_split_i in seq_len(length(pot_splits))) {
        win_split_val = pot_splits[[win_split_i]]
        #TODO: Refactor with continuous case
        if (verbosity>0 && !is.null(pb)) utils::setTxtProgressBar(pb, utils::getTxtProgressBar(pb)+1)
        if(!valid_breaks_k[[window_i]][win_split_i]) next
        
        new_split = partition_split(k, win_split_val)
        tent_partition = add_partition_split(partition, new_split)
        
        tent_split_fac_k = get_factors_from_splits_dim_m(X_d, X_range[[k]], tent_partition$s_by_dim[[k]], k)
        tent_cell_factor = interaction2_m(cell_factor_nk, tent_split_fac_k, m_mode==DS.MULTI_SAMPLE)
        if(!is.null(X_aux)) {
          tent_split_fac_k_aux = get_factors_from_splits_dim_m(X_aux, X_range[[k]], tent_partition$s_by_dim[[k]], k)
          tent_cell_factor_aux = interaction2_m(cell_factor_nk_aux, tent_split_fac_k_aux, is_sep_sample(X_aux))
        }
        
        # do_window_approach
        win_split_cond = gen_cat_win_split_cond_m(X_d, win_mask, k, win_split_val)
        win_cell_factor = interaction2_m(apply_mask_m(cell_factor_nk, win_mask, m_mode==DS.MULTI_SAMPLE), win_split_cond, m_mode==DS.MULTI_SAMPLE)
        win_d = apply_mask_m(d, win_mask, m_mode==DS.MULTI_SAMPLE)
        valid_ret = valid_partition_m(m_mode==DS.MULTI_SAMPLE, valid_fn, win_cell_factor, d=win_d, min_size=min_size)
        if(!valid_ret$fail) {
          if(!is.null(X_aux) && !allow_empty_aux) {
            win_split_cond_aux = factor(gen_cat_win_split_cond_m(X_aux, win_mask_aux, k, win_split_val), levels=c(FALSE, TRUE))
            win_cell_factor_aux = interaction2_m(apply_mask_m(cell_factor_nk_aux, win_mask_aux, is_sep_sample(X_aux)), 
                                                 win_split_cond_aux, is_sep_sample(X_aux), drop=allow_empty_aux)
            win_d_aux = if(!allow_est_errors_aux) apply_mask_m(d_aux, win_mask_aux, is_sep_sample(X_aux)) else NULL
            valid_ret = valid_partition_m(is_sep_sample(X_aux), valid_fn, win_cell_factor_aux, d=win_d_aux, min_size=min_size_aux)
          }
        }
        if(valid_ret$fail) {
          #cat("Invalid partition\n")
          valid_breaks_k[[window_i]][win_split_i] = FALSE
          next
        }
        if(debug) cat(paste("k", k, ". X_k", win_split_val, "\n"))
        obj_ret = obj_fn(y, X_d, d, N_est=N_est, cell_factor_tr = tent_cell_factor, debug=debug, est_plan=est_plan, 
                         sample="trtr")
        if(obj_ret[3]>0) { #don't need to check [2] (empty cells) as we already did that
          #cat("Estimation errors\n")
          valid_breaks_k[[window_i]][win_split_i] = FALSE
          next
        }
        val = obj_ret[1]
        stopifnot(is.finite(val))
        
        
        if(val<best_new_val) {
          #if(verbosity>0) print(paste("Testing split at ", X_k_cut, ". Val=", split_res$val))
          best_new_val = val
          new_factors_by_dim = replace_k_factor_m(factors_by_dim, k, tent_split_fac_k, is_sep_sample(X_d))
          if(!is.null(X_aux)) {
            new_factors_by_dim_aux = replace_k_factor_m(factors_by_dim_aux, k, tent_split_fac_k_aux, is_sep_sample(X_aux))
          }
          else new_factors_by_dim_aux = NULL
          search_ret = list(val=val, new_split=new_split, new_factors_by_dim=new_factors_by_dim, 
                            new_factors_by_dim_aux=new_factors_by_dim_aux)
        }
        
      }
    }
  }

  return(list(search_ret, valid_breaks_k))
}

# There are three general problems with a partition. 
# 1) Empty cells
# 2) Non-empty cells where objective can't be calculated
# 3) Cells where it can be calulcated but due to small sizes we don't want
# Main sample: We assume that a valid partition removes #1 and #2. Use min_size for 3
# For Aux: Use allow_empty_aux, allow_est_errors_aux, min_size_aux
# Include d_aux if you want to make sure that non-empty cells in aux have positive variance in d
# FOr CV:allow_empty_aux=TRUE, allow_est_errors_aux=FALSE, min_size_aux=1 (weaker check than removing estimation errors)
# Can set nsplits_k_warn_limit=Inf to disable
fit_partition_full <- function(y, X, d=NULL, X_aux=NULL, d_aux=NULL, X_range, max_splits=Inf, max_cells=Inf, 
                               min_size=2, verbosity=0, breaks_per_dim, N_est, obj_fn, allow_est_errors_aux=TRUE, 
                               min_size_aux=2, est_plan, partition=NULL, nsplits_k_warn_limit=200, pr_cl=NULL, 
                               ...) {
  assert_that(max_splits>=0, max_cells>=1, min_size>=1, msg="Need max_splits>=0, max_cells>=1, min_size>=1.")
  assert_that(is.flag(allow_est_errors_aux), msg="allow_est_errors_aux needs to be a flag")
  assert_that(is.na(nsplits_k_warn_limit) || nsplits_k_warn_limit>=1, msg="nsplits_k_warn_limit not understood")
  list[M, m_mode, N, K] = get_sample_type(y, X, d, checks=TRUE)
  est_min = ifelse(is.null(d), 2, 3) #If don't always need variance calc: ifelse(is.null(d), ifelse(honest, 2, 1), ifelse(honest, 3, 2))
  min_size = max(min_size, est_min)
  if(!allow_est_errors_aux)  min_size_aux = max(min_size_aux, est_min)
  debug = FALSE
  if(is.null(partition)) partition = grid_partition(X_range, colnames(X))
  breaks_per_dim = get_usable_break_points(breaks_per_dim, X, X_range, partition$dim_cat)
  valid_breaks = vector("list", length=K) #splits_by_dim(s_seq) #stores Xk_val's
  
  for(k in 1:K) {
    n_split_breaks_k = n_breaks_k(breaks_per_dim, k, partition, X_range)
    valid_breaks[[k]] = list(rep(TRUE, n_split_breaks_k))
    if(!is.na(nsplits_k_warn_limit) && n_split_breaks_k>nsplits_k_warn_limit) warning(paste("Warning: Many splits (", n_split_breaks_k, ") along dimension", k, "\n"))
  }
  factors_by_dim = get_factors_from_partition(partition, X)
  if(!is.null(X_aux)) {
    factors_by_dim_aux = get_factors_from_partition(partition, X_aux)
  }
  
  if(verbosity>0){
    cat("Grid > Fitting: Started.\n")
    t0 = Sys.time()
  } 
  split_i = 1
  seq_val = c()
  obj_ret = obj_fn(y, X, d, N_est=N_est, cell_factor_tr = interaction_m(factors_by_dim, is_sep_sample(X)), est_plan=est_plan, sample="trtr")
  if(obj_ret[3]>0 || !is.finite(obj_ret[1])) {
    stop("Estimation error with initial partition")
  }
  seq_val[1] = obj_ret[1]
  partition_seq = list()
  split_seq = list()
  partition_seq[[1]] = partition
  tent_cell_factor_aux = NULL
  style = if(summary(stdout())$class=="terminal") 3 else 1
  if(!is.null(pr_cl) & !requireNamespace("parallel", quietly = TRUE)) {
    stop("Package \"parallel\" needed for this function to work. Please install it.", call. = FALSE)
  }
  do_pbapply = requireNamespace("pbapply", quietly = TRUE) & (verbosity>0) & (is.null(pr_cl) || length(pr_cl)<K)
  
  while(TRUE) {
    if(split_i>max_splits) break
    if(num_cells(partition)==max_cells) break
    n_cuts_k = rep(0, K)
    for(k in 1:K) {
      n_cuts_k[k] = n_breaks_k(breaks_per_dim, k, partition, X_range)
    }
    n_cuts_total = sum(n_cuts_k)
    if(n_cuts_total==0) break
    if(verbosity>0) {
      cat(paste("Grid > Fitting > split ", split_i, ": Started\n"))
      t1 = Sys.time()
      if(is.null(pr_cl)) pb = utils::txtProgressBar(0, n_cuts_total, style = style)
    }
    
    params = c(list(y=y, X_d=X, d=d, X_range=X_range, pb=NULL, debug=debug, valid_breaks=valid_breaks, 
                  factors_by_dim=factors_by_dim, X_aux=X_aux, factors_by_dim_aux=factors_by_dim_aux, partition=partition, 
                  verbosity=verbosity, d_aux=d_aux, allow_est_errors_aux=allow_est_errors_aux, 
                  min_size=min_size, min_size_aux=min_size_aux, obj_fn=obj_fn, N_est=N_est, est_plan=est_plan, 
                  breaks_per_dim=breaks_per_dim), list(...))
    
    col_rets = my_apply(1:K, fit_partition_full_k, verbosity, pr_cl, params)
    
    best_new_val = Inf
    best_new_split = NULL
    for(k in 1:K) {
      col_ret = col_rets[[k]]
      search_ret = col_ret[[1]]
      valid_breaks[[k]] = col_ret[[2]]
      if(length(search_ret)>0 && search_ret$val<best_new_val) {
        #if(verbosity>0) print(paste("Testing split at ", X_k_cut, ". Val=", split_res$val))
        best_new_val = search_ret$val
        best_new_split = search_ret$new_split
        best_new_factors_by_dim = search_ret$new_factors_by_dim
        if(!is.null(X_aux)) {
          best_new_factors_by_dim_aux = search_ret$new_factors_by_dim_aux
        }
      }
    }
    
    if (verbosity>0) {
      t2 = Sys.time() #can us as.numeric(t1) to convert to seconds
      td = t2-t1
      if(is.null(pr_cl)) close(pb)
    }
    if(is.null(best_new_split)) {
      if (verbosity>0) cat(paste("Grid > Fitting > split ", split_i, ": Finished. Duration: ", format(as.numeric(td)), " ", attr(td, "units"), ". No valid splits\n"))
      break
    } 
    best_new_partition = add_partition_split(partition, best_new_split)
    if(num_cells(best_new_partition)>max_cells) {
      if (verbosity>0) cat(paste("Grid > Fitting > split ", split_i, ": Finished. Duration: ", format(as.numeric(td)), " ", attr(td, "units"), ". Best split has results in too many cells\n"))
      break
    } 
    partition = best_new_partition
    factors_by_dim = best_new_factors_by_dim
    if(!is.null(X_aux)) {
      factors_by_dim_aux = best_new_factors_by_dim_aux
    }
    split_i = split_i + 1
    seq_val[split_i] = best_new_val
    partition_seq[[split_i]] = partition
    split_seq[[split_i-1]] = best_new_split
    if(best_new_split[[1]] %in% partition$dim_cat) {
      k = best_new_split[[1]]
      windows = get_windows_cat(partition$s_by_dim[[k]], X_range[[k]])
      nwindows = length(windows)
      v_breaks = vector("list", length=nwindows)
      for(window_i in seq_len(nwindows)) {
        v_breaks[[window_i]] = rep(TRUE, n_cat_window_splits(length(windows[[window_i]])))
      }
      valid_breaks[[k]] = v_breaks
    }
    if (verbosity>0) { 
      cat(paste("Grid > Fitting > split ", split_i, ": Finished.",
                " Duration: ", format(as.numeric(td)), " ", attr(td, "units"), ".",
                " New split: k=", best_new_split[[1]], ", cut=", best_new_split[[2]], ", val=", best_new_val, "\n"))
    }
  }
  if (verbosity>0) {
    tn = Sys.time()
    td = tn-t0
    cat("Grid > Fitting: Finished.")
    cat(paste(" Entire Search Duration: ", format(as.numeric(td)), " ", attr(td, "units"), "\n"))
  }
  return(list(partition_seq=partition_seq, is_obj_val_seq=seq_val, split_seq=split_seq))
}

#Allows two lists or two datasets 
add_samples <- function (X, X_aux, M_mult) {
  if(is.null(X_aux)) return(X)
  if(M_mult) return(c(X, X_aux))
  return(list(X, X_aux))
}


fit_partition_bump_b <- function(b, samples, y, X_d, d=NULL, m_mode, X_aux, d_aux, verbosity, nsplits_k_warn_limit=NA, ...){
  if(verbosity>0) cat(paste("Grid > Bumping > b = ", b, "\n"))
  sample = samples[[b]]
  list[y_b, X_b, d_b] = subsample_m(y, X_d, d, sample)
  X_aux2 = add_samples(X_d, X_aux, is_sep_sample(X_d))
  d_aux2 = add_samples(d, d_aux, is_sep_sample(X_d))
  fit_partition_full(y=y_b, X=X_b, d=d_b, X_aux=X_aux2, d_aux=d_aux2, verbosity=verbosity, nsplits_k_warn_limit=NA, ...)
}

# These are bump wrappers
get_part_for_lambda <- function(obj, lambda, is_bumped=FALSE) {
  if(is_bumped) {
    is_obj_val_seq = unlist(lapply(obj, function(f) f$is_obj_val_seq))
    complexity_seq = unlist(lapply(obj, function(f) sapply(f$partition_seq, num_cells) - 1))
    partition_seq = unlist(lapply(obj, function(f) f$partition_seq ), recursive = FALSE)
  }
  else {
    is_obj_val_seq = obj$is_obj_val_seq
    complexity_seq = sapply(obj$partition_seq, num_cells) - 1
    partition_seq = obj$partition_seq
  }
  partition_i = which.min(is_obj_val_seq + lambda*complexity_seq)
  return(list(partition_i, partition_seq[[partition_i]]))
}
get_num_parts <- function(cvtr_fit, is_bumped=FALSE) {
  if(is_bumped)
    return(sum(sapply(cvtr_fit, function(part) length(part$is_obj_val_seq))))
  return(length(cvtr_fit$is_obj_val_seq))
}

get_all_lambda_ties <- function(cvtr_fit, is_bumped=FALSE) {
  if(is_bumped) {
    return(unlist(lapply(cvtr_fit, function(f) get_lambda_ties(f$is_obj_val_seq, sapply(f$partition_seq, num_cells) - 1))))
  }
  return(get_lambda_ties(cvtr_fit$is_obj_val_seq, sapply(cvtr_fit$partition_seq, num_cells) - 1))
}

# ... params sent to fit_partition_full()
cv_pick_lambda_f <- function(f, y, X_d, d, folds_ret, nfolds, potential_lambdas, N_est, 
                             verbosity, obj_fn, cv_tr_min_size, est_plan, cv_bump_samples, bump_ratio, 
                             nsplits_k_warn_limit=NA, min_size_aux=1, allow_empty_aux=TRUE, allow_est_errors_aux=FALSE, recal_is_obj_b=TRUE, ...) { #catch some of the params that might still be in ...
  if(verbosity>0) cat(paste("Grid > CV > Fold", f, "\n"))
  supplied_lambda = !is.null(potential_lambdas)
  if(supplied_lambda) n_lambda = length(potential_lambdas)
  
  list[y_f_tr, y_f_cv, X_f_tr, X_f_cv, d_f_tr, d_f_cv] = split_sample_folds_m(y, X_d, d, folds_ret, f)

  do_bump = (length(cv_bump_samples)>1 || cv_bump_samples>0)
  cvtr_fit = fit_partition_full(y_f_tr, X_f_tr, d_f_tr, X_f_cv, d_f_cv, 
                                min_size=cv_tr_min_size, verbosity=verbosity, 
                                N_est=N_est, obj_fn=obj_fn, allow_empty_aux=TRUE, 
                                allow_est_errors_aux=FALSE, min_size_aux=1, est_plan=est_plan, 
                                nsplits_k_warn_limit=NA, ...) #min_size_aux is weaker than removing est errors
  
  if(do_bump) {
    if(length(cv_bump_samples)>1) cv_bump_samples = cv_bump_samples[[f]]
    list[M, m_mode, N_tr, K] = get_sample_type(y_f_tr, X_f_tr, d_f_tr, checks=FALSE)
    cvtr_fit_bumps = gen_bumped_partitions(bump_samples=cv_bump_samples, bump_ratio, N_tr, m_mode, verbosity=verbosity, pr_cl=NULL, 
                                     min_size=cv_tr_min_size*bump_ratio, 
                                     y=y_f_tr, X_d=X_f_tr, d=d_f_tr, X_aux=X_f_cv, d_aux=d_f_cv, 
                                     N_est=N_est, obj_fn=obj_fn, 
                                     min_size_aux=1, est_plan=est_plan, nsplits_k_warn_limit=NA,
                                     allow_empty_aux=allow_empty_aux, allow_est_errors_aux=allow_est_errors_aux,
                                     ...)
    if(recal_is_obj_b) {
      #Use the updated values on the unbumped sample
      for(b in 1:length(cvtr_fit_bumps)) {
        #partition_seq=partition_seq, is_obj_val_seq
        cvtr_fit_bumps[[b]]$is_obj_val_seq = sapply(cvtr_fit_bumps[[b]]$partition_seq, function(p){
          obj_fn(y_f_tr, X_f_tr, d_f_tr, partition=p, est_plan=est_plan, sample="trtr")[1]
        })
      }
    }
    cvtr_fit = c(list(cvtr_fit), cvtr_fit_bumps)
  }
  
  if(!supplied_lambda) {
    return(cvtr_fit)
  } 
  #If we know the lambdas, eval data while we have it
  return(eval_lambdas(obj_fn, est_plan, potential_lambdas, cvtr_fit, y_f_tr, X_f_tr, d_f_tr, y_f_cv, X_f_cv, d_f_cv, N_est, do_bump))
}


eval_lambdas <- function(obj_fn, est_plan, potential_lambdas, cvtr_fit, y_f_tr, X_f_tr, d_f_tr, y_f_cv, X_f_cv, d_f_cv, N_est, is_bumped) {
  partition_oos_cache = rep(NA, get_num_parts(cvtr_fit, is_bumped))
  n_lambda = length(potential_lambdas)
  
  lambda_oos = rep(NA, n_lambda)
  for(lambda_i in seq_len(n_lambda)) {
    lambda = potential_lambdas[lambda_i]
    list[partition_i, part] = get_part_for_lambda(cvtr_fit, lambda, is_bumped)
    if(is.na(partition_oos_cache[partition_i])) {
      debug = FALSE 
      if(debug) cat(paste("s_by_dim", paste(part$s_by_dim, collapse=" "), "\n"))
      obj_ret = obj_fn(y_f_tr, X_f_tr, d_f_tr, y_f_cv, X_f_cv, d_f_cv, N_est=N_est, partition=part, debug=debug, 
                       est_plan=est_plan, sample="trcv")
      oos_obj_val = obj_ret[1]
      partition_oos_cache[partition_i] = oos_obj_val
    }
    lambda_oos[lambda_i] = partition_oos_cache[partition_i]
  }
  return(lambda_oos)
}

lambda_1se_selector <- function(potential_lambdas, lambda_oos, min_obs_1se, max_oos_err_allowed, verbosity) {
  n_lambda = ncol(lambda_oos)
  nfolds = nrow(lambda_oos)
  lambda_oos_means = colMeans(lambda_oos)
  lambda_oos_min_i = which.min(lambda_oos_means)
  #lambda_oos_sd = apply(X=lambda_oos, MARGIN=2, FUN=sd) #don't need all for now
  obs = lambda_oos[, lambda_oos_min_i]
  if(min_obs_1se>nfolds) {
    for(delta in 1:min(n_lambda-lambda_oos_min_i, lambda_oos_min_i-1)) {
      if(lambda_oos_min_i+delta<=n_lambda) {
        obs = c(obs, lambda_oos[, lambda_oos_min_i+delta])
      }
      if(lambda_oos_min_i-delta>=1) {
        obs = c(obs, lambda_oos[, lambda_oos_min_i-delta])
      }
      if(length(obs)>=min_obs_1se) break
    }
    print(obs)
  }
  max_oos_err_allowed = min(lambda_oos_means) + sd(obs)
  if(verbosity>0) cat(paste("max_oos_err_allowed:", paste(max_oos_err_allowed, collapse=" "), "\n"))
  lambda_star_i = min(which(lambda_oos_means <= max_oos_err_allowed))
  lambda_star = potential_lambdas[lambda_star_i]
  
  return(lambda_star)
}

# cv_tr_min_size: We don't want this too large as (since we have less data) otherwise we might not find the 
#                 MSE-min lambda and if the most detailed partition on the full data is best we might have a 
#                 lambda too large and choose one coarser. could choose 2
#                 On the other hand the types of partitions we generate when this param is too small will be different 
#                 and incomparable. Could choose (nfolds-2)/nfolds
#                 Therefore I take the average of the above two approaches.
# Used to warn if best lamda was the smallest, but since there's not much to do about it (we already scale 
# cv_tr_min_size), stopped reporting
# Note: We do not want to first fit the full grid and then take potential lambdas as one from each segment that 
#       picks another grid. Those lambdas aren't gauranteed to include the true lambda min. We basically 
#       roughly sampling the true  CV lambda function (which is is a step-function) and we might miss it and 
#       wrongly evaluate the benefit of each subgrid and therefore  pick the wrong one.
# ... params sent to cv_pick_lambda_f
# use cv_obj_fn if you want to a different obj fun for cv eval (rather than tr,tr training)
cv_pick_lambda <- function(y, X, d, folds_ret, nfolds, potential_lambdas=NULL, N_est=NA, min_size=5, verbosity=0, lambda_1se=FALSE, 
                           min_obs_1se=5, obj_fn, cv_tr_min_size=NA, est_plan, pr_cl=NULL, cv_bump_samples=0, bump_ratio=1, cv_obj_fn=NULL, ...) {
  #If potential_lambdas is NULL, then only have to iterate through lambda values that change partition_i (for any fold)
  #If is_obj_val_seq is monotonic then this is easy and can do sequentially, but not sure if this is the case
  supplied_lambda = !is.null(potential_lambdas)
  if(supplied_lambda) {
    n_lambda = length(potential_lambdas)
    lambda_oos = matrix(NA, nrow=nfolds, ncol=n_lambda)
  }
  else {
    lambda_ties = list()
    cvtr_fits = list()
  }
  if(verbosity>0) cat("Grid > CV: Started.\n")
  if(is.na(cv_tr_min_size)) cv_tr_min_size = as.integer(ifelse(nfolds==2, (2+min_size/2)/2, (nfolds-2)/nfolds)*min_size)
  
  params = c(list(y=y, X_d=X, d=d, folds_ret=folds_ret, nfolds=nfolds, potential_lambdas=potential_lambdas, 
                N_est=N_est, verbosity=verbosity-1,
                obj_fn=obj_fn, cv_tr_min_size=cv_tr_min_size, 
                est_plan=est_plan, cv_bump_samples=cv_bump_samples, bump_ratio=bump_ratio), list(...))
  
  col_rets = my_apply(1:nfolds, cv_pick_lambda_f, verbosity==1 || !is.null(pr_cl), pr_cl, params)
  do_bump = (length(cv_bump_samples)>1 || cv_bump_samples>0)
  
  if(is.null(cv_obj_fn)) cv_obj_fn = obj_fn
  
  # Process nfolds loop
  if(!supplied_lambda) {
    for(f in 1:nfolds) {
      cvtr_fits[[f]] = col_rets[[f]]
      lambda_ties[[f]] = get_all_lambda_ties(cvtr_fits[[f]], do_bump)  #build lambdas. Assuming no slope ties
    }
  }
  else {
    for(f in 1:nfolds) {
      lambda_oos[f, ] = col_rets[[f]]
    }
    n_cell_table = NULL
  }
  
  if(!supplied_lambda) {
    union_lambda_ties = sort(unlist(lambda_ties), decreasing=TRUE)
    mid_points = union_lambda_ties[-length(union_lambda_ties)] + diff(union_lambda_ties)/2
    potential_lambdas = c(union_lambda_ties[1]+1, mid_points, mid_points[length(mid_points)]/2)
    if(length(potential_lambdas)==0) {
      if(verbosity>0) cat("Note: CV folds consistently picked initial model (complexity didn't improve in-sample objective). Defaulting to lambda=0.\n")
      potential_lambdas=c(0)
    }
    n_lambda = length(potential_lambdas)
    lambda_oos = matrix(NA, nrow=nfolds, ncol=n_lambda)
    n_cell_table = matrix(NA, nrow=nfolds, ncol=n_lambda)
    
    for(f in 1:nfolds) {
      list[y_f_tr, y_f_cv, X_f_tr, X_f_cv, d_f_tr, d_f_cv] = split_sample_folds_m(y, X, d, folds_ret, f)
      lambda_oos[f,] = eval_lambdas(cv_obj_fn, est_plan, potential_lambdas, cvtr_fits[[f]], y_f_tr, X_f_tr, d_f_tr, y_f_cv, X_f_cv, d_f_cv, N_est, do_bump)
      if(FALSE) {
        #ns = get_num_parts(cvtr_fits[[f]])
        
        for(lambda_i in 1:length(potential_lambdas)) {
          lambda = potential_lambdas[lambda_i]
          list[partition_i, part] = get_part_for_lambda(cvtr_fits[[f]], lambda)
          n_cell_table[f, lambda_i] = num_cells(part)
          #print(paste("fit=good. lambda: ", lambda))
          #good_obj_ret = cv_obj_fn(y_f_tr, X_f_tr, d_f_tr, y_f_cv, X_f_cv, d_f_cv, N_est=N_est, partition=part, debug=TRUE, 
          #                 est_plan=est_plan, sample="trcv")
        }
      }
    }
  }
  lambda_oos_means = colMeans(lambda_oos)
  lambda_oos_min_i = which.min(lambda_oos_means)
  #if(lambda_oos_min_i==length(lambda_oos_means)) cat(paste("Warning: MSE-min lambda is the smallest (of",length(lambda_oos_means),"potential lambdas)\n"))
  if(lambda_1se) {
    lambda_star = lambda_1se_selector(potential_lambdas, lambda_oos, min_obs_1se, max_oos_err_allowed, verbosity)
  }
  else {
    lambda_star = potential_lambdas[lambda_oos_min_i]
  }
  if(verbosity>0){ 
    cat("Grid > CV: Finished.")
    cat(paste(" lambda_oos_means=[", paste(lambda_oos_means, collapse=" "), "]."))
    if(length(lambda_star)==0) cat(" Couldn't find any suitable lambdas, returning 1.\n")
    else {
      cat(paste(" potential_lambdas=[", paste(potential_lambdas, collapse=" "), "]."))
      cat(paste(" lambda_star=", lambda_star, ".\n"))
    }
  }
  if(length(lambda_star)==0) {
    lambda_star=1
  }
  
  return(list(lambda_star,lambda_oos, n_cell_table))
}
