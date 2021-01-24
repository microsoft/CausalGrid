

#' Create 2D plots of parameter estimates
#' 
#' Creates a 2D plot of parameter estimates or a series of such slices if partition is across >2 features.
#'
#' @param grid_fit grid_fit
#' @param X_names_2D X_names_2D
#'
#' @return ggplot2 object or list of such objects
#' @export
plot.estimated_partition <- function(grid_fit, X_names_2D=NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  split_dims = (grid_fit$partition$nsplits_by_dim > 0)
  n_split_dims = sum(split_dims)
  if(n_split_dims==0) {
    print("Nothing to graph as no heterogeneity")
    return(NULL)
  }
  desc_range_df = get_desc_df.grid_partition(grid_fit$partition, drop_unsplit=TRUE, cont_bounds_inf=FALSE)
  if(n_split_dims==1) {
    desc_range_df = do.call(cbind, lapply(desc_range_df, function(c) as.data.frame(t(matrix(unlist(c), nrow=2)))))
    desc_range_df['ymin'] = 0
    desc_range_df['ymax'] = 1
    colnames(desc_range_df)<-c("xmin", "xmax", "ymin", "ymax")
    desc_range_df["estimate"] = grid_fit$cell_stats$param_ests
    xname = if(!is.null(X_names_2D)) X_names_2D[1]  else grid_fit$partition$varnames[split_dims]
    plt = ggplot2::ggplot() + 
      ggplot2::scale_x_continuous(name=xname) + 
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) + xlab(xname) +
      ggplot2::geom_rect(data=desc_range_df, mapping=ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=estimate), color="black")
    return(plt)
  }
  if(n_split_dims==2){
    if(is.null(X_names_2D)) X_names_2D = grid_fit$partition$varnames[split_dims]
    return(gen_one_plt(desc_range_df, grid_fit$cell_stats$param_ests, X_names_2D))
  } 
  
  desc_range_df_fact = data.frame(lapply(get_desc_df.grid_partition(grid_fit$partition, drop_unsplit=TRUE, do_str=TRUE), unclass))
  if(is.null(X_names_2D)){ 
    if(is.null(grid_fit$importance_weights)) {
      X_names_2D = grid_fit$partition$varnames[split_dims][1:2]
    }
    else {
      X_names_2D = grid_fit$partition$varnames[order(imp_weights, decreasing=FALSE)]
    }
  }
  other_idx = !(names(desc_range_df) %in% X_names_2D)
  n_segs_other = (grid_fit$partition$nsplits_by_dim+1)[other_idx]
  names_other = names(desc_range_df)[other_idx]
  size_other = cumprod(n_segs_other)
  test_row_equals_vec <- function(M, v) {
    rowSums(M == rep(v, each = nrow(M))) == ncol(M)
  }
  plts = list()
  for(slice_i in 1:size_other) {
    segment_indexes = segment_indexes_from_cell_i(slice_i, n_segs_other)
    row_idx = test_row_equals_vec(desc_range_df_fact[,other_idx,drop=FALSE], segment_indexes)
    #levels_desc = segment_indexes
    levels_desc = c()
    for(k in 1:length(segment_indexes)){
      levels_desc[k] = levels(desc_range_df_fact[,which(other_idx)[k]])[segment_indexes[k]]
    }
    plts[[slice_i]] = gen_one_plt(desc_range_df[row_idx,X_names_2D], grid_fit$cell_stats$param_ests[row_idx], X_names_2D) +
      ggtitle(paste(paste(names_other, levels_desc), collapse=", "))
  }
  
  return(plts)
}


gen_one_plt <- function(desc_range_df, param_ests, X_names_2D) {
  desc_range_df = do.call(cbind, lapply(desc_range_df, function(c) as.data.frame(t(matrix(unlist(c), nrow=2)))))
  
  colnames(desc_range_df)<-c("xmin", "xmax", "ymin", "ymax")
  desc_range_df["estimate"] = param_ests
  
  plt = ggplot2::ggplot() + 
    ggplot2::scale_x_continuous(name=X_names_2D[1]) +ggplot2::scale_y_continuous(name=X_names_2D[2]) +
    ggplot2::geom_rect(data=desc_range_df, mapping=ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=estimate), color="black")
  return(plt)
}
