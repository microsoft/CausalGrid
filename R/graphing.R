# TODO:
# - Add marginal plots: https://www.r-graph-gallery.com/277-marginal-histogram-for-ggplot2.html
# - When more than 2-d, have the 2d graphs be the most important wones and split on the least

gen_one_plt <- function(desc_range_df, param_ests, X_names_2D) {
  desc_range_df = do.call(cbind, lapply(desc_range_df, function(c) as.data.frame(t(matrix(unlist(c), nrow=2)))))
  
  colnames(desc_range_df)<-c("xmin", "xmax", "ymin", "ymax")
  desc_range_df["fill"] = param_ests
  
  plt = ggplot2::ggplot() + 
    ggplot2::scale_x_continuous(name=X_names_2D[1]) +ggplot2::scale_y_continuous(name=X_names_2D[2]) +
    ggplot2::geom_rect(data=desc_range_df, mapping=ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=fill), color="black")
  return(plt)
}

#' Create 2D plots of parameter estimates
#' 
#' Creates a 2D plot of parameter estimates or a series of such slices if partition is across >2 features.
#'
#' @param grid_fit grid_fit
#' @param X_names_2D X_names_2D
#'
#' @return ggplot2 object or list of such objects
#' @export
plot_2D_partition.estimated_partition <- function(grid_fit, X_names_2D) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  split_dims = (grid_fit$partition$nsplits_by_dim > 0)
  n_split_dims = sum(split_dims)
  if(n_split_dims<2) {
    warning("Less than 2 dimensions of heterogeneity")
  }
  desc_range_df = get_desc_df.grid_partition(grid_fit$partition, drop_unsplit=TRUE)
  if(n_split_dims==2) return(gen_one_plt(desc_range_df, grid_fit$cell_stats$stats$param_ests))
  
  desc_range_df_fact = data.frame(lapply(get_desc_df.grid_partition(grid_fit$partition, drop_unsplit=TRUE, do_str=TRUE), unclass))
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
    plts[[slice_i]] = gen_one_plt(desc_range_df[row_idx,X_names_2D], grid_fit$cell_stats$stats$param_ests[row_idx], X_names_2D) +
      ggtitle(paste(paste(names_other, "= level", segment_indexes), collapse=", "))
  }
  
  return(plts)
}
