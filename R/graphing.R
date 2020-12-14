# TODO:
# - Add marginal plots: https://www.r-graph-gallery.com/277-marginal-histogram-for-ggplot2.html
# - When more than 2-d, have the 2d graphs be the most important wones and split on the least

#' plot_2D_partition.estimated_partition
#'
#' @param grid_fit grid_fit
#' @param X_names_2D X_names_2D
#'
#' @return ggplot2 object
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
  desc_range_df = get_desc_df(grid_fit$partition, drop_unsplit=T)
  desc_range_df = do.call(cbind, lapply(desc_range_df, function(c) as.data.frame(t(matrix(unlist(c), nrow=2)))))
  
  colnames(desc_range_df)<-c("xmin", "xmax", "ymin", "ymax")
  desc_range_df["fill"] = grid_fit$cell_stats$stats$param_ests
  
  plt = ggplot2::ggplot() + 
    ggplot2::scale_x_continuous(name=X_names_2D[1]) +ggplot2::scale_y_continuous(name=X_names_2D[2]) +
    ggplot2::geom_rect(data=desc_range_df, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=fill), color="black")
  return(plt)
}
