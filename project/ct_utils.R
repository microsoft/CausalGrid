library(ggplot2)
library(causalTree)

rpart_label_component <- function(object) {
  #TODO: Didn't copy the part with categorical variables, so this might not work then. Could do that.
  #Copied from print.rpart
  ff <- object$frame
  n <- nrow(ff)
  is.leaf <- (ff$var == "<leaf>")
  whichrow <- !is.leaf
  index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + !is.leaf))
  irow <- index[c(whichrow, FALSE)]
  ncat <- object$splits[irow, 2L]
  jrow <- irow[ncat < 2L]
  cutpoint <- object$splits[jrow, 4L]
  lsplit <- rsplit <- numeric(length(irow))
  lsplit[ncat<2L] <- cutpoint
  rsplit[ncat<2L] <- cutpoint
  
  
  vnames <- ff$var[whichrow]
  varname <- (as.character(vnames))
  node <- as.numeric(row.names(ff))
  parent <- match(node %/% 2L, node[whichrow])
  odd <- (as.logical(node %% 2L))
  labels_var <- character(n)
  labels_num <- numeric(n)
  labels_num[odd] <- rsplit[parent[odd]]
  labels_num[!odd] <- lsplit[parent[!odd]]
  labels_num[1L] <- NA
  labels_var[odd] <- varname[parent[odd]]
  labels_var[!odd] <- varname[parent[!odd]]
  labels_var[1L] <- "root"
  list(labels_var, labels_num)
}

plot_partition.rpart <- function(yvals, varnames, x_min, x_max, y_min, y_max, labels_var, labels_num, depth) {
  nnode = length(depth)
  if(nnode<=1) {
    return(data.frame(xmin=c(x_min), xmax=c(x_max), ymin=c(y_min), ymax=c(y_max), fill=c(yvals[1])))
  }
  yvals = yvals[2:nnode]
  labels_var = labels_var[2:nnode]
  labels_num = labels_num[2:nnode]
  depth = depth[2:nnode]
  nnode = length(depth)
  i1 = which(depth==0)[1]
  i2 = which(depth==0)[2]
  varname = labels_var[1]
  dim = which(varnames==varname)[1]
  cutoff = labels_num[1]
  if(dim==1) {
    x_max1 = cutoff
    x_min2 = cutoff
    y_max1 = y_max
    y_min2 = y_min
  }
  else {
    x_max1 = x_max
    x_min2 = x_min
    y_max1 = cutoff
    y_min2 = cutoff
  }
  
  ret1 = plot_partition.rpart(yvals[1:(i2-1)], varnames, x_min, x_max1, y_min, y_max1, labels_var[1:(i2-1)], labels_num[1:(i2-1)], depth[1:(i2-1)]-1)
  ret2 = plot_partition.rpart(yvals[i2:nnode], varnames, x_min2, x_max, y_min2, y_max, labels_var[i2:nnode], labels_num[i2:nnode], depth[i2:nnode]-1)
  
  return(rbind(ret1, ret2))
}


#' plot_2D_partition.rpart
#'
#' @param cart_fit rpart fit object
#' @param X_range X_range
#' @param cs color_list
#'
#' @return ggplot fig
plot_2D_partition.rpart <- function(cart_fit, X_range) {
  #Note: plotmo doesn't work well because it's just a grid and doesn't find the boundaries
  varnames = names(cart_fit$ordered)
  node <- as.numeric(row.names(cart_fit$frame))
  yvals = cart_fit$frame$yval
  depth <- rpart:::tree.depth(node)
  list[labels_var, labels_num] <- rpart_label_component(cart_fit)
  nnode = length(depth)
  rects = plot_partition.rpart(yvals, varnames, X_range[[1]][1], X_range[[1]][2], X_range[[2]][1], X_range[[2]][2], labels_var, labels_num, depth-1)
  
  plt = ggplot() + 
    scale_x_continuous(name=varnames[1]) +scale_y_continuous(name=varnames[2]) +
    geom_rect(data=rects, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=fill), color="black")
  return(plt)
}

ct_cv_tree <- function(form, data, treatment, index_tr=NULL, tr_split=NA, split.Honest=TRUE, cv.Honest=TRUE, 
                       minsize=2L, split.Bucket=FALSE, bucketNum=5, xval=10) {
  N = nrow(data)
  if(is.null(index_tr)) {
    if(is.na(tr_split)) tr_split=0.5
    index_tr = sample(N, tr_split*N)
  }
  #could've done causalTree() and then estimate.causalTree
  #fn_ret <- capture.output(ctree<-causalTree(form, data = data, treatment = treatment,
  #                    split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, 
  #                    xval = 2, cp = 0, minsize=minsize),
  #                    type="output") #does some random output
  #print(fn_ret)
  #opcp <- ctree$cptable[,1][which.min(ctree$cptable[,4])]
  #ct_opfit <- prune(ctree, opcp)
  
  split.alpha = if(split.Honest) 0.5 else 1
  
  fn_ret <- capture.output(honestTree <- honest.causalTree(form, data = data[index_tr,], treatment = treatment[index_tr],
                                                           est_data = data[-index_tr,],
                                                           est_treatment = treatment[-index_tr],
                                                           split.Rule = "CT", split.Honest = split.Honest,
                                                           HonestSampleSize = nrow(data[-index_tr,]),
                                                           split.Bucket = split.Bucket, bucketNum=bucketNum,
                                                           cv.option = "CT", cv.Honest = cv.Honest, minsize=minsize, 
                                                           split.alpha=split.alpha, xval=xval))
  #print(fn_ret)
  opcp <- honestTree$cptable[,1][which.min(honestTree$cptable[,4])]
  opTree <- prune(honestTree, opcp)
  
  return(opTree)
}

num_cells.rpart <- function(obj){
  sum(obj$frame[["var"]]=='<leaf>')
}

ct_nsplits_by_dim <- function(obj, ndim) {
  library(stringr)
  strs = paste(obj$frame$var[obj$frame$var!="<leaf>"])
  int_tbl = table(as.integer(str_sub(strs, start=2)))
  ret = rep(0, ndim)
  for(k in 1:ndim) {
    k_str = as.character(k)
    if(k_str %in% names(int_tbl))
      ret[k] = int_tbl[[k_str]]
  }
  return(ret)
}

#Just nodes and treatment effects
ct_desc <- function(ct_m, tex_table=TRUE, digits=3) {
  ct_m_desc <- capture.output(print(ct_m))
  ct_m_desc = ct_m_desc[-c(1:5)]
  new_str = c()
  for(i in 1:length(ct_m_desc)) {
    non_num_init = str_extract(ct_m_desc[i], paste0("^[ ]*[:digit:]+[)] (root|[:alnum:]*(< |>=))[ ]*"))
    nums = as.numeric(str_split(str_sub(ct_m_desc[i], start=str_length(non_num_init)+1, end=str_length(ct_m_desc[i])-2), " ")[[1]])
    node_path = if(i==1) non_num_init else paste(non_num_init, format(nums[1], digits=digits))
    str_effect = format(nums[length(nums)], digits=digits)
    is_leaf = str_sub(ct_m_desc[i],start=str_length(ct_m_desc[i]))=="*"
    if(tex_table) {
      n_spaces = str_length(str_extract(node_path, "^[ ]*"))
      node_path = paste0(paste(replicate(n_spaces, "~"), collapse = ""), str_sub(node_path, start=n_spaces))
      if(is_leaf)
        new_str[i] = paste0(node_path, " & ", str_effect, " \\\\")
      else
        new_str[i] = paste0(node_path, " &  \\\\")
    }
    else {
      new_str[i] = paste(node_path, str_effect)
      if(is_leaf) new_str[i]= paste(new_str[i], "*")
    }
  }
  #cat(file_cont, file=o_fname)
  if(tex_table) {
    new_str = c("\\begin{tabular}{lr}", "  \\hline", "Node & Est. \\\\", "  \\hline", new_str, "  \\hline", "\\end{tabular}")
  }
  return(new_str)
}

