# To run in the command-line with load_all: change do_load_all=T, then run the code in the first if(FALSE), subsequent runs just run that last line of the False block
# Undo for building project

library(testthat)
library(rprojroot)
testthat_root_dir <- rprojroot::find_testthat_root_file() #R cmd check doesn't copy over git and RStudio proj file

if(FALSE) { #Run manually to debug
  library(rprojroot)
  testthat_root_dir <- rprojroot::find_testthat_root_file()
  debugSource(paste0(testthat_root_dir,"/testrun.R"))
}

do_load_all=F
if(!do_load_all){
  library(CausalGrid)
} else {
  library(devtools)
  #devtools::load_all(".", export_all=FALSE, helpers=FALSE)
}

set.seed(1337)

source(paste0(testthat_root_dir,"/../dgps.R"))

data <- mix_data_d(n=1000)
pot_break_points = list(c(0.5), c(0))

# Just y ---------------

ret1 <- fit_estimate_partition(data$y, data$X, cv_folds=2, verbosity=0, pot_break_points=pot_break_points)
print(ret1$partition)
test_that("We get OK results (OOS)", {
  expect_equal(ret1$partition$nsplits_by_dim, c(1,1))
})


# Include d ---------------

ret1d <- fit_estimate_partition(data$y, data$X, data$d, cv_folds=2, verbosity=0, pot_break_points=pot_break_points)
print(ret1d$partition)
test_that("We get OK results (OOS)", {
  expect_equal(ret1d$partition$nsplits_by_dim, c(1,1))
})
any_sign_effect(ret1d, check_negative=T, method="fdr") #
#any_sign_effect(ret1d, check_negative=T, method="sim_mom_ineq") #the sim produces treatment effect with 0 std err, so causes problems

ret2d <- fit_estimate_partition(data$y, data$X, data$d, cv_folds=2, verbosity=0, pot_break_points=pot_break_points, ctrl_method="all")
print(ret2d$partition)
#TODO: Should I check this?
#test_that("We get OK results (OOS)", {
#  expect_equal(ret2d$partition$nsplits_by_dim, c(1,1))
#})

ret3d <- fit_estimate_partition(data$y, data$X, data$d, cv_folds=3, verbosity=0, pot_break_points=pot_break_points, ctrl_method="LassoCV")
print(ret3d$partition)
#TODO: Should I check this?
#test_that("We get OK results (OOS)", {
#  expect_equal(ret3d$partition$nsplits_by_dim, c(1,1))
#})

ret4d <- fit_estimate_partition(data$y, data$X, data$d, cv_folds=2, verbosity=0, pot_break_points=pot_break_points, ctrl_method="rf")
print(ret4d$partition)
#TODO: Should I check this?
#test_that("We get OK results (OOS)", {
#  expect_equal(ret4d$partition$nsplits_by_dim, c(1,1))
#})

ret1db <- fit_estimate_partition(data$y, data$X, data$d, cv_folds=2, verbosity=0, pot_break_points=pot_break_points, bump_B=2)


ret1dc <- fit_estimate_partition(data$y, data$X, data$d, cv_folds=2, verbosity=0, pot_break_points=pot_break_points, importance_type="single")

X_3 = data$X
X_3$X3 = data$X$X2
pot_break_points_3 = pot_break_points
pot_break_points_3[[3]] = pot_break_points[[2]]
print("---------------")
ret1dd <- fit_estimate_partition(data$y, X_3, data$d, cv_folds=2, verbosity=2, pot_break_points=pot_break_points_3, importance_type="interaction", bump_B=3)
print("---------------")
ret1dd <- fit_estimate_partition(data$y, X_3, data$d, cv_folds=2, verbosity=1, pot_break_points=pot_break_points_3, importance_type="interaction", bump_B=3)
print("---------------")
ret1dd <- fit_estimate_partition(data$y, X_3, data$d, cv_folds=2, verbosity=0, pot_break_points=pot_break_points_3, importance_type="interaction", bump_B=3)

# Old test
if(FALSE) {
  dim_D = 1
  n_4 = 100
  data <- exp_data(n_4=n_4, dim_D, err_sd = 1e-7)
  K = ncol(data$X)
  # Limit break points?
  # pot_break_points = list()
  # g=4
  # for(k in 1:ncol(data$X)) {
  #   pot_break_points[[k]] = quantile(data$X[,k], seq(0,1,length.out=g+1))[-c(g+1,1)]
  # }
  pot_break_points = NULL
  tr_index = c(sample(n_4, n_4/2), sample(n_4, n_4/2)+n_4, sample(n_4, n_4/2) + 2*n_4, sample(n_4, n_4/2) + 3*n_4)
  X = as.data.frame(data$X)
  X[[2]] = factor(c("a", "a", "b", "c"))
  ret2 <- fit_estimate_partition(data$y, X, tr_split = tr_index, cv_folds=5, max_splits=Inf, verbosity=1, pot_break_points=pot_break_points, d=data$d, ctrl_method="LassoCV") #bucket_min_d_var, bucket_min_n
  print(ret2)
  cat(paste("s_by_dim", paste(ret2$partition$s_by_dim, collapse=" "),"\n"))
  cat(paste("lambda", paste(ret2$lambda, collapse=" "),"\n"))
  cat(paste("param_ests", paste(ret2$cell_stats$param_ests, collapse=" "),"\n"))
  cat(paste("var_ests", paste(ret2$cell_stats$var_ests, collapse=" "),"\n"))
  cat(paste("cell_sizes", paste(ret2$cell_stats$cell_sizes, collapse=" "),"\n"))
  
  #View implied model
  est_df = data.frame(y=data$y, f = get_factor_from_partition(ret2$partition, data$X)) 
  if (dim_D) est_df = cbind(est_df, data$d)
  if(dim_D==0) {
    ols_fit = lm(y~0+f, data=est_df)
  }
  if (dim_D==1) {
    ols_fit = lm(y~0+f+d:f, data=est_df)
  }
  if (dim_D==2) {
    ols_fit = lm(y~0+f+(d1+d2):f, data=est_df)
  }
  print(summary(ols_fit))
  #0.5003187, 0.5000464
  
  
  #Compare to manually-specified split
  my_partition = add_split.grid_partition(add_split.grid_partition(grid_partition(get_X_range(data$X)),partition_split(1, .5)), partition_split(2, .5))
  y_tr = data$y[ret2$index_tr]
  X_tr = data$X[ret2$index_tr, , drop=FALSE]
  d_tr = data$d[ret2$index_tr, , drop=FALSE]
  X_es = data$X[-ret2$index_tr, , drop=FALSE]
  N_est = nrow(X_es)
  my_part_mse = mse_hat_obj(y_tr,X_tr,d_tr, N_est=N_est, partition=my_partition)
  print(paste("emse:",my_part_mse,". +pen", my_part_mse + ret2$lambda*(num_cells(my_partition)-1)))
  est_df[['f2']] = interaction(get_factors_from_partition(my_partition, data$X))
  if(dim_D==0) {
    ols_fit = lm(y~0+f2, data=est_df)
  }
  if(dim_D==1) {
    ols_fit = lm(y~0+f2+d:f2, data=est_df)
  }
  if(dim_D==2) {
    ols_fit = lm(y~0+f2+(d1+d2):f2, data=est_df)
  }
  print(summary(ols_fit))
  
}
