# To run in the command-line with load_all: then run the code in the first if(FALSE), subsequent runs just run that last line of the False block

library(testthat)
loadNamespace("rprojroot")
root_dir <- rprojroot::find_package_root_file() #R cmd check doesn't copy over git and RStudio proj file

if(FALSE) { #Run manually to debug
  debugSource(paste0(rprojroot::find_testthat_root_file(),"/testres.R"))
}

library(CausalGrid)

set.seed(1337)

context("Test result")

if(getwd()==root_dir) {
  box::use(tests/dgps)
} else{ #testing
  box::use(../dgps)
}


# Mean outcome
data <- dgps$exp_data(n_4=100, dim_D=0, err_sd = 0.00)
ret1 <- fit_estimate_partition(data$y, data$X, cv_folds=2, verbosity=0)
print(ret1$splits$s_by_dim)
test_that("We get OK results (OOS)", {
  expect_lt(ret1$partition$s_by_dim[[1]][1], .6)
  expect_gt(ret1$partition$s_by_dim[[1]][1], .4)
  expect_lt(ret1$partition$s_by_dim[[2]][1], .6)
  expect_gt(ret1$partition$s_by_dim[[2]][1], .4)
})

# 
# # Treatment effect (1)
# set.seed(1337)
# data <- exp_data(n_4=100, dim_D=1, err_sd = 0.01)
# ret2 <- fit_partition(data$y, data$X, d = data$d, cv_folds=2, verbosity=0)
# #print(ret2$splits$s_by_dim)
# test_that("We get OK results (OOS)", {
#   expect_lt(ret2$partition$s_by_dim[[1]][1], .6)
#   expect_gt(ret2$partition$s_by_dim[[1]][1], .4)
#   expect_lt(ret2$partition$s_by_dim[[2]][1], .6)
#   expect_gt(ret2$partition$s_by_dim[[2]][1], .4)
# })
# 
# # Treatment effect (multiple)
# set.seed(1337)
# data <- exp_data(n_4=100, dim_D=2, err_sd = 0.01)
# ret3 <- fit_partition(data$y, data$X, d = data$d, cv_folds=2, verbosity=0)
# #print(ret2$partition$s_by_dim)
# test_that("We get OK results (OOS)", {
#   expect_lt(ret3$partition$s_by_dim[[1]][1], .6)
#   expect_gt(ret3$partition$s_by_dim[[1]][1], .4)
#   expect_lt(ret3$partition$s_by_dim[[2]][1], .6)
#   expect_gt(ret3$partition$s_by_dim[[2]][1], .4)
# })

