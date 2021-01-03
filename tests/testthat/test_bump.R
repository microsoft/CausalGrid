# To run in the command-line with load_all: then run the code in the first if(FALSE), subsequent runs just run that last line of the False block

library(testthat)
library(rprojroot)
testthat_root_dir <- rprojroot::find_testthat_root_file() #R cmd check doesn't copy over git and RStudio proj file

if(FALSE) { #Run manually to debug
  library(rprojroot)
  testthat_root_dir <- rprojroot::find_testthat_root_file()
  debugSource(paste0(testthat_root_dir,"/test_bump.R"))
}

library(CausalGrid)

set.seed(1337)

context("Test bumping")

source(paste0(testthat_root_dir,"/../dgps.R"))

data_xor = XOR_sim(n=500)

#Show how bumping helps with fixed depth
ret_nobmp_fixed <- fit_estimate_partition(data_xor$y, data_xor$X, data_xor$w, verbosity=0, 
                                          partition_i=3)
print(ret_nobmp_fixed)
ret_bmpa_fixed <- fit_estimate_partition(data_xor$y, data_xor$X, data_xor$w, verbosity=0, 
                                         partition_i=3, bump_samples=20, bump_complexity=list(doCV=FALSE, incl_comp_in_pick=FALSE))
print(ret_bmpa_fixed)
ret_bmpb_fixed <- fit_estimate_partition(data_xor$y, data_xor$X, data_xor$w, verbosity=0, 
                                         partition_i=3, bump_samples=20, bump_complexity=list(doCV=FALSE, incl_comp_in_pick=TRUE))
print(ret_bmpb_fixed)

#Then add in CV complexity
ret_nobmp_cv <- fit_estimate_partition(data_xor$y, data_xor$X, data_xor$w, verbosity=0, 
                                       cv_folds=2)
print(ret_nobmp_cv)
ret_bmpa_cv <- fit_estimate_partition(data_xor$y, data_xor$X, data_xor$w, verbosity=0, 
                                   cv_folds=2, bump_samples=20, bump_complexity=list(doCV=FALSE, incl_comp_in_pick=FALSE)) 
print(ret_bmpa_cv)
ret_bmpb_cv <- fit_estimate_partition(data_xor$y, data_xor$X, data_xor$w, verbosity=0, 
                                   cv_folds=2, bump_samples=20, bump_complexity=list(doCV=FALSE, incl_comp_in_pick=TRUE)) 
print(ret_bmpb_cv)
ret_bmpc_cv <- fit_estimate_partition(data_xor$y, data_xor$X, data_xor$w, verbosity=0, 
                                   cv_folds=2, bump_samples=20, bump_complexity=list(doCV=TRUE, incl_comp_in_pick=FALSE)) 
print(ret_bmpc_cv)
ret_bmpd_cv <- fit_estimate_partition(data_xor$y, data_xor$X, data_xor$w, verbosity=0, 
                                   cv_folds=2, bump_samples=20, bump_complexity=list(doCV=TRUE, incl_comp_in_pick=TRUE)) 
print(ret_bmpd_cv)
