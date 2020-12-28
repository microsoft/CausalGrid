# Notes:
# - To change parallel/single-threaded, search for PARALLEL comments
# - There is a "suffix" for output files so that one can test w/o overwritting

library(gsubfn)
library(devtools)
suppressPackageStartupMessages(library(data.table))
library(CausalGrid)
library(causalTree)
library(doParallel)
library(foreach)
library(xtable)
library(stringr)
library(glmnet)
library(ranger)
library(gridExtra)
library(ggplot2)
source("project/ct_utils.R")
source("tests/dgps.R")

#Paths
export_dir = "C:/Users/bquist/Dropbox/Causal Grid/writeup/" #can be turned off below
sim_rdata_fname = "project/sim.RData"
log_file = "project/log.txt"
tbl_export_path = paste0(export_dir, "tables/")
fig_export_path = paste0(export_dir, "figs/")
suffix = "-temp"

#Estimation config
b = 4
nfolds = 5
minsize=25

#Sim parameters
S = 3 #3 100, TODO: Is this just 1!?
Ns = c(500, 1000) #c(500, 1000)
D = 4 #3
Ks = c(2, 2, 10, 20)
N_test = 8000
NN = length(Ns)
NIters = NN*D*S

good_features = list(c(T,T),c(T, F), c(T, T, rep(F, 8)), c(rep(T, 4), rep(F, 16)))

#Execution config
n_parallel = 1 #1 5; PARALLEL
my_seed = 1337
set.seed(my_seed)
rf.num.threads = 1 #NULL will multi-treatd, doesn't seem to help much with small data



# Helper functions --------
sim_data <- function(n=500, design=1) {
  if(design<=3) return(AI_sim(n, design))
  return(XOR_sim(n))
}

get_iter <- function(d, N_i, s) {
  (d-1)*NN*S + (N_i-1)*S + (s-1) + 1
}

yX_data <- function(y, X) {
  yX = cbind(y, X)
  colnames(yX) = c("Y", paste("X", 1:ncol(X), sep=""))
  yX = as.data.frame(yX)
}

sim_ct_fit <- function(y, X, w, tr_sample, honest=FALSE) {
  yX = yX_data(y, X)
  set.seed(my_seed)
  fit = ct_cv_tree("Y~.", data=yX, treatment=w, index_tr=tr_sample, minsize=minsize, split.Bucket=TRUE, bucketNum=b, xval=nfolds, split.Honest=honest, cv.Honest=honest)
  attr(fit$terms, ".Environment") <- NULL #save space other captures environment
  return(fit)
}

sim_ct_predict_te <- function(obj, y_te, X_te) {
  yX_te = yX_data(y_te, X_te)
  return(predict(obj, newdata=yX_te, type="vector"))
}

sim_eval_ct <- function(data1, data2, good_mask, honest=FALSE) {
  list[y, X, w, tau] = data1
  N = nrow(X)/2
  tr_sample = c(rep(TRUE, N), rep(FALSE, N))
  ct_fit = sim_ct_fit(y, X, w, tr_sample, honest=honest)
  nl = num_cells(ct_fit)
  nsplits_by_dim = ct_nsplits_by_dim(ct_fit, ncol(X))
  ngood = sum(nsplits_by_dim[good_mask])
  ntot = sum(nsplits_by_dim)
  list[y_te, X_te, w_te, tau_te] = data2
  mse = mean((tau_te - sim_ct_predict_te(ct_fit, y_te, X_te))^2)
  return(list(ct_fit, nl, mse, ngood, ntot))
}

sim_cg_fit <- function(y, X, w, tr_sample, verbosity=0, honest=FALSE, do_rf=FALSE, num.threads=rf.num.threads, ...) {
  set.seed(my_seed)
  if(do_rf) {
    return(fit_estimate_partition(y, X, d=w, tr_split=tr_sample, cv_folds=nfolds, verbosity=verbosity, min_size=2*minsize, max_splits=10, bucket_min_d_var=TRUE, bucket_min_n=2*b, honest=honest, ctrl_method=grid_rf(num.threads=num.threads), ...))
  }
  else {
    return(fit_estimate_partition(y, X, d=w, tr_split=tr_sample, cv_folds=nfolds, verbosity=verbosity, min_size=2*minsize, max_splits=10, bucket_min_d_var=TRUE, bucket_min_n=2*b, honest=honest, ...))
  }
}

sim_cg_vectors <- function(grid_fit, mask, X_te, tau_te) {
  nsplits = grid_fit$partition$nsplits_by_dim
  preds = predict_te.estimated_partition(grid_fit, new_X=X_te)
  cg_mse = mean((preds - tau_te)^2)
  return(c(num_cells(grid_fit), cg_mse, sum(nsplits[mask]), sum(nsplits)))
}

ct_ndim_splits <- function(obj, K) {
  return(sum(ct_nsplits_by_dim(obj, K)>0))
}

cg_ndim_splits <- function(obj) {
  return(sum(obj$partition$nsplits_by_dim>0))
}

# Generate Data --------

cat("Generating Data\n")
data1 = list()
data2 = list()
for(d in 1:D) {
  for(N_i in 1:NN) {
    N = Ns[N_i]
    for(s in 1:S){
      iter = get_iter(d, N_i, s)
      data1[[iter]] = sim_data(n=2*N, design=d)
      data2[[iter]] = sim_data(n=N_test, design=d)
    }
  }
}

# Eval CT ------
cat("Eval CT\n")

results_ct_h = matrix(0,nrow=0, ncol=7)
results_ct_a = matrix(0,nrow=0, ncol=7)
ct_h_fit_models = list()
ct_a_fit_models = list()
for(d in 1:D) {
  for(N_i in 1:NN) {
    for(s in 1:S){
      iter = get_iter(d, N_i, s)
      list[ct_h_fit, nl_h, mse_h, ngood_h, ntot_h] = sim_eval_ct(data1[[iter]], data2[[iter]], good_features[[d]], honest=T)
      results_ct_h = rbind(results_ct_h, c(d, N_i, s, nl_h, mse_h, ngood_h, ntot_h))
      ct_h_fit_models[[iter]] = ct_h_fit
      list[ct_a_fit, nl, mse, ngood, ntot] = sim_eval_ct(data1[[iter]], data2[[iter]], good_features[[d]])
      results_ct_a = rbind(results_ct_a, c(d, N_i, s, nl, mse, ngood, ntot))
      ct_a_fit_models[[iter]] = ct_a_fit
    }
  }
}
ct_a_nl = results_ct_a[,4]
ct_h_nl = results_ct_h[,4]

# Eval CG -----
cat("Eval CG\n")


if(n_parallel>1) {
  if(file.exists(log_file)) file.remove(log_file)
  cl <- makeCluster(n_parallel, outfile=log_file)
  registerDoParallel(cl)
}

bar_length = if(n_parallel>1) NN*D else S*NN*D
t1 = Sys.time()
cat(paste("Start time: ",t1,"\n"))
pb = utils::txtProgressBar(0, bar_length, style = 3)
run=1
N_bumps = 20
outer_results = list()
for(d in 1:D) {
  for(N_i in 1:NN) {
    if(n_parallel>1) {
      utils::setTxtProgressBar(pb, run)
      run = run+1
    }
    results_s = list() #PARALLEL: comment-out
    for(s in 1:S){  #PARALLEL: comment-out, uncomment next
    #results_s = foreach(s=1:S, .packages=c("proto","gsubfn","rpart", "rpart.plot", "data.table","causalTree", "ranger",  "lattice", "ggplot2", "caret", "Matrix", "foreach", "CausalGrid"), .errorhandling = "pass") %dopar% { #, .combine=rbind
      if(n_parallel==1) {
        utils::setTxtProgressBar(pb, run)
        run = run+1
      }
      res = c(s)
      N = Ns[N_i]
      iter = get_iter(d, N_i, s)
      list[y, X, w, tau] = data1[[iter]]
      list[y_te, X_te, w_te, tau_te] = data2[[iter]]
      tr_sample = c(rep(TRUE, N), rep(FALSE, N))
      
      
      grid_a_fit <- sim_cg_fit(y, X, w, tr_sample, honest=FALSE)
      grid_a_LassoCV_fit <- sim_cg_fit(y, X, w, tr_sample, honest=FALSE, ctrl_method="LassoCV")
      grid_a_RF_fit <- sim_cg_fit(y, X, w, tr_sample, honest=FALSE, do_rf=TRUE)
      grid_a_LassoCV_b_fit <- sim_cg_fit(y, X, w, tr_sample, honest=FALSE, ctrl_method="LassoCV", bump_samples=N_bumps)
      grid_a_RF_b_fit <- sim_cg_fit(y, X, w, tr_sample, honest=FALSE, do_rf=TRUE, bump_samples=N_bumps)

      grid_a_am_fit <- change_complexity(grid_a_fit, y, X, d=w, which.min(abs(ct_a_nl[iter] - (grid_a_fit$complexity_seq + 1))))
      grid_a_LassoCV_b_hm_fit <- change_complexity(grid_a_LassoCV_b_fit, y, X, d=w, which.min(abs(ct_h_nl[iter] - (grid_a_LassoCV_b_fit$complexity_seq + 1))))
      grid_a_RF_b_hm_fit <- change_complexity(grid_a_RF_b_fit, y, X, d=w, which.min(abs(ct_h_nl[iter] - (grid_a_RF_b_fit$complexity_seq + 1))))
      
      res = c(res, sim_cg_vectors(grid_a_fit, good_features[[d]], X_te, tau_te))
      res = c(res, sim_cg_vectors(grid_a_LassoCV_fit, good_features[[d]], X_te, tau_te))
      res = c(res, sim_cg_vectors(grid_a_RF_fit, good_features[[d]], X_te, tau_te))
      res = c(res, sim_cg_vectors(grid_a_LassoCV_b_fit, good_features[[d]], X_te, tau_te))
      res = c(res, sim_cg_vectors(grid_a_RF_b_fit, good_features[[d]], X_te, tau_te))
      res = c(res, sim_cg_vectors(grid_a_am_fit, good_features[[d]], X_te, tau_te))
      res = c(res, sim_cg_vectors(grid_a_LassoCV_b_hm_fit, good_features[[d]], X_te, tau_te))
      res = c(res, sim_cg_vectors(grid_a_RF_b_hm_fit, good_features[[d]], X_te, tau_te))
      
      #Save space
      grid_a_RF_fit$est_plan$rf_y_fit <- grid_a_RF_fit$est_plan$rf_d_fit <- NULL
      grid_a_RF_b_fit$est_plan$rf_y_fit <- grid_a_RF_b_fit$est_plan$rf_d_fit <- NULL
      res = list(grid_a_fit, grid_a_LassoCV_fit, grid_a_RF_fit, grid_a_LassoCV_b_fit, grid_a_RF_b_fit, res)
      
      results_s[[s]] = res #PARALLEL: comment-out, uncomment next
      #res
    }
    outer_results[[(d-1)*NN + (N_i-1) + 1]] = results_s
  }
}

t2 = Sys.time() #can us as.numeric(t1) to convert to seconds
td = t2-t1
close(pb)
cat(paste("Total time: ",format(as.numeric(td))," ", attr(td,"units"),"\n"))

if(n_parallel>1) stopCluster(cl)


# Collect results ----

cg_a_fit_models = list()
cg_a_LassoCV_fit_models = list()
cg_a_RF_fit_models = list()
cg_a_LassoCV_b_fit_models = list()
cg_a_RF_b_fit_models = list()
results_cg_a = matrix(0,nrow=0, ncol=7)
results_cg_a_LassoCV = matrix(0,nrow=0, ncol=7)
results_cg_a_RF = matrix(0,nrow=0, ncol=7)
results_cg_a_LassoCV_b = matrix(0,nrow=0, ncol=7)
results_cg_a_RF_b = matrix(0,nrow=0, ncol=7)
results_cg_a_am = matrix(0,nrow=0, ncol=7)
results_cg_a_LassoCV_b_hm = matrix(0,nrow=0, ncol=7)
results_cg_a_RF_b_hm = matrix(0,nrow=0, ncol=7)
n_errors = matrix(0, nrow=D, ncol=NN)
for(d in 1:D) {
  for(N_i in 1:NN) {
    results_s = outer_results[[(d-1)*NN + (N_i-1) + 1]]
    for(s in 1:S) {
      res = results_s[[s]]
      if(inherits(res, "error")) {
        cat(paste("Error: d=", d, "N_i=", N_i, "s=", s, "\n"))
        n_errors[d,N_i] = n_errors[d,N_i] + 1
        next
      }
      iter = get_iter(d, N_i, s)
      cg_a_fit_models[[iter]] = res[[1]]
      cg_a_LassoCV_fit_models[[iter]] = res[[2]]
      cg_a_RF_fit_models[[iter]] = res[[3]]
      cg_a_LassoCV_b_fit_models[[iter]] = res[[4]]
      cg_a_RF_b_fit_models[[iter]] = res[[5]]
      res = res[[6]]
      #s = res[1]
      results_cg_a = rbind(results_cg_a, c(d, N_i, s, res[2:5]))
      results_cg_a_LassoCV = rbind(results_cg_a_LassoCV, c(d, N_i, s, res[6:9]))
      results_cg_a_RF = rbind(results_cg_a_RF, c(d, N_i, s, res[10:13]))
      results_cg_a_LassoCV_b = rbind(results_cg_a_LassoCV, c(d, N_i, s, res[14:17]))
      results_cg_a_RF_b = rbind(results_cg_a_RF, c(d, N_i, s, res[18:21]))
      results_cg_a_am = rbind(results_cg_a_am, c(d, N_i, s, res[22:25]))
      results_cg_a_LassoCV_b_hm = rbind(results_cg_a_LassoCV_b_hm, c(d, N_i, s, res[26:29]))
      results_cg_a_RF_b_hm = rbind(results_cg_a_RF_b_hm, c(d, N_i, s, res[30:33]))
    }
  }
}
if(sum(n_errors)>0){
  cat("N errors")
  print(n_errors)
}


if(F){ #Output raw results
  save(S, Ns, D, data1, data2, 
       ct_h_fit_models, ct_a_fit_models, 
       cg_a_fit_models, cg_a_LassoCV_fit_models, cg_a_RF_fit_models, cg_a_LassoCV_b_fit_models, cg_a_RF_b_fit_models,
       results_ct_h, results_ct_a, 
       results_cg_a, results_cg_a_LassoCV, results_cg_a_RF, results_cg_a_LassoCV_b, results_cg_a_RF_b, results_cg_a_am, results_cg_a_LassoCV_b_hm, results_cg_a_RF_b_hm,
       file = sim_rdata_fname)
}
if(F){
  load(sim_rdata_fname)
  ct_a_nl = results_ct_a[,4]
  ct_h_nl = results_ct_h[,4]
}

# Output Results ------


sum_res = function(full_res) {
  nl = matrix(0, nrow=D, ncol=NN)
  mse = matrix(0, nrow=D, ncol=NN)
  pct_good = matrix(0, nrow=D, ncol=NN)
  for(d in 1:D) {
    for(N_i in 1:NN) {
      start = get_iter(d, N_i, 1)
      end = get_iter(d, N_i, S)
      nl[d,N_i] = mean(full_res[start:end,4])
      mse[d,N_i] = mean(full_res[start:end,5])
      pct_good[d,N_i] = sum(full_res[start:end,6])/sum(full_res[start:end,7])
    }
  }
  return(list(nl, mse, pct_good))
}

list[nl_CT_h, mse_CT_h, pct_good_CT_h] = sum_res(results_ct_h)
list[nl_CT_a, mse_CT_a, pct_good_CT_a] = sum_res(results_ct_a)
list[nl_CG_a, mse_CG_a, pct_good_CG_a] = sum_res(results_cg_a)
list[nl_CG_a_LassoCV, mse_CG_a_LassoCV, pct_good_CG_a_LassoCV] = sum_res(results_cg_a_LassoCV)
list[nl_CG_a_RF, mse_CG_a_RF, pct_good_CG_a_RF] = sum_res(results_cg_a_RF)
list[nl_CG_a_LassoCV_b, mse_CG_a_LassoCV_b, pct_good_CG_a_LassoCV_b] = sum_res(results_cg_a_LassoCV_b)
list[nl_CG_a_RF_b, mse_CG_a_RF_b, pct_good_CG_a_RF_b] = sum_res(results_cg_a_RF_b)
list[nl_CG_a_am, mse_CG_a_am, pct_good_CG_a_am] = sum_res(results_cg_a_am)
list[nl_CG_a_LassoCV_b_hm, mse_CG_a_LassoCV_b_hm, pct_good_CG_a_LassoCV_b_hm] = sum_res(results_cg_a_LassoCV_b_hm)
list[nl_CG_a_RF_b_hm, mse_CG_a_RF_b_hm, pct_good_CG_a_RF_b_hm] = sum_res(results_cg_a_RF_b_hm)


flatten_table <- function(mat) {
  new_mat = cbind(mat[1,, drop=F], mat[2,, drop=F], mat[3,, drop=F], mat[4,, drop=F])
  colnames(new_mat) = rep(c("N=500", "N=1000"), D)
  new_mat
}
compose_table <- function(mat_CT_h, mat_CT_a, mat_CG_a_am, mat_CG_a, mat_CG_a_LassoCV, mat_CG_a_RF,
                          mat_CG_a_LassoCV_b, mat_CG_a_RF_b, mat_CG_a_LassoCV_b_hm, mat_CG_a_RF_b_hm) {
  new_mat = rbind(flatten_table(mat_CT_h), flatten_table(mat_CT_a), 
                  flatten_table(mat_CG_a_am), flatten_table(mat_CG_a), 
                  flatten_table(mat_CG_a_LassoCV), flatten_table(mat_CG_a_RF),
                  flatten_table(mat_CG_a_LassoCV_b), flatten_table(mat_CG_a_RF_b),
                  flatten_table(mat_CG_a_LassoCV_b_hm), flatten_table(mat_CG_a_RF_b_hm))
  rownames(new_mat) = c("Causal Tree (CT)", "Causal Tree - Adaptive (CT-A)", 
                        "Causal Grid - Matching CT-A complexity (CG(M))", "Causal Grid (CG)", 
                        "CG w/ Linear Controls (CG-X)", "CG w/ RF Controls (CG-RF)", 
                        "CG-X w/ Bumping (CG-X-B)", "CG-RF w/ Bumping (CG-RF-B)", 
                        "CG-X-B - Matching CT complexity (CG-X-B-HM)", "CG-RF-B - Matching CT complexity (CG-RF-B-HM)")
  new_mat
}

fmt_table <- function(xtbl, o_fname) {
  capt_ret <- capture.output(file_cont <- print(xtbl, floating=F, comment = F))
  file_cont = paste0(str_sub(file_cont, end=35), " &\\multicolumn{2}{c}{Design 1}&\\multicolumn{2}{c}{Design 2}&\\multicolumn{2}{c}{Design 3}&\\multicolumn{2}{c}{Design 4}\\\\ \n", str_sub(file_cont, start=36))
  cat(file_cont, file=o_fname)
}

n_cells_comp = compose_table(nl_CT_h, nl_CT_a, nl_CG_a_am, nl_CG_a, nl_CG_a_LassoCV, nl_CG_a_RF,
                             nl_CG_a_LassoCV_b, nl_CG_a_RF_b, nl_CG_a_LassoCV_b_hm, nl_CG_a_RF_b_hm)
mse_comp = compose_table(mse_CT_h, mse_CT_a, mse_CG_a_am, mse_CG_a, mse_CG_a_LassoCV, mse_CG_a_RF,
                         mse_CG_a_LassoCV_b, mse_CG_a_RF_b,mse_CG_a_LassoCV_b_hm, mse_CG_a_RF_b_hm)
ratio_good_comp = compose_table(pct_good_CT_h, pct_good_CT_a, pct_good_CG_a_am, pct_good_CG_a, pct_good_CG_a_LassoCV, pct_good_CG_a_RF,
                                pct_good_CG_a_LassoCV_b, pct_good_CG_a_RF_b,pct_good_CG_a_LassoCV_b_hm, pct_good_CG_a_RF_b_hm)

if(T){ #Output tables
  fmt_table(xtable(n_cells_comp, digits=2), paste0(tbl_export_path, "n_cells",suffix,".tex"))
  fmt_table(xtable(mse_comp, digits=3), paste0(tbl_export_path, "mse",suffix,".tex"))
  fmt_table(xtable(ratio_good_comp), paste0(tbl_export_path, "ratio_good",suffix,".tex"))
}

# Output examples ------------
ct_sim_plot <- function(d, N_i, s, honest=TRUE) {
  iter = get_iter(d, N_i, s)
  list[y, X, w, tau] = data1[[iter]]
  X_range = get_X_range(X)
  if(honest)
    plt = plot_2D_partition.rpart(ct_h_fit_models[[iter]], X_range=X_range)
  else
    plt = plot_2D_partition.rpart(ct_a_fit_models[[iter]], X_range=X_range)
  return(plt + ggtitle("Causal Tree") + labs(fill = "tau(X)"))
}
cg_sim_plot <- function(d, N_i, s, lasso=FALSE) {
  iter = get_iter(d, N_i, s)
  list[y, X, w, tau] = data1[[iter]]
  N = Ns[N_i]
  X_range = get_X_range(X)
  if(lasso)
    grid_fit = cg_a_LassoCV_fit_models[[iter]]
  else
    grid_fit = cg_a_fit_models[[iter]]
  grid_a_m_fit <- change_complexity(grid_fit, y, X, d=w, which.min(abs(ct_h_nl[iter] - (grid_fit$complexity_seq + 1))))
  plt = plot_2D_partition.estimated_partition(grid_a_m_fit, c("X1", "X2"))
  return(plt + ggtitle("Causal Tree") + labs(fill = "tau(X)"))
}
ct_cg_plot <- function(d, N_i, s, ct_honest=TRUE, cg_lasso=FALSE) {
  ct_plt = ct_sim_plot(d, N_i, s, honest=ct_honest)
  cg_plt = cg_sim_plot(d, N_i, s, lasso=cg_lasso)
  grid.arrange(ct_plt + ggtitle("Causal Tree"), cg_plt + ggtitle("Causal Grid"), ncol=2)
}

if(F) {
  #Pick which one to show. 
  ct_a_nl = results_ct_a[,4]
  ct_h_nl = results_ct_h[,4]
  n_dim_splits = matrix(0,nrow=D*NN*S, ncol=6)
  for(d in 1:D) {
    k = Ks[d]
    for(N_i in 1:NN) {
      for(s in 1:S) {
        iter = get_iter(d, N_i, s)
        N = Ns[N_i]
        list[y, X, w, tau] = data1[[iter]]
        #list[y_te, X_te, w_te, tau_te] = data2[[iter]]
        tr_sample = c(rep(TRUE, N), rep(FALSE, N))
        
        grid_a_fit <- cg_a_fit_models[[iter]]
        grid_a_LassoCV_fit <- cg_a_LassoCV_fit_models[[iter]]
        grid_a_m_fit <- change_complexity(grid_a_fit, y, X, d=w, which.min(abs(ct_h_nl[iter] - (grid_a_fit$complexity_seq + 1))))
        grid_a_LassoCV_m_fit <- change_complexity(grid_a_LassoCV_fit, y, X, d=w, which.min(abs(ct_h_nl[iter] - (grid_a_LassoCV_fit$complexity_seq + 1))))
        
        n_dim_splits[iter, 1] = ct_ndim_splits(ct_h_fit_models[[iter]], k)
        n_dim_splits[iter, 2] = ct_ndim_splits(ct_a_fit_models[[iter]], k)
        n_dim_splits[iter, 3] = cg_ndim_splits(grid_a_fit)
        n_dim_splits[iter, 4] = cg_ndim_splits(grid_a_LassoCV_fit)
        n_dim_splits[iter, 5] = cg_ndim_splits(grid_a_m_fit)
        n_dim_splits[iter, 6] = cg_ndim_splits(grid_a_LassoCV_m_fit)
      }
    }
  }
  #const_mask = results_ct_h[,1]==2 & results_ct_h[,4]>=4 & n_dim_splits[,1]==2 & n_dim_splits[,5]==2 #cg: normal
  #grph_pick = cbind(results_ct_h[const_mask,c(1,2,3, 4)],results_cg_a_m[const_mask,c(4)])
  const_mask = results_ct_h[,1]==2 & results_ct_h[,4]>=4 & n_dim_splits[,1]==2 & n_dim_splits[,6]==2 #cg: Lasso
  grph_pick = cbind(results_ct_h[const_mask,c(1,2,3, 4)],results_cg_a_LassoCV_m[const_mask,c(4)])
  grph_pick
  ct_cg_plot(2, 1, 57, ct_honest=TRUE, cg_lasso=TRUE)
}

cg_table <- function(obj, digits=3){
  stats = obj$cell_stats$stats[c("param_est")]
  colnames(stats) <- c("Est.")
  tbl = cbind(get_desc_df(obj$partition, do_str=TRUE, drop_unsplit=TRUE, digits=digits), stats)
}

ggplot_to_pdf <-function(plt_obj, filename) {
  pdf(file=filename)
  print(plt_obj)
  dev.off()
}

if(F){
  d=2 #so we can hve a 2d model
  N_i=1
  s=57
  iter = get_iter(d, N_i, s)
  list[y, X, w, tau] = data1[[iter]]
  X_range = get_X_range(X)
  
  ct_m = ct_h_fit_models[[iter]]
  ct_m_desc = ct_desc(ct_m)
  cat(ct_m_desc, file=paste0(tbl_export_path, "ct_ex_2d",suffix,".tex"), sep="\n")
  ct_plt = plot_2D_partition.rpart(ct_m, X_range=X_range) + labs(fill = "tau(X)")
  print(ct_plt + ggtitle("Causal Tree"))
  
  grid_fit = cg_a_LassoCV_fit_models[[iter]] #cg_a_fit_models[[iter]]
  cg_m <- change_complexity(grid_fit, y, X, d=w, which.min(abs(ct_h_nl[iter] - (grid_fit$complexity_seq + 1))))
  print(cg_m)
  cg_m_tbl = cg_table(cg_m)
  temp <- capture.output(cg_tbl_file_cont <- print(xtable(cg_m_tbl, digits=3), floating=F, comment = F))
  cat(cg_tbl_file_cont, file=paste0(tbl_export_path, "cg_ex_2d",suffix,".tex"))
  cg_plt = plot_2D_partition.estimated_partition(cg_m, c("X1", "X2")) + labs(fill = "tau(X)")
  print(cg_plt + ggtitle("Causal Grid"))
  
  ggplot_to_pdf(ct_plt + ggtitle("Causal Tree"), paste0(fig_export_path, "ct_ex_2d",suffix,".pdf"))
  ggplot_to_pdf(cg_plt + ggtitle("Causal Grid"), paste0(fig_export_path, "cg_ex_2d",suffix,".pdf"))
  
  pdf(file=paste0(fig_export_path, "cg_ct_ex_2d",suffix,".pdf"), width=8, height=4)
  grid.arrange(ct_plt + ggtitle("Causal Tree"), cg_plt + ggtitle("Causal Grid"), ncol=2)
  dev.off()
}

