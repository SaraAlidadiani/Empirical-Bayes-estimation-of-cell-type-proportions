rm(list=ls())
work_dir       = "G:\\My Drive\\SOP\\snRNAseq_brain\\paper\\panel\\github\\"

# specify info needed to simulate bulk data
panel           = read.csv(paste0(work_dir,'panel_scaled.csv'),stringsAsFactors=F,row.names=1)
cellprop_stats  = read.csv(paste0(work_dir,'cell_prop_stats.csv'),row.names=1) # mean, variance, min, max of cell type proportions
n_subjects      = 10000   # number of bulk samples to simulate
error_ratio     = 0.5  # error_ratio * VAR(gene) is extra error variance in bulk data when bulk data cannot be perfectly predicted from the the panel means and the subjects cell type proportions 
shrinkSD_factor = 0.5  # SD / shrink_factor give SD of the prior
n_iter          = 5000 # Number of iterations in Empirical Bayes estimation
# begin
source(paste0(work_dir,"functions.R"))
library(penalized)
library(rstanarm) 
library(openxlsx)


results_list = list()
temp       = sim_bulk_data(n_subjects,cellprop_stats,panel,error_ratio)
bulk_data  = temp[[1]]
cell_props = temp[[2]]

temp       = constraint_estimation(bulk_data,panel) 
estimates  = temp[[1]]
qc         = temp[[2]]
results_list[[length(results_list)+1]] = estimates 
results_list[[length(results_list)+1]] = qc 

prior_mean = apply( estimates, 2, mean, na.rm=T)
prior_sd   = apply( estimates, 2, sd, na.rm=T)
temp       = eb_estimation(bulk_data,panel,prior_mean,prior_sd,shrinkSD_factor,n_iter)
estimates  = temp[[1]]
qc         = temp[[2]]
results_list[[length(results_list)+1]] = estimates 
results_list[[length(results_list)+1]] = qc 

results     = matrix(NA,ncol(estimates),5)
colnames(results) = c('Sim. mean','Est. mean','Est. Bias','Est. RMSE','Est. zero')
rownames(results) = colnames(estimates)
results[,1]  = apply(cell_props, 2, mean, na.rm=T) 
results[,2]  = apply(estimates , 2, mean, na.rm=T) 
results[,3]  = results[,1] - results[,2]
results[,4]  = sqrt( apply((estimates-cell_props )^2, 2, mean, na.rm=T) )
results[,5]  = apply(estimates==0, 2, sum, na.rm=T) 

Overall_mean = apply(results, 2, mean, na.rm=T) 
results_list[[length(results_list)+1]] = rbind(results,Overall_mean)

names(results_list) = c('estimates_constraint','qc_constraint','estimates_eb','qc_eb','simulation_results')
write.xlsx(results_list,paste0(work_dir,'estimation_results.xlsx'),rowNames=T)




