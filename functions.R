
eb_estimation = function(bulk_data,panel,prior_mean,prior_sd,shrinkSD_factor,n_iter) { 
  
  log_dir = getwd()
  panel   = as.matrix(panel)
  
  n_samples       = nrow(bulk_data)
  n_sites         = ncol(bulk_data)
  n_cellTypes     = ncol(panel)
  
  props           = matrix(NA,n_samples,n_cellTypes)
  colnames(props) = colnames(panel)
  r2              = rep(0,n_samples)
  
  fmla = bulk ~ -1 + panel

  print( "Start estimation: Empirical bayes")
#  sink(paste0(log_dir,"/R.log"))  
  for(i in 1:n_samples) {  # i=1
    
    print( paste0("Processing sample ",i) )  
    bulk = as.vector( t(bulk_data[i,]) )
    data = data.frame(bulk, panel )        
    
    mod_stan = stan_glm(fmla, data = data,
                        chains = 1, iter = n_iter, 
                        prior = normal(location = prior_mean, scale=prior_sd/shrinkSD_factor),
                        refresh = 0)
    props[i,1:n_cellTypes] =  coef(mod_stan)    #  posterior median estimates      
    r2[i] = median( bayes_R2(  mod_stan  ) )
    
  }
#  closeAllConnections() 
  
  props[is.na(props)] = 0
  nprops              = data.frame(props)
  nprops[nprops < 0]  = 0
  estimates           = nprops/rowSums(nprops) # standardize
  rownames(estimates) = rownames(bulk_data)
  colnames(estimates) = colnames(props)
  
  colnames(nprops)    = paste0(colnames(props),"raw")
  nprops$nprops_max1  = apply(nprops[,colnames(nprops)], 1, max) 
  nprops$nprops_min0  = apply(nprops[,colnames(nprops)], 1, min) 
  nprops$nprops_max1  = round(nprops$nprops_max1 == 1) 
  nprops$nprops_min0  = round(nprops$nprops_min0 == 0) 
  
  results             = list()
  results$estimates   = estimates
  results$qc          = data.frame(nprops,r2,shrinkSD_factor) 
  rownames(results$qc) = rownames(bulk_data)
  results
  
  
  results
  
} 

constraint_estimation = function(bulk_data,panel) { 
  
  bulk_data  = as.matrix(bulk_data)
  panel      = as.matrix(panel)
  n_samples  = nrow(bulk_data)
  n_celltype = ncol(panel)
  
  cell_props           = matrix(NA,n_samples,n_celltype)
  colnames(cell_props) = colnames(panel)
  r2                   = rep(0,n_samples)
  
  print( "Start estimation: Use boundary contraint >= 0")
  bol_positive    = rep(T,n_celltype)
  for(i in 1:n_samples) { # i =1 
    print( paste0("Processing sample ",i) )  
    mod = penalized(bulk_data[i,], ~ panel, ~-1,lambda1=0, lambda2=0, positive=bol_positive,trace=F)
    r2[i] = cor( fitted.values(mod), bulk_data[i,],use="pairwise.complete.obs")^2
    cell_props[i,1:n_celltype] = coef(mod)[1:n_celltype]           
  }  
  
  cell_props[is.na(cell_props)] = 0
  ncell_props              = data.frame(cell_props)
  ncell_props[ncell_props < 0]  = 0
  estimates           = ncell_props/rowSums(ncell_props) # standardize
  rownames(estimates) = rownames(bulk_data)
  colnames(estimates) = colnames(cell_props)
  
  colnames(ncell_props)    = paste0(colnames(cell_props),"raw")
  ncell_props$ncell_props_max1  = apply(ncell_props[,colnames(ncell_props)], 1, max) 
  ncell_props$ncell_props_min0  = apply(ncell_props[,colnames(ncell_props)], 1, min) 
  ncell_props$ncell_props_max1  = round(ncell_props$ncell_props_max1 == 1) 
  ncell_props$ncell_props_min0  = round(ncell_props$ncell_props_min0 == 0) 
  
  results             = list()
  results$estimates   = estimates
  results$qc          = data.frame(ncell_props,r2) 
  rownames(results$qc) =  rownames(bulk_data)
  results
  
} 



sim_bulk_data = function(n_subjects,cellprop_stats,panel,error_ratio) {
  
  panel          = as.matrix(panel)
  cellprop_stats = as.matrix(cellprop_stats)
  
  n_celltype     = ncol(panel) 
  cell_props     = matrix(NA,n_subjects,n_celltype)
  rownames(cell_props) = paste0('Sample_',1:n_subjects)
  colnames(cell_props) = colnames(panel) 
  for (i in 1:n_celltype) cell_props[,i]  = rgbeta(n_subjects, cellprop_stats[1,i], cellprop_stats[2,i], cellprop_stats[3,i], cellprop_stats[4,i])
  cell_props     = cell_props / apply(cell_props,1,sum)
  
  bulk_data = cell_props %*% t(panel)
  
  if ( error_ratio > 0 ) {
    gene_var  = apply(bulk_data,2,var) 
    error_var = error_ratio * gene_var 
    bulk_data = bulk_data + rnorm(nrow(bulk_data)*ncol(bulk_data), mean=0, sd=sqrt(error_var))
  }

  list(bulk_data,cell_props) 
}

rgbeta = function(n, mean, var, min = 0, max = 1) {
  
  dmin = mean - min
  dmax = max - mean
  
  if (dmin <= 0 || dmax <= 0) stop(paste("mean must be between min =", min, "and max =", max)) 
  if (var >= dmin * dmax) stop(paste("var must be less than (mean - min) * (max - mean) =", dmin * dmax))
  
  # mean and variance of the standard beta distributed variable
  mx = (mean - min) / (max - min)
  vx = var / (max - min)^2
  
  # find the corresponding alpha-beta parameterization
  a = ((1 - mx) / vx - 1 / mx) * mx^2
  b = a * (1 / mx - 1)
  
  # generate standard beta observations and transform
  x = rbeta(n, a, b)
  y = (max - min) * x + min
  
  y
}

