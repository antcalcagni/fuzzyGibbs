
rm_source <- FALSE
dm_conversion <- 1 # turn-off the c++ conversion for skew-normal based approximation of pi(theta|y,..)


##### Converting db_internal() #####

# Converting fY (no append)
hds <- list(y="double",X="vector<double>",Z="vector<double>",lb="double",ub="double",pars="vector<double>")
R2Cpp_funConverter(funIN = fY,namefun = "fy",typefun="double",fun_heading=hds,rcpp_export = TRUE,fileout = fileout,fileout_append = FALSE,filehead = filehead,filetail = filetail)

# Converting d1lnfy (yes append)
R2Cpp_funConverter(funIN = d1lnFy,namefun = "d1lnfy",typefun="double",fun_heading=hds,rcpp_export = TRUE,fileout = fileout,fileout_append = TRUE,filehead = filehead,filetail = filetail)

# Converting d2lnfy (yes append)
R2Cpp_funConverter(cpp_finalize = TRUE,funIN = d2lnFy,namefun = "d2lnfy",typefun="double",fun_heading=hds,rcpp_export = TRUE,fileout = fileout,fileout_append = TRUE,filehead = filehead,filetail = filetail)

# Compile C++ file
# Rcpp::sourceCpp(file = fileout,verbose = FALSE) --compile the file just at the end of the conversion process!


##### Converting functions for dm_internal() #####

if(dm_conversion){

# Converting fY_internal (no append)
xb <- paste0("b",1:(J+H)); pars_list <- as.list(rep("double",J+H)); names(pars_list) <- xb
hds <- c(list(y="double",X="vector<double>",Z="vector<double>",lb="double",ub="double"),pars_list)
R2Cpp_funConverter(funIN = fY_internal,namefun = "fy_internal",typefun="double",fun_heading=hds,rcpp_export = TRUE,fileout = fileout,fileout_append = TRUE,filehead = filehead)

# Converting grad_lnFy_internal
for(j in 1:length(grad_lnFy_internal)){
  R2Cpp_funConverter(funIN = grad_lnFy_internal[[j]],namefun = paste0("grad_lnFy_internal",j),typefun="auto",fun_heading=hds,rcpp_export = TRUE,fileout = fileout,fileout_append = TRUE,filehead = filehead)
}
R2Cpp_grad_lnfy_exec(JH = J+H,fileout = fileout,fileout_append = TRUE)
R2Cpp_grad_lnfy_exec_loop(JH = J+H,fileout = fileout,fileout_append = TRUE)

# Converting hess_lnFy_internal
for(j in 1:length(hess_lnFy_internal)){
  R2Cpp_funConverter(funIN = hess_lnFy_internal[[j]],namefun = paste0("hess_lnfy_internal",j),typefun="auto",fun_heading=hds,rcpp_export = TRUE,fileout = fileout,fileout_append = TRUE,filehead = filehead)
}
R2Cpp_hess_lnfy_exec(K = K,JH = J+H,fileout = fileout,fileout_append = TRUE)
R2Cpp_hess_lnfy_exec_loop(JH = J+H,fileout = fileout,fileout_append = TRUE)

}

# Compile C++ file
Rcpp::sourceCpp(file = fileout,verbose = verboseCpp,cacheDir = cacheDir_rcpp,cleanupCacheDir = FALSE)

# Remove original files (now converted in C++)
if(rm_source && dm_conversion){rm(grad_lnFy_internal,grad_lnFy_exec,hess_lnFy_internal,hessian_lnFy_exec)}
if(rm_source){rm(lnFy,lnfyA,d1lnFy,d2lnFy,dpi_y_approx,conditional_sampling_y)}
#if(rm_source){file.remove(paste0(getwd(),"/",fileout))}
