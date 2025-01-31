args <- commandArgs(trailingOnly = TRUE); for (arg in args) {eval(parse(text = arg))} #to retrieve input args

source("utilities/AGS_utilities.R")
load(paste0("cache/",modelname,"/routines",".rds"))
source(list.files(path = paste0("cache/",modelname,"/"),recursive = TRUE,pattern = "\\.cpp\\.R$",full.names = TRUE)) #retrieve and load the .cpp.R file

load(paste0("data/datain_",modelname,".rds"))

X_list <- lapply(1:nrow(datain$X), FUN = function(i) datain$X[i,,drop=FALSE]) #it is required by the cpp functions
Z_list <- lapply(1:nrow(datain$X), FUN = function(i) datain$Z[i,,drop=FALSE])
J <- NCOL(datain$X); H <- NCOL(datain$Z)

set.seed(seedx + chainId) 

parx_init <- runif(J+H,-0.5,0.5)
parx_init <- tryCatch(optim(par = parx_init,fn = function(x){-sum(log(fY(datain$m,datain$X,datain$Z,datain$lb,datain$ub,x)))},method = "BFGS")$par,error=function(e){parx_init})

#X = datain$X;Z = datain$Z; m = datain$m;s = datain$s;lb = datain$lb;ub = datain$ub; mu.theta = datain$mu_theta;sigma.2.theta = datain$sigma_theta;print.out=TRUE;sigma0_type=-1;sigma0delta=2.5
#out <- AGS_alg(B,J,H,datain$X,datain$Z,X_list,Z_list,datain$m,datain$s,datain$lb,datain$ub,parx_init,datain$mu_theta,datain$sigma_theta,print.out=TRUE,sigma0_type=-1,sigma0delta=2.5)

out <- NULL; kkMax <- 25; kk <- 0
while(is.null(out) && kk<=kkMax){
  print(kk)
  out <- tryCatch(AGS_alg(B,J,H,datain$X,datain$Z,X_list,Z_list,datain$m,datain$s,datain$lb,datain$ub,parx_init,datain$mu_theta,datain$sigma_theta,print.out=FALSE,sigma0_type=-1,sigma0delta=2.5),error=function(e){NULL})
  kk <- kk+1}

saveRDS(out,file = paste0("results/dataout_",paste0(modelname,"_",chainId,".rds")))
