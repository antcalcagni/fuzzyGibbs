
# External validation: Dataset 1 ------------------------------------------
## SOURCE - Colubi, A.: Statistical inference about the means of fuzzy random variables, Fuzzy Sets and Systems, 160(3), pp. 344-356 (2009)

rm(list=ls());graphics.off(); source("utilities/utilities.R")
J <- 3; H <- 1 #number of predictors for mu and phi

##### Model 1 [fy = Beta] ####
modelName <- "dataset1_beta"

### (Step 1) Model building
fY <- function(y,X,Z,lb,ub,pars){
  b <- pars[1:ncol(X)]
  g <- pars[(1+ncol(X)):(ncol(X)+ncol(Z))]
  eta <- X%*%b
  theta <- exp(eta)/(1+exp(eta))
  phi <- exp(Z%*%g)
  alpha <- theta*phi
  beta <- phi-theta*phi
  ((y - lb)^(alpha - 1)) * ((ub - y)^(beta - 1)) / ((ub - lb)^(alpha + beta - 1) * (gamma(alpha) * gamma(beta) / gamma(alpha + beta)))
}
db_internal() #create internal R functions

# Preparing for c++ conversion of R functions
filehead <- "utilities/R2Cpp_internal.head"; filetail <- "utilities/R2Cpp_internal.tail"; fileout <- "utilities/R2Cpp_internal.cpp"
verboseCpp <- TRUE
source("utilities/R2Cpp_utilities.R")
cacheDir_rcpp <-  paste0(getwd(),"/cache/",modelName)

dm_internal(fY,J,H) #create internal R functions
source("utilities/R2Cpp_funConverter_exec.R") #convert to c++

# Saving compiled files
save(fY,fY_internal,tud_lnFy_exec,tud_lnFy_internal,grad_fun,hess_fun,file = paste0("cache/",modelName,"/routines.rds"))
rm(list=ls())


### (Step 2) Preparing for model fitting
modelName <- "dataset1_beta"

# Load data
datax <- readRDS("data/colubi_2009.rds") #load data
str(datax)

datax <- datax[sample(1:nrow(datax),200,replace = FALSE),]

n <- NROW(datax)
datax$tree <- as.factor(datax$tree)
X <- model.matrix(~tree,data = datax)
Z <- matrix(1,nrow = n,ncol = 1)
m <- datax$m01; s <- datax$s;

lb <- 0; ub <- 1
J <- NCOL(X); H <- 1

# Defining priors
mu_theta <- c(rep(0,J),rep(-0.5,H)); sigma_theta <- c(rep(3.5,J),rep(13.5,H))

# Saving data
datain <- list(X=X,Z=Z,m=m,s=s,lb=lb,ub=ub,mu_theta=mu_theta,sigma_theta=sigma_theta) #input for the AGS algorithm
save(datain,file = paste0("data/datain_",modelName,".rds")) #saving inputs

### (Step 3) Running multiple chains via GNU Parallel
nchains <- 5 
B <- 4e3 #samples

cmds <- sapply(1:nchains,function(u){
  cmd <- paste0(
    "R CMD BATCH --no-restore --no-save '--args B=", B, 
    " chainId=", u,
    " seedx=", 240125,
    " modelname=\"", modelName, "\"", 
    "' AGS_exec.R ", 
    "cache/",modelName,"/out_chainId_",u,".stdout"
  )
})

ncores <- nchains #ncores=nchains
btcs <- split(cmds, ceiling(seq_along(cmds)/ncores)) 
for(jobs in btcs){ 
  cat("\n Running:");for(x in jobs){cat("\n\t",strsplit(x,split = "--args")[[1]][2])}
  
  batch_file <- tempfile("job_list_"); writeLines(jobs, batch_file)
  system(paste0("parallel -j ", nchains, " < ", batch_file)); unlink(batch_file)
}

### (Step 4) Save final mcmc as array (B x Nchains x Params)
propBurIn <- 0.5
load(paste0("cache/",modelName,"/routines.rds"))

# MCMCs
burnin <- floor(B*propBurIn)
out <- lapply(list.files(path = "results/",pattern = paste0("dataout_",modelName),full.names = TRUE),readRDS)
out <- lapply(out,function(x)coda::as.mcmc(x[(burnin+1):B,]))
out_list <- coda::as.mcmc.list(out)
out_array <- posterior::as_draws_array(out_list)

# Relevant MCMCs statistics
out <- matrix(out_array,nrow = dim(out_array)[1]*dim(out_array)[2],ncol = dim(out_array)[3])

Log_lik <- matrix(NA,nrow = nrow(out),ncol = nrow(datain$X))
for(b in 1:nrow(out)){
  muy <- plogis(datain$X%*%out[b,1:J]); phiy <- exp(datain$Z%*%out[b,J+H])
  ystar <- mapply(function(i)extraDistr::rnsbeta(n = 1,shape1 = muy[i]*phiy[i],shape2 = phiy[i]-phiy[i]*muy[i],min = datain$lb,max = datain$ub),1:nrow(datain$X))
  Log_lik[b,] <- log(fY(y = ystar,X = datain$X,Z = datain$Z,lb = datain$lb,ub = datain$ub,pars = out[b,]))
}

lppd <- sum(log(rowMeans(exp(Log_lik))))  # Log pointwise predictive density
p_waic <- sum(apply(Log_lik, 2, var))     # Pointwise variance of log-likelihood
waic <- -2 * (lppd - p_waic)              # WAIC measure
pointwise_waic <- -2 * (log(rowMeans(exp(Log_lik))) - apply(Log_lik, 1, var))
waic_se <- sqrt(length(pointwise_waic) * var(pointwise_waic)) # WAIC Standard Error

deviance_samples <- -2 * rowSums(Log_lik)  # Deviance for each posterior sample
D_bar <- mean(deviance_samples)
mean_log_lik <- colMeans(exp(Log_lik))  # Average likelihood across samples
D_theta_hat <- -2 * sum(log(mean_log_lik))     # Deviance at posterior mean
dic <- 2 * D_bar - D_theta_hat

posterior_table <- summarise_draws(posterior::as_draws_array(out))

# Saving all together
saveRDS(list(mcmc=out_array,stats=list(ptable=posterior_table,waic=waic,waic_se=waic_se,dic=dic)),file = paste0("results/mcmc_",modelName,".rds"))


rm(list=ls());graphics.off(); source("utilities/utilities.R")
J <- 3; H <- 1 #number of predictors for mu and phi


##### Model 2 [fy = LogitNorm] ####
modelName <- "dataset1_logitnorm"

### (Step 1) Model building
fY <- function(y,X,Z,lb,ub,pars){
  b <- pars[1:ncol(X)]
  g <- pars[(1+ncol(X)):(ncol(X)+ncol(Z))]
  eta <- X%*%b
  theta <- eta
  phi <- exp(Z%*%g)
  (1/(phi*sqrt(2*pi)))*(1/(y*(1-y)))*exp(-(log(y/(1-y))-theta)^2/(2*phi^2))
}
db_internal() #create internal R functions

# Preparing for c++ conversion of R functions
filehead <- "utilities/R2Cpp_internal.head"; filetail <- "utilities/R2Cpp_internal.tail"; fileout <- "utilities/R2Cpp_internal.cpp"
verboseCpp <- TRUE
source("utilities/R2Cpp_utilities.R")
cacheDir_rcpp <-  paste0(getwd(),"/cache/",modelName)

dm_internal(fY,J,H) #create internal R functions
source("utilities/R2Cpp_funConverter_exec.R") #convert to c++

# Saving compiled files
save(fY,fY_internal,tud_lnFy_exec,tud_lnFy_internal,grad_fun,hess_fun,file = paste0("cache/",modelName,"/routines.rds"))
rm(list=ls())

### (Step 2) Preparing for model fitting
modelName <- "dataset1_logitnorm"

# Load data
datax <- readRDS("data/colubi_2009.rds") #load data
str(datax)

datax <- datax[sample(1:nrow(datax),200,replace = FALSE),]

n <- NROW(datax)
datax$tree <- as.factor(datax$tree)
X <- model.matrix(~tree,data = datax)
Z <- matrix(1,nrow = n,ncol = 1)
m <- datax$m01; s <- datax$s;

lb <- 0; ub <- 1
J <- NCOL(X); H <- 1

# Defining priors
mu_theta <- c(rep(0,NCOL(X)),rep(-0.5,H)); sigma_theta <- c(rep(1,NCOL(X)),rep(0.35,H))

# Saving data
datain <- list(X=X,Z=Z,m=m,s=s,lb=lb,ub=ub,mu_theta=mu_theta,sigma_theta=sigma_theta) #input for the AGS algorithm
save(datain,file = paste0("data/datain_",modelName,".rds")) #saving inputs

### (Step 3) Running multiple chains via GNU Parallel
nchains <- 5 
B <- 4e3 #samples

cmds <- sapply(1:nchains,function(u){
  cmd <- paste0(
    "R CMD BATCH --no-restore --no-save '--args B=", B, 
    " chainId=", u,
    " seedx=", 240125,
    " modelname=\"", modelName, "\"", 
    "' AGS_exec.R ", 
    "cache/",modelName,"/out_chainId_",u,".stdout"
  )
})

ncores <- nchains #ncores=nchains
btcs <- split(cmds, ceiling(seq_along(cmds)/ncores)) 
for(jobs in btcs){ 
  cat("\n Running:");for(x in jobs){cat("\n\t",strsplit(x,split = "--args")[[1]][2])}
  
  batch_file <- tempfile("job_list_"); writeLines(jobs, batch_file)
  system(paste0("parallel -j ", nchains, " < ", batch_file)); unlink(batch_file)
}


### (Step 4) Save final mcmc as array (B x Nchains x Params)

propBurIn <- 0.5
load(paste0("cache/",modelName,"/routines.rds"))

# MCMCs
burnin <- floor(B*propBurIn)
out <- lapply(list.files(path = "results/",pattern = paste0("dataout_",modelName),full.names = TRUE),readRDS)
out <- lapply(out,function(x)coda::as.mcmc(x[(burnin+1):B,]))
out_list <- coda::as.mcmc.list(out)
out_array <- posterior::as_draws_array(out_list)

# Relevant MCMCs statistics
out <- matrix(out_array,nrow = dim(out_array)[1]*dim(out_array)[2],ncol = dim(out_array)[3])

Log_lik <- matrix(NA,nrow = nrow(out),ncol = nrow(datain$X))
for(b in 1:nrow(out)){
  muy <- datain$X%*%out[b,1:J]; phiy <- exp(datain$Z%*%out[b,J+H])
  ystar <- mapply(function(i)greybox::rlogitnorm(1,muy[i],phiy[i]),1:nrow(datain$X))
  Log_lik[b,] <- log(fY(y = ystar,X = datain$X,Z = datain$Z,lb = datain$lb,ub = datain$ub,pars = out[b,]))
}

lppd <- sum(log(rowMeans(exp(Log_lik))))  # Log pointwise predictive density
p_waic <- sum(apply(Log_lik, 2, var))     # Pointwise variance of log-likelihood
waic <- -2 * (lppd - p_waic)              # WAIC measure
pointwise_waic <- -2 * (log(rowMeans(exp(Log_lik))) - apply(Log_lik, 1, var))
waic_se <- sqrt(length(pointwise_waic) * var(pointwise_waic)) # WAIC Standard Error

deviance_samples <- -2 * rowSums(Log_lik)  # Deviance for each posterior sample
D_bar <- mean(deviance_samples)
mean_log_lik <- colMeans(exp(Log_lik))  # Average likelihood across samples
D_theta_hat <- -2 * sum(log(mean_log_lik))     # Deviance at posterior mean
dic <- 2 * D_bar - D_theta_hat

posterior_table <- summarise_draws(posterior::as_draws_array(out))

# Saving all together
saveRDS(list(mcmc=out_array,stats=list(ptable=posterior_table,waic=waic,waic_se=waic_se,dic=dic)),file = paste0("results/mcmc_",modelName,".rds"))



# External validation: Dataset 2 ------------------------------------------
## SOURCE - Kim, B., & Bishu, R. R. (1998). Evaluation of fuzzy linear regression models by comparing membership functions. Fuzzy sets and systems, 100(1-3), 343-352

rm(list=ls());graphics.off(); source("utilities/utilities.R")
J <- 4; H <- 1 #number of predictors for mu and phi

##### Model 1 [fy = Lognormal] ####
modelName <- "dataset2_lognormal"

### (Step 1) Model building
fY <- function(y,X,Z,lb,ub,pars){
  b <- pars[1:ncol(X)]
  g <- pars[(1+ncol(X)):(ncol(X)+ncol(Z))]
  eta <- X%*%b
  theta <- eta
  phi <- exp(Z%*%g)
  (1 / (y*phi * sqrt(2 * pi))) * exp(-((log(y) - theta)^2) / (2 * phi^2))
}
db_internal() #create internal R functions

# Preparing for c++ conversion of R functions
filehead <- "utilities/R2Cpp_internal.head"; filetail <- "utilities/R2Cpp_internal.tail"; fileout <- "utilities/R2Cpp_internal.cpp"
verboseCpp <- TRUE
source("utilities/R2Cpp_utilities.R")
cacheDir_rcpp <-  paste0(getwd(),"/cache/",modelName)

dm_internal(fY,J,H) #create internal R functions
source("utilities/R2Cpp_funConverter_exec.R") #convert to c++

# Saving compiled files
save(fY,fY_internal,tud_lnFy_exec,tud_lnFy_internal,grad_fun,hess_fun,file = paste0("cache/",modelName,"/routines.rds"))
rm(list=ls())


### (Step 2) Preparing for model fitting
modelName <- "dataset2_lognormal"

# Load data
datax <- readRDS("data/kim_bishu_1998.rds") #load data
str(datax)

n <- NROW(datax)
X <- model.matrix(~in_CRexp+out_CRexp+education,data = datax)
X[,2:4] <- scale(X[,2:4])
Z <- matrix(1,nrow = n,ncol = 1)
m <- datax$RT_m_beta; s <- datax$RT_s;

# Approximate bounds as Y is on R in this case!
source("utils.R")
rescale <- function(x,lb=0,ub=1){(x-lb)/(ub-lb)}   #linearly rescale into [lb,ub]
rescale_inv <- function(x,lb=0,ub=1){lb+(ub-lb)*x} #inverse of rescale

alphax <- 1e-9; xsup <- seq(0,1,length.out = 1e3)
Yobs_a0 <- t(mapply(function(i){
  mi_current <- rescale(m[i],min(m),max(m))
  fy <- beta_fn(xsup,mi_current,s[i])
  xout <- c(min(xsup[fy>alphax]),max(xsup[fy>alphax]))
  return(xout)
},1:n))
Yobs_a0 <- apply(Yobs_a0,2,function(x)rescale_inv(x,min(m),max(m)))
lb <- min(Yobs_a0[,1]); ub <- max(Yobs_a0[,2]) #proxies for the true bnds based on 0-alpha cuts
lb <- 0.01; ub <- ub+0.1

J <- NCOL(X); H <- 1

# Defining priors
mu_theta <- c(rep(0,J),rep(0,H)); sigma_theta <- c(rep(3.5,J),rep(0.25,H))

# Saving data
datain <- list(X=X,Z=Z,m=m,s=s,lb=lb,ub=ub,mu_theta=mu_theta,sigma_theta=sigma_theta) #input for the AGS algorithm
save(datain,file = paste0("data/datain_",modelName,".rds")) #saving inputs

### (Step 3) Running multiple chains via GNU Parallel
nchains <- 5 
B <- 4e3 #samples

cmds <- sapply(1:nchains,function(u){
  cmd <- paste0(
    "R CMD BATCH --no-restore --no-save '--args B=", B, 
    " chainId=", u,
    " seedx=", 240125,
    " modelname=\"", modelName, "\"", 
    "' AGS_exec.R ", 
    "cache/",modelName,"/out_chainId_",u,".stdout"
  )
})

ncores <- nchains #ncores=nchains
btcs <- split(cmds, ceiling(seq_along(cmds)/ncores)) 
for(jobs in btcs){ 
  cat("\n Running:");for(x in jobs){cat("\n\t",strsplit(x,split = "--args")[[1]][2])}
  
  batch_file <- tempfile("job_list_"); writeLines(jobs, batch_file)
  system(paste0("parallel -j ", nchains, " < ", batch_file)); unlink(batch_file)
}

### (Step 4) Save final mcmc as array (B x Nchains x Params)
propBurIn <- 0.5
load(paste0("cache/",modelName,"/routines.rds"))

# MCMCs
burnin <- floor(B*propBurIn)
out <- lapply(list.files(path = "results/",pattern = paste0("dataout_",modelName),full.names = TRUE),readRDS)
out <- lapply(out,function(x)coda::as.mcmc(x[(burnin+1):B,]))
out_list <- coda::as.mcmc.list(out)
out_array <- posterior::as_draws_array(out_list)

# Relevant MCMCs statistics
# out <- matrix(out_array,nrow = dim(out_array)[1]*dim(out_array)[2],ncol = dim(out_array)[3])

# Log_lik <- matrix(NA,nrow = nrow(out),ncol = nrow(datain$X))
# for(b in 1:nrow(out)){
#   muy <- plogis(datain$X%*%out[b,1:J]); phiy <- plogis(datain$Z%*%out[b,J+H]); qhi <- (log(1-0.5^phiy)/log(muy)) 
#   ystar <- mapply(function(i)extraDistr::rkumar(1,qhi[i],1/phi[i]),1:nrow(datain$X))  
#   Log_lik[b,] <- log(fY(y = ystar,X = datain$X,Z = datain$Z,lb = datain$lb,ub = datain$ub,pars = out[b,]))
# }

# lppd <- sum(log(rowMeans(exp(Log_lik))))  # Log pointwise predictive density
# p_waic <- sum(apply(Log_lik, 2, var))     # Pointwise variance of log-likelihood
# waic <- -2 * (lppd - p_waic)              # WAIC measure
# pointwise_waic <- -2 * (log(rowMeans(exp(Log_lik))) - apply(Log_lik, 1, var))
# waic_se <- sqrt(length(pointwise_waic) * var(pointwise_waic)) # WAIC Standard Error

# deviance_samples <- -2 * rowSums(Log_lik)  # Deviance for each posterior sample
# D_bar <- mean(deviance_samples)
# mean_log_lik <- colMeans(exp(Log_lik))  # Average likelihood across samples
# D_theta_hat <- -2 * sum(log(mean_log_lik))     # Deviance at posterior mean
# dic <- 2 * D_bar - D_theta_hat

posterior_table <- posterior::summarise_draws(posterior::as_draws_array(out))

# Saving all together
saveRDS(list(mcmc=out_array,stats=list(ptable=posterior_table,waic=waic,waic_se=waic_se,dic=dic)),file = paste0("results/mcmc_",modelName,".rds"))


# External validation: Dataset 3 ------------------------------------------
## SOURCE - ## Jiang, H., Kwong, C. K., Chan, C. Y., & Yung, K. L. (2019). A multi-objective evolutionary approach for fuzzy regression analysis. Expert Systems with Applications, 130, 225-235

rm(list=ls());graphics.off(); source("utilities/utilities.R")
J <- 2; H <- 1 #number of predictors for mu and phi

##### Model 1 [fy = Kumaraswamy] ####
modelName <- "dataset3_kumar"

### (Step 1) Model building
fY <- function(y,X,Z,lb,ub,pars){
  b <- pars[1:ncol(X)]
  g <- pars[(1+ncol(X)):(ncol(X)+ncol(Z))]
  eta <- X%*%b
  theta <- exp(eta)/(1+exp(eta))
  phi <- (Z%*%g)
  qhi <- exp(phi)/(1+exp(phi))
  alpha <- (log(1-0.5^qhi)/log(theta))
  alpha*(1/qhi)*y^(alpha - 1)*(1-y^alpha)^((1/qhi)-1)
}
db_internal() #create internal R functions

# Preparing for c++ conversion of R functions
filehead <- "utilities/R2Cpp_internal.head"; filetail <- "utilities/R2Cpp_internal.tail"; fileout <- "utilities/R2Cpp_internal.cpp"
verboseCpp <- TRUE
source("utilities/R2Cpp_utilities.R")
cacheDir_rcpp <-  paste0(getwd(),"/cache/",modelName)

dm_internal(fY,J,H) #create internal R functions
source("utilities/R2Cpp_funConverter_exec.R") #convert to c++

# Saving compiled files
save(fY,fY_internal,tud_lnFy_exec,tud_lnFy_internal,grad_fun,hess_fun,file = paste0("cache/",modelName,"/routines.rds"))
rm(list=ls())


### (Step 2) Preparing for model fitting
modelName <- "dataset3_kumar"

# Load data
datax <- readRDS("data/gong_etal_2018.rds") #load data
str(datax)

n <- NROW(datax)
X <- model.matrix(~X1,data = datax); X[,2] <- scale(X[,2])
Z <- matrix(1,nrow = n,ncol = 1)
m <- datax$ym01; s <- datax$ys;

lb <- 0.01; ub <- 0.99
J <- NCOL(X); H <- 1

# Defining priors
mu_theta <- c(rep(0,J),rep(0,H)); sigma_theta <- c(rep(3.5,J),rep(3.5,H))

# Saving data
datain <- list(X=X,Z=Z,m=m,s=s,lb=lb,ub=ub,mu_theta=mu_theta,sigma_theta=sigma_theta) #input for the AGS algorithm
save(datain,file = paste0("data/datain_",modelName,".rds")) #saving inputs

### (Step 3) Running multiple chains via GNU Parallel
nchains <- 2
B <- 1e3 #samples

seedx <- as.numeric(format(Sys.time(), format = "%d%M%S"))
cmds <- sapply(1:nchains,function(u){
  cmd <- paste0(
    "R CMD BATCH --no-restore --no-save '--args B=", B, 
    " chainId=", u,
    " seedx=", seedx,
    " modelname=\"", modelName, "\"", 
    "' AGS_exec.R ", 
    "cache/",modelName,"/out_chainId_",u,".stdout"
  )
})

ncores <- nchains #ncores=nchains
btcs <- split(cmds, ceiling(seq_along(cmds)/ncores)) 
for(jobs in btcs){ 
  cat("\n Running:");for(x in jobs){cat("\n\t",strsplit(x,split = "--args")[[1]][2])}
  
  batch_file <- tempfile("job_list_"); writeLines(jobs, batch_file)
  system(paste0("parallel -j ", nchains, " < ", batch_file)); unlink(batch_file)
}

### (Step 4) Save final mcmc as array (B x Nchains x Params)
propBurIn <- 0.5
load(paste0("cache/",modelName,"/routines.rds"))

# MCMCs
burnin <- floor(B*propBurIn)
out <- lapply(list.files(path = "results/",pattern = paste0("dataout_",modelName),full.names = TRUE),readRDS)
out <- lapply(out,function(x)coda::as.mcmc(x[(burnin+1):B,]))
out_list <- coda::as.mcmc.list(out)
out_array <- posterior::as_draws_array(out_list)

# Relevant MCMCs statistics
out <- matrix(out_array,nrow = dim(out_array)[1]*dim(out_array)[2],ncol = dim(out_array)[3])

Log_lik <- matrix(NA,nrow = nrow(out),ncol = nrow(datain$X))
for(b in 1:nrow(out)){
  muy <- plogis(datain$X%*%out[b,1:J]); phiy <- plogis(datain$Z%*%out[b,J+H]); qhi <- (log(1-0.5^phiy)/log(muy)) 
  ystar <- mapply(function(i)extraDistr::rkumar(1,qhi[i],1/phi[i]),1:nrow(datain$X))  
  Log_lik[b,] <- log(fY(y = ystar,X = datain$X,Z = datain$Z,lb = datain$lb,ub = datain$ub,pars = out[b,]))
}


lppd <- sum(log(rowMeans(exp(Log_lik))))  # Log pointwise predictive density
p_waic <- sum(apply(Log_lik, 2, var))     # Pointwise variance of log-likelihood
waic <- -2 * (lppd - p_waic)              # WAIC measure
pointwise_waic <- -2 * (log(rowMeans(exp(Log_lik))) - apply(Log_lik, 1, var))
waic_se <- sqrt(length(pointwise_waic) * var(pointwise_waic)) # WAIC Standard Error

deviance_samples <- -2 * rowSums(Log_lik)  # Deviance for each posterior sample
D_bar <- mean(deviance_samples)
mean_log_lik <- colMeans(exp(Log_lik))  # Average likelihood across samples
D_theta_hat <- -2 * sum(log(mean_log_lik))     # Deviance at posterior mean
dic <- 2 * D_bar - D_theta_hat

posterior_table <- summarise_draws(posterior::as_draws_array(out))

# Saving all together
saveRDS(list(mcmc=out_array,stats=list(ptable=posterior_table,waic=waic,waic_se=waic_se,dic=dic)),file = paste0("results/mcmc_",modelName,".rds"))



# Case study --------------------------------------------------------------
## SOURCE - Calcagni et al., A Bayesian beta linear model to analyze fuzzy rating responses. Book of short papers of the SIS Conference 2022

rm(list=ls());graphics.off(); source("utilities/utilities.R")
J <- 6; H <- 1 #number of predictors for mu and phi

##### Model 1 [fy = Beta] ####
modelName <- "casestudy_beta1"

### (Step 1) Model building
fY <- function(y,X,Z,lb,ub,pars){
  b <- pars[1:ncol(X)]
  g <- pars[(1+ncol(X)):(ncol(X)+ncol(Z))]
  eta <- X%*%b
  theta <- exp(eta)/(1+exp(eta))
  phi <- exp(Z%*%g)
  alpha <- theta*phi
  beta <- phi-theta*phi
  ((y - lb)^(alpha - 1)) * ((ub - y)^(beta - 1)) / ((ub - lb)^(alpha + beta - 1) * (gamma(alpha) * gamma(beta) / gamma(alpha + beta)))
}
db_internal() #create internal R functions

# Preparing for c++ conversion of R functions
filehead <- "utilities/R2Cpp_internal.head"; filetail <- "utilities/R2Cpp_internal.tail"; fileout <- "utilities/R2Cpp_internal.cpp"
verboseCpp <- TRUE
source("utilities/R2Cpp_utilities.R")
cacheDir_rcpp <-  paste0(getwd(),"/cache/",modelName)

dm_internal(fY,J,H) #create internal 1R functions
source("utilities/R2Cpp_funConverter_exec.R") #convert to c++

# Saving compiled files
save(fY,fY_internal,tud_lnFy_exec,tud_lnFy_internal,grad_fun,hess_fun,file = paste0("cache/",modelName,"/routines.rds"))
rm(list=ls())


### (Step 2) Preparing for model fitting
modelName <- "casestudy_beta1"

# Load data
load("data/VanLankveld_etal_2021.Rdata") #load data

n <- NROW(datax)
X = model.matrix(~ Age + Rel_length + desire + partner_respo + Gender_of_partner,data = datax); X[,2:NCOL(X)] <- scale(X[,2:NCOL(X)])
Z <- matrix(1,nrow = n,ncol = 1)
m <- datax$intimacy_m; s <- datax$intimacy_s;

lb <- 1; ub <- 5
J <- NCOL(X); H <- 1

# Defining priors
mu_theta <- c(rep(0,J),rep(0,H)); sigma_theta <- c(rep(3.5,J),rep(1.25,H))

# Saving data
datain <- list(X=X,Z=Z,m=m,s=s,lb=lb,ub=ub,mu_theta=mu_theta,sigma_theta=sigma_theta) #input for the AGS algorithm
save(datain,file = paste0("data/datain_",modelName,".rds")) #saving inputs

### (Step 3) Running multiple chains via GNU Parallel
nchains <- 5
B <- 4e3 #samples

seedx <- as.numeric(format(Sys.time(), format = "%d%M%S"))
cmds <- sapply(1:nchains,function(u){
  cmd <- paste0(
    "R CMD BATCH --no-restore --no-save '--args B=", B, 
    " chainId=", u,
    " seedx=", seedx,
    " modelname=\"", modelName, "\"", 
    "' AGS_exec.R ", 
    "cache/",modelName,"/out_chainId_",u,".stdout"
  )
})

ncores <- nchains #ncores=nchains
btcs <- split(cmds, ceiling(seq_along(cmds)/ncores)) 
for(jobs in btcs){ 
  cat("\n Running:");for(x in jobs){cat("\n\t",strsplit(x,split = "--args")[[1]][2])}
  
  batch_file <- tempfile("job_list_"); writeLines(jobs, batch_file)
  system(paste0("parallel -j ", nchains, " < ", batch_file)); unlink(batch_file)
}

### (Step 4) Save final mcmc as array (B x Nchains x Params)
propBurIn <- 0.5
load(paste0("cache/",modelName,"/routines.rds"))

# MCMCs
xfls <- list.files(path = "results/",pattern = paste0("dataout_",modelName),full.names = TRUE)
c(length(xfls),nchains) #converged vs run

burnin <- floor(B*propBurIn)
out <- lapply(xfls,readRDS)
out <- lapply(out,function(x)coda::as.mcmc(x[(burnin+1):B,]))
out_list <- coda::as.mcmc.list(out)
out_array <- posterior::as_draws_array(out_list)

# Relevant MCMCs statistics
out <- matrix(out_array,nrow = dim(out_array)[1]*dim(out_array)[2],ncol = dim(out_array)[3])

Log_lik <- matrix(NA,nrow = nrow(out),ncol = nrow(datain$X))
for(b in 1:nrow(out)){
  muy <- plogis(datain$X%*%out[b,1:J]); phiy <- exp(datain$Z%*%out[b,J+H])
  ystar <- mapply(function(i)extraDistr::rnsbeta(n = 1,shape1 = muy[i]*phiy[i],shape2 = phiy[i]-phiy[i]*muy[i],min = datain$lb,max = datain$ub),1:nrow(datain$X))
  Log_lik[b,] <- log(fY(y = ystar,X = datain$X,Z = datain$Z,lb = datain$lb,ub = datain$ub,pars = out[b,]))
}

lppd <- sum(log(rowMeans(exp(Log_lik))))  # Log pointwise predictive density
p_waic <- sum(apply(Log_lik, 2, var))     # Pointwise variance of log-likelihood
waic <- -2 * (lppd - p_waic)              # WAIC measure
pointwise_waic <- -2 * (log(rowMeans(exp(Log_lik))) - apply(Log_lik, 1, var))
waic_se <- sqrt(length(pointwise_waic) * var(pointwise_waic)) # WAIC Standard Error
waic2 <- loo::waic(x = Log_lik)

deviance_samples <- -2 * rowSums(Log_lik)  # Deviance for each posterior sample
D_bar <- mean(deviance_samples)
mean_log_lik <- colMeans(exp(Log_lik))  # Average likelihood across samples
D_theta_hat <- -2 * sum(log(mean_log_lik))     # Deviance at posterior mean
dic <- 2 * D_bar - D_theta_hat

posterior_table <- summarise_draws(posterior::as_draws_array(out))
hdpis <- coda::HPDinterval(obj = coda::as.mcmc(out))

# Saving all together
saveRDS(list(mcmc=out_array,stats=list(ptable=posterior_table,HDPIs=hdpis,waic=waic,waic2=waic2,waic_se=waic_se,dic=dic)),file = paste0("results/mcmc_",modelName,".rds"))


rm(list=ls());graphics.off(); source("utilities/utilities.R")
J <- 7; H <- 1 #number of predictors for mu and phi


##### Model 2 [fy = Beta] ####
modelName <- "casestudy_beta2"

### (Step 1) Model building
fY <- function(y,X,Z,lb,ub,pars){
  b <- pars[1:ncol(X)]
  g <- pars[(1+ncol(X)):(ncol(X)+ncol(Z))]
  eta <- X%*%b
  theta <- exp(eta)/(1+exp(eta))
  phi <- exp(Z%*%g)
  alpha <- theta*phi
  beta <- phi-theta*phi
  ((y - lb)^(alpha - 1)) * ((ub - y)^(beta - 1)) / ((ub - lb)^(alpha + beta - 1) * (gamma(alpha) * gamma(beta) / gamma(alpha + beta)))
}
db_internal() #create internal R functions

# Preparing for c++ conversion of R functions
filehead <- "utilities/R2Cpp_internal.head"; filetail <- "utilities/R2Cpp_internal.tail"; fileout <- "utilities/R2Cpp_internal.cpp"
verboseCpp <- TRUE
source("utilities/R2Cpp_utilities.R")
cacheDir_rcpp <-  paste0(getwd(),"/cache/",modelName)

dm_internal(fY,J,H) #create internal 1R functions
source("utilities/R2Cpp_funConverter_exec.R") #convert to c++

# Saving compiled files
save(fY,fY_internal,tud_lnFy_exec,tud_lnFy_internal,grad_fun,hess_fun,file = paste0("cache/",modelName,"/routines.rds"))
rm(list=ls())


### (Step 2) Preparing for model fitting
modelName <- "casestudy_beta2"

# Load data
load("data/VanLankveld_etal_2021.Rdata") #load data

n <- NROW(datax)
X = model.matrix(~ Age + Rel_length + desire + partner_respo + Gender_of_partner + partner_respo:Gender_of_partner,data = datax); X[,2:NCOL(X)] <- scale(X[,2:NCOL(X)])
Z <- matrix(1,nrow = n,ncol = 1)
m <- datax$intimacy_m; s <- datax$intimacy_s;

lb <- 1; ub <- 5
J <- NCOL(X); H <- 1

# Defining priors
mu_theta <- c(rep(0,J),rep(0,H)); sigma_theta <- c(rep(3.5,J),rep(1.25,H))

# Saving data
datain <- list(X=X,Z=Z,m=m,s=s,lb=lb,ub=ub,mu_theta=mu_theta,sigma_theta=sigma_theta) #input for the AGS algorithm
save(datain,file = paste0("data/datain_",modelName,".rds")) #saving inputs

### (Step 3) Running multiple chains via GNU Parallel
nchains <- 5
B <- 4e3 #samples

seedx <- as.numeric(format(Sys.time(), format = "%d%M%S"))
cmds <- sapply(1:nchains,function(u){
  cmd <- paste0(
    "R CMD BATCH --no-restore --no-save '--args B=", B, 
    " chainId=", u,
    " seedx=", seedx,
    " modelname=\"", modelName, "\"", 
    "' AGS_exec.R ", 
    "cache/",modelName,"/out_chainId_",u,".stdout"
  )
})

ncores <- nchains #ncores=nchains
btcs <- split(cmds, ceiling(seq_along(cmds)/ncores)) 
for(jobs in btcs){ 
  cat("\n Running:");for(x in jobs){cat("\n\t",strsplit(x,split = "--args")[[1]][2])}
  
  batch_file <- tempfile("job_list_"); writeLines(jobs, batch_file)
  system(paste0("parallel -j ", nchains, " < ", batch_file)); unlink(batch_file)
}

### (Step 4) Save final mcmc as array (B x Nchains x Params)
propBurIn <- 0.5
load(paste0("cache/",modelName,"/routines.rds"))

# MCMCs
xfls <- list.files(path = "results/",pattern = paste0("dataout_",modelName),full.names = TRUE)
c(length(xfls),nchains) #converged vs run

burnin <- floor(B*propBurIn)
out <- lapply(xfls,readRDS)
out <- lapply(out,function(x)coda::as.mcmc(x[(burnin+1):B,]))
out_list <- coda::as.mcmc.list(out)
out_array <- posterior::as_draws_array(out_list)

# Relevant MCMCs statistics
out <- matrix(out_array,nrow = dim(out_array)[1]*dim(out_array)[2],ncol = dim(out_array)[3])

Log_lik <- matrix(NA,nrow = nrow(out),ncol = nrow(datain$X))
for(b in 1:nrow(out)){
  muy <- plogis(datain$X%*%out[b,1:J]); phiy <- exp(datain$Z%*%out[b,J+H])
  ystar <- mapply(function(i)extraDistr::rnsbeta(n = 1,shape1 = muy[i]*phiy[i],shape2 = phiy[i]-phiy[i]*muy[i],min = datain$lb,max = datain$ub),1:nrow(datain$X))
  Log_lik[b,] <- log(fY(y = ystar,X = datain$X,Z = datain$Z,lb = datain$lb,ub = datain$ub,pars = out[b,]))
}

lppd <- sum(log(rowMeans(exp(Log_lik))))  # Log pointwise predictive density
p_waic <- sum(apply(Log_lik, 2, var))     # Pointwise variance of log-likelihood
waic <- -2 * (lppd - p_waic)              # WAIC measure
pointwise_waic <- -2 * (log(rowMeans(exp(Log_lik))) - apply(Log_lik, 1, var))
waic_se <- sqrt(length(pointwise_waic) * var(pointwise_waic)) # WAIC Standard Error
waic2 <- loo::waic(x = Log_lik)

deviance_samples <- -2 * rowSums(Log_lik)  # Deviance for each posterior sample
D_bar <- mean(deviance_samples)
mean_log_lik <- colMeans(exp(Log_lik))  # Average likelihood across samples
D_theta_hat <- -2 * sum(log(mean_log_lik))     # Deviance at posterior mean
dic <- 2 * D_bar - D_theta_hat

posterior_table <- summarise_draws(posterior::as_draws_array(out))
hdpis <- coda::HPDinterval(obj = coda::as.mcmc(out))

# Saving all together
saveRDS(list(mcmc=out_array,stats=list(ptable=posterior_table,HDPIs=hdpis,waic=waic,waic2=waic2,waic_se=waic_se,dic=dic)),file = paste0("results/mcmc_",modelName,".rds"))
