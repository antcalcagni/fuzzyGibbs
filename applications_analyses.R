#########
rm(list=ls());graphics.off(); source("utilities/utilities.R")

modelName <- "casestudy_beta2"
draws <- readRDS(paste0("results/mcmc_",modelName,".rds"))$mcmc

B <- dim(draws)[1]; JH <- dim(draws)[3]; nChains <- dim(draws)[2]

library(posterior)
options(cli.unicode = FALSE, cli.num_colors = 1)
summarise_draws(draws)
out <- matrix(draws,B*nChains,JH)

source("utilities/utilities.R")
x11(); auto.layout(JH*2);for(j in 1:(JH)){plot(out[,j],type="l",main=paste0("param_",j));hist(out[,j],main=paste0("param_",j))}
#######

defuzzify_centroid = function(mi,si,what=c("mean","mode")){
  a1 = 1+si*mi; b1 = 1+si*(1-mi) 
  if(what=="mean"){(a1)/(a1+b1)}
  else if(what=="mode"){(a1-1)/(a1+b1-2)}
}

xsup_fn = function(xsup,m,s,alpha){
  x=diff(range(xsup[beta_fn(xsup,m,s)>=alpha]))
  if(is.infinite(x)){x=0}
  return(x)
}

kaufmann_index = function(mux){
  return(2*sum(abs(mux-(mux>=0.5)))/length(mux))
}

rescale <- function(x,lb=0,ub=1){(x-lb)/(ub-lb)}   #linearly rescale into [lb,ub]
rescale_inv <- function(x,lb=0,ub=1){lb+(ub-lb)*x} #inverse of rescale

add_legend = function(...) {
  #From: https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


# External validation: Dataset 1 ------------------------------------------
## SOURCE - Colubi, A.: Statistical inference about the means of fuzzy random variables, Fuzzy Sets and Systems, 160(3), pp. 344-356 (2009)

## Model comparison
mod1 <- readRDS("results/mcmc_dataset1_beta.rds")
mod2 <- readRDS("results/mcmc_dataset1_logitnorm.rds")

c(mod1$stats$dic, mod2$stats$dic)
c(mod1$stats$waic, mod2$stats$waic)
mod1$stats$ptable

modFinal <- mod1 #final model to be chosen

## Posterior analysis (Simulation-based diagnostics)
load("data/datain_dataset1_beta.rds")
seedx <- format(Sys.time(), format = "%g%m%y"); set.seed(seedx)

s <- datain$s; m <- datain$m; n <- nrow(datain$X)

res <- optim(par = c(1,1),fn = function(x)-sum(dgamma(s,x[1],scale = x[2]/x[1],log = TRUE)),method = "L-BFGS-B",lower = c(0.01,0.01),upper = c(Inf,Inf))
alpha_s <- res$par[1]
beta_s <- res$par[2]/res$par[1] #the scale parameter is large enough!

pars_hat <- as.numeric(unlist(modFinal$stats$ptable[,2]))
muy <- plogis(datain$X%*%pars_hat[1:3])
phiy <- exp(datain$Z%*%pars_hat[4]) #dispersion parameters

z <- apply(datain$X,1,function(x)max(which(x==1)))

B <- 5e2
xsup = seq(1e-2,1-1e-2,length=1001)
Kfi1 <- matrix(NA,n,B); kfi1 = matrix(NA,1,n); Out_bxp1 = matrix(NA,n,5)
Kfi2 <- matrix(NA,n,B); kfi2 = matrix(NA,1,n); Out_bxp2 = matrix(NA,n,5)
Kfi3 <- matrix(NA,n,B); kfi3 = matrix(NA,1,n); Out_bxp3 = matrix(NA,n,5)
for(i in 1:n){
  s_hat <- rgamma(B,shape = alpha_s,scale = beta_s)
  y_hat <- mapply(function(b)rbeta(1,muy[i]*phiy[i],phiy[i]-phiy[i]*muy[i]),1:B)
  m_hat = mapply(function(b)rbeta(1,y_hat[b]*s_hat[b],s_hat[b]-s_hat[b]*y_hat[b]),1:B)
  
  # centroids
  Kfi1[i,] = mapply(function(b)defuzzify_centroid(m_hat[b],s_hat[b],"mean"),1:B)
  kfi1[i] = defuzzify_centroid(m[i],s[i],"mean")
  Out_bxp1[i,] = boxplot(Kfi1[i,],plot = FALSE)$stats
  
  # 0-cut (supports)
  Kfi2[i,] = mapply(function(b)xsup_fn(xsup,m_hat[b],s_hat[b],0.001),1:B)
  kfi2[i] = xsup_fn(xsup,m[i],s[i],0.001)
  Out_bxp2[i,] = boxplot(Kfi2[i,],plot = FALSE)$stats
  
  # Fuzziness
  Kfi3[i,] = mapply(function(b)kaufmann_index(beta_fn(xsup,m_hat[b],s_hat[b])),1:B)
  kfi3[i] = kaufmann_index(beta_fn(xsup,m[i],s[i]))
  Out_bxp3[i,] = boxplot(Kfi3[i,],plot = FALSE)$stats
}


## Figure 4
#x11(); 
tikzDevice::tikz(file='../IJAR/fig4app1.tex',width=7,height=3,sanitize = TRUE)
par(mfcol=c(1,3))
par(mar = c(2, 2, 2, 2))  # Smaller margins for each plot (bottom, left, top, and right)
par(oma = c(4, 1, 1, 1))  # Small outer margins

# Subplot1
Kfi <- Kfi1; kfi <- kfi1; Out_bxp <- Out_bxp1
ymin=min(c(min(kfi),min(Kfi)));ymax=max(c(max(kfi),max(Kfi))); fl=0.1
plot(lowess(1:n,Out_bxp[,1],f=fl)$y,1:n,type="l",col="gray",bty="n",xlim=c(ymin,ymax),xlab="",ylab=""); points(lowess(1:n,Out_bxp[,5],f=fl)$y,1:n,type="l",col="gray")
polygon(c(lowess(1:n,Out_bxp[,1],f=fl)$y, rev(lowess(1:n,Out_bxp[,5],f=fl)$y)), c(1:n, rev(1:n)), col = "#BBBDBD", border = NA)
polygon(c(lowess(1:n,Out_bxp[,2],f=fl)$y, rev(lowess(1:n,Out_bxp[,4],f=fl)$y)), c(1:n, rev(1:n)), col = "#DADADA", border = NA)
points(kfi,1:n,col="#307EB7",pch=18,cex=1.25,lty=3,type="p",lwd=2)
title("(A) Centroid", line = 1,adj=0,cex=2.45)

# Subplot2
Kfi <- Kfi2; kfi <- kfi2; Out_bxp <- Out_bxp2
ymin=min(c(min(kfi),min(Kfi)));ymax=max(c(max(kfi),max(Kfi))); fl=0.1
plot(lowess(1:n,Out_bxp[,1],f=fl)$y,1:n,type="l",col="gray",bty="n",xlim=c(ymin,ymax),xlab="",ylab=""); points(lowess(1:n,Out_bxp[,5],f=fl)$y,1:n,type="l",col="gray")
polygon(c(lowess(1:n,Out_bxp[,1],f=fl)$y, rev(lowess(1:n,Out_bxp[,5],f=fl)$y)), c(1:n, rev(1:n)), col = "#BBBDBD", border = NA)
polygon(c(lowess(1:n,Out_bxp[,2],f=fl)$y, rev(lowess(1:n,Out_bxp[,4],f=fl)$y)), c(1:n, rev(1:n)), col = "#DADADA", border = NA)
points(kfi,1:n,col="#307EB7",pch=18,cex=1.25,lty=3,type="p",lwd=2)
title("(B) Support", line = 1,adj=0,cex=2.45)

# Subplot3
Kfi <- Kfi3; kfi <- kfi3; Out_bxp <- Out_bxp3
ymin=min(c(min(kfi),min(Kfi)));ymax=max(c(max(kfi),max(Kfi))); fl=0.1
plot(lowess(1:n,Out_bxp[,1],f=fl)$y,1:n,type="l",col="gray",bty="n",xlim=c(ymin,0.3),xlab="",ylab=""); points(lowess(1:n,Out_bxp[,5],f=fl)$y,1:n,type="l",col="gray")
polygon(c(lowess(1:n,Out_bxp[,1],f=fl)$y, rev(lowess(1:n,Out_bxp[,5],f=fl)$y)), c(1:n, rev(1:n)), col = "#BBBDBD", border = NA)
polygon(c(lowess(1:n,Out_bxp[,2],f=fl)$y, rev(lowess(1:n,Out_bxp[,4],f=fl)$y)), c(1:n, rev(1:n)), col = "#DADADA", border = NA)
points(kfi,1:n,col="#307EB7",pch=18,cex=1.25,lty=3,type="p",lwd=2)
title("(C) Fuzziness", line = 1,adj=0,cex=2.45)

dev.off()

# Posterior stats
Kfi <- Kfi1; kfi <- kfi1; Out_bxp <- Out_bxp1
Kfi_qtls <- cbind(t(apply(Kfi,1,quantile,probs=c(0.05,0.95))),t(kfi))
covg <- sum(apply(Kfi_qtls,1,function(x)x[3]>=x[1]&x[3]<=x[2]))/length(kfi) #95% coverage
bayes_pv <-abs(0.5-mean(mapply(function(i)sum(Kfi[i,]>=kfi[i])/B,1:n))) #Bayes p-value
post_intervalW <- mean(apply(Kfi_qtls[,1:2],1,diff))
c(covg,bayes_pv,post_intervalW)

Kfi <- Kfi2; kfi <- kfi2; Out_bxp <- Out_bxp2
Kfi_qtls <- cbind(t(apply(Kfi,1,quantile,probs=c(0.05,0.95))),t(kfi))
covg <- sum(apply(Kfi_qtls,1,function(x)x[3]>=x[1]&x[3]<=x[2]))/length(kfi) #95% coverage
bayes_pv <-abs(0.5-mean(mapply(function(i)sum(Kfi[i,]>=kfi[i])/B,1:n))) #Bayes p-value
post_intervalW <- mean(apply(Kfi_qtls[,1:2],1,diff))
c(covg,bayes_pv,post_intervalW)

Kfi <- Kfi3; kfi <- kfi3; Out_bxp <- Out_bxp3
Kfi_qtls <- cbind(t(apply(Kfi,1,quantile,probs=c(0.05,0.95))),t(kfi))
covg <- sum(apply(Kfi_qtls,1,function(x)x[3]>=x[1]&x[3]<=x[2]))/length(kfi) #95% coverage
bayes_pv <-abs(0.5-mean(mapply(function(i)sum(Kfi[i,]>=kfi[i])/B,1:n))) #Bayes p-value
post_intervalW <- mean(apply(Kfi_qtls[,1:2],1,diff))
c(covg,bayes_pv,post_intervalW)


# External validation: Dataset 2 ------------------------------------------
## SOURCE - Kim, B., & Bishu, R. R. (1998). Evaluation of fuzzy linear regression models by comparing membership functions. Fuzzy sets and systems, 100(1-3), 343-352

modFinal <- readRDS("results/mcmc_dataset2_lognormal.rds")

## Posterior analysis (Simulation-based diagnostics)
load("data/datain_dataset2_lognormal.rds")
seedx <- format(Sys.time(), format = "%g%m%y"); set.seed(seedx)

s <- datain$s; m <- datain$m; n <- nrow(datain$X)

res <- optim(par = c(1,1),fn = function(x)-sum(dgamma(s,x[1],scale = x[2]/x[1],log = TRUE)),method = "L-BFGS-B",lower = c(0.01,0.01),upper = c(Inf,Inf))
alpha_s <- res$par[1]
beta_s <- res$par[2]/res$par[1] #the scale parameter is large enough!

pars_hat <- as.numeric(unlist(modFinal$stats$ptable[,2]))
muy <- datain$X%*%pars_hat[1:4]
phiy <- exp(datain$Z%*%pars_hat[5]) #dispersion parameters

B <- 5e2
xsup = seq(1e-2,1-1e-2,length=1e4+1)
Kfi1 <- matrix(NA,n,B); kfi1 = matrix(NA,1,n); Out_bxp1 = matrix(NA,n,5)
Kfi2 <- matrix(NA,n,B); kfi2 = matrix(NA,1,n); Out_bxp2 = matrix(NA,n,5)
Kfi3 <- matrix(NA,n,B); kfi3 = matrix(NA,1,n); Out_bxp3 = matrix(NA,n,5)
for(i in 1:n){
  s_hat <- rgamma(B,shape = alpha_s,scale = beta_s)
  y_hat <- mapply(function(b)rlnorm(1,muy[i],phiy[i]),1:B)
  lb0 <- min(y_hat)-0.1; ub0 <- max(y_hat)+0.1
  m_hat01 <- mapply(function(b)rbeta(1,rescale(y_hat[i],lb0,ub0)*s[i],s[i]-s[i]*rescale(y_hat[i],lb0,ub0)),1:B)
  m_hat <- rescale_inv(m_hat01,lb0,ub0)
  
  # centroids
  Kfi1[i,] = mapply(function(b)defuzzify_centroid(m_hat[b],s_hat[b],"mean"),1:B)
  kfi1[i] = defuzzify_centroid(m[i],s[i],"mean")
  Out_bxp1[i,] = boxplot(Kfi1[i,],plot = FALSE)$stats
  
  # 0-cut (supports)
  Kfi2[i,] = mapply(function(b)rescale_inv(xsup_fn(xsup,m_hat01[b],s_hat[b],0.0001),lb0,ub0),1:B)
  kfi2[i] = rescale_inv(xsup_fn(xsup,rescale(m[i],datain$lb,datain$ub),s[i],0.0001),datain$lb,datain$ub)
  Out_bxp2[i,] = boxplot(Kfi2[i,],plot = FALSE)$stats
  
  # Fuzziness
  Kfi3[i,] = mapply(function(b)kaufmann_index(beta_fn(xsup,m_hat01[b],s_hat[b])),1:B)
  kfi3[i] = kaufmann_index(beta_fn(xsup,rescale(m[i],datain$lb,datain$ub),s[i]))
  Out_bxp3[i,] = boxplot(Kfi3[i,],plot = FALSE)$stats
}


## Figure 5
#x11(); 
tikzDevice::tikz(file='../IJAR/fig5app1.tex',width=7,height=3,sanitize = TRUE)
par(mfcol=c(1,3))
par(mar = c(2, 2, 2, 2))  # Smaller margins for each plot (bottom, left, top, and right)
par(oma = c(4, 1, 1, 1))  # Small outer margins

# Subplot1
Kfi <- Kfi1; kfi <- kfi1; Out_bxp <- Out_bxp1
ymin=min(c(min(kfi),min(Kfi)));ymax=max(c(max(kfi),max(Kfi))); fl=0.25
plot(lowess(1:n,Out_bxp[,1],f=fl)$y,1:n,type="l",col="gray",bty="n",xlim=c(-2,15),xlab="",ylab=""); points(lowess(1:n,Out_bxp[,5],f=fl)$y,1:n,type="l",col="gray")
polygon(c(lowess(1:n,Out_bxp[,1],f=fl)$y, rev(lowess(1:n,Out_bxp[,5],f=fl)$y)), c(1:n, rev(1:n)), col = "#BBBDBD", border = NA)
polygon(c(lowess(1:n,Out_bxp[,2],f=fl)$y, rev(lowess(1:n,Out_bxp[,4],f=fl)$y)), c(1:n, rev(1:n)), col = "#DADADA", border = NA)
points(kfi,1:n,col="#307EB7",pch=18,cex=1.25,lty=3,type="p",lwd=2)
title("(A) Centroid", line = 1,adj=0,cex=2.45)

# Subplot2
Kfi <- Kfi2; kfi <- kfi2; Out_bxp <- Out_bxp2
ymin=min(c(min(kfi,na.rm=TRUE),min(Kfi,na.rm=TRUE)));ymax=max(c(max(kfi,na.rm=TRUE),max(Kfi,na.rm=TRUE))); fl=0.25
plot(lowess(1:n,Out_bxp[,1],f=fl)$y,1:n,type="l",col="gray",bty="n",xlim=c(-5,50),xlab="",ylab=""); points(lowess(1:n,Out_bxp[,5],f=fl)$y,1:n,type="l",col="gray")
polygon(c(lowess(1:n,Out_bxp[,1],f=fl)$y, rev(lowess(1:n,Out_bxp[,5],f=fl)$y)), c(1:n, rev(1:n)), col = "#BBBDBD", border = NA)
polygon(c(lowess(1:n,Out_bxp[,2],f=fl)$y, rev(lowess(1:n,Out_bxp[,4],f=fl)$y)), c(1:n, rev(1:n)), col = "#DADADA", border = NA)
points(kfi,1:n,col="#307EB7",pch=18,cex=1.25,lty=3,type="p",lwd=2)
title("(B) Support", line = 1,adj=0,cex=2.45)

# Subplot3
Kfi <- Kfi3; kfi <- kfi3; Out_bxp <- Out_bxp3
ymin=min(c(min(kfi,na.rm=TRUE),min(Kfi,na.rm=TRUE)));ymax=max(c(max(kfi,na.rm=TRUE),max(Kfi,na.rm=TRUE))); fl=0.25
plot(lowess(1:n,Out_bxp[,1],f=fl)$y,1:n,type="l",col="gray",bty="n",xlim=c(-0.05,0.6),xlab="",ylab=""); points(lowess(1:n,Out_bxp[,5],f=fl)$y,1:n,type="l",col="gray")
polygon(c(lowess(1:n,Out_bxp[,1],f=fl)$y, rev(lowess(1:n,Out_bxp[,5],f=fl)$y)), c(1:n, rev(1:n)), col = "#BBBDBD", border = NA)
polygon(c(lowess(1:n,Out_bxp[,2],f=fl)$y, rev(lowess(1:n,Out_bxp[,4],f=fl)$y)), c(1:n, rev(1:n)), col = "#DADADA", border = NA)
points(kfi,1:n,col="#307EB7",pch=18,cex=1.25,lty=3,type="p",lwd=2)
title("(C) Fuzziness", line = 1,adj=0,cex=2.45)

dev.off()

# Posterior stats
Kfi <- Kfi1; kfi <- kfi1; Out_bxp <- Out_bxp1
Kfi_qtls <- cbind(t(apply(Kfi,1,quantile,probs=c(0.05,0.95))),t(kfi))
covg <- sum(apply(Kfi_qtls,1,function(x)x[3]>=x[1]&x[3]<=x[2]))/length(kfi) #95% coverage
bayes_pv <-abs(0.5-mean(mapply(function(i)sum(Kfi[i,]>=kfi[i],na.rm=TRUE)/B,1:n))) #Bayes p-value
post_intervalW <- mean(apply(Kfi_qtls[,1:2],1,diff))
c(covg,bayes_pv,post_intervalW)

Kfi <- Kfi2; kfi <- kfi2; Out_bxp <- Out_bxp2
Kfi_qtls <- cbind(t(apply(Kfi,1,quantile,probs=c(0.05,0.95),na.rm=TRUE)),t(kfi))
covg <- sum(apply(Kfi_qtls,1,function(x)x[3]>=x[1]&x[3]<=x[2]))/length(kfi) #95% coverage
bayes_pv <-abs(0.5-mean(mapply(function(i)sum(Kfi[i,]>=kfi[i],na.rm=TRUE)/B,1:n),na.rm=TRUE)) #Bayes p-value
post_intervalW <- mean(apply(Kfi_qtls[,1:2],1,diff))
c(covg,bayes_pv,post_intervalW)

Kfi <- Kfi3; kfi <- kfi3; Out_bxp <- Out_bxp3
Kfi_qtls <- cbind(t(apply(Kfi,1,quantile,probs=c(0.05,0.95),na.rm=TRUE)),t(kfi))
covg <- sum(apply(Kfi_qtls,1,function(x)x[3]>=x[1]&x[3]<=x[2]))/length(kfi) #95% coverage
bayes_pv <-abs(0.5-mean(mapply(function(i)sum(Kfi[i,]>=kfi[i],na.rm=TRUE)/B,1:n),na.rm=TRUE)) #Bayes p-value
post_intervalW <- mean(apply(Kfi_qtls[,1:2],1,diff))
c(covg,bayes_pv,post_intervalW)


# External validation: Dataset 3 ------------------------------------------
## SOURCE - Jiang, H., Kwong, C. K., Chan, C. Y., & Yung, K. L. (2019). A multi-objective evolutionary approach for fuzzy regression analysis. Expert Systems with Applications, 130, 225-235

modFinal <- readRDS("results/mcmc_dataset3_kumar.rds")

## Posterior analysis (Simulation-based diagnostics)
load("data/datain_dataset3_kumar.rds")
seedx <- format(Sys.time(), format = "%g%m%y"); set.seed(as.numeric(seedx))

s <- datain$s; m <- datain$m; n <- nrow(datain$X)

res <- optim(par = c(1,1),fn = function(x)-sum(dgamma(s,x[1],scale = x[2]/x[1],log = TRUE)),method = "L-BFGS-B",lower = c(0.01,0.01),upper = c(Inf,Inf))
alpha_s <- res$par[1]
beta_s <- res$par[2]/res$par[1] #the scale parameter is large enough!

iid <- sample(1:n,30,FALSE)

pars_hat <- as.numeric(unlist(modFinal$stats$ptable[,2]))
muy <- plogis(datain$X[iid,]%*%pars_hat[1:2])
phiy <- plogis(datain$Z%*%pars_hat[3]) #dispersion parameters
qhiy <- (log(1-0.5^phiy)/log(muy)) 

B <- 5e2
xsup = seq(1e-2,1-1e-2,length=1e4+1)
Kfi1 <- matrix(NA,n,B); kfi1 = matrix(NA,1,n); Out_bxp1 = matrix(NA,n,5)
Kfi2 <- matrix(NA,n,B); kfi2 = matrix(NA,1,n); Out_bxp2 = matrix(NA,n,5)
Kfi3 <- matrix(NA,n,B); kfi3 = matrix(NA,1,n); Out_bxp3 = matrix(NA,n,5)
for(i in 1:n){
  s_hat <- rgamma(B,shape = alpha_s,scale = beta_s)
  y_hat <- mapply(function(b)extraDistr::rkumar(1,qhiy[i],1/phiy[i]),1:B)
  lb0 <- 0; ub0 <- 1
  m_hat01 <- mapply(function(b)rbeta(1,rescale(y_hat[i],lb0,ub0)*s[i],s[i]-s[i]*rescale(y_hat[i],lb0,ub0)),1:B)
  m_hat <- rescale_inv(m_hat01,lb0,ub0)
  
  # centroids
  Kfi1[i,] = mapply(function(b)defuzzify_centroid(m_hat[b],s_hat[b],"mean"),1:B)
  kfi1[i] = defuzzify_centroid(m[iid][i],s[iid][i],"mean")
  Out_bxp1[i,] = boxplot(Kfi1[i,],plot = FALSE)$stats
  
  # 0-cut (supports)
  Kfi2[i,] = mapply(function(b)rescale_inv(xsup_fn(xsup,m_hat01[b],s_hat[b],0.0001),lb0,ub0),1:B)
  kfi2[i] = rescale_inv(xsup_fn(xsup,rescale(m[iid][i],datain$lb,datain$ub),s[iid][i],0.0001),datain$lb,datain$ub)
  Out_bxp2[i,] = boxplot(Kfi2[i,],plot = FALSE)$stats
  
  # Fuzziness
  Kfi3[i,] = mapply(function(b)kaufmann_index(beta_fn(xsup,m_hat01[b],s_hat[b])),1:B)
  kfi3[i] = kaufmann_index(beta_fn(xsup,rescale(m[iid][i],datain$lb,datain$ub),s[iid][i]))
  Out_bxp3[i,] = boxplot(Kfi3[i,],plot = FALSE)$stats
}


## Figure 6
#x11(); 
tikzDevice::tikz(file='../IJAR/fig6app1.tex',width=7,height=3,sanitize = TRUE)
par(mfcol=c(1,3))
par(mar = c(2, 2, 2, 2))  # Smaller margins for each plot (bottom, left, top, and right)
par(oma = c(4, 1, 1, 1))  # Small outer margins

# Subplot1
Kfi <- Kfi1; kfi <- kfi1; Out_bxp <- Out_bxp1
ymin=min(c(min(kfi),min(Kfi)));ymax=max(c(max(kfi),max(Kfi))); fl=0.5
plot(lowess(1:n,Out_bxp[,1],f=fl)$y,1:n,type="l",col="gray",bty="n",xlim=c(ymin,ymax),xlab="",ylab=""); points(lowess(1:n,Out_bxp[,5],f=fl)$y,1:n,type="l",col="gray")
polygon(c(lowess(1:n,Out_bxp[,1],f=fl)$y, rev(lowess(1:n,Out_bxp[,5],f=fl)$y)), c(1:n, rev(1:n)), col = "#BBBDBD", border = NA)
polygon(c(lowess(1:n,Out_bxp[,2],f=fl)$y, rev(lowess(1:n,Out_bxp[,4],f=fl)$y)), c(1:n, rev(1:n)), col = "#DADADA", border = NA)
points(kfi,1:n,col="#307EB7",pch=18,cex=1.25,lty=3,type="p",lwd=2)
title("(A) Centroid", line = 1,adj=0,cex=2.45)

# Subplot2
Kfi <- Kfi2; kfi <- kfi2; Out_bxp <- Out_bxp2
ymin=min(c(min(kfi,na.rm=TRUE),min(Kfi,na.rm=TRUE)));ymax=max(c(max(kfi,na.rm=TRUE),max(Kfi,na.rm=TRUE))); fl=0.25
plot(lowess(1:n,Out_bxp[,1],f=fl)$y,1:n,type="l",col="gray",bty="n",xlim=c(ymin,ymax),xlab="",ylab=""); points(lowess(1:n,Out_bxp[,5],f=fl)$y,1:n,type="l",col="gray")
polygon(c(lowess(1:n,Out_bxp[,1],f=fl)$y, rev(lowess(1:n,Out_bxp[,5],f=fl)$y)), c(1:n, rev(1:n)), col = "#BBBDBD", border = NA)
polygon(c(lowess(1:n,Out_bxp[,2],f=fl)$y, rev(lowess(1:n,Out_bxp[,4],f=fl)$y)), c(1:n, rev(1:n)), col = "#DADADA", border = NA)
points(kfi,1:n,col="#307EB7",pch=18,cex=1.25,lty=3,type="p",lwd=2)
title("(B) Support", line = 1,adj=0,cex=2.45)

# Subplot3
Kfi <- Kfi3; kfi <- kfi3; Out_bxp <- Out_bxp3
ymin=min(c(min(kfi,na.rm=TRUE),min(Kfi,na.rm=TRUE)));ymax=max(c(max(kfi,na.rm=TRUE),max(Kfi,na.rm=TRUE))); fl=0.25
plot(lowess(1:n,Out_bxp[,1],f=fl)$y,1:n,type="l",col="gray",bty="n",xlim=c(ymin,ymax),xlab="",ylab=""); points(lowess(1:n,Out_bxp[,5],f=fl)$y,1:n,type="l",col="gray")
polygon(c(lowess(1:n,Out_bxp[,1],f=fl)$y, rev(lowess(1:n,Out_bxp[,5],f=fl)$y)), c(1:n, rev(1:n)), col = "#BBBDBD", border = NA)
polygon(c(lowess(1:n,Out_bxp[,2],f=fl)$y, rev(lowess(1:n,Out_bxp[,4],f=fl)$y)), c(1:n, rev(1:n)), col = "#DADADA", border = NA)
points(kfi,1:n,col="#307EB7",pch=18,cex=1.25,lty=3,type="p",lwd=2)
title("(C) Fuzziness", line = 1,adj=0,cex=2.45)

dev.off()

# Posterior stats
Kfi <- Kfi1; kfi <- kfi1; Out_bxp <- Out_bxp1
Kfi_qtls <- cbind(t(apply(Kfi,1,quantile,probs=c(0.05,0.95))),t(kfi))
covg <- sum(apply(Kfi_qtls,1,function(x)x[3]>=x[1]&x[3]<=x[2]))/length(kfi) #95% coverage
bayes_pv <-abs(0.5-mean(mapply(function(i)sum(Kfi[i,]>=kfi[i],na.rm=TRUE)/B,1:n))) #Bayes p-value
post_intervalW <- mean(apply(Kfi_qtls[,1:2],1,diff))
c(covg,bayes_pv,post_intervalW)

Kfi <- Kfi2; kfi <- kfi2; Out_bxp <- Out_bxp2
Kfi_qtls <- cbind(t(apply(Kfi,1,quantile,probs=c(0.05,0.95),na.rm=TRUE)),t(kfi))
covg <- sum(apply(Kfi_qtls,1,function(x)x[3]>=x[1]&x[3]<=x[2]))/length(kfi) #95% coverage
bayes_pv <-abs(0.5-mean(mapply(function(i)sum(Kfi[i,]>=kfi[i],na.rm=TRUE)/B,1:n),na.rm=TRUE)) #Bayes p-value
post_intervalW <- mean(apply(Kfi_qtls[,1:2],1,diff))
c(covg,bayes_pv,post_intervalW)

Kfi <- Kfi3; kfi <- kfi3; Out_bxp <- Out_bxp3
Kfi_qtls <- cbind(t(apply(Kfi,1,quantile,probs=c(0.05,0.95),na.rm=TRUE)),t(kfi))
covg <- sum(apply(Kfi_qtls,1,function(x)x[3]>=x[1]&x[3]<=x[2]))/length(kfi) #95% coverage
bayes_pv <-abs(0.5-mean(mapply(function(i)sum(Kfi[i,]>=kfi[i],na.rm=TRUE)/B,1:n),na.rm=TRUE)) #Bayes p-value
post_intervalW <- mean(apply(Kfi_qtls[,1:2],1,diff))
c(covg,bayes_pv,post_intervalW)


# Case Study --------------------------------------------------------------
## SOURCE - Calcagnì et al., A Bayesian beta linear model to analyze fuzzy rating responses. Book of short papers of the SIS Conference 2022

## Descriptive stats 
load("data/VanLankveld_etal_2021.Rdata") #load data
table(datax$Gender_of_partner)
aggregate(Age~Gender_of_partner,data=datax,FUN = function(x)c(mean(x),sd(x)))
aggregate(Rel_length~Gender_of_partner,data=datax,FUN = function(x)c(mean(x),sd(x)))

## Model comparison
mod2 <- readRDS("results/mcmc_casestudy_beta2.rds")
mod1 <- readRDS("results/mcmc_casestudy_beta1.rds")
c(mod1$stats$waic,mod2$stats$waic)

## Posterior analysis
A <- as.matrix(mod1$stats$ptable)[,-1]
A <- apply(A,2,as.numeric)

# Table 3
A <- round(cbind(A[,c(1,3)],mod1$stats$HDPIs,A[,c(8:9)]),3); A[7,1] <- exp(A[7,1])
colnames(A) = c("mean","sd","95% HDI lb","95% HDI ub","ESS-bulk","ESS-tail")
Xtab_tex = xtable::xtable(A)
attributes(Xtab_tex)$caption = ""
attributes(Xtab_tex)$label = "tab3app"
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})

# Figure 7
tikzDevice::tikz(file='../IJAR/fig7app1.tex',width=6.5,height=4)
Y_posts <- matrix(mod1$mcmc,nrow = dim(mod1$mcmc)[1]*dim(mod1$mcmc)[2],ncol = dim(mod1$mcmc)[3])
cols = c("orangered3","orchid4","palegreen4","peru","darkorange4","gray15","black")
par(mfrow=c(1,2),mai=c(1.35, 0.85, 0.45, 0.15))
plot(density(Y_posts[,1],bw = 0.25),bty="n",xlab="",ylab="",main="",col=cols[1],xlim=c(-1.85,2.25),ylim=c(0,1.8),lwd=1.35); title("(A)",line=-1,adj=0)
for(j in 2:6){lines(density(Y_posts[,j],bw = 0.25),col=cols[j],lwd=1.35)}
abline(v = 0,lty=2,lwd=2,col="gray")
plot(density(exp(Y_posts[,7]),bw = 0.45),bty="n",xlab="",ylab="",main="",col=cols[7],xlim=c(10,30),ylim=c(0,0.22),lwd=1.35); title("(B)",line=-1,adj=0)
add_legend("bottom",fill = cols[1:6],
           legend = c("gender\\_p=M","age","rel\\_len","sex\\_desire","gender\\_p=F"),border = FALSE,bty = "n",ncol = 3)
dev.off()




