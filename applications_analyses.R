#########
rm(list=ls());graphics.off(); source("utilities/utilities.R")

modelName <- "casestudy_beta2"
draws <- readRDS(paste0("results/mcmc_",modelName,".rds"))$mcmc

B <- dim(draws)[1]; JH <- dim(draws)[3]; nChains <- dim(draws)[2]

library(posterior)
options(cli.unicode = FALSE, cli.num_colors = 1)
summarise_draws(draws)
out <- matrix(draws,B*nChains,JH)

source("utils.R")

x11(); auto.layout(JH*2);for(j in 1:(JH)){plot(out[,j],type="l",main=paste0("param_",j));hist(out[,j],main=paste0("param_",j))}
#######



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
tikzDevice::tikz(file='fig4app1.tex',width=7,height=3,sanitize = TRUE)
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
  m_hat01 <- mapply(function(b)rbeta(1,rescale(y_hat[b],lb0,ub0)*s_hat[b],s_hat[b]-s_hat[b]*rescale(y_hat[b],lb0,ub0)),1:B)
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
tikzDevice::tikz(file='fig5app1.tex',width=7,height=3,sanitize = TRUE)
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
  m_hat01 <- mapply(function(b)rbeta(1,rescale(y_hat[b],lb0,ub0)*s_hat[i],s_hat[i]-s_hat[i]*rescale(y_hat[b],lb0,ub0)),1:B)
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
tikzDevice::tikz(file='../all_versions_rev/fig6app1.tex',width=7,height=3,sanitize = TRUE)
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
tikzDevice::tikz(file='fig7app1.tex',width=6.5,height=4)
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


## PPC analysis
m <- datax$intimacy_m 
s <- datax$intimacy_s;
res <- optim(par = c(1,1),fn = function(x)-sum(dgamma(s,x[1],scale = x[2]/x[1],log = TRUE)),method = "L-BFGS-B",lower = c(0.01,0.01),upper = c(Inf,Inf))
alpha_s <- res$par[1]
beta_s <- res$par[2]/res$par[1] #the scale parameter is large enough!

n <- NROW(datax)
X = model.matrix(~ Age + Rel_length + desire + partner_respo + Gender_of_partner,data = datax); X[,2:NCOL(X)] <- scale(X[,2:NCOL(X)])
Z <- matrix(1,nrow = n,ncol = 1)

pars_hat <- as.numeric(unlist(mod1$stats$ptable[,2]))
muy <- plogis(X%*%pars_hat[1:6])
phiy <- exp(Z%*%pars_hat[7]) #dispersion parameters

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
  kfi1[i] = defuzzify_centroid(rescale(m[i],1,5),s[i],"mean")
  Out_bxp1[i,] = boxplot(Kfi1[i,],plot = FALSE)$stats
  
  # 0-cut (supports)
  Kfi2[i,] = mapply(function(b)xsup_fn(xsup,m_hat[b],s_hat[b],0.001),1:B)
  kfi2[i] = xsup_fn(xsup,rescale(m[i],1,5),s[i],0.001)
  Out_bxp2[i,] = boxplot(Kfi2[i,],plot = FALSE)$stats
  
  # Fuzziness
  Kfi3[i,] = mapply(function(b)kaufmann_index(beta_fn(xsup,m_hat[b],s_hat[b])),1:B)
  kfi3[i] = kaufmann_index(beta_fn(xsup,rescale(m[i],1,5),s[i]))
  Out_bxp3[i,] = boxplot(Kfi3[i,],plot = FALSE)$stats
}

## Figure S.21 (Supplementary Materials)
#x11(); 
tikzDevice::tikz(file='casestudy_fig10.tex',width=7,height=3,sanitize = TRUE)
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


# Case Study - Supplementary Materials ------------------------------------

source("utils.R")
load("data/VanLankveld_etal_2021.Rdata") #load data

n <- NROW(datax)
X = model.matrix(~ Age + Rel_length + desire + partner_respo + Gender_of_partner + partner_respo + Gender_of_partner,data = datax); X[,2:NCOL(X)] <- scale(X[,2:NCOL(X)])
Z <- matrix(1,nrow = n,ncol = 1)
m <- datax$intimacy_m; s <- datax$intimacy_s;

pars_hat <- as.numeric(unlist(mod1$stats$ptable[,2]))
muy <- plogis(X%*%pars_hat[1:6])
phiy <- exp(Z%*%pars_hat[7]) #dispersion parameters

res <- optim(par = c(1,1),fn = function(x)-sum(dgamma(s,x[1],scale = x[2]/x[1],log = TRUE)),method = "L-BFGS-B",lower = c(0.01,0.01),upper = c(Inf,Inf))
alpha_s <- res$par[1]
beta_s <- res$par[2]/res$par[1] #the scale parameter is large enough!

lb <- 1; ub <- 5
J <- NCOL(X); H <- 1

### Convert beta fuzzy data into trg fuzzy data
Ytrg <- t(mapply(function(i)beta2trg(rescale(m[i],1,5),s[i]),1:n)) #from beta to trg fuzzy numbers (output: lb,ub,m)
Ytrg <- apply(Ytrg,2,function(x){1+unlist(x)*4})

### Generate data according to our model
set.seed(112)
Y_hat <- S_hat <- M_hat <- matrix(nrow = n,ncol = B)
for(i in 1:n){
  cat("\t",i)
  S_hat[i,] <- rgamma(B,shape = alpha_s,scale = beta_s)
  Y_hat[i,] <- mapply(function(b)rbeta(1,muy[i]*phiy[i],phiy[i]-phiy[i]*muy[i]),1:B)
  M_hat[i,] = mapply(function(b)rbeta(1,Y_hat[i,b]*S_hat[i,b],S_hat[i,b]-S_hat[i,b]*Y_hat[i,b]),1:B)
}  

### fuzzy-OLS (D'Urso, 2003)

## estimate parameters of the crisp-input/fuzzy-output (interactive) linear regression model
out_fols <- optim(fn = flr1,par = rep(1,NCOL(X)+2),X,m=Ytrg[,3],l=Ytrg[,3]-Ytrg[,1],r=Ytrg[,2]-Ytrg[,3],X=X,J=NCOL(X)-1,control=list(trace=1,maxit=1e3,reltol=1e-4,fnscale=1e3),method="Nelder-Mead")
beta_fols <- out_fols$par

## estimate var(beta_est) and CI(beta) via standard bootstrap
set.seed(131)
B=1000
Beta_boot = matrix(data = NA,nrow = B,ncol = length(beta_fols)) #each column is a regression parameter (in this case, four)
for(b in 1:B){
  iid = sample(x = 1:n,size = n,replace = TRUE)
  Y_boot = Ytrg[iid,]
  Beta_boot[b,] = optim(fn = flr1,par = rep(1,J+2),X,m=Y_boot[,2],l=Y_boot[,2]-Y_boot[,1],r=Y_boot[,2]+Y_boot[,3],X=X,J=NCOL(X)-1,control=list(trace=0,maxit=1e3,reltol=1e-4,fnscale=1e3),method="Nelder-Mead")$par
  if(b%%100 == 0){print(b)}
}

beta_sd = apply(Beta_boot,2,sd) #approximated
Beta_ci = apply(Beta_boot,2,function(x)quantile(x,probs = c(0.05/2, 1 - 0.05/2), type = 1)) #approximated
Beta_fols_tab <- cbind(beta_fols,beta_sd,t(Beta_ci)); rownames(Beta_fols_tab) <- c(colnames(X),"beta_l","beta_r")

## evaluate the goodness-of-fit
y_hat <- X%*%beta_fols[1:J]
l_hat <- y_hat*beta_fols[J+1]
r_hat <- y_hat*beta_fols[J+2]

constraints_violation_1 <- c(sum(y_hat<1),sum(y_hat>5),sum(l_hat<0),sum(r_hat<0))
l_hat[l_hat<0] <- 1e-4; r_hat[r_hat<0] <- 1e-4; y_hat[y_hat<1] <- 1; y_hat[y_hat>5] <- 5

# evaluate estimated fuzzy numbers using Grzegorzewski and Romaniuk's metrics
# See: Grzegorzewski, P., Romaniuk, M. (2022) Bootstrap methods for fuzzy data Uncertainty and Imprecision in Decision Making and Decision Support: New Advances, Challenges, and Perspectives, pp. 28-47 Springer
B <- 1e2
BD_fgibbs <- matrix(nrow = n,ncol = B); bd_fols <- matrix(nrow = n,ncol = 1)
FZ_fgibbs <- matrix(nrow = n,ncol = B); fz_fols <- matrix(nrow = n,ncol = 1)
xsup <- seq(1,5,length.out=1001)
for(i in 1:n){
  cat("\t",i)
  
  # convert estimated beta fuzzy numbers into trg fuzzy numbers
  U <- t(mapply(function(b)unlist(beta2trg(m = M_hat[i,b],s = S_hat[i,b])),1:B)) #lb,ub,m
  Uresc <- rescale_inv(U,1,5)
  
  # dist Bertol
  BD_fgibbs[i,] <- mapply(function(b){
    s1 <- cbind(Uresc[b,1],Uresc[b,3],Uresc[b,3],Uresc[b,2]); s2 <- cbind(Ytrg[i,1],Ytrg[i,3],Ytrg[i,3],Ytrg[i,2])
    FuzzyResampling::BertoluzzaDistance(s1,s2)
  },1:B)
  s1 <- cbind(y_hat[i]-l_hat[i],y_hat[i],y_hat[i],y_hat[i]+l_hat[i]); s2 <- cbind(Ytrg[i,1],Ytrg[i,3],Ytrg[i,3],Ytrg[i,2])
  bd_fols[i] <- FuzzyResampling::BertoluzzaDistance(s1,s2)
  
  # Fuzz
  S1 <- cbind(Uresc[,1],Uresc[,3],Uresc[,3],Uresc[,2])
  FZ_fgibbs[i,] <- FuzzyResampling::CalculateFuzziness(S1)/FuzzyResampling::CalculateFuzziness(s2)
  fz_fols[i] <- FuzzyResampling::CalculateFuzziness(s1)/FuzzyResampling::CalculateFuzziness(s2)
}

cls <- c(adjustcolor("#FFB90F", alpha.f = 0.7),adjustcolor("#A4D3EE", alpha.f = 0.7))

## Figure S.22 (Supplementary materials)
#x11(); 
tikzDevice::tikz(file='casestudy_add_fig1.tex',width=7,height=8,sanitize = TRUE)
par(mfcol=c(2,1))
par(mar = c(2, 2, 3, 2))  # Smaller margins for each plot (bottom, left, top, and right)
par(oma = c(4, 1, 1, 1))  # Small outer margins

boxplot(BD_fgibbs,frame=FALSE,col=cls[1],outline=FALSE,border=FALSE,medcol="#CD6600",medlwd=6,axes=FALSE,ylim=c(0,quantile(BD_fgibbs,0.95)))
points(1:B,rep(quantile(bd_fols,0.50),B),lwd=4,col=cls[2],type="l"); points(1:B,rep(quantile(bd_fols,0.25),B),lwd=2,col="gray",type="l",lty=2);points(1:B,rep(quantile(bd_fols,0.75),B),lwd=2,col="gray",type="l",lty=2)    
axis(side = 2,at = round(seq(0,quantile(BD_fgibbs,0.95),length.out=11),2))
title(main = "(A) Bertoluzza distance",line = 1,adj=0)

boxplot(1-FZ_fgibbs,frame=FALSE,col=cls[1],outline=FALSE,border=FALSE,medcol="#CD6600",medlwd=6,axes=FALSE,ylim=c(quantile(1-FZ_fgibbs,0.05),quantile(1-FZ_fgibbs,0.95)))
points(1:B,rep(quantile(1-fz_fols,0.50),B),lwd=4,col=cls[2],type="l"); points(1:B,rep(quantile(1-fz_fols,0.25),B),lwd=2,col="gray",type="l",lty=2);points(1:B,rep(quantile(1-fz_fols,0.75),B),lwd=2,col="gray",type="l",lty=2)    
axis(side = 2,at = round(seq(quantile(1-FZ_fgibbs,0.05),quantile(1-FZ_fgibbs,0.95),length.out=11),2))
abline(h = 0.0,lty=1,lwd=3,col="#696969")
title(main = "(B) Fuzziness ratio",line = 1,adj=0)

legend("bottom",fill = cls,legend = c("fCM","fOLS"),border = FALSE,bty = "n",ncol = 2,cex=1.25)
dev.off()

U1 <- cbind(summary(apply(BD_fgibbs,1,mean)),summary(apply(bd_fols,1,mean)),
           summary(apply(1-FZ_fgibbs,1,mean)),summary(apply(1-fz_fols,1,mean)))



### Tanaka's approach (Possibilistic linear regression)
yy <- Ytrg[,3]
lr <- apply(cbind(Ytrg[,3]-Ytrg[,1],Ytrg[,2]-Ytrg[,3]),1,sum)
out_plr <- fuzzyreg::plr(x = X,y = cbind(yy,lr),h = 0)
beta_plr <- out_plr$coef[,1:2]

## estimate var(beta_est) and CI(beta) via standard bootstrap
set.seed(231)
B=1000
Beta_boot_c = matrix(data = NA,nrow = B,ncol = nrow(beta_plr))
Beta_boot_l = matrix(data = NA,nrow = B,ncol = nrow(beta_plr))
for(b in 1:B){
  iid = sample(x = 1:n,size = n,replace = TRUE)
  Y_boot = cbind(yy,lr)[iid,]
  out <- fuzzyreg::plr(x = X,y = Y_boot,h = 0)
  Beta_boot_c[b,] <- out$coef[,1]
  Beta_boot_l[b,] <- out$coef[,2]
  if(b%%100 == 0){print(b)}
}

beta_sd_c = apply(Beta_boot_c,2,sd) #approximated
beta_sd_l = apply(Beta_boot_l,2,sd) #approximated
Beta_ci_c = apply(Beta_boot_c,2,function(x)quantile(x,probs = c(0.05/2, 1 - 0.05/2), type = 1)) #approximated
Beta_ci_l = apply(Beta_boot_l,2,function(x)quantile(x,probs = c(0.05/2, 1 - 0.05/2), type = 1)) #approximated
Beta_plr_tab <- cbind(beta_plr[,1],beta_sd_c,t(Beta_ci_c),
                      beta_plr[,2],beta_sd_l,t(Beta_ci_l));


## evaluate the goodness-of-fit
y_hat <- X%*%beta_plr[,1]
l_hat <- X%*%beta_plr[,2]
r_hat <- l_hat

constraints_violation_2 <- c(sum(y_hat<1),sum(y_hat>5),sum(l_hat<0),sum(r_hat<0))
#l_hat[l_hat<0] <- 1e-4; r_hat[r_hat<0] <- 1e-4; y_hat[y_hat<1] <- 1; y_hat[y_hat>5] <- 5

# evaluate estimated fuzzy numbers using Grzegorzewski and Roaniuk's metrics
B <- 1e2
BD_fgibbs <- matrix(nrow = n,ncol = B); bd_fols <- matrix(nrow = n,ncol = 1)
FZ_fgibbs <- matrix(nrow = n,ncol = B); fz_fols <- matrix(nrow = n,ncol = 1)
xsup <- seq(1,5,length.out=1001)
for(i in 1:n){
  cat("\t",i)
  
  # convert estimated beta fuzzy numbers into trg fuzzy numbers
  U <- t(mapply(function(b)unlist(beta2trg(m = M_hat[i,b],s = S_hat[i,b])),1:B)) #lb,ub,m
  Uresc <- rescale_inv(U,1,5)
  
  # dist Bertol
  BD_fgibbs[i,] <- mapply(function(b){
    s1 <- cbind(Uresc[b,1],Uresc[b,3],Uresc[b,3],Uresc[b,2]); s2 <- cbind(Ytrg[i,1],Ytrg[i,3],Ytrg[i,3],Ytrg[i,2])
    FuzzyResampling::BertoluzzaDistance(s1,s2)
  },1:B)
  s1 <- cbind(y_hat[i]-l_hat[i],y_hat[i],y_hat[i],y_hat[i]+l_hat[i]); s2 <- cbind(Ytrg[i,1],Ytrg[i,3],Ytrg[i,3],Ytrg[i,2])
  bd_fols[i] <- FuzzyResampling::BertoluzzaDistance(s1,s2)
  
  # Fuzz
  S1 <- cbind(Uresc[,1],Uresc[,3],Uresc[,3],Uresc[,2])
  FZ_fgibbs[i,] <- FuzzyResampling::CalculateFuzziness(S1)/FuzzyResampling::CalculateFuzziness(s2)
  fz_fols[i] <- FuzzyResampling::CalculateFuzziness(s1)/FuzzyResampling::CalculateFuzziness(s2)
}

## Figure S.22 (Supplementary materials)
#x11(); 
tikzDevice::tikz(file='casestudy_add_fig2.tex',width=7,height=8,sanitize = TRUE)
par(mfcol=c(2,1))
par(mar = c(2, 2, 3, 2))  # Smaller margins for each plot (bottom, left, top, and right)
par(oma = c(4, 1, 1, 1))  # Small outer margins

boxplot(BD_fgibbs,frame=FALSE,col=cls[1],outline=FALSE,border=FALSE,medcol="#CD6600",medlwd=6,axes=FALSE,ylim=c(0,quantile(BD_fgibbs,0.95)))
points(1:B,rep(quantile(bd_fols,0.50),B),lwd=4,col=cls[2],type="l"); points(1:B,rep(quantile(bd_fols,0.25),B),lwd=2,col="gray",type="l",lty=2);points(1:B,rep(quantile(bd_fols,0.75),B),lwd=2,col="gray",type="l",lty=2)    
axis(side = 2,at = round(seq(0,quantile(BD_fgibbs,0.95),length.out=11),2))
title(main = "(A) Bertoluzza distance",line = 1,adj=0)

boxplot(1-FZ_fgibbs,frame=FALSE,col=cls[1],outline=FALSE,border=FALSE,medcol="#CD6600",medlwd=6,axes=FALSE,ylim=c(quantile(1-fz_fols,0.05),quantile(1-FZ_fgibbs,0.95)))
points(1:B,rep(quantile(1-fz_fols,0.50),B),lwd=4,col=cls[2],type="l"); points(1:B,rep(quantile(1-fz_fols,0.25),B),lwd=2,col="gray",type="l",lty=2);points(1:B,rep(quantile(1-fz_fols,0.75),B),lwd=2,col="gray",type="l",lty=2)    
axis(side = 2,at = round(seq(quantile(1-fz_fols,0.05),quantile(1-FZ_fgibbs,0.95),length.out=11),2))
abline(h = 0.0,lty=1,lwd=3,col="#696969")
title(main = "(B) Fuzziness ratio",line = 1,adj=0)

add_legend("bottom",fill = cls,legend = c("fCM","PLR"),border = FALSE,bty = "n",ncol = 2,cex=1.25)
dev.off()

U2 <- cbind(summary(apply(BD_fgibbs,1,mean)),summary(apply(bd_fols,1,mean)),
           summary(apply(1-FZ_fgibbs,1,mean)),summary(apply(1-fz_fols,1,mean)))



### Nasrabadi's approach (Multi-objective fuzzy linear regression)
out_moflr <- fuzzyreg::moflr(x = cbind(X,0,0,0,0,0),y = cbind(yy,lr))
beta_moflr <- out_moflr$coef[,1:2]

## estimate var(beta_est) and CI(beta) via standard bootstrap
set.seed(231)
B=1000
Beta_boot_c = matrix(data = NA,nrow = B,ncol = nrow(beta_plr))
Beta_boot_l = matrix(data = NA,nrow = B,ncol = nrow(beta_plr))
for(b in 1:B){
  iid = sample(x = 1:n,size = n,replace = TRUE)
  Y_boot <- cbind(yy,lr)[iid,]
  out <- fuzzyreg::moflr(x = cbind(X,0,0,0,0,0),y = Y_boot)
  Beta_boot_c[b,] <- out$coef[,1]
  Beta_boot_l[b,] <- out$coef[,2]
  if(b%%10 == 0){print(b)}
}

beta_sd_c = apply(Beta_boot_c,2,sd) #approximated
beta_sd_l = apply(Beta_boot_l,2,sd) #approximated
Beta_ci_c = apply(Beta_boot_c,2,function(x)quantile(x,probs = c(0.05/2, 1 - 0.05/2), type = 1)) #approximated
Beta_ci_l = apply(Beta_boot_l,2,function(x)quantile(x,probs = c(0.05/2, 1 - 0.05/2), type = 1)) #approximated
Beta_moflr_tab <- cbind(beta_moflr[,1],beta_sd_c,t(Beta_ci_c),
                      beta_moflr[,2],beta_sd_l,t(Beta_ci_l));


## evaluate the goodness-of-fit
y_hat <- X%*%beta_moflr[,1]
l_hat <- X%*%beta_moflr[,2]
r_hat <- l_hat

constraints_violation_3 <- c(sum(y_hat<1),sum(y_hat>5),sum(l_hat<0),sum(r_hat<0))
#l_hat[l_hat<0] <- 1e-4; r_hat[r_hat<0] <- 1e-4; y_hat[y_hat<1] <- 1; y_hat[y_hat>5] <- 5

# evaluate estimated fuzzy numbers using Grzegorzewski and Roaniuk's metrics
B <- 1e2
BD_fgibbs <- matrix(nrow = n,ncol = B); bd_fols <- matrix(nrow = n,ncol = 1)
FZ_fgibbs <- matrix(nrow = n,ncol = B); fz_fols <- matrix(nrow = n,ncol = 1)
xsup <- seq(1,5,length.out=1001)
for(i in 1:n){
  cat("\t",i)
  
  # convert estimated beta fuzzy numbers into trg fuzzy numbers
  U <- t(mapply(function(b)unlist(beta2trg(m = M_hat[i,b],s = S_hat[i,b])),1:B)) #lb,ub,m
  Uresc <- rescale_inv(U,1,5)
  
  # dist Bertol
  BD_fgibbs[i,] <- mapply(function(b){
    s1 <- cbind(Uresc[b,1],Uresc[b,3],Uresc[b,3],Uresc[b,2]); s2 <- cbind(Ytrg[i,1],Ytrg[i,3],Ytrg[i,3],Ytrg[i,2])
    FuzzyResampling::BertoluzzaDistance(s1,s2)
  },1:B)
  s1 <- cbind(y_hat[i]-l_hat[i],y_hat[i],y_hat[i],y_hat[i]+l_hat[i]); s2 <- cbind(Ytrg[i,1],Ytrg[i,3],Ytrg[i,3],Ytrg[i,2])
  bd_fols[i] <- FuzzyResampling::BertoluzzaDistance(s1,s2)
  
  # Fuzz
  S1 <- cbind(Uresc[,1],Uresc[,3],Uresc[,3],Uresc[,2])
  FZ_fgibbs[i,] <- FuzzyResampling::CalculateFuzziness(S1)/FuzzyResampling::CalculateFuzziness(s2)
  fz_fols[i] <- FuzzyResampling::CalculateFuzziness(s1)/FuzzyResampling::CalculateFuzziness(s2)
}

## Figure S.23 (Supplementary materials)
#x11(); 
tikzDevice::tikz(file='casestudy_add_fig3.tex',width=7,height=8,sanitize = TRUE)
par(mfcol=c(2,1))
par(mar = c(2, 2, 3, 2))  # Smaller margins for each plot (bottom, left, top, and right)
par(oma = c(4, 1, 1, 1))  # Small outer margins

boxplot(BD_fgibbs,frame=FALSE,col=cls[1],outline=FALSE,border=FALSE,medcol="#CD6600",medlwd=6,axes=FALSE,ylim=c(0,quantile(BD_fgibbs,0.95)))
points(1:B,rep(quantile(bd_fols,0.50),B),lwd=4,col=cls[2],type="l"); points(1:B,rep(quantile(bd_fols,0.25),B),lwd=2,col="gray",type="l",lty=2);points(1:B,rep(quantile(bd_fols,0.75),B),lwd=2,col="gray",type="l",lty=2)    
axis(side = 2,at = round(seq(0,quantile(BD_fgibbs,0.95),length.out=11),2))
title(main = "(A) Bertoluzza distance",line = 1,adj=0)

boxplot(1-FZ_fgibbs,frame=FALSE,col=cls[1],outline=FALSE,border=FALSE,medcol="#CD6600",medlwd=6,axes=FALSE,ylim=c(quantile(1-fz_fols,0.05),quantile(1-FZ_fgibbs,0.95)))
points(1:B,rep(quantile(1-fz_fols,0.50),B),lwd=4,col=cls[2],type="l"); points(1:B,rep(quantile(1-fz_fols,0.25),B),lwd=2,col="gray",type="l",lty=2);points(1:B,rep(quantile(1-fz_fols,0.75),B),lwd=2,col="gray",type="l",lty=2)    
axis(side = 2,at = round(seq(quantile(1-fz_fols,0.05),quantile(1-FZ_fgibbs,0.95),length.out=11),2))
abline(h = 0.0,lty=1,lwd=3,col="#696969")
title(main = "(B) Fuzziness ratio",line = 1,adj=0)

add_legend("bottom",fill = cls,legend = c("fCM","MOFLR"),border = FALSE,bty = "n",ncol = 2,cex=1.25)

dev.off()

U3 <- cbind(summary(apply(BD_fgibbs,1,mean)),summary(apply(bd_fols,1,mean)),
           summary(apply(1-FZ_fgibbs,1,mean)),summary(apply(1-fz_fols,1,mean)))


### Fuzzy rule-based regression (Riza et al, 2015)

data_train_y <- cbind(X[,-1],Ytrg[,3])
data_train_l <- cbind(X[,-1],Ytrg[,3]-Ytrg[,1])
data_train_r <- cbind(X[,-1],Ytrg[,2]-Ytrg[,3])
range_data <- apply(rbind(data_train_y, data_train_l, data_train_r), 2, range)

out_frbs_y <- frbs::frbs.learn(data.train = data_train_y, range.data = range_data, method.type = "FIR.DM")
out_frbs_l <- frbs::frbs.learn(data.train = data_train_l, range.data = range_data, method.type = "FIR.DM")
out_frbs_r <- frbs::frbs.learn(data.train = data_train_r, range.data = range_data, method.type = "FIR.DM")

#Note: frbs does not output coefficients as it infers the functional form of the linear relationship among y and X using fuzzy logical rules


## evaluate the goodness-of-fit
y_hat <- predict(out_frbs_y, X[,-1])
l_hat <- predict(out_frbs_l, X[,-1])
r_hat <- predict(out_frbs_r, X[,-1])

constraints_violation_4 <- c(sum(y_hat<1),sum(y_hat>5),sum(l_hat<0),sum(r_hat<0))
l_hat[l_hat<0] <- 1e-4; r_hat[r_hat<0] <- 1e-4; y_hat[y_hat<1] <- 1; y_hat[y_hat>5] <- 5

# evaluate estimated fuzzy numbers using Grzegorzewski and Roaniuk's metrics
B <- 1e2; n <- NROW(Ytrg)
BD_fgibbs <- matrix(nrow = n,ncol = B); bd_fols <- matrix(nrow = n,ncol = 1)
FZ_fgibbs <- matrix(nrow = n,ncol = B); fz_fols <- matrix(nrow = n,ncol = 1)
xsup <- seq(1,5,length.out=1001)
for(i in 1:n){
  cat("\t",i)
  
  # convert estimated beta fuzzy numbers into trg fuzzy numbers
  U <- t(mapply(function(b)unlist(beta2trg(m = M_hat[i,b],s = S_hat[i,b])),1:B)) #lb,ub,m
  Uresc <- rescale_inv(U,1,5)
  
  # dist Bertol
  BD_fgibbs[i,] <- mapply(function(b){
    s1 <- cbind(Uresc[b,1],Uresc[b,3],Uresc[b,3],Uresc[b,2]); s2 <- cbind(Ytrg[i,1],Ytrg[i,3],Ytrg[i,3],Ytrg[i,2])
    FuzzyResampling::BertoluzzaDistance(s1,s2)
  },1:B)
  s1 <- cbind(y_hat[i]-l_hat[i],y_hat[i],y_hat[i],y_hat[i]+l_hat[i]); s2 <- cbind(Ytrg[i,1],Ytrg[i,3],Ytrg[i,3],Ytrg[i,2])
  bd_fols[i] <- FuzzyResampling::BertoluzzaDistance(s1,s2)
  
  # Fuzz
  S1 <- cbind(Uresc[,1],Uresc[,3],Uresc[,3],Uresc[,2])
  FZ_fgibbs[i,] <- FuzzyResampling::CalculateFuzziness(S1)/FuzzyResampling::CalculateFuzziness(s2)
  fz_fols[i] <- FuzzyResampling::CalculateFuzziness(s1)/FuzzyResampling::CalculateFuzziness(s2)
}

## Figure S.24 (Supplementary materials)
#x11(); 
tikzDevice::tikz(file='casestudy_add_fig4.tex',width=7,height=8,sanitize = TRUE)
par(mfcol=c(2,1))
par(mar = c(2, 2, 3, 2))  # Smaller margins for each plot (bottom, left, top, and right)
par(oma = c(4, 1, 1, 1))  # Small outer margins

boxplot(BD_fgibbs,frame=FALSE,col=cls[1],outline=FALSE,border=FALSE,medcol="#CD6600",medlwd=6,axes=FALSE,ylim=c(0,quantile(bd_fols,0.95)))
points(1:B,rep(quantile(bd_fols,0.50),B),lwd=4,col=cls[2],type="l"); points(1:B,rep(quantile(bd_fols,0.25),B),lwd=2,col="gray",type="l",lty=2);points(1:B,rep(quantile(bd_fols,0.75),B),lwd=2,col="gray",type="l",lty=2)    
axis(side = 2,at = round(seq(0,quantile(bd_fols,0.95),length.out=11),2))
title(main = "(A) Bertoluzza distance",line = 1,adj=0)

boxplot(1-FZ_fgibbs,frame=FALSE,col=cls[1],outline=FALSE,border=FALSE,medcol="#CD6600",medlwd=6,axes=FALSE,ylim=c(-1.45,quantile(1-fz_fols,0.95)))
points(1:B,rep(quantile(1-fz_fols,0.50),B),lwd=4,col=cls[2],type="l"); points(1:B,rep(quantile(1-fz_fols,0.25),B),lwd=2,col="gray",type="l",lty=2);points(1:B,rep(quantile(1-fz_fols,0.75),B),lwd=2,col="gray",type="l",lty=2)    
axis(side = 2,at = round(seq(-1.45,quantile(1-fz_fols,0.95),length.out=11),2))
abline(h = 0.0,lty=1,lwd=3,col="#696969")
title(main = "(B) Fuzziness ratio",line = 1,adj=0)

add_legend("bottom",fill = cls,legend = c("fCM","FRBS"),border = FALSE,bty = "n",ncol = 2,cex=1.25)

dev.off()

U4 <- cbind(summary(apply(BD_fgibbs,1,mean)),summary(apply(bd_fols,1,mean)),
           summary(apply(1-FZ_fgibbs,1,mean)),summary(apply(1-fz_fols,1,mean)))


## Tables

# Table S.1 (Supplementary materials)
B <- rbind(A[,1:4],
           Beta_fols_tab[c(1:4,6,5,7:8),],
           rbind(Beta_plr_tab[c(1:4,6,5),c(1:4)],Beta_plr_tab[c(1:4,6,5),c(5:8)]),
           rbind(Beta_moflr_tab[c(1:4,6,5),c(1:4)],Beta_moflr_tab[c(1:4,6,5),c(5:8)])
)

Xtab_tex = xtable::xtable(B)
attributes(Xtab_tex)$caption = ""
attributes(Xtab_tex)$label = ""
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})

# Table S.2 (Supplementary materials)
B <- rbind(t(cbind(U1[,1],U1[,2],U2[,2],U3[,2],U4[,2])),
           t(cbind(U1[,3],U1[,4],U2[,4],U3[,4],U4[,4])) )

Xtab_tex = xtable::xtable(B)
attributes(Xtab_tex)$caption = ""
attributes(Xtab_tex)$label = ""
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})

# Table S.3 (Supplementary materials)

B <- rbind(constraints_violation_1,
           constraints_violation_2,
           constraints_violation_3,
           constraints_violation_4)

Xtab_tex = xtable::xtable(B)
attributes(Xtab_tex)$caption = ""
attributes(Xtab_tex)$label = ""
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})



# Assessing model's mispecifications (Supplementary materials) ------------

### Case fy(y;theta_y) = Normal
rm(list=ls()); source("utils.R"); source("utilities/utilities.R")
graphics.off()

## 1. Data generation
set.seed(90210)
n <- 2.5e2; p <- 3; q <- 0
X <- cbind(1,matrix(runif(n*p,min = -1,max = 1),n,p))
X[,4] <- model.matrix(~sample(c("A","B"),size = n,replace = TRUE))[,-1]
X <- model.matrix(~X[,-1]); colnames(X)=paste0("X",1:NCOL(X))
Z <- cbind(1,scale(matrix(runif(n*q,min = -10,max = 10),n,q)))
#Z <- model.matrix(~sample(c("A","B"),size = n,replace = TRUE))
bx <- c(runif(1+p,-5,5)); gx <- c(runif(1+q,0.5,6.5)); c(bx,exp(gx))

# 1a. Generation of the true unobserved signal (0-level)
mu <- X%*%bx;  phi <- exp(Z%*%gx)
y0 <- mapply(function(i)rnorm(1,mu[i],phi[i]),1:n)
lb0 <- min(y0)-rnorm(1,sd = 0.1); ub0 <- max(y0)+rnorm(1,sd = 0.1) #true lower/upper bnds

# 1b. Generation of the shifted unobserved signal (1-level)
d0 <- quantile(y0,0.94)
y1 <- mapply(function(i)truncnorm::rtruncnorm(n = 1,a = lb0,mean = y0[i]+d0,b = ub0,sd = 1.5),1:n)
lb1 <- min(y1)-0.01; ub1 <- max(y1)+0.01 #true lower/upper bnds

# 1c. Generation of the (independent) observed spread 
s <- rgamma(n,25.0,25.0/50.0)

# 1d. Generation of the observed modes (which depends upon the 1-level blurred true signal)
m0 <- rescale_inv(mapply(function(i)rbeta(1,rescale(y0[i],lb0,ub0)*s[i],s[i]-s[i]*rescale(y0[i],lb0,ub0)),1:n),lb0,ub0) #no d-shifted observed fuzzy data
m1 <- rescale_inv(mapply(function(i)rbeta(1,rescale(y1[i],lb1,ub1)*s[i],s[i]-s[i]*rescale(y1[i],lb1,ub1)),1:n),lb1,ub1) #d-shifted observed fuzzy data

cbind(summary(y0),summary(y1),summary(m1))
boxplot(cbind(y0,y1,m1)); abline(h = mean(y0))
# Note: This is a hierarchical model. Actually, the true unobserved random vector Y is obscured by Y1, which is itself blurred by S.
# If not properly represented into the model, our model recovers Y1 instead of Y.

# 1e. Computing support of the observed fuzzy data
alphax=1e-9; xsup=seq(0,1,length.out = 1e3)
Yobs_a0 <- t(mapply(function(i){
  mi_current <- rescale(m1[i],min(m1),max(m1))
  fy <- beta_fn(x = xsup,mi = mi_current,gi = s[i])
  xout <- c(min(xsup[fy>alphax]),max(xsup[fy>alphax]))
  return(xout)
},1:n))
Yobs_a0_rescaled <- apply(Yobs_a0,2,rescale_inv,min(m1),max(m1))

# 1f. Evaluate 0-level signal wrt observed fuzzy data
#  how many fuzzy observations are above/below the true unobserved rv y0 ?
(1-sum(y0>Yobs_a0_rescaled[,1] & y0<Yobs_a0_rescaled[,2])/n) * 100

# 1g. Proxy for the bounds
lb_est <- min(Yobs_a0_rescaled[,1]); ub_est <- max(Yobs_a0_rescaled[,2]) #proxies for the true bnds based on 0-alpha cuts
lb_est <- lb_est-0.1; ub_est <- ub_est+0.1

## 2. Model building
fY <- function(y,X,Z,lb,ub,pars){
  b <- pars[1:ncol(X)]
  g <- pars[(1+ncol(X)):(ncol(X)+ncol(Z))]
  eta <- X%*%b
  theta <- eta
  phi <- exp(Z%*%g)
  (1 / (phi * sqrt(2 * pi))) * exp(-((y - theta)^2) / (2 * phi^2))
}
db_internal() 

# 2.1 Preparing for c++ conversion of R functions
filehead <- "utilities/R2Cpp_internal.head"; filetail <- "utilities/R2Cpp_internal.tail"; fileout <- "utilities/R2Cpp_internal.cpp"
verboseCpp <- TRUE
source("utilities/R2Cpp_utilities.R")
modelName <- "casestudy_addMat_Normal_0"
cacheDir_rcpp <-  paste0(getwd(),"/cache/",modelName)

dm_internal(fY,NCOL(X),NCOL(Z)) #create internal R functions
source("utilities/R2Cpp_funConverter_exec.R") #convert to c++

dir.create(path = paste0(getwd(),"/cache/casestudy_addMat_Normal_1"))
system(command = paste0("scp -r ",paste0("cache/casestudy_addMat_Normal_0/*"), " ",paste0("cache/casestudy_addMat_Normal_1")))

save(fY,fY_internal,tud_lnFy_exec,tud_lnFy_internal,grad_fun,hess_fun,file = paste0("cache/","casestudy_addMat_Normal_0","/routines.rds")) #for not d-shifted
save(fY,fY_internal,tud_lnFy_exec,tud_lnFy_internal,grad_fun,hess_fun,file = paste0("cache/","casestudy_addMat_Normal_1","/routines.rds")) #for d-shifted

## 3. Defining priors
J <- NCOL(X); H <- 1
mu_theta <- c(rep(0,J),rep(0,H)); sigma_theta <- c(rep(3.5,J),rep(3.5,H))

## 4. Saving data for GNU Parallel
# d-shifted fuzzy data
datain <- list(X=X,Z=Z,m=m1,s=s,lb=lb_est,ub=ub_est,mu_theta=mu_theta,sigma_theta=sigma_theta,y=y1) 
save(datain,file = paste0("data/datain_",gsub(x = modelName,pattern = "_0",replacement = "_1"),".rds")) 

# no d-shifted fuzzy data
datain <- list(X=X,Z=Z,m=m0,s=s,lb=lb0-0.1,ub=ub0+0.1,mu_theta=mu_theta,sigma_theta=sigma_theta,y=y0)
save(datain,file = paste0("data/datain_",modelName,".rds"))

## 5. Running multiple chains via GNU Parallel
currentModel <- "casestudy_addMat_Normal_1"  #whether to use the data for d-shifted (1) or not d-shifted (0) model
nchains <- 10
B <- 2e3 #samples

seedx <- as.numeric(format(Sys.time(), format = "%d%M%S"))
cmds <- sapply(1:nchains,function(u){
  cmd <- paste0(
    "R CMD BATCH --no-restore --no-save '--args B=", B, 
    " chainId=", u,
    " seedx=", seedx,
    " modelname=\"", currentModel, "\"", 
    "' AGS_exec.R ", 
    "cache/",currentModel,"/out_chainId_",u,".stdout"
  )
})

ncores <- nchains #ncores=nchains
btcs <- split(cmds, ceiling(seq_along(cmds)/ncores)) 
for(jobs in btcs){ 
  cat("\n Running:");for(x in jobs){cat("\n\t",strsplit(x,split = "--args")[[1]][2])}
  
  batch_file <- tempfile("job_list_"); writeLines(jobs, batch_file)
  system(paste0("parallel -j ", nchains, " < ", batch_file)); unlink(batch_file)
}


## 6. Save final mcmc as array (B x Nchains x Params)
currentModel <- "casestudy_addMat_Normal_1"  #whether to use the data for d-shifted (1) or not d-shifted (0) model
propBurIn <- 0.5
load(paste0("cache/",currentModel,"/routines.rds"))

xfls <- list.files(path = "results/",pattern = paste0("dataout_",currentModel),full.names = TRUE)
c(length(xfls),nchains) #converged vs run

burnin <- floor(B*propBurIn)
out <- lapply(xfls,readRDS)
out <- lapply(out,function(x)coda::as.mcmc(x[(burnin+1):B,]))
out_list <- coda::as.mcmc.list(out)
out_array <- posterior::as_draws_array(out_list)
out <- matrix(out_array,nrow = dim(out_array)[1]*dim(out_array)[2],ncol = dim(out_array)[3])
str(out_array)

posterior_table <- summarise_draws(posterior::as_draws_array(out))
hdpis <- coda::HPDinterval(obj = coda::as.mcmc(out))
saveRDS(list(mcmc=out_array,stats=list(ptable=posterior_table,HDPIs=hdpis,waic=NULL,waic2=NULL,waic_se=NULL,dic=NULL)),file = paste0("results/mcmc_",currentModel,".rds"))


## 7. Posterior analysis
# currentModel <- "casestudy_addMat_Normal_0" 
# mod0 <- readRDS(paste0("results/mcmc_",currentModel,".rds"))

currentModel <- "casestudy_addMat_Normal_1"  
mod1 <- readRDS(paste0("results/mcmc_",currentModel,".rds"))

## Posterior analysis
# A <- as.matrix(mod0$stats$ptable)[,-1]; A <- apply(A,2,as.numeric)
# A <- round(cbind(A[,c(1,3)],mod0$stats$HDPIs,A[,c(8:9)]),3)
A <- NULL

B <- as.matrix(mod1$stats$ptable)[,-1]; B <- apply(B,2,as.numeric)
B <- round(cbind(B[,c(1,3)],mod1$stats$HDPIs,B[,c(8:9)]),3)

# Table 4 (Supplementary materials)
U <- cbind(c(bx,gx),A[,1:4],B[,1:4])
#colnames(U) = c("true",rep(c("mean","sd","95% HDI lb","95% HDI ub"),2))
colnames(U) = c("true",rep(c("mean","sd","95% HDI lb","95% HDI ub"),1))
Xtab_tex = xtable::xtable(U)
attributes(Xtab_tex)$caption = ""
attributes(Xtab_tex)$label = ""
attributes(Xtab_tex)$align = rep("c",length(attributes(Xtab_tex)$align))
xtable::print.xtable(Xtab_tex,table.placement = "h!",sanitize.text.function = function(x){x})


# Figure S.25 (Supplementary materials)
load("data/datain_casestudy_addMat_Normal_1.rds")
#seedx <- format(Sys.time(), format = "%g%m%y"); set.seed(as.numeric(seedx))

s <- datain$s; m <- datain$m; n <- nrow(datain$X)
lb <- datain$lb; ub <- datain$ub

res <- optim(par = c(1,1),fn = function(x)-sum(dgamma(s,x[1],scale = x[2]/x[1],log = TRUE)),method = "L-BFGS-B",lower = c(0.01,0.01),upper = c(Inf,Inf))
alpha_s <- res$par[1]
beta_s <- res$par[2]/res$par[1] #the scale parameter is large enough!

pars_hat <- as.numeric(unlist(mod1$stats$ptable[,2]))
muy <- datain$X%*%pars_hat[1:4]
phiy <- exp(datain$Z%*%pars_hat[5]) 

B <- 5e2
xsup = seq(1e-2,1-1e-2,length=1e4+1)
Kfi1 <- matrix(NA,n,B); kfi1 = matrix(NA,1,n); Out_bxp1 = matrix(NA,n,5)
Kfi2 <- matrix(NA,n,B); kfi2 = matrix(NA,1,n); Out_bxp2 = matrix(NA,n,5)
Kfi3 <- matrix(NA,n,B); kfi3 = matrix(NA,1,n); Out_bxp3 = matrix(NA,n,5)
for(i in 1:n){
  s_hat <- rgamma(B,shape = alpha_s,scale = beta_s)
  y_hat <- mapply(function(b)rnorm(1,muy[i],phiy[i]),1:B)
  
  lb <- min(y_hat)-0.01; ub <- max(y_hat)+0.01
  m_hat01 <- mapply(function(b)rbeta(1,rescale(y_hat[b],lb,ub)*s_hat[b],s_hat[b]-s_hat[b]*rescale(y_hat[b],lb,ub)),1:B)
  m_hat <- rescale_inv(m_hat01,lb,ub)
  
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

#x11(); 
tikzDevice::tikz(file='violations_fig1.tex',width=7,height=3,sanitize = TRUE)
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

## Figure S.26 (Supplementary materials)
b_est <- pars_hat[1:4]; phi_est <- pars_hat[5]
Y_pred <- t(mvtnorm::rmvnorm(n = 5e2,mean = X%*%b_est,sigma = diag(n)*phi_est))

#x11()
tikzDevice::tikz(file='violations_fig2.tex',width=16,height=8,sanitize = TRUE)
par(mfrow=c(1,2))
boxplot(t(Y_pred),frame=FALSE,col="lightgray",outline=FALSE,border=FALSE,medcol="#CD6600",medlwd=6,axes=FALSE,ylim=c(quantile(as.numeric(y0),0.25),quantile(as.numeric(y0),0.99)))
axis(side = 2,at = round(seq(quantile(as.numeric(y0),0.05),quantile(as.numeric(y0),0.99),length.out=11),2))
points(y0,col="#698B69",pch=15); abline(h = mean(y0),col="#698B69",lwd=1.25,lty=2)
title(main = "(A)",line=1,adj=0,cex.main=2)

plot(density(Y_pred[,1],bw = 0.8),col="lightgray",bty="n",main="",xlab="",ylab="",ylim=c(0,0.2),xlim=c(min(y0),max(y0)))
set.seed(112); bbd <- sample(1:B,100,replace = FALSE)
for(b in bbd){lines(density(Y_pred[,b],bw = 0.8),col="lightgray")}
lines(density(y0,bw = 3),col="#698B69",lwd=4)
title(main = "(B)",line=1,adj=0,cex.main=2)

add_legend("bottomleft",fill = c("lightgray","#698B69"),legend = c("predicted y","true y"),border = FALSE,bty = "n",ncol = 1,cex=2)
dev.off()
