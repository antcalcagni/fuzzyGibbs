beta_fn <- function(x,mi,gi){
  a <- 1+gi*mi; b= 1+gi*(1-mi)
  C <- (((a-1)/(a+b-2))^(a-1) * (1-((a-1)/(a+b-2)))^(b-1))
  fx <- 1/C * x^(a-1) * (1-x)^(b-1)
  return(fx)
}

polygamma <- function(n,z){pracma::psi(n,z)} #polygamma function
power = function(x,a){x^a} #power function
rescale <- function(x,lb=0,ub=1){(x-lb)/(ub-lb)}   #linearly rescale into [lb,ub]
rescale_inv <- function(x,lb=0,ub=1){lb+(ub-lb)*x} #inverse of rescale

fyA <- function(m,y,s,lb=0,ub=1){ #h() function --see SMPS paper
  y <- rescale(y,lb,ub); m <- rescale(m,lb,ub)
  return(1/beta(y*s,s-s*y) * (m/(1-m))^(s*y))
}

lnfyA <- function(m,y,s,lb=0,ub=1){ #h() function --see SMPS paper
  y <- rescale(y,lb,ub); m <- rescale(m,lb,ub)
  return(log(1/beta(y*s,s-s*y) * (m/(1-m))^(s*y)))
}

dpi_y <- function(fY,y,m,s,X,Z,pars,lb=0,ub=1){ #conditional posterior y|thetay,..
  z <- integrate(function(x)fyA(m,x,s,lb,ub)*fY(x,X,Z,pars),lb+sign(lb)*1e-3,ub-sign(ub)*1e-3)$val
  fx <- (fyA(m,y,s,lb,ub)*fY(y,X,Z,pars))/z
  return(fx)
}

dpi_y_approx <- function(m,s,lambda0,sigma0,X,Z,pars,lb=0,ub=1,eps=1e-9,print.out=FALSE){ #derivative approximation of pi(y|thetay,..)
  lambda <- lambda0; sigma <- sigma0; y <- 1.99; k <- 0
  m=rescale(m,lb,ub)
  while((abs(y/lambda-1))>eps){
    k <- k+1
    y <- lambda
    K1 <- d1lnFy(rescale_inv(y,lb,ub),X,Z,pars) +                                   s*log(m/(1 - m)) + s*(-polygamma(0,s*y) + polygamma(0,s - s*y))
    K2 <- d2lnFy(rescale_inv(y,lb,ub),X,Z,pars) -                                    power(s,2)*(polygamma(1,s*y) + polygamma(1,s - s*y))
    lambda <- (1 + y*(-2 + K1 + sigma - K1*y))/sigma
    sigma <- (1 - (-1 + y)*y*(-2 + K2*(-1 + y)*y))/(lambda - 2*lambda*y + power(y,2))
    if(print.out){print(c(k,lambda,rescale_inv(lambda,lb,ub),sigma))}
    if(is.na(lambda)||is.na(sigma)){break}
  }
  if(lambda<0||lambda>1||sigma<0||is.na(lambda)||is.na(sigma)){
    return(c(NA,NA))
  }else{
    return(c(rescale_inv(lambda,lb,ub),sigma))
  }
}

conditional_sampling_y <- function(m,s,X,Z,pars,lb=0,ub=1,lambda0=NULL,sigma0=NULL,maxIter=150,print.out=FALSE,seedx=NULL){ #sampling from pi(y|thetay,..)
  if(is.null(lambda0)){lambda0 = rescale(m,lb,ub)} #starting point for lambda
  
  out <- rep(NA,2); k <- 0; cond <- 0; yhat <- NA
  while(cond<=0 && k<=maxIter){
    k <- k+1
    if(is.null(sigma0)){sigma0 <- runif(1,4.99,300.1)} #starting point for sigma
    out <- dpi_y_approx(m,s,lambda0=lambda0,sigma0=sigma0,X,Z,pars,lb,ub,print.out=print.out)
    cond <- 1-sum(is.na(out))
  }
  if(!is.null(seedx)){set.seed(seedx)}
  if(cond>0 && k<=maxIter){
    out[1] <- rescale(out[1],lb,ub)
    yhat   <- rescale_inv(rbeta(1,out[1]*out[2],out[2]-out[2]*out[1]),lb,ub)
  }else{
    warning("conditional_sampling_y: No solution found. Returning observed values!")
    yhat   <- rescale_inv(rbeta(1,m*s,s-m*s),lb,ub)
  }
  return(list(y=yhat,lambda=out[1],sigma=out[2],conv=cond,runs=k))
}

db_internal <- function(){ #internal functions for computing the 4p-Beta approximation of pi(y|theta,..) 
  out <- c()
  x <- strsplit(x = deparse(fY)[1],split = " ")[[1]]
  out <- c(out,paste0("lnFy <- function",paste(x[2:length(x)],collapse = ""),"{log(",paste(c("fY",x[2:length(x)]),collapse=""),")}"))
  out <- c(out,"d1lnFy <- Deriv::Deriv(lnFy,'y',nderiv=1)")
  out <- c(out,"d2lnFy <- Deriv::Deriv(lnFy,'y',nderiv=2)")
  
  writeLines(text = out,con = paste0(getwd(),"/dfY_internal.temp"))
  source(paste0(getwd(),"/dfY_internal.temp"))
  file.remove(paste0(getwd(),"/dfY_internal.temp"))
  
  return(invisible())
}

dm_internal <- function(fY_internal,J,H){ #internal functions for computing the Skew-Normal approximation of pi(theta|y,..)

  #fY_internal <- fY
  out <- deparse(fY_internal)
  
  ii <- grep(x = out,pattern = " <- pars")
  if(length(ii)>0){out[ii][1] <- paste(c("#",out[ii][1]),collapse=""); out[ii][2] <- paste(c("#",out[ii][2]),collapse="")}
  
  ii <- grep(x = out,pattern = "eta <- ")[1]
  x <- strsplit(x = out[ii],split = " ")[[1]]
  x[length(x)] <- paste0("c(",paste(paste0("b",1:J),collapse = ","),")")
  out[ii] <- paste(x,collapse = " ")
  
  ii <- grep(x = out,pattern = "phi <- ")[1]
  x <- strsplit(x = out[ii],split = " ")[[1]]
  x[length(x)] <- paste0("c(",paste(paste0("b",(J+1):(J+H)),collapse = ","),"))")
  out[ii] <- paste(x,collapse = " ")
  
  x <- strsplit(x = out[1],split = " ")[[1]]
  x[length(x)] <- paste0(paste(c(paste0("b",1:J),paste0("b",(J+1):(J+H))),collapse=","),")")
  out[1] <- paste(c("fY_internal <-",x),collapse = " ")
  
  # Define function lnFy_internal()
  out <- c(out,paste(c("lnFy_internal <- ",x,"{","log(fY_internal",paste(x[2:length(x)],collapse = ""),")}"),collapse = ""))
  
  # Define symbolic grad function
  out <- c(out,paste0("grad_lnFy_internal <- mapply(function(j)Deriv::Deriv(lnFy_internal,paste0('b',j)),1:",(J+H),",SIMPLIFY = FALSE)"))
  
  # Define executable grad function
  #out <- c(out,paste0("grad_lnFy_exec <- function(J,H,y,X,Z,pars){grad <- matrix(NA,J+H,1);for(j in 1:(J+H)){grad[j] <- sum(grad_lnFy_internal[[j]](y,X,Z,",paste(paste0("pars[",1:(J+H),"]"),collapse = ","),"))};return(grad)}"))
  out <- c(out,paste0("grad_lnFy_exec <- function(J,H,y,X,Z,lb,ub,pars){grad <- matrix(NA,J+H,1);for(j in 1:(J+H)){grad[j] <- sum(grad_lnFy_internal[[j]](y,X,Z,lb,ub,",paste(paste0("pars[",1:(J+H),"]"),collapse = ","),"))};return(grad)}")) 
  
  # Define symbolic hess function
  out <- c(out,paste0("IJ <- expand.grid(1:",(J+H),",1:",(J+H),"); IJ <- IJ[apply(IJ,1,function(x){x[1]>=x[2]}),]; K <- nrow(IJ)"))
  out <- c(out,paste0("hess_lnFy_internal <- mapply(function(k)Deriv::Deriv(Deriv::Deriv(lnFy_internal,paste0('b',IJ[k,1])),paste0('b',IJ[k,2])),1:K,SIMPLIFY = FALSE)"))

  # Define executable hess function
  #out <- c(out,paste0("hessian_lnFy_exec <- function(IJ,K,y,X,Z,pars){Hess <- matrix(NA,max(IJ),max(IJ));for(k in 1:K){Hess[IJ[k,1],IJ[k,2]] <- sum(hess_lnFy_internal[[k]](y,X,Z,",paste(paste0("pars[",1:(J+H),"]"),collapse = ","),"))}; Hess <- as.matrix(Matrix::forceSymmetric(Hess,uplo='L')); return(Hess)}"))
  out <- c(out,paste0("hessian_lnFy_exec <- function(IJ,K,y,X,Z,lb,ub,pars){Hess <- matrix(NA,max(IJ),max(IJ));for(k in 1:K){Hess[IJ[k,1],IJ[k,2]] <- sum(hess_lnFy_internal[[k]](y,X,Z,lb,ub,",paste(paste0("pars[",1:(J+H),"]"),collapse = ","),"))}; Hess <- as.matrix(Matrix::forceSymmetric(Hess,uplo='L')); return(Hess)}"))
  
  # Define symbolic third-order unmixed derivative (for SN-based dm schema)
  out <- c(out,paste0("tud_lnFy_internal <- mapply(function(j)Deriv::Deriv(lnFy_internal,paste0('b',j),nderiv=3),1:",(J+H),",SIMPLIFY = FALSE)"))
  
  # Define executable tud function
  out <- c(out,paste0("tud_lnFy_exec <- function(J,H,y,X,Z,lb,ub,pars){tud <- matrix(NA,(J+H),1);for(j in 1:(J+H)){tud[j] <- sum(tud_lnFy_internal[[j]](y,X,Z,lb,ub,",paste(paste0("pars[",1:(J+H),"]"),collapse = ","),"))};return(tud)}"))
  
  # Define grad_fun (for NR algorithm)
  out <- c(out,paste0("grad_fun <- function(y,X_list,Z_list,lb,ub,theta){",
                      "x <- grad_lnfy_exec_loop(y,X_list,Z_list,lb,ub,",paste(paste0("theta[",1:(J+H),"]"),collapse = ","),");",
                      "return(matrix(x,ncol=1))}"))
  
  # Define hess_fun (for NR algorithm)
  out <- c(out,paste0("hess_fun <- function(y,X_list,Z_list,lb,ub,theta){",
                      "x <- hess_lnfy_exec_loop(y,X_list,Z_list,lb,ub,",paste(paste0("theta[",1:(J+H),"]"),collapse = ","),");",
                      "return(v2m(x,length(theta)))}"))
  
  writeLines(text = out,con = paste0(getwd(),"/fY_internal.temp"))
  source(paste0(getwd(),"/fY_internal.temp"))
  file.remove(paste0(getwd(),"/fY_internal.temp"))
  
  return(invisible())
}

## Syntax to call the c++ functions to evaluate gradient and Hessian:
#   grad_to_call <- paste0("grad_lnfy_exec_loop", "(",paste("y","X_list","Z_list",paste(paste0("cme[",1:jh,"]"),collapse=","),sep = ","),")")
#   cg <- matrix(eval(parse(text = grad_to_call)),ncol = 1) - solve(diag(sigma.2.theta))%*%(cme-matrix(mu.theta,ncol=1)) 
#   hess_to_call <- paste0("hess_lnfy_exec_loop", "(",paste("y","X_list","Z_list",paste(paste0("cme[",1:jh,"]"),collapse=","),sep = ","),")")
#   cH <- vec2symMat(x=eval(parse(text = hess_to_call)),diag = TRUE) - diag(sigma.2.theta);
#
# Microbenchmark results show that native R functions are faster than compiled ones. Probably, this is due to the eval(parse(..)) call used to run the c++ functions.
#
#   Results (ms)
#   -- Hessian 
#      expr      min       lq     mean   median       uq    neval
# R   2.010016 2.112602 2.772321 2.219633 2.367367 173.4951 10000
# C++ 7.626380 7.800421 7.952061 7.854125 7.989010  14.5210 10000
#
#   -- Gradient
#       expr      min       lq     mean    median       uq     neval
# R    626.826  648.764  821.547  672.5255  710.492 152617.100 10000
# C++ 2916.641 3000.856 3067.846 3031.5325 3098.847   5700.362 10000 
#
# Note that to use c++ functions for gradient and Hessian evaluations, the raw data need to be converted as lists (instead of matrices):
#  X_list <- lapply(1:nrow(X), FUN = function(i) X[i,,drop=FALSE]) 
#  Z_list <- lapply(1:nrow(Z), FUN = function(i) Z[i,,drop=FALSE])
#

laplace_approx_0ld <- function(X,Z,y,J,H,IJ,K,mu.theta=NULL,sigma.2.theta=NULL, maxit = 250, tol = 1.0E-5, startx = NULL,lowerx = NULL,upperx = NULL) {
  ## Laplace approximation of the posterior distribution, using line-search augmented Newton-Raphson
  ## Note: J,IJ,K are not user-defined input, rather they are internal variables created by the call dm_internal()
  if(is.null(sigma.2.theta)){sigma.2.theta <- rep(1e4,J+H)}
  if(is.null(mu.theta)){mu.theta <- rep(0,J+H)}
  if(is.null(lowerx)){lowerx <- rep(-Inf,J+H)}; if(is.null(upperx)){upperx <- rep(Inf,J+H)}
  
  if(is.null(startx)){
    startx <- rep(0,J+H)
    cme <- tryCatch(optim(par = startx,fn = function(x){-sum(log(fY(y,X,Z,x)))},method = "L-BFGS-B",lower = lowerx,upper = upperx)$par,error=function(e){startx})
  }else{
    cme <- startx
  }
  
  mx <- NA
  cme <- matrix(cme,ncol = 1)
  conv <- -1
  for (i in 1:maxit) {
    #print(i)
    cg <- grad_lnFy_exec(J,H,y,X,Z,cme) - solve(diag(sigma.2.theta))%*%(cme-matrix(mu.theta,ncol=1))
    if (max(abs(cg)) < tol) {
      mx <- cme
      conv <- 1
      break
    }
    
    cH <- hessian_lnFy_exec(IJ,K,y,X,Z,cme) - diag(sigma.2.theta);
    delta <- solve(cH)%*%cg
    
    # Line search for optimizing lambda
    lambda <- 1
    f.old <- sum(log(fY(y,X,Z,cme))) + sum(mvtnorm::dmvnorm(as.numeric(cme),(mu.theta),diag(sigma.2.theta),log = TRUE))
    for (k in 1:10) {
      new.point <- cme - lambda*delta
      f.curr <- sum(log(fY(y,X,Z,new.point))) + sum(mvtnorm::dmvnorm(as.numeric(new.point),(mu.theta),diag(sigma.2.theta),log = TRUE))
      if (f.curr < f.old) {
        lambda <- 0.5*lambda
      } else {
        break
      }
    }
    cme <- cme - lambda*delta
  }
  
  if(i>maxit || conv<1){
    res <- list(mu=startx*NA,Sigma=diag(J+H)*NA)
    cat("laplace_approx: Laplace Approximation failed to converged!")
  }else{
    JJ <- -1*hessian_lnFy_exec(IJ,K,y,X,Z,mx) + solve(diag(sigma.2.theta))  #Negative Hessian at the mode (the last term is about the Gaussian prior)
    res <- list(mu = mx, Sigma = solve(JJ))
  }
  return(res)
}

dm_approx <- function(X, Z, y, J, H, IJ, K, mu.theta=NULL,sigma.2.theta=NULL,startx=NULL,lowerx=NULL,upperx=NULL,maxit=500,tol=1e-5) {
  ## Multivariate skew-normal approximation - derivative matching
  ## Note: J,IJ,K are not user-defined input, rather they are internal variables created by the call dm_internal()
  
  laplace.res <- laplace_approx(X, Z, y, J, H, IJ, K, mu.theta,sigma.2.theta,maxit,tol,startx,lowerx,upperx)
  laplace.mu <- laplace.res$mu
  laplace.Sigma <- laplace.res$Sigma
  
  tx <- tud_lnFy_exec(J,H,y,X,Z,laplace.mu) #Third-order unmixed derivative of the log-posterior density function
  matching.values <- solve_matching_dm(m = laplace.mu, J = solve(laplace.Sigma), t = tx) #Run the matching algorithm
  
  return(list(mu = matching.values$mu,
              Sigma = matching.values$Sigma,
              d = matching.values$d,
              conv = matching.values$conv))
}



auto.layout = function(n, layout=T){
  ### figure out how many rows
  sq = sqrt(n)
  rws = round(sq)
  
  #### if it's a perfect square, fill the matrix
  if (sqrt(n) == round(sqrt(n))){
    numbs = sort(rep(1:n, times=2))
    m = matrix(numbs, nrow=sq, byrow=T)
  } else {
    
    #### repeat twice the numbers that fit nicely
    topNum = trunc(n/rws)*rws
    numbs = sort(rep(1:topNum, times=2))
    if (topNum==n){
      m = matrix(numbs, nrow=rws, byrow=T)
    } else {
      #### get the rest figured out
      rem = n-topNum  ### remaining numbers
      rest = sort(rep((topNum+1):n, times=2))
      cols = (topNum/rws)*2
      rest = c(rep(0, times=(cols-length(rest))/2), rest, rep(0, times=(cols-length(rest))/2))
      m = matrix(c(numbs, rest), nrow=rws+1, byrow=T)
    }
  }
  
  if (layout){
    layout(m)
  } else {
    m
  }
}
