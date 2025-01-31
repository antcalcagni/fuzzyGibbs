AGS_alg <- function(B=2.5e3,J,H,X,Z,X_list,Z_list,m,s,lb,ub,parx_init,mu.theta,sigma.2.theta,sigma0_type=NA,sigma0delta=1.0,print.out=TRUE,NR_backtrack=0){
  # Note: sigma0delta > 0 works only if sigma0_type = -2
  
  parx <- matrix(NA,nrow=B,ncol = J+H); parx[1,] <- parx_init
  Inv_Diag_sigma_theta <- solve(diag(sigma.2.theta))
  
  for(b in 2:B){
    # Sampling from pi_y|theta
    out <- sampling_y_theta_loop(m,s,X_list,Z_list,parx[b-1,],lb = lb,ub = ub,sigma0 = sigma0_type,sigma0delta = sigma0delta,print_out = FALSE,maxiter = 1e3,eps = 1e-3); 
    ystar <- unlist(lapply(out,function(x)x[1]))
    
    # Sampling from pitheta|y
    #startx <- parx_init;y <- ystar; backtrack <- 0
    out <- NR_routine(parx_init,ystar,X,Z,X_list,Z_list,lb,ub,mu.theta,sigma.2.theta,stabil = TRUE,backtrack = NR_backtrack)
    
    if(out$conv==0){
      JJ <- solve(-hess_fun(ystar,X_list,Z_list,lb,ub,out$par) + Inv_Diag_sigma_theta)
      tx <- tud_lnFy_exec(J,H,ystar,X,Z,lb,ub,out$par) #Third-order unmixed derivative of the log-posterior density function
      out_dm <- solve_matching_dm(m = out$par, J = solve(JJ), t = tx) #Run the matching algorithm
      if(out_dm$conv>0){
        parx[b,] <- as.numeric(sn::rmsn(n=1, as.numeric(out_dm$mu), matrix(as.numeric(out_dm$Sigma),J+H,J+H), as.numeric(out_dm$d)))
      }else{
        parx[b,] <- as.numeric(mvtnorm::rmvnorm(1,out_dm$mu,out_dm$Sigma))
      }
    }else{
      parx[b,] <- rep(NA,J+H)
    }
    if(b%%100 && print.out){cat("\rProgress: ", paste0(round(b/B*100), "%"))}
  }
  
  return(parx)
}

NR_routine <- function(startx, y, X, Z, X_list, Z_list,lb,ub, mu.theta, sigma.2.theta,stabil=FALSE,stepsize=1,lambda=0,gamma2=0,backtrack=0,armijo_c=1e-4,wolfe_c=0.9,iter.max=250,tol=1e-6,epsilon=1e-9) {
  
  ## Netwon-Raphson algorithm with step-size optimized via backtracking ##
  
  # backtrack = {0: gradient-based, 1: Armijo, 2: Armijo + Wolfe}
  # Armijo constant c: Represents the "sufficient decrease" threshold in the objective function.
  #     Typical Range: from 1e-4 to 1e-1
  #     Behavior: Smaller c: Stricter decrease requirement
  #               Larger  c: More permissive, potentially accepting poor steps
  # Wolfe curvature r: Controls the curvature condition, ensuring the step size balances stability and sufficient descent.
  #     Typical Range: from 0.1 to 0.9
  #     Behavior: Smaller r: Easier to satisfy, might allow overshooting the optimum.
  #               Larger  c: Stricter condition, better convergence control but might reject good steps.
  
  # Initialize variables
  conv <- 0
  Diag_sigma_theta <- diag(sigma.2.theta)
  Inv_Diag_sigma_theta <- solve(diag(sigma.2.theta))
  theta <- startx
  D_prev <- rep(Inf, length(startx))
  obj_prev <- sum(log(fY(y, X, Z,lb,ub, theta))) + sum(mvtnorm::dmvnorm(as.numeric(theta), mu.theta, Diag_sigma_theta, log = TRUE))
  
  for (iter in seq_len(iter.max)) {
    
    # Compute gradient and Hessian
    #grad_eval <- grad_lnFy_exec(J, H, y, X, Z,lb,ub, theta) - solve(diag(sigma.2.theta)) %*% (theta - matrix(mu.theta, ncol = 1))
    grad_eval <- grad_fun(y,X_list,Z_list,lb,ub,theta) - Inv_Diag_sigma_theta %*% (theta - matrix(mu.theta, ncol = 1))
    #hess_eval <- -hessian_lnFy_exec(IJ, K, y, X, Z,lb,ub, theta) + diag(sigma.2.theta)
    hess_eval <- -hess_fun(y,X_list,Z_list,lb,ub,theta) + Diag_sigma_theta
    
    # Stabilization (if enabled)
    if (stabil) {
      sigma <- max(lambda, t(grad_eval) %*% grad_eval)
      hess_eval <- hess_eval + gamma2 * sigma * diag(nrow = nrow(hess_eval))
    }
    
    # Matrix inversion or pseudo-inverse fallback
    iH <- tryCatch(solve(hess_eval), error = function(e) {
      svd_decomp <- svd(hess_eval)
      inv_d <- ifelse(svd_decomp$d > epsilon, 1 / svd_decomp$d, 0)
      svd_decomp$v %*% diag(inv_d) %*% t(svd_decomp$u)
    })
    
    # Compute Newton step
    Delta <- stepsize * iH %*% grad_eval
    Lambda <- 1 # Initial step size
    theta_new <- theta + Lambda * as.vector(Delta)
    
    # Backtracking logic for step size adjustment
    if (backtrack == 0) {
      # Case 1: Gradient-based backtracking
      while (mean(grad_eval^2) >= mean(D_prev^2)) {
        Lambda <- Lambda / 2
        theta_new <- theta + Lambda * as.vector(Delta)
        #grad_eval <- grad_lnFy_exec(J, H, y, X, Z,lb,ub, theta_new) - solve(diag(sigma.2.theta)) %*% (theta - matrix(mu.theta, ncol = 1))
        grad_eval <- grad_fun(y,X_list,Z_list,lb,ub,theta_new) - Inv_Diag_sigma_theta %*% (theta_new - matrix(mu.theta, ncol = 1))
        if (Lambda < 1e-4) break
      } 
      # Note: more complex backtracking is disabled for the moment
      # } else if (backtrack >= 1) { 
      #   # Case 2: Objective-based backtracking (e.g., Armijo and Wolfe conditions)
      #   grad_new <- grad_eval
      #   k <- 0
      #   max_backtrack <- 50
      #   while (TRUE) {
      #     # Compute objective value and gradient at new parameters
      #     obj_new <- sum(log(fY(y, X, Z, theta_new))) + sum(mvtnorm::dmvnorm(as.numeric(theta_new), mu.theta, Diag_sigma_theta, log = TRUE))
      #     #grad_eval <- grad_lnFy_exec(J, H, y, X, Z, theta_new) - solve(diag(sigma.2.theta)) %*% (theta - matrix(mu.theta, ncol = 1))
      #     grad_new <- grad_fun(y,X_list,Z_list,theta_new) - Inv_Diag_sigma_theta %*% (theta_new - matrix(mu.theta, ncol = 1))
      #     
      #     # Armijo condition
      #     armijo <- obj_new <= obj_prev - armijo_c * Lambda * sum(grad_eval * Delta)
      #     
      #     # Wolfe curvature condition
      #     if (backtrack > 1) {
      #       wolfe <- abs(sum(grad_new * Delta)) <= wolfe_c * abs(sum(grad_eval * Delta))
      #     } else {
      #       wolfe <- TRUE
      #     }
      #     
      #     # Check if both conditions are satisfied
      #     if (armijo && wolfe) break
      #     
      #     # Reduce step size
      #     Lambda <- Lambda / 2
      #     theta_new <- theta + Lambda * as.vector(Delta)
      #     k <- k+1
      #     
      #     # Break if step size becomes too small
      #     if (Lambda < 1e-4 || k > max_backtrack) {
      #       warning("Backtracking failed to satisfy conditions after maximum iterations.")
      #       break
      #     }
      #   }
      #  obj_prev <- obj_new
    }
    
    # Update parameters
    D_prev <- grad_eval
    theta <- theta_new
    
    # Check convergence
    if (sqrt(mean(D_prev^2)) < tol) break
  }
  
  # Warn if maximum iterations reached
  if (iter == iter.max) {
    conv <- -1
    #warning("NR2: Maximum iterations reached without convergence.")
  }
  
  # Output results
  list(
    par = as.vector(theta),
    iterations = iter,
    gradient = D_prev,
    iH = iH,
    objective = obj_prev,
    conv = conv
  )
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

solve_matching_dm <- function(m, J, t) {
  # Solves the system of derivative matching equations analytically
  # Matched quantities are as follows (input):
  # m: mode
  # J: the negative Hessian at the mode 
  # t: the third-order unmixed derivatives
  
  u <- cbrt(t)
  R <- as.numeric(t(u)%*%solve(J)%*%u)
  
  final.fun <- function(kappa) {
    ret.val <- kappa*sn::zeta(3, kappa)^(2/3)/sn::zeta(1, kappa) + sn::zeta(2, kappa)*R^2/(sn::zeta(3, kappa)^(2/3) + sn::zeta(2, kappa)*R) - R
    return(ret.val)
  }
  
  aux.fun <- function(kappa) {
    ret.val <- sn::zeta(3, kappa)^(2/3) + sn::zeta(2, kappa)*R
    return(ret.val)
  }
  
  discont <- tryCatch(uniroot(aux.fun, lower = -1, upper = 1, extendInt = "yes", tol = .Machine$double.eps)$root,error=function(e){0})
  kappa.final <- tryCatch(uniroot(final.fun, lower = discont + 0.0001, upper = discont + 1, extendInt = "upX", tol = .Machine$double.eps)$root, error=function(e){NA})
  
  if (!is.na(kappa.final)) {
    d.final <- u/as.numeric(sn::zeta(3, kappa.final)^(1/3))
    Sigma.final <- solve(J + sn::zeta(2, kappa.final)*d.final%*%t(d.final))
    Sigma.final[lower.tri(Sigma.final)] <- t(Sigma.final)[lower.tri(Sigma.final)] # Making Sigma.final symmetric
    mu.final <- m - sn::zeta(1, kappa.final)*Sigma.final%*%d.final
    conv <- 1
  }else{
    cat("solve_matching_tm: Error (NA in final.fun)\n Return the Laplace Approximation solutions..\n")
    d.final <- u
    Sigma.final <- solve(J)
    Sigma.final[lower.tri(Sigma.final)] <- t(Sigma.final)[lower.tri(Sigma.final)] # Making Sigma.final symmetric
    mu.final <- m 
    conv <- -1
  }
  
  return(list(mu = mu.final, 
              Sigma = Sigma.final,
              d = d.final,
              kappa = kappa.final,
              conv=conv))
}

cbrt <- function(x) {
  # Find the real cube root of a number (used in skew-normal matching method)
  return(sign(x)*abs(x)^(1/3))
}


vec2symMat <- function(vec, diag = TRUE) {
  # Determine the size of the symmetric matrix
  n <- if (diag) {
    (-1 + sqrt(1 + 8 * length(vec))) / 2  # Solves n(n + 1)/2 = length(vec)
  } else {
    (1 + sqrt(1 + 8 * length(vec))) / 2  # Solves n(n - 1)/2 = length(vec)
  }
  
  if (n != floor(n)) stop("The vector length is incompatible with a symmetric matrix.")
  n <- as.integer(n)
  
  # Initialize an empty matrix
  mat <- matrix(0, n, n)
  
  # Fill the lower triangular part
  if (diag) {
    mat[lower.tri(mat, diag = TRUE)] <- vec
  } else {
    mat[lower.tri(mat, diag = FALSE)] <- vec
  }
  
  # Copy the lower triangular part to the upper triangular part
  mat <- mat + t(mat)
  if (diag) diag(mat) <- diag(mat) / 2  # Avoid doubling the diagonal
  
  return(mat)
}

v2m <- function (vec,K=NULL) {
  ltm <- diag(K)
  ltm[lower.tri(ltm,diag = TRUE)] <- vec
  ltm[upper.tri(ltm, diag = FALSE)] <- t(ltm)[upper.tri(ltm, diag = FALSE)]
  return(ltm)
}
