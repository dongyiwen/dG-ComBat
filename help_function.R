est.multi <- function(n.batches,n.array,cov.data,penalty_para){
  
  J <- length(cov.data)
  p <- ncol(cov.data[[1]])
  
  
  site.weights <- n.batches/n.array
  Sigma.comb <- site.weights[1]*cov.data[[1]]
  for(j in 2:J){
    Sigma.comb <- Sigma.comb+site.weights[j]*cov.data[[j]]
  }
  

  # detection of network structure
  
  
  prctile_vec <- seq(30, 99, length.out = 20)
  lam_vec <- seq(0.4, 1.2, length.out = 20)
  
  
  
  param <- param_tuning_new(Sigma = Sigma.comb, prctile_vec, lam_vec, norm.type = "F", penalty = penalty_para)## 0.8 for pet
  results <- dense(W_original = Sigma.comb, threshold = param$cut_out, lambda = param$lambda_out)
  
  CID <- results$CID[results$CID>=9] #>10
  if(length(CID)==0){
    CID <- p
    Clist <- 1:p
  }else{
    if(sum(CID) < p){
      Clist <- c(results$Clist[1: sum(CID)],sort(results$Clist[(sum(CID)+1):p]))
      CID <- c(CID, p-sum(CID))
    }else{
      Clist <- results$Clist[1:sum(CID)]
    }
  }
  
  K.est <- length(CID)
  L <- matrix(0, p, K.est)
  indices <- rep(1:K.est, CID)[order(Clist)]
  L[cbind(1:p, indices)] <- 1
  
  orig.CID <- results$CID[results$CID>=9]
  
  ## add new code
  
  K= dim(L)[2]
  G= dim(L)[1]

  block_genes <- lapply(1:K, function(k) which(L[,k] != 0))
  CID_L <- sapply(block_genes, length)
  LtL        <- diag(CID_L )           # K x K
  LtL_inv    <- diag(1/CID_L )


  Sigma_fhat <- LtL_inv %*% t(L) %*% Sigma.comb %*% L %*% LtL_inv

  if( sum(orig.CID) != G ){
    Sigma_fhat[,K] <- rep(0)
    Sigma_fhat[K,] <- rep(0)
  }
  updated.Sigma.comb =  Sigma_fhat
  
  return(list(L=L,Sigma.comb = Sigma.comb,W_dense=results$W_dense,CID=orig.CID,updated.Sigma.comb=updated.Sigma.comb))
}

estimate_hyperparams <- function(n.arrays, mu_hat, L,central_out_CID,cov.data) {
  # Z: n x G data matrix
  # L: G x K loading matrix (block structure)
  n <- n.arrays
  G <- ncol(cov.data)
  K <- ncol(L)
  
  # 1) Per-gene means
  mu_hat <- mu_hat
  block_genes <- lapply(1:K, function(k) which(L[,k] != 0))
  CID <- sapply(block_genes, length)
  
  # 3) Compute block-wise mu_bar and tau2_bar
  mu_bar   <- sapply(block_genes, function(idx) mean(mu_hat[idx]))
  tau2_bar <- sapply(block_genes, function(idx) {
    m <- mu_hat[idx]
    var(m)})
  
  
  # 4) Sample covariance matrix
  V_tilde <- cov.data
  #V_tilde  <- crossprod(Z_tilde) / (n - 1)  # G x G
  
  # 6) Factor scores and their covariance
  LtL        <- diag(CID)           # K x K
  LtL_inv    <- diag(1/CID)
  # SIGMA_F <- cov(t(f_hat))
  Sigma_fhat <- LtL_inv %*% t(L) %*% V_tilde %*% L %*% LtL_inv  # K x K
  
  if( sum(central_out_CID) != G ){
    Sigma_fhat[,K] <- rep(0)
    Sigma_fhat[K,] <- rep(0)
    
  }
  
  
  # 7) Residual covariance and per-gene variances
  Sigma_hat  <- V_tilde - L %*% Sigma_fhat %*% t(L)  # G x G
  delta2_hat <- diag(Sigma_hat)                     # length G
  
  V_bar_wins  <- numeric(length(block_genes))
  S2_bar_wins <- numeric(length(block_genes))
  
  for (k in seq_along(block_genes)) {
    idx   <- block_genes[[k]]
    vals  <- delta2_hat[idx]
    
    # 1) compute the 95th percentile
    p95   <- quantile(vals, 0.95)
    
    # 2) cap any value > p95 down to p95
    vals_w <- pmin(vals, p95)
    
    # 3) compute mean & variance on the winsorized vector
    V_bar_wins[k]  <- mean(vals_w)
    S2_bar_wins[k] <- var(vals_w) 
  }
  
  # 9') Method‐of‐moments on the winsorized moments
  lambda_bar <- V_bar_wins^2 / S2_bar_wins + 2
  theta_bar  <- V_bar_wins^3 / S2_bar_wins + V_bar_wins
  
  # 10) Return
  list(
    mu_bar     = mu_bar,      # length K
    tau2_bar   = tau2_bar,    # length K
    lambda_bar = lambda_bar,  # length K
    theta_bar  = theta_bar,    # length K
    V_tilde    = V_tilde,
    Sigma_fhat = Sigma_fhat,    #SIGMA_F = SIGMA_F,
    Sigma_hat  = Sigma_hat,
    block_genes = block_genes,
    mu_hat = mu_hat,
    prev_delta= delta2_hat
  )
}

update_combat_params <- function(n.arrays,cov.data, L, hyper, central_out_CID, delta_prev = NULL) {
  # Z          : n x G data matrix
  # L          : G x K loading matrix
  # hyper      : list from estimate_hyperparams(Z, L)
  # contains mu_hat, mu_bar, tau2_bar, Sigma_fhat, Sigma_hat, block_genes
  #
  # Returns a list with
  #   gamma_star : length-G posterior means
  #   delta_star : length-G posterior means
  
  n <- n.arrays
  G <- ncol(cov.data)
  K <- ncol(L)
  
  psi2_vec <- diag(L %*% hyper$Sigma_fhat %*% t(L))
  
  # 2) Determine v2_g = psi2_g + delta_prev_g (or fallback to var(Z[,g]))
  if (is.null(delta_prev)) {
    # first pass: use the sample variance of each column of Z
    #v2_vec <- psi2_vec + hyper$prev_delta
    v2_vec <- diag(cov.data)
  } else {
    v2_vec <- psi2_vec + delta_prev
  }
  
  # 3) Update gamma^*_g for each block
  gamma_star <- numeric(G)
  for (k in seq_along(hyper$block_genes)) {
    idx       <- hyper$block_genes[[k]]
    mu_bar_k  <- hyper$mu_bar[k]
    tau2_k    <- hyper$tau2_bar[k]
    
    mu_hat_k  <- hyper$mu_hat[idx]
    v2_k      <- v2_vec[idx]
    
    gamma_star[idx] <- (n * tau2_k * mu_hat_k + v2_k * mu_bar_k) /
      (n * tau2_k + v2_k)
  }
  
  
  # 4) Placeholder for delta^2 update
  #      delta_star[idx] = E[delta^2 | data, gamma_star, ...]
  #    according to the formula,
  
  I <- diag(nrow(L))
  block_genes <- lapply(1:K, function(k) which(L[,k] != 0))
  CID <- sapply(block_genes, length)
  LtL        <- diag(CID)           # K x K
  LtL_inv    <- diag(1/CID)
  
  if(sum(central_out_CID) != G){
    idx_raw <- which(L[, ncol(L)] == 1)
    idx_proj <- which(L[, ncol(L)] != 1)
    
    H.org <- L %*% LtL_inv  %*% t(L)
    H.org[idx_raw, idx_raw] <- 0
    H=I-H.org
    
  }else{
    H <- I-L %*% LtL_inv  %*% t(L)
  }
  part1<- (n-1)*H %*% cov.data %*% H
  part2<-I%*%matrix(hyper$mu_hat, nrow = G, ncol = n)-t(matrix(gamma_star, n,G, byrow = TRUE))
  
  #test part 2 with <- hyper$mu_hat - gamma_star
  final= part1+part2%*%t(part2)
  ssq=diag(final)
  
  
  
  delta_star <- numeric(G)
  for (k in seq_along(hyper$block_genes)) {
    idx    <- hyper$block_genes[[k]]
    lam_k  <- hyper$lambda_bar[k]
    the_k  <- hyper$theta_bar[k]
    delta_star[idx] <- (the_k + 0.5 * ssq[idx]) / (n/2 + lam_k - 1)
  }
  list(
    gamma_star = gamma_star,
    delta_star = delta_star
  )
}

run_combat_iterative <- function(n.arrays, mu_hat, L,cov.data,central_out_CID, tol = 1e-7, max_iter = 30, verbose = TRUE) {
  # Initial hyperparameters and holders
  hyper_params <- estimate_hyperparams(n.arrays, mu_hat, L,central_out_CID,cov.data) 
  delta_prev   <- NULL
  gamma_old    <- NULL
  delta_old    <- NULL
  
  for (iter in seq_len(max_iter)) {
    # Update ComBat parameters
    upd        <- update_combat_params(n.arrays,cov.data, L, hyper_params, central_out_CID, delta_prev)
    gamma_new  <- upd$gamma_star
    delta_new  <- upd$delta_star
    
    # Compute convergence metric
    if (!is.null(gamma_old)) {
      d_gamma <- sum(abs(gamma_new - gamma_old))
      d_delta <- sum(abs(delta_new - delta_old))
    } else {
      d_gamma <- NA_real_
      d_delta <- NA_real_
    }
    
    # Print progress
    if (verbose) {
      message(sprintf("[Iter %d]: Δγ = %s, Δδ^2 = %s",
                      iter,
                      ifelse(is.na(d_gamma), "NA", format(d_gamma, digits=4)),
                      ifelse(is.na(d_delta), "NA", format(d_delta, digits=4))))
    }
    
    # Check convergence
    if (!is.na(d_gamma) && max(d_gamma, d_delta) < tol) {
      if (verbose) message(sprintf("Converged in %d iterations (tol=%g)", iter, tol))
      break
    }
    
    # Prepare next iteration
    gamma_old  <- gamma_new
    delta_old  <- delta_new
    delta_prev <- delta_new
  }
  
  # Return results
  list(
    gamma       = gamma_new,
    delta2      = delta_new,
    iterations  = iter,
    converged   = (!is.na(d_gamma) && max(d_gamma, d_delta) < tol)
  )
}