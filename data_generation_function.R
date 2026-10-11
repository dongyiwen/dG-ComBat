Spd.proj<-function(SigA.hat, eps=NULL){
  p=ncol(SigA.hat)
  if(is.null(eps)){
    eps<-5/p
  }
  feasible=1
  SigA.t<-SigA.hat
  if(min(eigen(SigA.t)$values) <=eps ){
    feasible=2
    SigA.t<-ADMM_proj(SigA.hat, epsilon=eps)$mat
  }
  SigA.t<-lavaSearch2:::symmetrize(SigA.t, update.upper = TRUE)
  list(mat=SigA.t, conv=feasible)
}


# ─── Sub‐function: generate one block’s data ──────────────────────────────────
# Inputs:
#   n        – number of subjects
#   latent   – numeric vector of length n (the latent factor for this block)
#   Gk       – number of voxels/genes in this block
#   batch_m  – batch‐effect mean for rnorm
#   batch_sd – batch‐effect sd   for rnorm
#   lambda   – shape parameter for rinvgamma
#   theta    – scale parameter for rinvgamma
#
# Returns:
#   an n×Gk matrix: batch + latent + noise

generate_block <- function(n, latent, Gk, batch_m, batch_sd, lambda, theta) {
  #  one batch‐effect per feature
  gamma_vec <- rnorm(Gk, mean = batch_m, sd = batch_sd)
  batch_mat <- matrix(gamma_vec, nrow = n, ncol = Gk, byrow = TRUE)
  
  #  noise variances from inverse‐gamma
  noise_var <- MCMCpack::rinvgamma(Gk, shape = lambda, scale = theta)
  noise_mat <- sapply(noise_var, function(v) rnorm(n, mean = 0, sd = sqrt(v)))
  
  #  latent effect broadcasted across features
  latent_mat <- matrix(latent, nrow = n, ncol = Gk, byrow = FALSE)
  
  # data matrix
  data_mat <-batch_mat + latent_mat + noise_mat
  
  list(
    data   = data_mat,
    gamma  = gamma_vec,
    delta2 = noise_var
  )
}

# ─── Main simulator: multiple blocks ─────────────────────────────────────────
# Inputs:
#   n         – number of subjects/rows
#   G_vec     – integer vector of length K: G_1, …, G_K
#   Sigma_f   – K×K covariance matrix for latent factors
#   batch_means – numeric vector length K: means for each block’s batch effect
#   batch_sds   – numeric vector length K: sds   for each block’s batch effect
#   lambdas     – numeric vector length K: shape parameters for noise
#   thetas      – numeric vector length K: scale parameters for noise
#
# Returns:
#   list(
#     F = n×K matrix of latent factor scores,
#     Z = n×(sum G_k) matrix of simulated observations
#   )

simulate_multiblock_sgt <- function(n, G_vec, Sigma_f, batch_means, batch_sds, lambdas, thetas) {
  K <- length(G_vec)
  stopifnot(
    is.matrix(Sigma_f), ncol(Sigma_f) == K-1, nrow(Sigma_f) == K-1,
    length(batch_means) == K,
    length(batch_sds)   == K,
    length(lambdas)     == K,
    length(thetas)      == K
  )
  
  #  draw latent factors: n × K
  F_mat <- MASS::mvrnorm(n, mu = rep(0, K-1), Sigma = Sigma_f)
  F_mat <- cbind(F_mat,rep(0))
  
  # containers
  Z_blocks     <- vector("list", K)
  gamma_blocks <- vector("list", K)
  delta2_blocks<- vector("list", K)
  
  # simulate each block
  for (k in seq_len(K)) {
    out <- generate_block(
      n         = n,
      latent    = F_mat[, k],
      Gk        = G_vec[k],
      batch_m   = batch_means[k],
      batch_sd  = batch_sds[k],
      lambda    = lambdas[k],
      theta     = thetas[k]
    )
    Z_blocks[[k]]      <- out$data
    gamma_blocks[[k]]  <- out$gamma
    delta2_blocks[[k]] <- out$delta2
  }
  
  # combine
  Z          <- do.call(cbind, Z_blocks)
  gamma_true <- unlist(gamma_blocks)
  delta2_true<- unlist(delta2_blocks)
  
  list(
    Z           = Z,
    gamma_true  = gamma_true,
    delta2_true = delta2_true
  )
}


simulate_multiblock <- function(n, G_vec, Sigma_f, batch_means, batch_sds, lambdas, thetas) {
  K <- length(G_vec)
  stopifnot(
    is.matrix(Sigma_f), ncol(Sigma_f) == K, nrow(Sigma_f) == K,
    length(batch_means) == K,
    length(batch_sds)   == K,
    length(lambdas)     == K,
    length(thetas)      == K
  )
  
  #  draw latent factors: n × K
  F_mat <- MASS::mvrnorm(n, mu = rep(0, K), Sigma = Sigma_f)
  F_mat <- cbind(F_mat,rep(0))
  
  # containers
  Z_blocks     <- vector("list", K)
  gamma_blocks <- vector("list", K)
  delta2_blocks<- vector("list", K)
  
  # simulate each block
  for (k in seq_len(K)) {
    out <- generate_block(
      n         = n,
      latent    = F_mat[, k],
      Gk        = G_vec[k],
      batch_m   = batch_means[k],
      batch_sd  = batch_sds[k],
      lambda    = lambdas[k],
      theta     = thetas[k]
    )
    Z_blocks[[k]]      <- out$data
    gamma_blocks[[k]]  <- out$gamma
    delta2_blocks[[k]] <- out$delta2
  }
  
  # combine
  Z          <- do.call(cbind, Z_blocks)
  gamma_true <- unlist(gamma_blocks)
  delta2_true<- unlist(delta2_blocks)
  
  list(
    Z           = Z,
    gamma_true  = gamma_true,
    delta2_true = delta2_true
  )
}


