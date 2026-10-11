
# ==============================================================================
# 1. SETUP
# ==============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(tidyverse)
})


# ==============================================================================
# 2. SIMULATION PARAMETERS
# ==============================================================================

## ---- 2.1 Per-site covariance structures ------------------------------------
grid <- expand.grid(
  nSite = c(2, 4, 6, 8), nBlock = 2, p = 100,
  min.block = 0.8, max.block = 1.1,
  min.background = 0.2, max.background = 0.6
)
grid$seed <- runif(nrow(grid))
row <- grid[3, ]

J <- row$nSite
K <- row$nBlock
prop.vec <- rep(1 / K, K)
err.sigma.vec <- runif(J, min = 0, max = 0.1)

generate.Sigma.f <- function(K, seed) {
  set.seed(seed)
  diag  <- runif(K, min = row$min.block, max = row$max.block)
  rho_1 <- runif(1, min = row$min.background, max = row$max.background)
  rho_2 <- runif(1, min = row$min.background, max = row$max.background)
  rho_3 <- runif(1, min = row$min.background, max = row$max.background)
  p.mat <- matrix(c(diag[1], -rho_1, rho_2,
                    -rho_1, diag[2], rho_3,
                    rho_2, rho_3, diag[3]), nrow = 3)
  Spd.proj(SigA.hat = p.mat, eps = 0.001)$mat
}

JSite.pmat <- lapply(1:J, function(x) generate.Sigma.f(K = 3, seed = 16 + x))

## ---- 2.2 Batch (site) nuisance-parameter generators ------------------------
generate_param_set <- function(set_id, n, sig = TRUE) {
  set.seed(set_id)
  para <- list(
    batch_means = runif(n, min = 1, max = 2),
    batch_sds   = runif(n, min = 0.1, max = 0.4),
    lambdas     = runif(n, min = 5, max = 8),
    thetas      = runif(n, min = 0.2, max = 0.4)
  )
  if (sig) {
    para$batch_means <- c(para$batch_means, runif(1, min = 1, max = 2))
    para$batch_sds   <- c(para$batch_sds,   runif(1, min = 0.3, max = 0.4))
    para$lambdas     <- c(para$lambdas,     runif(1, min = 5, max = 8))
    para$thetas      <- c(para$thetas,      runif(1, min = 0.3, max = 0.7))
  }
  para
}

generate_param_set2 <- function(set_id, n, sig = TRUE) {
  set.seed(set_id)
  para <- list(
    batch_means = runif(n, min = 1, max = 2),
    batch_sds   = runif(n, min = 0.1, max = 0.4),
    lambdas     = runif(n, min = 5, max = 8),
    thetas      = runif(n, min = 0.2, max = 0.4)
  )
  if (sig) {
    para$batch_means <- c(para$batch_means, runif(1, min = 1, max = 2))
    para$batch_sds   <- c(para$batch_sds,   runif(1, min = 0.3, max = 0.4))
    para$lambdas     <- c(para$lambdas,     runif(1, min = 5, max = 8))
    para$thetas      <- c(para$thetas,      runif(1, min = 0.3, max = 0.7))
  }
  para
}

known_para <- list()
for (i in 1:2) known_para[[i]] <- generate_param_set(set_id = i + 12, n = 3, sig = TRUE)
for (i in 3:4) known_para[[i]] <- generate_param_set2(set_id = i + 123, n = 3, sig = TRUE)

## ---- 2.3 Covariate (AD / gender) effect structure --------------------------
m <- 200
overlap_frac_site <- 0

V_AD_full     <- 1:50
V_gender_full <- 51:100
V_shared      <- 150:200

set.seed(12)
AD_prob <- runif(4, 0.3, 0.7)

n_AD_shared     <- round(0.1 * length(V_shared))
AD_effect_idx   <- c(V_AD_full, sample(V_shared, n_AD_shared))
n_gender_shared <- round(0.1 * length(V_shared))
gender_effect_idx <- c(V_gender_full, sample(V_shared, n_gender_shared))

current_overlap <- intersect(AD_effect_idx, gender_effect_idx)
desired_overlap <- round(overlap_frac_site * length(V_shared))
if (length(current_overlap) < desired_overlap) {
  extra_overlap <- sample(setdiff(V_shared, current_overlap),
                          desired_overlap - length(current_overlap))
  AD_effect_idx     <- unique(c(AD_effect_idx, extra_overlap))
  gender_effect_idx <- unique(c(gender_effect_idx, extra_overlap))
}

ad_effect_sig_min_max <- c(0.1, 0.3)
AD_effect_base <- numeric(m)
AD_beta <- 0.8
delta   <- 0.01
AD_effect_base[intersect(AD_effect_idx, V_AD_full)] <-
  AD_beta + runif(length(intersect(AD_effect_idx, V_AD_full)), -delta, delta)

AD_shared_idx <- setdiff(AD_effect_idx, V_AD_full)
AD_effect_base[AD_shared_idx] <-
  runif(length(AD_shared_idx), ad_effect_sig_min_max[1], ad_effect_sig_min_max[2])


# ==============================================================================
# 3. DATA GENERATION
# ==============================================================================


p_site    <- c(300, 300, 300, 300)   # per-site N
G_vec     <- c(50, 50, 50, 50)
n_repeats <- 1                        

rep <- 1
set.seed(100 + rep)

df_list         <- vector("list", 4)
gamma_list      <- vector("list", 4)
delta_list      <- vector("list", 4)
true_AD_effects <- vector("list", 4)

for (i in 1:4) {
  set.seed(rep * 100 + i + 13)
  
  Sigma_f   <- JSite.pmat[[i]]
  para_list <- known_para[[i]]
  n <- p_site[i]
  AD_new <- rbinom(n, 1, AD_prob[i])
  true_AD_effects[[i]] <- AD_effect_base
  
  sim <- simulate_multiblock_sgt(
    n, G_vec, Sigma_f,
    para_list$batch_means, para_list$batch_sds,
    para_list$lambdas, para_list$thetas
  )
  df <- as.data.frame(sim$Z)
  
  if (length(AD_effect_idx) > 0) {
    subject_random_AD <- rnorm(n, 1, 0)
    df[, AD_effect_idx] <- df[, AD_effect_idx, drop = FALSE] +
      sweep(matrix(AD_new[1:n] * subject_random_AD, n, length(AD_effect_idx)),
            2, AD_effect_base[AD_effect_idx], `*`)
  }
  
  df$batch <- factor(i)
  df$AD    <- factor(AD_new)
  
  df_list[[i]]    <- df
  gamma_list[[i]] <- sim$gamma_true
  delta_list[[i]] <- sim$delta2_true
}

combined_df <- dplyr::bind_rows(df_list)

# ---- Package everything from this one-time data-generation run ------------
demo_data <- list(
  combined_df     = combined_df,      # the generated multi-site dataset
  gamma_list      = gamma_list,       # true gamma per batch
  delta_list      = delta_list,       # true delta^2 per batch
  true_AD_effects = true_AD_effects,  # true AD effect vector per batch
  AD_effect_idx   = AD_effect_idx,
  gender_effect_idx = gender_effect_idx,
  JSite.pmat      = JSite.pmat,
  known_para      = known_para,
  p_site          = p_site,
  G_vec           = G_vec
)

saveRDS(demo_data, file = "demo_data.rds")

# ==============================================================================
# 4. HARMONIZATION 
# ==============================================================================


mod   <- model.matrix(~ combined_df$AD)
dat   <- t(as.matrix(combined_df[, 1:200]))
batch <- as.factor(combined_df$batch)


site.outs <- lapply(unique(batch), function(b) {
  s <- batch == b
  fed_scombat_site_step1(dat[, s], batch[s], mod[s, ])
})
central <- fed_scombat_central_step1(site.outs, 0.5, 1)

site.outs2 <- lapply(unique(batch), function(b) {
  s <- batch == b
  fed_scombat_site_step1(dat[, s], batch[s], mod[s, ], central.out = central)
})
central2 <- fed_scombat_central_step1(site.outs2, 0.5, 1, assign_CID = 3)

site.outs3 <- Map(function(b, site_id) {
  s <- batch == b
  fed_scombat_site_step1(dat[, s], batch[s], mod[s, ],
                         central.out = central2, site_id = site_id)
}, unique(batch), seq_along(unique(batch)))

Y_new <- dplyr::bind_rows(lapply(site.outs3, function(x) as.data.frame(t(x$dat.combat))))



