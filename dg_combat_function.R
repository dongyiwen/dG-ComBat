
# ==============================================================================
# 1. MAIN ENTRY POINTS —  (site <-> central)
# ==============================================================================

## ------------------------------------------------------------------------
## 1.1 fed_scombat_site_step1()
##     Site-side step of the function. Behavior depends on what
##     `central.out` currently holds:
##       - central.out = NULL              -> Round 1: return local LS pieces
##       - central.out$var.pooled = NULL   -> Round 2: return local covariance
##                                            summary (cov.data, s.data, ...)
##       - central.out fully populated     -> Round 3: return corrected/
##                                            harmonized data (dat.combat)
##
## Args:
##   dat          p x n matrix (features x subjects) for this site
##   batch        factor of batch/site labels for this site's subjects
##   mod          model matrix of covariates (e.g. diagnosis), or NULL
##   central.out  output of fed_scombat_central_step1() from the previous
##                round, or NULL on the very first round
##   eb           unused flag retained from original signature
##   parametric   unused flag retained from original signature
##   mean.only    if TRUE, allows single-subject batches / mean-only ComBat
##   verbose      print progress messages
##   site_id      integer index of this site into central.out$hyper_all /
##                central.out$eb.estimate_all (required once central.out is
##                fully populated, i.e. Round 3)
##
## Returns:
##   Round 1: list(ls.site, dataDict, sigma.site = NULL)
##   Round 2: list(ls.site, dataDict, cov.data, s.data, mu.hat.unstand,
##                 sigma.site)
##   Round 3: list(dat.combat, estimates)
## ------------------------------------------------------------------------
fed_scombat_site_step1 <- function(dat, batch, 
                                   mod=NULL,
                                   central.out=NULL,
                                   eb=TRUE, 
                                   parametric=TRUE, mean.only = FALSE,verbose=TRUE,site_id = NULL){
  dat <- as.matrix(dat)
  #  .checkConstantRows(dat)
  # .checkNARows(dat)
  ## Check for missing values
  hasNAs <- any(is.na(dat))
  if (hasNAs & verbose){
    cat(paste0("Found ", sum(is.na(dat)), " missing data values. \n"))
  }
  if(mean.only){
    if (verbose) cat("Performing ComBat with mean only\n")
  }
  
  ##################### Getting design ############################
  dataDict <- getDataDictDC(batch, mod, verbose=verbose, mean.only=mean.only)
  design <- dataDict[["design"]]
  n.array<-dataDict[["n.array"]]
  ####################################################################
  
  
  ############### Site matrices for standardization #################
  # W^T W used in LS estimation
  ls_site <- NULL
  ls_site[[1]] <- crossprod(design)
  ls_site[[2]] <- tcrossprod(t(design), dat)
  
  
  # new dataDict with batches within current site
  inclBat <- dataDict$n.batches > 0
  dataDictSite <- dataDict
  dataDictSite$batch <- droplevels(dataDict$batch)
  dataDictSite$batches <- dataDict$batches[inclBat]
  dataDictSite$n.batch <- sum(inclBat)
  dataDictSite$n.batches <- dataDict$n.batches[inclBat]
  dataDictSite$batch.design <- as.matrix(dataDict$batch.design[,inclBat])
  
  if (is.null(central.out)) {
    site_out <- list(
      ls.site = ls_site,
      dataDict = dataDict,
      sigma.site = NULL
    )
    return(site_out)
  }
  # If beta.estimates given, get summary statistics for cov matrix
  
  stdObjects <- getStandardizedDataDC(dat=dat, 
                                      dataDict=dataDict,
                                      design=design,
                                      hasNAs=hasNAs,
                                      central.out=central.out
  )
  s.data <- stdObjects[["s.data"]]
  
  if (is.null(central.out$var.pooled)) {
    sigma_site<- getSigmaSummary(dat, dataDict, design, hasNAs, central.out)
    site_out <- list(
      ls.site = ls_site,
      dataDict = dataDict,
      cov.data = stdObjects[["cov.data"]],
      s.data = stdObjects[["s.data"]],
      mu.hat.unstand = stdObjects[["mu.hat.unstand"]],
      sigma.site = sigma_site
    )
    return(site_out)
  }
  
  var.pooled = unname(central.out$var.pooled)
  update.s.data<- s.data/(tcrossprod(var.pooled, rep(1,n.array)))
  
  p=nrow(update.s.data)

  if (is.null(site_id)) stop("site_id must be provided when central.out is not NULL.")
  hyper <- central.out$hyper_all[[site_id]]
  eb.estimate <- central.out$eb.estimate_all[[site_id]]
  correct.data<-getCorrectedData(dat,
                                 update.s.data, 
                                 dataDictSite, 
                                 eb.estimate, 
                                 central.out,
                                 stdObjects,
                                 hyper
  )
  
  estimates <- list(gamma.star=eb.estimate[["gamma"]],
                    delta.star=eb.estimate[["delta2"]], 
                    stand.mean=stdObjects[["stand.mean"]], 
                    mod.mean=stdObjects[["mod.mean"]], 
                    mod=mod,
                    mu = hyper[["mu_bar"]],
                    tau2 = hyper[["tau2_bar"]],
                    lambda = hyper[["lambda_bar"]],
                    theta = hyper[["theta_bar"]],
                    s.data = update.s.data
  )
  
  site_out <- list(dat.combat=correct.data, estimates=estimates)
  return(site_out )
  
}


## ------------------------------------------------------------------------
## 1.2 fed_scombat_central_step1()
##     Central-side aggregation step. 
##
## Args:
##   site.outs     list of per-site outputs from fed_scombat_site_step1(),
##                 all from the same round
##   penalty_para  penalty parameter passed to est.multi() for structure
##                 estimation
##   reference_site index of the site whose covariance is used as the
##                 alignment reference
##   assign_CID    FALSE to use est.multi()'s own CID, or an integer to
##                 keep only the first N components of it
##   known_CID/known_L  optionally supply a precomputed CID/L instead of
##                 estimating them
##   verbose       unused flag retained from original signature
##
## Returns (depends on round, same as before):
##   Round 1 (step1 == TRUE):  list(B.hat, stand.mean, var.pooled = NULL)
##   Round 2 (step1 == FALSE): list(B.hat, stand.mean, L_hat, var.pooled,
##                                 CID, reference, hyper_all, eb.estimate_all,
##                                 cov_full_all, new_cov.data, [Sigma.weight,
##                                 W.dense])
## ------------------------------------------------------------------------
fed_scombat_central_step1  <- function(site.outs,penalty_para,reference_site,
                                       assign_CID = FALSE,known_CID = FALSE,known_L = FALSE, verbose = FALSE) {
  J <- length(site.outs) # number of sites
  
  # get n.batches and n.array from sites
  batch_levels <- levels(site.outs[[1]]$dataDict$batch)
  n.batches <- Reduce("+", lapply(site.outs, function(x) x$dataDict$n.batches))
  n.batch <- length(n.batches)
  n.array <- sum(n.batches)
  n.arrays <- lapply(site.outs, function(x) x$dataDict$n.array)
  
  step1s <- sapply(site.outs, function(x) is.null(x$cov.data))
  if (length(unique(step1s)) > 1) {
    stop("Not all sites are at the same step, please confirm with each site.")
  }
  step1 <- all(step1s)
  
  #### Step 1: Get LS estimate across sites ####
  ls1 <- Reduce("+", lapply(site.outs, function(x) x$ls.site[[1]]))
  ls2 <- Reduce("+", lapply(site.outs, function(x) x$ls.site[[2]]))
  #lambda <- 1e-7
  #ls1_reg <- ls1 + diag(lambda, ncol(ls1))
  #B.hat <-crossprod(solve(ls1_reg), ls2)
  B.hat <-crossprod(ginv(ls1), ls2)
  #B.hat <- crossprod(solve(ls1), ls2)
  grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
  stand.mean <- crossprod(grand.mean, t(rep(1,n.array)))
  
  if (step1) {
    central_out <- list(
      B.hat = B.hat,
      stand.mean = stand.mean,
      var.pooled = NULL)
    return(central_out)
  }
  
  # --- Step 2: Standardization Parameters from Covariance Matrices ---
  vars <- lapply(site.outs, function(x) diag(x$cov.data))
  
  var.pooled = rep(0, length(vars[[1]]))
  for (i in 1:J) {
    var.pooled = var.pooled + n.arrays[[i]]*vars[[i]]
  }
  var.pooled = var.pooled/n.array
  
  stand.var.pooled<- sqrt(var.pooled)
  
  var.matrix <- outer(stand.var.pooled, stand.var.pooled, "*")
  
  # --- Standardize cov.data and mu.hat.unstand ---
  stand_cov.data <- lapply(site.outs, function(mat) {
    p<- mat$dataDict$n.array
    mat$cov.data*p/(p-1) / var.matrix})
  
  
  new_mu.hat.unstand <- lapply(site.outs, function(mat) {
    p <- ncol(mat$cov.data)
    matrix(mat$mu.hat.unstand / stand.var.pooled, nrow = p, ncol = 1)
  })
  
  # --- Step 3: Estimate Low-rank Structure (L) ---
  if (!isFALSE(known_CID) && !isFALSE(known_L)) {
    # if user provided both
    CID <- known_CID
    L_hat <- known_L
  }else{
    structure.res <- est.multi(n.batches, n.array, stand_cov.data, penalty_para)
    L_hat <- structure.res$L
    if (isFALSE(assign_CID)) {
      CID <- structure.res$CID
    } else if (is.numeric(assign_CID) && length(assign_CID) == 1) {
      CID <- structure.res$CID[seq_len(assign_CID)]
    } else {
      stop("assign_CID must be FALSE or a single numeric value.")
    }
  }
  
  
  
  # --- Step 4: Empirical Bayes Estimation ---
  hyper_all <- list()
  eb.estimate_all <- list()
  cov_full_all <- list()
  
  for (i in 1:J) {
    hyper_all[[i]] <- estimate_hyperparams(n.arrays[[i]], new_mu.hat.unstand[[i]], L_hat, CID, stand_cov.data[[i]])
    eb.estimate_all[[i]] <- run_combat_iterative(n.arrays[[i]], new_mu.hat.unstand[[i]], L_hat, stand_cov.data[[i]], CID)
    
    Sigma_fhat <- hyper_all[[i]]$Sigma_fhat
    delta2_diag <- diag(eb.estimate_all[[i]]$delta2)
    cov_full_all[[i]] <- L_hat %*% Sigma_fhat %*% t(L_hat) + delta2_diag
  }
  
  reference_matrix <- cov_full_all[[reference_site]]
  
  # --- Final Output ---
  # Base list of common elements
  central_out <- list(
    B.hat = B.hat,
    stand.mean = stand.mean,
    L_hat = L_hat,
    var.pooled = stand.var.pooled,
    CID = CID,
    reference = reference_matrix,
    hyper_all = hyper_all,
    eb.estimate_all = eb.estimate_all,
    cov_full_all = cov_full_all,
    new_cov.data = stand_cov.data
  )
  
  # Add Sigma.weight and W.dense if known_CID or known_L are not provided
  if (isFALSE(known_CID) || isFALSE(known_L)) {
    central_out$Sigma.weight <- structure.res$updated.Sigma.comb
    central_out$W.dense <- structure.res$W_dense
  }
  
  
  return(central_out)
}


# ==============================================================================
# 2. DATA DICTIONARY / DESIGN HELPERS
# ==============================================================================

## ------------------------------------------------------------------------
## 2.1 getDataDictDC()
##     Builds the batch/design bookkeeping used throughout the protocol:
##     which subject indices belong to which batch, the combined
##     batch + covariate design matrix, and basic counts/checks (e.g.
##     refusing single-subject batches unless mean.only = TRUE).
##
## Args:
##   batch      factor of batch/site labels
##   mod        covariate model matrix, or NULL
##   verbose    print progress messages
##   mean.only  allow single-subject batches (mean-only ComBat)
##
## Returns: list with batch, batches (indices per batch), n.batch,
##   n.batches, n.array, n.covariates, design, batch.design
## ------------------------------------------------------------------------
getDataDictDC <- function(batch, mod, verbose, mean.only){
  batch <- as.factor(batch)
  n.batch <- nlevels(batch)
  batches <- lapply(levels(batch), function(x)which(batch==x))
  n.batches <- sapply(batches, length)
  n.array  <- sum(n.batches)
  batchmod <- model.matrix(~-1+batch)  
  if (verbose) cat("[combat] Found",nlevels(batch),'batches\n')
  if(any(n.batches==1) & mean.only==FALSE){
    stop("Found one site with only one sample. Consider using the mean.only=TRUE option")
  }
  
  #combine batch variable and covariates
  design <- cbind(batchmod,mod)
  n.covariates <- ncol(design)-ncol(batchmod)
  if (verbose) cat("[combat] Adjusting for ",n.covariates,' covariate(s) or covariate level(s)\n')
  out <- list()
  #Making sure to keep track of names:
  names(batches)   <- names(n.batches) <- levels(batch)
  colnames(design) <- gsub("batch", "", colnames(design))
  out[["batch"]] <- batch
  out[["batches"]] <- batches
  out[["n.batch"]] <- n.batch
  out[["n.batches"]] <- n.batches
  out[["n.array"]] <- n.array
  out[["n.covariates"]] <- n.covariates
  out[["design"]] <- design
  out[["batch.design"]] <- design[,1:n.batch]
  return(out)
}


# ==============================================================================
# 3. STANDARDIZATION HELPERS
# ==============================================================================

## ------------------------------------------------------------------------
## 3.1 getStandardizedDataDC()
##     Centers the raw data using the pooled B.hat/stand.mean from
##     fed_scombat_central_step1(), removing the covariate ("mod") mean, and
##     computes the resulting per-site mean/covariance summary
##     (mu.hat.unstand, cov.data) that gets sent back to the central step.
##
## Args:
##   dat          p x n matrix for this site
##   dataDict     output of getDataDictDC()
##   design       design matrix (batch + covariates)
##   hasNAs       whether dat contains missing values (unused directly here,
##                kept for interface consistency)
##   central.out  central aggregation output containing B.hat, stand.mean
##
## Returns: list(s.data, stand.mean, mod.mean, beta.hat, cov.data,
##   mu.hat.unstand)
## ------------------------------------------------------------------------
getStandardizedDataDC <- function(dat, dataDict, design, hasNAs, central.out){
  batches=dataDict$batches
  n.batches=dataDict$n.batches
  n.batch=dataDict$n.batch
  n.array=dataDict$n.array
  
  
  B.hat <- central.out$B.hat
  stand.mean <- central.out$stand.mean[,1:n.array]
  
  ##################### standardization #######################
  
  if(!is.null(design)){
    tmp <- design
    tmp[,c(1:n.batch)] <- 0
    mod.mean <- t(tmp%*%B.hat)
    #stand.mean <- stand.mean+t(tmp%*%B.hat)
  } else {
    mod.mean <- 0
  }
  s.data <- (dat-stand.mean-mod.mean)
  
  rownames(stand.mean) <- rownames(mod.mean) <- rownames(dat)
  colnames(stand.mean) <- colnames(mod.mean) <- colnames(dat)
  
  ##################### mu hat #######################
  
  mu_hat_unstand<-colMeans(t(s.data))
  
  ##################### covariance structure #######################
  
  
  #cov.data<-cov(t(s.data))* (nrow(t(s.data))-1)/(nrow(t(s.data)))## or
  cov.data<-cov(t(s.data))* (nrow(t(s.data))-1)/(nrow(t(s.data)))
  
  
  return(list(s.data=s.data, 
              stand.mean=stand.mean,
              mod.mean=mod.mean, 
              beta.hat=B.hat,
              cov.data=cov.data,
              mu.hat.unstand = mu_hat_unstand
              
  )
  )
}


## ------------------------------------------------------------------------
## 3.2 getSigmaSummary()
##     Alternative, simpler per-site pooled-variance summary computed
##     directly from residuals of the design fit (dat - design %*% B.hat),
##     used as `sigma.site` in fed_scombat_site_step1()'s Round 2 return value.
##
## Args:
##   dat          p x n matrix for this site
##   dataDict     output of getDataDictDC()
##   design       design matrix (batch + covariates)
##   hasNAs       whether dat contains missing values
##   central.out  central aggregation output containing B.hat, stand.mean
##
## Returns: numeric vector of per-feature pooled variance estimates
## ------------------------------------------------------------------------
getSigmaSummary <- function(dat, dataDict, design, hasNAs, central.out){
  
  batches=dataDict$batches
  n.batches=dataDict$n.batches
  n.array=dataDict$n.array
  n.batch=dataDict$n.batch
  B.hat <- central.out$B.hat
  stand.mean <- central.out$stand.mean[,1:n.array]
  factors <- (n.array/(n.array-1))
  var.pooled <- rowVars(dat-t(design %*% B.hat), na.rm=TRUE)/factors
  
  return(var.pooled)
}


# ==============================================================================
# 4. EMPIRICAL BAYES CORRECTION (final harmonization step)
# ==============================================================================

## ------------------------------------------------------------------------
## 4.1 getCorrectedData()
##
## Args:
##   dat          original p x n matrix for this site (used only for
##                rownames on the output)
##   s.data       standardized data (post var.pooled scaling)
##   dataDict     site-restricted data dictionary (dataDictSite)
##   eb.estimate  this site's empirical Bayes estimate (gamma, delta2)
##   central.out  full central aggregation output (L_hat, reference, CID, ...)
##   stdObjects   output of getStandardizedDataDC() (stand.mean, mod.mean)
##   hyper        this site's hyperparameter estimate (contains Sigma_fhat)
##
## Returns: p x n matrix of harmonized ("ComBat-corrected") data, with
##   rownames restored from `dat`
## ------------------------------------------------------------------------
getCorrectedData <- function(dat, 
                             s.data, 
                             dataDict, 
                             eb.estimate, 
                             central.out, 
                             stdObjects, 
                             hyper) {
  # Extract needed components
  stand.mean <- stdObjects$stand.mean
  mod.mean <- stdObjects$mod.mean
  batches <- dataDict$batches
  batch.design <- dataDict$batch.design
  n.batches <- dataDict$n.batches
  n.array <- dataDict$n.array
  
  var.pooled <- central.out$var.pooled
  gamma.star <- eb.estimate[["gamma"]]
  delta.star <- eb.estimate[["delta2"]]
  L_hat <- central.out$L_hat
  Sigma_fhat <- hyper$Sigma_fhat
  G=nrow(s.data)
  
  
  Sigma.comb <- central.out$reference
  ##
  
  
  # --- Helper for standard group correction ---
  correct_standard <- function(rows) {
    print("have block structures")
    sm <- stand.mean[rows, ]
    mm <- mod.mean[rows, ]
    sd <- s.data[rows, ]
    vp <- var.pooled[rows]
    gs <- gamma.star[rows]
    ds <- delta.star[rows]
    
    
    new_sigma <- L_hat %*% Sigma_fhat %*% t(L_hat) + diag(delta.star)
    new_sigma_sub <- new_sigma[rows, rows]
    Sigmaiv <- Spd.proj(SigA.hat = new_sigma_sub)$mat
    U_list <- chol(Sigmaiv)
    U_inv<- backsolve(U_list, diag(ncol(U_list)))
    
    
    Sigma_ref <- Spd.proj(SigA.hat = Sigma.comb[rows, rows])$mat
    U_ref <- chol(Sigma_ref)
    
    U_inv<-as.matrix(U_inv)
    U_ref<-as.matrix(U_ref)
    
    
    for (i in batches) {
      sd[, i] <- sd[, i] - t(batch.design[i, ] %*% t(gs))
    }
    
    ##based on correlation align
    
    
    stand_sd <- t(sd) %*% U_inv %*% U_ref
    final <- (t(stand_sd) * tcrossprod(vp, rep(1, n.array))) + sm + mm
    return(final)
  }
  
  # --- Helper for signal group correction ---
  correct_signal <- function(rows) {
    print("have singleton variables")
    sm <- stand.mean[rows, ]
    mm <- mod.mean[rows, ]
    sd <- s.data[rows, ]
    vp <- var.pooled[rows]
    gs <- gamma.star[rows ]
    ds <- delta.star[rows]
    
    sc <- diag(Sigma.comb[rows, rows])
    #sc<- rep(1,length(diag_sc))
    
    j <- 1
    for (i in batches) {
      top <- sd[, i] - t(batch.design[i, ] %*% t(gs))
      bottom <- tcrossprod(sqrt(ds), rep(1, n.batches[j]))
      sd[, i] <- top / bottom
      j <- j + 1
    }
    #inv_d <- 1 / sc
    #stand_sd <- sd * tcrossprod(sc, rep(1, n.array))
    stand_sd<-sd* tcrossprod(sqrt(sc), rep(1, n.array))
    final <- (stand_sd * tcrossprod(vp, rep(1, n.array))) + sm + mm
    return(final)
  }
  
  if(sum(central.out$CID) != G){
    rows_to_keep <- L_hat[, ncol(L_hat)] != 1
    rows_to_remove <- !rows_to_keep
    
    corrected_str <- correct_standard(rows_to_keep)
    corrected_sig <- correct_signal(rows_to_remove)
    
    # Recombine and restore order
    combined <- rbind(corrected_str, corrected_sig)
    final_data <- combined[order(c(which(rows_to_keep), which(rows_to_remove))), ]
    
  }else{
    rows_to_keep<- which(rowSums(L_hat) == 1)
    final_data <- correct_standard(rows_to_keep)
  }
  
  rownames(final_data) <- rownames(dat)
  return(final_data)
}


# ==============================================================================
# 5. SPD PROJECTION / ADMM UTILITIES
# ==============================================================================
## ------------------------------------------------------------------------
## 5.1 Spd.proj()
##     Projects a symmetric matrix onto the (numerically) positive-definite
##     cone: if its smallest eigenvalue is already above `eps`, returns it
##     unchanged (after symmetrizing); otherwise projects via ADMM_proj().
##
## Args:
##   SigA.hat  symmetric matrix to project
##   eps       minimum acceptable eigenvalue; defaults to 5/p if NULL
##
## Returns: list(mat = projected matrix, conv = 1 if no projection was
##   needed, 2 if ADMM_proj() was used)
## ------------------------------------------------------------------------
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

## ------------------------------------------------------------------------
## 5.2 l1proj()
##     Efficient Euclidean projection of a vector onto an L1 ball of radius
##     `b`. Used inside ADMM_proj()'s S-update step.
##     Ref: Duchi et al. (2008), "Efficient Projections onto the L1-Ball
##     for Learning in High Dimensions", ICML.
##
## Args:
##   v  numeric vector to project
##   b  L1 ball radius (must be > 0)
##
## Returns: numeric vector, the projection of v onto {w : ||w||_1 <= b}
## ------------------------------------------------------------------------
l1proj<-function(v, b){
  
  # Efficient projection onto L1 ball of specified radius (i.e. b), used by the admm algo
  # Ref. Duchi et al. (2008). Efficient Projections onto the L1-Ball for Learning in High Dimensions, ICML
  
  stopifnot(b>0)
  
  u <- sort(abs(v),decreasing=TRUE)
  sv <- cumsum(u)
  rho <- max(which(u>(sv-b)/1:length(u)))
  theta <- max(0, (sv[rho]-b)/rho)
  w <-sign(v) * pmax(abs(v)-theta,0)
  
  return(w)
}

## ------------------------------------------------------------------------
## 5.3 ADMM_proj()
##     ADMM algorithm that projects a symmetric matrix onto the PSD cone
##     (with eigenvalues >= epsilon) while sparsely regularizing the
##     off-diagonal deviation (S) between the projected matrix (R) and the
##     input (mat), via l1proj(). Called by Spd.proj() when a matrix's
##     smallest eigenvalue is too small to use directly.
##
## Args:
##   mat            symmetric matrix to project
##   epsilon        minimum eigenvalue floor for the projected matrix
##   mu             ADMM penalty parameter (halved every 20 iterations)
##   it.max         maximum number of ADMM iterations
##   etol           convergence tolerance on R/S/primal-residual changes
##   etol_distance  convergence tolerance on ||R - mat|| stabilizing
##
## Returns: list(mat = projected matrix R, df_ADMM = per-iteration
##   diagnostics data frame: eps_R, eps_S, eps_primal, time, distance)
## ------------------------------------------------------------------------
ADMM_proj<-function(mat,
                    epsilon=1e-4,
                    mu=10,
                    it.max=1e3,
                    etol=1e-4,
                    etol_distance = 1e-4){
  
  
  
  p<-nrow(mat)
  
  # Initialization
  R<-diag(mat)
  S<-matrix(0,p,p)
  L<-matrix(0,p,p)
  
  itr<-0
  iteration <- eps_R <- eps_S <- eps_primal <- time <- distance <- NULL
  while (itr<it.max) {
    #print(itr)
    Rp<-R
    Sp<-S
    start <- Sys.time()
    # Subproblem I: R step
    W<-mat+S+mu*L
    W.eigdec<-eigen(W, symmetric=TRUE) 
    W.V<-W.eigdec$vectors
    W.D<-W.eigdec$values
    R<-W.V%*%diag(pmax(W.D,epsilon))%*%t(W.V)
    
    # Subproblem II: S step
    M<-R-mat-mu*L     
    S[lower.tri(S, diag = TRUE)]<-M[lower.tri(M, diag = TRUE)]-l1proj(v=M[lower.tri(M, diag = TRUE)],b=mu/2)    
    for (i in 2:p){
      for (j in 1:(i-1)){
        S[j,i]<-S[i,j]
      }
    }
    
    # L step: update the Lagrange parameter
    L<-L-(R-S-mat)/mu
    end <- Sys.time()
    #Stocking the values of different parameters with the number of iterations
    iteration <- c(iteration, itr)
    eps_R <- c(eps_R,max(abs(R-Rp)))
    eps_S <- c(eps_S,max(abs(S-Sp)))
    eps_primal <- c(eps_primal, max(abs(R-S-mat)))
    time <- c(time, end - start)
    distance <- c(distance,max(abs(R-mat)))
    
    # Stopping Rule                        
    #cat("check the stopping criterion:",max(abs(R-S-mat)),"\n")
    if (((max(abs(R-Rp))<etol) && (max(abs(S-Sp))<etol) && (max(abs(R-S-mat))<etol)) || (abs(max(abs(Rp-mat)) - max(abs(R-mat)))<etol_distance)){
      itr<-it.max
    } else {
      itr<-itr+1
    }
    
    if (itr%%20==0) {
      mu<-mu/2
    }
  }
  df_ADMM <- data.frame(iteration = iteration, eps_R = eps_R, eps_S=eps_S, eps_primal=eps_primal, time=time, distance=distance)
  return(list(mat=R,df_ADMM=df_ADMM))
  
}