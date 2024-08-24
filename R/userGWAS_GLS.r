.userGWAS_GLS = function (covstruc = covstruc, SNPs = SNPs,
                          estimator = "iterative",model = model,
                          diagGLS=diagGLS,fix_measurement = TRUE) {

 list.of.packages <- c("data.table", "GenomicSEM","dplyr","stringr","stringr","simsalapar","gdata","Matrix","lavaan","progress")
 new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
 if(length(new.packages)) stop("Missing package(s) ", paste0(new.packages, collapse=" and "))
 lapply(list.of.packages, library, character.only = TRUE)

  if (fix_measurement != TRUE & estimator == "analytic") {
    stop("Analytic estimator is only available when the argument fix_measurement = TRUE. Please set the fix_measurement argument to TRUE")
  }

  sumstatsGLS = SNPs
  nosnpmod = Model1   
  # Extract lambdas
  factors <- unique(nosnpmod$lhs[nosnpmod$op == "=~"])
  traits <- colnames(covstruc$S)
  num_traits <- ncol(covstruc$S)
  num_factors <- length(factors)
  # Generate all combinations of traits and factors
  combinations <- expand.grid(traits = traits, factors = factors)
  column_names <- paste0("lambda.", combinations$traits, "_", combinations$factors)
  # Store lambdas
  lambdas <- .extract_lambdas(nosnpmod)
  ### Create data structure for output ###
  GLS_mGWAS_results <- sumstatsGLS[, 1:6]
  # Extract necessary columns from sumstatsGLS
  betas <- sumstatsGLS %>% select(contains("beta."))
  SEs <- sumstatsGLS %>% select(contains("se."))

  ### Create data structure for output ###
  GLS_mGWAS_results <- sumstatsGLS[, 1:6]
  # Extract necessary columns from sumstatsGLS
  betas <- sumstatsGLS %>% select(contains("beta."))
  SEs <- sumstatsGLS %>% select(contains("se."))

  # Add factor-specific beta, SE, Z_beta, and p_val columns
  # Loop through each factor to add the necessary columns
  for (j in factors) {
    # Column names
    beta_col <- paste0("beta_", j)
    SE_col <- paste0("SE_", j)
    Z_beta_col <- paste0("Z_beta_", j)
    p_val_col <- paste0("p_val_", j)
    Q_col <- paste0("Q_", j)
    Q_df_col <- paste0("Q_df_", j)
    Qpval_col <- paste0("Qpval_", j)
    
   # Initialize columns with NA
    GLS_mGWAS_results[[beta_col]] <- NA
    GLS_mGWAS_results[[SE_col]] <- NA
    GLS_mGWAS_results[[Z_beta_col]] <- NA
    GLS_mGWAS_results[[p_val_col]] <- NA
    GLS_mGWAS_results[[Q_col]] <- NA
    GLS_mGWAS_results[[Q_df_col]] <- NA
    GLS_mGWAS_results[[Qpval_col]] <- NA
  } 

  # Add omnibus columns
  GLS_mGWAS_results <- GLS_mGWAS_results %>%
  mutate(Q_omnibus = NA,Q_omnibus_df = NA,Q_omnibus_pval = NA)

  # Initialize objects required for GLS calculation
  # Store lambdas
  lambdas_snp <- as.numeric(lambdas)
  # Create sampling correlation matrix
  R_SNP=covstruc$I 
  diag(R_SNP)[diag(R_SNP)<1]=1
  X = matrix(NA,num_traits,num_factors) 
  X[] = lambdas_snp
  colnames(X) = factors
  
  #Iterate across SNPs in sumstatsGLS input file
  pb <- progress_bar$new(
  format = paste("Running GLSmGWAS on", nrow(sumstatsGLS),"SNPs and",num_factors,"factors","[:bar]",":percent in :elapsed"),
  total = nrow(sumstatsGLS), clear = FALSE, width = 60)
  
  for (i in 1:nrow(sumstatsGLS)){
  SNP = sumstatsGLS[i,]
  beta_snp <- as.numeric(betas[i,])
  se_snp <- as.numeric(SEs[i,])  
  #Create V_SNP  
  V_SNP=cor2cov(R=as.matrix(R_SNP),sds=se_snp) #turn to sampling covariance matrix
  #GLS to conduct userGWAS and get betas and SEs for SNP effects on factors
    V_d = diag(diag(V_SNP))
    V = V_SNP
    y = beta_snp
  
  if(diagGLS==TRUE){
    Beta=solve(t(X) %*% solve(V_d) %*% X) %*% t(X) %*% solve(V_d) %*% y #beta (https://en.wikipedia.org/wiki/Generalized_least_squares and http://dx.doi.org/10.1080/10705511.2013.824793)
    bread=solve(t(X) %*% solve(V_d) %*% X)
    meat=t(X) %*% solve(V_d) %*% V %*% solve(V_d) %*% X
    sandwich=bread %*% meat %*% bread
    SE=sqrt(diag(sandwich))
  } else {
    Beta=solve(t(X) %*% solve(V) %*% X) %*% t(X) %*% solve(V) %*% y #beta (https://en.wikipedia.org/wiki/Generalized_least_squares and http://dx.doi.org/10.1080/10705511.2013.824793)
    rownames(Beta) = factors
    SE=sqrt(diag(solve(t(X) %*% solve(V) %*% X))) #SE of beta
  }
  names(SE) = factors

  # Model-implied SNP-phenotype betas and residuals of implied minus observed
  # Initialize beta_hats to zero
  beta_hats <- 0
  # Calculate beta_hats using all relevant lambda factors and beta values
  for (j in factors) {
    factor_col <- grep(paste0("lambda\\..*_", j, "$"), colnames(lambdas[1,]), value = TRUE)
    beta_hats <- as.numeric(beta_hats + lambdas[1, factor_col] * Beta[j,])
  }
  Resid=beta_snp-beta_hats

  ################# CALCULATE OMNIBUS AND FACTOR-SPECIFIC QSNP ##################### 
  # Omnibus Qsnp using Ronald's magic (see "Model fit statistics" section of Method in original Grotzinger et al. Genomic SEM paper)
  # this is a Qsnp for the entire factor model, not one factor at a time (i.e. compare a common pathways model with SNP effects only on the two factors to an independent pathways model for the SNP effects only on all of the indicators and not the factors)
  # Eigen decomposition of V_SNP
  E <- diag(eigen(V_SNP)$values)
  P1 <- eigen(V_SNP)$vectors
  # Initialize a list to store Q-statistics for each factor
  QsnpF_list <- list()
  Indicators_list <- list()
  # Loop through each factor to compute Q-statistics
  for (j in factors) {
    # Adjust Resid for the current factor
    lambda_cols <- grep(paste0("lambda\\..*_", j, "$"), colnames(lambdas), value = TRUE)
    factor_col <- grep(paste0("lambda\\..*_", j, "$"), colnames(lambdas[1,]), value = TRUE)
    Resid_adjusted <- Resid
    Resid_adjusted <- Resid_adjusted[lambdas[1, factor_col] != 0]
    Indicators_list[[j]]=length(Resid_adjusted)
    # Compute Q-statistic for the current factor
    E_factor=diag(eigen(V_SNP[lambdas[1, factor_col] != 0,lambdas[1, factor_col] != 0])$values)
    P1_factor=eigen(V_SNP[lambdas[1, factor_col] != 0,lambdas[1, factor_col] != 0])$vectors
    QsnpF <- Resid_adjusted %*% P1_factor %*% solve(E_factor) %*% t(P1_factor) %*% Resid_adjusted
    # Store the Q-statistic in the list
    QsnpF_list[[paste0("Qsnp_", j)]] <- QsnpF
  }

  # Calculate the overall Q-statistic
  Qsnp <- Resid %*% P1 %*% solve(E) %*% t(P1) %*% Resid 
  # Populate the GLS_mGWAS_results for factor-specific variables
    for (j in factors) {
      GLS_mGWAS_results[i, paste0("beta_", j)] <- Beta[j, 1]
      GLS_mGWAS_results[i, paste0("SE_", j)] <- SE[j]
      GLS_mGWAS_results[i, paste0("Q_", j)] <- QsnpF_list[[paste("Qsnp_",j,sep="")]][1]
    }
    # Add the Qmnibus column
  GLS_mGWAS_results[i, "Q_omnibus"] <- Qsnp
  pb$tick()  # Update the progress bar
  }

  # Calculate df and pval for factor-specfic Q statistics
  for (j in factors){
  GLS_mGWAS_results[, paste0("Z_beta_", j)] <- GLS_mGWAS_results[, paste0("beta_", j)]/GLS_mGWAS_results[, paste0("SE_", j)]
  GLS_mGWAS_results[, paste0("p_val_", j)] <- 2 * pnorm(-abs(GLS_mGWAS_results[, paste0("Z_beta_", j)]))
  GLS_mGWAS_results[, paste0("Q_df_", j)] <- Indicators_list[[j]] - 1
  GLS_mGWAS_results[, paste0("Qpval_", j)] <- pchisq(GLS_mGWAS_results[, paste0("Q_", j)],
                                                      GLS_mGWAS_results[, paste0("Q_df_", j)],lower.tail=FALSE)
  }
  # df and pval for Q omnibus
  GLS_mGWAS_results[, "Q_omnibus_df"] <- length(colnames(betas)) - length(Beta)
  GLS_mGWAS_results[, "Q_omnibus_pval"] <- pchisq(GLS_mGWAS_results[, "Q_omnibus"],
                                                   df=GLS_mGWAS_results[, "Q_omnibus_df"] , lower.tail = FALSE)
    
  #######################################################################################
  return(GLS_mGWAS_results) 
}
