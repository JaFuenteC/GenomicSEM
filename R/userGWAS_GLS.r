# Define the userGWAS2 function with GLS integration
userGWAS2 = function (covstruc = NULL, SNPs = NULL, estimation = "DWLS",
  estimator = "iterative",model = "", usermod = NULL,diagGLS=TRUE,
  printwarn = TRUE, sub = FALSE, cores = NULL, 
  toler = FALSE, SNPSE = FALSE, parallel = TRUE, GC = "standard", 
  MPI = FALSE, smooth_check = FALSE, TWAS = FALSE, std.lv = FALSE, 
  fix_measurement = TRUE, Q_SNP = FALSE) {
 list.of.packages <- c("data.table", "GenomicSEM","dplyr","stringr","stringr","simsalapar","gdata","Matrix","lavaan","progress")
 new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
 if(length(new.packages)) stop("Missing package(s) ", paste0(new.packages, collapse=" and "))
 lapply(list.of.packages, library, character.only = TRUE)

  if (!toler) 
  toler <- .Machine$double.eps
  GenomicSEM:::.check_one_of(estimation, c("DWLS", "ML"))
  GenomicSEM:::.check_one_of(estimator, c("iterative", "analytic"))
  GenomicSEM:::.check_boolean(printwarn)
  if (!is.null(cores)) 
  GenomicSEM:::.check_range(cores, min = 0, max = Inf)
  GenomicSEM:::.check_range(toler, min = 0, max = Inf)
  GenomicSEM:::.check_boolean(parallel)
  GenomicSEM:::.check_one_of(GC, c("standard", "conserv", "none"))
  GenomicSEM:::.check_boolean(MPI)
  GenomicSEM:::.check_boolean(smooth_check)
  GenomicSEM:::.check_boolean(TWAS)
  GenomicSEM:::.check_boolean(std.lv)
  time <- proc.time()
  Operating <- Sys.info()[["sysname"]]
  if (MPI == TRUE & Operating == "Windows") {
    stop("MPI is not currently available for Windows operating systems. Please set the MPI argument to FALSE, or switch to a Linux or Mac operating system.")
  }
  if (fix_measurement != TRUE & estimator == "analytic") {
    stop("Analytic estimator is only available when the argument fix_measurement = TRUE. Please set the fix_measurement argument to TRUE")
  }
  if (exists("Output")) {
    stop("Please note that an update was made to commonfactorGWAS on 4/1/21 so that addSNPs output CANNOT be fed directly to the function. It now expects the\n            output from ldsc (using covstruc = ...)  followed by the output from sumstats (using SNPs = ... ) as the first two arguments.")
  }
  test <- c(str_detect(model, "~"), str_detect(model, "="), 
    str_detect(model, "\\+"))
  if (!all(test)) {
    warning("Your model name may be listed in quotes; please remove the quotes and try re-running if the function has returned stopped running after returning an error.")
  }
  print("Please note that an update was made to userGWAS on Sept 1 2023  so that the default behavior is to fix the measurement model using the fix_measurement argument.")
  if (class(SNPs)[1] == "character") {
    print("You are likely listing arguments in the order of a previous version of userGWAS, if you have yur results stored after running addSNPs you can still explicitly call Output = ... to provide them to userGWAS. The current version of the function is faster and saves memory. It expects the\n          output from ldsc followed by the output from sumstats (using SNPs = ... ) as the first two arguments. See ?userGWAS for help on propper usag")
    warning("You are likely listing arguments (e.g. Output = ...) in the order of a previous version of userGWAS, if you have yur results stored after running addSNPs you can still explicitly call Output = ... to provide them to userGWAS. The current version of the function is faster and saves memory. It expects the\n            output from ldsc (using covstruc = ...)  followed by the output from sumstats (using SNPs = ... ) as the first two arguments. See ?userGWAS for help on propper usage")
  }  else {
    if (is.null(SNPs) | is.null(covstruc)) {
      print("You may be listing arguments in the order of a previous version of userGWAS, if you have yur results stored after running addSNPs you can still explicitly call Output = ... to provide them to userGWAS;if you already did this and the function ran then you can disregard this warning. The current version of the function is faster and saves memory. It expects the\n            output from ldsc followed by the output from sumstats (using SNPs = ... ) as the first two arguments. See ?userGWAS for help on propper usag")
      warning("You may be listing arguments (e.g. Output = ...) in the order of a previous version of userGWAS, if you have yur results stored after running addSNPs you can still explicitly call Output = ... to provide them to userGWAS; ; if you already did this and the function ran then you can disregard this warning. The current version of the function is faster and saves memory. It expects the\n              output from ldsc (using covstruc = ...)  followed by the output from sumstats (using SNPs = ... ) as the first two arguments. See ?userGWAS for help on propper usage")
    }
  }
  if (sub[[1]] != FALSE) {
    sub <- str_replace_all(sub, fixed(" "), "")
  }
  SNPs <- data.frame(SNPs)
  if (TWAS) {
    SNPs$Gene <- as.character(SNPs$Gene)
    SNPs$Panel <- as.character(SNPs$Panel)
    varSNP <- SNPs$HSQ
  }  else {
    SNPs$A1 <- as.character(SNPs$A1)
    SNPs$A2 <- as.character(SNPs$A2)
    SNPs$SNP <- as.character(SNPs$SNP)
    varSNP <- 2 * SNPs$MAF * (1 - SNPs$MAF)
  }
  if (SNPSE == FALSE) {
    varSNPSE2 <- (5e-04)^2
  }  else {
    varSNPSE2 <- SNPSE^2
  }
  V_LD <- as.matrix(covstruc[[1]])
  S_LD <- as.matrix(covstruc[[2]])
  I_LD <- as.matrix(covstruc[[3]])
  Model1 <- model
  

  ############################## GLS ANALYTIC ESTIMATOR ############################
  if(estimator=="analytic" & fix_measurement == TRUE) {
  # Reorder sumstats columns
  sumstats_colnames <- c("SNP", "CHR", "BP", "MAF", "A1", "A2",
                        colnames(SNPs)[grep(colnames(SNPs),patter="beta.*")],
                          colnames(SNPs)[grep(colnames(SNPs),patter="se.*")])
  
  # Select the desired columns from the dataframe
  sumstatsGLS = SNPs[, sumstats_colnames]
  
  # If usermod is NULL
  if(is.null(usermod)){
  nosnp_model <- paste(grep("~\\s*SNP$", strsplit(model, "\n")[[1]], value = TRUE, invert = TRUE), collapse = "\n")
  # Extract the factor names from these lines
  snp_lines <- grep("~\\s*SNP$", strsplit(model, "\n")[[1]], value = TRUE)
  factors <- unique(gsub("~.*", "", snp_lines))
  # Fit no-SNP model
  nosnpmod <- usermodel(covstruc,estimation = "DWLS", 
                        model = nosnp_model, 
                        CFIcalc = FALSE, std.lv = FALSE, imp_cov = FALSE)  
  usermod = nosnpmod
  }   
  
  nosnpmod = usermod   
  # Extract lambdas from usermodel output
  factors <- unique(nosnpmod$results$lhs[nosnpmod$results$op == "=~"])
  traits <- colnames(covstruc$S)
  num_traits <- ncol(covstruc$S)
  num_factors <- length(factors)
  # Generate all combinations of traits and factors
  combinations <- expand.grid(traits = traits, factors = factors)
  column_names <- paste0("lambda.", combinations$traits, "_", combinations$factors)
  extract_lambdas <- function(df) {
   # Define the traits and factors
   traits <- colnames(covstruc$S)
   # Initialize lambdas vector dynamically
   num_traits <- length(traits)
   num_factors <- num_factors
   lambdas <- rep(0, num_traits * num_factors)
   # Fill the lambdas vector
   for (factor_idx in seq_along(factors)) {
     factor <- factors[factor_idx]
     for (trait_idx in seq_along(traits)) {
       trait <- traits[trait_idx]
       row <- df[df$lhs == factor & df$rhs == trait & df$op == "=~", ]
       if (nrow(row) > 0) {
         lambda_value <- row$Unstand_Est
       } else {
         lambda_value <- 0
       }
       lambdas[(factor_idx - 1) * num_traits + trait_idx] <- lambda_value
     }
   }
  # Convert lambdas to a data frame and assign column names
  lambdas_df <- data.frame(matrix(lambdas, nrow=1, byrow=TRUE))
  colnames(lambdas_df) <- column_names
    return(lambdas_df)
  }
  # Store lambdas
  lambdas <- extract_lambdas(nosnpmod$results)
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

  } else {

  if (fix_measurement) {
    rownames(S_LD) <- colnames(S_LD)
    lines <- strsplit(model, "\n")[[1]]
    if (TWAS) {
      filtered_lines <- lines[!grepl(c("Gene"), lines)]
    } else {
      filtered_lines <- lines[!grepl("SNP", lines)]
    }
    filtered_lines <- filtered_lines[!grepl(c(":="), filtered_lines)]
    noSNPmodel <- paste(filtered_lines, collapse = "\n")
    smoothS <- ifelse(eigen(S_LD)$values[nrow(S_LD)] <= 
      0, S_LD <- as.matrix((nearPD(S_LD, corr = FALSE))$mat), 
      S_LD <- S_LD)
    smoothV <- ifelse(eigen(V_LD)$values[nrow(V_LD)] <= 
      0, V_LD <- as.matrix((nearPD(V_LD, corr = FALSE))$mat), 
      V_LD <- V_LD)
    W <- solve(V_LD, tol = toler)
    testnoSNP <- GenomicSEM:::.tryCatch.W.E(ReorderModelnoSNP <- sem(noSNPmodel, 
      sample.cov = S_LD, estimator = "DWLS", WLS.V = W, 
      sample.nobs = 2, optim.dx.tol = +Inf, optim.force.converged = TRUE, 
      control = list(iter.max = 1), std.lv = std.lv))
    order <- GenomicSEM:::.rearrange(k = ncol(S_LD), fit = ReorderModelnoSNP, 
      names = colnames(S_LD))
    V_Reorder <- V_LD[order, order]
    u <- nrow(V_Reorder)
    W_Reorder <- diag(u)
    diag(W_Reorder) <- diag(V_Reorder)
    W_Reorder <- solve(W_Reorder, tol = toler)
    if (estimation == "DWLS") {
      emptynoSNP <- GenomicSEM:::.tryCatch.W.E(Model1_Results <- sem(noSNPmodel, 
        sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, 
        sample.nobs = 2, optim.dx.tol = +Inf, std.lv = std.lv))
    }
    if (estimation == "ML") {
      emptynoSNP <- GenomicSEM:::.tryCatch.W.E(Model1_Results <- sem(noSNPmodel, 
        sample.cov = S_LD, estimator = "ML", sample.nobs = 200, 
        optim.dx.tol = +Inf, sample.cov.rescale = FALSE, 
        std.lv = std.lv))
    }
    Model1 <- parTable(Model1_Results)
    for (p in 1:nrow(Model1)) {
      Model1$free[p] <- ifelse(Model1$lhs[p] != Model1$rhs[p], 
        0, Model1$free[p])
    }
  }
  beta_SNP <- SNPs[, grep("beta.", fixed = TRUE, colnames(SNPs))]
  SE_SNP <- SNPs[, grep("se.", fixed = TRUE, colnames(SNPs))]
  n_phenotypes <- ncol(beta_SNP)
  diag(I_LD) <- ifelse(diag(I_LD) <= 1, 1, diag(I_LD))
  coords <- which(I_LD != "NA", arr.ind = T)
  V_SNP <- GenomicSEM:::.get_V_SNP(SE_SNP, I_LD, varSNP, GC, coords, n_phenotypes, 
    1)
  V_full <- GenomicSEM:::.get_V_full(n_phenotypes, V_LD, varSNPSE2, V_SNP)
  kv <- nrow(V_full)
  smooth2 <- ifelse(eigen(V_full)$values[kv] <= 0, V_full <- as.matrix((nearPD(V_full, 
    corr = FALSE))$mat), V_full <- V_full)
  S_Full <- GenomicSEM:::.get_S_Full(n_phenotypes, S_LD, varSNP, beta_SNP, 
    TWAS, 1)
  ks <- nrow(S_Full)
  smooth1 <- ifelse(eigen(S_Full)$values[ks] <= 0, S_Full <- as.matrix((nearPD(S_Full, 
    corr = FALSE))$mat), S_Full <- S_Full)
  k2 <- ncol(S_Full)
  for (i in 1) {
    W <- solve(V_full, tol = toler)
    test2 <- GenomicSEM:::.tryCatch.W.E(ReorderModel <- sem(model, sample.cov = S_Full, 
      estimator = "DWLS", WLS.V = W, sample.nobs = 2, 
      optim.dx.tol = +Inf, optim.force.converged = TRUE, 
      control = list(iter.max = 1), std.lv = std.lv))
    if (fix_measurement) {
      withSNP <- parTable(ReorderModel)
      if (TWAS) {
        for (p in 1:nrow(withSNP)) {
          if (withSNP$rhs[p] == "Gene" | withSNP$lhs[p] == 
            "Gene") {
            Model1 <- rbind(Model1, withSNP[p, ])
          }
        }
      } else {
        for (p in 1:nrow(withSNP)) {
          if (withSNP$rhs[p] == "SNP" | withSNP$lhs[p] == 
            "SNP") {
            Model1 <- rbind(Model1, withSNP[p, ])
          }
        }
      }
      for (p in 1:nrow(withSNP)) {
        if (withSNP$op[p] == ":=") {
          Model1 <- rbind(Model1, withSNP[p, ])
        }
      }
      test3 <- GenomicSEM:::.tryCatch.W.E(ReorderModel <- sem(Model1, 
        sample.cov = S_Full, estimator = "DWLS", WLS.V = W, 
        sample.nobs = 2, optim.dx.tol = +Inf, optim.force.converged = TRUE, 
        control = list(iter.max = 1), std.lv = std.lv))
    }
    order <- GenomicSEM:::.rearrange(k = k2, fit = ReorderModel, names = rownames(S_Full))
    suppressWarnings(df <- lavInspect(ReorderModel, "fit")["df"])
    suppressWarnings(npar <- lavInspect(ReorderModel, "fit")["npar"])
  }
  if (TWAS) {
    SNPs <- SNPs[, 1:3]
  }  else {
    SNPs <- SNPs[, 1:6]
  }
  f <- nrow(beta_SNP)
  LavModel1 <- GenomicSEM:::.userGWAS_main(i = 1, cores = 1, n_phenotypes, 
    n = 1, I_LD, V_LD, S_LD, std.lv, varSNPSE2, order, SNPs, 
    beta_SNP, SE_SNP, varSNP, GC, coords, smooth_check, 
    TWAS, printwarn, toler, estimation, sub, Model1, df, 
    npar, returnlavmodel = TRUE, Q_SNP = Q_SNP, model = model)
  if (!parallel) {
    if (sub[[1]] == FALSE) {
      Results_List <- vector(mode = "list", length = nrow(beta_SNP))
    }
    if (TWAS) {
      print("Starting TWAS Estimation")
    } else {
      print("Starting GWAS Estimation")
    }
    for (i in 1:nrow(beta_SNP)) {
      if (i == 1) {
        cat(paste0("Running Model: ", i, "\n"))
      }    else {
        if (i%%1000 == 0) {
          cat(paste0("Running Model: ", i, "\n"))
        }
      }
      final2 <- GenomicSEM:::.userGWAS_main(i, cores = 1, n_phenotypes, 
        1, I_LD, V_LD, S_LD, std.lv, varSNPSE2, order, 
        SNPs, beta_SNP, SE_SNP, varSNP, GC, coords, 
        smooth_check, TWAS, printwarn, toler, estimation, 
        sub, Model1, df, npar, basemodel = LavModel1, 
        Q_SNP = Q_SNP, model = model)
      final2$i <- NULL
      if (sub[[1]] != FALSE) {
        final3 <- as.data.frame(matrix(NA, ncol = ncol(final2), 
          nrow = length(sub)))
        final3[1:length(sub), ] <- final2[1:length(sub), 
          ]
        if (i == 1) {
          Results_List <- vector(mode = "list", length = length(sub))
          for (y in 1:length(sub)) {
            Results_List[[y]] <- as.data.frame(matrix(NA, 
              ncol = ncol(final3), nrow = f))
            colnames(Results_List[[y]]) <- colnames(final2)
            Results_List[[y]][1, ] <- final3[y, ]
          }
        } else {
          for (y in 1:nrow(final3)) {
            Results_List[[y]][i, ] <- final3[y, ]
          }
        }
      }  else {
        Results_List[[i]] <- final2
      }
    }
    time_all <- proc.time() - time
    print(time_all[3])
    return(Results_List)
  }  else {
    if (is.null(cores)) {
      int <- min(c(nrow(SNPs2), detectCores() - 1))
    } else {
      if (cores > nrow(SNPs)) 
        warning(paste0("Provided number of cores was greater than number of SNPs, reverting to cores=", 
          nrow(SNPs2)))
      int <- min(c(cores, nrow(SNPs)))
    }
    if (MPI) {
      cl <- getMPIcluster()
      registerDoParallel(cl)
    } else {
      if (Operating != "Windows") {
        cl <- makeCluster(int, type = "FORK")
      } else {
        cl <- makeCluster(int, type = "PSOCK")
      }
      registerDoParallel(cl)
      on.exit(stopCluster(cl))
    }
    SNPs <- suppressWarnings(split(SNPs, 1:int))
    beta_SNP <- suppressWarnings(split(beta_SNP, 1:int))
    SE_SNP <- suppressWarnings(split(SE_SNP, 1:int))
    varSNP <- suppressWarnings(split(varSNP, 1:int))
    if (TWAS) {
      print("Starting TWAS Estimation")
    }  else {
      print("Starting GWAS Estimation")
    }
    if (Operating != "Windows") {
    results <- foreach(n = icount(int), .combine = 'rbind', .packages = 'lavaan') %dopar% {
    # Process each n
    temp_results <- foreach(i = 1:nrow(beta_SNP[[n]]), .combine = 'rbind', .packages = 'lavaan') %do% {
    GenomicSEM:::.userGWAS_main(
      i, int, n_phenotypes, n, I_LD, V_LD, S_LD, std.lv, 
      varSNPSE2, order, SNPs[[n]], beta_SNP[[n]], 
      SE_SNP[[n]], varSNP[[n]], GC, coords, smooth_check, 
      TWAS, printwarn, toler, estimation, sub, Model1, 
      df, npar, basemodel = LavModel1, Q_SNP = Q_SNP, 
      model = model)
            }
      temp_results
        }
    }  else {
      utilfuncs <- list()
    utilfuncs[[".tryCatch.W.E"]] <- GenomicSEM:::.tryCatch.W.E
    utilfuncs[[".get_V_SNP"]] <- GenomicSEM:::.get_V_SNP
    utilfuncs[[".get_Z_pre"]] <- GenomicSEM:::.get_Z_pre
    utilfuncs[[".get_V_full"]] <- GenomicSEM:::.get_V_full

    results <- foreach(n = icount(int), .combine = "rbind") %:% 
        foreach(i = 1:nrow(beta_SNP[[n]]), .combine = "rbind", 
            .packages = c("lavaan", "gdata"), .export = c(".userGWAS_main")) %dopar% 
    {
        GenomicSEM:::.userGWAS_main(
            i, int, n_phenotypes, n, I_LD, V_LD, S_LD, std.lv, 
            varSNPSE2, order, SNPs[[n]], beta_SNP[[n]], 
            SE_SNP[[n]], varSNP[[n]], GC, coords, smooth_check, 
            TWAS, printwarn, toler, estimation, sub, Model1, 
            df, npar, utilfuncs, basemodel = LavModel1, Q_SNP = Q_SNP, 
            model = model)
        }
    }
    results <- results[order(results$i), ]
    results$i <- NULL
    if (sub[[1]] != FALSE) {
      Results_List <- vector(mode = "list", length = length(sub))
      for (y in 1:length(sub)) {
        Results_List[[y]] <- as.data.frame(matrix(NA, 
          ncol = ncol(results), nrow = nrow(results)/length(sub)))
        colnames(Results_List[[y]]) <- colnames(results)
        Results_List[[y]] <- subset(results, paste0(results$lhs, 
          results$op, results$rhs, sep = "") %in% sub[[y]] | 
          is.na(results$lhs))
      }
      rm(results)
    }
    if (sub[[1]] == FALSE) {
      if (TWAS) {
        names <- unique(results$Panel)
        Results_List <- vector(mode = "list", length = length(names))
        for (y in 1:length(names)) {
          Results_List[[y]] <- subset(results, results$Panel == 
            names[[y]])
          Results_List[[y]]$Model_Number <- NULL
        }
      } else {
        names <- unique(results$SNP)
        Results_List <- vector(mode = "list", length = length(names))
        for (y in 1:length(names)) {
          Results_List[[y]] <- subset(results, results$SNP == 
            names[[y]])
          Results_List[[y]]$Model_Number <- NULL
        }
      }
      rm(results)
      rm(names)
    }
    time_all <- proc.time() - time
    print(time_all[3])
    return(Results_List)

    }
  }
}

# RUN EXAMPLE WITH PSYCH TRAITS
# First load required packages, ldsc() and sumstats() output
library(data.table)
library(GenomicSEM)
load("LDSC_PSYCH.RData")
sumstats = fread("Psych_sumstats_4GLS.txt",data.table=FALSE)
sumstats = sumstats[1:20,] # first 20 SNPs for this example

# USING CONVENTIONAL ITERATIVE APPROACH
# Here we store the userGWAS model syntax (model with SNP effects)
model = "Comp=~OCD+TS+AN
Psych=~a*SCZ+a*BIP
Neuro=~ADHD+MDD+TS+ASD+PTSD
Int=~MDD+ANX+PTSD
SUD=~ALCH+CUD+OUD
Comp~~Psych+Neuro+Int+SUD
Psych~~Neuro+Int+SUD
Neuro~~Int+SUD
Int~~SUD
Comp~SNP
Psych~SNP
Neuro~SNP
Int~SNP
SUD~SNP"

# Define parameters of interest (SNP - factor associations)
sub = c("Comp~SNP","Psych~SNP","Neuro~SNP","Int~SNP","SUD~SNP")

# Run the function using conventional iterative approach
GWASoutput_iterative = userGWAS2(covstruc = LDSC_P, SNPs = sumstats, estimation = "DWLS",
                    estimator = "iterative",model = model,diagGLS=NULL, 
                    usermod = NULL, printwarn = TRUE, sub = sub, cores = NULL, 
                    toler = FALSE, SNPSE = FALSE, parallel = FALSE, GC = "standard", 
                    MPI = FALSE, smooth_check = FALSE, TWAS = FALSE, std.lv = FALSE, 
                    fix_measurement = TRUE, Q_SNP = TRUE) 
# Check output
GWASoutput_iterative

# USING ANALYTIC GLS-BASED APPROACH WITHOUT PROVIDING THE NO-SNP MODEL OUTPUT
# In this case we use the GLS approach WITHOUT providing the no-SNP model output (the argument usermod = NULL).
# Since the no-SNP model output is not provided, the function automatically extracts the factor names from the 
# userGWAS model syntax (as specified in the model argument) and estimates the no-SNP model to extract the (unstandardized) lambda
# parameters, and fix the measurement model across SNPs

# userGWAS model syntax (model with SNP effects)
model = "Comp=~OCD+TS+AN
Psych=~a*SCZ+a*BIP
Neuro=~ADHD+MDD+TS+ASD+PTSD
Int=~MDD+ANX+PTSD
SUD=~ALCH+CUD+OUD
Comp~~Psych+Neuro+Int+SUD
Psych~~Neuro+Int+SUD
Neuro~~Int+SUD
Int~~SUD
Comp~SNP
Psych~SNP
Neuro~SNP
Int~SNP
SUD~SNP"

# Define parameters of interest (SNP - factor associations)
sub = c("Comp~SNP","Psych~SNP","Neuro~SNP","Int~SNP","SUD~SNP")

GWASoutput_analytic = userGWAS2(covstruc = LDSC_P, SNPs = sumstats, estimation = "DWLS",
                    estimator = "analytic",model = model,diagGLS=TRUE, 
                    usermod = NULL, printwarn = TRUE, sub = sub, cores = NULL, 
                    toler = FALSE, SNPSE = FALSE, parallel = FALSE, GC = "standard", 
                    MPI = FALSE, smooth_check = FALSE, TWAS = FALSE, std.lv = FALSE, 
                    fix_measurement = TRUE, Q_SNP = TRUE) 
# Check output
GWASoutput_analytic

# USING ANALYTIC GLS-BASED APPROACH WITHOUT PROVIDING THE NO-SNP MODEL OUTPUT
# In this case we use the GLS approach PROVIDING the no-SNP model output (the argument usermod = CFAmod).
# The function automatically extracts the (unstandardized) lambda parameters from the output of the CFA model,
# fix the measurement model across SNPs

# Define the no-SNP model
noSNPmodel <- "Comp=~OCD+TS+AN
Psych=~a*SCZ+a*BIP
Neuro=~ADHD+MDD+TS+ASD+PTSD
Int=~MDD+ANX+PTSD
SUD=~ALCH+CUD+OUD
Comp~~Psych+Neuro+Int+SUD
Psych~~Neuro+Int+SUD
Neuro~~Int+SUD
Int~~SUD
"

# Estimate the no-SNP model using the usermodel() function
CFAmodel <- usermodel(LDSC_P,estimation = "DWLS", 
                        model = noSNPmodel, 
                        CFIcalc = FALSE, std.lv = FALSE, imp_cov = FALSE)

# Run the userGWAS function using diagonalized GLS and providing the no-SNP model output
GWASoutput_analytic_2 = userGWAS2(covstruc = LDSC_P, SNPs = sumstats, estimation = "DWLS",
                    estimator = "analytic",model = model,diagGLS=TRUE, 
                    usermod = CFAmodel, printwarn = TRUE, sub = sub, cores = NULL, 
                    toler = FALSE, SNPSE = FALSE, parallel = FALSE, GC = "standard", 
                    MPI = FALSE, smooth_check = FALSE, TWAS = FALSE, std.lv = FALSE, 
                    fix_measurement = TRUE, Q_SNP = TRUE) 
# Check output
GWASoutput_analytic_2


# CHECK CORRESPONDENCE BETWEEN ITERATIVE AND ANALYTIC APPROACHES
# COMP factor
cor(GWASoutput_iterative[[1]]$est,GWASoutput_analytic$beta_Comp)#SNP betas
cor(GWASoutput_iterative[[1]]$Z_Estimate,GWASoutput_analytic$Z_beta_Comp)#SNP betas
cor(GWASoutput_iterative[[1]]$ Q_SNP,GWASoutput_analytic$Q_Comp)#Q stat
# Psych factor
cor(GWASoutput_iterative[[2]]$est,GWASoutput_analytic$beta_Psych)#SNP betas
cor(GWASoutput_iterative[[2]]$Z_Estimate,GWASoutput_analytic$Z_beta_Psych)#SNP betas
cor(GWASoutput_iterative[[2]]$ Q_SNP,GWASoutput_analytic$Q_Psych)#Q stat
# Neuro factor
cor(GWASoutput_iterative[[3]]$est,GWASoutput_analytic$beta_Neuro)#SNP betas
cor(GWASoutput_iterative[[3]]$Z_Estimate,GWASoutput_analytic$Z_beta_Neuro)#SNP betas
cor(GWASoutput_iterative[[3]]$ Q_SNP,GWASoutput_analytic$Q_Neuro)#Q stat
# Int factor
cor(GWASoutput_iterative[[4]]$est,GWASoutput_analytic$beta_Int)#SNP betas
cor(GWASoutput_iterative[[4]]$Z_Estimate,GWASoutput_analytic$Z_beta_Int)#SNP betas
cor(GWASoutput_iterative[[4]]$ Q_SNP,GWASoutput_analytic$Q_Int)#Q stat
# SUD factor
cor(GWASoutput_iterative[[5]]$est,GWASoutput_analytic$beta_SUD)#SNP betas
cor(GWASoutput_iterative[[5]]$Z_Estimate,GWASoutput_analytic$Z_beta_SUD)#SNP betas
cor(GWASoutput_iterative[[5]]$ Q_SNP,GWASoutput_analytic$Q_SUD)#Q stat


