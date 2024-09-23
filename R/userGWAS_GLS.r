.userGWAS_GLS = function (covstruc = covstruc, SNPs = SNPs,
                          nosnpmod = Model1) {

 list.of.packages <- c("data.table", "GenomicSEM","dplyr","stringr","stringr","simsalapar","gdata","Matrix","lavaan","progress")
 new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
 if(length(new.packages)) stop("Missing package(s) ", paste0(new.packages, collapse=" and "))
 lapply(list.of.packages, library, character.only = TRUE)
 VMV = function (V1,M,V2){V1 %*% M %*% V2}
  sumstatsGLS = SNPs 
  # Extract lambdas from usermodel output
  factors <- unique(nosnpmod$lhs[nosnpmod$op == "=~"])
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
         lambda_value <- row$est
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
  lambdas <- extract_lambdas(nosnpmod)
  
  ### Create data structure for output ###
  GLS_mGWAS_results <- sumstatsGLS[, 1:6]
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
  
  # Extract necessary columns from sumstatsGLS
  betas <- sumstatsGLS %>% select(contains("beta."))
  SEs <- sumstatsGLS %>% select(contains("se."))

  # Initialize objects required for GLS calculation
  # Store lambdas
  lambdas_snp <- as.numeric(lambdas)
  # Create sampling correlation matrix
  R_SNP=covstruc$I 
  diag(R_SNP)[diag(R_SNP)<1]=1
  X = matrix(NA,num_traits,num_factors) 
  X[] = lambdas_snp
  colnames(X) = factors
  
  # Create lists of V matrices and SNP-trait betas vectors
  # List with vectors of SNP-trait betas
  beta_l <- lapply(split(betas, seq(nrow(betas))), function(x) as.numeric(unlist(x)))
  # List with vectors of standard errors of SNP-trait betas
  se_snp <- lapply(split(SEs, seq(nrow(SEs))), function(x) as.numeric(unlist(x)))
  #Create V_SNP list  
  V_SNP_list <- lapply(se_snp, function(se) {
                  cor2cov(R = as.matrix(R_SNP), sds = as.numeric(se))
  })
  #Create diagonalized V_SNP list
  V_d_list <- lapply(V_SNP_list, function(V_SNP) {
                diag(diag(V_SNP))  # Extract the diagonal elements and create a diagonal matrix
  })
    
  # Compute SNP-factor betas
  Beta_list <- unname(Map(function(V_d, beta) {
      # Compute factor Beta for each SNP using diagonalized SNP-specific V_SNP matrices
      solve(t(X) %*% solve(V_d) %*% X) %*% t(X) %*% solve(V_d) %*% beta
  }, V_d_list, beta_l))
  # SNP-factor betas SEs with sandwich correction
  # Initialize an empty list to store SE_parallel for each SNP
  SE_parallel_list <- list()
  # Map across the SNP-specific matrices
    SE_parallel_list <- Map(function(V_d, V_SNP) {
      bread <- solve(t(X) %*% solve(V_d) %*% X)     
      meat <- t(X) %*% solve(V_d) %*% V_SNP %*% solve(V_d) %*% X    
      sandwich <- bread %*% meat %*% bread  
      sqrt(diag(sandwich))    
    }, V_d_list, V_SNP_list)
  # Unname the list
  SE_parallel_list <- unname(SE_parallel_list)

  # Convert SE_parallel_list to a data frame
  SE_parallel_df <- do.call(rbind, SE_parallel_list)
  SE_parallel_df <- as.data.frame(SE_parallel_df)
  colnames(SE_parallel_df) <- paste("SE_",factors,sep="")

  # Convert Beta_list to a data frame
  Beta_parallel_df <- do.call(rbind, lapply(Beta_list, as.numeric))
  Beta_parallel_df <- as.data.frame(Beta_parallel_df)
  colnames(Beta_parallel_df) <- paste("Beta_",factors,sep="")
    
  # Z_Statistics
  Z_df_parallel = Beta_parallel_df/SE_parallel_df
  colnames(Z_df_parallel) <- paste("Z_",factors,sep="")
  
  ################################### QSNP CALCULATION #####################################################
  solveI=solve(covstruc$I) #calculate this just once (constant acrossSNPs)  
  ########################## 1. QSNP OMNIBUS ##########################
    # Initialize an empty list to store the beta_hats for each SNP
    # First, ensure each element in Beta_list is a numeric vector
    beta_hats_list <- lapply(Beta_list, function(beta) {
      # Convert the named column vector to a numeric matrix if necessary
      beta_matrix <- as.matrix(beta) 
      # Perform matrix multiplication 
      beta_hat <- beta_matrix[,1] %*% t(X) # Extract column vector and multiply by t(X)
      return(beta_hat) 
    })

    # Calculate residuals for each SNP and store them as numeric vectors without names
    Resid_parallel_list <- mapply(function(beta, beta_hat) {
      # Ensure that beta is converted to a numeric vector
      beta_vector <- as.vector(beta) 
      # Compute residuals
      residuals <- as.numeric(beta_vector - beta_hat)
      # Remove any names from the resulting vector
      unname(residuals)
    }, beta_l, beta_hats_list, SIMPLIFY = FALSE)

    compute_inside_list_omnibus <- function(solveI, se_snp_list) {
        # Loop through each SNP in the se_snp_list
        inside_list <- lapply(se_snp_list, function(SEs_snp) {
          # Ensure SEs_snp is a numeric vector
          SEs_snp <- as.numeric(SEs_snp)
          # Calculate the inside matrix for this SNP
          solveI / (SEs_snp %*% t(SEs_snp))
        })
        
        return(inside_list)  # Return a list of inside matrices (one per SNP)
      }

    inside_list = compute_inside_list_omnibus(solveI, se_snp)

    #OMNIBUS Q#
    Q_Omnibus_parallel = mapply(Resid_parallel_list,inside_list,Resid_parallel_list,FUN=VMV) 
  #################################### 2. FACTOR-SPECIFIC QSNP ########################################## 
  # Precompute lambdas_factor for each factor and store in a list
    lambdas_factor_list <- lapply(factors, function(factor_name) {
      factor_col <- grep(paste0("lambda\\..*_", factor_name, "$"), colnames(lambdas), value = TRUE)
      lambdas_factor <- lambdas[1, factor_col]
      lambdas_factor != 0  # Return TRUE where lambdas_factor is not equal to 0
    })
    names(lambdas_factor_list) <- factors

  # Create list with solveI subset for each factor  
    # Subset SolveI based on lambdas_factor_list
    solveI_list <- lapply(lambdas_factor_list, function(lambdas_factor) {
      # Subset SolveI by selecting rows and columns where lambdas_factor is TRUE
      solveI[lambdas_factor, lambdas_factor, drop = FALSE]
    })
    # Name the list elements based on the factor names
    names(solveI_list) <- factors
       
    # Subset SEs based on lambdas_factor_list
    SEs_factor_list <- lapply(seq_along(factors), function(i) {
        # Get the logical vector for the current factor
        lambdas_factor <- lambdas_factor_list[[i]]
        # Subset SEs by selecting the traits (columns) with non-zero loadings
        SEs[, lambdas_factor, drop = FALSE]
      })
    names(SEs_factor_list) <- factors
  
  # Convert the list of residuals into a data frame
  Resid_parallel_df <- do.call(rbind, Resid_parallel_list)
  # Subset residuals based on lambdas_factor_list
  Resid_factor_list <- lapply(seq_along(factors), function(i) {
    # Get the logical vector for the current factor
    lambdas_factor <- lambdas_factor_list[[i]]
    # Subset Resid_parallel_df by selecting the traits (columns) with non-zero loadings
    Resid_parallel_df[, lambdas_factor, drop = FALSE]
  })
  # Name the list elements based on the factors
  names(Resid_factor_list) <- factors
  # Store number of indicators per factor
  Indicators_list = lapply(Resid_factor_list, ncol)
   
  # Inside Matrix Calculation (using precomputed SolveI_list and already subset SEs_factor_list)
  compute_inside_list <- function(solveI_list, SEs_factor_list) {
    # Loop through factors (top-level list)
    inside_list_factor <- lapply(seq_along(solveI_list), function(i) {
      # Get the corresponding SolveI for this factor
      SolveI <- solveI_list[[i]]
      # Get the corresponding SEs_factor_list for this factor
      SEs_factor <- SEs_factor_list[[i]]
      # Now, for each SNP in this factor, compute the inside matrix
      lapply(seq_len(nrow(SEs_factor)), function(snp_index) {
        # Extract SEs for the current SNP as a numeric vector (row-wise)
        SEs_snp <- as.numeric(SEs_factor[snp_index, ])
        # Calculate the inside matrix for this SNP
        SolveI / (SEs_snp %*% t(SEs_snp))
            })  # End SNP loop
      })  # End factor loop
    return(inside_list_factor)
    }
  inside_list_factor <- compute_inside_list(solveI_list, SEs_factor_list)    
  
  # Calculate factor-specific QSNP statistics
  QsnpF_list <- lapply(seq_along(factors), function(i) {
                  factor <- factors[i]
                  # Retrieve the inside_list_factor for the current factor (i-th element)
                  inside_matrix_list <- inside_list_factor[[i]]            
                  # Extract residuals for the current factor
                  residuals_matrix <- Resid_factor_list[[i]]  # Matrix of residuals for this factor
                  # Compute QsnpF for the current factor, iterating over each SNP
                  QsnpF <- mapply(function(snp_idx, inside_matrix) {
                    # Extract the residuals for the current SNP
                    residuals <- residuals_matrix[snp_idx, , drop = FALSE]  # Residuals for a single SNP
                    # Apply the VMV function for this SNP
                    VMV(V1 = as.numeric(residuals), M = inside_matrix, V2 = as.numeric(residuals))
                    }, snp_idx = 1:nrow(residuals_matrix), inside_matrix = inside_matrix_list, SIMPLIFY = TRUE)
                  return(QsnpF)
    })
    # Combine results into a data frame
    QsnpF_df <- do.call(cbind, QsnpF_list)
    colnames(QsnpF_df) <- paste0("Q_", factors)
    
    # Populate the GLS_mGWAS_results for factor-specific variables
      for (j in factors) {
        GLS_mGWAS_results[, paste0("beta_", j)] <- Beta_parallel_df[, paste("Beta_",j,sep="")]
        GLS_mGWAS_results[, paste0("SE_", j)] <- SE_parallel_df[, paste("SE_",j,sep="")]
        GLS_mGWAS_results[, paste0("Z_beta_", j)] <- Z_df_parallel[, paste("Z_",j,sep="")]
        GLS_mGWAS_results[, paste0("Q_", j)] <- QsnpF_df[, paste("Q_",j,sep="")]
      }
      # Add the Qmnibus column
      GLS_mGWAS_results[, "Q_omnibus"] <- Q_Omnibus_parallel   

    # Calculate df and pval for factor-specfic Q statistics
        for (j in factors){
        GLS_mGWAS_results[, paste0("p_val_", j)] <- 2 * pnorm(-abs(GLS_mGWAS_results[, paste0("Z_beta_", j)]))
        GLS_mGWAS_results[, paste0("Q_df_", j)] <- Indicators_list[[j]] - 1
        GLS_mGWAS_results[, paste0("Qpval_", j)] <- pchisq(GLS_mGWAS_results[, paste0("Q_", j)],
                                                            GLS_mGWAS_results[, paste0("Q_df_", j)],lower.tail=FALSE)
        }
    # df and pval for Q omnibus
        GLS_mGWAS_results[, "Q_omnibus_df"] <- length(colnames(betas)) - ncol(Beta_parallel_df)
        GLS_mGWAS_results[, "Q_omnibus_pval"] <- pchisq(GLS_mGWAS_results[, "Q_omnibus"],
                                                        df=GLS_mGWAS_results[, "Q_omnibus_df"] , lower.tail = FALSE)  
  
    #pb$finish()   # Update the progress bar
    return(GLS_mGWAS_results)

}
