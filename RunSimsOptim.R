## define variables ##

library(purrr)

nModels = 7
nReps = 3
nGen = 10
nVar = 10

## establish empty matrices to hold outputs for Selfing and Recombination Population ##

geneticvalues <- matrix(nrow=nGen, ncol=nReps)
correlations <- matrix(nrow=nModels, ncol=nReps)
variances <- matrix(nrow=nVar,ncol=nReps)
alleles <- vector("list", length = nReps)


## Run repeat loop to run reps ##

i = 1
repeat{
  source("rrblup_rd_snp_yield.R") ##Source the SCript for the SCenario you would like to run##
  
  geneticvalues[,i] <-gvMat
  
  correlations[,i] <- corMat
  
  variances[,i] <- varMat
  
  alleles[[i]] <- allelesMat
  
  i <- i + 1
  
  if (i > nReps){ ##break at number of desired reps##
    break
  }
  
  
  
  ##create data frames and label##
  geneticvalues <- as.data.frame(geneticvalues)
  colnames(geneticvalues) <- c(1:nReps)
  rownames(geneticvalues) <- c("PrevCycPYT","NewParents","F1","F2","F3","F4","F5","PYT","AYT","Variety")
  
  correlations <- as.data.frame(correlations)
  colnames(correlations) <- c(1:nReps)
  rownames(correlations) <- c("PrevCycPYT", "F2","F3","F4","F5","PYT","AYT")
  
  variances <- as.data.frame(variances)
  colnames(variances) <- c(1:nReps)
  rownames(variances) <- c("PrevCycPYT", "newParents","F1","F2", "F3","F4", "F5", "PYT","AYT",
                           "Variety")
  
  
  ##write files
  write.csv(geneticvalues, "rrblup_rd_gvs_snp_yield.csv")
  write.csv(correlations, "rrblup_rd_cors_snp_yield.csv")
  write.csv(variances,"rrblup_rd_vars_snp_yield.csv")
  saveRDS(alleles, file="rrblup_rd_alleles_snp_yield.rds")
  
}
