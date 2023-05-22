library(purrr)

nModels = 7
nReps = 2
nGen = 10
nVar = 9

## establish empty matrices to hold outputs for Selfing and Recombination Population ##

geneticvalues <- matrix(nrow=nGen, ncol=nReps)
correlations <- matrix(nrow=nModels, ncol=nReps)
variances <- matrix(nrow=nVar,ncol=nReps)
alleles <- vector("list", length = nReps)
bv_ebv <- vector("list", length = nReps)



## Run repeat loop to run reps ##

i = 1
repeat{
  source("rrblup_rd_snp_yield.R") ##Source the SCript for the SCenario you would like to run##
  
  geneticvalues[,i] <-gvMat
  
  correlations[,i] <- corMat
  
  variances[,i] <- varMat
  
  alleles[[i]] <- allelesMat
  
  bv_ebvC1[[i]] <- bv_ebv

  
  i <- i + 1
  
  if (i > nReps){ ##break at number of desired reps##
    break
  }
  
  
  
  ##create data frames and label##
  geneticvalues <- as.data.frame(geneticvalues)
  gain <- as.data.frame(geneticvaluesC1[10,] - geneticvaluesC1[2,])
  Allgeneticvalues <- as.data.frame(rbind(geneticvaluesC1, C1gain))
  colnames(geneticvalues) <- c(1:nReps)
  rownames(geneticvalues) <- c("PrevCycPYT","NewParents","F1","F2","F3","F4","F5","PYT","AYT","Variety","meanGV")
  
  correlations <- as.data.frame(correlations)
  mean <- as.data.frame(rowMeans(correlations))
  Allcorrelations <- as.data.frame(cbind(mean,correlationsC3))
  rownames(correlationsC3) <- c("NewParents","F2","F3","F4","F5","PYT","AYT")
  colnames(correlationsC3) <- c(1:nReps+1)
  
  variances <- as.data.frame(variances)
  mean <- as.data.frame(rowMeans(variances))
  Allvariances <- as.data.frame(cbind(mean,variances))
  colnames(variancesC1) <- c(1:nReps+1)
  rownames(variancesC1) <- c("PrevCycPYT", "newParents","F1","F2", "F3","F4", "F5", "PYT","AYT")
  
  ##write files
  write.csv(Allgeneticvalues, "rrblup_rd_gvs_snp_yield.csv")
  write.csv(Allcorrelations, "rrblup_rd_cors_snp_yield.csv")
  write.csv(Allvariances,"rrblup_rd_vars_snp_yield.csv")
  saveRDS(alleles, file="rrblup_rd_alleles_snp_yield.rds")
  saveRDS(bv_ebv, file="rrblup_rd_bvebv_snp_yield.rds")

  
}
