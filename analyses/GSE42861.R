# Libraries
library(dplyr)
library(EpiDISH)

# DATA LOCALLY STORED - change accordingly
path_beta <- "<your-path-beta>"
path_pheno <- "<your-path-pheno>"

# Load EpiDISH reference matrix
load(file.path(getwd(),"data/cent12CT.m.rda"))

# Load pre-processed Beta and Pheno (processed with EUTOPS pipeline)
load(path_beta) # beta_merged
load(path_pheno) # pheno

# Check if pheno and beta are in the same order
if(!identical(pheno$basename, colnames(beta_merged))){
  cat("WARNING: Basename pheno and beta are not in the same order")
}else{cat("Basename pheno and beta in the same order")}

# EpiDISH
load(file.path(getwd(),"0-data/cent12CT.m.rda"))
frac.m <- hepidish(beta.m = beta_merged,
                   ref1.m = centEpiFibIC.m,
                   ref2.m = cent12CT.m,
                   h.CT.idx = 3,
                   method = 'RPC',
                   maxit = 500)

# Selected CpGs (all)
cgs <- c("cg14507403",
         "cg09821790",
         "cg15694422",
         "cg27527887",
         "cg11754974",
         "cg16652347", "cg13828440","cg18637238")

present <- cgs[cgs %in% rownames(beta_merged)] # 4/8 are present (450k data)

# Subset the dataframe based on the matching row names
t <- beta_merged[present,] 

data_GSE42861 <- cbind(pheno, frac.m,t(t))
#score.sum <- apply(data_GSE42861[,c("cg11754974","cg16652347","cg13828440","cg18637238")], 1, sum) # Composite score of 4 CpGs (not all are available)
score.sum <- apply(data_GSE42861[,c("cg16652347","cg13828440","cg18637238")], 1, sum) # Composite score of 3 CpGs
data_GSE42861 <- cbind(data_GSE42861, score.sum)
data_GSE42861$dataset <- 'GSE42861'

# Save DF for plotting
save(data_GSE42861, file=file.path(getwd(),"data/GSE42861.Rdata"))