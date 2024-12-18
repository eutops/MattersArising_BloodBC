# Libraries
library(dplyr)
library(EpiDISH)

# DATA LOCALLY STORED - change accordingly
path_beta <- "<your-path-beta>"
path_pheno <- "<your-path-pheno>"

# Load EpiDISH reference matrix
load(file.path(getwd(),"data/cent12CT.m.rda"))

# Load pre-processed beta_merged_merged and Pheno (processed with EUTOPS pipeline)
load(path_beta) # beta_merged
load(path_pheno) # pheno

# Check if pheno and beta_merged are in the same order
if(!identical(pheno$basename, colnames(beta))){
  cat("WARNING: Basename pheno and beta_merged are not in the same order")
}else{cat("Basename pheno and beta_merged in the same order")}

# EpiDISH 
data(cent12CT.m)
frac.m <- hepidish(beta.m = beta,
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

present <- cgs[cgs %in% rownames(beta)] # 450k so only "cg15694422" "cg16652347" "cg13828440" "cg18637238"

# Subset the dataframe based on the matching row names
t <- beta[present,] # 4 x 31

data_GSE117929 <- cbind(pheno, frac.m,t(t))
#score.sum <- apply(data_GSE117929[,c("cg11754974","cg16652347","cg13828440","cg18637238")], 1, sum) # Composite score of 4 CpGs (in case of 850k data)
score.sum <- apply(data_GSE117929[,c("cg16652347","cg13828440","cg18637238")], 1, sum) # Composite score of 3 CpGs  (in case of 450k data)
data_GSE117929 <- cbind(data_GSE117929, score.sum)
data_GSE117929$dataset <- 'GSE117929'

# Save DF for plotting
save(data_GSE117929, file=file.path(getwd(),"data/GSE117929.Rdata"))