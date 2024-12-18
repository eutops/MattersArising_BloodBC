# Methylation QTL analysis
library(MatrixEQTL)

# Helper function
source(file.path(getwd(),"src/coerce_numeric.R"))

# Libraries
if(!require("vcfR")){
  BiocManager::install("vcfR")
}

if(!require("MatrixEQTL")){
  install.packages("MatrixEQTL")
}


# path_vcf <- file.path(eutopsPath,"data/raw-data/snp/20211104_mqtl_forecee_snp/PLINK/merged/3c_after_qc_fixref.vcf.gz")
path_vcf <- "path to VCF file" # This file contains genetic information and is available on request.

# Read in VCF file and create a matrix from it --------
x <- vcfR::read.vcfR(path_vcf)
y <- vcfR::extract.gt(x) # extract genotype (0/1, 0/0 or 1/1)
z <- as.matrix(y) 
# Genotypes currently coded as 0/0, 0/1, ... -> need to code to 0, 1, 2:
# 0/0 = 0
# 0/1 or 1/0 = 1
# 1/1 = 2
zz <- gsub("0/0", "0", z, perl = T)
zz <- gsub("0/1|1/0", "1", zz, perl = T)
zz <- gsub("1/1", "2", zz, perl = T)
z <- coerce_numeric(zz)
rm(zz, y, x);gc()

# Find samples that are overlapping between SNPs (FR/ST) and methylation in mdat
# path_pheno_file <- file.path(eutopsPath,"data/raw-data/snp/20211104_mqtl_forecee_snp/3C_External_Validation_SNPs_051121.xlsx")
path_pheno_file <- "path to pheno file"
# path_mdat <- file.path(eutopsPath,"eca/0-master-dataset/mdat.Rdata")
path_mdat <- "path to master dataset"

pheno_file <- readxl::read_xlsx(path_pheno_file) |> 
  dplyr::filter(!is.na(`sample ID`))# this has a column for ST number

load(path_mdat)
mdat <- mdat[mdat$experiment_name=="3C_VALIDATION_DATA_BLOOD" & mdat$sampletype == "blood",]
intersect1 <- intersect(mdat$ST_number, pheno_file$ST_number)
pheno_file <- pheno_file[match(intersect1, pheno_file$ST_number),]
mdat <- mdat[match(intersect1, pheno_file$ST_number),]

# Samples overlapping between pheno file and SNP matrix
# need to remove first item from colnames of the dataframe e.g. like this:
colnames(z) <- stringr::str_split(colnames(z), "_", simplify = T)[,2]
intersect <- intersect(colnames(z), pheno_file$`sample ID`) # remove any NAs and keep those that arepresent, append ST number, ...
z <- z[,intersect]
pheno_file <- pheno_file[match(intersect, pheno_file$`sample ID`),]
identical(colnames(z), pheno_file$`sample ID`)

# lastly, keep overlapping methylation
intersect <- intersect(pheno_file$ST_number, mdat$ST_number)
mdat <- mdat[match(intersect, mdat$ST_number),]

# Check overlaps
identical(pheno_file$ST_number, mdat$ST_number)
# Rename column names in snp matrix with basename for simplicity
colnames(z) <- rownames(mdat)

# Load beta and keep only relevant cg and samples
load("<your-path-forecee-beta>") # beta_merged
beta <- coerce_numeric(beta_merged["cg14507403",colnames(z)])

# Check colnames
identical(colnames(z), colnames(beta))

# MQTL analysis


# SNP data
# Create a SlicedData variable
snp = SlicedData$new();
# Show the details of the empty object
show(snp)
# Create a matrix of values and assign to sd
snp$CreateFromMatrix(z);
# Show the detail of the object (one slice)
show(snp);


# Methylation data
# Create a SlicedData variable
methyl = SlicedData$new();
# Show the details of the empty object
show(methyl)
# Create a matrix of values and assign to sd
methyl$CreateFromMatrix(beta);
# Show the detail of the object (one slice)
show(methyl);

out <- Matrix_eQTL_main(snp, 
                        methyl, 
                        cvrt = SlicedData$new(), 
                        output_file_name = "", 
                        pvOutputThreshold = 1e-5,
                        useModel = modelLINEAR, 
                        errorCovariance = numeric(), 
                        verbose = TRUE, 
                        output_file_name.cis = "", 
                        pvOutputThreshold.cis = 0,
                        snpspos = NULL, 
                        genepos = NULL,
                        cisDist = 1e6,
                        pvalue.hist = FALSE,
                        min.pv.by.genesnp = FALSE,
                        noFDRsaveMemory = FALSE)
pheno <- mdat

# Get mqtl IDs
mqtls <- out$all$eqtls$snps
pheno <- cbind(mdat, t(beta),
               t(z[mqtls,]))

# Plotting
pheno_long <- pheno |>
  tidyr::pivot_longer(cols = any_of(mqtls),
                      names_to = "mqtl",
                      values_to = "genotype")

save(pheno_long, file = file.path(getwd(),"data/snp.Rdata"))

