# GSE237036 ---------------------------------------------------------------------
# https://www.nature.com/articles/s41467-023-40389-5#Abs1
# Genome wide DNA methylation profiling of peripheral blood mononuclear cells (PBMCs) 
# in normal and BC samples. The Illumina Infinium 850k Human DNA methylation Beadchip 
# was used to obtain DNA methylation profiles across approximately 820,000 CpGs in PBMC samples. 
# Samples included 50 newly diagnosed BC patients and 30 normal controls.

GSE <- "GSE237036"

# Path where data will be downloaded to, locally
dataDir = "<path-to-directory>" 

# Libraries
library(GEOquery)
library(minfi)
library(EpiDISH)
library(ChAMP)
library(dplyr)
library(ggplot2)
library(eutopsQC)

# Functions
source(file.path(getwd(),"src/unpack_GEO.R"))

## Get files from GEO ---------------------------------------------------------------------
gse <- getGEO(GSE)

# Pheno ----------------

# Pheno original
pheno <- pData(gse$GSE237036_series_matrix.txt.gz) # Pheno file
fileName <- file.path(dataDir, paste0(GSE, "_pheno.Rdata"))
save(pheno, file=fileName)

# Create 'pheno_trimmed' dataframe for QC
pheno_trimmed <- pheno %>%
  # Selects only the case/control column
  select('disease state:ch1') %>%
  # Rename 'disease state:ch1' to 'type'
  rename(type = 'disease state:ch1') %>%
  # Extract the desired part from 'supplementary_file.1' aka basename
  mutate(basename = sub('.*?/(GSM\\d+_\\d+_R\\d+C\\d+).*', '\\1', pheno$supplementary_file.1))

fileName <- file.path(dataDir, paste0(GSE, "_pheno_trimmed.Rdata"))
save(pheno_trimmed, file=fileName)


# Retrieve files from GEO ----------------
files <- getGEOSuppFiles(GSE, makeDirectory = FALSE, baseDir = dataDir)

# Create the 'idats' directory
dir.create(file.path(dataDir, "idats"), showWarnings = FALSE)

setwd(file.path(dataDir, "idats"))

# Untar the file
untar(file.path(dataDir, "GSE237036_RAW.tar"))

gz_files <- list.files(path = file.path(dataDir, "idats"), pattern = "\\.gz$", full.names = TRUE)

# Iteratively unpack the files  ----------------
for (fileNr in 1:length(gz_files)){
  gunzip(gz_files[fileNr], remove=TRUE)
}


# EUTOPS QC  ----------------

# Create the 'qc' directory
dir.create(file.path(dataDir, "qc"), showWarnings = FALSE)

phenoPath = file.path(dataDir, "GSE237036_pheno_trimmed.Rdata") # Pheno file
inputDir = file.path(dataDir,"idats")
outputDir = file.path(dataDir, "qc")


preprocessData(input=inputDir, output=outputDir, report=outputDir,pheno=phenoPath)

writeLines(capture.output(sessionInfo()), "~/Desktop/sessionInfo.txt")
