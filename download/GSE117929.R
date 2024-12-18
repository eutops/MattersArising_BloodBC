
# GSE117929 ---------------------------------------------------------------------
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117929
# Integration of genome-wide DNA methylation and transcription uncovered aberrant methylation-regulated genes and pathways in the peripheral blood mononuclear cells of systemic sclerosis [methylation]

GSE <- "GSE117929"

# Path where data will be downloaded to, locally
dataDir = "<path-to-directory>" 

# Libraries
library(GEOquery)
library(minfi)
library(EpiDISH)
library(ChAMP)
library(dplyr)
library(ggplot2)

# Functions
source(file.path(getwd(),"src/unpack_GEO.R"))

## Get files from GEO ---------------------------------------------------------------------

gse <- getGEO(GSE)
getGEOSuppFiles(GSE,makeDirectory = F,baseDir = file.path(dataDir,GSE)) # Get files from GEO
unpack_GEO(dataDir, GSE) # Unpack 

# Pheno 
pheno <- pData(gse$GSE117929_series_matrix.txt.gz)
save(pheno,file=file.path(dataDir,GSE,"pheno_raw.Rdata")) # Save locally 

# Check name of methylated / unmethylated signal
x <- readLines(file.path(dataDir,GSE,"GSE117929_methylated_and_unmethylated_signal.txt"),n=4)

# Filled in name of file, Uname, Mname, and sep manually after inspecting x
Mset <- minfi::readGEORawFile(filename = file.path(dataDir,GSE,"GSE117929_methylated_and_unmethylated_signal.txt"),
                              sep = "\t",
                              Uname = "Unmethylated signal",
                              Mname = "Methylated signal",
                              array = "IlluminaHumanMethylation450k",
                              annotation = "ilmn12.hg19",
                              showProgress = TRUE)

## Edit and trim pheno file --------------------------------------------------------------------- 
load(file.path(dataDir,GSE,"pheno_raw.Rdata"))

# 1. Rename columns
pheno$type <- NULL
pheno <- pheno %>%
  rename(type = description.1,
         diagnosis = `diagnosis:ch1`,
         gender = `gender:ch1`)
pheno$basename <- gsub(" \\[methylation\\]", "", pheno$title)

# 2. Select specific columns
pheno <- select(pheno, type, diagnosis, gender,description, geo_accession, basename)


save(pheno,file=file.path(dataDir,GSE,"pheno.Rdata")) # Save locally 


## QC & Pipeline --------------------------------------------------------------------- 


# Define global thresholds
INTENSITY_THRESHOLD <- 9.5     # minimum median intensity required
DETECTION_P_THRESHOLD <- 0.01  # maximum detection p-value
FAILED_PROBE_THESHOLD <- 0.1   # maximum proportion of failed probes per sample

qc <- getQC(Mset)
plotQC(qc) # All above the `good` threshold

# Filter any samples with median (un)methylated intensity less than threshold
low_intensity_samples <- rownames(qc)[qc$mMed<INTENSITY_THRESHOLD | qc$uMed<INTENSITY_THRESHOLD] # no samples below threshold

# Read in detP
colnames <- strsplit(readLines(file.path(dataDir,GSE,"GSE117929_methylated_and_unmethylated_signal.txt"), n = 1), "\t")[[1]]
select <- sort(grep("Detection Pval",colnames))
detP <- data.table::fread(file.path(dataDir,GSE,"GSE117929_methylated_and_unmethylated_signal.txt"),
                          sep = "\t",
                          select = select)

# Filter samples with too many failed probes
failed_samples <- colnames(detP)[colSums(detP>DETECTION_P_THRESHOLD) > (nrow(detP) * FAILED_PROBE_THESHOLD)] # No samples with low intensity fails
rm(detP, colnames, select, qc);gc()

# Convert to RatioSet and then beta
RSet <- ratioConvert(Mset, what = "both", keepCN = TRUE)
rm(Mset);gc()

beta <- getBeta(RSet)
beta <- na.omit(beta)
rm(RSet);gc()

densityPlot(beta) # Plot densities


### ChAMP Normalisation ---------------------------------------------------------------------
norm <- champ.norm(beta, arraytype = "450k", method = "BMIQ", cores = 8)
beta <- norm
densityPlot(beta)

## Rename pheno ---------------------------------------------------------------------
rownames(pheno) <- pheno$basename
ind <- match(rownames(pheno), colnames(beta)) # NOT all samples present


# Find entries that occur in both beta and pheno, since there are entries not in pheno and not in beta so we need to find their intersection.
common_samples <- intersect(colnames(beta), rownames(pheno)) # 31 samples

# Subset both beta and pheno
beta2 <- beta[,colnames(beta) %in% common_samples]
pheno2 <- pheno[rownames(pheno) %in% common_samples,]
pheno2 <- pheno2[match(colnames(beta2),rownames(pheno2)),]


identical(colnames(beta2), pheno2$basename) # TRUE

beta <- beta2
pheno <- pheno2

## Save subsetted beta and pheno ---------------------------------------------------------------------
save(beta, file = file.path(dataDir,GSE,"beta.Rdata"))
save(pheno, file = file.path(dataDir,GSE,"pheno.Rdata"))
identical(rownames(pheno), colnames(beta)) # TRUE

# Look at batch effects ~ UMAP ; prep pheno further - OPTIONAL --------
b <- beta[1:100000,]
umap <- uwot::umap(t(b))

ind <- match(rownames(umap), rownames(pheno))
pheno$umap1 <- umap[,1]
pheno$umap2 <- umap[,2]


## Plot ---------------------------------------------------------------------
pheno |>
  ggplot(aes(x = umap1,
             y = umap2,
             colour = type)) +
  geom_point()
