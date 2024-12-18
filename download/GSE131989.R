
# GSE131989 ---------------------------------------------------------------------
# Cell-type-specific differential methylation in rheumatoid arthritis peripheral blood samples
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131989
# Different cell types 
# Note: The code below should work, however since this concerns a very large
# _RAW.tar file (> 3 gb), download was done via browser instead of the code below.

GSE <- "GSE131989"

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
getGEOSuppFiles(GSE,makeDirectory = TRUE,baseDir = dataDir) # Get files from GEO
unpack_GEO(dataDir, GSE) # Unpack 

# Pheno 
pheno <- pData(gse$GSE131989_series_matrix.txt.gz)
save(pheno,file=file.path(dataDir,GSE,"pheno_raw.Rdata")) # Save locally 

## Rename pheno --------------------------------------------------------------------- 
load(file.path(dataDir,GSE,"pheno_raw.Rdata"))  # pheno: 371 x 53

# Extracting basename
pheno$basename <- sub('.*/([^/]+)_\\w+\\.idat\\.gz', '\\1', pheno$supplementary_file.1)
rownames(pheno) <- pheno$basename

# Get the column names where all entries are equal (to remove)
to_remove <- names(pheno)[apply(pheno, 2, function(x) length(unique(x))) == 1]
pheno <- pheno %>% select(-all_of(to_remove)) 

# Remove columns starting with 'characteristics' as these are double
pheno <- pheno %>% select(-all_of(grep("^characteristics_", colnames(pheno), value = TRUE)))

names(pheno)[names(pheno) == "subjecttype:ch1"] <- "type"

# Recode the type column 
pheno <- pheno %>%
  mutate(type = case_when(
    type == "1" ~ "RA",
    type == "2" ~ "Control",
    TRUE ~ as.character(type)
  ))

# Extracting identifier numbers from title column
pheno$id <- as.numeric(gsub(".*\\s(\\d+)(?:\\s\\(.*\\))?", "\\1", pheno$title))

# Create a new column 'excluded' in pheno indicating whether the title contains '(excluded from processed data)'
pheno$excluded <- ifelse(grepl("\\(excluded from processed data\\)", pheno$title), 1, 0)

# Rename
names(pheno)[names(pheno) == "ageatblooddraw:ch1"] <- "age"
names(pheno)[names(pheno) == "sample_type:ch1"] <- "celltype"
names(pheno)[names(pheno) == "smokerever:ch1"] <- "smoker_ever"

# Change datatype to numeric if needed
pheno$age <- as.numeric(pheno$age)

# Remove some columns we don't need, manually defined
pheno <- pheno %>% select(-all_of(c("title", "description", "supplementary_file.1", "supplementary_file", "geo_accession","sentrix_id:ch1")))

# Save the trimmed pheno for QC
save(pheno,file=file.path(dataDir,GSE,"pheno.Rdata")) # Save locally


## EUTOPSQC ---------------------------------------------------------------------

preprocessData(input=file.path(dataDir, GSE),output=file.path(dataDir, GSE),report=file.path(dataDir, GSE),pheno=file.path(dataDir, GSE,"pheno.Rdata"))


## Check beta and pheno ---------------------------------------------------------------------
load(file.path(dataDir,GSE,"beta_merged.Rdata"))  # 434327  x  371
load(file.path(dataDir,GSE,"pheno.Rdata"))

ind <- match(rownames(pheno), colnames(beta_merged)) # all samples present
pheno <- pheno[ind ,]
identical(colnames(beta_merged), pheno$basename) # TRUE 