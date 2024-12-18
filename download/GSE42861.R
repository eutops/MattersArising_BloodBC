
# GSE42861 ---------------------------------------------------------------------
# Differential DNA methylation in Rheumatoid arthritis
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42861
# Genome-wide DNA methylation level was studied to determine whether Rheumatoid arthritis patients (cases) 
# has methylation differences comparing to normal controls in peripheral blood leukocytes (PBLs).
# We used Illumina HumanMethylation450 BeadChip array to determine the genome-wide DNA methylation 
# difference in PBLs from Rheumatoid arthritis patients (cases) and normal controls

GSE <- "GSE42861"

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
library(data.table)
library(stringr)

# Functions
source(file.path(getwd(),"src/unpack_GEO.R"))

## Get files from GEO ---------------------------------------------------------------------

# Manually downloaded the series matrix and put it in dataDir/GSE42861; pheno is extremely big so API fails.
unpack_GEO(dataDir, GSE) # Unpack 

# Define path to series matrix (pheno)
path.pheno <- file.path(dataDir,GSE,"GSE42861_series_matrix.txt")

# Unzips the file to .txt file
gunzip(paste0("path.pheno",".gz"), remove=FALSE)

# Now open the .txt file in Notepad++ (download from internet first) and select ONLY rows 36 up to and including 69
# Paste this into a new .txt file that you will save in the same folder, as "GSE42861_series_matrix_trimmed.txt". 

path.pheno <- file.path(dataDir,GSE,"GSE42861_series_matrix_trimmed.txt")

# Read in the trimmed .txt file
pheno <- fread(path.pheno)

# Transpose the dataframe
pheno <- t(pheno)

# Step 1: Remove the prefix "!Sample_" from the entries in the first row
pheno[1,] <- sub("^!Sample_", "", pheno[1,])

# Step 2: Set the first row as column names
colnames(pheno) <- pheno[1,]

# Step 3: Remove the first row from the dataframe
pheno <- pheno[-1,]

# Setp 4: convert to dataframe
pheno <- as.data.frame(pheno)


# Get unique column names (. suffix)
colnames(pheno) <- make.unique(names(pheno), sep = ".")

# Pheno
save(pheno,file=file.path(dataDir,GSE,"pheno_raw.Rdata")) # Save locally


## Rename and trim pheno ---------------------------------------------------------------------
load(file.path(dataDir,GSE,"pheno_raw.Rdata"))  # pheno: 689 x 33

# Extracting basename
pheno$basename <- sub('.*/([^/]+)_\\w+\\.idat\\.gz', '\\1', pheno$supplementary_file)
rownames(pheno) <- pheno$basename

# Get the column names where all entries are equal (to remove)
to_remove <- names(pheno)[apply(pheno, 2, function(x) length(unique(x))) == 1]
pheno <- pheno %>% select(-all_of(to_remove)) # pheno: 689 x 9

# Rename columns
names(pheno)[names(pheno) == "characteristics_ch1.1"] <- "type"
names(pheno)[names(pheno) == "characteristics_ch1.3"] <- "age"
names(pheno)[names(pheno) == "characteristics_ch1.4"] <- "gender"
names(pheno)[names(pheno) == "characteristics_ch1.5"] <- "smoking"

# Recode the type column
pheno <- pheno %>%
  mutate(type = case_when(
    type == "disease state: Normal" ~ "Control",
    type == "disease state: rheumatoid arthritis" ~ "RA",
    TRUE ~ as.character(type)
  ))

# Extract numbers from the 'age' column
pheno$age <- as.numeric(str_extract(pheno$age, "\\d+"))

# Extract gender from the 'gender' column
pheno$gender <- str_extract(pheno$gender, "[fm]")

# Extract smoking status from the 'smoking' column
pheno$smoking <- str_extract(pheno$smoking, "(?<=smoking status: )\\w+")

# Remove some columns we don't need, manually defined
pheno <- pheno %>% select(-all_of(c("characteristics_ch1.2", "supplementary_file.1", "supplementary_file", "geo_accession")))

# Save the trimmed pheno for QC
save(pheno,file=file.path(dataDir,GSE,"pheno.Rdata")) # Save locally


## EUTOPSQC ---------------------------------------------------------------------

preprocessData(input=file.path(dataDir, GSE),output=file.path(dataDir, GSE),report=file.path(dataDir, GSE),pheno=file.path(dataDir, GSE,"pheno.Rdata"))


## Check beta and pheno ---------------------------------------------------------------------
load(file.path(dataDir,GSE,"beta_merged.Rdata"))  # 434175  x  689
load(file.path(dataDir,GSE,"pheno.Rdata")) # 689 x  5

identical(colnames(beta_merged), pheno$basename) # TRUE


