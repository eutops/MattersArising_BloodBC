# Download TCGA-BRCA data
suppressPackageStartupMessages(library(TCGAbiolinks))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(here))

here::i_am("TCGA.R") ## replace this with final filename of the tcga download script
id <- 'TCGA-BRCA'
outDir <- '<your output file directory>/'

# Get & save pheno/query files ----
suppressMessages(query <- GDCquery(project = id,
                                   data.category = "DNA Methylation",
                                   platform = "Illumina Human Methylation 450",
                                   sample.type = c("Primary Tumor",
                                                   "Solid Tissue Normal"),
                                   data.type = "Methylation Beta Value"))
dat <- query$results[[1]]
query_pheno <- GDCquery_clinic(project = id,
                               type = "clinical")
ind <- match(dat$cases.submitter_id, query_pheno$submitter_id)
tcga_pheno <- query_pheno[ind,]
save(tcga_pheno, dat, file = paste0(outDir, substr(id, 6, nchar(id)), ".Rdata"))

# Get & save methylation files ----
cat("Starting download (", as.character(Sys.time()), ").\n", sep = "")
setwd(paste0(outDir, "/data/"))
GDCdownload(query, method = "client", files.per.chunk = 10)
suppressMessages(data <- GDCprepare(query))
beta <- assays(data)[[1]]
beta <- na.omit(beta)
cat("Done downloading, saving file.\n")
save(beta, file = paste0(outDir, "/data/beta.Rdata"))

# delete intermediate files
    unlink("GDCdata/", recursive = T)
    suppressMessages(file.remove(c("gdc-client_v1.6.1_OSX_x64.zip",
                                   "gdc_manifest.txt",
                                   "gdc-client",
                                   "gdc_client_configuration.dtt")))
    
# return wd to normal
setwd(here())
