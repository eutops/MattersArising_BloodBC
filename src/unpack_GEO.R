unpack_GEO <- function(dataDir, GSE, grab="_RAW.tar$") {
  # Set path to data
  setwd(file.path(dataDir, GSE))
  
  # Get the list of files in the directory
  files <- list.files()
  
  # Get the file that ends with "_RAW.tar"
  tarFile <- files[grepl(grab, files)]
  
  # Unzip the raw data
  untar(tarFile)
  
  # Get the updated list of files
  files <- list.files()
  
  # Filter out file names that end with "_RAW.tar"
  filtered_files <- files[!grepl("_RAW.tar$", files)]
  
  # Filter out file names that do not end with ".gz"
  gz_files <- filtered_files[grep("\\.gz$", filtered_files)]
  
  # Print the filtered list of .gz file names for visual check
  print(gz_files)
  
  # Iteratively unpack the files
  for (fileNr in 1:length(gz_files)) {
    gunzip(gz_files[fileNr], remove = TRUE)
  }
}


