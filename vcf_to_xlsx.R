################################################################################
##                                 Name
##                                 Date
##                               Institute
##        Extracting Onco Info from VCF files analyzed using Ion Reporter
################################################################################

# 1. load libraries ------------------------------------------------------------
library(tidyverse)
library(xlsx)

# 2. make a clean code for extraction------------------------------------

# 2.1 create an empty dataframe to store the combined data and an empty character vector
OCA_plus_analysis_data <- matrix(nrow = 0, ncol=28, byrow=TRUE)

# 2.2 create a vector of onco info to match from 1st file
colnames_vcf <- read_lines("/path/to/example_vcf_file.vcf", n_max = 57)
colnames_vcf <- gsub("##", "#", colnames_vcf[c(8,26,10:12,15:21,39:40,42:54,57)])
colnames_vcf <- gsub("^.|=.*", "", colnames_vcf)

# 3.3 define list of files
files <- list.files(path="./Input/", pattern="*.vcf", full.names=TRUE, recursive=FALSE) 

# 3.4 initialize the for loop
for (i in 1:length(files)) {
  # extract important onco info in a vector
  all_onco_info <- readLines(files[i], n = 57) # load file
  
  # check progress
  message('Processing file ', i, ' of ', length(files))
  
  # create an empty character vector
  onco_info <- character(0)
  
  # Loop through patterns in colnames_vcf
  for (pattern in colnames_vcf) {
    # Find lines that match the pattern
    matching_lines <- grep(paste0("##", pattern, "="), all_onco_info, value = TRUE, fixed = TRUE)
    
    # Check if matching lines were found
    if (length(matching_lines) > 0) {
      # Extract and append matching lines
      onco_info <- c(onco_info, matching_lines)
    } else {
      # If no matching lines found, append NA
      onco_info <- c(onco_info, paste0("##", pattern, "=NA"))
    }
  }
  
  # Remove special characters and make each iteration a data frame
  table <- matrix(gsub(".*=", "", onco_info), ncol = 28, byrow = TRUE)
  
  # rbind each table of onco info to the empty dataframe
  OCA_plus_analysis_data <- rbind(OCA_plus_analysis_data, table)
  
  # rename cols
  colnames(OCA_plus_analysis_data) <- colnames_vcf
}


# 3.6 write to excel file
write.xlsx(OCA_plus_analysis_data, "/path/to/outputfile/OCA_plus_analysis_data.xlsx", row.names=FALSE, col.names=TRUE)
write.csv(OCA_plus_analysis_data, "/path/to/outputfile/OCA_plus_analysis_data.csv", row.names=FALSE, col.names=TRUE)

