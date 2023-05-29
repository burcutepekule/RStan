rm(list=ls())
setwd('/Users/burcutepekule/Library/CloudStorage/Dropbox/criticalwindow/code/R/RStan')
library(httr)

# Specify the ENA accession number
accession_number <- "PRJEB38984"

# Make an HTTP request to retrieve the data
response <- GET(paste0("https://www.ebi.ac.uk/ena/browser/api/fasta/", accession_number))


