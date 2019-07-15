# R script to work on Graeme Lloyd's local computer to get age data for each
# filtering surviving data set.

# Load metatree and Claddis libraries:
library(metatree)
library(Claddis)

# Set working directory to NEXUS folder:
setwd("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/NEXUS")

# Get N states folder names (contains data sets used for Bayesian searches):
NStatesFolders <- unique(unlist(lapply(as.list(list.dirs()), function(x) {y <- strsplit(x, "/")[[1]][2]; y <- y[!is.na(y)]; y[grep("[:0-9:]{1}", y)]})))

# Create empty data sets vector:
DataSets <- vector(mode = "character")

# For each folder:
for(i in NStatesFolders) {
  
  # Set working directory to ith folder:
  setwd(paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/NEXUS", "/", i, sep = ""))
  
  # Add data sets found to vector:
  DataSets <- sort(c(DataSets, unlist(lapply(as.list(list.files()), function(x) x[grep(".nex", x, fixed = TRUE)]))))
  
}

# For each independent data set:
for(i in DataSets) {
  
  # Read in XML file:
  XML <- readLines(paste("~/Documents/Homepage/www.graemetlloyd.com/xml/", gsub(".nex", "", i, fixed = TRUE), ".xml", sep =""))
  
  # Isolate reconciliation numbers:
  ReconNumbers <- unlist(strsplit(unlist(lapply(as.list(XML[(grep("<Taxa number", XML) + 1):(grep("</Taxa>", XML) - 1)]), function(x) strsplit(x, "recon_no=\"|\">")[[1]][2])), ";"))
  
  # Get OTU names:
  OTUNames <- unlist(lapply(as.list(XML[(grep("<Taxa number", XML) + 1):(grep("</Taxa>", XML) - 1)]), function(x) strsplit(x, ">|<")[[1]][3]))
  
  # Get updated (reconciled) reconciliation numbers:
  UpdatedReconNumbers <- unname(unlist(lapply(as.list(ReconNumbers), function(x) {y <- PaleobiologyDBTaxaQuerier(x, original = FALSE)[, c("OriginalTaxonNo", "ResolvedTaxonNo")]; gsub("txn:|var:", "", y[!is.na(y)][1])})))
  
  # Build ages matrix:
  AgesMatrix <- do.call(rbind, lapply(lapply(apply(cbind(UpdatedReconNumbers, OTUNames), 1, as.list), unlist), function(x) {y <- PaleobiologyDBOccurrenceQuerier(x[1]); z <- matrix(nrow = 0, ncol = 3); if(sum(!is.na(y[, "MaxMa"])) > 0) z <- cbind(rep(x[2], sum(!is.na(y[, "MaxMa"]))), y[!is.na(y[, "MaxMa"]), c("MaxMa", "MinMa"), drop = FALSE]); z}))
  
  # For each OTU name:
  for(j in OTUNames) {
    
    # Get matrix of firsts and lasts for jth taxon:
    TaxonMatrix <- AgesMatrix[AgesMatrix[, 1] == j, c("MaxMa", "MinMa"), drop = FALSE]
    
    # Make numeric:
    TaxonMatrix <- cbind(as.numeric(TaxonMatrix[, 1]), as.numeric(TaxonMatrix[, 2]))
    
    # Prune out anything that is too young to possibly be the true FAD (younger or equal in age to top of oldest occurrence):
    TaxonMatrix <- TaxonMatrix[!max(TaxonMatrix[, 2]) > TaxonMatrix[, 1], , drop = FALSE]
    
    # Add abck to ages matrix with original jth taxon occurrences removed:
    AgesMatrix <- rbind(AgesMatrix[-which(AgesMatrix[, 1] == j), ], cbind(rep(j, nrow(TaxonMatrix)), TaxonMatrix))
    
  }
  
  # Check for missing taxa:
  if(length(setdiff(OTUNames, unique(AgesMatrix[, 1]))) > 0) stop("Taxa missing!")
  
  # Add column name for taxon:
  colnames(AgesMatrix)[1] <- "Taxon"
  
  # Write out to ages folder:
  write.table(AgesMatrix, paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/Ages/", gsub(".nex", "", i, fixed = TRUE), ".txt", sep = ""), row.names = FALSE)

  # Output loop position:
  cat(i, " ")

}

