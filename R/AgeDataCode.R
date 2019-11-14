#' Get Age data from PBDB to pair to the tips of trees
#'
#' For each nexus dataset, locate an XML file containing ages. If any are missing, query PBDB to see if they can be found. Returns a TSV-file of ages per tip.
#' @import Claddis
#' @import metatree
#' @param data_path A path to nexus datafiles
#' @param xml_path A path to xml metadata for each dataset
#' @param outpath Where you'd like to write the data to. 
#' @return AgesMatrix A dataframe of taxa and their minimum and maximum occurrence
#' @export



build_ages <- function(data_path, xml_path, out_path){
# Build data sets from all NEXUS files instead of just the April ones:
  DataSets <- list.files(data_path, pattern = '.nex')
  print(DataSets)
# For each independent data set:
  for(i in DataSets) {
  print(i)
  # Read in XML file:
  XML <- readLines(paste(xml_path, gsub(".nex", "", i, fixed = TRUE), ".xml", sep =""))
  print(XML)

  # Isolate reconciliation numbers as list:
  ReconNumbers <- strsplit(unlist(lapply(as.list(XML[(grep("<Taxa number", XML) + 1):(grep("</Taxa>", XML) - 1)]), function(x) strsplit(x, "recon_no=\"|\">")[[1]][2])), ";")
  
  # Isolate reconciliation numbers as list:
  ReconNames <- strsplit(unlist(lapply(as.list(XML[(grep("<Taxa number", XML) + 1):(grep("</Taxa>", XML) - 1)]), function(x) strsplit(x, "recon_name=\"|\"")[[1]][2])), ",")
  
  # Get OTU names:
  OTUNames <- unlist(lapply(as.list(XML[(grep("<Taxa number", XML) + 1):(grep("</Taxa>", XML) - 1)]), function(x) strsplit(x, ">|<")[[1]][3]))
  
  # Add any specimen-level OTUs found to vector:
  SpecimenLevelOTUs <- unlist(ReconNames)[unlist(lapply(strsplit(unlist(ReconNames), split = "_"), length)) > 2]
  
  # If specimen-level OTUs were found:
  if(length(SpecimenLevelOTUs) > 0) {
    
    # Find any missing from the database:
    MissingSpecimenLevelOTUs <- sort(unlist(lapply(as.list(SpecimenLevelOTUs), function(x) ifelse(!any(SpecimenLevelOccurrences[, "TaxonName"] == x), x, NA))))
    
    # If found stop and warn user:
    if(length(MissingSpecimenLevelOTUs) > 0) stop(paste("The following specimen-level OTUs are missing from the database: ", paste(MissingSpecimenLevelOTUs, collapse = ", "), ".", sep = ""))
    
  }
  
  # Get updated (reconciled) reconciliation numbers:
  UpdatedReconNumbers <- unname(unlist(lapply(as.list(unlist(ReconNumbers)), function(x) {y <- metatree::PaleobiologyDBTaxaQuerier(x, original = FALSE)[, c("OriginalTaxonNo", "ResolvedTaxonNo")]; gsub("txn:|var:", "", y[!is.na(y)][1])})))
  
  # Build new lists from recon numbers:
  UpdatedReconNumbersList <- OTUNamesList <- ReconNumbers
  
  # For each value in order:
  for(j in 1:length(UpdatedReconNumbers)) {
    
    # Find list position for jth value:
    ListPosition <- ListPosition2ListCoordinates(j, unlist(lapply(ReconNumbers, length)))
    
    # Fill ou updated recon numbers list with updated recon numbers:
    UpdatedReconNumbersList[[ListPosition$ListNumber]][ListPosition$ListPosition] <- UpdatedReconNumbers[j]
    
    # Fill out OTU list with OTU names:
    OTUNamesList[[ListPosition$ListNumber]][ListPosition$ListPosition] <- OTUNames[ListPosition$ListNumber]
    
  }
  
  # Build ages matrix:
  AgesMatrix <- do.call(rbind, lapply(lapply(apply(cbind(unlist(UpdatedReconNumbersList), unlist(OTUNamesList), unlist(ReconNames)), 1, as.list), unlist), function(x) {y <- PaleobiologyDBOccurrenceQuerier(x[1]); z <- matrix(nrow = 0, ncol = 3); if(sum(!is.na(y[, "MaxMa"])) > 0) z <- cbind(rep(x[2], sum(!is.na(y[, "MaxMa"]))), rep(x[3], sum(!is.na(y[, "MaxMa"]))), y[!is.na(y[, "MaxMa"]), c("MaxMa", "MinMa"), drop = FALSE]); z}))
  
  # Remove any extraneous row names:
  rownames(AgesMatrix) <- NULL
  
  # Add column name for taxon:
  colnames(AgesMatrix)[1:2] <- c("Taxon", "ReconName")
  
  # If specimen-level OTUs were found:
  if(length(SpecimenLevelOTUs) > 0) {
    
    # Find any missing from the database:
    MissingSpecimenLevelOTUs <- sort(unlist(lapply(as.list(SpecimenLevelOTUs), function(x) ifelse(!any(SpecimenLevelOccurrences[, "TaxonName"] == x), x, NA))))
    
    # If found stop and warn user:
    if(length(MissingSpecimenLevelOTUs) > 0) stop(paste("The following specimen-level OTUs are missing from the database: ", paste(MissingSpecimenLevelOTUs, collapse = ", "), ".", sep = ""))
    
  }
  
  # Remove recon name column:
  AgesMatrix <- AgesMatrix[, c("Taxon", "MaxMa", "MinMa")]
  
  # For each OTU name:
  for(j in OTUNames) {
    
    # Get matrix of firsts and lasts for jth taxon:
    TaxonMatrix <- AgesMatrix[AgesMatrix[, 1] == j, c("MaxMa", "MinMa"), drop = FALSE]
    
    # Make numeric:
    TaxonMatrix <- cbind(as.numeric(TaxonMatrix[, 1]), as.numeric(TaxonMatrix[, 2]))
    
    # Prune out anything that is too young to possibly be the true FAD (younger in age to top of oldest occurrence):
    TaxonMatrix <- TaxonMatrix[!max(TaxonMatrix[, 2]) > TaxonMatrix[, 1], , drop = FALSE]
    
    # Add abck to ages matrix with original jth taxon occurrences removed:
    AgesMatrix <- rbind(AgesMatrix[-which(AgesMatrix[, 1] == j), ], cbind(rep(j, nrow(TaxonMatrix)), TaxonMatrix))
    
  }
  
  # Check for missing taxa:
  if(length(setdiff(OTUNames, unique(AgesMatrix[, 1]))) > 0) stop("Taxa missing!")
  
  # Write out to ages folder:
  write.table(AgesMatrix, paste( out_path,"/", gsub(".nex", "", i, fixed = TRUE), ".txt", sep = ""), row.names = FALSE)

  # Output loop position:
  cat(i, " ")
  }
  return(AgesMatrix)
}
