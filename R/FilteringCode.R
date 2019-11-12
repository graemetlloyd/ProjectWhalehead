##' Custom functions to filter datasets down to those that do not have unreconciled taxa
#'
#' Some datasets in PBDB have unreconciled taxa names. We remove these before 
#' proceding.
#' @param xml_path A directory of datasets that you would like to clean of unreconiled taxa. These must be xml metadata files.
#' @return vector A vector datasets with the unreconciled taxa removed
#'
#' @export remove_unreconciled

remove_unreconciled <- function(xml_path){
# For each independent data set:
  datasets <- list.files(xml_path, pattern = '.xml')
  for(i in datasets) {
  # Read in XML file:
    
    XML <- readLines(paste(xml_path, '/', i, sep =""))
  # If there are unreconciled taxa then remove these from the list:
    if(length(grep("recon_no=\"-1\"", XML)) > 0) datasets <- setdiff(datasets, i)
  }
  return(datasets)
}

##' Custom functions to filter datasets by completeness of tip sampling
#'
#' Some datasets in PBDB have unreconciled lack full tip sampling, or have fewer than five tips. These functions allow us to remove those datasets.
#' @param datasets A vector of datasets that you would like to clean of too-small datasets or datasets with tip ages missing. These must be xml metadata files.
#' @return vector A vector datasets with the datasets with too few tips or missing tips removed
#'
#' @export remove_missing_tips
#' @examples

remove_missing_tips <- function(xml_path){
# For each independent data set:
  datasets <- list.files(xml_path, pattern = '.xml')
  for(i in datasets) {
    # Read in XML file:
    
    XML <- readLines(paste(xml_path, '/', i, sep =""))
  # Isolate reconciliation numbers:
    ReconNumbers <- unlist(strsplit(unlist(lapply(as.list(XML[(grep("<Taxa number", XML) + 1):(grep("</Taxa>", XML) - 1)]), function(x) strsplit(x, "recon_no=\"|\">")[[1]][2])), ";"))
  
  # Get updated (reconciled) reconciliation numbers:
  UpdatedReconNumbers <- unname(unlist(lapply(as.list(ReconNumbers), function(x) {y <- PaleobiologyDBTaxaQuerier(x, original = FALSE)[, c("OriginalTaxonNo", "ResolvedTaxonNo")]; gsub("txn:|var:", "", y[!is.na(y)][1])})))
  
  # Build ages matrix:
  AgesMatrix <- do.call(rbind, lapply(as.list(ReconNumbers), function(x) {y <- PaleobiologyDBOccurrenceQuerier(x); z <- matrix(nrow = 0, ncol = 3); if(sum(!is.na(y[, "MaxMa"])) > 0) z <- cbind(rep(x, sum(!is.na(y[, "MaxMa"]))), y[!is.na(y[, "MaxMa"]), c("MaxMa", "MinMa"), drop = FALSE]); z}))
  
  # Remove any data sets where there is not at least one age for each tip:
  if(length(setdiff(UpdatedReconNumbers, unique(AgesMatrix[, 1]))) > 0) datasets <- setdiff(datasets, i)
  
  # If there are fewer than five ages in total remove those data sets from the pool:
  if(length(unique((as.numeric(AgesMatrix[, "MaxMa"]) + as.numeric(AgesMatrix[, "MinMa"])) / 2)) < 5) datasets <- setdiff(datasets, i)
  
  # Output loop position:
  cat(i, " ")
  }
  return(datasets)
}


