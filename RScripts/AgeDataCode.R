# R script to work on Graeme Lloyd's local computer to get age data for each
# filtering surviving data set.

# Load metatree and Claddis libraries:
library(metatree)
library(Claddis)

# Ad hoc database of specimen-level OTUs:
SpecimenLevelOccurrences <- matrix(c("Hadrosaurinae_indet_LACM_CIT_2852", "Maastrichtian", "Maastrichtian", "72.1", "66", "Herpetocetus_sp_NMNS_PV19540", "Tabianian", "Tabianian", "5.332", "3.6", "Thalassotherii_indet_USNM_187416", "Langhian", "Langhian", "15.97", "13.63", "Xenophoridae_indet_ChM_PV4834", "Late Oligocene", "Late Oligocene", "28.4", "23.03", "Squalodontidae_indet_ChM_PV4991", "Chattian", "Chattian", "27.82", "23.03", "Burnetiamorpha_indet_TM_4305", "Guadalupian", "Guadalupian", "270.6", "260.4", "Burnetiamorpha_indet_NHMUK_R871", "Lopingian", "Lopingian", "260.4", "251.0", "Iguanodontia_indet_NHMUK_R1831", "Valanginian", "Valanginian", "139", "134", "cf_Mantellisaurus_sp_NHMUK_R3741", "Hauterivian", "Aptian", "134", "113", "Ceratosauria_indet_MSMN_V6235", "Bathonian", "Bathonian", "168", "166", "Abelisauroidea_indet_MPCM_13573", "Cenomanian", "Cenomanian", "100", "93.9", "Abelisauridae_indet_MPCA_56", "Maastrichtian", "Maastrichtian", "72.1", "66", "Rhynchosauria_indet_unreferred_Nova_Scotia", "Otischalkian", "Otischalkian", "228", "216.5", "Hyperodapedon_sp_GSI_Unreferred", "Norian", "Rhaetian", "228", "199.6", "Hyperodapedon_sp_Z_53_slash_3", "Norian", "Rhaetian", "228", "199.6", "Hyperodapedon_sp_USNM_494329", "Carnian", "Carnian", "228", "216.5", "Paratypothoracisini_indet_SMNS_19003", "Norian", "Norian", "228", "209", "Tetanurae_indet_MM_2and8to9and11to21", "Valanginian", "Valanginian", "139", "134", "Diplodocus_sp_CM_11161", "Oxfordian", "Tithonian", "161.2", "145.5", "Diplodocus_sp_CMC_VP14128", "Oxfordian", "Tithonian", "161.2", "145.5", "Plateosaurus_sp_GPIT_18392", "Bathonian", "Bathonian", "168", "166", "Sauropodomorpha_indet_SMNS_12216", "Norian", "Norian", "228", "209", "Abelisauridae_indet_MCF_dash_PVPH_dash_237", "Turonian", "Turonian", "93.5", "89.3", "Abelisauridae_indet_unreferred_La_Boucharde", "Campanian", "Campanian", "83.5", "70.6", "Abelisauridae_indet_unreferred_Pourcieux", "Campanian", "Campanian", "83.5", "70.6", "Abelisauroidea_indet_MCT_1783_dash_R", "Maastrichtian", "Maastrichtian", "72.1", "66", "Alvarezsauridae_indet_MPC_100_99and120", "Campanian", "Campanian", "83.5", "70.6", "Alvarezsauridae_indet_YPM_1049", "Maastrichtian", "Maastrichtian", "72.1", "66", "Alvarezsauridae_indet_MPC_100_99and120", "Campanian", "Campanian", "83.5", "70.6", "Alvarezsauridae_indet_YPM_1049", "Maastrichtian", "Maastrichtian", "72.1", "66", "Burnetiamorpha_indet_BP_1_7098", "Guadalupian", "Guadalupian", "270.6", "260.4", "Ceratosauria_indet_CCG_20011", "Aalenian", "Bajocian", "175.6", "167.7", "Ceratosauria_indet_MNN_tig6", "Bathonian", "Oxfordian", "168", "155.7", "Ceratosauria_indet_USNM_8415", "Kimmeridgian", "Kimmeridgian", "155.7", "150.8", "Enantiornithes_indet_FRDC_05_CM_004", "Aptian", "Aptian", "125", "112", "Enantiornithes_indet_FRDC_06_CM_012", "Aptian", "Aptian", "125", "112", "Enantiornithes_indet_FRDC_07_CM_001", "Aptian", "Aptian", "125", "112", "Enantiornithes_indet_RAM_14306", "Campanian", "Campanian", "83.5", "70.6", "Heterodontosauridae_indet_NHMUK_RU_A100", "Hettangian", "Sinemurian", "199.6", "189.6", "Ornithomimidae_indet_CMN_12068_and_12069_and_12070", "Maastrichtian", "Maastrichtian", "72.1", "66", "Pelagiarctos_sp_SDNHM_131041", "Middle Miocene", "Middle Miocene", "17.5", "15", "Sauropodomorpha_indet_ISI_R260", "Sinemurian", "Sinemurian", "196.5", "189.6", "Tetanurae_indet_MNN_GAD1_dash_2", "Aptian", "Albian", "125", "99.6", "Titanosauriformes_indet_MPEF_dash_PV_10606", "Toarcian", "Toarcian", "183", "175.6", "Tyrannosauroidea_indet_FMNH_PR_2750", "Aptian", "Aptian", "125", "112"), ncol = 5, byrow = TRUE, dimnames = list(c(), c("TaxonName", "FAD", "LAD", "Max", "Min")))

# Set working directory to NEXUS folder:
setwd("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/NEXUS")

# Get N states folder names (contains data sets used for Bayesian searches):
#NStatesFolders <- unique(unlist(lapply(as.list(list.dirs()), function(x) {y <- strsplit(x, "/")[[1]][2]; y <- y[!is.na(y)]; y[grep("[:0-9:]{1}", y)]})))

# Create empty data sets vector:
#DataSets <- vector(mode = "character")

# For each folder:
#for(i in NStatesFolders) {
  
  # Set working directory to ith folder:
  #setwd(paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/NEXUS", "/", i, sep = ""))
  
  # Add data sets found to vector:
  #DataSets <- sort(c(DataSets, unlist(lapply(as.list(list.files()), function(x) x[grep(".nex", x, fixed = TRUE)]))))
  
#}

# Build data sets from all NEXUS files instead of just the April ones:
DataSets <- gsub(".nex", "", setdiff(list.files(), gsub("\\.|/", "", list.dirs())), fixed = TRUE)

# For each independent data set:
for(i in DataSets) {
  
  # Read in XML file:
  XML <- readLines(paste("~/Documents/Homepage/www.graemetlloyd.com/xml/", gsub(".nex", "", i, fixed = TRUE), ".xml", sep =""))
  
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
  UpdatedReconNumbers <- unname(unlist(lapply(as.list(unlist(ReconNumbers)), function(x) {y <- PaleobiologyDBTaxaQuerier(x, original = FALSE)[, c("OriginalTaxonNo", "ResolvedTaxonNo")]; gsub("txn:|var:", "", y[!is.na(y)][1])})))
  
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
  write.table(AgesMatrix, paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/Ages/", gsub(".nex", "", i, fixed = TRUE), ".txt", sep = ""), row.names = FALSE)

  # Output loop position:
  cat(i, " ")

}
