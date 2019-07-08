# R script to work on Graeme Lloyd's local computer to perform initial filter of
# data sets to include in the analysis.

# Load metatree and Claddis libraries:
library(metatree)
library(Claddis)

# Get independent data sets (excludes MRP, trace fossil etc. data sets):
IndependentDatasets <- GetPseudoIndependentPhyloFiles(xmlwd = "~/Documents/Homepage/www.graemetlloyd.com/xml", exclude.list = paste(c("Averianov_inpressa", "Bravo_et_Gaete_2015a", "Brocklehurst_etal_2013a", "Brocklehurst_etal_2015aa", "Brocklehurst_etal_2015ab", "Brocklehurst_etal_2015ac", "Brocklehurst_etal_2015ad", "Brocklehurst_etal_2015ae", "Brocklehurst_etal_2015af", "Bronzati_etal_2012a", "Bronzati_etal_2015ab", "Brusatte_etal_2009ba", "Campbell_etal_2016ab", "Carr_et_Williamson_2004a", "Carr_etal_2017ab", "Frederickson_et_Tumarkin-Deratzian_2014aa", "Frederickson_et_Tumarkin-Deratzian_2014ab", "Frederickson_et_Tumarkin-Deratzian_2014ac", "Frederickson_et_Tumarkin-Deratzian_2014ad", "Garcia_etal_2006a", "Gatesy_etal_2004ab", "Grellet-Tinner_2006a", "Grellet-Tinner_et_Chiappe_2004a", "Grellet-Tinner_et_Makovicky_2006a", "Knoll_2008a", "Kurochkin_1996a", "Lopez-Martinez_et_Vicens_2012a", "Lu_etal_2014aa", "Norden_etal_inpressa", "Pisani_etal_2002a", "Ruiz-Omenaca_etal_1997a", "Ruta_etal_2003ba", "Ruta_etal_2003bb", "Ruta_etal_2007a", "Selles_et_Galobart_2016a", "Sereno_1993a", "Sidor_2001a", "Skutschas_etal_inpressa", "Tanaka_etal_2011a", "Toljagic_et_Butler_2013a", "Tsuihiji_etal_2011aa", "Varricchio_et_Jackson_2004a", "Vila_etal_2017a", "Wilson_2005aa", "Wilson_2005ab", "Zelenitsky_et_Therrien_2008a"), ".xml", sep = ""))

# For each independent data set:
for(i in IndependentDatasets) {
  
  # Read in XML file:
  XML <- readLines(paste("~/Documents/Homepage/www.graemetlloyd.com/xml/", i, ".xml", sep =""))
  
  # If there are unreconciled taxa then remove these from the list:
  if(length(grep("recon_no=\"-1\"", XML)) > 0) IndependentDatasets <- setdiff(IndependentDatasets, i)
  
}

# For each independent data set:
for(i in IndependentDatasets) {
  
  # Read in XML file:
  XML <- readLines(paste("~/Documents/Homepage/www.graemetlloyd.com/xml/", i, ".xml", sep =""))
  
  # Isolate reconciliation numbers:
  ReconNumbers <- unlist(strsplit(unlist(lapply(as.list(XML[(grep("<Taxa number", XML) + 1):(grep("</Taxa>", XML) - 1)]), function(x) strsplit(x, "recon_no=\"|\">")[[1]][2])), ";"))
  
  # Get updated (reconciled) reconciliation numbers:
  UpdatedReconNumbers <- unname(unlist(lapply(as.list(ReconNumbers), function(x) {y <- PaleobiologyDBTaxaQuerier(x, original = FALSE)[, c("OriginalTaxonNo", "ResolvedTaxonNo")]; gsub("txn:|var:", "", y[!is.na(y)][1])})))
  
  # Build ages matrix:
  AgesMatrix <- do.call(rbind, lapply(as.list(ReconNumbers), function(x) {y <- PaleobiologyDBOccurrenceQuerier(x); z <- matrix(nrow = 0, ncol = 3); if(sum(!is.na(y[, "MaxMa"])) > 0) z <- cbind(rep(x, sum(!is.na(y[, "MaxMa"]))), y[!is.na(y[, "MaxMa"]), c("MaxMa", "MinMa"), drop = FALSE]); z}))
  
  # Remove any data sets where there is not at least one age for each tip:
  if(length(setdiff(UpdatedReconNumbers, unique(AgesMatrix[, 1]))) > 0) IndependentDatasets <- setdiff(IndependentDatasets, i)

  # Output loop position:
  cat(i, " ")

}

# For each independent data set:
for(i in IndependentDatasets) {
  
  # Read in XML file:
  XML <- readLines(paste("~/Documents/Homepage/www.graemetlloyd.com/xml/", i, ".xml", sep =""))
  
  # Isolate reconciliation numbers:
  ReconNumbers <- unlist(strsplit(unlist(lapply(as.list(XML[(grep("<Taxa number", XML) + 1):(grep("</Taxa>", XML) - 1)]), function(x) strsplit(x, "recon_no=\"|\">")[[1]][2])), ";"))
  
  # Get updated (reconciled) reconciliation numbers:
  UpdatedReconNumbers <- unname(unlist(lapply(as.list(ReconNumbers), function(x) {y <- PaleobiologyDBTaxaQuerier(x, original = FALSE)[, c("OriginalTaxonNo", "ResolvedTaxonNo")]; gsub("txn:|var:", "", y[!is.na(y)][1])})))
  
  # Build ages matrix:
  AgesMatrix <- do.call(rbind, lapply(as.list(ReconNumbers), function(x) {y <- PaleobiologyDBOccurrenceQuerier(x); z <- matrix(nrow = 0, ncol = 3); if(sum(!is.na(y[, "MaxMa"])) > 0) z <- cbind(rep(x, sum(!is.na(y[, "MaxMa"]))), y[!is.na(y[, "MaxMa"]), c("MaxMa", "MinMa"), drop = FALSE]); z}))
  
  # If there are fewer than five ages in total remove those data sets from the pool:
  if(length(unique((as.numeric(AgesMatrix[, "MaxMa"]) + as.numeric(AgesMatrix[, "MinMa"])) / 2)) < 5) IndependentDatasets <- setdiff(IndependentDatasets, i)
  
  # Output loop position:
  cat(i, " ")
  
}

# Independent data sets from above as found on 8/7/19:
IndependentDatasets <- c("Aguirre-Fernandez_etal_2009a", "Amson_et_de_Muizon_2014a", "Arnold_etal_2005a", "Bandeira_etal_2016a", "Barbosa_etal_2008a", "Bell_et_Evans_2010aa", "Berta_1991aa", "Bhullar_2011a", "Bianucci_2013a", "Bisconti_inpressa", "Boessenecker_et_Churchill_2013aa", "Bolotsky_et_Godefroit_2004a", "Bona_et_de_la_Fuente_2005a", "Boy_1989a", "Brinkman_et_Nicholls_1991a", "Brochu_1997aa", "Burns_et_Currie_2014ac", "Burns_et_Currie_2014ad", "Burns_etal_2011a", "Butler_etal_2007a", "Carballido_etal_2017bb", "Carpenter_2001aa", "Carpenter_2001ab", "Carpenter_2001ac", "Carpenter_2001ad", "Carpenter_etal_1998a", "Casanovas_etal_2001a", "Charig_et_Milner_1997a", "Chatterjee_1998a", "Chatterjee_2002a", "Chiappe_1995a", "Chiappe_1996a", "Chiappe_et_Calvo_1994a", "Chiappe_etal_1996a", "Cisneros_etal_2012a", "Clos_1995a", "Cullen_etal_2013a", "Dalla_Vecchia_2009a", "Day_etal_2016a", "Delcourt_et_Iori_inpressa", "Demar_2013a", "Druckenmiller_et_Maxwell_2010a", "Ekdale_etal_2011a", "Elzanowski_1999a", "Evander_1989a", "Ezcurra_etal_2016a", "Farke_et_Patel_2012a", "Farke_etal_2014ba", "Farke_etal_2014bb", "Fernandez_et_Talevi_2014a", "Ferreira_etal_2015a", "Fischer_etal_2016a", "Gallina_et_Otero_2015a", "Gasulla_etal_2015ab", "Gatesy_et_Amato_1992a", "Gatesy_etal_1993a", "Gay_2010a", "Geisler_et_Luo_1996a", "Gerlach_2001a", "Godefroit_etal_1998a", "Godefroit_etal_2004b", "Gower_et_Sennikov_1997a", "Han_etal_2015a", "Harper_etal_2000a", "Hill_etal_2003a", "Kammerer_2016b", "Kimura_et_Hasegawa_2010a", "Kutty_etal_2007ab", "Lambert_2008a", "Lee_1996a", "Liu_et_Rieppel_2005a", "Liu_etal_2013a", "Lu_etal_2002a", "Madzia_etal_inpressab", "Maganuco_etal_2007a", "Maidment_etal_2006a", "Maisch_et_Gebauer_2005a", "Makovicky_et_Sues_1998a", "Marek_etal_2015a", "Martinez_etal_2016aa", "Martinez_etal_2016ab", "Milner_etal_2009a", "Monks_et_Owen_2000a", "Mukherjee_et_Ray_2014a", "Norell_etal_2006ab", "Norman_2014a", "Novas_1996a", "Novas_1997a", "Parker_2016a", "Parker_2018a", "Parsons_et_Parsons_2009a", "Penkalski_2014a", "Penkalski_et_Dodson_1999a", "Pereda-Suberbiola_etal_2009a", "Perez-Moreno_etal_1993a", "Perez-Moreno_etal_1994a", "Rauhut_et_Carrano_inpressa", "Rauhut_et_Xu_2005a", "Russell_et_Dong_1993a", "Russell_et_Zheng_1993a", "Ryan_2007ab", "Sallam_etal_inpressaa", "Sallam_etal_inpressab", "Sanz_et_Buscalioni_1992a", "Schoch_2006a", "Schoch_et_Milner_2008a", "Schoch_etal_2007a", "Schott_et_Evans_inpressa", "Sereno_1999ah", "Sereno_2010a", "Sereno_etal_1996a", "Sereno_etal_1998a", "Shibata_etal_2015a", "Stocker_etal_inpressa", "Sues_1997a", "Takasaki_etal_2018a", "Thewissen_1992a", "Thomas_2015a", "Trueman_1998aa", "Tsogtbaatar_etal_2014a", "Tsuihiji_etal_2011ab", "Tsuji_etal_2013b", "Vickaryous_etal_2001a", "Vincent_etal_2012a", "Winkler_etal_1997a", "Woodruff_etal_2018ac", "Wu_etal_1996a", "Xu_etal_1999b", "Xu_etal_2011d", "Xu_etal_2012bb", "Xu_etal_2013a", "Xu_etal_inpressa", "Yates_2003a", "You_2003a", "You_et_Dodson_2003a", "You_et_Dodson_2004a", "Zanno_et_Makovicky_2011ab", "Zhang_et_Zhou_2000a")

# For each independent data set:
for(i in IndependentDatasets) {
  
  # Read in XML file:
  XML <- readLines(paste("/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/xml/", i, ".xml", sep = ""))
  
  # Write XML file to github folder:
  write(XML, paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/XML/", i, ".xml", sep = ""))
  
  # Read in NEXUS file:
  NEXUS <- Claddis::ReadMorphNexus(paste("/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/nexus/", i, ".nex", sep = ""))
  
  # Write NEXUS file to github folder:
  Claddis::WriteMorphNexus(NEXUS, paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/NEXUS/", i, ".nex", sep = ""))
  
  # Form MPTs file name:
  MPTsFileName <- paste(i, ".tre", sep = "")
  
  # Form MPTS file name if zipped:
  MPTsZipFileName <- paste(i, ".tre.zip", sep = "")
  
  # Set working directory to mpts folder:
  setwd("/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/mpts/")
  
  # Update i to MPTs file name:
  i <- ifelse(any(list.files() == MPTsFileName), MPTsFileName, MPTsZipFileName)
  
  # Copy file to github folder:
  file.copy(from = paste("/Users/eargtl/Documents/Homepage/www.graemetlloyd.com/mpts/", i, sep = ""), to = paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/MPTs/", i, sep = ""))
  
}

# Set working directory to github mpts:
setwd("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/MPTs")

# Find any zipped tree files:
ZipFiles <- list.files()[grep(".tre.zip", list.files(), fixed = TRUE)]

# If zip files found:
if(length(ZipFiles) > 0) {
  
  # Unzip zipped files:
  lapply(as.list(ZipFiles), function(x) unzip(x))
  
  # Remove zipped version:
  lapply(as.list(ZipFiles), function(x) file.remove(x))
  
}

# For each MPTs file:
for(i in list.files()) {
  
  # Read in trees:
  Trees <- ape::read.tree(i)
  
  # If more than 1000 trees then randomly sample only 1000 of them:
  if(length(Trees) > 1000) Trees <- Trees[sample(1:1000)]
  
  # Write trees back to MPTs folder:
  write.tree(Trees, i)
  
  # Output loop position:
  cat(i, " ")
  
}

# SPECIES-LEVEL OTU DATA SETS ONLY? HOW TO CHECK?
