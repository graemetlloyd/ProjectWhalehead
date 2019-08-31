# R script to perform strat congruence tests. All very slow!

# Load metatree and strap libraries:
library(hypRspace) # Install from github: devtools::install_github("graemetlloyd/hypRspace")
library(metatree) # Install from github: devtools::install_github("graemetlloyd/metatree")
library(strap)

# Local directory (modify for your local machine and hopefully rest should just work assuming the github directory structure matches):
TopDirectory <- "~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead"

# Set as Rev Bayes directory:
setwd(paste(TopDirectory, "/Data/RevBayesOutput", sep = ""))

# Find all Rev bayes files there:
RevBayesFiles <- list.files()

# For each Rev Bayes file:
for(i in RevBayesFiles) {
  
  # (Re)set to Rev Bayes directory:
  setwd(paste(TopDirectory, "/Data/RevBayesOutput", sep = ""))
  
  # Read in Rev bayes output:
  RevBayesOutput <- read.table(i, header = TRUE, stringsAsFactors = FALSE)
  
  # Set working directory as NEXUS data:
  setwd(paste(TopDirectory, "/Data/NEXUS", sep = ""))
  
  # Read in outgroup taxon:
  OutgroupTaxon <- rownames(Claddis::ReadMorphNexus(gsub(".trees", "", i, fixed = TRUE))$Matrix_1$Matrix)[1]
  
  # Get just the trees in ape format:
  BayesianTrees <- ape::read.tree(text = RevBayesOutput[, "tree"])
  
  # Remove branch lengths (allows collapsing to unique topologies below):
  BayesianTrees <- lapply(BayesianTrees, function(x) {x$edge.length <- NULL; x})
  
  # Reformat as a multiPhylo object:
  class(BayesianTrees) <- "multiPhylo"
  
  # Collapse to just unique topologies:
  BayesianTrees <- ape::read.tree(text = hypRspace::UniqueNewickStrings(ape::write.tree(BayesianTrees)))
  
  # Reroot trees on outgroup taxon:
  BayesianTrees <- root(BayesianTrees, outgroup = OutgroupTaxon)
  
  # Now reroot properly because ape is useless:
  BayesianTrees <- ape::read.tree(text = unlist(lapply(BayesianTrees, function(x) paste("(", OutgroupTaxon, ",", gsub(";", ");", ape::write.tree(drop.tip(x, tip = OutgroupTaxon))), sep = ""))))
  
  # Set working directory as Bayesian topologies:
  setwd(paste(TopDirectory, "/Data/BayesianTopologies", sep = ""))
  
  # If over 1000 trees randomly downsample to 1000:
  if(length(BayesianTrees) > 1000) BayesianTrees <- BayesianTrees[sample(1:length(BayesianTrees), size = 1000)]
  
  # Write out Bayesian topologies to folder:
  ape::write.tree(BayesianTrees, gsub(".nex.trees", ".tre", i, fixed = TRUE))
  
}

# Set working directory to bayesian topologies:
setwd(paste(TopDirectory, "/Data/BayesianTopologies", sep = ""))

# Make vector of file names found there:
BayesianTopologiesProcessed <- gsub(".tre", "", list.files(), fixed = TRUE)

# Show any unsampled data sets:
paste(setdiff(gsub(".nex.trees", "", RevBayesFiles, fixed = TRUE), BayesianTopologiesProcessed), collapse = " ")
