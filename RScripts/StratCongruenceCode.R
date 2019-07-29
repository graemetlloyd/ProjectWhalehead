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
  RevBayesOutput <- read.table(i, header = TRUE, stringsAsFactors=FALSE)
  
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
  BayesianTrees <- ape::read.tree(text = paste("(", OutgroupTaxon, ",", gsub(";", ");", ape::write.tree(drop.tip(BayesianTrees, tip = OutgroupTaxon))), sep = ""))
  
  # Set working directory as NEXUS data:
  setwd(paste(TopDirectory, "/Data/MPTs", sep = ""))
  
  # Read in parsimony trees:
  ParsimonyTrees <- ape::read.tree(gsub(".nex.trees", ".tre", i, fixed = TRUE))
  
  # Set working directory as NEXUS data:
  setwd(paste(TopDirectory, "/Data/Ages", sep = ""))
  
  # Read in parsimony trees:
  RawAges <- read.table(gsub(".nex.trees", ".txt", i, fixed = TRUE), header = TRUE, stringsAsFactors = FALSE)
  
  # Collapse raw ages to a set of tip ages using hard lower bound (max LAD) and smallest window (lowest FAD older than LAD):
  TipAges <- do.call(rbind, lapply(as.list(unique(RawAges[, "Taxon"])), function(x) {y <- unique(apply(RawAges[RawAges[, "Taxon"] == x, c("MaxMa", "MinMa"), drop = FALSE], 1, paste, collapse = "-")); y <- do.call(rbind, strsplit(y, "-")); LAD <- max(as.numeric(y[, 2])); FADs <- as.numeric(y[, 1]); FAD <- min(FADs[FADs >= LAD]); c(FAD, LAD)}))
  
  # Add rownames to tip ages:
  rownames(TipAges) <- unique(RawAges[, "Taxon"])
  
  # Add column names to tip ages:
  colnames(TipAges) <- c("FAD", "LAD")
  
  # If parsimony trees are a phylo object (one tree):
  if(class(ParsimonyTrees) == "phylo") {
    
    # Form into list:
    ParsimonyTrees <- list(ParsimonyTrees)
    
    # Make a multiphylo object (makes subsequent behaviour consistent):
    class(ParsimonyTrees) <- "multiPhylo"
    
  }

  # If parsimont tip names and age names do not match (happens because TNT often truncates names):
  if(length(setdiff(sort(rownames(TipAges)), sort(ParsimonyTrees[[1]]$tip.label))) > 0) {
    
    # Sort names independently (should form matches):
    NameMatches <- do.call(rbind, lapply(apply(cbind(sort(rownames(TipAges)), sort(ParsimonyTrees[[1]]$tip.label)), 1, as.list), function(x) unlist(x)))
    
    # Find just names that do not match:
    NameMatches <- NameMatches[NameMatches[, 1] != NameMatches [, 2], , drop = FALSE]
    
    # For each non-matching name:
    for(j in 1:nrow(NameMatches)) {
      
      # Update parsimony trees with longer tip age version:
      ParsimonyTrees <- lapply(ParsimonyTrees, function(x) {x$tip.label[x$tip.label == NameMatches [j, 2]] <- NameMatches [j, 1]; x})
      
      # Ensure list is still a multiphylo object:
      class(ParsimonyTrees) <- "multiPhylo"
      
    }
    
  }
  
  # Perform strat congruence metrics on input trees plus 1000 randomly generated trees:
  StratPhyloCongruenceResults <- StratPhyloCongruence(trees = c(BayesianTrees, ParsimonyTrees), ages = TipAges, rlen = 0, method = "basic", samp.perm = 1000, rand.perm = 1000, hard = TRUE, randomly.sample.ages = FALSE, fix.topology = FALSE, fix.outgroup = FALSE)
  
  # Pool sampled and randomly generated trees together:
  TreePool <- c(StratPhyloCongruenceResults$input.trees, StratPhyloCongruenceResults$rand.trees)
  
  # Pool minimum implied gaps for these trees together:
  MIGs <- c(StratPhyloCongruenceResults$input.tree.results[, "MIG"], StratPhyloCongruenceResults$rand.permutations[, "MIG"])
  
  # Get tree RF distances:
  RFMatrix <- MultiTreeDistance(trees = TreePool, distance = "RF", scale = TRUE)
  
  # Build PCoA from RF matrix:
  PCoA <- ape::pcoa(asin(sqrt(as.dist(RFMatrix))))
  
  # Plot as a circular Vornoi to file (probably a bad idea, but gives a basic plot - fancier contour plots would be better):
  #pdf(paste(TopDirectory, "/Data/Sample Plots/", gsub(".nex.trees", "", i, fixed = TRUE), ".pdf", sep = ""))
  #CircularVoronoi(x = PCoA$vectors[, 1], y = PCoA$vectors[, 2], heat = (MIGs - min(MIGs)) / max(MIGs - min(MIGs)), EdgeFactor = 1.1, NColours = 100)
  #dev.off()

}
