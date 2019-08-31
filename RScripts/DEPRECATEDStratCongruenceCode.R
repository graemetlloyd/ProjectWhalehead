# R script to perform strat congruence tests. All very slow!

# Load metatree and strap libraries:
library(hypRspace) # Install from github: devtools::install_github("graemetlloyd/hypRspace")
library(metatree) # Install from github: devtools::install_github("graemetlloyd/metatree")
library(strap)

# Local directory (modify for your local machine and hopefully rest should just work assuming the github directory structure matches):
TopDirectory <- "~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead"

# Set as Rev Bayes directory:
setwd(paste(TopDirectory, "/Data/BayesianTopologies", sep = ""))

# Find all Rev bayes files there:
BayesianFiles <- list.files()

# For each Rev Bayes file:
for(i in BayesianFiles) {
  
  # Print loop position to screen:
  cat(i, " ")
  
  # Set workibg directory to Bayesian topology directory:
  setwd(paste(TopDirectory, "/Data/BayesianTopologies", sep = ""))
  
  # Read in Bayesian trees:
  BayesianTrees <- ape::read.tree(i)
  
  # If Bayesian trees are a phylo object (one tree):
  if(class(BayesianTrees) == "phylo") {
    
    # Form into list:
    BayesianTrees <- list(BayesianTrees)
    
    # Make a multiphylo object (makes subsequent behaviour consistent):
    class(BayesianTrees) <- "multiPhylo"
    
  }
  
  # Set working directory as NEXUS data:
  setwd(paste(TopDirectory, "/Data/MPTs", sep = ""))
  
  # Read in parsimony trees:
  ParsimonyTrees <- ape::read.tree(i)
  
  # Set working directory as NEXUS data:
  setwd(paste(TopDirectory, "/Data/Ages", sep = ""))
  
  # Read in parsimony trees:
  RawAges <- read.table(gsub(".tre", ".txt", i, fixed = TRUE), header = TRUE, stringsAsFactors = FALSE)
  
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
  StratPhyloCongruenceResults <- StratPhyloCongruence(trees = c(BayesianTrees, ParsimonyTrees), ages = TipAges, rlen = 0, method = "basic", samp.perm = 1000, rand.perm = 1000, hard = TRUE, randomly.sample.ages = TRUE, fix.topology = TRUE, fix.outgroup = TRUE)
  
  # Build lists of each congruence metric for boxplots:
  SCIList <- list(Bayesian = unname(StratPhyloCongruenceResults$input.tree.results[1:length(BayesianTrees), "SCI"]), Parsimony = unname(StratPhyloCongruenceResults$input.tree.results[(length(BayesianTrees) + 1):(length(BayesianTrees) + length(ParsimonyTrees)), "SCI"]))
  RCIList <- list(Bayesian = unname(StratPhyloCongruenceResults$input.tree.results[1:length(BayesianTrees), "RCI"]), Parsimony = unname(StratPhyloCongruenceResults$input.tree.results[(length(BayesianTrees) + 1):(length(BayesianTrees) + length(ParsimonyTrees)), "RCI"]))
  GERList <- list(Bayesian = unname(StratPhyloCongruenceResults$input.tree.results[1:length(BayesianTrees), "GER"]), Parsimony = unname(StratPhyloCongruenceResults$input.tree.results[(length(BayesianTrees) + 1):(length(BayesianTrees) + length(ParsimonyTrees)), "GER"]))
  MSMStarList <- list(Bayesian = unname(StratPhyloCongruenceResults$input.tree.results[1:length(BayesianTrees), "MSM*"]), Parsimony = unname(StratPhyloCongruenceResults$input.tree.results[(length(BayesianTrees) + 1):(length(BayesianTrees) + length(ParsimonyTrees)), "MSM*"]))
  GERStarList <- list(Bayesian = unname(StratPhyloCongruenceResults$input.tree.results[1:length(BayesianTrees), "GER*"]), Parsimony = unname(StratPhyloCongruenceResults$input.tree.results[(length(BayesianTrees) + 1):(length(BayesianTrees) + length(ParsimonyTrees)), "GER*"]))
  GERtList <- list(Bayesian = unname(StratPhyloCongruenceResults$input.tree.results[1:length(BayesianTrees), "GERt"]), Parsimony = unname(StratPhyloCongruenceResults$input.tree.results[(length(BayesianTrees) + 1):(length(BayesianTrees) + length(ParsimonyTrees)), "GERt"]))
  MIGList <- list(Bayesian = unname(StratPhyloCongruenceResults$input.tree.results[1:length(BayesianTrees), "MIG"]), Parsimony = unname(StratPhyloCongruenceResults$input.tree.results[(length(BayesianTrees) + 1):(length(BayesianTrees) + length(ParsimonyTrees)), "MIG"]))
  
  # Make boxplots comparing Bayesian and parsimony trees:
  pdf(paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/Boxplots/SCI/", gsub(".tre", "", i), ".pdf", sep = ""))
  boxplot(SCIList, ylab = "Stratigraphic Consistency Index")
  dev.off()
  pdf(paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/Boxplots/RCI/", gsub(".tre", "", i), ".pdf", sep = ""))
  boxplot(RCIList, ylab = "Relative Completeness Index")
  dev.off()
  pdf(paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/Boxplots/MSMStar/", gsub(".tre", "", i), ".pdf", sep = ""))
  boxplot(MSMStarList, ylab = "Manhattan Stratigraphic Measure*")
  dev.off()
  pdf(paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/Boxplots/GER/", gsub(".tre", "", i), ".pdf", sep = ""))
  boxplot(GERList, ylab = "Gap Excess Ratio (GER)")
  dev.off()
  pdf(paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/Boxplots/GERStar/", gsub(".tre", "", i), ".pdf", sep = ""))
  boxplot(GERStarList, ylab = "Gap Excess Ratio (GER*)")
  dev.off()
  pdf(paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/Boxplots/GERt/", gsub(".tre", "", i), ".pdf", sep = ""))
  boxplot(GERtList, ylab = "Gap Excess Ratio (GERt)")
  dev.off()
  pdf(paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/Boxplots/MIG/", gsub(".tre", "", i), ".pdf", sep = ""))
  boxplot(MIGList, ylab = "Minimum Implied Gap")
  dev.off()
  
  # Write out results tables:
  write.table(cbind(StratPhyloCongruenceResults$input.tree.results, Bayes = c(rep(1, length(BayesianTrees)), rep(0, length(ParsimonyTrees)))), paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/StratCongruenceResults/InputTreeResults/", gsub(".tre", "", i), ".txt", sep = ""))
  write.table(StratPhyloCongruenceResults$rand.permutations, paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/StratCongruenceResults/RandomTreeResults/", gsub(".tre", "", i), ".txt", sep = ""))
  
  # Write out trees:
  ape::write.tree(StratPhyloCongruenceResults$input.trees[1:length(BayesianTrees)], paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/StratCongruenceResults/TimeScaledTrees/BayesianTrees/", gsub(".tre", "", i), ".tre", sep = ""))
  ape::write.tree(StratPhyloCongruenceResults$input.trees[(length(BayesianTrees) + 1):(length(BayesianTrees) + length(ParsimonyTrees))], paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/StratCongruenceResults/TimeScaledTrees/ParsimonyTrees/", gsub(".tre", "", i), ".tre", sep = ""))
  ape::write.tree(StratPhyloCongruenceResults$rand.trees, paste("~/Documents/Publications/in prep/Strat congruence - April/ProjectWhalehead/Data/StratCongruenceResults/TimeScaledTrees/RandomTrees/", gsub(".tre", "", i), ".tre", sep = ""))

}
