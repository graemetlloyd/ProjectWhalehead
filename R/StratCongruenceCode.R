#' R script to perform strat congruence tests.

#' This R function runs statigraphic congruence tests on both a Bayesian posterior sample and a set of most parsimonious trees. Nothing is returned, but boxplots of these figures are saved to a results directory.

#' @import hypRspace 
#' @import metatree 
#' @import strap 
#' @export make_strat_congruence

make_strat_congruence <- function(RevBayesFiles, NexusFiles, AgeFiles, ParsimonyFiles, RB_output, Parsimony_output) {
  
  RevBayesFiles <- "../Data/RevBayesOutput/"
  NexusFiles <- "../Data/NEXUS/"
  AgeFiles <- "../Data/Ages/"
  ParsimonyFiles <- "../Data/MPTS/"
  RB_output <- "revbayes_output"
  Parsimony_output <- "parsimony_out"
  
  StartDirectory <- getwd()
  setwd("../")
  TopDirectory <- getwd()
  setwd(StartDirectory)
  
  # Find all Rev bayes files there:
  RBTrees <- list.files(RevBayesFiles, full.names = TRUE)

  # For each Rev Bayes file:
  for(i in  RBTrees) {
  
  # (Re)set to Rev Bayes directory:

  # Read in Rev bayes output:
  RevBayesOutput <- read.table(i, header = TRUE, stringsAsFactors = FALSE)
  
  # Set working directory as NEXUS data:
  setwd(paste(TopDirectory, gsub("../", "/", NexusFiles, fixed = TRUE), sep = ""))
  
  # Read in outgroup taxon:
  OutgroupTaxon <- rownames(Claddis::ReadMorphNexus(gsub(".trees", "", i, fixed = TRUE))$Matrix_1$Matrix)[1]
  
  # Get just the trees in ape format:
  BayesianTrees <- ape::read.tree(text = RevBayesOutput[, "tree"])
  
  # Remove branch lengths (allows collapsing to unique topologies below):
  BayesianTrees <- lapply(BayesianTrees, function(x) {x$edge.length <- NULL; x})
  
  # Reformat as a multiPhylo object:
  class(BayesianTrees) <- "multiPhylo"
  
  # Reroot trees on outgroup taxon:
  BayesianTrees <- root(BayesianTrees, outgroup = OutgroupTaxon)
  
  # Now reroot properly because ape is useless:
  BayesianTrees <- ape::read.tree(text = unlist(lapply(BayesianTrees, function(x) paste("(", OutgroupTaxon, ",", gsub(";", ");", ape::write.tree(drop.tip(x, tip = OutgroupTaxon))), sep = ""))))
  
  # Ste working directory to ages folder:
  setwd(AgeFiles)
  
  # Read in parsimony trees:
  RawAges <- read.table(gsub(".nex.trees", ".txt", i, fixed = TRUE), header = TRUE, stringsAsFactors = FALSE)
  
  # Collapse raw ages to a set of tip ages using hard lower bound (max LAD) and smallest window (lowest FAD older than LAD):
  TipAges <- do.call(rbind, lapply(as.list(unique(RawAges[, "Taxon"])), function(x) {y <- unique(apply(RawAges[RawAges[, "Taxon"] == x, c("MaxMa", "MinMa"), drop = FALSE], 1, paste, collapse = "-")); y <- do.call(rbind, strsplit(y, "-")); LAD <- max(as.numeric(y[, 2])); FADs <- as.numeric(y[, 1]); FAD <- min(FADs[FADs >= LAD]); c(FAD, LAD)}))
  
  # Add rownames to tip ages:
  rownames(TipAges) <- unique(RawAges[, "Taxon"])
  
  # Add column names to tip ages:
  colnames(TipAges) <- c("FAD", "LAD")
  
  # Set working directory as NEXUS data:
  setwd(paste(TopDirectory, ParsimonyFiles, sep = ""))
  
  # Read in parsimony trees:
  ParsimonyTrees <- ape::read.tree(gsub(".nex.trees", ".tre", i))
  
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
  StratPhyloCongruenceResults <- StratPhyloCongruence(trees = c(BayesianTrees, ParsimonyTrees), ages = TipAges, outgroup.taxon = OutgroupTaxon, rlen = 0, method = "basic", samp.perm = 1000, rand.perm = 1000, hard = TRUE, randomly.sample.ages = FALSE, fix.topology = TRUE, fix.outgroup = TRUE)

  # Build lists of each congruence metric for boxplots:
  SCIList <- list(Bayesian = unname(StratPhyloCongruenceResults$input.tree.results[1:length(BayesianTrees), "SCI"]), Parsimony = unname(StratPhyloCongruenceResults$input.tree.results[(length(BayesianTrees) + 1):(length(BayesianTrees) + length(ParsimonyTrees)), "SCI"]))
  RCIList <- list(Bayesian = unname(StratPhyloCongruenceResults$input.tree.results[1:length(BayesianTrees), "RCI"]), Parsimony = unname(StratPhyloCongruenceResults$input.tree.results[(length(BayesianTrees) + 1):(length(BayesianTrees) + length(ParsimonyTrees)), "RCI"]))
  GERList <- list(Bayesian = unname(StratPhyloCongruenceResults$input.tree.results[1:length(BayesianTrees), "GER"]), Parsimony = unname(StratPhyloCongruenceResults$input.tree.results[(length(BayesianTrees) + 1):(length(BayesianTrees) + length(ParsimonyTrees)), "GER"]))
  MSMStarList <- list(Bayesian = unname(StratPhyloCongruenceResults$input.tree.results[1:length(BayesianTrees), "MSM*"]), Parsimony = unname(StratPhyloCongruenceResults$input.tree.results[(length(BayesianTrees) + 1):(length(BayesianTrees) + length(ParsimonyTrees)), "MSM*"]))
  GERStarList <- list(Bayesian = unname(StratPhyloCongruenceResults$input.tree.results[1:length(BayesianTrees), "GER*"]), Parsimony = unname(StratPhyloCongruenceResults$input.tree.results[(length(BayesianTrees) + 1):(length(BayesianTrees) + length(ParsimonyTrees)), "GER*"]))
  GERtList <- list(Bayesian = unname(StratPhyloCongruenceResults$input.tree.results[1:length(BayesianTrees), "GERt"]), Parsimony = unname(StratPhyloCongruenceResults$input.tree.results[(length(BayesianTrees) + 1):(length(BayesianTrees) + length(ParsimonyTrees)), "GERt"]))
  MIGList <- list(Bayesian = unname(StratPhyloCongruenceResults$input.tree.results[1:length(BayesianTrees), "MIG"]), Parsimony = unname(StratPhyloCongruenceResults$input.tree.results[(length(BayesianTrees) + 1):(length(BayesianTrees) + length(ParsimonyTrees)), "MIG"]))
  
  # Make boxplots comparing Bayesian and parsimony trees:
  pdf(paste("Data/Boxplots/SCI/", gsub(".nex.trees", "", i), ".pdf", sep = ""))
  boxplot(SCIList, ylab = "Stratigraphic Consistency Index")
  dev.off()
  pdf(paste("Data/Boxplots/RCI/", gsub(".nex.trees", "", i), ".pdf", sep = ""))
  boxplot(RCIList, ylab = "Relative Completeness Index")
  dev.off()
  pdf(paste("Data/Boxplots/MSMStar/", gsub(".nex.trees", "", i), ".pdf", sep = ""))
  boxplot(MSMStarList, ylab = "Manhattan Stratigraphic Measure*")
  dev.off()
  pdf(paste("Data/Boxplots/GER/", gsub(".nex.trees", "", i), ".pdf", sep = ""))
  boxplot(GERList, ylab = "Gap Excess Ratio (GER)")
  dev.off()
  pdf(paste("Boxplots/GERStar/", gsub(".nex.trees", "", i), ".pdf", sep = ""))
  boxplot(GERStarList, ylab = "Gap Excess Ratio (GER*)")
  dev.off()
  pdf(paste("Data/Boxplots/GERt/", gsub(".nex.trees", "", i), ".pdf", sep = ""))
  boxplot(GERtList, ylab = "Gap Excess Ratio (GERt)")
  dev.off()
  pdf(paste("Data/Boxplots/MIG/", gsub(".nex.trees", "", i), ".pdf", sep = ""))
  boxplot(MIGList, ylab = "Minimum Implied Gap")
  dev.off()
  
  # Write out results tables:
  write.table(cbind(StratPhyloCongruenceResults$input.tree.results, Bayes = c(rep(1, length(BayesianTrees)), rep(0, length(ParsimonyTrees)))), paste("Data/StratCongruenceResults/InputTreeResults/", gsub(".nex.trees", "", i), ".txt", sep = ""))
  write.table(StratPhyloCongruenceResults$rand.permutations, paste("Data/StratCongruenceResults/RandomTreeResults/", gsub(".nex.trees", "", i), ".txt", sep = ""))
  
  # Write out trees:
  ape::write.tree(StratPhyloCongruenceResults$input.trees[1:length(BayesianTrees)], paste("Data/StratCongruenceResults/TimeScaledTrees/BayesianTrees/", gsub(".nex.trees", "", i), ".tre", sep = ""))
  ape::write.tree(StratPhyloCongruenceResults$input.trees[(length(BayesianTrees) + 1):(length(BayesianTrees) + length(ParsimonyTrees))], paste("StratCongruenceResults/TimeScaledTrees/ParsimonyTrees/", gsub(".nex.trees", "", i), ".tre", sep = ""))
  ape::write.tree(StratPhyloCongruenceResults$rand.trees, paste("~Data/StratCongruenceResults/TimeScaledTrees/RandomTrees/", gsub(".nex.trees", "", i), ".tre", sep = ""))

  }
return(0)
}
