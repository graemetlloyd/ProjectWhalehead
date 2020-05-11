# Set working directory to combined trees folder:
setwd("~/Documents/Publications/in review/Strat congruence - April/ProjectWhalehead/Data/CombinedTrees")

# Get list of all data sets:
DataSets <- gsub(".nex", "", list.files())

# Set up empty vectors to store N axes covering 95% of the varance and the tip count for each data set:
Var95Axes <- NTips <- vector(mode = "numeric")

# For each data set:
for(i in DataSets) {
  
  # Set data set name:
  DataSet <- i
  
  # Read in trees from GitHub repo:
  Trees <- ape::read.nexus(paste("https://raw.githubusercontent.com/graemetlloyd/ProjectWhalehead/master/Data/CombinedTrees/", DataSet, ".nex", sep = ""))
  
  # Count and store number of tips:
  NTips <- c(NTips, ape::Ntip(Trees[[1]]))
  
  # Downsample to just 1000 trees (so we don't have to wait until the inevitable heat death of the universe):
  Trees <- Trees[sample(x = 1:length(Trees), size = 1000)]
  
  # Note outgroup taxon:
  OutgroupTaxon <- rownames(Claddis::ReadMorphNexus(paste("https://raw.githubusercontent.com/graemetlloyd/ProjectWhalehead/master/Data/NEXUS/", DataSet, ".nex", sep =""))$Matrix_1$Matrix)[1]
  
  # Read in age data:
  RawAges <- read.table(paste("https://raw.githubusercontent.com/graemetlloyd/ProjectWhalehead/master/Data/Ages/", DataSet, ".txt", sep = ""), header = TRUE, stringsAsFactors = FALSE)
  
  # Build tip ages matrix:
  TipAges <- do.call(rbind, lapply(as.list(unique(RawAges[, "Taxon"])), function(x) {y <- unique(apply(RawAges[RawAges[, "Taxon"] == x, c("MaxMa", "MinMa"), drop = FALSE], 1, paste, collapse = "-")); y <- do.call(rbind, strsplit(y, "-")); LAD <- max(as.numeric(y[, 2])); FADs <- as.numeric(y[, 1]); FAD <- min(FADs[FADs >= LAD]); c(FAD, LAD)}))
  
  # Add rownames (taxa) to tip ages:
  rownames(TipAges) <- unique(RawAges[, "Taxon"])
  
  # Add column names to tip ages:
  colnames(TipAges) <- c("FAD", "LAD")
  
  # Perform stratigraphic congruence metrics and store them:
  StratPhyloCongruenceResults <- strap::StratPhyloCongruence(trees = Trees, ages = TipAges, outgroup.taxon = OutgroupTaxon, rlen = 0, method = "basic", samp.perm = 1000, rand.perm = 1000, hard = TRUE, randomly.sample.ages = FALSE, fix.topology = TRUE, fix.outgroup = TRUE)
  
  # Get Robinson-Foulds distances:
  RFD <- phangorn::RF.dist(Trees, normalize = TRUE)
  
  # Produce 2-D MDS:
  MDSk2 <- cmdscale(RFD, k = 2)
  
  # Produce maximum k MDS:
  MDSmaxk <- cmdscale(RFD, k = length(Trees) - 1)
  
  # Get 2-D MDS pairwse distances:
  MDSDk2 <- dist(MDSk2)
  
  # Get k-D MDS pairwse distances:
  MDSDmaxk <- dist(MDSmaxk)
  
  # Find out how many MDS axes are required to cover 95% of the variance and store:
  Var95Axes <- c(Var95Axes, min(which(cumsum(apply(MDSmaxk, 2, var) / sum(apply(MDSmaxk, 2, var))) >= 0.95)))
  
  # Get MIG differences:
  MIGDistances <- dist(unname(StratPhyloCongruenceResults$input.tree.results[, "MIG"]))
  
  # Start plotting to PDF:
  pdf(paste("~/Documents/Publications/in review/Strat congruence - April/ProjectWhalehead/Data/MIGDifferences/", i, ".pdf", sep = ""), width = 10, height = 5)
  
  # Create one-by-two plot set up:
  par(mfrow = c(1, 2))
  
  # Plot MDS distance versus MIG difference for first two axes (2-D MDS):
  plot(MDSDk2, MIGDistances, pch = 20, cex = 0.5, col = "red", xlab = "MDS distance (2-dimensional)", ylab = "MIG difference")
  
  # Plot MDS distance versus MIG difference for all axes (k-D MDS):
  plot(MDSDmaxk, MIGDistances, pch = 20, cex = 0.5, col = "red", xlab = "MDS distance (multi-dimensional)", ylab = "MIG difference")
  
  # Stop plotting to PDF:
  dev.off()
  
}

# Add data set names to N axes and N tips vectors:
names(Var95Axes) <- names(NTips) <- DataSets

# Start plotting to PDF:
pdf("~/Documents/Publications/in review/Strat congruence - April/ProjectWhalehead/Figures/FigS1.pdf", width = 10, height = 5)

# Set plotting layout:
graphics::layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))

# Show histogram of N axes required for 95% variance:
hist(Var95Axes, breaks = seq(from = 0, to = 400, length.out = 33), col = "grey", border = 0, xlab = "N MDS axes required to cover 95% of the variance", main = "")

# Show relationship between N tips and N axes:
plot(NTips, Var95Axes, log = "xy", pch = 20, xlab = "Number of tips", ylab = "N MDS axes required to capture 95% of the variance"); abline(a = summary(lm(Var95Axes ~ NTips))$coefficients["(Intercept)", "Estimate"], b = summary(lm(Var95Axes ~ NTips))$coefficients["NTips", "Estimate"], untf = TRUE)

# Stop plotting to PDF:
dev.off()
