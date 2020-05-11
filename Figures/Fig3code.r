# Install required packages (run once):
#install.packages(c("ape", "devtools", "phangorn", "plot3D"), dependencies = TRUE)
#devtools::install_github("graemetlloyd/hypRspace")

# Set number of tips:
Ntips <- 6

# Generate all birfurcating trees with N tips and label A, B, C...[N]:
AllLabelledFourTipBifurcatingTrees <- ape::read.tree(text = unlist(hypRspace::AllLabelledNewicksFromNewickNumberStrings(hypRspace::AllNewickNumbersByPartition(Ntips)[[(Ntips - 1)]], LETTERS[1:Ntips])))

# Generate fake ages for each tip:
Ages1 <- matrix(c(1:Ntips, 1:Ntips), ncol = 2, dimnames = list(LETTERS[1:Ntips], c("FAD", "LAD")))

# Calculate minimum Implied Gaps for each tree:
MIGs1 <- unlist(lapply(AllLabelledFourTipBifurcatingTrees, function(x) sum(strap::DatePhylo(tree = x, ages = Ages1)$edge.length)))

# Get pairwise rescaled Robinson-Fould distances between each tree:
RFD <- phangorn::RF.dist(AllLabelledFourTipBifurcatingTrees, normalize = TRUE)

# Generate 2D MDS from RFDs:
MDS2k <- cmdscale(RFD, k = 2)

# Set Cartesian coordinates from 2D MDS and MIG:
x <- MDS2k[, 1]
y <- MDS2k[, 2]
z1 <- MIGs1

# Function to sample space as a grid:
SampleAsGrid <- function(x, y, z, density, xlim, ylim, weight = 1) {
  
  # Set x values as density-squared many points between xlimits:
  x_out <- matrix(rep(x = seq(from = min(xlim), to = max(xlim), length.out = density), times = density), ncol = density, byrow = TRUE)
  
  # Set y values as density-squared many points between ylimits:
  y_out <- matrix(rep(x = seq(from = min(ylim), to = max(ylim), length.out = density), times = density), ncol = density, byrow = FALSE)
  
  # Crete empty z values matrix:
  z_out <- matrix(NA, nrow = density, ncol = density)
  
  # For each x-value:
  for(i in 1:density) {
    
    # For each y-value:
    for(j in 1:density) {
      
      # Get current x-value:
      x_value <- x_out[i, j]
      
      # Get current y-value:
      y_value <- y_out[i, j]
      
      # Find absolute differences between current x-value and each input x value:
      x_diffs <- abs(x - x_value)
      
      # Find absolute differences between current y-value and each input y value:
      y_diffs <- abs(y - y_value)
      
      # Set initial z weights as distances between current point and all input points:
      z_weights <- sqrt((x_diffs ^ 2) + (y_diffs ^ 2))
      
      # Invert weights (want nearest points to be highest weight):
      z_weights <- 1 / (z_weights ^ weight)
      
      # If any weights are infinity (sampling exactly a data point):
      if(any(z_weights == Inf)) {
        
        # Set z_out as exactly the sampled data point (or mean of multiple data points):
        z_out[i, j] <- mean(z[z_weights == Inf])
        
      # If no weights are zero (sampling exactly a data point)
      } else {
        
        # Set z_out as weighted mean:
        z_out[i, j] <- sum(z * z_weights) / sum(z_weights)
        
      }
      
    }
    
  }
  
  # Return x y and z matrices (grid of values):
  return(list(x = x_out, y = y_out, z = z_out))
  
}

# Sample MIG landscape as a grid (interpolating interstitial values):
Output1 <- SampleAsGrid(x, y, z1, density = 250, xlim = c(min(c(x, y)), max(c(x, y))), ylim = c(min(c(x, y)), max(c(x, y))), weight = 3)

# Start plotting to file:
pdf("~/Documents/Publications/in review/Strat congruence - April/ProjectWhalehead/Figures/Fig3.pdf", width = 10, height = 5, useDingbats = FALSE)

# Set up two-panel plot:
par(mfrow = c(1, 2), mai = c(0.3, 0.3, 0.3, 0.3))

# Plot 3D MIG landscape:
plot3D::surf3D(Output1$x, Output1$y, Output1$z * -1, colvar = Output1$z * -1, colkey = TRUE, zlim = c(floor(min(Output1$z * -1)), ceiling(max(Output1$z * -1))))

# Plot top-down plot of space:
plot3D::surf3D(Output1$x, Output1$y, Output1$z * -1, colvar = Output1$z * -1, colkey = FALSE, zlim = c(floor(min(Output1$z * -1)), ceiling(max(Output1$z * -1))), theta = 0, phi = 90, d = 1000)

# Add x-axis:
lines(x = c(-0.5, 0.5), y = c(0, 0))

# Add y-axis:
lines(x = c(0, 0), y = c(-0.5, 0.5))

# Add grid lines:
for(i in -5:5) {lines(x = c(-0.5, 0.5), y = c(i/10, i/10), lty = 3); lines(x = c(i/10, i/10), y = c(-0.5, 0.5), lty = 3)}

# Highlight three points:
points(x = MDS2k[c(223, 266, 942), 1], y = MDS2k[c(223, 266, 942), 2], pch = 21, cex = 2.5, col = "black", bg = "white")
text(x = MDS2k[c(223, 266, 942), 1], y = MDS2k[c(223, 266, 942), 2], labels = c("A", "B", "C"))

# Stop plotting to file:
dev.off()
