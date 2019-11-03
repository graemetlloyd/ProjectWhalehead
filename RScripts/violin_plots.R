# Violin Plots

files <- list.files(path = "Data/StratCongruenceResults/InputTreeResults/Yates_2003a.txt", pattern = "*.txt", full.names = TRUE, recursive = FALSE)

makeViolins <- function(x) {
  fi <- read.table(x, header = TRUE, sep = " ")
  p <- ggplot(fi, aes(factor(Bayes), MIG))
  q <- p + geom_violin(scale = "count") + geom_jitter(height = 0, width = 0.2, alpha = .4)
  new_filename <- gsub("StratCongruenceResults/InputTreeResults/", "ViolinPlots/", gsub(".txt", ".pdf", x))
  ggsave(new_filename, q)
}

lapply(files, makeViolins)
