# Violin Plots

files <- list.files(path="Data/StratCongruenceResults/InputTreeResults", pattern="*.txt", full.names=TRUE, recursive=FALSE)

lapply(files, makeViolins)

  
makeViolins <- function(x){
  fi <- read.table(x, header=TRUE, sep=",")
  p <- ggplot(fi, aes(factor(Bayes), MIG))
  q <- p + geom_violin(scale = "count")+ geom_jitter(height = 0, width = 0.2, alpha = .4)
  new_filename <- paste0(x, '.pdf')
  ggsave(new_filename, q)
}
  