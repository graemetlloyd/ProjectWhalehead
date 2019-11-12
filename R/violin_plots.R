#' Create violin plots of our data
#' 
#' This function will take a file of stratigraphic congrunce metrics and make a violin plot. Returns nothing.
#' @import ggplot2
#' @param x file containing stratigraphic congruence metrics
#'
#' @export makeViolins


makeViolins <- function(x) {
  fi <- read.table(x, header = TRUE, sep = " ")
  p <- ggplot2::ggplot(fi, aes(factor(Bayes), MIG))
  q <- p + geom_violin(scale = "count") + geom_jitter(height = 0, width = 0.2, alpha = .4)
  just_names <- strsplit(x, split = '/')[[1]][-1][3]
  n <- just_names[[1]][length(just_names[[1]])] 
  print(n)
  new_filename <- paste0("ViolinPlots/",n, '.pdf')
  ggsave(new_filename, q)
}

