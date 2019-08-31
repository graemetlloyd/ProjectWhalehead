library(tools)

files <- list.files(path = "Data/StratCongruenceResults/InputTreeResults", pattern = "*.txt", full.names = TRUE, recursive = FALSE)
tre_files <- list.files(path = "Data/CombinedTrees", pattern = "*.nex", full.names = TRUE, recursive = FALSE)

custom_plots <- function (x, y) {
  #trees <- ape::read.nexus(y)
  #new_filename <- paste0(y, ".nex")
  #ape::write.nexus(trees, file = new_filename)
  rwty_obj <- rwty::load.trees(file = y, logfile = x)
  
  # THE ABOVE DOESN'T SEEM TO WORK CORRECTLY AS THE HEADER ROW IS BEING IGRNORED AND BEING SET AS THE FIRST ROW INSTEAD
  
  points = custom_treespace(rwty_obj, length(trees))
  pp <- ggplot(points %>% arrange(bayes), aes(x = x, y = y, color = rwty_obj$ptable$X32))
  final <- pp + geom_point(aes(size=bayes)) + scale_size(range = c(5, 1)) + guides(size = guide_legend(reverse=TRUE)) + scale_colour_gradient(low = "blue", high = "yellow")
  just_file <- file_path_sans_ext(y)
  new_filename <- paste0(just_file, '.pdf')
  ggsave(new_filename, final)
}
mapply(custom_plots, files, tre_files)
