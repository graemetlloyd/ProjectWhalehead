#' Create treespace ordination including extra columns.
#' 
#' This function will take list of rwty.chains objects and produce a matrix of ordinations in treespace based on Robinson-Foulds Distance. The points can be colored and sized by extra columns listed in the `ext_columns` argument. This function is modified from the `RWTY` package.
#'
#' @param CombinedTrees A list of all the parsimony and Bayesian trees you'd like to plot
#' @param StratFile The logfile containing your stratigraphic congrunence values. Should contain a column to indicate if the tree is Bayesian or not.
#' @param n.points The number of points on each plot
#' @param fill.color The name of any column in your parameter file that you would like to use as a fill colour for the points of the plot.
#' @param ext_columns Any extra columns you'd like to color your plot by
#' @return A matrix of ordinations in treespace.
#'
#' @export custom_plots
#'  
custom_plots <- function (CombinedTrees, StratFile, ext_columns = c("MIG", "Bayes")) {
  rwty_obj <- rwty::load.trees(CombinedTrees, logfile = StratFile, skip=0) 
  points <- custom_treespace(rwty_obj, ext_columns = ext_columns)
  color_col <- ext_columns[1]
  size_column <- ext_columns[2]
  final<- ggplot(data = points, aes(x, y, colour=points$MIG, alpha = points$Bayes)) + geom_point(aes(size = points$Bayes)) +  scale_size(range = c(5,1))  + scale_alpha(range = c(1, .5))+ scale_colour_gradient(low = "blue", high = "yellow")
  just_file <- tools::file_path_sans_ext(CombinedTrees)
   new_filename <- gsub("CombinedTrees", "TreeSpace3D", paste0(just_file, '.pdf'))
  ggsave(new_filename, final)
  return(final)
} 