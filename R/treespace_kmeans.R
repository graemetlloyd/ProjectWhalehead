Create space ordination including extra columns.
#' 
#' This function will take list of rwty.chains objects and produce a matrix of ordinations in treespace based on Robinson-Foulds Distance. The points can be colored and sized by extra columns listed in the `ext_columns` argument. This function is from Dan Warren and Rob Lanfear's `RWTY` package.
#'
#' @param chains A list of one or more rwty.chain objects
#' @param burnin The number of samples to remove from the start of the chain as burnin
#' @param n.points The number of points on each plot
#' @param fill.color The name of any column in your parameter file that you would like to use as a fill colour for the points of the plot.
#' @param ext_columns Any extra columns you'd like to color your plot by
#' @return A matrix of ordinations in treespace.
#'
#' @export custom_treespace

custom_treespace_kmeans <- function(name, chains, n.points = 1000, burnin = 0, fill.color = NA, ext_columns = NULL)
{
  chains = check.chains(chains)
  labels = names(chains)
  ptable = combine.ptables(chains, burnin = 0)
  if (!is.na(fill.color)) {
    if (fill.color %in% names(ptable)) {
    }
    else stop(sprintf("The fill.color name you supplied ('%s') wasn't found in your parameter table",fill.color))
  }
  chain = chains[[1]]
  step = as.integer((length(chain$trees) - burnin)/n.points)
  indices = seq(from = burnin + 1, to = length(chain$trees),
                by = step)
  indices = c(indices, range(length(chain$trees) -5, length(chain$trees)))
  
  trees = lapply(chains, function(x) x[["trees"]][indices])
  additional = unlist(lapply(length(chain$trees) * (0:(length(chains) -
                                                         1)), function(x) rep(x, length(indices))))
  ptable = ptable[(indices + additional), ]
  dimensions = 2
  alltrees = trees[[1]]
  if (length(trees) > 1) {
    for (i in 2:length(trees)) {
      alltrees = c(alltrees, trees[[i]])
    }
  }
  d <- tree.dist.matrix(alltrees)
  kclusts <- tibble(k = 1:9) %>%
    mutate(
      kclust = map(k, ~kmeans(d, .x)),
      tidied = map(kclust, tidy),
      glanced = map(kclust, glance),
      augmented = map(kclust, augment, d)
    )
  clusterings <- kclusts %>%
    unnest(glanced, .drop = TRUE)
  q <- ggplot(clusterings, aes(k, tot.withinss)) +
    geom_line()
  
  just_names <- strsplit(name, split = '/')[[1]][-1][2]
  n <- just_names[[1]][length(just_names[[1]])] 
  print(n)
  new_filename <- paste0("kMeans/",n, '.pdf')
  ggsave(new_filename, q)
}