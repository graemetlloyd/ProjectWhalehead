custom_treespace <- function (chains, n.points = 100, burnin = 0, fill.color = NA) 
{
  chains = check.chains(chains)
  labels = names(chains)
  ptable = combine.ptables(chains, burnin = 0)
  if (!is.na(fill.color)) {
    if (fill.color %in% names(ptable)) {
    }
    else stop(sprintf("The fill.color name you supplied ('%s') wasn't found in your parameter table", 
                      fill.color))
  }
  chain = chains[[1]]
  step = as.integer((length(chain$trees) - burnin)/n.points)
  indices = seq(from = burnin + 1, to = length(chain$trees), 
                by = step)
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
  if (sum(d) == 0) {
    x <- rep(0, length(alltrees))
    y <- rep(0, length(alltrees))
    mds <- data.frame(x = x, y = y)
  }
  else {
    mds <- cmdscale(d, k = dimensions)
  }
  points <- as.data.frame(mds)
  row.names(points) <- seq(nrow(points))
  names(points) <- c("x", "y")
  points$chain = ptable$chain
  points$sample = ptable$sample
  points$generation = ptable$generation
  points$bayes = ptable$X1
  if (!is.na(fill.color)) 
    points[, fill.color] = as.numeric(ptable[[fill.color]])
  return(points)
}