#' Create k-means dimension plots for data
#' 
#' This function will take a file of stratigraphic congrunce metrics, perform K-Means clustering
#' and plot the results
#' @param x file containing stratigraphic congruence metrics
#'
#' @export makeBoxplot


makeKmeansPlot <- function(x) {
  fi <- read.csv(x, header = TRUE, sep = " ")
  kclusts <- tibble(k = 1:9) %>%
  mutate(
        kclust = map(k, ~kmeans(data, .x)),
        tidied = map(kclust, tidy),
        glanced = map(kclust, glance),
        augmented = map(kclust, augment, data)
        )
  clusterings <- kclusts %>%
                 unnest(glanced, .drop = TRUE)
  q <- ggplot(clusterings, aes(k, tot.withinss)) +
          geom_line()
  
  just_names <- strsplit(x, split = '/')[[1]][-1][3]
  n <- just_names[[1]][length(just_names[[1]])] 
  print(n)
  new_filename <- paste0("kMeans/",n, '.pdf')
  ggsave(new_filename, q)
}