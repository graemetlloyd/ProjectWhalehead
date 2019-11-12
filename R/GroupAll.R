#' Process MIG score for all the files used in this study.
#'
#' This function reads in the output of Strap's MIG claculations. It can read in one file, or many, provided a 
#' relative path to your files and the extension of the files. This function then prcoesses the files to retrieve 
#' the median MIG scores for both Bayesian and parsimony analyses.
#' @importFrom purrr map_dfr
#' @importFrom dplyr %>% 
#' @param path relative path to your data
#' @param extension Extension of the files to process
#' @return median_tbl A table of median values for MIG across datasets.
#' @export


summarize_medians <- function(user.path = NULL, extension = NULL){ 
  files <- list.files(path = user.path, pattern = '.txt', full.names = TRUE)
  print(files)
  map_data <- purrr::map_dfr(files, function(filename) {
    read.csv(filename, sep = ' ') %>% 
    mutate(fname = filename)
           }
        )
median_tbl <-  map_data %>% 
  group_by(Bayes) %>% 
  summarise_all(median)
return(median_tbl)
}