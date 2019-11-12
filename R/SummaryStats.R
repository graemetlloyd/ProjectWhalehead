#' Generate the S1 table of Wright and Lloyd 2019
#' 
#' This function takes as its arguments what straitgraphic congruence measure you want
#' from a file of dataset of stratigraphic congruence metrics. Assumes by default that you want to group by the column that says whether the analysis was Bayesian or not.
#' 
#' @param dataf The file you'd like to examine
#' @param measure Which stratigraphic congruence measure you'd like to look at
#' @param by_col If you'd like the results grouped. Defaults to a column called "Bayes"
#' @return summary a data frame showing mean scores for your metric of choice, grouped by the by_col col
#' @export summary_stat


summary_stat <- function(dataf){
  table <- read.csv(dataf, sep = ' ')
  tib <- tibble::as_tibble(table)
  print(tib)
  outt <- tib %>%
    dplyr::group_by(Bayes) %>% 
    summarize(mean_val =mean(MIG),
              max_val = max(MIG),
              min_val = min(MIG))
  return(outt)  
}