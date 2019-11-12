#' Generate the S1 table of Wright and Lloyd 2019
#' 
#' This function takes as its arguments what straitgraphic congruence measure you want
#' from a file of dataset of stratigraphic congruence metrics. Assumes by default that you want to group by the column that says whether the analysis was Bayesian or not.
#' 
#' @param files a list of datafiles you'd like the process
#' @param measure Which stratigraphic congruence measure you'd like to look at
#' @param by_col If you'd like the results grouped. Defaults to a column called "Bayes"
#' @return summary a data frame showing mean scores for your metric of choice, grouped by the by_col col
#' @export generate_fig_s1


generate_fig_s1 <- function(files, measure, by_col = "Bayes"){
summary <- do.call(rbind, 
                      lapply(as.list(list.files()), function(x) 
                      {y <- read.table(x); c(BayesMeanMIG = mean(y[y[, by_col] == 1, measure]), ParsimonyMeanMIG = mean(y[y[, by_col] == 0, measure]), 
BayesMinMIG = min(y[y[, by_col] == 1, measure]), 
ParsimonyMinMIG = min(y[y[, by_col] == 0, measure]), 
BayesMaxMIG = max(y[y[, by_col] == 1, measure]), 
ParsimonyMaxMIG = max(y[y[, by_col] == 0, measure]))
                      }))
return(summary)
}