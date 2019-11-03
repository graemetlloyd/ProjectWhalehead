MIGSummary <- do.call(rbind, 
                      lapply(as.list(list.files()), function(x) 
                      {y <- read.table(x); c(BayesMeanMIG = mean(y[y[, "Bayes"] == 1, "MIG"]), ParsimonyMeanMIG = mean(y[y[, "Bayes"] == 0, "MIG"]), 
BayesMinMIG = min(y[y[, "Bayes"] == 1, "MIG"]), 
ParsimonyMinMIG = min(y[y[, "Bayes"] == 0, "MIG"]), 
BayesMaxMIG = max(y[y[, "Bayes"] == 1, "MIG"]), 
ParsimonyMaxMIG = max(y[y[, "Bayes"] == 0, "MIG"]))
                      }))