demar <- read.csv("Data/StratCongruenceResults/InputTreeResults/Yates_2003a.txt", sep = ' ')
demar_tibble <- as_tibble(demar)
demar_tibble %>%
  group_by(Bayes) %>% 
  summarize(mean_MIG = mean(MIG),
            max_MIG = max(MIG),
            min_MIG = min(MIG))