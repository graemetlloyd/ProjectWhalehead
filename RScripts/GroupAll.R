 
files <- list.files(pattern = "*.txt")

map_data <- map_dfr(files, function(filename) {
  read.csv(filename, sep = ' ') %>% 
  mutate(fname = filename)
           }
                    )
median_tbl <-  map_data %>% 
  group_by(Bayes, date) %>% 
  summarise_all(median)