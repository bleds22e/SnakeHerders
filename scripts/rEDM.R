# Attempt rEDM #
# EKB; June 2022 #

library(tidyverse)

# load prairie dog and snake ethograms
seq <- read_csv("data/EKP_matrix_snake.csv")

# PDOG ------------------------------------------------------------####

# Clean up ethogram tags and states #

# get required columns, replace blanks with NA, and fill in videoID values
sequence <- seq %>% dplyr::select(videoID = `Video ID`, # will need later
                                  ethogram_tag_pdog = `Ethogram Tag...10`) %>% 
  mutate(across(where(is.character), ~ na_if(.,""))) %>% 
  fill(videoID, .direction = "down")

# standardize ethogram_tags
sequence <- sequence %>% 
  mutate(ethogram_tag_pdog = replace(ethogram_tag_pdog, ethogram_tag_pdog == "Rocking +Sniffing", "Rocking + Sniffing"),
         ethogram_tag_pdog = replace(ethogram_tag_pdog, ethogram_tag_pdog == "Jump Yip", "Jump yip"),
         ethogram_tag_pdog = replace(ethogram_tag_pdog, ethogram_tag_pdog == "Moving away", "Retreat"),
         ethogram_tag_pdog = replace(ethogram_tag_pdog, ethogram_tag_pdog == "Moving toward", "Approach"))

# create df with numbers for each ethogram_tag
states_pdog <- tibble(number_pdog = seq(1:22),
                 ethogram_tag_pdog = c("Posting", "Jump yip", "Approach", "Scratching",
                                  "Pushing dirt", "Rocking", "Sniffing", 
                                  "Rocking + Sniffing", "Approach + Sniffing",
                                  "Retreat", "Jump back", "Escorting", "Foraging",
                                  "Investigating", "Touching face", "Strike",
                                  "Vocalizing", "Head to Head", "Fake foraging",
                                  "Observing", "Greet-kiss", "Running"))

# merge df and remove NAs from ethogram_tag or number columns
sequence <- full_join(sequence, states_pdog)

# SNAKE --------------------------------------------------------------------####

# Clean up ethogram tags and states #

# get required columns, replace blanks with NA, and fill in videoID values
snake_seq <- seq %>% dplyr::select(videoID = `Video ID`, # will need later
                                         number_snake = `Ethogram Tag...14`) %>% 
  mutate(across(where(is.character), ~ na_if(.,""))) %>% 
  fill(videoID, .direction = "down")

states_snake <- seq %>% 
  dplyr::select(states = `Numerical Category...13`) %>% 
  drop_na() %>% 
  separate(states, c("number_snake", "ethogram_tag_snake"), sep = " = ") %>% 
  mutate(number_snake = as.numeric(number_snake))

# merge df and remove NAs from ethogram_tag or number columns
snake_seq <- full_join(snake_seq, states_snake)


# COMBINE ------------------------------------------------------------------####

all_dat <- bind_cols(sequence, dplyr::select(snake_seq, -videoID)) %>% 
  slice(1:(n()-20)) %>%         # remove last 20 rows for juvnile pdogs
  group_split(videoID) %>%      # create multiple dfs based on videoID column
  map_dfr(~ .x %>% add_row())

