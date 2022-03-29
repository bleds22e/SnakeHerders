# Prairie Dog Behavior MCMC #
# EKB; March 16, 2022 #


# PACKAGES and DATA #

# packages
library(tidyverse)
library(markovchain)

# load matrix
matrix <- read_csv("data/JLV_matrix.csv", na = "")
seq <- read_csv("data/JLV_seq.csv")

# PREP #

# rename first column and convert the rest to numeric
# turn NAs into 0
matrix <- matrix %>% 
  rename(Behavior = `...1`) %>% 
  mutate(across(.cols = `1= Posting`:ncol(matrix), ~as.numeric(.))) %>% 
  replace(is.na(.), 0)

# calculate row sums and create proportion table
matrix <- matrix %>% 
  mutate(row_sums = rowSums(.[2:ncol(matrix)])) %>% 
  mutate(across(.cols = c(-Behavior, -row_sums), ~(./row_sums))) %>% 
  select(-row_sums) %>% 
  column_to_rownames("Behavior")
  
# save column names in a vector
states <- colnames(matrix)

## Sequence ##

# clean up ethogram tags and states
sequence <- seq %>% select(ethogram_tag = `Ethogram Tag...10`) %>% 
  mutate(across(where(is.character), ~ na_if(.,"")))

sequence <- sequence %>% 
  mutate(ethogram_tag = replace(ethogram_tag, ethogram_tag == "Rocking +Sniffing", "Rocking + Sniffing"),
         ethogram_tag = replace(ethogram_tag, ethogram_tag == "Jump Yip", "Jump yip"),
         ethogram_tag = replace(ethogram_tag, ethogram_tag == "Moving away", "Retreat"),
         ethogram_tag = replace(ethogram_tag, ethogram_tag == "Moving toward", "Approach"))

states <- tibble(number = seq(1:22),
                 ethogram_tag = c("Posting", "Jump yip", "Approach", "Scratching",
                                  "Pushing dirt", "Rocking", "Sniffing", 
                                  "Rocking + Sniffing", "Approach + Sniffing",
                                  "Retreat", "Jump back", "Escorting", "Foraging",
                                  "Investigating", "Touching face", "Strike",
                                  "Vocalizing", "Head to Head", "Fake foraging",
                                  "Observing", "Greet-kiss", "Running"))

sequence <- full_join(sequence, states)

# matrix with all 22 states
seq22_matrix <- createSequenceMatrix(stringchar = as.vector(sequence$ethogram_tag),
                                     toRowProbs = TRUE)

# reduce to 13 states
sequence13 <- sequence %>% 
  filter(number <= 13)

seq13_matrix <- createSequenceMatrix(stringchar = as.vector(sequence13$ethogram_tag),
                                     toRowProbs = TRUE)
  
# ATTEMPT MARKOV CHAIN #

# need to convert to character matrix??
seq22_matrix <- as.character(seq22_matrix)
mc_fit_22 <- markovchainFit(data = seq22_matrix)
verifyMarkovProperty(seq22_matrix)

seq13_matrix <- as.character(seq13_matrix)
mc_fit_13 <- markovchainFit(seq13_matrix)

verifyMarkovProperty(seq13_matrix)
plot(mc_fit_13)
