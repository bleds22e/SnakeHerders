# Prairie Dog Behavior MCMC #
# EKB; March 16, 2022 #


# PREP #

# packages
library(tidyverse)
library(markovchain)

# load matrix
matrix <- read_csv("data/Snake Whisperers_Pdog Data_JLV.csv", na = "")

# rename first column and convert the rest to numeric
matrix <- matrix %>% 
  rename(Behavior = `...1`) %>% 
  mutate(across(.cols = `1= Posting`:ncol(matrix), ~as.numeric(.)))
matrix <- as.matrix(matrix, nrow = 22, ncol = 22)

# ATTEMPT MARKOV CHAIN #

# need to convert to probability matrix
createSequenceMatrix(matrix,
                     toRowProbs = TRUE,
                     sanitize = TRUE)
