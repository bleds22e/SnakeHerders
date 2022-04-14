# Prairie Dog Behavior MCMC #
# EKB; March 16, 2022 #


# PACKAGES and DATA --------------------------------------------------------####

# packages
library(tidyverse)
library(markovchain)
library(diagram)
library(kequate)

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
  #mutate(row_sums = rowSums(.[2:ncol(matrix)])) %>% 
  #mutate(across(.cols = c(-Behavior, -row_sums), ~(./row_sums))) %>% 
  #dplyr::select(-c(row_sums)) %>% 
  column_to_rownames("Behavior")
  
# save column names in a vector
states <- colnames(matrix)

## Sequence ##

# clean up ethogram tags and states
sequence <- seq %>% dplyr::select(ethogram_tag = `Ethogram Tag...10`) %>% 
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

# TEST MARKOV PROPERTY ------------------------------------------------------####

## CREATE MATRICES ##

# matrix with all 22 states
seq22_matrix <- createSequenceMatrix(stringchar = as.vector(sequence$ethogram_tag),
                                     toRowProbs = FALSE)

# reduce to 13 states
sequence13 <- sequence %>% 
  filter(number <= 13)

seq13_matrix <- createSequenceMatrix(stringchar = as.vector(sequence13$ethogram_tag),
                                     toRowProbs = FALSE)
## ATTEMPT PLOTTING ##
plotmat(seq13_matrix,
        box.size = .05,
        box.cex = .5,
        cex.txt = .5)

## TEST MARKOV PROPERTY ##

# need to convert to character matrix??
seq22_matrix_char <- as.character(seq22_matrix)
mc_fit_22 <- markovchainFit(data = seq22_matrix_char)
verifyMarkovProperty(seq22_matrix_char)

# is significant, meaning is non-random
seq13_matrix_char <- as.character(seq13_matrix)
mc_fit_13 <- markovchainFit(seq13_matrix_char)
verifyMarkovProperty(seq13_matrix_char)


# FREEMAN-TUKEY RESIDUALS --------------------------------------------------####

# run chi-sq test to get expected values matrix
chi_sq <- chisq.test(seq13_matrix)
chi_sq$expected

# attempt FTres() from kequate
FT_residuals <- FTres(seq13_matrix, chi_sq$expected)
FT_res_matrix <- matrix(FT_residuals, nrow = 13, ncol = 13)

# plot matrix
plot.matrix:::plot.matrix(FT_res_matrix)

# chi.sq pairwise? why do these not work?
rmngb::pairwise.chisq.test(x = seq13_matrix)
RVAideMemoire::chisq.multcomp(x = seq13_matrix)
RVAideMemoire::fisher.multcomp(x = seq13_matrix)
rstatix::fisher_test(seq13_matrix, simulate.p.value = TRUE)
rstatix::row_wise_fisher_test(seq13_matrix, simulate.p.values=TRUE)
