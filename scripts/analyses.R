#################################
# Prairie Dog Behavior Analysis #
#       EKB; May 2022           #
#################################


# PACKAGES and DATA --------------------------------------------------------####

# Packages #

# install `pacman` package if not installed;
# then load (and install if needed) all packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, markovchain, kequate, Hmisc, corrplot)

# Correlation Table Functions #
source("scripts/correlation_matrix_fxn.R")

# Data #

# load ethogram sequence data from JLV
seq <- read_csv("data/JLV_seq.csv")

# DATA CLEANING ------------------------------------------------------------####

# Clean up ethogram tags and states #

# get required columns, replace blanks with NA, and fill in videoID values
sequence <- seq %>% dplyr::select(videoID = `Video ID`, # will need later
                                  ethogram_tag = `Ethogram Tag...10`) %>% 
  mutate(across(where(is.character), ~ na_if(.,""))) %>% 
  fill(videoID, .direction = "down") %>% 
  filter(videoID != 4795, videoID != 4793, videoID != 4778)

# standardize ethogram_tags
sequence <- sequence %>% 
  mutate(ethogram_tag = replace(ethogram_tag, ethogram_tag == "Rocking +Sniffing", "Rocking + Sniffing"),
         ethogram_tag = replace(ethogram_tag, ethogram_tag == "Jump Yip", "Jump yip"),
         ethogram_tag = replace(ethogram_tag, ethogram_tag == "Moving away", "Retreat"),
         ethogram_tag = replace(ethogram_tag, ethogram_tag == "Moving toward", "Approach"),
         ethogram_tag = replace(ethogram_tag, ethogram_tag == "Pushing dirt", "Foreleg drumming"))

# create df with numbers for each ethogram_tag
states <- tibble(number = seq(1:22),
                 ethogram_tag = c("Posting", "Jump yip", "Approach", "Scratching",
                                  "Foreleg drumming", "Rocking", "Sniffing", 
                                  "Rocking + Sniffing", "Approach + Sniffing",
                                  "Retreat", "Jump back", "Escorting", "Foraging",
                                  "Investigating", "Touching face", "Strike",
                                  "Vocalizing", "Head to Head", "Fake foraging",
                                  "Observing", "Greet-kiss", "Running"))

# merge df and remove NAs from ethogram_tag or number columns
sequence <- full_join(sequence, states) %>% 
  drop_na(ethogram_tag, number)
  

# Given the way the sequence data were entered, the last action from one obs and
# the first action from the next would be included in the matrix. By adding an
# NA back in at the start of every video, those action pairs will not be included
# in the the matrix.

sequence <- sequence %>% 
  group_split(videoID) %>%      # create multiple dfs based on videoID column
  map_dfr(~ .x %>% add_row())   # add a row with NAs to the end of each dataframe 
                                # the map_dfr function then recombines into one df

# remove any doubled behaviors
sequence <- sequence %>%  
  filter(ethogram_tag != lag(ethogram_tag))

# reduce to 13 key states and remove any new doubles created
sequence12 <- sequence %>% 
  filter(number <= 12) %>% 
  filter(ethogram_tag != lag(ethogram_tag))


# TEST MARKOV PROPERTY ------------------------------------------------------####

# Create Matrices #

# matrix with all 22 states
seq22_matrix <- createSequenceMatrix(stringchar = as.vector(sequence$ethogram_tag),
                                     toRowProbs = FALSE)
# matrix with core 12 states
seq12_matrix <- createSequenceMatrix(stringchar = as.vector(sequence12$ethogram_tag),
                                     toRowProbs = FALSE)


# Test Markov Property #

# full matrix
seq22_matrix_char <- as.character(seq22_matrix) # convert to character matrix
markov_prop22 <- (verifyMarkovProperty(seq22_matrix_char))

# reduced matrix
seq12_matrix_char <- as.character(seq12_matrix)
markov_prop12 <- (verifyMarkovProperty(seq12_matrix_char))


# TESTING RESIDUALS --------------------------------------------------------####

# Focusing on reduced matrix #

# run chi-sq test 
chi_sq <- chisq.test(seq12_matrix)
chi_sq

chi_sq_psim <- chisq.test(seq12_matrix, simulate.p.value = TRUE)
chi_sq_psim

# get expected values to compare to observed
expected <- chi_sq$expected

# get Freeman-Tukey residuals (observed vs. expected values)
FT_res <- kequate::FTres(seq12_matrix, expected)

# turn into a matrix of FT values
FT_res_matrix <- matrix(FT_res, nrow = 12, ncol = 12)
rownames(FT_res_matrix) <- rownames(seq12_matrix)
colnames(FT_res_matrix) <- colnames(seq12_matrix)

# get correlations
rcorr_FTres <- rcorr(FT_res_matrix, type = "pearson")

# save significant values
rcorr_FT_sig <- as.data.frame(rcorr_FTres$P, row.names = colnames(FT_res_matrix))
write.csv(rcorr_FT_sig, "data/rcorr_FTres_sig.csv")

# plot correlation
png("plots/FT_residuals.png", res = 300, width = 200, height = 200, units = "mm")
corrplot(rcorr_FTres$r, is.corr = FALSE, diag = FALSE, 
         type = "lower", tl.col = "black",)
dev.off()

# correlation value table with significance values marked with **
corr_matrix <- correlation_matrix(as.data.frame(FT_res_matrix), 
                                  replace_diagonal = TRUE, use = "lower")
corr_matrix <- as.data.frame(corr_matrix, row.names = colnames(corr_matrix))
write.csv(corr_matrix, "data/corr_matrix_FTres.csv")
