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
  fill(videoID, .direction = "down")

# standardize ethogram_tags
sequence <- sequence %>% 
  mutate(ethogram_tag = replace(ethogram_tag, ethogram_tag == "Rocking +Sniffing", "Rocking + Sniffing"),
         ethogram_tag = replace(ethogram_tag, ethogram_tag == "Jump Yip", "Jump yip"),
         ethogram_tag = replace(ethogram_tag, ethogram_tag == "Moving away", "Retreat"),
         ethogram_tag = replace(ethogram_tag, ethogram_tag == "Moving toward", "Approach"))

# create df with numbers for each ethogram_tag
states <- tibble(number = seq(1:22),
                 ethogram_tag = c("Posting", "Jump yip", "Approach", "Scratching",
                                  "Pushing dirt", "Rocking", "Sniffing", 
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
sequence13 <- sequence %>% 
  filter(number <= 13) %>% 
  filter(ethogram_tag != lag(ethogram_tag))


# TEST MARKOV PROPERTY ------------------------------------------------------####

# Create Matrices #

# matrix with all 22 states
seq22_matrix <- createSequenceMatrix(stringchar = as.vector(sequence$ethogram_tag),
                                     toRowProbs = FALSE)
# matrix with core 13 states
seq13_matrix <- createSequenceMatrix(stringchar = as.vector(sequence13$ethogram_tag),
                                     toRowProbs = FALSE)


# Test Markov Property #

# full matrix
seq22_matrix_char <- as.character(seq22_matrix) # convert to character matrix
markov_prop22 <- (verifyMarkovProperty(seq22_matrix_char))

# reduced matrix
seq13_matrix_char <- as.character(seq13_matrix)
markov_prop13 <- (verifyMarkovProperty(seq13_matrix_char))


# TESTING RESIDUALS --------------------------------------------------------####

# Focusing on reduced matrix #

# run chi-sq test 
chi_sq <- chisq.test(seq13_matrix, simulate.p.value = TRUE)
chi_sq

# get expected values to compare to observed
expected <- chi_sq$expected

# get Freeman-Tukey residuals (observed vs. expected values)
FT_res <- kequate::FTres(seq13_matrix, expected)

# turn into a matrix of FT values
FT_res_matrix <- matrix(FT_res, nrow = 13, ncol = 13)
rownames(FT_res_matrix) <- rownames(seq13_matrix)
colnames(FT_res_matrix) <- colnames(seq13_matrix)

# get correlations
rcorr_FTres <- rcorr(FT_res_matrix, type = "pearson")

# plot correlation
dev.copy(png, "plots/FT_residuals.png")
corrplot(rcorr_FTres$r, is.corr = TRUE, diag = FALSE, 
         type = "upper", tl.col = "black")
dev.off()

# pull out significance values
corr_matrix <- correlation_matrix(as.data.frame(FT_res_matrix), replace_diagonal = )
save_correlation_matrix(as.data.frame(FT_res_matrix),
                        "data/correlation_matrix_FTres.csv",
                        digits = 3, use = "lower", replace_diagonal = TRUE)
