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
rmngb::pairwise.chisq.test(x = FT_res_matrix)
RVAideMemoire::chisq.multcomp(x = seq13_matrix)
RVAideMemoire::fisher.multcomp(x = seq13_matrix)
rstatix::fisher_test(seq13_matrix, simulate.p.value = TRUE)
rstatix::row_wise_fisher_test(seq13_matrix, simulate.p.values=TRUE)

library(corrplot)
corrplot(FT_res_matrix, is.corr = FALSE)

correlation_matrix

#' correlation_matrix
#' Creates a publication-ready / formatted correlation matrix, using `Hmisc::rcorr` in the backend.
#'
#' @param df dataframe; containing numeric and/or logical columns to calculate correlations for
#' @param type character; specifies the type of correlations to compute; gets passed to `Hmisc::rcorr`; options are `"pearson"` or `"spearman"`; defaults to `"pearson"`
#' @param digits integer/double; number of decimals to show in the correlation matrix; gets passed to `formatC`; defaults to `3`
#' @param decimal.mark character; which decimal.mark to use; gets passed to `formatC`; defaults to `.`
#' @param use character; which part of the correlation matrix to display; options are `"all"`, `"upper"`, `"lower"`; defaults to `"all"`
#' @param show_significance boolean; whether to add `*` to represent the significance levels for the correlations; defaults to `TRUE`
#' @param replace_diagonal boolean; whether to replace the correlations on the diagonal; defaults to `FALSE`
#' @param replacement character; what to replace the diagonal and/or upper/lower triangles with; defaults to `""` (empty string)
#'
#' @return a correlation matrix
#' @export
#'
#' @examples
#' `correlation_matrix(iris)`
#' `correlation_matrix(mtcars)`
correlation_matrix <- function(df, 
                               type = "pearson",
                               digits = 3, 
                               decimal.mark = ".",
                               use = "all", 
                               show_significance = TRUE, 
                               replace_diagonal = FALSE, 
                               replacement = ""){
  
  # check arguments
  stopifnot({
    is.numeric(digits)
    digits >= 0
    use %in% c("all", "upper", "lower")
    is.logical(replace_diagonal)
    is.logical(show_significance)
    is.character(replacement)
  })
  # we need the Hmisc package for this
  require(Hmisc)
  
  # retain only numeric and boolean columns
  isNumericOrBoolean = vapply(df, function(x) is.numeric(x) | is.logical(x), logical(1))
  if (sum(!isNumericOrBoolean) > 0) {
    cat('Dropping non-numeric/-boolean column(s):', paste(names(isNumericOrBoolean)[!isNumericOrBoolean], collapse = ', '), '\n\n')
  }
  df = df[isNumericOrBoolean]
  
  # transform input data frame to matrix
  x <- as.matrix(df)
  
  # run correlation analysis using Hmisc package
  correlation_matrix <- Hmisc::rcorr(x, type = type)
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value 
  
  # transform correlations to specific character format
  Rformatted = formatC(R, format = 'f', digits = digits, decimal.mark = decimal.mark)
  
  # if there are any negative numbers, we want to put a space before the positives to align all
  if (sum(!is.na(R) & R < 0) > 0) {
    Rformatted = ifelse(!is.na(R) & R > 0, paste0(" ", Rformatted), Rformatted)
  }
  
  # add significance levels if desired
  if (show_significance) {
    # define notions for significance levels; spacing is important.
    stars <- ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", ""))))
    Rformatted = paste0(Rformatted, stars)
  }
  
  # make all character strings equally long
  max_length = max(nchar(Rformatted))
  Rformatted = vapply(Rformatted, function(x) {
    current_length = nchar(x)
    difference = max_length - current_length
    return(paste0(x, paste(rep(" ", difference), collapse = ''), sep = ''))
  }, FUN.VALUE = character(1))
  
  # build a new matrix that includes the formatted correlations and their significance stars
  Rnew <- matrix(Rformatted, ncol = ncol(x))
  rownames(Rnew) <- colnames(Rnew) <- colnames(x)
  
  # replace undesired values
  if (use == 'upper') {
    Rnew[lower.tri(Rnew, diag = replace_diagonal)] <- replacement
  } else if (use == 'lower') {
    Rnew[upper.tri(Rnew, diag = replace_diagonal)] <- replacement
  } else if (replace_diagonal) {
    diag(Rnew) <- replacement
  }
  
  return(Rnew)
}

correlation_matrix(as.data.frame(seq13_matrix))
