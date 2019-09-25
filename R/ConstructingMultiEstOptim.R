### #
### # Constructing multivariate estimates from analyses by parts
### # 2019-09-23 (skn)
### # ------------------------------------------------------------


### # Library
library(magrittr)
library(dplyr)
library(tidyr)

### # Read all VCE results
 s_vce_result <- "../../work/VCE_resluts.csv"
#s_data_dir <- file.path(here::here(), "inst","extdata")
#s_vce_result <- file.path(s_data_dir, "VCE_resluts.csv")
tbl_vce <- readr::read_delim(file = s_vce_result, delim = ";")
tbl_vce$estimate[tbl_vce$estimate == "---"] <- "0"
tbl_vce$estimate <- as.numeric(as.character(tbl_vce$estimate))

### # 'Matrices: NATRUAL' are the inputs for mix99 in the routine evaluation containing variance and covariance
tbl_varCovar <- tbl_vce %>% filter(type == "variance" | type == "covariance") %>% select(type,traits,random_effect,estimate)

#' Split traits into trait 1 and trait 2
#+ trait-split
tbl_varCovar <- tbl_varCovar %>% separate(traits, c('trait', 'surrogate'), remove = FALSE)
tbl_varCovar[is.na(tbl_varCovar$surrogate),'surrogate'] <- tbl_varCovar[is.na(tbl_varCovar$surrogate),'trait']
tbl_varCovar

#' change order of trait and surrogate based on alphabetic order of them
vec_tmp_trait <- tbl_varCovar[tbl_varCovar[,'trait'] > tbl_varCovar[,'surrogate'],'trait']
tbl_varCovar[tbl_varCovar[,'trait'] > tbl_varCovar[,'surrogate'],'trait'] <- tbl_varCovar[tbl_varCovar[,'trait'] > tbl_varCovar[,'surrogate'],'surrogate']
tbl_varCovar[tbl_varCovar[,'trait'] > tbl_varCovar[,'surrogate'],'surrogate'] <- vec_tmp_trait


### # ###############################################################
### # Build Matrice with for example Ccc Cca Cfc Cfa Cwc Cwa Cac Caa
### # ###############################################################

### # For each traits&random_effect get the mean value of all variants
by_grp <- tbl_varCovar %>% group_by(trait, surrogate, random_effect)
smry <- summarise(by_grp,
                  meanEstimate = mean(estimate, na.rm = TRUE))

#' take only values for random effect animal
smry_animal <- smry %>% filter(random_effect == 'animal')

#' Prepare matrix to be able to read from tbl
vec_trait_name <- unique(smry$trait)
n_nr_trait <- length(vec_trait_name)
mat_animal <- matrix(0, nrow = n_nr_trait, ncol = n_nr_trait)
rownames(mat_animal) <- vec_trait_name
colnames(mat_animal) <- vec_trait_name

#' loop over rows of tbl and write elements to matrix
for (i in 1:nrow(smry_animal)){
  mat_animal[smry_animal$trait[i], smry_animal$surrogate[i]] <- smry_animal$meanEstimate[i]
  mat_animal[smry_animal$surrogate[i], smry_animal$trait[i]] <- smry_animal$meanEstimate[i]
}

#' Get random effects
vec_randomEffect_name <- unique(smry$random_effect)
