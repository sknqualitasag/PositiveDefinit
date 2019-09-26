#' ---
#' title: Constructing multivariate estimates from analyses by parts
#' date:  "`r Sys.Date()`"
#' ---

#' Library
library(magrittr)
library(dplyr)
library(tidyr)
#rm(list=ls())


#' Read all VCE results
s_vce_result <- "../../work/VCE_results.csv"
#s_data_dir <- file.path(here::here(), "inst","extdata")
#s_vce_result <- file.path(s_data_dir, "VCE_results.csv")
tbl_vce <- readr::read_delim(file = s_vce_result, delim = ";")
tbl_vce$estimate[tbl_vce$estimate == "---"] <- "0"
tbl_vce$estimate <- as.numeric(as.character(tbl_vce$estimate))


#' Build Matrix
#' #############################################################
#' 'Matrices: NATRUAL' are the inputs for mix99 in the routine evaluation containing variance and covariance
#tbl_varCovar <- tbl_vce %>% filter(type == "variance" | type == "covariance") %>% select(type,traits,random_effect,estimate)
# Mit Heterosis 2 Kombination fehlen
#tbl_varCovar <- tbl_vce %>% filter(type == "variance" | type == "covariance") %>% filter(model_name == "Het") %>% select(type,traits,random_effect,estimate)
# Ohne Heterosis ! spezifisch zu Sophie
tbl_vce[is.na(tbl_vce[,"model_name"]),"model_name"] <- "Without_Het"
tbl_varCovar <- tbl_vce %>% filter(type == "variance" | type == "covariance") %>% filter(model_name == "Without_Het") %>% select(type,traits,random_effect,estimate)

#' Split traits into trait 1 and trait 2
#+ trait-split
tbl_varCovar <- tbl_varCovar %>% separate(traits, c('trait', 'surrogate'), remove = FALSE)
tbl_varCovar[is.na(tbl_varCovar$surrogate),'surrogate'] <- tbl_varCovar[is.na(tbl_varCovar$surrogate),'trait']

#' change order of trait and surrogate based on alphabetic order of them
idx <- tbl_varCovar[,'trait'] > tbl_varCovar[,'surrogate']
vec_tmp_trait <- tbl_varCovar[idx,'trait']
tbl_varCovar[idx,'trait'] <- tbl_varCovar[idx,'surrogate']
tbl_varCovar[idx,'surrogate'] <- vec_tmp_trait$trait

#' For each traits&random_effect get the mean value of all variants
by_grp <- tbl_varCovar %>% group_by(trait, surrogate, random_effect)
smry <- summarise(by_grp,
                  meanEstimate = mean(estimate, na.rm = TRUE))

#' Prepare matrix to be able to read from tbl
vec_trait_name <- unique(smry$trait)
n_nr_trait <- length(vec_trait_name)
mat_randomEffect <- matrix(0, nrow = n_nr_trait, ncol = n_nr_trait)
rownames(mat_randomEffect) <- vec_trait_name
colnames(mat_randomEffect) <- vec_trait_name

#' Get random effects
vec_randomEffect_name <- unique(smry$random_effect)

resultList <- NULL
for(Z in vec_randomEffect_name){
  # take only values for a random effect
  smry_Z <- smry %>% filter(random_effect == Z)

  # loop over rows of tbl and write elements to matrix
  for (i in 1:nrow(smry_Z)){
    mat_randomEffect[smry_Z$trait[i], smry_Z$surrogate[i]] <- smry_Z$meanEstimate[i]
    mat_randomEffect[smry_Z$surrogate[i], smry_Z$trait[i]] <- smry_Z$meanEstimate[i]
  }
  resultList[[Z]] <- mat_randomEffect
}

#' Check or Transfrom Matrix to insure beeing Positive Definit
#' #############################################################

PDresultList <- NULL
for(Z in vec_randomEffect_name){
  # Optimized function of Schaeffer
  PDresultList[[Z]] <- makePD2(resultList[[Z]])
}

#' Build Variance/Covariance Parameter File for Mix99
#' #############################################################

#' Prepare the different input to build the parameter file
n_nr_randomEffect <- length(vec_randomEffect_name)

#for(Z in vec_randomEffect_name){
  Z <- "animal"
  for(i in 1:n_nr_trait){
    for(j in 1:(n_nr_trait-1)){
      PDresultList[[Z]][[i,j]]
    }
  }
#}

