#' ---
#' title: Bending and Production of Mix99 Paramater-File
#' date:  "`r Sys.Date()`"
#' ---

#' Get the function makePD2() in Rscript MakePD.R
source("MakePD.R")

#' Run rscript via console
#+arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=2) stop("didn't receive 2 arguments")
s_vce_result <- args[1]
s_result_file <- args[2]
if(!file.exists(s_vce_result)) stop("first argument isn't an existing file")

#' Library for this script
#+ library
suppressPackageStartupMessages(if(! require("magrittr")) {
  install.packages("magrittr", repos="https://stat.ethz.ch/CRAN/")
  require("magrittr")
})
suppressPackageStartupMessages(if(! require("dplyr")) {
  install.packages("dplyr", repos="https://stat.ethz.ch/CRAN/")
  require("dplyr")
})
suppressPackageStartupMessages(if(! require("tidyr")) {
  install.packages("tidyr", repos="https://stat.ethz.ch/CRAN/")
  require("tidyr")
})
#rm(list=ls())

#' Read all VCE results
#+ read_transform_misssing_o_type
#s_vce_result <- '/Volumes/data_projekte/projekte/singularity_data_zws_gslim/muku_CarcassVK/work/VCE_results.csv'
tbl_vce <- readr::read_delim(file = s_vce_result, delim = ";")
tbl_vce$estimate[tbl_vce$estimate == "---"] <- "0"
tbl_vce$estimate <- as.numeric(as.character(tbl_vce$estimate))

#' Build Matrix
#' ##########################################################

#' 'Matrices: NATRUAL' are the inputs for mix99 in the routine evaluation containing variance and covariance
#+ varianceCovariance
#tbl_varCovar <- tbl_vce %>% filter(type == "variance" | type == "covariance") %>% select(type,traits,random_effect,estimate)
# Mit Heterosis 2 Kombination fehlen
#tbl_varCovar <- tbl_vce %>% filter(type == "variance" | type == "covariance") %>% filter(model_name == "Het") %>% select(type,traits,random_effect,estimate)
# Ohne Heterosis ! spezifisch zu Sophie
tbl_vce[is.na(tbl_vce[,"model_name"]),"model_name"] <- "Without_Het"
tbl_varCovar <- tbl_vce %>% filter(type == "variance" | type == "covariance") %>% filter(model_name == "Without_Het") %>% select(type,traits,random_effect,estimate)
#tbl_varCovar <- tbl_vce %>% filter(type == "heritability" | type == "genetic_correlation") %>% filter(model_name == "Without_Het") %>% select(type,traits,random_effect,estimate)

#' Split traits into trait 1 and trait 2
#+ trait-split
tbl_varCovar <- tbl_varCovar %>% separate(traits, c('trait', 'surrogate'), remove = FALSE)
tbl_varCovar[is.na(tbl_varCovar$surrogate),'surrogate'] <- tbl_varCovar[is.na(tbl_varCovar$surrogate),'trait']

#' change order of trait and surrogate based on alphabetic order of them
#+ switch_trait_surrogate
idx <- tbl_varCovar[,'trait'] > tbl_varCovar[,'surrogate']
vec_tmp_trait <- tbl_varCovar[idx,'trait']
tbl_varCovar[idx,'trait'] <- tbl_varCovar[idx,'surrogate']
tbl_varCovar[idx,'surrogate'] <- vec_tmp_trait$trait

#' For each traits&random_effect get the mean value of all variants
#+ mean_per_random_effect
by_grp <- tbl_varCovar %>% group_by(trait, surrogate, random_effect)
smry <- summarise(by_grp,
                  meanEstimate = mean(estimate, na.rm = TRUE))

#' Prepare matrix to be able to read from tbl
#+ prepare_matrix
vec_trait_name <- unique(smry$trait)
n_nr_trait <- length(vec_trait_name)
mat_randomEffect <- matrix(0, nrow = n_nr_trait, ncol = n_nr_trait)
rownames(mat_randomEffect) <- vec_trait_name
colnames(mat_randomEffect) <- vec_trait_name

#' Get random effects
#+ random_list
vec_randomEffect_name <- unique(smry$random_effect)

resultList <- NULL
for(Z in vec_randomEffect_name){
  # take only values for a random effect
  smry_Z <- smry %>% filter(random_effect == Z)

  # loop over rows of tbl and write elements to matrix
  for (i in 1:nrow(smry_Z)){
    mat_randomEffect[smry_Z$trait[i], smry_Z$surrogate[i]] <- round(smry_Z$meanEstimate[i],3)
    mat_randomEffect[smry_Z$surrogate[i], smry_Z$trait[i]] <- round(smry_Z$meanEstimate[i],3)

#    mat_randomEffect[smry_Z$trait[i], smry_Z$surrogate[i]] <- smry_Z$meanEstimate[i]
#    mat_randomEffect[smry_Z$surrogate[i], smry_Z$trait[i]] <- smry_Z$meanEstimate[i]
  }
  resultList[[Z]] <- mat_randomEffect
}

#' Check or Transfrom Matrix to insure beeing Positive Definit
#' ##########################################################"
#+ positive_definite_function_makePD2

PDresultList <- NULL
for(Z in vec_randomEffect_name){
  # Optimized function of Schaeffer
  PDresultList[[Z]] <- round(makePD2(resultList[[Z]]),digits = 3)
#  PDresultList[[Z]] <- makePD2(resultList[[Z]])
}

#' Build Variance/Covariance Parameter File for Mix99
#' #############################################################

#' Prepare the different input to build the parameter file
#+ length_random_effect
n_nr_randomEffect <- length(vec_randomEffect_name)

#' Animal and residual should have a specific order in mix99
#+ order_random_effect
vec_random_effect_req <- c("animal", "residual")
if (!all(vec_random_effect_req %in% vec_randomEffect_name))
  stop(" * ERROR: Required random effects animal and residual are not both in list of random effects")
vec_random_effects_mand <- setdiff(vec_randomEffect_name, vec_random_effect_req)
vec_random_effect_order <- c(vec_random_effects_mand, vec_random_effect_req)

#' Build Variance/Covariance Parameter-File for Mix99
#+ create_var_file
#s_result_file <- 'mix99_carcass.var'
if (file.exists(s_result_file))
  file.remove(s_result_file)
idx_rand_eff <- 1
for(Z in vec_random_effect_order){
  for(i in 1:n_nr_trait){
    for(j in i:n_nr_trait){
      cat(idx_rand_eff, i, j, format(PDresultList[[Z]][[i,j]], scientific = FALSE), "\n", file = s_result_file, append = TRUE)
#      cat(idx_rand_eff, i, j, round(PDresultList[[Z]][[i,j]], digits = 3), "\n", file = s_result_file, append = TRUE)
    }
  }
  idx_rand_eff <- idx_rand_eff + 1
}

