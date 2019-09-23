### #
### # Constructing multivariate estimates from analyses by parts
### # 2019-09-23 (skn)
### # ----------------------------------------------


### # Library
library(magrittr)
library(dplyr)

### # Read all VCE results
s_vce_result <- "../../work/VCE_resluts.csv"
tbl_vce <- readr::read_delim(file = s_vce_result, delim = ";")
### # Only Sample 1 da meistens status 1
tbl_vce <- tbl_vce %>% filter(data_subset_number == "Sample_1")

### # 'Matrices: NATRUAL' are the inputs for mix99 in the routine evaluation containing variance and covariance
tbl_varCovar <- tbl_vce %>% filter(type == "variance" | type == "covariance") %>% select(type,traits,random_effect,estimate)


### # Build Matrice with Ccc Cca Cfc Cfa Cwc Cwa Cac Caa
by_grp <- tbl_varCovar %>% group_by(traits, random_effect)
smry <- summarise(by_grp,
                  meanEstimate = mean(estimate))
