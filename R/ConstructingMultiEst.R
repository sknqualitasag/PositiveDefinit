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
tbl_vce$estimate <- as.numeric(as.character(tbl_vce$estimate))

### # 'Matrices: NATRUAL' are the inputs for mix99 in the routine evaluation containing variance and covariance
tbl_varCovar <- tbl_vce %>% filter(type == "variance" | type == "covariance") %>% select(type,traits,random_effect,estimate,trait_combination)


### # Build Matrice with Ccc Cca Cfc Cfa Cwc Cwa Cac Caa
### # Step 1: For each traits&random_effect get the mean value of all variants
by_grp <- tbl_varCovar %>% group_by(traits, random_effect)
smry <- summarise(by_grp,
                  meanEstimate = mean(estimate, na.rm = TRUE))
### # Step 2: Split in 3 dataset: animal, herdyear, residual
tbl_animal <- smry %>% filter(random_effect == "animal")
tbl_herdyear <- smry %>% filter(random_effect == "herdyear")
tbl_residual <- smry %>% filter(random_effect == "residual")
### # Step 3:
