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
tbl_varCovar <- tbl_vce %>% filter(type == "variance" | type == "covariance") %>% select(type,traits,random_effect,estimate)


### # Build Matrice with Ccc Cca Cfc Cfa Cwc Cwa Cac Caa
### # Step 1: For each traits&random_effect get the mean value of all variants
by_grp <- tbl_varCovar %>% group_by(traits, random_effect)
smry <- summarise(by_grp,
                  meanEstimate = mean(estimate, na.rm = TRUE))
### # Step 2: Split in 3 dataset: animal, herdyear, residual
tbl_animal <- smry %>% filter(random_effect == "animal")
tbl_herdyear <- smry %>% filter(random_effect == "herdyear")
tbl_residual <- smry %>% filter(random_effect == "residual")


### # Step 3: Build matrice
### # For animal
animalVarCovar <- matrix(0, nrow = 8, ncol = 8)
#1. Zeile
animalVarCovar[1,1] <- tbl_animal%>%filter(traits == "ccc")%>%magrittr::extract2("meanEstimate")
animalVarCovar[1,2] <- tbl_animal%>%filter(traits == "ccc+cca")%>%magrittr::extract2("meanEstimate")
#animalVarCovar[1,3] <- tbl_animal%>%filter(traits == "ccc+cfc")%>%magrittr::extract2("meanEstimate") #leer
animalVarCovar[1,3] <- 0.01613 #Value of Sample_1 CacCccCfcCwc
animalVarCovar[1,4] <- tbl_animal%>%filter(traits == "ccc+cfa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[1,5] <- tbl_animal%>%filter(traits == "ccc+cwc")%>%magrittr::extract2("meanEstimate")
animalVarCovar[1,6] <- tbl_animal%>%filter(traits == "ccc+cwa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[1,7] <- tbl_animal%>%filter(traits == "ccc+cac")%>%magrittr::extract2("meanEstimate")
animalVarCovar[1,8] <- tbl_animal%>%filter(traits == "ccc+caa")%>%magrittr::extract2("meanEstimate")
#2. Zeile
animalVarCovar[2,1] <- tbl_animal%>%filter(traits == "ccc+cca")%>%magrittr::extract2("meanEstimate")
animalVarCovar[2,2] <- tbl_animal%>%filter(traits == "cca")%>%magrittr::extract2("meanEstimate")
animalVarCovar[2,3] <- tbl_animal%>%filter(traits == "cca+cfc")%>%magrittr::extract2("meanEstimate")
animalVarCovar[2,4] <- tbl_animal%>%filter(traits == "cca+cfa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[2,5] <- tbl_animal%>%filter(traits == "cca+cwc")%>%magrittr::extract2("meanEstimate")
animalVarCovar[2,6] <- tbl_animal%>%filter(traits == "cca+cwa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[2,7] <- tbl_animal%>%filter(traits == "cca+cac")%>%magrittr::extract2("meanEstimate")
animalVarCovar[2,8] <- tbl_animal%>%filter(traits == "cca+caa")%>%magrittr::extract2("meanEstimate")
#3. Zeile
#animalVarCovar[3,1] <- tbl_animal%>%filter(traits == "cfc+ccc")%>%magrittr::extract2("meanEstimate") #leer
animalVarCovar[3,1] <- 0.01613 #Value of Sample_1 CacCccCfcCwc
animalVarCovar[3,2] <- tbl_animal%>%filter(traits == "cca+cfc")%>%magrittr::extract2("meanEstimate")
animalVarCovar[3,3] <- tbl_animal%>%filter(traits == "cfc")%>%magrittr::extract2("meanEstimate")
animalVarCovar[3,4] <- tbl_animal%>%filter(traits == "cfc+cfa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[3,5] <- tbl_animal%>%filter(traits == "cfc+cwc")%>%magrittr::extract2("meanEstimate")
animalVarCovar[3,6] <- tbl_animal%>%filter(traits == "cfc+cwa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[3,7] <- tbl_animal%>%filter(traits == "cfc+cac")%>%magrittr::extract2("meanEstimate")
animalVarCovar[3,8] <- tbl_animal%>%filter(traits == "cfc+caa")%>%magrittr::extract2("meanEstimate")
#4. Zeile
animalVarCovar[4,1] <- tbl_animal%>%filter(traits == "ccc+cfa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[4,2] <- tbl_animal%>%filter(traits == "cca+cfa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[4,3] <- tbl_animal%>%filter(traits == "cfc+cfa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[4,4] <- tbl_animal%>%filter(traits == "cfa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[4,5] <- tbl_animal%>%filter(traits == "cfa+cwc")%>%magrittr::extract2("meanEstimate")
animalVarCovar[4,6] <- tbl_animal%>%filter(traits == "cfa+cwa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[4,7] <- tbl_animal%>%filter(traits == "cfa+cac")%>%magrittr::extract2("meanEstimate")
animalVarCovar[4,8] <- tbl_animal%>%filter(traits == "cfa+caa")%>%magrittr::extract2("meanEstimate")
#5. Zeile
animalVarCovar[5,1] <- tbl_animal%>%filter(traits == "ccc+cwc")%>%magrittr::extract2("meanEstimate")
animalVarCovar[5,2] <- tbl_animal%>%filter(traits == "cca+cwc")%>%magrittr::extract2("meanEstimate")
animalVarCovar[5,3] <- tbl_animal%>%filter(traits == "cfc+cwc")%>%magrittr::extract2("meanEstimate")
animalVarCovar[5,4] <- tbl_animal%>%filter(traits == "cfa+cwc")%>%magrittr::extract2("meanEstimate")
animalVarCovar[5,5] <- tbl_animal%>%filter(traits == "cwc")%>%magrittr::extract2("meanEstimate")
animalVarCovar[5,6] <- tbl_animal%>%filter(traits == "cwc+cwa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[5,7] <- tbl_animal%>%filter(traits == "cac+cwc")%>%magrittr::extract2("meanEstimate")
animalVarCovar[5,8] <- tbl_animal%>%filter(traits == "caa+cwc")%>%magrittr::extract2("meanEstimate")
#6. Zeile
animalVarCovar[6,1] <- tbl_animal%>%filter(traits == "ccc+cwa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[6,2] <- tbl_animal%>%filter(traits == "cca+cwa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[6,3] <- tbl_animal%>%filter(traits == "cfc+cwa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[6,4] <- tbl_animal%>%filter(traits == "cfa+cwa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[6,5] <- tbl_animal%>%filter(traits == "cwc+cwa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[6,6] <- tbl_animal%>%filter(traits == "cwa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[6,7] <- tbl_animal%>%filter(traits == "cac+cwa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[6,8] <- tbl_animal%>%filter(traits == "caa+cwa")%>%magrittr::extract2("meanEstimate")
#7. Zeile
animalVarCovar[7,1] <- tbl_animal%>%filter(traits == "ccc+cac")%>%magrittr::extract2("meanEstimate")
animalVarCovar[7,2] <- tbl_animal%>%filter(traits == "cca+cac")%>%magrittr::extract2("meanEstimate")
animalVarCovar[7,3] <- tbl_animal%>%filter(traits == "cfc+cac")%>%magrittr::extract2("meanEstimate")
animalVarCovar[7,4] <- tbl_animal%>%filter(traits == "cfa+cac")%>%magrittr::extract2("meanEstimate")
animalVarCovar[7,5] <- tbl_animal%>%filter(traits == "cac+cwc")%>%magrittr::extract2("meanEstimate")
animalVarCovar[7,6] <- tbl_animal%>%filter(traits == "cac+cwa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[7,7] <- tbl_animal%>%filter(traits == "cac")%>%magrittr::extract2("meanEstimate")
animalVarCovar[7,8] <- tbl_animal%>%filter(traits == "cac+caa")%>%magrittr::extract2("meanEstimate")
#8. Zeile
animalVarCovar[8,1] <- tbl_animal%>%filter(traits == "ccc+caa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[8,2] <- tbl_animal%>%filter(traits == "cca+caa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[8,3] <- tbl_animal%>%filter(traits == "cfc+caa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[8,4] <- tbl_animal%>%filter(traits == "cfa+caa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[8,5] <- tbl_animal%>%filter(traits == "caa+cwc")%>%magrittr::extract2("meanEstimate")
animalVarCovar[8,6] <- tbl_animal%>%filter(traits == "caa+cwa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[8,7] <- tbl_animal%>%filter(traits == "cac+caa")%>%magrittr::extract2("meanEstimate")
animalVarCovar[8,8] <- tbl_animal%>%filter(traits == "caa")%>%magrittr::extract2("meanEstimate")

animalVarCovar
# not positive definite
matrixcalc::is.positive.definite(animalVarCovar)
eigen(animalVarCovar)
# transform to positive definite
PD_animalVarCovar <- makePD2(animalVarCovar)
eigen(PD_animalVarCovar)


### # For herdyear
herdyearVarCovar <- matrix(0, nrow = 8, ncol = 8)
#1. Zeile
herdyearVarCovar[1,1] <- tbl_herdyear%>%filter(traits == "ccc")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[1,2] <- tbl_herdyear%>%filter(traits == "ccc+cca")%>%magrittr::extract2("meanEstimate")
#animalVarCovar[1,3] <- tbl_animal%>%filter(traits == "ccc+cfc")%>%magrittr::extract2("meanEstimate") #leer
herdyearVarCovar[1,3] <- 0.019362 #Value of Sample_1 CacCccCfcCwc
herdyearVarCovar[1,4] <- tbl_herdyear%>%filter(traits == "ccc+cfa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[1,5] <- tbl_herdyear%>%filter(traits == "ccc+cwc")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[1,6] <- tbl_herdyear%>%filter(traits == "ccc+cwa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[1,7] <- tbl_herdyear%>%filter(traits == "ccc+cac")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[1,8] <- tbl_herdyear%>%filter(traits == "ccc+caa")%>%magrittr::extract2("meanEstimate")
#2. Zeile
herdyearVarCovar[2,1] <- tbl_herdyear%>%filter(traits == "ccc+cca")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[2,2] <- tbl_herdyear%>%filter(traits == "cca")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[2,3] <- tbl_herdyear%>%filter(traits == "cca+cfc")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[2,4] <- tbl_herdyear%>%filter(traits == "cca+cfa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[2,5] <- tbl_herdyear%>%filter(traits == "cca+cwc")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[2,6] <- tbl_herdyear%>%filter(traits == "cca+cwa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[2,7] <- tbl_herdyear%>%filter(traits == "cca+cac")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[2,8] <- tbl_herdyear%>%filter(traits == "cca+caa")%>%magrittr::extract2("meanEstimate")
#3. Zeile
#animalVarCovar[3,1] <- tbl_animal%>%filter(traits == "cfc+ccc")%>%magrittr::extract2("meanEstimate") #leer
herdyearVarCovar[3,1] <- 0.019362 #Value of Sample_1 CacCccCfcCwc
herdyearVarCovar[3,2] <- tbl_herdyear%>%filter(traits == "cca+cfc")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[3,3] <- tbl_herdyear%>%filter(traits == "cfc")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[3,4] <- tbl_herdyear%>%filter(traits == "cfc+cfa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[3,5] <- tbl_herdyear%>%filter(traits == "cfc+cwc")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[3,6] <- tbl_herdyear%>%filter(traits == "cfc+cwa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[3,7] <- tbl_herdyear%>%filter(traits == "cfc+cac")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[3,8] <- tbl_herdyear%>%filter(traits == "cfc+caa")%>%magrittr::extract2("meanEstimate")
#4. Zeile
herdyearVarCovar[4,1] <- tbl_herdyear%>%filter(traits == "ccc+cfa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[4,2] <- tbl_herdyear%>%filter(traits == "cca+cfa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[4,3] <- tbl_herdyear%>%filter(traits == "cfc+cfa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[4,4] <- tbl_herdyear%>%filter(traits == "cfa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[4,5] <- tbl_herdyear%>%filter(traits == "cfa+cwc")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[4,6] <- tbl_herdyear%>%filter(traits == "cfa+cwa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[4,7] <- tbl_herdyear%>%filter(traits == "cfa+cac")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[4,8] <- tbl_herdyear%>%filter(traits == "cfa+caa")%>%magrittr::extract2("meanEstimate")
#5. Zeile
herdyearVarCovar[5,1] <- tbl_herdyear%>%filter(traits == "ccc+cwc")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[5,2] <- tbl_herdyear%>%filter(traits == "cca+cwc")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[5,3] <- tbl_herdyear%>%filter(traits == "cfc+cwc")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[5,4] <- tbl_herdyear%>%filter(traits == "cfa+cwc")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[5,5] <- tbl_herdyear%>%filter(traits == "cwc")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[5,6] <- tbl_herdyear%>%filter(traits == "cwc+cwa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[5,7] <- tbl_herdyear%>%filter(traits == "cac+cwc")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[5,8] <- tbl_herdyear%>%filter(traits == "caa+cwc")%>%magrittr::extract2("meanEstimate")
#6. Zeile
herdyearVarCovar[6,1] <- tbl_herdyear%>%filter(traits == "ccc+cwa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[6,2] <- tbl_herdyear%>%filter(traits == "cca+cwa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[6,3] <- tbl_herdyear%>%filter(traits == "cfc+cwa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[6,4] <- tbl_herdyear%>%filter(traits == "cfa+cwa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[6,5] <- tbl_herdyear%>%filter(traits == "cwc+cwa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[6,6] <- tbl_herdyear%>%filter(traits == "cwa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[6,7] <- tbl_herdyear%>%filter(traits == "cac+cwa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[6,8] <- tbl_herdyear%>%filter(traits == "caa+cwa")%>%magrittr::extract2("meanEstimate")
#7. Zeile
herdyearVarCovar[7,1] <- tbl_herdyear%>%filter(traits == "ccc+cac")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[7,2] <- tbl_herdyear%>%filter(traits == "cca+cac")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[7,3] <- tbl_herdyear%>%filter(traits == "cfc+cac")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[7,4] <- tbl_herdyear%>%filter(traits == "cfa+cac")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[7,5] <- tbl_herdyear%>%filter(traits == "cac+cwc")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[7,6] <- tbl_herdyear%>%filter(traits == "cac+cwa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[7,7] <- tbl_herdyear%>%filter(traits == "cac")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[7,8] <- tbl_herdyear%>%filter(traits == "cac+caa")%>%magrittr::extract2("meanEstimate")
#8. Zeile
herdyearVarCovar[8,1] <- tbl_herdyear%>%filter(traits == "ccc+caa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[8,2] <- tbl_herdyear%>%filter(traits == "cca+caa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[8,3] <- tbl_herdyear%>%filter(traits == "cfc+caa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[8,4] <- tbl_herdyear%>%filter(traits == "cfa+caa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[8,5] <- tbl_herdyear%>%filter(traits == "caa+cwc")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[8,6] <- tbl_herdyear%>%filter(traits == "caa+cwa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[8,7] <- tbl_herdyear%>%filter(traits == "cac+caa")%>%magrittr::extract2("meanEstimate")
herdyearVarCovar[8,8] <- tbl_herdyear%>%filter(traits == "caa")%>%magrittr::extract2("meanEstimate")

herdyearVarCovar
# not positive definite
matrixcalc::is.positive.definite(herdyearVarCovar)
eigen(herdyearVarCovar)

### # For residual
residualVarCovar <- matrix(0, nrow = 8, ncol = 8)
# Prob with NaN deswegen auskommentiert
#1. Zeile
residualVarCovar[1,1] <- tbl_residual%>%filter(traits == "ccc")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[1,2] <- tbl_residual%>%filter(traits == "ccc+cca")%>%magrittr::extract2("meanEstimate")
#animalVarCovar[1,3] <- tbl_animal%>%filter(traits == "ccc+cfc")%>%magrittr::extract2("meanEstimate") #leer
residualVarCovar[1,3] <- 0.15656 #Value of Sample_1 CacCccCfcCwc
#residualVarCovar[1,4] <- tbl_residual%>%filter(traits == "ccc+cfa")%>%magrittr::extract2("meanEstimate")
residualVarCovar[1,5] <- tbl_residual%>%filter(traits == "ccc+cwc")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[1,6] <- tbl_residual%>%filter(traits == "ccc+cwa")%>%magrittr::extract2("meanEstimate")
residualVarCovar[1,7] <- tbl_residual%>%filter(traits == "ccc+cac")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[1,8] <- tbl_residual%>%filter(traits == "ccc+caa")%>%magrittr::extract2("meanEstimate")
#2. Zeile
#residualVarCovar[2,1] <- tbl_residual%>%filter(traits == "ccc+cca")%>%magrittr::extract2("meanEstimate")
residualVarCovar[2,2] <- tbl_residual%>%filter(traits == "cca")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[2,3] <- tbl_residual%>%filter(traits == "cca+cfc")%>%magrittr::extract2("meanEstimate")
residualVarCovar[2,4] <- tbl_residual%>%filter(traits == "cca+cfa")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[2,5] <- tbl_residual%>%filter(traits == "cca+cwc")%>%magrittr::extract2("meanEstimate")
residualVarCovar[2,6] <- tbl_residual%>%filter(traits == "cca+cwa")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[2,7] <- tbl_residual%>%filter(traits == "cca+cac")%>%magrittr::extract2("meanEstimate")
residualVarCovar[2,8] <- tbl_residual%>%filter(traits == "cca+caa")%>%magrittr::extract2("meanEstimate")
#3. Zeile
#animalVarCovar[3,1] <- tbl_animal%>%filter(traits == "cfc+ccc")%>%magrittr::extract2("meanEstimate") #leer
residualVarCovar[3,1] <- 0.15656 #Value of Sample_1 CacCccCfcCwc
#residualVarCovar[3,2] <- tbl_residual%>%filter(traits == "cca+cfc")%>%magrittr::extract2("meanEstimate")
residualVarCovar[3,3] <- tbl_residual%>%filter(traits == "cfc")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[3,4] <- tbl_residual%>%filter(traits == "cfc+cfa")%>%magrittr::extract2("meanEstimate")
residualVarCovar[3,5] <- tbl_residual%>%filter(traits == "cfc+cwc")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[3,6] <- tbl_residual%>%filter(traits == "cfc+cwa")%>%magrittr::extract2("meanEstimate")
residualVarCovar[3,7] <- tbl_residual%>%filter(traits == "cfc+cac")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[3,8] <- tbl_residual%>%filter(traits == "cfc+caa")%>%magrittr::extract2("meanEstimate")
#4. Zeile
#residualVarCovar[4,1] <- tbl_residual%>%filter(traits == "ccc+cfa")%>%magrittr::extract2("meanEstimate")
residualVarCovar[4,2] <- tbl_residual%>%filter(traits == "cca+cfa")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[4,3] <- tbl_residual%>%filter(traits == "cfc+cfa")%>%magrittr::extract2("meanEstimate")
residualVarCovar[4,4] <- tbl_residual%>%filter(traits == "cfa")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[4,5] <- tbl_residual%>%filter(traits == "cfa+cwc")%>%magrittr::extract2("meanEstimate")
residualVarCovar[4,6] <- tbl_residual%>%filter(traits == "cfa+cwa")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[4,7] <- tbl_residual%>%filter(traits == "cfa+cac")%>%magrittr::extract2("meanEstimate")
residualVarCovar[4,8] <- tbl_residual%>%filter(traits == "cfa+caa")%>%magrittr::extract2("meanEstimate")
#5. Zeile
residualVarCovar[5,1] <- tbl_residual%>%filter(traits == "ccc+cwc")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[5,2] <- tbl_residual%>%filter(traits == "cca+cwc")%>%magrittr::extract2("meanEstimate")
residualVarCovar[5,3] <- tbl_residual%>%filter(traits == "cfc+cwc")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[5,4] <- tbl_residual%>%filter(traits == "cfa+cwc")%>%magrittr::extract2("meanEstimate")
residualVarCovar[5,5] <- tbl_residual%>%filter(traits == "cwc")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[5,6] <- tbl_residual%>%filter(traits == "cwc+cwa")%>%magrittr::extract2("meanEstimate")
residualVarCovar[5,7] <- tbl_residual%>%filter(traits == "cac+cwc")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[5,8] <- tbl_residual%>%filter(traits == "caa+cwc")%>%magrittr::extract2("meanEstimate")
#6. Zeile
#residualVarCovar[6,1] <- tbl_residual%>%filter(traits == "ccc+cwa")%>%magrittr::extract2("meanEstimate")
residualVarCovar[6,2] <- tbl_residual%>%filter(traits == "cca+cwa")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[6,3] <- tbl_residual%>%filter(traits == "cfc+cwa")%>%magrittr::extract2("meanEstimate")
residualVarCovar[6,4] <- tbl_residual%>%filter(traits == "cfa+cwa")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[6,5] <- tbl_residual%>%filter(traits == "cwc+cwa")%>%magrittr::extract2("meanEstimate")
residualVarCovar[6,6] <- tbl_residual%>%filter(traits == "cwa")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[6,7] <- tbl_residual%>%filter(traits == "cac+cwa")%>%magrittr::extract2("meanEstimate")
residualVarCovar[6,8] <- tbl_residual%>%filter(traits == "caa+cwa")%>%magrittr::extract2("meanEstimate")
#7. Zeile
residualVarCovar[7,1] <- tbl_residual%>%filter(traits == "ccc+cac")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[7,2] <- tbl_residual%>%filter(traits == "cca+cac")%>%magrittr::extract2("meanEstimate")
residualVarCovar[7,3] <- tbl_residual%>%filter(traits == "cfc+cac")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[7,4] <- tbl_residual%>%filter(traits == "cfa+cac")%>%magrittr::extract2("meanEstimate")
residualVarCovar[7,5] <- tbl_residual%>%filter(traits == "cac+cwc")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[7,6] <- tbl_residual%>%filter(traits == "cac+cwa")%>%magrittr::extract2("meanEstimate")
residualVarCovar[7,7] <- tbl_residual%>%filter(traits == "cac")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[7,8] <- tbl_residual%>%filter(traits == "cac+caa")%>%magrittr::extract2("meanEstimate")
#8. Zeile
#residualVarCovar[8,1] <- tbl_residual%>%filter(traits == "ccc+caa")%>%magrittr::extract2("meanEstimate")
residualVarCovar[8,2] <- tbl_residual%>%filter(traits == "cca+caa")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[8,3] <- tbl_residual%>%filter(traits == "cfc+caa")%>%magrittr::extract2("meanEstimate")
residualVarCovar[8,4] <- tbl_residual%>%filter(traits == "cfa+caa")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[8,5] <- tbl_residual%>%filter(traits == "caa+cwc")%>%magrittr::extract2("meanEstimate")
residualVarCovar[8,6] <- tbl_residual%>%filter(traits == "caa+cwa")%>%magrittr::extract2("meanEstimate")
#residualVarCovar[8,7] <- tbl_residual%>%filter(traits == "cac+caa")%>%magrittr::extract2("meanEstimate")
residualVarCovar[8,8] <- tbl_residual%>%filter(traits == "caa")%>%magrittr::extract2("meanEstimate")

residualVarCovar
# not positive definite
matrixcalc::is.positive.definite(residualVarCovar)
eigen(residualVarCovar)
