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
tbl_vce <- readr::read_delim(file = s_vce_result, delim = ";")
tbl_vce$estimate[tbl_vce$estimate == "---"] <- "0"
tbl_vce$estimate <- as.numeric(as.character(tbl_vce$estimate))

### # 'Matrices: NATRUAL' are the inputs for mix99 in the routine evaluation containing variance and covariance
tbl_varCovar <- tbl_vce %>% filter(type == "variance" | type == "covariance") %>% select(type,traits,random_effect,estimate)

### # ###############################################################
### # Build Matrice with for example Ccc Cca Cfc Cfa Cwc Cwa Cac Caa
### # ###############################################################

### # For each traits&random_effect get the mean value of all variants
by_grp <- tbl_varCovar %>% group_by(traits, random_effect)
smry <- summarise(by_grp,
                  meanEstimate = mean(estimate, na.rm = TRUE))

### # Split the column `traits` in order to facilitate the matrice building
smry <- smry %>% separate(traits, c("trait","surrogate"))

### # For variance, NA's are removed
for(i in 1:dim(smry)[1]){
     if(is.na(smry$surrogate[i])){
       smry$surrogate[i] <- smry$trait[i]
     }
}

### # Get random effects
randomEffectsName <- unique(smry %>% select(random_effect) %>% unique() %>% use_series(random_effect))

#### # Get trait
TraitName <- unique(smry %>% select(trait) %>% unique() %>% use_series(trait))

##for(Z in randomEffectsName){
  Z <- "animal"
  # Tibble per random effect
  tbl_Z <- smry %>% filter(random_effect == Z)

  # Check that all Traitname are as trait and as surrogate
  for(g in 1:length(TraitName)){
    g <- 1
    trtg <- TraitName[g]
    for(h in 1:length(TraitName)){
      h <- 2
      trth <- TraitName[h]

      #Combinaison trtg & trth is missing
      if(dim(tbl_Z %>% filter(trait == trtg & surrogate == trth))[1] == 0){
        #Combinaison trth & trtg is not missing, so get the info of the tibble for the combinaison trtg & trth
        if(dim(tbl_Z %>% filter(trait == trth & surrogate == trtg))[1] != 0){
          missing_re <- tbl_Z %>% filter(trait == trth & surrogate == trtg) %>% magrittr::extract2("random_effect")
          missing_me <- tbl_Z %>% filter(trait == trth & surrogate == trtg) %>% magrittr::extract2("meanEstimate")
          tbl_missing_trtg_trth <- tibble(trtg,trth,missing_re,missing_me)
          tbl_Z <- bind_rows(tbl_Z, tbl_missing_trtg_trth)


        }

      }
    }
  }





#  # Matrix with NULL values
#  matrix_VarCovar <- matrix(0, nrow = length(TraitName), ncol = length(TraitName))
#
##  for(i in 1:length(TraitName)){
#    i <- 1
##    for(j in 1:length(TraitName)){
#    j <- 1
#    #Variance and Covariance
#    matrix_VarCovar[i,j] <- tbl_Z %>% filter(trait == TraitName[[i]] && surrogate == TraitName[[j]]) %>% magrittr::extract2("meanEstimate") # NA's für Variance, nicht alle Variante in `trait` und ``surrogate!!!
#    j <- j + 1
##    }
#    i <- i + 1
##  }
#
##}






  ##  for(i in 1:nrow(VarCovar)){
  ##    #i <- 1
  ##    for(j in 1:ncol(VarCovar)){
  ##      #j <- 1
  ##      # ! error mit ccc+cfc deswegen [1,3] leer
  ##      VarCovar[i,j] <- tbl_Z %>% filter(traits == traits4Matrix[i,j]) %>% magrittr::extract2("meanEstimate")
  ##      j <- j + 1
  ##    }
  ##    i <- i + 1
  ##  }

## alle combination? wie, damit z.B. Madeleine mit ccv verwenden kann -> automatisch
#traits4Matrix <- rbind(c("ccc","ccc+cca","ccc+cfc","ccc+cfa","ccc+cwc","ccc+cwa","ccc+cac","ccc+caa"),
#                       c("ccc+cca","cca","cca+cfc","cca+cfa","cca+cwc","cca+cwa","cca+cac","cca+caa"),
#                       c("cfc+ccc","cca+cfc","cfc","cfc+cfa","cfc+cwc","cfc+cwa","cfc+cac","cfc+caa"),
#                       c("ccc+cfa","cca+cfa","cfc+cfa","cfa","cfa+cwc","cfa+cwa","cfa+cac","cfa+caa"),
#                       c("ccc+cwc","cca+cwc","cfc+cwc","cfa+cwc","cwc","cwc+cwa","cac+cwc","caa+cwc"),
#                       c("ccc+cwa","cca+cwa","cfc+cwa","cfa+cwa","cwc+cwa","cwa","cac+cwa","caa+cwa"),
#                       c("ccc+cac","cca+cac","cfc+cac","cfa+cac","cac+cwc","cac+cwa","cac","cac+caa"),
#                       c("ccc+caa","cca+caa","cfc+caa","cfa+caa","caa+cwc","caa+cwa","cac+caa","caa"))
