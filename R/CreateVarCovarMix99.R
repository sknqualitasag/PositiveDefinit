#' ---
#' title: Build Parameter with Variances and Covariances for Mix99
#' date:  "`r Sys.Date()`"
#' ---

#' @title Build Parameter with Variances and Covariances for Mix99
#'
#' @export create_parameter_varCovar_mix99
create_parameter_varCovar_mix99 <- function(psInputFile,
                                            psOutputFile,
                                            pbLog = FALSE){

  # Prepare the different input to build the parameter file
  n_nr_randomEffect <- length(vec_randomEffect_name) #### ! Sophie : vec_randomEffect_name wie macht man damit es von Funktion hinweg gebraucht werden kann?
  # Check if in inputFile the random effects animal and residual are present
  vec_random_effect_req <- c("animal", "residual")
  if (!all(vec_random_effect_req %in% vec_randomEffect_name))
    stop(" * ERROR: Required random effects animal and residual are not both in list of random effects")
  # Get the other random effects #### ! Sophie : if not other? Implement this check Sophie!
  vec_random_effects_mand <- setdiff(vec_randomEffect_name, vec_random_effect_req)
  # Animal and residual should have a specific order in mix99
  vec_random_effect_order <- c(vec_random_effects_mand, vec_random_effect_req)

  # Check if psOutputFile is existing
  if (file.exists(psOutputFile))
    file.remove(psOutputFile)

  # Build Variance/Covariance Parameter-File for Mix99
  idx_rand_eff <- 1
  for(Z in vec_random_effect_order){ #### ! Sophie : vec_randomEffect_name and Z wie macht man damit es von Funktion hinweg gebraucht werden kann?
    for(i in 1:n_nr_trait){ #### ! Sophie : n_nr_trait wie macht man damit es von Funktion hinweg gebraucht werden kann?
      for(j in i:n_nr_trait){
        cat(idx_rand_eff, i, j, format(psInputFile[[Z]][[i,j]], scientific = FALSE), "\n", file = psOutputFile, append = TRUE)
      }
    }
    idx_rand_eff <- idx_rand_eff + 1
  }

}
