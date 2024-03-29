###
###
###
###   Purpose:   Build Parameter with Variances and Covariances for Mix99
###   started:   2019/09/30 (skn)
###
### ######################################### ###


#' @title Build Parameter with Variances and Covariances for Mix99
#'
#' @export create_parameter_varCovar_mix99
create_parameter_varCovar_mix99 <- function(psInputFile,
                                            psOutputFile,
                                            pbLog = FALSE){

  # Prepare the different input to build the parameter file
  vec_randomEffect_name <- names(psInputFile)
  n_nr_randomEffect <- length(vec_randomEffect_name)
  # Check if in inputFile the random effects animal and residual are present
  vec_random_effect_req <- c("animal", "residual")
  if (!all(vec_random_effect_req %in% vec_randomEffect_name))
    stop(" * ERROR: Required random effects animal and residual are not both in list of random effects")
  # Get the other random effects
  vec_random_effects_mand <- setdiff(vec_randomEffect_name, vec_random_effect_req)
  # Animal and residual should have a specific order in mix99
  vec_random_effect_order <- c(vec_random_effects_mand, vec_random_effect_req)

  # Check if psOutputFile is existing
  if (file.exists(psOutputFile))
    file.remove(psOutputFile)

  # Build Variance/Covariance Parameter-File for Mix99
  n_nr_trait <- dim(psInputFile[[1]])[1]
  idx_rand_eff <- 1
  for(Z in vec_random_effect_order){
    for(i in 1:n_nr_trait){
      for(j in i:n_nr_trait){
        cat(idx_rand_eff, i, j, format(psInputFile[[Z]][[i,j]], scientific = FALSE), "\n", file = psOutputFile, append = TRUE)
      }
    }
    idx_rand_eff <- idx_rand_eff + 1
  }

}



#' @title Parameter File for Mix99
#'
#' @description Read csv-file of the variance component estimates and construct
#' a matrix with this estimates by parts. This matrix is checked if she is
#' positive definite or not. If this matrix is not positive definite, 2
#' functions may be called. If you put the paramter psOptionRatio = FALSE,
#' makePD2() is called . For the parameter psOptionRatio = TRUE the function
#' make_pd_rat_ev() is called. The difference is, that with make_pd_rat_ev() you
#' have the option to whish a maximum ratio between largest and smallest
#' eigenvalue.
#'
#' The software mix99 need a paramter-File with the variance and covariance. So
#' this function is building this paramter-file.
#'
#' @param psInputFile input csv-file
#' @param psOptionRatio TRUE or FALSE
#' @param psRatio number
#' @param psOutputFile output txt-file
#'
#' @examples
#' \dontrun{
#' PositiveDefinit::parameter_varCovar_mix99(psInputFile = '/qualstorzws01/data_projekte/projekte/singularity_data_zws_gslim/muku_CarcassVK/work/VCE_results.csv',
#' psOutputFile = 'par_varCovar_mix99.txt')
#' }
#' @export parameter_varCovar_mix99
parameter_varCovar_mix99 <- function(psInputFile = psInputFile,
                                     psOptionRatio = FALSE,
                                     psRatio = 100,
                                     psOutputFile = psOutputFile){

  ### # Check or Transfrom Matrix if necessary to insure beeing Positive Definit
  ResultPD <- PositiveDefinit::positivedefinit(psInputFile,
                                              psOptionRatio,
                                              psRatio)

#  ### # Build Parameter-File in txt-Format with Variances for Mix99
  PositiveDefinit::create_parameter_varCovar_mix99(psInputFile = ResultPD,
                                                   psOutputFile = psOutputFile)

}
