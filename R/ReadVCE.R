#' ---
#' title: Read VCE results
#' date:  "`r Sys.Date()`"
#' ---

#' @title Reading VCE results
#'
#' @description
#' First attempt to read variance component estimation (VCE) from a csv-file
#' and storing the input in tibble
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tidyr
#' @importFrom readr read_delim
#' @export read_vce
read_vce <- function(psInputFile,
                     pbLog = FALSE){

  ### # Read all VCE results
  tbl_vce <- readr::read_delim(file = psInputFile, delim = ";")
  # Transform "---" to 0 coming from VCE software some times
  tbl_vce$estimate[tbl_vce$estimate == "---"] <- "0"
  # Transform estimates to numeric
  tbl_vce$estimate <- as.numeric(as.character(tbl_vce$estimate))

  ### # Resulting tibble
  return(tbl_vce)

}



#' @title Reading VCE results for graphics
#'
#' @export read_vce4grafics
read_vce4grafics <- function(psInputFile,
                     pbLog = FALSE){

  ### # Read all VCE results
  d.vr <- read.table(file = psInputFile, sep =";", header = TRUE)

  orderTraitNames <- function(traits) {
    traits <- as.character(traits)
    trtvec <- strsplit(traits, split = "+", fixed = TRUE)[[1]]
    paste(trtvec[order(trtvec)], collapse = "+")
  }
  d.vr$traits <- as.factor(apply(d.vr[,"traits",drop=FALSE], 1, orderTraitNames))
  d.vr$trait_combination <- as.factor(apply(d.vr[,"trait_combination",drop=FALSE], 1, orderTraitNames))
  d.vr$estimate <- as.numeric(as.character(d.vr$estimate))
  d.vr$STD_ERR_estimate <- as.numeric(as.character(d.vr$STD_ERR_estimate))

  ### # Resulting dataframe
  return(d.vr)

}
