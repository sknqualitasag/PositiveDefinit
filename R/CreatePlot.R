#' ---
#' title: Create Plot to Overview VCE results
#' date:  "`r Sys.Date()`"
#' ---

#' @title Create Plot to Overview VCE results for genetic correlation
#'
#' @export plot_gencorr
plot_gencorr <- function(psInputFile){

  t.corr <- psInputFile[psInputFile$type=="genetic_correlation",]
  gg <- ggplot(data = t.corr, aes(traits, estimate)) +
    geom_point(aes(colour = data_subset_number)) +
    theme(axis.text.x = element_text(angle = 90)) +
    facet_grid(rows= vars(model_name)) +
    ggtitle("Genetic correlation")

  return(print(gg))

}



#' @title Create Plot to Overview VCE results for heritability
#'
#' @export plot_h2
plot_h2 <- function(psInputFile){

  t.h2 <- psInputFile[psInputFile$type=="heritability",]
  gg <- ggplot(data = t.h2, aes(traits, estimate)) +
    geom_jitter(aes(colour = data_subset_number), width = 0.15) +
    facet_grid(rows= vars(model_name)) +
    ggtitle("Heritability")

  return(print(gg))

}



plot_var
#' @title Create Plot to Overview VCE results for variances
#'
#' @export plot_var
plot_var <- function(psInputFile){

  t.variance <- psInputFile[psInputFile$type=="variance",]

  gg <- ggplot(data = t.variance, aes(traits, estimate)) +
    geom_point(aes(colour = factor(data_subset_number))) +
    facet_grid(rows= vars(model_name), cols = vars(random_effect)) +
    ggtitle("Variance")

  return(print(gg))

}
