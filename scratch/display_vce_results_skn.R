setwd("/Volumes/data_projekte/projekte/singularity_data_zws_gslim/muku_CarcassVK/work/")

suppressPackageStartupMessages(if(! require("ggplot2")) {
  install.packages("ggplot2", repos="https://stat.ethz.ch/CRAN/")
  require("ggplot2")
})

d.vr <- read.table(file = "VCE_results.csv", sep =";", header = TRUE)

orderTraitNames <- function(traits) {
  traits <- as.character(traits)
  trtvec <- strsplit(traits, split = "+", fixed = TRUE)[[1]]
  paste(trtvec[order(trtvec)], collapse = "+")
}
d.vr$traits <- as.factor(apply(d.vr[,"traits",drop=FALSE], 1, orderTraitNames))
d.vr$trait_combination <- as.factor(apply(d.vr[,"trait_combination",drop=FALSE], 1, orderTraitNames))
d.vr$estimate <- as.numeric(as.character(d.vr$estimate))
d.vr$STD_ERR_estimate <- as.numeric(as.character(d.vr$STD_ERR_estimate))

### # Genetische Korrelation ----------------------
t.corr <- d.vr[d.vr$type=="genetic_correlation",]
gg <- ggplot(data = t.corr, aes(traits, estimate)) +
  geom_point(aes(colour = data_subset_number)) +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(rows= vars(model_name)) +
  ggtitle("Genetic correlation")
gg



### # Heritabilitaet ----------------------
t.h2 <- d.vr[d.vr$type=="heritability",]

gg <- ggplot(data = t.h2, aes(traits, estimate)) +
  geom_jitter(aes(colour = data_subset_number), width = 0.15) +
  facet_grid(rows= vars(model_name)) +
  ggtitle("Heritability")
gg

gg <- ggplot(data = t.h2, aes(traits, estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin=estimate-STD_ERR_estimate, ymax=estimate+STD_ERR_estimate), size=0.3, width=0.2)
gg


### # Varianz ----------------------
t.variance <- d.vr[d.vr$type=="variance",]

gg <- ggplot(data = t.variance, aes(traits, estimate)) +
  geom_point(aes(colour = factor(data_subset_number))) +
  facet_grid(rows= vars(model_name), cols = vars(random_effect)) +
  ggtitle("Variance")

gg
