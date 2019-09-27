###
###
###
###   Purpose:   Run Graphics as Overview for VCE
###   started:   2019/09/27 (skn)
###
### ######################################### ###

### # argument parsing from commandline
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=2) stop("didn't receive 2 arguments")
s_vce_result <- args[1]
if(!file.exists(s_vce_result)) stop("first argument isn't an existing file")

### # Check that files exist
if (!file.exists(s_vce_result))
  stop("Cannot find input file: ",s_vce_result)

### # Run function read_vce
ResultDF <- PositiveDefinit::read_vce4grafics(psInputFile = s_vce_result)

### # Plot genetic correlation
PositiveDefinit::plot_gencorr(psInputFile = ResultDF) ### Sophie: How make direct ppt oder pdf with plots

### # Plot heritability
PositiveDefinit::plot_h2(psInputFile = ResultDF) ### Sophie: How make direct ppt oder pdf with plots

### # Plot variance
PositiveDefinit::plot_var(psInputFile = ResultDF) ### Sophie: How make direct ppt oder pdf with plots
