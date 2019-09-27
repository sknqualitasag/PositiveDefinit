###
###
###
###   Purpose:   Run Bending and Production of Mix99 Paramater-File
###   started:   2019/09/27 (skn)
###
### ######################################### ###

### # argument parsing from commandline
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=2) stop("didn't receive 2 arguments")
s_vce_result <- args[1]
s_result_file <- args[2]
if(!file.exists(s_vce_result)) stop("first argument isn't an existing file")

### # Check that files exist
if (!file.exists(s_vce_result))
  stop("Cannot find input file: ",s_vce_result)

### # Run function read_vce
ResultTibble <- PositiveDefinit::read_vce(psInputFile = s_vce_result)

### # Run function build_matrix
ResultMatrixAsList <- PositiveDefinit::build_matrix(psInputFile = ResultTibble)


### # Check or Transfrom Matrix if necessary to insure beeing Positive Definit
ResultPD <- PositiveDefinit::check_transform_positivedefinit(psInputFile = ResultMatrixAsList,
                                                             psOptionRatio = FALSE,
                                                             psRatio = 100)


### # Build Parameter-File in txt-Format with Variances for Mix99
PositiveDefinit::create_parameter_varCovar_mix99(psInputFile = ResultPD,
                                                 psOutputFile = s_result_file)
