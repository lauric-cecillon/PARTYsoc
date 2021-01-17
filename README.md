# PARTYsoc
## A statistical model partitioning soil organic carbon (SOC) into its centennially active and stable fractions based on Rock-Eval thermal analysis

### USER MANUAL of the second version: PARTYsocv2.0 and PARTYsocv2.0eu
[![DOI](https://zenodo.org/badge/319097495.svg)](https://zenodo.org/badge/latestdoi/319097495)

#### This short manual and the R scripts, R data file and csv file provided in this repository accompany the draft: "Partitioning soil organic carbon into its centennially active and stable fractions with statistical models based on Rock-Eval® thermal analysis (PARTYsocv2.0 and PARTYsocv2.0eu)"
#### Submitted to Geoscientific Model Development by Lauric Cécillon et al.
#### Contact e-mail: lauric.cecillon@inrae.fr

#### The part 1 of the R script entitled "importing_rock-eval_data_computing_rock-eval_parameters_GMD" can be used to extract Rock-Eval data from raw Rock-Eval 6 files (R.00 and S.00).

#### The part 2 of the R script entitled "importing_rock-eval_data_computing_rock-eval_parameters_GMD" can be used to calculate Rock-Eval parameters that are not provided by the software of Vinci Technologies.

#### The R script entitled "partysocv2.0_partysocv2.0eu_Rscript_GMD" can be used to reproduce the models (PARTYsocv2.0 and PARTYsoc2.0eu statistical models) and the results presented in the above-mentionned draft, using the data (csv) file entitled "PARTYsocv2.0_GMD". This csv file contains the 40 Rock-Eval parameters of the 105 topsoil samples from seven reference sites described in the above-mentionned draft.

#### The R script entitled "PARTYsocv2.0eu_prediction_Rscript_GMD" can be used to run directly the PARTYsocv2.0eu statistical model (i.e. predicting the size of the centennially active and stable SOC fraction in new topsoils), using 18 Rock-Eval parameters (measured on these new topsoils) as predictor variables. This R scripts uses the R data file entitled "PARTYsocv2.0eu_data_for_prediction_GMD".
