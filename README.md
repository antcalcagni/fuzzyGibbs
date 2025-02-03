# fuzzyGibbs

This repository contains the complete scripts to reproduce the Sections "Simulation Study" and "Applications" of the research article https://arxiv.org/abs/2501.18974

The folder tree is as follow:
- cache:        it contains temporary files (e.g., compiled c++ algorithms) used during the estimation procedure via the Approximate Gibbs Sampler (AGS)
- data:         it contains the observed fuzzy data for the proposed applications
- results:      it contains the mcmc samples (dataout_[applicationName_modelTested_chainId].rds) for each model being tested in the applications;
                it also contains the aggregated chains along with some posterior statistics (mcmc_[applicationName_modelTested].rds)
- utilities:    it contains algorithms and routines used internally by AGS;
                Note: AGS_utilities.R contains the workhorse AGS algorithm together with subroutines internally used by AGS
                    R2Cpp_funConverter_exec.R is the script used to convert R functions into executable cpp routines (the output is: R2Cpp_internal.cpp)
                    The compiled cpp files are stored into the folder cache/[modelName]/sourceCpp..

applications.R: It is the main R file used to define and fit the hypothesised models on each dataset. The file is structured into large sections, 
                each of them containing the 4 steps used to prepare the data, compile the model, run the AGS algorithm, compute posterior summary stats.

applications_analyses.R: It contains the script used for the PPC analysis on the three datasets used in applications.R
simulation_study.R: It contains the script used for running the simulation study along with the associated graphical results.

AGS_exec.R:  It is the R script internally used by applications.R to distribute AGS across cores via GNU Parallel
