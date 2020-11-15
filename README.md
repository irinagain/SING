# SING

SING method implementation and code to reproduce results in corresponding paper

Risk, B. and Gaynanova, I. *"Simultaneous Non-Gaussian Component Analysis (SING) for Data Integration in Neuroimaging"* [arXiv link](https://github.com/irinagain/SING)

## Main functions

**jngcaFunctions.R** - all the R functions needed for method's implementation and corresponding analyses

**mCCAjointICA.R** - supporting R functions implementing Joint ICA, mCCA + jICA, and SING-averaged variant via orthogonal Procrustes with scaling

**CfunctionsCurvilinear.cpp** - additional C++ functions for implementation of curvilinear algorithm used for SING optimization; equivalent R functions can be found in **jngcaFunctions.R**

**makecifti.R** - R function to make cifti files, modified version from [Mandy Mejia](https://mandymejia.com), see her original version with instructions [here](https://mandymejia.com/2016/10/28/r-function-to-write-cifti-files/)

## Supporting data

**ComponentsForSim_setting2.Rda** - Shared and individual components used to simulate large-scale data for Simulation Setting 2. Also used as input to generate sparse components for Simulation Setting 3.


**community_affiliation_mmplus.csv** - File to support visualization of resting state correlations, our version is modified from [this paper](https://www.nature.com/articles/s41598-019-55738-y). A similar version can be also found [here](https://github.com/emergelab/hierarchical-brain-networks/blob/master/brainmaps/node_affiliations/AA_all_maps.csv).

## Simulations scripts

### Simulation Setting 1

**s1_rank_simulations_small.R** - Simulation setting 1, rank estimation. (WARNING: this should only be run on a cluster as it takes significant time due to a total of 400 evaluations, and 1000 permutations for rank selection)

**s2_estimation_simulations_small.R** - Simulation setting 1, component and mixing matrices estimation performance. Can be run in parallel with s1.

**s3_results_analysis_small.R** - Simulation setting 1, uses the saved results from s1 and s2 to generate figures in the paper

### Simulation Setting 2

**s4_generate_data_sim_setting2.R** - Simulation setting 2, uses the components in supplied .Rda file as true components to generate simulated data.

**s5_estimation_sim_setting2.R** - Simulation setting 2, apply joint ICA, mCCA + jICA and SING variants, also apply rank estimation approach.

**s6_results_sim_setting2.R** - Simulation setting 2, use the results from s5 to calculate component and mixing matrix estimation errors, create figures in the paper.

### Simulation Setting 3

**s7_generate_data_sim_setting3.R** - Simulation setting 3, uses the components in supplied .Rda file with additional truncation as true **sparse** components to generate simulated data.

**s8_estimation_sim_setting3.R** - Simulation setting 3, apply joint ICA, mCCA + jICA and SING variants, also apply rank estimation approach.

**s9_results_sim_setting3.R** - Simulation setting 3, use the results from s8 to calculate component and mixing matrix estimation errors, create figures in the paper.

