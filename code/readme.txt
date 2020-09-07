THIS SOURCE CODE IS SUPPLIED “AS IS” WITHOUT WARRANTY OF ANY KIND, AND ITS AUTHOR AND THE JOURNAL OF
MACHINE LEARNING RESEARCH (JMLR) AND JMLR’S PUBLISHERS AND DISTRIBUTORS, DISCLAIM ANY AND ALL WARRANTIES,
INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, AND ANY WARRANTIES OR NON INFRINGEMENT. THE USER
ASSUMES ALL LIABILITY AND RESPONSIBILITY FOR USE OF THIS
SOURCE CODE, AND NEITHER THE AUTHOR NOR JMLR, NOR
JMLR’S PUBLISHERS AND DISTRIBUTORS, WILL BE LIABLE FOR
DAMAGES OF ANY KIND RESULTING FROM ITS USE. Without limiting the generality of the foregoing, neither the author, nor JMLR, nor
JMLR’s publishers and distributors, warrant that the Source Code will be
error-free, will operate without interruption, or will meet the needs of the
user.

Readme file as Supplementary Material to the paper 
Nicholson, Wilms, Bien, Matteson (2020), High Dimensional Forecasting via Interpretable Vector Autoregression, Journal of Machine Learning Research 21, 1-52.

Below, we describe how to reproduce all results from the paper using the code provided in the zip-folder "code"

####################
# Simulation Study #
####################

# Section Forecast Comparisons &  Robustness on HLag as pmax increases #
- Go to code/simulation/forecast and set this as your working directory in R
- Open the R-script "forecastscript.R"
- This R-script contains all the code to reproduce the forecast results from Simulation Scenarios 1-5 for all estimators except for the method of Giannone et al. (2015), the factor models, AR and VAR(1)
- Open the R-script "forecast_other_benchmarks.R"
- This R-script contains all the code to reproduce the forecast results from Simulation Scenarios 1-4 for the factor models, AR and VAR(1)


# Section Lag Order Selection #
- Go to Code/Simulation/Lagselection and set this as your working directory in R
- Open the R-script "lagselectionscript.R"
- This R-script contains all the code to reproduce the lag order selection results for all estimators 


#################
# Applications  #
#################

##############
# macro data #
##############
- Go to code/applications/macro 
- Open the R-script "smallmedium.R". This R-script contains all the code to reproduce the results from the smallmedium VAR for all estimators except for the method of Giannone et al. (2015) and Carriero et al. (2017)
- Open the R-script "medium.R". This R-script contains all the code to reproduce the results from the medium VAR for all estimators except for the method of Giannone et al. (2015) and Carriero et al. (2017)
- Open the R-script "mediumlarge.R". This R-script contains all the code to reproduce the results from the mediumlarge VAR for all estimators except for the method of Giannone et al. (2015) and Carriero et al. (2017)
- Open the R-script "large.R". This R-script contains all the code to reproduce the results from the large VAR for all estimators except for the method of Giannone et al. (2015) and Carriero et al. (2017)

################
# finance data #
################
- Go to code/applications/finance
- Open the R-script "finance.R". This R-script contains all the code to reproduce the results from the finance application

################
# energy data #
################
- Go to code/applications/energy
- Open the R-script "energy.R". This R-script contains all the code to reproduce the results from the energy application


####################
# Bayesian methods #
####################
Note that the github page of the first author contains in addition the following files, to replicate the results of the Bayesian methods:
# Simulation Study #
- Go to Code/Simulation/Forecast/GLPmatlab
- Open the matlab-script "SimX.m" with X=1,2,3,4,5
- Each script contains all code to reproduce the forecast results for the method of Giannone et al. (2015) in  Simulation Scenario X. 
  We use the matlab code as provided by the authors to obtain the GLP BVAR estimates

- Go to Code/Simulation/Forecast/CCMmatlab
- Open the matlab-script "simX_CCM.m" with X=1,2,3,4,5
- Each script contains all code to reproduce the forecast results for the method of Carriero et al. (2017) in  Simulation Scenario X. 
  We combine the matlab code as provided by the authors with the Matlab code for BVAR using Gibbs sampling available at  https://sites.google.com/site/dimitriskorobilis/matlab/code-for-vars to obtain the CCM BVAR estimates

# Applications  #
- Go to Code/Application/SW/GLPmatlab
- Open the matlab-script "smallmedium.m" This script contains all code to reproduce the results for the method of Giannone et al. (2015) on the smallmedium VAR data.
- Open the matlab-script "medium.m" This script contains all code to reproduce the results for the method of Giannone et al. (2015) on the medium VAR data.
- Open the matlab-script "mediumlarge.m" This script contains all code to reproduce the results for the method of Giannone et al. (2015) on the mediumlarge VAR data.
- Open the matlab-script "large.m" This script contains all code to reproduce the results for the method of Giannone et al. (2015) on the large VAR data.
  We use the matlab code as provided by the authors to obtain the GLP BVAR estimates

- Go to Code/Application/SW/CCMmatlab
- Open the matlab-script "smallmedium.m" This script contains all code to reproduce the results for the method of Carriero et al. (2017) on the smallmedium VAR data.
- Open the matlab-script "medium.m" This script contains all code to reproduce the results for the method of Carriero et al. (2017) on the medium VAR data.
- Open the matlab-script "mediumlarge.m" This script contains all code to reproduce the results for the method of Carriero et al. (2017) on the mediumlarge VAR data.
- Open the matlab-script "large.m" This script contains all code to reproduce the results for the method of Carriero et al. (2017) on the large VAR data.
  We combine the matlab code as provided by the authors with the Matlab code for BVAR using Gibbs sampling available at  https://sites.google.com/site/dimitriskorobilis/matlab/code-for-vars to obtain the CCM BVAR estimates
