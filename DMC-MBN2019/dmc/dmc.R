# System functions for the DMC (Dynamic Models of Choice)
#    Usually user does not need to edit
#
# This script loads in all required files for DMC
#

require(parallel) # parallel processing
require(msm)  # For truncated normal priors 
require(truncdist) # truncated t etc.
require(coda) # Sampling output analysis
require(loo) # For WAIC and looaic calculation
require(hypergeo) # For population plausible values
require(statmod) # Wald model
require(pracma)  # For gng and stop signal robust integration
require(numDeriv) # Prior transformations
require(vioplot) # Stop signal graphs
require(ggplot2) # For fancy graphs
require("mvtnorm") # For Bayes Factors
require("Matrix") # For Bayes Factors
require("Brobdingnag") # For Bayes Factors
require("stringr") # For Bayes Factors

# Pull in the file utilities so load_model can be called from dmc_hierarchical 
temp_wd <- getwd (); setwd(file.path (temp_wd, "dmc")) 

source ("file_utils.R")

# Load in all the DMC modules
source ("dmc_model.R")
source ("dmc_sampling.R")
source ("dmc_hierarchical.R")
source ("dmc_plotting.R")
source ("dmc_analysis.R")
source ("bridge_sampling_functions_warp3.R") # Written by Quentin F. Gronau

setwd(temp_wd)


