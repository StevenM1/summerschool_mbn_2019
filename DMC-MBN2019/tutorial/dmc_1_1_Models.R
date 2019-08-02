##################  DMC Lesson 1 Models and Data
#
# Note: Before running these scripts, the current working directory must be set 
# to the top-level folder containing the dmc and tutorial subfolders 

### Lesson 1.1 Models

rm(list=ls()) 
source ("dmc/dmc.R")

# For this lesson we're using the LBA model. Only one model can be active at
# a time. Models are specified in the "models" directory under the "dmc" 
# directory. It also contains further sub-directories, the names of which are 
# the first argument to the load_model function. These subdirectories contain 
# a file "dists.R" which defines the random and likelihood functions that 
# define the model in its core parameterization. The second argument to
# load_model specifies a particular model parameterization. For the following 
# LBA model the core parameterization is in terms of the threshold "b" parameter
# but we will use a re-parameterization B=b-A (where A is start point 
# variability) defined in "lba_B.R".  
load_model (dir_name="LBA",model_name="lba_B.R")

# NB1: See ?LBA for rtdists package help on the LBA parameters.
# NB2: More advanced users can create their own reparameterizations using an
#      existing model file as a template to alter the three functions it contains
#      transform.dmc, random.dmc and likelihood.dmc (See the lessons in 
#      dmc_2_3_AddModel1.R and dmc_2_4_AddModel2.R, but it is recommended you
#      work through to the end of lesson 4 first.

# SETTING UP A MODEL requires 6 steps, in three groups. 
# 
# Steps 1-5 create a set of list and vector objects, and can be done in any 
# order. Step 6 uses these objects as arguments to the "model.dmc" function, 
# which creates an array object ("model") with a set of attributes specifying a
# particular model and parameterization.

#### Example 1: a single subject with parameters varying with response (R) 
####   and correctness (M, for match) in a minimal design with one factor (S, 
####   a two level stimulus factor) that is used to score two possible 
####   responses named "r1" and "r2".
#
#  The first three steps specify the factors, possible responses, and how they 
#  are scored for accuracy (this mapping creates a new factor "M", with levels
#  "true" and "false"). 
#
#  The next two steps specify details about the model parameters and how they 
#  relate to factors (including M).
# 
#  1) Create a named list ( e.g., "factors" ) specifying factor names 
#     ( = list names ) and factor levels ( = character vectors assigned to 
#     corresponding list elements ).
#
factors <- list(S=c("s1","s2"))  # Stimulus 1 and Stimulus 2
#
#     NB: Factor levels cannot contain a "." as this is used by DMC to 
#         concantanate factor levels from multiple factors.

#
#  2) Create a character vector ( e.g., "responses" ) giving the names of the
#     possible responses ( usually levels (data$R) in the data frame )
#
responses <- c("r1","r2")        # Response 1 and Response 2

#
#  3) Create "match.map", a list of lists, where the first entry is a list 
#     named "M" that specifies which accumulator (or boundary, for the DDM
#     model) matches the correct response for each factor level. 
#     Accumulators/boundaries can be specified by the corresponding response 
#     name or an integer 1..n for n response possibilities (see dmc_1_6 and 
#     dmc_1_7 for more advanced scoring when there are more than two responses).
#    
match.map <- list(M=list(s1="r1",s2="r2"))
#
#     NB: Typically the levels all come from one factor, but can be a 
#         combination of levels from different factors. In this case the
#         levels must be specified in the same order as factor names 
#         (see 1 above and dmc_5_1 and dmc_5_2) 

####  For this example we will setup an LBA model where the B parameter varies 
####  with the response (R), mean_v and sd_v vary with match (M), and there are
####  constants for st0 and for the sd_v of the mismatching accumulator. 
#
#  4) Create "p.map", a list which tells DMC how the factors specified above map
#     onto the parameters of the model.
#     These "EXTERNAL PARAMETERS TYPES" are defined in a model parameterization 
#     file ("lba_B.R" in this example), and "FACTORS" are specified in 
#     steps (1) & (3) (usually M for the latter). 
#
p.map <- list(A="1",B="R",t0="1",mean_v="M",sd_v="M",st0="1")
#
#     NB1: The R and M factors MUST be given last when used with a parameter
#          influenced by more than one factor
#     NB2: The "1" indicates that the same estimated value will be assumed for
#          all cells in the design (i.e., an "intercept" parameter).
#         

#
#  5) Create "const", a vector which specifies the parameters to  
#     be held constant (i.e., not varied in fitting), and their values.
#
constants <- c(sd_v.false=1,st0=0)
#
#     NB: DMC names parameters by the parameter-type name and 
#         (possibly) concatenated factor levels to which they apply, 
#           e.g., "parameter_type_name.level_factor_A.level_factor_B ..."
#         Again these levels must be in the same order as names in step (1)
#         Once a model is created you can view these names (assuming you
#         called the model "model") using : attr(model,"p.vector")
#

####  Now we create the model object, specifying the "norm" distribution 
####  function (one of the options for the LBA model in the "rtdists" package 
####  "n1PDF" function, passed through to the "likelihood.dmc" function in the
####  "lba_B.R" script.)

#
#  6) Use the objects created in steps 1 to 5 as arguments to the function 
#     "model.dmc" in order to create a Boolean array (with a set of attributes, 
#      e.g., "model") that specifies the model to be fit.
#
model <- model.dmc(p.map,match.map=match.map,constants=constants, 
                   factors=factors,responses=responses,type="norm")
#
#    NB1: The following things are checked and will cause a STOP
#         a) Factors names "M", "R", and "1" (intecept factor) are reserved
#         b) Factor levels must be unique within AND BETWEEN factors
#         c) List "match.map" must have values in responses  
#         d) Where match.map contains more than one list their names cannot
#            be the same as factor names from (1)

# This produces the following output:

# Parameter vector names are: ( see attr(,"p.vector") )
# [1] "A"            "B.r1"         "B.r2"         "t0"           "mean_v.true" 
# [6] "mean_v.false" "sd_v.true"   
# 
# Constants are (see attr(,"constants") ):
# sd_v.false        st0 
#          1          0 
# 
# Model type = norm (posdrift= TRUE ) 


# This information is also available as attributes of the model object, e.g.,
attr(model,"p.vector")
attr(model,"constants")
attr(model,"type")
attr(model,"posdrift")   

# NB: The "posdrift" boolean is a modifier specific to the norm type of the LBA
#     model (see "rtdists") 

# THE MODEL ARRAY is
model

# and produces the following output:
#
# , , r1
# 
#          A B.r1  B.r2   t0 mean_v.true mean_v.false sd_v.true sd_v.false  st0
# s1.r1 TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# s2.r1 TRUE TRUE FALSE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# s1.r2 TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# s2.r2 TRUE TRUE FALSE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# 
# , , r2
# 
#          A  B.r1 B.r2   t0 mean_v.true mean_v.false sd_v.true sd_v.false  st0
# s1.r1 TRUE FALSE TRUE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# s2.r1 TRUE FALSE TRUE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# s1.r2 TRUE FALSE TRUE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# s2.r2 TRUE FALSE TRUE TRUE        TRUE        FALSE      TRUE      FALSE TRUE

# There is one matrix for each accumulator, with row names showing the 
# level.response naming and the columns the parameter-type.level naming for
# parameters. The booleans show which parameters map to which design+response
# cell.

# "model" is used by the function "p.df.dmc" (in dmc_model.R) to perform 
# parameter assignment for each row (i) of the model array, and transforms 
# and creates parameters as specified in "transform.dmc" (in "lba_B.R") 

# To see this, first create a parameter vector to map. These are not very 
# sensible values, chosen just so we can see differences in the following 
# printout
p.vector  <- c(A=3,B.r1=4,B.r2=5,
               mean_v.true=2,mean_v.false=-1,
               sd_v.true=0.5,
               t0=.2)

# The "model" array above is computationally efficient but it is not easy to
# see how parameters map to design cells. The function "print.cell.p" makes it
# easy to loop through cells and print the mapping.
print.cell.p(p.vector,model)

# [1] "s1.r1"
#    A b  t0 mean_v sd_v st0
# r1 3 7 0.2      2  0.5   0
# r2 3 8 0.2     -1  1.0   0
# 
# [1] "s2.r1"
#    A b  t0 mean_v sd_v st0
# r1 3 7 0.2     -1  1.0   0
# r2 3 8 0.2      2  0.5   0
# 
# [1] "s1.r2"
#    A b  t0 mean_v sd_v st0
# r2 3 8 0.2     -1  1.0   0
# r1 3 7 0.2      2  0.5   0
# 
# [1] "s2.r2"
#    A b  t0 mean_v sd_v st0
# r2 3 8 0.2      2  0.5   0
# r1 3 7 0.2     -1  1.0   0

# NB1: the rows are in "n1" order, that is, for ".r1" cells the first
#      row is "r1", for ".r2" cells the first row is "r2". This is because 
#      "n1PDF" (in the "rtdists" package) returns the likelihood for first unit 
#      or node in a set of accumulators (the term n1 mean "node one") with 
#      parameters corresponding to the first entries in each parameter type's
#      parameter vector.

# NB2: "transform.dmc" creates b = A+B, so the columns are the "internal" 
#      parameter types for the LBA (i.e., the name expected by "n1PDF" given 
#      the type designator passed to likelihood.dmc). "transform.dmc" can also
#      be used to change the scale on which parameters are sampled. For example, 
#      you might sample sd_v on a log scale to enforce positivity, so 
#      "transform.dmc" would have to calculate 
#      par.df$sd_v <- exp(par.df$log_sd_v) 
# save_data(model,p.vector,file="dmc_1_1.Rdata")

# Most of the DMC lessons use very simple designs. If you want to know how to
# setup more complex designs see the advanced lesson dmc_5_1_Factorial.R. If you 
# want to see how to setup more complex scoring see dmc_5_2_Scoring.R. 

