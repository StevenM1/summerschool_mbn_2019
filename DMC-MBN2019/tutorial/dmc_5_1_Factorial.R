##################  DMC Lesson 5: Advanced Models and Data
#
# Note: Before running these scripts, the current working directory must be set 
# to the top-level folder containing the dmc and tutorial subfolders 

### Lesson 5.1:  Factorial Designs

rm(list=ls()); par(mfrow = c(1,1)) 
source ("dmc/dmc.R")

load_model ("LBA","lba_B.R") # LBA model with B=b-A parameterization 

# EXAMPLE 1: 2 x 2 design, S x F

factors <- list(S=c("s1","s2"),F=c("f1","f2"))
responses <- c("r1","r2")
match.map <- list(M=list(s1="r1",s2="r2"))

# Use same simple model as for Lesson 1 (i.e., factor F present but not used)
p.map <- list(A="1",B="R",t0="1",mean_v="M",sd_v="M",st0="1")
const <- c(sd_v.false=1,st0=0)

model.dmc(p.map,match.map=match.map,constants=const,responses=responses,
          factors=factors)

# , , r1
# 
#             A B.r1  B.r2   t0 mean_v.true mean_v.false sd_v.true sd_v.false  st0
# s1.f1.r1 TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# s2.f1.r1 TRUE TRUE FALSE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# s1.f2.r1 TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# s2.f2.r1 TRUE TRUE FALSE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# s1.f1.r2 TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# s2.f1.r2 TRUE TRUE FALSE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# s1.f2.r2 TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# s2.f2.r2 TRUE TRUE FALSE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# ...


# EXAMPLE 2: 2x2x2 (Stimulus and factors F and G) specified in model.dmc call

model.dmc(p.map,match.map=match.map,constants=const,responses=responses,
          factors=list(S=c("s1","s2"),F=c("f1","f2"),G=c("g1","g2")))
# , , r1
# 
#             A    B.r1 B.r2  t0   mean_v.true mean_v.false sd_v.true  sd_v.false st0
# s1.f1.g1.r1 TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# s2.f1.g1.r1 TRUE TRUE FALSE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# s1.f2.g1.r1 TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# s2.f2.g1.r1 TRUE TRUE FALSE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# s1.f1.g2.r1 TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# s2.f1.g2.r1 TRUE TRUE FALSE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# s1.f2.g2.r1 TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# s2.f2.g2.r1 TRUE TRUE FALSE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# s1.f1.g1.r2 TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# s2.f1.g1.r2 TRUE TRUE FALSE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# s1.f2.g1.r2 TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# s2.f2.g1.r2 TRUE TRUE FALSE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# s1.f1.g2.r2 TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# s2.f1.g2.r2 TRUE TRUE FALSE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# s1.f2.g2.r2 TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# s2.f2.g2.r2 TRUE TRUE FALSE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# ...


# EXAMPLE 3: Lexical decision example where the single S factor has different 
# levels TO response. This demonstrates how to make a many-to-one mapping from 
# stimulus to response in match.map. Also shows we can use names (of response 
# levels) rather than numbers in match.map

model.dmc(p.map,constants=const,
          factors=list(S=c("lf","hf","nw")),
          responses=c("WW","NW"),
          match.map=list(M=list(lf="WW",hf="WW",nw="NW")))

# , , WW
# 
#         A. B.WW  B.NW   t0 mean_v.true mean_v.false sd_v.true sd_v.false  st0
# lf.WW TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# hf.WW TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# nw.WW TRUE TRUE FALSE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# lf.NW TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# hf.NW TRUE TRUE FALSE TRUE        TRUE        FALSE      TRUE      FALSE TRUE
# nw.NW TRUE TRUE FALSE TRUE       FALSE         TRUE     FALSE       TRUE TRUE
# ...


# EXAMPLE 4: In this example a factor in the experimental design is named factor 
# A, which we assume varies with the start-point parameter of the model (A).
# That is we now use an extra factor in p.map. The design has the 3-level 
# stimulus and factor A has two levels (i.e., 3x2), As parameter A differs as a 
# function of factor A there are now parameters A.a1 and A.a2, one for each 
# level of factor A (i.e., a1 and a2)

p.map <- list(A="A",B="R",t0="1",mean_v="M",sd_v="M",st0="1")

model <- model.dmc(p.map,constants=const,
                   factors=list(S=c("lf","hf","nw"),A=c("a1","a2")),
                   responses=c("WW","NW"),
                   match.map=list(M=list(lf=1,hf=1,nw=2)))

# Specifying a long p.vector can be error prone. "check.p.vector" checks if the
# parameter vector is compatible with a model, printing warnings if there are
# any problems.

p.vector  <- c(A.a1=3,A.a2=2,
               B.WW=4,B.NW=5,
               t0=.2,mean_v.true=2,mean_v.false=-1,
               sd_v.true=0.5)

# If everything is OK nothing is returned.
check.p.vector(p.vector,model)

# NB: order of names in p.vector doesn't matter
check.p.vector(p.vector[length(p.vector):1],model)

# Here is what happens if you get it wrong
check.p.vector(c(p.vector,A.a3=1),model)

# Here is how to check the way elements of p.vector are mapped to cells
print.cell.p(p.vector,model)
#  [1] "lf.a1.WW"
#    A b  t0 mean_v sd_v st0
# WW 3 7 0.2      2  0.5   0
# NW 3 8 0.2     -1  1.0   0
# 
# [1] "hf.a1.WW"
#    A b  t0 mean_v sd_v st0
# WW 3 7 0.2      2  0.5   0
# NW 3 8 0.2     -1  1.0   0
# 
# [1] "nw.a1.WW"
#    A b  t0 mean_v sd_v st0
# WW 3 7 0.2     -1  1.0   0
# NW 3 8 0.2      2  0.5   0
# 
# [1] "lf.a2.WW"
#    A b  t0 mean_v sd_v st0
# WW 2 6 0.2      2  0.5   0
# NW 2 7 0.2     -1  1.0   0
# 
# ...
# 
# [1] "nw.a2.NW"
#    A b  t0 mean_v sd_v st0
# NW 2 7 0.2      2  0.5   0
# WW 2 6 0.2     -1  1.0   0


# EXAMPLE 5: 3-way design with stimulus, factor A & factor B. Map both factors 
# A and B to parameter A, (so parameter names are A.a1.b1 etc.)

p.map <- list(A=c("A","B"),B="R",t0="1",mean_v="M",sd_v="M",st0="1")

model <- model.dmc(p.map,constants=const,
                   factors=list(S=c("lf","hf","nw"),A=c("a1","a2"),B=c("b1","b2")),
                   match.map=list(M=list(lf=1,hf=1,nw=2)),
                   responses=c("WW","NW"))

p.vector  <- c(A.a1.b1=3,A.a2.b1=4,A.a1.b2=1,A.a2.b2=2,B.WW=4,B.NW=5,t0=.2,
               mean_v.true=2,mean_v.false=-1,sd_v.true=0.5)

check.p.vector(p.vector,model)

print.cell.p(p.vector,model)

# # 24 CELLS = S=3 X A=2 X B=2 X R=2 
# 
# [1] "lf.a1.b1.WW"
#    A b  t0 mean_v sd_v st0
# WW 3 7 0.2      2  0.5   0
# NW 3 8 0.2     -1  1.0   0
# 
# [1] "hf.a1.b1.WW"
#    A b  t0 mean_v sd_v st0
# WW 3 7 0.2      2  0.5   0
# NW 3 8 0.2     -1  1.0   0
# 
# [1] "nw.a1.b1.WW"
#    A b  t0 mean_v sd_v st0
# WW 3 7 0.2     -1  1.0   0
# NW 3 8 0.2      2  0.5   0
# 
# [1] "lf.a2.b1.WW"
#    A b  t0 mean_v sd_v st0
# WW 4 8 0.2      2  0.5   0
# NW 4 9 0.2     -1  1.0   0
#  ...
# [1] "lf.a1.b2.WW"
#    A b  t0 mean_v sd_v st0
# WW 1 5 0.2      2  0.5   0
# NW 1 6 0.2     -1  1.0   0
# # ...
# [1] "nw.a2.b2.WW"
#    A b  t0 mean_v sd_v st0
# WW 2 6 0.2     -1  1.0   0
# NW 2 7 0.2      2  0.5   0
#  ...
# [1] "lf.a1.b1.NW"
#    A b  t0 mean_v sd_v st0
# NW 3 8 0.2     -1  1.0   0
# WW 3 7 0.2      2  0.5   0
# ...
# [1] "nw.a2.b2.NW"
#    A b  t0 mean_v sd_v st0
# NW 2 7 0.2      2  0.5   0
# WW 2 6 0.2     -1  1.0   0


# EXAMPLE 6: 3 x 2 X 2 Design with stimulus, factor A (with A.a1 and A.a2), 
# parameter B with factor B and R (so B.b1.w etc.)

# NB: The "R" factor MUST be put last in p.map!!
p.map <- list(A="A",B=c("B","R"),t0="1",mean_v="M",sd_v="M",st0="1")

model <- model.dmc(p.map,constants=const,
                   factors=list(S=c("lf","hf","nw"),A=c("a1","a2"),B=c("b1","b2")),
                   match.map=list(M=list(lf=1,hf=1,nw=2)),
                   responses=c("WW","NW"))


p.vector  <- c(A.a1=3,A.a2=4,B.b1.WW=4,B.b2.WW=5,B.b1.NW=6,B.b2.NW=7,
               t0=.2,mean_v.true=2,mean_v.false=-1,sd_v.true=0.5)

check.p.vector(p.vector,model)

# 24 cells
print.cell.p(p.vector,model)
# 
# [1] "lf.a1.b1.WW"
#    A b  t0 mean_v sd_v st0
# WW 3 7 0.2      2  0.5   0
# NW 3 9 0.2     -1  1.0   0
#
# ...
# 
# [1] "nw.a2.b2.NW"
#    A  b  t0 mean_v sd_v st0
# NW 4 11 0.2      2  0.5   0
# WW 4  9 0.2     -1  1.0   0

# EXAMPLE 7: 3 X 2 X 2 Design with stimulus, factor A (with A.a1 and A.a2), 
# parameter mean_v with mean_v.b1.true etc.

# NB: The "M" factor MUST be put last in p.map!!
p.map <- list(A="A",B="R",t0="1",mean_v=c("B","M"),sd_v="M",st0="1")

model <- model.dmc(p.map,constants=const,
                   factors=list(S=c("lf","hf","nw"),A=c("a1","a2"),B=c("b1","b2")),
                   responses=c("WW","NW"),match.map=list(M=list(lf=1,hf=1,nw=2)))


p.vector  <- c(A.a1=3,A.a2=4,B.WW=4,B.NW=5,t0=.2,sd_v.true=0.5,
               mean_v.b1.true=2,mean_v.b2.true=-1,mean_v.b1.false=3,mean_v.b2.false=-2)

check.p.vector(p.vector,model)

print.cell.p(p.vector,model)

# [1] "lf.a1.b1.WW"
#    A b  t0 mean_v sd_v st0
# WW 3 7 0.2      2  0.5   0
# NW 3 8 0.2      3  1.0   0
# ...
# [1] "nw.a2.b2.NW"
#    A b  t0 mean_v sd_v st0
# NW 4 9 0.2     -1  0.5   0
# WW 4 8 0.2     -2  1.0   0


