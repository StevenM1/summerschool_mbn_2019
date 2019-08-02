##################  DMC Lesson 5: Advanced Models and Data
#
# Note: Before running these scripts, the current working directory must be set 
# to the top-level folder containing the dmc and tutorial subfolders 

### Lesson 5.2: Advanced Scoring 

rm(list=ls()); par(mfrow = c(1,1)) 
source ("dmc/dmc.R")

load_model ("LBA","lba_B.R") # LBA model with B=b-A parameterization 

######## Section 1: Scoring with more than one factor

# Sometimes scoring what is a correct response is depends on more than one 
# factor. For example, where there is a "compatible" response condition where 
# you respond with a left button press to a stimulus on the left and a right 
# button press to a stimulus on the right, but also an "incompatible" response 
# condition where you have to respond with the hand opposite the stimulus 
# position. 

# Example 1: SxA design, scoring depends on both S and A

# NB: The order of names in "factors" requires S and A to be together and 
#     ordered so rows in model contain contiguous strings "s1.a1" etc.

model <- model.dmc(
  p.map=list(A="1",B="1",t0="1",mean_v="M",sd_v="1",st0="1"),
  factors=list(S=c("s1","s2"),A=c("a1","a2")),
  responses=c("r1","r2"),constants=c(sd_v=1,st0=0),
  match.map=list(M=list(s1.a1="r1",s1.a2="r2",s2.a1="r2",s2.a2="r1"))
)

# Recall e.g., s1.a1.r1 is cell S=1,A=1,response=r1 so e,g., 
# s1.a1 and s2.a2 cells for accumulator r1 (, , r1) and accumulator 2 (, , r2) 
# are reversed.
model
# , , r1
# 
#             A    B   t0 mean_v.true mean_v.false sd_v  st0
# s1.a1.r1 TRUE TRUE TRUE        TRUE        FALSE TRUE TRUE
# s2.a1.r1 TRUE TRUE TRUE       FALSE         TRUE TRUE TRUE
# s1.a2.r1 TRUE TRUE TRUE       FALSE         TRUE TRUE TRUE
# s2.a2.r1 TRUE TRUE TRUE        TRUE        FALSE TRUE TRUE
# s1.a1.r2 TRUE TRUE TRUE        TRUE        FALSE TRUE TRUE
# s2.a1.r2 TRUE TRUE TRUE       FALSE         TRUE TRUE TRUE
# s1.a2.r2 TRUE TRUE TRUE       FALSE         TRUE TRUE TRUE
# s2.a2.r2 TRUE TRUE TRUE        TRUE        FALSE TRUE TRUE
# 
# , , r2
# 
#             A    B   t0 mean_v.true mean_v.false sd_v  st0
# s1.a1.r1 TRUE TRUE TRUE       FALSE         TRUE TRUE TRUE
# s2.a1.r1 TRUE TRUE TRUE        TRUE        FALSE TRUE TRUE
# s1.a2.r1 TRUE TRUE TRUE        TRUE        FALSE TRUE TRUE
# s2.a2.r1 TRUE TRUE TRUE       FALSE         TRUE TRUE TRUE
# s1.a1.r2 TRUE TRUE TRUE       FALSE         TRUE TRUE TRUE
# s2.a1.r2 TRUE TRUE TRUE        TRUE        FALSE TRUE TRUE
# s1.a2.r2 TRUE TRUE TRUE        TRUE        FALSE TRUE TRUE
# s2.a2.r2 TRUE TRUE TRUE       FALSE         TRUE TRUE TRUE


p.vector  <- c(A=1,B=1,t0=.2,mean_v.true=2,mean_v.false=0)
check.p.vector(p.vector,model)
print.cell.p(p.vector,model)
# [1] "s1.a1.r1"
#    A b  t0 mean_v sd_v st0
# r1 1 2 0.2      2    1   0
# r2 1 2 0.2      0    1   0
# 
# [1] "s2.a1.r1"
#    A b  t0 mean_v sd_v st0
# r1 1 2 0.2      0    1   0
# r2 1 2 0.2      2    1   0
# ...


# Example 2: 2 x 2 x 2 (S x A x B) design, scoring depends on S and B 

model <- model.dmc(
  p.map=list(A="1",B="1",t0="1",mean_v="M",sd_v="1",st0="1"),
  factors=list(S=c("s1","s2"),B=c("b1","b2"),A=c("a1","a2")),
  responses=c("r1","r2"),constants=c(sd_v=1,st0=0),
  match.map=list(M=list(s1.b1="r1",s1.b2="r2",s2.b1="r2",s2.b2="r1"))
)

# NB: "model.dmc" checks that all cells are scored, returns an error if not
#     e.g., s2.a2="r1" missing in match.map as in the following example.

# model <- model.dmc(
#   p.map=list(A="1",B="1",t0="1",mean_v="M",sd_v="1",st0="1"),
#   factors=list(S=c("s1","s2"),A=c("a1","a2")),
#   responses=c("r1","r2"),constants=c(sd_v=1,st0=0),
#   match.map=list(M=list(s1.a1="r1",s1.a2="r2",s2.a1="r2"))
# )


######## Section 2:  Scoring with more than one choice

# 3 stimulus and 3 choice example, with binary scoring

model.dmc(
  factors=list(S=c("s1","s2","s3")),
  responses=c("r1","r2","r3"),
  p.map=list(A="1",B="R",t0="1",mean_v="M",sd_v="1",st0="1"),
  match.map=list(M=list(s1="r1",s2="r2",s3="r3")),
  constants=c(sd_v=1,st0=0))

# e.g. for accumulator r1 (, , r1) only s1 cells are TRUE 

# , , r1
# 
#          A B.r1  B.r2  B.r3   t0 mean_v.true mean_v.false sd_v  st0
# s1.r1 TRUE TRUE FALSE FALSE TRUE        TRUE        FALSE TRUE TRUE
# s2.r1 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s3.r1 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s1.r2 TRUE TRUE FALSE FALSE TRUE        TRUE        FALSE TRUE TRUE
# s2.r2 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s3.r2 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s1.r3 TRUE TRUE FALSE FALSE TRUE        TRUE        FALSE TRUE TRUE
# s2.r3 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s3.r3 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# ...


# 3 x 2 (S x a) design, reverse scoring for second level of A 

model.dmc(
  factors=list(S=c("s1","s2","s3"),A=c("a1","a2")),
  responses=c("r1","r2","r3"),
  p.map=list(A="1",B="R",t0="1",mean_v="M",sd_v="1",st0="1"),
  match.map=list(M=list(s1.a1="r1",s2.a1="r2",s3.a1="r3",
                        s1.a2="r3",s2.a2="r2",s3.a2="r1")),
  constants=c(sd_v=1,st0=0))

# e.g. for response r1 (, , r1) both s1.a1 and s3.a2 cells are TRUE 

# , , r1
# 
#             A B.r1  B.r2  B.r3   t0 mean_v.true mean_v.false sd_v  st0
# s1.a1.r1 TRUE TRUE FALSE FALSE TRUE        TRUE        FALSE TRUE TRUE
# s2.a1.r1 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s3.a1.r1 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s1.a2.r1 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s2.a2.r1 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s3.a2.r1 TRUE TRUE FALSE FALSE TRUE        TRUE        FALSE TRUE TRUE
# s1.a1.r2 TRUE TRUE FALSE FALSE TRUE        TRUE        FALSE TRUE TRUE
# s2.a1.r2 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s3.a1.r2 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s1.a2.r2 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s2.a2.r2 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s3.a2.r2 TRUE TRUE FALSE FALSE TRUE        TRUE        FALSE TRUE TRUE
# s1.a1.r3 TRUE TRUE FALSE FALSE TRUE        TRUE        FALSE TRUE TRUE
# s2.a1.r3 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s3.a1.r3 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s1.a2.r3 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s2.a2.r3 TRUE TRUE FALSE FALSE TRUE       FALSE         TRUE TRUE TRUE
# s3.a2.r3 TRUE TRUE FALSE FALSE TRUE        TRUE        FALSE TRUE TRUE
# ...


######## Section 3: Non-factorial scoring with 3 choices: 27 CASES in match.map!
#
# This section demonstrates scoring in a very complex case, demonstrating
# how to use dmc when very flexible scoring is required. 
# 
#
# Stroop task example with 3 choices with different rates for two types of 
# correct response and two types of error response
#
# Factors: relevant stimulus (ie., ink colour) RS = R/G/B (red/green/blue) 
#     and: irrelevant stimulus (i.e., word) IS = r/g/b
# Accumulators specified by response factor R = R/G/B
#
# NB: Have to use different level labels for RS and IS. However, it is 
#     convenient to use upper vs. lower case so you can still test equality 
#     by ignoring case as detailed below. More generally you can specify your
#     own function to make this sort of equality test.
#
#       Suppose you wish to use the following scoring, two types of correct:
#              R=RS=IS     (congruent)    mean_v = dm (double match)
#              R=RS, R<>IS (incongruent)  mean_v = rm (relevant match)
#       and two types of error:
#              R=IS, R<>RS (IS error)     mean_v = im (irrelevant match)
#              R<>RS,R<>IS (error)        mean_v = nm (no match)
#
#     You could then estimate corresponding mean_v for dm,rm,im and nm (note
#     that none of these levels can be called "true" or "false")
#
#
# First create a factor with all 27 cases (3x3x3) that you have to map
# either dm,rm,im or nm to, using the convenience function "empty.map",
# where the first entry is a list of the experimental factors and response 
# factor (and their levels) over which the scoring is done, and the output
# is a factor vector with no values yet assigned

map <- empty.map(list(RS=c("R","G","B"),IS=c("r","g","b"),
                      R=c("R","G","B")),levels=c("dm","rm","im","nm"))

map
# R.r.R G.r.R B.r.R R.g.R G.g.R B.g.R R.b.R G.b.R B.b.R 
# <NA>  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>  <NA> 
#   R.r.G G.r.G B.r.G R.g.G G.g.G B.g.G R.b.G G.b.G B.b.G 
# <NA>  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>  <NA> 
#   R.r.B G.r.B B.r.B R.g.B G.g.B B.g.B R.b.B G.b.B B.b.B 
# <NA>  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>  <NA>
#   Levels: dm rm im nm

# NB: You must include ALL factors in the design in the empty map in the SAME
#     ORDER as in the models factor = argument. This is true even for factors 
#     that are IRRELEVANT to the mapping (to accommodate these, values have to
#     be repeated appropriately).

# You could now assign appropriate values by hand. For example the first entry
# in map is "R.r.R, being RS = red, IS = red for accumulator R, which is a dm
#
# You could also do it programmatically using the elements of the "map" names 
# between the periods. In some cases this requires making a wrong assignment
# for a subset (e.g, the dm cases in this example) and then fixing later by 
# overwriting. In general do all of the single equality conditions first, then
# do double equality conditons (like dm), then triple etc.
#
# For example, "rm" is appropriate where the first and third entries (RS and R)
# are equal. This can be done with the convenience function "assign.map". The
# eq.list parameter specifies a list containing vectors specifying the position
# of the elements of the name to compare.
#
map <- assign.map(map,value="rm",eq.list=list(c(1,3)))
#
# R.r.R G.r.R B.r.R R.g.R G.g.R B.g.R R.b.R G.b.R B.b.R 
#    rm  <NA>  <NA>    rm  <NA>  <NA>    rm  <NA>  <NA> 
# R.r.G G.r.G B.r.G R.g.G G.g.G B.g.G R.b.G G.b.G B.b.G 
#  <NA>    rm  <NA>  <NA>    rm  <NA>  <NA>    rm  <NA> 
# R.r.B G.r.B B.r.B R.g.B G.g.B B.g.B R.b.B G.b.B B.b.B 
#  <NA>  <NA>    rm  <NA>  <NA>    rm  <NA>  <NA>    rm 
#
# Note that this is wrong for when all entries are equal (e.g., the first) but 
# we will deal with that later by overwriting.

# For "im" equality ignoring case is required for the second vs. third elements. 
# "assign.map" has a "funs" argument specifying a list of lists that give
# pairs of functions to transform the elements before checking equality.
#
# Note: the position of arguments for funs=list(list("toupper","") corresponds 
#       to the eq.list=list(c(2,3))
map <- assign.map(map,value="im",eq.list=list(c(2,3)),
                  funs=list(list("toupper","")))
#
# R.r.R G.r.R B.r.R R.g.R G.g.R B.g.R R.b.R G.b.R B.b.R 
#    im    im    im    rm  <NA>  <NA>    rm  <NA>  <NA> 
# R.r.G G.r.G B.r.G R.g.G G.g.G B.g.G R.b.G G.b.G B.b.G 
#  <NA>    rm  <NA>    im    im    im  <NA>    rm  <NA> 
# R.r.B G.r.B B.r.B R.g.B G.g.B B.g.B R.b.B G.b.B B.b.B 
#  <NA>  <NA>    rm  <NA>  <NA>    rm    im    im    im 
#
# NB1: A "" in funs is automatically expanded to "identity", meaning
#      no transformation, so funs=list(list("toupper","") is equivalent to
#      funs=list(list("toupper","identity"). If no funs argument is specified
#      funs is created automatically with all "identity" entries.
# NB2: funs=list(list("","tolower") would also work.
#
# For "dm" both equalities must hold. "assign.map" allows two or more equality
# tests to be made and the result to depend on all being true by adding 
# extra terms to the eq.list and funs arguments
#
eq.list <- list(c(1,3),c(2,3))
funs <- list(list("",""),list("toupper",""))
map <- assign.map(map,value="dm",eq.list=eq.list,funs=funs)
#
# R.r.R G.r.R B.r.R R.g.R G.g.R B.g.R R.b.R G.b.R B.b.R 
#    dm    im    im    rm  <NA>  <NA>    rm  <NA>  <NA> 
# R.r.G G.r.G B.r.G R.g.G G.g.G B.g.G R.b.G G.b.G B.b.G 
#  <NA>    rm  <NA>    im    dm    im  <NA>    rm  <NA> 
# R.r.B G.r.B B.r.B R.g.B G.g.B B.g.B R.b.B G.b.B B.b.B 
#  <NA>  <NA>    rm  <NA>  <NA>    rm    im    im    dm 
#
# NB: This example requires dm to be run last, but with appropriate (more
#     complicated) arguments to funs for im and rm (so they do not make
#     assignments to the dm cases) any order can be used.
#
# Finally nm can be assigned to all of the entries of map that remain empty. 
# This is the default behaviour of "assign.map" so is accomplished by
map <- assign.map(map,value="nm")
#
# R.r.R G.r.R B.r.R R.g.R G.g.R B.g.R R.b.R G.b.R B.b.R 
#    dm    im    im    rm    nm    nm    rm    nm    nm 
# R.r.G G.r.G B.r.G R.g.G G.g.G B.g.G R.b.G G.b.G B.b.G 
#    nm    rm    nm    im    dm    im    nm    rm    nm 
# R.r.B G.r.B B.r.B R.g.B G.g.B B.g.B R.b.B G.b.B B.b.B 
#    nm    nm    rm    nm    nm    rm    im    im    dm 


# Unfortunatley the names used for factor levels and responses above will cause
# problems with the model.dmc function, as some names are not unique. Here
# is a version that avoids these problems but requires a more sophisticated way
# of setting up the map, using an inline function.

# Here an inline function can be used to pick out first letters and set them to 
# the same case for equality testing
map <- empty.map(list(RS=c("Rr","Gr","Br"),IS=c("ri","gi","bi"),
                      R=c("red","green","blue")),levels=c("dm","rm","im","nm"))
funs=list(list(function(x){tolower(substr(x,1,1))},function(x){tolower(substr(x,1,1))}),
          list(function(x){tolower(substr(x,1,1))},function(x){tolower(substr(x,1,1))}))
map <- assign.map(map,value="rm",eq.list=list(c(1,3)),funs=funs[-2])
map <- assign.map(map,value="im",eq.list=list(c(2,3)),funs=funs[-2])
map <- assign.map(map,value="dm",eq.list=list(c(1,3),c(2,3)),funs=funs)
map <- assign.map(map,value="nm")
#   Rr.ri.red   Gr.ri.red   Br.ri.red   Rr.gi.red   Gr.gi.red   Br.gi.red 
#          dm          im          im          rm          nm          nm 
#   Rr.bi.red   Gr.bi.red   Br.bi.red Rr.ri.green Gr.ri.green Br.ri.green 
#          rm          nm          nm          nm          rm          nm 
# Rr.gi.green Gr.gi.green Br.gi.green Rr.bi.green Gr.bi.green Br.bi.green 
#          im          dm          im          nm          rm          nm 
#  Rr.ri.blue  Gr.ri.blue  Br.ri.blue  Rr.gi.blue  Gr.gi.blue  Br.gi.blue 
#          nm          nm          rm          nm          nm          rm 
#  Rr.bi.blue  Gr.bi.blue  Br.bi.blue 
#          im          im          dm 
# Levels: dm rm im nm

# The map can then be passed to "model.dmc" as a named entry in the match.map 
# list argument, and the name (RI below) used to parameterize mean_v in the  
# p.map argument. M must be passed in match.map in the usual way (making an
# assignment for every response to a level of a factor in the design) but in 
# this example is never used.
model <- model.dmc(
  factors=list(RS=c("Rr","Gr","Br"),IS=c("ri","gi","bi")),
  responses=c("red","green","blue"),
  p.map=list(A="1",B="1",t0="1",mean_v="RI",sd_v="1",st0="1"),
  match.map=list(M=list(Rr="red",Gr="green",Br="blue"),RI=map),
  constants=c(sd_v=1,st0=0) )
# Note: we could also just use list(R=1,G=2,B=3), i.e. in responses order, for 
#       the M parameter.

# To see this, first create a parameter vector to map. 
p.vector  <- c(A=1,B=1,t0=.2,
               mean_v.dm=4,mean_v.rm=3,mean_v.im=2,mean_v.nm=0)

check.p.vector(p.vector,model)

print.cell.p(p.vector,model)
# [1] "Rr.ri.red"
#       A b  t0 mean_v sd_v st0
# red   1 2 0.2      4    1   0
# green 1 2 0.2      0    1   0
# blue  1 2 0.2      0    1   0
# 
# [1] "Gr.ri.red"
#       A b  t0 mean_v sd_v st0
# red   1 2 0.2      2    1   0
# green 1 2 0.2      3    1   0
# blue  1 2 0.2      0    1   0
# ...

# Given this is so complicated, it's a good idea to simulate some data and check 
# that profiles make sense

data.model <- data.model.dmc(simulate.dmc(p.vector,model,n=1e3),model)

# Can be slow ...
par(mfrow=c(2,4))
profile.dmc("A",.1,2,p.vector,data.model)
profile.dmc("B",.1,2,p.vector,data.model)
profile.dmc("mean_v.dm",3,5,p.vector,data.model)
profile.dmc("mean_v.rm",2,4,p.vector,data.model)
profile.dmc("mean_v.im",1,3,p.vector,data.model)
profile.dmc("mean_v.nm",-1,1,p.vector,data.model)
profile.dmc("t0",.01,.4,p.vector,data.model)

