##################  DMC Lesson 5: Advanced Multiple subjects
#for this lesson we will also use another package, gridExtra, to stack ggplots together
install.packages("gridExtra") 



### Lesson 5.5:  Using ggplot2 with DMC samples.

rm(list=ls())
setwd(your_directory_here)
source("dmc/dmc.R")
require(gridExtra)

# Please send any questions about the lesson to luke.strickland@utas.edu.au

# The posterior prediction functions, post.predict.dmc and h.post.predict.dmc, 
# can save an attribute called "gglist", to be used with ggplot2.
# This attribute is a list with several elements, each containing a data frame.
# The first df contains data and posterior predictions for response proportion, 
# that is, the proportion of times a certain response is made to a stimulus. 
# Response proportion for the 'correct' response to a stimulus is the accuracy.
# The second element contains proportions of `non-responses', for designs in which
# responses can be omitted.
# The third element contains RT quantiles, for example the median, which is the 
# middle value in the RT distribution, the 0.1 quantile, which is the value 
# greater than or equal to 10% of observed responses, and the 0.9 quantile,
# which is the value greater or equal to 90% of observed responses.
# The RT quantiles that are calculated can be changed with the argument
# probs.gglist 
#The fourth and fifth element contain posterior predictive for mean RT,
# and SD of RTs, respectively.

# We use a simple LNR example, from lesson 4.4,
# but many of the techniques demonstrated are most useful
# for more complex designs.

load_data("dmc_5_5.RData") # samples3 from lesson 4.4
load_model ("LNR","lnrPP.R")

# Below is an example getting the gglist for a single participant.
# post.predict.dmc can get a list for ggplot by turning on the argument
# gglist = TRUE. It is saved as an attribute - "gglist"
#

lnrppredicts <- post.predict.dmc(samples3[[1]], gglist=TRUE)
attr(lnrppredicts, "gglist")
# $pps
#    S  R    mean  lower  upper  data
# 1 s1 r1 0.56580 0.4998 0.6280 0.516
# 2 s2 r1 0.43392 0.3737 0.4981 0.424
# 3 s1 r2 0.43420 0.3720 0.5002 0.484
# 4 s2 r2 0.56608 0.5019 0.6263 0.576
# 
# $RTs
#     S  R      mean     lower     upper      data quantile
# 1  s1 r1 0.3024461 0.2777056 0.3328993 0.2931947      0.1
# ...
# 5  s1 r1 0.5563816 0.4963841 0.6240004 0.5225163      0.5
# ...
# 12 s2 r2 1.5502984 1.2327825 2.0540991 1.5485011      0.9

#different number of RT quantiles - 5 not 3

lnrppredicts <- post.predict.dmc(samples3[[1]], gglist=T, probs.gglist=
                                   seq(0.1, 0.9, 0.2))
attr(lnrppredicts, "gglist")
# ...
# $RTs
#     S  R      mean     lower     upper      data quantile
# 1  s1 r1 0.3036429 0.2757735 0.3352187 0.2931947      0.1
# ...
# 8  s2 r2 0.4098602 0.3720858 0.4585516 0.4070878      0.3
# 9  s1 r1 0.5599552 0.5043002 0.6485499 0.5225163      0.5
# ...
# 16 s2 r2 0.8002950 0.7027513 0.9071369 0.8885420      0.7
# ...
# 20 s2 r2 1.5162413 1.1731450 1.9423548 1.5485011      0.9

# We will now run through an example with predictions for multiple participants 
# We use h.post.predict.dmc, which can save gglist as an attribute of the 
# "av" attribute of the output.
#
# The list contains means and posterior predictions for the average of all
# participants. The function gets these averages by stacking all the data 
# together into one big data frame then calculating all statistics, that is the 
# RT quantiles and response proportions, for the whole data frame. The same is 
# done for the simulated data (sim). The "mean" value of the sim is the average 
# of all the posterior predictions. The CI values provide some lower and upper 
# credible interval, which is the quantile of the posterior predictions. The 
# range of the credible intervals is by default 95%, but can be changed with the 
# argument CI.gglist.

lnr.group.ppredicts <- h.post.predict.dmc(samples3, gglist=TRUE)
attr(attr(lnr.group.ppredicts, "av"), "gglist")
# $pps
#    S  R     mean     lower     upper   data
# 1 s1 r1 0.651963 0.6393750 0.6644675 0.6502
# 2 s2 r1 0.348307 0.3379950 0.3578625 0.3503
# 3 s1 r2 0.348037 0.3355325 0.3606250 0.3498
# 4 s2 r2 0.651693 0.6421375 0.6620050 0.6497
# 
# $RTs
#     S  R      mean     lower     upper      data quantile
# 1  s1 r1 0.2865556 0.2826444 0.2909831 0.2852836      0.1
# ...
# 12 s2 r2 1.1871971 1.1477369 1.2234497 1.1768945      0.9


# The ggplot2 package provides researchers with intuitive graphics
# syntax to quickly produce plots. Some researchers may already be very familiar
# with ggplot, in which case they may wish to grab the "gglist" attribute from
# their PP object and do their own thing. We include three convenience functions 
# in DMC that generate quick ggplots of RT and response proportion data with 
# sensible settings. These functions can read the output of post.predict.dmc 
# and h.post.predict.dmc directly, provided that a gglist attribute was saved. 

# In the Bayesian manner "the data are the data" so are shown with no error 
# bars (as open circles joined by broken lines). The model fits are shown 
# as credible intervals, with the median indicated by a solid circle.

# The first convenience function plots the proportion of each response. 
ggplot.RP.dmc(lnr.group.ppredicts)

# The layout of these functions can be adjusted with panels.ppage (if there's 
# more than one page, the output is a list of ggplots element by element). Also 
# with nrow and ncol, e.g.:

# Set minimalist theme (theme defined in dmc code ) 

theme_set(theme_simple())

ggplot.RP.dmc(lnr.group.ppredicts, panels.ppage=1)
ggplot.RP.dmc(lnr.group.ppredicts, ncol=1)


# The second convenience function plots the RTs quantiles.
ggplot.RT.dmc(lnr.group.ppredicts)

# To bind the plots together using the gridExtra package:
grid.arrange(ggplot.RP.dmc(lnr.group.ppredicts), 
             ggplot.RT.dmc(lnr.group.ppredicts),
             layout_matrix = cbind(c(1,2), c(1,2)))


# The output of the dmc convenience functions are either a single ggplot object,
# (if one page), or a list of ggplot objects (one page per element). 
# Thus the output can be modified with ggplot syntax.
# For example add a different theme and line at chance response accuracy:
ggplot.RP.dmc(lnr.group.ppredicts) +  geom_hline(yintercept=0.5, linetype=3) +
  geom_text(aes(x= 0.8, y= 0.52, label="Chance"))+ theme_bw()

# For multiple page graphics (a list of ggplots) do the same thing but use lapply:
lapply(ggplot.RP.dmc(lnr.group.ppredicts, panels.ppage=1), 
       function(x) x +  geom_hline(yintercept=0.5, linetype=3) +
         geom_text(aes(x= 0.8, y= 0.52, label="Chance"))+ theme_bw())


# The title, x axis title, y axis title, x labels and y labels, 
# can all be modified easily with ggplot syntax, e.g.:
ggplot.RT.dmc(lnr.group.ppredicts) + 
  scale_x_discrete(labels=c("Response One","Response Two")) +
  xlab("Responses") +ylab ("Response Time") +
  ggtitle("LNR quantile RT predictions") +
  theme(plot.title = element_text(hjust = 0.5))


# Although the default output of the above functions may sometimes suffice,
# in some cases the user may want to modify the data in some way.
# For example, they may wish to subset the data, or modify the factor labels and
# levels.

# To do this, extract the attribute from the ppobject -
# for a post.predict.dmc object, stored as attribute "gglist"
# for a h.post.predict.dmc object, attribute "gglist" of attribute "av"

# For example, the user may only want to plot response accuracy,
# that is the proportion of correct responses to a stimulus. 
# A quick way to achieve this is to modify the data being plotted.

lnr.group.gglist <- attr(attr(lnr.group.ppredicts, "av"), "gglist")

# Thus, extract the rp df, take only the correct responses, and drop the 
# now redundant "R" column that denotes the observed response. The response
# column is redundant because only the correct response for each stimulus will 
# be plotted.
acc.df <- lnr.group.gglist[[1]]
acc.df <- acc.df[(acc.df$S=="s1" & acc.df$R=="r1") | 
                   (acc.df$S=="s2" & acc.df$R=="r2"),]
acc.df<- acc.df[-2]
acc.df
#    S     mean     lower     upper   data
# 1 s1 0.651963 0.6393750 0.6644675 0.6502
# 4 s2 0.651693 0.6421375 0.6620050 0.6497

# The ggplot dmc functions can accept data frames. 
# One modification necessary in this example is to change the x axis.
# Because response, the default x axis, was removed from the data frame, 
# something else needs to be specified. In this case, "S", i.e., the stimulus
# column.
ggplot.RP.dmc(acc.df, xaxis="S") + ylab("accuracy") + ylim(c(0.6, 0.7))

# Similarly, if you wish to examine only the response times of correct responses. 
correctonly.df <- lnr.group.gglist[[3]]
correctonly.df<- correctonly.df[(correctonly.df$S=="s1" &
                                   correctonly.df$R=="r1")| 
                                  (correctonly.df$S=="s2" & correctonly.df$R=="r2"),]
correctonly.df <- correctonly.df[-2]

# As response was dropped, again we need a different x axis. 
ggplot.RT.dmc(correctonly.df, xaxis="S") + ylab("Correct RT")

# The user may also wish to modify the labels of the panels, a
# quick way to do this is to modify the factor labels in the data frame e.g.,:
df.newlabels <- lnr.group.gglist[[3]]
levels(df.newlabels$S)<- c("Stimulus One", "Stimulus Two")
ggplot.RT.dmc(df.newlabels)

# Note this can be done more flexibly, without changing the data, 
# with labellers, but this approach can get a bit tricky: 
ggplot.RT.dmc(lnr.group.ppredicts) + 
  facet_wrap(c("S"),labeller=
               as_labeller(c(s1 = "Stimulus One", s2= "Stimulus Two"), multi_line=F))

# One way to change the order in which plots are presented is to change the 
# order in the data frame. 
df.neworder <- lnr.group.gglist[[3]]
df.neworder$S <- factor(df.neworder$S, levels= c("s2", "s1"))
ggplot.RT.dmc(df.neworder)

# The gglist attribute also has the mean RT and the sd of RT - elements
# [[3]] and [[4]]. They will work with the ggplot.RT.dmc
# function with do.quantiles=FALSE
ggplot.RT.dmc(lnr.group.gglist[[4]], do.quantiles=FALSE) +ylab("Mean RT")
ggplot.RT.dmc(lnr.group.gglist[[5]], do.quantiles=FALSE) +ylab("SD RT")

# The user may wish to aggregate the data at different 
# levels to the default. 
# A good way to do this is to first save the simulation 
# with post.predict, e.g.,
lnr.group.sim <- h.post.predict.dmc(samples3,save.sim=TRUE)

# Then make 'sim' and 'data' data frames.
# Stack each subject's similation into one data frame
sim <- do.call(rbind, lnr.group.sim)
# Do the same for the data
data <- lapply(lnr.group.sim, function(x) attr(x, "data"))
data <- do.call(rbind, data)

# To get a data frame of summary stats for ggplot, run get.fitgglist.dmc on this df.
# The advantage of this approach is the user can aggregate at different levels
# without re-simulating. That is, to treat all levels of a factor as one when 
# first calculating response proportion and RT quantiles for the ggplot df.

# For example, the stimulus factor could be ignored by setting factors = NULL
# In designs with more factors, to ignore a factor out then set factors to a 
# vector that contains the name of all remaining factors not to be ignored.

noS <- get.fitgglist.dmc(sim, data, factors=NULL, noR=F)

ggplot.RP.dmc(noS[[1]]) 
ggplot.RT.dmc(noS[[3]]) 

# It is also possible to drop the "R" argument, corresponding to the observed
# response that participants made. In other words, to aggregate the observed RT
# for all responses on a type of trial together.

noRs <- get.fitgglist.dmc(sim, data, noR=TRUE)

# Dropping the "R" argument will cause the get.fitgglist.dmc function to get 
# response accuracy rather than "response proportion" by default. 

ggplot.RP.dmc(noRs[[1]])

# Response accuacy is scored by checking whether the levels of S and R are equal, 
# but this will not work in some cases. However, the user is able to specify any 
# scoring function they want with the acc.fun argument. The default form is:
acc.fun=function(x){as.numeric(x$S)==as.numeric(x$R)}
# Note the "as.numeric" here means that the corresponding elements of S and R
# need only be in the same order, but may have different labels.

# Suppose the response factor had levels "word" and "nonword" and a factor WF
# (word frequency) was used for scoring with three levels "common" "rare" and 
# "never". Then one might use
#
acc.fun=function(x){(x$WF=="nonword" & x$R=="never") | 
    (x$WF!="nonword" & x$R!="never")}
# 
# If the acc.fun argument is NULL, response proportion will be returned. 

# The convenience function gplot.RA.dmc drops the proportion error, which is 
# redundant, and adds a better ylab. 

ggplot.RA.dmc(noRs[[1]]) 


# The noR argument drops out the response argument for RTs (i.e., aggregates 
# response times over all responses rather than separating them).
#
# By default ggplot.RT.dmc uses the response factor 'R' as the xaxis. The 
# default will casue an error here as there is no reponse factor in noRs. 
# Instead we use the S(timulus) factor.
ggplot.RT.dmc(noRs[[3]], xaxis='S')





