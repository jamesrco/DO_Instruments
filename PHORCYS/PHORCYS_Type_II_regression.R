# Type II regression of PHORCYS data and Winklers/shipboard incubations

library(smatr)
library(lmodel2)

setwd("/Users/jrcollins/Code/DO_Instruments/PHORCYS/data/processed/")

# read in data; create variables

rd <- read.csv("PHORCYS_other_methods.csv")

PHOR <- rd[,1]
PHOR_err <- rd[,2]
SB_inc <- rd[,5]
SB_inc_err <- rd[,6]
Wink <- rd[,3]
Wink_err <- rd[,4]

# fit incubation data, Winkler data to PHORCYS data using Type II major axis regression via lmodel2

fitMA.SB_inc <- lmodel2(SB_inc~PHOR)
intfitMA.SB_inc <- fitMA.SB_inc$regression.results[2,2]
slopefitMA.SB_inc <- fitMA.SB_inc$regression.results[2,3]

fitMA.Wink <- lmodel2(Wink~PHOR)
intfitMA.Wink <- fitMA.Wink$regression.results[2,2]
slopefitMA.Wink <- fitMA.Wink$regression.results[2,3]

# evaulate the two regression models at the first and last values of PHOR (needed for plotting of the two regressions in Igor)

PHOR_pred.SB_inc <- intfitMA.SB_inc + PHOR[c(1,6)]*slopefitMA.SB_inc
PHOR_pred.Wink <- intfitMA.Wink + PHOR[c(1,6)]*slopefitMA.Wink

