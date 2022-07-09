########################################################################################################
## ANOVA / DiffEx -- create ANOVAout dataframe of tests for differential expression/abundance
##                -- create PDF and HTML Volcanoe Plots
##
## By Eric Dammer, Duc Duong, and Qiudong Deng
#########################################
## NOTES
##
## - can repeat for different subgroup comparisons, if necessary
## - parallelized function. Requires doParallel, parallel, and dependencies
## - writes DEX table to pipeline variable ANOVAout and .csv, optionally
## - 0 Tukey values cannot be resuced with precision if <1e-9 or -10;
## - Option to fallback from 0 Tukey p values to Bonferroni-corrected T test (unequal variance)
## - if only 2 comparison groups are present, ANOVA overall p value is equivalent to T test; FDR correction for all proteinwide comparisons provided in that special case (default="BH").
## - avoid dashes ("-") in group names (set in vector of character strings to Grouping parameter)
##
##
#########################################
## Required Loaded Data and Parameters ##
#########################################

rootdir="e:/SysBioPipeline/"
setwd(rootdir)
load("SeyfriedPipelineOutput-test.RData")


# Parameters (default for function if unspecified here are generally the same as below in this template)
                         Grouping=numericMeta$Group                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.
                         parallelThreads=30                             # number of CPU threads to speed calculation (recommended, 2 or more):
                         NETcolors=net$colors                           # list net with slot/vector containing module color assignments; length of vector must be equal to number of rows in cleanDat.
                         twoGroupCorrMethod="BH"                        # default method for full FDR correction when only 2 groups present is Benjamini-Hochberg; see p.adjust(..., methods= ) options.
                         outputCSV=TRUE                                 # Output Full Table of statistics?  TRUE/FALSE
                         outFilePrefix="4"                              # typically "4", or step # in pipeline being run
                         outFileSuffix=FileBaseName                     # A description of the project, used as a filename suffix
                         fallbackIfZeroTukeyP=TRUE 


## Create parallel backend before running the function -- what options you use depends on your computer and OS.
  library(parallel)
  library(doParallel)
  stopCluster(clusterLocal)
  clusterLocal <- snow::makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
  registerDoParallel(clusterLocal)


source("./parANOVA.dex.R")
ANOVAout <- parANOVA.dex(cleanDat)

## Before plotting volcanoes, choose comparison (column #s) to flip denominator in log2(ratios);
## - any column numbers specified flip the volcano X axis for that pairwise comparison.
head(ANOVAout)
#                  F-Value       Pr(>F)      CT-ADAD    LOAD-ADAD   sEOAD-ADAD      LOAD-CT     sEOAD-CT  sEOAD-LOAD diff CT-ADAD diff LOAD-ADAD diff sEOAD-ADAD diff LOAD-CT diff sEOAD-CT diff sEOAD-LOAD NETcolors
#DYNC1H1|Q14204   6.842056 1.662086e-04 8.395118e-04 9.989343e-01 6.355620e-02 5.475599e-03 1.544454e-01 0.133368891 -0.033488813    0.001707582     -0.02036040  0.035196394   0.013128415    -0.022067980     green
#SPTAN1|Q13813    4.803079 2.676945e-03 2.985694e-03 1.132990e-01 2.334785e-03 9.832338e-01 9.964918e-01 0.994515348  0.023474420    0.020485155      0.02244106 -0.002989265  -0.001033357     0.001955908 lightcyan
#ANK2|Q01484     30.578501 7.963266e-18 9.321317e-12 2.661617e-04 1.374872e-09 7.035131e-03 1.296163e-05 0.973271732  0.128902746    0.075660044      0.08243078 -0.053242702  -0.046471967     0.006770735    purple
#SPTBN1|Q01082    1.301731 2.734578e-01 7.060692e-01 9.931906e-01 2.802188e-01 9.372470e-01 8.333703e-01 0.673174551  0.007426268    0.002526550      0.01160397 -0.004899718   0.004177701     0.009077419      grey
#PLEC|Q15149    116.386227 2.931614e-54 2.255820e-23 3.634339e-05 6.752541e-02 2.246769e-10 4.576154e-55 0.004697824 -0.668760920   -0.287562376     -0.10669039  0.381198545   0.562070535     0.180871990 turquoise
#SYNE1|Q8NF91    33.090858 3.909341e-19 8.525109e-10 1.189322e-02 3.295228e-03 6.598376e-04 3.003174e-15 0.854027054  0.114212304    0.053553888      0.04160811 -0.060658416  -0.072604192    -0.011945776      blue



## Generate Volcano plots, PDF and possibly interactive (mouseover) HTML
#########################################################################

########################
##   Volcano Plots    ##
########################
## Default Parameters ##
## - set in fn call - ##
########################

FCmin=0                    # 0.25 for 25%, 0 for no threshold (vertical minimum FC threshold dashed lines)
selectComps="ALL"          # "ALL" for volcano output(s) on all pairwise comparisons in ANOVAout
flip=c()                   # p value column numbers in which to swap denominator of pair for x axis range (gene products high in denominator, will be on left)
signifP=0.05               # p value threshold for counting Differential Expression points
useNETcolors=TRUE          # use module colors saved to ANOVAout, if available; otherwise when FALSE, specify downColor upColor, and NCcolor (must be valid R color specifications in quotes)
downColor="royalblue"      # significant points above/beyond thresholds on the upper left are this color if useNETcolors=FALSE
upColor="red"              # significant points above/beyond thresholds on the upper right are this color if useNETcolors=FALSE
NCcolor="grey"             # points not significant are this color if useNETcolors=FALSE
splitColors=FALSE          # create a separate volcano plot(s) for each color in an outputfigs/splitVolcanoes subfolder (folder created if it does not exist)
highlightGeneProducts=c()  # c("APP|P05067","MAPT|P10636","APOE|P02649") ; a list of uniqueID rownames to highlight as larger gold points. If symbolsOnly=TRUE, this can be a list of symbols, like c("APP","SMOC1","MAPT")
symbolsOnly=FALSE          # for mouse-over HTML plots and the above highlight callouts, consider only displaying and using official gene symbol from first part of UniqueID rownames of ANOVAout.
HTMLout=TRUE               # output interactive HTML copies that can be opened in browser. Requires plotly package.
outFilePrefix="4"          # typically "4", or step # in pipeline being run
outFileSuffix=FileBaseName # A description of the project, used as a filename suffix
outputFigs=getwd()         # Location to save figure file output(s)

plotVolc(ANOVAout, flip=c(3,4,5))

# not run: highlight gene products of interest:
# plotVolc(ANOVAout, flip=c(3,4,5), highlightGeneProducts=c("APP","SMOC1","MAPT"), symbolsOnly=TRUE)

# not run:  highlight gene products of interest by symbol, and do not use WGCNA module colors:
# plotVolc(ANOVAout, flip=c(3,4,5), useNETcolors=FALSE, highlightGeneProducts=c("APP","SMOC1","MAPT"), symbolsOnly=TRUE)






#save.image(paste0("4.saved.image.",FileBaseName,".RData"))  #overwrites if present
