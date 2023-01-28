########################################################################################################
## ANOVA / DiffEx -- 3 core user functions:
##
##   parANOVA.dex      -- create ANOVAout dataframe of tests for differential expression/abundance
##   plotVolc          -- create PDF and HTML Volcano Plots, output volcano settings to variables used later
##   DEXpercentStacked -- create PDF 
##
## By Eric Dammer, Duc Duong, and Qiudong Deng
#########################################
## NOTES
##
## - Output ANOVA+Tukey pairwise stats for volcano and downstream analyses
## - if only 2 comparison groups are present, ANOVA overall p value is equivalent to T test;
##   FDR correction for all proteinwide comparisons provided in that special case (default twoGroupCorrMethod="BH")
## - can be repeat for different subgroup comparisons, if necessary
## - parallelized function. Requires R packages doParallel, parallel, and dependencies
## - writes DEX table typically to pipeline variable ANOVAout and .csv table
## - 0 Tukey values are inaccurate when <1e-9 or -10; Estimates of very small values become 0 in R base Tukey post hoc calculations.
## - So, we implement an option to fallback from Tukey p values less than 1e-10 to Bonferroni-corrected T test (unequal variance)
## - avoid dashes ("-") in strings representing comparison groups (Grouping vector); they will be substituted with '.'
##
##
#########################################
## Required Loaded Data and Parameters ##
#########################################

rootdir="e:/SysBioPipeline/"
setwd(rootdir)
load("SeyfriedPipelineOutput-test.RData")


# Parameters (function may fallback to defaults if unspecified)
Grouping=numericMeta$Group                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat; same # of values as samples.
parallelThreads=30                             # number of CPU threads to speed calculation (recommended, 2 or more):
NETcolors=net$colors                           # list net with slot/vector containing module color assignments; length of vector must be equal to number of rows in cleanDat.
twoGroupCorrMethod="BH"                        # default method for full FDR correction when only 2 groups present is Benjamini-Hochberg; see p.adjust(..., methods= ) options.
outputCSV=TRUE                                 # Output Full Table of statistics?  TRUE/FALSE
outFilePrefix="4"                              # typically "4", or step # in pipeline being run; for output file sorting by filename.
outFileSuffix="DeepADproteome"                 # A description of the project, used as a filename suffix
fallbackIfSmallTukeyP=TRUE                     # Inaccurate Tukey p values < 1e-10 will not replaced with reliable Bonferroni FDR from T test for the pairwise comparison.


source("./parANOVA.dex.R")
ANOVAout <- parANOVA.dex()                     # runs on cleanDat and Grouping variables as required input.
#...Tukey p<1e-10 Fallback calculations using Bonferroni corrected T test: 8167 [15%]


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

#-------------- ---------- ------------ THESE COLUMNS ARE TUKEY POST-HOC P VALUE COLUMNS FOR PAIRWISE COMPARISONS---  THESE COLUMNS ARE LOG2 MEAN GROUP DIFFERENCES FOR THE SAME PAIRWISE COMPARISONS------

########################################
## Volcano Plots: PDF and HTML        ##
########################################
## Parameters set as variables here.  ##
########################################
## Most will be set to default values ##
## if not explicitly specified.       ##
########################################

FCmin=0                    # 0.25 for 25%, 0 for no threshold (vertical minimum FC threshold dashed lines)
selectComps="ALL"          # "ALL" for volcano output(s) on all pairwise comparisons in ANOVAout
flip=c(3,4,5)              # ANOVAout column index numbers for p values in which to swap denominator of pair for x axis range (gene products high in denominator, will be on left)
                           # As a general rule, the group with less severe effects is usually set to be the denominator (represented by what is right of '-' in ANOVAout column names)
signifP=0.05               # p value threshold for counting Differential Expression points
useNETcolors=TRUE          # use module colors saved to ANOVAout, if available; otherwise when FALSE, specify downColor upColor, and NCcolor (must be valid R color specifications in quotes)
downColor="royalblue"      # significant points above/beyond thresholds on the upper left are this color if useNETcolors=FALSE
upColor="red"              # significant points above/beyond thresholds on the upper right are this color if useNETcolors=FALSE
NCcolor="grey"             # points not significant are this color if useNETcolors=FALSE
splitColors=FALSE          # create a separate volcano plot(s) for each color in an outputfigs/splitVolcanoes subfolder (folder created if it does not exist)
highlightGeneProducts=c()  # c("APP|P05067","MAPT|P10636","APOE|P02649") ; a list of uniqueID rownames to highlight as larger gold points. If symbolsOnly=TRUE, this can be a list of symbols, like c("APP","SMOC1","MAPT")
labelHighlighted=FALSE     # if true, highlighted spots get text labels with their rownames from ANOVAout
symbolsOnly=FALSE          # for mouse-over HTML plots and the above highlight callouts, consider only displaying and using official gene symbol from first part of UniqueID rownames of ANOVAout.
labelTop=0                 # maximum p below which to label all points in the PDF output; OR an integer number of top ranked most significant points to label
labelSize=4.5              # text label font size, if any labels are found (when labelHighlighted=TRUE or labelTop>0)
sameScale=FALSE            # When multiple plots are drawn, should they all be the same scale with min and max x and y ranges?
HTMLout=TRUE               # output interactive HTML copies that can be opened in browser. Requires plotly package.
outFilePrefix="4"          # typically the step # in the pipeline being run
outFileSuffix="DeepADproteome"
                           # A description of the project, used as a filename suffix
outputfigs=getwd()         # Location to save figure file output(s)

plotVolc()                 # runs on ANOVAout as input (need not be specified).

## not run: highlight gene products of interest:
# highlightGeneProducts=c("APP","SMOC1","MAPT")
# symbolsOnly=TRUE
# plotVolc(ANOVAout)

## not run:  highlight gene products of interest by symbol, and do not use WGCNA module colors:
# flip=c(3,4,5)
# useNETcolors=FALSE
# highlightGeneProducts=c("APP","SMOC1","MAPT")
# symbolsOnly=TRUE
# plotVolc(ANOVAout)



#################################################################
## DEx Stacked Bar Plots: PDF(s)                               ##
#################################################################
## Parameters set as variables exported by plotVolc() before.  ##
#################################################################
## Use this function if you have modules in memory created     ##
## using the Seyfried Systems Biology Pipeline, after running  ##
## parANOVA.dex() and plotVolc() functions.                    ##
#################################################################

DEXpercentStacked()        # runs on prior function outputs as input; writes stacked bar plot(s) to PDF.



## End of examples using 3 core functions of this repository.
##############################################################################



########################################################################################################
## Correlation Alternative to ANOVA / DiffEx --
##
##   trait.corStat     -- create CORout dataframe of correlation p and R values which can be used in volcano
##                        and other steps of the pipeline
##
#########################################
## Required Loaded Data and Parameters ##
#########################################

rootdir="e:/SysBioPipeline/"
setwd(rootdir)
load("SeyfriedPipelineOutput-test.RData")


# Parameters (function may fallback to defaults if unspecified)
NETcolors=net$colors                           # list net with slot/vector containing module color assignments; length of vector must be equal to number of rows in cleanDat.
outputCSV=TRUE                                 # Output Full Table of statistics?  TRUE/FALSE
outFilePrefix="4"                              # typically "4", or step # in pipeline being run; for output file sorting by filename.
outFileSuffix="DeepADproteome"                 # A description of the project, used as a filename suffix

# These parameters are specific to trait correlation statistics generation; traits are provided as columns of the data frame stored in the provided example RData as the variable numericMeta. 
cor.traits=c("Astro","Endo","Micro","Neuro","Oligo")                 # Brain (DLPFC) tissue sample cell type proportion estimates from an abundance x sample matrix deconvolution algorithm
filter.trait="Center"                                                # Center at which the tissue samples were collected during autopsy
filter.trait.subsets=c("ALL","Emory","Hopkins","IU","UCSF","WashU")  # Subsets of case samples will be used for correlation to the cell type proportion estimates
                                                                     # (5 separate cor.traits x 6 sample subsets = 30 total p and R value columns to generate)
corFn="bicor"                                                        #'bicor' or anything else will cause Pearson (cor) to be used


source("./parANOVA.dex.R")
CORout <- trait.corStat()                      # runs on cleanDat and Grouping variables as required input.
# Correlation p + R table calculations complete. If you want to use the table with plotVolc(), set the variable corVolc=TRUE and use variable CORout to store the table generated.


#head(CORout)
dim(CORout)
#[1] 9379   63   #column count is for 2 unused columns, 30 p value columns, 30 R value (bicor) columns, and the network module color assignments as the last column





## Find protein rows with desired correlation stat profile(s) - add TRUE/FALSE filter for criteria met in first 2 (unused) columns of CORout data frame
## These boolean columns can be used as filters later. (optional)
colnames(CORout)[1:2]<-c("Filter.ALLsamples.NeuronCor.Sig",
                         "Filter.OnlyOneCenterSamples.NeuronCor.Sig")

# Which proteins are signigicantly correlated with neuronal proportions across all samples, and NOT significantly correlated with any of the other 4 cell type proportions?
this.column=colnames(CORout)[1]
CORout[,this.column]<-FALSE
CORout[which(CORout[,"p Neuro.ALL"]<0.01 & CORout[,"p Astro.ALL"]>0.05 & CORout[,"p Endo.ALL"]>0.05 & CORout[,"p Micro.ALL"]>0.05 & CORout[,"p Oligo.ALL"]>0.05), this.column] <- TRUE
# Which proteins are significantly correlated with neuronal propotions only in samples from a single center?
this.column=colnames(CORout)[2]
CORout[,this.column]<-FALSE
Neuron.ALLsample.ColumnIndex<-which(grepl("p Neuro.ALL",colnames(CORout)))
Neuron.SeparateCenter.ColumnIndices<-which(grepl("p Neuro.",colnames(CORout)) & !grepl("p Neuro.ALL",colnames(CORout)))
CORout[unlist(sapply(1:nrow(CORout),function(x) { row=CORout[x,]; if(row[Neuron.ALLsample.ColumnIndex]>0.05 & length(which(row[Neuron.SeparateCenter.ColumnIndices]<0.01))==1) x; } )), this.column] <- TRUE

table(CORout[,1])
#FALSE  TRUE 
# 9364    15

table(CORout[,2])
#FALSE  TRUE 
# 9129   250


corVolc=TRUE        # changes the behavior of plotVolc, DEXpercentStacked, and GOparallel functions later in the pipeline, to use CORout
useNETcolors=TRUE
sameScale=TRUE
highlightGeneProducts=rownames(CORout)[which(CORout$Filter.ALLsamples.NeuronCor.Sig)]  # Specifies which spots should be large
labelHighlighted=TRUE
plotVolc()          # Plots PDFs and HTMLs with large spots for proteins significant (p<0.01) for neuron proportion correlation in all samples,
                    # and not significant (p>0.05) for all the other 4 cell types



DEXpercentStacked(CORout)
#- corVolc=TRUE so plotting volcanoes using correlation statistics in the table stored in variable CORout.
#- Found plotVolc function output variable testIndexMasterList to recall comparisons selected for selectComps. Using the following comparisons:
#...
#- Plotting using correlation volcano(es) stats. Variable flip will be ignored so positive correlations remain positive.
#- Plotting significant correlation hit counts. Variable FCmin will be ignored.
#- Variable outputfigs not specified. Saving volcano plots to e:/OneDrive/SystemsBioPipeline .
#
#Performing calculations for mean differential abundance and fractions of differential gene products per module:
#- Generating PDF file(s)...


## Gene Ongology Gene Set Enrichment Analysis of significantly trait-correlated proteins:
source("./GOparallel-FET.R")  #Available from the repository file https://github.com/edammer/GOparallel/blob/main/GOparallel-FET.R
ANOVAgroups=TRUE
parallelThreads=30

GOparallel(CORout)
#- ANOVAout not in memory, trying to use input provided to this function (could be ANOVAout or CORout).
#- corVolc=TRUE so getting gene lists from significant correlations in statistics table stored in variable CORout.
#- Significant groups were defined by trait corrrelation. Variable flip will be ignored so positive correlations remain positive.
#- Extracting names of groups being compared in ANOVAout table from p value column headers...
#
#[1] "...Applying a minimum p value cutoff of 0.05 for correlation statistics lists..."
#- Using correlation significance for defining groups. Variable FCmin will be ignored.
#
#Running FET overlap statistics in parallel for  60  symbol lists using up to  30  threads...
