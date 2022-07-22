########################################################################################################
## ANOVA / DiffEx -- list of ANOVAout dataframes for different subgroup comparisons, if necessary
##
## By Eric Dammer and Duc Duong
#########################################


#########################################
## Required Loaded Data and Parameters ##
#########################################


#parANOVA.dex <- function(cleanDat.=cleanDat, ..., env=.GlobalEnv) {     # cleanDat, data table of 'cleaned' relative abundance with rows (genes or proteins) and columns (samples)
parANOVA.dex <- function(cleanDat.=cleanDat,                             # cleanDat, data table of 'cleaned' relative abundance with rows (genes or proteins) and columns (samples)
                         Grouping=numericMeta$Group,                     # Named groups (N>=2) for comparison of difference of means, in sample (column) order of cleanDat
                         parallelThreads=30,                             # number of CPU threads to speed calculation (recommended, 2 or more):
                         NETcolors=net$colors,                           #list net with slot/vector containing module color assignments; length of vector must be equal to number of rows in cleanDat.
                         twoGroupCorrMethod="BH",                        # default method for full FDR correction when only 2 groups present is Benjamini-Hochberg
                         outputCSV=TRUE,                                 #Output Full Table of statistics?  TRUE/FALSE
                         outFilePrefix="4",                              # typically "4", or step # in pipeline being run
                         outFileSuffix=FileBaseName,                     # A description of the project, used as a filename suffix
                         fallbackIfZeroTukeyP=TRUE, env=.GlobalEnv) {

############################
## END PARAMETERS SECTION ##
############################
  if(!exists("cleanDat.")) stop("cleanDat must be supplied, in form of log2(relative abundance) data with rows as genes/proteins and columns as samples.")

  ## Set up parallel backend as a local cluster using n (parallelThreads) CPU threads (in the .GlobalEnv outside of this function!)
#  require(doParallel, quietly=TRUE)
  require(parallel, quietly=TRUE)
  #if(exists("clusterLocal")) stopCluster(clusterLocal) #in case already started.
#  clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
#  registerDoParallel(clusterLocal)

  if (nchar(outFilePrefix)>0) outFilePrefix=paste0(outFilePrefix,".")
  if (nchar(outFileSuffix)>0) outFileSuffix=paste0("-",outFileSuffix)

  ANOVAoutList<-list()
  caseSubset="ALL"

  data = as.data.frame(cbind(as.character(colnames(cleanDat.)), as.character(Grouping), t(cleanDat.)))
  colnames(data)[1:2]<-c("CODE","SampleType")
  #test run gets column headers for output
  i=3
  aov<-aov(as.double(rnorm(n=nrow(data)))~SampleType, data=data)
  anovaresult<-anova(aov)
  tuk <- TukeyHSD(aov)
  tukresult1<-data.frame(tuk$SampleType) #With random (rnorm fn output), NO PAIRWISE COMPARISONS ARE MISSING--this is the template for comparison columns in ANOVAoutList[[caseSubset]]
  j=length(rownames(tukresult1))
  comparisonList<-rownames(tukresult1)

  line = c(paste("Protein", "F-Value", "Pr(>F)", sep=","))
  for (a in 1:length(comparisonList)) {
    line=c(paste(line,comparisonList[a],sep=","))
  }
  for (a in 1:length(comparisonList)) {
    line=c(paste(line,paste0("diff ",comparisonList[a]),sep=","))
  }

  ## Fast code with parApply
  dataclipped <- data[, 3:ncol(data)] # as.matrix(as.numeric(data[,3:ncol(data)]))
  SampleType <- as.character(data$SampleType)
  parallel::clusterExport(cl=clusterLocal, list("SampleType","tukresult1","fallbackIfZeroTukeyP"), envir=environment())  # below helperfn not running in .GlobalEnv
  parts <- splitIndices(ncol(dataclipped), length(clusterLocal))
  dataclippedParts <- lapply(parts, function(i) dataclipped[,i,drop=FALSE])
#  resParts <- parallel::parApply(cl=clusterLocal, dataclippedParts, 2, function(x) {
  resParts <- parallel::clusterApply(cl=clusterLocal, dataclippedParts, helperfn)
  ANOVAoutList[[caseSubset]] <- do.call(cbind,resParts)
  ANOVAoutList[[caseSubset]] <- t(ANOVAoutList[[caseSubset]])
  countZeroFallbacks=sum(ANOVAoutList[[caseSubset]][,ncol(ANOVAoutList[[caseSubset]])])
  print(paste0("...Tukey p=0 Fallback calculations using Bonferroni corrected T test: ", countZeroFallbacks, " [",signif(countZeroFallbacks/(nrow(tukresult1)*nrow(ANOVAoutList[[caseSubset]]))*100,2),"%]"))
  ANOVAoutList[[caseSubset]] <- ANOVAoutList[[caseSubset]][,-ncol(ANOVAoutList[[caseSubset]])]
  ANOVAcols <- as.vector(data.frame(do.call("rbind", strsplit(as.character(line), "[,]"))))
  ANOVAcols <- ANOVAcols[2:length(ANOVAcols)]
  if (length(unique(SampleType))==2) { #for single pairwise comparison, essentially Pr(>F) from ANOVA is equivalent to T-test result (except 1 vs 2 non-NA measurements are allowed)
    ANOVAoutList[[caseSubset]][,3] <- p.adjust(ANOVAoutList[[caseSubset]][,2],method=twoGroupCorrMethod) #get BH FDR for ANOVA/T-test (equivalent) p values
    ANOVAoutList[[caseSubset]][,2:3]<-ANOVAoutList[[caseSubset]][,c(3,2)]
    ANOVAcols[2] <- paste0("FDR (",twoGroupCorrMethod,")")
  }
  colnames(ANOVAoutList[[caseSubset]]) <- ANOVAcols
  ANOVAoutList[[caseSubset]] <- as.data.frame(ANOVAoutList[[caseSubset]])

  numComp=(ncol(ANOVAoutList[[caseSubset]])-2)/2

  ## Add network module colors
  if(exists("NETcolors")) if (length(NETcolors)==nrow(ANOVAoutList[[caseSubset]])) { ANOVAoutList[[caseSubset]]$NETcolors <- NETcolors } else { warning("Network color assignment vector not of length in rows of cleanDat supplied -- not included in output table.") }

  #*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
  if(outputCSV)  write.csv(ANOVAoutList[[caseSubset]],file=paste0(outFilePrefix,"ANOVA_diffEx-",caseSubset,outFileSuffix,".csv"))
  #*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

#  stopCluster(clusterLocal)
  return(ANOVAoutList[[caseSubset]])
} # end function parANOVA.dex


## Per solution to spur parallel use from within function (worked with parApply above only outside of function):  https://stackoverflow.com/questions/43269142/r-parlapply-not-parallel/45447830#45447830
helperfn <- function(xx) {
    apply(xx,2,function(x) {
    if(length(unique(SampleType[which(!is.na(x))]))<2) { #<2, handles "only one contrast level..." 
      x=data.frame(x=rep(1,length(x)), SampleType = SampleType)  #handles too many missing values
    } else {
      x <- data.frame(x = as.double(x), SampleType = SampleType) #as.double(x) instead of x corrects:   Error in lm.fit(x, y,... NA/NaN/Inf in 'y'
    }
    aov <- aov(as.double(x) ~ SampleType, data = x)
    anovaresult <- anova(aov)
    tuk <- TukeyHSD(aov)
    tukresult <- data.frame(tuk$SampleType)
    if (!length(rownames(tukresult)) == length(rownames(tukresult1))) tukresult <- tukresult[match(rownames(tukresult1), rownames(tukresult)), ]
    if (fallbackIfZeroTukeyP) {
      # Fallback to Bonferroni correction of t.test output if Tukey p returned as 0:
      zeroFallbacks=length(which(tukresult[,"p.adj"]==0))
      if (zeroFallbacks>0) for (comp in which(tukresult[,"p.adj"]==0)) {
        this.comp=rownames(tukresult)[comp]
        grp1=gsub("^(.*)\\-.*","\\1",this.comp)
        grp2=gsub("^.*\\-(.*)$","\\1",this.comp)
        lowTuk.p.estimate <- p.adjust(t.test(x[SampleType==grp1,"x"],x[SampleType==grp2,"x"],alternative="two.sided",var.equal=FALSE)$p.value, method="bonferroni",n=nrow(tukresult))
        tukresult[comp,"p.adj"] <- lowTuk.p.estimate
      }
    }

    c(anovaresult$F[1], anovaresult$Pr[1], as.vector(tukresult[, "p.adj"]), as.vector(tukresult[, "diff"]), zeroFallbacks)

    })  #close apply()

#  }, chunk.size = round(ncol(dataclipped)/parallelThreads,0)) #end parApply
  } # end function helperfn






##############################
## Volcano Plotter function ##
##############################
## 
## requires ggplot2, plotly packages.
##

plotVolc<- function(ANOVAout=ANOVAout,
                    FCmin=0,                     # 0.25 for 25%, 0 for no threshold (vertical minimum FC threshold dashed lines)
                    selectComps="ALL",           # "ALL" for volcano output(s) on all pairwise comparisons in ANOVAout
                    flip=c(),                    # p value column numbers in which to swap denominator of pair for x axis range (gene products high in denominator, will be on left)
                    signifP=0.05,                # p value threshold for counting Differential Expression points
                    useNETcolors=TRUE,           # use module colors saved to ANOVAout, if available; otherwise specify downColor upColor, and NCcolor
                    downColor="royalblue",       # significant points above/beyond thresholds on the upper left are this color if useNETcolors=FALSE
                    upColor="red",               # significant points above/beyond thresholds on the upper right are this color if useNETcolors=FALSE
                    NCcolor="grey",              # points not significant are this color if useNETcolors=FALSE
                    splitColors=FALSE,           # create a separate volcano plot(s) for each color in an outputfigs/splitVolcanoes subfolder (folder created if it does not exist)
                    highlightGeneProducts=c(),   # c("APP|P05067","MAPT|P10636","APOE|P02649") ; a list of uniqueID rownames to highlight as larger gold points. If symbolsOnly=TRUE, this can be a list of symbols, like c("APP","SMOC1","MAPT")
                    symbolsOnly=FALSE,           # for mouse-over HTML plots and the above highlight callouts, consider only displaying and using official gene symbol from first part of UniqueID rownames of ANOVAout.
                    HTMLout=TRUE,                # output interactive HTML copies that can be opened in browser. Requires plotly package.
                    outFilePrefix="4",           # typically "4", or step # in pipeline being run
                    outFileSuffix=FileBaseName,  # A description of the project, used as a filename suffix
                    outputfigs=getwd(),          # "drive:/folder/to/location" to save output PDFs and possible HTML files
                    env=.GlobalEnv) {

############################
## END PARAMETERS SECTION ##
############################

require(ggplot2,quietly=TRUE)

#baseNameVolcanoes <- paste0(FileBaseName) #** paste0(FileBaseName,".",caseSubset)

numberOfNonComparisonColumns=length(colnames(ANOVAout)) - length(which(grepl("diff ",colnames(ANOVAout))))*2
numComp <- (length(colnames(ANOVAout)) - numberOfNonComparisonColumns) / 2 # of columns separating comparisons from matched column of log2(diffs), i.e. # of comparisons

# What columns (column numbers) of ANOVAout (T-test or Tukey p values) do we plot volcanoes for? 
testIndexMasterList <- if (selectComps=="ALL" | selectComps=="all" | selectComps=="All") { c(3:(numComp+2)) } else { if(is.integer(as.numeric(selectComps))) { as.integer(selectComps) } else { error("selectComps must be set to 'all' or valid integer column indexes of ANOVAout.") }}
# log2(1) means NO change minimum to be counted in the volcano bookends; log2(1.25) for 25% FC min.
cutoff <- log2(1+FCmin)
# p value cutoff for Volcano significant hit counting; dashed line at -log10(sigCutoff)
sigCutoff <- signifP

# Gene products to highlight
BIGspots <- highlightGeneProducts

cutoff=log2(1+FCmin)
# shows what your cutoff for log2(FC) calculates as
print(paste0("...Applying a ", FCmin*100,"% minimum fold change threshold at + and - x=", signif(cutoff,2)," ..."))


n <- nrow(ANOVAout)
dexComps <- list()
iter <- length(testIndexMasterList) + 1
comparisonIDs <- data.frame(dfVariable = rep(NA, length(testIndexMasterList)), Comparison = rep(NA, length(testIndexMasterList)))
for (i in testIndexMasterList) {
  iter <- iter - 1
  # dexRows<-which(ANOVAout[,i]<sigCutoff) #choose rows where the DEX p<sigCutoff
  comparisonIDs[iter, ] <- as.vector(c(paste0("dexTargets.", gsub("-", ".", colnames(ANOVAout)[i])), paste0(as.character(gsub("-", " vs ", colnames(ANOVAout)[i])))))
  dexComps[[comparisonIDs[iter, 1]]] <- ANOVAout
  if (!is.na(match(i, flip))) {
    dexComps[[comparisonIDs[iter, 1]]][, i + numComp] <- -1 * as.numeric(dexComps[[comparisonIDs[iter, 1]]][, i + numComp])
    comparisonIDs[iter, 2] <- gsub("(*.*) vs (*.*)", "\\2 vs \\1", comparisonIDs[iter, 2]) # flip label "vs" in comParisonIDs$Comparison[iter]
  }
}
comparisonIDs # list element names and Logical comparisons for those retrievable Dex measurements in the list elements
ls(dexComps) # list elements are dataframes with the DEX entries for that comparison


## volcano plots with (module or other) colors, SPLITTABLE on COLOR

pointColorsVectorListForPlots <- list()
volcListModColorsWeb <- list()
volcListModColors <- list()
dfListModColors <- list()

iter <- length(testIndexMasterList) + 1
for (testIndex in testIndexMasterList) {
  iter <- iter - 1
  df <- eval(parse(text = "dexComps[[comparisonIDs$dfVariable[iter]]]"))
  print(paste0("Processing ANOVA column ", testIndex, " (", comparisonIDs$Comparison[iter], ") for volcano..."))

  # correct 0 Tukey pValues to ANOVA p (in column 2); it's better than taking -log10 of 0 in the next step
  if (length(which(df[, testIndex] == 0))>0) print(paste0("Found ",length(which(df[, testIndex] == 0))," p=0 values for this comparison! Substituted -log10(p) using one-way ANOVA overall p value (Consider fallback=TRUE in ANOVAout calculation.)"))
  df[which(df[, testIndex] == 0), testIndex] <- as.numeric(df[which(df[, testIndex] == 0), 2])
 
  ## Check if ANOVA pVal is Significant and above FC cutoff defined above. Thresholds are used to set volcano point colors
  df[, testIndex][is.na(df[, testIndex])] <- 1 # p=0.9999 instead of NA
  df[, testIndex + numComp][is.na(df[, testIndex + numComp])] <- 0 # log2(difference)=0 instead of NA
  
  df$negLogP <- -log10(as.numeric(df[, testIndex]))
  
  df$threshold1 <- as.numeric(rep(0, n))
  ## Any COMPARISON SIGNIFICANT (uses ANOVA p in column 2 of df instead of Tukey p): # for (i in 1:n) { if (abs(as.numeric(df[i,testIndex+numComp]))<cutoff | df[i,2]>sigCutoff ) {df$threshold1[i]=3} else { if (df[i,testIndex+numComp]<cutoff) {df$threshold1[i]=2} else {df$threshold1[i]=1}} }
  for (i in 1:n) {
    if (abs(as.numeric(df[i, testIndex + numComp])) < cutoff | as.numeric(df[i, testIndex]) > sigCutoff) {
      df$threshold1[i] <- 3
    } else {
      if (as.numeric(df[i, testIndex + numComp]) < cutoff) {
        df$threshold1[i] <- 2
      } else {
        df$threshold1[i] <- 1
      }
    }
  }
  df$threshold1 <- as.factor(df$threshold1)
  
  df$Symbol <- rownames(df)
  if (symbolsOnly) df$Symbol <- as.data.frame(do.call("rbind", strsplit( as.data.frame(do.call("rbind", strsplit(as.character(rownames(df)), "[|]")))[, 1], "[;]")))[,1]
  
  ## Color Interesting Gene Product Spots DIFFERENTLY as 4th color if doing blue/red/green (no module colors) -- (4=gold1 below)
  if(useNETcolors) { 
    df$color1<-df$NETcolors
    df$threshold2 <- as.numeric(df$threshold1)
    df$threshold2[match(intersect(df$Symbol, BIGspots), df$Symbol)] <- 4 
  } else {
    df$color1 <- as.numeric(df$threshold1)
    df$color1[match(intersect(df$Symbol, BIGspots), df$Symbol)] <- 4
    df$threshold2 <- as.numeric(df$threshold1)
    df$threshold2[match(intersect(df$Symbol, BIGspots), df$Symbol)] <- 4 
  }
  df$color1 <- as.factor(df$color1)

  if(useNETcolors) {
    df$size1 <- as.numeric(df$threshold2)
  } else {
    df$size1 <- as.numeric(df$color1)
  }
  df$size1[df$size1 < 4] <- 2.5
  df$size1[df$size1 == 4] <- 6

  #df$color2: actual spot color as.character() ; df$color1 is a factorable dummy
  if(useNETcolors) {
    df$color2<-as.character(df$color1)
    df$color2[df$threshold2 == 4] <- "gold1" #for BIGspots
  } else {
    df$color2 <- as.numeric(df$color1)
    df$color2[df$color2 == 1] <- upColor
    df$color2[df$color2 == 2] <- downColor
    df$color2[df$color2 == 3] <- NCcolor
    df$color2[df$color2 == 4] <- "gold1" #for BIGspots
  }                              #gold1 is also the 435th, last WGCNA unique color

  #df$color3 is outline color, where outlined pch symbols (21) used
  df$color3 <- df$color2
  df$color3[df$color3 == "gold1"] <- "black" #for BIGspots outline
  df$pch <- as.numeric(df$threshold2)
  df$pch[df$pch < 4] <- 16 # unfilled circles (use color2)
  df$pch[df$pch == 4] <- 21 # filled, outlined circles (border uses color3)

  #put gold1 back to module color for fill of BIGspots if (useNETcolors)
  if (useNETcolors) { df$color2[df$color2=="gold1"] <- as.character(df$color1[df$color2=="gold1"]) }

  df <- df[order(df$size1, decreasing = FALSE), ] # puts larger dots on top (at bottom of df)
  
  
  #splitColors TRUE/FALSE: make one volcano with all colors (FALSE), or make volcanoes for each color (TRUE)
  # SPLIT DATA FRAME FOR VOLCANO PLOT BY COLORS (if multiple eachColorSplit items)
  df.AllColors <- df
  eachColorSplit <- if (splitColors) {
    unique(df.AllColors$NETcolors)
  } else {
    c("allcolors")
  }
  for (eachColor in eachColorSplit) {
    if (splitColors) {
      df.oneColor <- df.AllColors[which(df.AllColors$NETcolors == eachColor), ]
      df.oneColor$color1 <- factor(df.oneColor$NETcolors)
    }
    else {
      df.oneColor <- df.AllColors
    } # end if (splitColors)
    
    names(df.oneColor)[testIndex + numComp] <- "xdata" # x=as.numeric(df.oneColor[,testIndex+numComp])
    list_element <- paste(comparisonIDs$dfVariable[iter], eachColor, sep = ".") # colnames(df)[testIndex]
    pointColorsVectorListForPlots[[list_element]] <- data.frame(color1 = factor(as.integer(df.oneColor$color1)), color2 = as.character(df.oneColor$color2), color3 = as.character(df.oneColor$color3), size = as.numeric(df.oneColor$size), pch = as.numeric(df.oneColor$pch)) #*** df.oneColor$NETcolors
    
    
    volcano1 <- ggplot(data = df.oneColor, aes(x = xdata, y = negLogP, color = color1, text = Symbol)) +
      # scale_colour_manual(values = unique(data.frame(col1=df.oneColor$color1,col2=df.oneColor$color2))[order(unique(data.frame(col1=df.oneColor$color1,col2=df.oneColor$color2))[,1]),2] )+ #THIS COLOR(S) IS LOOKED UP ACTIVELY DURING GRAPHICS OUTPUT IN THE VARIABLE, SO WE'VE USED A LIST ELEMENT THAT IS NEVER CHANGED
      geom_point(aes(fill = pointColorsVectorListForPlots[[list_element]][, "color3"]), alpha = 0.66, size = pointColorsVectorListForPlots[[list_element]]$size, pch = pointColorsVectorListForPlots[[list_element]][, "pch"], color = pointColorsVectorListForPlots[[list_element]][, "color3"]) +
      theme(legend.position = "none") +
      xlim(c(min(as.numeric(df.oneColor[, testIndex + numComp])), max(as.numeric(df.oneColor[, testIndex + numComp])))) + ylim(c(0, max(df.oneColor$negLogP))) +
      xlab(as.expression(bquote("Difference, log"[2] ~ .(comparisonIDs$Comparison[iter])))) + # colnames(df.oneColor)[testIndex]
      ylab(as.expression(bquote("-log"[10] ~ "p value"))) +
      theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
      theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
      
      geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black", size = 1.2) +
      # geom_text(aes(0,1.30103,label = 1.30103, vjust = -1))+
      geom_vline(xintercept = cutoff, linetype = "dashed", color = "black", size = 1.2) +
      geom_vline(xintercept = -cutoff, linetype = "dashed", color = "black", size = 1.2) +
      annotate("text", x = min(as.numeric(df.oneColor[, testIndex + numComp])) / 2, y = max(df.oneColor$negLogP) * .95, size = 5, label = paste0("Downregulated: ", bquote(.(length(which(as.numeric(df.oneColor$threshold1) == 2)))))) +
      annotate("text", x = max(as.numeric(df.oneColor[, testIndex + numComp])) / 2, y = max(df.oneColor$negLogP) * .95, size = 5, label = paste0("Upregulated: ", bquote(.(length(which(as.numeric(df.oneColor$threshold1) == 1)))))) +
      
      theme(
        # axis.text = element_text(size = 14),
        # legend.key = element_rect(fill = "navy"),
        # legend.background = element_rect(fill = "white"),
        # legend.position = c(0.14, 0.80),
        panel.grid.major = element_line(color = "darkgrey", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white")
      )
    
    
    # web version doesn't use as.expression! plotly fails with those, so we rebuild the volcano for the web.
    volcanoweb <- ggplot(data = df.oneColor, aes(x = xdata, y = negLogP, color = color1, text = Symbol)) +
      scale_colour_manual(values = unique(data.frame(col1 = df.oneColor$color1, col2 = df.oneColor$color2))[order(unique(data.frame(col1 = df.oneColor$color1, col2 = df.oneColor$color2))[, 1]), 2]) + # THIS COLOR(S) IS LOOKED UP ACTIVELY BY PLOTLY IN THE VARIABLE, SO WE'VE USED A LIST ELEMENT THAT IS NEVER CHANGED
      # scale_y_continuous(breaks = seq(0, 8, by = 1))+
      # geom_point(aes(fill=pointColorsVectorListForPlots[[list_element]][,"color2"]), alpha=0.66, size=pointColorsVectorListForPlots[[list_element]]$size, pch=16, color=pointColorsVectorListForPlots[[list_element]][,"color3"]) + #pch=pointColorsVectorListForPlots[[list_element]][,"pch"] not using variable symbol types (21 for outlined circle only)
      geom_point(alpha = 0.66, size = pointColorsVectorListForPlots[[list_element]]$size, pch = 16) + # pch=pointColorsVectorListForPlots[[list_element]][,"pch"] just uses the higher pch code in the web render.
      theme(legend.position = "none") +
      xlim(c(min(as.numeric(df.oneColor[, testIndex + numComp])), max(as.numeric(df.oneColor[, testIndex + numComp])))) + ylim(c(0, max(df.oneColor$negLogP))) +
      xlab(paste0("Difference, log2 ", comparisonIDs$Comparison[iter])) +
      ylab(paste0("-log10 p value")) +
      theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
      theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
      
      geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black", size = 1.2) +
      # geom_text(aes(0,1.30103,label = 1.30103, vjust = -1))+
      geom_vline(xintercept = cutoff, linetype = "dashed", color = "black", size = 1.2) +
      geom_vline(xintercept = -cutoff, linetype = "dashed", color = "black", size = 1.2) +
      annotate("text", x = min(as.numeric(df.oneColor[, testIndex + numComp])) / 2, y = max(df.oneColor$negLogP) * .95, size = 5, label = paste0("Downregulated: ", bquote(.(length(which(as.numeric(df.oneColor$threshold1) == 2)))))) +
      annotate("text", x = max(as.numeric(df.oneColor[, testIndex + numComp])) / 2, y = max(df.oneColor$negLogP) * .95, size = 5, label = paste0("Upregulated: ", bquote(.(length(which(as.numeric(df.oneColor$threshold1) == 1)))))) +
      
      theme(
        # axis.text = element_text(size = 14),
        # legend.key = element_rect(fill = "navy"),
        # legend.background = element_rect(fill = "white"),
        # legend.position = c(0.14, 0.80),
        panel.grid.major = element_line(color = "darkgrey", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white")
      )
    
    
    volcListModColors[[list_element]] <- volcano1
    volcListModColorsWeb[[list_element]] <- volcanoweb
    
    print(volcano1) # prints to active output (separate page)
    rm(volcano1)
    rm(volcanoweb)
    dfListModColors[[list_element]] <- df.oneColor
  } # closes for(eachColor...
} # closes for(testIndex...


## Write Files
if(splitColors) { dir.create(file.path(outputfigs, "/SplitVolcano/")) }

if (nchar(outFilePrefix)>0) outFilePrefix=paste0(outFilePrefix,".")
if (nchar(outFileSuffix)>0) outFileSuffix=paste0("-",outFileSuffix)


# Print to PDFs, one per color (per comparison, if multiple)
iter <- length(testIndexMasterList) + 1
for (testIndex in testIndexMasterList) {
  iter <- iter - 1
  for (eachColor in eachColorSplit) {
    list_element <- paste(comparisonIDs$dfVariable[iter], eachColor, sep = ".")
    df.oneColor <- dfListModColors[[list_element]]
    if(splitColors) {
      file <- paste0(outputfigs,"/SplitVolcano/",outFilePrefix,"Volcano_", eachColor,"-",gsub(" ", "_", comparisonIDs$Comparison[iter]),outFileSuffix,".pdf")
    } else {
      file <- paste0(outputfigs,"/",outFilePrefix,"Volcano_", eachColor,"-",gsub(" ", "_", comparisonIDs$Comparison[iter]),"-",outFileSuffix,".pdf")
    }
    pdf(file = file, colormodel="srgb", height = 8, width = 8)
    par(mfrow = c(1, 1))
    par(mar = c(6, 8.5, 3, 3))
    print(volcListModColorsWeb[[list_element]])
    #     volcListModColors[[list_element]])    #bigspots colors not handled correctly for PDF output of modColored spots:
    dev.off()
  }
}

## Print to html volcano plots (one per module colors and per comparison, if applicable)
if (HTMLout) {
require(plotly, quietly=TRUE)

iter <- length(testIndexMasterList) + 1
for (testIndex in testIndexMasterList) {
  iter <- iter - 1
  for (eachColor in eachColorSplit) {
    list_element <- paste(comparisonIDs$dfVariable[iter], eachColor, sep = ".")
    df.oneColor <- dfListModColors[[list_element]]
    webPlot <- ggplotly(volcListModColorsWeb[[list_element]])
    if(splitColors) {
      tempfilename=paste0(outputfigs,"/SplitVolcano/",outFilePrefix,"HTMLvolcano_",eachColor,"-",gsub(" ", "_", comparisonIDs$Comparison[iter]),"-",outFileSuffix,".html")
    } else {
      tempfilename=paste0(outputfigs,"/",outFilePrefix,"HTMLvolcano_",eachColor,"-",gsub(" ", "_", comparisonIDs$Comparison[iter]),"-",outFileSuffix,".html")
    }
    htmlwidgets::saveWidget(webPlot, tempfilename, selfcontained = TRUE, libdir = "delete.me")
  }
}

} #end if (HTMLout)

}
