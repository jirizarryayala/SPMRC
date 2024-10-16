###Sample Differential Expression Analysis
###Author: jirizarry
###Date: 10/11/2024

###Note: This code is mean to be run on Tulane University's High Performance Computing cluster Cypress using R version 4.3.2
###This requires use of the centos7 partition!
###Example: idev --time=0-05:00:00 --partition=centos7

##Load packages, will install themselves if necessary
if(!require("data.table")){
  install.packages("data.table")
}
if(!require("BiocManager")){
  install.packages("BiocManager")
}
if(!require("edgeR")){
  BiocManager::install("edgeR")
}

##Example with simple comparison
DGE=fread("/path/to/sample") #Read in data
sampleGroups=factor(c("wt", "wt", "s", "s", "wt", "s")) #Define grouping structure here,
#example if samples alternate between male and female group=c("male", "female", "male", "female")
DGE1=DGEList(DGE, group = sampleGroups)
expressedDGE=filterByExpr(DGE1) #Find unexpressed genes
#table(expressedDGE) can show how many genes are expressed or not
DGE2=DGE1[expressedDGE,,keep.lib.size=FALSE] #Remove unexpressed genes
DGE3=normLibSizes(DGE2)
DGE3=estimateDisp(DGE3)
modelDesign=model.matrix(~sampleGroups)
fit = glmQLFit(DGE3, modelDesign) #Fit inverse binomial model
qlf=glmQLFTest(fit) #Significance test
hits=topTags(qlf) #Selects top hits
hits$table$Gene #Returns top hits

##Example with multiple comparisons
DGE=fread("/path/to/sample")
sampleGroups=factor(c("wt", "s1", "s2", "s1", "wt", "s2"))
DGE1=DGEList(DGE, group = sampleGroups)
expressedDGE=filterByExpr(DGE1)
DGE2=DGE1[expressedDGE,,keep.lib.size=FALSE]
DGE3=normLibSizes(DGE2)
DGE3=estimateDisp(DGE3)
modelDesign=model.matrix(~0+sampleGroups, data=DGE3$samples)
colnames(modelDesign)=levels(DGE3$samples$group)
fit = glmQLFit(DGE3, modelDesign)
contrasts=makeContrasts(A=s1-wt, B=s2-wt, C=s1-s2,levels=modelDesign)#Here we define what comparisons we are able to do later on
qlf.A=glmQLFTest(fit,contrast=contrasts[,"A"])
resultsA=topTags(qlf.A)
qlf.B=glmQLFTest(fit,contrast=contrasts[,"B"])
resultsB=topTags(qlf.B)
qlf.C=glmQLFTest(fit,contrast=contrasts[,"C"])
