###Sample Differential Expression Analysis
###Author: jirizarry
###Date: 10/11/2024

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

##Read in data
DEG=fread("/path/to/sample")
sampleGroups=factor(c("wt", "wt", "s", "s", "wt", "s")) #Define grouping structure here,
#exmaple if samples alternate between male and female group=c("male", "female", "male", "female")
DEG1=DGEList(DEG, group = sampleGroups)
expressedDEG=filterByExpr(DEG1) #Find unexpressed genes
#table(expressedDEG) can show how many genes are expressed or not
DEG2=DEG1[expressedDEG,,keep.lib.size=FALSE] #Remove unexpressed genes
DEG3=normLibSizes(DEG2)
modelDesign=model.matrix(~sampleGroups)
fit = glmQLFit(DEG3, modelDesign) #Fit inverse binomial model
qlf=glmQLFTest(fit) #Significance test
hits=topTags(qlf) #Selects top hits
hits$table$gene #Returns top hits
