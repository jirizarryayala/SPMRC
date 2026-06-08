###RNA Seq Data Analysis & Viz w/ no replicates
###author: jirizarry
###date: 5/19/2026
##File paths
mondrinos_dir="Censored/File/Path/" #change with your own!
file_name="PQSPTX-expression-matrix.tsv"

##Libraries
if(!require("data.table")){
  install.packages("data.table")
}
if(!require("edgeR")){
  BiocManager::install("edgeR")
}
if(!require("tidyverse")){
  install.packages("tidyverse")
}
#remotes::install_version("systemfonts", version = "< 1.3.0") Newest version of systemfonts was breaking package installation
#remotes::install_version("ggplot2", version = "< 4.0.0") New ggplot also breaks this
if (!require("clusterProfiler")){
  BiocManager::install("clusterProfiler")
}
if (!require("org.Hs.eg.db")){
  BiocManager::install("org.Hs.eg.db")
}
if (!require("enrichplot")){
  BiocManager::install("enrichplot")
}
if (!require("ComplexHeatmap")){
  BiocManager::install("ComplexHeatmap")
}
if (!require("circlize")){
  library(circlize)
}

##Define DGE analysis function

DGE_Analysis=function(count_matrix, sample_groups, gene_annotation, bcv=0.2){  #the parameter bcv manually sets dispersion since we can't calculate it without replicates
  DGE_List_1=DGEList(counts=count_matrix, group=sample_groups, genes=gene_annotation)
  expressedDGE=filterByExpr(DGE_List_1) #Find unexpressed genes
  DGE_List_2=DGE_List_1[expressedDGE,,keep.lib.size=FALSE] #Remove unexpressed genes
  DGE_List_3=normLibSizes(DGE_List_2)
  fit=glmFit(DGE_List_3, dispersion=bcv^2) #Fit model
  return(fit)
}

##Read in expression data
expression_matrix=fread(paste(mondrinos_dir, file_name, sep=""))

names(expression_matrix); dim(expression_matrix)
# [1] "gene_id"         "gene_name"       "gene_biotype"    "PQSPTX_10_cpm"
# [5] "PQSPTX_10_count" "PQSPTX_11_cpm"   "PQSPTX_11_count" "PQSPTX_12_cpm"
# [9] "PQSPTX_12_count" "PQSPTX_13_cpm"   "PQSPTX_13_count" "PQSPTX_14_cpm"
# [13] "PQSPTX_14_count" "PQSPTX_15_cpm"   "PQSPTX_15_count" "PQSPTX_16_cpm"
# [17] "PQSPTX_16_count" "PQSPTX_1_cpm"    "PQSPTX_1_count"  "PQSPTX_2_cpm"
# [21] "PQSPTX_2_count"  "PQSPTX_3_cpm"    "PQSPTX_3_count"  "PQSPTX_4_cpm"
# [25] "PQSPTX_4_count"  "PQSPTX_5_cpm"    "PQSPTX_5_count"  "PQSPTX_6_cpm"
# [29] "PQSPTX_6_count"  "PQSPTX_7_cpm"    "PQSPTX_7_count"  "PQSPTX_8_cpm"
# [33] "PQSPTX_8_count"  "PQSPTX_9_cpm"    "PQSPTX_9_count"
# [1] 78986    35


##Grouping structure as defined by expert
# So, for the HRMVEC, it would be comparing samples 
# 1-3 (effect of hormones in male cells), 
# comparing samples 4-6 (effect of hormones in female cells), 
# then comparing sample 1 vs. 4 (male vs. female). 
# Same for HOF, compare 7-9, 10-12, and 7 vs. 10.
##
#Sample numbers not zero-padded on the left; need to keep that in mind
#Grouping structures are similar, we can make general ones and apply them as necessary
g1=c("A", "B", "C")
g2=c("A", "B")


#Keep only raw counts
keep_list=grep("count", names(expression_matrix))


expression_matrix_1=as.matrix(expression_matrix[, ..keep_list])
names(expression_matrix[, ..keep_list])
# [1] "PQSPTX_10_count" "PQSPTX_11_count" "PQSPTX_12_count" "PQSPTX_13_count"
# [5] "PQSPTX_14_count" "PQSPTX_15_count" "PQSPTX_16_count" "PQSPTX_1_count"
# [9] "PQSPTX_2_count"  "PQSPTX_3_count"  "PQSPTX_4_count"  "PQSPTX_5_count"
# [13] "PQSPTX_6_count"  "PQSPTX_7_count"  "PQSPTX_8_count"  "PQSPTX_9_count"
#head(expression_matrix_1)

gene_annotation=expression_matrix$gene_id
gene_names=expression_matrix$gene_name

##Comp 1: 1-3 (effect of hormones in male cells)
##Compile DGE List
expression_matrix_1_comp1=expression_matrix_1[,8:10]
comp_1_res=DGE_Analysis(expression_matrix_1_comp1, g1, gene_annotation)
comp_1=glmLRT(comp_1_res, coef=2:3)
topTags(comp_1) #Nothing of note

##Comp 2: 4-6 (effect of hormones in female cells)
expression_matrix_1_comp2=expression_matrix_1[,11:13]
comp_2_res=DGE_Analysis(expression_matrix_1_comp2, g1, gene_annotation)
comp_2=glmLRT(comp_2_res, coef=2:3)
topTags(comp_2) #Nothing of note

##Comp 3: 1 vs. 4 
expression_matrix_1_comp3=expression_matrix_1[,c(8,11)]
comp_3_res=DGE_Analysis(expression_matrix_1_comp3, g2, gene_annotation)
comp_3=glmLRT(comp_3_res)
topTags(comp_3) #Nothing of note

##Comp 4:  7-9
expression_matrix_1_comp4=expression_matrix_1[,14:16]
comp_4_res=DGE_Analysis(expression_matrix_1_comp4, g1, gene_annotation)
comp_4=glmLRT(comp_4_res, coef=2:3)
topTags(comp_4) #Stuff!

res_tab_4=data.frame(topTags(comp_4))
res_tab_4=res_tab_4[res_tab_4$FDR<.05,]
head(res_tab_4)

go_enrich_4 = enrichGO(gene = res_tab_4$genes, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'BP')

png(filename = "~/dotplot_4.png", type="cairo")
dotplot(go_enrich_4)
dev.off()

##Comp 5:  10-12
expression_matrix_1_comp5=expression_matrix_1[,1:3]
comp_5_res=DGE_Analysis(expression_matrix_1_comp5, g1, gene_annotation)
comp_5=glmLRT(comp_5_res, coef=2:3)
topTags(comp_5) #Nothing of note

##Comp 6: 7 vs. 12
expression_matrix_1_comp6=expression_matrix_1[,c(1,14)]
comp_6_res=DGE_Analysis(expression_matrix_1_comp6, g2, gene_annotation)
comp_6=glmLRT(comp_6_res)
topTags(comp_6) # Lots of stuff!

res_tab_6=data.frame(topTags(comp_6, n=Inf))
res_tab_6=res_tab_6[res_tab_6$FDR<.05,]
head(res_tab_6)

go_enrich_6 = enrichGO(gene = res_tab_6$genes, OrgDb = 'org.Hs.eg.db', keyType = 'ENSEMBL', ont = 'BP')
go_enrich_6

##Dotplot
png(filename = "~/dotplot_6.png", type="cairo")
dotplot(go_enrich_6)
dev.off()

#Select out nonesense pathways
go_enrich_6_f=go_enrich_6
hits_bool = grep("pregnancy|embryo", go_enrich_6@result$Description, ignore.case = TRUE)
go_enrich_6_f@result = go_enrich_6@result[-hits_bool, ]

png(filename = "~/dotplot_6_2.png", type="cairo", height = 600)
dotplot(go_enrich_6_f, showCategory=Inf)
dev.off()


#Heatmap
res_tab_6_1=res_tab_6[res_tab_6$logFC>-7,]
sig_genes=res_tab_6_1$genes
sig_gene_counts=expression_matrix_1_comp6[gene_annotation %in% sig_genes,]

gene_names[gene_names==""]=gene_annotation[gene_names==""]


row.names(sig_gene_counts)=gene_names[gene_annotation %in% sig_genes]

sig_gene_counts=sig_gene_counts[!grepl(pattern = "ENSG000",row.names(sig_gene_counts)),] #Remove noncoding RNAs

sig_gene_counts=as.matrix(sig_gene_counts)
sig_gene_counts_1=log2(sig_gene_counts+0.00000001) # Correct 0 coutn
colnames(sig_gene_counts_1)=c("Female", "Male")


png(filename = "~/heatmap_comp_6_final.png", type="cairo", height= 300, width = 600)
  hm=Heatmap(t(scale(sig_gene_counts_1)), 
          heatmap_legend_param=list(title=""),
          cluster_rows=F, show_heatmap_legend = FALSE)
  col_fun=colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  lgd=Legend(col_fun = col_fun, title = "")
  draw(hm)
  draw(lgd, x = unit(20, "cm"), y = unit(6.1, "cm"))
dev.off()

##Barplots

sig_gene_counts_sub=sig_gene_counts[grepl(pattern = "RAMP2|PAPPA|MMP3|KDR|ITGB3",row.names(sig_gene_counts)),] #Select expert-indicated genes
sig_gene_counts_sub=sig_gene_counts_sub[!grepl(pattern = "AS1",row.names(sig_gene_counts_sub)),] #Remove non-coding antisense
sig_gene_counts_sub_log2=log2(sig_gene_counts_sub)

#Try log2 transformed and raw.
#log2 first
png(filename = "~/count_barplot_g6_final_v1.png", type="cairo", width = 600, height = 300)
par(mar = c(4, 4.1, 2, 1))  
barplot(t(sig_gene_counts_sub_log2), beside = TRUE, col = c("red", "blue"),
        names.arg = row.names(sig_gene_counts_sub_log2), las=3, ylab=expression("Log2 Transcript Count") )
legend("topright", 
       legend = c("Female", "Male"), 
       fill = c("red", "blue"))
dev.off()

#Now raw
png(filename = "~/count_barplot_g6_final_v2.png", type="cairo", width = 600, height = 300)
par(mar = c(4, 4.1, 2, 1))  
barplot(t(sig_gene_counts_sub), beside = TRUE, col = c("red", "blue"),
        names.arg = row.names(sig_gene_counts_sub), las=3, ylab=expression("Transcript Count") )
legend("topright", 
       legend = c("Female", "Male"), 
       fill = c("red", "blue"))
dev.off()

library(edgeR)
packageVersion("edgeR")

#Original
# png(filename = "~/count_barplot_g6.png", type="cairo", width = 1200, height = 720)
#   par(mar = c(10, 4, 4, 2))  
#   barplot(t(sig_gene_counts), beside = TRUE, col = c("red", "blue"),
#         names.arg = row.names(sig_gene_counts), las=3)
#   legend("topright", 
#          legend = c("Female", "Male"), 
#          fill = c("red", "blue"))
# dev.off()

##Per pathway barplots
gene_ids=go_enrich_6_f@result$geneID[go_enrich_6_f@result$p.adjust<.05]
path_descs=go_enrich_6_f@result$Description[go_enrich_6_f@result$p.adjust<.05]


gene_id_lists=c()
for(i in 1:length(gene_ids)){
  gene_id_lists[i]=str_split(gene_ids[i], pattern="/")
}

for(i in 1:length(gene_id_lists)){
  sig_gene_counts=expression_matrix_1_comp6[gene_annotation %in% gene_id_lists[[i]],]
  row.names(sig_gene_counts)=gene_names[gene_annotation %in% gene_id_lists[[i]]]
  sig_gene_counts=as.matrix(sig_gene_counts)
  png(filename = paste("~/count_barplot_g6", as.character(i),".png", sep=""), type="cairo", width = 1200, height = 720)
  par(mar = c(10, 4, 4, 2))  
  barplot(t(sig_gene_counts), beside = TRUE, col = c("red", "blue"),
          names.arg = row.names(sig_gene_counts), las=3, main = path_descs[i])
  legend("topright", 
         legend = c("PQSPTX_10", "PQSPTX_7"), 
         fill = c("red", "blue"))
  dev.off()
}
