###Analysis of predicted metagenomes generated in PICRUSt2
library(lme4)
library(lmtest)
library(ggplot2)
library(reshape2)
library(vegan)
library(gtools)
library(ggpubr)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)

##Re run
reRun<- FALSE
##Load data
##1) Sample data
if(!exists("sample")){
  if(isTRUE(reRun)){
    source("Dada2_Pipeline.R") ## Run the script at base directory of repository!   
  } else {
    sample<- readRDS(file = "/SAN/Victors_playground/Ascaris_Microbiome/output/sample.Rds") ##
  }
}

##2)Predicted metagenomes
setwd("/SAN/Victors_playground/Ascaris_Microbiome/output/picrust2_out_pipeline/")
##How the ASVs contribute to gene family abundances in each sample
PredMet<- read.table("EC_metagenome_out/pred_metagenome_contrib.tsv", header = T, sep = "\t")
##Enzyme classification abundances per sample (Count table)
PredDes<- read.table("EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv", header = T, sep = "\t")

##3)Predicted pathway functions
##How the ASVs contribute to metabolic pathways abundances in each sample
PredPath<- read.table("pathways_out/path_abun_contrib.tsv", header = T, sep = "\t")
##Metabolic pathway abundances per sample
PredPathDes<- read.table("pathways_out/path_abun_unstrat_descrip.tsv", header = T, sep = "\t")

##Adjust tables
##Matrix EC sample
PredDesM<- PredDes
rownames(PredDesM)<-PredDesM$function.
PredDesM$function.<- NULL
PredDesM$description<- NULL
PredDesM<- as.matrix(PredDesM)
coldata<- as.tibble(sample)
rownames(coldata) <- paste0("Sample", coldata$BeGenDiv_Name)

keep <- data.frame(name = colnames(PredDesM))
coldata<- coldata[keep$name,]
rownames(coldata) <- paste0("Sample", coldata$BeGenDiv_Name)

ddsTable <- DESeqDataSetFromMatrix(
  countData = round(PredDesM),
  colData = coldata,
  design = ~Compartment)

ddsTable<- DESeq(ddsTable)

results(ddsTable, cooksCutoff = FALSE)
##Select the top 25 more abundant genes 
## Order results by padj values
PredDes%>%
  replace(is.na(.), 0)%>%
  mutate(Total = rowSums(select(., contains("Sample"))))%>%
  select(function., description, Total)-> tmp

PredMet_ordered <- PredMet[order(PredMet$), ]

top20_sigOE_genes <- rownames(res_tableOE_ordered[1:20, ])


##Heatmap
heatmap(PredDesM, scale = "row")

######### Bray-Curtis dissimilarity estimation (Figure 2) ########
foo.matrix<- as.matrix(foo)
foo.braycurt<- vegdist(foo.matrix, method = "bray")
as.matrix(foo.braycurt)

###Using pheatmap to include annotations 
foo.clust <- hclust(dist(foo.braycurt), method = "complete") ##Dendogram
require(dendextend)
as.dendrogram(foo.clust) %>%
  plot(horiz = TRUE)

foo.col <- cutree(tree = foo.clust, k = 2)
foo.col  <- data.frame(cluster = ifelse(test = foo.col  == 1, yes = "cluster 1", no = "cluster 2"))
foo.col$Primer_comb_ID <- rownames(foo.col)

foo.col <- merge(foo.col, primerInput, by="Primer_comb_ID", sort= F)

col_groups <- foo.col %>%
  select("Primer_comb_ID", "Gen", "Region") ##Here It is possible to add the expected size 

row.names(col_groups)<- col_groups$Primer_comb_ID

col_groups$Primer_comb_ID<- NULL

colour_groups <- list( Gen= c("12S"= "#E3DAC9", "16S"= "pink","18S"= "#440154FF", "28S"= "#21908CFF", 
                              "COI"= "#FDE725FF", "ITS"= "#C46210", "rbcL"= "#D0FF14"),
                       Region= c("V1-V2"= "#8DD3C7","V1-V3"= "#009999" ,"V3-V4"= "#FFFFB3", "V4"= "#BEBADA", "V4-V5"= "#FB8072", "V6-V7"= "#80B1D3",
                                 "V6-V8"="#FDB462", "V7-V8"= "#B3DE69", "V7-V9"="#FC4E07","V8-V9"= "#FCCDE5", "V9"= "#D9D9D9",
                                 "D2"="#D95F02", "D3"= "#7570B3", "Folmer"= "#E7298A", "ITS1"= "#A6761D", "Chloroplast"= "#66A61E", "Mitochondrial"= "#E6AB02"))
require(pheatmap)
require(viridis)
BCheatmap <- pheatmap(foo.braycurt, 
                      color = plasma(100),
                      border_color = NA,
                      annotation_col = col_groups, 
                      #annotation_row = col_groups,
                      annotation_colors = colour_groups,
                      #cutree_rows = 2,
                      #cutree_cols = 2,
                      show_rownames = F,
                      show_colnames = F,
                      main= "Bray-Curtis dissimilarity among primers")

#pdf(file = "~/AA_Primer_evaluation/Figures/Manuscript/Figure_2.pdf", width = 10, height = 8)
#BCheatmap
#dev.off()
