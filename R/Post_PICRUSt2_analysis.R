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
library(ALDEx2)
library(tidyverse)
library(viridis)

##Re run
reRun<- FALSE
##Load data
##1) Sample data
if(!exists("sample")){
  if(isTRUE(reRun)){
    source("Dada2_Pipeline.R") ## Run the script at base directory of repository!   
  } else {
    sample <- read.csv("/SAN/Victors_playground/Ascaris_Microbiome/Pig_Ascaris_16S_Samples_barcode.csv", dec=",", stringsAsFactors=FALSE)
    ##Add sample names used by Ulrike
    sample[ , "MDC_Names"] <- samplesMDC ##
  }
}
##Define factor variables
fac.vars <- c("Barcode_Plate", "Barcode_Well", "Sample_ID", 
  "Barcode_ID", "Barcode_Seq","Compartment", "System", "DPI", "InfectionStatus",
  "AnimalSpecies", "WormSex", "Live", "MDC_Names") 
 
sample[, fac.vars] <- apply(sample[, fac.vars], 2, as.factor)

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
rownames(sample) <- paste0("Sample", sample$BeGenDiv_Name)

keep <- data.frame(name = colnames(PredDesM))
coldata<- sample[keep$name,]
rownames(coldata) <- paste0("Sample", coldata$BeGenDiv_Name)
rm(keep)
##Create DESeq table 
ddsTable <- DESeqDataSetFromMatrix(
  countData = round(PredDesM),
  colData = coldata,
  design = ~Compartment)

##Select the top 25 more abundant genes 
PredDes%>%
  replace(is.na(.), 0)%>%
  mutate(Total = rowSums(select(., contains("Sample"))))%>%
  select(function., description, Total)%>%
  column_to_rownames(var = "function.")-> tmp

tmp <- tmp[order(-tmp$Total), ]
tmp <- rownames(tmp[1:25, ])
PredDestop25 <- data.frame(PredDesM[tmp, ])

##Modify data frame to canonical format
mPredDestop25 <- data.frame(Gene = rownames(PredDestop25), PredDestop25)
mPredDestop25 <- melt(mPredDestop25, id.vars = "Gene")
colnames(mPredDestop25) <- c("Gene", "Sample_Name", "Counts")
##Add metadata 
coldata%>%
  rownames_to_column(var= "Sample_Name")->coldata
mPredDestop25 <- plyr::join(mPredDestop25, coldata, by="Sample_Name")
mPredDestop25[, fac.vars] <- apply(mPredDestop25[, fac.vars], 2, as.factor)

##Check genes more samples
## plot using ggplot2
mPredDestop25%>%
ggplot(aes(x = reorder(Gene, -Counts), y = Counts)) +
  geom_point() +
  scale_y_log10(name = "log10 Gene Count", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  annotation_logticks(sides = "l")+
  labs(tag= "A)")+
  xlab("Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

mPredDestop25%>%
  filter(Compartment%in%c("Ascaris"))%>%
  ggplot(aes(y= Counts, x= reorder(Gene, -Counts)))+
  scale_y_log10("log10 Gene counts", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot(aes(color= WormSex))+
  xlab("Enzyme Classification ID")+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))+
  annotation_logticks(sides = "l")

mPredDestop25%>%
  filter(Compartment%in%c("Faeces", "Duodenum","Jejunum", "Colon", "Cecum", "Ileum"))%>%
  ggplot(aes(y= Counts, x= reorder(Gene, -Counts)))+
  scale_y_log10("log10 Gene counts", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot(aes(color= Compartment))+
  xlab("Enzyme Classification ID")+
  scale_color_brewer(palette = "Set3")+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))+
  annotation_logticks(sides = "l")

mPredDestop25%>%
  filter(Compartment%in%c("Faeces"))%>%
  ggplot(aes(y= Counts, x= reorder(Gene, -Counts)))+
  scale_y_log10("log10 Gene counts", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot(aes(color= DPI))+
  xlab("Enzyme Classification ID")+
  labs(tag= "D)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))+
  annotation_logticks(sides = "l")

mPredDestop25%>%
  filter(Compartment%in%c("Faeces"))%>%
  filter(Gene%in%c("EC:1.6.5.3"))%>%
  ggplot(aes(y= Counts, x= as.factor(DPI)))+
  scale_y_log10("log10 Gene counts", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot(color= "black")+
  geom_jitter(shape=21, position=position_jitter(0.2), size=3, aes(fill= DPI), color= "black")+
  xlab("DPI")+
  labs(tag= "D)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  annotation_logticks(sides = "l")+
  stat_compare_means(label= "p.signif", method = "t.test", ref.group = "2", paired = F, na.rm = TRUE)+
  stat_compare_means(method =  "anova", label.y = 2.5, label.x = 1)

mPredDestop25%>%
  filter(AnimalSpecies%in%c("Pig", "Ascaris"))%>%
  ggplot(aes(y= Counts, x= as.factor(AnimalSpecies)))+
  scale_y_log10("log10 Gene counts", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot(aes(color= AnimalSpecies))+
  xlab("Animal species")+
  labs(tag= "D)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))+
  annotation_logticks(sides = "l")

mPredDestop25%>%
  filter(AnimalSpecies%in%c("Pig"))%>%
  ggplot(aes(y= Counts, x= as.factor(System)))+
  scale_y_log10("log10 Gene counts", 
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  geom_boxplot(aes(color= System))+
  xlab("Pig individual")+
  labs(tag= "D)")+
  theme_bw()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))+
  annotation_logticks(sides = "l")

######### Bray-Curtis dissimilarity estimation ########
heatmap(as.matrix(PredDestop25))

PredDestop25.matrix<- as.matrix(PredDestop25)
PredDestop25.matrix<- t(PredDestop25.matrix)
PredDestop25.braycurt<- vegdist(PredDestop25.matrix, method = "bray")
PredDestop25.braycurt<-as.matrix(PredDestop25.braycurt)

###Using pheatmap to include annotations 
PredDestop25.clust <- hclust(dist(PredDestop25.braycurt), method = "complete") ##Dendogram

require(dendextend)
as.dendrogram(PredDestop25.clust) %>%
  plot(horiz = TRUE)

PredDestop25.col <- cutree(tree = PredDestop25.clust, k = 2)
PredDestop25.col  <- data.frame(cluster = ifelse(test = PredDestop25.col  == 1, yes = "cluster 1", no = "cluster 2"))
PredDestop25.col$Sample_Name <- rownames(PredDestop25.col)

PredDestop25.col <- merge(PredDestop25.col, coldata, by="Sample_Name", sort= F)

col_groups <- PredDestop25.col %>%
  select("Sample_Name", "AnimalSpecies", "Compartment") ##Here It is possible to add the expected size 

row.names(col_groups)<- col_groups$Sample_Name

col_groups$Sample_Name<- NULL

colour_groups <- list(AnimalSpecies= c("Pig"= "pink", "Ascaris"= "#C46210"),
                      Region= c("Faeces"= "#8DD3C7","Colon"= "#009999" ,"Duodenum"= "#FFFFB3", "Jejunum"= "#BEBADA", 
                                "Negative"= "#FB8072", "Ascaris"= "#80B1D3"))

BCPredDestop25 <- pheatmap(PredDestop25.braycurt, 
                      color = viridis(100),
                      border_color = NA,
                      annotation_col = col_groups, 
                     #annotation_row = col_groups,
                      annotation_colors = colour_groups,
                      #cutree_rows = 2,
                      #cutree_cols = 2,
                      show_rownames = F,
                      show_colnames = F,
                      main= "Bray-Curtis dissimilarity among samples")

#pdf(file = "/SAN/Victors_playground/Ascaris_Microbiome/output/BC_Predicted_genes.pdf", width = 10, height = 8)
#BCPredDestop25
#dev.off()

##DESeq analysis
ddsTable<- DESeq(ddsTable)
