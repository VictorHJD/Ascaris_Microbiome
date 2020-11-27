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


##Scaling abundances function
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
##Load data
##1) Sample data
if(!exists("sample")){
    sample <- read.csv("/SAN/Victors_playground/Ascaris_Microbiome/Pig_Ascaris_16S_Samples_barcode.csv", dec=",", stringsAsFactors=FALSE)
    ##Add sample names used by Ulrike
    path <- "/SAN/Victors_playground/Ascaris_Microbiome/2018_22_Nem1/" # CHANGE ME to the directory containing the fastq files after unzipping.
    list.files(path)
    fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
    fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
    fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 
    samplesMDC<- gsub("-\\d+_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
    sample[ , "MDC_Names"] <- samplesMDC ##
}
rm(fastqF, fastqFiles, fastqR, path, samplesMDC)

##Define factor variables
fac.vars <- c("Barcode_Plate", "Barcode_Well", "Sample_ID", 
  "Barcode_ID", "Barcode_Seq","Compartment", "System", "DPI", "InfectionStatus",
  "AnimalSpecies", "WormSex", "Live", "MDC_Names") 
 
sample[, fac.vars] <- apply(sample[, fac.vars], 2, as.factor)

##2)Predicted metagenomes
setwd("/SAN/Victors_playground/Ascaris_Microbiome/output/picrust2_out_pipeline/")
##How the ASVs contribute to gene family abundances in each sample
PredMet<- read.table("EC_metagenome_out/pred_metagenome_contrib.tsv", header = T, sep = "\t")
##Enzyme classification (EC) abundances per sample (Count table)
PredDes<- read.table("EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv", header = T, sep = "\t")

##3)How the ASVs contribute to KEGG onthology (KO) abundances in each sample
PredKO<- read.table("KO_metagenome_out/pred_metagenome_contrib.tsv", header = T, sep = "\t")
##KO abundances per sample (Count table)
PredKODes<- read.table("KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz", header = T, sep = "\t")

##4)Predicted pathway functions
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
pheatmap(as.matrix(PredDestop25))

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

##Add annotation for samples
col_groups <- coldata %>%
  select("Sample_Name", "AnimalSpecies", "Compartment", "Barcode_Plate") ##Here It is possible to add the expected size 

row.names(col_groups)<- col_groups$Sample_Name

col_groups$Sample_Name<- NULL

colour_groups <- list(AnimalSpecies= c("Pig"= "pink", "Ascaris"= "#C46210"),
                      Region= c("Faeces"= "#8DD3C7","Colon"= "#009999" ,"Duodenum"= "#FFFFB3", "Jejunum"= "#BEBADA", 
                                "Negative"= "#FB8072", "Ascaris"= "#80B1D3"),
                      Barcode_Plate= c(A1= "red", A2= "blue"))

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

####Heat map pathways#####
##Matrix Pathways sample
PredPathDesM<- PredPathDes
rownames(PredPathDesM)<-PredPathDesM$pathway
PredPathDesM$pathway<- NULL
PredPathDesM$description<- NULL
PredPathDesM%>%
  replace(is.na(.), 0)-> PredPathDesM
  
PredPathDesM<- as.matrix(PredPathDesM)
pheatmap(PredPathDesM)

PredPathDesM.norm <- t(apply(PredPathDesM, 1, cal_z_score))
as.data.frame(PredPathDesM.norm)%>%
  replace(is.na(.), 0)-> PredPathDesM.norm
PredPathDesM.norm<- as.matrix(PredPathDesM.norm)
pheatmap(PredPathDesM.norm)

PredPathDesM.clust <- hclust(dist(PredPathDesM), method = "complete") ##Dendogram

require(dendextend)
as.dendrogram(PredPathDesM.clust) %>%
  plot(horiz = TRUE)
##Detect cluster of pathways (apparently 2)
PredPathDesM.col <- cutree(tree = PredPathDesM.clust, k = 2)
PredPathDesM.col  <- data.frame(cluster = ifelse(test = PredPathDesM.col  == 1, yes = "cluster 1", no = "cluster 2"))
PredPathDesM.col%>%
  filter(cluster== "cluster 1")-> PathC1

PathC1<- rownames(PathC1)
PathC1<-gsub("-", ".", PathC1)

##Select the top 25 more abundant pathways out of cluster 1 pathways 
PredPathDes%>%
  replace(is.na(.), 0)%>%
  mutate(Total = rowSums(select(., contains("Sample"))))%>%
  select(pathway, description, Total)%>%
  column_to_rownames(var = "pathway")-> tmp

tmp <- tmp[order(-tmp$Total), ]
tmp <- rownames(tmp[1:25, ])
PredPathDestop25 <- data.frame(PredPathDesM[tmp, ])

rownames(PredPathDes)<-PredPathDes$pathway
Pathdes<- PredPathDes[tmp,]
Pathdes<- as.data.frame(Pathdes[,1:2])
#write.csv(Pathdes, "/SAN/Victors_playground/Ascaris_Microbiome/output/Top_predicted_pathways.csv")

df<- data.frame(PredPathDestop25)
df<- t(df)
df<- data.frame(df)
df%>% 
  mutate(Total= rowSums(select(., contains("PWY"))))%>%
  mutate_all(.funs = ~./Total)->df

rownames(df)<-rownames(t(data.frame(PredPathDestop25)))
df$Total<- NULL
df<- t(df)

mPredPathDestop25 <- data.frame(Pathway = rownames(df), df)
mPredPathDestop25 <- melt(mPredPathDestop25, id.vars = "Pathway")
colnames(mPredPathDestop25) <- c("Pathway", "Sample_Name", "Relative_abundance")

##Add metadata 
mPredPathDestop25<- plyr::join(mPredPathDestop25, coldata, by="Sample_Name")
Path25<- unique(mPredPathDestop25$Pathway)
#pdf(file = "/SAN/Victors_playground/Ascaris_Microbiome/output/Top_predicted_pathways.pdf", width = 10, height = 8)
mPredPathDestop25%>%
  ggplot(aes(x = reorder(Pathway, -Relative_abundance), y = Relative_abundance, color= Pathway)) +
  scale_color_viridis_d()+
  geom_jitter(size = 2, alpha = 0.05, width = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 5)+
  coord_flip() +
  scale_y_continuous(name = "Pathway relative abundance", limits = c(0, 0.075))+
  labs(tag= "A)")+
  xlab("Pathway") +
  theme_bw()+
  theme(legend.position = "none")+
  geom_segment(
    aes(x = Pathway, xend = Pathway,
        y = 0.04, yend = mean(Relative_abundance)),
    size = 0.8) +
  geom_hline(aes(yintercept = 0.04), color = "gray70", size = 0.6)
#dev.off()

mPredPathDestop25%>%
  filter(Compartment%in%c("Ascaris"))%>%
  ggplot(aes(y= Relative_abundance, x= reorder(Pathway, -Relative_abundance)))+
  scale_y_continuous("Pathway relative abundance")+
  geom_boxplot(aes(color= WormSex))+
  xlab("Pathway")+
  labs(tag= "B)")+
  coord_flip() +
  theme_bw()+
  theme(text = element_text(size=16))

###Using pheatmap to include annotations 
df.matrix<-as.matrix(PredPathDestop25)

##Scaling abundances
df.matrix.norm <- t(apply(df.matrix, 1, cal_z_score))
pheatmap(df.matrix.norm)

df.clust <- hclust(dist(df.matrix), method = "complete") ##Dendogram

require(dendextend)
as.dendrogram(df.clust) %>%
  plot(horiz = TRUE)
##Detect cluster of pathways (apparently none :S)
df.col <- cutree(tree = df.clust, k = 2)
df.col  <- data.frame(cluster = ifelse(test = df.col  == 1, yes = "cluster 1", no = "cluster 2"))

##Add annotation for samples
col_groups <- coldata %>%
  select("Sample_Name", "AnimalSpecies", "Compartment", "Barcode_Plate") ##Here It is possible to add the expected size 

row.names(col_groups)<- col_groups$Sample_Name

col_groups$Sample_Name<- NULL

colour_groups <- list(AnimalSpecies= c("Pig"= "pink", "Ascaris"= "#C46210"),
                      Region= c("Faeces"= "#8DD3C7","Colon"= "#009999" ,"Duodenum"= "#FFFFB3", "Jejunum"= "#BEBADA", 
                                "Negative"= "#FB8072", "Ascaris"= "#80B1D3"),
                      Barcode_Plate= c(A1= "red", A2= "blue"))

HMPredPathDestop25 <- pheatmap(df.matrix.norm, 
                           #color = viridis(100),
                           border_color = NA,
                           annotation_col = col_groups, 
                           #annotation_row = col_groups,
                           annotation_colors = colour_groups,
                           #cutree_rows = 2,
                           cutree_cols = 2,
                           show_rownames = ,
                           show_colnames = F,
                           main= "Pathway relative abundance among samples")

#pdf(file = "/SAN/Victors_playground/Ascaris_Microbiome/output/Predicted_pathways.pdf", width = 10, height = 8)
#HMPredPathDestop25
#dev.off()

##Select jusr Ascaris samples to compare between sexes 
coldata%>%
  filter(AnimalSpecies == "Ascaris")-> Asc.df
keep<- Asc.df$Sample_Name

tmp<- as.data.frame(t(PredPathDestop25))
tmp$Sample_Name<- rownames(tmp)
tmp<- tmp[tmp$Sample_Name %in% keep,]
tmp$Sample_Name<- NULL
tmp<- t(tmp)
Asc.matrix<-as.matrix(tmp)

##Scaling abundances
Asc.matrix.norm <- t(apply(Asc.matrix, 1, cal_z_score))
##Add annotation for samples
Asc_groups <- Asc.df %>%
  select("Sample_Name", "Barcode_Plate", "System", "WormSex") ##Here It is possible to add the expected size 

row.names(Asc_groups)<- Asc_groups$Sample_Name

Asc_groups$Sample_Name<- NULL

colour_Asc_groups <- list(WormSex= c("Male"= "green", "Female"= "#C46210"),
                      System= c("Pig1"= "#8DD3C7","Pig2"= "#009999" ,"Pig3"= "#FFFFB3", "Pig5"= "#BEBADA", 
                                "SH"= "#FB8072"),
                      Barcode_Plate= c(A1= "red", A2= "blue"))

AscPredPathDestop25 <- pheatmap(Asc.matrix.norm, 
                               border_color = NA,
                               annotation_col = Asc_groups,
                               annotation_colors = colour_Asc_groups,
                               cutree_cols = 2,
                               show_rownames = T,
                               show_colnames = F,
                               main= "Pathway scaled abundance among Ascaris samples")

#pdf(file = "/SAN/Victors_playground/Ascaris_Microbiome/output/Predicted_pathways_Ascaris.pdf", width = 10, height = 8)
#AscPredPathDestop25
#dev.off()

##Select jusr Ascaris samples to compare between sexes 
coldata%>%
  filter(AnimalSpecies == "Pig")-> Pig.df
keep<- Pig.df$Sample_Name

tmp<- as.data.frame(t(PredPathDestop25))
tmp$Sample_Name<- rownames(tmp)
tmp<- tmp[tmp$Sample_Name %in% keep,]
tmp$Sample_Name<- NULL
tmp<- t(tmp)
Pig.matrix<-as.matrix(tmp)

##Scaling abundances
Pig.matrix.norm <- t(apply(Pig.matrix, 1, cal_z_score))
##Add annotation for samples
Pig_groups <- Pig.df %>%
  select("Sample_Name", "Barcode_Plate", "System", "InfectionStatus") ##Here It is possible to add the expected size 

row.names(Pig_groups)<- Pig_groups$Sample_Name

Pig_groups$Sample_Name<- NULL

colour_Pig_groups <- list(InfectionStatus= c("Infected"= "green", "Non_infected"= "#C46210"),
                          System= c("Pig1"= "#8DD3C7","Pig2"= "#009999" ,"Pig3"= "#FFFFB3", "Pig4"= "#E6AB02","Pig5"= "#BEBADA", 
                                     "Pig6"= "#80B1D3", "Pig7"= "#FDB462", "Pig8"= "#B3DE69", "Pig9"= "#FC4E07"),
                          Barcode_Plate= c(A1= "red", A2= "blue"))

PigPredPathDestop25 <- pheatmap(Pig.matrix.norm, 
                                border_color = NA,
                                annotation_col = Pig_groups,
                                annotation_colors = colour_Pig_groups,
                                cutree_cols = 2,
                                show_rownames = T,
                                show_colnames = F,
                                main= "Pathway scaled abundance among Pig samples")

#pdf(file = "/SAN/Victors_playground/Ascaris_Microbiome/output/Predicted_pathways_Pig.pdf", width = 10, height = 8)
#PigPredPathDestop25
#dev.off()
##DESeq analysis
ddsTable<- DESeq(ddsTable)
