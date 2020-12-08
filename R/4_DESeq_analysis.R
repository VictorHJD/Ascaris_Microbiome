##DESeq analysis
library(DESeq2)
library(phyloseq)
library(microbiome)
library(tidyverse)
library(data.table)
library(scales)
library(gridExtra)
library(grid)

reRun <- FALSE
##Load data (not rarefied or cleaned)
if(!exists("PS")){
  if(isTRUE(reRun)){
    source("1_Dada2_Pipeline.R") ## Run the script at base directory of repository!   
  } else {
    PS<- readRDS(file = "/SAN/Victors_playground/Ascaris_Microbiome/output/PhyloSeqComp.Rds") ##New annotation SILVA
    sample<- data.table(as(sample_data(PS), "data.frame"), keep.rownames = T)
    ##Eliminate "empty" samples 
    PS <- prune_samples(sample_sums(PS)>0, PS)
  }
}

##Glom to genus level 
PS.raw.genus <- tax_glom(PS, "Genus")
##Eliminate "Negative" controls an faeces from the analysis 
PS.raw.genus<- subset_samples(PS.raw.genus, !(Compartment%in%c("Faeces", "Negative")))
##Select just Jejunum and Ascaris from infected individuals 
PS.PA.genus<- subset_samples(PS.raw.genus, !(Compartment%in%c("Duodenum", "Colon", "Cecum", "Ileum")))
PS.PA.genus<- subset_samples(PS.PA.genus, !(System%in%c("SH", "Pig6", "Pig7", "Pig8", "Pig9")))
##Select just Ascaris samples
PS.Asc.genus<- subset_samples(PS.raw.genus, (Compartment%in%c("Ascaris")))
##Select just Jejunum and Duodenum from Infected
PS.PJD.genus<- subset_samples(PS.raw.genus, !(Compartment%in%c("Ascaris", "Colon", "Cecum", "Ileum")))
PS.PJD.genus<- subset_samples(PS.PJD.genus, !(System%in%c("Pig6", "Pig7", "Pig8", "Pig9")))

### Comparison Pig Jejunum vs Ascaris from Jejunum
##Differences by phyla
Bac.diagdds <- phyloseq_to_deseq2(subset_taxa(PS.PA.genus, Kingdom%in%"Bacteria"), ~ AnimalSpecies)
Bac.diagdds <- DESeq(Bac.diagdds, test="LRT", fitType="parametric", reduced= ~ 1)

Bac.res <- results(Bac.diagdds, cooksCutoff = FALSE)
alpha <- 0.05
Bac.sigtab <- Bac.res[which(Bac.res$padj < alpha), ]
Bac.sigtab <- cbind(as(Bac.sigtab, "data.frame"), as(tax_table(PS.PA.genus)[rownames(Bac.sigtab), ], "matrix"))
rownames(Bac.sigtab) <- NULL
Bac.sigtab

# Phylum order
x <- tapply(Bac.sigtab$log2FoldChange, Bac.sigtab$Phylum, function(x) max(x))
x <- sort(x, TRUE)
Bac.sigtab$Phylum <- factor(as.character(Bac.sigtab$Phylum), levels=names(x))
# Genus order
x <- tapply(Bac.sigtab$log2FoldChange, Bac.sigtab$Genus, function(x) max(x))
x <- sort(x, TRUE)
Bac.sigtab$Genus <- factor(as.character(Bac.sigtab$Genus), levels=names(x))

ggplot(Bac.sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= Phylum), color= "black") + 
  scale_fill_brewer(palette = "Set1")+
  coord_flip()+
  geom_hline(aes(yintercept = 0), color = "gray70", size = 0.6)+
  xlab("Bacteria Genus")+
  ylab("Ascaris <-- Log-2-Fold-Change --> Jejunum")+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5), text = element_text(size=16))+
  theme_bw()-> A
  
png("Figures/Q6_Diffferential_Genus_Jejunum_Worms.png", units = 'in', res = 300, width=14, height=14)
grid.arrange(A)
dev.off()

### Comparison Ascaris Pig vs SH
##Differences by phyla
Bac.diagdds <- phyloseq_to_deseq2(subset_taxa(PS.Asc.genus, Kingdom%in%"Bacteria"), ~ Live)
Bac.diagdds <- DESeq(Bac.diagdds, test="LRT", fitType="parametric", reduced= ~ 1)

Bac.res <- results(Bac.diagdds, cooksCutoff = FALSE)
alpha <- 0.05
Bac.sigtab <- Bac.res[which(Bac.res$padj < alpha), ]
Bac.sigtab <- cbind(as(Bac.sigtab, "data.frame"), as(tax_table(PS.Asc.genus)[rownames(Bac.sigtab), ], "matrix"))
rownames(Bac.sigtab) <- NULL
Bac.sigtab

# Phylum order
x <- tapply(Bac.sigtab$log2FoldChange, Bac.sigtab$Phylum, function(x) max(x))
x <- sort(x, TRUE)
Bac.sigtab$Phylum <- factor(as.character(Bac.sigtab$Phylum), levels=names(x))
# Genus order
x <- tapply(Bac.sigtab$log2FoldChange, Bac.sigtab$Genus, function(x) max(x))
x <- sort(x, TRUE)
Bac.sigtab$Genus <- factor(as.character(Bac.sigtab$Genus), levels=names(x))

ggplot(Bac.sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) +
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= Phylum), color= "black") + 
  scale_fill_brewer(palette = "Set1")+
  coord_flip()+
  geom_hline(aes(yintercept = 0), color = "gray70", size = 0.6)+
  xlab("Bacteria Genus")+
  ylab("Ascaris(SH) <-- Log-2-Fold-Change --> Ascaris(FU)")+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5), text = element_text(size=16))+
  theme_bw()-> B

png("Figures/Q6_Diffferential_Genus_FU_SH.png", units = 'in', res = 300, width=14, height=14)
grid.arrange(B)
dev.off()

### Comparison Ascaris Pig Sex
##Differences by phyla
Bac.diagdds <- phyloseq_to_deseq2(subset_taxa(PS.Asc.genus, Kingdom%in%"Bacteria"), ~ WormSex)
Bac.diagdds <- DESeq(Bac.diagdds, test="LRT", fitType="parametric", reduced= ~ 1)

Bac.res <- results(Bac.diagdds, cooksCutoff = FALSE)
alpha <- 0.05
Bac.sigtab <- Bac.res[which(Bac.res$padj < alpha), ]
Bac.sigtab <- cbind(as(Bac.sigtab, "data.frame"), as(tax_table(PS.Asc.genus)[rownames(Bac.sigtab), ], "matrix"))
rownames(Bac.sigtab) <- NULL
Bac.sigtab

##No Genus is significative "expressed" in either male or female worms 

### Comparison Jejunum Pig vs Duodenum (Infected Pigs)
##Differences by phyla
Bac.diagdds <- phyloseq_to_deseq2(subset_taxa(PS.PJD.genus, Kingdom%in%"Bacteria"), ~ Compartment)
Bac.diagdds <- DESeq(Bac.diagdds, test="LRT", fitType="parametric", reduced= ~ 1)

Bac.res <- results(Bac.diagdds, cooksCutoff = FALSE)
alpha <- 0.05
Bac.sigtab <- Bac.res[which(Bac.res$padj < alpha), ]
Bac.sigtab <- cbind(as(Bac.sigtab, "data.frame"), as(tax_table(PS.PJD.genus)[rownames(Bac.sigtab), ], "matrix"))
rownames(Bac.sigtab) <- NULL
Bac.sigtab

##Just a genus from Succinivibrionaceae was significative over represented in duodenum 