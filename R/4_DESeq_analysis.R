##DESeq analysis
library(DESeq2)
library(phyloseq)
library(microbiome)
library(tidyverse)
library(data.table)
library(scales)

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

##Differences by phyla

Bac.diagdds <- phyloseq_to_deseq2(subset_taxa(PS.raw.genus, Kingdom%in%"Bacteria"), ~ AnimalSpecies)
Bac.diagdds <- DESeq(Bac.diagdds, test="LRT", fitType="parametric", reduced= ~ 1)

Bac.res <- results(Bac.diagdds, cooksCutoff = FALSE)
alpha <- 0.05
Bac.sigtab <- Bac.res[which(Bac.res$padj < alpha), ]
Bac.sigtab <- cbind(as(Bac.sigtab, "data.frame"), as(tax_table(PS.genus)[rownames(Bac.sigtab), ], "matrix"))
rownames(Bac.sigtab) <- NULL
Bac.sigtab

get.sigtab.plot <- function (sigtab){
  sigtab$Phylum <- factor(as.character(sigtab$Phylum), levels=main_phyla)
  ## Genus order
  x <- tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
  x <- sort(x, TRUE)
  sigtab$Genus <- factor(as.character(sigtab$Genus), levels=names(x))
  ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum, label=scientific(padj))) +
    geom_point(size=6) +
    geom_text(vjust=1.6, color="black")+
    theme_bw()+
    scale_color_brewer(palette = "Set1") +
    theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust=1))
}

plot.bacterial.diff <- get.sigtab.plot(Bac.sigtab)

devSVG("figures/figures_Hyena/Figure3d_single_diff.svg")
plot.bacterial.diff
dev.off()

pdf("figures/figures_Hyena/Figure3d_single_diff.pdf")
plot.bacterial.diff
dev.off()
