library(phyloseq)
library(microbiome)
library(tidyverse)

reRun <- FALSE

##Load data 
if(!exists("PS")){
    if(isTRUE(reRun)){
        source("Dada2_Pipeline.R") ## Run the script at base directory of repository!   
    } else {
        #PS<- readRDS(file="/SAN/Victors_playground/Ascaris_Microbiome/PhyloSeqCombi.Rds") #Old annotation BLAST
      PS<- readRDS(file = "/SAN/Victors_playground/Ascaris_Microbiome/output/PhyloSeqComp.Rds") ##New annotation SILVA
    }
}

### Inspect
summarize_phyloseq(PS)
rank_names(PS)

## HOW many ASVs for off-target eukaryotes and archaea
table(tax_table(PS)[, "Kingdom"], exclude = NULL) ## ---> README results summary

## HOW many reads for off-target eukaryotes and archaea
by(t(otu_table(PS)), tax_table(PS)[, "Kingdom"], sum) ## --->  README results summary

### which different phyla for each sample
## barplot.PS <- plot_bar(PS, fill = "Phylum") ##Without normalization

##Rarefaction curves TODO: COLOR BY sample type (best make rather a ggplot)
rarcurv <- vegan::rarecurve(otu_table(PS),
                            label = T, xlab = "Read number", ylab = "ASVs")

pdf("Figures/Rarefaction_curves.pdf", 
    width=8, height=10, onefile=T)
rarcurv 
dev.off()

##Eliminate "empty" samples 
PS <- prune_samples(sample_sums(PS)>0, PS)

#Visualize alpha-diversity (raw data)
alphaDiv.PS <- plot_richness(PS, x= "Compartment", color = "Compartment" , measures = c("Observed","Chao1", "Shannon")) +
  #geom_boxplot()+
  geom_jitter(alpha= 0.005)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))

png("Figures/Alphadiv_PCoA.png", units = 'in', res = 300, width=4, height=3)
alphaDiv.PS
dev.off()

##Beta diversity (raw data)
ord.bray <- ordinate(PS, method="PCoA", distance="bray")

betadiv.PS <- plot_ordination(PS, ord.bray, color = "AnimalSpecies")+
  theme_bw()+
  geom_point(size=5, alpha= 0.75)

png("Figures/Betadiv_PCoA.png", units = 'in', res = 300, width=8, height=4)
betadiv.PS 
dev.off()

betadiv.PS.label <- plot_ordination(PS, ord.bray, color = "AnimalSpecies",
                                    label= "Barcode_Well", shape="Barcode_Plate")+
  theme_bw()

png("Figures/Betadiv_PCoA_label.png", units = 'in', res = 300, width=8, height=4)
betadiv.PS.label
dev.off()

###Bray-Curtis dissimilarity

PS1 <- subset_samples(PS, Barcode_Plate%in%"A1")
PS2 <- subset_samples(PS, Barcode_Plate%in%"A2")

dis.bray.1<- phyloseq::distance(PS1, method="bray")
dis.bray.2<- phyloseq::distance(PS2, method="bray")

sdat.1 <- sample_data(PS1)
sdat.2 <- sample_data(PS2)

## plate distance between all samples
as_tibble(sdat.1) %>%
    mutate(Ypos = as.numeric(as.factor(substr(Barcode_Well, 1, 1)))) %>%
    mutate(Xpos = as.numeric(substr(Barcode_Well, 2, 2))) %>%
    select(BeGenDiv_Name, Barcode_Plate, Barcode_Well, Ypos, Xpos) -> sdatPos.1

as_tibble(sdat.2) %>%
    mutate(Ypos = as.numeric(as.factor(substr(Barcode_Well, 1, 1)))) %>%
    mutate(Xpos = as.numeric(substr(Barcode_Well, 2, 2))) %>%
    select(BeGenDiv_Name, Barcode_Plate, Barcode_Well, Ypos, Xpos) -> sdatPos.2

XposDist.1 <- dist(sdatPos.1$Xpos)
XposDist.2 <- dist(sdatPos.2$Xpos)

YposDist.1 <- dist(sdatPos.1$Ypos)
YposDist.2 <- dist(sdatPos.2$Ypos)

CombiDist.1 <- XposDist.1+YposDist.1
CombiDist.2 <- XposDist.2+YposDist.2

CombiDist <- data.frame(plateX=c(as.vector(XposDist.1), as.vector(XposDist.2)),
                        plateY=c(as.vector(YposDist.1), as.vector(XposDist.2)),
                        plateCombi=c(as.vector(CombiDist.1), as.vector(CombiDist.2)),
                        bray=c(as.vector(dis.bray.1), as.vector(dis.bray.2)),
                        Plate=c(rep("A", times=length(XposDist.1)),
                                rep("B", times=length(XposDist.2))))

png("Figures/RowPCR_vs_bray.png", units = 'in', res = 300, width=6, height=4)
ggplot(CombiDist, aes(plateY, bray, color=Plate)) +
    geom_jitter(alpha=0.3, width=0.3, height=0) + 
    stat_smooth(se=FALSE, method="lm") +
    scale_x_continuous("Physical row-distance on PCR plate") +
    scale_y_continuous("Curtis-Bray dissimilarity between samples") +
    theme_bw()
dev.off()

png("Figures/ColumnPCR_vs_bray.png", units = 'in', res = 300, width=6, height=4)
ggplot(CombiDist, aes(plateX, bray, color=Plate)) +
    geom_jitter(alpha=0.3, width=0.3, height=0) + 
    stat_smooth(se=FALSE, method="lm") +
    scale_x_continuous("Physical column-distance on PCR plate") +
    scale_y_continuous("Curtis-Bray dissimilarity between samples") +
    theme_bw()
dev.off()

png("Figures/CombiPCR_vs_bray.png", units = 'in', res = 300, width=6, height=4)
ggplot(CombiDist, aes(plateCombi, bray, color=Plate)) +
    geom_jitter(alpha=0.3, width=0.3, height=0) + 
    stat_smooth(se=FALSE, method="lm") +
    scale_x_continuous("Combined physical distance on plate") +
    scale_y_continuous("Curtis-Bray dissimilarity between samples") +
    theme_bw()
dev.off()

##Filtering 
##1) Sample filtering: Filtering samples with low counts  
PS1 <- prune_samples(sample_sums(PS)>=2000, PS)
summarize_phyloseq(PS1)
##Filter low-occurrence, poorly-represented ASVs from this data
##Remove ASVs that do not show appear more than 2 times in more than half the samples
#asvkeep<-genefilter_sample(PS1, filterfun_sample(function(x) x > 2), A=0.5*nsamples(PS1))
#PS1<- prune_taxa(asvkeep, PS1)
##Transform to even sampling depth
## Rarefy without replacement
PS1<- rarefy_even_depth(PS1, rngseed=1, sample.size=0.9*min(sample_sums(PS1)), replace=F)
readcount(PS1)

## Merge ASVs that have the same taxonomy at a certain taxonomic rank (in this case Phylum and Family)
PS.Fam<-  tax_glom(PS1, "Family", NArm = F)
summarize_phyloseq(PS.Fam)

PS.Phy<-  tax_glom(PS1, "Phylum", NArm = F)
summarize_phyloseq(PS.Phy)

plot_bar(PS.Phy, fill="Phylum") + facet_wrap(~Compartment, scales= "free_x", nrow=1)

##Alpha diversity (rarefied)
plot_richness(PS1, x= "Compartment", color = "Compartment" , measures = c("Observed","Chao1", "Shannon")) +
  #geom_boxplot()+
  geom_jitter(alpha= 0.005)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))

alphadiv<- estimate_richness(PS1)
pairwise.wilcox.test(alphadiv$Observed, sample_data(PS1)$Compartment)

alphadiv%>%
  rownames_to_column()->tmp1

as.tibble(sample)%>%
  mutate(rowname= paste0("Sample", 1:nrow(sample)))->tmp2

tmp1<-inner_join(tmp1, tmp2, by="rowname")
rownames(tmp1)<- tmp1$rowname
tmp1$rowname<- NULL
alphadiv<- tmp1
rm(tmp1,tmp2)

require(ggpubr)
alphadiv%>%
  ggplot(aes(x= Compartment, y= Observed))+
  geom_boxplot(aes(color= Compartment))+
  geom_point(aes(color=Compartment))+
  xlab("Sample type")+
  scale_color_brewer(palette = "Set3")+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_compare_means(label= "p.signif", method = "t.test", ref.group = "Ascaris", paired = F, na.rm = TRUE)+
  stat_compare_means(method =  "anova", label.y = 10.5, label.x = 2)
