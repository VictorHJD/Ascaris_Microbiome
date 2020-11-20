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
by((otu_table(PS)), tax_table(PS)[, "Kingdom"], sum) ## --->  README results summary

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
    scale_y_continuous("Bray-Curtis dissimilarity between samples") +
    theme_bw()
dev.off()

png("Figures/ColumnPCR_vs_bray.png", units = 'in', res = 300, width=6, height=4)
ggplot(CombiDist, aes(plateX, bray, color=Plate)) +
    geom_jitter(alpha=0.3, width=0.3, height=0) + 
    stat_smooth(se=FALSE, method="lm") +
    scale_x_continuous("Physical column-distance on PCR plate") +
    scale_y_continuous("Bray-Curtis dissimilarity between samples") +
    theme_bw()
dev.off()

png("Figures/CombiPCR_vs_bray.png", units = 'in', res = 300, width=6, height=4)
ggplot(CombiDist, aes(plateCombi, bray, color=Plate)) +
    geom_jitter(alpha=0.3, width=0.3, height=0) + 
    stat_smooth(se=FALSE, method="lm") +
    scale_x_continuous("Combined physical distance on plate") +
    scale_y_continuous("Bray-Curtis dissimilarity between samples") +
    theme_bw()
dev.off()

###Model 1: Can the Bray-Curtis dissimilarity between samples be predicted by the position in the plate 
summary(lm(bray~plateCombi, data = subset(CombiDist, Plate== "A"))) 
summary(lm(bray~plateCombi, data = subset(CombiDist, Plate== "B")))

##Filtering 
##1) Sample filtering: Filtering samples with low counts  
PS3 <- prune_samples(sample_sums(PS)>=2000, PS)
summarize_phyloseq(PS3)
##2) Taxa filtering: Remove "uncharachterized" ASVs
PS3<- subset_taxa(PS3, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

##3) Taxa filtering: Remove low prevalent taxa
##Create a prevalence dataframe 
Prevdf<- apply(X = otu_table(PS3),
                     MARGIN = 1,
                     FUN = function(x){sum(x > 0)})

##Add taxonomy and total read counts to this data.frame
Prevdf<- data.frame(Prevalence = Prevdf,
                          TotalAbundance = taxa_sums(PS3),
                          tax_table(PS3))

plyr::ddply(Prevdf, "Phylum", function(df1){
  data.frame(mean_prevalence=mean(df1$Prevalence),total_abundance=sum(df1$TotalAbundance,na.rm = T),stringsAsFactors = F)
})

phyla2Filter<- c("Deinococcota", "Dependentiae")
PS3<- subset_taxa(PS3, !Phylum %in% phyla2Filter)

Prevdf<- apply(X = otu_table(PS3),
               MARGIN = 1,
               FUN = function(x){sum(x > 0)})
Prevdf<- data.frame(Prevalence = Prevdf,
                    TotalAbundance = taxa_sums(PS3),
                    tax_table(PS3))

ggplot(Prevdf, aes(TotalAbundance, Prevalence / nsamples(PS3),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Log 10 Total Reads") + ylab("Prevalence [Prop. of Samples]") +
  theme_bw()+
  facet_wrap(~Phylum) + theme(legend.position="none")
  
##4) Transform to even sampling depth
## Rarefy without replacement
vegan::rarecurve(t(otu_table(PS3)), step=50, cex=0.5)
set.seed(2020)
PS3<- rarefy_even_depth(PS3, rngseed=1, sample.size=0.99*min(sample_sums(PS3)), replace=F)
readcount(PS3)

## Merge ASVs that have the same taxonomy at a certain taxonomic rank (in this case Phylum and Family)
PS.Fam<-  tax_glom(PS3, "Family", NArm = F)
summarize_phyloseq(PS.Fam)

PS.Phy<-  tax_glom(PS3, "Phylum", NArm = F)
summarize_phyloseq(PS.Phy)

plot_bar(PS.Phy, fill="Phylum") + facet_wrap(~Compartment, scales= "free_x", nrow=1)

##Alpha diversity (rarefied)
plot_richness(PS3, x= "Compartment", color = "Compartment" , measures = c("Observed","Chao1", "Shannon")) +
  #geom_boxplot()+
  geom_jitter(alpha= 0.005)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))

alphadiv<- estimate_richness(PS3)

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
  geom_boxplot(color= "black", alpha= 0.5)+
  geom_jitter(shape=21, position=position_jitter(0.2), size=3, aes(fill= Compartment), color= "black")+
  xlab("Sample type")+
  scale_color_brewer(palette = "Set3")+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_compare_means(label= "p.signif", method = "t.test", ref.group = "Ascaris", paired = F, na.rm = TRUE)+
  stat_compare_means(method =  "anova", label.y = 10.5, label.x = 2)

pairwise.wilcox.test(alphadiv$Observed, alphadiv$Compartment, p.adjust.method = "bonferroni")

##Beta diversity (rarefied)
# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist<- phyloseq::distance(PS3, method="unifrac", weighted=F)
ordination<- ordinate(PS3, method="PCoA", distance=wunifrac_dist)
plot_ordination(PS3, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(shape=21, size=3, aes(fill= Compartment), color= "black")+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))

vegan::adonis(wunifrac_dist ~ sample_data(PS3)$Compartment)

# PCoA plot using the weighted UniFrac as distance
unifrac_dist<- phyloseq::distance(PS3, method="unifrac", weighted=T)
ordination<- ordinate(PS3, method="PCoA", distance=unifrac_dist)
plot_ordination(PS3, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(shape=21, size=3, aes(fill= Compartment), color= "black")+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))

# PCoA plot using Bray-Curtis as distance
bray_dist<- phyloseq::distance(PS3, method="bray", weighted=T)
ordination<- ordinate(PS3, method="PCoA", distance=bray_dist)
plot_ordination(PS3, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(shape=21, size=3, aes(fill= Compartment), color= "black")+
  labs(tag= "D)")+
  theme_bw()+
  theme(text = element_text(size=16))

vegan::adonis(bray_dist ~ sample_data(PS3)$Compartment)

##Create biom format object for PICRUSt2
require("biomformat")
asvmat.rare<- as.matrix(PS3@otu_table)
taxmat.rare<- as.matrix(PS3@tax_table)
sample.rare<- as.data.frame(PS3@sam_data)

biom.rare<- make_biom(asvmat.rare, sample_metadata = sample.rare, observation_metadata = taxmat.rare,
                      id = NULL, matrix_element_type = "int")

write_biom(biom.rare,"/SAN/Victors_playground/Ascaris_Microbiome/output/rare.biom") ##This biom file is not working

biom.tmp<- make_biom(asvmat.rare, matrix_element_type = "int")
write_biom(biom.tmp,"/SAN/Victors_playground/Ascaris_Microbiome/output/biom_tmp.biom") ##Temporal biom for test

##Select sequences from the ASV in PS3
keep <- data.frame(name = rownames(asvmat.rare))
names(dna)
dna.rare<- dna[keep$name]
writeXStringSet(dna.rare, "/SAN/Victors_playground/Ascaris_Microbiome/output/Rare_ASV.fasta") #-> For Picrust2

#write.table(asvmat, "/SAN/Victors_playground/Ascaris_Microbiome/output/Rare_ASV_matrix.txt", sep = "\t") #-> For Picrust2
#write.csv(taxamat, "/SAN/Victors_playground/Ascaris_Microbiome/output/Rare_Taxa_matrix.csv") #-> For Picrust2
#write.csv(as.matrix(bray_dist), "/SAN/Victors_playground/Ascaris_Microbiome/output/Rare_Bray_Curtis.csv")

#Transform to relative abundance matix summarized to genus level
PS4<- transform_sample_counts(PS3, function(x) x / sum(x) )

#asvmat.ra<- as.matrix(PS4@otu_table)
#taxmat.ra<- as.matrix(PS4@tax_table)
#write.table(asvmat.ra, "/SAN/Victors_playground/Ascaris_Microbiome/output/RelAb_ASV_matrix.txt", sep = "\t") 
#write.csv(taxmat.ra, "/SAN/Victors_playground/Ascaris_Microbiome/output/RelAb_Taxa_matrix.csv") 


PS.Gen<-  tax_glom(PS4, "Genus", NArm = T)
summarize_phyloseq(PS.Gen)

asvmat.gen<- as.matrix(PS.Gen@otu_table)
taxmat.gen<- as.data.frame(PS.Gen@tax_table)
taxmat.gen<- taxmat.gen[,"Genus", drop=FALSE]

genus.relab<- cbind(asvmat.gen, taxmat.gen)
genus.relab<- genus.relab[, c(167, 1:166)]
gsub(" ", "_", genus.relab$Genus)
write.table(genus.relab, "/SAN/Victors_playground/Ascaris_Microbiome/output/Genus_relative_abundance.txt", 
            sep = "\t", row.names = FALSE, quote=F)