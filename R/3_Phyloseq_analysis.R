###Filtering and rarefaction
library(phyloseq)
library(microbiome)
library(tidyverse)
library(data.table)

reRun <- FALSE

##Load data 
if(!exists("PS")){
  if(isTRUE(reRun)){
    source("1_Dada2_Pipeline.R") ## Run the script at base directory of repository!   
  } else {
    PS<- readRDS(file = "/SAN/Victors_playground/Ascaris_Microbiome/output/PhyloSeqComp.Rds") ##New annotation SILVA
    ##Eliminate "empty" samples 
    PS <- prune_samples(sample_sums(PS)>0, PS)
  }
}

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

phyla2Filter<- c("Aquificota", "Dependentiae")
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

PS3<- rarefy_even_depth(PS3, rngseed=2020, sample.size=min(sample_sums(PS3)), replace=F)
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
  dplyr::filter(Compartment%in%c("Duodenum"))%>%
  ggplot(aes(x= InfectionStatus, y= Observed))+
  geom_boxplot(color= "black", alpha= 0.5)+
  geom_jitter(shape=21, position=position_jitter(0.2), size=3, aes(fill= System), color= "black")+
  xlab("Sample type")+
  scale_color_brewer(palette = "Set3")+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_compare_means(label= "p.signif", method = "t.test", ref.group = "Non_infected", paired = F, na.rm = TRUE)+
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
biom.tmp<- make_biom(asvmat.rare, matrix_element_type = "int")
write_biom(biom.tmp,"/SAN/Victors_playground/Ascaris_Microbiome/output/biom_tmp.biom") ##Good biom for test

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

##Prune samples 
##Pig samples (not faeces)
PS3.Pig<- subset_samples(PS3, AnimalSpecies=="Pig"&Compartment!="Faeces")
sdt.pig <- data.table(as(sample_data(PS3.Pig), "data.frame"), keep.rownames = T)
##Just faeces samples
PS3.Fec<- subset_samples(PS3, Compartment=="Faeces")
sdt.fec <- data.table(as(sample_data(PS3.Fec), "data.frame"), keep.rownames = T)
##Pig samples (duodenum) and Ascaris
PS3.PA<- subset_samples(PS3, Compartment%in%c("Duodenum", "Ascaris"))
sdt.PA <- data.table(as(sample_data(PS3.PA), "data.frame"), keep.rownames = T)

##Pig samples (no faeces) and Ascaris
PS3.PA2<- subset_samples(PS3, Compartment!="Faeces")
sdt.PA2 <- data.table(as(sample_data(PS3.PA2), "data.frame"), keep.rownames = T)
sample_data(PS3.PA2)$Replicates<- paste(sdt.PA2$System, sdt.PA2$Compartment, sep = ".")
sdt.PA2$Replicate<- paste(sdt.PA2$System, sdt.PA2$Compartment, sep = ".")
sdt.PA2%>%
  select(InfectionStatus,AnimalSpecies,WormSex,Live,Compartment,System, Replicate)%>%
  distinct()->sdt.PA2
sdt.PA2<- sample_data(sdt.PA2)
sample_names(sdt.PA2) <- sdt.PA2$Replicate

PS3.PA2<-merge_samples(PS3.PA2, "Replicates")


bray_dist<- phyloseq::distance(PS3.Fec, method="bray", weighted=T)
ordination<- ordinate(PS3.Fec, method="PCoA", distance=bray_dist)
plot_ordination(PS3.Fec, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(shape=21, size=3, aes(fill= DPI), color= "black")+
  labs(tag= "D)")+
  theme_bw()+
  theme(text = element_text(size=16))

bray_dist<- phyloseq::distance(PS3.Pig, method="bray", weighted=T)
ordination<- ordinate(PS3.Pig, method="PCoA", distance=bray_dist)
plot_ordination(PS3.Pig, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(shape=21, size=3, aes(fill= Compartment), color= "black")+
  labs(tag= "D)")+
  theme_bw()+
  theme(text = element_text(size=16))

bray_dist<- phyloseq::distance(PS3.PA, method="bray", weighted=T)
ordination<- ordinate(PS3.PA, method="PCoA", distance=bray_dist)
plot_ordination(PS3.PA, ordination, color= "System", shape="Compartment")+ 
  theme(aspect.ratio=1)+
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "D)")+
  theme_bw()+
  scale_colour_brewer(type="qual", palette="Paired")+
  theme(text = element_text(size=16))

vegan::adonis(bray_dist ~ System+Compartment+Barcode_Plate+Barcode_Well, data = sdt.PA)

bray_dist<- phyloseq::distance(PS3.PA2, method="bray", weighted=T)
ordination<- ordinate(PS3.PA2, method="PCoA", distance=bray_dist)
plot_ordination(PS3.PA2, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "D)")+
  theme_bw()+
  scale_colour_brewer(type="qual", palette="Paired")+
  theme(text = element_text(size=16))