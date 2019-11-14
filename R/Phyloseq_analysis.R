library(ggplot2)
library(reshape)
library(phyloseq)
library(data.table)
library(parallel)
library(microbiome)

##Load data 
if(!exists("PS")){
  source("Dada2_Pipeline.R")
}

#PS<- readRDS(file="/SAN/Victors_playground/Ascaris_Microbiome/PhyloSeqCombi.Rds")

summarize_phyloseq(PS)
rarcurv <- vegan::rarecurve(otu_table(PS),
                            label = T)

PSRare <- rarefy_even_depth(PSHigh, rngseed=1, sample.size= min(sample_sums(PSHigh)), replace=F)

#Visualize alpha-diversity
plot_richness(PS, x= "Compartment", color = "Compartment" , measures = c("Observed","Chao1", "Shannon")) +
  #geom_boxplot()+
  geom_jitter(alpha= 0.005)+
  theme_bw()

##Beta diversity
ord.nmds.bray <- ordinate(PS, method="PCoA", distance="bray")
plot_ordination(PS, ord.nmds.bray, title="Bray NMDS")

##Barplot
top20 <- names(sort(taxa_sums(PS), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(PS, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="BeGenDiv_Name", fill="Phylum") #+ facet_wrap(~When, scales="free_x")


# Create table, number of features for each phyla
table(tax_table(PS)[, "Phylum"], exclude = NULL)
PS2 <- subset_taxa(PS, !is.na(Family) & !Family %in% c("", "uncharacterized"))

# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(PS2),
                MARGIN = ifelse(taxa_are_rows(PS2), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(PS2),
                     tax_table(PS2))

prev <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
#write.csv(prev, "~/Ascaris_Microbiome/Prevalence_Phylum.csv")
