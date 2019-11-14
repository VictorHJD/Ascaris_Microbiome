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
rank_names(PS)
hist(rowSums(otu_table(PS)))
table(tax_table(PS)[, "Kingdom"], exclude = NULL) ## Mainly Bacteria, a bit of Eukaryota and Archaea
table(tax_table(PS)[, "Phylum"], exclude = NULL)
barplot.PS <- plot_bar(PS, fill = "Phylum") ##Without normalization

##Rarefaction curves
rarcurv <- vegan::rarecurve(otu_table(PS),
                            label = T, xlab = "Read number", ylab = "ASVs")

#pdf("~/GitProjects/Ascaris_Microbiome/Figures/Rarefaction_curves.pdf", 
#    width=8, height=10, onefile=T)
#rarcurv 
#dev.off()

##Eliminate empty samples 
PS <- prune_samples(sample_sums(PS)>0, PS)

##Separate by plate
PS.plate1<- subset_samples(PS, Barcode_Plate=="A1")
PS.plate2<- subset_samples(PS, Barcode_Plate=="A2")

#Visualize alpha-diversity (raw data)
alphaDiv.PS <- plot_richness(PS, x= "Compartment", color = "Compartment" , measures = c("Observed","Chao1", "Shannon")) +
  #geom_boxplot()+
  geom_jitter(alpha= 0.005)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))+
  labs(tag= "A)")
  
##Beta diversity (raw data)
ord.bray <- ordinate(PS, method="PCoA", distance="bray")

betadiv.PS <- plot_ordination(PS, ord.bray, color = "AnimalSpecies")+
  theme_bw()+
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "A)")

pdf("~/GitProjects/Ascaris_Microbiome/Figures/Betadiv_PCoA.pdf", 
    width=10, height=8, onefile=T)
betadiv.PS 
dev.off()

##Other ordination plots
plot_ordination(PS, ord.bray, color = "Barcode_Plate")+
  theme_bw()+
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "B)")

plot_ordination(PS, ord.bray, color = "DPI")+
  theme_bw()+
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "B)")

plot_ordination(PS, ord.bray, color = "Compartment")+
  theme_bw()+
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "B)")

plot_ordination(PS, ord.bray, color = "System")+
  theme_bw()+
  facet_wrap(~Compartment) +
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "B)")

plot_ordination(PS, ord.bray, color = "Compartment")+
  theme_bw()+
  facet_wrap(~Compartment) +
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "C)")

plot_ordination(PS, ord.bray, color = "Compartment")+
  theme_bw()+
  facet_wrap(~Compartment) +
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "C)")

###Bray-Curtis dissimilarity
dis.bray<- phyloseq::distance(PS, method="bray")
dis.bray.mat <- as.matrix(dis.bray) ##distance matrix all samples

##By plate
dis.bray1<- phyloseq::distance(PS.plate1, method="bray")
dis.bray.mat1 <- as.matrix(dis.bray1) ##distance matrix samples in PCR plate 1

dis.bray2<- phyloseq::distance(PS.plate2, method="bray")
dis.bray.mat2 <- as.matrix(dis.bray2) ##distance matrix samples in PCR plate 2

###Plate distance estimation (based on characters)
samnam <- sample_data(PS)$BeGenDiv_Name
position<- paste(sample_data(PS)$Barcode_Plate,sample_data(PS)$Barcode_Well)
plat.dist.mat <- stringdist::stringdistmatrix(position, position) ##Not perfect! Just counts character differences 
dimnames(plat.dist.mat)<- list(samnam, samnam)

matplot(plat.dist.mat, dis.bray.mat)

###Customized plate distance 
platcol.dist<-as.matrix(dist(as.numeric(as.factor(letters[1:12]))))

platrow.dist<-as.matrix(dist(as.numeric(as.factor(letters[1:8]))))

plate.dis<- matrix(nrow = 8, ncol= 12)
plate.dis[1,]<- c(0:11)
plate.dis[2,]<- c(1,2:12)
plate.dis[3,]<- c(2,3:13)
plate.dis[4,]<- c(3,4:14)
plate.dis[5,]<- c(4,5:15)
plate.dis[6,]<- c(5,6:16)
plate.dis[7,]<- c(6,7:17)
plate.dis[8,]<- c(7,8:18)


# Create table, number of features for each phyla
table(tax_table(PS)[, "Phylum"], exclude = NULL)
PS2 <- subset_taxa(PS, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(PS2),
                MARGIN = ifelse(taxa_are_rows(PS2), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(PS2),
                     tax_table(PS2))

prev <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

##Rarefaction
PSRare <- rarefy_even_depth(PS, rngseed=1, sample.size=0.9*min(sample_sums(PS)), replace=F)
