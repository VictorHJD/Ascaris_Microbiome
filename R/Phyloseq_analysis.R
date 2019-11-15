## library(ggplot2)
## library(reshape)
library(phyloseq)
library(microbiome)

## library(data.table)
## library(parallel)

reRun <- FALSE

##Load data 
if(!exists("PS")){
    if(isTRUE(reRun)){
        source("Dada2_Pipeline.R") ## Run the script at base directory of repository!   
    } else {
        PS<- readRDS(file="/SAN/Victors_playground/Ascaris_Microbiome/PhyloSeqCombi.Rds")
    }
}

### ## inspect
summarize_phyloseq(PS)
rank_names(PS)

## hist(rowSums(otu_table(PS)), xlab="number of reads", ylab="number of samples")

## ## HOW many ASVs for off-target eukaryotes and archaea
table(tax_table(PS)[, "Kingdom"], exclude = NULL) ## ---> README results summary

## ## HOW many reads for off-target eukaryotes and archaea
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

pdf("Figures/Alphadiv_PCoA.pdf", width=10, height=8, onefile=T)
alphaDiv.PS
dev.off()

##Beta diversity (raw data)
ord.bray <- ordinate(PS, method="PCoA", distance="bray")

betadiv.PS <- plot_ordination(PS, ord.bray, color = "AnimalSpecies")+
  theme_bw()+
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "A)")

pdf("Figures/Betadiv_PCoA.pdf", 
    width=10, height=8, onefile=T)
betadiv.PS 
dev.off()

##Other possible ordination plots
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
#samnam <- sample_data(PS)$BeGenDiv_Name
#position<- paste(sample_data(PS)$Barcode_Plate,sample_data(PS)$Barcode_Well)
#plat.dist.mat <- stringdist::stringdistmatrix(position, position) ##Not perfect! Just counts character differences 
#dimnames(plat.dist.mat)<- list(samnam, samnam)

#matplot(plat.dist.mat, dis.bray.mat)

###Customized plate distance 
#platcol.dist<-as.matrix(dist(as.numeric(as.factor(letters[1:12]))))
#platrow.dist<-as.matrix(dist(as.numeric(as.factor(letters[1:8]))))

##Create a plate with the distances 
plate.dis<- matrix(nrow = 8, ncol= 12)
plate.dis[1,]<- c(0:11)
plate.dis[2,]<- c(1,2:12)
plate.dis[3,]<- c(2,3:13)
plate.dis[4,]<- c(3,4:14)
plate.dis[5,]<- c(4,5:15)
plate.dis[6,]<- c(5,6:16)
plate.dis[7,]<- c(6,7:17)
plate.dis[8,]<- c(7,8:18)

colnames(plate.dis)<- as.numeric(as.factor(letters[1:12]))
rownames(plate.dis)<- as.factor(letters[1:8])

plate.dis<- t(plate.dis)
plate.dis<- as.data.frame(c(plate.dis))

##Link the distances whithin the plate to the positions 
plate.dis[,2]<- samdata$Barcode_Well[1:96]
plate.dis[64,2]<- "H8"
plate.dis[65:72,2]<- c("A9", "B9", "C9", "D9", "E9", "F9", "G9", "H9")
plate.dis[73:80,2]<- c("A10", "B10", "C10", "D10", "E10", "F10", "G10", "H10")
plate.dis[81:88,2]<- c("A11", "B11", "C11", "D11", "E11", "F11", "G11", "H11")
plate.dis[89:96,2]<- c("A12", "B12", "C12", "D12", "E12", "F12", "G12", "H12")
colnames(plate.dis)<- c("Plate_dis", "Plate_pos")

###Assign distances per plate position (PCR plate 1)
position.P1<- sample_data(PS.plate1)$Barcode_Well

position.P1<- as.data.frame(position.P1)
colnames(position.P1)<- "Plate_pos"

plate.1.dis<- plyr::join(position.P1, plate.dis, by= "Plate_pos")
rownames(plate.1.dis)<- sample_data(PS.plate1)$BeGenDiv_Name

###Create the distance matrix for PCR plate 1
plate.1.dis.mat<- as.matrix(dist(plate.1.dis$Plate_dis))
dimnames(plate.1.dis.mat)<- list(rownames(plate.1.dis), rownames(plate.1.dis)) 

##plot Bray-Curtis disimilarity vs plate position
plate.dist.bray.plot <- matplot(plate.1.dis.mat, dis.bray.mat1, pch = c(16,1), xlab = "Plate distance", ylab = "Bray-Curtis dissimilarity") ###It work!!! 

#pdf("~/GitProjects/Ascaris_Microbiome/Figures/Plate_distances_BC_dissimilarity.pdf", 
#    width=10, height=8, onefile=T)
#plate.dist.bray.plot
#dev.off()


# Create table, number of features for each phyla
#table(tax_table(PS)[, "Phylum"], exclude = NULL)
#PS2 <- subset_taxa(PS, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Compute prevalence of each feature, store as data.frame
#prevdf <- apply(X = otu_table(PS2),
#                MARGIN = ifelse(taxa_are_rows(PS2), yes = 1, no = 2),
#                FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
#prevdf <- data.frame(Prevalence = prevdf,
#                     TotalAbundance = taxa_sums(PS2),
#                     tax_table(PS2))

#prev <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

##Rarefaction
#PSRare <- rarefy_even_depth(PS, rngseed=1, sample.size=0.9*min(sample_sums(PS)), replace=F)
