###Ascaris-Microbiome project (Ankur)
library("dada2") ##Following pipeline for v1.16

Refil<- FALSE ##Turn to TRUE just if filter parameters are changed 
Taxanot<- FALSE ##Turn to TRUE just in case new taxonomic assignment is required
Phylobj<- TRUE ##Start from ASV matrix, Tax matrix 

if(Refil){
###Load files 
path <- "/SAN/Victors_playground/Ascaris_Microbiome/2018_22_Nem1/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 

samples <- gsub("S\\d+-", "\\1", basename(fastqF))
samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(samples))

samplesMDC<- gsub("-\\d+_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))

###Check quality of the reads 
plotQualityProfile(fastqF[10:11])
plotQualityProfile(fastqR[10:11])

#filt_path <- "/SAN/Victors_playground/Ascaris_Microbiome/filtered_old" ##old filtering

filt_path <- "/SAN/Victors_playground/Ascaris_Microbiome/filtered2"  
#Pipeline filtration 
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

##Old filtering conditions
#out <- for(i in seq_along(fastqF)) {
#        fastqPairedFilter(c(fastqF[i], fastqR[i]), c(filtFs[i], filtRs[i]),
#                    truncLen=c(240,240),  
#                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, 
#                    compress=TRUE, verbose=TRUE)}

out <- filterAndTrim(fastqF, filtFs, fastqR, filtRs, truncLen=c(240, 240), ##This gives just 10bp overlap, Prev. conditions were 240,250
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft = c(17, 21), ##Remove primers
                     compress=TRUE, multithread=TRUE)
head(out)

###Learning errors
errF <- learnErrors(filtFs, multithread=TRUE)
#105311973 total bases in 472251 reads from 15 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
#105311973 total bases in 472251 reads from 15 samples will be used for learning the error rates.

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

##Sample inference 
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

##Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

##Construction of sequence table 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab))) ##Everything looks good but we have some short fragments 
plot(table(nchar(getSequences(seqtab))))

##Since our expected amplicon size is 440bp with primers, Let's make an in silico cut
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 400:430] ##New filtering gives amplicons ~426bp expected as amplicon
plot(table(nchar(getSequences(seqtab2))))

##Remove chimeras 
##Use just the data with the expected size
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
##With new filtering conditions 53.5% of the reads passed 
saveRDS(seqtab.nochim, "/SAN/Victors_playground/Ascaris_Microbiome/seqtab_final.rds") ##Final sequence table without chimeras
}

if(Taxanot){
  seqtab.nochim<- readRDS("/SAN/Victors_playground/Ascaris_Microbiome/seqtab_final.rds")
  
##Create an ASV matrix with ASV as rows and samples as columns (for Alessio)
asvmat <- t(seqtab.nochim) #Removing sequence rownames for display only
asv.names <- as.data.frame(rownames(asvmat))
rownames(asv.names) <- paste0("ASV", 1:nrow(asv.names))
rownames(asvmat)<-NULL
rownames(asvmat) <- paste0("ASV", 1:nrow(asvmat))
colnames(asvmat) <- paste0("Sample", 1:ncol(asvmat))
head(asvmat)
write.csv(asvmat, "/SAN/Victors_playground/Ascaris_Microbiome/output/ASV_matrix.csv")

##Get count of ASVs detected by sample
asv.sample<- as.data.frame(asvmat)
test<- data.frame()
for (i in 1:ncol(asv.sample)) {
  asv<- data.frame()
  asv[1,1]<- sum(asv.sample[,i]!=0)
  rownames(asv)<- paste0("Sample", i)
  test <- rbind(test, asv) ### Join all the "individual" data frames into the final data frame 
}
asv.sample<- as.matrix(test)
colnames(asv.sample)<- "ASVs_dada2"
rm(test,asv, i)

###Track reads through the pipeline 
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(fastqF, getN), sapply(fastqR, getN), sapply(filtFs, getN), sapply(filtRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track<-track[,c(1,3,5:8)]
colnames(track) <- c("input_dada2", "filtered_dada2", "denoisedF_dada2", "denoisedR_dada2", "merged_dada2", "nonchim_dada2")
rownames(track) <- samples
rownames(track) <- paste0("Sample", 1:nrow(track))
track<-cbind(track, asv.sample)
head(track)
saveRDS(track, "/SAN/Victors_playground/Ascaris_Microbiome/Track_DADA2.rds") ##Final track of dada2 pipeline

##Taxonomic annotation using naive Bayesian classifier from dada2 with SILVA db version 138
taxa <- assignTaxonomy(seqtab.nochim, "/SAN/db/RDP_Silva/Silva_138.1/dada2format/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/SAN/db/RDP_Silva/Silva_138.1/dada2format/silva_species_assignment_v138.fa.gz")
saveRDS(taxa, "/SAN/Victors_playground/Ascaris_Microbiome/tax_final.rds") ##Final taxonomy

##Create an Taxa matrix with ASV as rows and taxonomic level as columns (for Alessio)
taxamat <- taxa # Removing sequence rownames for display only
rownames(taxamat) <- NULL
rownames(taxamat) <- paste0("ASV", 1:nrow(taxamat))
head(taxamat)
write.csv(taxamat, "/SAN/Victors_playground/Ascaris_Microbiome/output/Taxa_matrix.csv")

##Transform seqtab.nochim to fasta file 
library(DECIPHER); packageVersion("DECIPHER")
##Create a DNAString set from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))
names(dna)<- paste0("ASV", 1:length(dna)) ##Give short names to each sequence
writeXStringSet(dna, "/SAN/Victors_playground/Ascaris_Microbiome/output/ASV.fasta") ##Export fasta seq

##Create a Phylogenetic tree
Align16S<- AlignSeqs(dna, anchor= NA, verbose= FALSE) ##Alignment

phangAlign16S <- phyDat(as(Align16S, "matrix"), type="DNA")
dm16S <- dist.ml(phangAlign16S) ## Distance matrix
treeNJ16S <- NJ(dm16S) # Note, tip order != sequence order
plot(treeNJ16S, type= "unrooted", use.edge.length= FALSE, no.margin= TRUE, show.tip.label= TRUE) ##Neighbour-Joining tree
fit1 <- pml(treeNJ16S, data=phangAlign16S)
fitGTR16S <- update(fit1, k=4, inv=0.2)
fitGTR16S <- optim.pml(fitGTR16S, model="GTR", optInv=TRUE, optGamma=TRUE,
                       rearrangement = "stochastic", control = pml.control(trace = 0))
plot(fitGTR16S, type= "unrooted", use.edge.length= FALSE, no.margin= TRUE, show.tip.label= FALSE)

#tree<- fitGTR16S$tree
#save tree (to add)

## Using IdTaxa taxonomic classification method (To be modified!)
##Load SILVA db version 138
#load("/SAN/db/RDP_Silva/Silva_138.1/dada2format/silva_nr99_v138_train_set.fa.gz") #
#ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
#ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
#taxid <- t(sapply(ids, function(x) {
#  m <- match(ranks, x$rank)
#  taxa <- x$taxon[m]
#  taxa[startsWith(taxa, "unclassified_")] <- NA
#  taxa
#}))
#colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
}

if(Phylobj){
  ##Load matrices
  asvmat<- read.csv("/SAN/Victors_playground/Ascaris_Microbiome/output/ASV_matrix.csv")
  rownames(asvmat)<-asvmat$X
  asvmat$X<- NULL
  taxamat<- read.csv("/SAN/Victors_playground/Ascaris_Microbiome/output/Taxa_matrix.csv")
  rownames(taxamat)<-taxamat$X
  taxamat$X<- NULL
  dna<- readDNAStringSet("/SAN/Victors_playground/Ascaris_Microbiome/output/ASV.fasta")
  ##Load sample data
  sample <- read.csv("/SAN/Victors_playground/Ascaris_Microbiome/Pig_Ascaris_16S_Samples_barcode.csv", dec=",", stringsAsFactors=FALSE)
 ##Add sample names used by Ulrike
  sample[ , "MDC_Names"] <- samplesMDC
  
##To phyloseq
library(phyloseq)
library(ggplot2)
library(dplyr)

#keep<- rownames(seqtab.nochim)
### Sample data includes those that didn't worked, so let's eliminate them 
#samdata <- samdata[samdata$BeGenDiv_Name %in% keep, ]
#rownames(samdata) <- samdata$sample_names

##To make Phyloseq object
##1) Use the ASV matrix and transform it to "OTU table" format
asv<- otu_table(asvmat, taxa_are_rows = T)
sample_names(asv)
##2) Use sample dataframe and transform it to "sample data" format
sample<- sample_data(sample)
sample_names(sample) <- sample_names(asv)
##3) Use taxa matrix and transform it to "tax table" format
tax<-tax_table(as.matrix(taxamat))
sample_names(tax)

PS <- merge_phyloseq(asv, tax)

###Add phylogenetic tree
require(ape)
phylotree<- rtree(ntaxa(PS), rooted=TRUE, tip.label=taxa_names(PS))

PS <- merge_phyloseq(asv, sample, tax, phylotree)

table(sample$System, sample$Compartment) ## ---> README sample overview (previous filtering)

saveRDS(PS, file="/SAN/Victors_playground/Ascaris_Microbiome/output/PhyloSeqComp.Rds")
saveRDS(sample, file="/SAN/Victors_playground/Ascaris_Microbiome/output/sample.Rds")
}

