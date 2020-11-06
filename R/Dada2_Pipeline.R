###Ascaris-Microbiome project (Ankur)
##library("MultiAmplicon") ##Not necessary for taxonomic assignment  
library("dada2") ##Following pipeline for v1.16
###Load files 
path <- "/SAN/Victors_playground/Ascaris_Microbiome/2018_22_Nem1/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 

samples <- gsub("S\\d+-", "\\1", basename(fastqF))
samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(samples))

###Check quality of the reads 
plotQualityProfile(fastqF[10:11])
plotQualityProfile(fastqR[10:11])

###They look okish!!!!

filt_path <- "/SAN/Victors_playground/Ascaris_Microbiome/filtered"

#Pipeline filtration 
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

out <- for(i in seq_along(fastqF)) {
        fastqPairedFilter(c(fastqF[i], fastqR[i]), c(filtFs[i], filtRs[i]),
                    truncLen=c(240,240), ##This gives just 10bp overlap, Prev. conditions were 240,250 
                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE)}

###Learning errors
errF <- learnErrors(filtFs, multithread=TRUE)
##100930560 total bases in 420544 reads from 15 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
##100930560 total bases in 420544 reads from 15 samples will be used for learning the error rates.

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
##Since our expected amplicon size is 460bp, Let's make an in silico cut: everything below 439 bp will be cutted 
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 439:468]
plot(table(nchar(getSequences(seqtab2))))
##Remove chimeras 
#seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#dim(seqtab.nochim)

##Use just the data with the expected size
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) ###Uuuu just 41% of the read pass... let's continue 

###Track reads through the pipeline 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "merged", "nonchim")
rownames(track) <- samples
head(track)

##Taxonomic annotation using naive Bayesian classifier from dada2 with SILVA db version 138
taxa <- assignTaxonomy(seqtab.nochim, "/SAN/db/RDP_Silva/Silva_138.1/dada2format/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/SAN/db/RDP_Silva/Silva_138.1/dada2format/silva_species_assignment_v138.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


## Using IdTaxa taxonomic classification method (To be modified!)
library(DECIPHER); packageVersion("DECIPHER")

##Create a DNAString set from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))
##Load SILVA db version 138
load("/SAN/db/RDP_Silva/Silva_138.1/dada2format/silva_nr99_v138_train_set.fa.gz") #
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)


##To phyloseq
library(phyloseq)
library(ggplot2)
library(dplyr)

samdata <- read.csv("/SAN/Victors_playground/Ascaris_Microbiome/Pig_Ascaris_16S_Samples_barcode.csv", dec=",", stringsAsFactors=FALSE)

#keep<- rownames(seqtab.nochim)
### Sample data includes those that didn't worked, so let's eliminate them 
#samdata <- samdata[samdata$BeGenDiv_Name %in% keep, ]
#rownames(samdata) <- samdata$sample_names

##To make Phyloseq object

asv<- otu_table(seqtab.nochim, taxa_are_rows=FALSE)
sample_names(asv)
sample<- sample_data(samdata)
sample_names(sample) <- sample_names(asv)
tax<-tax_table(taxa)
sample_names(tax)

###Add phylogenetic tree
require(ape)
phylotree<- rtree(ntaxa(PS), rooted=TRUE, tip.label=taxa_names(PS))

PS <- merge_phyloseq(asv, sample, tax, phylotree)

saveRDS(PS, file="/SAN/Victors_playground/Ascaris_Microbiome/PhyloSeqCombi.Rds") 


