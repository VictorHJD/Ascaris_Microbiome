###Ascaris-Microbiome project (Ankur)
library("MultiAmplicon")
library("dada2")
###Load files 
path <- "/SAN/Victors_playground/Ascaris_Microbiome/2018_22_Nem1/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 

samples <- gsub("S\\d+-", "\\1", basename(fastqF))
samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(samples))

###Check quality of the reads 
plotQualityProfile(fastqF[[177]])
plotQualityProfile(fastqR[[177]])

###They look okish!!!! Trimming at 250 should be ok!! 

filt_path <- "/SAN/Victors_playground/Ascaris_Microbiome/filtered"

#Pipeline filtration 
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

#out<- filterAndTrim(fastqF, filtFs, fastqR, filtRs,
#                truncLen=c(250,250), minLen=c(250,250), 
#                maxN=0, maxEE=c(2,2), truncQ=2, rm.phix = TRUE,
#                compress=TRUE, verbose=TRUE)

out <- for(i in seq_along(fastqF)) {
        fastqPairedFilter(c(fastqF[i], fastqR[i]), c(filtFs[i], filtRs[i]),
                    truncLen=c(240,250), 
                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE)}

#head(out)

##In order to not eliminate samples that did not pass the filter 
not.lost <- file.exists(filtFs) 
filtFs <- filtFs[not.lost]

not.lostR <- file.exists(filtRs) 
filtRs <- filtRs[not.lostR]

###Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- samples[not.lost]
names(derepRs) <- samples[not.lostR]

###Learning errors
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errF <- dadaFs.lrn[[1]]$err_out
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errR <- dadaRs.lrn[[1]]$err_out
#errF <- learnErrors(filtFs, multithread=TRUE)
#errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

##Sample inference 
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

##Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[2]])

##Construction of sequence table 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab))) ##Everything looks good but we have some short fragments 

##Let's cut everything above 400 bp 
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 400:473]

##Remove chimeras 
#seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#dim(seqtab.nochim)

##Use just the data with the expected size
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) ###Uuuu just 44% of the read pass... let's continue 

###Track reads through the pipeline 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "merged", "nonchim")
rownames(track) <- samples
head(track)

##Taxonomic annotation
taxa <- assignTaxonomy(seqtab.nochim, "/SAN/db/RDP/Silva_132/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/SAN/db/RDP/silva_species_assignment_v123.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

##Using Multiamplicon function
#taxa.MA <- blastTaxAnnot(seqtab.nochim,
#                    db = "/SAN/db/blastdb/nt/nt",
#                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
#                    infasta = "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/HMHZ2_2_in.fasta",
#                    outblast = "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/blast2_2_out.fasta",
#                    taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
#                    num_threads = 20)



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


