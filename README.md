# Ascaris_Microbiome
Analysis of microbiome data from Ascaris infections in pigs

## Goal 1: How much "eukaryote conamination" is present in the sequencing 

### How many ASVs and reads

|    ||Archaea   | Bacteria  |Eukaryota   |  NA |
|---|---|---|---|---|
|number ASVs   | 25  | 35588  |101   | 4  |
|number Sequences   | 329  | 1542188  |2602   | 105  |
                      
-> Eukaryote "contamination" seems neglegible I'd say

## But another problem:  PCR chimeras 

The overall number of reads counted towards ASVs is very low: Starting
from 5233526 filtered reads means less than 30% (1545224) made it
through the pipeline: dada2 detected loads of PCR chimreras (more than
50% of the reads)!!??

## Goal 2: Figure out whether (even how) cross-contamination between PCR wells is present 

### Do negative samples with reads/ASVs reported taxonomically diverse composition?

### YES!

![alt text](https://raw.githubusercontent.com/)

### Does proximity of samples on PCR plates predict higher compositonal similarity

![alt text](https://raw.githubusercontent.com/)




