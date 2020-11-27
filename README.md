# *Ascaris*-Pig-Microbiome project
Analysis of microbiome data from *Ascaris* infections in pigs with the dada2 pipeline

## Sample overview

### How many samples were included in the sequence library:
In total **177** samples

|         |Ascaris| Cecum| Colon| Duodenum| Faeces| Ileum| Jejunum| Negative|Total|
|---|---|---|---|---|---|---|---|---|---|
|DNA_kit  |     0 |    0 |   0  |      0  |    0  |   0  |      0 |       2 | 2   |  
|H2O      |     0 |    0 |   0  |      0  |    0  |   0  |      0 |       2 | 2   |
|Pig1     |     5 |    3 |   2  |      3  |    5  |   1  |      1 |       0 |20   |
|Pig2     |     7 |    4 |   2  |      2  |    5  |   0  |      3 |       0 |23   |
|Pig3     |    10 |    2 |   2  |      2  |    5  |   1  |      3 |       0 |25   |
|Pig4     |     0 |    2 |   2  |      2  |    5  |   0  |      4 |       0 |15   |
|Pig5     |     1 |    2 |   2  |      2  |    5  |   0  |      3 |       0 |15   |
|Pig6     |     0 |    2 |   2  |      2  |    5  |   1  |      3 |       0 |15   |
|Pig7     |     0 |    2 |   1  |      2  |    5  |   1  |      3 |       0 |14   |
|Pig8     |     0 |    2 |   2  |      2  |    5  |   0  |      3 |       0 |14   |
|Pig9     |     0 |    2 |   2  |      2  |    5  |   1  |      3 |       0 |15   |
|SH       |    15 |    0 |   0  |      0  |    0  |   0  |      0 |       0 |15   |
|Tris     |     0 |    0 |   0  |      0  |    0  |   0  |      0 |       2 | 2   |

### How many samples passed the filtering and sequencing processing:
In total **167** samples

|         |Ascaris| Cecum| Colon| Duodenum| Faeces| Ileum| Jejunum| Negative|Total|
|---|---|---|---|---|---|---|---|---|---|
|DNA_kit  |     0 |    0 |   0  |      0  |    0  |   0  |      0 |       2 | 2   |
|H2O      |     0 |    0 |   0  |      0  |    0  |   0  |      0 |       1 | 1   |
|Pig1     |     5 |    3 |   2  |      3  |    5  |   1  |      1 |       0 |20   |
|Pig2     |     7 |    4 |   2  |      2  |    4  |   0  |      3 |       0 |22   |
|Pig3     |    10 |    2 |   2  |      2  |    5  |   1  |      2 |       0 |24   |
|Pig4     |     0 |    2 |   2  |      2  |    5  |   0  |      4 |       0 |15   |
|Pig5     |     1 |    2 |   2  |      2  |    4  |   0  |      3 |       0 |14   |
|Pig6     |     0 |    2 |   2  |      2  |    4  |   1  |      3 |       0 |14   |
|Pig7     |     0 |    2 |   1  |      1  |    4  |   1  |      3 |       0 |12   |
|Pig8     |     0 |    2 |   2  |      2  |    4  |   0  |      3 |       0 |13   |
|Pig9     |     0 |    2 |   2  |      2  |    4  |   1  |      3 |       0 |14   |
|SH       |    15 |    0 |   0  |      0  |    0  |   0  |      0 |       0 |15   |
|Tris     |     0 |    0 |   0  |      0  |    0  |   0  |      0 |       1 | 1   |

## Goal 1: How much "eukaryote contamination" is present in the sequencing 

### How many ASVs and reads

|    |Archaea   | Bacteria  |Eukaryota   |  NA |
|---|---|---|---|---|
|number ASVs   | 0  | 16268  |0   |1  |
|number Sequences   | 0  | 1542188  |0   | 17  |
                      
-> NO Eukaryote contamination

## But another problem:  PCR chimeras 

The overall number of reads counted towards ASVs is very low: Starting
from 5935573 filtered reads means less than 40% (2374869) made it
through the pipeline: dada2 detected loads of *PCR chimeras* (more than
*50% of the reads*)!!??

## Goal 2: Figure out whether (even how) cross-contamination between PCR wells is present 

### Do negative samples with reads/ASVs reported taxonomically diverse composition?

YES! The tree problematic samples with a high number of reads have
also a diversity of ASVs.

![See this figure](https://github.com/VictorHJD/Ascaris_Microbiome/blob/master/Figures/Alphadiv_PCoA.png)

### Are contaminated negative samples similar in their composition to samples in their physical proximity on PCR plates

No! Negative samples are not particularly close in their composition to the samples surrounding them on PCR plates. 
 
![See this figure](https://github.com/VictorHJD/Ascaris_Microbiome/blob/master/Figures/Betadiv_PCoA_label.png)

NA are the negative controls.

### Does proximity of samples on PCR plates predict higher compositional similarity

Not sure. 

![This seems to suggest,](https://github.com/VictorHJD/Ascaris_Microbiome/blob/master/Figures/CombiPCR_vs_bray.png) but the samples close to each other have are also from the same experimental group. This could be an effect of either. 

(The combined distance is simply the sum of the row and column distance on the PCR plate, I also made plots for each sperately [see Figures]. 







