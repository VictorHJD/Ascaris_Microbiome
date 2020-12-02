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
    sample<- data.table(as(sample_data(PS), "data.frame"), keep.rownames = T)
    ##Eliminate "empty" samples 
    PS <- prune_samples(sample_sums(PS)>0, PS)
  }
}
##Extract sample data

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

table(alphadiv$System, alphadiv$Compartment) ## ---> README sample overview (post filtering)

##Prune samples for different questions
##Pig samples (Just GI compartment)
PS3.Pig<- subset_samples(PS3, !(Compartment%in%c("Faeces", "Negative", "Ascaris")))
sdt.pig <- data.table(as(sample_data(PS3.Pig), "data.frame"), keep.rownames = T)
alphadiv.pig<- estimate_richness(PS3.Pig) ###Estimate alpha diversity values 
alphadiv.pig%>%
  rownames_to_column()->tmp1

row.names(sdt.pig)<- sdt.pig$rn
names(sdt.pig)[names(sdt.pig) == "rn"] <- "rowname"

tmp1<-inner_join(tmp1, sdt.pig, by="rowname")
rownames(tmp1)<- tmp1$rowname
tmp1$rowname<- NULL
sdt.pig<- tmp1
rm(tmp1,alphadiv.pig)

##Merge compartment data
##Pig samples (no faeces) 
sample_data(PS3.Pig)$Replicates<- paste(sdt.pig$System, sdt.pig$Compartment, sep = ".")
sdt.pig$Replicate<- paste(sdt.pig$System, sdt.pig$Compartment, sep = ".")
sdt.pig%>%
  select(InfectionStatus,AnimalSpecies,WormSex,Live,Compartment,System, Replicate)%>%
  distinct()->sdt.pig2
sdt.pig2<- sample_data(sdt.pig2)
sample_names(sdt.pig2) <- sdt.pig2$Replicate

PS3.pig2<-merge_samples(PS3.Pig, "Replicates")
sample_data(PS3.pig2)<- sdt.pig2

alphadiv.pig2<- estimate_richness(PS3.pig2) ###Estimate alpha diversity values
alphadiv.pig2%>%
  rownames_to_column()->tmp1

row.names(sdt.pig2)<- sdt.pig2$Replicate
names(sdt.pig2)[names(sdt.pig2) == "Replicate"] <- "rowname"

tmp1<-inner_join(tmp1, sdt.pig2, by="rowname")
rownames(tmp1)<- tmp1$rowname
tmp1$rowname<- NULL
sdt.pig2<- tmp1
rm(tmp1,alphadiv.pig2)

##Just faeces samples (difference accross time)
PS3.Fec<- subset_samples(PS3, Compartment=="Faeces")
sdt.fec <- data.table(as(sample_data(PS3.Fec), "data.frame"), keep.rownames = T)
alphadiv.fec<- estimate_richness(PS3.Fec) ###Estimate alpha diversity values
alphadiv.fec%>%
  rownames_to_column()->tmp1

row.names(sdt.fec)<- sdt.fec$rn
names(sdt.fec)[names(sdt.fec) == "rn"] <- "rowname"

tmp1<-inner_join(tmp1, sdt.fec, by="rowname")
rownames(tmp1)<- tmp1$rowname
tmp1$rowname<- NULL
sdt.fec<- tmp1
rm(tmp1,alphadiv.fec)

##Fecal samples have no technical replicates NOT merge

##Pig samples compartment and Ascaris (not SH)
PS3.PA<- subset_samples(PS3, !(Compartment%in%c("Negative", "Faeces")))
PS3.PA<- subset_samples(PS3.PA, !(System%in%c("SH")))
sdt.PA <- data.table(as(sample_data(PS3.PA), "data.frame"), keep.rownames = T)
alphadiv.PA<- estimate_richness(PS3.PA) ###Estimate alpha diversity values
alphadiv.PA%>%
  rownames_to_column()->tmp1

row.names(sdt.PA)<- sdt.PA$rn
names(sdt.PA)[names(sdt.PA) == "rn"] <- "rowname"

tmp1<-inner_join(tmp1, sdt.PA, by="rowname")
rownames(tmp1)<- tmp1$rowname
tmp1$rowname<- NULL
sdt.PA<- tmp1
rm(tmp1,alphadiv.PA)

##Merge compartment and Ascaris data
sample_data(PS3.PA)$Replicates<- paste(sdt.PA$System, sdt.PA$Compartment, sep = ".")
sdt.PA$Replicate<- paste(sdt.PA$System, sdt.PA$Compartment, sep = ".")
sdt.PA%>%
  select(InfectionStatus,AnimalSpecies,Live,Compartment,System, Replicate)%>%
  distinct()->sdt.PA2
sdt.PA2<- sample_data(sdt.PA2)
sample_names(sdt.PA2) <- sdt.PA2$Replicate

PS3.PA2<-merge_samples(PS3.PA, "Replicates")
sample_data(PS3.PA2)<- sdt.PA2

alphadiv.PA2<- estimate_richness(PS3.PA2) ###Estimate alpha diversity values
alphadiv.PA2%>%
  rownames_to_column()->tmp1

row.names(sdt.PA2)<- sdt.PA2$Replicate
names(sdt.PA2)[names(sdt.PA2) == "Replicate"] <- "rowname"

tmp1<-inner_join(tmp1, sdt.PA2, by="rowname")
rownames(tmp1)<- tmp1$rowname
tmp1$rowname<- NULL
sdt.PA2<- tmp1
rm(tmp1,alphadiv.PA2)

##Ascaris samples
PS3.Asc<- subset_samples(PS3, Compartment%in%c("Ascaris"))
sdt.Asc <- data.table(as(sample_data(PS3.Asc), "data.frame"), keep.rownames = T)
alphadiv.Asc<- estimate_richness(PS3.Asc) ###Estimate alpha diversity values
alphadiv.Asc%>%
  rownames_to_column()->tmp1

row.names(sdt.Asc)<- sdt.Asc$rn
names(sdt.Asc)[names(sdt.Asc) == "rn"] <- "rowname"

tmp1<-inner_join(tmp1, sdt.Asc, by="rowname")
rownames(tmp1)<- tmp1$rowname
tmp1$rowname<- NULL
sdt.Asc<- tmp1
rm(tmp1,alphadiv.Asc)

###Question 1:
###How does Ascaris impact the porcine microbiome and does this differ in different gut regions?
###General comparison between infected and non infected pigs (all compartments and not merged replicates)
require(ggpubr)
require(RColorBrewer)
require(rstatix)
###General comparison between infected and non infected pigs by compartments (merged replicates)
###Group comparisons
###Change to sdt.pig for not merged samples
### Infected vs Non Infected
sdt.pig2%>% 
  mutate(Compartment = fct_relevel(Compartment, 
                                             "Duodenum", "Jejunum", "Ileum", 
                                             "Cecum", "Colon"))%>%
  dplyr::group_by(Compartment)%>%
  wilcox_test(Chao1 ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Infected_NonInfected_Compartment.csv")

sdt.pig2%>% 
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  dplyr::group_by(Compartment)%>%
  wilcox_effsize(Chao1 ~ InfectionStatus)

##Plot 
sdt.pig2%>%
  dplyr::filter(!(Compartment%in%c("Faeces", "Negative", "Ascaris")))%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  ggplot(aes(x= Compartment, y= Chao1))+
  geom_boxplot(aes(color= InfectionStatus), alpha= 0.5)+
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= InfectionStatus), color= "black")+
  xlab("GI compartment")+
  ylab("Diversity (Chao1 Index)")+
  labs(tag= "A)", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -2, hide.ns = F,label = "{p.adj}{p.adj.signif}")-> A

###Infected pigs, non infected compartment
sdt.pig2%>% 
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  dplyr::group_by(InfectionStatus)%>%
  wilcox_test(Chao1 ~ Compartment)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Compartment.csv")

sdt.pig2%>% 
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  dplyr::group_by(InfectionStatus)%>%
  wilcox_effsize(Chao1 ~ Compartment)
  
##PLot 
sdt.pig2%>%
  dplyr::filter(!(Compartment%in%c("Faeces", "Negative", "Ascaris")))%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  ggplot(aes(x= Compartment, y= Chao1))+
  geom_boxplot(color= "black", alpha= 0.5)+
  geom_jitter(shape=21, position=position_jitter(0.2), size=3, aes(fill= System), color= "black")+
  xlab("GI compartment")+
  ylab("Diversity (Chao1 Index)")+
  scale_color_brewer(palette = "Set3")+
  labs(tag= "B)", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 45))+
  facet_wrap(~InfectionStatus)+ 
  stat_pvalue_manual(stats.test, bracket.nudge.y = -2, hide.ns = T,label = "{p.adj}{p.adj.signif}")-> B

##Beta diversity (rarefied)
# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist<- phyloseq::distance(PS3.pig2,
                                   method="unifrac", weighted=F)
ordination<- ordinate(PS3.pig2,
                      method="PCoA", distance=wunifrac_dist)
plot_ordination(PS3.pig2, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(shape=21, size=3, aes(fill= Compartment), color= "black")+
  labs(title = "Unweighted UniFrac",tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))

# PCoA plot using the weighted UniFrac as distance
unifrac_dist<- phyloseq::distance(PS3.pig2,
                                  method="unifrac", weighted=T)
ordination<- ordinate(PS3.pig2,
                      method="PCoA", distance=unifrac_dist)
plot_ordination(PS3.pig2, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(shape=21, size=3, aes(fill= Compartment), color= "black")+
  labs(title = "Weighted UniFrac", tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))

# PCoA plot using Bray-Curtis as distance
bray_dist<- phyloseq::distance(PS3.pig2, 
                               method="bray", weighted=T)
ordination<- ordinate(PS3.pig2,
                      method="PCoA", distance=bray_dist)
plot_ordination(PS3.pig2, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(shape=21, size=3, aes(fill= Compartment), color= "black")+
  labs(title = "Bray-Curtis dissimilarity", tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_text (x = 0.00, y = -0.35, label = paste ("Bray-Curtis~ Compartment+Pig+Infection, \n PERMANOVA, Compartment p= 0.001; R-squared= 0.2665"))-> C

##Adonis PERMANOVA
#Is beta diversity vary function of the compartment, pig, infection status or technical predictors
bd.pig<- vegan::adonis(bray_dist~ Compartment + System + InfectionStatus,
              permutations = 999, data = sdt.pig2) ##Technical factors can't be included for merged samples
bd.pig ##Manually added to the plot 

write.csv(bd.pig[[1]], "Tables/Q1_Permanova.csv")

###Yes, Compartment and Plate are the most significant predictors
## Calculate multivariate dispersion (aka distance to the centroid)
mvd.pig<- vegan::betadisper(bray_dist, sdt.pig2$Compartment, type = "centroid")
vegan::permutest(mvd.pig, permutations = 999)
anova(mvd.pig)
plot(mvd.pig) ##Same than plot C
boxplot(mvd.pig)
plot(TukeyHSD(mvd.pig))
###Add sample to the centroid from each sample to sdt.pig
tmp<- as.data.frame(mvd.pig$distances)
colnames(tmp)<- c("distances")
cbind(sdt.pig, tmp)-> sdt.pig2
###Linear model 
summary(lm(sdt.pig2, formula = distances~ Compartment, na.action = na.exclude))

###Group comparisons 
### Infected vs Non Infected
sdt.pig2%>% 
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  dplyr::group_by(Compartment)%>%
  wilcox_test(distances ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Infected_NonInfected_distances.csv")

sdt.pig2%>% 
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  dplyr::group_by(Compartment)%>%
  wilcox_effsize(distances ~ InfectionStatus)

##Plot 
sdt.pig2%>%
  dplyr::filter(!(Compartment%in%c("Faeces", "Negative", "Ascaris")))%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  ggplot(aes(x= Compartment, y= distances))+
  geom_boxplot(aes(color= InfectionStatus), alpha= 0.5)+
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= InfectionStatus), color= "black")+
  xlab("GI compartment")+
  ylab("Beta diversity (distance to coentroid)")+
  labs(tag= "A)", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_pvalue_manual(stats.test, bracket.nudge.y = 0.05, hide.ns = T,label = "{p.adj}{p.adj.signif}")-> D

###Infected pigs, non infected by compartment
sdt.pig2%>% 
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  dplyr::group_by(InfectionStatus)%>%
  wilcox_test(distances ~ Compartment)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q1_Compartment_distances.csv")

sdt.pig2%>% 
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  dplyr::group_by(InfectionStatus)%>%
  wilcox_effsize(distances ~ Compartment)

##PLot 
sdt.pig2%>%
  dplyr::filter(!(Compartment%in%c("Faeces", "Negative", "Ascaris")))%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  ggplot(aes(x= Compartment, y= distances))+
  geom_boxplot(color= "black", alpha= 0.5)+
  geom_jitter(shape=21, position=position_jitter(0.2), size=3, aes(fill= System), color= "black")+
  xlab("GI compartment")+
  ylab("Beta diversity (distance to centroid)")+
  labs(tag= "B)", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 45))+
  facet_wrap(~InfectionStatus)+ 
  stat_pvalue_manual(stats.test, bracket.nudge.y = -.095, hide.ns = TRUE,label = "{p.adj}{p.adj.signif}")->E

require(grid)
require(gridExtra)
png("Figures/Q1_Alphadiv_Compartment.png", units = 'in', res = 300, width=14, height=14)
grid.arrange(A, B)
dev.off()

png("Figures/Q1_PCoA_Betadiv_Compartment.png", units = 'in', res = 300, width=10, height=8)
grid.arrange(C)
dev.off()

png("Figures/Q1_Betadiv_Distances_Compartment.png", units = 'in', res = 300, width=14, height=14)
grid.arrange(D,E)
dev.off()

###Question 2:
###How does the Ascaris microbiome compare to that of the porcine intestinal microbiome?
###General comparison between infected pigs (Jejunum and Ascaris and not merged replicates)
### Infected vs Non Infected
sdt.PA2$InfectionStatus[is.na(sdt.PA2$InfectionStatus)] <- "Worm" ##Remove NAs
sdt.PA2%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon", "Ascaris"))%>%
  #dplyr::group_by(Compartment)%>%
  wilcox_test(Chao1 ~ InfectionStatus)%>%
  adjust_pvalue(method = "bonferroni")%>%
  add_significance()%>%
  add_xy_position(x = "InfectionStatus")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q2_Infected_NonInfected_Worm.csv")

sdt.PA2%>% 
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  #dplyr::group_by(Compartment)%>%
  wilcox_effsize(Chao1 ~ InfectionStatus)

sdt.PA2%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon", "Ascaris"))%>%
  ggplot(aes(x= Compartment, y= Chao1))+
  geom_boxplot(aes(color= InfectionStatus), alpha= 0.5)+
  scale_color_manual(values = c("#F8766D", "#619CFF", "#00BA38"),
                     aesthetics = c("colour", "fill"))+
  geom_point(shape=21, position=position_jitter(0.2), size=3, 
             aes(fill= InfectionStatus, color= InfectionStatus), color= "black")+
  xlab("GI compartment")+
  ylab("Diversity (Chao1)")+
  theme_bw()+
  labs(tag= "A)", caption = get_pwc_label(stats.test))+
  theme(text = element_text(size=16))-> f

###Infected pigs, non infected compartment
sdt.PA2%>% 
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  #dplyr::group_by(InfectionStatus)%>%
  wilcox_test(Chao1 ~ Compartment)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q2_Compartment.csv")

sdt.pig2%>% 
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  #dplyr::group_by(InfectionStatus)%>%
  wilcox_effsize(Chao1 ~ Compartment)

##Plot 
sdt.PA2%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  ggplot(aes(x= Compartment, y= Chao1))+
  geom_boxplot(color= "black", alpha= 0.5)+
  geom_point(shape=21, position=position_jitter(0.2), size=3, 
             aes(fill= System), color= "black")+
  xlab("GI compartment")+
  ylab("Diversity (Chao1)")+
  labs(tag= "B)", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 45))+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -2, hide.ns = T,label = "{p.adj}{p.adj.signif}")-> G

##Beta diversity (rarefied)
# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist<- phyloseq::distance(PS3.PA2,
                                   method="unifrac", weighted=F)
ordination<- ordinate(PS3.PA2,
                      method="PCoA", distance=wunifrac_dist)
plot_ordination(PS3.PA2, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(shape= Compartment, color= System))+
  scale_color_brewer(palette = "Paired")+
  labs(title = "Unweighted UniFrac",tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))

# PCoA plot using the weighted UniFrac as distance
unifrac_dist<- phyloseq::distance(PS3.PA2,
                                  method="unifrac", weighted=T)
ordination<- ordinate(PS3.PA2,
                      method="PCoA", distance=unifrac_dist)
plot_ordination(PS3.PA2, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(size=3, aes(shape= Compartment, color= System))+
  scale_color_brewer(palette = "Paired")+
  labs(title = "Weighted UniFrac",tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))

# PCoA plot using Bray-Curtis as distance
bray_dist<- phyloseq::distance(PS3.PA2, 
                               method="bray", weighted=T)
ordination<- ordinate(PS3.PA2,
                      method="PCoA", distance=bray_dist)
plot_ordination(PS3.PA2, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(shape=21, size=3, aes(fill= Compartment), color= "black")+
  labs(title = "Bray-Curtis dissimilarity",tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_text (x = 0, y = 0.35, 
             label = paste ("Bray-Curtis~ Compartment+Pig+Infection, \n PERMANOVA, Compartment p= 0.001; R^2= 0.3009"))-> H

##Adonis PERMANOVA
#Is beta diversity vary function of the compartment, pig, infection status or technical predictors
bd.PA<- vegan::adonis(bray_dist~ Compartment + System + InfectionStatus,
                       permutations = 999, data = sdt.PA2)
bd.PA ##Manually added to the plot 
write.csv(bd.PA[[1]], "Tables/Q2_Permanova.csv")

###Yes, Compartment, System and Plate are the most significant predictors
## Calculate multivariate dispersion (aka distance to the centroid)
mvd.PA<- vegan::betadisper(bray_dist, sdt.PA2$Compartment, type = "centroid")
vegan::permutest(mvd.PA, permutations = 999)
anova(mvd.PA)
plot(mvd.PA)
boxplot(mvd.PA)
plot(TukeyHSD(mvd.PA))
###Add sample to the centroid from each sample to sdt.PA
tmp<- as.data.frame(mvd.PA$distances)
colnames(tmp)<- c("distances")
cbind(sdt.PA2, tmp)-> sdt.PA2
###Linear model 
summary(lm(sdt.PA2, formula = distances~ Compartment, na.action = na.exclude))

###Infected pigs, non infected by compartment
sdt.PA2%>% 
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  #dplyr::group_by(InfectionStatus)%>%
  wilcox_test(distances ~ Compartment)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "Compartment")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q2_Compartment_distances.csv")

sdt.PA2%>% 
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  #dplyr::group_by(InfectionStatus)%>%
  wilcox_effsize(distances ~ Compartment)

##PLot 
sdt.PA2%>%
  mutate(Compartment = fct_relevel(Compartment, 
                                   "Duodenum", "Jejunum", "Ileum", 
                                   "Cecum", "Colon"))%>%
  ggplot(aes(x= Compartment, y= distances))+
  geom_boxplot(color= "black", alpha= 0.5)+
  geom_jitter(shape=21, position=position_jitter(0.2), size=3, aes(fill= System), color= "black")+
  xlab("GI compartment")+
  ylab("Beta diversity (distance to centroid)")+
  labs(tag= "B)", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 45))+
  #facet_wrap(~InfectionStatus)+ 
  stat_pvalue_manual(stats.test, bracket.nudge.y = -.095, hide.ns = TRUE,label = "{p.adj}{p.adj.signif}")->I

png("Figures/Q2_Alphadiv_Compartment.png", units = 'in', res = 300, width=14, height=14)
grid.arrange(f, G)
dev.off()

png("Figures/Q2_PCoA_Betadiv_Compartment.png", units = 'in', res = 300, width=10, height=8)
grid.arrange(H)
dev.off()

png("Figures/Q2_Betadiv_Distances_Compartment.png", units = 'in', res = 300, width=14, height=14)
grid.arrange(I)
dev.off()

###Question 3: Does the Ascaris microbiome differ between male and female worms?
##
### Local Ascaris vs SH
sdt.Asc%>% 
  #dplyr::group_by(System)%>%
  wilcox_test(Chao1 ~ Live)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "System")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q3_Pig_SH.csv")

sdt.Asc%>% 
  wilcox_effsize(Chao1 ~ Live)

##Plot 
sdt.Asc%>%
  ggplot(aes(x= Live, y= Chao1))+
  geom_boxplot(color= "black", alpha= 0.5)+
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= System), color= "black")+
  xlab("Worm origin")+
  ylab("Diversity (Chao1 Index)")+
  labs(tag= "A)", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -2, hide.ns = F,label = "{p.adj}{p.adj.signif}")-> J

# PCoA plot using Bray-Curtis as distance
bray_dist<- phyloseq::distance(PS3.Asc, 
                               method="bray", weighted=T)
ordination<- ordinate(PS3.Asc,
                      method="PCoA", distance=bray_dist)
plot_ordination(PS3.Asc, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(shape=21, size=3, aes(fill= System), color= "black")+
  labs(title = "Bray-Curtis dissimilarity",tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_text (x = -.20, y = 0.0, 
             label = paste ("Bray-Curtis~ System+Origin+Sex, \n PERMANOVA, System p= 0.001; R^2= 0.3189"))-> K

##Adonis PERMANOVA
#Is beta diversity vary function of the compartment, pig, infection status or technical predictors
bd.Asc<- vegan::adonis(bray_dist~ System+Live+WormSex,
                      permutations = 999, data = sdt.Asc)
bd.Asc ##Manually added to the plot 
write.csv(bd.Asc[[1]], "Tables/Q3_Permanova.csv")

###Yes, System is the most significant predictors
## Calculate multivariate dispersion (aka distance to the centroid)
mvd.Asc<- vegan::betadisper(bray_dist, sdt.Asc$System, type = "centroid")
vegan::permutest(mvd.Asc, permutations = 999)
anova(mvd.Asc)
plot(mvd.Asc)
boxplot(mvd.Asc)
plot(TukeyHSD(mvd.Asc))
###Add sample to the centroid from each sample to sdt.PA
tmp<- as.data.frame(mvd.Asc$distances)
colnames(tmp)<- c("distances")
cbind(sdt.Asc, tmp)-> sdt.Asc
###Linear model 
summary(lm(sdt.Asc, formula = distances~ System, na.action = na.exclude))

###Distances by system
sdt.Asc%>%
  wilcox_test(distances ~ System)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "System")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q3_System_distances.csv")

sdt.Asc%>%
  wilcox_effsize(distances ~ System)

##PLot 
sdt.Asc%>%
  ggplot(aes(x= System, y= distances))+
  geom_boxplot(color= "black", alpha= 0.5)+
  geom_jitter(shape=21, position=position_jitter(0.2), size=3, aes(fill= System), color= "black")+
  xlab("Worm Origin")+
  ylab("Beta diversity (distance to centroid)")+
  labs(tag= "A)", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 45))+
  stat_pvalue_manual(stats.test, hide.ns = TRUE,label = "{p.adj}{p.adj.signif}")->L

### Female vs Sex Ascaris No SH
sdt.Asc%>% 
  dplyr::filter(!System%in%("SH"))%>%
  wilcox_test(Chao1 ~ WormSex)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "WormSex")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q3_Worm_Sex.csv")

sdt.Asc%>%  
  dplyr::filter(!System%in%("SH"))%>%
  wilcox_effsize(Chao1 ~ WormSex)

##Plot 
sdt.Asc%>%
  dplyr::filter(!System%in%("SH"))%>%
  ggplot(aes(x= WormSex, y= Chao1))+
  geom_boxplot(color= "black", alpha= 0.5)+
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= System), color= "black")+
  xlab("Worm sex")+
  ylab("Diversity (Chao1 Index)")+
  labs(tag= "A)", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_pvalue_manual(stats.test, bracket.nudge.y = -2, hide.ns = F,label = "{p.adj}{p.adj.signif}")-> M

# PCoA plot using Bray-Curtis as distance
bray_dist<- phyloseq::distance(subset_samples(PS3.Asc, System!="SH"), 
                               method="bray", weighted=T)
ordination<- ordinate(subset_samples(PS3.Asc, System!="SH"),
                      method="PCoA", distance=bray_dist)
plot_ordination(subset_samples(PS3.Asc, System!="SH"), ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(shape=21, size=3, aes(fill= WormSex), color= "black")+
  labs(title = "Bray-Curtis dissimilarity",tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_text (x = 0, y = 0.2, 
             label = paste ("Bray-Curtis~ System+Sex, \n PERMANOVA, System p= 0.008; R^2= 0.2404"))-> N

##Adonis PERMANOVA
#Is beta diversity vary function of the compartment, pig, infection status or technical predictors
bd.Asc<- vegan::adonis(bray_dist~ System+WormSex,
                       permutations = 999, data = subset(sdt.Asc, System!="SH"))
bd.Asc ##Manually added to the plot 
write.csv(bd.Asc[[1]], "Tables/Q3_Permanova2.csv")

###Yes, System is the most significant predictors
## Calculate multivariate dispersion (aka distance to the centroid)
mvd.Asc<- vegan::betadisper(bray_dist, subset(sdt.Asc, System!="SH")$WormSex, type = "centroid")
vegan::permutest(mvd.Asc, permutations = 999)
anova(mvd.Asc)
plot(mvd.Asc)
boxplot(mvd.Asc)
plot(TukeyHSD(mvd.Asc))
###Add sample to the centroid from each sample to sdt.PA
tmp<- as.data.frame(mvd.Asc$distances)
colnames(tmp)<- c("distances2")
sdt.Asc%>%
  dplyr::filter(!System%in%("SH"))->sdt.Asc2
cbind(sdt.Asc2, tmp)-> sdt.Asc2
###Linear model 
summary(lm(sdt.Asc2, formula = distances~ System+WormSex, na.action = na.exclude))

###Distances by system
sdt.Asc2%>%
  wilcox_test(distances ~ WormSex)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "WormSex")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q3_Sex_distances_No_SH.csv")

sdt.Asc2%>%
  wilcox_effsize(distances ~ WormSex)

##PLot 
sdt.Asc2%>%
  ggplot(aes(x= WormSex, y= distances))+
  geom_boxplot(color= "black", alpha= 0.5)+
  geom_jitter(shape=21, position=position_jitter(0.2), size=3, aes(fill= System), color= "black")+
  xlab("System")+
  ylab("Beta diversity (distance to centroid)")+
  labs(tag= "B)", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16), axis.text.x = element_text(angle = 45))+
  stat_pvalue_manual(stats.test, hide.ns = TRUE,label = "{p.adj}{p.adj.signif}")->O


png("Figures/Q3_Alphadiv_Worm_Origin.png", units = 'in', res = 300, width=14, height=14)
grid.arrange(J, K)
dev.off()

png("Figures/Q3_Betadiv_Distances_Origin.png", units = 'in', res = 300, width=10, height=8)
grid.arrange(L)
dev.off()

png("Figures/Q3_Alphadiv_Worm_Sex.png", units = 'in', res = 300, width=14, height=14)
grid.arrange(M, N)
dev.off()

png("Figures/Q3_Betadiv_Distances_Sex.png", units = 'in', res = 300, width=10, height=8)
grid.arrange(O)
dev.off()

###Additional data --> Track microbiome in fecal samples
sdt.fec%>%
  mutate(DPI = fct_relevel(DPI, "2", "14", "21", "42", "49"))%>%
  wilcox_test(Chao1 ~ DPI)%>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()%>%
  add_xy_position(x = "DPI")-> stats.test

##Save statistical analysis
x <- stats.test
x$groups<- NULL
write.csv(x, "Tables/Q4_Track_infection_Diversity.csv")

sdt.fec%>%
  mutate(DPI = fct_relevel(DPI, "2", "14", "21", "42", "49"))%>%
  wilcox_effsize(Chao1 ~ DPI)

sdt.fec%>%
  mutate(DPI = fct_relevel(DPI, "2", "14", "21", "42", "49"))%>%
  ggplot(aes(x= DPI, y= Chao1))+
  geom_boxplot(color= "black")+
  geom_point(shape=21, position=position_jitter(0.2), size=3, aes(fill= System), color= "black")+
  xlab("Day post infection")+
  ylab("Diversity (Chao1 Index)")+
  geom_line(aes(group = System), color= "gray")+
  labs(tag= "A)", caption = get_pwc_label(stats.test))+
  theme_bw()+
  theme(text = element_text(size=16))+
  stat_pvalue_manual(stats.test, hide.ns = TRUE,label = "{p.adj}{p.adj.signif}")->P

# PCoA plot using Bray-Curtis as distance
bray_dist<- phyloseq::distance(PS3.Fec, 
                               method="bray", weighted=T)
ordination<- ordinate(PS3.Fec,
                      method="PCoA", distance=bray_dist)
plot_ordination(PS3.Fec, ordination)+ 
  theme(aspect.ratio=1)+
  geom_point(shape=21, size=3, aes(fill= DPI), color= "black")+
  labs(title = "Bray-Curtis dissimilarity",tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))+
  geom_text (x = 0, y = 0.4, 
             label = paste ("Bray-Curtis~ System+DPI, \n PERMANOVA, DPI p= 0.001; R^2= 0.3601"))->Q

##Adonis PERMANOVA
#Is beta diversity vary function of the compartment, pig, infection status or technical predictors
bd.Fec<- vegan::adonis(bray_dist~ System+DPI,
                       permutations = 999, data = sdt.fec)
bd.Fec ##Manually added to the plot 
write.csv(bd.Fec[[1]], "Tables/Q4_Permanova.csv")


png("Figures/Q4_Alphadiv_Fecal.png", units = 'in', res = 300, width=14, height=14)
grid.arrange(P, Q)
dev.off()



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

