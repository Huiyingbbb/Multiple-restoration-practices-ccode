#### community experiment ####
library(knitr)
library(tidyverse)
library(ggplot2)
library(ggdendro)
library(dendextend)
library(splitstackshape)
library(vegan)
library(caret)
library(ggpubr)

library(ggepi)
library(ggridges)
library(patchwork)
library(party)
library(dplyr)
library(randomForest)

setwd('C:/huiying data/3rd chapter/a_community')

#### first to prepare a data file for pcoa- calculate evenness and richness ####
df<- read.csv('df.csv')
df<- na.omit(df) # remove missing cells
shannon_index <- function(abundance) {
  proportions <- abundance / sum(abundance)
  proportions <- proportions[proportions > 0]  # Remove zero values
  -sum(proportions * log(proportions))
}
df$X<- NULL
species_data <- df[, 4:12] # to cover all species columns
species_richness <- apply(species_data, 1, function(row) {
  sum(row > 0)  # Count non-zero values
})

# Calculate species evenness
species_evenness <- apply(species_data, 1, function(row) {
  H <- shannon_index(row)
  S <- sum(row > 0)  # Species richness
  H / log(S)         # Evenness
})
df$evenness<- species_evenness
df$richness <- species_richness

df$plant <- rowSums(df[, 4:12]) #to calculate total biomass (plant)
df$FG_G<- rowSums(df[, 4:6])
df$FG_H<- rowSums(df[, 7:9])
df$FG_L<- rowSums(df[, 10:12])

df$G1<- NULL
df$G2<- NULL
df$G3<- NULL
df$L1<- NULL
df$L2<- NULL
df$L3<- NULL
df$H1<- NULL
df$H2<- NULL
df$H3<- NULL

####  pcoa analysis- single factors ####
df_single<- df %>%
  filter(Lv == "1")
df_single$Lv<- NULL
df_single$MWD<- NULL
df_single$respiration<- NULL
colnames(df_single)[1]<- 'group'
df_single$richness<- as.numeric(df_single$richness)

## need to change here while switching between soil, plant or all response variables ##
columns_to_summarize <- colnames(df_single)[2:11]
standardized_sig_facts <- df_single %>%
  group_by(group) %>%
  summarise(across(all_of(columns_to_summarize), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}")) %>%
  as.data.frame() %>%
  column_to_rownames(var = "group") %>%
  scale()

single2<- data.frame(standardized_sig_facts)
a_single2<- data.frame(t(single2[-1]))

dist<- vegdist(t(a_single2),method='euclidean',diag=T,upper=T)
dist<- as.matrix(dist)

pcoa<- cmdscale(dist,eig=T)
eig<- summary(eigenvals(pcoa))
aixs<- paste0('PCoA',1:ncol(eig))
eig<- data.frame(aixs,t(eig)[,-3])
pco1<- round(eig[1,3]*100,3)
pco2<- round(eig[2,3]*100,3)
xlab<- paste0('PCoA1(',pco1,'%)')
ylab<- paste0('PCoA2(',pco2,'%)')
pcoa_points<- as.data.frame(pcoa$points)
pcoa_points<- data.frame(pcoa_points,group=colnames(a_single2))

pall<-ggplot(pcoa_points,aes(V1,V2)) +
  geom_point(size=4,aes(color=group))+
  geom_text(aes(label=rownames(pcoa_points)), size=4, vjust=1, hjust=0)+
  labs(x=xlab,y=ylab)+ 
  theme(plot.title = element_text(size = 12,hjust = 0.5))+  
  geom_hline(yintercept=0, linetype=4) +        
  geom_vline(xintercept=0 ,linetype=4)+ 
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

## only consider soil ##
columns_to_summarize <- colnames(df_single)[5:11]
standardized_sig_facts <- df_single %>%
  group_by(group) %>%
  summarise(across(all_of(columns_to_summarize), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}")) %>%
  as.data.frame() %>%
  column_to_rownames(var = "group") %>%
  scale()

single2<- data.frame(standardized_sig_facts)
a_single2<- data.frame(t(single2[-1]))

#df_single2<- data.frame(t(df_single[-1]))
#colnames(df_single2)<- df_single[,1]

dist<- vegdist(t(a_single2),method='euclidean',diag=T,upper=T)
dist<- as.matrix(dist)

pcoa<- cmdscale(dist,eig=T)
eig<- summary(eigenvals(pcoa))
aixs<- paste0('PCoA',1:ncol(eig))
eig<- data.frame(aixs,t(eig)[,-3])
pco1<- round(eig[1,3]*100,3)
pco2<- round(eig[2,3]*100,3)
xlab<- paste0('PCoA1(',pco1,'%)')
ylab<- paste0('PCoA2(',pco2,'%)')
pcoa_points<- as.data.frame(pcoa$points)
pcoa_points<- data.frame(pcoa_points,group=colnames(a_single2))

ggplot(pcoa_points,aes(V1,V2)) +
  geom_point(size=4,aes(color=group))+
  geom_text(aes(label=rownames(pcoa_points)), size=4, vjust=1, hjust=0)+
  labs(x=xlab,y=ylab)+ 
  theme(plot.title = element_text(size = 12,hjust = 0.5))+  
  geom_hline(yintercept=0, linetype=4) +        
  geom_vline(xintercept=0 ,linetype=4)+ 
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
