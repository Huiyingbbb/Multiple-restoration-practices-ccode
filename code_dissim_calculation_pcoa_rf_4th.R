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

#### first to prepare a clean data file fpr pcoa ####
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
  #geom_polygon(data = group_border, aes(fill = group),alpha=0.6)+
  #填充颜色，设置颜色透明度  
  geom_text(aes(label=rownames(pcoa_points)), size=4, vjust=1, hjust=0)+
  #添加样本名称,名称位置大小  
  labs(x=xlab,y=ylab)+ 
  #添加x、y轴标题和图标题 
  #stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
  #添加置信椭圆  
  theme(plot.title = element_text(size = 12,hjust = 0.5))+  
  geom_hline(yintercept=0, linetype=4) +        
  #添加原点垂直辅助线  
  geom_vline(xintercept=0 ,linetype=4)+ 
  #添加原点水平辅助线  
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
  #geom_polygon(data = group_border, aes(fill = group),alpha=0.6)+
  #填充颜色，设置颜色透明度  
  geom_text(aes(label=rownames(pcoa_points)), size=4, vjust=1, hjust=0)+
  #添加样本名称,名称位置大小  
  labs(x=xlab,y=ylab)+ 
  #添加x、y轴标题和图标题 
  #stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
  #添加置信椭圆  
  theme(plot.title = element_text(size = 12,hjust = 0.5))+  
  geom_hline(yintercept=0, linetype=4) +        
  #添加原点垂直辅助线  
  geom_vline(xintercept=0 ,linetype=4)+ 
  #添加原点水平辅助线  
  theme_bw()+  
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())


#### dissimilarity calculation ####
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
df_single<- df %>%
  filter(lv == "1")
df_single$lv<- NULL
df_single$MWD<- NULL
df_single$respiration<- NULL

columns_to_summarize <- colnames(df_single)[2:11]

standardized_sig_facts <- df_single %>%
  group_by(treatment) %>%
  summarise(across(all_of(columns_to_summarize), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}")) %>%
  as.data.frame() %>%
  column_to_rownames(var = "treatment") %>%
  scale()


bray_dis<- vegdist(standardized_sig_facts,method = "euclidean")
dis<- as.matrix(bray_dis)

df_2<-df%>%
  filter(lv=='2')
f3<-df%>%
  filter(lv=='3')
f4<-df%>%
  filter(lv=='4')
f5<-df%>%
  filter(lv=='5')

# calculate for 2 factors #
f2_mut<-df_2%>%
  separate(treatment, into = c("label_1", "label_2"), sep = "-")

a<-f2_mut$label_1
b<-f2_mut$label_2

N<-1:84 #you have 84 2 factor level units
c = vector()
for (i in N) {
  c[i]<-dis[a[i],b[i]]
}
f2_mut<-f2_mut%>%
  mutate(dissim = c)

# center/descale the dissimilarity index # 
f2_dis<-as.data.frame(f2_mut$dissim)
colnames(f2_dis)<-c("dissim")
pro_f2<-caret::preProcess(f2_dis, method = c("range"))
f2_dis<-predict(pro_f2, f2_dis)
c<-vector()
c<-f2_dis$dissim
df_2<-df_2%>%
  mutate(dissim = c)

# calculate for 3 factors #

f3_mut<-f3%>%
  separate(treatment,c("label_1", "label_2",'label_3'), sep = "-")

a<-f3_mut$label_1
b<-f3_mut$label_2
c<-f3_mut$label_3

X_1 = vector()
X_2 = vector()
X_3 = vector()

Y = vector()
N<-1:10
for (i in N) {
  X_1[i]<-dis[a[i],b[i]]
  X_2[i]<-dis[a[i],c[i]]
  X_3[i]<-dis[b[i],c[i]]
  Y[i]<-X_1[i]+X_2[i]+X_3[i]
}
f3_mut<-f3_mut%>%
  mutate(dissim = Y)

# center/descale the dissimilarity index # 
f3_dis<-as.data.frame(f3_mut$dissim)
colnames(f3_dis)<-c("dissim")
pro_f3<-caret::preProcess(f3_dis, method = c("range"))
f3_dis<-predict(pro_f3, f3_dis)
c<-vector()
c<-f3_dis$dissim
f3<-f3%>%
  mutate(dissim = c)

# calculate for 4 factor #
f4_mut<-f4%>%
  separate(treatment,c("label_1", "label_2",'label_3','label_4'), sep = "-")

a<-f4_mut$label_1
b<-f4_mut$label_2
c<-f4_mut$label_3
d<-f4_mut$label_4

X_1 = vector()
X_2 = vector()
X_3 = vector()
X_4 = vector()
X_5 = vector()
X_6 = vector()

Y = vector()
N<-1:10
for (i in N) {
  X_1[i]<-dis[a[i],b[i]]
  X_2[i]<-dis[a[i],c[i]]
  X_3[i]<-dis[a[i],d[i]]
  X_4[i]<-dis[b[i],c[i]]
  X_5[i]<-dis[b[i],d[i]]
  X_6[i]<-dis[c[i],d[i]]
  Y[i]<-X_1[i]+X_2[i]+X_3[i]+X_4[i]+X_5[i]+X_6[i]
}
f4_mut<-f4_mut%>%
  mutate(dissim = Y)

# center/descale the dissimilarity index # 
f4_dis<-as.data.frame(f4_mut$dissim)
colnames(f4_dis)<-c("dissim")
pro_f4<-caret::preProcess(f4_dis, method = c("range"))
f4_dis<-predict(pro_f4, f4_dis)
c<-vector()
c<-f4_dis$dissim
f4<-f4%>%
  mutate(dissim = c)

# calculate for 5 factor #
f5_mut<-f5%>%
  separate(treatment,c("label_1", "label_2",'label_3','label_4','label_5'), sep = "-")

a<-f5_mut$label_1
b<-f5_mut$label_2
c<-f5_mut$label_3
d<-f5_mut$label_4
e<-f5_mut$label_5

X_1 = vector()
X_2 = vector()
X_3 = vector()
X_4 = vector()
X_5 = vector()
X_6 = vector()
X_7 = vector()
X_8 = vector()
X_9 = vector()
X_10 = vector()
Y = vector()
N<-1:10
for (i in N) {
  X_1[i]<-dis[a[i],b[i]]
  X_2[i]<-dis[a[i],c[i]]
  X_3[i]<-dis[a[i],d[i]]
  X_4[i]<-dis[a[i],e[i]]
  X_5[i]<-dis[b[i],c[i]]
  X_6[i]<-dis[b[i],d[i]]
  X_7[i]<-dis[b[i],e[i]]
  X_8[i]<-dis[c[i],d[i]]
  X_9[i]<-dis[c[i],e[i]]
  X_10[i]<-dis[d[i],e[i]]
  Y[i]<-X_1[i]+X_2[i]+X_3[i]+X_4[i]+X_5[i]+X_6[i]+X_7[i]+X_8[i]+X_9[i]+X_10[i]
}
f5_mut<-f5_mut%>%
  mutate(dissim = Y)

# center/descale the dissimilarity index # 
f5_dis<-as.data.frame(f5_mut$dissim)
colnames(f5_dis)<-c("dissim")
pro_f5<-caret::preProcess(f5_dis, method = c("range"))
f5_dis<-predict(pro_f5, f5_dis)
c<-vector()
c<-f5_dis$dissim
f5<-f5%>%
  mutate(dissim = c)

# data cleaning #
df_2$label_1<- NULL
df_2$label_2<- NULL
f3$label_1<- NULL
f3$label_2<- NULL
f3$label_3<- NULL
f4$label_1<- NULL
f4$label_2<- NULL
f4$label_3<- NULL
f4$label_4<- NULL
f5$label_1<- NULL
f5$label_2<- NULL
f5$label_3<- NULL
f5$label_4<- NULL
f5$label_5<- NULL
all_dis<- rbind(df_2,f3)
all_dis<- rbind(all_dis,f4)
all_dis<- rbind(all_dis,f5)
write.csv(all_dis,'all_dissimilarity.csv')

#### regression on level and dissimilarity ####
df<- read.csv('df_rf.csv')
f2<- df %>%
  filter(df$lv=='2')
df_lv<-df %>%
  filter(remark!='1')
df_all<- df %>%
  filter(lv!='1'&lv!='0')
dis1<- ggplot(f2,aes(x=dissim,y=FG_G))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  theme_bw()+
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank() # Remove x-axis title
  )
all_dis1<- ggplot(df_all,aes(x=dissim,y=FG_G))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  theme_bw()+
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank() # Remove x-axis title
  )
  
  
lv1<- ggplot(df_lv,aes(x=lv,y=FG_G))+
  geom_point(size = 1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = '#5c7896', fill = '#5c7896', alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )

dis2<- ggplot(f2,aes(x=dissim,y=FG_H))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )

all_dis2<- ggplot(df_all,aes(x=dissim,y=FG_H))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )
lv2<- ggplot(df_lv,aes(x=lv,y=FG_H))+
  geom_point(size = 1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = '#5c7896', fill = '#5c7896', alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )

dis3<- ggplot(f2,aes(x=dissim,y=FG_L))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )

all_dis3<- ggplot(df_all,aes(x=dissim,y=FG_L))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )
lv3<- ggplot(df_lv,aes(x=lv,y=FG_L))+
  geom_point(size = 1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = '#5c7896', fill = '#5c7896', alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )

dis4<- ggplot(f2,aes(x=dissim,y=evenness))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )

all_dis4<- ggplot(df_all,aes(x=dissim,y=evenness))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )
lv4<- ggplot(df_lv,aes(x=lv,y=evenness))+
  geom_point(size = 1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = '#5c7896', fill = '#5c7896', alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )
dis5<- ggplot(f2,aes(x=dissim,y=richness))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )
all_dis5<- ggplot(df_all,aes(x=dissim,y=richness))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )

lv5<- ggplot(df_lv,aes(x=lv,y=richness))+
  geom_point(size = 1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = '#5c7896', fill = '#5c7896', alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )
dis6<- ggplot(f2,aes(x=dissim,y=root))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )
all_dis6<- ggplot(df_all,aes(x=dissim,y=root))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )
lv6<- ggplot(df_lv,aes(x=lv,y=root))+
  geom_point(size = 1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = '#5c7896', fill = '#5c7896', alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )

dis7<- ggplot(f2,aes(x=dissim,y=plant))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )
all_dis7<- ggplot(df_all,aes(x=dissim,y=plant))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )
lv7<- ggplot(df_lv,aes(x=lv,y=plant))+
  geom_point(size = 1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = '#5c7896', fill = '#5c7896', alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )


dis8<- ggplot(f2,aes(x=dissim,y=WSA))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )
all_dis8<- ggplot(df_all,aes(x=dissim,y=WSA))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )
lv8<- ggplot(df_lv,aes(x=lv,y=WSA))+
  geom_point(size = 1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = '#5c7896', fill = '#5c7896', alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )

dis9<- ggplot(f2,aes(x=dissim,y=WHC))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )

all_dis9<-ggplot(df_all,aes(x=dissim,y=WHC))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )
lv9<- ggplot(df_lv,aes(x=lv,y=WHC))+
  geom_point(size = 1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = '#5c7896', fill = '#5c7896', alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )

dis10<- ggplot(f2,aes(x=dissim,y=pH))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )
all_dis10<- ggplot(df_all,aes(x=dissim,y=pH))+
  geom_point(size =1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = "#BAA43F", fill = "#BAA43F", alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )
lv10<- ggplot(df_lv,aes(x=lv,y=pH))+
  geom_point(size = 1.6,alpha=0.1)+
  stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
           method = "spearman", 
           label.x = Inf, label.y = Inf,
           hjust = 1.1, vjust = 1.5,size=3.5) +
  geom_smooth(method = "lm",color = '#5c7896', fill = '#5c7896', alpha = 0.3)+
  theme_bw()+
  theme(
    axis.title.y = element_blank(), # Remove y-axis title
    axis.title.x = element_blank()
  )
#### random forest analysis ####
df<- read.csv('df_rf.csv')

df$remark=factor(df$remark, levels = unique(df$remark))
levels = c("1", "2", "3", "4",'5')
stressors = c("A","Be","Bi","P","St","Si","V")
responses = c("FG_G","FG_H",'FG_L',"evenness","richness","root","plant",'WSA','WHC','pH')
treatment = as.vector(unique(df$remark))
n_iter = 100

# Estimating mean and its 95% confidence interval

BootStrap_mean = function(response, data=df, target = treatment, n_perm = n_iter){
  summary = list()
  
  for(treatment in target){
    bs = numeric(0)
    if(treatment=="1") population = data[data$remark%in%stressors, response]
    if(treatment!="1") population = data[data$remark==treatment, response]
    size = length(population)-sum(is.na(population))
    
    for(id in c(1:n_perm)){
      k = mean(sample(population, size, replace = T), na.rm = TRUE)
      bs = append(bs, k)
    }
    summary[[treatment]] = c(quantile(bs, .025,na.rm = TRUE), mean(bs,na.rm = TRUE), quantile(bs, .975,na.rm = TRUE))
    names(summary[[treatment]]) = c("2.5%", "mean", "97.5%")
  }
  summary = t(data.frame(summary))
  summary = data.frame("target" = target, summary); row.names(summary) = c()
  return(summary)
}
BootStrap_ES_rep = function(response, data=df, target = treatment, n_perm = n_iter){
  resampled = list()
  
  population_CT = data[data$remark=="CT", response]
  
  for(treatment in target){
    bs = numeric(0)
    if(treatment=="1") population_TR = data[data$remark%in%stressors, response]
    if(treatment!="1") population_TR = data[data$remark==treatment, response]
    size_CT = length(population_CT)-sum(is.na(population_CT))
    size_TR = length(population_TR)-sum(is.na(population_TR))
    
    for(id in c(1:n_perm)){
      k_CT = mean(sample(population_CT, size_CT, replace = T), na.rm = TRUE)
      k_TR = mean(sample(population_TR, size_TR, replace = T), na.rm = TRUE)
      bs = append(bs, k_TR - k_CT)
    }
    resampled[[treatment]] = bs
  }
  resampled[["CT"]] = rep(0, n_perm)
  return(resampled)
}
BootStrap_ES_summary = function(data){
  summary = list()
  p = 0
  summary[["CT"]] = c(0,0,0,1)
  target = names(data)
  
  for(treatment in target[-1]){
    bs = data[[treatment]]
    p = length(which(bs>0))/length(bs)
    p = min(p, 1-p)
    summary[[treatment]] = c(quantile(bs, .025,na.rm = TRUE), mean(bs,na.rm = TRUE), quantile(bs, .975,na.rm = TRUE), p)
  }
  summary = t(data.frame(summary))
  colnames(summary) = c("2.5%", "mean", "97.5%", "p_value")
  summary = data.frame(target, summary); row.names(summary) = c()
  
  return(summary)
}
Null_distribution_rep = function(response, data=df, n_perm=n_iter){
  
  output = list()
  for(Lv in levels){
    
    resampled = list()
    
    # Checking which stressor combinations were jointly tested
    if(Lv=="1") combination = data[data$remark%in%stressors,c(1:7)]##changed the selected columns
    if(Lv!="1") combination = data[data["remark"]==Lv,c(1:7)]##changed the selected columns
    Level = sum(combination[1,]) 
    
    
    # Null distributions can be taken based on three different assumptions
    for(type in c("Additive", "Multiplicative", "Dominative")){
      
      population_CT = df[df$remark=="CT", response]
      size_CT = length(population_CT)-sum(is.na(population_CT)) ##subtract NA value
      
      # For each combination, bootstrap resampling is conducted
      for(j in c(1:nrow(combination))){
        bs = numeric(0)
        selected_stressors = stressors[which(combination[j,]==1)]
        sub_n_perm = ceiling(n_perm/nrow(combination))*5 #*5 increase the permutation number for sub-sampling
        
        # bootstrap resampling
        for(id in c(1:sub_n_perm)){
          each_effect = numeric(0)
          k_CT = mean(sample(population_CT, size_CT, replace = T),na.rm = TRUE) #ignore NA
          
          for(treatment in selected_stressors){
            population_TR = df[df$remark==treatment, response]
            size_TR = length(population_TR)-sum(is.na(population_TR))
            k_TR = mean(sample(population_TR, size_TR, replace = T),na.rm = TRUE)#ignore NA
            
            # ES estimate depending on the type of null hypotheses
            if(type=="Additive")       each_effect = append(each_effect, (k_TR - k_CT))
            if(type=="Multiplicative") each_effect = append(each_effect, (k_TR - k_CT)/k_CT)
            if(type=="Dominative")      each_effect = append(each_effect, (k_TR - k_CT))
          }
          
          # Calculating an expected ES after collecting the ESs of all relevant single stressors
          if(type=="Additive")       joint_effect = sum(each_effect)
          if(type=="Multiplicative"){
            z = 1
            for(m in c(1:Level)) z = z * (1 + each_effect[m])
            joint_effect = (z - 1)*k_CT
          }
          if(type=="Dominative")      joint_effect = each_effect[which(max(abs(each_effect))==abs(each_effect))]
          
          bs = append(bs, joint_effect)
        }
        resampled[[type]][[j]] = bs
      }
      
    }
    output[[Lv]] = resampled
  }  
  return(output)
}
Null_distribution_rep_transform = function(data){
  output = list()
  for(Lv in levels){
    for(type in c("Additive", "Multiplicative", "Dominative")){
      output[[Lv]][[type]] = sample(unlist(data[[Lv]][[type]]), n_iter, replace=F)
    }
  }
  return(output)
}
NHST_summary = function(null_data, Actual_data){
  output = list()
  for(Lv in levels){
    summary = list()
    summary[["Actual"]] = c(quantile(Actual_data[[Lv]], .025,na.rm = TRUE), mean(Actual_data[[Lv]],na.rm = TRUE), quantile(Actual_data[[Lv]], .975,na.rm = TRUE), 1)
    p = 0
    assumptions = c("Additive", "Multiplicative", "Dominative")
    
    for(i_assumption in assumptions){
      bs   = (Actual_data[[Lv]] - null_data[[Lv]][[i_assumption]])
      p = length(which(bs>0))/length(bs)
      p = min(p, 1-p)
      summary[[i_assumption]] = c(quantile(null_data[[Lv]][[i_assumption]], .025,na.rm = TRUE), mean(null_data[[Lv]][[i_assumption]],na.rm = TRUE), quantile(null_data[[Lv]][[i_assumption]], .975,na.rm = TRUE), p)
    }
    summary = t(data.frame(summary))
    colnames(summary) = c("2.5%", "mean", "97.5%", "p_value")
    summary = data.frame(ES = c("Actual", "Additive","Multiplicative","Dominative"), summary); row.names(summary) = c()
    
    output[[Lv]] = summary
  }
  
  return(output)
}

NHST_summary_transform = function(data){
  output = list()
  for(i in 1:4){
    summary = rbind(data[["1"]][i, 2:4], data[["2"]][i, 2:4], data[["3"]][i, 2:4],
                    data[["4"]][i, 2:4], data[["5"]][i, 2:4])
    summary = cbind(levels, summary)
    colnames(summary) = c("Lv", "Low", "Mean", "High")
    output[[c("Actual", "Additive", "Multiplicative", "Dominative")[i]]] = summary
  }
  return(output)
}
Expected_ES_for_each = function(data){
  output = numeric(0)
  for(type in c("Additive", "Multiplicative", "Dominative")){
    tmp = numeric(0)
    for(Lv in levels){
      n_len = length(data[[Lv]][[type]])
      for(i in 1:n_len){
        tmp = append(tmp, mean(data[[Lv]][[type]][[i]]))
      }
    }
    output = cbind(output,tmp)
  }
  colnames(output)= c("E1", "E2", "E3")
  return(output)
}
response_mean_all = list()
response_ES_all   = list()
joint_ES_null_all = list()
ES_for_each_all   = list()
g_rawdata_all     = list()
g_ms_all          = list()
for (i_response in responses) {
  #bootstrap mean
  response_mean = BootStrap_mean(i_response)
  response_mean_all[[i_response]] = response_mean
  #effect size
  response_ES_bs  = BootStrap_ES_rep(i_response)
  response_ES     = BootStrap_ES_summary(response_ES_bs)
  response_ES_all[[i_response]]   = response_ES
  #Null model
  Null_ES_bs0     = Null_distribution_rep(i_response)
  Null_ES_bs      = Null_distribution_rep_transform(Null_ES_bs0)
  
  joint_ES_null   = NHST_summary(Null_ES_bs, response_ES_bs)
  joint_ES_null_all[[i_response]] = bind_rows(joint_ES_null, .id = "column_label")
  
  ES_plot         = NHST_summary_transform(joint_ES_null)
  ES_for_each     = Expected_ES_for_each(Null_ES_bs0)
  ES_for_each_all[[i_response]]   = ES_for_each
  
  g_rawdata_all[[i_response]] =  local({
    i_response = i_response
    response_mean = response_mean
    Mycolor=c("#9EA9A1","#7a8e8f","#7a8e8f","#7a8e8f","#7a8e8f","#7a8e8f","#7a8e8f","#7a8e8f","#ACBDCF","#8ca5c0","#7990A9",'#5c7896',"#4e5b6d")
    
    ggplot() +
      xlab(i_response) +
      coord_flip() +
      theme_bw()+
      theme(legend.position = 'none', axis.title.x=element_blank())+
     stat_density_ridges(data=df, aes_string(x = i_response, y = "remark",color= "remark",fill= "remark"),
                          geom = "density_ridges_gradient",
                          rel_min_height = 0.01, 
                          jittered_points = TRUE, 
                          alpha = .5,
                          position = position_points_jitter(height = .2, yoffset = .15),
                          point_size = 1, point_alpha = .3, 
                          scale = .5) +
      scale_fill_manual(values  = Mycolor)+
      scale_alpha_manual(values = c(rep(0.5,17)))+
      geom_estci(data=response_mean, aes(x = mean, y = target, xmin=X2.5., xmax=X97.5., 
                                         xintercept=response_mean[1,"mean"]), center.linecolour = "black",
                 size=0.4, ci.linesize = 0.5, position=position_nudge(y = -0.15)) +
      scale_color_manual(values  = Mycolor)
    
    
  })
  
  # Ploting the number of stressors and ES relationship
  
  g_ms_all[[i_response]] = local({
    i_response = i_response
    ct_value = response_mean[1,"mean"]
    ES_plot = ES_plot
    for(i in 1:4) ES_plot[[i]][,2:4] = ES_plot[[i]][,2:4] + ct_value
    
    ggplot()+
      theme_bw()+
      theme(legend.position = 'none', axis.title.x=element_blank(), axis.title.y=element_blank())+
      coord_flip() +
      scale_y_discrete(limits = factor(c("1","2","3","4","5"), levels=c("1","2","3","4","5"))) +
      
      ### Mean & CI ###
      # 3 assumptions
      geom_estci(data=ES_plot[["Additive"]], aes(x = Mean, y = Lv, xmin=Low, xmax=High, xintercept=ct_value), 
                 color="#9BB9B9", size=0.4, ci.linesize = 0.4, position=position_nudge(y = +0.1)) +
      geom_estci(data=ES_plot[["Multiplicative"]], aes(x = Mean, y = Lv, xmin=Low, xmax=High, xintercept=ct_value), 
                 color="#B5CBE2", size=0.4, ci.linesize = 0.4, position=position_nudge(y = +0.2)) +
      geom_estci(data=ES_plot[["Dominative"]], aes(x = Mean, y = Lv, xmin=Low, xmax=High, xintercept=ct_value), 
                 color="#455D99", size=0.4, ci.linesize = 0.4, position=position_nudge(y = +0.3)) +
      
      # Actual ES
      geom_estci(data=ES_plot[["Actual"]], aes(x = Mean, y = Lv, xmin=Low, xmax=High, xintercept=ct_value), 
                 size=0.6, ci.linesize = 0.6, position=position_nudge(y = 0))
  })
  
}
postResample <- function(pred, obs)
{
  isNA <- is.na(pred)
  pred <- pred[!isNA]
  obs <- obs[!isNA]
  if (!is.factor(obs) && is.numeric(obs))
  {
    if(length(obs) + length(pred) == 0)
    {
      out <- rep(NA, 3)
    } else {
      if(length(unique(pred)) < 2 || length(unique(obs)) < 2)
      {
        resamplCor <- NA
      } else {
        resamplCor <- try(cor(pred, obs, use = "pairwise.complete.obs"), silent = TRUE)
        if (inherits(resamplCor, "try-error")) resamplCor <- NA
      }
      mse <- mean((pred - obs)^2)
      mae <- mean(abs(pred - obs))
      out <- c(sqrt(mse), resamplCor^2, mae)
    }
    names(out) <- c("RMSE", "Rsquared", "MAE")
  } else {
    if(length(obs) + length(pred) == 0)
    {
      out <- rep(NA, 2)
    } else {
      pred <- factor(pred, levels = levels(obs))
      requireNamespaceQuietStop("e1071")
      out <- unlist(e1071::classAgreement(table(obs, pred)))[c("diag", "kappa")]
    }
    names(out) <- c("Accuracy", "Kappa")
  }
  if(any(is.nan(out))) out[is.nan(out)] <- NA
  out
}


df.rf = df[df[, "remark"] %in% levels,]

lv_list  = unique(df.rf[,"lv"])
id_lv1   = which(df.rf[,"lv"]==1)
id_lv2   = which(df.rf[,"lv"]==2)
id_lvh   = which(df.rf[,"lv"]>2)

n_data   = nrow(df.rf)
n_lv1    = sum(df.rf[,"lv"]==1)
n_lvh    = sum(df.rf[,"lv"]!=1)
n_lv2<- 10
n_eachlv = 10 #changed it to 10 for balanced resampling. Our experiment has 10 replicates for high levels.
n_tree   = 50
n_iter2  = 50

rf.r2 = data.frame(matrix(NA, ncol=3, nrow=4*n_iter2*length(responses)))
rf.r2[,1] = rep(responses,each=4*n_iter2)
rf.r2[,2] = rep(c("ES","lv+ES","lv+ES+Dis","All"), n_iter2*length(responses))
rf.r2[,2] = factor(rf.r2[,2], levels=c("ES","lv+ES","lv+ES+Dis","All"))
colnames(rf.r2) = c("Response", "Model", "R2")

rf.pred = data.frame(matrix(NA, ncol=6, nrow=n_iter2*length(responses)))
rf.pred[,1] = rep(responses,each=n_iter2)
colnames(rf.pred) = c("Response", lv_list)

rf.vimp = data.frame(matrix(NA, ncol=13, nrow=n_iter2*length(responses)))
rf.vimp[,1] = rep(responses,each=n_iter2) 
colnames(rf.vimp) = c("Response", "Lv", stressors, 'dissim',"E1","E2","E3")

rf.prediction.all = list()

j = 0
for(i_response in responses){
  
  df.rf.tmp = cbind(df.rf, ES_for_each_all[[i_response]])
  #define the formulas for three random forest models
  
  eval(parse(text=(paste("fml      = formula(",i_response,"~", paste(c("lv","E1","dissim", stressors), collapse=" + "),")", sep=""))))
  eval(parse(text=(paste("fml.ES   = formula(",i_response,"~E1)", sep=""))))
  eval(parse(text=(paste("fml.LvES = formula(",i_response,"~", paste(c("lv","E1"), collapse=" + "),")", sep=""))))
  eval(parse(text=(paste("fml.LvESDis = formula(",i_response,"~", paste(c("lv","E1","dissim"), collapse=" + "),")", sep=""))))
  for(i in 1:n_iter2){
    j = j + 1
    # bootstrap resampling
    set.seed(j)
    # take 10 sample from Lv1, 10 from Lv2, and take 30 sample from the other levels for balanced resampling
    rid = c(sample(id_lv1, n_eachlv, replace=T), sample(id_lv2, n_lv2, replace=T),sample(id_lvh, n_lvh, replace=T))
    rdf = df.rf.tmp[rid,]
    rf_model      = tryCatch({
      cforest(fml,      data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))},
      error = function(e) {cforest(fml.ES,   data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))})
    rf_model.ES   = cforest(fml.ES,   data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
    rf_model.LvES = cforest(fml.LvES, data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
    rf_model.LvESDis = cforest(fml.LvESDis, data = rdf, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
    
    # shuffling corresponding variables for evaluating R2 reduction
    # rdf.test.LvID = rdf
    # rdf.test.LvID[,colnames(ES_for_each_all[[i_response]])] = apply(rdf.test.LvID[,colnames(ES_for_each_all[[i_response]])],2,sample)  
    
    # evaluating fitting performance
    
    rf_prdct.ES   = tryCatch({predict(rf_model.ES, OOB=T)}, error=function(e){rep(0,length(rid))})
    rf_prdct.LvES = tryCatch({predict(rf_model.LvES, OOB=T)}, error=function(e){rep(0,length(rid))})
    rf_prdct.LvESDis = tryCatch({predict(rf_model.LvESDis, OOB=T)}, error=function(e){rep(0,length(rid))})
    rf_prdct.All  = tryCatch({predict(rf_model, OOB=T)}, error=function(e){rep(0,length(rid))})
    
    rf.r2[4*j-3,3]    = postResample(rf_prdct.ES,rdf[,i_response])[2]
    rf.r2[4*j-2,3]    = postResample(rf_prdct.LvES,rdf[,i_response])[2]
    rf.r2[4*j-1,3]    = postResample(rf_prdct.LvESDis,rdf[,i_response])[2]
    rf.r2[4*j-0,3]    = postResample(rf_prdct.All,rdf[,i_response])[2]
    
    
    # variable importance
    #set.seed(j)
    # tmp.vimp = numeric(0)
    # for(itmp in 1:5) tmp.vimp = rbind(tmp.vimp,varimp(rf_model))
    # tmp.vimp = apply(tmp.vimp,2,mean)
    #tmp.vimp = tmp.vimp/sum(tmp.vimp)
    # rf.vimp[j, 2:(length(tmp.vimp)+1)] = tmp.vimp*rf.r2[2*j-0,3]*100
    
    # fitting curve
    #tmp.curve = c()
    #for(i_lv in lv_list){
    #  tmp.curve = append(tmp.curve, mean(rf_prdct.All[rdf[,"lv"]==i_lv]))
    #}
    #rf.pred[j,2:5]    =  tmp.curve
  }
  
  rf_model      = cforest(fml,      data = df.rf.tmp, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
  rf_model.ES   = cforest(fml.ES,   data = df.rf.tmp, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
  rf_model.LvES = cforest(fml.LvES, data = df.rf.tmp, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
  rf_model.LvESDis = cforest(fml.LvESDis, data = df.rf.tmp, control=cforest_control(ntree=n_tree, minsplit = 5, minbucket=2))
  
  rf.prediction.all[[i_response]] = data.frame(
    Model     = rep(c("ES","lv+ES","lv+ES+Dis","All"), each = nrow(df.rf.tmp)),
    Predicted = c(predict(rf_model.ES, OOB=F),predict(rf_model.LvES, OOB=F),predict(rf_model.LvESDis, OOB=F),predict(rf_model, OOB=F)),
    Observed  = rep(df.rf.tmp[,i_response], 4))
}

rf.r2.summary = data.frame(matrix(NA,ncol=5,nrow=4*length(responses)))
colnames(rf.r2.summary) = c("Response","Model", "CI.low", "Mean", "CI.high")
rf.r2.summary[,1] = rep(responses,each=4)
rf.r2.summary[,2] = rep(c("ES","lv+ES","lv+ES+Dis","All"),length(responses))
rf.r2.summary[,2] = factor(rf.r2.summary[,2],levels=c("ES","lv+ES","lv+ES+Dis","All"))


rf.pred.summary = data.frame(matrix(NA,ncol=5,nrow=5*length(responses)))#5 is level here
colnames(rf.pred.summary) = c("Response","lv","CI.low", "Mean", "CI.high")
rf.pred.summary[,1] = rep(responses,each=length(lv_list))
rf.pred.summary[,2] = rep(lv_list, length(responses))

#rf.vimp.summary = data.frame(matrix(NA,ncol=5,nrow=13*length(responses))) #change here nr*length here, see rf.vimp col before
#colnames(rf.vimp.summary) = c("Response","Variable","CI.low", "Mean", "CI.high")
#rf.vimp.summary[,1] = rep(responses,each=13)
#rf.vimp.summary[,2] = rep(colnames(rf.vimp)[2:13], length(responses))
#rf.vimp.summary[,2] = factor(rf.vimp.summary[,2], levels = unique(rf.vimp.summary[,2]))
j=0
for(i_response in responses){
  jj = 0
  jjj = 0
  rf.r2.summary[4*j+1,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="ES"),   3],c(.025,.50,.975), na.rm=T)
  rf.r2.summary[4*j+2,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="lv+ES"),   3],c(.025,.50,.975), na.rm=T)
  rf.r2.summary[4*j+3,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="lv+ES+Dis"),3],c(.025,.50,.975), na.rm=T)
  rf.r2.summary[4*j+4,3:5] = quantile(rf.r2[which(rf.r2[,1]==i_response&rf.r2[,2]=="All"),  3],c(.025,.50,.975), na.rm=T)
  for(i_lv in lv_list){
    jj = jj + 1
    rf.pred.summary[5*j+jj,3:5]     = quantile(rf.pred[which(rf.pred[,1]==i_response),as.character(i_lv)],c(.025,.50,.975), na.rm=T)
  }
  
  
  # for(i_var in colnames(rf.vimp)[2:11]){
  #  jjj = jjj + 1
  #   rf.vimp.summary[10*j+jjj,3:5]     = quantile(rf.vimp[which(rf.vimp[,1]==i_response),as.character(i_var)],c(.025,.50,.975), na.rm=T)
  # }#change 14-7*j+jjj,3:5
  
  j = j + 1
  
}
g_rf_all     = list()
g_vimp_all   = list()
g_r2_all     = list()

for(i in 1:(length(responses))){
  
  # g_vimp_all[[i]] = local({
  #  i = i
  # ggplot(data=rf.vimp.summary[rf.vimp.summary[,2]==responses[i],],aes(x=factor(Variable, level = c("Lv","ES","Dis","Bi","Co","Cl","S","Ba","M","PD",'OM')), y=Mean))+
  #  geom_bar(stat="identity", fill="#999999", alpha=0.5) +
  # xlab("Variability explained [%]")+
  #geom_errorbar(aes(ymax=CI.high, ymin=CI.low), width=.2, position=position_dodge(width=0.0)) +
  #  theme_bw()+
  # theme(legend.position = 'none', axis.title.x=element_blank(),axis.title.y=element_blank())
  #  })
  
  g_r2_all[[i]] = local({
    i = i
    rf.r2 = rf.r2
    rf.r2.summary = rf.r2.summary 
    level_order=c("ES","lv+ES","lv+ES+Dis","All")
    ggplot(data=rf.r2[rf.r2[,1]==responses[i],],aes(x=Model,y=R2, fill=Model))+
      geom_violin( color="#00000000",alpha=.5,position=position_dodge(width=0.3),trim=F)+
      geom_pointrange(data=rf.r2.summary[rf.r2.summary[,1]==responses[i],], aes(y=Mean, ymax=CI.high, ymin=CI.low,color=Model), position=position_dodge(width=0.2),size=0.3) +
      scale_fill_manual(values  = c("#5c7896","#495B46","#BAA43F","#A9796A"))+ 
      scale_color_manual(values = c("#5c7896","#495B46","#BAA43F","#A9796A"))+ 
      theme_bw() + ylim(c(0,1.0))+
      theme(axis.title.y=element_text(size=8))+
      ylab('Variability explained (R2%)')+
      theme(legend.position = 'none', axis.title.x=element_blank())
  })
  
}

## port, 11*12 plot size ##
g_rawdata_all[[1]]+all_dis1+lv1+g_r2_all[[1]]+
  g_rawdata_all[[2]]+all_dis2+lv2+g_r2_all[[2]]+
  g_rawdata_all[[3]]+all_dis3+lv3+g_r2_all[[3]]+
  g_rawdata_all[[4]]+all_dis4+lv4+g_r2_all[[4]]+
  g_rawdata_all[[5]]+all_dis5+lv5+g_r2_all[[5]]+
  plot_layout(ncol=4, widths=c(2.5,1,1,1))


  g_rawdata_all[[6]]+all_dis6+lv6+g_r2_all[[6]]+
  g_rawdata_all[[7]]+all_dis7+lv7+g_r2_all[[7]]+
 plot_layout(ncol=4, widths=c(2.5,1,1,1))
  
  g_rawdata_all[[8]]+all_dis8+lv8+g_r2_all[[8]]+
    g_rawdata_all[[9]]+all_dis9+lv9+g_r2_all[[9]]+
    g_rawdata_all[[10]]+all_dis10+lv10+g_r2_all[[10]]+
    plot_layout(ncol=4, widths=c(2.5,1,1,1))
  
  g_ms_all[[1]]+  g_ms_all[[2]]+  g_ms_all[[3]]+  g_ms_all[[4]]+
    g_ms_all[[5]]+  g_ms_all[[6]]+  g_ms_all[[7]]+  g_ms_all[[8]]+
    g_ms_all[[9]]+g_ms_all[[10]]+plot_layout(ncol=3, widths=c(1,1,1))
  
  
#### singnificant calculation ####
  df_p<-df %>%
    filter(remark!='1')
  
  df_p$remark=factor(df_p$remark, levels = c('CT','A','Bi','Be','P','St','Si','V','2','3','4','5'))
  

  root_p<- lm(WSA~remark,data=df_p)  
anova(root_p) 
qqnorm(residuals(root_p)) 
shapiro.test(residuals(root_p)) # if p>0.05, fail to reject H0
## no normal distribution ##
library(dunn.test)
#root_p<-kruskal.test(WSA~remark,data=df_p)

root_after<- dunnTest(evenness ~ remark, data = df_p, method = "bh")

# Extract the comparison results
comparison_results <- root_after$res

# Filter to only include comparisons involving 'CT'
control_comparisons <- comparison_results[grep("CT", comparison_results$Comparison), ]

# View the filtered results
print(control_comparisons)
write.csv(control_comparisons,'evenness_p.csv')


#### multifunctionality calculation ####
#z score normalization
setwd('C:/huiying data/3rd chapter/a_community')
df<- read.csv('df_rf.csv')
df_tr<-df %>%
  filter(remark!='1')


df_tr<- df_tr%>%
  mutate(WSA = ((WSA-mean(df$WSA))/sd(df[,"WSA"])),WHC = (WHC-mean(df$WHC))/sd(df[,"WHC"]),root = (root-mean(df$root))/sd(df[,"root"]),plant =(plant-mean(df$plant))/sd(df[,"plant"]),richness =(richness-mean(df$richness))/sd(df[,"richness"]),)


n_iter   = 100      # permutaion size

treatments = c("A","Be","Bi","P","St","Si","V",'2','3','4','5')
responses = c("root","plant",'WSA','WHC','richness')


###                    Functions                              ##
Bootstrap_ES_rep = function(response, data=df_tr, target=treatments, n_perm = n_iter){
  resampled = list()
  for (treatment in target) {
    bs = numeric(0)
    population_TR =data[data["remark"]==treatment,response]
    population_CK =data[data["remark"]=="CT",response]
    size_CK = length(population_CK)
    size_TR = length(population_TR)
    
    for (id in 1:n_perm) {
      k_CK = mean(sample(population_CK, size_CK, replace = T), na.rm = TRUE)
      k_TR = mean(sample(population_TR, size_TR, replace = T), na.rm = TRUE)
      bs = append(bs, k_TR - k_CK)
    }
    resampled[[treatment]] = mean(bs, na.rm = TRUE)
  }
  return(resampled)
}

##  mean effect size    ##

response_ES_bs =list()
for (i_response in responses) {
  response_ES_bs[[i_response]] = Bootstrap_ES_rep(i_response)
}
mean_sig_ES=data.frame(matrix(data = NA, nrow = 11, ncol =5))
colnames(mean_sig_ES)=responses
rownames(mean_sig_ES)=treatments

for (i_response in responses) {
  for (treatment in treatments) {
    mean_sig_ES[treatment,i_response]=response_ES_bs[[i_response]][[treatment]]
  }
}
library(dplyr)
library(caret)
#preproc <- preProcess(mean_sig_ES, method = c("range"))
#df_scaled <- predict(preproc, mean_sig_ES)
#df_scaled$Multifunctionality <- rowMeans(df_scaled, na.rm = TRUE)

mean_sig_ES$Multifunctionality <- rowMeans(mean_sig_ES, na.rm = TRUE)
## plot heatmap effect size ##
library(pheatmap)
clrsp <- colorRampPalette(c("#D54846", "white", "#619DBB"))   
clrs <- clrsp(200) 

breaks1 <- seq(-1.15,1.15, length.out = 200)
pheatmap(mean_sig_ES,border_color='black',color =  clrs, breaks = breaks1,number_color='black',gaps_row=c(7),cluster_rows=F,legend=T,cluster_cols=F,display_numbers = T,frontsize_number=10,angle_col = 45,)


#### two factor treatments figure (supplementary information) ####
setwd('C:/huiying data/3rd chapter/a_community')
df<- read.csv('df_rf.csv')
df_f2<-df %>%
  filter(remark=='2'|remark=='CT')
library(Rmisc)
df_ci<- summarySE(df_f2, measurevar="WHC", groupvars=c("treatment"))
df_ct <- df_ci[df_ci$treatment == 'CT', ]
df_others <- df_ci[df_ci$treatment != 'CT', ]

df_others_sorted <- df_others[order(df_others$WHC), ]

df_ci_sorted <- rbind(df_ct, df_others_sorted)
df_ci_sorted$treatment <- factor(df_ci_sorted$treatment, levels = df_ci_sorted$treatment)
df_f2$treatment <- factor(df_f2$treatment, levels = df_ci_sorted$treatment)


df_ci_sorted$order <-  1:nrow(df_ci_sorted)
df_f2 <- df_f2 %>%
  left_join(df_ci_sorted %>% select(treatment, order), by = "treatment")

mean_ct<- mean(df_f2$WHC[df_f2$remark == 'CT'], na.rm = TRUE)
p_whc<- ggplot(df_f2,aes(x=treatment, y=WHC))+
  geom_errorbar(data=df_ci_sorted,aes(x=as.numeric(factor(treatment)),ymin=WHC-se,ymax=WHC+se),width=0.14,size=1,color='black')+
  geom_point(data=df_f2,position=position_jitterdodge(jitter.width=0.8, dodge.width = 0.12), 
             aes(colour = order), size = 3.1,alpha=0.3,show.legend = FALSE) +
  scale_color_gradient(low = "coral3", high = "springgreen4") + 
  stat_summary(aes(group = treatment), fun = mean, geom = "point", size = 2.6, color = 'black', shape = 21, fill = 'black') +  # Add mean as dot
  scale_x_discrete(guide = guide_axis(angle = 45))+  geom_hline(yintercept=mean_ct,linetype='dashed')+
  xlab('') +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12),plot.title = element_text(size = 12))+
  ylab('Water holding capacity (%)')+
   theme_bw()

## wsa ##
df<- read.csv('df_rf.csv')
df_f2<-df %>%
  filter(remark=='2'|remark=='CT')
df_ci<- summarySE(df_f2, measurevar="WSA", groupvars=c("treatment"))
df_ct <- df_ci[df_ci$treatment == 'CT', ]
df_others <- df_ci[df_ci$treatment != 'CT', ]

df_others_sorted <- df_others[order(df_others$WSA), ]

df_ci_sorted <- rbind(df_ct, df_others_sorted)
df_ci_sorted$treatment <- factor(df_ci_sorted$treatment, levels = df_ci_sorted$treatment)
df_f2$treatment <- factor(df_f2$treatment, levels = df_ci_sorted$treatment)


df_ci_sorted$order <-  1:nrow(df_ci_sorted)
df_f2 <- df_f2 %>%
  left_join(df_ci_sorted %>% select(treatment, order), by = "treatment")

mean_ct<- mean(df_f2$WSA[df_f2$remark == 'CT'], na.rm = TRUE)
p_WSA<- ggplot(df_f2,aes(x=treatment, y=WSA))+
  geom_errorbar(data=df_ci_sorted,aes(x=as.numeric(factor(treatment)),ymin=WSA-se,ymax=WSA+se),width=0.14,size=1,color='black')+
  geom_point(data=df_f2,position=position_jitterdodge(jitter.width=0.8, dodge.width = 0.12), 
             aes(colour = order), size = 3.1,alpha=0.3,show.legend = FALSE) +
  scale_color_gradient(low = "coral3", high = "springgreen4") + 
  stat_summary(aes(group = treatment), fun = mean, geom = "point", size = 2.6, color = 'black', shape = 21, fill = 'black') +  # Add mean as dot
  scale_x_discrete(guide = guide_axis(angle = 45))+  geom_hline(yintercept=mean_ct,linetype='dashed')+
  xlab('') +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12),plot.title = element_text(size = 12))+
  ylab('Water stable aggregates (%)')+
  theme_bw()
## root ###
df<- read.csv('df_rf.csv')
df_f2<-df %>%
  filter(remark=='2'|remark=='CT')
df_ci<- summarySE(df_f2, measurevar="root", groupvars=c("treatment"))
df_ct <- df_ci[df_ci$treatment == 'CT', ]
df_others <- df_ci[df_ci$treatment != 'CT', ]

df_others_sorted <- df_others[order(df_others$root), ]

df_ci_sorted <- rbind(df_ct, df_others_sorted)
df_ci_sorted$treatment <- factor(df_ci_sorted$treatment, levels = df_ci_sorted$treatment)
df_f2$treatment <- factor(df_f2$treatment, levels = df_ci_sorted$treatment)


df_ci_sorted$order <-  1:nrow(df_ci_sorted)
df_f2 <- df_f2 %>%
  left_join(df_ci_sorted %>% select(treatment, order), by = "treatment")

mean_ct<- mean(df_f2$root[df_f2$remark == 'CT'], na.rm = TRUE)
p_root<- ggplot(df_f2,aes(x=treatment, y=root))+
  geom_errorbar(data=df_ci_sorted,aes(x=as.numeric(factor(treatment)),ymin=root-se,ymax=root+se),width=0.14,size=1,color='black')+
  geom_point(data=df_f2,position=position_jitterdodge(jitter.width=0.8, dodge.width = 0.12), 
             aes(colour = order), size = 3.1,alpha=0.3,show.legend = FALSE) +
  scale_color_gradient(low = "coral3", high = "springgreen4") + 
  stat_summary(aes(group = treatment), fun = mean, geom = "point", size = 2.6, color = 'black', shape = 21, fill = 'black') +  # Add mean as dot
  scale_x_discrete(guide = guide_axis(angle = 45))+  geom_hline(yintercept=mean_ct,linetype='dashed')+
  xlab('') +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12),plot.title = element_text(size = 12))+
  ylab('Belowground biomass (g)')+
  theme_bw()

## plant ###
df<- read.csv('df_rf.csv')
df_f2<-df %>%
  filter(remark=='2'|remark=='CT')
df_ci<- summarySE(df_f2, measurevar="plant", groupvars=c("treatment"))
df_ct <- df_ci[df_ci$treatment == 'CT', ]
df_others <- df_ci[df_ci$treatment != 'CT', ]

df_others_sorted <- df_others[order(df_others$plant), ]

df_ci_sorted <- rbind(df_ct, df_others_sorted)
df_ci_sorted$treatment <- factor(df_ci_sorted$treatment, levels = df_ci_sorted$treatment)
df_f2$treatment <- factor(df_f2$treatment, levels = df_ci_sorted$treatment)


df_ci_sorted$order <-  1:nrow(df_ci_sorted)
df_f2 <- df_f2 %>%
  left_join(df_ci_sorted %>% select(treatment, order), by = "treatment")

mean_ct<- mean(df_f2$plant[df_f2$remark == 'CT'], na.rm = TRUE)
p_plant<- ggplot(df_f2,aes(x=treatment, y=plant))+
  geom_errorbar(data=df_ci_sorted,aes(x=as.numeric(factor(treatment)),ymin=plant-se,ymax=plant+se),width=0.14,size=1,color='black')+
  geom_point(data=df_f2,position=position_jitterdodge(jitter.width=0.8, dodge.width = 0.12), 
             aes(colour = order), size = 3.1,alpha=0.3,show.legend = FALSE) +
  scale_color_gradient(low = "coral3", high = "springgreen4") + 
  stat_summary(aes(group = treatment), fun = mean, geom = "point", size = 2.6, color = 'black', shape = 21, fill = 'black') +  # Add mean as dot
  scale_x_discrete(guide = guide_axis(angle = 45))+  geom_hline(yintercept=mean_ct,linetype='dashed')+
  xlab('') +
  theme(axis.text = element_text(size = 12),axis.title = element_text(size = 12),plot.title = element_text(size = 12))+
  ylab('Aboveground biomass (g)')+
  theme_bw()
library(patchwork)
p_whc + p_WSA + p_root + p_plant + 
  plot_layout(ncol = 1) & theme(legend.position = "bottom")
