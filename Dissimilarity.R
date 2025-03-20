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