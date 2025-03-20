## list all two factor combinations ##
## supplementary information figure ##

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
