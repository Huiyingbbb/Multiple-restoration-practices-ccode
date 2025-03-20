
#1. Import data

library(ggplot2)
library(ggridges)
library(patchwork)
library(party)
library(caret)
library(dplyr)
library(randomForest)
library(dabestr)

setwd('C:/huiying data/3rd chapter/a_community')
df<- read.csv('df_rf.csv')

df$treatment=factor(df$treatment, levels = unique(df$treatment))
levels = c('A-Si','A-St','A-V','A-P','A-Bi','A-Be','Be-Bi','Be-P','Be-St','Be-Si','Be-V','Bi-V','Bi-P','Bi-Si','Bi-St','Si-St','St-V','P-St','P-Si','Si-V','P-V')
stressors = c("A","Be","Bi","P","St","Si","V")
responses = c("FG_G","FG_H",'FG_L',"evenness","richness","root","plant",'WSA','WHC','pH')


df_2<- df%>%
  filter(lv!=1&lv!=3&lv!=4&lv!=5)
df_2$treatment=factor(df_2$treatment, levels = unique(df_2$treatment))
levels = c('A-Si','A-St','A-V','A-P','A-Bi','A-Be','Be-Bi','Be-P','Be-St','Be-Si','Be-V','Bi-V','Bi-P','Bi-Si','Bi-St','Si-St','St-V','P-St','P-Si','Si-V','P-V')

df_H=df_2[df_2[,"treatment"]%in%levels,]
df_CK=df_2[df_2[,"treatment"]=="CT",]


# 2. Loading functions

# Functions for calculating null model predictions for specific factor combinations
# Tutorial: https://mohanb96.github.io/Nullmodels.html

NullModel = function(object, selected_factors=vector(), n_perm = 50){
  
  output = list()
  CT = object[["summary"]][["mean"]][[1]]
  input=object[["data"]]
  population_CT= input[input[,as.character(object[["x"]][[2]])] == object[["result"]][["control_group"]][1], object[["result"]][["variable"]][1]]
  size_CT=length(population_CT)
  
  for (type in c("additive","multiplicative", "dominative")) {
    bs = numeric(0)
    for (id in c(1:n_perm)) {
      each_effect = numeric(0)
      k_CT = mean(sample(population_CT, size_CT, replace = T))
      for (treatment in selected_factors) {
        population_TR = input[input[,as.character(object[["x"]][[2]])] == object[["result"]][["test_group"]][which(object[["result"]][["test_group"]]==treatment)], object[["result"]][["variable"]][1]]
        size_TR = length(population_TR)
        k_TR = mean(sample(population_TR,size_TR,replace = T))
        
        # ES estimate depending on the type of null hypothesis
        if(type == "additive")        each_effect = append(each_effect,(k_TR - k_CT))
        if(type == "multiplicative")  each_effect = append(each_effect, (k_TR-k_CT)/k_CT)
        if(type == "dominative")      each_effect = append(each_effect,(k_TR - k_CT))
      }
      
      if(type == "additive") {
        joint_effect = sum(each_effect)
        pre_response = joint_effect+CT
      }
      
      if(type=="multiplicative"){
        z = 1
        for(m in c(1:length(selected_factors))) {
          z = z * (1 + each_effect[m])
          joint_effect = (z - 1)*k_CT
        }
        pre_response =joint_effect+CT
      }
      
      if(type=="dominative")  {
        joint_effect = each_effect[which(max(abs(each_effect))==abs(each_effect))]
        pre_response =joint_effect+CT
      }
      
      bs = append(bs, pre_response)
    }
    output[[type]] = bs
  }
  return(output)
}



NullModel_summary = function(Null_modle_i_Treatment, actual_data){
  output = list()

    output[["actual"]] = c(quantile(actual_value[[1]], .025), mean(actual_value[[1]]),quantile(actual_value[[1]], .975),1)
    for (type in c("additive","multiplicative", "dominative")) { 
      bs = numeric(0)
      actu_bs = numeric(0)
      for (i in 1:50) {
        k_CK = mean(sample(df_CK[,i_response], length(df_CK[,i_response]), replace = T), na.rm = TRUE)
        k_TR = mean(sample(Null_modle_i_Treatment[[type]], length(Null_modle_i_Treatment[[type]]), replace = T), na.rm = TRUE)
        bs = append(bs, k_TR - k_CK)
        k_actu = mean(sample(actual_data[,1], length(actual_data[,1]), replace = T), na.rm = TRUE)
        actu_bs = append(actu_bs, k_actu - k_CK)
      }
      compare<- data.frame(actu_bs,bs)
      colnames(compare)<- c('actual','predict')
      an<- aov(actual~predict,data=compare)
      p = summary(an)[[1]][["Pr(>F)"]][1]
      output[[type]] = c(quantile(Null_modle_i_Treatment[[type]], .025), mean(Null_modle_i_Treatment[[type]]), quantile(Null_modle_i_Treatment[[type]], .975), p)
  }
  output = t(data.frame(output))
  colnames(output) = c("X2.5%", "mean", "X97.5%", "p_value")
  output = data.frame(ES = c("actual", "additive","multiplicative","dominative"), output)
  row.names(output) = c()
  
  return(output)
}






# 3. Effect size calculation



ES.plot_FG_H<-df_2%>%
  dabest(treatment, FG_H,
         idx = c("CT",'A-Si','A-St','A-V','A-P','A-Bi','A-Be','Be-Bi','Be-P','Be-St','Be-Si','Be-V','Bi-V','Bi-P','Bi-Si','Bi-St','Si-St','St-V','P-St','P-Si','Si-V','P-V'),
         paired = FALSE)
ES.plot.meandiff_FG_H<-mean_diff(ES.plot_FG_H,reps = 500)

ES.plot_FG_G<-df_2%>%
  dabest(treatment, FG_G,
         idx = c("CT",'A-Si','A-St','A-V','A-P','A-Bi','A-Be','Be-Bi','Be-P','Be-St','Be-Si','Be-V','Bi-V','Bi-P','Bi-Si','Bi-St','Si-St','St-V','P-St','P-Si','Si-V','P-V'),
         paired = FALSE)
ES.plot.meandiff_FG_G<-mean_diff(ES.plot_FG_G,reps = 500)

ES.plot_FG_L<-df_2%>%
  dabest(treatment, FG_L,
         idx = c("CT",'A-Si','A-St','A-V','A-P','A-Bi','A-Be','Be-Bi','Be-P','Be-St','Be-Si','Be-V','Bi-V','Bi-P','Bi-Si','Bi-St','Si-St','St-V','P-St','P-Si','Si-V','P-V'),
         paired = FALSE)
ES.plot.meandiff_FG_L<-mean_diff(ES.plot_FG_L,reps = 500)

ES.plot_richness<-df_2%>%
  dabest(treatment, richness,
         idx = c("CT",'A-Si','A-St','A-V','A-P','A-Bi','A-Be','Be-Bi','Be-P','Be-St','Be-Si','Be-V','Bi-V','Bi-P','Bi-Si','Bi-St','Si-St','St-V','P-St','P-Si','Si-V','P-V'),
         paired = FALSE)
ES.plot.meandiff_richness<-mean_diff(ES.plot_richness,reps = 500)

ES.plot_evenness<-df_2%>%
  dabest(treatment, evenness,
         idx = c("CT",'A-Si','A-St','A-V','A-P','A-Bi','A-Be','Be-Bi','Be-P','Be-St','Be-Si','Be-V','Bi-V','Bi-P','Bi-Si','Bi-St','Si-St','St-V','P-St','P-Si','Si-V','P-V'),
         paired = FALSE)
ES.plot.meandiff_evenness<-mean_diff(ES.plot_evenness,reps = 500)


ES.plot_plant<-df_2%>%
  dabest(treatment, plant,
         idx = c("CT",'A-Si','A-St','A-V','A-P','A-Bi','A-Be','Be-Bi','Be-P','Be-St','Be-Si','Be-V','Bi-V','Bi-P','Bi-Si','Bi-St','Si-St','St-V','P-St','P-Si','Si-V','P-V'),
         paired = FALSE)
ES.plot.meandiff_plant<-mean_diff(ES.plot_plant,reps = 500)

ES.plot_root<-df_2%>%
  dabest(treatment, root,
         idx = c("CT",'A-Si','A-St','A-V','A-P','A-Bi','A-Be','Be-Bi','Be-P','Be-St','Be-Si','Be-V','Bi-V','Bi-P','Bi-Si','Bi-St','Si-St','St-V','P-St','P-Si','Si-V','P-V'),
         paired = FALSE)
ES.plot.meandiff_root<-mean_diff(ES.plot_root,reps = 500)

ES.plot_WSA<-df_2%>%
  dabest(treatment, WSA,
         idx = c("CT",'A-Si','A-St','A-V','A-P','A-Bi','A-Be','Be-Bi','Be-P','Be-St','Be-Si','Be-V','Bi-V','Bi-P','Bi-Si','Bi-St','Si-St','St-V','P-St','P-Si','Si-V','P-V'),
         paired = FALSE)
ES.plot.meandiff_WSA<-mean_diff(ES.plot_WSA,reps = 500)

ES.plot_pH<-df_2%>%
  dabest(treatment, pH,
         idx = c("CT",'A-Si','A-St','A-V','A-P','A-Bi','A-Be','Be-Bi','Be-P','Be-St','Be-Si','Be-V','Bi-V','Bi-P','Bi-Si','Bi-St','Si-St','St-V','P-St','P-Si','Si-V','P-V'),
         paired = FALSE)
ES.plot.meandiff_pH<-mean_diff(ES.plot_pH,reps = 500)

ES.plot_WHC<-df_2%>%
  dabest(treatment, WHC,
         idx = c("CT",'A-Si','A-St','A-V','A-P','A-Bi','A-Be','Be-Bi','Be-P','Be-St','Be-Si','Be-V','Bi-V','Bi-P','Bi-Si','Bi-St','Si-St','St-V','P-St','P-Si','Si-V','P-V'),
         paired = FALSE)
ES.plot.meandiff_WHC<-mean_diff(ES.plot_WHC,reps = 500)

plot(ES.plot.meandiff_FG_G)
#p

# 4. Main analysis (Null model)

ES.plot.meandiff = list()

ES.plot.meandiff[["FG_G"]]=ES.plot.meandiff_FG_G
ES.plot.meandiff[["FG_H"]]=ES.plot.meandiff_FG_H
ES.plot.meandiff[["FG_L"]]=ES.plot.meandiff_FG_L
ES.plot.meandiff[["evenness"]]=ES.plot.meandiff_FG_H
ES.plot.meandiff[["richness"]]=ES.plot.meandiff_FG_L

ES.plot.meandiff[["plant"]]=ES.plot.meandiff_plant
ES.plot.meandiff[["root"]]=ES.plot.meandiff_root
ES.plot.meandiff[["pH"]]=ES.plot.meandiff_pH
ES.plot.meandiff[["WSA"]]=ES.plot.meandiff_WSA
ES.plot.meandiff[["WHC"]]=ES.plot.meandiff_WHC

df_H=df_2[df_2[,"treatment"]%in%levels,]
df_CK=df_2[df_2[,"treatment"]=="CT",]
n_perm = 50

Diviation_additive = data.frame(matrix(data = NA, nrow = 21, ncol = 10))
colnames(Diviation_additive)=c("Lv","diviation","P","dissim","actual_ES","null_ES",'low_actual','high_actual','low_model','high_model')

Diviation_multiplicative = data.frame(matrix(data = NA, nrow = 21, ncol = 10))
colnames(Diviation_multiplicative)=c("Lv","diviation","P","dissim","actual_ES","null_ES",'low_actual','high_actual','low_model','high_model')

Diviation_dominative = data.frame(matrix(data = NA, nrow = 21, ncol = 10))
colnames(Diviation_dominative)=c("Lv","diviation","P","dissim","actual_ES","null_ES",'low_actual','high_actual','low_model','high_model')

# add function


Diviation_3_models = list()
for (i_response in responses) {
  i=1
  for (combination_i in levels) {
    print(paste("Combination:", combination_i))
      actual_value <- df_H %>%
      filter(treatment ==combination_i) %>%
      select(i_response)

    print(paste("Actual data:", actual_value))
    print(paste("i=", i))
    print(paste("processed measurement", i_response))
    # If it fails here, you will see which row causes the issue
    Null_modle_i_Treatment = NullModel(ES.plot.meandiff[[i_response]], selected_factors = combination_i, n_perm = n_perm)
    
    Null_modle_summary_i_Treatment = NullModel_summary(Null_modle_i_Treatment, actual_data = actual_value)
    
    Diviation_additive[i,2]=(Null_modle_summary_i_Treatment[1,3]-Null_modle_summary_i_Treatment[2,3])/Null_modle_summary_i_Treatment[2,3]  #changed
    Diviation_additive[i,3]=Null_modle_summary_i_Treatment[2,5]
    Diviation_additive[i,4]= df_H %>% filter(treatment == combination_i) %>% pull(22) %>% unique() %>% .[1]
    Diviation_additive[i,1]=combination_i
    Diviation_additive[i,5]=Null_modle_summary_i_Treatment[1,3]-mean(df_CK[,i_response])
    Diviation_additive[i,6]=Null_modle_summary_i_Treatment[2,3]-mean(df_CK[,i_response])
    Diviation_additive[i,7]=Null_modle_summary_i_Treatment[1,2]-mean(df_CK[,i_response])
    Diviation_additive[i,8]=Null_modle_summary_i_Treatment[1,4]-mean(df_CK[,i_response])
    Diviation_additive[i,9]=Null_modle_summary_i_Treatment[2,2]-mean(df_CK[,i_response])
    Diviation_additive[i,10]=Null_modle_summary_i_Treatment[2,4]-mean(df_CK[,i_response])
    
    Diviation_multiplicative[i,2]=(Null_modle_summary_i_Treatment[1,3]-Null_modle_summary_i_Treatment[3,3])/Null_modle_summary_i_Treatment[3,3]
    Diviation_multiplicative[i,3]=Null_modle_summary_i_Treatment[3,5]
    Diviation_multiplicative[i,4]= df_H %>% filter(treatment == combination_i) %>% pull(22) %>% unique() %>% .[1]
    Diviation_multiplicative[i,1]=combination_i
    Diviation_multiplicative[i,5]=Null_modle_summary_i_Treatment[1,3]-mean(df_CK[,i_response])
    Diviation_multiplicative[i,6]=Null_modle_summary_i_Treatment[3,3]-mean(df_CK[,i_response])
    Diviation_multiplicative[i,7]=Null_modle_summary_i_Treatment[1,2]-mean(df_CK[,i_response])
    Diviation_multiplicative[i,8]=Null_modle_summary_i_Treatment[1,4]-mean(df_CK[,i_response])
    Diviation_multiplicative[i,9]=Null_modle_summary_i_Treatment[3,2]-mean(df_CK[,i_response])
    Diviation_multiplicative[i,10]=Null_modle_summary_i_Treatment[3,4]-mean(df_CK[,i_response])
    
    Diviation_dominative[i,2]=(Null_modle_summary_i_Treatment[1,3]-Null_modle_summary_i_Treatment[4,3])/Null_modle_summary_i_Treatment[4,3]
    Diviation_dominative[i,3]=Null_modle_summary_i_Treatment[4,5]
    Diviation_dominative[i,4]= df_H %>% filter(treatment == combination_i) %>% pull(22) %>% unique() %>% .[1]
    Diviation_dominative[i,1]=combination_i
    Diviation_dominative[i,5]=Null_modle_summary_i_Treatment[1,3]-mean(df_CK[,i_response])
    Diviation_dominative[i,6]=Null_modle_summary_i_Treatment[4,3]-mean(df_CK[,i_response])
    Diviation_dominative[i,7]=Null_modle_summary_i_Treatment[1,2]-mean(df_CK[,i_response])
    Diviation_dominative[i,8]=Null_modle_summary_i_Treatment[1,4]-mean(df_CK[,i_response])
    Diviation_dominative[i,9]=Null_modle_summary_i_Treatment[4,2]-mean(df_CK[,i_response])
    Diviation_dominative[i,10]=Null_modle_summary_i_Treatment[4,4]-mean(df_CK[,i_response])
    
    i=i+1
    }
  
  Diviation_3_models[[i_response]][["additive"]]=Diviation_additive
  Diviation_3_models[[i_response]][["multiplicative"]]=Diviation_multiplicative
  Diviation_3_models[[i_response]][["dominative"]]=Diviation_dominative
  
}



# 5.1 Identifying interaction type according to p values

null_models = c("additive", "multiplicative", "dominative")

for (i_response in responses) {
  for (n_model in null_models) {
    for (i in 1:21) {
      # Check for Synergistic interaction
      if (Diviation_3_models[[i_response]][[n_model]][["P"]][i]<0.05) {
        if (Diviation_3_models[[i_response]][[n_model]][["actual_ES"]][i]>Diviation_3_models[[i_response]][[n_model]][["null_ES"]][i]){
        Diviation_3_models[[i_response]][[n_model]][["interaction_type"]][i] <- "synergistic"
        } else {Diviation_3_models[[i_response]][[n_model]][["interaction_type"]][i] <- "antagonistic"}
      }  else {
        Diviation_3_models[[i_response]][[n_model]][["interaction_type"]][i] <- "Null"
      }
    }
  }
}

# Ploting (point)

p_diviation = list()

for (i_response in responses) {
  
  p_diviation[[i_response]][["additive"]] = ggplot(Diviation_3_models[[i_response]][["additive"]], aes(x=dissim, y=actual_ES, color = interaction_type))+
    geom_point(size =2.2)+
    scale_color_manual(values = c("synergistic"='#5c7896',"Null"="grey70",'antagonistic'='#A9796A'))+    theme_bw()+
    geom_hline(yintercept=0,linetype='dashed')+
    #geom_errorbar(data=Diviation_3_models[[i_response]][["additive"]],aes(x=dissim,ymin=low_actual,ymax=high_actual),width=0.02,size=0.5,color='black')+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  ylim(-5,5)
  
  p_diviation[[i_response]][["multiplicative"]] = ggplot(Diviation_3_models[[i_response]][["multiplicative"]], aes(x=dissim, y=actual_ES, color = interaction_type))+
    geom_point(size =2.2)+
    scale_color_manual(values = c("synergistic"='#5c7896',"Null"="grey70",'antagonistic'='#A9796A'))+    theme_bw()+
    geom_hline(yintercept=0,linetype='dashed')+
     theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  ylim(-5,5)
  
  p_diviation[[i_response]][["dominative"]] = ggplot(Diviation_3_models[[i_response]][["dominative"]], aes(x=dissim, y=actual_ES, color = interaction_type))+
    geom_point(size =2.2)+ 
    scale_color_manual(values = c("synergistic"='#5c7896',"Null"="grey70",'antagonistic'='#A9796A'))+
    geom_hline(yintercept=0,linetype='dashed')+
    theme_bw()+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  
}

p_diviation[["WHC"]][["additive"]]+p_diviation[["WHC"]][["multiplicative"]]+p_diviation[["WHC"]][["dominative"]]+
  p_diviation[["pH"]][["additive"]]+p_diviation[["pH"]][["multiplicative"]]+p_diviation[["pH"]][["dominative"]]+
  p_diviation[["WSA"]][["additive"]]+p_diviation[["WSA"]][["multiplicative"]]+p_diviation[["WSA"]][["dominative"]]+
  plot_layout(ncol=3, widths=c(2,2,2))

p_diviation[["FG_G"]][["additive"]]+p_diviation[["FG_G"]][["multiplicative"]]+p_diviation[["FG_G"]][["dominative"]]+
  p_diviation[["FG_H"]][["additive"]]+p_diviation[["FG_H"]][["multiplicative"]]+p_diviation[["FG_H"]][["dominative"]]+
  p_diviation[["FG_L"]][["additive"]]+p_diviation[["FG_L"]][["multiplicative"]]+p_diviation[["FG_L"]][["dominative"]]+
  p_diviation[["richness"]][["additive"]]+p_diviation[["richness"]][["multiplicative"]]+p_diviation[["richness"]][["dominative"]]+
  p_diviation[["evenness"]][["additive"]]+p_diviation[["evenness"]][["multiplicative"]]+p_diviation[["evenness"]][["dominative"]]+
  plot_layout(ncol=3, widths=c(2,2,2))

p_diviation[["plant"]][["additive"]]+p_diviation[["plant"]][["multiplicative"]]+p_diviation[["plant"]][["dominative"]]+
  p_diviation[["root"]][["additive"]]+p_diviation[["root"]][["multiplicative"]]+p_diviation[["root"]][["dominative"]]+
  plot_layout(ncol=3, widths=c(2,2,2))




#### till here ####
# 5.2 Summarize proportions of occurence of factor net interaction types by three null model assumptions

# 5.2.1 Summarize proportion of occurence of factor net interaction types by factor level


### 3 dissimilarity index ###

null_models = c("additive", "multiplicative", "dominative")


count_additive = data.frame(matrix(data = NA, nrow = 9, ncol = 3))
colnames(count_additive)=c("interaction_type","dissimilarity_range",'count')
count_additive[,1]<- rep(c('synergistic','antagonistic','no'),3)
count_additive[,2]<- rep(c('low','mid','high'),each=3)

count_dominative = data.frame(matrix(data = NA, nrow = 9, ncol = 3))
colnames(count_dominative)=c("interaction_type","dissimilarity_range",'count')
count_dominative[,1]<- rep(c('synergistic','antagonistic','no'),3)
count_dominative[,2]<- rep(c('low','mid','high'),each=3)

count_multiplicative= data.frame(matrix(data = NA, nrow = 9, ncol = 3))
colnames(count_multiplicative)=c("interaction_type","dissimilarity_range",'count')
count_multiplicative[,1]<- rep(c('synergistic','antagonistic','no'),3)
count_multiplicative[,2]<- rep(c('low','mid','high'),each=3)

Overall_count = list()
count<- list()

for (i_response in responses) {
 
    a_i_s=0
    a_i_a=0
    a_i_n=0
    b_i_s=0
    b_i_a=0
    b_i_n=0
 
    d_i_s=0
    d_i_a=0
    d_i_n=0
    for(i in 1:21){
    if(Diviation_3_models[[i_response]][['additive']][["dissim"]][i]<0.29){
      if(Diviation_3_models[[i_response]][['additive']][["interaction_type"]][i]=='synergistic'){
        a_i_s<- a_i_s+1
      }else if(Diviation_3_models[[i_response]][['additive']][["interaction_type"]][i]=='antagonistic'){
        a_i_a<- a_i_a+1
      } else(a_i_n<- a_i_n+1)
    }
      else if (Diviation_3_models[[i_response]][['additive']][["dissim"]][i]>0.29 &Diviation_3_models[[i_response]][['additive']][["dissim"]][i]<0.53 ){
        if(Diviation_3_models[[i_response]][['additive']][["interaction_type"]][i]=='synergistic'){
          b_i_s<- b_i_s+1
        }else if(Diviation_3_models[[i_response]][['additive']][["interaction_type"]][i]=='antagonistic'){
          b_i_a<- b_i_a+1
        } else(b_i_n<- b_i_n+1)
      }
      else{
        if(Diviation_3_models[[i_response]][['additive']][["interaction_type"]][i]=='synergistic'){
          d_i_s<- d_i_s+1
        }else if(Diviation_3_models[[i_response]][['additive']][["interaction_type"]][i]=='antagonistic'){
          d_i_a<- d_i_a+1
        } else(d_i_n<- d_i_n+1)
      }
    }
    count_additive[,3]<- c(a_i_s,a_i_a,a_i_n,b_i_s,b_i_a,b_i_n,d_i_s,d_i_a,d_i_n)
  
    a_i_s=0
    a_i_a=0
    a_i_n=0
    b_i_s=0
    b_i_a=0
    b_i_n=0
  
    d_i_s=0
    d_i_a=0
    d_i_n=0
    
    for(i in 1:21){
      if(Diviation_3_models[[i_response]][['dominative']][["dissim"]][i]<0.29){
        if(Diviation_3_models[[i_response]][['dominative']][["interaction_type"]][i]=='synergistic'){
          a_i_s<- a_i_s+1
        }else if(Diviation_3_models[[i_response]][['dominative']][["interaction_type"]][i]=='antagonistic'){
          a_i_a<- a_i_a+1
        } else(a_i_n<- a_i_n+1)
      }
      else if (Diviation_3_models[[i_response]][['dominative']][["dissim"]][i]>0.29 &Diviation_3_models[[i_response]][['dominative']][["dissim"]][i]<0.53 ){
        if(Diviation_3_models[[i_response]][['dominative']][["interaction_type"]][i]=='synergistic'){
          b_i_s<- b_i_s+1
        }else if(Diviation_3_models[[i_response]][['dominative']][["interaction_type"]][i]=='antagonistic'){
          b_i_a<- b_i_a+1
        } else(b_i_n<- b_i_n+1)
      }
      else{
        if(Diviation_3_models[[i_response]][['dominative']][["interaction_type"]][i]=='synergistic'){
          d_i_s<- d_i_s+1
        }else if(Diviation_3_models[[i_response]][['dominative']][["interaction_type"]][i]=='antagonistic'){
          d_i_a<- d_i_a+1
        } else(d_i_n<- d_i_n+1)
      }
    }
    count_dominative[,3]<- c(a_i_s,a_i_a,a_i_n,b_i_s,b_i_a,b_i_n,d_i_s,d_i_a,d_i_n)
    
    a_i_s=0
    a_i_a=0
    a_i_n=0
    b_i_s=0
    b_i_a=0
    b_i_n=0

    d_i_s=0
    d_i_a=0
    d_i_n=0
    
    for(i in 1:21){
      if(Diviation_3_models[[i_response]][['multiplicative']][["dissim"]][i]<0.29){
        if(Diviation_3_models[[i_response]][['multiplicative']][["interaction_type"]][i]=='synergistic'){
          a_i_s<- a_i_s+1
        }else if(Diviation_3_models[[i_response]][['multiplicative']][["interaction_type"]][i]=='antagonistic'){
          a_i_a<- a_i_a+1
        } else(a_i_n<- a_i_n+1)
      }
      else if (Diviation_3_models[[i_response]][['multiplicative']][["dissim"]][i]>0.29 &Diviation_3_models[[i_response]][['multiplicative']][["dissim"]][i]<0.53 ){
        if(Diviation_3_models[[i_response]][['multiplicative']][["interaction_type"]][i]=='synergistic'){
          b_i_s<- b_i_s+1
        }else if(Diviation_3_models[[i_response]][['multiplicative']][["interaction_type"]][i]=='antagonistic'){
          b_i_a<- b_i_a+1
        } else(b_i_n<- b_i_n+1)
      }
     else{
        if(Diviation_3_models[[i_response]][['multiplicative']][["interaction_type"]][i]=='synergistic'){
          d_i_s<- d_i_s+1
        }else if(Diviation_3_models[[i_response]][['multiplicative']][["interaction_type"]][i]=='antagonistic'){
          d_i_a<- d_i_a+1
        } else(d_i_n<- d_i_n+1)
      }
    }
    count_multiplicative[,3]<- c(a_i_s,a_i_a,a_i_n,b_i_s,b_i_a,b_i_n,d_i_s,d_i_a,d_i_n)
    
  Overall_count[[i_response]][["additive"]]=count_additive
  Overall_count[[i_response]][["multiplicative"]]=count_multiplicative
  Overall_count[[i_response]][["dominative"]]=count_dominative
  
  
  }


### 2 dissimilarity range ###
null_models = c("additive", "multiplicative", "dominative")


count_additive = data.frame(matrix(data = NA, nrow = 6, ncol = 3))
colnames(count_additive)=c("interaction_type","dissimilarity_range",'count')
count_additive[,1]<- rep(c('synergistic','antagonistic','no'),2)
count_additive[,2]<- rep(c('low','high'),each=3)

count_dominative = data.frame(matrix(data = NA, nrow = 6, ncol = 3))
colnames(count_dominative)=c("interaction_type","dissimilarity_range",'count')
count_dominative[,1]<- rep(c('synergistic','antagonistic','no'),2)
count_dominative[,2]<- rep(c('low','high'),each=3)

count_multiplicative= data.frame(matrix(data = NA, nrow = 6, ncol = 3))
colnames(count_multiplicative)=c("interaction_type","dissimilarity_range",'count')
count_multiplicative[,1]<- rep(c('synergistic','antagonistic','no'),2)
count_multiplicative[,2]<- rep(c('low','high'),each=3)

Overall_count = list()
count<- list()

for (i_response in responses) {
  
  a_i_s=0
  a_i_a=0
  a_i_n=0

  
  d_i_s=0
  d_i_a=0
  d_i_n=0
  for(i in 1:21){
    if(Diviation_3_models[[i_response]][['additive']][["dissim"]][i]<0.34){
      if(Diviation_3_models[[i_response]][['additive']][["interaction_type"]][i]=='synergistic'){
        a_i_s<- a_i_s+1
      }else if(Diviation_3_models[[i_response]][['additive']][["interaction_type"]][i]=='antagonistic'){
        a_i_a<- a_i_a+1
      } else(a_i_n<- a_i_n+1)
    }

    else{
      if(Diviation_3_models[[i_response]][['additive']][["interaction_type"]][i]=='synergistic'){
        d_i_s<- d_i_s+1
      }else if(Diviation_3_models[[i_response]][['additive']][["interaction_type"]][i]=='antagonistic'){
        d_i_a<- d_i_a+1
      } else(d_i_n<- d_i_n+1)
    }
  }
  count_additive[,3]<- c(a_i_s,a_i_a,a_i_n,d_i_s,d_i_a,d_i_n)
  
  a_i_s=0
  a_i_a=0
  a_i_n=0
  
  d_i_s=0
  d_i_a=0
  d_i_n=0
  
  for(i in 1:21){
    if(Diviation_3_models[[i_response]][['dominative']][["dissim"]][i]<0.34){
      if(Diviation_3_models[[i_response]][['dominative']][["interaction_type"]][i]=='synergistic'){
        a_i_s<- a_i_s+1
      }else if(Diviation_3_models[[i_response]][['dominative']][["interaction_type"]][i]=='antagonistic'){
        a_i_a<- a_i_a+1
      } else(a_i_n<- a_i_n+1)
    }
    else{
      if(Diviation_3_models[[i_response]][['dominative']][["interaction_type"]][i]=='synergistic'){
        d_i_s<- d_i_s+1
      }else if(Diviation_3_models[[i_response]][['dominative']][["interaction_type"]][i]=='antagonistic'){
        d_i_a<- d_i_a+1
      } else(d_i_n<- d_i_n+1)
    }
  }
  count_dominative[,3]<- c(a_i_s,a_i_a,a_i_n,d_i_s,d_i_a,d_i_n)
  
  a_i_s=0
  a_i_a=0
  a_i_n=0

  d_i_s=0
  d_i_a=0
  d_i_n=0
  
  for(i in 1:21){
    if(Diviation_3_models[[i_response]][['multiplicative']][["dissim"]][i]<0.34){
      if(Diviation_3_models[[i_response]][['multiplicative']][["interaction_type"]][i]=='synergistic'){
        a_i_s<- a_i_s+1
      }else if(Diviation_3_models[[i_response]][['multiplicative']][["interaction_type"]][i]=='antagonistic'){
        a_i_a<- a_i_a+1
      } else(a_i_n<- a_i_n+1)
    }
    else{
      if(Diviation_3_models[[i_response]][['multiplicative']][["interaction_type"]][i]=='synergistic'){
        d_i_s<- d_i_s+1
      }else if(Diviation_3_models[[i_response]][['multiplicative']][["interaction_type"]][i]=='antagonistic'){
        d_i_a<- d_i_a+1
      } else(d_i_n<- d_i_n+1)
    }
  }
  count_multiplicative[,3]<- c(a_i_s,a_i_a,a_i_n,d_i_s,d_i_a,d_i_n)
  
  Overall_count[[i_response]][["additive"]]=count_additive
  Overall_count[[i_response]][["multiplicative"]]=count_multiplicative
  Overall_count[[i_response]][["dominative"]]=count_dominative
  
  
}

# Ploting percentage barplot

p_percent_level = list()

for (i_response in responses) {
  for (n_model in null_models) {
    Overall_count[[i_response]][[n_model]][['dissimilarity_range']]<- factor(Overall_count[[i_response]][[n_model]][['dissimilarity_range']],levels=c('low','high'))
}}

for (i_response in responses) {
  p_percent_level[[i_response]][["additive"]]=ggplot(Overall_count[[i_response]][["additive"]], aes(fill = interaction_type, y=count, x=dissimilarity_range))+
    geom_bar(position = "fill", stat = "identity",alpha=0.52)+
    scale_fill_manual(values = c("synergistic"='#5c7896',"no"="grey70",'antagonistic'='#A9796A'))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  p_percent_level[[i_response]][["multiplicative"]]=ggplot(Overall_count[[i_response]][["multiplicative"]], aes(fill = interaction_type, y=count, x=dissimilarity_range))+
    geom_bar(position = "fill", stat = "identity",alpha=0.52)+
    scale_fill_manual(values = c("synergistic"='#5c7896',"no"="grey70",'antagonistic'='#A9796A'))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  p_percent_level[[i_response]][["dominative"]]=ggplot(Overall_count[[i_response]][["dominative"]], aes(fill = interaction_type, y=count, x=dissimilarity_range))+
    geom_bar(position = "fill", stat = "identity",alpha=0.52)+
    scale_fill_manual(values = c("synergistic"='#5c7896',"no"="grey70",'antagonistic'='#A9796A'))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    theme(legend.position = 'none',axis.title.x = element_blank())+
    xlab(i_response)
  
}

  p_percent_level[["WSA"]][["additive"]]+p_percent_level[["WSA"]][["multiplicative"]]+p_percent_level[["WSA"]][["dominative"]]+
    p_percent_level[["WHC"]][["additive"]]+p_percent_level[["WHC"]][["multiplicative"]]+p_percent_level[["WHC"]][["dominative"]]+
     p_percent_level[["pH"]][["additive"]]+p_percent_level[["pH"]][["multiplicative"]]+p_percent_level[["pH"]][["dominative"]]+
     plot_layout(ncol=3, widths=c(2,2,2))

p_percent_level[["FG_G"]][["additive"]]+p_percent_level[["FG_G"]][["multiplicative"]]+p_percent_level[["FG_G"]][["dominative"]]+
  p_percent_level[["FG_H"]][["additive"]]+p_percent_level[["FG_H"]][["multiplicative"]]+p_percent_level[["FG_H"]][["dominative"]]+
  p_percent_level[["FG_L"]][["additive"]]+p_percent_level[["FG_L"]][["multiplicative"]]+p_percent_level[["FG_L"]][["dominative"]]+
  p_percent_level[['evenness']][["additive"]]+p_percent_level[["evenness"]][["multiplicative"]]+p_percent_level[["evenness"]][["dominative"]]+
  p_percent_level[['richness']][["additive"]]+p_percent_level[["richness"]][["multiplicative"]]+p_percent_level[["richness"]][["dominative"]]+
  plot_layout(ncol=3, widths=c(2,2,2))

p_percent_level[["plant"]][["additive"]]+p_percent_level[["plant"]][["multiplicative"]]+p_percent_level[["plant"]][["dominative"]]+
  p_percent_level[["root"]][["additive"]]+p_percent_level[['root']][["multiplicative"]]+p_percent_level[["root"]][["dominative"]]+
  plot_layout(ncol=3, widths=c(2,2,2))

top_plots <- p_percent_level[["WSA"]][["additive"]] + 
  p_percent_level[["WHC"]][["additive"]] + 
  p_percent_level[["pH"]][["additive"]] + 
  p_percent_level[["FG_G"]][["additive"]] + 
  p_percent_level[["FG_H"]][["additive"]] + 
  p_percent_level[["FG_L"]][["additive"]] + 
  plot_layout(ncol=3, widths=c(2,2,2))

# Define the second row (evenness & richness), aligned to the right
middle_plots <- 
  plot_spacer() + 
  p_percent_level[["plant"]][["additive"]] + 
  p_percent_level[["root"]][["additive"]] + 
  plot_spacer() + 
  p_percent_level[['evenness']][["additive"]] + 
  p_percent_level[['richness']][["additive"]] +  
  plot_layout(ncol=3, widths=c(2,2,2))


# Combine all plots
top_plots / middle_plots 
