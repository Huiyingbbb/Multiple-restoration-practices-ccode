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
