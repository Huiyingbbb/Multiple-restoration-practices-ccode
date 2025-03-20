df<- read.csv('df_rf.csv')
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

