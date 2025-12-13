library(readxl)
library(tidyverse)
library(ggpubr)
library(rstatix)  

metadata_biomarker<-read_excel("Final_Altered_metadata_biomarker_only.xlsx")

for(x in names(metadata_biomarker)){
  c<-trimws(metadata_biomarker[[x]])
  c_num<-as.numeric(str_remove(metadata_biomarker[[x]],"<"))
  metadata_biomarker[[x]]<-c_num
  test_result<- wilcox.test(metadata_biomarker[[x]]~metadata_biomarker$`Opioid=Substance`)
  if (test_result$p.value<0.05){
    print(x)
  }
}

colnames(metadata_biomarker)[colnames(metadata_biomarker) == "Opioid=Substance"] <- "Opioid"

bxp6 <- ggboxplot(
  metadata_biomarker, x = "Opioid", y = "AUC insulin", 
  ylab = "AUC insulin", xlab = "Opioid=Substance",add = "jitter")+theme(
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 18)  
  )
bxp6

stat.test <- metadata_biomarker %>% 
  wilcox_test(`AUC insulin` ~ Opioid) %>%
  add_significance() %>% 
  add_xy_position(x = "Opioid")
stat.test

plot6 <- bxp6 + 
  stat_pvalue_manual(stat.test, tip.length = 0) 

plot6

ggsave("plot6.png", width = 6, height = 4, dpi = 300)

#Significant Biomarkers:
#[1] "FPG"
#[1] "Fasting insulin"
#[1] "Fasting Proinsulin"
#[1] "OGIS"
#[1] "AUC glucose"
#[1] "AUC insulin"
