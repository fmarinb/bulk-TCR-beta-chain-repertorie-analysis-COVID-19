###################### PACKAGE LODAING ########################

library(immunarch)
library(tidyverse)
library(rstatix)
library(easystats)
require(nortest)
library(car)
library(PMCMRplus)
library(ggpubr)
library(cowplot)

###################### DATA LOADING ############################
# Repertorie data
file_path<-"C:/Users/Usuario/OneDrive/Escritorio/TCR/001_imput_vdjtools"
tcr_data<-repLoad(file_path, .mode='single', .coding=T)

################## EXPLORATORY METRICS #########################

# Number of unique clonotypes
uniq_clonotypes<-repExplore(tcr_data$data, .method = 'volume')

# Number of clones (=counts)
prof_seq<-repExplore(tcr_data$data, .method='clones')

# Clonotype abundance distribution
clonotypes_abund<-repExplore(tcr_data$data, .method='count')

##################### DIVERSITY METRICS #######################

# Chao1 abundance estimator
chao1<-repDiversity(tcr_data$data, 'chao1', .col='aa')%>%
  as.data.frame()%>%
  tibble::rownames_to_column(., "sample_id")

# Alpha-1  diversity: normalized Shannon-Wiener
# Imput the results table from VDJtools
sw_df<-read.delim("C:/Users/Usuario/OneDrive/Escritorio/TCR/TCR/diversity/001_shannon.txt",
                  header=T)

# Alpha-2 diversity: Gini-Simpson Index
gs<-repDiversity(tcr_data$data, 'gini.simp', .col='aa')%>%
  rename(sample_id=Sample)

# Div50
div_50<-repDiversity(tcr_data$data, 'd50', .col='aa')%>%
  as.data.frame()%>%
  tibble::rownames_to_column(., "sample_id")

# Gini's inequality index
gini<-repDiversity(tcr_data$data, .method='gini', .col='aa' )%>%
  as.data.frame()%>%
  tibble::rownames_to_column(., "sample_id")

# Diversity metrics per sample

chao1_2<-chao1%>%
  select(.,sample_id, Estimator)%>%
  rename(CHAO1=Estimator)

gs2<-gs%>%
  rename(GS=Value)

div_50_2<-div_50%>%
  select(.,sample_id, Clones)%>%
  rename(D50=Clones)

gini2<-gini%>%
  select(.,sample_id, V1)%>%
  rename(GINI=V1)  

diversity_df<-sw_df%>%
  select(sample_id, Morethan55, Sex, Severity, 
         sex_sev, sex_55, sev_55, 
         normalizedShannonWienerIndex_mean)%>%
  full_join(.,chao1_2, by='sample_id')%>%
  full_join(.,gs2, by='sample_id')%>%
  full_join(.,div_50_2, by='sample_id')%>%
  full_join(.,gini2, by='sample_id')

###################### STATISTICAL TESTS #########################

source('stats_functions.R')

# sex, morethan 55 and severity
test_metrics<-diversity_df%>%
  select(normalizedShannonWienerIndex_mean, CHAO1, GS, D50, GINI)%>%
  colnames()
test_groups<-diversity_df%>%
  select(Morethan55, Severity, Sex)%>%
  colnames()

norm<-check_normality(df = diversity_df, test_groups = test_groups, 
                      test_metrics= test_metrics, n = 60, g = 2)

test2<-test2(df = diversity_df, test_groups = test_groups, 
             test_metrics= test_metrics, var = T)

# Combined groups
test_metrics<-diversity_df%>%
  select(normalizedShannonWienerIndex_mean, CHAO1, GS, D50, GINI)%>%
  colnames()
test_groups<-diversity_df%>%
  select(sev_55, sex_55, sex_sev)%>%
  colnames()

norm4<-check_normality(df = diversity_df, test_groups = test_groups, 
                       test_metrics= test_metrics, n = 60, g = 4)

test_4_div<-test_more(df = diversity_df, test_groups = test_groups, 
                      test_metrics= test_metrics, var=T)

# Save the results
library(openxlsx)

statistics_results_4 <- createWorkbook()
for (i in 1:length(test_4_div)) {
  addWorksheet(statistics_results_4, sheetName = paste0("Sheet", i))
  writeData(statistics_results_4, sheet = i, x = test_4_div[[i]])
}
saveWorkbook(statistics_results_4, file = "div_stats4.xlsx")


statistics_results_2 <- createWorkbook()
for (i in 1:length(test2)) {
  addWorksheet(statistics_results_2, sheetName = paste0("Sheet", i))
  writeData(statistics_results_2, sheet = i, x = test2[[i]])
}
saveWorkbook(statistics_results_2, file = "div_stats2.xlsx")


#################### PLOTTING SIGNIFICATIVE RESULTS DIVERSITY ##################

# Morethan 55
CHAO1_55<-ggplot(diversity_df, aes(y=CHAO1, x=Morethan55))+
  geom_violin(aes(fill=Morethan55), color='white')+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  scale_fill_material() +
  ylab("Chao1 Index") +
  theme_modern()+
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='wilcox.test', label.x=1, label.y=6.5e+05, size=8)



SHANNON_55<-ggplot(diversity_df, aes(y=normalizedShannonWienerIndex_mean, 
                                     x=Morethan55))+
  geom_violin(aes(fill=Morethan55), color='white')+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  scale_fill_material() +
  ylab("Normalized Shannon-Wiener Index") +
  theme_modern()+
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='wilcox.test', label.x=1, label.y=1.01, size=8)

GS_55<-ggplot(diversity_df, aes(y=GS, x=Morethan55))+
  geom_violin(aes(fill=Morethan55), color='white')+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  scale_fill_material() +
  ylab("Gini-Simpson Index") +
  theme_modern()+
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='wilcox.test', label.x=1, label.y=1.01, size=8)

D50_55<-ggplot(diversity_df, aes(y=D50, x=Morethan55))+
  geom_violin(aes(fill=Morethan55), color='white')+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  scale_fill_material() +
  ylab("D50 Index") +
  theme_modern()+
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='wilcox.test', label.x=1, label.y=70000, size=8)

GINI_55<-ggplot(diversity_df, aes(y=GINI, x=Morethan55))+
  geom_violin(aes(fill=Morethan55), color='white')+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  scale_fill_material() +
  ylab("Gini's Inequality Index") +
  theme_modern()+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='t.test', label.x=1, label.y=1.01, size=8)

# SEv_55

stat.test<-dunn_test(CHAO1 ~ sev_55 ,
                     data = diversity_df,
                     p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "sev_55")%>%filter(., p.adj<0.05)
CHAO1_55_SEV<-ggplot(diversity_df, aes(y=CHAO1, x=sev_55))+
  geom_violin(aes(fill=sev_55), color='white')+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ylab("Chao1 Index") +
  theme_modern()+
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method = 'kruskal.test', label.y = 7.2e+05, size=8)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, size=6)



stat.test<-dunn_test(normalizedShannonWienerIndex_mean ~ sev_55 ,
                     data = diversity_df,
                     p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "sev_55")%>%filter(., p.adj<0.05)
SHANNON_55_SEV<-ggplot(diversity_df, aes(y=normalizedShannonWienerIndex_mean, x=sev_55))+
  geom_violin(aes(fill=sev_55), color='white')+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ylab("Normalized Shannon-Wiener Index Index") +
  theme_modern()+
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method = 'kruskal.test', label.y = 1.1, size=8 )+
  stat_pvalue_manual(stat.test, label = "p.adj", y.position=c(1.01, 1.035, 1.06), size=6)


stat.test<-dunn_test(GS~ sev_55 ,
                     data = diversity_df,
                     p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "sev_55")%>%filter(., p.adj<0.05)
GS_55_SEV<-ggplot(diversity_df, aes(y=GS, x=sev_55))+
  geom_violin(aes(fill=sev_55), color='white')+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ylab("Gini-Simpson Index") +
  theme_modern()+
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method = 'kruskal.test', label.y = 1.03,size=8)+
  stat_pvalue_manual(stat.test, label = "p.adj", y.position=c(1.002, 1.01, 1.018), size=6)


stat.test<-dunn_test(GS~ sev_55 ,
                     data = diversity_df,
                     p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "sev_55")%>%filter(., p.adj<0.05)
D50_55_SEV<-ggplot(diversity_df, aes(y=D50, x=sev_55))+
  geom_violin(aes(fill=sev_55), color='white')+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ylab("D50 Index") +
  theme_modern()+
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method = 'kruskal.test', label.y = 60000, size=8 )+
  stat_pvalue_manual(stat.test, label = "p.adj", y.position=c(50000, 53000, 56000),  size=6)


stat.test<-dunn_test(GINI~ sev_55 ,
                     data = diversity_df,
                     p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "sev_55")%>%filter(., p.adj<0.05)
GINI_55_SEV<-ggplot(diversity_df, aes(y=GINI, x=sev_55))+
  geom_violin(aes(fill=sev_55), color='white')+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ylab("Gini`s Inequality Index") +
  theme_modern()+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method = 'kruskal.test', label.y = 1, size=8)+
  stat_pvalue_manual(stat.test, label = "p.adj", y.position=c(0.9, 0.93, 0.96), size=6)


#################### CLONALITY FREQUENCIES #######################

############### TABLE OF CLONALITY FREQUENCIES ###############

clonality<-repClonality(.data=tcr_data$data, .method='homeo')%>%
  as.data.frame()%>%
  tibble::rownames_to_column(., "Sample")%>%
  rename(Rare='Rare (0 < X <= 1e-05)',
         Small='Small (1e-05 < X <= 1e-04)',
         Medium="Medium (1e-04 < X <= 0.001)",
         Large="Large (0.001 < X <= 0.01)",
         Hyperexpanded="Hyperexpanded (0.01 < X <= 1)")

df1<-tcr_data$meta%>%
  select(Sample, Morethan55, Sex, Severity, sev_55, sex_55, sex_sev)

clonality<-df1%>%full_join(., clonality, by='Sample')

###################### STATISTICAL TESTS #########################

source('stats_functions.R')

test_groups<-clonality%>%
  select( Morethan55, Sex, Severity)%>%
  colnames()
test_metrics<-clonality%>%
  select(Rare, Small, Medium, Large, Hyperexpanded)%>%
  colnames()

norm_clon<-check_normality(df = clonality, test_groups = test_groups, 
                           test_metrics= test_metrics, n = 60, g = 2)
test2_clon<-test2(df = clonality, test_groups = test_groups, 
                  test_metrics= test_metrics, var=T)


test_groups<-clonality%>%
  select( sev_55, sex_55, sex_sev)%>%
  colnames()
test_metrics<-clonality%>%
  select(Rare, Small, Medium, Large, Hyperexpanded)%>%
  colnames()

norm_clon4<-check_normality(df = clonality, test_groups = test_groups, 
                            test_metrics= test_metrics, n = 30, g = 4)
test4_clon<-test_more(df = clonality, test_groups = test_groups, 
                      test_metrics= test_metrics, var=T)


# Save the results
library(openxlsx)

statistics_results_clon4 <- createWorkbook()

for (i in 1:length(test4_clon)) {
  addWorksheet(statistics_results_clon4, sheetName = paste0("Sheet", i))
  writeData(statistics_results_clon4, sheet = i, x = test4_clon[[i]])
}
saveWorkbook(statistics_results_clon4, file = "clon4_stats.xlsx")


statistics_results_clon2 <- createWorkbook()

for (i in 1:length(test2_clon)) {
  addWorksheet(statistics_results_clon2, sheetName = paste0("Sheet", i))
  writeData(statistics_results_clon2, sheet = i, x = test2_clon[[i]])
}
saveWorkbook(statistics_results_clon2, file = "clon2_stats.xlsx")

############ PLOTTING SIGNIFICATIVE RESULTS CLONALITY FREQUENCIES###############

grouping_var <- "Morethan55"
quant_vars <- c("Hyperexpanded")
df_melted_hiper <- reshape2::melt(clonality, id = grouping_var, measure.vars = quant_vars)
df_melted_hiper<-df_melted_hiper%>%mutate(., sd=sd(value))

head(df_melted_hiper)

M55_HIPER<-ggbarplot(df_melted_hiper, x='Morethan55', y='value', fill='Morethan55', 
                     add = "mean_se") + 
  scale_fill_material() +
  theme_modern()+
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  stat_compare_means(method='wilcox.test', label.x=1, label.y=0.16, size=8)


grouping_var <- "Severity"
quant_vars <- c("Hyperexpanded")
df_melted_hiper <- reshape2::melt(clonality, id = grouping_var, measure.vars = quant_vars)
df_melted_hiper<-df_melted_hiper%>%mutate(., sd=sd(value))

SEV_HIPER<-ggbarplot(df_melted_hiper, x='Severity', y='value', fill='Severity', 
                     add = "mean_se") + 
  scale_fill_material() +
  theme_modern()+
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  ylab("Frequency of hiperexpanded clonotypes") +
  stat_compare_means(method='wilcox.test', label.x=1, label.y=0.16, size=8)



stat.test<-dunn_test(Hyperexpanded~ sev_55 ,
                    data = clonality,
                    p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "sev_55")%>%filter(., p.adj<0.05)

M55_SEV_HIPER<-ggbarplot(clonality, x='sev_55', y='Hyperexpanded', fill='sev_55', 
                         add = "mean_se") + 
  theme_modern()+
  ylab("Frequency") +
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ggtitle(('Hyperexpanded'))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  theme(plot.title = element_text(size = 24))+
  stat_compare_means(method = 'kruskal.test', label.y = 0.65, size=8)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01,
                     y.position=c(0.3, 0.35, 0.4), size=6)



stat.test<-dunn_test(Rare~ sev_55 ,
                     data = clonality,
                     p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "sev_55")%>%filter(., p.adj<0.05)

M55_SEV_RARE<-ggbarplot(clonality, x='sev_55', y='Rare', fill='sev_55', 
                        add = "mean_se") +
  theme_modern()+
  ylab("Frequency") +
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ggtitle('Rare')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  theme(plot.title = element_text(size = 24))+
  stat_compare_means(method = 'kruskal.test', label.y = 0.65, size=8)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01,
                     y.position=c(0.6), size=6)



M55_SEV_MEDIUM<-ggbarplot(clonality, x='sev_55', y='Medium', fill='sev_55', 
                          add = "mean_se") + 
  theme_modern()+
  ylab("Frequency") +
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ggtitle('Medium')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  theme(plot.title = element_text(size = 24))+
  stat_compare_means(method = 'anova', label.y = 0.65, size=8)




M55_SEV_SMALL<-ggbarplot(clonality, x='sev_55', y='Small', fill='sev_55', 
                         add = "mean_se") + 
  theme_modern()+
  ylab("Frequency") +
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ggtitle('Small')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  theme(plot.title = element_text(size = 24))+
  stat_compare_means(method = 'kruskal.test', label.y = 0.65, size=8)



M55_SEV_LARGE<-ggbarplot(clonality, x='sev_55', y='Large', fill='sev_55', 
                         add = "mean_se") + 
  theme_modern()+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ggtitle('Large')+
  ylab("Frequency") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  theme(plot.title = element_text(size = 24))+
  stat_compare_means(method = 'kruskal.test', label.y = 0.65, size=8)

