###################### PACKAGE LODAING ########################

library(immunarch)
library(ggstatsplot)
library(tidyverse)
library(rstatix)
library(easystats)
require(nortest)
library(car)
library(PMCMRplus)
library(ggpubr)
library(cowplot)
library(pheatmap)

###################### DATA LOADING ############################
# Repertorie data
file_path<-"C:/Users/Usuario/OneDrive/Escritorio/TCR/001_imput_vdjtools"
tcr_data<-repLoad(file_path, .mode='single', .coding=T)

##################### V ALLELES FREQUENCY TABLE ################
#  Excluding ambiguous alleles and weight them by the clonotype count
imm_gu2_w<-geneUsage(tcr_data$data, "hs.trbv", .quant = 'count', .norm=T, .ambig = "exc")

imm_gu2_w<-t(imm_gu2_w)
colnames(imm_gu2_w)<-imm_gu2_w[1,]
imm_gu2_w<-imm_gu2_w[-1,]%>%as.data.frame%>%
  tibble::rownames_to_column(., "Sample")%>%
  mutate_at(vars(-1), as.numeric)

df1<-tcr_data$meta%>%
  select(Sample, Morethan55, Sex, Severity, sev_55, sex_55, sex_sev)

v_usage_prop_w<-df1%>%full_join(., imm_gu2_w, by='Sample')%>% 
  rename_all(~str_replace_all(., "-", "_"))

###################### STATISTICAL TESTS FOR V ALLELE ###########
source('TCR/stats_functions.R')

groups<-v_usage_prop_w%>%select(Sex, Severity, Morethan55)%>%colnames()
metrics<-v_usage_prop_w%>%select(-Sample, -Sex, -Severity, -Morethan55,
                               -sex_sev, -sev_55, -sex_55)%>%colnames()


v_norm2_w<-check_normality(df = v_usage_prop_w, test_groups = groups, test_metrics = metrics, 
                           n=60, g=2)
v_homced2_w<-check_homced(df = v_usage_prop_w, 
                          test_groups = groups, 
                          test_metrics = metrics)
test2_v_w<-wilcox_groups_test(df = v_usage_prop_w, 
                 test_groups = groups, 
                 test_metrics = metrics)

# combined groups

groups<-v_usage_prop_w%>%select(sev_55, sex_55, sex_sev)%>%colnames()
metrics<-v_usage_prop_w%>%select(-Sample, -Sex, -Severity, -Morethan55,
                               -sex_sev, -sev_55, -sex_55)%>%colnames()


v_norm4_w<-check_normality(df = v_usage_prop_w, 
                           test_groups = groups, 
                           test_metrics = metrics, 
                           n=60, g=4)
v_homced4_w<-check_homced(df = v_usage_prop_w, 
                          test_groups = groups, 
                          test_metrics = metrics)
test4_v_w<-test_more(df = v_usage_prop_w, 
                     test_groups = groups, 
                     test_metrics = metrics, 
                     var = T)

#################### PLOTTING SIGNIFICATIVE V ALLELES ##############

# severity
fun<-function(x)(x*100)
v_usage_prop_perc<-v_usage_prop_w%>%
  mutate_at(.,vars(matches('TRBV')), .funs = fun )


grouping_var <- "Severity"
quant_vars <- c("TRBV15", "TRBV19", "TRBV14", "TRBV6_4")
df_melted <- reshape2::melt(v_usage_prop_w, id = grouping_var, measure.vars = quant_vars)
df_melted


v_usage_severity<-ggplot(df_melted, aes(x = Severity, y = value))+
  geom_violin(aes(fill=Severity), color='white')+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  scale_fill_material() +
  theme_modern()+
  ylab("Frequency") +
  theme(legend.title = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(axis.title.y=element_text(size=14))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x=element_blank())+
  facet_grid(. ~ variable)+
  theme(strip.text.x = element_text(size = 12))+
  stat_compare_means(method='wilcox.test', label.x=1, label.y=0.35, label = "p.format")

# Morethan 55

grouping_var <- "Morethan55"
quant_vars <- c("TRBV12_3", "TRBV4_2", "TRBV10_3", "TRBV6_4", "TRBV14")
df_melted <- reshape2::melt(v_usage_prop_w, id = grouping_var, measure.vars = quant_vars)


v_usage_55<-ggplot(df_melted, aes(x = Morethan55 , y = value))+
  geom_violin(aes(fill=Morethan55), color='white')+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  scale_fill_material() +
  theme_modern()+
  ylab("Frequency") +
  theme(legend.title = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(axis.title.y=element_text(size=14))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x=element_blank())+
  facet_grid(. ~ variable)+
  theme(strip.text.x = element_text(size = 16))+
  stat_compare_means(method='wilcox.test', label.x=1, label.y=0.15, label = "p.format")

# Sev_55

stat.test<-dunn_test(TRBV12_3 ~ sev_55 ,data = v_usage_prop_w, p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "sev_55")%>%filter(., p.adj<0.05)
stat.test

TRBV12_3_sev_55<-ggplot(v_usage_prop_w, aes(x = sev_55 , y = TRBV12_3))+
  geom_violin(aes(fill=sev_55), color='white')+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  theme_modern()+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ylab("Frequency") +
  ggtitle('TRBV12_3')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")+
  theme(plot.title = element_text(face="bold", size = 20))+
  theme(axis.ticks.x = element_blank())+
  theme(axis.title.y=element_text(size=20))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x=element_blank())+
  theme(axis.text.y = element_text(size=16))+
  stat_compare_means(method='kruskal.test', label.x=1, label.y=0.125, size=6)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, 
                     y.position=c(0.105, 0.115), size=6)



stat.test<-dunn_test(TRBV15 ~ sev_55 ,data = v_usage_prop_w, p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "sev_55")%>%filter(., p.adj<0.05)
stat.test

TRBV15_sev_55<-ggplot(v_usage_prop_w, aes(x = sev_55 , y = TRBV15))+
  geom_violin(aes(fill=sev_55), color='white')+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  theme_modern()+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ylab("Frequency") +
  ggtitle('TRBV15')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")+
  theme(plot.title = element_text(face="bold", size = 20))+
  theme(axis.ticks.x = element_blank())+
  theme(axis.title.y=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x=element_blank())+
  theme(axis.text.y = element_text(size=16))+
  stat_compare_means(method='kruskal.test', label.x=1, label.y=0.135, size=6)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, 
                     y.position=c(0.11, 0.117, 0.124), size=6)




stat.test<-dunn_test(TRBV5_5 ~ sev_55 ,data = v_usage_prop_w, p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "sev_55")%>%filter(., p.adj<0.05)
stat.test

TRBV5_5_sev_55<-ggplot(v_usage_prop_w, aes(x = sev_55 , y = TRBV5_5))+
  geom_violin(aes(fill=sev_55), color='white')+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  theme_modern()+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ylab("Frequency") +
  ggtitle('TRBV5_5')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")+
  theme(plot.title = element_text(face="bold", size = 20))+
  theme(axis.ticks.x = element_blank())+
  theme(axis.title.y=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x=element_blank())+
  theme(axis.text.y = element_text(size=16))+
  stat_compare_means(method='kruskal.test', label.x=1, label.y=0.062, size=6)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, 
                     y.position=c(0.04), size=6)




stat.test<-dunn_test(TRBV6_4 ~ sev_55 ,data = v_usage_prop_w, p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "sev_55")%>%filter(., p.adj<0.05)
stat.test


TRBV6_4_sev_55<-ggplot(v_usage_prop_w, aes(x = sev_55 , y = TRBV6_4))+
  geom_violin(aes(fill=sev_55), color='white')+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  theme_modern()+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ylab("Frequency") +
  ggtitle('TRBV6_4')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")+
  theme(plot.title = element_text(face="bold", size = 20))+
  theme(axis.ticks.x = element_blank())+
  theme(axis.title.y=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=12))+
  theme(axis.title.x=element_blank())+
  theme(strip.text.x = element_text(size = 16))+
  stat_compare_means(method='kruskal.test', label.x=1, label.y=0.028, size=6)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, 
                     bracket.nudge.y = +0.01,  y.position=c(0.012, 0.0135, 0.015), size=6)





stat.test<-dunn_test(TRBV15 ~ sev_55 ,data = v_usage_prop_w, p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "sev_55")%>%filter(., p.adj<0.05)
stat.test

TRBV14_sev_55<-ggplot(v_usage_prop_w, aes(x = sev_55 , y = TRBV14))+
  geom_violin(aes(fill=sev_55), color='white')+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  theme_modern()+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ylab("Frequency") +
  ggtitle('TRBV14')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title = element_blank())+
  theme(plot.title = element_text(face="bold", size = 20))+
  theme(legend.position = 'none')+
  theme(axis.title.y=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=16))+
  theme(axis.title.x=element_blank())+
  stat_compare_means(method='kruskal.test', label.x=1, label.y=0.127, size=6)+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, 
                     y.position=c(0.095, 0.105, 0.115), sieze=6)



##################### J ALLELES FREQUENCY TABLE ################

imm_gu2_w_j<-geneUsage(tcr_data$data, "hs.trbj", .quant = 'count', .norm = T, .ambig = "exc")

imm_gu2_w_j<-t(imm_gu2_w_j)
colnames(imm_gu2_w_j)<-imm_gu2_w_j[1,]
imm_gu2_w_j<-imm_gu2_w_j[-1,]%>%as.data.frame%>%
  tibble::rownames_to_column(., "Sample")%>%
  mutate_at(vars(-1), as.numeric)



df1<-tcr_data$meta%>%
  select(Sample, Morethan55, Sex, Severity, sev_55, sex_55, sex_sev)

j_usage_prop_w<-df1%>%full_join(., imm_gu2_w_j, by='Sample')%>% 
  rename_all(~str_replace_all(., "-", "_"))


###################### STATISTICAL TESTS FOR V ALLELE ###########

groups<-j_usage_prop_w%>%select(Sex, Severity, Morethan55)%>%colnames()
metrics<-j_usage_prop_w%>%select(-Sample, -Sex, -Severity, -Morethan55,
                                 -sex_sev, -sev_55, -sex_55)%>%colnames()


j_norm2_w<-check_normality(df = j_usage_prop_w,
                           test_groups = groups,
                           test_metrics = metrics, 
                           n=60, g=2)
j_homced2_w<-check_homced(df = j_usage_prop_w, test_groups = groups, test_metrics = metrics)
test2_j_w<-wilcox_groups_test(df = j_usage_prop_w,
                              test_groups = groups,
                              test_metrics = metrics,
                              var = T)



groups<-j_usage_prop_w%>%select(sev_55, sex_55, sex_sev)%>%colnames()
metrics<-j_usage_prop_w%>%select(-Sample, -Sex, -Severity, -Morethan55,
                                 -sex_sev, -sev_55, -sex_55)%>%colnames()

j_norm4_w<-check_normality(df = j_usage_prop_w,
                           test_groups = groups,
                           test_metrics = metrics, 
                           n=60, g=4)
j_homced4_w<-check_homced(df = j_usage_prop_w,
                          test_groups = groups,
                          test_metrics = metrics)
test4_j_w<-test_more(df = j_usage_prop_w,
                     test_groups = groups,
                     test_metrics = metrics,
                     var = T)

# No significative differences in J usage were found related to severity and age


############################# HEATMAPS ###########################

# V usage count table and clonotype weight

imm_gu2_w_count<-geneUsage(tcr_data$data, "hs.trbv", .quant = 'count', .ambig = "exc")
View(imm_gu2_w_count)

imm_gu2_w_count<-t(imm_gu2_w_count)
colnames(imm_gu2_w_count)<-imm_gu2_w_count[1,]
imm_gu2_w_count<-imm_gu2_w_count[-1,]%>%as.data.frame%>%
  tibble::rownames_to_column(., "Sample")%>%
  mutate_at(vars(-1), as.numeric)

df1<-tcr_data$meta%>%
  select(Sample, Morethan55, Sex, Severity, sev_55, sex_55, sex_sev)

v_usage_prop_w_count<-df1%>%full_join(., imm_gu2_w_count, by='Sample')%>% 
  rename_all(~str_replace_all(., "-", "_"))

    # Create the matrix and represent measures with z-score

cal_z_score <- function(x)((x - mean(x)) / sd(x))
v_usage_matrix<-v_usage_prop_w_count%>%
  select(-Sample, -Sex, -Severity, -Morethan55,
         -sex_sev, -sev_55, -sex_55)%>% 
  replace(is.na(.), 0)%>%
  as.matrix(.)
v_usage_matrix<-t(apply(v_usage_matrix, 1, cal_z_score))

    # Add annotations 

annotations<-v_usage_prop_w_count%>%
  select(Sample, Severity, Morethan55, Sex)%>%
  column_to_rownames(var = 'Sample')

row.names(v_usage_matrix)<-row.names(annotations)

my_colour = list(
  Sex = c( Man= '#54BEBA', Woman= '#833D98'),
  Morethan55 = c('<55' = '#EFBD6D', '>55' = '#682714'),
  Severity = c(Severe = "#F44336", Mild = "#2196F3"))


  # Heatmap
v_usage_count<-pheatmap(v_usage_matrix, 
                        annotation_row = annotations,  
                        annotation_colors=my_colour,
                        annotation_names_row=F,
                        show_rownames=F,
                        clustering_distance_rows="euclidean")



# J usage count table and clonotype weight

imm_gu2_w_count<-geneUsage(tcr_data$data, "hs.trbj",
                           .quant = 'count', .ambig = "exc")

imm_gu2_w_count<-t(imm_gu2_w_count)
colnames(imm_gu2_w_count)<-imm_gu2_w_count[1,]
imm_gu2_w_count<-imm_gu2_w_count[-1,]%>%as.data.frame%>%
  tibble::rownames_to_column(., "Sample")%>%
  mutate_at(vars(-1), as.numeric)

df1<-tcr_data$meta%>%
  select(Sample, Morethan55, Sex, Severity, sev_55, sex_55, sex_sev)

j_usage_prop_w_count<-df1%>%full_join(., imm_gu2_w_count, by='Sample')%>% 
  rename_all(~str_replace_all(., "-", "_"))

    # Create the matrix and represent measures with z-score

cal_z_score <- function(x)((x - mean(x)) / sd(x))
j_usage_matrix<-j_usage_prop_w_count%>%
  select(-Sample, -Sex, -Severity, -Morethan55,
         -sex_sev, -sev_55, -sex_55)%>% 
  replace(is.na(.), 0)%>%
  as.matrix(.)
j_usage_matrix<-t(apply(j_usage_matrix, 1, cal_z_score))

  # Add annotations 

annotations<-j_usage_prop_w_count%>%
  select(Sample, Severity, Morethan55, Sex)%>%
  column_to_rownames(var = 'Sample')
row.names(j_usage_matrix)<-row.names(annotations)


    # Heatmap

j_usage_count<-pheatmap(j_usage_matrix, 
                        annotation_row = annotations, 
                        annotation_names_row=F,
                        show_rownames=F,
                        annotation_colors=my_colour,
                        clustering_distance_rows="euclidean")


