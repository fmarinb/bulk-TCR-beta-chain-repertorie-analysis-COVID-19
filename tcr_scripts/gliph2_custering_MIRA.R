###################### PACKAGE LODAING ########################

library(immunarch)
library(turboGliph)
library(tidyverse)
library(arsenal)
library(ggpubr)
library(easystats)
library(effectsize)
library(rstatix)

###################### DATA LOADING ############################
# Repertorie data
setwd("C:/Users/Usuario/OneDrive/Escritorio/TCR")
file_path<-"C:/Users/Usuario/OneDrive/Escritorio/TCR/001_imput_vdjtools"
tcr_data<-repLoad(file_path, .mode='single', .coding=T)
filtered_data<-top(.data = tcr_data$data, .n = 100)
#Turbogliph example data
data("gliph_input_data")

###################### DATA SETUP IN TRBOGLIPH FORMAT ###########
# Add a new column with sample name in each repertorie file
names_dfs<-list()
for (i in seq_along(filtered_data)) {
  df_name<-names(filtered_data)[i]
  filtered_data[[i]]<-filtered_data[[i]]%>%mutate(., patient=df_name)
}
# Concatenate by row each repertorie file
combined_df<-bind_rows(filtered_data)

# Rename and reorder variables 
combined_df2<-combined_df%>%
  select(., Clones, CDR3.aa, V.name, J.name, patient)%>%
  mutate(TRAV=NA, TRAJ=NA, CDR3a=NA)%>%
  rename(CDR3b=CDR3.aa, counts=Clones, TRBV=V.name, TRBJ=J.name)%>%
  select(all_of(names(gliph_input_data[, c(1:8)])))

combined_df2$TRAV<-as.character(combined_df2$TRAV)
combined_df2$TRAJ<-as.character(combined_df2$TRAJ)
combined_df2$CDR3a<-as.character(combined_df2$CDR3a)

imput_gliph2_covid<-as.data.frame(combined_df2)


######################## GLIPH2 REFERENCE DATA LOADING ##########
CD48<-read.delim('TCR/gliph_ref/ref_CD48_v2.0.txt', header = T)
names(CD48)<-c('CDR3b', 'TRBV', 'TRBJ' )

######################## APPLY GLIPH2 CLUSTERING ################
# gliph2 results
res4<-gliph2(cdr3_sequences=imput_gliph2_covid, 
             refdb_beta = CD48,
             kmer_mindepth = 5,
             sim_depth=1000,
             n_cores= parallel::detectCores(),
             cluster_min_size = 5,
             motif_length = c(3,4,5),
             result_folder = 'gliph2_cmin5_d5_345')
# Network building and representation
network_covid<-plot_network(clustering_output = res4,
                            n_cores=parallel::detectCores(),
                            show_additional_columns = c('patient', 
                                                        'clonal.expansion.score', 
                                                        'fisher.score'))

# Change sample names per sev_55 class name 
new_dataframe<-imput_gliph2_covid%>%
  rename(., 'Sample'='patient')
new_dataframe2<-new_dataframe%>%
  left_join(., tcr_data$meta[,c(1:12)])%>%
  select(., CDR3b, TRBV, TRBJ, CDR3a, TRAV, TRAJ, -Sample, sev_55, counts)%>%
  rename(., 'patient'='sev_55')

# Re-run GLIPH2 with previous changes
res4_sev55<-gliph2(cdr3_sequences=new_dataframe2, 
                   refdb_beta = CD48,
                   kmer_mindepth = 5,
                   sim_depth=1000,
                   n_cores= parallel::detectCores(),
                   cluster_min_size = 5,
                   motif_length = c(3,4,5),
                   result_folder = 'gliph2_cmin5_d5_345_sev55')

# Network visualization
network_covidsev55<-plot_network(clustering_output = res4_sev55,
                                 n_cores=parallel::detectCores(),
                                 cluster_min_size =5,
                                 color_info = c("patient"),
                                 show_additional_columns = c('patient', 
                                                             'clonal.expansion.score', 
                                                             'fisher.score', 
                                                             'OvE'))

######################## MIRA FILTERING ########################
# GLIPH2 cluster results loading
all_cdr3_results<-read.csv('TCR/gliph2_cmin5_d5_345_sev55/cluster_member_details.txt', sep='\t')

# MIRA data loading
mira<-read.csv('TCR/MIRA.csv', header = T , stringsAsFactors = F, sep = ";")

# MIRA dataset filtering of CDR3b sequences present in GLIPH2 clusters (MIRA positive CDR3b sequences)
common_mira3<-mira%>%
  filter(., TCR.BioIdentity %in% unique(all_cdr3_results$CDR3b))

# Table of counts of each CDR3b, reported in GLIPH2 clustering, in each sample
cdr3_sequences<-map(tcr_data$data, ~select(., c(1,4)))%>%
  map(., ~filter(., CDR3.aa %in% unique(all_cdr3_results$CDR3b)))%>%
  map(., ~group_by(., CDR3.aa))%>%
  map(., ~summarise(., Clones=sum(Clones)))

cdr3_table<-reduce(cdr3_sequences, function(x, y) full_join(x, y, by = "CDR3.aa"))%>%
  rename_at(vars(-CDR3.aa), ~ names(cdr3_sequences))%>%
  mutate_all(~ ifelse(is.na(.), 0, .))%>%
  rename(., 'CDR3b'='CDR3.aa')

################# CLUSTER TAG COUNT TABLE (MIRA POSITIVE) #####

# Table filtering to obtain GLIPH2 MIRA positive CDR3b count table
# Change the name of the CDR3b to its respective cluster tag. 
# Sum of counts to obtain the final custer tag count table.
cdr3_table3 <- cdr3_table %>%
  full_join(., all_cdr3_results[, c(1:3)], by = 'CDR3b') %>%
  filter(., CDR3b %in% common_mira3$TCR.BioIdentity) %>%
  select(., -CDR3b, -seq_ID) %>%
  select(tag, everything()) %>%
  group_by(., tag) %>%
  summarize(across(everything(), sum, na.rm = TRUE)) %>%
  column_to_rownames(., var = 'tag') %>%
  mutate_all(~ ifelse(is.na(.), 0, .))

# Transposition and metadata addition

t_cdr3_table3_2<-t(cdr3_table3)%>%
  as.data.frame()%>%
  rownames_to_column(., var='Sample')%>%
  mutate_at(vars(-Sample),
            .funs = ~ ifelse(. == 0, 0, 1))%>%
  full_join(., tcr_data$meta[, c(1:12)], by='Sample')%>%
  select(.,Sample, -Code_BD, -Platform_code, 
         -Age, -Morethan50, -Morethan60, 
         Morethan55, Sex, Severity, 
         sex_sev, sex_55, sev_55, 
         `G%YE_DEGNS`:YSSGE_4_22)

t_cdr3_table3_2<-t_cdr3_table3_2%>%
  mutate_at(vars( `G%YE_DEGNS`:YSSGE_4_22), factor)


################ CLUSTER TAG FREQUENCY TABLE (MIRA POSITIVE) ####

# Table of frequencies of each CDR3b, reported in GLIPH2 clustering, in each sample
cdr3_sequences2<-map(tcr_data$data, ~select(., c(2,4)))%>%
  map(., ~filter(., CDR3.aa %in% unique(all_cdr3_results$CDR3b)))%>%
  map(., ~group_by(., CDR3.aa))%>%
  map(., ~summarise(., Freq=sum(Proportion)))


cdr3_table_freq<-reduce(cdr3_sequences2, function(x, y) full_join(x, y, by = "CDR3.aa"))%>%
  rename_at(vars(-CDR3.aa), ~ names(cdr3_sequences))%>%
  mutate_all(~ ifelse(is.na(.), 0, .))%>%
  rename(., 'CDR3b'='CDR3.aa')

# Table filtering to obtain GLIPH2 MIRA positive CDR3b frequency table
# Change the name of the CDR3b to its respective cluster tag. 
# Sum of frequencies to obtain the final custer tag count table.

t_cdr3_table3_freq_mira<-t(cdr3_table3_freq_mira)%>%
  as.data.frame()%>%
  rownames_to_column(., var='Sample')%>%
  full_join(., tcr_data$meta[, c(1:12)], by='Sample')%>%
  select(.,Sample, -Code_BD, -Platform_code, 
         -Age, -Morethan50, -Morethan60,
         Morethan55, Sex, Severity,
         sex_sev, sex_55, sev_55, 
         `G%YE_DEGNS`:YSSGE_4_22)


########### STATISTICAL TESTS OF MIRA POSITIVE CLUSTER TAGS BETWEEN SEV_55 GROUPS ############

source('stats_functions.R')
kw_motif_mira_prop<-kw_test(df=t_cdr3_table3_3,
                            test_groups=select(t_cdr3_table3_3, sev_55)%>%colnames(), 
                            test_metrics = select(t_cdr3_table3_3, 
                                                  `G%YE_DEGNS`:YSSGE_4_22 )%>%colnames())


################## PLOTTING SIGNIFICATIVE CLUSTERS ###############

# SL%SYE_DGNST
stat.test<-dunn_test(`SP%YE_GHNRST` ~ sev_55 ,data = t_cdr3_table3_freq_mira, p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "sev_55")%>%filter(., p.adj<0.05)
stat.test2<-kruskal_test(data=t_cdr3_table3_freq_mira, formula = `SP%YE_GHNRST` ~sev_55)


GHNRST_plot<-ggplot(t_cdr3_table3_freq_mira, aes(x = sev_55 , y = log(`SP%YE_GHNRST`)))+
  geom_violin(aes(fill=sev_55), color='white')+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  theme_modern()+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ylab("log10(Frequency)") +
  ggtitle('SP%YE_GHNRST')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")+
  theme(plot.title = element_text(face="bold", size = 12))+
  theme(axis.ticks.x = element_blank())+
  theme(axis.title.y=element_text(size=20))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x=element_blank())+
  theme(strip.text.x = element_text(size = 24))+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, 
                     bracket.nudge.y = +0.01,  y.position=c(-4, -3.5))+
  annotate("text", x = 1.2, y = -2, label = paste("Kruskal-Wallis, p =", format(stat.test2$p, digits = 4)))


# SL%SYE_DGNST
stat.test<-dunn_test(`SL%SYE_DGNST` ~ sev_55 ,data =t_cdr3_table3_freq_mira, p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "sev_55")%>%filter(., p.adj<0.05)
stat.test2<-kruskal_test(data=t_cdr3_table3_freq_mira, formula =`SL%SYE_DGNST` ~sev_55)


DGNST_plot<-ggplot(t_cdr3_table3_freq_mira, aes(x = sev_55 , y = log(`SL%SYE_DGNST`)))+
  geom_violin(aes(fill=sev_55), color='white')+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  theme_modern()+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ylab("log10(Frequency)") +
  ggtitle('SL%SYE_DGNST')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")+
  theme(plot.title = element_text(face="bold", size = 12))+
  theme(axis.ticks.x = element_blank())+
  theme(axis.title.y=element_text(size=20))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x=element_blank())+
  theme(strip.text.x = element_text(size = 24))+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, 
                     bracket.nudge.y = +0.01,  y.position=c(-3))+
  annotate("text", x = 1.2, y = -2, label = paste("Kruskal-Wallis, p =", format(stat.test2$p, digits = 4)))

# SS%YE_AGST
stat.test<-dunn_test(`SS%YE_AGST` ~ sev_55 ,data = t_cdr3_table3_freq_mira, p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "sev_55")%>%filter(., p.adj<0.05)
stat.test2<-kruskal_test(data=t_cdr3_table3_freq_mira, formula =`SS%YE_AGST`~sev_55)


SSAGST_plot<-ggplot(t_cdr3_table3_freq_mira, aes(x = sev_55 , y = log(`SS%YE_AGST`)))+
  geom_violin(aes(fill=sev_55), color='white')+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  theme_modern()+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ylab("log10(Frequency)") +
  ggtitle('SS%YE_AGST')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")+
  theme(plot.title = element_text(face="bold", size = 12))+
  theme(axis.ticks.x = element_blank())+
  theme(axis.title.y=element_text(size=20))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x=element_blank())+
  theme(strip.text.x = element_text(size = 24))+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, 
                     bracket.nudge.y = +0.01,  y.position=c(-4))+
  annotate("text", x = 1.2, y = -3, label = paste("Kruskal-Wallis, p =", format(stat.test2$p, digits = 4)))

# S%GYE_AEFGHLRSVWY
stat.test<-dunn_test(`S%GYE_AEFGHLRSVWY` ~ sev_55 ,data = t_cdr3_table3_freq_mira , p.adjust.method = 'BH')
stat.test$p.adj<-round(stat.test$p.adj, 4)
stat.test <- stat.test %>% add_xy_position(x = "sev_55")%>%filter(., p.adj<0.05)
stat.test2<-kruskal_test(data=t_cdr3_table3_freq_mira, formula =`S%GYE_AEFGHLRSVWY`~sev_55)


SGESVWY_plot<-ggplot(t_cdr3_table3_freq_mira, aes(x = sev_55 , y = log(`S%GYE_AEFGHLRSVWY`)))+
  geom_violin(aes(fill=sev_55), color='white')+
  geom_boxplot(width=0.2, lwd=1, outlier.shape = NA, fill="transparent")+
  geom_jitter2(width = 0.05, alpha = 0.5)+
  theme_modern()+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ylab("log10(Frequency)") +
  ggtitle('S%GYE_AEFGHLRSVWY')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position="none")+
  theme(plot.title = element_text(face="bold", size = 12))+
  theme(axis.ticks.x = element_blank())+
  theme(axis.title.y=element_text(size=20))+
  theme(axis.text.x = element_blank())+
  theme(axis.title.x=element_blank())+
  theme(strip.text.x = element_text(size = 24))+
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01, 
                     bracket.nudge.y = +0.01,  y.position=c(-4))+
  annotate("text", x = 1.2, y = -3, label = paste("Kruskal-Wallis, p =", format(stat.test2$p, digits = 4)))

######## FREQUENCIES OF SIGNIFICATIVE TAG CDR3b SEQUENCES IN EACH SAMPLE ##########

# SP%YE_GHNRST
target_GHNRST <- c("CASSPGYEQYF",
                   "CASSPSYEQYF",
                   "CASSPTYEQYF")
tc_GHNRST <- trackClonotypes(tcr_data$data, target_GHNRST , .col = "aa")

vis(tc_GHNRST, .plot = 'area')+
  ylab('Frequency')+
  theme(axis.title.y=element_text(size=20))+
  geom_vline(xintercept= 58.5, linetype="dashed") + 
  geom_vline(xintercept= 98.5, linetype="dashed")+ 
  geom_vline(xintercept = 133.5, linetype="dashed")+
  theme(axis.text.x = element_text(size=6))

# SL%SYE_DGNST
target_DGNST<-c("CASSLDSYEQYF",
                "CASSLGSYEQYF",
                "CASSLGSYEQYF",
                "CASSLNSYEQYF",
                "CASSLSSYEQYF",
                "CASSLSSYEQYF",
                "CASSLSSYEQYF",
                "CASSLSSYEQYF",
                "CASSLSSYEQYF",
                "CASSLTSYEQYF")

tc_DGNST<-trackClonotypes(tcr_data$data, unique(target_DGNST), .col = "aa")
vis(tc_DGNST, .plot = 'area')+
  ylab('Frequency')+
  theme(axis.title.y=element_text(size=20))+
  geom_vline(xintercept= 58.5, linetype="dashed") + 
  geom_vline(xintercept= 98.5, linetype="dashed")+ 
  geom_vline(xintercept = 133.5, linetype="dashed")+
  theme(axis.text.x = element_text(size=6))

# SS%YE_AGST
target_AGST<-c("CASSLAYEQYF",
               "CASSLGYEQYF",
               "CASSLGYEQYF",
               "CASSLSYEQYF",
               "CASSLSYEQYF",
               "CASSLTYEQYF")

tc_AGST<-trackClonotypes(tcr_data$data, unique(target_AGST), .col = "aa",.norm=T)
vis(tc_AGST, .plot='area')+
  ylab('Frequency')+
  theme(axis.title.y=element_text(size=20))+
  geom_vline(xintercept= 58.5, linetype="dashed") + 
  geom_vline(xintercept= 98.5, linetype="dashed")+ 
  geom_vline(xintercept = 133.5, linetype="dashed")+
  theme(axis.text.x = element_text(size=6))

# S%GYE_AEFGHLRSVWY
target_AEFGHLRSVWY<-c("CASSAGYEQYF",
                      "CASSFGYEQYF",
                      "CASSGGYEQYF",
                      "CASSHGYEQYF",
                      "CASSLGYEQYF",
                      "CASSLGYEQYF",
                      "CASSSGYEQYF",
                      "CASSVGYEQYF",
                      "CASSWGYEQYF",
                      "CASSYGYEQYF",
                      "CATSRGYEQYF",
                      "CATSRGYEQYF")

tc_AEFGHLRSVWY<-trackClonotypes(tcr_data$data, unique(target_AEFGHLRSVWY), .col = "aa")
vis(tc_AEFGHLRSVWY, .plot = 'area')+
  ylab('Frequency')+
  theme(axis.title.y=element_text(size=20))+
  geom_vline(xintercept= 58.5, linetype="dashed") + 
  geom_vline(xintercept= 98.5, linetype="dashed")+ 
  geom_vline(xintercept = 133.5, linetype="dashed")+
  theme(axis.text.x = element_text(size=6))
