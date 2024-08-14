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
library(igraph)

###################### DATA LOADING ############################
file_path<-"C:/Users/Usuario/OneDrive/Escritorio/TCR/001_imput_vdjtools"
tcr_data<-repLoad(file_path, .mode='single', .coding=T)

##################### OVERLAP ALL SAMPLES #####################
#### ONLY AMINOACIDS #####
# JACCARD
accard_total_aa<-repOverlap(.data=tcr_data$data,
                            .method='jaccard',
                            .col='aa',
                            .verbose=T,
                            .force.matrix= T)
# MORISITA

morisita_total_aa<-repOverlap(.data=tcr_data$data,
                              .method='morisita',
                              .col='aa',
                              .verbose=T,
                              .force.matrix= T)
# PUBLIC

public_total_aa<-repOverlap(.data=tcr_data$data,
                            .method='public',
                            .col='aa',
                            .verbose=T,
                            .force.matrix= T)

save(jaccard_total_aa, file='jaccard_total_aa.RData')
save(morisita_total_aa, file='morisita_total_aa.RData')
save(public_total_aa, file='pub_total_aa.RData')

jc_tot_aa<-load('jaccard_total_aa.RData')
ms_tot_aa<-load('morisita_total_aa.RData')
pub_tot_aa<-load('pub_total_aa.RData')

#### AMINOACID AND V-J ####

# JACCARD

jaccard_total_aavj<-repOverlap(.data=tcr_data$data,
                               .method='jaccard',
                               .col='aa+v+j',
                               .verbose=T,
                               .force.matrix= T)

# MORISITA

morisita_total_aavj<-repOverlap(.data=tcr_data$data,
                                .method='morisita',
                                .col='aa+v+j',
                                .verbose=T,
                                .force.matrix= T)
# PUBLIC

public_total_aavj<-repOverlap(.data=tcr_data$data,
                              .method='public',
                              .col='aa+v+j',
                              .verbose=T,
                              .force.matrix= T)

save(jaccard_total_aavj, file='jaccard_total_aavj.RData')
save(morisita_total_aavj, file='morisita_total_aavj.RData')
save(public_total_aavj, file='pub_total_aavj.rds')

jc_tot_aavj<-load('jaccard_total_aavj.RData')
ms_tot_aavj<-load('morisita_total_aavj.RData')
pub_tot_aavj<-load('pub_total_aavj.RData')


################ FUNCTIONS TO COMPUTE OVERLAPS EASILY ########

overlaps<-function(data, type){
  if (type=='aa'){
    jaccard_aa<-repOverlap(.data=data,
                           .method='jaccard',
                           .col='aa',
                           .verbose=T,
                           .force.matrix= T)
    
    morisita_aa<-repOverlap(.data=data,
                            .method='morisita',
                            .col='aa',
                            .verbose=T,
                            .force.matrix= T)
    
    public_aa<-repOverlap(.data=data,
                          .method='public',
                          .col='aa',
                          .verbose=T,
                          .force.matrix= T)
    
    overlap_results<-list(jaccard_aa, morisita_aa, public_aa)
    save(overlap_results, file='overlap_aavj.RData')
    return(overlap_results)
    
  }else if (type=='aavj'){
    jaccard_aavj<-repOverlap(.data=data,
                             .method='jaccard',
                             .col='aa+v+j',
                             .verbose=T,
                             .force.matrix= T)
    
    morisita_aavj<-repOverlap(.data=data,
                              .method='morisita',
                              .col='aa+v+j',
                              .verbose=T,
                              .force.matrix= T)
    
    public_aavj<-repOverlap(.data=data,
                            .method='public',
                            .col='aa+v+j',
                            .verbose=T,
                            .force.matrix= T)
    
    
    overlap_results<-list(jaccard_aavj, morisita_aavj, public_aavj)
    save(overlap_results, file='overlap_aavj.RData')
    return(overlap_results)
  }else{
    print("Argument error, values are 'aa' or 'aavj'")
    
  }
}



combine_dataframes <- function(matrix_list, group_names) {
  df_list <- list()
  for (i in 1:length(matrix_list)) {
    upper_tri_values <- matrix_list[[i]][upper.tri(matrix_list[[i]])] %>% as.vector()
    group <- rep(group_names[i], length(upper_tri_values))
    df <- data.frame(Values = upper_tri_values, Groups = group)
    df_list[[i]] <- df
  }
  df <- bind_rows(df_list)
  return(df)
}




plot_inc_ov <- function(tablelist){
  # Is empty my list?
  if (length(tablelist) == 0){
    stop("La lista de tablas est? vac?a.")
  }
  # Is data.frame object checking
  if (!all(map_lgl(tablelist, is.data.frame))){
    stop("No todas las listas son tablas.")
  }
  # Join tables
  df <- bind_rows(tablelist)
  inc_plot<-ggplot(df, aes(x=depth, y=Mean))+
    geom_point(aes(color=source), size=3)+
    geom_line(aes(group = source, color=source))+
    facet_grid(. ~ source)+
    xlab("Depth") +
    theme_bw()+
    theme(legend.title = element_blank())+
    theme(axis.title.x = element_text(size=18))+
    theme(axis.title.y = element_text(size=18))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10))+
    theme(strip.text = element_blank(), strip.background = element_blank())
  
  results<-list(df=df, plot=inc_plot)
  
  
  return(results)
}


## Cummulative overlap


inc_overlap_sum<-function(inc_ov, group){
  
  # Upper triangle
  inc_ut <- lapply(inc_ov, function(m) m[upper.tri(m)])
  # Mean overlap
  inc_avg <- sapply(inc_ut, function(m) mean(m))
  # SD
  inc_sd <- sapply(inc_ut, function(m) sd(m))
  
  inc_sum<-data.frame(Mean=inc_avg, 
                      SD=inc_sd,
                      depth=names(inc_avg),
                      source=group)
  inc_sum$depth <- factor(inc_sum$depth, 
                          levels = c(50000, 10000, 5000, 1000, 500, 100))
  
  
  
  return(inc_sum)
  
}



#################### GROUP-BASED OVERLAPS ####################

################### AMINOACIDS AND V-J #######################

#### DATA RETRIEVAL ####
# Severity

severe<-repFilter(.data = tcr_data, .method = 'by.meta', 
                  list(Severity=include('Severe')))
mild<-repFilter(.data = tcr_data, .method = 'by.meta', 
                list(Severity=include('Mild')))


sev_overlaps<-overlaps(data = severe$data, type = 'aavj')
mild_overlaps<-overlaps(data = mild$data, type = 'aavj')

sev_aavj<-save(sev_overlaps,  file ='sev_aavj.RData')
mild_aavj<-save(mild_overlaps, file = 'mild_aavj.RData')


# Morethan 55

morethan55<-repFilter(.data = tcr_data, .method = 'by.meta', 
                      list(Morethan55=include('>55')))
lessthan55<-repFilter(.data = tcr_data, .method = 'by.meta', 
                      list(Morethan55=include('<55')))


m55_overlaps<-overlaps(data = morethan55$data, type = 'aavj')
l55_overlaps<-overlaps(data = lessthan55$data, type = 'aavj')

m55_aavj<-save(m55_overlaps,  file ='m55_aavj.RData')
l55_aavj<-save(l55_overlaps, file = '55_aavj.RData')

# Sex

male<-repFilter(.data = tcr_data, .method = 'by.meta', 
                list(Sex=include('Man')))
female<-repFilter(.data = tcr_data, .method = 'by.meta', 
                  list(Sex=include('Woman')))


male_overlaps<-overlaps(data = male$data, type = 'aavj')
female_overlaps<-overlaps(data = female$data, type = 'aavj')

male_aavj<-save(sev_overlaps,  file ='male_aavj.RData')
female_aavj<-save(mild_overlaps, file = 'female_aavj.RData')


# Sev_55

mild_l55<-repFilter(.data = tcr_data, .method = 'by.meta', 
                    list(sev_55=include('Mild;<55')))
mild_m55<-repFilter(.data = tcr_data, .method = 'by.meta', 
                    list(sev_55=include('Mild;>55')))
sev_l55<-repFilter(.data = tcr_data, .method = 'by.meta', 
                   list(sev_55=include('Severe;<55')))
sev_m55<-repFilter(.data = tcr_data, .method = 'by.meta', 
                   list(sev_55=include('Severe;>55')))


mild_l55_overlaps<-overlaps(data = mild_l55$data, type = 'aavj')
mild_m55_overlaps<-overlaps(data = mild_m55$data, type = 'aavj')
sev_l55_overlaps<-overlaps(data = sev_l55$data, type = 'aavj')
sev_m55_overlaps<-overlaps(data = sev_m55$data, type = 'aavj')


mild_l55_aavj<-save(mild_l55_overlaps,  file ='mild_l55_aavj.RData')
mild_m55_aavj<-save(mild_m55_overlaps, file = 'mild_m55_aavj.RData')
sev_l55_aavj<-save(sev_l55_overlaps,  file ='sev_l55_aavj.RData')
sev_m55_aavj<-save(sev_m55_overlaps, file = 'sev_m55_aavj.RData')

#### TESTING DIFERENCES AND PLOTTING ####

# Dataframes

df_jc_severity<-combine_dataframes(matrix_list = list(sev_overlaps[[1]], mild_overlaps[[1]]),
                                   group_names = c('Severe', 'Mild'))
df_mor_severity<-combine_dataframes(matrix_list = list(sev_overlaps[[2]], mild_overlaps[[2]]),
                                    group_names = c('Severe', 'Mild'))
df_pub_severity<-combine_dataframes(matrix_list = list(sev_overlaps[[3]], mild_overlaps[[3]]),
                                    group_names = c('Severe', 'Mild'))


df_jc_55<-combine_dataframes(matrix_list = list(m55_overlaps[[1]], l55_overlaps[[1]]),
                             group_names = c('>=55', '<55'))
df_mor_55<-combine_dataframes(matrix_list = list(m55_overlaps[[2]], l55_overlaps[[2]]),
                              group_names = c('>=55', '<55'))
df_pub_55<-combine_dataframes(matrix_list = list(m55_overlaps[[3]], l55_overlaps[[3]]),
                              group_names = c('>=55', '<55'))

df_jc_sex<-combine_dataframes(matrix_list = list(male_overlaps[[1]], female_overlaps[[1]]),
                              group_names = c('Male', 'Female'))
df_mor_sex<-combine_dataframes(matrix_list = list(male_overlaps[[2]], female_overlaps[[2]]),
                               group_names = c('Male', 'Female'))
df_pub_sex<-combine_dataframes(matrix_list = list(male_overlaps[[3]], female_overlaps[[3]]),
                               group_names = c('Male', 'Female'))

df_overlaps_jc<-bind_rows(df_jc_severity, df_jc_55, df_jc_sex)%>%
  mutate(Groups2 = case_when(
    Groups %in% c('Mild', 'Severe') ~ 'Severidad',
    Groups %in% c('Male', 'Female') ~ 'Sexo',
    Groups %in% c('>=55', '<55') ~ 'Edad',
    TRUE ~ NA_character_
  ))


df_overlaps_mor<-bind_rows(df_mor_severity, df_mor_55, df_mor_sex)%>%
  mutate(Groups2 = case_when(
    Groups %in% c('Mild', 'Severe') ~ 'Severidad',
    Groups %in% c('Male', 'Female') ~ 'Sexo',
    Groups %in% c('>=55', '<55') ~ 'Edad',
    TRUE ~ NA_character_
  ))


df_overlaps_pub<-bind_rows(df_pub_severity, df_pub_55, df_pub_sex)%>%
  mutate(Groups2 = case_when(
    Groups %in% c('Mild', 'Severe') ~ 'Severidad',
    Groups %in% c('Male', 'Female') ~ 'Sexo',
    Groups %in% c('>=55', '<55') ~ 'Edad',
  ))



# Normality test

table(df_overlaps_jc$Groups)
groups<-unique(df_overlaps_jc$Groups)


for (group in groups) {
  subset_data <- df_overlaps_jc%>%filter(., Groups==group)
  result <- lillie.test(subset_data$Values)
  
  cat("Group:", group, "\n")
  cat("p-value:", result$p.value, "\n")
  cat("Test statistic:", result$statistic, "\n")
  cat("-----\n")
}

for (group in groups) {
  subset_data <- df_overlaps_mor%>%filter(., Groups==group)
  result <- lillie.test(subset_data$Values)
  
  cat("Group:", group, "\n")
  cat("p-value:", result$p.value, "\n")
  cat("Test statistic:", result$statistic, "\n")
  cat("-----\n")
}

for (group in groups) {
  subset_data <- df_overlaps_pub%>%filter(., Groups==group)
  result <- lillie.test(subset_data$Values)
  
  cat("Group:", group, "\n")
  cat("p-value:", result$p.value, "\n")
  cat("Test statistic:", result$statistic, "\n")
  cat("-----\n")
}

# Test overlaps in severity groups

result_jc_sev <- wilcox_test(df_jc_severity, formula = Values~Groups, p.adjust.method = 'BH',
                             conf.level =.95, detailed = T )
effect_size_jc_sev <- wilcox_effsize(data = df_jc_severity, formula = Values~Groups )

ggplot(df_jc_severity, aes(y=log10(Values), x=Groups))+
  geom_violin(aes(fill=Groups), color='white')+
  geom_boxplot(width=0.2, lwd=1, fill="transparent")+
  scale_fill_material() +
  geom_jitter2(width = 0.05, alpha = 0.2, size=1)+
  ylab("log(Jaccard Index)") +
  theme_modern()+
  theme(legend.title = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='wilcox.test', label.y=0, label.x = 1.2, size=6 )



result_mor_sev <- wilcox_test(df_mor_severity, formula = Values~Groups, p.adjust.method = 'BH',
                              conf.level =.95, detailed = T )
effect_size_mor_sev <- wilcox_effsize(data = df_mor_severity, formula = Values~Groups )

ggplot(df_mor_severity, aes(y=log10(Values), x=Groups))+
  geom_violin(aes(fill=Groups), color='white')+
  geom_boxplot(width=0.2, lwd=1, fill="transparent")+
  scale_fill_material() +
  geom_jitter2(width = 0.05, alpha = 0.2, size=1)+
  ylab("log(Morisita's Index)") +
  theme_modern()+
  theme(legend.title = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='wilcox.test', label.y=0.3, label.x = 1.2, size=6)



result_pub_sev <- wilcox_test(df_pub_severity, formula = Values~Groups, p.adjust.method = 'BH',
                              conf.level =.95, detailed = T )
effect_size_pub_sev <- wilcox_effsize(data = df_pub_severity, formula = Values~Groups )


ggplot(df_pub_severity, aes(y=log10(Values), x=Groups))+
  geom_violin(aes(fill=Groups), color='white')+
  geom_boxplot(width=0.2, lwd=1, fill="transparent")+
  scale_fill_material() +
  geom_jitter2(width = 0.05, alpha = 0.2, size=1)+
  ylab("log(Public repertories)") +
  theme_modern()+
  theme(legend.title = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='wilcox.test', label.y=5.2, label.x = 1.2, size=6)

# Test overlaps in morethan55 groups

result_jc_55 <- wilcox_test(df_jc_55, formula = Values~Groups, p.adjust.method = 'BH',
                            conf.level =.95, detailed = T )
effect_size_jc_55 <- wilcox_effsize(data = df_jc_55, formula = Values~Groups )

ggplot(df_jc_55, aes(y=log10(Values), x=Groups))+
  geom_violin(aes(fill=Groups), color='white')+
  geom_boxplot(width=0.2, lwd=1, fill="transparent")+
  scale_fill_material() +
  geom_jitter2(width = 0.05, alpha = 0.2, size=1)+
  ylab("log(Jaccard Index)") +
  theme_modern()+
  theme(legend.title = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='wilcox.test', label.y=0, label.x = 1.2, size=6 )



result_mor_55 <- wilcox_test(df_mor_55, formula = Values~Groups, p.adjust.method = 'BH',
                             conf.level =.95, detailed = T )
effect_size_mor_55 <- wilcox_effsize(data = df_mor_55, formula = Values~Groups )

ggplot(df_mor_55, aes(y=log10(Values), x=Groups))+
  geom_violin(aes(fill=Groups), color='white')+
  geom_boxplot(width=0.2, lwd=1, fill="transparent")+
  scale_fill_material() +
  geom_jitter2(width = 0.05, alpha = 0.2, size=1)+
  ylab("log(Morisita's Index)") +
  theme_modern()+
  theme(legend.title = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='wilcox.test', label.y=0.3, label.x = 1.2, size=6)



result_pub_55 <- wilcox_test(df_pub_55, formula = Values~Groups, p.adjust.method = 'BH',
                             conf.level =.95, detailed = T )
effect_size_pub_55 <- wilcox_effsize(data = df_pub_55, formula = Values~Groups )

ggplot(df_mor_55, aes(y=log10(Values), x=Groups))+
  geom_violin(aes(fill=Groups), color='white')+
  geom_boxplot(width=0.2, lwd=1, fill="transparent")+
  scale_fill_material() +
  geom_jitter2(width = 0.05, alpha = 0.2, size=1)+
  ylab("log10(Public repertories)") +
  theme_modern()+
  theme(legend.title = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='wilcox.test', label.y=5.2, label.x = 1.2, size=6)

# Test overlaps in sex groups

result_jc_sex <- wilcox_test(df_jc_sex, formula = Values~Groups, p.adjust.method = 'BH',
                             conf.level =.95, detailed = T )
effect_size_jc_sex <- wilcox_effsize(data = df_jc_sex, formula = Values~Groups )

ggplot(df_jc_sex, aes(y=log10(Values), x=Groups))+
  geom_violin(aes(fill=Groups), color='white')+
  geom_boxplot(width=0.2, lwd=1, fill="transparent")+
  scale_fill_material() +
  geom_jitter2(width = 0.05, alpha = 0.2, size=1)+
  ylab("log10(Jaccard's overlap)") +
  theme_modern()+
  theme(legend.title = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='wilcox.test', label.y=0, label.x = 1.2, size=6 )



result_mor_sex <- wilcox_test(df_mor_sex, formula = Values~Groups, p.adjust.method = 'BH',
                              conf.level =.95, detailed = T )
effect_size_mor_sex <- wilcox_effsize(data = df_mor_sex, formula = Values~Groups )

ggplot(df_mor_sex, aes(y=log10(Values), x=Groups))+
  geom_violin(aes(fill=Groups), color='white')+
  geom_boxplot(width=0.2, lwd=1, fill="transparent")+
  scale_fill_material() +
  geom_jitter2(width = 0.05, alpha = 0.2, size=1)+
  ylab("log10(Morisita's Index)") +
  theme_modern()+
  theme(legend.title = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='wilcox.test', label.y=0.3, label.x = 1.2, size=6)



result_pub_sex <- wilcox_test(df_pub_sex, formula = Values~Groups, p.adjust.method = 'BH',
                              conf.level =.95, detailed = T )
effect_size_pub_sex <- wilcox_effsize(data = df_pub_sex, formula = Values~Groups )


ggplot(df_pub_sex, aes(y=log10(Values), x=Groups))+
  geom_violin(aes(fill=Groups), color='white')+
  geom_boxplot(width=0.2, lwd=1, fill="transparent")+
  scale_fill_material() +
  geom_jitter2(width = 0.05, alpha = 0.2, size=1)+
  ylab("log10(Public repertories)") +
  theme_modern()+
  theme(legend.title = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='wilcox.test', label.y=5.2, label.x = 1.2, size=6)

# Test overlaps in sev_55 groups (combined severity and morethan55)

# Dataframes

df_jc_sev55<-combine_dataframes(matrix_list = list(mild_l55_overlaps[[1]], 
                                                   mild_m55_overlaps[[1]],
                                                   sev_l55_overlaps[[1]],
                                                   sev_m55_overlaps[[1]]),
                                group_names = c('Mild;<55', 'Mild;>=55',
                                                'Severe;<55', 'Severe;>=55'))


df_mor_sev55<-combine_dataframes(matrix_list = list(mild_l55_overlaps[[2]], 
                                                    mild_m55_overlaps[[2]],
                                                    sev_l55_overlaps[[2]],
                                                    sev_m55_overlaps[[2]]),
                                 group_names = c('Mild;<55', 'Mild;>=55',
                                                 'Severe;<55', 'Severe;>=55'))




df_pub_sev55<-combine_dataframes(matrix_list = list(mild_l55_overlaps[[3]], 
                                                    mild_m55_overlaps[[3]],
                                                    sev_l55_overlaps[[3]],
                                                    sev_m55_overlaps[[3]]),
                                 group_names = c('Mild;<55', 'Mild;>=55',
                                                 'Severe;<55', 'Severe;>=55'))
# Normality

groups<-unique(df_jc_sev55$Groups)

for (group in groups) {
  subset_data <- df_jc_sev55%>%filter(., Groups==group)
  result <- lillie.test(subset_data$Values)
  
  cat("Group:", group, "\n")
  cat("p-value:", result$p.value, "\n")
  cat("Test statistic:", result$statistic, "\n")
  cat("-----\n")
}

for (group in groups) {
  subset_data <- df_mor_sev55%>%filter(., Groups==group)
  result <- lillie.test(subset_data$Values)
  
  cat("Group:", group, "\n")
  cat("p-value:", result$p.value, "\n")
  cat("Test statistic:", result$statistic, "\n")
  cat("-----\n")
}

for (group in groups) {
  subset_data <- df_pub_sev55%>%filter(., Groups==group)
  result <- lillie.test(subset_data$Values)
  
  cat("Group:", group, "\n")
  cat("p-value:", result$p.value, "\n")
  cat("Test statistic:", result$statistic, "\n")
  cat("-----\n")
}


# Statistical tests and plotting

  # Jaccard

kw_jc_sev55 <- kruskal_test(Values ~ Groups, data = df_jc_sev55)
ef_jc_sev55<-kruskal_effsize(Values ~ Groups, data = df_jc_sev55)

stat.test<-dunn_test(Values~ Groups,
                     data = df_jc_sev55,
                     p.adjust.method = 'BH')
stat_test_jc_sev55 <- stat.test %>% add_xy_position(x = "Groups")%>%filter(., p.adj<0.05)


ggplot(df_jc_sev55, aes(y=log10(Values), x=Groups))+
  geom_violin(aes(fill=Groups), color='white')+
  geom_boxplot(width=0.2, lwd=1, fill="transparent")+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  geom_jitter2(width = 0.05, alpha = 0.2, size=1)+
  ylab("log10(Jaccard's overlap)") +
  theme_modern()+
  theme(legend.title = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='kruskal.test', label.y=-1.1, label.x = 1.2, size=6)+
  stat_pvalue_manual(stat_test_jc_sev55, label = "p.adj.signif" , 
                     y.position=c(-1.8, -1.7, -1.6, -1.5, -1.4, -1.3 ))

  
  # Morisita

kw_mor_sev55 <- kruskal_test(Values ~ Groups, data = df_mor_sev55)
ef_mor_sev55<-kruskal_effsize(Values ~ Groups, data = df_mor_sev55)

stat.test<-dunn_test(Values~ Groups,
                     data = df_mor_sev55,
                     p.adjust.method = 'BH')
stat_test_mor_sev55 <- stat.test %>% add_xy_position(x = "Groups")%>%filter(., p.adj<0.05)

ggplot(df_mor_sev55, aes(y=log10(Values), x=Groups))+
  geom_violin(aes(fill=Groups), color='white')+
  geom_boxplot(width=0.2, lwd=1, fill="transparent")+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  geom_jitter2(width = 0.05, alpha = 0.2, size=1)+
  ylab("log10(Morisita's Overlap)") +
  theme_modern()+
  theme(legend.title = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='kruskal.test', label.y=2.5, label.x = 1.2, size=6)+
  stat_pvalue_manual(stat_test_mor_sev55, label = "p.adj.signif" , 
                     y.position=c(0, 0.4, 0.8, 1.2, 1.6))

  
  # Public

kw_pub_sev55 <- kruskal_test(Values ~ Groups, data = df_pub_sev55)
ef_pub_sev55<-kruskal_effsize(Values ~ Groups, data = df_pub_sev55)


stat.test<-dunn_test(Values~ Groups,
                     data = df_pub_sev55,
                     p.adjust.method = 'BH')
stat_test_pub_sev55 <- stat.test %>% add_xy_position(x = "Groups")%>%filter(., p.adj<0.05)

ggplot(df_pub_sev55, aes(y=log10(Values), x=Groups))+
  geom_violin(aes(fill=Groups), color='white')+
  geom_boxplot(width=0.2, lwd=1, fill="transparent")+
  scale_fill_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  geom_jitter2(width = 0.05, alpha = 0.2, size=1)+
  ylab("log10(Public repertories)") +
  theme_modern()+
  theme(legend.title = element_blank())+
  theme(axis.text.y = element_text(size=20))+
  theme(axis.title.x=element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.title.y=element_text(size=24))+
  stat_compare_means(method='kruskal.test', label.y=5.25, label.x = 1.2, size=6)+
  stat_pvalue_manual(stat_test_pub_sev55, label = "p.adj.signif" , 
                     y.position=c(4, 4.2, 4.4, 4.6, 4.8, 5))


# Summarise statistics results

statistics_results_list<-list(result_jc_sev,  result_jc_55, result_jc_sex,
                              result_mor_sev, result_mor_55, result_mor_sex,
                              result_pub_sev, result_pub_55, result_pub_sex)
statistics_overlaps_aavj<-bind_rows(statistics_results_list)

effect_size_results_list<-list(effect_size_jc_sev, effect_size_jc_55, effect_size_jc_sex,
                               effect_size_mor_sev, effect_size_mor_55, effect_size_mor_sex,
                               effect_size_pub_sev, effect_size_pub_55, effect_size_pub_sex)
effect_size_overlaps_aavj<-bind_rows(effect_size_results_list)

write.csv(statistics_overlaps_aavj, 'statistics_overlaps_aavj.csv')
write.csv(effect_size_overlaps_aavj, 'effect_size_overlaps_aavj.csv')


statistics_results_list_sev55<-list(stat_test_jc_sev55,
                                    stat_test_mor_sev55,
                                    stat_test_pub_sev55)
statistics_overlaps_aavj_sev55<-bind_rows(statistics_results_list_sev55)
statistics_overlaps_aavj_sev55$groups<-as.character(statistics_overlaps_aavj_sev55$groups)

effect_size_results_list_sev55<-list(ef_jc_sev55,
                                     ef_mor_sev55,
                                     ef_pub_sev55)
effect_size_overlaps_aavj_sev55<-bind_rows(effect_size_results_list_sev55)


statistics_results_list_kw_sev55<-list(kw_jc_sev55,
                                       kw_mor_sev55,
                                       kw_pub_sev55)
statistics_overlaps_aavj_kw_sev5<-bind_rows(statistics_results_list_kw_sev55)

write.table(statistics_overlaps_aavj_sev55, 'statistics_overlaps_aavj_sev55.csv')
write.csv(statistics_overlaps_aavj_kw_sev5, 'statistics_overlaps_aavj_kw_sev5.csv')
write.csv(effect_size_overlaps_aavj_sev55, 'effect_size_overlaps_aavj_sev55.csv')


###################### INCREMENTAL OVERLAP ###########################

# Dataframes

  # Severity

inc_sev_jc<-repOverlap(severe$data, 
                       'inc+jaccard', 
                       .col = "aa+v+j",
                       .step=c(100, 500, 1000, 5000, 10000, 50000),
                       .verbose = TRUE)
inc_mild_jc<-repOverlap(mild$data, 
                        'inc+jaccard',
                        .col = "aa+v+j",
                        .step=c(100, 500, 1000, 5000, 10000, 50000),
                        .verbose = TRUE)


inc_sev_pub<-repOverlap(severe$data, 
                        'inc+public', 
                        .col = "aa+v+j",
                        .step=c(100, 500, 1000, 5000, 10000, 50000),
                        .verbose = TRUE)
inc_mild_pub<-repOverlap(mild$data, 
                         'inc+public', 
                         .col = "aa+v+j",
                         .step=c(100, 500, 1000, 5000, 10000, 50000),
                         .verbose = TRUE)


inc_sev_morisita<-repOverlap(severe$data, 
                             'inc+morisita', 
                             .col = "aa+v+j",
                             .step=c(100, 500, 1000, 5000, 10000, 50000),
                             .verbose = TRUE)
inc_mild_morisita<-repOverlap(mild$data, 
                              'inc+morisita', 
                              .col = "aa+v+j",
                              .step=c(100, 500, 1000, 5000, 10000, 50000),
                              .verbose = TRUE)




inc_jc_sum_sev<-inc_overlap_sum(inc_sev_jc, group='Severe')
inc_jc_sum_mild<-inc_overlap_sum(inc_mild_jc, group='Mild')

inc_pub_sum_sev<-inc_overlap_sum(inc_sev_pub, group='Severe')
inc_pub_sum_mild<-inc_overlap_sum(inc_mild_pub, group='Mild')

inc_mor_sum_sev<-inc_overlap_sum(inc_sev_morisita, group='Severe')
inc_mor_sum_mild<-inc_overlap_sum(inc_mild_morisita, group='Mild')

inc_jc_sum_sev<-save(inc_jc_sum_sev,  file ='inc_jc_sum_sev.RData')
inc_jc_sum_mild<-save(inc_jc_sum_mild, file = 'inc_jc_sum_mild.RData')

inc_pub_sum_sev<-save(inc_pub_sum_sev,  file ='inc_pub_sum_sev.RData')
inc_pub_sum_mild<-save(inc_pub_sum_mild, file = 'inc_pub_sum_mild.RData')

inc_mor_sum_sev<-save(inc_mor_sum_sev,  file ='inc_mor_sum_sev.RData')
inc_mor_sum_mild<-save(inc_mor_sum_mild, file = 'inc_mor_sum_mild.RData')

  

# Morethan 55

inc_m55_jc<-repOverlap(morethan55$data, 
                       'inc+jaccard', 
                       .col = "aa+v+j",
                       .step=c(100, 500, 1000, 5000, 10000, 50000),
                       .verbose = TRUE)
inc_l55_jc<-repOverlap(lessthan55$data, 
                       'inc+jaccard', 
                       .col = "aa+v+j",
                       .step=c(100, 500, 1000, 5000, 10000, 50000),
                       .verbose = TRUE)


inc_m55_pub<-repOverlap(morethan55$data, 
                        'inc+public',
                        .col = "aa+v+j",
                        .step=c(100, 500, 1000, 5000, 10000, 50000),
                        .verbose = TRUE)
inc_l55_pub<-repOverlap(lessthan55$data, 
                        'inc+public', 
                        .col = "aa+v+j",
                        .step=c(100, 500, 1000, 5000, 10000, 50000),
                        .verbose = TRUE)


inc_m55_morisita<-repOverlap(morethan55$data, 
                             'inc+morisita',
                             .col = "aa+v+j",
                             .step=c(100, 500, 1000, 5000, 10000, 50000),
                             .verbose = TRUE)
inc_l55_morisita<-repOverlap(lessthan55$data,
                             'inc+morisita',
                             .col = "aa+v+j",
                             .step=c(100, 500, 1000, 5000, 10000, 50000),
                             .verbose = TRUE)


inc_jc_sum_m55<-inc_overlap_sum(inc_m55_jc, group='>=55')
inc_jc_sum_l55<-inc_overlap_sum(inc_l55_jc, group='<55')

inc_pub_sum_m55<-inc_overlap_sum(inc_m55_pub, group='>=55')
inc_pub_sum_l55<-inc_overlap_sum(inc_l55_pub, group='<55')

inc_mor_sum_m55<-inc_overlap_sum(inc_m55_morisita, group='>=55')
inc_mor_sum_l55<-inc_overlap_sum(inc_l55_morisita, group='<55')


inc_jc_sum_m55<-save(inc_jc_sum_m55,  file ='inc_jc_sum_m55.RData')
inc_jc_sum_l55<-save(inc_jc_sum_l55, file = 'inc_jc_sum_l55.RData')

inc_pub_sum_m55<-save(inc_pub_sum_m55,  file ='inc_pub_sum_m55.RData')
inc_pub_sum_l55<-save(inc_pub_sum_l55, file = 'inc_pub_sum_l55.RData')

inc_mor_sum_m55<-save(inc_mor_sum_m55,  file ='inc_mor_sum_m55.RData')
inc_mor_sum_l55<-save(inc_mor_sum_l55, file = 'inc_mor_sum_l55.RData')


  # Sex
inc_male_jc<-repOverlap(male$data, 
                        'inc+jaccard', 
                        .col = "aa+v+j",
                        .step=c(100, 500, 1000, 5000, 10000, 50000),
                        .verbose = TRUE)
inc_female_jc<-repOverlap(female$data, 
                          'inc+jaccard',
                          .col = "aa+v+j",
                          .step=c(100, 500, 1000, 5000, 10000, 50000),
                          .verbose = TRUE)


inc_male_pub<-repOverlap(male$data, 
                         'inc+public',
                         .col = "aa+v+j",
                         .step=c(100, 500, 1000, 5000, 10000, 50000),
                         .verbose = TRUE)
inc_female_pub<-repOverlap(female$data, 
                           'inc+public', 
                           .col = "aa+v+j",
                           .step=c(100, 500, 1000, 5000, 10000, 50000),
                           .verbose = TRUE)


inc_male_morisita<-repOverlap(male$data, 
                              'inc+morisita',
                              .col = "aa+v+j",
                              .step=c(100, 500, 1000, 5000, 10000, 50000),
                              .verbose = TRUE)
inc_femal_morisita<-repOverlap(female$data, 
                               'inc+morisita', 
                               .col = "aa+v+j",
                               .step=c(100, 500, 1000, 5000, 10000, 50000),
                               .verbose = TRUE)


inc_jc_sum_male<-inc_overlap_sum(inc_male_jc, group='Male')
inc_jc_sum_female<-inc_overlap_sum(inc_female_jc, group='Female')

inc_pub_sum_male<-inc_overlap_sum(inc_male_pub, group='Male')
inc_pub_sum_female<-inc_overlap_sum(inc_female_pub, group='Female')

inc_mor_sum_male<-inc_overlap_sum(inc_male_morisita, group='Male')
inc_mor_sum_female<-inc_overlap_sum(inc_femal_morisita, group='Female')



inc_jc_sum_male<-save(inc_jc_sum_male,  file ='inc_jc_sum_male.RData')
inc_jc_sum_female<-save(inc_jc_sum_female, file = 'inc_jc_sum_female.RData')

inc_pub_sum_male<-save(inc_pub_sum_male,  file ='inc_pub_sum_male.RData')
inc_pub_sum_female<-save(inc_pub_sum_female, file = 'inc_pub_sum_female.RData')

inc_mor_sum_male<-save(inc_mor_sum_male,  file ='inc_mor_sum_male.RData')
inc_mor_sum_female<-save(inc_mor_sum_female, file = 'inc_mor_sum_female.RData')



  #Sev_55

inc_ml55_jc<-repOverlap(mild_l55$data, 
                        'inc+jaccard', 
                        .col = "aa+v+j",
                        .step=c(100, 500, 1000, 5000, 10000, 50000),
                        .verbose = TRUE)
inc_mm55_jc<-repOverlap(mild_m55$data, 
                        'inc+jaccard',
                        .col = "aa+v+j",
                        .step=c(100, 500, 1000, 5000, 10000, 50000),
                        .verbose = TRUE)

inc_sl55_jc<-repOverlap(sev_l55$data, 
                        'inc+jaccard', 
                        .col = "aa+v+j",
                        .step=c(100, 500, 1000, 5000, 10000, 50000),
                        .verbose = TRUE)
inc_sm55_jc<-repOverlap(sev_m55$data, 
                        'inc+jaccard',
                        .col = "aa+v+j",
                        .step=c(100, 500, 1000, 5000, 10000, 50000),
                        .verbose = TRUE)



inc_ml55_pub<-repOverlap(mild_l55$data, 
                         'inc+public', 
                         .col = "aa+v+j",
                         .step=c(100, 500, 1000, 5000, 10000, 50000),
                         .verbose = TRUE)
inc_mm55_pub<-repOverlap(mild_m55$data, 
                         'inc+public', 
                         .col = "aa+v+j",
                         .step=c(100, 500, 1000, 5000, 10000, 50000),
                         .verbose = TRUE)
inc_sl55_pub<-repOverlap(sev_l55$data, 
                         'inc+public', 
                         .col = "aa+v+j",
                         .step=c(100, 500, 1000, 5000, 10000, 50000),
                         .verbose = TRUE)
inc_sm55_pub<-repOverlap(sev_m55$data, 
                         'inc+public',
                         .col = "aa+v+j",
                         .step=c(100, 500, 1000, 5000, 10000, 50000),
                         .verbose = TRUE)




inc_ml55_morisita<-repOverlap(mild_l55$data, 
                              'inc+morisita', 
                              .col = "aa+v+j",
                              .step=c(100, 500, 1000, 5000, 10000, 50000),
                              .verbose = TRUE)
inc_mm55_morisita<-repOverlap(mild_m55$data, 
                              'inc+morisita', 
                              .col = "aa+v+j",
                              .step=c(100, 500, 1000, 5000, 10000, 50000),
                              .verbose = TRUE)
inc_sl55_morisita<-repOverlap(sev_l55$data, 
                              'inc+morisita', 
                              .col = "aa+v+j",
                              .step=c(100, 500, 1000, 5000, 10000, 50000),
                              .verbose = TRUE)
inc_sm55_morisita<-repOverlap(sev_m55$data, 
                              'inc+morisita',
                              .col = "aa+v+j",
                              .step=c(100, 500, 1000, 5000, 10000, 50000),
                              .verbose = TRUE)





inc_jc_sum_ml55<-inc_overlap_sum(inc_ml55_jc, group='Mild;<55')
inc_jc_sum_mm55<-inc_overlap_sum(inc_mm55_jc, group='Mild;>=55')
inc_jc_sum_sl55<-inc_overlap_sum(inc_sl55_jc, group='Severe;<55')
inc_jc_sum_sm55<-inc_overlap_sum(inc_sm55_jc, group='Severe;>=55')

inc_pub_sum_ml55<-inc_overlap_sum(inc_ml55_pub, group='Mild;<55')
inc_pub_sum_mm55<-inc_overlap_sum(inc_mm55_pub, group='Mild;>=55')
inc_pub_sum_sl55<-inc_overlap_sum(inc_sl55_pub, group='Severe;<55')
inc_pub_sum_sm55<-inc_overlap_sum(inc_sm55_pub, group='Severe;>=55')

inc_mor_sum_ml55<-inc_overlap_sum(inc_ml55_morisita, group='Mild;<55')
inc_mor_sum_mm55<-inc_overlap_sum(inc_mm55_morisita, group='Mild;>=55')
inc_mor_sum_sl55<-inc_overlap_sum(inc_sl55_morisita, group='Severe;<55')
inc_mor_sum_sm55<-inc_overlap_sum(inc_sm55_morisita, group='Severe;>=55')






inc_jc_sev_55<-save(inc_jc_sev_55, file = 'inc_jc_sev_55.RData')
inc_pub_sev_55<-save(inc_pub_sev_55, file = 'inc_pub_sev_55.RData')
inc_mor_sev_55<-save(inc_mor_sev_55, file = 'inc_mor_sev_55.RData')


#### PLOTTING MEAN INCREMENTAL OVERLAPS PER GROUP ####

# Severity, sex and morethan55

inc_jc<-plot_inc_ov(tablelist=list(inc_jc_sum_sev, inc_jc_sum_mild, 
                                   inc_jc_sum_m55, inc_jc_sum_l55,
                                   inc_jc_sum_female, inc_jc_sum_male))

inc_pub<-plot_inc_ov(tablelist=list(inc_pub_sum_sev, inc_pub_sum_mild, 
                                    inc_pub_sum_m55, inc_pub_sum_l55,
                                    inc_pub_sum_female, inc_pub_sum_male))

inc_mor<-plot_inc_ov(tablelist=list(inc_mor_sum_sev, inc_mor_sum_mild, 
                                    inc_mor_sum_m55, inc_mor_sum_l55,
                                    inc_mor_sum_female, inc_mor_sum_male))

inc_jc<-inc_jc$plot+ggtitle('Jaccard Index')+
  theme(plot.title = element_text(size = 18, hjust = 0.5))
inc_pub<-inc_pub$plot+ggtitle('Public repertories')+
  theme(plot.title = element_text(size = 18, hjust = 0.5))
inc_mor<-inc_mor$plot+ggtitle('Morisita index')+
  theme(plot.title = element_text(size = 18, hjust = 0.5))

# sev_55
inc_jc_sev_55<-plot_inc_ov(tablelist=list(inc_jc_sum_ml55, inc_jc_sum_mm55, 
                                          inc_jc_sum_sl55, inc_jc_sum_sm55))

inc_pub_sev_55<-plot_inc_ov(tablelist=list(inc_pub_sum_ml55, inc_pub_sum_mm55, 
                                           inc_pub_sum_sl55, inc_pub_sum_sm55))

inc_mor_sev_55<-plot_inc_ov(tablelist=list(inc_mor_sum_ml55, inc_mor_sum_mm55, 
                                           inc_mor_sum_sl55, inc_mor_sum_sm55))

inc_jc_sev_55$plot<-inc_jc_sev_55$plot+
  scale_color_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ggtitle('Jaccard Index')+
  theme(plot.title = element_text(size = 18, hjust = 0.5))

inc_pub_sev_55$plot<-inc_pub_sev_55$plot+
  scale_color_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ggtitle('Public Repertories')+
  theme(plot.title = element_text(size = 18, hjust = 0.5))

inc_mor_sev_55$plot<-inc_mor_sev_55$plot+
  scale_color_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ggtitle('Morisita Index')+
  theme(plot.title = element_text(size = 18, hjust = 0.5))

library(cowplot)
p<-plot_grid(inc_jc_sev_55$plot, inc_pub_sev_55$plot, inc_mor_sev_55$plot, nrow=3)

######################### AMINOACIDS  ###############################
#### DATA RETRIEVAL ####

  # Severity

sev_overlaps_aa<-overlaps(data = severe$data, type = 'aa')
mild_overlaps_aa-overlaps(data = mild$data, type = 'aa')

sev_aa<-save(sev_overlaps_aa,  file ='sev_aa.RData')
mild_aa<-save(mild_overlaps_aa, file = 'mild_aa.RData')

  # Morethan 55
m55_overlaps_aa<-overlaps(data = severe$data, type = 'aa')
l55_overlaps_aa<-overlaps(data = mild$data, type = 'aa')

m55_aa<-save(m55_overlaps_aa,  file ='m55_aa.RData')
l55_aa<-save(l55_overlaps_aa, file = '55_aavj.RData')

  # Sex
male_overlaps_aa<-overlaps(data = severe$data, type = 'aa')
female_overlaps_aa<-overlaps(data = mild$data, type = 'aa')

male_aa<-save(male_overlaps_aa,  file ='male_aa.RData')
female_aa<-save(female_overlaps_aa, file = 'female_aa.RData')


###################### INCREMENTAL OVERLAP ###########################

# Dataframes

  # Severity

inc_sev_jc_aa<-repOverlap(severe$data, 
                          'inc+jaccard', 
                          .col = "aa",
                          .step=c(100, 500, 1000, 5000, 10000, 50000),
                          .verbose = TRUE)
inc_mild_jc_aa<-repOverlap(mild$data, 
                           'inc+jaccard',
                           .step=c(100, 500, 1000, 5000, 10000, 50000),
                           .verbose = TRUE)


inc_sev_pub_aa<-repOverlap(severe$data, 
                           'inc+public', 
                           .col = "aa+",
                           .step=c(100, 500, 1000, 5000, 10000, 50000),
                           .verbose = TRUE)
inc_mild_pub_aa<-repOverlap(mild$data, 
                            'inc+public', 
                            .col = "aa+",
                            .step=c(100, 500, 1000, 5000, 10000, 50000),
                            .verbose = TRUE)


inc_sev_morisita_aa<-repOverlap(severe$data, 
                                'inc+morisita', 
                                .col = "aa",
                                .step=c(100, 500, 1000, 5000, 10000, 50000),
                                .verbose = TRUE)
inc_mild_morisita_aa<-repOverlap(mild$data, 
                                 'inc+morisita', 
                                 .col = "aa",
                                 .step=c(100, 500, 1000, 5000, 10000, 50000),
                                 .verbose = TRUE)




inc_jc_sum_sev_aa<-inc_overlap_sum(inc_sev_jc_aa, group='Severe')
inc_jc_sum_mild_aa<-inc_overlap_sum(inc_mild_jc_aa, group='Mild')

inc_pub_sum_sev_aa<-inc_overlap_sum(inc_sev_pub_aa, group='Severe')
inc_pub_sum_mild_aa<-inc_overlap_sum(inc_mild_pub_aa, group='Mild')

inc_mor_sum_sev_aa<-inc_overlap_sum(inc_sev_morisita_aa, group='Severe')
inc_mor_sum_mild_aa<-inc_overlap_sum(inc_mild_morisita_aa, group='Mild')

inc_jc_sum_sev_aa<-save(inc_jc_sum_sev_aa,  file ='inc_jc_sum_sev_aa.RData')
inc_jc_sum_mild_aa<-save(inc_jc_sum_mild_aa, file = 'inc_jc_sum_mild_aa.RData')

inc_pub_sum_sev_aa<-save(inc_pub_sum_sev_aa,  file ='inc_pub_sum_sev_aa.RData')
inc_pub_sum_mild_aa<-save(inc_pub_sum_mild_aa, file = 'inc_pub_sum_mild_aa.RData')

inc_mor_sum_sev_aa<-save(inc_mor_sum_sev_aa,  file ='inc_mor_sum_sev_aa.RData')
inc_mor_sum_mild_aa<-save(inc_mor_sum_mild_aa, file = 'inc_mor_sum_mild_aa.RData')


  # Morethan 55

inc_m55_jc_aa<-repOverlap(morethan55$data, 
                          'inc+jaccard',
                          .col = "aa",
                          .step=c(100, 500, 1000, 5000, 10000, 50000),
                          .verbose = TRUE)
inc_l55_jc_aa<-repOverlap(lessthan55$data, 
                          'inc+jaccard',
                          .col = "aa",
                          .step=c(100, 500, 1000, 5000, 10000, 50000),
                          .verbose = TRUE)


inc_m55_pub_aa<-repOverlap(morethan55$data, 
                           'inc+public', 
                           .col = "aa",
                           .step=c(100, 500, 1000, 5000, 10000, 50000),
                           .verbose = TRUE)
inc_l55_pub_aa<-repOverlap(lessthan55$data, 
                           'inc+public', 
                           .col = "aa",
                           .step=c(100, 500, 1000, 5000, 10000, 50000),
                           .verbose = TRUE)


inc_m55_morisita_aa<-repOverlap(morethan55$data, 
                                'inc+morisita', 
                                .col = "aa",
                                .step=c(100, 500, 1000, 5000, 10000, 50000),
                                .verbose = TRUE)
inc_l55_morisita_aa<-repOverlap(lessthan55$data,
                                'inc+morisita', 
                                .col = "aa",
                                .step=c(100, 500, 1000, 5000, 10000, 50000),
                                .verbose = TRUE)


inc_jc_sum_m55_aa<-inc_overlap_sum(inc_m55_jc_aa, group='>=55')
inc_jc_sum_l55_aa<-inc_overlap_sum(inc_l55_jc_aa, group='<55')

inc_pub_sum_m55_aa<-inc_overlap_sum(inc_m55_pub_aa, group='>=55')
inc_pub_sum_l55_aa<-inc_overlap_sum(inc_l55_pub_aa, group='<55')

inc_mor_sum_m55_aa<-inc_overlap_sum(inc_m55_morisita_aa, group='>=55')
inc_mor_sum_l55_aa<-inc_overlap_sum(inc_l55_morisita_aa, group='<55')

inc_jc_sum_m55_aa<-save(inc_jc_sum_m55_aa,  file ='inc_jc_sum_m55_aa.RData')
inc_jc_sum_l55_aa<-save(inc_jc_sum_l55_aa, file = 'inc_jc_sum_l55_aa.RData')

inc_pub_sum_m55_aa<-save(inc_pub_sum_m55_aa,  file ='inc_pub_sum_m55_aa.RData')
inc_pub_sum_l55_aa<-save(inc_pub_sum_l55_aa, file = 'inc_pub_sum_l55_aa.RData')

inc_mor_sum_m55_aa<-save(inc_mor_sum_m55_aa,  file ='inc_mor_sum_m55_aa.RData')
inc_mor_sum_l55_aa<-save(inc_mor_sum_l55_aa, file = 'inc_mor_sum_l55_aa.RData')


  # Sex

inc_male_jc_aa<-repOverlap(male$data, 
                           'inc+jaccard', 
                           .col = "aa",
                           .step=c(100, 500, 1000, 5000, 10000, 50000),
                           .verbose = TRUE)
inc_female_jc_aa<-repOverlap(female$data, 
                             'inc+jaccard',
                             .col = "aa",
                             .step=c(100, 500, 1000, 5000, 10000, 50000),
                             .verbose = TRUE)


inc_male_pub_aa<-repOverlap(male$data, 
                            'inc+public', 
                            .col = "aa",
                            .step=c(100, 500, 1000, 5000, 10000, 50000),
                            .verbose = TRUE)
inc_female_pub_aa<-repOverlap(female$data, 
                              'inc+public', 
                              .col = "aa",
                              .step=c(100, 500, 1000, 5000, 10000, 50000),
                              .verbose = TRUE)


inc_male_morisita_aa<-repOverlap(male$data, 
                                 'inc+morisita', 
                                 .col = "aa",
                                 .step=c(100, 500, 1000, 5000, 10000, 50000),
                                 .verbose = TRUE)
inc_femal_morisita_aa<-repOverlap(female$data, 
                                  'inc+morisita', 
                                  .col = "aa",
                                  .step=c(100, 500, 1000, 5000, 10000, 50000),
                                  .verbose = TRUE)


inc_jc_sum_male_aa<-inc_overlap_sum(inc_male_jc_aa, group='Male')
inc_jc_sum_female_aa<-inc_overlap_sum(inc_female_jc_aa, group='Female')

inc_pub_sum_male_aa<-inc_overlap_sum(inc_male_pub_aa, group='Male')
inc_pub_sum_female_aa<-inc_overlap_sum(inc_female_pub_aa, group='Female')

inc_mor_sum_male_aa<-inc_overlap_sum(inc_male_morisita_aa, group='Male')
inc_mor_sum_female_aa<-inc_overlap_sum(inc_femal_morisita_aa, group='Female')


inc_jc_sum_male_aa<-save(inc_jc_sum_male_aa,  file ='inc_jc_sum_male_aa.RData')
inc_jc_sum_female_aa<-save(inc_jc_sum_female_aa, file = 'inc_jc_sum_female_aa.RData')

inc_pub_sum_male_aa<-save(inc_pub_sum_male_aa,  file ='inc_pub_sum_male_aa.RData')
inc_pub_sum_female_aa<-save(inc_pub_sum_female_aa, file = 'inc_pub_sum_female_aa.RData')

inc_mor_sum_male_aa<-save(inc_mor_sum_male_aa,  file ='inc_mor_sum_male_aa.RData')
inc_mor_sum_female_aa<-save(inc_mor_sum_female_aa, file = 'inc_mor_sum_female_aa.RData')





  # sev_55
inc_ml55_jc_aa<-repOverlap(mild_l55$data, 
                           'inc+jaccard', 
                           .col = "aa",
                           .step=c(100, 500, 1000, 5000, 10000, 50000),
                           .verbose = TRUE)
inc_mm55_jc_aa<-repOverlap(mild_m55$data, 
                           'inc+jaccard',
                           .col = "aa",
                           .step=c(100, 500, 1000, 5000, 10000, 50000),
                           .verbose = TRUE)

inc_sl55_jc_aa<-repOverlap(sev_l55$data, 
                           'inc+jaccard', 
                           .col = "aa",
                           .step=c(100, 500, 1000, 5000, 10000, 50000),
                           .verbose = TRUE)
inc_sm55_jc_aa<-repOverlap(sev_m55$data, 
                           'inc+jaccard',
                           .col = "aa",
                           .step=c(100, 500, 1000, 5000, 10000, 50000),
                           .verbose = TRUE)



inc_ml55_pub_aa<-repOverlap(mild_l55$data, 
                            'inc+public', 
                            .col = "aa",
                            .step=c(100, 500, 1000, 5000, 10000, 50000),
                            .verbose = TRUE)
inc_mm55_pub_aa<-repOverlap(mild_m55$data, 
                            'inc+public', 
                            .col = "aa",
                            .step=c(100, 500, 1000, 5000, 10000, 50000),
                            .verbose = TRUE)
inc_sl55_pub_aa<-repOverlap(sev_l55$data, 
                            'inc+public', 
                            .col = "aa",
                            .step=c(100, 500, 1000, 5000, 10000, 50000),
                            .verbose = TRUE)
inc_sm55_pub_aa<-repOverlap(sev_m55$data, 
                            'inc+public',
                            .col = "aa",
                            .step=c(100, 500, 1000, 5000, 10000, 50000),
                            .verbose = TRUE)




inc_ml55_morisita_aa<-repOverlap(mild_l55$data, 
                                 'inc+morisita', 
                                 .col = "aa",
                                 .step=c(100, 500, 1000, 5000, 10000, 50000),
                                 .verbose = TRUE)
inc_mm55_morisita_aa<-repOverlap(mild_m55$data, 
                                 'inc+morisita', 
                                 .col = "aa",
                                 .step=c(100, 500, 1000, 5000, 10000, 50000),
                                 .verbose = TRUE)
inc_sl55_morisita_aa<-repOverlap(sev_l55$data, 
                                 'inc+morisita', 
                                 .col = "aa",
                                 .step=c(100, 500, 1000, 5000, 10000, 50000),
                                 .verbose = TRUE)
inc_sm55_morisita_aa<-repOverlap(sev_m55$data, 
                                 'inc+morisita',
                                 .col = "aa",
                                 .step=c(100, 500, 1000, 5000, 10000, 50000),
                                 .verbose = TRUE)





inc_jc_sum_ml55_aa<-inc_overlap_sum(inc_ml55_jc_aa, group='Mild;<55')
inc_jc_sum_mm55_aa<-inc_overlap_sum(inc_mm55_jc_aa, group='Mild;>=55')
inc_jc_sum_sl55_aa<-inc_overlap_sum(inc_sl55_jc_aa, group='Severe;<55')
inc_jc_sum_sm55_aa<-inc_overlap_sum(inc_sm55_jc_aa, group='Severe;>=55')

inc_pub_sum_ml55_aa<-inc_overlap_sum(inc_ml55_pub_aa, group='Mild;<55')
inc_pub_sum_mm55_aa<-inc_overlap_sum(inc_mm55_pub_aa, group='Mild;>=55')
inc_pub_sum_sl55_aa<-inc_overlap_sum(inc_sl55_pub_aa, group='Severe;<55')
inc_pub_sum_sm55_aa<-inc_overlap_sum(inc_sm55_pub_aa, group='Severe;>=55')

inc_mor_sum_ml55_aa<-inc_overlap_sum(inc_ml55_morisita_aa, group='Mild;<55')
inc_mor_sum_mm55_aa<-inc_overlap_sum(inc_mm55_morisita_aa, group='Mild;>=55')
inc_mor_sum_sl55_aa<-inc_overlap_sum(inc_sl55_morisita_aa, group='Severe;<55')
inc_mor_sum_sm55_aa<-inc_overlap_sum(inc_sm55_morisita_aa, group='Severe;>=55')



#### PLOTTING MEAN INCREMENTAL OVERLAPS PER GROUP ####

inc_jc_aa<-plot_inc_ov(tablelist=list(inc_jc_sum_sev_aa, inc_jc_sum_mild_aa, 
                                      inc_jc_sum_m55_aa, inc_jc_sum_l55_aa,
                                      inc_jc_sum_female_aa, inc_jc_sum_male_aa))

inc_pub_aa<-plot_inc_ov(tablelist=list(inc_pub_sum_sev_aa, inc_pub_sum_mild_aa, 
                                       inc_pub_sum_m55_aa, inc_pub_sum_l55_aa,
                                       inc_pub_sum_female_aa, inc_pub_sum_male_aa))

inc_mor_aa<-plot_inc_ov(tablelist=list(inc_mor_sum_sev_aa, inc_mor_sum_mild_aa, 
                                       inc_mor_sum_m55_aa, inc_mor_sum_l55_aa,
                                       inc_mor_sum_female_aa, inc_mor_sum_male_aa))


inc_jc_aa<-inc_jc_aa+ggtitle('Jaccard Index')+
  theme(plot.title = element_text(size = 18, hjust = 0.5))
inc_pub_aa<-inc_pub_aa+ggtitle('Public repertories')+
  theme(plot.title = element_text(size = 18, hjust = 0.5))
inc_mor_aa$plot+ggtitle('Morisita Index')+
  theme(plot.title = element_text(size = 18, hjust = 0.5))



inc_jc_sev_55_aa<-plot_inc_ov(tablelist=list(inc_jc_sum_ml55_aa, inc_jc_sum_mm55_aa, 
                                             inc_jc_sum_sl55_aa, inc_jc_sum_sm55_aa))

inc_pub_sev_55_aa<-plot_inc_ov(tablelist=list(inc_pub_sum_ml55_aa, inc_pub_sum_mm55_aa, 
                                              inc_pub_sum_sl55_aa, inc_pub_sum_sm55_aa))

inc_mor_sev_55_aa<-plot_inc_ov(tablelist=list(inc_mor_sum_ml55_aa, inc_mor_sum_mm55_aa, 
                                              inc_mor_sum_sl55_aa, inc_mor_sum_sm55_aa))


inc_jc_sev_55_aa$plot<-inc_jc_sev_55_aa$plot+
  scale_color_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ggtitle('Jaccard Index')+
  theme(plot.title = element_text(size = 18, hjust = 0.5))

inc_pub_sev_55_aa$plot<-inc_pub_sev_55_aa$plot+
  scale_color_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ggtitle('Public Repertories')+
  theme(plot.title = element_text(size = 18, hjust = 0.5))

inc_mor_sev_55_aa$plot<-inc_mor_sev_55_aa$plot+
  scale_color_manual(values=c("#1C5F9E", "#6EB6D9","#F77964","#C5240E"))+
  ggtitle('Morisita Index')+
  theme(plot.title = element_text(size = 18, hjust = 0.5))


p<-plot_grid(inc_jc_sev_55_aa$plot, inc_pub_sev_55_aa$plot, inc_mor_sev_55_aa$plot, nrow=3)
