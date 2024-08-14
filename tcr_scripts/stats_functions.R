###### Functions to check the normality and homocedasticity####################
library(nortest)
library(easystats)
library(effectsize)
library(rstatix)
check_normality<-function(df, test_groups, test_metrics, n, g){
  
  df_total<-df
  
  groups<-df_total%>%
    select(., all_of(test_groups))%>%
    colnames()
  
  div_names<-df_total%>%
    select(., all_of(test_metrics))%>%
    colnames()
  
  group<-c()
  div<-c()
  pval<-c()
  tests<-c()
  
  # Conditional to choose the apropiate test
  for (j in 1:length(groups)){ 
    for (i in 1:length(div_names)){
      
      # Create the table of values
      test_table<-df_total%>%
        select(groups[j], div_names[i])%>%
        as.data.frame()%>%
        na.omit() 
      
      # Select what test is convenient depending on the subsample's size
      if (n<50){
        print('N<50, perform Shapiro-Wilk test')
        var<-paste('Variable name is',div_names[i])
        sw<-by(data=test_table, INDICES = test_table[,1], FUN=function(x){shapiro.test(x[,2])})
        tests<-append(tests, sw)
        
        
      }else if(n>=50){
        print('N<50, perform Lilliefors KS test')
        print(paste('Variable name is',div_names[i]))
        ks<-by(data=test_table, INDICES = test_table[,1], FUN=function(x){lillie.test(x[,2])})
        tests<-append(tests, ks)
        

        
        
      }else{print('Error, no test avaliable')}
      
      # Do qqplot
      qqp<-by(data=test_table, INDICES = test_table[,1], FUN=function(x){qqPlot(x[,2])})
    }
    
  }
  
  # Make the results table
  div_names2<-rep(div_names, each=g,  times=length(groups))
  cat_names<-names(tests)
  
  
  
  results<-tibble(Index=div_names2,
                  Group=cat_names,
                  Test=tests)
 

  results <- results %>% separate(Test, into = c("STAT", "PVAL", "METHOD", "DATA_NAME"), sep = ",")%>%
    select(-DATA_NAME)%>%
    separate(PVAL, into=c('X', 'p.val'), sep='=')%>%
    select(-X)%>%
    arrange(Index)
  results$p.val<-as.numeric(results$p.val)
  
  return(results)
  
}





check_homced<-function(df, test_groups, test_metrics){
  
  df_total<-df
  
  groups<-df_total%>%
    select(., all_of(test_groups))%>%
    colnames()
  
  div_names<-df_total%>%
    select(., all_of(test_metrics))%>%
    colnames()
  
  lv_p_val<-c()
  group<-c()
  div<-c()
  
  # Do the Levene's test
  for (j in 1:length(groups)){ 
    for (i in 1:length(div_names)){
      index<-as.character(div_names[i])
      fo<-as.formula(paste(index, '~',groups[j]))
      lv<-leveneTest(fo,data=df_total)
      p_val<-lv$`Pr(>F)`[1]
      
      group<-append(group, groups[j])
      div<-append(div, index)
      lv_p_val<-c(lv_p_val, p_val)
    }
  }
  #Create the results table
  results<-tibble(Group=group,
                  Index=div,
                  Levenne_p.val=lv_p_val)
  
  return(results)
}




#### WILCOXON RANK SUM TEST FUNCTION ############################################

wilcox_groups_test<-function(df, test_groups, test_metrics){
  df_total<-df
  
  groups<-df_total%>%
    select(., all_of(test_groups))%>%
    colnames()
  
  div_names<-df_total%>%
    select(., all_of(test_metrics))%>%
    colnames()
  
  p_val<-c()
  ci_l<-c()
  ci_h<-c()
  mw_u_val<-c()
  rankbis<-c()
  rankbis_ci_l<-c()
  rankbis_ci_h<-c()
  rankbis_int<-c()
  group<-c()
  div<-c()
  
  for (j in 1:length(groups)){ 
    for (i in 1:length(div_names)){
      print(paste('Wilcoxon rank-sum (Mann-Whitneys U) test between' ,
                  groups[j], 'and', div_names[i]))
      
      # Create the table of values
      test_table<-df_total%>%
        select(groups[j], div_names[i])%>%
        as.data.frame()%>%
        na.omit()
      print(test_table)
      
      group<-append(group, groups[j])
      div<-append(div, div_names[i])
      
      # Agregate them to show the median value by group  
      test_table_param<-aggregate(test_table[,2], 
                                  by = list(group=test_table[,1]), 
                                  median)
      
      colnames(test_table_param)<-c('Group', 'Median')
      print(test_table_param)
      
      
      # Perform the 'Wilcoxon rank-sum (Mann-Whitney?s U) test
      mw_u<-wilcox.test(test_table[,2]~test_table[,1], data = test_table,
                        exact = FALSE, conf.int = TRUE)
      print(mw_u)
      
      # Get the p-value, statistic and CI. Append then while iterating
      p_val<-append(p_val, mw_u$p.value)
      mw_u_val<-append(mw_u_val, mw_u$statistic)
      ci_l<-append(ci_l, mw_u$conf.int[1])
      ci_h<-append(ci_h, mw_u$conf.int[2])
      
      # Get the rank biserial size effect. Calculate the value, the CI and 
      # its interpretation. Append then inside their respective vectors while 
      # ?iterating
      rb<-rank_biserial(test_table[,2]~test_table[,1], data = test_table)
      rankbis<-append(rankbis, rb$r_rank_biserial)
      rankbis_ci_l<-append(rankbis_ci_l, rb$CI_low)
      rankbis_ci_h<-append(rankbis_ci_h, rb$CI_high)
      
      # The interpretation is done following the most recent criteria for it
      rankbis_int<-append(rankbis_int, 
                          interpret_rank_biserial( rb$r_rank_biserial,
                                                   rules = 'funder2019'))
      # Save the grouping variable name and the sympthom variable name in their
      # respective vectors while iterating
      
    }
  }
  results<-tibble(Group=group,
                  Index=div,
                  MW_Utest.p.value = p_val,
                  Statistic=mw_u_val,
                  MW_95_low_CI=ci_l,
                  MW_95_high_CI=ci_h,
                  Rank_biserial=rankbis,
                  Rank_biserial_95_low_CI=rankbis_ci_l,
                  Rank_biserial_95_high_CI=rankbis_ci_h,
                  Rank_biserial_meaning=rankbis_int)
  
  return(results)
}

####T TEST FUNCTION #############################################################

t_groups_test<-function(df, test_groups, test_metrics, var){
  df_total<-df
  
  groups<-df_total%>%
    select(., all_of(test_groups))%>%
    colnames()
  
  div_names<-df_total%>%
    select(., all_of(test_metrics))%>%
    colnames()
  
  ttestp_val<-c()
  ci_l_t<-c()
  ci_h_t<-c()
  ttest_val<-c()
  cohens_d<-c()
  cohens_d_ci_l<-c()
  cohens_d_ci_h<-c()
  cohensd_int<-c()
  group<-c()
  div<-c()
  
  
  for (j in 1:length(groups)){ 
    for (i in 1:length(div_names)){
      print(paste('T-test for mean comparison test between' ,
                  groups[j], 'and', div_names[i]))
      
      # Create the table of values
      test_table<-df_total%>%
        select(groups[j], div_names[i])%>%
        as.data.frame()%>%
        na.omit()
      
    print(test_table[, 1])
      
      # Agregate them to show the mean value by group  
      test_table_param<-aggregate(test_table[,2], 
                                  by = list(group=test_table[,1]), 
                                  mean)
      
      colnames(test_table_param)<-c('Group', 'Mean')
      print(test_table_param)
      
      
      # Perform the t test
      ttest<-t.test(test_table[,2]~test_table[,1], data = test_table,
                    exact = FALSE, conf.int = TRUE, var.equal=var)
      
      # Get the p-value, statistic and CI. Append then while iterating
      ttestp_val<-append(ttestp_val, ttest$p.value)
      ttest_val<-append(ttest_val, ttest$statistic)
      ci_l_t<-append(ci_l_t, ttest$conf.int[1])
      ci_h_t<-append(ci_h_t, ttest$conf.int[2])
      
      
      cd<-cohens_d(test_table[,2]~test_table[,1], data = test_table)
      cohens_d<-append(cohens_d, cd$Cohens_d)
      cohens_d_ci_l<-append(cohens_d_ci_l, cd$CI_low)
      cohens_d_ci_h<-append(cohens_d_ci_h, cd$CI_high)
      
      # The interpretation is done following the criteria for it
      cohensd_int<-append(cohensd_int, 
                          interpret_cohens_d(cd$Cohens_d,
                                             rules = 'cohen1988'))
      # Save the grouping variable name and the sympthom variable name in their
      # respective vectors while iterating
      group<-append(group, groups[j])
      div<-append(div, div_names[i])
      
    }
  }
  results<-tibble(Group=group,
                  Index=div,
                  ttestp_val=ttestp_val,
                  Statistic=ttest_val,
                  Ttest_95_low_CI=ci_l_t,
                  Test_95_high_CI=ci_h_t,
                  Cohen_d=cohens_d,
                  Cohen_d_95_low_CI=cohens_d_ci_l,
                  Cohen_d_95_high_CI=cohens_d_ci_h,
                  Cohen_d_meaning=cohensd_int)
  
  return(results)
}



#### KRUSKALL-WALLIS TEST ######################################################

kw_test<-function(df, test_groups, test_metrics){
  df_total<-df
  
  groups<-df_total%>%
    select(., all_of(test_groups))%>%
    colnames()
  
  div_names<-df_total%>%
    select(., all_of(test_metrics))%>%
    colnames()
  
  p_val<-c()
  ks_val<-c()
  res_v<-c()
  res_ci_l<-c()
  res_ci_h<-c()
  res_int<-c()
  post_hoc<-c()
  group<-c()
  div<-c()
  dt<-c()
  
  for (j in 1:length(groups)){ 
    for (i in 1:length(div_names)){
      print(paste('Kruskall-Wallis test between' ,
                  groups[j], 'and', div_names[i]))
      
      # Create the table of values
      test_table<-df_total%>%
        select(groups[j], div_names[i])%>%
        as.data.frame()%>%
        na.omit()
      
      # Agregate them to show the median value by group  
      test_table_param<-aggregate(test_table[,2], 
                                  by = list(group=test_table[,1]), 
                                  median)
      
      colnames(test_table_param)<-c('Group', 'Median')
      print(test_table_param)
      
      # Perform the Kruskall-Wallis test
      ks<-kruskal.test(test_table[,2]~test_table[,1])
      p_val<-append(p_val, ks$p.value)
      ks_val<-append(ks_val, ks$statistic)
      
      # Calculate size effect
      res<-rank_epsilon_squared(test_table[,2]~test_table[,1], data = test_table)
      res_v<-append(res_v, res$rank_epsilon_squared)
      res_ci_l<-append(res_ci_l, res$CI_low)
      res_ci_h<-append(res_ci_h, res$CI_high)
      
      # The interpretation is done following the criteria for it
      res_int<-append(res_int, 
                      interpret_epsilon_squared(res$rank_epsilon_squared))
      
      # Post hoc analysis with Dunn`s test and Bonferroni-Holm correction
      print("Dunn's test for post-hoc analysis")
      values<-test_table[,2]
      gr<-test_table[,1]
      df_temp<-data.frame(values=values, gr=gr)
      pt<-dunn_test(data=df_temp, values ~ gr,  p.adjust.method = 'BH')
      pt<-pt%>%mutate(., Metric=rep(div_names[i], 6), .before=group1)
      post_hoc<-bind_rows(post_hoc, pt, .id = NULL)
      
      
      group<-append(group, groups[j])
      div<-append(div, div_names[i])
      
    }
  }
  results_kw<-tibble(Group=group,
                     Index=div,
                     Kruskall_Wallis_p.val=p_val,
                     Statistic=ks_val,
                     `Rank_epsilon_squared(RES)`=res_v,
                     RES_95_low_CI=res_ci_l,
                     RES_d_95_high_CI=res_ci_h,
                     RES_meaning=res_int)
  
  results_posthoc<-compact(post_hoc)
  
  results<-list(Kruskall_Wallis=results_kw, Post_hoc=results_posthoc)
  
  return(results)
} 



#### ANOVA ######################################################


anova_test<-function(df, test_groups, test_metrics, .var){
  df_total<-df
  
  groups<-df_total%>%
    select(., all_of(test_groups))%>%
    colnames()
  
  div_names<-df_total%>%
    select(., all_of(test_metrics))%>%
    colnames()
  
  p_val<-c()
  f_val<-c()
  omega_v<-c()
  omega_ci_l<-c()
  omega_ci_h<-c()
  omega_int<-c()
  group<-c()
  div<-c()
  dt<-c()
  post_hoc<-tibble()
  
  
  for (j in 1:length(groups)){ 
    for (i in 1:length(div_names)){
      print(paste('ANOVA test between' ,
                  groups[j], 'and', div_names[i]))
      
      # Create the table of values
      test_table<-df_total%>%
        select(groups[j], div_names[i])%>%
        as.data.frame()%>%
        na.omit()
      
      # Agregate them to show the median value by group  
      test_table_param<-aggregate(test_table[,2], 
                                  by = list(group=test_table[,1]), 
                                  mean)
      
      colnames(test_table_param)<-c('Group', 'Mean')
      print(test_table_param)
      
      # Perform the ANOVA)
      res.aov<-oneway.test(test_table[,2]~test_table[,1], var.equal = .var)
      p_val<-append(p_val, res.aov$p.value)
      f_val<-append(f_val, res.aov$statistic)
      
      # Calculate size effect
      om<-omega_squared(res.aov)
      omega_v<-append(omega_v, om$Omega2)
      omega_ci_l<-append(omega_ci_l, om$CI_low)
      omega_ci_h<-append(omega_ci_h, om$CI_high)
      
      # The interpretation is done following the criteria for it
      omega_int<-append(omega_int, 
                        interpret_omega_squared(om$Omega2))
      
      
      
      # Post hoc analysis with pairwise t-test or Games Howell test
      
      if(.var==T){
        print("Pairwise t-test for post-hoc analysis for equal var ANOVA")
        values<-test_table[,2]
        gr<-test_table[,1]
        df_temp<-data.frame(values=values, gr=gr)
        pt<-pairwise_t_test(data=df_temp, values ~ gr,  p.adjust.method = 'BH')
        pt<-pt%>%mutate(., Metric_PTT=rep(div_names[i], 6), .before=group1)
        post_hoc<-bind_rows(post_hoc, pt, .id = NULL)
        
        
        
        
      }else{
        print("Games Howell Test for not equal var ANOVA")
        values<-test_table[,2]
        gr<-test_table[,1]
        df_temp<-data.frame(values=values, gr=gr)
        pt<-games_howell_test(data=df_temp, values ~ gr)
        pt<-pt%>%mutate(., Metric_GHT=rep(div_names[i], 6), .before=group1)
        post_hoc<-bind_rows(post_hoc, pt, .id = NULL)
        
        
      }
      
      group<-append(group, groups[j])
      div<-append(div, div_names[i])
      
    }
  }
  
  results_aov<-tibble(Group=group,
                      Index=div,
                      ANOVA_p.val=p_val,
                      Statistic=f_val,
                      `Omega-Squared (OS)`=omega_v,
                      OS_95_low_CI=omega_ci_l,
                      OS_d_95_high_CI=omega_ci_h,
                      OS_meaning=omega_int)
  
  results_posthoc<-post_hoc
  
  results<-list(ANOVA=results_aov, Post_hoc=results_posthoc)
  
  return(results)
}



##### Joinning functions #######################################################
test2<-function(df, test_groups, test_metrics, var=TRUE){
  ttest<-t_groups_test(df=df, test_groups=test_groups , test_metrics=test_metrics, var = var)
  mwu<-wilcox_groups_test(df=df, test_groups=test_groups , test_metrics=test_metrics)
  return(list(T_test=ttest, is.var.equal=var,  MWU=mwu))
  
}


test_more<-function(df, test_groups, test_metrics, var=TRUE){
  anova<-anova_test(df=df, test_groups=test_groups , test_metrics=test_metrics, .var = var)
  kw<-kw_test(df=df, test_groups=test_groups , test_metrics=test_metrics)
  return(list(ANOVA=anova, is.var.equal=var,  Kruskall_wallis=kw))
}
################################################################################