library(tidyverse)
library(FactoMineR)
library(factoextra)
library(forcats)
library(ggpubr)
library(gridExtra)
library(smacof)
library(NbClust)


div<-read.csv("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/diversity.csv",
              header=T, sep=',')[,-c(1, c(3:8))]
clon<-read.csv("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/clonality.csv",
               header=T, sep=',')[,-c(1, c(3:8))]
vus<-read.csv("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/v_usage_prop_w.csv",
              header=T, sep=',')[,-c(1, c(3:8))]%>%replace(is.na(.), 0)

jus<-read.csv("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/j_usage_prop_w.csv",
              header=T, sep=',')[,-c(1, c(3:8))]
cdr3<-read.csv("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/t_cdr3_table3_freq_mira.csv",
               header=T, sep=',')[,-c(1, c(3:8))]
meta<-read.csv("C:/Users/Usuario/OneDrive/Escritorio/covid_integration/data/metadata.csv",
               header = T, sep=';')[, c(1:12)]

dfs <- list(meta, div, clon, vus, jus, cdr3)
tcr_df <- reduce(dfs, full_join, by = "sample")

tcr_scale_df <- tcr_df %>%
  mutate_if(is.numeric, function(x) (x - median(x)) / IQR(x))%>%
  dplyr::select(-SPT.AGE_ASV)

tcr_matrix<-tcr_scale_df%>%select_if(., is.numeric)%>%select(, -Age)
rownames(tcr_matrix)<-tcr_scale_df$sample





###### Multiple logistic regression #####

library(MASS)
library(nnet)
### Multinomial

df_regression<-tcr_df%>%
  select(., Severity, sev_55, Morethan55, Sex, sex_sev, which(colnames(tcr_df) %in% c(colnames(cdr3))))
df_regression_multinom <- df_regression[ , !(names(df_regression) %in% c("Severity", "Morethan55", "Sex", "sex_sev", "sample"))]
as.factor(df_regression_multinom$sev_55)
complete_multinom_model<-multinom(sev_55 ~ ., data=df_regression_multinom)
selected_multinom_model<-stepAIC(complete_multinom_model, direction='both')
summary(selected_multinom_model)





#### Binomial: Severe;<55 againts each sev_55 subgroup ####


### vs Mild:<55 

s55_1.df_regression_binom <- df_regression_multinom %>%
  filter(sev_55 %in% c('Severe;<55', 'Mild;<55'))
s55_1.df_regression_binom$sev_55 <- ifelse(s55_1.df_regression_binom$sev_55 == 'Severe;<55', 1, 0)
s55_1.df_regression_binom$sev_55<-as.factor(s55_1.df_regression_binom$sev_55)


set.seed(1234)

fit_model_and_metrics <- function(data, test_data) {
  s55_1.complete_binom_model <- glm(sev_55 ~ ., data = data, family = binomial)
  s55_1.selected_binom_model <- stepAIC(s55_1.complete_binom_model, direction = "both")
  
  predicted <- stats::predict(s55_1.selected_binom_model, test_data, type = "response")
  rocobj <- pROC::roc(test_data$sev_55, predicted)
  auc_value <- pROC::auc(rocobj)
  selected_variables <- names(coef(s55_1.selected_binom_model)[-1])
  
  return(list(roc = rocobj, auc = auc_value, selected_variables = selected_variables))
}

num_simulations <- 1000
auc_values <- numeric(num_simulations)
roc_list <- vector("list", num_simulations)
selected_variables_list <- vector("list", num_simulations)

for (i in 1:num_simulations) {
  sample <- sample(c(TRUE, FALSE), nrow(s55_2.df_regression_binom), replace = TRUE, prob = c(0.75, 0.25))
  train <- s55_2.df_regression_binom[sample, ]
  test <- s55_2.df_regression_binom[!sample, ]
  
  result <- fit_model_and_metrics(train, test)
  auc_values[i] <- result$auc
  roc_list[[i]] <- result$roc
  selected_variables_list[[i]] <- result$selected_variables
}

mean_auc <- mean(auc_values)
cat("Mean AUC:", mean_auc, "\n")

plot(roc_list[[1]], col = rgb(0, 0, 1, alpha = 0.2), main = "ROC curve for Severe;<55 vs Mild;<55")
for (i in 2:num_simulations) {
  lines(roc_list[[i]], col = rgb(0, 0, 1, alpha = 0.05))
}

abline(a = 1, b = -1, col = "red", lwd = 2)

all_selected_variables <- unique(unlist(selected_variables_list))
variable_counts <- sapply(all_selected_variables, function(x) sum(sapply(selected_variables_list, function(y) x %in% y)))
sorted_variables <- names(sort(variable_counts, decreasing = TRUE))
for (variable in sorted_variables) {
  cat(variable, ":", variable_counts[variable], "\n")
}





### vs Mild;>55



s55_2.df_regression_binom <- df_regression_multinom %>%
  filter(sev_55 %in% c('Severe;<55', 'Mild;>55'))
s55_2.df_regression_binom$sev_55 <- ifelse(s55_2.df_regression_binom$sev_55 == 'Severe;<55', 1, 0)
s55_2.df_regression_binom$sev_55<-as.factor(s55_2.df_regression_binom$sev_55)


set.seed(1234)

fit_model_and_metrics <- function(data, test_data) {
  s55_2.complete_binom_model <- glm(sev_55 ~ ., data = data, family = binomial)
  s55_2.selected_binom_model <- stepAIC(s55_2.complete_binom_model, direction = "both")
  
  predicted <- stats::predict(s55_2.selected_binom_model, test_data, type = "response")
  rocobj <- pROC::roc(test_data$sev_55, predicted)
  auc_value <- pROC::auc(rocobj)
  selected_variables <- names(coef(s55_2.selected_binom_model)[-1])
  
  return(list(roc = rocobj, auc = auc_value, selected_variables = selected_variables))
}

num_simulations <- 100
auc_values <- numeric(num_simulations)
roc_list <- vector("list", num_simulations)
selected_variables_list <- vector("list", num_simulations)

for (i in 1:num_simulations) {
  sample <- sample(c(TRUE, FALSE), nrow(s55_2.df_regression_binom), replace = TRUE, prob = c(0.75, 0.25))
  train <- s55_2.df_regression_binom[sample, ]
  test <- s55_2.df_regression_binom[!sample, ]
  
  result <- fit_model_and_metrics(train, test)
  auc_values[i] <- result$auc
  roc_list[[i]] <- result$roc
  selected_variables_list[[i]] <- result$selected_variables
}

mean_auc <- mean(auc_values)
cat("Mean AUC:", mean_auc, "\n")

plot(roc_list[[1]], col = rgb(0, 0, 1, alpha = 0.2), main = "ROC curve for Severe;<55 vs Mild;>55")
for (i in 2:num_simulations) {
  lines(roc_list[[i]], col = rgb(0, 0, 1, alpha = 0.05))
}

abline(a = 1, b = -1, col = "red", lwd = 2)

all_selected_variables <- unique(unlist(selected_variables_list))
variable_counts <- sapply(all_selected_variables, function(x) sum(sapply(selected_variables_list, function(y) x %in% y)))
sorted_variables <- names(sort(variable_counts, decreasing = TRUE))
for (variable in sorted_variables) {
  cat(variable, ":", variable_counts[variable], "\n")
}





# vs Severe;>55


s55_3.df_regression_binom <- df_regression_multinom %>%
  filter(sev_55 %in% c('Severe;<55', 'Severe;>55'))
s55_3.df_regression_binom$sev_55 <- ifelse(s55_3.df_regression_binom$sev_55 == 'Severe;<55', 1, 0)
s55_3.df_regression_binom$sev_55<-as.factor(s55_3.df_regression_binom$sev_55)



set.seed(1234)

fit_model_and_metrics <- function(data, test_data) {
  s55_3.complete_binom_model <- glm(sev_55 ~ ., data = data, family = binomial)
  s55_3.selected_binom_model <- stepAIC(s55_3.complete_binom_model, direction = "both")
  
  predicted <- stats::predict(s55_3.selected_binom_model, test_data, type = "response")
  rocobj <- pROC::roc(test_data$sev_55, predicted)
  auc_value <- pROC::auc(rocobj)
  selected_variables <- names(coef(s55_3.selected_binom_model)[-1])
  
  return(list(roc = rocobj, auc = auc_value, selected_variables = selected_variables))
}

num_simulations <- 100
auc_values <- numeric(num_simulations)
roc_list <- vector("list", num_simulations)
selected_variables_list <- vector("list", num_simulations)

for (i in 1:num_simulations) {
  sample <- sample(c(TRUE, FALSE), nrow(s55_3.df_regression_binom), replace = TRUE, prob = c(0.75, 0.25))
  train <- s55_3.df_regression_binom[sample, ]
  test <- s55_3.df_regression_binom[!sample, ]
  
  result <- fit_model_and_metrics(train, test)
  auc_values[i] <- result$auc
  roc_list[[i]] <- result$roc
  selected_variables_list[[i]] <- result$selected_variables
}

mean_auc <- mean(auc_values)
cat("Mean AUC:", mean_auc, "\n")

plot(roc_list[[1]], col = rgb(0, 0, 1, alpha = 0.2), main = "ROC curve for Severe;<55 vs Severe;>55 variable")
for (i in 2:num_simulations) {
  lines(roc_list[[i]], col = rgb(0, 0, 1, alpha = 0.05))
}

abline(a = 1, b = -1, col = "red", lwd = 2)

all_selected_variables <- unique(unlist(selected_variables_list))
variable_counts <- sapply(all_selected_variables, function(x) sum(sapply(selected_variables_list, function(y) x %in% y)))
sorted_variables <- names(sort(variable_counts, decreasing = TRUE))
for (variable in sorted_variables) {
  cat(variable, ":", variable_counts[variable], "\n")
}








### Binomial of Severity and Morethan55
set.seed(1234)

df_regression_binom <- df_regression[ , !(names(df_regression) %in% c("sev_55", "Morethan55", "Sex", "sex_sev", "sample"))]
df_regression_binom$Severity <- ifelse(df_regression_binom$Severity == "Severe", 1, 0)
df_regression_binom$Severity<-as.factor(df_regression_binom$Severity)


fit_model_and_metrics <- function(data, test_data) {
  complete_binom_model <- glm(Severity ~ ., data = data, family = binomial)
  selected_binom_model <- stepAIC(complete_binom_model, direction = "both")
  
  predicted <- stats::predict(selected_binom_model, test_data, type = "response")
  rocobj <- pROC::roc(test_data$Severity, predicted)
  auc_value <- pROC::auc(rocobj)
  selected_variables <- names(coef(selected_binom_model)[-1])
  
  return(list(roc = rocobj, auc = auc_value, selected_variables = selected_variables))
}

num_simulations <- 100
auc_values <- numeric(num_simulations)
roc_list <- vector("list", num_simulations)
selected_variables_list <- vector("list", num_simulations)

for (i in 1:num_simulations) {
  sample <- sample(c(TRUE, FALSE), nrow(df_regression_binom), replace = TRUE, prob = c(0.75, 0.25))
  train <- df_regression_binom[sample, ] 
  test <- df_regression_binom[!sample, ]
  
  result <- fit_model_and_metrics(train, test)
  auc_values[i] <- result$auc
  roc_list[[i]] <- result$roc
  selected_variables_list[[i]] <- result$selected_variables
}

mean_auc <- mean(auc_values)
cat("Mean AUC:", mean_auc, "\n")

plot(roc_list[[1]], col = rgb(0, 0, 1, alpha = 0.2), main = "ROC curve for Mild vs Severe")
for (i in 2:num_simulations) {
  lines(roc_list[[i]], col = rgb(0, 0, 1, alpha = 0.05))
}

abline(a = 1, b = -1, col = "red", lwd = 2)

all_selected_variables <- unique(unlist(selected_variables_list))
variable_counts <- sapply(all_selected_variables, function(x) sum(sapply(selected_variables_list, function(y) x %in% y)))
sorted_variables <- names(sort(variable_counts, decreasing = TRUE))
for (variable in sorted_variables) {
  cat(variable, ":", variable_counts[variable], "\n")
}










df_regression_binom <- df_regression[ , !(names(df_regression) %in% c("sev_55", "Severity", "Sex", "sex_sev", "sample"))]
df_regression_binom$Severity <- ifelse(df_regression_binom$Severity == ">55", 1, 0)
df_regression_binom$Severity<-as.factor(df_regression_binom$Severity)

set.seed(1234)

fit_model_and_metrics <- function(data, test_data) {
  m55.complete_binom_model <- glm(Morethan55 ~ ., data = data, family = binomial)
  m55.selected_binom_model <- stepAIC(m55.complete_binom_model, direction = "both")
  
  predicted <- stats::predict(m55.selected_binom_model, test_data, type = "response")
  rocobj <- pROC::roc(test_data$Morethan55, predicted)
  auc_value <- pROC::auc(rocobj)
  selected_variables <- names(coef(m55.selected_binom_model)[-1])
  
  return(list(roc = rocobj, auc = auc_value, selected_variables = selected_variables))
}

num_simulations <- 100
auc_values <- numeric(num_simulations)
roc_list <- vector("list", num_simulations)
selected_variables_list <- vector("list", num_simulations)

for (i in 1:num_simulations) {
  sample <- sample(c(TRUE, FALSE), nrow(m55.df_regression_binom), replace = TRUE, prob = c(0.75, 0.25))
  train <- m55.df_regression_binom[sample, ]
  test <- m55.df_regression_binom[!sample, ]
  
  result <- fit_model_and_metrics(train, test)
  auc_values[i] <- result$auc
  roc_list[[i]] <- result$roc
  selected_variables_list[[i]] <- result$selected_variables
}

mean_auc <- mean(auc_values)
cat("Mean AUC:", mean_auc, "\n")

plot(roc_list[[1]], col = rgb(0, 0, 1, alpha = 0.2), main = "ROC curve for <55 vs >=55 age groups")
for (i in 2:num_simulations) {
  lines(roc_list[[i]], col = rgb(0, 0, 1, alpha = 0.05))
}

abline(a = 1, b = -1, col = "red", lwd = 2)

all_selected_variables <- unique(unlist(selected_variables_list))
variable_counts <- sapply(all_selected_variables, function(x) sum(sapply(selected_variables_list, function(y) x %in% y)))
sorted_variables <- names(sort(variable_counts, decreasing = TRUE))
for (variable in sorted_variables) {
  cat(variable, ":", variable_counts[variable], "\n")
}















