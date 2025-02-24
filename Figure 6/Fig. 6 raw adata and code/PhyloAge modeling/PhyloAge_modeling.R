################################################################################
#Startup

rm(list = ls())

library(plotrix)
library(metafor)
library(investr)

setwd("C:/Users/tuo47278/Desktop")


################################################################################
#Train the model (healthy people)

df1 <- read.csv("Train_data_30X_coding.csv", header = T, sep = ",")

ft <- function(x, a, b) {
  (1/b)*(log(x) - a)
} 


#fit models

fit1 <- nls(Age ~ ft(gamma, a, b),  start = list(a =(1/1000), b = .3), 
            trace = T, control = list(maxiter = 50), data = df1)

pat <- paste("pa", 1:length(df1$Individual), sep = "")
df_meta <- rbind(cbind(as.numeric(predFit(fit1, se.fit = T)[[1]]), as.numeric(predFit(fit1, se.fit = T)[[2]]), pat))


colnames(df_meta) <- c("d", "sei", "p")
df_meta <- data.frame(df_meta)


#Meta-analysis

train_out <- NULL

for(i in 1:length(pat)){
  a <- rma(yi = as.numeric(d), sei = as.numeric(sei), method="REML", weighted=T, data = df_meta[df_meta$p == pat[i], ], level=95, intercept=T)
  train_out <- rbind(train_out, c(a$beta, a$se, a$zval, a$pval, a$ci.lb, a$ci.ub))
  print(c(a$beta, a$se, a$ci.lb, a$ci.ub))
}

colnames(train_out) <- c("estimate_Age", "SE", "z_val", "p_val", "low_CI", "high_CI")

write.csv(train_out, file = "training_result.csv")


################################################################################
#Test the model (MPN people)

test <- read.csv("Test_data_TCGA.csv", header = T, sep = ",")
test1 <- as.data.frame(test$gamma)
colnames(test1) <- "gamma"


#PhyloAge prediction
pat <- paste("pa", 1:length(test$Individual), sep = "")
df_meta2 <- rbind(cbind(as.numeric(predFit(fit1, test1, se.fit = T)[[1]]), as.numeric(predFit(fit1, test1, se.fit = T)[[2]]), pat))

colnames(df_meta2) <- c("d", "sei", "p")
df_meta2 <- data.frame(df_meta2)


#Meta-analysis

test_out <- NULL
for(i in 1:length(pat)){
  a <- rma(yi = as.numeric(d), sei = as.numeric(sei), method="REML", weighted=T, data = df_meta2[df_meta2$p == pat[i], ], level=95, intercept=T)
  test_out <- rbind(test_out, c(a$beta, a$se, a$zval, a$pval, a$ci.lb, a$ci.ub))
  print(c(a$beta, a$se, a$ci.lb, a$ci.ub))
}

colnames(test_out) <- c("estimate_Age", "SE", "z_val", "p_val", "low_CI", "high_CI")

write.csv(test_out, file = "testing_result.csv")


################################################################################
#Validate the model by leave-one-out analysis

loo_out <- NULL

for(j in 1:1:length(df1$Individual)){
  df_loo <- df1[-j,]
  
    fit1 <- nls(Age ~ ft(gamma, a, b),  start = list(a= 0.001, b = 0.03), 
               trace = T, control = list(maxiter = 500), data = df_loo)
  
  test1 <- as.data.frame(df1[j,3])
  colnames(test1) <- "gamma"
  
  pat <- paste("pa", 1, sep = "")
  df_meta2 <- rbind(cbind(as.numeric(predFit(fit1, test1, se.fit = T)[[1]]), as.numeric(predFit(fit1, test1, se.fit = T)[[2]]), pat))
  
  colnames(df_meta2) <- c("d", "sei", "p")
  df_meta2 <- data.frame(df_meta2)
  
  a <- rma(yi = as.numeric(d), sei = as.numeric(sei), method="REML", weighted=T, data = df_meta2, level=95, intercept=T)
  print(c(df1[j,2],a$beta, a$se, a$ci.lb, a$ci.ub))
  loo_out <- rbind(loo_out, c(df1[j,2],a$beta, a$se, a$ci.lb, a$ci.ub))
  
}

write.csv(loo_out, file = "LOO_result.csv")
