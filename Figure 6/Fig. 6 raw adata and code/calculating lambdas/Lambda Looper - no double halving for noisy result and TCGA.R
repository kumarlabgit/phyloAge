################################################################################
#Startup

rm(list = ls())

setwd("C:/Users/tuo47278/Desktop/mut_tables")

files <- list.files(path="C:/Users/tuo47278/Desktop/mut_tables", pattern="*.csv", full.names=TRUE, recursive=FALSE)
fileNames <- Sys.glob("*.csv")

lambdaTable <- data.frame(matrix(nrow = 0, ncol = 2))

for (fileName in fileNames) {
  
  #fileName <- fileNames[1]
  
  df <- read.csv(file = fileName)
  
  
  #filter out sites with only one mutant read (accounts for erroneous reads)
  df <- df[df$Mut > 1, ] 
  
  #filter out sites with less than 30 total reads (accounts for low coverage sites)
  #this must be skipped when average rd <30, as in several Williams individuals and true bulk
  #df <- df[df$Mut+df$Wild > 30, ]  
  
  #calculate VAF
  #does not include the extra halving here, which we need for Mitchell/Williams, but not noisy or TCGA
  VAF <- (df$Mut/(df$Mut + df$Wild))
  
  #make output df
  out <- data.frame(cbind(df$Wild,df$Mut,VAF))
  colnames(out) <- c("ref_counts", "var_counts", "VAF")
  
  #filter out very low-freq sites (accounts for ancestral lineages)
  out <- out[out$VAF > 0.01, ]
  
  
  #calculate lambda
  lambda <- sum(out$VAF)
  
  new_row <- data.frame(matrix(nrow = 1, ncol = 2))
  new_row[,1] <- fileName
  new_row[,2] <- lambda  
  
  lambdaTable <- rbind(lambdaTable, new_row)
  
} #closes "for (fileName in fileNames)"


colnames(lambdaTable) <- c("File", "lambda")


write.csv(lambdaTable, "lambdas.csv", row.names = FALSE)

