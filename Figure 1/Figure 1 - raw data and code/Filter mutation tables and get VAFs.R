################################################################################
#Startup

rm(list = ls())

setwd("C:/Users/tuo47278/Desktop/")

raw_df <- read.csv(file = "./PD5117_WildMutCount.csv")

#define the subset size (0.01 = 1% of the data will be kept, 99% dropped)
#can be used to test robustness of lambdas to sparse data
#if not subsetting, leave at 1
subset <- 1

################################################################################

#create a new df to work with
df <- raw_df


#Subset (optional; won't do anything if subset is set to 1)

#select rows to keep
rows_to_keep <- round(subset * nrow(df))
rows_to_keep_indices <- sample(1:nrow(df), size = rows_to_keep)

#subset the df
df <- df[rows_to_keep_indices, ]


#Filter

#filter out sites with only one mutant read (accounts for erroneous reads)
df <- df[df$Mut > 1, ] 

#filter out sites with less than 30 total reads (accounts for low coverage sites)
#this must be skipped when average rd <30, as in several Williams individuals and true bulk
df <- df[df$Mut+df$Wild > 30, ]  

#calculate VAF
VAF <- (df$Mut/(df$Mut + df$Wild))/2

#make output df
out <- data.frame(cbind(df$Wild,df$Mut,VAF))
colnames(out) <- c("ref_counts", "var_counts", "VAF")

#filter out very low-freq sites (accounts for ancestral lineages)
out <- out[out$VAF > 0.01, ]


#write the filtered mutation table with VAFs
write.csv(out, file = "filtered_mutation_table.csv", row.names=FALSE,  quote = FALSE)

################################################################################