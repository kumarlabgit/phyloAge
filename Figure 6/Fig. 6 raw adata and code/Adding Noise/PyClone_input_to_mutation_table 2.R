################################################################################
#Startup

rm(list = ls())

setwd("C:/Users/tuo47278/Desktop/make_mut_tables")

files <- list.files(path="C:/Users/tuo47278/Desktop/make_mut_tables", pattern="*.tsv", full.names=TRUE, recursive=FALSE)
fileNames <- Sys.glob("*.tsv")



for (fileName in fileNames) {
  
  #fileName <- fileNames[1]
  
  df <- read.table(file = fileName, header = TRUE)

  Position <- df$mutation_id
  Wild <- df$ref_counts
  Mut <- df$var_counts
  Miss <- 0  
  
  mut_table <- data.frame(cbind(Position, Wild, Mut, Miss))
  colnames(mut_table) <- c("Position", "Wild", "Mut", "Miss")
  
  filename_to_write <- paste0(fileName, "_mut_table.csv")
  write.csv(mut_table, file = filename_to_write, row.names=FALSE)

} #closes "for (fileName in fileNames)"


################################################################################
