#!/usr/bin/env Rscript
#####Section to CHange ##########
setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/APOC3")
mysource <- "tosource.R"

myArray=list.files("/mnt/sda/gagelo01/Projects/Pipelines/Robust_drug_pipeline/Analysis")
myArray=myArray[grepl("^[1-9]", myArray) & grepl("\\.R$", myArray)]
myArray <- myArray[c(1,2,3, 4,5,9)]
myArray <- paste0("/mnt/sda/gagelo01/Projects/Pipelines/Robust_drug_pipeline/Analysis/", myArray)
##################################
bloc1<- paste0("myArray=(", paste(myArray, collapse = " "), ")")

bloc2 <- "for str in ${myArray[@]}; do
chmod u+x $str
done"

bloc3 <- paste(paste(paste0("echo 'Initializing ", myArray, "' && ", myArray, " ", mysource," &&"), collapse = " "), "echo 'The master script finished without errors'")

fileConn<-file("Analysis/drugpipeline.sh")
writeLines(paste(c("#!/bin/bash",bloc1, bloc2, bloc3), collapse = "\n"), fileConn)
close(fileConn)
