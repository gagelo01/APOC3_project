#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(openxlsx)
wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/APOC3"
setwd(wd)
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")
wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/ANGPTL4"
setwd(wd)

all_outcome <- fread("Data/Modified/all_outcome.txt" )
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")

####align on the good reference allele
k <- all_outcome[id.outcome == "trait-16-4", ] %>% GagnonMR::convert_outcome_to_exposure(.)
harm <- TwoSampleMR::harmonise_data(k, all_outcome)
tosel <- colnames(harm)[!grepl("exposure", colnames(harm))]
data_res <- harm[,tosel] %>% as.data.table(.)
cox_res <- fread("Data/Modified/cox_res.R")
cox_dat <- fread("Data/Modified/cox_data.R")
res_cophescan <- fread("Data/Modified/res_cophescan_hyperpriors.txt")
####
dataset <- df_index[id %in% unique(c(data_res$id.outcome)), ]
dataset[, c("group_name", "author", "consortium", "unit", "nsnp", "initial_build", "category", "sd", "note" ) := NULL]

suptab2 <- data_res[SNP=="rs116843064",]
suptab2[,c("remove", "palindromic", "ambiguous", "mr_keep.outcome", "pval_origin.outcome", "action", "mr_keep") := NULL]
suptab2 <- merge(suptab2, df_index[,.(id, clean_variable_name)], by.x = "id.outcome", by.y = "id")

suptab3 <- merge(res_cophescan, df_index[,.(id, trait, clean_variable_name)], by.x = "querytrait", by.y = "id")
setnames(suptab3, "trait", "outcome")


list_supdat <- list( "Supplementary Table 1" = dataset,
                     "Supplementary Table 2" = suptab2,
                     "Supplementary Table 3" = suptab3)


dt_title <- data.table(title = paste0("ST", c(1:3)),
                       caption = c( "Description of the datasets used for Two-Sample Mendelian randomization.",
                                    "Effect of the E40K on the human phenome",
                                    "Cophescan results for all tested traits"))


###
col_description<- vector(mode = "list", length = length(list_supdat))
col_description[[1]] <- tribble(
  ~x, ~y,  
  "id", "a unique id", 
  "trait", "a unique name",
  "year", "the year the GWAS was published", 
  "trait", "The author of the GWAS",
  "consortium", "the name of the consortium of the GWAS",
  "sex", "sex included",
  "population", "ancestry",
  "sample_size", "the sample_size",
  "pmid", "the pubmed ID",
  "ncase", "the number of cases",
  "ncontrol", "the number of controls",
  "adjustments", "what variables were included in the GWAS adjustment model",
  "clean_variable_name", "A publication ready name that can be used to plot figures."
) %>% as.data.table(.)

col_description[[2]] <- tribble(
  ~x, ~y,  
  "id.outcome", "a unique id", 
  "chr.outcome", "chromosome",
  "pos.outcome", "Position build 37", 
  "other_allele.outcome", "The non-effect allele",
  "effect_allele.outcome", "The effect allele",
  "beta.outcome", "beta effect estimate",
  "se.outcome", "standard error of the estimate",
  "pval.outcome", "p-valueor of the estimate",
  "eaf.outcome", "effect allele frequency",
  "samplesize.outcome", "sample size",
  "ncase.outcome", "number of cases",
  "SNP", "rsid",
  "ncontrol.outcome", "number of controls",
  "outcome", "A unique name for the outcome",
  "clean_variable_name", "A publication ready name that can be used to plot figures."
) %>% as.data.table(.)

col_description[[3]] <- tribble(
  ~x, ~y,  
  "PP.Hn ", "posterior probability of no association with the query trait", 
  "PP.Ha", "posterior probability of association with a different variant than the query variant",
  "PP.Hc", "posterior probabilty of association with querry variant with the query trait", 
  "nsnps", "the number of snps included",
  "querysnp", "The rsid for the query snp",
  "querytrait", "the id for the qeury trait",
  "typeBF", "which algorithm was used. Either SuSiE coloc (multiple causal variant may exist in that locus) or coloc (single variant assumption)",
  "cophe.hyp.call", "Which hypothesis was prioritised by CoPheScan",
  "outcome", "A unique name for the outcome",
  "clean_variable_name", "A publication ready name that can be used to plot figures."
) %>% as.data.table(.)

bold_st <- createStyle(textDecoration = "Bold")
wb <- createWorkbook()
for(i in 1:length(list_supdat)) {
  addWorksheet(wb, sheetName =  dt_title[i, title])
  
  title <- dt_title[i,paste0(title, " : ", caption)]
  writeData(wb, sheet = i, x = title, startCol = 1, startRow = 1)
  addStyle(wb, sheet = i, style = bold_st, rows = 1, cols = 1:2)
  writeData(wb, sheet = i, x = col_description[[i]], startCol = 1, startRow = 2)
  addStyle(wb, sheet = i, style = bold_st, rows = 2:col_description[[i]][,.N+2], cols = 1)
  deleteData(wb, sheet = i, rows = 2, cols = 1:2, gridExpand = TRUE)
  writeData(wb, sheet = i, x = list_supdat[[i]], startCol = 1, startRow = col_description[[i]][,.N+4])
  addStyle(wb, sheet = i, style = bold_st, rows = col_description[[i]][,.N+4], cols = 1:length(colnames(list_supdat[[i]])), gridExpand = TRUE, stack = TRUE)
}
saveWorkbook(wb, file = "Results/supplementary_tables.xlsx", overwrite = TRUE)

message("This script finished without errors")

