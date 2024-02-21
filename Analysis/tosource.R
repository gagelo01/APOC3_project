#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
library(tictoc)
library(furrr)

wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/APOC3"
setwd(wd)
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref<-"/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
x <- paste0("awk -F ' ' '{print $2}' ","/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs", ".bim")
dt_ref_snp <- fread(cmd = x, header = FALSE) #Those are SNP in 1000G with MAF>0.01
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao_small <- ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/"), ]
dt_gene_region <- fread("/mnt/sda/gagelo01/Projects/small_MR_exploration/Test_potential_project/PWMR_anything/Data/Modified/dt_gene_region.txt")

options(future.globals.maxSize= 5e9)
plan(multicore, workers = 40, gc = TRUE)

###chose instrument
study_to_selectinst <- c("deCODE", "FENLAND", "ARIC")
study_noinst_butexposure <- NULL
# study_noinst_butexposure <- c("IUCPQ Biobank", "GTEX")
should_skip_homogenous = TRUE
typeof_sel = c("lead_snp", "multicis_independent_clumping")
######change parameters#########
######
parameters <- GagnonMR::default_param()
parameters$path <- c(GagnonMR::default_param()$path, paste0(wd, "/Data/Modified/Gtex_vcf/"))
parameters$uni_cis_minpval <- 1
parameters$ldref <- ldref
parameters$snp_bim <- dt_ref_snp$V1 #I only include SNP in 1000G with maf > 0.01
parameters$multicis_clumping["clump_r2"]<-0.6



############Choose the gene you wish to include and the outcome you wish to include  ###########
gene_toinclude <-  c("APOC3", "APOA1", "APOA4", "APOA5") 
vec_tissue_gtex <- "Liver" 
############Choose the outcome you wish to include  ###########
# ID_server_out <- c(df_index[grepl("dis-15-",id) & ncase > 1000, id]) %>% unique(.)
ID_server_out <- c("dis-13-1", 
                   # df_index[grepl("dis-15-",id) & ncase > 1000, id],
                   df_index[grepl("trait-16", id) & grepl("^HDL|^LDL|^logTG", trait) & population %in% c("European") & sex == "Males and Females", id],
                    "trait-19-5",
                   "trait-7-1", "trait-6-1","dis-14-6", "dis-18-1", "dis-19-1",
                   "dis-2-2", "trait-14-8", "trait-14-10", paste0("trait-27-", 1:2), "trait-2-2", "trait-2-4", "dis-23-2",
                   "trait-12-2", "dis-7-1", "trait-25-16", "trait-25-1", "trait-29-20", "dis-12-1",
                  "trait-13-1", "trait-13-2", "dis-4-1", "dis-3-1")
out_server <- paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", ID_server_out, "/", ID_server_out, ".vcf.gz")
if(file.exists("Data/Modified/dt_gene_region.txt")) {
  dt_gene_region <- fread("Data/Modified/dt_gene_region.txt")
  out_server <- c(out_server, dt_gene_region[hgnc%in% gene_toinclude, vcffile])}
out_mrbase<-NULL

split_outcome<-FALSE
all_mr_methods_short = FALSE
all_mr_methods_skip_presso = FALSE