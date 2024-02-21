#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(furrr)
library(gassocplot)
library(MetBrewer)
library("RColorBrewer")
library(ggh4x)
library(forestploter)
wd<-"/mnt/sda/gagelo01/Projects/small_MR_exploration/APOC3/"
setwd(wd)
gwasvcf::set_bcftools()
gwasvcf::set_plink()
dt_gene_region <- fread("Data/Modified/dt_gene_region.txt")
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
inst <- fread("Data/Modified/inst.txt")
inst <- inst[typeof_sel == "multicis_independent_clumping", ]
inst <- merge(inst, dt_gene_region[,.(id,  hgnc, study)], by.x = "id.exposure", by.y = "id")
all_outcome <- fread("Data/Modified/all_outcome.txt")
resmvmr <- fread( "Data/Modified/res_mvmr.txt")
resmvmr <- merge(resmvmr, distinct(res[,.(id.exposure, id.outcome, clean_variable_name, mechanism)]), by=c("id.exposure", "id.outcome"))

res <-fread("Data/Modified/res_multicis_independent.txt")
res <- merge(res, dt_gene_region[,.(id,  hgnc, study)], by.x = "id.exposure", by.y = "id")
res <- merge(res, df_index[,.(id,  clean_variable_name)], by.x = "id.outcome", by.y = "id")
res[, mechanism := ifelse(study %in% c("ARIC", "FENLAND", "deCODE"),
                             "Blood protein levels",
                             "Liver gene levels")]
dt<- rbindlist(list(
  data.table(nom = c(paste0("trait-2-", 1:6), paste0("dis-6-", 1:3), "dis-23-2"), value = "Glucose homeostasis"),
  data.table(nom = c(paste0("dis-13-", 1:4), "dis-5-1", df_index[grepl("dis-14-", id), id], paste0("trait-13-", 1:2), "dis-18-1", "dis-19-1", "dis-1-1", "trait-26-1"), value = "Vascular"),
  data.table(nom = c("dis-2-2", "trait-14-8", paste0("trait-27-", 1:2)), value = "Liver"),
  data.table(nom = c("dis-7-1", "trait-12-2"), value = "Kidney"),
  data.table(nom = c(unique(df_index[grepl("trait-16-", id), id]), "trait-19-5","trait-19-2", "trait-29-20") , value = "Lipids"),
  data.table(nom =   paste0("trait-28-", 1:2), value = "Bone"),
  data.table(nom = c(unique(df_index[grepl("met-", id), id]), "trait-20-8"), value = "Metabolites"),
  data.table(nom = c("trait-1-1", "trait-10-1", "trait-10-3", "trait-25-16", "trait-25-1"), value = "Anthropometric"),
  data.table(nom = c("trait-6-1", "trait-7-1"), value = "Lifespan"),
  data.table(nom = c("dis-12-1",  "trait-14-10"), value = "Pancreas"),
  data.table(nom = c("dis-4-1"), value = "Brain"),
  data.table(nom = res[grepl("prot-", id.outcome), id.outcome], value = "Protein")
))
res[id.outcome %in% setdiff(unique(res$id.outcome), dt$nom), .(id.outcome, outcome, clean_variable_name)] %>% distinct(.)

res <- merge(res, distinct(dt), by.x = "id.outcome", by.y = "nom")

#####forest plot########
dt_forest <-  data.table(res)
dt_forest <- dt_forest[method=="Inverse variance weighted", ]
dt_forest[,b:=b*-1]
dt_forest[,lci:=b-1.96*se]
dt_forest[,uci:=b+1.96*se]
dt_forest <- dt_forest[order(clean_variable_name, mechanism),]
dt_forest <- dt_forest[method == "Inverse variance weighted" & grepl("eqtl-74|prot-5", id.exposure),] #GTEX gives way cleaner results
setnames(dt_forest, "mechanism", "panel")

df_index[id %in% dt_forest$id.outcome, .(clean_variable_name, trait, consortium, id)][order(clean_variable_name)] %>% distinct(.)

dt_forest <- dt_forest[!(id.outcome%in%"trait-19-2"), ]
p <- my_forestplotter_fancy(data = dt_forest[grepl("prot-", id.exposure) & hgnc == "APOC3" & (grepl("dis-", id.outcome) | id.outcome == "trait-6-1"), ], 
                            col.below.header = "clean_variable_name", 
                            col.header = "value", 
                            col.header.heading = "", 
                            effect.name = "OR (95% CI)", 
                            col.right = NULL,
                            exponentiate = TRUE,
                            xlab = "Effect of 1-SD reduction of APOC3")

ggsave(filename = paste0("Results/Figure1_APOC3.png"),plot = p, dpi =300,
       width = 520/72,height = 320/72,units="in",scale=1, device = "png")

p <- my_forestplotter_fancy(data = dt_forest[grepl("prot-", id.exposure) & hgnc == "APOC3" & (grepl("trait-", id.outcome) & id.outcome != "trait-6-1"), ], 
                            col.below.header = "clean_variable_name", 
                            col.header = "value", 
                            col.header.heading = "", 
                            effect.name = "Beta (95% CI)", 
                            col.right = NULL,
                            exponentiate = FALSE,
                            xlab = "Effect of 1-SD reduction of APOC3")

ggsave(filename = paste0("Results/Figure2_APOC3.png"),plot = p, dpi =300,
       width = 500/72,height = 459/72,units="in",scale=1, device = "png")

# p <- my_forestplotter_fancy(data = dt_forest[grepl("eqtl-", id.exposure) & hgnc == "APOC3" & (grepl("dis-", id.outcome) | id.outcome == "trait-6-1"), ], 
#                             col.below.header = "clean_variable_name", 
#                             col.header = "value", 
#                             col.header.heading = "", 
#                             effect.name = "OR (95% CI)", 
#                             col.right = NULL,
#                             exponentiate = TRUE,
#                             xlab = "Effect of 1-SD reduction of APOC3")
# 
# ggsave(filename = paste0("Results/Figure5.png"),plot = p, dpi =300,
#        width = 520/72,height = 320/72,units="in",scale=1, device = "png")
# 
# p <- my_forestplotter_fancy(data = dt_forest[grepl("eqtl-", id.exposure) & hgnc == "APOC3" & (grepl("trait-", id.outcome) & id.outcome != "trait-6-1"), ], 
#                             col.below.header = "clean_variable_name", 
#                             col.header = "value", 
#                             col.header.heading = "", 
#                             effect.name = "Beta (95% CI)", 
#                             col.right = NULL,
#                             exponentiate = FALSE,
#                             xlab = "Effect of 1-SD reduction of APOC3")
# 
# ggsave(filename = paste0("Results/Figure6.png"),plot = p, dpi =300,
#        width = 558/72,height = 509/72,units="in",scale=1, device = "png")

#####Multivariable MR#####
mvmr_forest <- rbind(res, resmvmr, fill = TRUE)
mvmr_forest2 <- mvmr_forest[method %in% c("Multivariable IVW", "Inverse variance weighted") & study == "deCODE", ]
mvmr_forest2[,colheader := ifelse(method == "Inverse variance weighted", "Univariable", "Multivariable")]
mvmr_forest2[,colheader:=factor(colheader, levels = c("Univariable", "Multivariable"))]
mvmr_forest2[,b:=b*-1]
mvmr_forest2[,lci:=b-1.96*se]
mvmr_forest2[,uci:=b+1.96*se]
mvmr_forest2[,panel := paste0("Effect on ", tolower(clean_variable_name))]
mvmr_forest2 <- mvmr_forest2[order(colheader, hgnc, clean_variable_name)]
data <- data.table(mvmr_forest2)
p <- my_forestplotter_fancy(data = mvmr_forest2[outcome=="cad_aragam"], 
                            col.below.header = "hgnc", 
                            col.header = "colheader", 
                            col.header.heading = "Exposures", 
                            effect.name = "OR (95% CI)", 
                            # col.left = "clean_variable_name",
                            # col.left.toformat = "clean_variable_name",
                            # col.left.heading = "Outcomes",
                            col.left = "nsnp",
                            col.left.heading = "n SNPs",
                            col.right = NULL,
                            exponentiate = TRUE,
                            xlab = "Effect of 1-SD reduction of blood protein levels")

ggsave(filename = paste0("Results/Figure3_APOC3.png"),plot = p, dpi =300,
       width = 530/72,height = 230/72,units="in",scale=1, device = "png")  

p <- my_forestplotter_fancy(data = mvmr_forest2[outcome=="logTG_GLGC2022_European_MalesandFemales",], 
                            col.below.header = "hgnc", 
                            col.header = "colheader", 
                            col.header.heading = "Exposures", 
                            effect.name = "Effect (95% CI)", 
                            # col.left = "clean_variable_name",
                            # col.left.toformat = "clean_variable_name",
                            # col.left.heading = "Outcomes",
                            col.left = "nsnp",
                            col.left.heading = "n SNPs",
                            col.right = NULL,
                            exponentiate = FALSE,
                            xlab = "Effect of 1-SD reduction of blood protein levels")

ggsave(filename = paste0("Results/Figure4_APOC3.png"),plot = p, dpi =300,
       width = 530/72,height = 230/72,units="in",scale=1, device = "png")  

#######mediation analyses
res_mediation <- fread( "Data/Modified/res_mediation.txt")
k1<-res_mediation[method == "Multivariable IVW" & correctfor %in% c("trait-16-4", "trait-19-5") & id.exposure == "prot-5-3713", ]
k1[,id_exposure_outcome := paste0(id.exposure, "_", id.outcome)]
res[,id_exposure_outcome := paste0(id.exposure, "_", id.outcome)]
k2<-res[method == "Inverse variance weighted" & id_exposure_outcome %in% k1$id_exposure_outcome, ]
mediation_forest <- rbind(k1, k2, fill = TRUE)
mediation_forest[,b:=b*-1]
mediation_forest[,lci:=b-1.96*se]
mediation_forest[,uci:=b+1.96*se]
mediation_forest[,clean_variable_name := NULL]
mediation_forest <-  merge(mediation_forest, df_index[,.(id, clean_variable_name)], by.x = "correctfor", by.y = "id", all.x = TRUE)
setnames(mediation_forest, "clean_variable_name", "correct_for_clean")
mediation_forest <-  merge(mediation_forest, df_index[,.(id, clean_variable_name)], by.x = "id.outcome", by.y = "id", all.x = TRUE)
setnames(mediation_forest, "clean_variable_name", "outcome_name")
mediation_forest[is.na(correct_for_clean), correct_for_clean :="Unadjusted"]
mediation_forest[,panel := ""]
mediation_forest <- mediation_forest[order(outcome_name), ]
mediation_forest[,dummy:=""]
data <- data.table(mediation_forest)
p <- my_forestplotter_fancy(data = data, 
                            col.below.header = "dummy", 
                            col.header = "outcome_name", 
                            col.header.heading = "Outcomes", 
                            effect.name = "OR (95% CI)", 
                            col.left = "correct_for_clean",
                            col.left.toformat = NULL,
                            col.left.heading = "Adjusted for",
                            col.right = NULL,
                            exponentiate = TRUE,
                            xlab = "Effect of 1-SD reduction of APOC3")

ggsave(filename = paste0("Results/Figure5_APOC3.png"),plot = p, dpi =300,
       width = 472/72,height = 256/72,units="in",scale=1, device = "png")  

message("This script finished without errors")
