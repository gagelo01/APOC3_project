#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(DT)
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
res <-fread("Data/Modified/res_multicis_independent.txt")
res <- merge(res, dt_gene_region[,.(id,  hgnc, study)], by.x = "id.exposure", by.y = "id")
res <- merge(res, df_index[,.(id,  clean_variable_name)], by.x = "id.outcome", by.y = "id")
res[, mechanism := ifelse(study %in% c("ARIC", "FENLAND", "deCODE"),
                          "Blood protein levels",
                          "Liver gene levels")]
FandQ <- fread("Data/Modified/FandQ_univariate.txt")
egger_intercept <- fread( "Data/Modified/egger_intercept.txt")
egger_intercept <- merge(egger_intercept, dt_gene_region[,.(id,  hgnc, study)], by.x = "id.exposure", by.y = "id")
egger_intercept <- merge(egger_intercept, df_index[,.(id,  clean_variable_name)], by.x = "id.outcome", by.y = "id")
# dt<- rbindlist(list(
#   data.table(nom = c(paste0("trait-2-", 1:6), paste0("dis-6-", 1:3), "dis-23-2"), value = "Glucose homeostasis"),
#   data.table(nom = c(paste0("dis-13-", 1:4), "dis-5-1", df_index[grepl("dis-14-", id), id], paste0("trait-13-", 1:2), "dis-18-1", "dis-19-1", "dis-1-1"), value = "Vascular"),
#   data.table(nom = c("dis-2-2", "trait-14-8", paste0("trait-27-", 1:2)), value = "Liver"),
#   data.table(nom = c("dis-7-1", "trait-12-2"), value = "Kidney"),
#   data.table(nom = c(unique(df_index[grepl("trait-16-", id), id]), "trait-19-5","trait-19-2", "trait-29-20") , value = "Lipids"),
#   data.table(nom =   paste0("trait-28-", 1:2), value = "Bone"),
#   data.table(nom = c(unique(df_index[grepl("met-", id), id]), "trait-20-8"), value = "Metabolites"),
#   data.table(nom = c("trait-1-1", "trait-10-1", "trait-10-3", "trait-25-16", "trait-25-1"), value = "Anthropometric"),
#   data.table(nom = c("trait-6-1", "trait-7-1"), value = "Lifespan"),
#   data.table(nom = c("dis-12-1",  "trait-14-10"), value = "Pancreas")
# ))
# res[id.outcome %in% setdiff(unique(res$id.outcome), dt$nom), .(id.outcome, outcome, clean_variable_name)] %>% distinct(.)
# res <- merge(res, distinct(dt), by.x = "id.outcome", by.y = "nom")

######
return_format_data<-function(data) {
  k <- data[, paste0(format(round(exp(b), digits = 2), nsmall = 2), " 95% CI=", format(round(exp(lci), digits = 2), nsmall = 2), " to ",  format(round(exp(uci), digits = 2), nsmall = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))]
  names(k) <- data[,paste0(exposure," on ", outcome)]
  return(k)
}
return_format_data_noexp <-function(data) {
  k<-data[, paste0(format(round(b, digits = 2), nsmall = 2) , " 95% CI=", format(round(lci, digits = 2), nsmall = 2), " to ",  format(round(uci, digits = 2), nsmall = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))]
  names(k) <- data[,paste0(exposure," on ", outcome)]
  return(k)
}
return_format_fstat <-function(data) {
  return(data[, paste0(N, " SNPs (r2 = ", round(rsq*100, digits =2), "%; F-statistics = ",  round(fstat, digits = 0), ")")])
}
######
meth_tosel <- c("Weighted median", "Inverse variance weighted","MR-PRESSO (outlier-corrected)", "Contamination mixture")
k<-res[hgnc == "APOC3" & study == "deCODE" & method %in% meth_tosel, all(pval<0.05) & (all(b>0)|all(b<0)), by = c("id.outcome", "outcome", "study")]
idout_causal <- intersect(k[V1==TRUE,][order(id.outcome)]$id.outcome,
                          egger_intercept[hgnc == "APOC3" & study == "deCODE" & pval > 0.05,id.outcome ]) %>% unique

#####Abstract######
inst[hgnc == "APOC3", max(samplesize.exposure),by = "study"]
res[grepl("trait-", id.outcome), length(unique(id.outcome))]
res[grepl("dis-", id.outcome), length(unique(id.outcome))]
inst[hgnc == "APOC3", .N,by = "study"]
res_inhi <- data.table(res)
res_inhi[, b := b*-1]
res_inhi[, lci := b-se*1.96]
res_inhi[, uci := b+se*1.96]
k<-res_inhi[hgnc == "APOC3" & method == "Inverse variance weighted"  & id.outcome %in% idout_causal & study == "deCODE", ]
k[grepl("trait-", id.outcome)][order(b), ] %>% return_format_data_noexp
k[grepl("dis-", id.outcome)][order(b), ] %>% return_format_data
k[grepl("trait-7", id.outcome)][order(b), ] %>% return_format_data_noexp()

resmvmr <- fread( "Data/Modified/res_mvmr.txt")
resmvmr <- merge(resmvmr, distinct(res[,.(id.exposure, id.outcome, clean_variable_name, mechanism)]), by=c("id.exposure", "id.outcome"))
resmvmr_inhi <- data.table(resmvmr)
resmvmr_inhi[,b:=b*-1]
resmvmr_inhi[,lci:=b-1.96*se]
resmvmr_inhi[,uci:=b+1.96*se]

mvmr_forest <- rbind(res, resmvmr, fill = TRUE)
k<- mvmr_forest[method %in% c("Inverse variance weighted", "Multivariable IVW") & hgnc == "APOC3" & study == "deCODE", .(id.outcome, method, b)]
k_wide <- dcast(k, id.outcome ~ method, value.var = "b")
k_wide[, `Inverse variance weighted`- `Multivariable IVW`] %>% mean(.,na.rm = TRUE)
res_inhi[hgnc == "APOC3" & method == "Inverse variance weighted"  & pval > 0.05 & study == "deCODE" & grepl("dis-", id.outcome), ]$outcome %>% unique
#####Results######
#Para 1
maxsnp <- res_inhi[hgnc == "APOC3", max(nsnp)]
maxsnp
FandQ[method=="Inverse variance weighted" & id.exposure == "prot-5-3713" & N == maxsnp, rsq][1]
FandQ[method=="Inverse variance weighted" & id.exposure == "prot-5-3713" & N == maxsnp, ]
res_inhi[id.exposure=="prot-5-3713" & id.outcome == "trait-16-4" & method == "Inverse variance weighted", ] %>% return_format_data_noexp

#Para2
res_inhi[id.outcome %in% idout_causal & id.exposure == "prot-5-3713" & method == "Inverse variance weighted"
         & grepl("trait", id.outcome), ][order(b)] %>% return_format_data_noexp

res_inhi[!(id.outcome %in% idout_causal) & id.exposure == "prot-5-3713" & method == "Inverse variance weighted" & pval < 0.05 
         & grepl("trait", id.outcome), ][order(b)] %>% return_format_data_noexp


##para 3
res_inhi[id.outcome %in% idout_causal & id.exposure == "prot-5-3713" & method == "Inverse variance weighted"
         & grepl("dis", id.outcome), ][order(b)] %>% return_format_data

res_inhi[!(id.outcome %in% idout_causal) & id.exposure == "prot-5-3713" & method == "Inverse variance weighted" & pval < 0.05 
         & grepl("dis", id.outcome), ][order(b)] %>% return_format_data

res_inhi[!(id.outcome %in% idout_causal) & id.exposure == "prot-5-3713" & method == "Inverse variance weighted" & pval > 0.05 
         & grepl("dis", id.outcome), ][order(b)] %>% return_format_data

#para 4 
resmvmr_inhi[study=="deCODE" & method == "Multivariable IVW" & nsnp == max(nsnp), max(nsnp), by = "hgnc"]
resmvmr_inhi[study=="deCODE" & method == "Multivariable IVW" & nsnp == max(nsnp), max(F_stastistics), by = "hgnc"]
resmvmr_inhi[id.outcome%in%c("dis-13-1") & method == "Multivariable IVW" & id.exposure == "prot-5-3713", exp(b)]
resmvmr_inhi[id.outcome%in%c("trait-16-4") & method == "Multivariable IVW" & id.exposure == "prot-5-3713", b]
res_inhi[id.outcome%in%c("trait-16-4") & method == "Inverse variance weighted" & id.exposure == "prot-5-3713", b]

#Para 5
FandQ[method=="Inverse variance weighted" & id.exposure == "eqtl-74-3", ][N == max(N), rsq][1]
FandQ[method=="Inverse variance weighted" & id.exposure == "eqtl-74-3", ][N == max(N), fstat][1]

id_outcome_robust <- res_inhi[method %in% met_toinc & hgnc == "APOC3" & study == "GTEX", (all(b<0)|all(b>0)) & all(pval<0.05), by = "id.outcome"][V1==TRUE, id.outcome]
res_inhi[!(id.outcome %in% id_outcome_robust) & id.exposure == "prot-5-3713" & method == "Inverse variance weighted" & pval < 0.05 
         & grepl("dis", id.outcome), ][order(b)] %>% return_format_data



all_outcome[grepl("prot-", id.outcome), unique(samplesize.outcome), by = id.outcome]
inst_unicis$SNP %>% unique %>% length(.)
no_epitope <- setdiff(inst_unicis$SNP, df_VEP[LD_with!="", unique(SNP)])
epitope <- intersect(inst_unicis$SNP, df_VEP[LD_with!="", unique(SNP)])
length(epitope)
epitope
length(no_epitope)


