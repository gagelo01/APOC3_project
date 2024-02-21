#!/usr/bin/env Rscript
setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/APOC3/")
source("Analysis/tosource.R")
inst <- fread("Data/Modified/inst.txt")
inst <- inst[typeof_sel == "multicis_independent_clumping", ]
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")
inst <- merge(inst, dt_gene_region[,.(id,  hgnc, study)], by.x = "id.exposure", by.y = "id")
all_outcome <- fread("Data/Modified/all_outcome.txt")


#What is the best way to design a multivariable MR ?
#I suggest taking the full set of SNPs, clumping it at 0.1. Then doing every combination of 
#study. so deCODE-deCODE, Fenland-deCODE, and Fenland for every outcome studied.
ID <- c(inst[grepl("prot-5-|prot-7-", id.exposure), unique(id.exposure)])
all_inst_mvmr <- all_outcome[id.outcome %in% ID, ] %>% GagnonMR::convert_outcome_to_exposure(.) %>% as.data.table(.)

k <- dt_gene_region[id %in% ID, unique(id), by = "study"]
id_exposure_vec <- lapply(split(k, k$study), function(x) x$V1)
id_outcome_vec <- ID_server_out

arguments <- tidyr::expand_grid(id_exposure_vec, id_outcome_vec)

res_mvmr <- future_map(split(arguments, 1:nrow(arguments)), function(x) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  resmvmmr <- GagnonMR::performmvmr(id_exposure_vec = x$id_exposure_vec[[1]], 
                        id_outcome_vec = x$id_outcome_vec, 
                        inst_all_sign_clump = inst, 
                        all_inst_mvmr = all_inst_mvmr, 
                        all_outcome_mvmr = all_outcome,
                        clump_r2 = 0.1,
                        clump_kb = 1000,
                        parameters = parameters)
  return(resmvmmr)
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)

res_mvmr <- merge(res_mvmr, dt_gene_region[,.(id,  hgnc, study)], by.x = "id.exposure", by.y = "id")

fwrite(res_mvmr, "Data/Modified/res_mvmr.txt")

###Test MVMR APOC3 -> TG/APOB -> CAD/AP
k2 <- all_outcome[id.outcome %in% c("trait-16-4", "trait-19-5"), ] %>%  
  convert_outcome_to_exposure() %>% 
  as.data.table(.)
inst_mediation <- rbind(inst,k2, fill = TRUE)

adjustfor <- list(c("trait-16-4", "trait-19-5"), "trait-16-4", "trait-19-5")
arg1 <- tidyr::expand_grid(id_exposure =  dt_gene_region[id %in% ID & hgnc == "APOC3", unique(id)], adjustfor)
arguments <- tidyr::expand_grid(arg1, id_outcome= c("dis-13-1", "dis-12-1", "dis-19-1"))

res_mediation <- future_map(split(arguments, 1:nrow(arguments)), function(x) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  resmvmmr <- GagnonMR::performmvmr(id_exposure_vec = x$id_exposure, 
                          id_outcome_vec = x$id_outcome, 
                          correctfor = x$adjustfor[[1]],
                          inst_all_sign_clump = inst_mediation, 
                          all_inst_mvmr = inst_mediation, 
                          all_outcome_mvmr = all_outcome,
                          should_clump = FALSE,
                          clump_r2 = 0.1,
                          clump_kb = 1000,
                          parameters = parameters)
  return(resmvmmr)
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)

fwrite(res_mediation, "Data/Modified/res_mediation.txt")

message("This script finished without errors")
