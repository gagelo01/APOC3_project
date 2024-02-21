setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/APOC3/")
source("Analysis/tosource.R")
res_multicis <- fread("Data/Modified/res_multicis_independent.txt")

res_multicis <- merge(res_multicis, dt_gene_region[,.(hgnc, id, study, gene_region)], by.x = "id.exposure", by.y ="id")
res_multicis <- merge(res_multicis, df_index[,.(id,clean_variable_name, note)], by.x = "id.outcome", by.y = "id")
res_multicis[, id.association := paste0(id.exposure,"_", id.outcome)]
method_tosel <- c("Inverse variance weighted", "Weighted median",
                  "MR-PRESSO (outlier-corrected)", "Contamination mixture")
tosel <- res_multicis[method %in% method_tosel, all(pval<0.05) & (all(b<0)|all(b>0)), by = "id.association"][V1==TRUE,]$id.association
res_multicis[id.association%in%tosel,.N, by = "hgnc" ]

res_multicis[id.association%in%tosel & method == "Inverse variance weighted",.N, by = c("hgnc", "outcome") ][order(hgnc)]

res_multicis[id.association%in%tosel & method == "Inverse variance weighted" & b < 0,]
res_multicis[hgnc == "APOC3" & method == "Inverse variance weighted", ][order(-b), .(study, outcome,nsnp, b, se, pval)]

