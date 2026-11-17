# Paths are configured through config/paths.sh or environment variables.
.mdd_repo <- Sys.getenv("MDD_CODE_ROOT", unset = getwd())
if (!file.exists(file.path(.mdd_repo, "config", "paths.R")))
  stop("Set MDD_CODE_ROOT to the repository root, or use scripts/run_r.sh.")
source(file.path(.mdd_repo, "config", "paths.R"))

library(purrr)
library(bigsnpr)
library(data.table)
library(tidyverse)

print( commandArgs(trailingOnly=TRUE))
## Arguments Loaded in
LD_ref = commandArgs(trailingOnly=TRUE)[1] # 1kg or ukbb
group = commandArgs(trailingOnly=TRUE)[2] # 1kg or ukbb

prs_cs_dir = mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_output_redo/")
prs_cs_output=  paste(prs_cs_dir, "ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_",LD_ref,"_LD_ref_",group,"_pst_eff_a1_b0.5_phiauto_allchr.txt", sep = "")
prs_cs = fread(prs_cs_output)
prs_cs_unique_id1 = paste(prs_cs$V1,prs_cs$V3, prs_cs$V4, prs_cs$V5, sep = "_")
prs_cs_unique_id2 = paste(prs_cs$V1,prs_cs$V3, prs_cs$V5, prs_cs$V4, sep = "_")

## Load bgen file
backingfile = mdd_path("MDD_BACKING", "")
obj.bigsnp <- snp_attach(backingfile)
obj_unique_id = paste(as.numeric(obj.bigsnp$map$chromosome), obj.bigsnp$map$physical.pos, obj.bigsnp$map$allele1, obj.bigsnp$map$allele2, sep = "_")
G <- obj.bigsnp$genotypes


## Obtain Indices and Effect Sizes for each heritability structure
ind.col.1 = match(prs_cs_unique_id1,obj_unique_id)
ind.col.2 = match(prs_cs_unique_id2,obj_unique_id)

sum(!is.na(ind.col.1))
sum(!is.na(ind.col.2))


if(!(all.equal(obj_unique_id[ind.col.1], prs_cs_unique_id1))){
  print("IDs not match!")
  afdasdffdsfa
}

print("start prediction")

# Calculate Xb for each heritability structure
pheno.prs_cs = big_prodVec(
  G,
  prs_cs$V6,
  ind.col = ind.col.1,
  center = NULL,
  scale = NULL,
  ncores = 1
)

out_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_"),LD_ref,"_LD_ref_",group,"_pst_eff_a1_b0.5_phiauto_prediction.txt", sep = "")
write.table(pheno.prs_cs, out_file)

print("done")
