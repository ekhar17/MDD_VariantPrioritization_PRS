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
trait = commandArgs(trailingOnly=TRUE)[2] # BIP, SCZ, ADHD
group = commandArgs(trailingOnly=TRUE)[3] # logbfover4, cr95set

prs_cs_dir = mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_output_redo/")
prs_cs_output= paste(prs_cs_dir,"/MDD_GWAS_",LD_ref,"_LD_ref_",trait,"_",group,"_set_pst_eff_a1_b0.5_phiauto_allchr_redo.txt", sep = "")
prs_cs_output
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

out_file = paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_prediction_redo/MDD_GWAS_"),LD_ref,"_LD_ref_",trait,"_",group,"_set_pst_eff_a1_b0.5_phiauto_prediction_updated.txt", sep = "")
write.table(pheno.prs_cs, out_file)

print(out_file)
print("done")
