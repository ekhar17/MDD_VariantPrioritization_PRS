# Paths are configured through config/paths.sh or environment variables.
.mdd_repo <- Sys.getenv("MDD_CODE_ROOT", unset = getwd())
if (!file.exists(file.path(.mdd_repo, "config", "paths.R")))
  stop("Set MDD_CODE_ROOT to the repository root, or use scripts/run_r.sh.")
source(file.path(.mdd_repo, "config", "paths.R"))

library(purrr)
library(bigsnpr)
library(data.table)
library(tidyverse)


## Load bgen file
backingfile = mdd_path("MDD_BACKING", "")
obj.bigsnp <- snp_attach(backingfile)
obj_unique_id = paste(as.numeric(obj.bigsnp$map$chromosome), obj.bigsnp$map$physical.pos, obj.bigsnp$map$allele1, obj.bigsnp$map$allele2, sep = "_")
G <- obj.bigsnp$genotypes


## Load sample ID
sample = fread(mdd_path("MDD_SAMPLE", ""))
head(sample,3)
## remove first individual
sample = sample[-1,]
print("remove placeholder")
head(sample,3)
## Create object ID
object_id = paste(sample$ID_1, sample$ID_2, sep = "_")
head(object_id)

## Load other ancestries
SAS_QC_all = fread(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/SAS_SinglePlink0.2_Het3SD_Sex_Exclusion_Full_Covar_QC_MDD.txt"))
SAS_qc_unique_id = paste(SAS_QC_all$FID, SAS_QC_all$IID, sep = "_")
AFR_QC_all = fread(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/AFR_SinglePlink0.2_Het3SD_Sex_Exclusion_Full_Covar_QC_MDD.txt"))
AFR_qc_unique_id = paste(AFR_QC_all$FID, AFR_QC_all$IID, sep = "_")
EUR_QC_all = fread(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/EUR_Single_Het_Sex_Exclusion_Full_Covar_QC_MDD.txt"))
EUR_qc_unique_id = paste(EUR_QC_all$FID, EUR_QC_all$IID, sep = "_")

AFR_MAF = snp_MAF(G, ind.row = which(object_id %in% AFR_qc_unique_id))
head(AFR_MAF)

SAS_MAF = snp_MAF(G, ind.row = which(object_id %in% SAS_qc_unique_id))
head(SAS_MAF)


map_info = obj.bigsnp$map
map_info$AFR_MAF = AFR_MAF
map_info$SAS_MAF = SAS_MAF
head(map_info)
fwrite(map_info, mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/map_info.csv"))

EUR_MAF = snp_MAF(G, ind.row = which(object_id %in% EUR_qc_unique_id))
head(EUR_MAF)

map_info$EUR_MAF = EUR_MAF
head(map_info)


fwrite(map_info, mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/map_info_wEUR.csv"))
