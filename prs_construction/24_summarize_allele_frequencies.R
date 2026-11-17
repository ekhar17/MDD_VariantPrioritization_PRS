# Paths are configured through config/paths.sh or environment variables.
.mdd_repo <- Sys.getenv("MDD_CODE_ROOT", unset = getwd())
if (!file.exists(file.path(.mdd_repo, "config", "paths.R")))
  stop("Set MDD_CODE_ROOT to the repository root, or use scripts/run_r.sh.")
source(file.path(.mdd_repo, "config", "paths.R"))


## Load library
library(data.table)
library(bigsnpr)

map_info  = fread(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/map_info_wEUR.csv"))
map_info


## Load Summary Statistics
summary_stat = fread(mdd_path("MDD_SOURCE_DATA", "/original_data/ispych_MDD_GWAS_wo_UKB23ndme/daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_Wray_FinnGen_MVPaf_2_HRC_MAF01.gz"))

## Load SNPs used for other PRS
mdd_snp_info = fread(mdd_path("MDD_SOURCE_DATA", "/derived_data/procedure_data/ipsychMDDGWAS.sumstats.QC.maf01info09_5grp.txt"))

## Load hapmap list
hapmap = fread(mdd_path("MDD_HAPMAP", ""))
hapmap_id1 = paste(hapmap$CHR,hapmap$BP, hapmap$A1, hapmap$A2, sep = "_")
hapmap_id2 = paste(hapmap$CHR,hapmap$BP,  hapmap$A2, hapmap$A1, sep = "_")

## Create summary stat and mdd snp info unique id
summary_stat$SNP_id1 = paste(summary_stat$CHR, summary_stat$BP, summary_stat$A1, summary_stat$A2, sep = "_")
mdd_snp_info$SNP_id1 = paste(mdd_snp_info$CHR, mdd_snp_info$POShg19, mdd_snp_info$A1, mdd_snp_info$A2, sep = "_")


## Keep only those in hapmap
summary_stat_hapmap = summary_stat[summary_stat$SNP_id1 %in% hapmap_id1 | summary_stat$SNP_id1 %in% hapmap_id2,]
mdd_snp_info_hapmap = mdd_snp_info[mdd_snp_info$SNP_id1 %in% hapmap_id1 |mdd_snp_info$SNP_id1 %in% hapmap_id2 ,]

## Add SNP id 2 in case flipped in UK B
mdd_snp_info_hapmap$SNP_id2 = paste(mdd_snp_info_hapmap$CHR, mdd_snp_info_hapmap$POShg19, mdd_snp_info_hapmap$A2, mdd_snp_info_hapmap$A1, sep = "_")

## Keep only summary stat in mdd gwas info
summary_stat_hapmap_keep = summary_stat_hapmap[summary_stat_hapmap$SNP_id1 %in% mdd_snp_info_hapmap$SNP_id1 | summary_stat_hapmap$SNP_id1 %in% mdd_snp_info_hapmap$SNP_id2,]

## Repeat for 3 groups
summary_stat_hapmap_logbf_over4 = summary_stat_hapmap[summary_stat_hapmap$SNP_id1 %in% mdd_snp_info_hapmap$SNP_id1[as.logical(mdd_snp_info_hapmap$logbf_over4)]
                                                      | summary_stat_hapmap$SNP_id1 %in% mdd_snp_info_hapmap$SNP_id2[as.logical(mdd_snp_info_hapmap$logbf_over4)],]

summary_stat_hapmap_cr95set = summary_stat_hapmap[summary_stat_hapmap$SNP_id1 %in% mdd_snp_info_hapmap$SNP_id1[as.logical(mdd_snp_info_hapmap$cr95set)]
                                                  | summary_stat_hapmap$SNP_id1 %in% mdd_snp_info_hapmap$SNP_id2[as.logical(mdd_snp_info_hapmap$cr95set)],]

summary_stat_hapmap_twassig = summary_stat_hapmap[summary_stat_hapmap$SNP_id1 %in% mdd_snp_info_hapmap$SNP_id1[as.logical(mdd_snp_info_hapmap$twassig)]
                                                  | summary_stat_hapmap$SNP_id1 %in% mdd_snp_info_hapmap$SNP_id2[as.logical(mdd_snp_info_hapmap$twassig)],]


## Replace rest of summary stat keep rsid with name in hapmap file
## Add SNPID2
summary_stat_hapmap_keep$SNP_id2 = paste(summary_stat_hapmap_keep$CHR, summary_stat_hapmap_keep$BP, summary_stat_hapmap_keep$A2, summary_stat_hapmap_keep$A1, sep = "_")
summary_stat_hapmap_logbf_over4$SNP_id2 = paste(summary_stat_hapmap_logbf_over4$CHR, summary_stat_hapmap_logbf_over4$BP, summary_stat_hapmap_logbf_over4$A2, summary_stat_hapmap_logbf_over4$A1, sep = "_")
summary_stat_hapmap_cr95set$SNP_id2 = paste(summary_stat_hapmap_cr95set$CHR, summary_stat_hapmap_cr95set$BP, summary_stat_hapmap_cr95set$A2, summary_stat_hapmap_cr95set$A1, sep = "_")
summary_stat_hapmap_twassig$SNP_id2 = paste(summary_stat_hapmap_twassig$CHR, summary_stat_hapmap_twassig$BP, summary_stat_hapmap_twassig$A2, summary_stat_hapmap_twassig$A1, sep = "_")

map_info$SNP_id = paste(as.numeric(map_info$chromosome),
                        map_info$physical.pos,
                        map_info$allele1,
                        map_info$allele2,
                        sep = "_")

map_info_keep = map_info[map_info$SNP_id %in% summary_stat_hapmap_keep$SNP_id1 |
                           map_info$SNP_id %in% summary_stat_hapmap_keep$SNP_id2,]
map_info_logbf_over4_keep = map_info[map_info$SNP_id %in% summary_stat_hapmap_logbf_over4$SNP_id1 |
                           map_info$SNP_id %in% summary_stat_hapmap_logbf_over4$SNP_id2,]
map_info_cr95set = map_info[map_info$SNP_id %in% summary_stat_hapmap_cr95set$SNP_id1 |
                                       map_info$SNP_id %in% summary_stat_hapmap_cr95set$SNP_id2,]
map_info_twassig = map_info[map_info$SNP_id %in% summary_stat_hapmap_twassig$SNP_id1 |
                                    map_info$SNP_id %in% summary_stat_hapmap_twassig$SNP_id2,]


hist(map_info_keep$EUR_MAF)
hist(map_info_keep$AFR_MAF)
hist(map_info_keep$SAS_MAF)
sum(map_info_keep$EUR_MAF < 1e-3)
sum(map_info_keep$AFR_MAF < 1e-3)
sum(map_info_keep$SAS_MAF < 1e-3)


sum(map_info_logbf_over4_keep$EUR_MAF < 1e-3)
sum(map_info_logbf_over4_keep$AFR_MAF < 1e-3)
sum(map_info_logbf_over4_keep$SAS_MAF < 1e-3)


sum(map_info_cr95set$EUR_MAF < 1e-3)
sum(map_info_cr95set$AFR_MAF < 1e-3)
sum(map_info_cr95set$SAS_MAF < 1e-3)

sum(map_info_twassig$EUR_MAF < 1e-3)
sum(map_info_twassig$AFR_MAF < 1e-3)
sum(map_info_twassig$SAS_MAF < 1e-3)


QC_all = fread(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/EUR_Single_Het_Sex_Exclusion_Full_Covar_QC_MDD.txt"))
SAS_QC_all = fread(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/SAS_SinglePlink0.2_Het3SD_Sex_Exclusion_Full_Covar_QC_MDD.txt"))
AFR_QC_all = fread(mdd_path("MDD_PROJECT", "/01.TrainTestSplit/data/AFR_SinglePlink0.2_Het3SD_Sex_Exclusion_Full_Covar_QC_MDD.txt"))

table(QC_all$MDD)
table(SAS_QC_all$MDD)
table(AFR_QC_all$MDD)


table(QC_all$MDD)/nrow(QC_all)
table(SAS_QC_all$MDD)/nrow(SAS_QC_all)
table(AFR_QC_all$MDD)/nrow(AFR_QC_all)
