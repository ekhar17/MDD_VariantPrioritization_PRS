# Paths are configured through config/paths.sh or environment variables.
.mdd_repo <- Sys.getenv("MDD_CODE_ROOT", unset = getwd())
if (!file.exists(file.path(.mdd_repo, "config", "paths.R")))
  stop("Set MDD_CODE_ROOT to the repository root, or use scripts/run_r.sh.")
source(file.path(.mdd_repo, "config", "paths.R"))

## Load library
library(data.table)
library(bigsnpr)

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

## Make unique id that match hapmap id 1
summary_stat_hapmap_keep$hapmap_id_match = summary_stat_hapmap_keep$SNP_id1
summary_stat_hapmap_keep$hapmap_id_match[summary_stat_hapmap_keep$SNP_id2 %in% hapmap_id1] =  summary_stat_hapmap_keep$SNP_id2[summary_stat_hapmap_keep$SNP_id2 %in% hapmap_id1]
summary_stat_hapmap_logbf_over4$hapmap_id_match = summary_stat_hapmap_logbf_over4$SNP_id1
summary_stat_hapmap_logbf_over4$hapmap_id_match[summary_stat_hapmap_logbf_over4$SNP_id2 %in% hapmap_id1] =  summary_stat_hapmap_logbf_over4$SNP_id2[summary_stat_hapmap_logbf_over4$SNP_id2 %in% hapmap_id1]
summary_stat_hapmap_cr95set$hapmap_id_match = summary_stat_hapmap_cr95set$SNP_id1
summary_stat_hapmap_cr95set$hapmap_id_match[summary_stat_hapmap_cr95set$SNP_id2 %in% hapmap_id1] =  summary_stat_hapmap_cr95set$SNP_id2[summary_stat_hapmap_cr95set$SNP_id2 %in% hapmap_id1]
summary_stat_hapmap_twassig$hapmap_id_match = summary_stat_hapmap_twassig$SNP_id1
summary_stat_hapmap_twassig$hapmap_id_match[summary_stat_hapmap_twassig$SNP_id2 %in% hapmap_id1] =  summary_stat_hapmap_twassig$SNP_id2[summary_stat_hapmap_twassig$SNP_id2 %in% hapmap_id1]

## Order hapmap to match summary stat
hapmap_match = hapmap[match(summary_stat_hapmap_keep$hapmap_id_match, hapmap_id1),]
hapmap_match_logbf_over4 = hapmap[match(summary_stat_hapmap_logbf_over4$hapmap_id_match, hapmap_id1),]
hapmap_match_cr95set = hapmap[match(summary_stat_hapmap_cr95set$hapmap_id_match, hapmap_id1),]
hapmap_match_twassig = hapmap[match(summary_stat_hapmap_twassig$hapmap_id_match, hapmap_id1),]

## Make SNP ID match
summary_stat_hapmap_keep$SNP[!summary_stat_hapmap_keep$SNP %in% hapmap_match$SNP] = hapmap_match$SNP[!summary_stat_hapmap_keep$SNP %in% hapmap_match$SNP]
summary_stat_hapmap_logbf_over4$SNP[!summary_stat_hapmap_logbf_over4$SNP %in% hapmap_match_logbf_over4$SNP] = hapmap_match_logbf_over4$SNP[!summary_stat_hapmap_logbf_over4$SNP %in% hapmap_match_logbf_over4$SNP]
summary_stat_hapmap_cr95set$SNP[!summary_stat_hapmap_cr95set$SNP %in% hapmap_match_cr95set$SNP] = hapmap_match_cr95set$SNP[!summary_stat_hapmap_cr95set$SNP %in% hapmap_match_cr95set$SNP]
summary_stat_hapmap_twassig$SNP[!summary_stat_hapmap_twassig$SNP %in% hapmap_match_twassig$SNP] = hapmap_match_twassig$SNP[!summary_stat_hapmap_twassig$SNP %in% hapmap_match_twassig$SNP]

## Keep only columns needed for PRS CS
prs_cs_sumstat = data.frame(SNP = summary_stat_hapmap_keep$SNP,
                            A1 = summary_stat_hapmap_keep$A1,
                            A2 = summary_stat_hapmap_keep$A2,
                            OR = summary_stat_hapmap_keep$OR,
                            SE = summary_stat_hapmap_keep$SE)
prs_cs_sumstat_logbf_over4 = data.frame(SNP = summary_stat_hapmap_logbf_over4$SNP,
                            A1 = summary_stat_hapmap_logbf_over4$A1,
                            A2 = summary_stat_hapmap_logbf_over4$A2,
                            OR = summary_stat_hapmap_logbf_over4$OR,
                            SE = summary_stat_hapmap_logbf_over4$SE)
prs_cs_sumstat_cr95set = data.frame(SNP = summary_stat_hapmap_cr95set$SNP,
                                        A1 = summary_stat_hapmap_cr95set$A1,
                                        A2 = summary_stat_hapmap_cr95set$A2,
                                        OR = summary_stat_hapmap_cr95set$OR,
                                        SE = summary_stat_hapmap_cr95set$SE)
prs_cs_sumstat_twassig = data.frame(SNP = summary_stat_hapmap_twassig$SNP,
                                    A1 = summary_stat_hapmap_twassig$A1,
                                    A2 = summary_stat_hapmap_twassig$A2,
                                    OR = summary_stat_hapmap_twassig$OR,
                                    SE = summary_stat_hapmap_twassig$SE)

## Calculate N_gwas
N_gwas = max(summary_stat_hapmap_keep$Nca + summary_stat_hapmap_keep$Nco , na.rm = T)
max(summary_stat_hapmap_logbf_over4$Nca + summary_stat_hapmap_logbf_over4$Nco, na.rm = T)
max(summary_stat_hapmap_cr95set$Nca + summary_stat_hapmap_cr95set$Nco, na.rm = T)
max(summary_stat_hapmap_twassig$Nca + summary_stat_hapmap_twassig$Nco, na.rm = T)
## All are 674,452

## Create Bim for UK Biobank based on SNPs in
## Load Backing File
backingfile = mdd_path("MDD_BACKING", "")
obj.bigsnp <- snp_attach(backingfile)
str(obj.bigsnp)

## Obtain SNP ID's and geno information
obj.bigsnp_id = paste(as.numeric(obj.bigsnp$map$chromosome), obj.bigsnp$map$physical.pos,
                      obj.bigsnp$map$allele1, obj.bigsnp$map$allele2, sep = "_")

## Obj.bigsnp_id_hapmap and add unique id
UKB_in_hapmap = as.data.frame(obj.bigsnp$map[obj.bigsnp_id %in% hapmap_id1 |obj.bigsnp_id %in% hapmap_id2,])
UKB_in_hapmap$SNP_id1 = paste(as.numeric(UKB_in_hapmap$chromosome), UKB_in_hapmap$physical.pos,
                              UKB_in_hapmap$allele1, UKB_in_hapmap$allele2, sep = "_")
UKB_in_hapmap$SNP_id2 = paste(as.numeric(UKB_in_hapmap$chromosome), UKB_in_hapmap$physical.pos,
                              UKB_in_hapmap$allele2, UKB_in_hapmap$allele1, sep = "_")
## Make a hapmap id vector
UKB_in_hapmap$hapmap_id_match = UKB_in_hapmap$SNP_id1
UKB_in_hapmap$hapmap_id_match[UKB_in_hapmap$SNP_id2 %in% hapmap_id1] = UKB_in_hapmap$SNP_id2[UKB_in_hapmap$SNP_id2 %in% hapmap_id1]

## Match hapmap
hapmap_ukb_match = hapmap[match(UKB_in_hapmap$hapmap_id_match, hapmap_id1),]

## Replace RSID with those that match
UKB_in_hapmap$rsid[!UKB_in_hapmap$rsid %in% hapmap_ukb_match$SNP] = hapmap_ukb_match$SNP[!hapmap_ukb_match$SNP %in% UKB_in_hapmap$rsid]

## Make Bim Out
bim = data.frame(CHR = as.numeric(UKB_in_hapmap$chromosome),
                 ID = UKB_in_hapmap$rsid,
                 POS = 0,
                 BP = UKB_in_hapmap$physical.pos,
                 A1 = UKB_in_hapmap$allele1,
                 A2 = UKB_in_hapmap$allele2)

## Save Reformatted GWAS, BIM, and N_gwas
fwrite(prs_cs_sumstat,
       mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap.txt"),
       sep = "\t", quote = F)
fwrite(prs_cs_sumstat_logbf_over4,
       mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_logbf_over4.txt"),
       sep = "\t", quote = F)
fwrite(prs_cs_sumstat_cr95set,
       mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_cr95set.txt"),
       sep = "\t", quote = F)
fwrite(prs_cs_sumstat_twassig,
       mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_twassig.txt"),
       sep = "\t", quote = F)
fwrite(bim,
       mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/bim/UKB_hapmap_bim.bim"),
       sep = '\t', quote = F, col.names = F)
write.table(N_gwas,
            mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/ispych_MDD_GWAS_wo_UKB23ndme_NGWAS.txt"),
            col.names = F, row.names = F, quote = F)
