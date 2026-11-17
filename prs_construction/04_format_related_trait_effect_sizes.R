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


dir = mdd_path("MDD_SOURCE_DATA", "/derived_data/procedure_data/fgannosumstats/")
end = "_GWAS_sumstats.QC.maf01info09_95crlgbf4.txt"

disorders = c("ADHD", "BIP", "SCZ")

for(dis in disorders){

  ## Load summary stat
  file = paste(dir, dis, end, sep = "")
  alttrait_summary_stat = fread(file)
  print("loaded sum stats")


  ## Add Unique ID
  alttrait_summary_stat$SNP_id1 = paste(alttrait_summary_stat$CHR, alttrait_summary_stat$POShg19, alttrait_summary_stat$A1, alttrait_summary_stat$A2, sep = "_")


  ## Keep only those in hapmap
  alttrait_summary_stat_hapmap = alttrait_summary_stat[alttrait_summary_stat$SNP_id1 %in% hapmap_id1 | alttrait_summary_stat$SNP_id1 %in% hapmap_id2,]
  alttrait_summary_stat_hapmap$SNP_id2 = paste(alttrait_summary_stat_hapmap$CHR, alttrait_summary_stat_hapmap$POShg19, alttrait_summary_stat_hapmap$A2, alttrait_summary_stat_hapmap$A1, sep = "_")
  alttrait_summary_stat_hapmap$hapmap_id_match = alttrait_summary_stat_hapmap$SNP_id1
  alttrait_summary_stat_hapmap$hapmap_id_match[alttrait_summary_stat_hapmap$SNP_id2 %in% hapmap_id1] = alttrait_summary_stat_hapmap$SNP_id2[alttrait_summary_stat_hapmap$SNP_id2 %in% hapmap_id1]

  alttrait_MDD_logbfover4 = alttrait_summary_stat_hapmap[alttrait_summary_stat_hapmap$SNP_id1 %in% summary_stat_hapmap_logbf_over4$SNP_id1 | alttrait_summary_stat_hapmap$SNP_id2 %in% summary_stat_hapmap_logbf_over4$SNP_id1, ]
alttrait_MDD_cr95set =  alttrait_summary_stat_hapmap[alttrait_summary_stat_hapmap$SNP_id1 %in% summary_stat_hapmap_cr95set$SNP_id1 | alttrait_summary_stat_hapmap$SNP_id2 %in% summary_stat_hapmap_cr95set$SNP_id1, ]
alttrait_MDD_twassig =  alttrait_summary_stat_hapmap[alttrait_summary_stat_hapmap$SNP_id1 %in% summary_stat_hapmap_twassig$SNP_id1 | alttrait_summary_stat_hapmap$SNP_id2 %in% summary_stat_hapmap_twassig$SNP_id1, ]

print(all.equal(nrow(alttrait_MDD_logbfover4), sum(alttrait_summary_stat_hapmap$RSID %in% summary_stat_hapmap_logbf_over4$SNP)))
print(all.equal(nrow(alttrait_MDD_cr95set), sum(alttrait_summary_stat_hapmap$RSID %in% summary_stat_hapmap_cr95set$SNP)))
print(all.equal(nrow(alttrait_MDD_twassig), sum(alttrait_summary_stat_hapmap$RSID %in% summary_stat_hapmap_twassig$SNP)))

## Order hapmap to match summary stat
hapmap_match_alttrait_logbfover4 = hapmap[match(alttrait_MDD_logbfover4$hapmap_id_match, hapmap_id1),]
hapmap_match_alttrait_cr95set = hapmap[match(alttrait_MDD_cr95set$hapmap_id_match, hapmap_id1),]
hapmap_match_alttrait_twassig = hapmap[match(alttrait_MDD_twassig$hapmap_id_match, hapmap_id1),]

## Make SNP ID match
alttrait_MDD_logbfover4$SNP = hapmap_match_alttrait_logbfover4$SNP
alttrait_MDD_cr95set$SNP = hapmap_match_alttrait_cr95set$SNP
alttrait_MDD_twassig$SNP = hapmap_match_alttrait_twassig$SNP

print(all.equal(alttrait_MDD_logbfover4$RSID, alttrait_MDD_logbfover4$SNP))
      print(all.equal(alttrait_MDD_cr95set$RSID, alttrait_MDD_cr95set$SNP))
            print(all.equal(alttrait_MDD_twassig$RSID, alttrait_MDD_twassig$SNP))

prs_cs_sumstat_alttrait_MDD_logbfover4 = data.frame(SNP = alttrait_MDD_logbfover4$SNP,
                                                    A1 = alttrait_MDD_logbfover4$A1,
                                                    A2 = alttrait_MDD_logbfover4$A2,
                                                    BETA = alttrait_MDD_logbfover4$B,
                                                    SE = alttrait_MDD_logbfover4$SE)

prs_cs_sumstat_alttrait_MDD_cr95set = data.frame(SNP = alttrait_MDD_cr95set$SNP,
                                                 A1 = alttrait_MDD_cr95set$A1,
                                                 A2 = alttrait_MDD_cr95set$A2,
                                                 BETA = alttrait_MDD_cr95set$B,
                                                 SE = alttrait_MDD_cr95set$SE)

prs_cs_sumstat_alttrait_MDD_twassig = data.frame(SNP = alttrait_MDD_twassig$SNP,
                                                 A1 = alttrait_MDD_twassig$A1,
                                                 A2 = alttrait_MDD_twassig$A2,
                                                 BETA = alttrait_MDD_twassig$B,
                                                 SE = alttrait_MDD_twassig$SE)

## Save Reformatted GWAS, BIM, and N_gwas
fwrite(prs_cs_sumstat_alttrait_MDD_logbfover4,
       paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/hapmap_MDD_Set_logbfover4_"),dis,"_coef.txt",sep = ""),
       sep = "\t", quote = F)
fwrite(prs_cs_sumstat_alttrait_MDD_cr95set,
       paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/hapmap_MDD_Set_cr95set_"),dis,"_coef.txt",sep = ""),
       sep = "\t", quote = F)
fwrite(prs_cs_sumstat_alttrait_MDD_twassig,
       paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/hapmap_MDD_Set_twassig_"),dis,"_coef.txt",sep = ""),
       sep = "\t", quote = F)

print("saved")
}


###########################################################                  ANX Meta Subtract                  ################################################################################


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


## Load ANX sumstats
ANX_sumstats = fread(mdd_path("MDD_SOURCE_DATA", "/derived_data/procedure_data/fgannosumstats/ANX_EUR_ukbsubstra_GWAS_sumstats.QC.maf01info09_95crlgbf4.txt"))
natgen = fread(mdd_path("MDD_SOURCE_DATA", "/original_data/GWASsumstats/ANX_EUR.txt.gz"))

## Add hg19 position to ANX_sumstats
natgen_match = natgen[match(ANX_sumstats$RSID, natgen$SNP),]
all.equal(natgen_match$SNP, ANX_sumstats$RSID)
ANX_sumstats$POShg19 = natgen_match$BP
ANX_sumstats$A1_hg19 = natgen_match$A1
ANX_sumstats$A2_hg19 = natgen_match$A2


## Add unique ID
ANX_sumstats$SNP_id1 = paste(ANX_sumstats$CHR, ANX_sumstats$POShg19, ANX_sumstats$A1_hg19, ANX_sumstats$A2_hg19, sep = "_")


## Keep only those in hapmap
ANX_sumstats_hapmap = ANX_sumstats[ANX_sumstats$SNP_id1 %in% hapmap_id1 | ANX_sumstats$SNP_id1 %in% hapmap_id2,]
ANX_sumstats_hapmap$SNP_id2 = paste(ANX_sumstats_hapmap$CHR, ANX_sumstats_hapmap$POShg19, ANX_sumstats_hapmap$A2_hg19, ANX_sumstats_hapmap$A1_hg19, sep = "_")

alttrait_MDD_logbfover4 = ANX_sumstats_hapmap[ANX_sumstats_hapmap$SNP_id1 %in% summary_stat_hapmap_logbf_over4$SNP_id1 | ANX_sumstats_hapmap$SNP_id2 %in% summary_stat_hapmap_logbf_over4$SNP_id1, ]
alttrait_MDD_cr95set = ANX_sumstats_hapmap[ANX_sumstats_hapmap$SNP_id1 %in% summary_stat_hapmap_cr95set$SNP_id1 | ANX_sumstats_hapmap$SNP_id2 %in% summary_stat_hapmap_cr95set$SNP_id1, ]
alttrait_MDD_twassig = ANX_sumstats_hapmap[ANX_sumstats_hapmap$SNP_id1 %in% summary_stat_hapmap_twassig$SNP_id1 | ANX_sumstats_hapmap$SNP_id2 %in% summary_stat_hapmap_twassig$SNP_id1, ]


  ## Make unique id that match hapmap id 1
alttrait_MDD_logbfover4$hapmap_id_match = alttrait_MDD_logbfover4$SNP_id1
alttrait_MDD_logbfover4$hapmap_id_match[alttrait_MDD_logbfover4$SNP_id2 %in% hapmap_id1] = alttrait_MDD_logbfover4$SNP_id2[alttrait_MDD_logbfover4$SNP_id2 %in% hapmap_id1]

alttrait_MDD_cr95set$hapmap_id_match = alttrait_MDD_cr95set$SNP_id1
alttrait_MDD_cr95set$hapmap_id_match[alttrait_MDD_cr95set$SNP_id2 %in% hapmap_id1] = alttrait_MDD_cr95set$SNP_id2[alttrait_MDD_cr95set$SNP_id2 %in% hapmap_id1]

alttrait_MDD_twassig$hapmap_id_match = alttrait_MDD_twassig$SNP_id1
alttrait_MDD_twassig$hapmap_id_match[alttrait_MDD_twassig$SNP_id2 %in% hapmap_id1] = alttrait_MDD_twassig$SNP_id2[alttrait_MDD_twassig$SNP_id2 %in% hapmap_id1]


  ## Order hapmap to match summary stat
  hapmap_match_alttrait_logbfover4 = hapmap[match(alttrait_MDD_logbfover4$hapmap_id_match, hapmap_id1),]
  hapmap_match_alttrait_cr95set = hapmap[match(alttrait_MDD_cr95set$hapmap_id_match, hapmap_id1),]
  hapmap_match_alttrait_twassig = hapmap[match(alttrait_MDD_twassig$hapmap_id_match, hapmap_id1),]

  ## Make SNP ID match
  alttrait_MDD_logbfover4$SNP = hapmap_match_alttrait_logbfover4$SNP
  alttrait_MDD_cr95set$SNP = hapmap_match_alttrait_cr95set$SNP
  alttrait_MDD_twassig$SNP = hapmap_match_alttrait_twassig$SNP

  all.equal(alttrait_MDD_logbfover4$RSID, alttrait_MDD_logbfover4$SNP)
  all.equal(alttrait_MDD_cr95set$RSID, alttrait_MDD_cr95set$SNP)
  all.equal(alttrait_MDD_twassig$RSID, alttrait_MDD_twassig$SNP)

  prs_cs_sumstat_alttrait_MDD_logbfover4 = data.frame(SNP = alttrait_MDD_logbfover4$SNP,
                                                  A1 = alttrait_MDD_logbfover4$A1_hg19,
                                                  A2 = alttrait_MDD_logbfover4$A2_hg19,
                                                  BETA = alttrait_MDD_logbfover4$Z,
                                                  P = alttrait_MDD_logbfover4$P)

  prs_cs_sumstat_alttrait_MDD_cr95set = data.frame(SNP = alttrait_MDD_cr95set$SNP,
                                                      A1 = alttrait_MDD_cr95set$A1_hg19,
                                                      A2 = alttrait_MDD_cr95set$A2_hg19,
                                                      BETA = alttrait_MDD_cr95set$Z,
                                                      P = alttrait_MDD_cr95set$P)

  prs_cs_sumstat_alttrait_MDD_twassig = data.frame(SNP = alttrait_MDD_twassig$SNP,
                                                      A1 = alttrait_MDD_twassig$A1_hg19,
                                                      A2 = alttrait_MDD_twassig$A2_hg19,
                                                      BETA = alttrait_MDD_twassig$Z,
                                                      P = alttrait_MDD_twassig$P)

 ## Save Reformatted GWAS, BIM, and N_gwas
  fwrite(prs_cs_sumstat_alttrait_MDD_logbfover4,
         paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/hapmap_MDD_Set_logbfover4_ANXmetasubtract_coef.txt"),sep = ""),
         sep = "\t", quote = F)
  fwrite(prs_cs_sumstat_alttrait_MDD_cr95set,
         paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/hapmap_MDD_Set_cr95set_ANXmetasubtract_coef.txt"),sep = ""),
         sep = "\t", quote = F)
  fwrite(prs_cs_sumstat_alttrait_MDD_twassig,
         paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/hapmap_MDD_Set_twassig_ANXmetasubtract_coef.txt"),sep = ""),
         sep = "\t", quote = F)

  print("saved")
