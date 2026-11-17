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


## Load hapmap list
hapmap = fread(mdd_path("MDD_HAPMAP", ""))
hapmap_id1 = paste(hapmap$CHR,hapmap$BP, hapmap$A1, hapmap$A2, sep = "_")
hapmap_id2 = paste(hapmap$CHR,hapmap$BP,  hapmap$A2, hapmap$A1, sep = "_")

## Create summary stat and mdd snp info unique id
summary_stat$SNP_id1 = paste(summary_stat$CHR, summary_stat$BP, summary_stat$A1, summary_stat$A2, sep = "_")
## Keep only those in hapmap
summary_stat_hapmap = summary_stat[summary_stat$SNP_id1 %in% hapmap_id1 | summary_stat$SNP_id1 %in% hapmap_id2,]


dir = mdd_path("MDD_SOURCE_DATA", "/derived_data/procedure_data/fgannosumstats/")
end = "_GWAS_sumstats.QC.maf01info09_95crlgbf4.txt"


disorders = c("ADHD", "BIP", "SCZ")

for(dis in disorders){

  ## Load summary stat
  file = paste(dir, dis, end, sep = "")
  alttrait_snpinfo = fread(file)
  print("loaded alt trait sum stats")

  ## Add Unique ID
  alttrait_snpinfo$SNP_id1 = paste(alttrait_snpinfo$CHR, alttrait_snpinfo$POShg19, alttrait_snpinfo$A1, alttrait_snpinfo$A2, sep = "_")

  ## Keep only those in hapmap
  alttrait_snpinfo_hapmap = alttrait_snpinfo[alttrait_snpinfo$SNP_id1 %in% hapmap_id1 | alttrait_snpinfo$SNP_id1 %in% hapmap_id2,]
  alttrait_snpinfo_hapmap$SNP_id2 = paste(alttrait_snpinfo_hapmap$CHR, alttrait_snpinfo_hapmap$POShg19, alttrait_snpinfo_hapmap$A2, alttrait_snpinfo_hapmap$A1, sep = "_")

  ## Filter to logbf over 4 and cr95set
  alttrait_logbf_over4 = alttrait_snpinfo_hapmap[alttrait_snpinfo_hapmap$logbf_over4 == 1,]
  alttrait_cr95set = alttrait_snpinfo_hapmap[alttrait_snpinfo_hapmap$cr95set == 1,]


  ## Keep only summary stat in mdd gwas info with logbf over 4 and cr95set for alt trait
  summary_stat_hapmap_keep_alttrait_logbfover4 = summary_stat_hapmap[summary_stat_hapmap$SNP_id1 %in% alttrait_logbf_over4$SNP_id1 | summary_stat_hapmap$SNP_id1 %in% alttrait_logbf_over4$SNP_id2,]
  summary_stat_hapmap_keep_alttrait_cr95set = summary_stat_hapmap[summary_stat_hapmap$SNP_id1 %in% alttrait_cr95set$SNP_id1 | summary_stat_hapmap$SNP_id1 %in% alttrait_cr95set$SNP_id2,]


  ## Add SNP id2
  summary_stat_hapmap_keep_alttrait_logbfover4$SNP_id2 = paste(summary_stat_hapmap_keep_alttrait_logbfover4$CHR, summary_stat_hapmap_keep_alttrait_logbfover4$BP, summary_stat_hapmap_keep_alttrait_logbfover4$A2, summary_stat_hapmap_keep_alttrait_logbfover4$A1, sep = "_")
  summary_stat_hapmap_keep_alttrait_cr95set$SNP_id2 = paste(summary_stat_hapmap_keep_alttrait_cr95set$CHR, summary_stat_hapmap_keep_alttrait_cr95set$BP, summary_stat_hapmap_keep_alttrait_cr95set$A2, summary_stat_hapmap_keep_alttrait_cr95set$A1, sep = "_")


  ## Make unique id that match hapmap id 1
  summary_stat_hapmap_keep_alttrait_logbfover4$hapmap_id_match = summary_stat_hapmap_keep_alttrait_logbfover4$SNP_id1
  summary_stat_hapmap_keep_alttrait_logbfover4$hapmap_id_match[summary_stat_hapmap_keep_alttrait_logbfover4$SNP_id2 %in% hapmap_id1] =  summary_stat_hapmap_keep_alttrait_logbfover4$SNP_id2[summary_stat_hapmap_keep_alttrait_logbfover4$SNP_id2 %in% hapmap_id1]

  summary_stat_hapmap_keep_alttrait_cr95set$hapmap_id_match = summary_stat_hapmap_keep_alttrait_cr95set$SNP_id1
  summary_stat_hapmap_keep_alttrait_cr95set$hapmap_id_match[summary_stat_hapmap_keep_alttrait_cr95set$SNP_id2 %in% hapmap_id1] =  summary_stat_hapmap_keep_alttrait_cr95set$SNP_id2[summary_stat_hapmap_keep_alttrait_cr95set$SNP_id2 %in% hapmap_id1]


  ## Order hapmap to match summary stat
  hapmap_match_alttrait_logbfover4 = hapmap[match(summary_stat_hapmap_keep_alttrait_logbfover4$hapmap_id_match, hapmap_id1),]
  hapmap_match_alttrait_cr95set = hapmap[match(summary_stat_hapmap_keep_alttrait_cr95set$hapmap_id_match, hapmap_id1),]

  ## Make SNP ID match
  summary_stat_hapmap_keep_alttrait_logbfover4$SNP[!summary_stat_hapmap_keep_alttrait_logbfover4$SNP %in% hapmap_match_alttrait_logbfover4$SNP] = hapmap_match_alttrait_logbfover4$SNP[!summary_stat_hapmap_keep_alttrait_logbfover4$SNP %in% hapmap_match_alttrait_logbfover4$SNP]
  summary_stat_hapmap_keep_alttrait_cr95set$SNP[!summary_stat_hapmap_keep_alttrait_cr95set$SNP %in% hapmap_match_alttrait_cr95set$SNP] = hapmap_match_alttrait_cr95set$SNP[!summary_stat_hapmap_keep_alttrait_cr95set$SNP %in% hapmap_match_alttrait_cr95set$SNP]

  prs_cs_sumstat_alttrait_logbfover4 = data.frame(SNP = summary_stat_hapmap_keep_alttrait_logbfover4$SNP,
                                                  A1 = summary_stat_hapmap_keep_alttrait_logbfover4$A1,
                                                  A2 = summary_stat_hapmap_keep_alttrait_logbfover4$A2,
                                                  OR = summary_stat_hapmap_keep_alttrait_logbfover4$OR,
                                                  SE = summary_stat_hapmap_keep_alttrait_logbfover4$SE)
  prs_cs_sumstat_alttrait_cr95set = data.frame(SNP = summary_stat_hapmap_keep_alttrait_cr95set$SNP,
                                               A1 = summary_stat_hapmap_keep_alttrait_cr95set$A1,
                                               A2 = summary_stat_hapmap_keep_alttrait_cr95set$A2,
                                               OR = summary_stat_hapmap_keep_alttrait_cr95set$OR,
                                               SE = summary_stat_hapmap_keep_alttrait_cr95set$SE)


  ## Save Reformatted GWAS, BIM, and N_gwas
  fwrite(prs_cs_sumstat_alttrait_logbfover4,
         paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),dis,"_hapmap_logbfover4_MDD_OR.txt",sep = ""),
         sep = "\t", quote = F)
  fwrite(prs_cs_sumstat_alttrait_cr95set,
         paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),dis,"_hapmap_cr95set_MDD_OR.txt",sep = ""),
         sep = "\t", quote = F)


  print("saved")
}


###########################################################                  ANX Meta Subtract                  ################################################################################


## Load library
library(data.table)
library(bigsnpr)

## Load Summary Statistics
summary_stat = fread(mdd_path("MDD_SOURCE_DATA", "/original_data/ispych_MDD_GWAS_wo_UKB23ndme/daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_Wray_FinnGen_MVPaf_2_HRC_MAF01.gz"))


## Load hapmap list
hapmap = fread(mdd_path("MDD_HAPMAP", ""))
hapmap_id1 = paste(hapmap$CHR,hapmap$BP, hapmap$A1, hapmap$A2, sep = "_")
hapmap_id2 = paste(hapmap$CHR,hapmap$BP,  hapmap$A2, hapmap$A1, sep = "_")

## Create summary stat and mdd snp info unique id
summary_stat$SNP_id1 = paste(summary_stat$CHR, summary_stat$BP, summary_stat$A1, summary_stat$A2, sep = "_")
## Keep only those in hapmap
summary_stat_hapmap = summary_stat[summary_stat$SNP_id1 %in% hapmap_id1 | summary_stat$SNP_id1 %in% hapmap_id2,]


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

## Keep only cr95set and logbf_over4 set
ANX_sumstats_hapmap_cr95set = ANX_sumstats_hapmap[ANX_sumstats_hapmap$cr95set == 1,]
ANX_sumstats_hapmap_logbf_over4 = ANX_sumstats_hapmap[ANX_sumstats_hapmap$logbf_over4 == 1,]


  ## Keep only summary stat in mdd gwas info with logbf over 4 and cr95set for alt trait
  summary_stat_hapmap_keep_alttrait_logbfover4 = summary_stat_hapmap[summary_stat_hapmap$SNP_id1 %in% ANX_sumstats_hapmap_logbf_over4$SNP_id1 | summary_stat_hapmap$SNP_id1 %in% ANX_sumstats_hapmap_logbf_over4$SNP_id2,]
  summary_stat_hapmap_keep_alttrait_cr95set = summary_stat_hapmap[summary_stat_hapmap$SNP_id1 %in% ANX_sumstats_hapmap_cr95set$SNP_id1 | summary_stat_hapmap$SNP_id1 %in% ANX_sumstats_hapmap_cr95set$SNP_id2,]


  ## Add SNP id2
  summary_stat_hapmap_keep_alttrait_logbfover4$SNP_id2 = paste(summary_stat_hapmap_keep_alttrait_logbfover4$CHR, summary_stat_hapmap_keep_alttrait_logbfover4$BP, summary_stat_hapmap_keep_alttrait_logbfover4$A2, summary_stat_hapmap_keep_alttrait_logbfover4$A1, sep = "_")
  summary_stat_hapmap_keep_alttrait_cr95set$SNP_id2 = paste(summary_stat_hapmap_keep_alttrait_cr95set$CHR, summary_stat_hapmap_keep_alttrait_cr95set$BP, summary_stat_hapmap_keep_alttrait_cr95set$A2, summary_stat_hapmap_keep_alttrait_cr95set$A1, sep = "_")


  ## Make unique id that match hapmap id 1
  summary_stat_hapmap_keep_alttrait_logbfover4$hapmap_id_match = summary_stat_hapmap_keep_alttrait_logbfover4$SNP_id1
  summary_stat_hapmap_keep_alttrait_logbfover4$hapmap_id_match[summary_stat_hapmap_keep_alttrait_logbfover4$SNP_id2 %in% hapmap_id1] =  summary_stat_hapmap_keep_alttrait_logbfover4$SNP_id2[summary_stat_hapmap_keep_alttrait_logbfover4$SNP_id2 %in% hapmap_id1]

  summary_stat_hapmap_keep_alttrait_cr95set$hapmap_id_match = summary_stat_hapmap_keep_alttrait_cr95set$SNP_id1
  summary_stat_hapmap_keep_alttrait_cr95set$hapmap_id_match[summary_stat_hapmap_keep_alttrait_cr95set$SNP_id2 %in% hapmap_id1] =  summary_stat_hapmap_keep_alttrait_cr95set$SNP_id2[summary_stat_hapmap_keep_alttrait_cr95set$SNP_id2 %in% hapmap_id1]


  ## Order hapmap to match summary stat
  hapmap_match_alttrait_logbfover4 = hapmap[match(summary_stat_hapmap_keep_alttrait_logbfover4$hapmap_id_match, hapmap_id1),]
  hapmap_match_alttrait_cr95set = hapmap[match(summary_stat_hapmap_keep_alttrait_cr95set$hapmap_id_match, hapmap_id1),]

  ## Make SNP ID match
  summary_stat_hapmap_keep_alttrait_logbfover4$SNP[!summary_stat_hapmap_keep_alttrait_logbfover4$SNP %in% hapmap_match_alttrait_logbfover4$SNP] = hapmap_match_alttrait_logbfover4$SNP[!summary_stat_hapmap_keep_alttrait_logbfover4$SNP %in% hapmap_match_alttrait_logbfover4$SNP]
  summary_stat_hapmap_keep_alttrait_cr95set$SNP[!summary_stat_hapmap_keep_alttrait_cr95set$SNP %in% hapmap_match_alttrait_cr95set$SNP] = hapmap_match_alttrait_cr95set$SNP[!summary_stat_hapmap_keep_alttrait_cr95set$SNP %in% hapmap_match_alttrait_cr95set$SNP]

  prs_cs_sumstat_alttrait_logbfover4 = data.frame(SNP = summary_stat_hapmap_keep_alttrait_logbfover4$SNP,
                                                  A1 = summary_stat_hapmap_keep_alttrait_logbfover4$A1,
                                                  A2 = summary_stat_hapmap_keep_alttrait_logbfover4$A2,
                                                  OR = summary_stat_hapmap_keep_alttrait_logbfover4$OR,
                                                  SE = summary_stat_hapmap_keep_alttrait_logbfover4$SE)
  prs_cs_sumstat_alttrait_cr95set = data.frame(SNP = summary_stat_hapmap_keep_alttrait_cr95set$SNP,
                                               A1 = summary_stat_hapmap_keep_alttrait_cr95set$A1,
                                               A2 = summary_stat_hapmap_keep_alttrait_cr95set$A2,
                                               OR = summary_stat_hapmap_keep_alttrait_cr95set$OR,
                                               SE = summary_stat_hapmap_keep_alttrait_cr95set$SE)


  ## Save Reformatted GWAS, BIM, and N_gwas
  fwrite(prs_cs_sumstat_alttrait_logbfover4,
         paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/ANXmetasubtract_hapmap_logbfover4_MDD_OR.txt"),sep = ""),
         sep = "\t", quote = F)
  fwrite(prs_cs_sumstat_alttrait_cr95set,
         paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/ANXmetasubtract_hapmap_cr95set_MDD_OR.txt"),sep = ""),
         sep = "\t", quote = F)


  print("saved")
