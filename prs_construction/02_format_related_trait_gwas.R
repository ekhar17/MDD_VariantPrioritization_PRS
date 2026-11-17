# Paths are configured through config/paths.sh or environment variables.
.mdd_repo <- Sys.getenv("MDD_CODE_ROOT", unset = getwd())
if (!file.exists(file.path(.mdd_repo, "config", "paths.R")))
  stop("Set MDD_CODE_ROOT to the repository root, or use scripts/run_r.sh.")
source(file.path(.mdd_repo, "config", "paths.R"))

## Load library
library(data.table)
library(bigsnpr)

dir = mdd_path("MDD_SOURCE_DATA", "/derived_data/procedure_data/fgannosumstats/")
end = "_GWAS_sumstats.QC.maf01info09_95crlgbf4.txt"

## Load hapmap list
hapmap = fread(mdd_path("MDD_HAPMAP", ""))
hapmap_id1 = paste(hapmap$CHR,hapmap$BP, hapmap$A1, hapmap$A2, sep = "_")
hapmap_id2 = paste(hapmap$CHR,hapmap$BP,  hapmap$A2, hapmap$A1, sep = "_")


disorders = c("ADHD", "BIP", "SCZ")

for(dis in disorders){

  ## Load summary stat
  file = paste(dir, dis, end, sep = "")
  summary_stat = fread(file)
  print("loaded sum stats")


  ## Add Unique ID
  summary_stat$SNP_id1 = paste(summary_stat$CHR, summary_stat$POShg19, summary_stat$A1, summary_stat$A2, sep = "_")

  ## Keep only those in hapmap
  summary_stat_hapmap = summary_stat[summary_stat$SNP_id1 %in% hapmap_id1 | summary_stat$SNP_id1 %in% hapmap_id2,]
  summary_stat_hapmap$SNP_id2 = paste(summary_stat_hapmap$CHR, summary_stat_hapmap$POShg19, summary_stat_hapmap$A2, summary_stat_hapmap$A1, sep = "_")
  summary_stat_hapmap$hapmap_id_match = summary_stat_hapmap$SNP_id1
  summary_stat_hapmap$hapmap_id_match[summary_stat_hapmap$SNP_id2 %in% hapmap_id1] = summary_stat_hapmap$SNP_id2[summary_stat_hapmap$SNP_id2 %in% hapmap_id1]
  print("only hapmap")

  ## Keep only cr95set and logbf_over4 set
  summary_stat_hapmap_cr95set = summary_stat_hapmap[summary_stat_hapmap$cr95set == 1,]
  summary_stat_hapmap_logbf_over4 = summary_stat_hapmap[summary_stat_hapmap$logbf_over4 == 1,]

  ## Order hapmap to match summary stat
  hapmap_match = hapmap[match(summary_stat_hapmap$hapmap_id_match, hapmap_id1),]
  hapmap_match_logbf_over4 = hapmap[match(summary_stat_hapmap_logbf_over4$hapmap_id_match, hapmap_id1),]
  hapmap_match_cr95set = hapmap[match(summary_stat_hapmap_cr95set$hapmap_id_match, hapmap_id1),]

  # ADD RSID
  summary_stat_hapmap$SNP = hapmap_match$SNP
  summary_stat_hapmap_logbf_over4$SNP = hapmap_match_logbf_over4$SNP
  summary_stat_hapmap_cr95set$SNP = hapmap_match_cr95set$SNP
  print("added rsid")

  ## Keep only columns needed for PRS CS
  prs_cs_sumstat = data.frame(SNP = summary_stat_hapmap$SNP,
                              A1 = summary_stat_hapmap$A1,
                              A2 = summary_stat_hapmap$A2,
                              BETA = summary_stat_hapmap$B,
                              SE = summary_stat_hapmap$SE)
  prs_cs_sumstat_logbf_over4 = data.frame(SNP = summary_stat_hapmap_logbf_over4$SNP,
                                          A1 = summary_stat_hapmap_logbf_over4$A1,
                                          A2 = summary_stat_hapmap_logbf_over4$A2,
                                          BETA = summary_stat_hapmap_logbf_over4$B,
                                          SE = summary_stat_hapmap_logbf_over4$SE)
  prs_cs_sumstat_cr95set = data.frame(SNP = summary_stat_hapmap_cr95set$SNP,
                                      A1 = summary_stat_hapmap_cr95set$A1,
                                      A2 = summary_stat_hapmap_cr95set$A2,
                                      BETA = summary_stat_hapmap_cr95set$B,
                                      SE = summary_stat_hapmap_cr95set$SE)


  ## Calculate N_gwas
  N_gwas = max(summary_stat_hapmap$NCASE + summary_stat_hapmap$NCONTROL , na.rm = T)
  print(N_gwas)
  print(max(summary_stat_hapmap_logbf_over4$NCASE + summary_stat_hapmap_logbf_over4$NCONTROL, na.rm = T))
  print(max(summary_stat_hapmap_cr95set$NCASE + summary_stat_hapmap_cr95set$NCONTROL, na.rm = T))

  ## Save Reformatted GWAS, BIM, and N_gwas
  fwrite(prs_cs_sumstat,
         paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),dis,"_hapmap.txt",sep = ""),
         sep = "\t", quote = F)
  fwrite(prs_cs_sumstat_logbf_over4,
         paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),dis,"_hapmap_logbf_over4.txt",sep = ""),
         sep = "\t", quote = F)
  fwrite(prs_cs_sumstat_cr95set,
         paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),dis,"_hapmap_cr95set.txt",sep = ""),
         sep = "\t", quote = F)
  write.table(N_gwas,
              paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),dis,"_NGWAS.txt",sep = ""),
              col.names = F, row.names = F, quote = F)

  print("saved")
}


###########################################################                  ANX Meta Subtract                  ################################################################################


library(data.table)


## Load hapmap list
hapmap = fread(mdd_path("MDD_HAPMAP", ""))
hapmap_id1 = paste(hapmap$CHR,hapmap$BP, hapmap$A1, hapmap$A2, sep = "_")
hapmap_id2 = paste(hapmap$CHR,hapmap$BP,  hapmap$A2, hapmap$A1, sep = "_")

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
ANX_sumstats_hapmap$hapmap_id_match = ANX_sumstats_hapmap$SNP_id1
ANX_sumstats_hapmap$hapmap_id_match[ANX_sumstats_hapmap$SNP_id2 %in% hapmap_id1] = ANX_sumstats_hapmap$SNP_id2[ANX_sumstats_hapmap$SNP_id2 %in% hapmap_id1]
print("only hapmap")

## Keep only cr95set and logbf_over4 set
ANX_sumstats_hapmap_cr95set = ANX_sumstats_hapmap[ANX_sumstats_hapmap$cr95set == 1,]
ANX_sumstats_hapmap_logbf_over4 = ANX_sumstats_hapmap[ANX_sumstats_hapmap$logbf_over4 == 1,]

## Order hapmap to match summary stat
hapmap_match = hapmap[match(ANX_sumstats_hapmap$hapmap_id_match, hapmap_id1),]
hapmap_match_logbf_over4 = hapmap[match(ANX_sumstats_hapmap_logbf_over4$hapmap_id_match, hapmap_id1),]
hapmap_match_cr95set = hapmap[match(ANX_sumstats_hapmap_cr95set$hapmap_id_match, hapmap_id1),]

# ADD RSID
ANX_sumstats_hapmap$SNP = hapmap_match$SNP
ANX_sumstats_hapmap_logbf_over4$SNP = hapmap_match_logbf_over4$SNP
ANX_sumstats_hapmap_cr95set$SNP = hapmap_match_cr95set$SNP
all.equal(ANX_sumstats_hapmap$RSID, ANX_sumstats_hapmap$SNP)
all.equal(ANX_sumstats_hapmap_logbf_over4$RSID, ANX_sumstats_hapmap_logbf_over4$SNP)
all.equal(ANX_sumstats_hapmap_cr95set$RSID, ANX_sumstats_hapmap_cr95set$SNP)
print("added rsid")

## Keep only columns needed for PRS CS
prs_cs_sumstat = data.frame(SNP = ANX_sumstats_hapmap$SNP,
                            A1 = ANX_sumstats_hapmap$A1_hg19,
                            A2 = ANX_sumstats_hapmap$A2_hg19,
                            BETA = ANX_sumstats_hapmap$Z,
                            P = ANX_sumstats_hapmap$P)
prs_cs_sumstat_logbf_over4 = data.frame(SNP = ANX_sumstats_hapmap_logbf_over4$SNP,
                                        A1 = ANX_sumstats_hapmap_logbf_over4$A1_hg19,
                                        A2 = ANX_sumstats_hapmap_logbf_over4$A2_hg19,
                                        BETA = ANX_sumstats_hapmap_logbf_over4$Z,
                                        P = ANX_sumstats_hapmap_logbf_over4$P)
prs_cs_sumstat_cr95set = data.frame(SNP = ANX_sumstats_hapmap_cr95set$SNP,
                                    A1 = ANX_sumstats_hapmap_cr95set$A1_hg19,
                                    A2 = ANX_sumstats_hapmap_cr95set$A2_hg19,
                                    BETA = ANX_sumstats_hapmap_cr95set$Z,
                                    P = ANX_sumstats_hapmap_cr95set$P)


## Calculate N_gwas
N_gwas = max(ANX_sumstats_hapmap$Ncase + ANX_sumstats_hapmap$Ncontrol , na.rm = T)
print(N_gwas)
print(max(ANX_sumstats_hapmap_logbf_over4$Ncase + ANX_sumstats_hapmap_logbf_over4$Ncontrol, na.rm = T))
print(max(ANX_sumstats_hapmap_cr95set$Ncase + ANX_sumstats_hapmap_cr95set$Ncontrol, na.rm = T))

## Save Reformatted GWAS, BIM, and N_gwas
fwrite(prs_cs_sumstat,
       paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/ANXmetasubtract_hapmap.txt"),sep = ""),
       sep = "\t", quote = F)
fwrite(prs_cs_sumstat_logbf_over4,
       paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/ANXmetasubtract_hapmap_logbf_over4.txt"),sep = ""),
       sep = "\t", quote = F)
fwrite(prs_cs_sumstat_cr95set,
       paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/ANXmetasubtract_hapmap_cr95set.txt"),sep = ""),
       sep = "\t", quote = F)
write.table(N_gwas,
            paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/ANXmetasubtract_NGWAS.txt"),sep = ""),
            col.names = F, row.names = F, quote = F)

print("saved")
