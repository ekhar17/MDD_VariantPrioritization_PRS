# Paths are configured through config/paths.sh or environment variables.
.mdd_repo <- Sys.getenv("MDD_CODE_ROOT", unset = getwd())
if (!file.exists(file.path(.mdd_repo, "config", "paths.R")))
  stop("Set MDD_CODE_ROOT to the repository root, or use scripts/run_r.sh.")
source(file.path(.mdd_repo, "config", "paths.R"))

library(data.table)
prs_cs_dir = mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/prs-cs_output_redo")
all_files = list.files(prs_cs_dir)

for(LD_ref in c("1kg")){
  print(LD_ref)
  prs_cs_output_full = data.frame()
  for(chr in 1:22){
    file = paste(prs_cs_dir, "/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_",LD_ref,"_LD_ref_pst_eff_a1_b0.5_phiauto_chr", chr, ".txt", sep = "")
    prs_cs_chr = fread(file)
    prs_cs_output_full = rbind(prs_cs_output_full, prs_cs_chr)
  }
  file_out=  paste(prs_cs_dir, "/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_",LD_ref,"_LD_ref_pst_eff_a1_b0.5_phiauto_allchr.txt", sep = "")
  print(file_out)
  fwrite(prs_cs_output_full, file_out, sep = "\t")
}


## TWAS SIG  - No CHR 4, 5, no 7, 12,13,17,21,22
for(LD_ref in c("1kg")){
  print(LD_ref)
  for(group in c("logbf_over4", "cr95set", "twassig")){
    prs_cs_output_full = data.frame()
    for(chr in 1:22){
      file = paste(prs_cs_dir, "/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_",LD_ref,"_LD_ref_",group,"_pst_eff_a1_b0.5_phiauto_chr", chr, ".txt", sep = "")
      prs_cs_chr = fread(file)
      prs_cs_output_full = rbind(prs_cs_output_full, prs_cs_chr)
    }
    file_out=  paste(prs_cs_dir, "/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_",LD_ref,"_LD_ref_",group,"_pst_eff_a1_b0.5_phiauto_allchr.txt", sep = "")
    print(file_out)
    fwrite(prs_cs_output_full, file_out, sep = "\t")
  }
}


## ANX LOGBF4  - No CHR 19,21,22
for(LD_ref in c("1kg")){
  print(LD_ref)
  for(trait in c("BIP", "ADHD", "SCZ", "ANXmetasubtract")){
    for(group in c("hapmap_logbf_over4", "hapmap_cr95set", "hapamap")){
      prs_cs_output_full = data.frame()
      for(chr in 1:22){
        file = paste(prs_cs_dir, "/",trait,"_",LD_ref,"_LD_ref_",group,"_redo_pst_eff_a1_b0.5_phiauto_chr", chr, ".txt", sep = "")
        prs_cs_chr = fread(file)
        prs_cs_output_full = rbind(prs_cs_output_full, prs_cs_chr)
      }
      file_out=  paste(prs_cs_dir, "/",trait,"_",LD_ref,"_LD_ref_",group,"_pst_eff_a1_b0.5_phiauto_allchr_redo.txt", sep = "")
      print(file_out)
      fwrite(prs_cs_output_full, file_out, sep = "\t")
    }
  }
}


for(LD_ref in c("1kg")){
  print(LD_ref)
  for(trait in c("BIP", "ADHD", "SCZ", "ANXmetasubtract" )){
    for(group in c("logbfover4", "cr95set")){
      prs_cs_output_full = data.frame()
      for(chr in 1:22){
        file = paste(prs_cs_dir, "/MDD_GWAS_",LD_ref,"_LD_ref_",trait,"_",group,"_set_redo_pst_eff_a1_b0.5_phiauto_chr", chr, ".txt", sep = "")
        prs_cs_chr = fread(file)
        prs_cs_output_full = rbind(prs_cs_output_full, prs_cs_chr)
      }
      file_out=  paste(prs_cs_dir, "/MDD_GWAS_",LD_ref,"_LD_ref_",trait,"_",group,"_set_pst_eff_a1_b0.5_phiauto_allchr_redo.txt", sep = "")
      print(file_out)
      fwrite(prs_cs_output_full, file_out, sep = "\t")
    }
  }
}


## TWAS SIG  - No CHR 4, 5, no 7, 12,13,17,21,22
for(LD_ref in c("1kg")){
  print(LD_ref)
  for(trait in c("BIP", "ADHD", "SCZ", "ANXmetasubtract")){
    for(group in c("logbfover4", "cr95set", "twassig")){
      prs_cs_output_full = data.frame()
      for(chr in 1:22){
        file = paste(prs_cs_dir, "/MDD_GWAS_",trait,"set_",LD_ref,"_LD_ref_",group,"_coef_redo_pst_eff_a1_b0.5_phiauto_chr", chr, ".txt", sep = "")
        prs_cs_chr = fread(file)
        prs_cs_output_full = rbind(prs_cs_output_full, prs_cs_chr)
      }
      file_out=  paste(prs_cs_dir, "/MDD_GWAS_",trait,"set_",LD_ref,"_LD_ref_",group,"_coef_pst_eff_a1_b0.5_phiauto_allchr_redo.txt", sep = "")
      print(file_out)
      fwrite(prs_cs_output_full, file_out, sep = "\t")
    }
  }
}
