# Paths are configured through config/paths.sh or environment variables.
.mdd_repo <- Sys.getenv("MDD_CODE_ROOT", unset = getwd())
if (!file.exists(file.path(.mdd_repo, "config", "paths.R")))
  stop("Set MDD_CODE_ROOT to the repository root, or use scripts/run_r.sh.")
source(file.path(.mdd_repo, "config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggpattern)
  library(colorspace)
})

library(data.table)
library(tidyverse)
library(UpSetR)

## Read in MDD cr95set and logbf4 set
MDD_logbf4 = fread(
  mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_logbf_over4.txt"))
MDD_cr95set = fread       (mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap_cr95set.txt"))
MDD_hapmap = fread(
  mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/ispych_MDD_GWAS_wo_UKB23ndme_hapmap_mdd_group_overlap.txt"))


## Read in BIP cr95set and logbf4 set
BIP_logbf4 = fread(paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),"BIP","_hapmap_logbf_over4.txt",sep = ""))
BIP_cr95set = fread(paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),"BIP","_hapmap_cr95set.txt",sep = ""))
BIP_hapmap = fread(paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),"BIP","_hapmap.txt",sep = ""))

## Read in ADHD cr95set and logbf4 set
adhd_logbf4 = fread(paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),"ADHD","_hapmap_logbf_over4.txt",sep = ""))
adhd_cr95set = fread(paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),"ADHD","_hapmap_cr95set.txt",sep = ""))
adhd_hapmap = fread(paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),"ADHD","_hapmap.txt",sep = ""))

## Read in SCZ cr95set and logbf4 set
scz_logbf4 = fread(paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),"SCZ","_hapmap_logbf_over4.txt",sep = ""))
scz_cr95set = fread(paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),"SCZ","_hapmap_cr95set.txt",sep = ""))
scz_hapmap = fread(paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),"SCZ","_hapmap.txt",sep = ""))


## Read in SCZ cr95set and logbf4 set
anx_logbf4 = fread(paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),"ANXmetasubtract","_hapmap_logbf_over4.txt",sep = ""))
anx_cr95set = fread(paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),"ANXmetasubtract","_hapmap_cr95set.txt",sep = ""))
anx_hapmap = fread(paste(mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/clean_sumstat/"),"ANXmetasubtract","_hapmap.txt",sep = ""))

nrow(anx_hapmap)
nrow(scz_hapmap)
nrow(adhd_hapmap)
nrow(MDD_hapmap)
nrow(BIP_hapmap)

listInput_cr95set <- list( MDD_cr95set = MDD_cr95set$SNP,
                           BIP_cr95set = BIP_cr95set$SNP,
                           ADHD_cr95set = adhd_cr95set$SNP,
                           SCZ_cr95set = scz_cr95set$SNP,
                           ANX_cr95set = anx_cr95set$SNP)
matrixInput_cr95set =fromList(listInput_cr95set)


listInput_logbf4 <- list( MDD_logbf4 = MDD_logbf4$SNP,
                          BIP_logbf4 = BIP_logbf4$SNP,
                          ADHD_logbf4 = adhd_logbf4$SNP,
                          SCZ_logbf4 = scz_logbf4$SNP,
                          ANX_logbf4 = anx_logbf4$SNP)
matrixInput_logbf4 =fromList(listInput_logbf4)

png(mdd_path("MDD_PROJECT", "/08.Generate_Figures/fig_v1/cr95se_upset_tmp.png"), width = 1000, height = 500)
upset(matrixInput_cr95set,  intersect = c("SCZ_cr95set",
                                         "BIP_cr95set",
                                         "ADHD_cr95set",
                                         "ANX_cr95set",
                                         "MDD_cr95set"),
      name='cr95set SNP group',
      width_ratio=0.15, sort_sets=FALSE,height_ratio=0.3,
      themes=upset_modify_themes(
        list(
          'intersections_matrix'=theme(text=element_text(size=16))
        )),
      base_annotations=list(
        'Intersection size'=intersection_size(
          mapping=aes(fill='bars_color'),
          text=list(
            vjust=-0.5,
            hjust=0.2,
            angle=45,
            size=4
          )
        ) + scale_fill_manual(values=c('bars_color'='#F0A430FF'),guide='none')
      ),
      #annotations = list(
      #   'GWAS P-value'=list(
      #      aes=aes(x=intersection, y=P),
      #     geom=geom_boxplot(na.rm=TRUE)
      #),
      #'GWAS Beta'=list(
      #   aes=aes(x=intersection, y=BETA),
      #  geom=list(
      #     #geom_jitter(aes(color=log10(votes)), na.rm=TRUE),
      #    geom_violin(alpha=0.5, na.rm=TRUE)
      #)
      #)
      #),
      set_sizes=(
        upset_set_size(geom=geom_bar(aes(fill='groupcolor')))
        + scale_fill_manual(values=c('groupcolor'='#768048FF'), guide='none')
        # + geom_text(aes(label=..count..), hjust=1.1, vjust=-0.2, stat='count')
        + theme(axis.text.x=element_text(angle=90,size=8),
                #axis.title.y = element_text(size=14) # Increase font size for "SNP group size" label
                axis.title.x = element_text(size=14)
        )   # Increase font size for "SNP group" label)
        + ylab('SNP group size')
        # + scale_y_continuous(trans=reverse_log_trans())
        #  + scale_x_continuous(trans=reverse_log_trans())
      )
)
dev.off()

png(mdd_path("MDD_PROJECT", "/08.Generate_Figures/fig_v1/logbf4_upset_tmp.png"), width = 1000, height = 500)
upset(matrixInput_logbf4,  intersect = c("SCZ_logbf4",
  "BIP_logbf4",
  "ADHD_logbf4",
  "ANX_logbf4",
  "MDD_logbf4"),
      name='logBF4 SNP group',
width_ratio=0.15, sort_sets=FALSE,height_ratio=0.3,
      themes=upset_modify_themes(
        list(
          'intersections_matrix'=theme(text=element_text(size=16))
        )),
      base_annotations=list(
        'Intersection size'=intersection_size(
          mapping=aes(fill='bars_color'),
          text=list(
            vjust=-0.5,
            hjust=0.2,
            angle=45,
            size=4
          )
        ) + scale_fill_manual(values=c('bars_color'='#F0A430FF'),guide='none')
      ),
      #annotations = list(
      #   'GWAS P-value'=list(
      #      aes=aes(x=intersection, y=P),
      #     geom=geom_boxplot(na.rm=TRUE)
      #),
      #'GWAS Beta'=list(
      #   aes=aes(x=intersection, y=BETA),
      #  geom=list(
      #     #geom_jitter(aes(color=log10(votes)), na.rm=TRUE),
      #    geom_violin(alpha=0.5, na.rm=TRUE)
      #)
      #)
      #),
      set_sizes=(
        upset_set_size(geom=geom_bar(aes(fill='groupcolor')))
        + scale_fill_manual(values=c('groupcolor'='#768048FF'), guide='none')
        # + geom_text(aes(label=..count..), hjust=1.1, vjust=-0.2, stat='count')
        + theme(axis.text.x=element_text(angle=90,size=8),
                #axis.title.y = element_text(size=14) # Increase font size for "SNP group size" label
                axis.title.x = element_text(size=14)
        )   # Increase font size for "SNP group" label)
        + ylab('SNP group size')
        # + scale_y_continuous(trans=reverse_log_trans())
        #  + scale_x_continuous(trans=reverse_log_trans())
      )
)
dev.off()


listInput <- list(
  BIP_logbf4 = BIP_logbf4$SNP, BIP_cr95set = BIP_cr95set$SNP,
  ADHD_logbf4 =  adhd_logbf4$SNP,  ADHD_cr95set = adhd_cr95set$SNP,
  SCZ_logbf4 =  scz_logbf4$SNP,  SCZ_cr95set = scz_cr95set$SNP,
  ANX_logbf4 =    anx_logbf4$SNP,  ANX_cr95set = anx_cr95set$SNP,
  MDD_logbf4 = MDD_logbf4$SNP, MDD_cr95set = MDD_cr95set$SNP)
matrixInput =fromList(listInput)
upset(matrixInput, sets =  c(
  "BIP_logbf4",
  "BIP_cr95set",
  "ADHD_logbf4",
  "ADHD_cr95set",
  "SCZ_logbf4",
  "SCZ_cr95set",
  "ANX_logbf4",
  "ANX_cr95set",
  "MDD_logbf4",

  "MDD_cr95set"),
  keep.order = T)

upset(b1, set_vars,
      name='SNP group', width_ratio=0.15, sort_sets=FALSE,height_ratio=0.3,
      themes=upset_modify_themes(
        list(
          'intersections_matrix'=theme(text=element_text(size=16))
        )),
      base_annotations=list(
        'Intersection size'=intersection_size(
          mapping=aes(fill='bars_color'),
          text=list(
            vjust=-0.5,
            hjust=0.2,
            angle=45,
            size=4
          )
        ) + scale_fill_manual(values=c('bars_color'='#F0A430FF'),guide='none')
      ),
      #annotations = list(
      #   'GWAS P-value'=list(
      #      aes=aes(x=intersection, y=P),
      #     geom=geom_boxplot(na.rm=TRUE)
      #),
      #'GWAS Beta'=list(
      #   aes=aes(x=intersection, y=BETA),
      #  geom=list(
      #     #geom_jitter(aes(color=log10(votes)), na.rm=TRUE),
      #    geom_violin(alpha=0.5, na.rm=TRUE)
      #)
      #)
      #),
      set_sizes=(
        upset_set_size(geom=geom_bar(aes(fill='groupcolor')))
        + scale_fill_manual(values=c('groupcolor'='#768048FF'), guide='none')
        # + geom_text(aes(label=..count..), hjust=1.1, vjust=-0.2, stat='count')
        + theme(axis.text.x=element_text(angle=90,size=8),
                #axis.title.y = element_text(size=14) # Increase font size for "SNP group size" label
                axis.title.x = element_text(size=14)
        )   # Increase font size for "SNP group" label)
        + ylab('SNP group size')
        # + scale_y_continuous(trans=reverse_log_trans())
        #  + scale_x_continuous(trans=reverse_log_trans())
      )
)
matrixOutput = ComplexHeatmap::list_to_matrix(listInput)
matrixTest = ComplexHeatmap::make_comb_mat(matrixOutput)
ComplexHeatmap::UpSet(matrixTest,
                      set_order = c(
                        "MDD_cr95set",
                        "MDD_logbf4",
                        "ANX_cr95set",
                        "ANX_logbf4",
                        "ADHD_cr95set",
                        "ADHD_logbf4",
                        "BIP_cr95set",
                        "BIP_logbf4",
                        "SCZ_cr95set",
                        "SCZ_logbf4"
                      ))
ComplexHeatmap::comb_size(matrixTest)
ComplexHeatmap::UpSet(matrixTest[ComplexHeatmap::comb_degree(matrixTest) == 1|
                                   ComplexHeatmap::comb_degree(matrixTest) < 4],
                      set_order = c(
                        "MDD_cr95set",
                        "MDD_logbf4",
                        "ANX_cr95set",
                        "ANX_logbf4",
                        "ADHD_cr95set",
                        "ADHD_logbf4",
                        "BIP_cr95set",
                        "BIP_logbf4",
                        "SCZ_cr95set",
                        "SCZ_logbf4"
                      ))
ComplexHeatmap::UpSet(matrixTest[ComplexHeatmap::comb_size(matrixTest) > 4e2|
                                   ComplexHeatmap::comb_degree(matrixTest)==1],
                      set_order = c(
                        "MDD_cr95set",
                        "MDD_logbf4",
                        "ANX_cr95set",
                        "ANX_logbf4",
                        "ADHD_cr95set",
                        "ADHD_logbf4",
                        "BIP_cr95set",
                        "BIP_logbf4",
                        "SCZ_cr95set",
                        "SCZ_logbf4"
                      ))
ComplexHeatmap::comb_size(matrixTest) == 1
ComplexHeatmap::comb_degree(matrixTest) == 1


df_long <- bind_rows(lapply(names(listInput), function(setname) {
  data.frame(SNP = listInput[[setname]], set = setname)
}))

df_binary <- df_long %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = set, values_from = value, values_fill = list(value = 0))

# ✅ Filter to only rows (SNPs) that are present in MDD_cr95set
df_subset <- df_binary %>%
  filter(MDD_cr95set == 1|
           MDD_logbf4 == 1)

# Plot UpSet showing only intersections that include MDD_cr95set
library(ComplexUpset)
upset(df_subset, intersect = names(listInput), min_size = 100)

upset(matrixInput, nintersects = 100)  # or Inf to show all

upset(matrixInput, intersections  = c(
  "BIP_logbf4",
  "BIP_cr95set",
  "ADHD_logbf4",
  "ADHD_cr95set",
  "SCZ_logbf4",
  "SCZ_cr95set",
  "ANX_logbf4",
  "ANX_cr95set",
  "MDD_logbf4",

  "MDD_cr95set"),
  keep.order = T)


upset(matrixInput,
      nintersects = 4)

four_set_snps <- matrixInput[rowSums(matrixInput) == 10, ]

# Upset plot only for 4-set overlaps
upset(four_set_snps,
      sets = colnames(four_set_snps),
      keep.order = TRUE)


colSums(matrixOutput)
