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

#fwrite(results, "${MDD_PROJECT}/05.PRS_CS_construction/data/auc_eur_res.csv")
results = fread( mdd_path("MDD_PROJECT", "/05.PRS_CS_construction/data/auc_eur_res.csv"))

results$class = gsub("_.*", "", results$type)
results$prs_type = gsub(".*Coef_", "", results$type)
clean_trait = function(df, col = "column",column_input){
  df$model = df[,column_input]

  df[grepl("MDD", df$model), col] = "MDD"
  df[grepl("BIP", df$model), col] = "BIP"
  df[grepl("SCZ", df$model), col] = "SCZ"
  df[grepl("ADHD", df$model), col] = "ADHD"
  df[grepl("ANX", df$model), col] = "ANX"
  return(df)
}
clean_type = function(df, col = "column", column_input){
  df$model = df[,column_input]
  df[grepl("hapmap", df$model), col] = "hapmap"
  df[grepl("logbf4", df$model), col] = "logbf4"
  df[grepl("logbf_over4", df$model), col] = "logbf4"
  df[grepl("cr95set", df$model), col] = "cr95set"
  df[grepl("twassig", df$model), col] = "twassig"
  df[grepl("TWAS", df$model), col] = "twassig"
  return(df)
}
results = clean_type(results %>% as.data.frame(), col = "group_type", column_input = "type")
results = clean_trait(results %>% as.data.frame(), col = "coef", column_input = "type")
results = clean_trait(results %>% as.data.frame(), col = "set", column_input = "type")
results = clean_trait(results %>% as.data.frame(), col = "trait", column_input = "type")
results[results $class == "AltTraitSet-MDDCoef","coef"] = "MDD"
results[results $class == "MDDSet-AltTraitCoef","set"] = "MDD"

performance = results
performance$type_type = "Trait Set: Trait Coef"
performance$type_type[performance$class == "AltTraitSet-MDDCoef"] = "AltTrait Set : MDD Coef"
performance$type_type[performance$class == "MDDSet-AltTraitCoef"] = "MDD Set : AltTrait Coef"

performance$coef = factor(performance$coef, levels = c("MDD", "ADHD", "ANX", "BIP", "SCZ"))
performance_stats = performance %>% group_by(type,class,trait,prs_type, set, coef, group_type,type_type) %>% summarise(
  mean_auc = mean(diff_AUC), max_auc = max(diff_AUC)
)

performance$class = factor(performance$class,
                           levels = c("MDDSet-MDDCoef","AltTraitSet-AltTraitCoef",
                                      "AltTraitSet-MDDCoef","MDDSet-AltTraitCoef"   ))
library(ggpattern)


performance$trait2 = factor(performance$trait,
                            levels = c("MDD", "ADHD", "ANX", "BIP", "SCZ", "EXTRA"))
performance2 = performance
performance2[performance2$class == "MDDSet-MDDCoef",]$type_type = "MDD SNP Set : MDD Coef"
v1 = performance2[performance2$class == "MDDSet-MDDCoef",]
v1_1 = v1
v1_1$trait2 = "ADHD"
v1_2 = v1
v1_2$trait2 = "ANX"
v1_3 = v1
v1_3$trait2 = "BIP"
v1_4 = v1
v1_4$trait2 = "SCZ"
performance3 = rbind(performance2, v1_1, v1_2, v1_3, v1_4)
performance3$type_type = factor(performance3$type_type, levels= c("AltTrait Set : MDD Coef",
                                                                  "MDD SNP Set : MDD Coef",
                                                                  "MDD Set : AltTrait Coef",
                                                                  "Trait Set: Trait Coef"))


ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/new_fig/AltTraitCoef_orSNP_AllEUR_barplot_wpattern_nolabel_small.png"),
       performance3 %>% filter(class == "AltTraitSet-MDDCoef"  | class == "MDDSet-AltTraitCoef") %>%
         filter(group_type != "twassig") %>% filter(group_type != "hapmap") %>% filter(trait2 != "MDD") %>%
         ggplot( aes(y = diff_AUC*100, fill = group_type, x = type_type))+geom_bar_pattern(stat = "identity",position = "dodge",
                                                                                           aes(pattern = class)
                                                                                           ,
                                                                                           pattern_spacing = 0.1,   # smaller = tighter hatching
                                                                                           pattern_size = 0.2
         ) + facet_grid(~trait2,
                        scales = "free_x")+
         scale_fill_manual(values = c("pink1", "dodgerblue3"))+theme_bw()+
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),

               axis.ticks.x = element_blank(),

               text = element_text(size = 20, color = "black"),

               axis.text = element_text(colour="black"),

               strip.text.x = element_text(size = 22, colour = "black"),

               legend.key.size = unit(1.5, "cm"))+
         scale_pattern_manual(values = c(
           "AltTraitSet-MDDCoef" = "circle",
           "MDDSet-AltTraitCoef" = "crosshatch",
           "MDDSet-MDDCoef" = "none"
         )) + labs(y = "Difference in AUC (%)", x = "Class", fill = "Set", pattern = "Class", linetype = "")+
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.00, 0.1)))+
         scale_x_discrete(labels =
                            c("AltTrait Set",
                              "AltTrait Coef"
                            ))+
         geom_hline(aes(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_cr95set_prs"],
                        linetype = "MDD cr95set"),
                    color =  "maroon3", size =1.3)+
         geom_hline(aes(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_logbf_over4_prs"],
                        linetype = "MDD logbf4")  ,color = "darkblue", size =1.3) +
         scale_linetype_manual(values = c(
           "MDD cr95set" = "dashed",
           "MDD logbf4"  = "dashed"
         )),
       height = 8, width = 18)


############ Scratchwork


ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/fig_v4/AltTraitCoef_orSNP_AllEUR_barplot_wpattern_nolabel.png"),
       performance3 %>% filter(class == "AltTraitSet-MDDCoef"  | class == "MDDSet-AltTraitCoef") %>%
         filter(group_type != "twassig") %>% filter(group_type != "hapmap") %>% filter(trait2 != "MDD") %>%
         ggplot( aes(y = diff_AUC*100, fill = group_type, x = type_type))+geom_bar_pattern(stat = "identity",position = "dodge",
                                                                                           aes(pattern = class)
                                                                                           ,
                                                                                           pattern_spacing = 0.1,   # smaller = tighter hatching
                                                                                           pattern_size = 0.2
         ) + facet_grid(~trait2,
                        scales = "free_x")+
         scale_fill_manual(values = c("pink1", "dodgerblue3"))+theme_bw()+
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),

               axis.ticks.x = element_blank(),

               text = element_text(size = 20, color = "black"),

               axis.text = element_text(colour="black"),

               strip.text.x = element_text(size = 22, colour = "black"),

               legend.key.size = unit(1.5, "cm"))+
         scale_pattern_manual(values = c(
           "AltTraitSet-MDDCoef" = "circle",
           "MDDSet-AltTraitCoef" = "crosshatch",
           "MDDSet-MDDCoef" = "none"
         )) + labs(y = "Difference in AUC (%)", x = "Class", fill = "Set", pattern = "Class", linetype = "")+
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.00, 0.1)))+
         scale_x_discrete(labels =
                            c("AltTrait SNP Set",
                              "AltTrait Coef"
                            ))+
         geom_hline(aes(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_cr95set_prs"],
                        linetype = "MDD cr95set"),
                    color =  "maroon3")+
         geom_hline(aes(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_logbf_over4_prs"],
                        linetype = "MDD logbf4")  ,color = "darkblue") +
         scale_linetype_manual(values = c(
           "MDD cr95set" = "dashed",
           "MDD logbf4"  = "dashed"
         ))+
         guides(
           linetype = guide_legend(
             override.aes = list(size = 5)  # legend-only thickness
           )
         ),
       height = 8, width = 18)
