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
performance$type_type = "Trait SNP Set: Trait Coef"
performance$type_type[performance$class == "AltTraitSet-MDDCoef"] = "AltTrait SNP Set : MDD Coef"
performance$type_type[performance$class == "MDDSet-AltTraitCoef"] = "MDD SNP Set : AltTrait Coef"

performance$coef = factor(performance$coef, levels = c("MDD", "ADHD", "ANX", "BIP", "SCZ"))
performance_stats = performance %>% group_by(type,class,trait,prs_type, set, coef, group_type,type_type) %>% summarise(
  mean_auc = mean(diff_AUC), max_auc = max(diff_AUC)
)

performance$class = factor(performance$class,
                           levels = c("MDDSet-MDDCoef","AltTraitSet-AltTraitCoef",
                                      "AltTraitSet-MDDCoef","MDDSet-AltTraitCoef"   ))
library(ggpattern)

ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/new_fig/TraitSNPSet_TraitCoef_AllEUR_barplot_wpattern.png"),
       performance %>% filter(class == "AltTraitSet-AltTraitCoef" | class == "MDDSet-MDDCoef") %>%
         filter(group_type != "twassig") %>%
         ggplot( aes(y = diff_AUC*100, fill = group_type, x = coef))+geom_bar_pattern(stat = "identity",position = "dodge",
                                                                                      aes(pattern = class)) +
         #facet_wrap(~type_type)+
         scale_pattern_manual(values = c( "none","stripe"), labels = c( "MDD","Alt Trait"))+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_cr95set_prs"],
                    color =  "maroon3", lty = 2)+
         geom_text(aes(5.77,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_cr95set_prs"],
                       label = "MDD cr95set", vjust = -1), size = 6)+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_logbf_over4_prs"],
                    color = "darkblue", lty = 2)+
         geom_text(aes(5.77,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_logbf_over4_prs"],
                       label = "MDD logbf4", vjust = -1), size = 6)+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_hapmap_prs"],
                    color = "gold4", lty = 2)+
         geom_text(aes(5.77,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_hapmap_prs"],
                       label = "MDD hapmap", vjust = -1), size = 6)+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_twassig_prs"],
                    color = "darkgreen", lty = 2)+
         geom_text(aes(5.77,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_twassig_prs"],
                       label = "MDD twassig", vjust = -1), size = 6)+
         scale_fill_manual(values = c("pink1", "gold2", "dodgerblue3"))+theme_bw()+
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),

               axis.ticks.x = element_blank(),

               text = element_text(size = 20, color = "black"),

               axis.text = element_text(colour="black"),

               strip.text.x = element_text(size = 22, colour = "black"),

               legend.key.size = unit(1.5, "cm"))+
         geom_text(
           data = performance_stats[
             (performance_stats$class %in% c("AltTraitSet-AltTraitCoef" ,"MDDSet-MDDCoef") )
             &
               performance_stats$group_type != "twassig", ],
           aes(x = coef, y = max_auc*100, label = round(mean_auc*100, 1)),
           position = position_dodge(width = 0.9), # aligns with boxplot dodging
           vjust = -1,
           size = 6
         ) + labs(y = "Difference in AUC (%)", x = "Disorder", fill = "Set", pattern = "Class",
                  title = "Trait Set - Trait Coef")+
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.00, 0.1)))+expand_limits(x=6.1)+
         guides(
           fill = guide_legend(override.aes = list(pattern = "none")),
           pattern = guide_legend(override.aes = list(fill = "white",
                                                      color = "black"))
         ),
       height = 8, width = 18)


ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/new_fig/TraitSNPSet_TraitCoef_AllEUR_barplot_wpattern_small.png"),
       performance %>% filter(class == "AltTraitSet-AltTraitCoef" | class == "MDDSet-MDDCoef") %>%
         filter(group_type != "twassig") %>%
         ggplot( aes(y = diff_AUC*100, fill = group_type, x = coef))+geom_bar_pattern(stat = "identity",position = "dodge",
                                                                                      aes(pattern = class)) +
         #facet_wrap(~type_type)+
         scale_pattern_manual(values = c( "none","stripe"), labels = c( "MDD","Alt Trait"))+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_cr95set_prs"],
                    color =  "maroon3", lty = 2)+
         geom_text(aes(5.77,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_cr95set_prs"],
                       label = "MDD cr95set", vjust = -1), size = 6)+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_logbf_over4_prs"],
                    color = "darkblue", lty = 2)+
         geom_text(aes(5.77,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_logbf_over4_prs"],
                       label = "MDD logbf4", vjust = -1), size = 6)+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_hapmap_prs"],
                    color = "gold4", lty = 2)+
         geom_text(aes(5.77,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_hapmap_prs"],
                       label = "MDD hapmap", vjust = -1), size = 6)+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_twassig_prs"],
                    color = "darkgreen", lty = 2)+
         geom_text(aes(5.77,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_twassig_prs"],
                       label = "MDD twassig", vjust = -1), size = 6)+
         scale_fill_manual(values = c("pink1", "gold2", "dodgerblue3"))+theme_bw()+
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),

               axis.ticks.x = element_blank(),

               text = element_text(size = 20, color = "black"),

               axis.text = element_text(colour="black"),

               strip.text.x = element_text(size = 22, colour = "black"),

               legend.key.size = unit(1.5, "cm"))+
         geom_text(
           data = performance_stats[
             (performance_stats$class %in% c("AltTraitSet-AltTraitCoef" ,"MDDSet-MDDCoef") )
             &
               performance_stats$group_type != "twassig", ],
           aes(x = coef, y = max_auc*100, label = round(mean_auc*100, 1)),
           position = position_dodge(width = 0.9), # aligns with boxplot dodging
           vjust = -1,
           size = 6
         ) + labs(y = "Difference in AUC (%)", x = "Disorder", fill = "Set", pattern = "Class",
                  title = "TraitSet-TraitCoef")+
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.00, 0.1)))+expand_limits(x=6.1)+
         guides(
           fill = guide_legend(override.aes = list(pattern = "none")),
           pattern = guide_legend(override.aes = list(fill = "white",
                                                      color = "black"))
         ),

       height = 8, width = 8)


performance$trait2 = factor(performance$trait,
                            levels = c("MDD", "ADHD", "ANX", "BIP", "SCZ", "EXTRA"))

#ggsave("${MDD_PROJECT}/08.Generate_Figures/fig_v2/AltTraitSNPSet_MDDCoef_viceversa_ALLEUR_barplot_wpattern.png",

png(mdd_path("MDD_PROJECT", "/08.Generate_Figures/fig_v2/AltTraitSNPSet_MDDCoef_viceversa_ALLEUR_barplot_wpattern.png"),
    width = 19, height = 8, units = "in", res = 300)

performance %>% filter(class == "AltTraitSet-MDDCoef" | class == "MDDSet-AltTraitCoef") %>%
  filter(group_type != "twassig") %>%
  ggplot( aes(y = diff_AUC*100, fill = group_type, x = trait2))+geom_bar_pattern(stat = "identity",position = "dodge",
                                                                                 aes(pattern = class)) + facet_wrap(~type_type,
                                                                                                                    scales = "free_x")+
  geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_cr95set_prs"],
             color = "maroon3", lty = 2)+
  geom_text(data = data.frame(type_type = "MDD SNP Set : AltTrait Coef",
                              group_type = "cr95set"),aes(5.1,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_cr95set_prs"],
                                                          label = "MDD cr95set", vjust = -1), size = 6)+
  geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_logbf_over4_prs"],
             color = "darkblue", lty = 2)+
  geom_text(data = data.frame(type_type = "MDD SNP Set : AltTrait Coef",
                              group_type = "cr95set"),aes(5.1,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_logbf_over4_prs"],
                                                          label = "MDD logbf4", vjust = -1), size = 6)+
  geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_hapmap_prs"],
             color = "gold4", lty = 2)+
  geom_text(data = data.frame(type_type = "MDD SNP Set : AltTrait Coef",
                              group_type = "cr95set"),aes(5.1,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_hapmap_prs"],
                                                          label = "MDD hapmap", vjust = -1), size = 6)+
  geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_twassig_prs"],
             color = "darkgreen", lty = 2)+
  geom_text(data = data.frame(type_type = "MDD SNP Set : AltTrait Coef",
                              group_type = "cr95set"),
            aes(5.1,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_twassig_prs"],
                label = "MDD twassig", vjust = -1), size = 6)+
  scale_fill_manual(values = c("pink1", "dodgerblue3"))+theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),

        axis.ticks.x = element_blank(),

        text = element_text(size = 20, color = "black"),

        axis.text = element_text(colour="black"),

        strip.text.x = element_text(size = 22, colour = "black"),

        legend.key.size = unit(1.5, "cm"))+
  geom_text(
    data = performance_stats[
      (performance_stats$class %in% c( "AltTraitSet-MDDCoef" , "MDDSet-AltTraitCoef") )
      &
        performance_stats$group_type != "twassig", ],
    aes(x = trait, y = max_auc*100, label = round(mean_auc*100, 1)),
    position = position_dodge(width = 0.9), # aligns with boxplot dodging
    vjust = -1,
    size = 6
  ) + labs(y = "Difference in AUC (%)", x = "Disorder", fill = "Set", pattern = "Class")+
  scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.00, 0.1)))+
  scale_pattern_manual(values = c(
    "AltTraitSet-MDDCoef" = "circle",
    "MDDSet-AltTraitCoef" = "crosshatch"
  ))+
  geom_blank(data = data.frame(trait2 = "EXTRA", type_type =  "MDD SNP Set : AltTrait Coef",
                               diff_AUC = 0, group_type = "cr95set"),
             aes(x = trait2))+
  scale_x_discrete(labels =
                     c("ADHD",
                       "ANX",
                       "BIP",
                       "SCZ",
                       ""
                     ),
                   expand = expansion(mult = c(0.05, 0.18)) )+
  guides(
    fill = guide_legend(override.aes = list(pattern = "none")),
    pattern = guide_legend(override.aes = list(fill = "white",
                                               color = "black"))

  )
dev.off()


ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/fig_v2/AltTraitSNPSet_MDDCoef_viceversa_ALLEUR_barplot_wpattern_small.png"),


       performance %>% filter(class == "AltTraitSet-MDDCoef" | class == "MDDSet-AltTraitCoef") %>%
         filter(group_type != "twassig") %>%
         ggplot( aes(y = diff_AUC*100, fill = group_type, x = trait2))+geom_bar_pattern(stat = "identity",position = "dodge",
                                                                                        aes(pattern = class)) + facet_wrap(~type_type,
                                                                                                                           scales = "free_x")+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_cr95set_prs"],
                    color = "maroon3", lty = 2)+
         geom_text(data = data.frame(type_type = "MDD SNP Set : AltTrait Coef",
                                     group_type = "cr95set"),aes(5.1,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_cr95set_prs"],
                                                                 label = "MDD cr95set", vjust = -1), size = 6)+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_logbf_over4_prs"],
                    color = "darkblue", lty = 2)+
         geom_text(data = data.frame(type_type = "MDD SNP Set : AltTrait Coef",
                                     group_type = "cr95set"),aes(5.1,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_logbf_over4_prs"],
                                                                 label = "MDD logbf4", vjust = -1), size = 6)+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_hapmap_prs"],
                    color = "gold4", lty = 2)+
         geom_text(data = data.frame(type_type = "MDD SNP Set : AltTrait Coef",
                                     group_type = "cr95set"),aes(5.1,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_hapmap_prs"],
                                                                 label = "MDD hapmap", vjust = -1), size = 6)+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_twassig_prs"],
                    color = "darkgreen", lty = 2)+
         geom_text(data = data.frame(type_type = "MDD SNP Set : AltTrait Coef",
                                     group_type = "cr95set"),
                   aes(5.1,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_twassig_prs"],
                       label = "MDD twassig", vjust = -1), size = 6)+
         scale_fill_manual(values = c("pink1", "dodgerblue3"))+theme_bw()+
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),

               axis.ticks.x = element_blank(),

               text = element_text(size = 20, color = "black"),

               axis.text = element_text(colour="black"),

               strip.text.x = element_text(size = 22, colour = "black"),

               legend.key.size = unit(1.5, "cm"))+
         geom_text(
           data = performance_stats[
             (performance_stats$class %in% c( "AltTraitSet-MDDCoef" , "MDDSet-AltTraitCoef") )
             &
               performance_stats$group_type != "twassig", ],
           aes(x = trait, y = max_auc*100, label = round(mean_auc*100, 1)),
           position = position_dodge(width = 0.9), # aligns with boxplot dodging
           vjust = -1,
           size = 6
         ) + labs(y = "Difference in AUC (%)", x = "Disorder", fill = "Set", pattern = "Class")+
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.00, 0.1)))+
         scale_pattern_manual(values = c(
           "AltTraitSet-MDDCoef" = "circle",
           "MDDSet-AltTraitCoef" = "crosshatch"
         ))+
         geom_blank(data = data.frame(trait2 = "EXTRA", type_type =  "MDD SNP Set : AltTrait Coef",
                                      diff_AUC = 0, group_type = "cr95set"),
                    aes(x = trait2))+
         scale_x_discrete(labels =
                            c("ADHD",
                              "ANX",
                              "BIP",
                              "SCZ",
                              ""
                            ),
                          expand = expansion(mult = c(0.05, 0.18)) )+
         guides(
           fill = guide_legend(override.aes = list(pattern = "none")),
           pattern = guide_legend(override.aes = list(fill = "white",
                                                      color = "black"))

         ),
       height = 8, width = 8)


ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/fig_v2/TraitSNPSet_TraitCoef_AllEUR_barplot_nopattern.png"),
       performance %>% filter(class == "AltTraitSet-AltTraitCoef" | class == "MDDSet-MDDCoef") %>%
         filter(group_type != "twassig") %>%
         ggplot( aes(y = diff_AUC*100, fill = group_type, x = coef))+geom_bar(stat = "identity",position = "dodge"
         ) + facet_wrap(~type_type)+
         scale_pattern_manual(values = c( "none","stripe"), labels = c( "MDD","Alt Trait"))+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_cr95set_prs"],
                    color =  "maroon3", lty = 2)+
         geom_text(aes(5.77,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_cr95set_prs"],
                       label = "MDD cr95set", vjust = -1), size = 6)+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_logbf_over4_prs"],
                    color = "darkblue", lty = 2)+
         geom_text(aes(5.77,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_logbf_over4_prs"],
                       label = "MDD logbf4", vjust = -1), size = 6)+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_hapmap_prs"],
                    color = "gold4", lty = 2)+
         geom_text(aes(5.77,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_hapmap_prs"],
                       label = "MDD hapmap", vjust = -1), size = 6)+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_twassig_prs"],
                    color = "darkgreen", lty = 2)+
         geom_text(aes(5.77,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_twassig_prs"],
                       label = "MDD twassig", vjust = -1), size = 6)+
         scale_fill_manual(values = c("pink1", "gold2", "dodgerblue3"))+theme_bw()+
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),

               axis.ticks.x = element_blank(),

               text = element_text(size = 20, color = "black"),

               axis.text = element_text(colour="black"),

               strip.text.x = element_text(size = 22, colour = "black"),

               legend.key.size = unit(1.5, "cm"))+
         geom_text(
           data = performance_stats[
             (performance_stats$class %in% c("AltTraitSet-AltTraitCoef" ,"MDDSet-MDDCoef") )
             &
               performance_stats$group_type != "twassig", ],
           aes(x = coef, y = max_auc*100, label = round(mean_auc*100, 1)),
           position = position_dodge(width = 0.9), # aligns with boxplot dodging
           vjust = -1,
           size = 6
         ) + labs(y = "Difference in AUC (%)", x = "Disorder", fill = "Set", pattern = "Class")+
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.00, 0.1)))+expand_limits(x=6.1),
       height = 8, width = 18)

ggsave(mdd_path("MDD_PROJECT", "/08.Generate_Figures/fig_v2/AltTraitSNPSet_MDDCoef_viceversa_ALLEUR_barplot_nopattern.png"),

       performance %>% filter(class == "AltTraitSet-MDDCoef" | class == "MDDSet-AltTraitCoef") %>%
         filter(group_type != "twassig") %>%
         ggplot( aes(y = diff_AUC*100, fill = group_type, x = trait2))+geom_bar(stat = "identity",position = "dodge") + facet_wrap(~type_type,
                                                                                                                                   scales = "free_x")+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_cr95set_prs"],
                    color = "maroon3", lty = 2)+
         geom_text(data = data.frame(type_type = "MDD SNP Set : AltTrait Coef",
                                     group_type = "cr95set"),aes(5.1,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_cr95set_prs"],
                                                                 label = "MDD cr95set", vjust = -1), size = 6)+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_logbf_over4_prs"],
                    color = "darkblue", lty = 2)+
         geom_text(data = data.frame(type_type = "MDD SNP Set : AltTrait Coef",
                                     group_type = "cr95set"),aes(5.1,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_logbf_over4_prs"],
                                                                 label = "MDD logbf4", vjust = -1), size = 6)+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_hapmap_prs"],
                    color = "gold4", lty = 2)+
         geom_text(data = data.frame(type_type = "MDD SNP Set : AltTrait Coef",
                                     group_type = "cr95set"),aes(5.1,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_hapmap_prs"],
                                                                 label = "MDD hapmap", vjust = -1), size = 6)+
         geom_hline(yintercept = 100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_twassig_prs"],
                    color = "darkgreen", lty = 2)+
         geom_text(data = data.frame(type_type = "MDD SNP Set : AltTrait Coef",
                                     group_type = "cr95set"),
                   aes(5.1,100*performance_stats$mean_auc[performance_stats$type == "MDDSet-MDDCoef_twassig_prs"],
                       label = "MDD twassig", vjust = -1), size = 6)+
         scale_fill_manual(values = c("pink1",  "dodgerblue3"))+theme_bw()+
         theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),

               axis.ticks.x = element_blank(),

               text = element_text(size = 20, color = "black"),

               axis.text = element_text(colour="black"),

               strip.text.x = element_text(size = 22, colour = "black"),

               legend.key.size = unit(1.5, "cm"))+
         geom_text(
           data = performance_stats[
             (performance_stats$class %in% c( "AltTraitSet-MDDCoef" , "MDDSet-AltTraitCoef") )
             &
               performance_stats$group_type != "twassig", ],
           aes(x = trait, y = max_auc*100, label = round(mean_auc*100, 1)),
           position = position_dodge(width = 0.9), # aligns with boxplot dodging
           vjust = -1,
           size = 6
         ) + labs(y = "Difference in AUC (%)", x = "Disorder", fill = "Set", pattern = "Class")+
         scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.00, 0.1)))+
         scale_pattern_manual(values = c( "circle", "crosshatch"))+
         geom_blank(data = data.frame(trait2 = "EXTRA", type_type =  "MDD SNP Set : AltTrait Coef",
                                      diff_AUC = 0, group_type = "cr95set"),
                    aes(x = trait2))+
         scale_x_discrete(labels =
                            c("ADHD",
                              "ANX",
                              "BIP",
                              "SCZ",
                              ""
                            ),
                          expand = expansion(mult = c(0.05, 0.18)) ),
       height = 8, width = 19)
