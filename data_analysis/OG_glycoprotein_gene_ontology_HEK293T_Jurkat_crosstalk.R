#import packages
packages_names <- c("tidyverse", "rstatix", "ggpubr", "showtext", "readxl")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
#general from HEK293T down to Jurkat up 
downHEK293T_medianJurkat <- semi_join(OG_glycoprotein_down_HEK293T, OG_glycoprotein_median_Jurkat, by = "UniprotID")

downHEK293T_upJurkat <- semi_join(OG_glycoprotein_down_HEK293T, OG_glycoprotein_up_Jurkat, by = "UniprotID")

medianHEK293T_upJurkat <- semi_join(OG_glycoprotein_median_HEK293T, OG_glycoprotein_up_Jurkat, by = "UniprotID")

general_from_HEK293T_down_to_Jurkat_up <- bind_rows(
  downHEK293T_medianJurkat, downHEK293T_upJurkat, medianHEK293T_upJurkat
)

write_xlsx(general_from_HEK293T_down_to_Jurkat_up, path = paste0(file_path, "general_from_HEK293T_down_to_Jurkat_up.xlsx"))

#general from HEK293T up to Jurkat down
upHEK293T_medianJurkat <- semi_join(OG_glycoprotein_up_HEK293T, OG_glycoprotein_median_Jurkat, by = "UniprotID")

upHEK293T_downJurkat <- semi_join(OG_glycoprotein_up_HEK293T, OG_glycoprotein_down_Jurkat, by = "UniprotID")

medianHEK293T_downJurkat <- semi_join(OG_glycoprotein_median_HEK293T, OG_glycoprotein_down_Jurkat, by = "UniprotID")

general_from_HEK293T_up_Jurkat_down <- bind_rows(
  upHEK293T_medianJurkat, upHEK293T_downJurkat, medianHEK293T_downJurkat
)

write_xlsx(general_from_HEK293T_up_Jurkat_down, path = paste0(file_path, "general_from_HEK293T_up_Jurkat_down.xlsx"))

#import data
gene_ontology_HEK293T_Jurkat_crosstalk_from_HEK293T_up_to_Jurkat_down <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_HEK293T_Jurkat_crosstalk.xlsx",
  sheet = "from_HEK293T_up_to_Jurkat_down",
  col_names = TRUE,
  .name_repair = "universal"
)
gene_ontology_HEK293T_Jurkat_crosstalk_from_HEK293T_down_to_Jurkat_up <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_HEK293T_Jurkat_crosstalk.xlsx",
  sheet = "from_HEK293T_down_to_Jurkat_up",
  col_names = TRUE,
  .name_repair = "universal"
)

#plot from_HEK293T_up_to_Jurkat_down
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

gene_ontology_HEK293T_Jurkat_crosstalk_from_HEK293T_up_to_Jurkat_down_plot <- gene_ontology_HEK293T_Jurkat_crosstalk_from_HEK293T_up_to_Jurkat_down |> 
  ggplot(aes(x = Count, y = fct_reorder(Term, Count), size = Fold.Enrichment, fill = -log10(P.Value))) +
  geom_hline(aes(yintercept = Term), color = "gray", linetype = "dashed") +
  geom_point(alpha = 0.5, pch = 21) +
  labs(x = "Count", y = "", size = "Fold Enrichment", fill = expression(-log[10](paste(italic(P), " Value")))) +
  scale_size(range = c(8, 12), breaks = c(2, 7, 12), guide = guide_legend(order = 2)) +
  scale_fill_gradient(low = alpha(Color_2, 0.1), high = Color_2) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 100),
    axis.title.y = element_text(size = 100),
    axis.text.x = element_text(color = "black", size = 100),
    axis.text.y = element_text(color = "black", size = 100, lineheight = 0.12),
    legend.title = element_text(size = 100, color = "black"),
    legend.text = element_text(size = 100, color = "black"),
    legend.spacing.x = unit(0.1, units = "in"),
    legend.spacing.y = unit(0.1, units = "in")
  )

ggsave(filename = paste0(file_path, "gene_ontology_HEK293T_Jurkat_crosstalk_from_HEK293T_up_to_Jurkat_down_plot.png"), 
       plot = gene_ontology_HEK293T_Jurkat_crosstalk_from_HEK293T_up_to_Jurkat_down_plot, height = 4, width = 8, units = c("in"), dpi = 600)

#plot from_HEK293T_down_to_Jurkat_up
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

gene_ontology_HEK293T_Jurkat_crosstalk_from_HEK293T_down_to_Jurkat_up_plot <- gene_ontology_HEK293T_Jurkat_crosstalk_from_HEK293T_down_to_Jurkat_up |> 
  ggplot(aes(x = Count, y = fct_reorder(Term, Count), size = Fold.Enrichment, fill = -log10(P.Value))) +
  geom_hline(aes(yintercept = Term), color = "gray", linetype = "dashed") +
  geom_point(alpha = 0.5, pch = 21) +
  labs(x = "Count", y = "", size = "Fold Enrichment", fill = expression(-log[10](paste(italic(P), " Value")))) +
  scale_size(range = c(8, 12), breaks = c(2, 4, 6), guide = guide_legend(order = 2)) +
  scale_fill_gradient(low = alpha(Color_4, 0.1), high = Color_4) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 100),
    axis.title.y = element_text(size = 100),
    axis.text.x = element_text(color = "black", size = 100),
    axis.text.y = element_text(color = "black", size = 100, lineheight = 0.12),
    legend.title = element_text(size = 100, color = "black"),
    legend.text = element_text(size = 100, color = "black"),
    legend.spacing.x = unit(0.1, units = "in"),
    legend.spacing.y = unit(0.1, units = "in")
  )

ggsave(filename = paste0(file_path, "gene_ontology_HEK293T_Jurkat_crosstalk_from_HEK293T_down_to_Jurkat_up_plot.png"), plot = gene_ontology_HEK293T_Jurkat_crosstalk_from_HEK293T_down_to_Jurkat_up_plot, height = 4, width = 8, units = c("in"), dpi = 600)
