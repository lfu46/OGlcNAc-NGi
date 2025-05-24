#import packages
packages_names <- c("tidyverse", "rstatix", "ggpubr", "showtext", "readxl", "ggbeeswarm")
lapply(packages_names, require, character.only = TRUE)

#import data
#HEK293T
gene_ontology_OG_glycoprotein_protein_exosome_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein_exosome.xlsx",
  sheet = 'HEK293T',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_glycoprotein_protein_exosome_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(UniprotID %in% gene_ontology_OG_glycoprotein_protein_exosome_HEK293T$UNIPROT_ACCESSION) |> 
  select(UniprotID, logFC) |> mutate(cell = "HEK293T")

#HepG2
gene_ontology_OG_glycoprotein_protein_exosome_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein_exosome.xlsx",
  sheet = 'HepG2',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_glycoprotein_protein_exosome_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(UniprotID %in% gene_ontology_OG_glycoprotein_protein_exosome_HepG2$UNIPROT_ACCESSION) |> 
  select(UniprotID, logFC) |> mutate(cell = "HepG2")

#Jurkat
gene_ontology_OG_glycoprotein_protein_exosome_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein_exosome.xlsx",
  sheet = 'Jurkat',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_glycoprotein_protein_exosome_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(UniprotID %in% gene_ontology_OG_glycoprotein_protein_exosome_Jurkat$UNIPROT_ACCESSION) |> 
  select(UniprotID, logFC) |> mutate(cell = "Jurkat")

#combine
OG_glycoprotein_protein_exosome_combined <- bind_rows(
  OG_glycoprotein_protein_exosome_HEK293T,
  OG_glycoprotein_protein_exosome_HepG2,
  OG_glycoprotein_protein_exosome_Jurkat
)

#wilcox test
OG_glycoprotein_protein_exosome_wilcox_test <- OG_glycoprotein_protein_exosome_combined |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#beeswarm boxplot
font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

beeswarm_boxplot_OG_glycoprotein_protein_exosome <- OG_glycoprotein_protein_exosome_combined |> 
  ggplot(aes(x = cell, y = logFC)) +
  geom_beeswarm(aes(color = cell)) +
  geom_boxplot(width = 0.2, outliers = FALSE) +
  scale_color_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  stat_pvalue_manual(data = OG_glycoprotein_protein_exosome_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 5, 
                     hide.ns = "p", y.position = c(1.6)) +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "", y = expression(Δlog[2]*"FC"), fill = "") +
  ggtitle("Exosome") +
  theme_bw() +
  theme(
    title = element_text(size = 10, color = "black"),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 10),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.line.x = element_blank(),
    axis.ticks.length.x = unit(0, "in"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "beeswarm_boxplot_OG_glycoprotein_protein_exosome.eps"), 
       device = "eps",
       plot = beeswarm_boxplot_OG_glycoprotein_protein_exosome,
       height = 1.8, width = 1.8, 
       units = "in", 
       dpi = 1200)
