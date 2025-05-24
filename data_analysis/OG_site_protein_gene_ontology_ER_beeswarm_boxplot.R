#import packages
packages_names <- c("tidyverse", "rstatix", "ggpubr", "showtext", "readxl", "ggbeeswarm")
lapply(packages_names, require, character.only = TRUE)

#import data
#HEK293T
gene_ontology_OG_site_protein_ER_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_site_protein_ER.xlsx",
  sheet = 'HEK293T',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_site_protein_ER_HEK293T <- OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(UniprotID %in% gene_ontology_OG_site_protein_ER_HEK293T$UNIPROT_ACCESSION) |> 
  select(Index, logFC) |> mutate(cell = "HEK293T")

#HepG2
gene_ontology_OG_site_protein_ER_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_site_protein_ER.xlsx",
  sheet = 'HepG2',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_site_protein_ER_HepG2 <- OG_site_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(UniprotID %in% gene_ontology_OG_site_protein_ER_HepG2$UNIPROT_ACCESSION) |> 
  select(Index, logFC) |> mutate(cell = "HepG2")

#Jurkat
gene_ontology_OG_site_protein_ER_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_site_protein_ER.xlsx",
  sheet = 'Jurkat',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_site_protein_ER_Jurkat <- OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(UniprotID %in% gene_ontology_OG_site_protein_ER_Jurkat$UNIPROT_ACCESSION) |> 
  select(Index, logFC) |> mutate(cell = "Jurkat")

#combine
OG_site_protein_ER_combined <- bind_rows(
  OG_site_protein_ER_HEK293T,
  OG_site_protein_ER_HepG2,
  OG_site_protein_ER_Jurkat
)

#wilcox test
OG_site_protein_ER_wilcox_test <- OG_site_protein_ER_combined |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#beeswarm boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

beeswarm_boxplot_OG_site_protein_ER <- OG_site_protein_ER_combined |> 
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
  stat_pvalue_manual(data = OG_site_protein_ER_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 50, 
                     hide.ns = "p", y.position = c(1.6, 1.9)) +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "", y = expression(Δlog[2]*"FC"), fill = "") +
  ggtitle("Endoplasmic reticulum") +
  theme_bw() +
  theme(
    title = element_text(size = 80, color = "black"),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 80, color = "black"),
    axis.line.x = element_blank(),
    axis.ticks.length.x = unit(0, "in"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "beeswarm_boxplot_OG_site_protein_ER.png"), 
       plot = beeswarm_boxplot_OG_site_protein_ER,
       height = 3, width = 3, units = "in", dpi = 600)

