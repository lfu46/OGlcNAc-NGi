#import packages
packages_names <- c("tidyverse", "rstatix", "ggpubr", "showtext", "readxl", "ggbeeswarm")
lapply(packages_names, require, character.only = TRUE)

#metabolism
#import data
gene_ontology_OG_glycoprotein_pathway_metabolism_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_pathway_metabolism.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "HepG2"
) |> select(UniprotID = UNIPROT_ACCESSION) |> mutate(cell = "HepG2")

gene_ontology_OG_glycoprotein_pathway_metabolism_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_pathway_metabolism.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "HEK293T"
) |> select(UniprotID = UNIPROT_ACCESSION) |> mutate(cell = "HEK293T")

gene_ontology_OG_glycoprotein_pathway_metabolism_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_pathway_metabolism.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "Jurkat"
) |> select(UniprotID = UNIPROT_ACCESSION) |> mutate(cell = "Jurkat")

#generate data frame
gene_ontology_OG_glycoprotein_pathway_metabolism_logFC_HepG2 <- gene_ontology_OG_glycoprotein_pathway_metabolism_HepG2 |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  select(UniprotID, cell, logFC)

gene_ontology_OG_glycoprotein_pathway_metabolism_logFC_HEK293T <- gene_ontology_OG_glycoprotein_pathway_metabolism_HEK293T |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  select(UniprotID, cell, logFC)

gene_ontology_OG_glycoprotein_pathway_metabolism_logFC_Jurkat <- gene_ontology_OG_glycoprotein_pathway_metabolism_Jurkat |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  select(UniprotID, cell, logFC)

OG_glycoprotein_pathway_metabolism_combined <- bind_rows(
  gene_ontology_OG_glycoprotein_pathway_metabolism_logFC_HepG2,
  gene_ontology_OG_glycoprotein_pathway_metabolism_logFC_HEK293T,
  gene_ontology_OG_glycoprotein_pathway_metabolism_logFC_Jurkat
)

#wilcox test
OG_glycoprotein_pathway_metabolism_wilcox_test <- OG_glycoprotein_pathway_metabolism_combined |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#swarm plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

swarm_plot_OG_glycoprotein_pathway_metabolism <- OG_glycoprotein_pathway_metabolism_combined |> 
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
  stat_pvalue_manual(data = OG_glycoprotein_pathway_metabolism_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 40, 
                     hide.ns = "p", y.position = c(1.8)) +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)"), fill = "") +
  ggtitle("Metabolism") +
  theme_bw() +
  theme(
    title = element_text(size = 100, color = "black"),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 80, color = "black"),
    axis.line.x = element_blank(),
    axis.ticks.length.x = unit(0, "in"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "swarm_plot_OG_glycoprotein_pathway_metabolism.png"), plot = swarm_plot_OG_glycoprotein_pathway_metabolism,
       height = 3, width = 3, units = "in", dpi = 600)
