#import packages
packages_names <- c("UpSetR", "ggupset", "ggalluvial", "showtext")
lapply(packages_names, require, character.only = TRUE)

#Alluvial plot
#generate data frame
#HepG2
OG_glycoprotein_up_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> select(UniprotID) |> 
  mutate(category = "up" , group = "Raw", count = 1)

OG_glycoprotein_down_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> select(UniprotID) |> 
  mutate(category = "down", group = "Raw", count = 1)

OG_glycoprotein_median_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(! UniprotID %in% OG_glycoprotein_up_HepG2$UniprotID) |> 
  filter(! UniprotID %in% OG_glycoprotein_down_HepG2$UniprotID) |> 
  select(UniprotID) |> 
  mutate(category = "median", group = "Raw", count = 1)

OG_glycoprotein_protein_up_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> select(UniprotID) |> 
  mutate(category = "up", group = "Normalized", count = 1)

OG_glycoprotein_protein_down_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> select(UniprotID) |> 
  mutate(category = "down", group = "Normalized", count = 1)

OG_glycoprotein_protein_median_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |>
  filter(! UniprotID %in% OG_glycoprotein_protein_up_HepG2$UniprotID) |> 
  filter(! UniprotID %in% OG_glycoprotein_protein_down_HepG2$UniprotID) |> 
  select(UniprotID) |> 
  mutate(category = "median", group = "Normalized", count = 1)

OG_glycoprotein_protein_missing_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(!UniprotID %in% OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2$UniprotID) |> 
  select(UniprotID) |> 
  mutate(category = "missing", group = "Normalized", count = 1)

#alluvial diagram
#combine
OG_combine_alluvial <- bind_rows(
  OG_glycoprotein_up_HepG2,
  OG_glycoprotein_down_HepG2,
  OG_glycoprotein_median_HepG2,
  OG_glycoprotein_protein_up_HepG2,
  OG_glycoprotein_protein_down_HepG2,
  OG_glycoprotein_protein_median_HepG2,
  OG_glycoprotein_protein_missing_HepG2
) |> mutate(category = factor(category, levels = c("missing" ,"up", "median", "down")))

font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

alluvial_plot_OG_combine_HepG2 <- OG_combine_alluvial |> 
  ggplot(aes(x = factor(group, levels = c("Raw", "Normalized")), stratum = category, alluvium = UniprotID, y = count,
             fill = category, label = category)) +
  geom_flow(width = 0.5) +
  geom_stratum(width = 0.5) +
  geom_text(aes(color = category), stat = "stratum", size = 25) +
  labs(x = "", y = "# of glycoprotein") +
  scale_color_manual(
    values = c(
      "up" = "black",
      "median" = "black",
      "missing" = "black",
      "down" = "white"
    )
  ) +
  scale_fill_manual(
    values = c(
      "missing" = "white",
      "up" = "orange",
      "down" = Color_6,
      "median" = "gray"
    )
  ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.position = 'none'
  )

ggsave(filename = paste0(file_path, "alluvial_plot_OG_combine_HepG2.png"), 
       plot = alluvial_plot_OG_combine_HepG2, 
       height = 3, width = 4, units = "in", dpi = 600)

#Dot plot
#import data
OG_site_protein_up_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein.xlsx",
  sheet = 'up_HEK293T',
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(Cell = "HEK293T")

OG_site_protein_up_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein.xlsx",
  sheet = 'up_HepG2',
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(Cell = "HepG2")

OG_site_protein_up_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein.xlsx",
  sheet = 'up_Jurkat',
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(Cell = "Jurkat")

OG_site_protein_down_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein.xlsx",
  sheet = 'down_HEK293T',
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(Cell = "HEK293T")

OG_site_protein_down_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein.xlsx",
  sheet = 'down_HepG2',
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(Cell = "HepG2")

OG_site_protein_down_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein.xlsx",
  sheet = 'down_Jurkat',
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(Cell = "Jurkat")

#generate data frame
#up
OG_site_protein_up_combined <- bind_rows(
  OG_site_protein_up_HEK293T,
  OG_site_protein_up_HepG2,
  OG_site_protein_up_Jurkat
) |> select(Term, Count, P.Value, Cell)

font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_OG_site_protein_up <- OG_site_protein_up_combined |> 
  ggplot(aes(x = Cell, y = factor(Term, levels = Term))) +
  geom_point(aes(fill = Cell, size = Count, alpha = P.Value), shape = 21) +
  scale_x_discrete(labels = c("HEK293T", "HepG2", "Jurkat")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_size(range = c(6, 10)) +
  scale_alpha_binned(range = c(1, 0.2)) +
  scale_fill_manual(
    name = "Cell",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.2),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black", lineheight = 0.1),
    legend.text = element_text(size = 80),
    legend.title = element_text(size = 80, color = "black"),
    legend.frame = element_rect(color = "black", linewidth = 0.2),
    legend.ticks = element_line(color = "black", linewidth = 0.2)
  )

ggsave(filename = paste0(file_path, "dot_plot_OG_site_protein_up.png"), 
       plot = dot_plot_OG_site_protein_up, 
       height = 7, width = 6, units = c("in"), dpi = 600)

#down
OG_site_protein_down_combined <- bind_rows(
  OG_site_protein_down_HEK293T,
  OG_site_protein_down_HepG2,
  OG_site_protein_down_Jurkat
) |> select(Term, Count, P.Value, Cell)

font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_OG_site_protein_down <- OG_site_protein_down_combined |> 
  ggplot(aes(x = Cell, y = factor(Term, levels = OG_site_protein_down_combined |> distinct(Term) |> pull()))) +
  geom_point(aes(fill = Cell, size = Count, alpha = P.Value), shape = 21) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_size(range = c(6, 10)) +
  scale_alpha_binned(range = c(1, 0.2), breaks = c(0.001, 0.005, 0.1)) +
  scale_fill_manual(
    name = "Cell",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.2),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black", lineheight = 0.1),
    legend.text = element_text(size = 80),
    legend.title = element_text(size = 80, color = "black"),
    legend.frame = element_rect(color = "black", linewidth = 0.2),
    legend.ticks = element_line(color = "black", linewidth = 0.2)
  )

ggsave(filename = paste0(file_path, "dot_plot_OG_site_protein_down.png"), 
       plot = dot_plot_OG_site_protein_down, 
       height = 7, width = 6, units = c("in"), dpi = 600)

#Selected terms
#Mitochondrion
#import data
#HEK293T
gene_ontology_OG_glycoprotein_protein_mito_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein_mito.xlsx",
  sheet = 'HEK293T',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_glycoprotein_protein_mito_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(UniprotID %in% gene_ontology_OG_glycoprotein_protein_mito_HEK293T$UNIPROT_ACCESSION) |> 
  select(UniprotID, logFC) |> mutate(cell = "HEK293T")

#HepG2
gene_ontology_OG_glycoprotein_protein_mito_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein_mito.xlsx",
  sheet = 'HepG2',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_glycoprotein_protein_mito_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(UniprotID %in% gene_ontology_OG_glycoprotein_protein_mito_HepG2$UNIPROT_ACCESSION) |> 
  select(UniprotID, logFC) |> mutate(cell = "HepG2")

#Jurkat
gene_ontology_OG_glycoprotein_protein_mito_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein_mito.xlsx",
  sheet = 'Jurkat',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_glycoprotein_protein_mito_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(UniprotID %in% gene_ontology_OG_glycoprotein_protein_mito_Jurkat$UNIPROT_ACCESSION) |> 
  select(UniprotID, logFC) |> mutate(cell = "Jurkat")

#combine
OG_glycoprotein_protein_mito_combined <- bind_rows(
  OG_glycoprotein_protein_mito_HEK293T,
  OG_glycoprotein_protein_mito_HepG2,
  OG_glycoprotein_protein_mito_Jurkat
)

#wilcox test
OG_glycoprotein_protein_mito_wilcox_test <- OG_glycoprotein_protein_mito_combined |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#beeswarm boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

beeswarm_boxplot_OG_glycoprotein_protein_mito <- OG_glycoprotein_protein_mito_combined |> 
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
  stat_pvalue_manual(data = OG_glycoprotein_protein_mito_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 50, 
                     hide.ns = "p", y.position = c(1.6)) +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "", y = expression(Δlog[2]*"FC"), fill = "") +
  ggtitle("Mitochondrion") +
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

ggsave(filename = paste0(file_path, "beeswarm_boxplot_OG_glycoprotein_protein_mito.png"), 
       plot = beeswarm_boxplot_OG_glycoprotein_protein_mito,
       height = 3, width = 3, units = "in", dpi = 600)

#Extracellular exosome
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
font_add(family = "calibri", regular = "calibri.ttf")
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
                     tip.length = 0, size = 50, 
                     hide.ns = "p", y.position = c(1.6)) +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "", y = expression(Δlog[2]*"FC"), fill = "") +
  ggtitle("Extracellular exosome") +
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

ggsave(filename = paste0(file_path, "beeswarm_boxplot_OG_glycoprotein_protein_exosome.png"), 
       plot = beeswarm_boxplot_OG_glycoprotein_protein_exosome,
       height = 3, width = 3, units = "in", dpi = 600)

#Metabolic pathway
#import data
#HEK293T
gene_ontology_OG_glycoprotein_protein_metabolic_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein_metabolic.xlsx",
  sheet = 'HEK293T',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_glycoprotein_protein_metabolic_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(UniprotID %in% gene_ontology_OG_glycoprotein_protein_metabolic_HEK293T$UNIPROT_ACCESSION) |> 
  select(UniprotID, logFC) |> mutate(cell = "HEK293T")

#HepG2
gene_ontology_OG_glycoprotein_protein_metabolic_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein_metabolic.xlsx",
  sheet = 'HepG2',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_glycoprotein_protein_metabolic_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(UniprotID %in% gene_ontology_OG_glycoprotein_protein_metabolic_HepG2$UNIPROT_ACCESSION) |> 
  select(UniprotID, logFC) |> mutate(cell = "HepG2")

#Jurkat
gene_ontology_OG_glycoprotein_protein_metabolic_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein_metabolic.xlsx",
  sheet = 'Jurkat',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_glycoprotein_protein_metabolic_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(UniprotID %in% gene_ontology_OG_glycoprotein_protein_metabolic_Jurkat$UNIPROT_ACCESSION) |> 
  select(UniprotID, logFC) |> mutate(cell = "Jurkat")

#combine
OG_glycoprotein_protein_metabolic_combined <- bind_rows(
  OG_glycoprotein_protein_metabolic_HEK293T,
  OG_glycoprotein_protein_metabolic_HepG2,
  OG_glycoprotein_protein_metabolic_Jurkat
)

#wilcox test
OG_glycoprotein_protein_metabolic_wilcox_test <- OG_glycoprotein_protein_metabolic_combined |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#beeswarm boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

beeswarm_boxplot_OG_glycoprotein_protein_metabolic <- OG_glycoprotein_protein_metabolic_combined |> 
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
  stat_pvalue_manual(data = OG_glycoprotein_protein_metabolic_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 50, 
                     hide.ns = "p", y.position = c(1.6)) +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "", y = expression(Δlog[2]*"FC"), fill = "") +
  ggtitle("Metabolic pathway") +
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

ggsave(filename = paste0(file_path, "beeswarm_boxplot_OG_glycoprotein_protein_metabolic.png"), 
       plot = beeswarm_boxplot_OG_glycoprotein_protein_metabolic,
       height = 3, width = 3, units = "in", dpi = 600)

#Membrane
#import data
#HEK293T
gene_ontology_OG_glycoprotein_protein_membrane_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein_membrane.xlsx",
  sheet = 'HEK293T',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_glycoprotein_protein_membrane_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(UniprotID %in% gene_ontology_OG_glycoprotein_protein_membrane_HEK293T$UNIPROT_ACCESSION) |> 
  select(UniprotID, logFC) |> mutate(cell = "HEK293T")

#HepG2
gene_ontology_OG_glycoprotein_protein_membrane_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein_membrane.xlsx",
  sheet = 'HepG2',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_glycoprotein_protein_membrane_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(UniprotID %in% gene_ontology_OG_glycoprotein_protein_membrane_HepG2$UNIPROT_ACCESSION) |> 
  select(UniprotID, logFC) |> mutate(cell = "HepG2")

#Jurkat
gene_ontology_OG_glycoprotein_protein_membrane_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_protein_membrane.xlsx",
  sheet = 'Jurkat',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_glycoprotein_protein_membrane_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(UniprotID %in% gene_ontology_OG_glycoprotein_protein_membrane_Jurkat$UNIPROT_ACCESSION) |> 
  select(UniprotID, logFC) |> mutate(cell = "Jurkat")

#combine
OG_glycoprotein_protein_membrane_combined <- bind_rows(
  OG_glycoprotein_protein_membrane_HEK293T,
  OG_glycoprotein_protein_membrane_HepG2,
  OG_glycoprotein_protein_membrane_Jurkat
)

#wilcox test
OG_glycoprotein_protein_membrane_wilcox_test <- OG_glycoprotein_protein_membrane_combined |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#beeswarm boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

beeswarm_boxplot_OG_glycoprotein_protein_membrane <- OG_glycoprotein_protein_membrane_combined |> 
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
  stat_pvalue_manual(data = OG_glycoprotein_protein_membrane_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 50, 
                     hide.ns = "p", y.position = c(1.6, 1.9, 1.6)) +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "", y = expression(Δlog[2]*"FC"), fill = "") +
  ggtitle("Membrane") +
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

ggsave(filename = paste0(file_path, "beeswarm_boxplot_OG_glycoprotein_protein_membrane.png"), 
       plot = beeswarm_boxplot_OG_glycoprotein_protein_membrane,
       height = 3, width = 3, units = "in", dpi = 600)
