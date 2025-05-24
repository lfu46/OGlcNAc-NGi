#import packages
packages_names <- c("tidyverse", "writexl", "readxl", "ComplexHeatmap", "circlize", 
                    "ggpubr", "rstatix", "showtext", "VennDiagram")
lapply(packages_names, require, character.only = TRUE)

#import data
#HEK293T
OG_glycoprotein_exosome_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_Exosome_Nucleus_HEK293T.xlsx",
  sheet = 'Exosome',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  mutate(cell = "HEK293T", category = "Exosome")

OG_glycoprotein_nucleus_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_Exosome_Nucleus_HEK293T.xlsx",
  sheet = 'Nucleus',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  mutate(cell = "HEK293T", category = "Nucleus")

OG_glycoprotein_exosome_nucleus_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_Exosome_Nucleus_HEK293T.xlsx",
  sheet = 'Exosome_Nucleus',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  mutate(cell = "HEK293T", category = "Exosome, Nucleus")

OG_glycoprotein_nucleus_not_exosome_HEK293T <- OG_glycoprotein_nucleus_HEK293T |> 
  filter(!UniprotID %in% OG_glycoprotein_exosome_nucleus_HEK293T$UniprotID) |> 
  mutate(cell = "HEK293T", category = "Nucleus, not Exosome")

OG_glycoprotein_exosome_not_nucleus_HEK293T <- OG_glycoprotein_exosome_HEK293T |> 
  filter(!UniprotID %in% OG_glycoprotein_nucleus_HEK293T$UniprotID) |> 
  mutate(cell = "HEK293T", category = "Exosome, not Nucleus")

#combine
OG_glycoprotein_exosome_nucleus_HEK293T_combined <- bind_rows(
  OG_glycoprotein_exosome_HEK293T,
  OG_glycoprotein_nucleus_HEK293T,
  OG_glycoprotein_exosome_nucleus_HEK293T,
  OG_glycoprotein_nucleus_not_exosome_HEK293T,
  OG_glycoprotein_exosome_not_nucleus_HEK293T
)

#wilcox test
OG_glycoprotein_exosome_nucleus_HEK293T_combined |> 
  wilcox_test(logFC ~ category, p.adjust.method = "BH") |> 
  add_significance("p")

#import data
#HepG2
OG_glycoprotein_exosome_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_Exosome_Nucleus_HepG2.xlsx",
  sheet = 'Exosome',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  mutate(cell = "HepG2", category = "Exosome")

OG_glycoprotein_nucleus_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_Exosome_Nucleus_HepG2.xlsx",
  sheet = 'Nucleus',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  mutate(cell = "HepG2", category = "Nucleus")

OG_glycoprotein_exosome_nucleus_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_Exosome_Nucleus_HepG2.xlsx",
  sheet = 'Exosome_Nucleus',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  mutate(cell = "HepG2", category = "Exosome, Nucleus")

OG_glycoprotein_nucleus_not_exosome_HepG2 <- OG_glycoprotein_nucleus_HepG2 |> 
  filter(!UniprotID %in% OG_glycoprotein_exosome_nucleus_HepG2$UniprotID) |> 
  mutate(cell = "HepG2", category = "Nucleus, not Exosome")

OG_glycoprotein_exosome_not_nucleus_HepG2 <- OG_glycoprotein_exosome_HepG2 |> 
  filter(!UniprotID %in% OG_glycoprotein_nucleus_HepG2$UniprotID) |> 
  mutate(cell = "HepG2", category = "Exosome, not Nucleus")

#combine
OG_glycoprotein_exosome_nucleus_HepG2_combined <- bind_rows(
  OG_glycoprotein_exosome_HepG2,
  OG_glycoprotein_nucleus_HepG2,
  OG_glycoprotein_exosome_nucleus_HepG2,
  OG_glycoprotein_nucleus_not_exosome_HepG2,
  OG_glycoprotein_exosome_not_nucleus_HepG2
)

write_xlsx(OG_glycoprotein_exosome_nucleus_HepG2_combined, path = paste0(file_path, "OG_glycoprotein_exosome_nucleus_HepG2_combined.xlsx"))

#wilcox test
OG_glycoprotein_exosome_nucleus_HepG2_combined |> 
  wilcox_test(logFC ~ category, p.adjust.method = "BH") |> 
  add_significance("p")

#import data
#Jurkat
OG_glycoprotein_exosome_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_Exosome_Nucleus_Jurkat.xlsx",
  sheet = 'Exosome',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  mutate(cell = "Jurkat", category = "Exosome")

OG_glycoprotein_nucleus_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_Exosome_Nucleus_Jurkat.xlsx",
  sheet = 'Nucleus',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  mutate(cell = "Jurkat", category = "Nucleus")

OG_glycoprotein_exosome_nucleus_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_Exosome_Nucleus_Jurkat.xlsx",
  sheet = 'Exosome_Nucleus',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  mutate(cell = "Jurkat", category = "Exosome, Nucleus")

OG_glycoprotein_nucleus_not_exosome_Jurkat <- OG_glycoprotein_nucleus_Jurkat |> 
  filter(!UniprotID %in% OG_glycoprotein_exosome_nucleus_Jurkat$UniprotID) |> 
  mutate(cell = "Jurkat", category = "Nucleus, not Exosome")

OG_glycoprotein_exosome_not_nucleus_Jurkat <- OG_glycoprotein_exosome_Jurkat |> 
  filter(!UniprotID %in% OG_glycoprotein_nucleus_Jurkat$UniprotID) |> 
  mutate(cell = "Jurkat", category = "Exosome, not Nucleus")

#combine
OG_glycoprotein_exosome_nucleus_Jurkat_combined <- bind_rows(
  OG_glycoprotein_exosome_Jurkat,
  OG_glycoprotein_nucleus_Jurkat,
  OG_glycoprotein_exosome_nucleus_Jurkat,
  OG_glycoprotein_nucleus_not_exosome_Jurkat,
  OG_glycoprotein_exosome_not_nucleus_Jurkat
)

#wilcox test
OG_glycoprotein_exosome_nucleus_Jurkat_combined |> 
  wilcox_test(logFC ~ category, p.adjust.method = "BH") |> 
  add_significance("p")

#combine all
OG_glycoprotein_exosome_nucleus_combined <- bind_rows(
  #OG_glycoprotein_exosome_nucleus_HEK293T_combined,
  OG_glycoprotein_exosome_nucleus_HepG2_combined
  #OG_glycoprotein_exosome_nucleus_Jurkat_combined
)

#wilcox test
OG_glycoprotein_exosome_nucleus_wilcox_test <- OG_glycoprotein_exosome_nucleus_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ category, p.adjust.method = "BH") |> 
  add_significance("p") |> 
  filter(p.signif != "ns") |> 
  slice(c(1, 3, 5, 6, 7))

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_exosome_nucleus <- OG_glycoprotein_exosome_nucleus_combined |> 
  ggplot() +
  geom_violin(aes(x = factor(category, levels = c("Exosome", "Nucleus", "Exosome, Nucleus", "Nucleus, not Exosome", "Exosome, not Nucleus")), 
                  y = logFC, fill = cell), color = "transparent") +
  geom_boxplot(aes(x = factor(category, levels = c("Exosome", "Nucleus", "Exosome, Nucleus", "Nucleus, not Exosome", "Exosome, not Nucleus")), y = logFC), 
               color = "black", outliers = FALSE, width = 0.2) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  #facet_grid(~ cell, scales = "free_x") +
  stat_pvalue_manual(data = OG_glycoprotein_exosome_nucleus_wilcox_test, label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(2.1, 2.4, 2.0, 2.0, 1.8)) +
  coord_cartesian(ylim = c(-1.8, 2.5)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 100, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_exosome_nucleus.png"), 
       plot = violin_boxplot_OG_glycoprotein_exosome_nucleus, height = 4, width = 4, units = "in", dpi = 600)

#gene ontology
#import data
gene_ontology_exosome_nucleus_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\gene_ontology_exosome_nuclear_HepG2.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(Exp = "Exosome, Nucleus")

#dot plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_exosome_nucleus_HepG2 <- gene_ontology_exosome_nucleus_HepG2 |> 
  ggplot(aes(x = Exp, y = Term, size = Count, fill = P.Value)) +
  geom_point(shape = 21) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_size(range = c(5, 11), breaks = c(8, 12, 16)) +
  scale_fill_stepsn(limits = c(0, 1), 
                    breaks = c(0, 1E-3, 0.01, 0.05, 1), 
                    labels = c("0", "1E-3", "0.01", "0.05", "1"), 
                    n.breaks = 6, 
                    values = scales::rescale(c(0, 1E-3, 0.01, 0.05, 1)),
                    colours = c(Color_8, "transparent")) +
  scale_x_discrete(labels = c("Exosome,\nNucleus")) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.2),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 100, color = "black", lineheight = 0.15),
    axis.text.y = element_text(size = 100, color = "black", lineheight = 0.1),
    legend.text = element_text(size = 80),
    legend.title = element_text(size = 80, color = "black"),
    legend.frame = element_rect(color = "black", linewidth = 0.2),
    legend.ticks = element_line(color = "black", linewidth = 0.2)
  )

ggsave(filename = paste0(file_path, "dot_plot_exosome_nucleus_HepG2.png"), 
       plot = dot_plot_exosome_nucleus_HepG2, height = 4, width = 6, units = c("in"), dpi = 600)

#venn diagram
#import data
Exosome_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_Exosome_Nucleus_HepG2.xlsx",
  sheet = 'Exosome',
  col_names = TRUE,
  .name_repair = "universal"
) |> pull(UNIPROT_ACCESSION)

Nucleus_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_Exosome_Nucleus_HepG2.xlsx",
  sheet = 'Nucleus',
  col_names = TRUE,
  .name_repair = "universal"
) |> pull(UNIPROT_ACCESSION)

venn.diagram(
  x = list(Exosome_HepG2, Nucleus_HepG2),
  category.names = c("Exosome \n (73)", "Nucleus \n (239)"),
  filename = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\Exo_Nuc_HepG2_venn_diagram.png", 
  fill = c(Color_1, Color_2),
  output = TRUE,
  imagetype = "png",
  height = 2,
  width = 2,
  units = c("in"),
  resolution = 600,
  lty = 'blank',
  cex = 1,
  fontfamily = "sans",
  cat.cex = 1,
  cat.default.pos = "outer",
  cat.pos = c(-20, 150),
  cat.dist = c(0.02, 0.02),
  cat.fontfamily = "sans",
  margin = 0.05
)
