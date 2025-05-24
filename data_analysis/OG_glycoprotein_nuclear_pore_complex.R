#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "showtext", 
                    "ggpubr", "rstatix", "ComplexHeatmap", "circlize", "eulerr")
lapply(packages_names, require, character.only = TRUE)

#import database
nuclear_pore_complex_database <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\nuclear_pore_complex\\uniprot_nuclear_pore_complex.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)

#generate data frame
#HEK293T
OG_glycoprotein_npc_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(UniprotID %in% nuclear_pore_complex_database$Entry) |> 
  mutate(cell = "HEK293T", group = "OG_HEK293T")

WP_protein_npc_HEK293T <- WP_protein_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_npc_HEK293T$UniprotID) |> 
  mutate(cell = "HEK293T", group = "WP")

#HepG2
OG_glycoprotein_npc_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(UniprotID %in% nuclear_pore_complex_database$Entry) |> 
  mutate(cell = "HepG2", group = "OG_HepG2")

WP_protein_npc_HepG2 <- WP_protein_Top_tb_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_npc_HepG2$UniprotID) |> 
  mutate(cell = "HepG2", group = "WP")

#Jurkat
OG_glycoprotein_npc_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(UniprotID %in% nuclear_pore_complex_database$Entry) |> 
  mutate(cell = "Jurkat", group = "OG_Jurkat")

WP_protein_npc_Jurkat <- WP_protein_Top_tb_Jurkat |> 
  filter(UniprotID %in% OG_glycoprotein_npc_Jurkat$UniprotID) |> 
  mutate(cell = "Jurkat", group = "WP")

#combine OG
OG_glycoprotein_npc_combined <- bind_rows(
  OG_glycoprotein_npc_HEK293T,
  OG_glycoprotein_npc_HepG2,
  OG_glycoprotein_npc_Jurkat
)

write_xlsx(OG_glycoprotein_npc_combined, path = paste0(file_path, "OG_glycoprotein_npc_combined.xlsx"))

#wilcox test
OG_glycoprotein_npc_wilcox_test <- OG_glycoprotein_npc_combined |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_npc <- OG_glycoprotein_npc_combined |> 
  ggplot() +
  geom_violin(aes(x = cell, y = logFC, fill = cell), 
              color = "transparent") +
  geom_boxplot(aes(x = cell, y = logFC), color = "black",
               outliers = FALSE, width = 0.2) +
  labs(x = "", y = expression(Δlog[2]*"FC")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  stat_pvalue_manual(data = OG_glycoprotein_npc_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.8)) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_npc.png"), plot = violin_boxplot_OG_glycoprotein_npc, height = 4, width = 4, units = "in", dpi = 600)

#combine OG WP
OG_WP_npc_combined <- bind_rows(
  OG_glycoprotein_npc_HEK293T,
  WP_protein_npc_HEK293T,
  OG_glycoprotein_npc_HepG2,
  WP_protein_npc_HepG2,
  OG_glycoprotein_npc_Jurkat,
  WP_protein_npc_Jurkat
)

#wilcox test
OG_WP_npc_wilcox_test <- OG_WP_npc_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ group, p.adjust.method = "BH") |> 
  add_significance("p")

#point
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

point_plot_OG_WP_npc <- OG_WP_npc_combined |> 
  ggplot(aes(x = factor(group, levels = c("WP", "OG_HEK293T", "OG_HepG2", "OG_Jurkat")), y = logFC)) +
  geom_line(aes(group = UniprotID), color = "gray") +
  geom_boxplot() +
  geom_point(aes(fill = group), shape = 21,
                 color = "black") +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(
    values = c(
      "OG_HEK293T" = Color_2,
      "OG_HepG2" = Color_3,
      "OG_Jurkat" = Color_4,
      "WP" = "gray"
    )
  ) +
  scale_x_discrete(labels = c("WP", "OG")) +
  stat_pvalue_manual(data = OG_WP_npc_wilcox_test, label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.5)) +
  coord_cartesian(ylim = c(-2, 2)) +
  facet_grid(~ cell, scale = "free_x") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "point_plot_OG_WP_npc.png"), plot = point_plot_OG_WP_npc, height = 4, width = 4, units = "in", dpi = 600)

#heatmap
write_xlsx(OG_glycoprotein_npc_combined, path = paste0(file_path, "OG_glycoprotein_npc_combined.xlsx"))

OG_glycoprotein_npc_overlap_protein_list <- c(
  "P12270", "P52948", "P49790", "P55735", "P09651", "P35658", "O60318", "P52597", "P49792"
)

OG_glycoprotein_npc_overlap_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_npc_overlap_protein_list) |> 
  mutate(cell = "HEK293T")
OG_glycoprotein_npc_overlap_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_npc_overlap_protein_list) |> 
  mutate(cell = "HepG2")
OG_glycoprotein_npc_overlap_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(UniprotID %in% OG_glycoprotein_npc_overlap_protein_list) |> 
  mutate(cell = "Jurkat")

#import data from uniprot
npc_uniprot <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\nuclear_pore_complex\\npc_uniprot.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)

#combine
OG_glycoprotein_npc_overlap <- bind_rows(
  OG_glycoprotein_npc_overlap_HEK293T,
  OG_glycoprotein_npc_overlap_HepG2,
  OG_glycoprotein_npc_overlap_Jurkat
) |> left_join(npc_uniprot, by = join_by(UniprotID == From)) |> 
  select(logFC, cell, Gene.Names) |> 
  pivot_wider(names_from = Gene.Names, values_from = logFC)

#Heatmap
mat <- t(data.matrix(OG_glycoprotein_npc_overlap))
colnames(mat) <- OG_glycoprotein_npc_overlap$cell
mat_adj <- mat[-1,]

col_mat <- colorRamp2(breaks = c(-1, 0, 1), colors = c("blue", "white", "red"))

h1 <- Heatmap(mat_adj, col = col_mat, name = 'Δlog2FC')


