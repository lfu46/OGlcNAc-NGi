#import packages
packages_names <- c("tidyverse", "writexl", "readxl", "ComplexHeatmap", "circlize", 
                    "ggpubr", "rstatix", "showtext", "ggbeeswarm")
lapply(packages_names, require, character.only = TRUE)

#import database
HumanMitoCarta_3 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\mito\\Human.MitoCarta3.xlsx",
  sheet = 'A Human MitoCarta3.0',
  col_names = TRUE,
  .name_repair = "universal"
) 

Mito_localization <- HumanMitoCarta_3 |> 
  select(Symbol, UniprotID = UniProt, MitoCarta3.0_SubMitoLocalization)

Mito_pathway <- HumanMitoCarta_3 |> 
  select(Symbol, UniprotID = UniProt, MitoCarta3.0_MitoPathways) |> 
  separate_rows(MitoCarta3.0_MitoPathways, sep = ' \\| ') |> 
  distinct()

#MitoCop
MitoCop_all_protein_group <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\mito\\MitoCoP.xlsx",
  sheet = '(A) All protein groups',
  skip = 1,
  col_names = TRUE,
  .name_repair = "universal"
)

MitoCop_all_protein_group_adj <- MitoCop_all_protein_group |> 
  select(UniprotID = Simplified.protein.IDs, MitoCoP, 
         substractive = Classification...32, spatial = Classification...41,
         importomics = Classification...66, membrane_association = Category, 
         mito_copy = Mean.mito.copies.per.cell...2.3.Reps.,
         cell_copy = Mean.copy.numbers.per.cell...2.3.Reps., half_lives = Protein.half.lives..h.) |> 
  filter(MitoCoP == 1) |> 
  separate_rows(UniprotID, sep = ";") |> 
  distinct()

MitoCop_all_protein_group_adj_2 <- MitoCop_all_protein_group |> 
  select(UniprotID = Simplified.protein.IDs, MitoCoP,
         'Mitochondrial.outer.membrane..GO.0005741...all.evidence.sources.':'Fe.S.protein.biogenesis') |> 
  filter(MitoCoP == 1) |> 
  separate_rows(UniprotID, sep = ";") |> 
  pivot_longer('Mitochondrial.outer.membrane..GO.0005741...all.evidence.sources.':'Fe.S.protein.biogenesis',
               names_to = 'submitochondrial', values_to = 'value') |> 
  filter(!is.na(value))

MitoCop_mito <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\mito\\MitoCoP.xlsx",
  sheet = '(B) MitoCoP (1,134 genes)',
  skip = 1,
  col_names = TRUE,
  .name_repair = "universal"
)

MitoCop_mito_function <- MitoCop_mito |> 
  select(UniprotID = Simplified.protein.IDs, MitoCoP, 
         Morphology..dynamics...organization:Unknown) |> 
  filter(MitoCoP == 1) |> 
  separate_rows(UniprotID, sep = ";") |> 
  pivot_longer(Morphology..dynamics...organization:Unknown, 
               names_to = 'Function', values_to = 'value') |> 
  filter(!is.na(value))

MitoCop_mito_disease <- MitoCop_mito |> 
  select(UniprotID = Simplified.protein.IDs, MitoCoP.disease.gene, 
         Central.nervous.system:Tumors) |> 
  filter(MitoCoP.disease.gene == 1) |> 
  separate_rows(UniprotID, sep = ";") |> 
  pivot_longer(Central.nervous.system:Tumors, 
               names_to = 'Function', values_to = 'value') |> 
  filter(!is.na(value)) |> 
  distinct()

#generate data frame
#HEK293T
OG_glycoprotein_mito_1_HEK293T <- MitoCop_all_protein_group_adj |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HEK293T")

OG_glycoprotein_mito_2_HEK293T <- MitoCop_all_protein_group_adj_2 |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HEK293T")

OG_glycoprotein_mito_3_HEK293T <- MitoCop_mito_disease |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HEK293T")

#HepG2
OG_glycoprotein_mito_1_HepG2 <- MitoCop_all_protein_group_adj |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HepG2")

OG_glycoprotein_mito_2_HepG2 <- MitoCop_all_protein_group_adj_2 |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HepG2")

OG_glycoprotein_mito_3_HepG2 <- MitoCop_mito_disease |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HepG2")

#Jurkat
OG_glycoprotein_mito_1_Jurkat <- MitoCop_all_protein_group_adj |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "Jurkat")

OG_glycoprotein_mito_2_Jurkat <- MitoCop_all_protein_group_adj_2 |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "Jurkat")

OG_glycoprotein_mito_3_Jurkat <- MitoCop_mito_disease |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "Jurkat")

#combined
OG_glycoprotein_mitocop_mito_1_combined <- bind_rows(
  OG_glycoprotein_mito_1_HEK293T,
  OG_glycoprotein_mito_1_HepG2,
  OG_glycoprotein_mito_1_Jurkat
)

write_xlsx(OG_glycoprotein_mitocop_mito_1_combined, 
           path = paste0(file_path, "OG_glycoprotein_mitocop_mito_1_combined.xlsx"))

OG_glycoprotein_mitocop_mito_2_combined <- bind_rows(
  OG_glycoprotein_mito_2_HEK293T,
  OG_glycoprotein_mito_2_HepG2,
  OG_glycoprotein_mito_2_Jurkat
)

write_xlsx(OG_glycoprotein_mitocop_mito_2_combined,
           path = paste0(file_path, "OG_glycoprotein_mitocop_mito_2_combined.xlsx"))

OG_glycoprotein_mitocop_mito_3_combined <- bind_rows(
  OG_glycoprotein_mito_3_HEK293T,
  OG_glycoprotein_mito_3_HepG2,
  OG_glycoprotein_mito_3_Jurkat
)

write_xlsx(OG_glycoprotein_mitocop_mito_3_combined,
           path = paste0(file_path, "OG_glycoprotein_mitocop_mito_3_combined.xlsx"))

#spatial proteomic
#wilcox test
OG_glycoprotein_mitocop_spatial_wilcox_test <- OG_glycoprotein_mitocop_mito_1_combined |> 
  filter(
    spatial == "Cluster 1" | spatial == "Cluster 2"
    ) |> 
  group_by(spatial) |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#point boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

point_boxplot_OG_glycoprotein_mito_spatial <- OG_glycoprotein_mitocop_mito_1_combined |> 
  filter(spatial != "NaN") |> 
  ggplot() +
  geom_jitter(aes(x = cell, y = logFC, color = cell), position = position_jitter(width = 0.2)) +
  geom_boxplot(aes(x = cell, y = logFC), color = "black", outliers = FALSE, fill = "transparent") +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_color_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_grid(~ spatial, scales = "free_x") +
  stat_pvalue_manual(data = OG_glycoprotein_mitocop_spatial_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.8, 1.7, 1.8)) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent")
  )

ggsave(filename = paste0(file_path, "point_boxplot_OG_glycoprotein_mito_spatial.png"), 
       plot = point_boxplot_OG_glycoprotein_mito_spatial, height = 3, width = 4, units = "in", dpi = 600)

#membrane association
#wilcox test
OG_glycoprotein_mitocop_membrane_wilcox_test <- OG_glycoprotein_mitocop_mito_1_combined |> 
  filter(
    membrane_association %in% c("Integral", "Peripheral", "Soluble", "Ambiguous")
  ) |> 
  mutate(group = ifelse(membrane_association %in% c("Integral", "Peripheral"), "Integral & Peripheral", "Soluble & Ambiguous")) |> 
  group_by(group) |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#point boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()
  
point_boxplot_OG_glycoprotein_mito_membrane <- OG_glycoprotein_mitocop_mito_1_combined |> 
    filter(
      membrane_association %in% c("Integral", "Peripheral", "Soluble", "Ambiguous")
    ) |> 
    mutate(group = ifelse(membrane_association %in% c("Integral", "Peripheral"), "Integral & Peripheral", "Soluble & Ambiguous")) |> 
  filter(group == "Integral & Peripheral") |>   
  ggplot(aes(x = cell, y = logFC)) +
  geom_jitter(aes(color = cell), position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", outliers = FALSE, fill = "transparent") +
    labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
    scale_color_manual(
      values = c(
        "HEK293T" = Color_2,
        "HepG2" = Color_3,
        "Jurkat" = Color_4
      )
    ) +
    stat_pvalue_manual(data = OG_glycoprotein_mitocop_membrane_wilcox_test, label = "p.signif", 
                       tip.length = 0, size = 50,
                       hide.ns = "p", y.position = c(1.8, 1.6)) +
    coord_cartesian(ylim = c(-2, 2)) +
  ggtitle("Integral & Peripheral") +
    theme_bw() +
    theme(
      title = element_text(color = "black", size = 70, lineheight = 0.15),
      panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
      panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
      axis.title = element_text(size = 100),
      axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
      axis.text.y = element_text(size = 80, color = "black"),
      legend.position = "none",
      strip.text = element_text(size = 80, color = "black"),
      strip.background = element_rect(fill = "transparent", color = "transparent")
    )

ggsave(filename = paste0(file_path, "point_boxplot_OG_glycoprotein_mito_membrane.png"), 
         plot = point_boxplot_OG_glycoprotein_mito_membrane, height = 3, width = 2.5, units = "in", dpi = 600)

#half lives
#wilcox test
OG_glycoprotein_mitocop_half_lives_wilcox_test <- OG_glycoprotein_mitocop_mito_1_combined |> 
  filter(!is.na(half_lives)) |>
  filter(!half_lives == "NA") |> 
  mutate(half_lives = as.numeric(half_lives)) |> 
  mutate(group = ifelse(half_lives < 50, "short", "long")) |> 
  group_by(group) |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_mito_half_lives <- OG_glycoprotein_mitocop_mito_1_combined |> 
  filter(!is.na(half_lives)) |>
  filter(!half_lives == "NA") |> 
  mutate(half_lives = as.numeric(half_lives)) |> 
  mutate(group = ifelse(half_lives < 50, "short", "long")) |> 
  ggplot() +
  geom_violin(aes(x = cell, y = logFC, fill = cell), color = "transparent") +
  geom_boxplot(aes(x = cell, y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_grid(~ factor(group, levels = c("short", "long")), scales = "free_x") +
  stat_pvalue_manual(data = OG_glycoprotein_mitocop_half_lives_wilcox_test, label = "p.adj.signif", 
                     tip.length = 0, size = 50,
                     hide.ns = "p.adj", y.position = c(1.8)) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_mito_half_lives.png"), 
       plot = violin_boxplot_OG_glycoprotein_mito_half_lives, height = 3, width = 4, units = "in", dpi = 600)

#mito copy
#wilcox test
OG_glycoprotein_mitocop_mito_copy_wilcox_test <- OG_glycoprotein_mitocop_mito_1_combined |> 
  filter(!mito_copy == "NaN") |>
  mutate(mito_copy = as.numeric(mito_copy)) |> 
  mutate(group = ifelse(mito_copy < 5E5, "Mito.Abs.low", "Mito.Abs.high")) |> 
  group_by(group) |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_mito_mito_copy <- OG_glycoprotein_mitocop_mito_1_combined |> 
  filter(!mito_copy == "NaN") |>
  mutate(mito_copy = as.numeric(mito_copy)) |> 
  mutate(group = ifelse(mito_copy < 5E5, "Mito.Abs.low", "Mito.Abs.high")) |> 
  ggplot() +
  geom_violin(aes(x = cell, y = logFC, fill = cell), color = "transparent") +
  geom_boxplot(aes(x = cell, y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_grid(~ factor(group, levels = c("Mito.Abs.low", "Mito.Abs.high")), scales = "free_x") +
  stat_pvalue_manual(data = OG_glycoprotein_mitocop_mito_copy_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.5, 1.8)) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_mito_mito_copy.png"), 
       plot = violin_boxplot_OG_glycoprotein_mito_mito_copy, height = 3, width = 4, units = "in", dpi = 600)

#cell copy
#wilcox test
OG_glycoprotein_mitocop_cell_copy_wilcox_test <- OG_glycoprotein_mitocop_mito_1_combined |> 
  filter(!cell_copy == "NaN") |>
  mutate(cell_copy = as.numeric(cell_copy)) |> 
  mutate(group = ifelse(cell_copy < 6E5, "Cell.Abs.low", "Cell.Abs.high")) |> 
  group_by(group) |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_mito_cell_copy <- OG_glycoprotein_mitocop_mito_1_combined |> 
  filter(!cell_copy == "NaN") |>
  mutate(cell_copy = as.numeric(cell_copy)) |> 
  mutate(group = ifelse(cell_copy < 6E5, "Cell.Abs.low", "Cell.Abs.high")) |> 
  ggplot() +
  geom_violin(aes(x = cell, y = logFC, fill = cell), color = "transparent") +
  geom_boxplot(aes(x = cell, y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_grid(~ factor(group, levels = c("Cell.Abs.low", "Cell.Abs.high")), scales = "free_x") +
  stat_pvalue_manual(data = OG_glycoprotein_mitocop_cell_copy_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.5, 1.8)) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_mito_cell_copy.png"), 
       plot = violin_boxplot_OG_glycoprotein_mito_cell_copy, height = 3, width = 4, units = "in", dpi = 600)

#sub mito 
#wilcox test
OG_glycoprotein_mitocop_submito_wilcox_test <- OG_glycoprotein_mitocop_mito_2_combined |> 
  mutate(group = ifelse(submitochondrial == "Mitochondrial.inner.membrane..GO.0005743...all.evidence.sources.", 'Mito.inner', NA),
         group = ifelse(submitochondrial == "Mitochondrial.outer.membrane..GO.0005741...all.evidence.sources.", 'Mito.outer', group),
         group = ifelse(submitochondrial == "Mitochondrial.matrix..GO.0005759...all.evidence.sources.", 'Mito.matrix', group), 
         group = ifelse(is.na(group), 'Other', group)) |> 
  group_by(group) |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

point_boxplot_OG_glycoprotein_mito_sub_mito <- OG_glycoprotein_mitocop_mito_2_combined |> 
  mutate(group = ifelse(submitochondrial == "Mitochondrial.inner.membrane..GO.0005743...all.evidence.sources.", 'Mito.inner', NA),
         group = ifelse(submitochondrial == "Mitochondrial.outer.membrane..GO.0005741...all.evidence.sources.", 'Mito.outer', group),
         group = ifelse(submitochondrial == "Mitochondrial.matrix..GO.0005759...all.evidence.sources.", 'Mito.matrix', group), 
         group = ifelse(is.na(group), 'Other', group)) |> 
  ggplot(aes(x = cell, y = logFC)) +
  geom_jitter(aes(color = cell), position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", outliers = FALSE, fill = "transparent") +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_color_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_grid(~ factor(group, levels = c("Mito.outer", "Mito.inner", "Mito.matrix", "Other")), scales = "free_x") +
  stat_pvalue_manual(data = OG_glycoprotein_mitocop_submito_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 6,
                     hide.ns = "p", y.position = c(2.2, 2.1, 2.1, 1.8, 2.3, 1.7, 2.0)) +
  coord_cartesian(ylim = c(-1.2, 2.5)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(size = 10, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 12, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent")
  )

ggsave(filename = paste0(file_path, "point_boxplot_OG_glycoprotein_mito_sub_mito.eps"), 
       device = "eps",
       plot = point_boxplot_OG_glycoprotein_mito_sub_mito, 
       height = 2, width = 4, 
       units = "in", 
       dpi = 1200
       )

#disease
#wilcox test
OG_glycoprotein_mitocop_disease_wilcox_test <- OG_glycoprotein_mitocop_mito_3_combined |> 
  filter(Function == "Metabolism" | Function == "Liver") |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#point boxplot
font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

point_boxplot_OG_glycoprotein_mito_disease <- OG_glycoprotein_mitocop_mito_3_combined |> 
  ggplot(aes(x = cell, y = logFC)) +
  geom_jitter(aes(color = cell), position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", fill = "transparent", outliers = FALSE) +
  labs(x = "", y = "") +
  scale_color_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  stat_pvalue_manual(data = OG_glycoprotein_mitocop_disease_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 6,
                     hide.ns = "p", y.position = c(1.4, 1.7, 2.0)) +
  coord_cartesian(ylim = c(-1.8, 2.2)) +
  theme_bw() +
  theme(
    title = element_text(color = "black", size = 10, lineheight = 0.15),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "point_boxplot_OG_glycoprotein_mito_disease.eps"), 
       device = "eps",
       plot = point_boxplot_OG_glycoprotein_mito_disease, 
       height = 2, width = 2, 
       units = "in", 
       dpi = 1200
       )

#generate data frame
#HEK293T
OG_glycoprotein_mito_pathway_HEK293T <- Mito_pathway |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HEK293T")

#HepG2
OG_glycoprotein_mito_pathway_HepG2 <- Mito_pathway |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HepG2")

#Jurkat
OG_glycoprotein_mito_pathway_Jurkat <- Mito_pathway |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "Jurkat")

#combined
OG_glycoprotein_mito_pathway_combined <- bind_rows(
  OG_glycoprotein_mito_pathway_HEK293T,
  OG_glycoprotein_mito_pathway_HepG2,
  OG_glycoprotein_mito_pathway_Jurkat
)

write_xlsx(OG_glycoprotein_mito_pathway_combined, path = paste0(file_path, "OG_glycoprotein_mito_pathway_combined.xlsx"))

#pathway
#wilcox test
OG_glycoprotein_mito_pathway_wilcox_test <- OG_glycoprotein_mito_pathway_combined |> 
  filter(str_detect(MitoCarta3.0_MitoPathways, "Mitochondrial central dogma")) |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

point_boxplot_OG_glycoprotein_mito_pathway <- OG_glycoprotein_mito_pathway_combined |> 
  ggplot(aes(x = cell, y = logFC)) +
  geom_jitter(aes(color = cell), position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", fill = "transparent", outliers = FALSE) +
  labs(x = "", y = "") +
  scale_color_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  stat_pvalue_manual(data = OG_glycoprotein_mito_pathway_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 6,
                     hide.ns = "p", y.position = c(1.7, 2.2, 1.8)) +
  coord_cartesian(ylim = c(-1.8, 2.4)) +
  theme_bw() +
  theme(
    title = element_text(color = "black", size = 10, lineheight = 0.15),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "point_boxplot_OG_glycoprotein_mito_pathway.eps"), 
       device = "eps",
       plot = point_boxplot_OG_glycoprotein_mito_pathway, 
       height = 2, width = 2, 
       units = "in", 
       dpi = 1200
       )
