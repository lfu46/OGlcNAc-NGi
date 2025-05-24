#import packages
packages_names <- c("tidyverse", "writexl", "readxl", "ggpubr", "rstatix", "showtext",
                    "ggridges", "ComplexHeatmap", "circlize", "UpSetR")
lapply(packages_names, require, character.only = TRUE)

#import data
transcription_database <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\transcription_factor\\DatabaseExtract_v_101.xlsx",
  sheet = 'Sheet1',
  col_names = TRUE,
  .name_repair = "universal"
)
uniprot_transcription_factor <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\transcription_factor\\uniprot_transcription_factor.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)

#generate data frame
transcription_database_adj <- transcription_database |> 
  left_join(uniprot_transcription_factor, by = join_by(HGNC.symbol == From)) |> 
  select(UniprotID = Entry, Gene_Name = HGNC.symbol, DNA_binding_domain = DBD) |> 
  filter(!is.na(UniprotID))

#generate data frame
#HEK293T
OG_glycoprotein_TF_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  left_join(transcription_database_adj, by = 'UniprotID') |> 
  filter(!is.na(DNA_binding_domain)) |> 
  mutate(cell = "HEK293T", group = "OG_TF", group2 = "OG_TF_HEK293T")

write_xlsx(OG_glycoprotein_TF_HEK293T, path = paste0(file_path, "OG_glycoprotein_TF_HEK293T.xlsx"))

OG_glycoprotein_not_TF_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(! UniprotID %in% OG_glycoprotein_TF_HEK293T$UniprotID) |> 
  mutate(cell = "HEK293T", group = "OG_not_TF", group2 = "OG_not_TF")

WP_protein_TF_HEK293T <- WP_protein_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_TF_HEK293T$UniprotID) |> 
  mutate(cell = "HEK293T", group = "WP_TF", group2 = "WP_TF")

#HepG2
OG_glycoprotein_TF_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  left_join(transcription_database_adj, by = 'UniprotID') |> 
  filter(!is.na(DNA_binding_domain)) |> 
  mutate(cell = "HepG2", group = "OG_TF", group2 = "OG_TF_HepG2")

OG_glycoprotein_not_TF_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(! UniprotID %in% OG_glycoprotein_TF_HepG2$UniprotID) |> 
  mutate(cell = "HepG2", group = "OG_not_TF", group2 = "OG_not_TF")

WP_protein_TF_HepG2 <- WP_protein_Top_tb_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_TF_HepG2$UniprotID) |> 
  mutate(cell = "HepG2", group = "WP_TF", group2 = "WP_TF")

#Jurkat
OG_glycoprotein_TF_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  left_join(transcription_database_adj, by = 'UniprotID') |> 
  filter(!is.na(DNA_binding_domain)) |> 
  mutate(cell = "Jurkat", group = "OG_TF", group2 = "OG_TF_Jurkat")

OG_glycoprotein_not_TF_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(! UniprotID %in% OG_glycoprotein_TF_Jurkat$UniprotID) |> 
  mutate(cell = "Jurkat", group = "OG_not_TF", group2 = "OG_not_TF")

WP_protein_TF_Jurkat <- WP_protein_Top_tb_Jurkat |> 
  filter(UniprotID %in% OG_glycoprotein_TF_Jurkat$UniprotID) |> 
  mutate(cell = "Jurkat", group = "WP_TF", group2 = "WP_TF")

#combine
OG_glycoprotein_TF_combined <- bind_rows(
  OG_glycoprotein_TF_HEK293T,
  OG_glycoprotein_TF_HepG2,
  OG_glycoprotein_TF_Jurkat
)

write_xlsx(OG_glycoprotein_TF_combined, path = paste0(file_path, "OG_glycoprotein_TF_combined.xlsx"))

OG_glycoprotein_TF_total <- OG_glycoprotein_TF_combined |> distinct(UniprotID)

OG_glycoprotein_TF_upset_HEK293T <- OG_glycoprotein_TF_HEK293T |> select(UniprotID) |> mutate(HEK293T = 1)

OG_glycoprotein_TF_upset_HepG2 <- OG_glycoprotein_TF_HepG2 |> select(UniprotID) |> mutate(HepG2 = 1)

OG_glycoprotein_TF_upset_Jurkat <- OG_glycoprotein_TF_Jurkat |> select(UniprotID) |> mutate(Jurkat = 1)

OG_glycoprotein_TF_total_upset <- OG_glycoprotein_TF_total |> 
  left_join(OG_glycoprotein_TF_upset_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_TF_upset_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_TF_upset_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  mutate(
    HepG2 = ifelse(is.na(HepG2), 0, 1),
    HEK293T = ifelse(is.na(HEK293T), 0, 1),
    Jurkat = ifelse(is.na(Jurkat), 0, 1)
  )

OG_glycoprotein_TF_total_upset_dataframe <- as.data.frame(OG_glycoprotein_TF_total_upset)

upset(OG_glycoprotein_TF_total_upset_dataframe, sets = c("HepG2", "HEK293T", "Jurkat"),
      main.bar.color = "black", 
      order.by = "freq", point.size = 5,
      text.scale = c(2, 1.5, 1.5, 1.5, 2, 2),
      mainbar.y.label = "# of quantified \nO-GlcNAcylated TF")

#wilcox test
OG_glycoprotein_TF_wilcox_test <- OG_glycoprotein_TF_combined |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_TF <- OG_glycoprotein_TF_combined |> 
  ggplot() +
  geom_violin(aes(x = cell, y = logFC, fill = cell), color = "transparent") +
  geom_boxplot(aes(x = cell, y = logFC), outliers = FALSE, width = 0.2) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  labs(x = "", y = expression(Δlog[2]*"FC")) +
  stat_pvalue_manual(data = OG_glycoprotein_TF_wilcox_test, label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.5, 1.5)) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    axis.ticks.length.x = unit(0, "in"),
    strip.text = element_text(size = 80, color = "black", margin = margin(b = 0.01, t = 0, unit = "line")),
    strip.background = element_rect(fill = NA, color = NA),
    legend.position = "none",
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_TF.png"), plot = violin_boxplot_OG_glycoprotein_TF, height = 3, width = 3, units = c("in"), dpi = 600)

#combine
OG_glycoprotein_TF_not_TF_combined <- bind_rows(
  OG_glycoprotein_TF_HEK293T,
  OG_glycoprotein_not_TF_HEK293T,
  OG_glycoprotein_TF_HepG2,
  OG_glycoprotein_not_TF_HepG2,
  OG_glycoprotein_TF_Jurkat,
  OG_glycoprotein_not_TF_Jurkat
)

#wilcox test
OG_glycoprotein_TF_not_TF_wilcox_test <- OG_glycoprotein_TF_not_TF_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ group, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_TF_not_TF <- OG_glycoprotein_TF_not_TF_combined |> 
  ggplot() +
  geom_violin(aes(x = factor(group, levels = c("OG_TF", "OG_not_TF")), y = logFC, fill = group2), color = "transparent") +
  geom_boxplot(aes(x = factor(group, levels = c("OG_TF", "OG_not_TF")), y = logFC), outliers = FALSE, width = 0.2) +
  scale_fill_manual(
    values = c(
      "OG_TF_HEK293T" = Color_2,
      "OG_TF_HepG2" = Color_3,
      "OG_TF_Jurkat" = Color_4,
      "OG_not_TF" = "gray"
    )
  ) +
  scale_x_discrete(
    label = c("TF", "not TF")
  ) +
  labs(x = "", y = expression(Δlog[2]*"FC")) +
  stat_pvalue_manual(data = OG_glycoprotein_TF_not_TF_wilcox_test, label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.7)) +
  coord_cartesian(ylim = c(-2, 2)) +
  facet_grid(~ cell, scales = "free_x") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    axis.ticks.length.x = unit(0, "in"),
    strip.text = element_text(size = 80, color = "black", margin = margin(b = 0.01, t = 0, unit = "line")),
    strip.background = element_rect(fill = NA, color = NA),
    legend.position = "none",
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_TF_not_TF.png"), plot = violin_boxplot_OG_glycoprotein_TF_not_TF, height = 3, width = 4, units = c("in"), dpi = 600)

#combine
OG_WP_TF_combined <- bind_rows(
  OG_glycoprotein_TF_HEK293T,
  WP_protein_TF_HEK293T,
  OG_glycoprotein_TF_HepG2,
  WP_protein_TF_HepG2,
  OG_glycoprotein_TF_Jurkat,
  WP_protein_TF_Jurkat
)

write_xlsx(OG_WP_TF_combined, path = paste0(file_path, "OG_WP_TF_combined.xlsx"))

#wilcox test
OG_WP_TF_wilcox_test <- OG_WP_TF_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ group, p.adjust.method = "BH") |> 
  add_significance("p")

#point boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

point_boxplot_OG_WP_TF <- OG_WP_TF_combined |> 
  ggplot() +
  geom_line(aes(x = factor(group, levels = c("WP_TF", "OG_TF")), y = logFC, group = UniprotID), color = "gray") +
  geom_boxplot(aes(x = factor(group, levels = c("WP_TF", "OG_TF")), y = logFC), outliers = FALSE) +
  geom_point(aes(x = factor(group, levels = c("WP_TF", "OG_TF")), y = logFC, fill = group2), 
             shape = 21, size = 2, color = "black") +
  scale_fill_manual(values = c(
    "OG_TF_HEK293T" = Color_2,
    "OG_TF_HepG2" = Color_3,
    "OG_TF_Jurkat" = Color_4,
    "WP_TF" = "gray"
  )) +
  scale_x_discrete(
    label = c("WP", "OG")
  ) +
  stat_pvalue_manual(data = OG_WP_TF_wilcox_test, label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.7)) +
  facet_grid(~ cell, scale = "free_x") +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)"), fill = "") +
  theme_bw() +
  theme(
    title = element_text(size = 100, color = "black"),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    axis.ticks.length.x = unit(0, "in"),
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = NA, color = NA),
    legend.position = "none",
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "point_boxplot_OG_WP_TF.png"), plot = point_boxplot_OG_WP_TF, height = 3, width = 4, units = c("in"), dpi = 600)

#generate data frame
OG_glycoprotein_DBD_HEK293T <- OG_glycoprotein_TF_HEK293T |> 
  select(DNA_binding_domain, logFC, cell) |> 
  separate_rows(DNA_binding_domain, sep = "; ")

OG_glycoprotein_DBD_HepG2 <- OG_glycoprotein_TF_HepG2 |> 
  select(DNA_binding_domain, logFC, cell) |> 
  separate_rows(DNA_binding_domain, sep = "; ")

OG_glycoprotein_DBD_Jurkat <- OG_glycoprotein_TF_Jurkat |> 
  select(DNA_binding_domain, logFC, cell) |> 
  separate_rows(DNA_binding_domain, sep = "; ")

#combine
OG_glycoprotein_DBD_combined <- bind_rows(
  OG_glycoprotein_DBD_HEK293T,
  OG_glycoprotein_DBD_HepG2,
  OG_glycoprotein_DBD_Jurkat
)

write_xlsx(OG_glycoprotein_DBD_combined, path = paste0(file_path, "OG_glycoprotein_DBD_combined.xlsx"))

#wilcox test
OG_glycoprotein_DBD_HEK293T_wilcox_test <- OG_glycoprotein_DBD_HEK293T |> 
  wilcox_test(logFC ~ DNA_binding_domain, p.adjust.method = "BH") |> 
  add_significance("p") |> 
  filter(!p.signif == "ns") |> 
  filter(n1 > 4, n2 > 4) |> 
  mutate(cell = "HEK293T")

OG_glycoprotein_DBD_HepG2_wilcox_test <- OG_glycoprotein_DBD_HepG2 |> 
  wilcox_test(logFC ~ DNA_binding_domain, p.adjust.method = "BH") |> 
  add_significance("p") |> 
  filter(!p.signif == "ns") |> 
  mutate(cell = "HepG2")

OG_glycoprotein_DBD_Jurkat_wilcox_test <- OG_glycoprotein_DBD_Jurkat |> 
  wilcox_test(logFC ~ DNA_binding_domain, p.adjust.method = "BH") |> 
  add_significance("p") |> 
  filter(!p.signif == "ns") |> 
  mutate(cell = "Jurkat")

OG_glycoprotein_DBD_wilcox_test_combined <- bind_rows(
  OG_glycoprotein_DBD_HEK293T_wilcox_test,
  OG_glycoprotein_DBD_HepG2_wilcox_test,
  OG_glycoprotein_DBD_Jurkat_wilcox_test
)

write_xlsx(OG_glycoprotein_DBD_wilcox_test_combined, path = paste0(file_path, "OG_glycoprotein_DBD_wilcox_test_combined.xlsx"))

#combine
OG_glycoprotein_DBD_HEK293T_list <- c(
  "bHLH", "C2H2 ZF", "Ets", "Forkhead", "GATA", "Homeodomain", "Rel"
)
OG_glycoprotein_DBD_HepG2_list <- c("C2H2 ZF", "Homeodomain")
OG_glycoprotein_DBD_Jurkat_list <- c("C2H2 ZF")

OG_glycoprotein_DBD_combined_adj <- bind_rows(
  OG_glycoprotein_DBD_HEK293T |> filter(DNA_binding_domain %in% OG_glycoprotein_DBD_HEK293T_list),
  OG_glycoprotein_DBD_HepG2 |> filter(DNA_binding_domain %in% OG_glycoprotein_DBD_HepG2_list),
  OG_glycoprotein_DBD_Jurkat |> filter(DNA_binding_domain %in% OG_glycoprotein_DBD_Jurkat_list)
)

#label
OG_glycoprotein_DBD_wilcox_test_label <- bind_rows(
  OG_glycoprotein_DBD_wilcox_test_combined |> 
    filter(group1 == "Ets", group2 == "Homeodomain", cell == "HEK293T"),
  OG_glycoprotein_DBD_wilcox_test_combined |> 
    filter(group1 == "GATA", group2 == "Homeodomain", cell == "HEK293T"),
  OG_glycoprotein_DBD_wilcox_test_combined |> 
    filter(group1 == "Homeodomain", group2 == "Rel", cell == "HEK293T")
) |> 
  mutate(group1 = ifelse(group2 %in% c("Rel"), group2, group1))

#ggridges
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

ggridges_OG_glycoprotein_DBD <- OG_glycoprotein_DBD_combined_adj |> 
  ggplot() +
  geom_density_ridges(
    aes(x = logFC, y = DNA_binding_domain, fill = cell), color = "black",
    quantile_lines = TRUE, quantiles = 2, scale = 0.55,
    jittered_points = TRUE, position = position_points_jitter(width = 0.05, height = 0),
    point_shape = "|", point_size = 8, point_alpha = 1
  ) +
  labs(x = expression(Δlog[2]*"FC"), y = "") +
  geom_text(data = OG_glycoprotein_DBD_wilcox_test_label, aes(y = group1, x = 1.6, label = p.signif), 
            size = 50) +
  coord_cartesian(xlim = c(-2, 2)) +
  facet_grid(cell ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black"),
    strip.text = element_text(size = 100, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "ggridges_OG_glycoprotein_DBD.png"), 
       plot = ggridges_OG_glycoprotein_DBD, height = 6, width = 4, units = "in", dpi = 600)

#overlap heatmap
#import data
uniprot_OG_glycoprotein_TF_overlap <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\transcription_factor\\uniprot_OG_glycoprotein_TF_overlap_2.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
) |> separate(col = Gene.Names, sep = " ", into = 'Gene_Name')

#generate data frame
#HEK293T
OG_glycoprotein_TF_overlap_HEK293T <- uniprot_OG_glycoprotein_TF_overlap |> 
  left_join(OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T, by = join_by(From == UniprotID)) |> 
  select(Gene_Name, logFC, Annotation) |> mutate(cell = "HEK293T")

#HepG2
OG_glycoprotein_TF_overlap_HepG2 <- uniprot_OG_glycoprotein_TF_overlap |> 
  left_join(OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2, by = join_by(From == UniprotID)) |> 
  select(Gene_Name, logFC, Annotation) |> mutate(cell = "HepG2")

#Jurkat
OG_glycoprotein_TF_overlap_Jurkat <- uniprot_OG_glycoprotein_TF_overlap |> 
  left_join(OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat, by = join_by(From == UniprotID)) |> 
  select(Gene_Name, logFC, Annotation) |> mutate(cell = "Jurkat")

#combine
OG_glycoprotein_TF_overlap_combined <- bind_rows(
  OG_glycoprotein_TF_overlap_HEK293T,
  OG_glycoprotein_TF_overlap_HepG2,
  OG_glycoprotein_TF_overlap_Jurkat
) |> mutate(Annotation = ifelse(is.na(Annotation), "Other", Annotation))

OG_glycoprotein_TF_overlap_combined_adj <- OG_glycoprotein_TF_overlap_combined |> 
  arrange(Annotation)

write_xlsx(OG_glycoprotein_TF_overlap_combined_adj, path = paste0(file_path, "OG_glycoprotein_TF_overlap_combined_adj.xlsx"))

OG_glycoprotein_TF_overlap_combined_adj_2 <- OG_glycoprotein_TF_overlap_combined_adj |> 
  select(cell, Gene_Name, logFC) |> 
  pivot_wider(names_from = 'Gene_Name', values_from = 'logFC')

#Heatmap
mat <- data.matrix(OG_glycoprotein_TF_overlap_combined_adj_2 |> select(!cell))
rownames(mat) <- OG_glycoprotein_TF_overlap_combined_adj_2$cell
split_category <- OG_glycoprotein_TF_overlap_combined_adj |> filter(cell == "HEK293T") |> pull(Annotation)
mat_adj <- t(mat)

col_mat <- colorRamp2(breaks = c(-1.5, 0, 1.5), colors = c("blue", "white", "red"))

h1 <- Heatmap(matrix = mat_adj, col = col_mat, name = 'Δlog2FC',
              row_split = split_category, row_names_side = "left", row_title = NULL,
              show_row_dend = FALSE,
              border = TRUE)
