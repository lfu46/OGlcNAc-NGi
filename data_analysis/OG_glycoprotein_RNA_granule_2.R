#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "showtext", 
                    "ggpubr", "rstatix", "ComplexHeatmap", "circlize", "eulerr")
lapply(packages_names, require, character.only = TRUE)

#import database
SG_PB_combined <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\RNA_Granule_Database_Uniprot.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = ProteinID) |> distinct()

SG_PB_overlap <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\PB_only_homo_sapiens.csv",
  col_names = TRUE,
  name_repair = "universal",
  skip = 2
) |> select(Gene.Name) |> distinct()

SG_only <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\SG_only_homo_sapiens.csv",
  col_names = TRUE,
  name_repair = "universal",
  skip = 2
) |> select(Gene.Name) |> distinct() |> 
  filter(! Gene.Name %in% SG_PB_overlap$Gene.Name)

PB_only <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\PB_only_homo_sapiens.csv",
  col_names = TRUE,
  name_repair = "universal",
  skip = 2
) |> select(Gene.Name) |> distinct() |> 
  filter(! Gene.Name %in% SG_PB_overlap$Gene.Name)

#import uniprot
SG_uniprot <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\SG_uniprot.xlsx",
  col_names = TRUE
)
PB_uniprot <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\PB_uniprot.xlsx",
  col_names = TRUE
)
SG_or_PB_uniprot <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\SG_or_PB_uniprot.xlsx",
  col_names = TRUE
)

#generate data frame
SG_only_uniprot <- SG_only |> 
  left_join(SG_uniprot, by = join_by(Gene.Name == From)) |> 
  select(Gene.Name, UniprotID = Entry) |> 
  mutate(category = "SG")

PB_only_uniprot <- PB_only |> 
  left_join(PB_uniprot, by = join_by(Gene.Name == From)) |> 
  select(Gene.Name, UniprotID = Entry) |> 
  mutate(category = "PB")

SG_PB_overlap_uniprot <- SG_PB_overlap |> 
  left_join(SG_or_PB_uniprot, by = join_by(Gene.Name == From)) |> 
  select(Gene.Name, UniprotID = Entry) |> 
  mutate(category = "Both")

RNA_granule_combined <- bind_rows(
  SG_only_uniprot,
  PB_only_uniprot,
  SG_PB_overlap_uniprot
)

#generate data frame
#HEK293T
OG_glycoprotein_RNA_granule_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  mutate(cell = "HEK293T") |> 
  left_join(RNA_granule_combined, by = "UniprotID") |> 
  filter(!is.na(category))

OG_glycoprotein_not_RNA_granule_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(! UniprotID %in% OG_glycoprotein_RNA_granule_HEK293T$UniprotID) |> 
  mutate(cell = "HEK293T", group = "OG_HEK293T", category = "not_granule")

WP_protein_RNA_granule_HEK293T <- WP_protein_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_RNA_granule_HEK293T$UniprotID) |> 
  left_join(RNA_granule_combined, by = "UniprotID") |> 
  mutate(cell = "HEK293T", group = "WP_RBP", group2 = "WP")

#HepG2
OG_glycoprotein_RNA_granule_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  mutate(cell = "HepG2") |> 
  left_join(RNA_granule_combined, by = "UniprotID") |> 
  filter(!is.na(category))

OG_glycoprotein_not_RNA_granule_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(! UniprotID %in% OG_glycoprotein_RNA_granule_HepG2$UniprotID) |> 
  mutate(cell = "HepG2", group = "OG_HepG2", category = "not_granule")

WP_protein_RNA_granule_HepG2 <- WP_protein_Top_tb_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_RNA_granule_HepG2$UniprotID) |> 
  left_join(RNA_granule_combined, by = "UniprotID") |> 
  mutate(cell = "HepG2", group = "WP_RBP", group2 = "WP")

#Jurkat
OG_glycoprotein_RNA_granule_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  mutate(cell = "Jurkat") |> 
  left_join(RNA_granule_combined, by = "UniprotID") |> 
  filter(!is.na(category))

OG_glycoprotein_not_RNA_granule_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(! UniprotID %in% OG_glycoprotein_RNA_granule_Jurkat$UniprotID) |> 
  mutate(cell = "Jurkat", group = "OG_Jurkat", category = "not_granule")

WP_protein_RNA_granule_Jurkat <- WP_protein_Top_tb_Jurkat |> 
  filter(UniprotID %in% OG_glycoprotein_RNA_granule_Jurkat$UniprotID) |> 
  left_join(RNA_granule_combined, by = "UniprotID") |> 
  mutate(cell = "Jurkat", group = "WP_RBP", group2 = "WP")

#combine
OG_glycoprotein_RNA_granule_combined <- bind_rows(
  OG_glycoprotein_RNA_granule_HEK293T,
  OG_glycoprotein_RNA_granule_HepG2,
  OG_glycoprotein_RNA_granule_Jurkat
)

write_xlsx(OG_glycoprotein_RNA_granule_combined, path = paste0(file_path, "OG_glycoprotein_RNA_granule_combined.xlsx"))

#wilcox test
OG_glycoprotein_RNA_granule_wilcox_test <- OG_glycoprotein_RNA_granule_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ category, p.adjust.method = "BH") |> 
  add_significance("p")

#density plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

density_plot_OG_glycoprotein_RNA_granule <- OG_glycoprotein_RNA_granule_combined |> 
  ggplot() +
  geom_density(aes(x = logFC, color = cell), linewidth = 1) +
  #geom_vline(aes(xintercept = quantile_0.5, color = cell)) +
  scale_color_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  labs(x = expression(Δlog[2]*"FC")) +
  coord_cartesian(xlim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 80, color = "black"),
    legend.position = "right"
  )

ggsave(filename = paste0(file_path, "density_plot_OG_glycoprotein_RNA_granule.png"), 
       plot = density_plot_OG_glycoprotein_RNA_granule, height = 4, width = 5, units = "in", dpi = 600)


#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_RNA_granule <- OG_glycoprotein_RNA_granule_combined |> 
  ggplot() +
  geom_violin(aes(x = factor(category, levels = c("SG", "Both")), y = logFC, fill = cell), color = "transparent") +
  geom_boxplot(aes(x = factor(category, levels = c("SG", "Both")), y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_grid(~ cell, scales = "free_x") +
  stat_pvalue_manual(data = OG_glycoprotein_RNA_granule_wilcox_test, label = "p.signif", tip.length = 0, size = 50,
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
    strip.background = element_rect(fill = "transparent", color = "transparent")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_RNA_granule.png"), plot = violin_boxplot_OG_glycoprotein_RNA_granule, height = 4, width = 6, units = "in", dpi = 600)

#combine
OG_WP_RNA_granule_combined <- bind_rows(
  OG_glycoprotein_RNA_granule_HEK293T,
  WP_protein_RNA_granule_HEK293T,
  OG_glycoprotein_RNA_granule_HepG2,
  WP_protein_RNA_granule_HepG2,
  OG_glycoprotein_RNA_granule_Jurkat,
  WP_protein_RNA_granule_Jurkat
) |> mutate(group2 = ifelse(is.na(group2), "OG", group2))

write_xlsx(OG_WP_RNA_granule_combined, path = paste0(file_path, "OG_WP_RNA_granule_combined.xlsx"))

#wilcox test
OG_WP_RNA_granule_wilcox_test <- OG_WP_RNA_granule_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ group, p.adjust.method = "BH") |> 
  add_significance("p")

#combine
OG_glycoprotein_RNA_granule_not_combined <- bind_rows(
  OG_glycoprotein_RNA_granule_HEK293T,
  OG_glycoprotein_not_RNA_granule_HEK293T,
  OG_glycoprotein_RNA_granule_HepG2,
  OG_glycoprotein_not_RNA_granule_HepG2,
  OG_glycoprotein_RNA_granule_Jurkat,
  OG_glycoprotein_not_RNA_granule_Jurkat
)

#wilcox test
OG_RNA_granule_not_wilcox_test <- OG_glycoprotein_RNA_granule_not_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ category, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_RNA_granule_not <- OG_glycoprotein_RNA_granule_not_combined |> 
  ggplot() +
  geom_violin(aes(x = factor(category, levels = c("SG", "Both", "not_granule")), y = logFC, fill = cell), color = "transparent") +
  geom_boxplot(aes(x = factor(category, levels = c("SG", "Both", "not_granule")), y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_grid(~ cell, scales = "free_x") +
  stat_pvalue_manual(data = OG_RNA_granule_not_wilcox_test, label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.8)) +
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

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_RNA_granule_not.png"), plot = violin_boxplot_OG_glycoprotein_RNA_granule_not, height = 3, width = 6, units = "in", dpi = 600)

#combine
OG_glycoprotein_RNA_granule_not_HEK293T <- bind_rows(
  OG_glycoprotein_RNA_granule_HEK293T,
  OG_glycoprotein_not_RNA_granule_HEK293T
)

OG_site_RNA_granule_HEK293T <- OG_glycopeptide_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_RNA_granule_HEK293T$UniprotID) |> 
  filter(!is.na(combined_site))

write_xlsx(OG_site_RNA_granule_HEK293T, path = paste0(file_path, "OG_site_RNA_granule_HEK293T.xlsx"))

OG_site_not_RNA_granule_HEK293T <- OG_glycopeptide_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_not_RNA_granule_HEK293T$UniprotID) |> 
  filter(!is.na(combined_site))

write_xlsx(OG_site_not_RNA_granule_HEK293T, path = paste0(file_path, "OG_site_not_RNA_granule_HEK293T.xlsx"))

#improt data
uniprot_RNA_granule_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\RNA_granule_HEK293T_uniprot.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)
uniport_not_RNA_granule_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\Not_RNA_granule_HEK293T_uniprot.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)

#generate data frame
OG_site_RNA_granule_13mer_HEK293T <- OG_site_RNA_granule_HEK293T |> 
  left_join(uniprot_RNA_granule_HEK293T, by = join_by(UniprotID == From)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence)

write_xlsx(OG_site_RNA_granule_13mer_HEK293T, path = paste0(file_path, "OG_site_RNA_granule_13mer_HEK293T.xlsx"))

OG_site_not_RNA_granule_13mer_HEK293T <- OG_site_not_RNA_granule_HEK293T |> 
  left_join(uniport_not_RNA_granule_HEK293T, by = join_by(UniprotID == From)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence)

write_xlsx(OG_site_not_RNA_granule_13mer_HEK293T, path = paste0(file_path, "OG_site_not_RNA_granule_13mer_HEK293T.xlsx"))

#import data from plogo
plogo_OG_site_RNA_granule_HEK293T_fix_S <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\plogo_matrix_OG_site_RNA_granule_13mer_HEK293T_fixed_S.txt",
  col_names = TRUE,
  delim = ","
)
plogo_OG_site_RNA_granule_HEK293T_fix_T <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\plogo_matrix_OG_site_RNA_granule_13mer_HEK293T_fixed_T.txt",
  col_names = TRUE,
  delim = ","
)
plogo_OG_site_not_RNA_granule_HEK293T_fix_S <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\plogo_matrix_OG_site_not_RNA_granule_13mer_HEK293T_fixed_S.txt",
  col_names = TRUE,
  delim = ","
)
plogo_OG_site_not_RNA_granule_HEK293T_fix_T <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\plogo_matrix_OG_site_not_RNA_granule_13mer_HEK293T_fixed_T.txt",
  col_names = TRUE,
  delim = ","
)

#heatmap
#RNA granule fixed S
plogo_OG_site_RNA_granule_HEK293T_fix_S_adj <- plogo_OG_site_RNA_granule_HEK293T_fix_S |> 
  mutate(score = ifelse(site == "S_0", NA, score)) |> 
  separate(col = site, into = c('site', 'position'), convert = FALSE, sep = "_") |> 
  arrange(site) |> 
  filter(!site %in% c("B", "U", "X", "Z")) |> 
  pivot_wider(names_from = site, values_from = score) |> 
  filter(position != '0')

mat_RNA_granule_fixed_S <- data.matrix(plogo_OG_site_RNA_granule_HEK293T_fix_S_adj |> select(!position))
rownames(mat_RNA_granule_fixed_S) <- plogo_OG_site_RNA_granule_HEK293T_fix_S_adj$position

col_mat <- colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"))

h1 <- Heatmap(matrix = mat_RNA_granule_fixed_S, col = col_mat, 
        cluster_rows = FALSE, cluster_columns = FALSE,
        rect_gp = gpar(col = "white", lty = 2),
        row_order = plogo_OG_site_RNA_granule_HEK293T_fix_S_adj$position,
        name = 'Enrichment Score', row_names_side = "left")

#RNA granule fixed T
plogo_OG_site_RNA_granule_HEK293T_fix_T_adj <- plogo_OG_site_RNA_granule_HEK293T_fix_T |> 
  mutate(score = ifelse(site == "T_0", NA, score)) |> 
  separate(col = site, into = c('site', 'position'), convert = FALSE, sep = "_") |> 
  arrange(site) |> 
  filter(!site %in% c("B", "U", "X", "Z")) |> 
  pivot_wider(names_from = site, values_from = score) |> 
  filter(position != '0')

mat_RNA_granule_fixed_T <- data.matrix(plogo_OG_site_RNA_granule_HEK293T_fix_T_adj |> select(!position))
rownames(mat_RNA_granule_fixed_T) <- plogo_OG_site_RNA_granule_HEK293T_fix_T_adj$position

col_mat <- colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"))

h2 <- Heatmap(matrix = mat_RNA_granule_fixed_T, col = col_mat, 
              cluster_rows = FALSE, cluster_columns = FALSE,
              rect_gp = gpar(col = "white", lty = 2),
              name = 'Enrichment Score', row_names_side = "left")

#not RNA granule fixed S
plogo_OG_site_not_RNA_granule_HEK293T_fix_S_adj <- plogo_OG_site_not_RNA_granule_HEK293T_fix_S |> 
  mutate(score = ifelse(site == "S_0", NA, score)) |> 
  separate(col = site, into = c('site', 'position'), convert = FALSE, sep = "_") |> 
  arrange(site) |> 
  filter(!site %in% c("B", "U", "X", "Z")) |> 
  pivot_wider(names_from = site, values_from = score) |> 
  filter(position != '0')

mat_not_RNA_granule_fixed_S <- data.matrix(plogo_OG_site_not_RNA_granule_HEK293T_fix_S_adj |> select(!position))
rownames(mat_not_RNA_granule_fixed_S) <- plogo_OG_site_not_RNA_granule_HEK293T_fix_S_adj$position

col_mat <- colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"))

h3 <- Heatmap(matrix = mat_not_RNA_granule_fixed_S, col = col_mat, 
              cluster_rows = FALSE, cluster_columns = FALSE,
              rect_gp = gpar(col = "white", lty = 2),
              name = 'Enrichment Score', row_names_side = "left")

#not RNA granule fixed T
plogo_OG_site_not_RNA_granule_HEK293T_fix_T_adj <- plogo_OG_site_not_RNA_granule_HEK293T_fix_T |> 
  mutate(score = ifelse(site == "T_0", NA, score)) |> 
  separate(col = site, into = c('site', 'position'), convert = FALSE, sep = "_") |> 
  arrange(site) |> 
  filter(!site %in% c("B", "U", "X", "Z")) |> 
  pivot_wider(names_from = site, values_from = score) |> 
  filter(position != '0')

mat_not_RNA_granule_fixed_T <- data.matrix(plogo_OG_site_not_RNA_granule_HEK293T_fix_T_adj |> select(!position))
rownames(mat_not_RNA_granule_fixed_T) <- plogo_OG_site_not_RNA_granule_HEK293T_fix_T_adj$position

col_mat <- colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"))

h4 <- Heatmap(matrix = mat_not_RNA_granule_fixed_T, col = col_mat, 
              cluster_rows = FALSE, cluster_columns = FALSE,
              rect_gp = gpar(col = "white", lty = 2),
              name = 'Enrichment Score', row_names_side = "left")

#import data from plogo 2
plogo_OG_site_RNA_granule_HEK293T_fix_S_2 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\plogo_matrix_OG_site_RNA_granule_13mer_HEK293T_fixed_S_2.txt",
  col_names = TRUE,
  delim = ","
)
plogo_OG_site_RNA_granule_HEK293T_fix_T_2 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\plogo_matrix_OG_site_RNA_granule_13mer_HEK293T_fixed_T_2.txt",
  col_names = TRUE,
  delim = ","
)
plogo_OG_site_not_RNA_granule_HEK293T_fix_S_2 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\plogo_matrix_OG_site_not_RNA_granule_13mer_HEK293T_fixed_S_2.txt",
  col_names = TRUE,
  delim = ","
)
plogo_OG_site_not_RNA_granule_HEK293T_fix_T_2 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\plogo_matrix_OG_site_not_RNA_granule_13mer_HEK293T_fixed_T_2.txt",
  col_names = TRUE,
  delim = ","
)

#heatmap
#RNA granule fixed S_2
plogo_OG_site_RNA_granule_HEK293T_fix_S_2_adj <- plogo_OG_site_RNA_granule_HEK293T_fix_S_2 |> 
  mutate(score = ifelse(site == "S_0", NA, score)) |> 
  separate(col = site, into = c('site', 'position'), convert = FALSE, sep = "_") |> 
  arrange(site) |> 
  filter(!site %in% c("B", "U", "X", "Z")) |> 
  pivot_wider(names_from = site, values_from = score) |> 
  filter(position != '0')

mat_RNA_granule_fixed_S_2 <- data.matrix(plogo_OG_site_RNA_granule_HEK293T_fix_S_2_adj |> select(!position))
rownames(mat_RNA_granule_fixed_S_2) <- plogo_OG_site_RNA_granule_HEK293T_fix_S_2_adj$position

col_mat <- colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"))

h5 <- Heatmap(matrix = mat_RNA_granule_fixed_S_2, col = col_mat, 
              cluster_rows = FALSE, cluster_columns = FALSE,
              rect_gp = gpar(col = "white", lty = 2),
              name = 'Enrichment Score', row_names_side = "left")

#RNA granule fixed T_2
plogo_OG_site_RNA_granule_HEK293T_fix_T_2_adj <- plogo_OG_site_RNA_granule_HEK293T_fix_T_2 |> 
  mutate(score = ifelse(site == "T_0", NA, score)) |> 
  separate(col = site, into = c('site', 'position'), convert = FALSE, sep = "_") |> 
  arrange(site) |> 
  filter(!site %in% c("B", "U", "X", "Z")) |> 
  pivot_wider(names_from = site, values_from = score) |> 
  filter(position != '0')

mat_RNA_granule_fixed_T_2 <- data.matrix(plogo_OG_site_RNA_granule_HEK293T_fix_T_2_adj |> select(!position))
rownames(mat_RNA_granule_fixed_T_2) <- plogo_OG_site_RNA_granule_HEK293T_fix_T_2_adj$position

col_mat <- colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"))

h6 <- Heatmap(matrix = mat_RNA_granule_fixed_T_2, col = col_mat, 
              cluster_rows = FALSE, cluster_columns = FALSE,
              rect_gp = gpar(col = "white", lty = 2),
              name = 'Enrichment Score', row_names_side = "left")

#not RNA granule fixed S_2
plogo_OG_site_not_RNA_granule_HEK293T_fix_S_2_adj <- 
  plogo_OG_site_not_RNA_granule_HEK293T_fix_S_2 |> 
  mutate(score = ifelse(site == "S_0", NA, score)) |> 
  separate(col = site, into = c('site', 'position'), convert = FALSE, sep = "_") |> 
  arrange(site) |> 
  filter(!site %in% c("B", "U", "X", "Z")) |> 
  pivot_wider(names_from = site, values_from = score) |> 
  filter(position != '0')

mat_not_RNA_granule_fixed_S_2 <- data.matrix(plogo_OG_site_not_RNA_granule_HEK293T_fix_S_2_adj |> select(!position))
rownames(mat_not_RNA_granule_fixed_S_2) <- plogo_OG_site_not_RNA_granule_HEK293T_fix_S_2_adj$position

col_mat <- colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"))

h7 <- Heatmap(matrix = mat_not_RNA_granule_fixed_S_2, col = col_mat, 
              cluster_rows = FALSE, cluster_columns = FALSE,
              rect_gp = gpar(col = "white", lty = 2),
              name = 'Enrichment Score', row_names_side = "left")

#not RNA granule fixed T_2
plogo_OG_site_not_RNA_granule_HEK293T_fix_T_2_adj <- 
  plogo_OG_site_not_RNA_granule_HEK293T_fix_T_2 |> 
  mutate(score = ifelse(site == "T_0", NA, score)) |> 
  separate(col = site, into = c('site', 'position'), convert = FALSE, sep = "_") |> 
  arrange(site) |> 
  filter(!site %in% c("B", "U", "X", "Z")) |> 
  pivot_wider(names_from = site, values_from = score) |> 
  filter(position != '0')

mat_not_RNA_granule_fixed_T_2 <- data.matrix(plogo_OG_site_not_RNA_granule_HEK293T_fix_T_2_adj |> select(!position))
rownames(mat_not_RNA_granule_fixed_T_2) <- plogo_OG_site_not_RNA_granule_HEK293T_fix_T_2_adj$position

col_mat <- colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"))

h8 <- Heatmap(matrix = mat_not_RNA_granule_fixed_T_2, col = col_mat, 
              cluster_rows = FALSE, cluster_columns = FALSE,
              rect_gp = gpar(col = "white", lty = 2),
              name = 'Enrichment Score', row_names_side = "left")

#RNA granule glycoprotein overlap
df <- OG_glycoprotein_RNA_granule_combined |> 
  filter(UniprotID %in% uniprot_RNA_granule_overlap$From) |> distinct(UniprotID)

write_xlsx(df, path = paste0(file_path, "OG_glycoprotein_RNA_granule_overlap.xlsx"))

#import uniprot
uniprot_RNA_granule_overlap <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RNA_granule\\uniprot_OG_glycoprotein_RNA_granule_overlap.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
) |> separate(col = Gene.Names, into = "Gene_Name")

OG_glycoprotein_RNA_granule_overlap <- OG_glycoprotein_RNA_granule_combined |> 
  left_join(uniprot_RNA_granule_overlap, by = join_by(UniprotID == From)) |> 
  filter(!is.na(Gene_Name)) |> 
  select(Gene_Name, cell, logFC, Annotation) |> 
  filter(!is.na(Annotation)) |> 
  arrange(Annotation)

write_xlsx(OG_glycoprotein_RNA_granule_overlap, path = paste0(file_path, "OG_glycoprotein_RNA_granule_overlap.xlsx"))

OG_glycoprotein_RNA_granule_overlap_adj <- OG_glycoprotein_RNA_granule_overlap |> 
  select(Gene_Name, cell, logFC) |> 
  pivot_wider(names_from = cell, values_from = logFC)

#heatmap
mat <- data.matrix(OG_glycoprotein_RNA_granule_overlap_adj |> select(!Gene_Name))
rownames(mat) <- OG_glycoprotein_RNA_granule_overlap_adj$Gene_Name
row_split <- OG_glycoprotein_RNA_granule_overlap |> filter(cell == "HEK293T") |> 
  filter(Annotation != "Other") |> pull(Annotation)

mat_col <- colorRamp2(breaks = c(-2, 0, 2), color = c("blue", "white", "red"))

h1 <- Heatmap(matrix = mat, col = mat_col, name = 'Δlog2FC', 
              row_split = row_split, row_names_side = "left", 
              row_title = NULL, border = TRUE, show_row_dend = FALSE)

#circular proportion plot/radical plot
df <- tribble(
  ~ cell, ~ total, ~ Q1, ~ Q2, ~ Q3, ~ Q4,
  "HEK293T", 410, 98, 104, 94, 410/3,
  "HepG2", 193, 19, 27, 94, 193/3,
  "Jurkat", 295, 36, 91, 94, 295/3
)

df_adj <- df |> 
  mutate(Q0 = total - (Q1 + Q2 + Q3)) |> 
  pivot_longer(cols = Q1:Q0, names_to = 'group', values_to = 'num') |> 
  mutate(cell = factor(cell, levels = c('AAA', "HEK293T", "HepG2", "Jurkat")))


font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

proportion_plot <- df_adj |> 
  ggplot(aes(x = cell, y = num, fill = factor(group, levels = c("Q4", "Q0", "Q1", "Q2", "Q3")))) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(aes(label = num, color = factor(group, levels = c("Q4", "Q0", "Q1", "Q2", "Q3"))), 
            position = position_fill(vjust = 0.5), size = 30, show.legend = FALSE) +
  scale_color_manual(
    values = c(
      "Q4" = "transparent",
      "Q1" = "black",
      "Q2" = "black",
      "Q3" = "black",
      "Q0" = "black"
    )) +
  scale_fill_manual(
    values = c(
      "Q4" = "transparent",
      "Q1" = Color_2,
      "Q2" = Color_3,
      "Q3" = Color_4,
      "Q0" = "gray"
    ),
    labels = c(
      "Q1" = "1",
      "Q2" = "2",
      "Q3" = "3",
      "Q0" = "not RNA granule",
      "Q4" = ""
    )
  ) +
  scale_x_discrete(drop = FALSE, labels = c("", "HEK293T", "HepG2", "Jurkat")) +
  scale_y_continuous(breaks = c(0, 3/16, 3/8, 9/16, 3/4, 3/4), labels = c(0, 25, 50, 75, 100, 0)) +
  labs(x = "", y = "Proportion (%)", fill = "Quantified in # cell line") +
  coord_polar(theta = "y") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 70, color = "black"),
    legend.text = element_text(size = 70, color = "black"),
    legend.position = "right"
  )

ggsave(filename = paste0(file_path, "proportion_plot.png"), plot = proportion_plot, height = 6, width = 7, units = "in", dpi = 600)

#circular proportion plot/radical plot
df <- tribble(
  ~ cell, ~ total, ~ SG, ~ Both, ~ Q4,
  "HEK293T", 410, 257, 39, 410/3,
  "HepG2", 193, 123, 17, 193/3,
  "Jurkat", 295, 193, 28, 295/3
)

df_adj <- df |> 
  mutate(Q0 = total - (SG + Both)) |> 
  pivot_longer(cols = SG:Q0, names_to = 'group', values_to = 'num') |> 
  mutate(cell = factor(cell, levels = c('AAA', "HepG2", "Jurkat", "HEK293T")))


font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

proportion_plot <- df_adj |> 
  ggplot(aes(x = cell, y = num, fill = factor(group, levels = c("Q4", "Q0", "SG", "Both")))) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(aes(label = num, color = factor(group, levels = c("Q4", "Q0", "SG", "Both"))), 
            position = position_fill(vjust = 0.5), size = 35, show.legend = FALSE) +
  scale_color_manual(
    values = c(
      "Q4" = "transparent",
      "SG" = "black",
      "Both" = "black",
      "Q0" = "black"
    )) +
  scale_fill_manual(
    values = c(
      "Q4" = "transparent",
      "SG" = Color_1,
      "Both" = Color_2,
      "Q0" = "gray"
    ),
    labels = c(
      "SG" = "SG",
      "Both" = "SG & PB",
      "Q0" = "not RNA granule",
      "Q4" = ""
    )
  ) +
  scale_x_discrete(drop = FALSE, labels = c("", "HepG2", "Jurkat", "HEK293T")) +
  scale_y_continuous(breaks = c(0, 3/16, 3/8, 9/16, 3/4, 3/4), labels = c(0, 25, 50, 75, 100, 0)) +
  labs(x = "", y = "Proportion (%)", fill = "Category") +
  coord_polar(theta = "y") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 120, color = "black"),
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.position = "right"
  )

ggsave(filename = paste0(file_path, "proportion_plot.png"), plot = proportion_plot, height = 6, width = 7, units = "in", dpi = 600)
