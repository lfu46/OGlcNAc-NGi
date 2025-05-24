#import packages
packages_names <- c("tidyverse", "rstatix", "showtext", "ggpubr", "readxl", "writexl", "ComplexHeatmap", "circlize")
lapply(packages_names, require, character.only = TRUE)

#site level
#generate OG site glycoprotein
OG_site_glycoprotein_HepG2 <- OG_glycopeptide_Top_tb_HepG2_singlesite |> distinct(UniprotID)
OG_site_glycoprotein_HEK293T <- OG_glycopeptide_Top_tb_HEK293T_singlesite |> distinct(UniprotID)
OG_site_glycoprotein_Jurkat <- OG_glycopeptide_Top_tb_Jurkat_singlesite |> distinct(UniprotID)

#site level
#combine
OG_site_glycoprotein_combined <- bind_rows(
  OG_site_glycoprotein_HepG2,
  OG_site_glycoprotein_HEK293T,
  OG_site_glycoprotein_Jurkat
) |> distinct()

write_xlsx(OG_site_glycoprotein_combined, path = paste0(file_path, "OG_site_glycoprotein_combined.xlsx"))

#import uniprot sequence data
OG_glycoprotein_uniprot_sequence <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\uniprot_sequence\\OG_site_glycoprotein_sequence.tsv",
  delim = "\t",
  col_names = TRUE,
  name_repair = "universal"
)

#site level
#generate data frame for 13mer 
#HepG2
OG_site_Mer13_HepG2 <- OG_glycopeptide_Top_tb_HepG2_singlesite |> 
  mutate(category = ifelse(logFC > 0.5 & adj.P.Val < 0.05, "up", ifelse(logFC < -0.5 & adj.P.Val < 0.05, "down", "median"))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence, category) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence) |> 
  separate(Mer13, into = paste0("position_", -7:6), sep = "", remove = FALSE) |> 
  select(!'position_-7') |> 
  pivot_longer(cols = starts_with("position"), names_to = "position", values_to = "Amino_Acid") |> 
  mutate(cell = "HepG2")

write_xlsx(OG_site_Mer13_HepG2, path = paste0(file_path, "OG_site_Mer13_HepG2.xlsx"))

#site level
#generate data frame for 13mer 
#HEK293T
OG_site_Mer13_HEK293T <- OG_glycopeptide_Top_tb_HEK293T_singlesite |> 
  mutate(category = ifelse(logFC > 0.5 & adj.P.Val < 0.05, "up", ifelse(logFC < -0.5 & adj.P.Val < 0.05, "down", "median"))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence, category) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence) |> 
  separate(Mer13, into = paste0("position_", -7:6), sep = "", remove = FALSE) |> 
  select(!'position_-7') |> 
  pivot_longer(cols = starts_with("position"), names_to = "position", values_to = "Amino_Acid") |> 
  mutate(cell = "HEK293T")

write_xlsx(OG_site_Mer13_HEK293T, path = paste0(file_path, "OG_site_Mer13_HEK293T.xlsx"))

#site level
#generate data frame for 13mer 
#Jurkat
OG_site_Mer13_Jurkat <- OG_glycopeptide_Top_tb_Jurkat_singlesite |> 
  mutate(category = ifelse(logFC > 0.5 & adj.P.Val < 0.05, "up", ifelse(logFC < -0.5 & adj.P.Val < 0.05, "down", "median"))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence, category) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence) |> 
  separate(Mer13, into = paste0("position_", -7:6), sep = "", remove = FALSE) |> 
  select(!'position_-7') |> 
  pivot_longer(cols = starts_with("position"), names_to = "position", values_to = "Amino_Acid") |> 
  mutate(cell = "Jurkat")

write_xlsx(OG_site_Mer13_Jurkat, path = paste0(file_path, "OG_site_Mer13_Jurkat.xlsx"))

#combine
OG_site_Mer13_combined <- bind_rows(
  OG_site_Mer13_HepG2,
  OG_site_Mer13_HEK293T,
  OG_site_Mer13_Jurkat
)

#site level
#acidic residue D,E
OG_site_Mer13_combined_acidic <- OG_site_Mer13_combined |> 
  filter(Amino_Acid %in% c("D", "E")) |> 
  group_by(cell, Index, category) |> 
  count(Amino_Acid) |> 
  ungroup()

#site level
#wilcox test
OG_site_Mer13_combined_wilcox_test <- OG_site_Mer13_combined_acidic |> 
  group_by(cell) |> 
  wilcox_test(n ~ category, p.adjust.method = "BH") |> 
  add_significance("p")

#site level
#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_site_Mer13_combined_acidic <- OG_site_Mer13_combined_acidic |>
  ggplot() +
  geom_violin(aes(x = factor(category, levels = c("up", "median", "down")), y = n, fill = cell),
              draw_quantiles = c(0.5), adjust = 1.2) +
  #geom_boxplot(aes(x = factor(category, levels = c("up", "median", "down")), y = n), 
  #            color = "black", outliers = FALSE, width = 0.2) +
  scale_fill_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  labs(x = "", y = "# of acidic residue") +
  stat_pvalue_manual(data = OG_site_Mer13_combined_wilcox_test, label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(3.8)) +
  facet_grid(~ cell, scale = "free_x") +
  coord_cartesian(ylim = c(0, 4)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.background = element_rect(fill = "transparent", color = "transparent"),
    strip.text = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_site_Mer13_combined_acidic.png"), plot = violin_boxplot_OG_site_Mer13_combined_acidic, height = 3, width = 5, units = "in", dpi = 600)

#site level
#polar residue S, T, N & Q
OG_site_Mer13_combined_polar <- OG_site_Mer13_combined |> 
  filter(Amino_Acid %in% c("S", "T", "N", "Q")) |> 
  group_by(cell, Index, category) |> 
  count(Amino_Acid) |> 
  ungroup()

#site level
#wilcox test
OG_site_Mer13_combined_polar_wilcox_test <- OG_site_Mer13_combined_polar |> 
  group_by(cell) |> 
  wilcox_test(n ~ category, p.adjust.method = "BH") |> 
  add_significance("p")

#site level
#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_site_Mer13_combined_polar <- OG_site_Mer13_combined_polar |>
  ggplot() +
  geom_violin(aes(x = factor(category, levels = c("up", "median", "down")), y = n, fill = cell),
              draw_quantiles = c(0.5), adjust = 0.9) +
  #geom_boxplot(aes(x = factor(category, levels = c("up", "median", "down")), y = n), 
  #             color = "black", outliers = FALSE, width = 0.2) +
  scale_fill_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  labs(x = "", y = "# of polar residue") +
  stat_pvalue_manual(data = OG_site_Mer13_combined_polar_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(5.3, 5.8)) +
  facet_grid(~ cell, scale = "free_x") +
  coord_cartesian(ylim = c(0, 6)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.background = element_rect(fill = "transparent", color = "transparent"),
    strip.text = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_site_Mer13_combined_polar.png"), plot = violin_boxplot_OG_site_Mer13_combined_polar, height = 3, width = 5, units = "in", dpi = 600)

#site level
#percentage basic residue
OG_site_Mer13_combined_basic <- OG_site_Mer13_combined |> 
  filter(Amino_Acid %in% c("F", "W", "Y")) |> 
  group_by(cell, Index, category) |> 
  count() |> 
  mutate(percentage = n/13) |> 
  ungroup() |> 
  group_by(cell, category) |> 
  summarize(avg_percentage = mean(percentage))

check_df |> 
  group_by(cell) |> 
  wilcox_test(n ~ category)

#site level
Mer13_HepG2 <- OG_site_Mer13_combined |> 
  filter(cell == "HepG2") |> 
  distinct(Mer13, category)

write_xlsx(Mer13_HepG2, path = paste0(file_path, "Mer13_HepG2.xlsx"))

Mer13_HEK293T <- OG_site_Mer13_combined |> 
  filter(cell == "HEK293T") |> 
  distinct(Mer13, category)

write_xlsx(Mer13_HEK293T, path = paste0(file_path, "Mer13_HEK293T.xlsx"))

Mer13_Jurkat <- OG_site_Mer13_combined |> 
  filter(cell == "Jurkat") |> 
  distinct(Mer13, category)

write_xlsx(Mer13_Jurkat, path = paste0(file_path, "Mer13_Jurkat.xlsx"))

#import data from plogo
#HEK293T
plogo_OG_site_up_HEK293T_fixed_S <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\figure8\\plogo_OG_site_13mer_up_HEK293T_fixed_S.txt",
  col_names = TRUE,
  delim = ","
)
plogo_OG_site_up_HEK293T_fixed_T <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\figure8\\plogo_OG_site_13mer_up_HEK293T_fixed_T.txt",
  col_names = TRUE,
  delim = ","
)
#HepG2
plogo_OG_site_up_HepG2_fixed_S <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\figure8\\plogo_OG_site_13mer_up_HepG2_fixed_S.txt",
  col_names = TRUE,
  delim = ","
)
plogo_OG_site_up_HepG2_fixed_T <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\figure8\\plogo_OG_site_13mer_up_HepG2_fixed_T.txt",
  col_names = TRUE,
  delim = ","
)
#Jurkat
plogo_OG_site_up_Jurkat_fixed_S <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\figure8\\plogo_OG_site_13mer_up_Jurkat_fixed_S.txt",
  col_names = TRUE,
  delim = ","
)
plogo_OG_site_up_Jurkat_fixed_T <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\figure8\\plogo_OG_site_13mer_up_Jurkat_fixed_T.txt",
  col_names = TRUE,
  delim = ","
)

#heatmap
#HEK293T up fixed S
plogo_OG_site_up_HEK293T_fixed_S_adj <- plogo_OG_site_up_HEK293T_fixed_S |> 
  mutate(score = ifelse(site == "S_0", NA, score)) |> 
  separate(col = site, into = c('site', 'position'), convert = FALSE, sep = "_") |> 
  arrange(site) |> 
  filter(!site %in% c("B", "U", "X", "Z")) |> 
  pivot_wider(names_from = site, values_from = score) |> 
  filter(position != '0')

mat_OG_site_up_HEK293T_fixed_S <- data.matrix(plogo_OG_site_up_HEK293T_fixed_S_adj |> select(!position))
rownames(mat_OG_site_up_HEK293T_fixed_S) <- plogo_OG_site_up_HEK293T_fixed_S_adj$position

col_mat <- colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"))

h1 <- Heatmap(matrix = mat_OG_site_up_HEK293T_fixed_S, col = col_mat, 
              cluster_rows = FALSE, cluster_columns = FALSE,
              rect_gp = gpar(col = "white", lty = 2),
              name = 'Enrichment Score', row_names_side = "left")

#HEK293T up fixed T
plogo_OG_site_up_HEK293T_fixed_T_adj <- plogo_OG_site_up_HEK293T_fixed_T |> 
  mutate(score = ifelse(site == "T_0", NA, score)) |> 
  separate(col = site, into = c('site', 'position'), convert = FALSE, sep = "_") |> 
  arrange(site) |> 
  filter(!site %in% c("B", "U", "X", "Z")) |> 
  pivot_wider(names_from = site, values_from = score) |> 
  filter(position != '0')

mat_OG_site_up_HEK293T_fixed_T <- data.matrix(plogo_OG_site_up_HEK293T_fixed_T_adj |> select(!position))
rownames(mat_OG_site_up_HEK293T_fixed_T) <- plogo_OG_site_up_HEK293T_fixed_T_adj$position

col_mat <- colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"))

h2 <- Heatmap(matrix = mat_OG_site_up_HEK293T_fixed_T, col = col_mat, 
              cluster_rows = FALSE, cluster_columns = FALSE,
              rect_gp = gpar(col = "white", lty = 2),
              name = 'Enrichment Score', row_names_side = "left")

#HepG2 up fixed S
plogo_OG_site_up_HepG2_fixed_S_adj <- plogo_OG_site_up_HepG2_fixed_S |> 
  mutate(score = ifelse(site == "S_0", NA, score)) |> 
  separate(col = site, into = c('site', 'position'), convert = FALSE, sep = "_") |> 
  arrange(site) |> 
  filter(!site %in% c("B", "U", "X", "Z")) |> 
  pivot_wider(names_from = site, values_from = score) |> 
  filter(position != '0')

mat_OG_site_up_HepG2_fixed_S <- data.matrix(plogo_OG_site_up_HepG2_fixed_S_adj |> select(!position))
rownames(mat_OG_site_up_HepG2_fixed_S) <- plogo_OG_site_up_HepG2_fixed_S_adj$position

col_mat <- colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"))

h3 <- Heatmap(matrix = mat_OG_site_up_HepG2_fixed_S, col = col_mat, 
              cluster_rows = FALSE, cluster_columns = FALSE,
              rect_gp = gpar(col = "white", lty = 2),
              name = 'Enrichment Score', row_names_side = "left")

#HepG2 up fixed T
plogo_OG_site_up_HepG2_fixed_T_adj <- plogo_OG_site_up_HepG2_fixed_T |> 
  mutate(score = ifelse(site == "T_0", NA, score)) |> 
  separate(col = site, into = c('site', 'position'), convert = FALSE, sep = "_") |> 
  arrange(site) |> 
  filter(!site %in% c("B", "U", "X", "Z")) |> 
  pivot_wider(names_from = site, values_from = score) |> 
  filter(position != '0')

mat_OG_site_up_HepG2_fixed_T <- data.matrix(plogo_OG_site_up_HepG2_fixed_T_adj |> select(!position))
rownames(mat_OG_site_up_HepG2_fixed_T) <- plogo_OG_site_up_HepG2_fixed_T_adj$position

col_mat <- colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"))

h4 <- Heatmap(matrix = mat_OG_site_up_HepG2_fixed_T, col = col_mat, 
              cluster_rows = FALSE, cluster_columns = FALSE,
              rect_gp = gpar(col = "white", lty = 2),
              name = 'Enrichment Score', row_names_side = "left")

#Jurkat up fixed S
plogo_OG_site_up_Jurkat_fixed_S_adj <- plogo_OG_site_up_Jurkat_fixed_S |> 
  mutate(score = ifelse(site == "S_0", NA, score)) |> 
  separate(col = site, into = c('site', 'position'), convert = FALSE, sep = "_") |> 
  arrange(site) |> 
  filter(!site %in% c("B", "U", "X", "Z")) |> 
  pivot_wider(names_from = site, values_from = score) |> 
  filter(position != '0')

mat_OG_site_up_Jurkat_fixed_S <- data.matrix(plogo_OG_site_up_Jurkat_fixed_S_adj |> select(!position))
rownames(mat_OG_site_up_Jurkat_fixed_S) <- plogo_OG_site_up_Jurkat_fixed_S_adj$position

col_mat <- colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"))

h5 <- Heatmap(matrix = mat_OG_site_up_Jurkat_fixed_S, col = col_mat, 
              cluster_rows = FALSE, cluster_columns = FALSE,
              rect_gp = gpar(col = "white", lty = 2),
              name = 'Enrichment Score', row_names_side = "left")

#Jurkat up fixed T
plogo_OG_site_up_Jurkat_fixed_T_adj <- plogo_OG_site_up_Jurkat_fixed_T |> 
  mutate(score = ifelse(site == "T_0", NA, score)) |> 
  separate(col = site, into = c('site', 'position'), convert = FALSE, sep = "_") |> 
  arrange(site) |> 
  filter(!site %in% c("B", "U", "X", "Z")) |> 
  pivot_wider(names_from = site, values_from = score) |> 
  filter(position != '0')

mat_OG_site_up_Jurkat_fixed_T <- data.matrix(plogo_OG_site_up_Jurkat_fixed_T_adj |> select(!position))
rownames(mat_OG_site_up_Jurkat_fixed_T) <- plogo_OG_site_up_Jurkat_fixed_T_adj$position

col_mat <- colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"))

h6 <- Heatmap(matrix = mat_OG_site_up_Jurkat_fixed_T, col = col_mat, 
              cluster_rows = FALSE, cluster_columns = FALSE,
              rect_gp = gpar(col = "white", lty = 2),
              name = 'Enrichment Score', row_names_side = "left")
