#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "showtext", "ggpubr", "rstatix", "ComplexHeatmap", "circlize", "eulerr")
lapply(packages_names, require, character.only = TRUE)

#import human RBP database
human_RBP_database <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\Table_Hs_RBP.txt",
  col_names = TRUE,
  name_repair = "universal",
  skip = 6
)

human_RBP_database_adj <- human_RBP_database |> 
  select(Entry_Name, Uniprot_ID, RBP2GO_Score, Baltz_HEK293_2012, Castello_HeLa.S3_2012, Beckmann_HuH.7_2015, Castello_HeLa.S3_2016,
         Conrad_K562_2016, Milek_MCF7_2017, Perez.Perri_Jurkat_RIC_2018, Perez.Perri_Jurkat_eRIC_2018, Garcia.Moreno_HEK293_2019, Backlund_HuH.7_Cytoplasmic_2020,
         Backlund_HuH.7_Nuclear_2020, Kramer_HeLa_2014, Panhale_HEK293_2019, Mullari_HEK293_2017)

write_xlsx(human_RBP_database_adj, path = paste0(file_path, "human_RBP_database_adj.xlsx"))

#import human RBP database adj
human_RBP_database_adj <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\human_RBP_database_adj.xlsx",
  sheet = 'Sheet2'
)

#generate data frame for RNA binding proteins in each cell line
#HEK293T
OG_glycoprotein_RBP_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(UniprotID %in% human_RBP_database_adj$Uniprot_ID) |> 
  mutate(cell = "HEK293T", group = "OG_RBP_HEK293T")

write_xlsx(OG_glycoprotein_RBP_HEK293T, path = paste0(file_path, "OG_glycoprotein_RBP_HEK293T.xlsx"))

OG_glycoprotein_non_RBP_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(! UniprotID %in% OG_glycoprotein_RBP_HEK293T$UniprotID) |> 
  mutate(cell = "HEK293T", group = "OG_non_RBP_HEK293T")

WP_protein_RBP_HEK293T <- WP_protein_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_RBP_HEK293T$UniprotID) |> 
  mutate(cell = "HEK293T", group = "WP_RBP")

#HepG2
OG_glycoprotein_RBP_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(UniprotID %in% human_RBP_database_adj$Uniprot_ID) |> 
  mutate(cell = "HepG2", group = "OG_RBP_HepG2")

write_xlsx(OG_glycoprotein_RBP_HepG2, path = paste0(file_path, "OG_glycoprotein_RBP_HepG2.xlsx"))

OG_glycoprotein_non_RBP_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(! UniprotID %in% OG_glycoprotein_RBP_HepG2$UniprotID) |> 
  mutate(cell = "HepG2", group = "OG_non_RBP_HepG2")

WP_protein_RBP_HepG2 <- WP_protein_Top_tb_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_RBP_HepG2$UniprotID) |> 
  mutate(cell = "HepG2", group = "WP_RBP")

#Jurkat
OG_glycoprotein_RBP_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(UniprotID %in% human_RBP_database_adj$Uniprot_ID) |> 
  mutate(cell = "Jurkat", group = "OG_RBP_Jurkat")

write_xlsx(OG_glycoprotein_RBP_Jurkat, path = paste0(file_path, "OG_glycoprotein_RBP_Jurkat.xlsx"))

OG_glycoprotein_non_RBP_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(! UniprotID %in% OG_glycoprotein_RBP_Jurkat$UniprotID) |> 
  mutate(cell = "Jurkat", group = "OG_non_RBP_Jurkat")

WP_protein_RBP_Jurkat <- WP_protein_Top_tb_Jurkat |> 
  filter(UniprotID %in% OG_glycoprotein_RBP_Jurkat$UniprotID) |> 
  mutate(cell = "Jurkat", group = "WP_RBP")

#eruler plot for RBP in three cell lines
mat_eruler <- c(
  "HEK293T (412)" = 138,
  "HepG2 (193)" = 31,
  "Jurkat (302)" = 50,
  "HEK293T (412)&HepG2 (193)" = 154,
  "HEK293T (412)&Jurkat (302)" = 224,
  "HepG2 (193)&Jurkat (302)" = 130,
  "HepG2 (193)&HEK293T (412)&Jurkat (302)" = 124
)

fit_eruler <- euler(mat_eruler)

plot(fit_eruler, 
     fill = list(fill = c(Color_2, Color_3, Color_4), alpha = 0.8),
     edges = list(col = "white", lwd = 3),
     labels = list(cex = c(1.2, 1.2, 1.2)),
     legend = list(fontsize = 12),
     check.overlap = TRUE
     )

#proportion plot
df <- tribble(
  ~ cell, ~ total, ~ Q1, ~ Q2, ~ Q3,
  "HEK293T", 1234, 138, 151, 121,
  "HepG2", 432, 31, 41, 121,
  "Jurkat", 850, 48, 126, 121
)

df_adj <- df |> 
  mutate(Q0 = total - (Q1 + Q2 + Q3)) |> 
  pivot_longer(cols = Q1:Q0, names_to = 'group', values_to = 'num')

font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

proportion_plot <- df_adj |> 
  ggplot(aes(x = cell, y = num, fill = factor(group, levels = c("Q0", "Q1", "Q2", "Q3")))) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(aes(label = num), 
            position = position_fill(vjust = 0.5), size = 30) +
  scale_fill_manual(
    values = c(
      "Q1" = Color_2,
      "Q2" = Color_3,
      "Q3" = Color_4,
      "Q0" = "gray"
    ),
    labels = c(
      "Q1" = "1",
      "Q2" = "2",
      "Q3" = "3",
      "Q0" = "not RBP"
    )
  ) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 25, 50, 75, 100)) +
  labs(x = "", y = "Proportion (%)", fill = "Quantified in # cell line") +
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

ggsave(filename = paste0(file_path, "proportion_plot.png"), plot = proportion_plot, height = 4, width = 5, units = "in", dpi = 600)

#OG glycoprotein RBP
OG_glycoprotein_RBP_combined <- bind_rows(
  OG_glycoprotein_RBP_HEK293T,
  OG_glycoprotein_RBP_HepG2,
  OG_glycoprotein_RBP_Jurkat
) |> group_by(cell) |> 
  mutate(
  quantile_0.5 = quantile(logFC, probs = 0.5)
)

#density plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

density_plot_OG_RBP <- OG_glycoprotein_RBP_combined |> 
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
  labs(x = expression(log[2]*"(Tuni/Ctrl)")) +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
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

ggsave(filename = paste0(file_path, "density_plot_OG_RBP.png"), plot = density_plot_OG_RBP, height = 4, width = 5, units = "in", dpi = 600)

#wilcox test
OG_glycoprotein_RBP_wilcox_test <- OG_glycoprotein_RBP_combined |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_RBP <- OG_glycoprotein_RBP_combined |> 
  ggplot() +
  geom_violin(aes(x = cell, y = logFC, fill = cell), color = "transparent") +
  geom_boxplot(aes(x = cell, y = logFC), outliers = FALSE, width = 0.2) +
  labs(x = "", y = expression(Δlog[2]*"FC")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  stat_pvalue_manual(data = OG_glycoprotein_RBP_wilcox_test, label = "p.signif", 
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
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_RBP.png"), plot = violin_boxplot_OG_glycoprotein_RBP, height = 3, width = 3, units = "in", dpi = 600)

#OG glycoprotein RBP non RBP
OG_glycoprotein_RBP_non_RBP_combined <- bind_rows(
  OG_glycoprotein_RBP_HEK293T,
  OG_glycoprotein_non_RBP_HEK293T,
  OG_glycoprotein_RBP_HepG2,
  OG_glycoprotein_non_RBP_HepG2,
  OG_glycoprotein_RBP_Jurkat,
  OG_glycoprotein_non_RBP_Jurkat
)

#wilcox test
OG_glycoprotein_RBP_non_RBP_wilcox_test <- OG_glycoprotein_RBP_non_RBP_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ group, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_RBP_non_RBP <- OG_glycoprotein_RBP_non_RBP_combined |> 
  ggplot(aes(x = factor(group, levels = c("OG_RBP_HEK293T",
                                          "OG_non_RBP_HEK293T",
                                          "OG_RBP_HepG2",
                                          "OG_non_RBP_HepG2",
                                          "OG_RBP_Jurkat",
                                          "OG_non_RBP_Jurkat")), y = logFC)) +
  geom_violin(aes(fill = group), color = "transparent") +
  geom_boxplot(outliers = FALSE, width = 0.2) +
  labs(x = "", y = expression(Δlog[2]*"FC")) +
  scale_fill_manual(
    values = c(
      "OG_RBP_HEK293T" = Color_2,
      "OG_non_RBP_HEK293T" = "gray",
      "OG_RBP_HepG2" = Color_3,
      "OG_non_RBP_HepG2" = "gray",
      "OG_RBP_Jurkat" = Color_4,
      "OG_non_RBP_Jurkat" = "gray"
    )
  ) +
  scale_x_discrete(labels = c("RBP", "not RBP")) +
  stat_pvalue_manual(data = OG_glycoprotein_RBP_non_RBP_wilcox_test, label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.8)) +
  coord_cartesian(ylim = c(-2, 2)) +
  facet_grid(~ cell, scale = "free_x") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.background = element_rect(color = "transparent", fill = "transparent"),
    strip.text = element_text(color = "black", size = 80)
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_RBP_non_RBP.png"), plot = violin_boxplot_OG_glycoprotein_RBP_non_RBP, height = 3, width = 4, units = "in", dpi = 600)

#OG WP RBP
OG_WP_RBP_combined <- bind_rows(
  OG_glycoprotein_RBP_HEK293T,
  WP_protein_RBP_HEK293T,
  OG_glycoprotein_RBP_HepG2,
  WP_protein_RBP_HepG2,
  OG_glycoprotein_RBP_Jurkat,
  WP_protein_RBP_Jurkat
)

write_xlsx(OG_WP_RBP_combined, path = paste0(file_path, "OG_WP_RBP_combined.xlsx"))

#wilcox test
OG_WP_RBP_wilcox_test <- OG_WP_RBP_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ group, p.adjust.method = "BH") |> 
  add_significance("p")

#point plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

point_plot_OG_WP_RBP <- OG_WP_RBP_combined |> 
  ggplot(aes(x = group, y = logFC)) +
  geom_line(aes(group = UniprotID), color = "gray") +
  geom_point(aes(fill = group), shape = 21, size = 2, color = "black") +
  scale_fill_manual(values = c(
    "OG_RBP_HEK293T" = Color_2,
    "OG_RBP_HepG2" = Color_3,
    "OG_RBP_Jurkat" = Color_4,
    "WP_RBP" = "gray"
  )) +
  scale_x_discrete(labels = c("OG", "WP")) +
  stat_pvalue_manual(data = OG_WP_RBP_wilcox_test, label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(2.5)) +
  facet_grid(~ cell, scale = "free_x") +
  coord_cartesian(ylim = c(-3, 3)) +
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
    strip.background = element_rect(fill = "transparent", color = "transparent"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "point_plot_OG_WP_RBP.png"), plot = point_plot_OG_WP_RBP, height = 4, width = 4, units = c("in"), dpi = 600)

#heatmap
#import RBP overlap
OG_glycoprotein_RBP_overlap <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_WP_RBP_combined.xlsx",
  sheet = "overlap"
)

#generate data frame
OG_glycoprotein_RBP_overlap_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_RBP_overlap$UniprotID) |> 
  select(UniprotID, logFC_HEK293T = logFC)

OG_glycoprotein_RBP_overlap_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_RBP_overlap$UniprotID) |> 
  select(UniprotID, logFC_HepG2 = logFC)

OG_glycoprotein_RBP_overlap_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(UniprotID %in% OG_glycoprotein_RBP_overlap$UniprotID) |> 
  select(UniprotID, logFC_Jurkat = logFC)

OG_glycoprotein_RBP_overlap_logFC <- OG_glycoprotein_RBP_overlap_HEK293T |> 
  left_join(OG_glycoprotein_RBP_overlap_HepG2, by = "UniprotID") |> 
  left_join(OG_glycoprotein_RBP_overlap_Jurkat, by = "UniprotID")

#heatmap
mat_OG_RBP <- t(data.matrix(OG_glycoprotein_RBP_overlap_logFC |> select(starts_with("logFC_"))))
rownames(mat_OG_RBP) <- c("HEK293T", "HepG2", "Jurkat")

col_mat <- colorRamp2(breaks = c(-3, 0, 3), color = c("blue", "white", "red"))

Heatmap(mat_OG_RBP, col = col_mat, 
        name = 'log2(Tuni/Ctrl)')
