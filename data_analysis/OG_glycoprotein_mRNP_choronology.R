#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "showtext", "ggpubr", "rstatix")
lapply(packages_names, require, character.only = TRUE)

#import mRNA chronology data base
mRNP_chronology <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\mRNP_chronology\\mRNP_chronology_hsa_HeLa.csv",
  col_names = TRUE,
  name_repair = "universal"
)

mRNP_chronology_cluster <- mRNP_chronology |> 
  select(UniprotID = MasterAccession, cluster)

#generate data frame
#HEK293T
OG_glycoprotein_mRNP_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(UniprotID %in% mRNP_chronology_cluster$UniprotID) |> 
  left_join(mRNP_chronology_cluster, by = "UniprotID") |> 
  filter(!is.na(cluster)) |> 
  mutate(cell = "HEK293T", group = "OG")

WP_protein_mRNP_HEK293T <- WP_protein_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_mRNP_HEK293T$UniprotID) |> 
  left_join(mRNP_chronology_cluster, by = "UniprotID") |> 
  mutate(cell = "HEK293T", group = "WP")

#HepG2
OG_glycoprotein_mRNP_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(UniprotID %in% mRNP_chronology_cluster$UniprotID) |> 
  left_join(mRNP_chronology_cluster, by = "UniprotID") |> 
  filter(!is.na(cluster)) |> 
  mutate(cell = "HepG2", group = "OG")

WP_protein_mRNP_HepG2 <- WP_protein_Top_tb_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_mRNP_HepG2$UniprotID) |> 
  left_join(mRNP_chronology_cluster, by = "UniprotID") |> 
  mutate(cell = "HepG2", group = "WP")

#Jurkat
OG_glycoprotein_mRNP_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(UniprotID %in% mRNP_chronology_cluster$UniprotID) |> 
  left_join(mRNP_chronology_cluster, by = "UniprotID") |> 
  filter(!is.na(cluster)) |> 
  mutate(cell = "Jurkat", group = "OG")

WP_protein_mRNP_Jurkat <- WP_protein_Top_tb_Jurkat |> 
  filter(UniprotID %in% OG_glycoprotein_mRNP_Jurkat$UniprotID) |> 
  left_join(mRNP_chronology_cluster, by = "UniprotID") |> 
  mutate(cell = "Jurkat", group = "WP")

#combine
OG_glycoprotein_mRNP_combined <- bind_rows(
  OG_glycoprotein_mRNP_HEK293T,
  OG_glycoprotein_mRNP_HepG2,
  OG_glycoprotein_mRNP_Jurkat
)

write_xlsx(OG_glycoprotein_mRNP_combined, path = paste0(file_path, "OG_glycoprotein_mRNP_combined.xlsx"))

#wilcox test
OG_glycoprotein_mRNP_wilcox_test <- OG_glycoprotein_mRNP_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ cluster, p.adjust.method = "BH") |> 
  add_significance("p") |> 
  filter(p.signif != "ns") |> 
  filter(group2 == "VII")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_mRNP <- OG_glycoprotein_mRNP_combined |> 
  ggplot() +
  geom_violin(aes(x = cluster, y = logFC, fill = cell), color = "transparent") +
  geom_boxplot(aes(x = cluster, y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_grid(~ cell, scales = "free_x") +
  stat_pvalue_manual(data = OG_glycoprotein_mRNP_wilcox_test , label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(2.0, 2.4, 2.8)) +
  coord_cartesian(ylim = c(-3, 3)) +
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

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_mRNP.png"), plot = violin_boxplot_OG_glycoprotein_mRNP, height = 4, width = 8, units = "in", dpi = 600)

#combine
OG_WP_mRNP_combined <- bind_rows(
  OG_glycoprotein_mRNP_HEK293T,
  WP_protein_mRNP_HEK293T,
  OG_glycoprotein_mRNP_HepG2,
  WP_protein_mRNP_HepG2,
  OG_glycoprotein_mRNP_Jurkat,
  WP_protein_mRNP_Jurkat
)

#wilcox test
OG_WP_mRNP_wilcox_test <- OG_WP_mRNP_combined |> 
  group_by(cell, cluster) |>
  wilcox_test(logFC ~ group, p.adjust.method = "BH") |> 
  add_significance("p") |> 
  filter(!p.signif == "ns")

#point boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

point_boxplot_OG_glycoprotein_mRNP_chronology_HEK293T <- OG_WP_mRNP_combined |> 
  filter(cell == "HEK293T", cluster %in% c("IV", "V", "VI", "VII")) |> 
  ggplot(aes(x = factor(group, levels = c("WP", "OG")), y = logFC)) +
  geom_line(aes(group = UniprotID), color = "gray") +
  geom_boxplot(outliers = FALSE) +
  geom_point(aes(fill = group), shape = 21, size = 2, color = "black") +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(
    values = c(
      "WP" = "gray",
      "OG" = Color_2
    )
  ) +
  facet_grid(~ cluster, scales = "free_x") +
  stat_pvalue_manual(data = OG_WP_mRNP_wilcox_test , label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.0)) +
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

ggsave(filename = paste0(file_path, "point_boxplot_OG_glycoprotein_mRNP_chronology_HEK293T.png"), plot = point_boxplot_OG_glycoprotein_mRNP_chronology_HEK293T, height = 3, width = 6, units = "in", dpi = 600)
