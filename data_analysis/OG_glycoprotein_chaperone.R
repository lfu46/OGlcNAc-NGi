#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "ggpubr", "rstatix", "showtext", "ComplexHeatmap", "circlize")
lapply(packages_names, require, character.only = TRUE)

#import hsp90 clients database
hsp90_clients <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\hsp90_clients\\hsp90_clients.xlsx"
)

#generate data frame
#HEK293T
OG_glycoprotein_hsp90_clients_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(UniprotID %in% hsp90_clients$Entry) |> 
  mutate(cell = "HEK293T", group = "OG_client_HEK293T")

OG_glycoprotein_hsp90_non_clients_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(! UniprotID %in% hsp90_clients$Entry) |> 
  mutate(cell = "HEK293T", group = "non-client")

WP_protein_hsp90_clients_HEK293T <- WP_protein_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_hsp90_clients_HEK293T$UniprotID) |> 
  mutate(cell = "HEK293T", group = "WP_client")

#HepG2
OG_glycoprotein_hsp90_clients_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(UniprotID %in% hsp90_clients$Entry) |> 
  mutate(cell = "HepG2", group = "OG_client_HepG2")

OG_glycoprotein_hsp90_non_clients_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(! UniprotID %in% hsp90_clients$Entry) |> 
  mutate(cell = "HepG2", group = "non-client")

WP_protein_hsp90_clients_HepG2 <- WP_protein_Top_tb_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_hsp90_clients_HepG2$UniprotID) |> 
  mutate(cell = "HepG2", group = "WP_client")

#Jurkat
OG_glycoprotein_hsp90_clients_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(UniprotID %in% hsp90_clients$Entry) |> 
  mutate(cell = "Jurkat", group = "OG_client_Jurkat")

OG_glycoprotein_hsp90_non_clients_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(! UniprotID %in% hsp90_clients$Entry) |> 
  mutate(cell = "Jurkat", group = "non-client")

WP_protein_hsp90_clients_Jurkat <- WP_protein_Top_tb_Jurkat |> 
  filter(UniprotID %in% OG_glycoprotein_hsp90_clients_Jurkat$UniprotID) |> 
  mutate(cell = "Jurkat", group = "WP_client")


#OG glycoprotein
#hsp90 clients
OG_glycoprotein_hsp90_clients_combined <- bind_rows(
  OG_glycoprotein_hsp90_clients_HEK293T,
  OG_glycoprotein_hsp90_clients_HepG2,
  OG_glycoprotein_hsp90_clients_Jurkat
)

write_xlsx(OG_glycoprotein_hsp90_clients_combined, path = paste0(file_path, "OG_glycoprotein_hsp90_clients_combined.xlsx"))

#wilcox test
OG_glycoprotein_hsp90_clients_wilcox_test <- OG_glycoprotein_hsp90_clients_combined |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_hsp90_clients_combined <- OG_glycoprotein_hsp90_clients_combined |> 
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
  stat_pvalue_manual(data = OG_glycoprotein_hsp90_clients_wilcox_test, label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(2.5, 2.8)) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none"
    )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_hsp90_clients_combined.png"), plot = violin_boxplot_OG_glycoprotein_hsp90_clients_combined, height = 4, width = 4, units = "in", dpi = 600)

#OG glycoprotein
#hsp90 clients vs. non hsp90 clients
OG_glycoprotein_hsp90_non_hsp90_combined <- bind_rows(
  OG_glycoprotein_hsp90_clients_HEK293T,
  OG_glycoprotein_hsp90_non_clients_HEK293T,
  OG_glycoprotein_hsp90_clients_HepG2,
  OG_glycoprotein_hsp90_non_clients_HepG2,
  OG_glycoprotein_hsp90_clients_Jurkat,
  OG_glycoprotein_hsp90_non_clients_Jurkat
)

#wilcox test
OG_glycoprotein_hsp90_non_hsp90_wilcox_test <- OG_glycoprotein_hsp90_non_hsp90_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ group, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_hsp90_non_hsp90_combined <- OG_glycoprotein_hsp90_non_hsp90_combined |> 
  mutate(group_color = ifelse(group == "client", paste(cell, group, sep = "_"), group)) |> 
  ggplot() +
  geom_violin(aes(x = group, y = logFC, fill = group_color), color = "transparent") +
  geom_boxplot(aes(x = group, y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(
    values = c(
      "HEK293T_client" = Color_2,
      "HepG2_client" = Color_3,
      "Jurkat_client" = Color_4,
      "non-client" = "gray"
    )
  ) +
  stat_pvalue_manual(data = OG_glycoprotein_hsp90_non_hsp90_wilcox_test, label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(2.6)) +
  coord_cartesian(ylim = c(-3, 3)) +
  facet_grid(~ cell, scale = "free_x") +
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

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_hsp90_non_hsp90_combined.png"), plot = violin_boxplot_OG_glycoprotein_hsp90_non_hsp90_combined, height = 4, width = 4, units = c("in"), dpi = 600)

#OG vs. WP hsp90 clients
#generate data frame
OG_WP_hsp90_clients_combined <- bind_rows(
  OG_glycoprotein_hsp90_clients_HEK293T,
  WP_protein_hsp90_clients_HEK293T,
  OG_glycoprotein_hsp90_clients_HepG2,
  WP_protein_hsp90_clients_HepG2,
  OG_glycoprotein_hsp90_clients_Jurkat,
  WP_protein_hsp90_clients_Jurkat
)

#wilcox test
OG_WP_hsp90_clients_wilcox_test <- OG_WP_hsp90_clients_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ group, p.adjust.method = "BH") |> 
  add_significance("p")
  
#point plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

point_plot_OG_WP_hsp90_clients <- OG_WP_hsp90_clients_combined |> 
  ggplot(aes(x = group, y = logFC)) +
  geom_line(aes(group = UniprotID), color = "gray") +
  geom_point(aes(fill = group), shape = 21, color = "black") +
  scale_fill_manual(
    values = c(
      "OG_client_HEK293T" = Color_2,
      "OG_client_HepG2" = Color_3,
      "OG_client_Jurkat" = Color_4,
      "WP_client" = "gray"
    )
  ) +
  scale_x_discrete(
    labels = c("OG", "WP")
  ) +
  stat_pvalue_manual(data = OG_WP_hsp90_clients_wilcox_test, label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(2.6)) +
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
    legend.position = "none",
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "point_plot_OG_WP_hsp90_clients.png"), plot = point_plot_OG_WP_hsp90_clients, height = 4, width = 4, units = c("in"), dpi = 600)

#heatmap
hsp90_clients_overlap <- c("Q9UI10", "Q15046", "Q7L0Y3", "Q9BQC3", "P49327")

OG_glycoprotein_hsp90_clients_overlap_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(UniprotID %in% hsp90_clients_overlap) |> 
  select(UniprotID, logFC_HEK293T = logFC)

OG_glycoprotein_hsp90_clients_overlap_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(UniprotID %in% hsp90_clients_overlap) |> 
  select(UniprotID, logFC_HepG2 = logFC)

OG_glycoprotein_hsp90_clients_overlap_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(UniprotID %in% hsp90_clients_overlap) |> 
  select(UniprotID, logFC_Jurkat = logFC)

OG_glycoprotein_hsp90_clients_overlap <- OG_glycoprotein_hsp90_clients_overlap_HEK293T |> 
  left_join(OG_glycoprotein_hsp90_clients_overlap_HepG2, by = "UniprotID") |> 
  left_join(OG_glycoprotein_hsp90_clients_overlap_Jurkat, by = "UniprotID")

#heatmap
mat_col <- colorRamp2(breaks = c(-1.5, 0, 1.5), colors = c("blue", "white", "red"))
mat <- data.matrix(OG_glycoprotein_hsp90_clients_overlap |> select(starts_with("logFC_")))
rownames(mat) <- c("EIF2B4", "KARS1", "TRMT10C", "DPH2", "FASN")
colnames(mat) <- c("HEK293T", "HepG2", "Jurkat")

Heatmap(mat, col = mat_col, name = 'log2(Tuni/Ctrl)')

#import data
uniprot_protein_folding_chaperone <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\protein_folding_chaperone\\uniprot_protein_folding_chaperone.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)

#generate data frame
#HEK293T
OG_glycoprotein_chaperone_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(UniprotID %in% uniprot_protein_folding_chaperone$Entry) |> 
  mutate(cell = "HEK293T", group = "OG_chaperone")

OG_glycoprotein_not_chaperone_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(!UniprotID %in% OG_glycoprotein_chaperone_HEK293T$UniprotID) |> 
  mutate(cell = "HEK293T", group = "OG_not_chaperone")

WP_protein_chaperone_HEK293T <- WP_protein_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_chaperone_HEK293T$UniprotID) |> 
  mutate(cell = "HEK293T", group = "WP_chaperone")

#HepG2
OG_glycoprotein_chaperone_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(UniprotID %in% uniprot_protein_folding_chaperone$Entry) |> 
  mutate(cell = "HepG2", group = "OG_chaperone")

OG_glycoprotein_not_chaperone_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(!UniprotID %in% OG_glycoprotein_chaperone_HepG2$UniprotID) |> 
  mutate(cell = "HepG2", group = "OG_not_chaperone")

WP_protein_chaperone_HepG2 <- WP_protein_Top_tb_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_chaperone_HepG2$UniprotID) |> 
  mutate(cell = "HepG2", group = "WP_chaperone")

#Jurkat
OG_glycoprotein_chaperone_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(UniprotID %in% uniprot_protein_folding_chaperone$Entry) |> 
  mutate(cell = "Jurkat", group = "OG_chaperone")

OG_glycoprotein_not_chaperone_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(!UniprotID %in% OG_glycoprotein_chaperone_Jurkat$UniprotID) |> 
  mutate(cell = "Jurkat", group = "OG_not_chaperone")

WP_protein_chaperone_Jurkat <- WP_protein_Top_tb_Jurkat |> 
  filter(UniprotID %in% OG_glycoprotein_chaperone_Jurkat$UniprotID) |> 
  mutate(cell = "Jurkat", group = "WP_chaperone")


#combine
OG_glycoprotein_chaperone_combined <- bind_rows(
  OG_glycoprotein_chaperone_HEK293T,
  OG_glycoprotein_chaperone_HepG2,
  OG_glycoprotein_chaperone_Jurkat
)

#wilcox test
OG_glycoprotein_chaperone_wilcox_test <- OG_glycoprotein_chaperone_combined |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#combine
OG_WP_chaperone_combined <- bind_rows(
  OG_glycoprotein_chaperone_HEK293T,
  WP_protein_chaperone_HEK293T,
  OG_glycoprotein_chaperone_HepG2,
  WP_protein_chaperone_HepG2,
  OG_glycoprotein_chaperone_Jurkat,
  WP_protein_chaperone_Jurkat
)

#wilcox test
OG_WP_chaperone_wilcox_test <- OG_WP_chaperone_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ group, p.adjust.method = "BH") |> 
  add_significance("p")

#combine
OG_glycoprotein_chaperone_not_combined <- bind_rows(
  OG_glycoprotein_chaperone_HEK293T,
  OG_glycoprotein_not_chaperone_HEK293T,
  OG_glycoprotein_chaperone_HepG2,
  OG_glycoprotein_not_chaperone_HepG2,
  OG_glycoprotein_chaperone_Jurkat,
  OG_glycoprotein_not_chaperone_Jurkat
)

#wilcox test
OG_glycoprotein_chaperone_not_wilcox_test <- OG_glycoprotein_chaperone_not_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ group, p.adjust.method = "BH") |> 
  add_significance("p")
