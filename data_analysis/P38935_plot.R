#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "showtext", "ggpubr", "rstatix", 
                    "ComplexHeatmap", "circlize", "eulerr", "ggridges")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
OG_glycoprotein_RBD_logFC_AAA_HEK293T_adj <- OG_glycoprotein_RBD_logFC_AAA_HEK293T |> 
  select(UniprotID = UNIPROT_ACCESSION, Gene_Name = GENE.NAME) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_HEK293T, by = "UniprotID") |> 
  select(UniprotID, Gene_Name, Tuni_1_sl_tmm:Ctrl_6_sl_tmm) |> 
  mutate(
    log2_ratio1 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    log2_ratio2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    log2_ratio3 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  )

write_xlsx(OG_glycoprotein_RBD_logFC_AAA_HEK293T_adj, path = paste0(file_path, "OG_glycoprotein_RBD_logFC_AAA_HEK293T_adj.xlsx"))

#heatmap
mat <- data.matrix(OG_glycoprotein_RBD_logFC_AAA_HEK293T_adj |> select(starts_with("scaled")))
rownames(mat) <- OG_glycoprotein_RBD_logFC_AAA_HEK293T_adj$Gene_Name
colnames(mat) <- c("Tuni 1", "Tuni 2", "Tuni 3", "Ctrl 1", "Ctrl 2", "Ctrl 3")

mat_col <- colorRamp2(breaks = c(-1, 0, 1), color = c("blue", "white", "red"))

h1 <- Heatmap(mat, col = mat_col, 
              cluster_columns = FALSE,
              show_column_names = TRUE, name = 'z-scale intensity')

#protein plot
P38935_glycoprotein_data <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(UniprotID == "P38935")

P38935_glycopeptide_data <- OG_glycopeptide_Top_tb_HEK293T |> 
  filter(UniprotID == "P38935")

P38935_protein_data <- WP_protein_Top_tb_HEK293T |> 
  filter(UniprotID == "P38935")

P38935_domain_data <- tribble(
  ~ domain, ~ start_position, ~ end_position, 
  'Helicase SMUBP-2/HCS1, 1B', 10, 139,
  'AAA', 192, 410,
  'AAA', 419, 614,
  'R3H', 732, 784,
  'AN1-like Zinc finger', 897, 935
)

P38935_phospho_data <- tribble(
  ~ phospho_site,
  160,
  218,
  289,
  539
)

P38935_rect_data <- tribble(
  ~ start, ~ end, ~ top, ~ bottom,
  0, 993, 0.39, 0.42
)

P38935_protein_text_data <- tribble(
  ~ text, ~ x_position, ~ y_position,
  "Unmod.protein", 800, -0.06,
  "Glyco.protein", 800, 0.3
)

font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

P38935_plot <- ggplot() +
  xlim(c(0, 993)) +
  ylim(c(-0.1, 0.53)) +
  geom_point(data = P38935_glycopeptide_data, aes(x = site_position, y = logFC), size = 2) +
  geom_hline(data = P38935_glycoprotein_data, aes(yintercept = logFC), linetype = "dashed", color = Color_2, linewidth = 1) +
  geom_hline(data = P38935_protein_data, aes(yintercept = logFC), linetype = "dashed", color = "black", linewidth = 1) +
  geom_rect(data = P38935_rect_data, aes(xmin = start, xmax = end, ymin = bottom, ymax = top), fill = "gray", color = "black", alpha = 0.5, inherit.aes = FALSE) +
  geom_rect(data = P38935_domain_data, aes(xmin = start_position, xmax = end_position, ymin = P38935_rect_data$bottom, ymax = P38935_rect_data$top, fill = domain), alpha = 0.5) +
  geom_point(data = P38935_glycopeptide_data, aes(x = site_position, y = 0.44), shape = 23, size = 2, fill = Color_2, color = "black") +
  geom_point(data = P38935_phospho_data, aes(x = phospho_site, y = 0.405), shape = 21, size = 2, fill = "yellow", color = "black") +
  geom_text(data = P38935_glycopeptide_data, aes(x = site_position, y = 0.50, label = combined_site), size = 30, angle = 90) +
  geom_text(data = P38935_protein_text_data, aes(x = x_position, y = y_position, label = text, color = text), size = 30, show.legend = FALSE) +
  scale_color_manual(values = c(
    "Unmod.protein" = "black",
    "Glyco.protein" = Color_2
  )) +
  labs(x = "Amino acid position", y = expression(log[2]*"(Tuni/Ctrl)", fill = "Domain")) +
  scale_fill_manual(
    name = "Domain",
    values = c(
      "Helicase SMUBP-2/HCS1, 1B" = Color_1,
      "AAA" = Color_3,
      "R3H" = Color_4,
      "AN1-like Zinc finger" = Color_5
    )
  ) +
  theme_bw() +
  theme(
    title = element_text(size = 100, color = "black"),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    axis.title.x = element_text(size = 100),
    axis.title.y = element_text(size = 100),
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black")
  )

ggsave(paste0(file_path, "P38935_plot.png"), plot = P38935_plot, height = 4, width = 7, dpi = 600)

#protein plot sequence only
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

P38935_plot_sequence_only <- ggplot() +
  xlim(c(0, 993)) +
  ylim(c(-0.1, 0.53)) +
  geom_rect(data = P38935_rect_data, aes(xmin = start, xmax = end, ymin = bottom, ymax = top), fill = "gray", color = "black", alpha = 0.5, inherit.aes = FALSE) +
  geom_rect(data = P38935_domain_data, aes(xmin = start_position, xmax = end_position, ymin = P38935_rect_data$bottom, ymax = P38935_rect_data$top, fill = domain), alpha = 0.5) +
  geom_point(data = P38935_glycopeptide_data, aes(x = site_position, y = 0.44), shape = 23, size = 2, fill = Color_2, color = "black") +
  geom_point(data = P38935_phospho_data, aes(x = phospho_site, y = 0.405), shape = 21, size = 2, fill = "yellow", color = "black") +
  geom_text(data = P38935_glycopeptide_data, aes(x = site_position, y = 0.50, label = combined_site), size = 30, angle = 90) +
  labs(x = "Amino acid position", y = expression(log[2]*"(Tuni/Ctrl)", fill = "Domain")) +
  scale_fill_manual(
    name = "Domain",
    values = c(
      "Helicase SMUBP-2/HCS1, 1B" = Color_1,
      "AAA" = Color_3,
      "R3H" = Color_4,
      "AN1-like Zinc finger" = Color_5
    )
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  )

ggsave(paste0(file_path, "P38935_plot_sequence_only.png"), plot = P38935_plot_sequence_only, height = 3.8, width = 3.8, dpi = 600)



#NCPR
#import data
P38935_NCPR_data <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\P38935_plot\\NCPR.xlsx",
  col_names = c('position', 'NCPR')
)

#bar plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

bar_plot_P38935_NCPR <- P38935_NCPR_data |> 
  mutate(group = ifelse(NCPR > 0, "high", ifelse(NCPR < 0, "low", "none"))) |> 
  ggplot() +
  geom_bar(aes(x = position, y = NCPR, fill = group), stat = "identity", color = "transparent") +
  scale_fill_manual(
    values = c(
      "high" = Color_3,
      "low" = Color_4
    )
  ) +
  labs(x = "Amino acid position") +
  coord_cartesian(xlim = c(150, 350)) +
  theme_classic() +
  theme(
    title = element_text(size = 100, color = "black"),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    axis.title.x = element_text(size = 100),
    axis.title.y = element_text(size = 100),
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "bar_plot_P38935_NCPR.png"), plot = bar_plot_P38935_NCPR, height = 1.6, width = 4, units = "in", dpi = 600)

#IUPred
#import data
P38935_IUPred <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\P38935_plot\\P38935_IUPred.txt",
  col_names = TRUE,
  name_repair = "universal"
)

font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

P38935_IUPred_plot <- P38935_IUPred |> 
  ggplot() +
  geom_line(aes(x = Position, y = AIUpred), color = "black", linewidth = 1) +
  geom_hline(aes(yintercept = 0.5), color = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "", y = "IUPred") +
  coord_cartesian(xlim = c(150, 350)) +
  scale_x_continuous(labels = NULL) +
  theme_classic() +
  theme(
    title = element_text(size = 100, color = "black"),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    axis.title.x = element_text(size = 80),
    axis.title.y = element_text(size = 100),
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "P38935_IUPred_plot.png"), plot = P38935_IUPred_plot, height = 1.5, width = 4, units = "in", dpi = 600)
