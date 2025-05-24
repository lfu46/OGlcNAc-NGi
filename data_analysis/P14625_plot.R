#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "showtext", "ggpubr", "rstatix", 
                    "ComplexHeatmap", "circlize", "eulerr", "ggridges")
lapply(packages_names, require, character.only = TRUE)

#protein plot
P14625_glycoprotein_data <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(UniprotID == "P14625")

P14625_glycopeptide_data <- OG_glycopeptide_Top_tb_HEK293T |> 
  filter(!is.na(combined_site)) |> 
  filter(UniprotID == "P14625")

P14625_protein_data <- WP_protein_Top_tb_HEK293T |> 
  filter(UniprotID == "P14625")

P14625_domain_data <- tribble(
  ~ domain, ~ start_position, ~ end_position, 
  'HATPase_c_3', 101, 221,
  'HSP90', 257, 773
)

P14625_rect_data <- tribble(
  ~ start, ~ end, ~ top, ~ bottom,
  0, 803, 1.7, 2.4
)

P14625_protein_text_data <- tribble(
  ~ text, ~ x_position, ~ y_position,
  "Unmod.protein", 700, 1.46,
  "Glyco.protein", 700, -1.18
)

font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

P14625_plot <- ggplot() +
  xlim(c(0, 803)) +
  ylim(c(-3, 2.5)) +
  geom_point(data = P14625_glycopeptide_data, aes(x = site_position, y = logFC), size = 2) +
  geom_hline(data = P14625_glycoprotein_data, aes(yintercept = logFC), linetype = "dashed", color = Color_2, linewidth = 0.5) +
  geom_hline(data = P14625_protein_data, aes(yintercept = logFC), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_rect(data = P14625_rect_data, aes(xmin = start, xmax = end, ymin = bottom, ymax = top), fill = "gray", color = "black", alpha = 0.5, inherit.aes = FALSE) +
  geom_rect(data = P14625_domain_data, aes(xmin = start_position, xmax = end_position, ymin = P14625_rect_data$bottom, ymax = P14625_rect_data$top, fill = domain), alpha = 0.5) +
  geom_point(data = P14625_glycopeptide_data, aes(x = site_position, y = 2.4), shape = 23, size = 2, fill = Color_2, color = "black") +
  labs(x = "Amino acid position", y = expression(log[2]*"(Tuni/Ctrl)", fill = "Domain")) +
  scale_fill_manual(
    name = "Domain",
    values = c(
      "HATPase_c_3" = Color_1,
      "HSP90" = Color_3
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    legend.key.size = unit(0.1, "in")
  )

ggsave(
  filename = paste0(file_path, "P14625_plot.eps"),
  plot = P14625_plot,
  device = cairo_ps,
  height = 1.5, width = 4, units = 'in',
  fallback_resolution = 1200
)

font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

P14625_plot_protein_only <- ggplot() +
  xlim(c(0, 803)) +
  ylim(c(-1.3, 1.6)) +
  geom_rect(data = P14625_rect_data, aes(xmin = start, xmax = end, ymin = bottom, ymax = top), fill = "gray", color = "black", alpha = 0.5, inherit.aes = FALSE) +
  geom_rect(data = P14625_domain_data, aes(xmin = start_position, xmax = end_position, ymin = P14625_rect_data$bottom, ymax = P14625_rect_data$top, fill = domain), alpha = 0.5) +
  geom_point(data = P14625_glycopeptide_data, aes(x = site_position, y = 0.12), shape = 23, size = 2, fill = Color_2, color = "black") +
  geom_point(data = P14625_ubi_data, aes(x = ubi_site, y = 0.06), shape = 21, size = 2, fill = "orange", color = "black") +
  geom_text(data = P14625_glycopeptide_data, aes(x = site_position, y = 0.25, label = combined_site), size = 30, angle = 90) +
  scale_fill_manual(
    name = "Domain",
    values = c(
      "HATPase_c_3" = Color_1,
      "HSP90" = Color_3
    )
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = paste0(file_path, "P14625_plot_protein_only.eps"),
  plot = P14625_plot_protein_only,
  device = cairo_ps,
  height = 2, width = 4, units = 'in',
  fallback_resolution = 1200
)

#IUPred
#import data
P14625_AIUPred <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\P14625_plot\\AIUPred_P14625.txt",
  col_names = TRUE,
  name_repair = "universal"
)

font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

P14625_AIUPred_plot <- P14625_AIUPred |> 
  ggplot() +
  geom_line(aes(x = Position, y = AIUpred), color = "black", linewidth = 0.8) +
  geom_hline(aes(yintercept = 0.5), color = "red", linetype = "dashed", linewidth = 0.5) +
  labs(x = "", y = "IUPred") +
  coord_cartesian(xlim = c(0, 803)) +
  scale_x_continuous(labels = NULL) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 10),
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    legend.position = "none"
  )


ggsave(
  filename = paste0(file_path, "P14625_AIUPred_plot.eps"),
  plot = P14625_AIUPred_plot,
  device = 'eps',
  height = 1, width = 2.9, units = 'in'
)

