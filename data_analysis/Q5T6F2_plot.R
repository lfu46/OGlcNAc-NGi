#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "showtext", "ggpubr", "rstatix", 
                    "ComplexHeatmap", "circlize", "eulerr", "ggridges")
lapply(packages_names, require, character.only = TRUE)

#protein plot
Q5T6F2_glycoprotein_data <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(UniprotID == "Q5T6F2")

Q5T6F2_glycopeptide_data <- OG_glycopeptide_Top_tb_HEK293T |> 
  filter(!is.na(combined_site)) |> 
  filter(UniprotID == "Q5T6F2") |> 
  mutate(site_position_adj = c(372, 970, 1020, 1119, 1073, 427, 517, 637, 469, 200))

Q5T6F2_protein_data <- WP_protein_Top_tb_HEK293T |> 
  filter(UniprotID == "Q5T6F2")

Q5T6F2_domain_data <- tribble(
  ~ domain, ~ start_position, ~ end_position, 
  'UBA', 53, 91,
  'UBAP2/protein lingerer', 512, 544
)

Q5T6F2_ubi_data <- tribble(
  ~ ubi_site,
  16,
  48,
  84,
  242,
  248,
  868,
  981,
  993
)

Q5T6F2_rect_data <- tribble(
  ~ start, ~ end, ~ top, ~ bottom,
  0, 1119, 0.03, 0.09
)

Q5T6F2_protein_text_data <- tribble(
  ~ text, ~ x_position, ~ y_position,
  "Unmod.protein", 950, -0.22,
  "Glyco.protein", 950, -0.73
)

font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

Q5T6F2_plot <- ggplot() +
  xlim(c(0, 1119)) +
  ylim(c(-1.1, 0.32)) +
  geom_point(data = Q5T6F2_glycopeptide_data, aes(x = site_position, y = logFC), size = 2) +
  geom_hline(data = Q5T6F2_glycoprotein_data, aes(yintercept = logFC), linetype = "dashed", color = Color_2, linewidth = 1) +
  geom_hline(data = Q5T6F2_protein_data, aes(yintercept = logFC), linetype = "dashed", color = "black", linewidth = 1) +
  geom_rect(data = Q5T6F2_rect_data, aes(xmin = start, xmax = end, ymin = bottom, ymax = top), fill = "gray", color = "black", alpha = 0.5, inherit.aes = FALSE) +
  geom_rect(data = Q5T6F2_domain_data, aes(xmin = start_position, xmax = end_position, ymin = Q5T6F2_rect_data$bottom, ymax = Q5T6F2_rect_data$top, fill = domain), alpha = 0.5) +
  geom_point(data = Q5T6F2_glycopeptide_data, aes(x = site_position, y = 0.12), shape = 23, size = 2, fill = Color_2, color = "black") +
  geom_point(data = Q5T6F2_ubi_data, aes(x = ubi_site, y = 0.06), shape = 21, size = 2, fill = "orange", color = "black") +
  geom_text(data = Q5T6F2_glycopeptide_data, aes(x = site_position_adj, y = 0.25, label = combined_site), size = 30, angle = 90) +
  geom_text(data = Q5T6F2_protein_text_data, aes(x = x_position, y = y_position, label = text, color = text), size = 30, show.legend = FALSE) +
  scale_color_manual(values = c(
    "Unmod.protein" = "black",
    "Glyco.protein" = Color_2
  )) +
  labs(x = "Amino acid position", y = expression(log[2]*"(Tuni/Ctrl)", fill = "Domain")) +
  scale_fill_manual(
    name = "Domain",
    values = c(
      "UBA" = Color_1,
      "UBAP2/protein lingerer" = Color_3
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

ggsave(paste0(file_path, "Q5T6F2_plot.png"), plot = Q5T6F2_plot, height = 4, width = 7, dpi = 600)

font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

Q5T6F2_plot_protein_only <- ggplot() +
  xlim(c(0, 1119)) +
  ylim(c(-1.1, 0.32)) +
  geom_rect(data = Q5T6F2_rect_data, aes(xmin = start, xmax = end, ymin = bottom, ymax = top), fill = "gray", color = "black", alpha = 0.5, inherit.aes = FALSE) +
  geom_rect(data = Q5T6F2_domain_data, aes(xmin = start_position, xmax = end_position, ymin = Q5T6F2_rect_data$bottom, ymax = Q5T6F2_rect_data$top, fill = domain), alpha = 0.5) +
  geom_point(data = Q5T6F2_glycopeptide_data, aes(x = site_position, y = 0.12), shape = 23, size = 2, fill = Color_2, color = "black") +
  geom_point(data = Q5T6F2_ubi_data, aes(x = ubi_site, y = 0.06), shape = 21, size = 2, fill = "orange", color = "black") +
  geom_text(data = Q5T6F2_glycopeptide_data, aes(x = site_position_adj, y = 0.25, label = combined_site), size = 30, angle = 90) +
  scale_fill_manual(
    name = "Domain",
    values = c(
      "UBA" = Color_1,
      "UBAP2/protein lingerer" = Color_3
    )
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  )

ggsave(paste0(file_path, "Q5T6F2_plot_protein_only.png"), plot = Q5T6F2_plot_protein_only, height = 3.5, width = 3.8, dpi = 600)

#NCPR
#import data
Q5T6F2_NCPR_data <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\Q5T6F2_plot\\NCPR.xlsx",
  col_names = c('position', 'NCPR')
)

#bar plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

bar_plot_Q5T6F2_NCPR <- Q5T6F2_NCPR_data |> 
  mutate(group = ifelse(NCPR > 0, "high", ifelse(NCPR < 0, "low", "none"))) |> 
  ggplot() +
  geom_bar(aes(x = position, y = NCPR, fill = group), stat = "identity", color = "transparent", width = 3) +
  scale_fill_manual(
    values = c(
      "high" = Color_3,
      "low" = Color_4
    )
  ) +
  labs(x = "Amino acid position") +
  coord_cartesian(xlim = c(0, 1119)) +
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

ggsave(filename = paste0(file_path, "bar_plot_Q5T6F2_NCPR.png"), plot = bar_plot_Q5T6F2_NCPR, height = 1.6, width = 4, units = "in", dpi = 600)

#IUPred
#import data
Q5T6F2_AIUPred <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\Q5T6F2_plot\\Q5T6F2_AIUPred.txt",
  col_names = TRUE,
  name_repair = "universal"
)

font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

Q5T6F2_AIUPred_plot <- Q5T6F2_AIUPred |> 
  ggplot() +
  geom_line(aes(x = Position, y = AIUpred), color = "black", linewidth = 0.8) +
  geom_hline(aes(yintercept = 0.5), color = "red", linetype = "dashed", linewidth = 0.8) +
  labs(x = "", y = "IUPred") +
  coord_cartesian(xlim = c(0, 1119)) +
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

ggsave(filename = paste0(file_path, "Q5T6F2_AIUPred_plot.png"), plot = Q5T6F2_AIUPred_plot, height = 1.5, width = 4, units = "in", dpi = 600)
