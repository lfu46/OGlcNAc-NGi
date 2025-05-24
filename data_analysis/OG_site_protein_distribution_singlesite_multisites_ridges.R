#import packages
packages_names <- c("showtext", "ggpubr", "rstatix", "ggridges")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
#HEK293T
OG_multiplesites_glycoprotein_HEK293T_list <- OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(!is.na(combined_site)) |> 
  count(UniprotID) |> 
  filter(n > 1) |> 
  pull(UniprotID)

OG_multiplesites_glycoprotein_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T$UniprotID) |> 
  mutate(multiplesites = ifelse(UniprotID %in% OG_multiplesites_glycoprotein_HEK293T_list, "multiple sites", "single site")) |> 
  select(UniprotID, multiplesites, logFC) |> 
  mutate(Cell = "HEK293T")

#HepG2
OG_multiplesites_glycoprotein_HepG2_list <- OG_site_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(!is.na(combined_site)) |> 
  count(UniprotID) |> 
  filter(n > 1) |> 
  pull(UniprotID)

OG_multiplesites_glycoprotein_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(UniprotID %in% OG_site_protein_OGlcNAc_effect_Top_tb_HepG2$UniprotID) |> 
  mutate(multiplesites = ifelse(UniprotID %in% OG_multiplesites_glycoprotein_HepG2_list, "multiple sites", "single site")) |> 
  select(UniprotID, multiplesites, logFC) |> 
  mutate(Cell = "HepG2")

#Jurkat
OG_multiplesites_glycoprotein_Jurkat_list <- OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(!is.na(combined_site)) |> 
  count(UniprotID) |> 
  filter(n > 1) |> 
  pull(UniprotID)

OG_multiplesites_glycoprotein_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(UniprotID %in% OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat$UniprotID) |> 
  mutate(multiplesites = ifelse(UniprotID %in% OG_multiplesites_glycoprotein_Jurkat_list, "multiple sites", "single site")) |> 
  select(UniprotID, multiplesites, logFC) |> 
  mutate(Cell = "Jurkat")

#combine
OG_multiplesites_glycoprotein_combined <- bind_rows(
  OG_multiplesites_glycoprotein_HepG2,
  OG_multiplesites_glycoprotein_HEK293T,
  OG_multiplesites_glycoprotein_Jurkat
)

#wilcox test
OG_multisites_glycoprotein_wilcox_test <- OG_multiplesites_glycoprotein_combined |> 
  group_by(Cell) |> 
  wilcox_test(logFC ~ multiplesites, p.adjust.method = "BH") |> 
  add_significance("p")

#density ridges plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

density_ridges_plot_OG_multiplesites_glycoprotein <- OG_multiplesites_glycoprotein_combined |>
  ggplot() +
  geom_density_ridges(
    aes(x = logFC, y = multiplesites, fill = Cell), color = "black",
    quantile_lines = TRUE, quantiles = 2, 
    scale = 0.55,
    jittered_points = TRUE, position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 8, point_alpha = 1
  ) +
  scale_fill_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  scale_color_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  labs(x = expression(Δlog[2]*"FC"), y = "") +
  geom_text(data = OG_multisites_glycoprotein_wilcox_test, 
            aes(y = group2, x = 2.7, label = p.signif), size = 50, vjust = -0.1) +
  facet_grid(~ Cell, scales = "free_x") +
  coord_cartesian(xlim = c(-3, 3.5)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black"),
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent")
  )

ggsave(filename = paste0(file_path, "density_ridges_plot_OG_multiplesites_glycoprotein.png"), 
       plot = density_ridges_plot_OG_multiplesites_glycoprotein, height = 3, width = 8, units = "in", dpi = 600)
