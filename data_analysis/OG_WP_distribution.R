#import packages
packages_names <- c("tidyverse", "rstatix", "ggpubr", "showtext", "introdataviz")
lapply(packages_names, require, character.only = TRUE)

#generate data frame for HepG2
OG_WP_distribution_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  select(UniprotID, logFC_OG = logFC) |> 
  left_join(WP_protein_Top_tb_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC_OG, logFC_WP = logFC) |> 
  filter(!is.na(logFC_WP)) |> 
  pivot_longer(cols = logFC_OG:logFC_WP, names_to = "Exp", values_to = "logFC") |> 
  mutate(Cell = "HepG2")

#ks test
ks.test(
  OG_WP_distribution_HepG2 |> filter(Exp == "logFC_OG") |> pull(logFC),
  OG_WP_distribution_HepG2 |> filter(Exp == "logFC_WP") |> pull(logFC)
)

#wilcox test
OG_WP_distribution_HepG2_wilcox_test <- OG_WP_distribution_HepG2 |> 
  wilcox_test(logFC ~ Exp, p.adjust.method = "BH") |> 
  add_significance("p")

#generate data frame for HEK293T
OG_WP_distribution_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  select(UniprotID, logFC_OG = logFC) |> 
  left_join(WP_protein_Top_tb_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC_OG, logFC_WP = logFC) |> 
  filter(!is.na(logFC_WP)) |> 
  pivot_longer(cols = logFC_OG:logFC_WP, names_to = "Exp", values_to = "logFC") |> 
  mutate(Cell = "HEK293T")

#ks test
ks.test(
  OG_WP_distribution_HEK293T |> filter(Exp == "logFC_OG") |> pull(logFC),
  OG_WP_distribution_HEK293T |> filter(Exp == "logFC_WP") |> pull(logFC)
)

#wilcox test
OG_WP_distribution_HEK293T_wilcox_test <- OG_WP_distribution_HEK293T |> 
  wilcox_test(logFC ~ Exp, p.adjust.method = "BH") |> 
  add_significance("p")

#generate data frame for Jurkat
OG_WP_distribution_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  select(UniprotID, logFC_OG = logFC) |> 
  left_join(WP_protein_Top_tb_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC_OG, logFC_WP = logFC) |> 
  filter(!is.na(logFC_WP)) |> 
  pivot_longer(cols = logFC_OG:logFC_WP, names_to = "Exp", values_to = "logFC") |> 
  mutate(Cell = "Jurkat")

#ks test
ks.test(
  OG_WP_distribution_Jurkat |> filter(Exp == "logFC_OG") |> pull(logFC),
  OG_WP_distribution_Jurkat |> filter(Exp == "logFC_WP") |> pull(logFC)
)

#wilcox test
OG_WP_distribution_Jurkat_wilcox_test <- OG_WP_distribution_Jurkat |> 
  wilcox_test(logFC ~ Exp, p.adjust.method = "BH") |> 
  add_significance("p")

#combine
OG_WP_distribution_combined <- bind_rows(
  OG_WP_distribution_HepG2,
  OG_WP_distribution_HEK293T,
  OG_WP_distribution_Jurkat
)

#combine ks test

#wilcox test
OG_WP_distribution_combined_wilcox_test <- OG_WP_distribution_combined |> 
  group_by(Cell) |> 
  wilcox_test(logFC ~ Exp, p.adjust.method = "BH") |> 
  add_significance("p")

#split violin plot
ks_test_label <- tibble(
  Cell = c("HEK293T", "HEK293T", "HepG2", "HepG2", "Jurkat", "Jurkat"),
  p_signif = c("****", "****", "****", "****", "****", "****"),
  logFC = 1.5,
  Name = c("glycoprotein_HEK293T", "glycoprotein_HEK293T", "glycoprotein_HepG2", "glycoprotein_HepG2", "glycoprotein_Jurkat", "glycoprotein_Jurkat"),
  Exp = c("OG_HEK293T", "WP", "OG_HepG2", "WP", "OG_Jurkat", "WP")
)

font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

split_violin_plot_OG_WP_distribution <- OG_WP_distribution_combined |> 
  mutate(Name = paste("glycoprotein", Cell, sep = "_"), Exp = ifelse(Exp == "logFC_OG", paste("OG", Cell, sep = "_"), "WP")) |> 
  ggplot(aes(x = Name, y = logFC, fill = factor(Exp, c("OG_HEK293T", "OG_HepG2", "OG_Jurkat", "WP")))) +
  geom_split_violin(color = "transparent") +
  geom_text(data = ks_test_label, aes(y = logFC, label = p_signif), size = 40, hjust = -0.2) +
  scale_fill_manual(values = c(
    "OG_HEK293T" = Color_2,
    "OG_HepG2" = Color_3,
    "OG_Jurkat" = Color_4,
    "WP" = "gray"
  )) +
  facet_grid(~ Cell, scale = "free_x") +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)"), fill = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.line.x = element_blank(),
    axis.ticks.length.x = unit(0, "in"),
    legend.position = "none",
    strip.text = element_text(size = 15, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(
  filename = paste0(file_path, "split_violin_plot_OG_WP_distribution.eps"),
  device = "eps",
  plot = split_violin_plot_OG_WP_distribution, 
  height = 3, width = 4, units = c("in"), 
  dpi = 1200
  )
