#import packages
packages_names <- c("showtext", "ggpubr", "rstatix")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
OG_glycoprotein_logFC_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  mutate(cell = "HepG2") |> 
  select(UniprotID, logFC, cell)

OG_glycoprotein_logFC_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  mutate(cell = "HEK293T") |> 
  select(UniprotID, logFC, cell)

OG_glycoprotein_logFC_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  mutate(cell = "Jurkat") |> 
  select(UniprotID, logFC, cell)

OG_glycoprotein_logFC_combined <- bind_rows(OG_glycoprotein_logFC_HepG2, OG_glycoprotein_logFC_HEK293T, OG_glycoprotein_logFC_Jurkat)

#ks.test
ks.test(OG_glycoprotein_Top_tb_HepG2$logFC, OG_glycoprotein_Top_tb_HEK293T$logFC)
ks.test(OG_glycoprotein_Top_tb_HEK293T$logFC, OG_glycoprotein_Top_tb_Jurkat$logFC)
ks.test(OG_glycoprotein_Top_tb_Jurkat$logFC, OG_glycoprotein_Top_tb_HepG2$logFC)

OG_glycoprotein_ks_test <- tribble(
  ~ .y., ~ group1, ~ group2, ~ p, ~ p.signif,
  "logFC", "HEK293T", "HepG2", 2.153e-13, "****",
  "logFC", "Jurkat", "HEK293T", 0.04457, "*",
  "logFC", "HepG2", "Jurkat", 6.29e-15, "****"
)

#wilcox test
OG_glycoprotein_logFC_wilcox_test <- OG_glycoprotein_logFC_combined |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_logFC <- OG_glycoprotein_logFC_combined |>
  ggplot() +
  geom_violin(aes(x = cell, y = logFC, fill = cell), color = "transparent") +
  geom_boxplot(aes(x = cell, y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  scale_fill_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  stat_pvalue_manual(data = OG_glycoprotein_ks_test, label = "p.signif", tip.length = 0, 
                     size = 7, y.position = c(1.5, 1.9, 1.7)) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 12, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none"
  )

ggsave(
  filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_logFC.eps"),
  device = "eps",
  plot = violin_boxplot_OG_glycoprotein_logFC, 
  height = 2.5, width = 2, units = "in", 
  dpi = 1200
  )
