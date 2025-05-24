#import packages
packages_names <- c("tidyverse", "showtext", "writexl", "ggpubr", "rstatix")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
OG_glycoprotein_protein_logFC_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  select(UniprotID, logFC) |> 
  mutate(cell = "HEK293T")

OG_glycoprotein_protein_logFC_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  select(UniprotID, logFC) |> 
  mutate(cell = "HepG2")

OG_glycoprotein_protein_logFC_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  select(UniprotID, logFC) |> 
  mutate(cell = "Jurkat")

#combine
OG_glycoprotein_protein_logFC_combined <- bind_rows(
  OG_glycoprotein_protein_logFC_HEK293T,
  OG_glycoprotein_protein_logFC_HepG2,
  OG_glycoprotein_protein_logFC_Jurkat
)

#wilcox test
OG_glycoprotein_protein_wilcox_test <- OG_glycoprotein_protein_logFC_combined |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#ks test
ks.test(OG_glycoprotein_protein_logFC_HEK293T$logFC, OG_glycoprotein_protein_logFC_HepG2$logFC)
ks.test(OG_glycoprotein_protein_logFC_HepG2$logFC, OG_glycoprotein_protein_logFC_Jurkat$logFC)
ks.test(OG_glycoprotein_protein_logFC_Jurkat$logFC, OG_glycoprotein_protein_logFC_HEK293T$logFC)

OG_glycoprotein_protein_ks_test <- OG_glycoprotein_protein_wilcox_test |> 
  mutate(
    p = c(1.647e-15, 0.0002669, 7.718e-15),
    p.signif = c('****', '***', '****')
  )

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_protein_logFC <- OG_glycoprotein_protein_logFC_combined |> 
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
  labs(x = "", y = expression(Δlog[2]*"FC")) +
  stat_pvalue_manual(data = OG_glycoprotein_protein_ks_test, label = "p.signif", tip.length = 0, 
                     hide.ns = "p", size = 50, y.position = c(1.8, 2.1, 2.4)) +
  coord_cartesian(ylim = c(-2.5, 2.5)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_protein_logFC.png"), 
       plot = violin_boxplot_OG_glycoprotein_protein_logFC, height = 4, width = 4, units = "in", dpi = 600)
