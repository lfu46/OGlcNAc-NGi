#import packages
packages_names <- c("tidyverse", "showtext", "writexl", "ggpubr", "rstatix")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
OG_site_protein_logFC_HEK293T <- OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  select(Index, logFC) |> 
  mutate(cell = "HEK293T")

OG_site_protein_logFC_HepG2 <- OG_site_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  select(Index, logFC) |> 
  mutate(cell = "HepG2")

OG_site_protein_logFC_Jurkat <- OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  select(Index, logFC) |> 
  mutate(cell = "Jurkat")

#combine
OG_site_protein_logFC_combined <- bind_rows(
  OG_site_protein_logFC_HEK293T,
  OG_site_protein_logFC_HepG2,
  OG_site_protein_logFC_Jurkat
)

#wilcox test
OG_site_protein_wilcox_test <- OG_site_protein_logFC_combined |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#ks test
ks.test(OG_site_protein_logFC_HEK293T$logFC, OG_site_protein_logFC_HepG2$logFC)
ks.test(OG_site_protein_logFC_HepG2$logFC, OG_site_protein_logFC_Jurkat$logFC)
ks.test(OG_site_protein_logFC_Jurkat$logFC, OG_site_protein_logFC_HEK293T$logFC)

OG_site_protein_ks_test <- OG_site_protein_wilcox_test |> 
  mutate(
    p = c(0.0005235, 8.328e-09, 2.517e-05),
    p.signif = c('***', '****', '****')
  )

#boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_OG_site_protein_logFC <- OG_site_protein_logFC_combined |> 
  ggplot(aes(x = cell, y = logFC)) +
  geom_boxplot(aes(color = cell), fill = "transparent", outliers = FALSE, notch = TRUE, size = 1.5) +
  scale_color_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  labs(x = "", y = expression(Δlog[2]*"FC")) +
  stat_pvalue_manual(data = OG_site_protein_ks_test, label = "p.signif", tip.length = 0, 
                     hide.ns = "p", size = 40, y.position = c(1.8, 2.1, 2.4), coord.flip = TRUE) +
  coord_flip(ylim = c(-2.5, 2.5)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black"),
    legend.position = 'none'
  )

ggsave(filename = paste0(file_path, "boxplot_OG_site_protein_logFC.png"), 
       plot = boxplot_OG_site_protein_logFC, height = 2, width = 4, units = "in", dpi = 600)
