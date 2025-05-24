#import packages
packages_names <- c("showtext", "ggpubr", "rstatix", "ggridges", "ggbeeswarm")
lapply(packages_names, require, character.only = TRUE)

#glycopeptide level
#generate data frame
OG_glycopeptide_logFC_HepG2 <- OG_glycopeptide_Top_tb_HepG2 |> 
  mutate(Cell = "HepG2") |> 
  select(UniprotID, logFC, Cell)

OG_glycopeptide_logFC_HEK293T <- OG_glycopeptide_Top_tb_HEK293T |> 
  mutate(Cell = "HEK293T") |> 
  select(UniprotID, logFC, Cell)

OG_glycopeptide_logFC_Jurkat <- OG_glycopeptide_Top_tb_Jurkat |> 
  mutate(Cell = "Jurkat") |> 
  select(UniprotID, logFC, Cell)

OG_glycopeptide_logFC_total <- bind_rows(
  OG_glycopeptide_logFC_HepG2,
  OG_glycopeptide_logFC_HEK293T,
  OG_glycopeptide_logFC_Jurkat
)

#glycopeptide level
#wilcox test
OG_glycopeptide_logFC_wilcox_test <- OG_glycopeptide_logFC_total |> 
  wilcox_test(logFC ~ Cell, p.adjust.method = "BH") |> 
  add_significance("p")

#glycopeptide level
#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycopeptide_logFC_total <- OG_glycopeptide_logFC_total |>
  ggplot() +
  geom_violin(aes(x = Cell, y = logFC, fill = Cell)) +
  geom_boxplot(aes(x = Cell, y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  scale_fill_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  stat_pvalue_manual(data = OG_glycopeptide_logFC_wilcox_test, label = "p.signif", tip.length = 0, size = 30,
                     hide.ns = "p", y.position = c(3.0, 3.2)) +
  coord_cartesian(ylim = c(-3, 3.5)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycopeptide_logFC_total.png"), plot = violin_boxplot_OG_glycopeptide_logFC_total, height = 4, width = 5, units = "in", dpi = 600)

#site level
#generate data frame
OG_site_logFC_HepG2 <- OG_glycopeptide_Top_tb_HepG2 |> 
  filter(!is.na(combined_site)) |> 
  mutate(Cell = "HepG2") |> 
  select(UniprotID, logFC, Cell)

OG_site_logFC_HEK293T <- OG_glycopeptide_Top_tb_HEK293T |> 
  filter(!is.na(combined_site)) |> 
  mutate(Cell = "HEK293T") |> 
  select(UniprotID, logFC, Cell)

OG_site_logFC_Jurkat <- OG_glycopeptide_Top_tb_Jurkat |> 
  filter(!is.na(combined_site)) |> 
  mutate(Cell = "Jurkat") |> 
  select(UniprotID, logFC, Cell)

OG_site_logFC_total <- bind_rows(
  OG_site_logFC_HepG2,
  OG_site_logFC_HEK293T,
  OG_site_logFC_Jurkat
)

#site level
#wilcox test
OG_site_logFC_wilcox_test <- OG_site_logFC_total |> 
  wilcox_test(logFC ~ Cell, p.adjust.method = "BH") |> 
  add_significance("p")

OG_site_logFC_wilcox_test_adj <- OG_site_logFC_wilcox_test |> 
  mutate(group1 = ifelse(p.signif == "**", "Jurkat", group1),
         group2 = ifelse(p.signif == "**", "HepG2", group2))

#site level
#ridges plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_OG_site_logFC <- OG_site_logFC_total |>
  ggplot(aes(x = Cell, y = logFC)) +
  geom_boxplot(aes(color = Cell), fill = "transparent", outliers = FALSE,
               notch = TRUE, size = 1.5) +
  scale_color_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  stat_pvalue_manual(data = OG_site_logFC_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 50, 
                     hide.ns = "p", y.position = c(1.8), coord.flip = TRUE) +
  coord_flip(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "boxplot_OG_site_logFC.png"), plot = boxplot_OG_site_logFC, 
       height = 2, width = 4, units = "in", dpi = 600)
