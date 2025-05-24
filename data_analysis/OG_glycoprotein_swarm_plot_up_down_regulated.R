#import packages
packages_names <- c("tidyverse", "rstatix", "showtext", "ggpubr", "readxl", "ggbeeswarm")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
OG_glycoprotein_up_down_HepG2 <- bind_rows(
  OG_glycoprotein_Top_tb_HepG2 |> filter(UniprotID %in% OG_glycoprotein_up_total$UniprotID) |> select(UniprotID, logFC) |> mutate(Category = "up"),
  OG_glycoprotein_Top_tb_HepG2 |> filter(UniprotID %in% OG_glycoprotein_down_total$UniprotID) |> select(UniprotID, logFC) |> mutate(Category = "down")
) |> mutate(Cell = "HepG2")

OG_glycoprotein_up_down_HEK293T <- bind_rows(
  OG_glycoprotein_Top_tb_HEK293T |> filter(UniprotID %in% OG_glycoprotein_up_total$UniprotID) |> select(UniprotID, logFC) |> mutate(Category = "up"),
  OG_glycoprotein_Top_tb_HEK293T |> filter(UniprotID %in% OG_glycoprotein_down_total$UniprotID) |> select(UniprotID, logFC) |> mutate(Category = "down")
) |> mutate(Cell = "HEK293T")

OG_glycoprotein_up_down_Jurkat <- bind_rows(
  OG_glycoprotein_Top_tb_Jurkat |> filter(UniprotID %in% OG_glycoprotein_up_total$UniprotID) |> select(UniprotID, logFC) |> mutate(Category = "up"),
  OG_glycoprotein_Top_tb_Jurkat |> filter(UniprotID %in% OG_glycoprotein_down_total$UniprotID) |> select(UniprotID, logFC) |> mutate(Category = "down")
) |> mutate(Cell = "Jurkat")

OG_glycoprotein_up_down_total_logFC <- bind_rows(
  OG_glycoprotein_up_down_HepG2,
  OG_glycoprotein_up_down_HEK293T,
  OG_glycoprotein_up_down_Jurkat
)

#t test
OG_glycoprotein_up_down_total_logFC_t_test <- OG_glycoprotein_up_down_total_logFC |> 
  group_by(Cell) |> 
  t_test(logFC ~ Category, p.adjust.method = "BH") |> 
  add_significance("p")

#swarm plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

swarm_plot_OG_glycoprotein_up_down_total_logFC <- OG_glycoprotein_up_down_total_logFC |> 
  ggplot(aes(x = factor(Category, levels = c("up", "down")), y = logFC, color = factor(Category, levels = c("up", "down")))) +
  geom_beeswarm() +
  scale_color_manual(
    name = "",
    values = c(
    "up" = Color_5,
    "down" = Color_6
  )) +
  stat_pvalue_manual(data = OG_glycoprotein_up_down_total_logFC_t_test, label = "p.signif", tip.length = 0, size = 30, hide.ns = "p.signif",
                     y.position = c(2.8)) +
  facet_grid(~ Cell, scale = "free_x") +
  coord_cartesian(ylim = c(-3, 3)) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)"), fill = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 80, color = "black"),
    axis.line.x = element_blank(),
    axis.ticks.length.x = unit(0, "in"),
    strip.text = element_text(size = 80, color = "black", margin = margin(b = 0.01, t = 0, unit = "line")),
    strip.background = element_rect(fill = NA, color = NA),
    legend.position = "bottom",
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "swarm_plot_OG_glycoprotein_up_down_total_logFC.png"), plot = swarm_plot_OG_glycoprotein_up_down_total_logFC, height = 4, width = 6, units = c("in"), dpi = 600)

