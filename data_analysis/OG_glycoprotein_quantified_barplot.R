#import packages
packages_names <- c("showtext", "ggpubr", "rstatix")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
OG_glycoprotein_quantified_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  select(UniprotID) |> 
  mutate(
    HepG2 = 1,
    HEK293T = ifelse(UniprotID %in% OG_glycoprotein_Top_tb_HEK293T$UniprotID, 1, 0),
    Jurkat = ifelse(UniprotID %in% OG_glycoprotein_Top_tb_Jurkat$UniprotID, 1, 0)
  ) |> 
  mutate(
    count = HepG2 + HEK293T + Jurkat,
    Cell = "HepG2"
  )

OG_glycoprotein_quantified_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  select(UniprotID) |> 
  mutate(
    HEK293T = 1,
    HepG2 = ifelse(UniprotID %in% OG_glycoprotein_Top_tb_HepG2$UniprotID, 1, 0),
    Jurkat = ifelse(UniprotID %in% OG_glycoprotein_Top_tb_Jurkat$UniprotID, 1, 0)
  ) |> 
  mutate(
    count = HepG2 + HEK293T + Jurkat,
    Cell = "HEK293T"
  )

OG_glycoprotein_quantified_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  select(UniprotID) |> 
  mutate(
    Jurkat = 1,
    HepG2 = ifelse(UniprotID %in% OG_glycoprotein_Top_tb_HepG2$UniprotID, 1, 0),
    HEK293T = ifelse(UniprotID %in% OG_glycoprotein_Top_tb_HEK293T$UniprotID, 1, 0)
  ) |> 
  mutate(
    count = HepG2 + HEK293T + Jurkat,
    Cell = "Jurkat"
  )

OG_glycoprotein_quantified_total <- bind_rows(
  OG_glycoprotein_quantified_HepG2,
  OG_glycoprotein_quantified_HEK293T,
  OG_glycoprotein_quantified_Jurkat
) |> 
  mutate(
    Cell = factor(Cell, levels = c("HEK293T", "HepG2", "Jurkat")),
    identified = factor(count, levels = c(1, 2, 3))
  ) |> 
  count(Cell, identified)

#barplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

barplot_OG_glycoprotein_quantified_total <- OG_glycoprotein_quantified_total |> 
  ggplot(aes(y = Cell, x = n, fill = Cell, alpha = identified)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), 
            position = position_stack(vjust = 0.5),
            color = "black",
            size = 30
          ) +
  labs(y = "", x = "# of quantified glycoprotein") +
  scale_fill_manual(
    values = c(
    "HEK293T" = Color_2,
    "HepG2" = Color_3,
    "Jurkat" = Color_4
  ), guide = "none") +
  scale_alpha_manual(
    name = "Quantified in # cell line",
    values = c(
      "1" = 0.6,
      "2" = 0.8,
      "3" = 1.0
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black"),
    legend.position = "top"
  )

ggsave(filename = paste0(file_path, "barplot_OG_glycoprotein_quantified_total.png"), plot = barplot_OG_glycoprotein_quantified_total, height = 4, width = 5, units = "in", dpi = 600)



