#import packages
packages_names <- c("tidyverse", "rstatix", "showtext", "ggpubr", "readxl")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
OG_glycoprotein_higher_zero_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC > 0) |> select(UniprotID)
OG_glycoprotein_higher_zero_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC > 0) |> select(UniprotID)
OG_glycoprotein_higher_zero_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC > 0) |> select(UniprotID)

OG_glycoprotein_lower_zero_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC < 0) |> select(UniprotID)
OG_glycoprotein_lower_zero_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC < 0) |> select(UniprotID)
OG_glycoprotein_lower_zero_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC < 0) |> select(UniprotID)

#generate data frame
#up
OG_glycoprotein_up_combined <- bind_rows(
  OG_glycoprotein_up_HepG2,
  OG_glycoprotein_up_HEK293T,
  OG_glycoprotein_up_Jurkat
) |> distinct() |> 
  filter(! UniprotID %in% OG_glycoprotein_lower_zero_HEK293T$UniprotID) |> 
  filter(! UniprotID%in% OG_glycoprotein_lower_zero_HepG2$UniprotID) |> 
  filter(! UniprotID %in% OG_glycoprotein_lower_zero_Jurkat$UniprotID)

write_xlsx(OG_glycoprotein_up_combined, path = paste0(file_path, "OG_glycoprotein_up_combined.xlsx"))

#down
OG_glycoprotein_down_combined <- bind_rows(
  OG_glycoprotein_down_HepG2,
  OG_glycoprotein_down_HEK293T,
  OG_glycoprotein_down_Jurkat
) |> distinct() |> 
  filter(! UniprotID %in% OG_glycoprotein_higher_zero_HEK293T$UniprotID) |> 
  filter(! UniprotID %in% OG_glycoprotein_higher_zero_HepG2$UniprotID) |> 
  filter(! UniprotID %in% OG_glycoprotein_higher_zero_Jurkat$UniprotID)

write_xlsx(OG_glycoprotein_down_combined, path = paste0(file_path, "OG_glycoprotein_down_combined.xlsx"))

#generate data frame
OG_glycoprotein_up_logFC_combined <- OG_glycoprotein_up_combined |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC_HepG2 = logFC) |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC_HepG2, logFC_HEK293T = logFC) |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC_HEK293T, logFC_HepG2, logFC_Jurkat = logFC) |> 
  arrange(desc(logFC_HEK293T)) |> 
  mutate(category = "up")

OG_glycoprotein_down_logFC_combined <- OG_glycoprotein_down_combined |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC_HepG2 = logFC) |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC_HepG2, logFC_HEK293T = logFC) |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC_HEK293T, logFC_HepG2, logFC_Jurkat = logFC) |> 
  arrange(desc(logFC_HEK293T)) |> 
  mutate(category = "down")

OG_glycoprotein_up_down_logFC_combined <- bind_rows(
  OG_glycoprotein_up_logFC_combined, 
  OG_glycoprotein_down_logFC_combined
) |> 
  pivot_longer(cols = starts_with("logFC"), names_to = 'Cell', values_to = 'logFC') |> 
  mutate(Cell = ifelse(Cell == "logFC_HEK293T", "HEK293T", ifelse(Cell == "logFC_HepG2", "HepG2", "Jurkat")))

#heatmap
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

heatmap_OG_glycoprotein_up_down_logFC <- OG_glycoprotein_up_down_logFC_combined |> 
  ggplot(aes(y = Cell, x = factor(UniprotID, levels = OG_glycoprotein_up_down_logFC_combined |> distinct(UniprotID) |> pull()), 
             fill = logFC)) +
  geom_tile() +
  labs(fill = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_binned(type = "viridis", na.value = "transparent", limits = c(-4, 4), breaks = c(-4, -0.5, 0, 0.5, 4)) +
  facet_grid(~ factor(category, levels = c("up", "down")), scale = "free_x") +
  theme_void() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0),
    panel.grid.minor = element_line(color = "gray", linewidth = 0),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 80, color = "black", lineheight = 0.1),
    axis.ticks.x = element_line(linewidth = 0),
    legend.text = element_text(size = 80),
    legend.title = element_text(size = 80, color = "black"),
    legend.frame = element_rect(color = "black", linewidth = 0.2),
    legend.ticks = element_line(color = "black", linewidth = 0.2),
    strip.text = element_text(color = "black", size = 100),
    strip.background = element_rect(color = "transparent", fill = "transparent")
  )

ggsave(filename = paste0(file_path, "hheatmap_OG_glycoprotein_up_down_logFC.png"), plot = heatmap_OG_glycoprotein_up_down_logFC, height = 3, width = 5, units = c("in"), dpi = 600)

