#import packages
packages_names <- c("UpSetR", "ggupset", "ggalluvial", "showtext")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
#HEK293T
OG_glycoprotein_up_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> select(UniprotID) |> 
  mutate(category = "up" , group = "Raw", count = 1)

OG_glycoprotein_down_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> select(UniprotID) |> 
  mutate(category = "down", group = "Raw", count = 1)

OG_glycoprotein_median_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(! UniprotID %in% OG_glycoprotein_up_HEK293T$UniprotID) |> 
  filter(! UniprotID %in% OG_glycoprotein_down_HEK293T$UniprotID) |> 
  select(UniprotID) |> 
  mutate(category = "median", group = "Raw", count = 1)

OG_glycoprotein_protein_up_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> select(UniprotID) |> 
  mutate(category = "up", group = "Normalized", count = 1)

OG_glycoprotein_protein_down_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> select(UniprotID) |> 
  mutate(category = "down", group = "Normalized", count = 1)

OG_glycoprotein_protein_median_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |>
  filter(! UniprotID %in% OG_glycoprotein_protein_up_HEK293T$UniprotID) |> 
  filter(! UniprotID %in% OG_glycoprotein_protein_down_HEK293T$UniprotID) |> 
  select(UniprotID) |> 
  mutate(category = "median", group = "Normalized", count = 1)

OG_glycoprotein_protein_missing_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(!UniprotID %in% OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T$UniprotID) |> 
  select(UniprotID) |> 
  mutate(category = "missing", group = "Normalized", count = 1)

#combine
OG_glycoprotein_protein_total <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\glycoprotein_OGlcNAc_effect\\OG_glycoprotein_protein_total.xlsx"
)

OG_combine_upset <- OG_glycoprotein_protein_total |> 
  left_join(OG_glycoprotein_up_HEK293T, by = "UniprotID") |> 
  left_join(OG_glycoprotein_down_HEK293T, by = "UniprotID") |> 
  left_join(OG_glycoprotein_median_HEK293T, by = "UniprotID") |> 
  left_join(OG_glycoprotein_protein_up_HEK293T, by = "UniprotID") |> 
  left_join(OG_glycoprotein_protein_down_HEK293T, by = "UniprotID") |> 
  left_join(OG_glycoprotein_protein_median_HEK293T, by = "UniprotID") |> 
  mutate(
    up_glycoprotein_HEK293T = ifelse(is.na(up_glycoprotein_HEK293T), 0, 1),
    down_glycoprotein_HEK293T = ifelse(is.na(down_glycoprotein_HEK293T), 0, 1),
    median_glycoprotein_HEK293T = ifelse(is.na(median_glycoprotein_HEK293T), 0, 1),
    up_glycoprotein_protein_HEK293T = ifelse(is.na(up_glycoprotein_protein_HEK293T), 0, 1),
    down_glycoprotein_protein_HEK293T = ifelse(is.na(down_glycoprotein_protein_HEK293T), 0, 1),
    median_glycoprotein_protein_HEK293T = ifelse(is.na(median_glycoprotein_protein_HEK293T), 0, 1)
  )

OG_combine_upset_dataframe <- as.data.frame(OG_combine_upset)

upset(OG_combine_upset_dataframe, 
      sets = c("up_glycoprotein_HEK293T", "down_glycoprotein_HEK293T", "median_glycoprotein_HEK293T",
               "up_glycoprotein_protein_HEK293T", "down_glycoprotein_protein_HEK293T", "median_glycoprotein_protein_HEK293T"),
      main.bar.color = "black", 
      order.by = "freq", point.size = 5,
      text.scale = c(2, 1.5, 1.5, 1.5, 2, 2)
      )

#alluvial diagram
#combine
OG_combine_alluvial_HEK293T <- bind_rows(
  OG_glycoprotein_up_HEK293T,
  OG_glycoprotein_down_HEK293T,
  OG_glycoprotein_median_HEK293T,
  OG_glycoprotein_protein_up_HEK293T,
  OG_glycoprotein_protein_down_HEK293T,
  OG_glycoprotein_protein_median_HEK293T,
  OG_glycoprotein_protein_missing_HEK293T
) |> mutate(category = factor(category, levels = c("missing" ,"up", "median", "down")))

font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

alluvial_plot_OG_combine_HEK293T <- OG_combine_alluvial_HEK293T |> 
  ggplot(aes(x = factor(group, levels = c("Raw", "Normalized")), stratum = category, alluvium = UniprotID, y = count,
             fill = category, label = category)) +
  geom_flow(width = 0.5) +
  geom_stratum(width = 0.5) +
  geom_text(aes(color = category), stat = "stratum", size = 30) +
  labs(x = "", y = "# of glycoprotein") +
  scale_color_manual(
    values = c(
      "up" = "black",
      "median" = "black",
      "missing" = "black",
      "down" = "white"
    )
  ) +
  scale_fill_manual(
    values = c(
      "missing" = "white",
      "up" = "orange",
      "down" = Color_6,
      "median" = "gray"
    )
  ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 120),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.position = 'none'
  )

ggsave(filename = paste0(file_path, "alluvial_plot_OG_combine_HEK293T.tiff"), 
       device = "tiff",
       plot = alluvial_plot_OG_combine_HEK293T, 
       height = 2, width = 2, 
       units = "in", 
       dpi = 1200)

#generate data frame
#HepG2
OG_glycoprotein_up_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> select(UniprotID) |> 
  mutate(category = "up" , group = "Raw", count = 1)

OG_glycoprotein_down_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> select(UniprotID) |> 
  mutate(category = "down", group = "Raw", count = 1)

OG_glycoprotein_median_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(! UniprotID %in% OG_glycoprotein_up_HepG2$UniprotID) |> 
  filter(! UniprotID %in% OG_glycoprotein_down_HepG2$UniprotID) |> 
  select(UniprotID) |> 
  mutate(category = "median", group = "Raw", count = 1)

OG_glycoprotein_protein_up_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> select(UniprotID) |> 
  mutate(category = "up", group = "Normalized", count = 1)

OG_glycoprotein_protein_down_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> select(UniprotID) |> 
  mutate(category = "down", group = "Normalized", count = 1)

OG_glycoprotein_protein_median_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |>
  filter(! UniprotID %in% OG_glycoprotein_protein_up_HepG2$UniprotID) |> 
  filter(! UniprotID %in% OG_glycoprotein_protein_down_HepG2$UniprotID) |> 
  select(UniprotID) |> 
  mutate(category = "median", group = "Normalized", count = 1)

OG_glycoprotein_protein_missing_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(!UniprotID %in% OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2$UniprotID) |> 
  select(UniprotID) |> 
  mutate(category = "missing", group = "Normalized", count = 1)

#alluvial diagram
#combine
OG_combine_alluvial <- bind_rows(
  OG_glycoprotein_up_HepG2,
  OG_glycoprotein_down_HepG2,
  OG_glycoprotein_median_HepG2,
  OG_glycoprotein_protein_up_HepG2,
  OG_glycoprotein_protein_down_HepG2,
  OG_glycoprotein_protein_median_HepG2,
  OG_glycoprotein_protein_missing_HepG2
) |> mutate(category = factor(category, levels = c("missing" ,"up", "median", "down")))

font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

alluvial_plot_OG_combine_HepG2 <- OG_combine_alluvial |> 
  ggplot(aes(x = factor(group, levels = c("Raw", "Normalized")), stratum = category, alluvium = UniprotID, y = count,
             fill = category, label = category)) +
  geom_flow(width = 0.5) +
  geom_stratum(width = 0.5) +
  geom_text(aes(color = category), stat = "stratum", size = 30) +
  labs(x = "", y = "# of glycoprotein") +
  scale_color_manual(
    values = c(
      "up" = "black",
      "median" = "black",
      "missing" = "black",
      "down" = "white"
    )
  ) +
  scale_fill_manual(
    values = c(
      "missing" = "white",
      "up" = "orange",
      "down" = Color_6,
      "median" = "gray"
    )
  ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 120),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.position = 'none'
  )

ggsave(filename = paste0(file_path, "alluvial_plot_OG_combine_HepG2.tiff"), 
       device = "tiff",
       plot = alluvial_plot_OG_combine_HepG2, 
       height = 2, width = 2, 
       units = "in", 
       dpi = 1200)

#generate data frame
#Jurkat
OG_glycoprotein_up_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> select(UniprotID) |> 
  mutate(category = "up" , group = "Raw", count = 1)

OG_glycoprotein_down_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> select(UniprotID) |> 
  mutate(category = "down", group = "Raw", count = 1)

OG_glycoprotein_median_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(! UniprotID %in% OG_glycoprotein_up_Jurkat$UniprotID) |> 
  filter(! UniprotID %in% OG_glycoprotein_down_Jurkat$UniprotID) |> 
  select(UniprotID) |> 
  mutate(category = "median", group = "Raw", count = 1)

OG_glycoprotein_protein_up_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> select(UniprotID) |> 
  mutate(category = "up", group = "Normalized", count = 1)

OG_glycoprotein_protein_down_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> select(UniprotID) |> 
  mutate(category = "down", group = "Normalized", count = 1)

OG_glycoprotein_protein_median_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |>
  filter(! UniprotID %in% OG_glycoprotein_protein_up_Jurkat$UniprotID) |> 
  filter(! UniprotID %in% OG_glycoprotein_protein_down_Jurkat$UniprotID) |> 
  select(UniprotID) |> 
  mutate(category = "median", group = "Normalized", count = 1)

OG_glycoprotein_protein_missing_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(!UniprotID %in% OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat$UniprotID) |> 
  select(UniprotID) |> 
  mutate(category = "missing", group = "Normalized", count = 1)

#alluvial diagram
#combine
OG_combine_alluvial <- bind_rows(
  OG_glycoprotein_up_Jurkat,
  OG_glycoprotein_down_Jurkat,
  OG_glycoprotein_median_Jurkat,
  OG_glycoprotein_protein_up_Jurkat,
  OG_glycoprotein_protein_down_Jurkat,
  OG_glycoprotein_protein_median_Jurkat,
  OG_glycoprotein_protein_missing_Jurkat
) |> mutate(category = factor(category, levels = c("missing" ,"up", "median", "down")))

font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

alluvial_plot_OG_combine_Jurkat <- OG_combine_alluvial |> 
  ggplot(aes(x = factor(group, levels = c("Raw", "Normalized")), stratum = category, alluvium = UniprotID, y = count,
             fill = category, label = category)) +
  geom_flow(width = 0.5) +
  geom_stratum(width = 0.5) +
  geom_text(aes(color = category), stat = "stratum", size = 30) +
  labs(x = "", y = "# of glycoprotein") +
  scale_color_manual(
    values = c(
      "up" = "black",
      "median" = "black",
      "missing" = "black",
      "down" = "white"
    )
  ) +
  scale_fill_manual(
    values = c(
      "missing" = "white",
      "up" = "orange",
      "down" = Color_6,
      "median" = "gray"
    )
  ) +
  theme_classic() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 120),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.position = 'none'
  )

ggsave(filename = paste0(file_path, "alluvial_plot_OG_combine_Jurkat.tiff"), 
       device = "tiff",
       plot = alluvial_plot_OG_combine_Jurkat, 
       height = 2, width = 2, 
       units = "in", 
       dpi = 1200)
