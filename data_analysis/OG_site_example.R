#import packages
packages_names <- c("tidyverse", "rstatix", "showtext", "ggpubr", "UpSetR")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
#HepG2
OG_site_up_HepG2 <- OG_glycopeptide_Top_tb_HepG2 |> 
  filter(!is.na(combined_site)) |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(Index) |> mutate(up_HepG2 = 1)

write_xlsx(OG_site_up_HepG2, path = paste0(file_path, "OG_site_up_HepG2.xlsx"))
  
OG_site_down_HepG2 <- OG_glycopeptide_Top_tb_HepG2 |> 
  filter(!is.na(combined_site)) |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(Index) |> mutate(down_HepG2 = 1)

write_xlsx(OG_site_down_HepG2, path = paste0(file_path, "OG_site_down_HepG2.xlsx"))

OG_site_median_HepG2 <- OG_glycopeptide_Top_tb_HepG2 |> 
  filter(!is.na(combined_site)) |> 
  filter(! Index %in% OG_site_up_HepG2) |> 
  filter(! Index %in% OG_site_down_HepG2) |> 
  select(Index) |> mutate(median_HepG2 = 1)

#HEK293T
OG_site_up_HEK293T <- OG_glycopeptide_Top_tb_HEK293T |> 
  filter(!is.na(combined_site)) |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(Index) |> mutate(up_HEK293T = 1)

write_xlsx(OG_site_up_HEK293T, path = paste0(file_path, "OG_site_up_HEK293T.xlsx"))

OG_site_down_HEK293T <- OG_glycopeptide_Top_tb_HEK293T |> 
  filter(!is.na(combined_site)) |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(Index) |> mutate(down_HEK293T = 1)

write_xlsx(OG_site_down_HEK293T, path = paste0(file_path, "OG_site_down_HEK293T.xlsx"))

OG_site_median_HEK293T <- OG_glycopeptide_Top_tb_HEK293T |> 
  filter(!is.na(combined_site)) |> 
  filter(! Index %in% OG_site_up_HEK293T) |> 
  filter(! Index %in% OG_site_down_HEK293T) |> 
  select(Index) |> mutate(median_HEK293T = 1)

#Jurkat
OG_site_up_Jurkat <- OG_glycopeptide_Top_tb_Jurkat |> 
  filter(!is.na(combined_site)) |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(Index) |> mutate(up_Jurkat = 1)

write_xlsx(OG_site_up_Jurkat, path = paste0(file_path, "OG_site_up_Jurkat.xlsx"))

OG_site_down_Jurkat <- OG_glycopeptide_Top_tb_Jurkat |> 
  filter(!is.na(combined_site)) |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(Index) |> mutate(down_Jurkat = 1)

write_xlsx(OG_site_down_Jurkat, path = paste0(file_path, "OG_site_down_Jurkat.xlsx"))

OG_site_median_Jurkat <- OG_glycopeptide_Top_tb_Jurkat |> 
  filter(!is.na(combined_site)) |> 
  filter(! Index %in% OG_site_up_Jurkat) |> 
  filter(! Index %in% OG_site_down_Jurkat) |> 
  select(Index) |> mutate(median_Jurkat = 1)

#OG site total
OG_site_total <- bind_rows(
  OG_glycopeptide_Top_tb_HEK293T |> filter(!is.na(combined_site)) |> select(Index),
  OG_glycopeptide_Top_tb_HepG2 |> filter(!is.na(combined_site)) |> select(Index),
  OG_glycopeptide_Top_tb_Jurkat |> filter(!is.na(combined_site)) |> select(Index)
) |> distinct()

write_xlsx(OG_site_total, path = paste0(file_path, "OG_site_total.xlsx"))

OG_site_logFC_total <- bind_rows(
  OG_glycopeptide_Top_tb_HepG2 |> filter(!is.na(combined_site)) |> mutate(cell = "HepG2"),
  OG_glycopeptide_Top_tb_HEK293T |> filter(!is.na(combined_site)) |> mutate(cell = "HEK293T"),
  OG_glycopeptide_Top_tb_Jurkat |> filter(!is.na(combined_site)) |> mutate(cell = "Jurakt")
)

write_xlsx(OG_site_logFC_total, path = paste0(file_path, "OG_site_logFC_total.xlsx"))

#generate data frame
OG_site_up_down_median_upset <- OG_site_total |> 
  left_join(OG_site_up_HepG2, by = "Index") |> 
  left_join(OG_site_down_HepG2, by = "Index") |> 
  left_join(OG_site_median_HepG2, by = "Index") |> 
  left_join(OG_site_up_HEK293T, by = "Index") |> 
  left_join(OG_site_down_HEK293T, by = "Index") |> 
  left_join(OG_site_median_HEK293T, by = "Index") |> 
  left_join(OG_site_up_Jurkat, by = "Index") |> 
  left_join(OG_site_down_Jurkat, by = "Index") |> 
  left_join(OG_site_median_Jurkat, by = "Index") |> 
  mutate(
    up_HepG2 = ifelse(is.na(up_HepG2), 0, 1),
    down_HepG2 = ifelse(is.na(down_HepG2), 0, 1),
    median_HepG2 = ifelse(is.na(median_HepG2), 0, 1),
    up_HEK293T = ifelse(is.na(up_HEK293T), 0, 1),
    down_HEK293T = ifelse(is.na(down_HEK293T), 0, 1),
    median_HEK293T = ifelse(is.na(median_HEK293T), 0, 1),
    up_Jurkat = ifelse(is.na(up_Jurkat), 0, 1),
    down_Jurkat = ifelse(is.na(down_Jurkat), 0, 1),
    median_Jurkat = ifelse(is.na(median_Jurkat), 0, 1)
  )

OG_site_up_down_median_upset_dataframe <- as.data.frame(OG_site_up_down_median_upset)

upset(OG_site_up_down_median_upset_dataframe, 
      sets = c(
        "up_HepG2", "down_HepG2",
        "up_HEK293T", "down_HEK293T", 
        "up_Jurkat", "down_Jurkat"
      ),
      point.size = 5, keep.order = TRUE, text.scale = c(1.5, 1.5, 1.5, 1.5, 2, 2))

#example
OG_site_example <- c(
  "Q9HCE1_12_24_1_1_S19",
  "P61019_152_165_1_1_S154",
  "P52948_1053_1069_1_1_S1064"
)

OG_site_example_sl_tmm_logFC <- bind_rows(
  OG_glycopeptide_raw_sl_tmm_HepG2 |> filter(Index %in% OG_site_example) |> 
    left_join(OG_glycopeptide_Top_tb_HepG2, by = "Index") |> 
    select(Index, combined_site = combined_site.x, UniprotID = UniprotID.x, logFC, ends_with("sl_tmm"), ) |> 
    mutate(
      ratio1 = Tuni_1_sl_tmm / Ctrl_4_sl_tmm,
      ratio2 = Tuni_2_sl_tmm / Ctrl_5_sl_tmm,
      ratio3 = Tuni_3_sl_tmm / Ctrl_6_sl_tmm,
    ) |> 
    mutate(
      log2_ratio1 = log2(ratio1),
      log2_ratio2 = log2(ratio2),
      log2_ratio3 = log2(ratio3)
    ) |> 
    select(Index, combined_site, UniprotID, logFC_site = logFC, starts_with("log2_")) |> 
    left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
    select(Index, combined_site, UniprotID, logFC_site, starts_with("log2_"), logFC_glycoprotein = logFC) |> 
    mutate(cell = "HepG2"),
  
  OG_glycopeptide_raw_sl_tmm_HEK293T |> filter(Index %in% OG_site_example) |> 
    left_join(OG_glycopeptide_Top_tb_HEK293T, by = "Index") |> 
    select(Index, combined_site = combined_site.x, UniprotID = UniprotID.x, logFC, ends_with("sl_tmm")) |> 
    mutate(
      ratio1 = Tuni_1_sl_tmm / Ctrl_4_sl_tmm,
      ratio2 = Tuni_2_sl_tmm / Ctrl_5_sl_tmm,
      ratio3 = Tuni_3_sl_tmm / Ctrl_6_sl_tmm,
    ) |> 
    mutate(
      log2_ratio1 = log2(ratio1),
      log2_ratio2 = log2(ratio2),
      log2_ratio3 = log2(ratio3)
    ) |> 
    select(Index, combined_site, UniprotID, logFC_site = logFC, starts_with("log2_")) |> 
    left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
    select(Index, combined_site, UniprotID, logFC_site, starts_with("log2_"), logFC_glycoprotein = logFC) |> 
    mutate(cell = "HEK293T"),
  
  OG_glycopeptide_raw_sl_tmm_Jurkat |> filter(Index %in% OG_site_example) |> 
    left_join(OG_glycopeptide_Top_tb_Jurkat, by = "Index") |> 
    select(Index, combined_site = combined_site.x, UniprotID = UniprotID.x, logFC, ends_with("sl_tmm")) |> 
    mutate(
      ratio1 = Tuni_1_sl_tmm / Ctrl_4_sl_tmm,
      ratio2 = Tuni_2_sl_tmm / Ctrl_5_sl_tmm,
      ratio3 = Tuni_3_sl_tmm / Ctrl_6_sl_tmm,
    ) |> 
    mutate(
      log2_ratio1 = log2(ratio1),
      log2_ratio2 = log2(ratio2),
      log2_ratio3 = log2(ratio3)
    ) |> 
    select(Index, combined_site, UniprotID, logFC_site = logFC, starts_with("log2_")) |> 
    left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
    select(Index, combined_site, UniprotID, logFC_site, starts_with("log2_"), logFC_glycoprotein = logFC) |> 
    mutate(cell = "Jurkat")
) |> pivot_longer(cols = starts_with("log2"), names_to = 'Exp', values_to = 'triplicate')

#plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

OG_site_example_plot <- OG_site_example_sl_tmm_logFC |> 
  ggplot() +
  geom_bar(aes(x = combined_site, y = logFC_site, fill = cell), stat = "identity", position = position_dodge(width = 0.92))  +
  geom_hline(aes(yintercept = logFC_glycoprotein, color = cell), linetype = "dotdash", linewidth = 1) +
  geom_point(aes(x = combined_site, y = triplicate, fill = cell), size = 3, shape = 21, color = "black", position = position_dodge(width = 0.92)) +
  facet_wrap(vars(UniprotID), nrow = 3, scales = "free",
             labeller = labeller(UniprotID = c(
               "Q9HCE1" = "MOV10",
               "P61019" = "RAB2A",
               "P52948" = "NUP98"
             ))) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(values = c(Color_2, Color_3, Color_4)) +
  scale_color_manual(values = c(Color_2, Color_3, Color_4)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 90, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    strip.text = element_text(size = 100, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_blank()
  )

ggsave(filename = paste0(file_path, "OG_site_example_plot.png"), plot = OG_site_example_plot, height = 6, width = 4, dpi = 600)
