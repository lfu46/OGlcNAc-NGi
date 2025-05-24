#import packages
packages_names <- c("tidyverse", "showtext", "writexl")
lapply(packages_names, require, character.only = TRUE)

#generate up or down regulated OG site for each cell line
#HEK293T
OG_site_protein_up_HEK293T <- OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(Index, UniprotID, combined_site, site_position)

write_xlsx(OG_site_protein_up_HEK293T, path = paste0(file_path, "OG_site_protein_up_HEK293T.xlsx"))

OG_site_protein_down_HEK293T <- OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(Index, UniprotID, combined_site, site_position)

write_xlsx(OG_site_protein_down_HEK293T, path = paste0(file_path, "OG_site_protein_down_HEK293T.xlsx"))

OG_site_protein_median_HEK293T <- OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(! Index %in% OG_site_protein_up_HEK293T$Index) |> 
  filter(! Index %in% OG_site_protein_down_HEK293T$Index) |> 
  select(Index, UniprotID, combined_site, site_position)

write_xlsx(OG_site_protein_median_HEK293T, path = paste0(file_path, "OG_site_protein_median_HEK293T.xlsx"))

#HepG2
OG_site_protein_up_HepG2 <- OG_site_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(Index, UniprotID, combined_site, site_position)

write_xlsx(OG_site_protein_up_HepG2, path = paste0(file_path, "OG_site_protein_up_HepG2.xlsx"))

OG_site_protein_down_HepG2 <- OG_site_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(Index, UniprotID, combined_site, site_position)

write_xlsx(OG_site_protein_down_HepG2, path = paste0(file_path, "OG_site_protein_down_HepG2.xlsx"))

OG_site_protein_median_HepG2 <- OG_site_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(! Index %in% OG_site_protein_up_HepG2$Index) |> 
  filter(! Index %in% OG_site_protein_down_HepG2$Index) |> 
  select(Index, UniprotID, combined_site, site_position)

write_xlsx(OG_site_protein_median_HepG2, path = paste0(file_path, "OG_site_protein_median_HepG2.xlsx"))

#Jurkat
OG_site_protein_up_Jurkat <- OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(Index, UniprotID, combined_site, site_position)

write_xlsx(OG_site_protein_up_Jurkat, path = paste0(file_path, "OG_site_protein_up_Jurkat.xlsx"))

OG_site_protein_down_Jurkat <- OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(Index, UniprotID, combined_site, site_position)

write_xlsx(OG_site_protein_down_Jurkat, path = paste0(file_path, "OG_site_protein_down_Jurkat.xlsx"))

OG_site_protein_median_Jurkat <- OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(! Index %in% OG_site_protein_up_Jurkat$Index) |> 
  filter(! Index %in% OG_site_protein_down_Jurkat$Index) |> 
  select(Index, UniprotID, combined_site, site_position)

write_xlsx(OG_site_protein_median_Jurkat, path = paste0(file_path, "OG_site_protein_median_Jurkat.xlsx"))

#OG site protein total
OG_site_protein_total <- bind_rows(
  OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T |> distinct(UniprotID),
  OG_site_protein_OGlcNAc_effect_Top_tb_HepG2 |> distinct(UniprotID),
  OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat |> distinct(UniprotID)
) |> distinct()

write_xlsx(OG_site_protein_total, path = paste0(file_path, "OG_site_protein_total.xlsx"))

#Minus-Average plot
#HEK293T
#up
OG_site_protein_up_logFC_sl_tmm_HEK293T <- OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(Index %in% OG_site_protein_up_HEK293T$Index) |> 
  left_join(OG_glycopeptide_raw_sl_tmm_HEK293T, by = "Index") |> 
  select(Index, UniprotID = UniprotID.x, combined_site = combined_site.x, logFC, adj.P.Val, ends_with("sl_tmm")) |> 
  mutate(
    sum_SN = Tuni_1_sl_tmm + Tuni_2_sl_tmm + Tuni_3_sl_tmm + Ctrl_4_sl_tmm + Ctrl_5_sl_tmm + Ctrl_6_sl_tmm
  )

#down
OG_site_protein_down_logFC_sl_tmm_HEK293T <- OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(Index %in% OG_site_protein_down_HEK293T$Index) |> 
  left_join(OG_glycopeptide_raw_sl_tmm_HEK293T, by = "Index") |> 
  select(Index, UniprotID = UniprotID.x, combined_site = combined_site.x, logFC, adj.P.Val, ends_with("sl_tmm")) |> 
  mutate(
    sum_SN = Tuni_1_sl_tmm + Tuni_2_sl_tmm + Tuni_3_sl_tmm + Ctrl_4_sl_tmm + Ctrl_5_sl_tmm + Ctrl_6_sl_tmm
  )

#median
OG_site_protein_median_logFC_sl_tmm_HEK293T <- OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(Index %in% OG_site_protein_median_HEK293T$Index) |> 
  left_join(OG_glycopeptide_raw_sl_tmm_HEK293T, by = "Index") |> 
  select(Index, UniprotID = UniprotID.x, combined_site = combined_site.x, logFC, adj.P.Val, ends_with("sl_tmm")) |> 
  mutate(
    sum_SN = Tuni_1_sl_tmm + Tuni_2_sl_tmm + Tuni_3_sl_tmm + Ctrl_4_sl_tmm + Ctrl_5_sl_tmm + Ctrl_6_sl_tmm
  )

#plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

MA_plot_OG_site_protein_HEK293T <- ggplot() +
  geom_point(data = OG_site_protein_median_logFC_sl_tmm_HEK293T, color = "gray", size = 1, alpha = 0.2,
             aes(x = logFC, y = log10(sum_SN))) +
  geom_point(data = OG_site_protein_up_logFC_sl_tmm_HEK293T, 
             aes(x = logFC, y = log10(sum_SN), size = -log10(adj.P.Val), color = -log10(adj.P.Val))) +
  geom_point(data = OG_site_protein_down_logFC_sl_tmm_HEK293T, 
             aes(x = logFC, y = log10(sum_SN), size = -log10(adj.P.Val), color = -log10(adj.P.Val))) +
  geom_vline(aes(xintercept = c(-0.5, 0.5)),  linetype= "dashed", color = "black") +
  coord_cartesian(xlim = c(-3, 3)) +
  ggtitle("HEK293T") +
  labs(x = expression(Δlog[2]*"FC"), y = "glycopeptide S/N") +
  scale_y_continuous(breaks = c(2, 3, 4, 5), labels = c('1E2', '1E3', '1E4', "1E5")) +
  scale_size(range = c(1, 3), name = expression(-log[10]*"(adj.P.Val)"), 
             breaks = c(3, 4), labels = c('3', '4')) +
  scale_color_gradient(low = alpha(Color_2, 0.5), high = Color_2, 
                       name = expression(-log[10]*"(adj.P.Val)"),
                       breaks = c(3, 4, 5), labels = c('3', '4', '5')) +
  theme_bw() +
  theme(
    title = element_text(size = 100, color = "black"),
    panel.grid.major = element_line(color = "gray", linewidth = 0.1),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.05),
    axis.text = element_text(size = 80, color = "black"),
    axis.title.x = element_text(size = 100, color = "black"),
    axis.title.y = element_text(size = 100, color = "black"),
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black")
  )

ggsave(paste0(file_path, "MA_plot_OG_site_protein_HEK293T.png"), 
       plot = MA_plot_OG_site_protein_HEK293T, height = 4, width = 5.5, units = c("in"), dpi = 600) 

#Minus-Average plot
#HepG2
#up
OG_site_protein_up_logFC_sl_tmm_HepG2 <- OG_site_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(Index %in% OG_site_protein_up_HepG2$Index) |> 
  left_join(OG_glycopeptide_raw_sl_tmm_HepG2, by = "Index") |> 
  select(Index, UniprotID = UniprotID.x, combined_site = combined_site.x, logFC, adj.P.Val, ends_with("sl_tmm")) |> 
  mutate(
    sum_SN = Tuni_1_sl_tmm + Tuni_2_sl_tmm + Tuni_3_sl_tmm + Ctrl_4_sl_tmm + Ctrl_5_sl_tmm + Ctrl_6_sl_tmm
  )

#down
OG_site_protein_down_logFC_sl_tmm_HepG2 <- OG_site_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(Index %in% OG_site_protein_down_HepG2$Index) |> 
  left_join(OG_glycopeptide_raw_sl_tmm_HepG2, by = "Index") |> 
  select(Index, UniprotID = UniprotID.x, combined_site = combined_site.x, logFC, adj.P.Val, ends_with("sl_tmm")) |> 
  mutate(
    sum_SN = Tuni_1_sl_tmm + Tuni_2_sl_tmm + Tuni_3_sl_tmm + Ctrl_4_sl_tmm + Ctrl_5_sl_tmm + Ctrl_6_sl_tmm
  )

#median
OG_site_protein_median_logFC_sl_tmm_HepG2 <- OG_site_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(Index %in% OG_site_protein_median_HepG2$Index) |> 
  left_join(OG_glycopeptide_raw_sl_tmm_HepG2, by = "Index") |> 
  select(Index, UniprotID = UniprotID.x, combined_site = combined_site.x, logFC, adj.P.Val, ends_with("sl_tmm")) |> 
  mutate(
    sum_SN = Tuni_1_sl_tmm + Tuni_2_sl_tmm + Tuni_3_sl_tmm + Ctrl_4_sl_tmm + Ctrl_5_sl_tmm + Ctrl_6_sl_tmm
  )

#plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

MA_plot_OG_site_protein_HepG2 <- ggplot() +
  geom_point(data = OG_site_protein_median_logFC_sl_tmm_HepG2, color = "gray", size = 1, alpha = 0.2,
             aes(x = logFC, y = log10(sum_SN))) +
  geom_point(data = OG_site_protein_up_logFC_sl_tmm_HepG2, 
             aes(x = logFC, y = log10(sum_SN), size = -log10(adj.P.Val), color = -log10(adj.P.Val))) +
  geom_point(data = OG_site_protein_down_logFC_sl_tmm_HepG2, 
             aes(x = logFC, y = log10(sum_SN), size = -log10(adj.P.Val), color = -log10(adj.P.Val))) +
  geom_vline(aes(xintercept = c(-0.5, 0.5)),  linetype= "dashed", color = "black") +
  coord_cartesian(xlim = c(-3, 3)) +
  ggtitle("HepG2") +
  labs(x = expression(Δlog[2]*"FC"), y = "glycopeptide S/N") +
  scale_y_continuous(breaks = c(2, 3, 4, 5), labels = c('1E2', '1E3', '1E4', "1E5")) +
  scale_size(range = c(1, 3), name = expression(-log[10]*"(adj.P.Val)"), 
             breaks = c(3, 4), labels = c('3', '4')) +
  scale_color_gradient(low = alpha(Color_3, 0.5), high = Color_3, 
                       name = expression(-log[10]*"(adj.P.Val)"),
                       breaks = c(3, 4, 5), labels = c('3', '4', '5')) +
  theme_bw() +
  theme(
    title = element_text(size = 100, color = "black"),
    panel.grid.major = element_line(color = "gray", linewidth = 0.1),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.05),
    axis.text = element_text(size = 80, color = "black"),
    axis.title.x = element_text(size = 100, color = "black"),
    axis.title.y = element_text(size = 100, color = "black"),
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black")
  )

ggsave(paste0(file_path, "MA_plot_OG_site_protein_HepG2.png"), 
       plot = MA_plot_OG_site_protein_HepG2, height = 4, width = 5.5, units = c("in"), dpi = 600)

#Minus-Average plot
#Jurkat
#up
OG_site_protein_up_logFC_sl_tmm_Jurkat <- OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(Index %in% OG_site_protein_up_Jurkat$Index) |> 
  left_join(OG_glycopeptide_raw_sl_tmm_Jurkat, by = "Index") |> 
  select(Index, UniprotID = UniprotID.x, combined_site = combined_site.x, logFC, adj.P.Val, ends_with("sl_tmm")) |> 
  mutate(
    sum_SN = Tuni_1_sl_tmm + Tuni_2_sl_tmm + Tuni_3_sl_tmm + Ctrl_4_sl_tmm + Ctrl_5_sl_tmm + Ctrl_6_sl_tmm
  )

#down
OG_site_protein_down_logFC_sl_tmm_Jurkat <- OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(Index %in% OG_site_protein_down_Jurkat$Index) |> 
  left_join(OG_glycopeptide_raw_sl_tmm_Jurkat, by = "Index") |> 
  select(Index, UniprotID = UniprotID.x, combined_site = combined_site.x, logFC, adj.P.Val, ends_with("sl_tmm")) |> 
  mutate(
    sum_SN = Tuni_1_sl_tmm + Tuni_2_sl_tmm + Tuni_3_sl_tmm + Ctrl_4_sl_tmm + Ctrl_5_sl_tmm + Ctrl_6_sl_tmm
  )

#median
OG_site_protein_median_logFC_sl_tmm_Jurkat <- OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(Index %in% OG_site_protein_median_Jurkat$Index) |> 
  left_join(OG_glycopeptide_raw_sl_tmm_Jurkat, by = "Index") |> 
  select(Index, UniprotID = UniprotID.x, combined_site = combined_site.x, logFC, adj.P.Val, ends_with("sl_tmm")) |> 
  mutate(
    sum_SN = Tuni_1_sl_tmm + Tuni_2_sl_tmm + Tuni_3_sl_tmm + Ctrl_4_sl_tmm + Ctrl_5_sl_tmm + Ctrl_6_sl_tmm
  )

#plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

MA_plot_OG_site_protein_Jurkat <- ggplot() +
  geom_point(data = OG_site_protein_median_logFC_sl_tmm_Jurkat, color = "gray", size = 1, alpha = 0.2,
             aes(x = logFC, y = log10(sum_SN))) +
  geom_point(data = OG_site_protein_up_logFC_sl_tmm_Jurkat, 
             aes(x = logFC, y = log10(sum_SN), size = -log10(adj.P.Val), color = -log10(adj.P.Val))) +
  geom_point(data = OG_site_protein_down_logFC_sl_tmm_Jurkat, 
             aes(x = logFC, y = log10(sum_SN), size = -log10(adj.P.Val), color = -log10(adj.P.Val))) +
  geom_vline(aes(xintercept = c(-0.5, 0.5)),  linetype= "dashed", color = "black") +
  coord_cartesian(xlim = c(-3, 3)) +
  ggtitle("Jurkat") +
  labs(x = expression(Δlog[2]*"FC"), y = "glycopeptide S/N") +
  scale_y_continuous(breaks = c(2, 3, 4, 5), labels = c('1E2', '1E3', '1E4', "1E5")) +
  scale_size(range = c(1, 3), name = expression(-log[10]*"(adj.P.Val)"), 
             breaks = c(3, 4), labels = c('3', '4')) +
  scale_color_gradient(low = alpha(Color_4, 0.5), high = Color_4, 
                       name = expression(-log[10]*"(adj.P.Val)"),
                       breaks = c(2, 3, 4), labels = c('2', '3', '4')) +
  theme_bw() +
  theme(
    title = element_text(size = 100, color = "black"),
    panel.grid.major = element_line(color = "gray", linewidth = 0.1),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.05),
    axis.text = element_text(size = 80, color = "black"),
    axis.title.x = element_text(size = 100, color = "black"),
    axis.title.y = element_text(size = 100, color = "black"),
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black")
  )

ggsave(paste0(file_path, "MA_plot_OG_site_protein_Jurkat.png"), 
       plot = MA_plot_OG_site_protein_Jurkat, height = 4, width = 5.5, units = c("in"), dpi = 600)
