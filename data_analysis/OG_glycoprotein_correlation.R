#import packages
packages_names <- c("tidyverse", "rstatix", "showtext", "ggpubr", "corrplot")
lapply(packages_names, require, character.only = TRUE)

#generate OG glycopeptide logFC overlap
tb1 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_common_total$UniprotID) |> 
  select(UniprotID, HepG2_logFC = logFC)

tb2 <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_common_total$UniprotID) |> 
  select(UniprotID, HEK293T_logFC = logFC)

tb3 <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(UniprotID %in% OG_glycoprotein_common_total$UniprotID) |> 
  select(UniprotID, Jurkat_logFC = logFC)

OG_glycoprotein_logFC_common <- tb1 |> 
  left_join(tb2, by = join_by(UniprotID == UniprotID)) |> 
  left_join(tb3, by = join_by(UniprotID == UniprotID))

#calculate person's r between cell lines
correlation_HepG2_HEK293T <- cor(OG_glycoprotein_logFC_common$HepG2_logFC, OG_glycoprotein_logFC_common$HEK293T_logFC, method = "pearson")
correlation_HEK293T_Jurkat <- cor(OG_glycoprotein_logFC_common$HEK293T_logFC, OG_glycoprotein_logFC_common$Jurkat_logFC, method = "pearson")

#correlation plot between different cell lines
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

correlation_plot_OG_glycoprotein_common <- ggplot(data = OG_glycoprotein_logFC_common) +
  geom_point(aes(x = HEK293T_logFC, y = HepG2_logFC), shape = 21, color = Color_3, size = 3) +
  geom_point(aes(x = HEK293T_logFC, y = Jurkat_logFC), shape = 21, color = Color_4, size = 3) +
  labs(x = "log2(Fold change) in HEK293T", y = "log2(Fold change) in HepG2 or Jurkat") +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  annotate("point", x = -2, y = 2, shape = 21, size = 2, color = Color_3) +
  annotate("text", x = -2, y = 2, size = 25, hjust = -0.12, color = "black", label = "HepG2 (r = 0.44)") +
  annotate("point", x = -2, y = 1.8, shape = 21, size = 2, color = Color_4) +
  annotate("text", x = -2, y = 1.8, size = 25, hjust = -0.12, color = "black", label = "Jurkat (r = 0.11)") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 80),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black")
  )

ggsave(filename = paste0(file_path, "correlation_plot_OG_glycoprotein_common.png"), plot = correlation_plot_OG_glycoprotein_common, height = 4, width = 4, dpi = 600)

#correlation matrix plot
OG_glycoprotein_raw_sl_tmm_common <- OG_glycoprotein_raw_sl_tmm_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_common_total$UniprotID) |> 
  mutate(
    ratio_1_HepG2 = Tuni_1_sl_tmm/Ctrl_4_sl_tmm,
    ratio_2_HepG2 = Tuni_2_sl_tmm/Ctrl_5_sl_tmm,
    ratio_3_HepG2 = Tuni_3_sl_tmm/Ctrl_6_sl_tmm
  ) |> 
  select(UniprotID, ends_with("HepG2")) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  mutate(
    ratio_1_HEK293T = Tuni_1_sl_tmm/Ctrl_4_sl_tmm,
    ratio_2_HEK293T = Tuni_2_sl_tmm/Ctrl_5_sl_tmm,
    ratio_3_HEK293T = Tuni_3_sl_tmm/Ctrl_6_sl_tmm
  ) |> 
  select(UniprotID, starts_with("ratio")) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  mutate(
    ratio_1_Jurkat = Tuni_1_sl_tmm/Ctrl_4_sl_tmm,
    ratio_2_Jurkat = Tuni_2_sl_tmm/Ctrl_5_sl_tmm,
    ratio_3_Jurkat = Tuni_3_sl_tmm/Ctrl_6_sl_tmm
  ) |> 
  select(UniprotID, starts_with("ratio"))

OG_glycorprotein_common_cor_matrix <- cor(OG_glycoprotein_raw_sl_tmm_common |> select(!UniprotID), method = "pearson")
rownames(OG_glycorprotein_common_cor_matrix) <- c(
                          "HepG2 rep1", 
                          "HepG2 rep2",
                          "HepG2 rep3",
                          "HEK293T rep1",
                          "HEK293T rep2",
                          "HEK293T rep3",
                          "Jurkat rep1",
                          "Jurkat rep2",
                          "Jurkat rep3")

colnames(OG_glycorprotein_common_cor_matrix) <- c(
                          "HepG2 rep1", 
                          "HepG2 rep2",
                          "HepG2 rep3",
                          "HEK293T rep1",
                          "HEK293T rep2",
                          "HEK293T rep3",
                          "Jurkat rep1",
                          "Jurkat rep2",
                          "Jurkat rep3")

corrplot(OG_glycorprotein_common_cor_matrix, type = "lower", method = "square",
         tl.cex = 1.2, tl.col = "black", tl.srt = 90,
         col = COL1('YlGn', n = 10),
         col.lim = c(0, 1),
         cl.cex = 0.8,
         is.corr = FALSE)

