#import packages
packages_names <- c("tidyverse", "limma", "edgeR", "showtext")
lapply(packages_names, require, character.only = TRUE)

#generate OG glycopeptide logFC data frame
PCA_OG_glycoprotein_HepG2 <- OG_glycoprotein_raw_sl_tmm_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_common_total$UniprotID) |> 
  select(UniprotID, ends_with("sl_tmm")) |> 
  mutate(
    ratio1 = Tuni_1_sl_tmm / Ctrl_4_sl_tmm,
    ratio2 = Tuni_2_sl_tmm / Ctrl_5_sl_tmm,
    ratio3 = Tuni_3_sl_tmm / Ctrl_6_sl_tmm,
  ) |> 
  mutate(
    log2_ratio1_HepG2 = log2(ratio1),
    log2_ratio2_HepG2 = log2(ratio2),
    log2_ratio3_HepG2 = log2(ratio3)
  ) |>
  select(UniprotID, log2_ratio1_HepG2, log2_ratio2_HepG2, log2_ratio3_HepG2)

PCA_OG_glycoprotein_HEK293T <- OG_glycoprotein_raw_sl_tmm_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_common_total$UniprotID) |> 
  select(UniprotID, ends_with("sl_tmm")) |> 
  mutate(
    ratio1 = Tuni_1_sl_tmm / Ctrl_4_sl_tmm,
    ratio2 = Tuni_2_sl_tmm / Ctrl_5_sl_tmm,
    ratio3 = Tuni_3_sl_tmm / Ctrl_6_sl_tmm,
  ) |> 
  mutate(
    log2_ratio1_HEK293T = log2(ratio1),
    log2_ratio2_HEK293T = log2(ratio2),
    log2_ratio3_HEK293T = log2(ratio3)
  ) |>
  select(UniprotID, log2_ratio1_HEK293T, log2_ratio2_HEK293T, log2_ratio3_HEK293T)

PCA_OG_glycoprotein_Jurkat <- OG_glycoprotein_raw_sl_tmm_Jurkat |> 
  filter(UniprotID %in% OG_glycoprotein_common_total$UniprotID) |> 
  select(UniprotID, ends_with("sl_tmm")) |> 
  mutate(
    ratio1 = Tuni_1_sl_tmm / Ctrl_4_sl_tmm,
    ratio2 = Tuni_2_sl_tmm / Ctrl_5_sl_tmm,
    ratio3 = Tuni_3_sl_tmm / Ctrl_6_sl_tmm,
  ) |> 
  mutate(
    log2_ratio1_Jurkat = log2(ratio1),
    log2_ratio2_Jurkat = log2(ratio2),
    log2_ratio3_Jurkat = log2(ratio3)
  ) |>
  select(UniprotID, log2_ratio1_Jurkat, log2_ratio2_Jurkat, log2_ratio3_Jurkat)

PCA_OG_glycoprotein_sl_tmm_common <- PCA_OG_glycoprotein_HepG2 |> 
  left_join(PCA_OG_glycoprotein_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  left_join(PCA_OG_glycoprotein_Jurkat, by = join_by(UniprotID == UniprotID))

#principle component analysis
df_t_OG_glycoprotein_sl_tmm_common <- t(PCA_OG_glycoprotein_sl_tmm_common |> select(log2_ratio1_HepG2:log2_ratio3_Jurkat))

pca_result_OG_glycoprotein_common <- prcomp(df_t_OG_glycoprotein_sl_tmm_common, scale. = TRUE)

variance_explained_OG_glycoprotein_common <- pca_result_OG_glycoprotein_common$sdev^2 / sum(pca_result_OG_glycoprotein_common$sdev^2)
percent_variance_OG_glycoprotein_common <- round(variance_explained_OG_glycoprotein_common * 100, 2)

pca_data_OG_glycoprotein_common <- as.data.frame(pca_result_OG_glycoprotein_common$x)

pca_data_OG_glycoprotein_common$Exp<- factor(c("HepG2", "HepG2", "HepG2", "HEK293T", "HEK293T", "HEK293T", "Jurkat", "Jurkat", "Jurkat"))

font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

PCA_OG_glycoprotein_common <- ggplot(pca_data_OG_glycoprotein_common, aes(x = PC1, y = PC2, color = Exp)) +
  geom_point(shape = 21, size = 3, stroke = 1.5) +
  labs(
    x = paste0("PC1 (", percent_variance_OG_glycoprotein_common[1], "%)"),
    y = paste0("PC2 (", percent_variance_OG_glycoprotein_common[2], "%)")
  ) +
  scale_color_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.text = element_text(size = 80),
    legend.title = element_blank()
  )

ggsave(filename = paste0(file_path, "PCA_OG_glycoprotein_common.png"), plot = PCA_OG_glycoprotein_common, height = 4, width = 4, units = c("in"), dpi = 600)
