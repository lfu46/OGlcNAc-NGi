#import packages
packages_names <- c("circlize", "scales", "ComplexHeatmap", "gridBase")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
#HepG2
OG_glycoprotein_up_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

OG_glycoprotein_down_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

#HEK293T
OG_glycoprotein_up_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

OG_glycoprotein_down_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

#Jurkat
OG_glycoprotein_up_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

OG_glycoprotein_down_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

#from down HEK293T to up Jurkat
downHEK293T_upJurkat <- semi_join(
  OG_glycoprotein_down_HEK293T, OG_glycoprotein_up_Jurkat, by = "UniprotID"
)

#from up HEK293T to down Jurkat
upHEK293T_downJurkat <- semi_join(
  OG_glycoprotein_up_HEK293T, OG_glycoprotein_down_Jurkat, by = "UniprotID"
)

#from down HEK293T to up HepG2
downHEK293T_upHepG2 <- semi_join(
  OG_glycoprotein_down_HEK293T, OG_glycoprotein_up_HepG2, by = "UniprotID"
)

#from up HEK293T to down HepG2
upHEK293T_downHepG2 <- semi_join(
  OG_glycoprotein_up_HEK293T, OG_glycoprotein_down_HepG2, by = "UniprotID"
)

#from up Jurkat to down HepG2
upJurkat_downHepG2 <- semi_join(
  OG_glycoprotein_up_Jurkat, OG_glycoprotein_down_HepG2, by = "UniprotID"
)

#from up HepG2 to down Jurkat
upHepG2_downJurkat <- semi_join(
  OG_glycoprotein_up_HepG2, OG_glycoprotein_down_Jurkat, by = "UniprotID"
)

#generate data frame for HEK293T_Jurkat
HEK293T_Jurkat <- bind_rows(downHEK293T_upJurkat, upHEK293T_downJurkat)

write_xlsx(HEK293T_Jurkat, path = paste0(file_path, "HEK293T_Jurkat.xlsx"))

HEK293T_Jurkat_adj <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\HEK293T_Jurkat.xlsx"
)

HEK293T_Jurkat_logFC <- HEK293T_Jurkat_adj |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  select(UniprotID, Gene_Name, logFC_HEK293T = logFC) |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  select(UniprotID, Gene_Name, logFC_HEK293T, logFC_Jurkat = logFC)

HEK293T_Jurkat_sl_tmm_logFC <- HEK293T_Jurkat_logFC |> 
  left_join(OG_glycoprotein_raw_sl_tmm_HEK293T, by = "UniprotID") |> 
  select(UniprotID, Gene_Name, logFC_HEK293T, logFC_Jurkat, ends_with("sl_tmm")) |> 
  mutate(
    log2_ratio1_HEK293T = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    log2_ratio2_HEK293T = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    log2_ratio3_HEK293T = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(UniprotID, Gene_Name, logFC_HEK293T, log2_ratio1_HEK293T:log2_ratio3_HEK293T, logFC_Jurkat) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_Jurkat, by = "UniprotID") |> 
  mutate(
    log2_ratio1_Jurkat = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    log2_ratio2_Jurkat = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    log2_ratio3_Jurkat = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(Gene_Name, logFC_HEK293T, log2_ratio1_HEK293T:log2_ratio3_HEK293T, logFC_Jurkat, log2_ratio1_Jurkat:log2_ratio3_Jurkat)

#complex heatmap
mat <- data.matrix(HEK293T_Jurkat_sl_tmm_logFC |> select(starts_with("log2_ratio")))
rownames(mat) <- HEK293T_Jurkat_sl_tmm_logFC$Gene_Name
colnames(mat) <- c("HEK293T 1", "HEK293T 2", "HEK293T 3", "Jurkat 1", "Jurkat 2", "Jurkat 3")

anno_matrix <- data.matrix(HEK293T_Jurkat_sl_tmm_logFC |> select(starts_with("log2_ratio")))

ha1 = rowAnnotation("log2(Tuni/Ctrl)" = anno_points(anno_matrix,
                   pch = c(1, 1, 1, 2, 2, 2), 
                   gp = gpar(col = c(Color_2, Color_2, Color_2, Color_4, Color_4, Color_4), lwd = 2),
                   ylim = c(-3, 3),
                   size = unit(4, "mm"),
                   width = unit(4, "cm"),
                   axis_param = list(
                     at = c(-3, 0, 3)
                   )
                   ))

h1 <- Heatmap(mat, name = "log2(Tuni/Ctrl)", show_row_names = TRUE, show_column_names = TRUE,
        right_annotation = ha1)

#generate data frame for HEK293T_HepG2
HEK293T_HepG2 <- bind_rows(downHEK293T_upHepG2, upHEK293T_downHepG2)

write_xlsx(HEK293T_HepG2, path = paste0(file_path, "HEK293T_HepG2.xlsx"))

HEK293T_HepG2_adj <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\HEK293T_HepG2.xlsx"
)

HEK293T_HepG2_logFC <- HEK293T_HepG2_adj |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  select(UniprotID, Gene_Name, logFC_HEK293T = logFC) |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  select(UniprotID, Gene_Name, logFC_HEK293T, logFC_HepG2 = logFC)

HEK293T_HepG2_sl_tmm_logFC <- HEK293T_HepG2_logFC |> 
  left_join(OG_glycoprotein_raw_sl_tmm_HEK293T, by = "UniprotID") |> 
  select(UniprotID, Gene_Name, logFC_HEK293T, logFC_HepG2, ends_with("sl_tmm")) |> 
  mutate(
    log2_ratio1_HEK293T = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    log2_ratio2_HEK293T = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    log2_ratio3_HEK293T = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(UniprotID, Gene_Name, logFC_HEK293T, log2_ratio1_HEK293T:log2_ratio3_HEK293T, logFC_HepG2) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_HepG2, by = "UniprotID") |> 
  mutate(
    log2_ratio1_HepG2 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    log2_ratio2_HepG2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    log2_ratio3_HepG2 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(Gene_Name, logFC_HEK293T, log2_ratio1_HEK293T:log2_ratio3_HEK293T, logFC_HepG2, log2_ratio1_HepG2:log2_ratio3_HepG2)

#complex heatmap
mat <- data.matrix(HEK293T_HepG2_sl_tmm_logFC |> select(starts_with("log2_ratio")))
rownames(mat) <- HEK293T_HepG2_sl_tmm_logFC$Gene_Name
colnames(mat) <- c("HEK293T 1", "HEK293T 2", "HEK293T 3", "HepG2 1", "HepG2 2", "HepG2 3")

col_mat <- colorRamp2(breaks = c(-3, 0, 3), colors = c("blue", "white", "red"))

anno_matrix <- data.matrix(HEK293T_HepG2_sl_tmm_logFC |> select(starts_with("log2_ratio")))

ha2 = rowAnnotation("log2(Tuni/Ctrl)" = anno_points(anno_matrix,
                                                    pch = c(1, 1, 1, 2, 2, 2), 
                                                    gp = gpar(col = c(Color_2, Color_2, Color_2, Color_3, Color_3, Color_3), lwd = 2),
                                                    ylim = c(-3, 3),
                                                    size = unit(4, "mm"),
                                                    width = unit(4, "cm"),
                                                    axis_param = list(
                                                      at = c(-3, 0, 3)
                                                    )
))

h2 <- Heatmap(mat, col = col_mat, 
              name = "log2(Tuni/Ctrl)", show_row_names = TRUE, show_column_names = TRUE,
              right_annotation = ha2)

#generate data frame for Jurkat_HepG2
Jurkat_HepG2 <- bind_rows(upHepG2_downJurkat, upJurkat_downHepG2)

write_xlsx(Jurkat_HepG2, path = paste0(file_path, "Jurkat_HepG2.xlsx"))

Jurkat_HepG2_adj <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\Jurkat_HepG2.xlsx"
)

Jurkat_HepG2_logFC <- Jurkat_HepG2_adj |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  select(UniprotID, Gene_Name, logFC_Jurkat = logFC) |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  select(UniprotID, Gene_Name, logFC_Jurkat, logFC_HepG2 = logFC)

Jurkat_HepG2_sl_tmm_logFC <- Jurkat_HepG2_logFC |> 
  left_join(OG_glycoprotein_raw_sl_tmm_Jurkat, by = "UniprotID") |> 
  select(UniprotID, Gene_Name, logFC_Jurkat, logFC_HepG2, ends_with("sl_tmm")) |> 
  mutate(
    log2_ratio1_Jurkat = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    log2_ratio2_Jurkat = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    log2_ratio3_Jurkat = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(UniprotID, Gene_Name, logFC_Jurkat, log2_ratio1_Jurkat:log2_ratio3_Jurkat, logFC_HepG2) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_HepG2, by = "UniprotID") |> 
  mutate(
    log2_ratio1_HepG2 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    log2_ratio2_HepG2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    log2_ratio3_HepG2 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(Gene_Name, logFC_Jurkat, log2_ratio1_Jurkat:log2_ratio3_Jurkat, logFC_HepG2, log2_ratio1_HepG2:log2_ratio3_HepG2)

#complex heatmap
mat <- data.matrix(Jurkat_HepG2_sl_tmm_logFC |> select(starts_with("log2_ratio")))
rownames(mat) <- Jurkat_HepG2_sl_tmm_logFC$Gene_Name
colnames(mat) <- c("Jurkat 1", "Jurkat 2", "Jurkat 3", "HepG2 1", "HepG2 2", "HepG2 3")

col_mat <- colorRamp2(breaks = c(-2, 0, 2), colors = c("blue", "white", "red"))

anno_matrix <- data.matrix(Jurkat_HepG2_sl_tmm_logFC |> select(starts_with("log2_ratio")))

ha3 = rowAnnotation("log2(Tuni/Ctrl)" = anno_points(anno_matrix,
                                                    pch = c(1, 1, 1, 2, 2, 2), 
                                                    gp = gpar(col = c(Color_4, Color_4, Color_4, Color_3, Color_3, Color_3), lwd = 2),
                                                    ylim = c(-2, 2),
                                                    size = unit(4, "mm"),
                                                    width = unit(4, "cm"),
                                                    axis_param = list(
                                                      at = c(-2, 0, 2)
                                                    )
))

h3 <- Heatmap(mat, col = col_mat, 
              name = "log2(Tuni/Ctrl)", show_row_names = TRUE, show_column_names = TRUE,
              right_annotation = ha3)
