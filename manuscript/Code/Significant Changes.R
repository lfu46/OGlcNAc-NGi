#import packages
packages_names <- c("circlize", "scales", "ComplexHeatmap", "gridBase")
lapply(packages_names, require, character.only = TRUE)

#Circular Heatmap
#generate data frame
#HepG2
OG_glycoprotein_up_sl_tmm_logFC_adjp_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "up", cell = "HepG2") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

OG_glycoprotein_down_sl_tmm_logFC_adjp_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "down", cell = "HepG2") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

OG_glycoprotein_median_sl_tmm_logFC_adjp_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(!(UniprotID %in% OG_glycoprotein_up_sl_tmm_logFC_adjp_HepG2$UniprotID)) |> 
  filter(!(UniprotID %in% OG_glycoprotein_down_sl_tmm_logFC_adjp_HepG2$UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "median", cell = "HepG2") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

#HEK293T
OG_glycoprotein_up_sl_tmm_logFC_adjp_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "up", cell = "HEK293T") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

OG_glycoprotein_down_sl_tmm_logFC_adjp_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "down", cell = "HEK293T") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

OG_glycoprotein_median_sl_tmm_logFC_adjp_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(!(UniprotID %in% OG_glycoprotein_up_sl_tmm_logFC_adjp_HEK293T$UniprotID)) |> 
  filter(!(UniprotID %in% OG_glycoprotein_down_sl_tmm_logFC_adjp_HEK293T$UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "median", cell = "HEK293T") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

#Jurkat
OG_glycoprotein_up_sl_tmm_logFC_adjp_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "up", cell = "Jurkat") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

OG_glycoprotein_down_sl_tmm_logFC_adjp_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "down", cell = "Jurkat") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

OG_glycoprotein_median_sl_tmm_logFC_adjp_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(!(UniprotID %in% OG_glycoprotein_up_sl_tmm_logFC_adjp_Jurkat$UniprotID)) |> 
  filter(!(UniprotID %in% OG_glycoprotein_down_sl_tmm_logFC_adjp_Jurkat$UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val) |> 
  mutate(category = "median", cell = "Jurkat") |> 
  arrange(adj.P.Val) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC, adj.P.Val, category, cell, ends_with("sl_tmm"))

#combine
OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total <- bind_rows(
  OG_glycoprotein_up_sl_tmm_logFC_adjp_HepG2,
  OG_glycoprotein_median_sl_tmm_logFC_adjp_HepG2,
  OG_glycoprotein_down_sl_tmm_logFC_adjp_HepG2,
  OG_glycoprotein_up_sl_tmm_logFC_adjp_HEK293T,
  OG_glycoprotein_median_sl_tmm_logFC_adjp_HEK293T,
  OG_glycoprotein_down_sl_tmm_logFC_adjp_HEK293T,
  OG_glycoprotein_up_sl_tmm_logFC_adjp_Jurkat,
  OG_glycoprotein_median_sl_tmm_logFC_adjp_Jurkat,
  OG_glycoprotein_down_sl_tmm_logFC_adjp_Jurkat
) |> 
  rowwise() |> 
  mutate(
    max_value = max(across(ends_with("sl_tmm"))),
    min_value = min(across(ends_with("sl_tmm")))
  ) |> 
  ungroup() |> 
  mutate(
    scaled_Tuni_1_sl_tmm = (Tuni_1_sl_tmm - min_value) * 2/(max_value - min_value) - 1,
    scaled_Tuni_2_sl_tmm = (Tuni_2_sl_tmm - min_value) * 2/(max_value - min_value) - 1,
    scaled_Tuni_3_sl_tmm = (Tuni_3_sl_tmm - min_value) * 2/(max_value - min_value) - 1,
    scaled_Ctrl_4_sl_tmm = (Ctrl_4_sl_tmm - min_value) * 2/(max_value - min_value) - 1,
    scaled_Ctrl_5_sl_tmm = (Ctrl_5_sl_tmm - min_value) * 2/(max_value - min_value) - 1,
    scaled_Ctrl_6_sl_tmm = (Ctrl_6_sl_tmm - min_value) * 2/(max_value - min_value) - 1
  ) |> 
  mutate(
    cell_line = ifelse(cell == "HEK293T", 1, ifelse(cell == "HepG2", 2, 3))
  )

write_xlsx(OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total, path = paste0(file_path, "OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total.xlsx"))

#import data
OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total.xlsx"
)

#generate data matrix
mat_sl_tmm <- data.matrix(OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total |> select(scaled_Tuni_1_sl_tmm:scaled_Ctrl_6_sl_tmm))

cell <- OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total$cell

cell_line <- OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total$cell_line

adjpVal <- OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total$adj.P.Val

category <- OG_glycoprotein_up_down_median_sl_tmm_logFC_adjp_total$category

#generate circos link data frame
#from down jurkat to up HEK293T
donwJurkat_upHEK293T <- semi_join(
  OG_glycoprotein_down_sl_tmm_logFC_adjp_Jurkat, OG_glycoprotein_up_sl_tmm_logFC_adjp_HEK293T,
  by = "UniprotID"
)

indices_from_downJurkat_to_upHEK293T_downJurkat <- which(OG_glycoprotein_down_sl_tmm_logFC_adjp_Jurkat$UniprotID %in% donwJurkat_upHEK293T$UniprotID)
indices_from_downJurkat_to_upHEK293T_upHEK293T <- which(OG_glycoprotein_up_sl_tmm_logFC_adjp_HEK293T$UniprotID %in% donwJurkat_upHEK293T$UniprotID)

link_from_downJurkat_to_upHEK293T <- tribble(
  ~ group1, ~ group1_index, ~ group2, ~ group2_index,
  "Jurkat", indices_from_downJurkat_to_upHEK293T_downJurkat, "HKE293T", indices_from_downJurkat_to_upHEK293T_upHEK293T
) |> unnest(cols = c(group1_index, group2_index))

#from down HEK293T to up Jurkat
downHEK293T_upJurkat <- semi_join(
  OG_glycoprotein_down_sl_tmm_logFC_adjp_HEK293T, OG_glycoprotein_up_sl_tmm_logFC_adjp_Jurkat,
  by = "UniprotID"
)

indices_from_downHEK293T_to_upJurkat_downHEK293T <- which(OG_glycoprotein_down_sl_tmm_logFC_adjp_HEK293T$UniprotID %in% downHEK293T_upJurkat$UniprotID)
indices_from_downHEK293T_to_upJurkat_upJurkat <- which(OG_glycoprotein_up_sl_tmm_logFC_adjp_Jurkat$UniprotID %in% downHEK293T_upJurkat$UniprotID)

link_from_downHEK293T_to_upJurkat <- tribble(
  ~ group1, ~ group1_index, ~ group2, ~ group2_index,
  "HEK293T", indices_from_downHEK293T_to_upJurkat_downHEK293T, "Jurkat", indices_from_downHEK293T_to_upJurkat_upJurkat
) |> unnest(cols = c(group1_index, group2_index))

#from down HepG2 to up HEK293T
downHepG2_upHEK293T <- semi_join(
  OG_glycoprotein_down_sl_tmm_logFC_adjp_HepG2, OG_glycoprotein_up_sl_tmm_logFC_adjp_HEK293T,
  by = "UniprotID"
)

indices_from_downHepG2_to_upHEK293T_downHepG2 <- which(OG_glycoprotein_down_sl_tmm_logFC_adjp_HepG2$UniprotID %in% downHepG2_upHEK293T$UniprotID)
indices_from_downHepG2_to_upHEK293T_upHEK293T <- which(OG_glycoprotein_up_sl_tmm_logFC_adjp_HEK293T$UniprotID %in% downHepG2_upHEK293T$UniprotID)

link_from_downHepG2_to_upHEK293T <- tribble(
  ~ group1, ~ group1_index, ~ group2, ~ group2_index,
  "HepG2", indices_from_downHepG2_to_upHEK293T_downHepG2, "HEK293T", indices_from_downHepG2_to_upHEK293T_upHEK293T
) |> unnest(cols = c(group1_index, group2_index))

#from down HEK293T to up HepG2
downHEK293T_upHepG2 <- semi_join(
  OG_glycoprotein_down_sl_tmm_logFC_adjp_HEK293T, OG_glycoprotein_up_sl_tmm_logFC_adjp_HepG2,
  by = "UniprotID"
)

indices_from_downHEK293T_to_upHepG2_downHEK293T <- which(OG_glycoprotein_down_sl_tmm_logFC_adjp_HEK293T$UniprotID %in% downHEK293T_upHepG2$UniprotID)
indices_from_donwHEK293T_to_upHepG2_upHepG2 <- which(OG_glycoprotein_up_sl_tmm_logFC_adjp_HepG2$UniprotID %in% downHEK293T_upHepG2$UniprotID)

link_from_downHEK293T_to_upHepG2 <- tribble(
  ~ group1, ~ group1_index, ~ group2, ~ group2_index,
  "HEK293T", indices_from_downHEK293T_to_upHepG2_downHEK293T, "HepG2", indices_from_donwHEK293T_to_upHepG2_upHepG2
) |> unnest(cols = c(group1_index, group2_index))

#from down HepG2 to up Jurkat
downHepG2_upJurkat <- semi_join(
  OG_glycoprotein_down_sl_tmm_logFC_adjp_HepG2, OG_glycoprotein_up_sl_tmm_logFC_adjp_Jurkat,
  by = "UniprotID"
)

indices_from_downHepG2_to_upJurkat_downHepG2 <- which(OG_glycoprotein_down_sl_tmm_logFC_adjp_HepG2$UniprotID %in% downHepG2_upJurkat$UniprotID)
indices_from_downHepG2_to_upJurkat_upJurkat <- which(OG_glycoprotein_up_sl_tmm_logFC_adjp_Jurkat$UniprotID %in% downHepG2_upJurkat$UniprotID)

link_from_downHepG2_to_upJurkat <- tribble(
  ~ group1, ~ group1_index, ~ group2, ~ group2_index,
  "HepG2", indices_from_downHepG2_to_upJurkat_downHepG2, "Jurkat", indices_from_downHepG2_to_upJurkat_upJurkat
) |> unnest(cols = c(group1_index, group2_index))

#from down Jurkat to up HepG2
downJurkat_upHepG2 <- semi_join(
  OG_glycoprotein_down_sl_tmm_logFC_adjp_Jurkat, OG_glycoprotein_up_sl_tmm_logFC_adjp_HepG2,
  by = "UniprotID"
)

indices_from_downJurkat_to_upHepG2_downJurkat <- which(OG_glycoprotein_down_sl_tmm_logFC_adjp_Jurkat$UniprotID %in% downJurkat_upHepG2$UniprotID)
indices_from_downJurkat_to_upHepG2_upHepG2 <- which(OG_glycoprotein_up_sl_tmm_logFC_adjp_HepG2$UniprotID %in% downJurkat_upHepG2$UniprotID)

link_from_downJurkat_to_upHepG2 <- tribble(
  ~ group1, ~ group1_index, ~ group2, ~ group2_index,
  "Jurkat", indices_from_downJurkat_to_upHepG2_downJurkat, "HepG2", indices_from_downJurkat_to_upHepG2_upHepG2
) |> unnest(cols = c(group1_index, group2_index))

#heatmap
circos.clear()
col_category <- c("up" = Color_5, "median" = "gray", "down" = Color_6)
circos.heatmap(category, split = cell_line, col = col_category, track.height = 0.03)

col_mat <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
circos.heatmap(mat_sl_tmm, col = col_mat, track.height = 0.1)

col_cell = c("HEK293T" = Color_2, "HepG2" = Color_3, "Jurkat" = Color_4)
circos.heatmap(cell, col = col_cell, track.height = 0.03)

col_adjpval_2 = colorRamp2(c(0, 0.001, 0.01, 0.01+1e-6, 0.05, 1), 
                           c(Color_1, alpha(Color_1, 0.8), alpha(Color_1, 0.6), 
                             alpha(Color_1, 0.6), alpha(Color_1, 0.4), "white"))
circos.heatmap(adjpVal, col = col_adjpval, track.height = 0.03)

#heatmap link 
#from down jurkat to up HEK293T
for(i in seq_len(nrow(link_from_downJurkat_to_upHEK293T))) {
  group1 = 3
  group2 = 1
  
  x1 = link_from_downJurkat_to_upHEK293T$group1_index[i] + 857
  x2 = link_from_downJurkat_to_upHEK293T$group2_index[i]
  
  circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_4, lwd = 3)
}

#from down HEK293T to up Jurkat
for(i in seq_len(nrow(link_from_downHEK293T_to_upJurkat))) {
  group1 = 1
  group2 = 3
  
  x1 = link_from_downHEK293T_to_upJurkat$group1_index[i] + 1121
  x2 = link_from_downHEK293T_to_upJurkat$group2_index[i]
  
  circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_2, lwd = 3)
}

#from down HepG2 to up HEK293T
for(i in seq_len(nrow(link_from_downHepG2_to_upHEK293T))) {
  group1 = 2
  group2 = 1
  
  x1 = link_from_downHepG2_to_upHEK293T$group1_index[i] + 349
  x2 = link_from_downHepG2_to_upHEK293T$group2_index[i]
  
  circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_3, lwd =3)
}

#from down HEK293T to up HepG2
for(i in seq_len(nrow(link_from_downHEK293T_to_upHepG2))) {
  group1 = 1
  group2 = 2
  
  x1 = link_from_downHEK293T_to_upHepG2$group1_index[i] + 1121
  x2 = link_from_downHEK293T_to_upHepG2$group2_index[i]
  
  circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_2, lwd = 3)
}

#from down HepG2 to up Jurkat
for(i in seq_len(nrow(link_from_downHepG2_to_upJurkat))) {
  group1 = 2
  group2 = 3
  
  x1 = link_from_downHepG2_to_upJurkat$group1_index[i] + 349
  x2 = link_from_downHepG2_to_upJurkat$group2_index[i]
  
  circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_3, lwd = 3)
}

#from down Jurkat to up HepG2
for(i in seq_len(nrow(link_from_downJurkat_to_upHepG2))) {
  group1 = 3
  group2 = 2
  
  x1 = link_from_downJurkat_to_upHepG2$group1_index[i] + 857
  x2 = link_from_downJurkat_to_upHepG2$group2_index[i]
  
  circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_4, lwd =3)
}

#circular heatmap wrapper
circlize_plot = function() {
  circos.heatmap(category, split = cell_line, col = col_category, track.height = 0.03)
  circos.heatmap(mat_sl_tmm, col = col_mat, track.height = 0.1)
  circos.heatmap(cell, col = col_cell, track.height = 0.03)
  circos.heatmap(adjpVal, col = col_adjpval_2, track.height = 0.03)
  
  #heatmap link 
  #from down jurkat to up HEK293T
  for(i in seq_len(nrow(link_from_downJurkat_to_upHEK293T))) {
    group1 = 3
    group2 = 1
    
    x1 = link_from_downJurkat_to_upHEK293T$group1_index[i] + 857
    x2 = link_from_downJurkat_to_upHEK293T$group2_index[i]
    
    circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_4, lwd = 3)
  }
  
  #from down HEK293T to up Jurkat
  for(i in seq_len(nrow(link_from_downHEK293T_to_upJurkat))) {
    group1 = 1
    group2 = 3
    
    x1 = link_from_downHEK293T_to_upJurkat$group1_index[i] + 1121
    x2 = link_from_downHEK293T_to_upJurkat$group2_index[i]
    
    circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_2, lwd = 3)
  }
  
  #from down HepG2 to up HEK293T
  for(i in seq_len(nrow(link_from_downHepG2_to_upHEK293T))) {
    group1 = 2
    group2 = 1
    
    x1 = link_from_downHepG2_to_upHEK293T$group1_index[i] + 349
    x2 = link_from_downHepG2_to_upHEK293T$group2_index[i]
    
    circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_3, lwd =3)
  }
  
  #from down HEK293T to up HepG2
  for(i in seq_len(nrow(link_from_downHEK293T_to_upHepG2))) {
    group1 = 1
    group2 = 2
    
    x1 = link_from_downHEK293T_to_upHepG2$group1_index[i] + 1121
    x2 = link_from_downHEK293T_to_upHepG2$group2_index[i]
    
    circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_2, lwd = 3)
  }
  
  #from down HepG2 to up Jurkat
  for(i in seq_len(nrow(link_from_downHepG2_to_upJurkat))) {
    group1 = 2
    group2 = 3
    
    x1 = link_from_downHepG2_to_upJurkat$group1_index[i] + 349
    x2 = link_from_downHepG2_to_upJurkat$group2_index[i]
    
    circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_3, lwd = 3)
  }
  
  #from down Jurkat to up HepG2
  for(i in seq_len(nrow(link_from_downJurkat_to_upHepG2))) {
    group1 = 3
    group2 = 2
    
    x1 = link_from_downJurkat_to_upHepG2$group1_index[i] + 857
    x2 = link_from_downJurkat_to_upHepG2$group2_index[i]
    
    circos.link(group1, x1 - 0.5, group2, x2 - 0.5, col = Color_4, lwd =3)
  }
  
  circos.clear()
}

#generate legend for circular heatmap
lgd_category <- Legend(title = "Category", at = names(col_category), legend_gp = gpar(fill = col_category))
lgd_mat <- Legend(title = "Z-score", col_fun = col_mat)
lgd_cell <- Legend(title = "Cell", at = names(col_cell), legend_gp = gpar(fill = col_cell))
lgd_adjpval <- Legend(title = "adj.P.Value", col_fun = col_adjpval_2, at = c(0, 0.001, 0.01, 0.05, 1),
                      break_dist = c(1, 1, 3, 3))

#combine circular heatmap and legends
plot.new()
circle_size = unit(1, "snpc")

pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circlize_plot()
upViewport()

h = dev.size()[2]
lgd_list <- packLegend(lgd_mat, lgd_cell, lgd_category, lgd_adjpval, 
                       max_height = unit(0.9*h, "inch"))
draw(lgd_list, x = circle_size, just = "left")

#Heatmap HEK293T vs. Jurkat
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

#Bar plot
#generate data frame
#HepG2
OG_glycoprotein_up_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_up_HepG2, path = paste0(file_path, "OG_glycoprotein_up_HepG2.xlsx"))

OG_glycoprotein_down_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_down_HepG2, path = paste0(file_path, "OG_glycoprotein_down_HepG2.xlsx"))

#HEK293T
OG_glycoprotein_up_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_up_HEK293T, path = paste0(file_path, "OG_glycoprotein_up_HEK293T.xlsx"))

OG_glycoprotein_down_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_down_HEK293T, path = paste0(file_path, "OG_glycoprotein_down_HEK293T.xlsx"))

#Jurkat
OG_glycoprotein_up_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_up_Jurkat, path = paste0(file_path, "OG_glycoprotein_up_Jurkat.xlsx"))

OG_glycoprotein_down_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_down_Jurkat, path = paste0(file_path, "OG_glycoprotein_down_Jurkat.xlsx"))

#import data
gene_ontology_OG_glycoprotein_up_down_specific_cell_up_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_up_down_specific_cell.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "up_HepG2"
)
gene_ontology_OG_glycoprotein_up_down_specific_cell_down_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_up_down_specific_cell.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "down_Jurkat"
)
gene_ontology_OG_glycoprotein_up_down_specific_cell_up_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_up_down_specific_cell.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "up_HEK293T"
)

#bar plot
#up HepG2
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

barplot_OG_glycoprotein_up_gene_ontology_HepG2 <- 
  gene_ontology_OG_glycoprotein_up_down_specific_cell_up_HepG2 |> 
  ggplot(aes(x = -log10(P.Value), y = fct_reorder(Term, -log10(P.Value)))) +
  geom_bar(stat = "identity", fill = Color_3, width = 0.7, alpha = 0.5) +
  labs(x = expression(-log[10]*"("*paste(italic(P), " Value")*")"), y = "") +
  geom_text(aes(label = Term, x = 0), hjust = 0, size = 15) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 40),
    axis.text.x = element_text(color = "black", size = 40),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "in")
  )

ggsave(filename = paste0(file_path, "barplot_OG_glycoprotein_up_gene_ontology_HepG2.png"), 
       plot = barplot_OG_glycoprotein_up_gene_ontology_HepG2, height = 1.5, width = 2, units = c("in"), dpi = 600)

#down Jurkat
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

barplot_OG_glycoprotein_down_gene_ontology_Jurkat <- 
  gene_ontology_OG_glycoprotein_up_down_specific_cell_down_Jurkat |> 
  ggplot(aes(x = -log10(P.Value), y = fct_reorder(Term, -log10(P.Value)))) +
  geom_bar(stat = "identity", fill = Color_4, width = 0.7, alpha = 0.5) +
  labs(x = expression(-log[10]*"("*paste(italic(P), " Value")*")"), y = "") +
  geom_text(aes(label = Term, x = 0), hjust = 0, size = 15) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 40),
    axis.text.x = element_text(color = "black", size = 40),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "in")
  )

ggsave(filename = paste0(file_path, "barplot_OG_glycoprotein_down_gene_ontology_Jurkat.png"), 
       plot = barplot_OG_glycoprotein_down_gene_ontology_Jurkat, height = 1.5, width = 2, units = c("in"), dpi = 600)

#up HEK293T
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

barplot_OG_glycoprotein_up_gene_ontology_HEK293T <- 
  gene_ontology_OG_glycoprotein_up_down_specific_cell_up_HEK293T |> 
  ggplot(aes(x = -log10(P.Value), y = fct_reorder(Term, -log10(P.Value)))) +
  geom_bar(stat = "identity", fill = Color_2, width = 0.7, alpha = 0.5) +
  labs(x = expression(-log[10]*"("*paste(italic(P), " Value")*")"), y = "") +
  geom_text(aes(label = Term, x = 0), hjust = 0, size = 15) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 40),
    axis.text.x = element_text(color = "black", size = 40),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "in")
  )

ggsave(filename = paste0(file_path, "barplot_OG_glycoprotein_up_gene_ontology_HEK293T.png"), 
       plot = barplot_OG_glycoprotein_up_gene_ontology_HEK293T, height = 1.5, width = 2, units = c("in"), dpi = 600)
