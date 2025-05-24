#import packages
packages_names <- c("UpSetR")
lapply(packages_names, require, character.only = TRUE)

#generate glycorprotein list for each category
OG_glycoprotein_up_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

OG_glycoprotein_down_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

OG_glycoprotein_median_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(!(UniprotID %in% OG_glycoprotein_up_HepG2$UniprotID)) |> 
  filter(!(UniprotID %in% OG_glycoprotein_down_HepG2$UniprotID)) |> 
  select(UniprotID)

OG_glycoprotein_up_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

OG_glycoprotein_down_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

OG_glycoprotein_median_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(!(UniprotID %in% OG_glycoprotein_up_HEK293T$UniprotID)) |> 
  filter(!(UniprotID %in% OG_glycoprotein_down_HEK293T$UniprotID)) |> 
  select(UniprotID)

OG_glycoprotein_up_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

OG_glycoprotein_down_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

OG_glycoprotein_median_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(!(UniprotID %in% OG_glycoprotein_up_Jurkat$UniprotID)) |> 
  filter(!(UniprotID %in% OG_glycoprotein_down_Jurkat$UniprotID)) |> 
  select(UniprotID)

OG_glycoprotein_up_HepG2 <- OG_glycoprotein_up_HepG2 |> mutate(up_HepG2 = 1)
OG_glycoprotein_down_HepG2 <- OG_glycoprotein_down_HepG2 |> mutate(down_HepG2 = 1)
OG_glycoprotein_median_HepG2 <- OG_glycoprotein_median_HepG2 |> mutate(median_HepG2 = 1)
OG_glycoprotein_up_HEK293T <- OG_glycoprotein_up_HEK293T |> mutate(up_HEK293T = 1)
OG_glycoprotein_down_HEK293T <- OG_glycoprotein_down_HEK293T |> mutate(down_HEK293T = 1)
OG_glycoprotein_median_HEK293T <- OG_glycoprotein_median_HEK293T |> mutate(median_HEK293T = 1)
OG_glycoprotein_up_Jurkat <- OG_glycoprotein_up_Jurkat |> mutate(up_Jurkat = 1)
OG_glycoprotein_down_Jurkat <- OG_glycoprotein_down_Jurkat |> mutate(down_Jurkat = 1)
OG_glycoprotein_median_Jurkat <- OG_glycoprotein_median_Jurkat |> mutate(median_Jurkat = 1)

OG_glycoprotein_up_down_median_upset <- OG_glycoprotein_total |> 
  left_join(OG_glycoprotein_up_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_down_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_median_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_up_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_down_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_median_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_up_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_down_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_median_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
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

OG_glycoprotein_up_down_median_upset_dataframe <- as.data.frame(OG_glycoprotein_up_down_median_upset)

upset(OG_glycoprotein_up_down_median_upset_dataframe, 
      sets = c("up_HepG2", "up_HEK293T", "up_Jurkat",
               "down_HepG2", "down_HEK293T", "down_Jurkat"
               ),
      point.size = 5, keep.order = TRUE,
      text.scale = c(1.5, 1.5, 1.5, 1.5, 2, 2))



