#import packages
packages_names <- c("UpSetR", "ggupset", "grDevices")
lapply(packages_names, require, character.only = TRUE)

#generate total glycoprotein
OG_glycoprotein_total <- bind_rows(
  OG_glycoprotein_Top_tb_HepG2 |> distinct(UniprotID),
  OG_glycoprotein_Top_tb_HEK293T |> distinct(UniprotID),
  OG_glycoprotein_Top_tb_Jurkat |> distinct(UniprotID)
) |> distinct()

OG_glycoprotein_upset_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> select(UniprotID) |> mutate(HepG2 = 1)

OG_glycoprotein_upset_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> select(UniprotID) |> mutate(HEK293T = 1)

OG_glycoprotein_upset_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> select(UniprotID) |> mutate(Jurkat = 1)

OG_glycoprotein_total_upset <- OG_glycoprotein_total |> 
  left_join(OG_glycoprotein_upset_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_upset_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_upset_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  mutate(
    HepG2 = ifelse(is.na(HepG2), 0, 1),
    HEK293T = ifelse(is.na(HEK293T), 0, 1),
    Jurkat = ifelse(is.na(Jurkat), 0, 1)
  )

OG_glycoprotein_total_upset_dataframe <- as.data.frame(OG_glycoprotein_total_upset)

upset(OG_glycoprotein_total_upset_dataframe, sets = c("HepG2", "HEK293T", "Jurkat"),
      main.bar.color = Color_2,
      sets.bar.color = Color_2,
      order.by = "freq", point.size = 5,
      text.scale = c(2, 1.5, 1.5, 1.5, 2, 3),
      mainbar.y.label = "# of quantified \nO-GlcNAcylated proteins")


