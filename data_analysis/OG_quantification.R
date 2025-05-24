#import packages
packages_names <- c("tidyverse", "limma", "edgeR", "showtext")
lapply(packages_names, require, character.only = TRUE)

#glycopeptide level quantification
#HepG2
OG_glycopeptide_raw_HepG2 <- OG_peptide_raw_HepG2 |> 
  group_by(Index, UniprotID, combined_site, Start.Position, End.Position, identified, localized) |> 
  summarize(
    Tuni_1 = sum(Sn.126),
    Tuni_2 = sum(Sn.127n),
    Tuni_3 = sum(Sn.128c),
    Ctrl_4 = sum(Sn.129n),
    Ctrl_5 = sum(Sn.130c),
    Ctrl_6 = sum(Sn.131)
  ) |> 
  ungroup()

write_xlsx(OG_glycopeptide_raw_HepG2, path = paste0(file_path, "OG_glycopeptide_raw_HepG2.xlsx"))

#HEK293T
OG_glycopeptide_raw_HEK293T <- OG_peptide_raw_HEK293T |> 
  group_by(Index, UniprotID, combined_site, Start.Position, End.Position, identified, localized) |> 
  summarize(
    Tuni_1 = sum(Sn.126),
    Tuni_2 = sum(Sn.127n),
    Tuni_3 = sum(Sn.128c),
    Ctrl_4 = sum(Sn.129n),
    Ctrl_5 = sum(Sn.130c),
    Ctrl_6 = sum(Sn.131)
  ) |> 
  ungroup()

write_xlsx(OG_glycopeptide_raw_HEK293T, path = paste0(file_path, "OG_glycopeptide_raw_HEK293T.xlsx"))

#Jurkat
OG_glycopeptide_raw_Jurkat <- OG_peptide_raw_Jurkat |> 
  group_by(Index, UniprotID, combined_site, Start.Position, End.Position, identified, localized) |> 
  summarize(
    Tuni_1 = sum(Sn.126),
    Tuni_2 = sum(Sn.127n),
    Tuni_3 = sum(Sn.128c),
    Ctrl_4 = sum(Sn.129n),
    Ctrl_5 = sum(Sn.130c),
    Ctrl_6 = sum(Sn.131)
  ) |> 
  ungroup()

write_xlsx(OG_glycopeptide_raw_Jurkat, path = paste0(file_path, "OG_glycopeptide_raw_Jurkat.xlsx"))

#glycoprotein level quantification
#HepG2
OG_glycoprotein_raw_HepG2 <- OG_glycopeptide_raw_HepG2 |> 
  group_by(UniprotID) |> 
  summarize(
    Tuni_1 = sum(Tuni_1),
    Tuni_2 = sum(Tuni_2),
    Tuni_3 = sum(Tuni_3),
    Ctrl_4 = sum(Ctrl_4),
    Ctrl_5 = sum(Ctrl_5),
    Ctrl_6 = sum(Ctrl_6)
  )

write_xlsx(OG_glycoprotein_raw_HepG2, path = paste0(file_path, "OG_glycoprotein_raw_HepG2.xlsx"))

#HEK293T
OG_glycoprotein_raw_HEK293T <- OG_glycopeptide_raw_HEK293T |> 
  group_by(UniprotID) |> 
  summarize(
    Tuni_1 = sum(Tuni_1),
    Tuni_2 = sum(Tuni_2),
    Tuni_3 = sum(Tuni_3),
    Ctrl_4 = sum(Ctrl_4),
    Ctrl_5 = sum(Ctrl_5),
    Ctrl_6 = sum(Ctrl_6)
  )

write_xlsx(OG_glycoprotein_raw_HEK293T, path = paste0(file_path, "OG_glycoprotein_raw_HEK293T.xlsx"))

#Jurkat
OG_glycoprotein_raw_Jurkat <- OG_glycopeptide_raw_Jurkat |> 
  group_by(UniprotID) |> 
  summarize(
    Tuni_1 = sum(Tuni_1),
    Tuni_2 = sum(Tuni_2),
    Tuni_3 = sum(Tuni_3),
    Ctrl_4 = sum(Ctrl_4),
    Ctrl_5 = sum(Ctrl_5),
    Ctrl_6 = sum(Ctrl_6)
  )

write_xlsx(OG_glycoprotein_raw_Jurkat, path = paste0(file_path, "OG_glycoprotein_raw_Jurkat.xlsx"))

