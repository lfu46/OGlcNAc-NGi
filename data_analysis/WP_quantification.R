#import packages
packages_names <- c("tidyverse", "limma", "edgeR", "showtext")
lapply(packages_names, require, character.only = TRUE)

#protein level quantification
#HepG2
WP_protein_raw_HepG2 <- WP_peptide_raw_HepG2 |> 
  group_by(UniprotID) |> 
  summarize(
    Tuni_1 = sum(Sn.126),
    Tuni_2 = sum(Sn.127n),
    Tuni_3 = sum(Sn.128c),
    Ctrl_4 = sum(Sn.129n),
    Ctrl_5 = sum(Sn.130c),
    Ctrl_6 = sum(Sn.131)
  )

write_xlsx(WP_protein_raw_HepG2, path = paste0(file_path, "WP_protein_raw_HepG2.xlsx"))

#HEK293T
WP_protein_raw_HEK293T <- WP_peptide_raw_HEK293T |> 
  group_by(UniprotID) |> 
  summarize(
    Tuni_1 = sum(Sn.126),
    Tuni_2 = sum(Sn.127n),
    Tuni_3 = sum(Sn.128c),
    Ctrl_4 = sum(Sn.129n),
    Ctrl_5 = sum(Sn.130c),
    Ctrl_6 = sum(Sn.131)
  )

write_xlsx(WP_protein_raw_HEK293T, path = paste0(file_path, "WP_protein_raw_HEK293T.xlsx"))

#Jurkat
WP_protein_raw_Jurkat <- WP_peptide_raw_Jurkat |> 
  group_by(UniprotID) |> 
  summarize(
    Tuni_1 = sum(Sn.126),
    Tuni_2 = sum(Sn.127n),
    Tuni_3 = sum(Sn.128c),
    Ctrl_4 = sum(Sn.129n),
    Ctrl_5 = sum(Sn.130c),
    Ctrl_6 = sum(Sn.131)
  )

write_xlsx(WP_protein_raw_Jurkat, path = paste0(file_path, "WP_protein_raw_Jurkat.xlsx"))
