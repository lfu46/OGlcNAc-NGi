#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "showtext", "ggpubr", "rstatix", "ComplexHeatmap", "circlize", "eulerr")
lapply(packages_names, require, character.only = TRUE)

#import data
nucleolus_database <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\nucleolus\\nucleolus_database.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)

nucleolus_uniprot <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\nucleolus\\nucleolus_uniprot.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)

nucleolus_uniprot_database <- nucleolus_database |> 
  left_join(nucleolus_uniprot, by = join_by(Gene_Name == From)) |> 
  select(UniprotID = Entry, Location)

#generate database
#HEK293T
OG_glycoprotein_nucleolus_HEK293T <- OG_glycoprotein_RBP_HEK293T |> 
  left_join(nucleolus_uniprot_database, by = "UniprotID") |> 
  filter(!is.na(Location))

WP_protein_nucleolus_HEK293T <- WP_protein_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_nucleolus_HEK293T$UniprotID) |> 
  left_join(nucleolus_uniprot_database, by = "UniprotID") |> 
  mutate(cell = "HEK293T", group = "WP_RBP")

#HepG2
OG_glycoprotein_nucleolus_HepG2 <- OG_glycoprotein_RBP_HepG2 |> 
  left_join(nucleolus_uniprot_database, by = "UniprotID") |> 
  filter(!is.na(Location))

WP_protein_nucleolus_HepG2 <- WP_protein_Top_tb_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_nucleolus_HepG2$UniprotID) |> 
  left_join(nucleolus_uniprot_database, by = "UniprotID") |> 
  mutate(cell = "HepG2", group = "WP_RBP")

#Jurkat
OG_glycoprotein_nucleolus_Jurkat <- OG_glycoprotein_RBP_Jurkat |> 
  left_join(nucleolus_uniprot_database, by = "UniprotID") |> 
  filter(!is.na(Location))

WP_protein_nucleolus_Jurkat <- WP_protein_Top_tb_Jurkat |> 
  filter(UniprotID %in% OG_glycoprotein_nucleolus_Jurkat$UniprotID) |> 
  left_join(nucleolus_uniprot_database, by = "UniprotID") |> 
  mutate(cell = "Jurkat", group = "WP_RBP")

#combine
OG_glycoprotein_nucleolus_combined <- bind_rows(
  OG_glycoprotein_nucleolus_HEK293T,
  WP_protein_nucleolus_HEK293T,
  OG_glycoprotein_nucleolus_HepG2,
  WP_protein_nucleolus_HepG2,
  OG_glycoprotein_nucleolus_Jurkat,
  WP_protein_nucleolus_Jurkat
)

#wilcox test
OG_glycoprotein_nucleolus_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ group, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

