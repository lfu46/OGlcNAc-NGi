#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "showtext", "ggpubr", "rstatix", "ComplexHeatmap", "circlize", "eulerr")
lapply(packages_names, require, character.only = TRUE)

#import data
stress_granule_gene <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\stress_granule_core\\stress_gruanle_core.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)

stress_granule_uniprot <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\stress_granule_core\\stress_gruanle_core_uniprot.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)

SG_core_list <- stress_granule_gene |> filter(Type == "SG-core") |> pull(Gene.Name)

stress_granule_gene_adj <- stress_granule_gene |> 
  mutate(Type = ifelse(Gene.Name %in% SG_core_list, "SG-core", Type)) |> 
  distinct()

stress_granule_gene_uniprot <- stress_granule_gene_adj |> 
  left_join(stress_granule_uniprot, by = join_by(Gene.Name == From), ) |> 
  filter(!is.na(Entry)) |> 
  select(UniprotID = Entry, Type)

#generate data frame
#HEK293T
OG_glycoprotein_SG_type_HEK293T <- OG_glycoprotein_RBP_HEK293T |> 
  left_join(stress_granule_gene_uniprot, by = "UniprotID") |> 
  filter(!is.na(Type))

WP_protein_SG_type_HEK293T <- WP_protein_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_SG_type_HEK293T$UniprotID) |> 
  left_join(stress_granule_gene_uniprot, by = "UniprotID") |> 
  mutate(cell = "HEK293T", group = "WP_RBP")

#HepG2
OG_glycoprotein_SG_type_HepG2 <- OG_glycoprotein_RBP_HepG2 |> 
  left_join(stress_granule_gene_uniprot, by = "UniprotID") |> 
  filter(!is.na(Type))

WP_protein_SG_type_HepG2 <- WP_protein_Top_tb_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_SG_type_HepG2$UniprotID) |> 
  left_join(stress_granule_gene_uniprot, by = "UniprotID") |> 
  mutate(cell = "HepG2", group = "WP_RBP")

#Jurkat
OG_glycoprotein_SG_type_Jurkat <- OG_glycoprotein_RBP_Jurkat |> 
  left_join(stress_granule_gene_uniprot, by = "UniprotID") |> 
  filter(!is.na(Type))

WP_protein_SG_type_Jurkat <- WP_protein_Top_tb_Jurkat |> 
  filter(UniprotID %in% OG_glycoprotein_SG_type_Jurkat$UniprotID) |> 
  left_join(stress_granule_gene_uniprot, by = "UniprotID") |> 
  mutate(cell = "Jurkat", group = "WP_RBP")

#combine
OG_WP_SG_type_Jurkat <- bind_rows(
  OG_glycoprotein_SG_type_Jurkat,
  WP_protein_SG_type_Jurkat
)

OG_WP_SG_type_Jurkat |> 
  group_by(Type) |> 
  wilcox_test(logFC ~ group, p.adjust.method = "BH")


