#import packages
packages_names <- c("tidyverse", "rstatix", "ggpubr", "showtext", "readxl", "ggbeeswarm")
lapply(packages_names, require, character.only = TRUE)

#import data
#HEK293T
gene_ontology_OG_site_protein_nucleocytoplasmic_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_site_protein_nucleocytoplasmic_transport.xlsx",
  sheet = 'HEK293T',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_site_protein_nucleocytoplasmic_HEK293T <- OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  filter(UniprotID %in% gene_ontology_OG_site_protein_nucleocytoplasmic_HEK293T$UNIPROT_ACCESSION) |> 
  select(Index, logFC) |> mutate(cell = "HEK293T")

#HepG2
gene_ontology_OG_site_protein_nucleocytoplasmic_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_site_protein_nucleocytoplasmic_transport.xlsx",
  sheet = 'HepG2',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_site_protein_nucleocytoplasmic_HepG2 <- OG_site_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  filter(UniprotID %in% gene_ontology_OG_site_protein_nucleocytoplasmic_HepG2$UNIPROT_ACCESSION) |> 
  select(Index, logFC) |> mutate(cell = "HepG2")

#Jurkat
gene_ontology_OG_site_protein_nucleocytoplasmic_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_site_protein_nucleocytoplasmic_transport.xlsx",
  sheet = 'Jurkat',
  col_names = TRUE,
  .name_repair = "universal"
)

OG_site_protein_nucleocytoplasmic_Jurkat <- OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  filter(UniprotID %in% gene_ontology_OG_site_protein_nucleocytoplasmic_Jurkat$UNIPROT_ACCESSION) |> 
  select(Index, logFC) |> mutate(cell = "Jurkat")

#combine
OG_site_protein_nucleocytoplasmic_combined <- bind_rows(
  OG_site_protein_nucleocytoplasmic_HEK293T,
  OG_site_protein_nucleocytoplasmic_HepG2,
  OG_site_protein_nucleocytoplasmic_Jurkat
)

#wilcox test
OG_site_protein_nucleoplasmic_wilcox_test <- OG_site_protein_nucleocytoplasmic_combined |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

