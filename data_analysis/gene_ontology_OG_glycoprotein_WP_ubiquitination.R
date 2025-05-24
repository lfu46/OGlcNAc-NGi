#import packages
packages_names <- c("tidyverse", "showtext", "readxl", "ggpubr", "rstatix")
lapply(packages_names, require, character.only = TRUE)

#import data
#HepG2
gene_ontology_OG_glycoprotein_BP_ubiquitination_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_BP_ubiquitination.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "HepG2"
) |> select(UniprotID = UNIPROT_ACCESSION)

#HEK293T
gene_ontology_OG_glycoprotein_BP_ubiquitination_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_BP_ubiquitination.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "HEK293T"
) |> select(UniprotID = UNIPROT_ACCESSION)

#Jurkat
gene_ontology_OG_glycoprotein_BP_ubiquitination_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_BP_ubiquitination.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "Jurkat"
) |> select(UniprotID = UNIPROT_ACCESSION)

#generate data frame
#HepG2
gene_ontology_OG_glycoprotein_WP_BP_ubiquitination_logFC_HepG2 <- gene_ontology_OG_glycoprotein_BP_ubiquitination_HepG2 |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  select(UniprotID, OG = logFC) |> 
  left_join(WP_protein_Top_tb_HepG2, by = "UniprotID") |> 
  filter(!is.na(logFC)) |> 
  select(UniprotID, OG, WP = logFC) |> 
  pivot_longer(cols = OG:WP, names_to = 'Group', values_to = 'logFC') |> 
  mutate(Cell = "HepG2")

OG_WP_ks_test_HepG2 <- ks.test(
  gene_ontology_OG_glycoprotein_WP_BP_ubiquitination_logFC_HepG2 |> filter(Group == "OG") |> pull(logFC),
  gene_ontology_OG_glycoprotein_WP_BP_ubiquitination_logFC_HepG2 |> filter(Group == "WP") |> pull(logFC)
)

#HEK293T
gene_ontology_OG_glycoprotein_WP_BP_ubiquitination_logFC_HEK293T <- gene_ontology_OG_glycoprotein_BP_ubiquitination_HEK293T |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  select(UniprotID, OG = logFC) |> 
  left_join(WP_protein_Top_tb_HEK293T, by = "UniprotID") |> 
  filter(!is.na(logFC)) |> 
  select(UniprotID, OG, WP = logFC) |> 
  pivot_longer(cols = OG:WP, names_to = 'Group', values_to = 'logFC') |> 
  mutate(Cell = "HEK293T")

OG_WP_ks_test_HEK293T <- ks.test(
  gene_ontology_OG_glycoprotein_WP_BP_ubiquitination_logFC_HEK293T |> filter(Group == "OG") |> pull(logFC),
  gene_ontology_OG_glycoprotein_WP_BP_ubiquitination_logFC_HEK293T |> filter(Group == "WP") |> pull(logFC)
)

#Jurkat
gene_ontology_OG_glycoprotein_WP_BP_ubiquitination_logFC_Jurkat <- gene_ontology_OG_glycoprotein_BP_ubiquitination_Jurkat |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  select(UniprotID, OG = logFC) |> 
  left_join(WP_protein_Top_tb_Jurkat, by = "UniprotID") |> 
  filter(!is.na(logFC)) |> 
  select(UniprotID, OG, WP = logFC) |> 
  pivot_longer(cols = OG:WP, names_to = 'Group', values_to = 'logFC') |> 
  mutate(Cell = "Jurkat")

OG_WP_ks_test_Jurkat <- ks.test(
  gene_ontology_OG_glycoprotein_WP_BP_ubiquitination_logFC_Jurkat |> filter(Group == "OG") |> pull(logFC),
  gene_ontology_OG_glycoprotein_WP_BP_ubiquitination_logFC_Jurkat |> filter(Group == "WP") |> pull(logFC)
)

#combine
gene_ontology_OG_glycoprotein_WP_BP_ubiquitination_logFC_total <- bind_rows(
  gene_ontology_OG_glycoprotein_WP_BP_ubiquitination_logFC_HepG2,
  gene_ontology_OG_glycoprotein_WP_BP_ubiquitination_logFC_HEK293T,
  gene_ontology_OG_glycoprotein_WP_BP_ubiquitination_logFC_Jurkat
)

#t test
gene_ontology_OG_glycoprotein_WP_BP_ubiquitination_logFC_t_test <- gene_ontology_OG_glycoprotein_WP_BP_ubiquitination_logFC_total |> 
  group_by(Cell) |> 
  wilcox_test(logFC ~ Group, p.adjust.method = "BH") |> 
  add_significance("p")

