#import packages
packages_names <- c("tidyverse", "rstatix", "ggpubr", "showtext", "readxl", "ggbeeswarm")
lapply(packages_names, require, character.only = TRUE)

#protein folding
#import data
gene_ontology_OG_glycoprotein_BP_proteinfolding_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_BP_proteinfolding.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "HepG2"
) |> select(UniprotID = UNIPROT_ACCESSION) |> mutate(cell = "HepG2")

gene_ontology_OG_glycoprotein_BP_proteinfolding_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_BP_proteinfolding.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "HEK293T"
) |> select(UniprotID = UNIPROT_ACCESSION) |> mutate(cell = "HEK293T")

gene_ontology_OG_glycoprotein_BP_proteinfolding_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_BP_proteinfolding.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "Jurkat"
) |> select(UniprotID = UNIPROT_ACCESSION) |> mutate(cell = "Jurkat")

#generate data frame
gene_ontology_OG_glycoprotein_BP_proteinfolding_logFC_HepG2 <- gene_ontology_OG_glycoprotein_BP_proteinfolding_HepG2 |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  select(UniprotID, cell, logFC)

gene_ontology_OG_glycoprotein_BP_proteinfolding_logFC_HEK293T <- gene_ontology_OG_glycoprotein_BP_proteinfolding_HEK293T |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  select(UniprotID, cell, logFC)

gene_ontology_OG_glycoprotein_BP_proteinfolding_logFC_Jurkat <- gene_ontology_OG_glycoprotein_BP_proteinfolding_Jurkat |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  select(UniprotID, cell, logFC)

gene_ontology_OG_glycoprotein_BP_proteinfolding_total <- bind_rows(
  gene_ontology_OG_glycoprotein_BP_proteinfolding_logFC_HepG2,
  gene_ontology_OG_glycoprotein_BP_proteinfolding_logFC_HEK293T,
  gene_ontology_OG_glycoprotein_BP_proteinfolding_logFC_Jurkat
)

#wilcox test
gene_ontology_OG_glycoprotein_BP_proteinfolding_wilcox_test <- gene_ontology_OG_glycoprotein_BP_proteinfolding_total |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")





