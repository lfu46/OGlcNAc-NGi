#import packages
packages_names <- c("tidyverse", "writexl", "readxl", "ComplexHeatmap", "circlize", 
                    "ggpubr", "rstatix", "showtext")
lapply(packages_names, require, character.only = TRUE)

#pathway
pathway_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\exosome_signaling_pathway.xlsx",
  sheet = 'HEK293T',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION, category) |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  mutate(cell = "HEK293T")

pathway_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\exosome_signaling_pathway.xlsx",
  sheet = 'HepG2',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION, category) |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  mutate(cell = "HepG2")

pathway_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\exosome_signaling_pathway.xlsx",
  sheet = 'Jurkat',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION, category) |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  mutate(cell = "Jurkat")

#combine
pathway_combined <- bind_rows(
  pathway_HEK293T,
  pathway_HepG2,
  pathway_Jurkat
)

#wilcox test
pathway_combined |> 
  filter(category == "AMPK signaling pathway") |>  
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")
