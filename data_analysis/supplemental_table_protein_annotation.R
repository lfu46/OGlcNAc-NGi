#import packages
packages_names <- c("tidyverse", "readxl", "writexl")
lapply(packages_names, require, character.only = TRUE)

#protein annotation
#HEK293T OG
protein_annotation_OG_HEK293T <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data\\NOCT_HEK293T_08102024\\OG_sequest\\ronghuwulab_1734144424.csv",
  col_names = TRUE,
  name_repair = 'universal'
) |> separate(reference, sep = '\\|', into = c('sp', 'UniProt.Accession', 'Protein.Name')) |> 
  filter(UniProt.Accession %in% OG_glycoprotein_Top_tb_HEK293T$UniprotID) |> 
  select(UniProt.Accession, Annotation) |> 
  arrange(UniProt.Accession) |> filter(!duplicated(UniProt.Accession))

write_xlsx(protein_annotation_OG_HEK293T, path = paste0(file_path, 'protein_annotation_OG_HEK293T.xlsx'))

#HepG2 OG
protein_annotation_OG_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data\\NOCT_HepG2_07042024\\OG_sequest\\ronghuwulab_1734175894.csv",
  col_names = TRUE,
  name_repair = 'universal'
) |> separate(reference, sep = '\\|', into = c('sp', 'UniProt.Accession', 'Protein.Name')) |> 
  filter(UniProt.Accession %in% OG_glycoprotein_Top_tb_HepG2$UniprotID) |> 
  select(UniProt.Accession, Annotation) |> 
  arrange(UniProt.Accession) |> filter(!duplicated(UniProt.Accession))

write_xlsx(protein_annotation_OG_HepG2, path = paste0(file_path, 'protein_annotation_OG_HepG2.xlsx'))

#Jurkat OG
protein_annotation_OG_Jurkat <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data\\NOCT_Jurkat_08102024\\OG_sequest\\ronghuwulab_1734176339.csv",
  col_names = TRUE,
  name_repair = 'universal'
) |> separate(reference, sep = '\\|', into = c('sp', 'UniProt.Accession', 'Protein.Name')) |> 
  filter(UniProt.Accession %in% OG_glycoprotein_Top_tb_Jurkat$UniprotID) |> 
  select(UniProt.Accession, Annotation) |> 
  arrange(UniProt.Accession) |> filter(!duplicated(UniProt.Accession))

write_xlsx(protein_annotation_OG_Jurkat, path = paste0(file_path, 'protein_annotation_OG_Jurkat.xlsx'))

#HEK293T WP
protein_annotation_WP_HEK293T <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data\\NOCT_HEK293T_08102024\\WP_sequest\\ronghuwulab_1734180390.csv",
  col_names = TRUE,
  name_repair = 'universal'
) |> separate(reference, sep = '\\|', into = c('sp', 'UniProt.Accession', 'Protein.Name')) |> 
  filter(UniProt.Accession %in% WP_protein_Top_tb_HEK293T$UniprotID) |> 
  select(UniProt.Accession, Annotation) |> 
  arrange(UniProt.Accession) |> filter(!duplicated(UniProt.Accession))

write_xlsx(protein_annotation_WP_HEK293T, path = paste0(file_path, 'protein_annotation_WP_HEK293T.xlsx'))

#HepG2 WP
protein_annotation_WP_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data\\NOCT_HepG2_07042024\\WP_sequest\\ronghuwulab_1734180762.csv",
  col_names = TRUE,
  name_repair = 'universal'
) |> separate(reference, sep = '\\|', into = c('sp', 'UniProt.Accession', 'Protein.Name')) |> 
  filter(UniProt.Accession %in% WP_protein_Top_tb_HepG2$UniprotID) |> 
  select(UniProt.Accession, Annotation) |> 
  arrange(UniProt.Accession) |> filter(!duplicated(UniProt.Accession))

write_xlsx(protein_annotation_WP_HepG2, path = paste0(file_path, 'protein_annotation_WP_HepG2.xlsx'))

#Jurkat WP
protein_annotation_WP_Jurkat <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data\\NOCT_Jurkat_08102024\\WP_sequest\\ronghuwulab_1734180202.csv",
  col_names = TRUE,
  name_repair = 'universal'
) |> separate(reference, sep = '\\|', into = c('sp', 'UniProt.Accession', 'Protein.Name')) |> 
  filter(UniProt.Accession %in% WP_protein_Top_tb_Jurkat$UniprotID) |> 
  select(UniProt.Accession, Annotation) |> 
  arrange(UniProt.Accession) |> filter(!duplicated(UniProt.Accession))

write_xlsx(protein_annotation_WP_Jurkat, path = paste0(file_path, 'protein_annotation_WP_Jurkat.xlsx'))

#HEK293T OG site
protein_annotation_OG_HEK293T <- read_xlsx(
  'E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\protein_annotation\\protein_annotation_OG_HEK293T.xlsx'
)

site_protein_annotation_OG_HEK293T <- OG_psm_2_raw_HEK293T |> select(UniprotID, localized) |> 
  filter(localized != 0) |> 
  select(!localized) |> 
  left_join(protein_annotation_OG_HEK293T, by = join_by('UniprotID' == 'UniProt.Accession'))

write_xlsx(site_protein_annotation_OG_HEK293T, path = paste0(file_path, "site_protein_annotation_OG_HEK293T.xlsx"))

#HepG2 OG site
protein_annotation_OG_HepG2 <- read_xlsx(
  'E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\protein_annotation\\protein_annotation_OG_HepG2.xlsx'
)

site_protein_annotation_OG_HepG2 <- OG_psm_2_raw_HepG2 |> select(UniprotID, localized) |> 
  filter(localized != 0) |> 
  select(!localized) |> 
  left_join(protein_annotation_OG_HepG2, by = join_by('UniprotID' == 'UniProt.Accession'))

write_xlsx(site_protein_annotation_OG_HepG2, path = paste0(file_path, "site_protein_annotation_OG_HepG2.xlsx"))

#Jurkat OG site
protein_annotation_OG_Jurkat <- read_xlsx(
  'E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\protein_annotation\\protein_annotation_OG_Jurkat.xlsx'
)

site_protein_annotation_OG_Jurkat <- OG_psm_2_raw_Jurkat |> select(UniprotID, localized) |> 
  filter(localized != 0) |> 
  select(!localized) |> 
  left_join(protein_annotation_OG_Jurkat, by = join_by('UniprotID' == 'UniProt.Accession'))

write_xlsx(site_protein_annotation_OG_Jurkat, path = paste0(file_path, "site_protein_annotation_OG_Jurkat.xlsx"))
