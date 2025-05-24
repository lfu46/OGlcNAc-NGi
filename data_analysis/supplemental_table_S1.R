#import packages
packages_names <- c("tidyverse", "writexl", "readxl", "clusterProfiler", "org.Hs.eg.db")
lapply(packages_names, require, character.only = TRUE)

#import data from uniprot
#HEK293T
OG_glycoprotein_Gene_Name_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Supplemental Data\\OG_glycoprotein_Gene_Name_HEK293T.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
) |> separate(Gene.Names, into = c('Gene Name'))

#HepG2
OG_glycoprotein_Gene_Name_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Supplemental Data\\OG_glycoprotein_Gene_Name_HepG2.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
) |> separate(Gene.Names, into = c('Gene Name'))

#Jurkat
OG_glycoprotein_Gene_Name_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Supplemental Data\\OG_glycoprotein_Gene_Name_Jurkat.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
) |> separate(Gene.Names, into = c('Gene Name'))

#import data from sequest
#HEK293T
OG_glycoprotein_coverage_HEK293T <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data\\NOCT_HEK293T_08102024\\OG_sequest\\ronghuwulab_1729615197.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> 
  filter(!str_detect(reference, '##')) |> 
  separate(reference, into = c('sp', 'UniProt_Accession', 'HUMAN'), sep = '\\|')
  
#HepG2
OG_glycoprotein_coverage_HepG2 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data\\NOCT_HepG2_07042024\\OG_sequest\\ronghuwulab_1729609528.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> 
  filter(!str_detect(reference, '##')) |> 
  separate(reference, into = c('sp', 'UniProt_Accession', 'HUMAN'), sep = '\\|')

#Jurkat
OG_glycoprotein_coverage_Jurkat <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data\\NOCT_Jurkat_08102024\\OG_sequest\\ronghuwulab_1729615436.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> 
  filter(!str_detect(reference, '##')) |> 
  separate(reference, into = c('sp', 'UniProt_Accession', 'HUMAN'), sep = '\\|')

#glycoprotein abundance change
#HEK293T
OG_glycoprotein_supple1_HEK293T_1 <- OG_glycoprotein_raw_sl_tmm_HEK293T |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  dplyr::select(UniProt_Accession = UniprotID, 
         'Tuni 1 (126)' = Tuni_1_sl_tmm, 'Tuni 2 (127)' = Tuni_2_sl_tmm, 'Tuni 3 (128)' = Tuni_3_sl_tmm,
         'Ctrl 1 (129)' = Ctrl_4_sl_tmm, 'Ctrl 2 (130)' = Ctrl_5_sl_tmm, 'Ctrl 3 (131)' = Ctrl_6_sl_tmm, 
         'Avg(log2(Tuni/Ctrl))' = logFC, 'adjusted P value' = adj.P.Val)

OG_glycoprotein_supple1_HEK293T_Gene <- tibble(bitr(OG_glycoprotein_supple1_HEK293T_1$UniProt_Accession, 
                                                    fromType = 'UNIPROT', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')) %>% 
  filter(!duplicated(UNIPROT))

OG_glycoprotein_supple1_HEK293T_2 <- OG_glycoprotein_supple1_HEK293T_1 |> 
  left_join(OG_glycoprotein_supple1_HEK293T_Gene, by = join_by('UniProt_Accession' == 'UNIPROT')) |> 
  left_join(OG_glycoprotein_coverage_HEK293T, by = 'UniProt_Accession') |> 
  dplyr::select(
    UniProt_Accession, Gene = SYMBOL, 'Total Peptide' = 'Total', 'Unique Peptide' = 'Unique',
    Length, Coverage, 'Coverage Percentage' = Coverage..,
    'Tuni 1 (126)':'Ctrl 3 (131)', 'Avg(log2(Tuni/Ctrl))', 'adjusted P value'
  ) |> mutate(
    category = case_when(
      `Avg(log2(Tuni/Ctrl))` > 0.5 & `adjusted P value` < 0.05 ~ 'up regulated',
      `Avg(log2(Tuni/Ctrl))` < -0.5 & `adjusted P value` < 0.05 ~ 'down regulated',
      .default = 'not regulated'
    )
  )

#HepG2
OG_glycoprotein_supple1_HepG2_1 <- OG_glycoprotein_raw_sl_tmm_HepG2 |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  dplyr::select(UniProt_Accession = UniprotID, 
         'Tuni 1 (126)' = Tuni_1_sl_tmm, 'Tuni 2 (127)' = Tuni_2_sl_tmm, 'Tuni 3 (128)' = Tuni_3_sl_tmm,
         'Ctrl 1 (129)' = Ctrl_4_sl_tmm, 'Ctrl 2 (130)' = Ctrl_5_sl_tmm, 'Ctrl 3 (131)' = Ctrl_6_sl_tmm, 
         'Avg(log2(Tuni/Ctrl))' = logFC, 'adjusted P value' = adj.P.Val)

OG_glycoprotein_supple1_HepG2_Gene <- tibble(bitr(OG_glycoprotein_supple1_HepG2_1$UniProt_Accession, 
                                                    fromType = 'UNIPROT', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')) %>% 
  filter(!duplicated(UNIPROT))

OG_glycoprotein_supple1_HepG2_2 <- OG_glycoprotein_supple1_HepG2_1 |> 
  left_join(OG_glycoprotein_supple1_HepG2_Gene, by = join_by('UniProt_Accession' == 'UNIPROT')) |> 
  left_join(OG_glycoprotein_coverage_HepG2, by = 'UniProt_Accession') |> 
  dplyr::select(
    UniProt_Accession, Gene = SYMBOL, 'Total Peptide' = 'Total', 'Unique Peptide' = 'Unique',
    Length, Coverage, 'Coverage Percentage' = Coverage..,
    'Tuni 1 (126)':'Ctrl 3 (131)', 'Avg(log2(Tuni/Ctrl))', 'adjusted P value'
  ) |> mutate(
    category = case_when(
      `Avg(log2(Tuni/Ctrl))` > 0.5 & `adjusted P value` < 0.05 ~ 'up regulated',
      `Avg(log2(Tuni/Ctrl))` < -0.5 & `adjusted P value` < 0.05 ~ 'down regulated',
      .default = 'not regulated'
    )
  )

#Jurkat
OG_glycoprotein_supple1_Jurkat_1 <- OG_glycoprotein_raw_sl_tmm_Jurkat |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  dplyr::select(UniProt_Accession = UniprotID, 
         'Tuni 1 (126)' = Tuni_1_sl_tmm, 'Tuni 2 (127)' = Tuni_2_sl_tmm, 'Tuni 3 (128)' = Tuni_3_sl_tmm,
         'Ctrl 1 (129)' = Ctrl_4_sl_tmm, 'Ctrl 2 (130)' = Ctrl_5_sl_tmm, 'Ctrl 3 (131)' = Ctrl_6_sl_tmm, 
         'Avg(log2(Tuni/Ctrl))' = logFC, 'adjusted P value' = adj.P.Val)

OG_glycoprotein_supple1_Jurkat_Gene <- tibble(bitr(OG_glycoprotein_supple1_Jurkat_1$UniProt_Accession, 
                                                  fromType = 'UNIPROT', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')) %>% 
  filter(!duplicated(UNIPROT))

OG_glycoprotein_supple1_Jurkat_2 <- OG_glycoprotein_supple1_Jurkat_1 |> 
  left_join(OG_glycoprotein_supple1_Jurkat_Gene,  by = join_by('UniProt_Accession' == 'UNIPROT')) |> 
  left_join(OG_glycoprotein_coverage_Jurkat, by = 'UniProt_Accession') |> 
  dplyr::select(
    UniProt_Accession, Gene = SYMBOL, 'Total Peptide' = 'Total', 'Unique Peptide' = 'Unique',
    Length, Coverage, 'Coverage Percentage' = Coverage..,
    'Tuni 1 (126)':'Ctrl 3 (131)', 'Avg(log2(Tuni/Ctrl))', 'adjusted P value'
  ) |> mutate(
    category = case_when(
      `Avg(log2(Tuni/Ctrl))` > 0.5 & `adjusted P value` < 0.05 ~ 'up regulated',
      `Avg(log2(Tuni/Ctrl))` < -0.5 & `adjusted P value` < 0.05 ~ 'down regulated',
      .default = 'not regulated'
    )
  )

#list
supple1_list <- list(
  OG_glycoprotein_supple1_HEK293T_2 |> 
    dplyr::select(UniProt_Accession, Gene, 
           'Tuni 1 (126)':'Ctrl 3 (131)'),
  OG_glycoprotein_supple1_HEK293T_2 |> 
    dplyr::select(UniProt_Accession, Gene, 
           'Total Peptide':'Unique Peptide', Length, Coverage, 'Coverage Percentage',
           'Avg(log2(Tuni/Ctrl))', 'adjusted P value', 'category'),
  OG_glycoprotein_supple1_HepG2_2 |> 
    dplyr::select(UniProt_Accession, Gene, 
           'Tuni 1 (126)':'Ctrl 3 (131)'),
  OG_glycoprotein_supple1_HepG2_2 |> 
    dplyr::select(UniProt_Accession, Gene, 
           'Total Peptide':'Unique Peptide', Length, Coverage, 'Coverage Percentage',
           'Avg(log2(Tuni/Ctrl))', 'adjusted P value', 'category'),
  OG_glycoprotein_supple1_Jurkat_2 |> 
    dplyr::select(UniProt_Accession, Gene, 
           'Tuni 1 (126)':'Ctrl 3 (131)'),
  OG_glycoprotein_supple1_Jurkat_2 |> 
    dplyr::select(UniProt_Accession, Gene, 
           'Total Peptide':'Unique Peptide', Length, Coverage, 'Coverage Percentage',
           'Avg(log2(Tuni/Ctrl))', 'adjusted P value', 'category')
)

write_xlsx(supple1_list, path = paste0(file_path, "supplemental_table_S1.xlsx"))
