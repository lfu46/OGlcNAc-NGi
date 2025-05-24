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

#O-GlcNAcylation site abundance change
#HEK293T
OG_site_supple3_HEK293T_1 <- OG_psm_2_raw_HEK293T |> 
  dplyr::select(
    Index, UniProt_Accession = UniprotID, 'ModScore Peptide' = ModScore.Peptide, XCorr, PPM, 
    'Theo m/z' = Theo.m.z, 'Obs m/z' = Obs.m.z, 'Precursor Charge' = charge,
    'Site 1' = Site1, 'Site 1 position' = Site.1.Position, 'Site 1 ModScore' = Site.1.Score,
    'Site 2' = Site2, 'Site 2 position' = Site.2.Position, 'Site 2 ModScore' = Site.2.Score,
    localized
  ) |> 
  filter(localized != 0) |> 
  left_join(OG_glycopeptide_raw_sl_tmm_HEK293T, by = c("Index")) |> 
  dplyr::select(
    Index:'Site 2 ModScore', 'Tuni 1 (126)' = Tuni_1_sl_tmm, 'Tuni 2 (127)' = Tuni_2_sl_tmm, 
    'Tuni 3 (128)' = Tuni_3_sl_tmm, 'Ctrl 1 (129)' = Ctrl_4_sl_tmm, 
    'Ctrl 2 (130)' = Ctrl_5_sl_tmm, 'Ctrl 3 (131)' = Ctrl_6_sl_tmm
  ) |> 
  left_join(OG_glycopeptide_Top_tb_HEK293T, by = c("Index")) |> 
  dplyr::select(
    Index:'Ctrl 3 (131)', 'Avg(log2(Tuni/Ctrl))' = logFC, 'adjusted P value' = adj.P.Val
  )

OG_site_supple3_HEK293T_Gene <- tibble(
  bitr(OG_site_supple3_HEK293T_1$UniProt_Accession, fromType = 'UNIPROT', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')
) |> filter(!duplicated(UNIPROT))
  
OG_site_supple3_HEK293T_2 <- OG_site_supple3_HEK293T_1 |> 
  left_join(OG_site_supple3_HEK293T_Gene, by = join_by(UniProt_Accession == UNIPROT)) |> 
  mutate('Site Assignment' = 'Unambiguous') |> 
  dplyr::select(
    Index, UniProt_Accession, Gene = SYMBOL, 'ModScore Peptide', 'XCorr', 'PPM', 
    'Theo m/z', 'Obs m/z', 'Precursor Charge',
    'Site 1', 'Site 1 position', 'Site 1 ModScore',
    'Site 2', 'Site 2 position', 'Site 2 ModScore',
    'Site Assignment', 'Tuni 1 (126)':'Ctrl 3 (131)', 
    'Avg(log2(Tuni/Ctrl))', 'adjusted P value'
  ) |> mutate(
    category = case_when(
      `Avg(log2(Tuni/Ctrl))` > 0.5 & `adjusted P value` < 0.05 ~ 'up regulated',
      `Avg(log2(Tuni/Ctrl))` < -0.5 & `adjusted P value` < 0.05 ~ 'down regulated',
      .default = 'not regulated'
    )
  )

#HepG2
OG_site_supple3_HepG2_1 <- OG_psm_2_raw_HepG2 |> 
  dplyr::select(
    Index, UniProt_Accession = UniprotID, 'ModScore Peptide' = ModScore.Peptide, XCorr, PPM, 
    'Theo m/z' = Theo.m.z, 'Obs m/z' = Obs.m.z, 'Precursor Charge' = charge,
    'Site 1' = Site1, 'Site 1 position' = Site.1.Position, 'Site 1 ModScore' = Site.1.Score,
    'Site 2' = Site2, 'Site 2 position' = Site.2.Position, 'Site 2 ModScore' = Site.2.Score,
    localized
  ) |> 
  filter(localized != 0) |> 
  left_join(OG_glycopeptide_raw_sl_tmm_HepG2, by = c("Index")) |> 
  dplyr::select(
    Index:'Site 2 ModScore', 'Tuni 1 (126)' = Tuni_1_sl_tmm, 'Tuni 2 (127)' = Tuni_2_sl_tmm, 
    'Tuni 3 (128)' = Tuni_3_sl_tmm, 'Ctrl 1 (129)' = Ctrl_4_sl_tmm, 
    'Ctrl 2 (130)' = Ctrl_5_sl_tmm, 'Ctrl 3 (131)' = Ctrl_6_sl_tmm
  ) |> 
  left_join(OG_glycopeptide_Top_tb_HepG2, by = c("Index")) |> 
  dplyr::select(
    Index:'Ctrl 3 (131)', 'Avg(log2(Tuni/Ctrl))' = logFC, 'adjusted P value' = adj.P.Val
  )

OG_site_supple3_HepG2_Gene <- tibble(
  bitr(OG_site_supple3_HepG2_1$UniProt_Accession, fromType = 'UNIPROT', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')
) |> filter(!duplicated(UNIPROT))

OG_site_supple3_HepG2_2 <- OG_site_supple3_HepG2_1 |> 
  left_join(OG_site_supple3_HepG2_Gene, by = join_by(UniProt_Accession == UNIPROT)) |> 
  mutate('Site Assignment' = 'Unambiguous') |> 
  dplyr::select(
    Index, UniProt_Accession, Gene = SYMBOL, 'ModScore Peptide', 'XCorr', 'PPM', 
    'Theo m/z', 'Obs m/z', 'Precursor Charge',
    'Site 1', 'Site 1 position', 'Site 1 ModScore',
    'Site 2', 'Site 2 position', 'Site 2 ModScore', 'Site Assignment',
    'Tuni 1 (126)':'Ctrl 3 (131)', 
    'Avg(log2(Tuni/Ctrl))', 'adjusted P value'
  ) |> mutate(
    category = case_when(
      `Avg(log2(Tuni/Ctrl))` > 0.5 & `adjusted P value` < 0.05 ~ 'up regulated',
      `Avg(log2(Tuni/Ctrl))` < -0.5 & `adjusted P value` < 0.05 ~ 'down regulated',
      .default = 'not regulated'
    )
  )

#Jurkat
OG_site_supple3_Jurkat_1 <- OG_psm_2_raw_Jurkat |> 
  dplyr::select(
    Index, UniProt_Accession = UniprotID, 'ModScore Peptide' = ModScore.Peptide, XCorr, PPM, 
    'Theo m/z' = Theo.m.z, 'Obs m/z' = Obs.m.z, 'Precursor Charge' = charge,
    'Site 1' = Site1, 'Site 1 position' = Site.1.Position, 'Site 1 ModScore' = Site.1.Score,
    'Site 2' = Site2, 'Site 2 position' = Site.2.Position, 'Site 2 ModScore' = Site.2.Score,
    localized
  ) |> 
  filter(localized != 0) |> 
  left_join(OG_glycopeptide_raw_sl_tmm_Jurkat, by = c("Index")) |> 
  dplyr::select(
    Index:'Site 2 ModScore', 'Tuni 1 (126)' = Tuni_1_sl_tmm, 'Tuni 2 (127)' = Tuni_2_sl_tmm, 
    'Tuni 3 (128)' = Tuni_3_sl_tmm, 'Ctrl 1 (129)' = Ctrl_4_sl_tmm, 
    'Ctrl 2 (130)' = Ctrl_5_sl_tmm, 'Ctrl 3 (131)' = Ctrl_6_sl_tmm
  ) |> 
  left_join(OG_glycopeptide_Top_tb_Jurkat, by = c("Index")) |> 
  dplyr::select(
    Index:'Ctrl 3 (131)', 'Avg(log2(Tuni/Ctrl))' = logFC, 'adjusted P value' = adj.P.Val
  ) 

OG_site_supple3_Jurkat_Gene <- tibble(
  bitr(OG_site_supple3_Jurkat_1$UniProt_Accession, fromType = 'UNIPROT', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')
) |> filter(!duplicated(UNIPROT))

OG_site_supple3_Jurkat_2 <- OG_site_supple3_Jurkat_1 |> 
  left_join(OG_site_supple3_Jurkat_Gene, by = join_by(UniProt_Accession == UNIPROT)) |> 
  mutate('Site Assignment' = 'Unambiguous') |> 
  dplyr::select(
    Index, UniProt_Accession, Gene = SYMBOL, 'ModScore Peptide', 'XCorr', 'PPM', 
    'Theo m/z', 'Obs m/z', 'Precursor Charge',
    'Site 1', 'Site 1 position', 'Site 1 ModScore',
    'Site 2', 'Site 2 position', 'Site 2 ModScore', 'Site Assignment',
    'Tuni 1 (126)':'Ctrl 3 (131)', 
    'Avg(log2(Tuni/Ctrl))', 'adjusted P value'
  ) |> mutate(
    category = case_when(
      `Avg(log2(Tuni/Ctrl))` > 0.5 & `adjusted P value` < 0.05 ~ 'up regulated',
      `Avg(log2(Tuni/Ctrl))` < -0.5 & `adjusted P value` < 0.05 ~ 'down regulated',
      .default = 'not regulated'
    )
  )

#list
supple3_list <- list(
  OG_site_supple3_HEK293T_2 |> 
    dplyr::select(Index, UniProt_Accession, Gene, `ModScore Peptide`:`Site 2 ModScore`),
  OG_site_supple3_HEK293T_2 |> 
    dplyr::select(Index, UniProt_Accession, Gene, `Site 1`, `Tuni 1 (126)`:`category`) |> 
    distinct(),
  OG_site_supple3_HepG2_2 |> 
    dplyr::select(Index, UniProt_Accession, Gene, `ModScore Peptide`:`Site 2 ModScore`),
  OG_site_supple3_HepG2_2 |> 
    dplyr::select(Index, UniProt_Accession, Gene, `Site 1`, `Tuni 1 (126)`:`category`) |> 
    distinct(),
  OG_site_supple3_Jurkat_2 |> 
    dplyr::select(Index, UniProt_Accession, Gene, `ModScore Peptide`:`Site 2 ModScore`),
  OG_site_supple3_Jurkat_2 |> 
    dplyr::select(Index, UniProt_Accession, Gene, `Site 1`, `Tuni 1 (126)`:`category`) |> 
    distinct()
)

write_xlsx(supple3_list, path = paste0(file_path, "supplemental_table_S3_R.xlsx"))

#OG site annotation
site_annotation_list <- list(
  OG_site_supple3_HEK293T_2 |> 
    dplyr::select(Index, UniProt_Accession) |> 
    distinct() |> 
    left_join(protein_annotation_OG_HEK293T, by = join_by('UniProt_Accession' == 'UniProt.Accession')),
  OG_site_supple3_HepG2_2 |> 
    dplyr::select(Index, UniProt_Accession) |> 
    distinct() |> 
    left_join(protein_annotation_OG_HepG2, by = join_by('UniProt_Accession' == 'UniProt.Accession')),
  OG_site_supple3_Jurkat_2 |> 
    dplyr::select(Index, UniProt_Accession) |> 
    distinct() |> 
    left_join(protein_annotation_OG_Jurkat, by = join_by('UniProt_Accession' == 'UniProt.Accession'))
)

write_xlsx(site_annotation_list, paste0(file_path, 'site_annotation_list.xlsx'))


