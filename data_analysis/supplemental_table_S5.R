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

#abundance normalization
#HEK293T
OG_glycoprotein_protein_supple5_HEK293T_1 <- OG_glycoprotein_raw_sl_tmm_HEK293T |> 
  dplyr::select(UniprotID, ends_with("sl_tmm")) |> 
  mutate(
    log2_ratio1_OG_glycoprotein = log2(Tuni_1_sl_tmm / Ctrl_4_sl_tmm),
    log2_ratio2_OG_glycoprotein = log2(Tuni_2_sl_tmm / Ctrl_5_sl_tmm),
    log2_ratio3_OG_glycoprotein = log2(Tuni_3_sl_tmm / Ctrl_6_sl_tmm)
  ) |> 
  dplyr::select(UniprotID, starts_with("log2")) |> 
  left_join(WP_protein_raw_sl_tmm_HEK293T, by = "UniprotID") |> 
  dplyr::select(UniprotID, log2_ratio1_OG_glycoprotein:log2_ratio3_OG_glycoprotein, ends_with("sl_tmm")) |> 
  mutate(
    log2_ratio1_WP = log2(Tuni_1_sl_tmm / Ctrl_4_sl_tmm),
    log2_ratio2_WP = log2(Tuni_2_sl_tmm / Ctrl_5_sl_tmm),
    log2_ratio3_WP = log2(Tuni_3_sl_tmm / Ctrl_6_sl_tmm)
  ) |> 
  dplyr::select(UniprotID, starts_with("log2")) |> 
  left_join(OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T, by = "UniprotID")

OG_glycoprotein_protein_supple5_HEK293T_Gene <- tibble(
  bitr(OG_glycoprotein_protein_supple5_HEK293T_1$UniprotID, fromType = 'UNIPROT', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')
) |> filter(!duplicated(UNIPROT))

OG_glycoprotein_protein_supple5_HEK293T_2 <- OG_glycoprotein_protein_supple5_HEK293T_1 |> 
  left_join(OG_glycoprotein_protein_supple5_HEK293T_Gene, by = join_by(UniprotID == UNIPROT)) |> 
  dplyr::select(
    UniProt_Accession = UniprotID, 
    Gene = SYMBOL,
    #'log2(Tuni/Ctrl)OG 1' = log2_ratio1_OG_glycoprotein, 
    #'log2(Tuni/Ctrl)OG 2' = log2_ratio2_OG_glycoprotein, 
    #'log2(Tuni/Ctrl)OG 3' = log2_ratio3_OG_glycoprotein,
    #'log2(Tuni/Ctrl)WP 1' = log2_ratio1_WP,
    #'log2(Tuni/Ctrl)WP 2' = log2_ratio2_WP,
    #'log2(Tuni/Ctrl)WP 3' = log2_ratio3_WP,
    'Δlog2(Tuni/Ctrl)' = logFC,
    'adjusted P value' = adj.P.Val
  ) |> mutate(
    category = case_when(
      `Δlog2(Tuni/Ctrl)` > 0.5 & `adjusted P value` < 0.05 ~ 'up regulated',
      `Δlog2(Tuni/Ctrl)` < -0.5 & `adjusted P value` < 0.05 ~ 'down regulated',
      .default = 'not regulated'
    )
  )

#HepG2
OG_glycoprotein_protein_supple5_HepG2_1 <- OG_glycoprotein_raw_sl_tmm_HepG2 |> 
  dplyr::select(UniprotID, ends_with("sl_tmm")) |> 
  mutate(
    log2_ratio1_OG_glycoprotein = log2(Tuni_1_sl_tmm / Ctrl_4_sl_tmm),
    log2_ratio2_OG_glycoprotein = log2(Tuni_2_sl_tmm / Ctrl_5_sl_tmm),
    log2_ratio3_OG_glycoprotein = log2(Tuni_3_sl_tmm / Ctrl_6_sl_tmm)
  ) |> 
  dplyr::select(UniprotID, starts_with("log2")) |> 
  left_join(WP_protein_raw_sl_tmm_HepG2, by = "UniprotID") |> 
  dplyr::select(UniprotID, log2_ratio1_OG_glycoprotein:log2_ratio3_OG_glycoprotein, ends_with("sl_tmm")) |> 
  mutate(
    log2_ratio1_WP = log2(Tuni_1_sl_tmm / Ctrl_4_sl_tmm),
    log2_ratio2_WP = log2(Tuni_2_sl_tmm / Ctrl_5_sl_tmm),
    log2_ratio3_WP = log2(Tuni_3_sl_tmm / Ctrl_6_sl_tmm)
  ) |> 
  dplyr::select(UniprotID, starts_with("log2")) |> 
  left_join(OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2, by = "UniprotID")

OG_glycoprotein_protein_supple5_HepG2_Gene <- tibble(
  bitr(OG_glycoprotein_protein_supple5_HepG2_1$UniprotID, fromType = 'UNIPROT', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')
) |> filter(!duplicated(UNIPROT))

OG_glycoprotein_protein_supple5_HepG2_2 <- OG_glycoprotein_protein_supple5_HepG2_1 |> 
  left_join(OG_glycoprotein_protein_supple5_HepG2_Gene, by = join_by(UniprotID == UNIPROT)) |> 
  dplyr::select(
    UniProt_Accession = UniprotID, 
    Gene = SYMBOL,
    #'log2(Tuni/Ctrl)OG 1' = log2_ratio1_OG_glycoprotein, 
    #'log2(Tuni/Ctrl)OG 2' = log2_ratio2_OG_glycoprotein, 
    #'log2(Tuni/Ctrl)OG 3' = log2_ratio3_OG_glycoprotein,
    #'log2(Tuni/Ctrl)WP 1' = log2_ratio1_WP,
    #'log2(Tuni/Ctrl)WP 2' = log2_ratio2_WP,
    #'log2(Tuni/Ctrl)WP 3' = log2_ratio3_WP,
    'Δlog2(Tuni/Ctrl)' = logFC,
    'adjusted P value' = adj.P.Val
  ) |> mutate(
    category = case_when(
      `Δlog2(Tuni/Ctrl)` > 0.5 & `adjusted P value` < 0.05 ~ 'up regulated',
      `Δlog2(Tuni/Ctrl)` < -0.5 & `adjusted P value` < 0.05 ~ 'down regulated',
      .default = 'not regulated'
    )
  )

#Jurkat
OG_glycoprotein_protein_supple5_Jurkat_1 <- OG_glycoprotein_raw_sl_tmm_Jurkat |> 
  dplyr::select(UniprotID, ends_with("sl_tmm")) |> 
  mutate(
    log2_ratio1_OG_glycoprotein = log2(Tuni_1_sl_tmm / Ctrl_4_sl_tmm),
    log2_ratio2_OG_glycoprotein = log2(Tuni_2_sl_tmm / Ctrl_5_sl_tmm),
    log2_ratio3_OG_glycoprotein = log2(Tuni_3_sl_tmm / Ctrl_6_sl_tmm)
  ) |> 
  dplyr::select(UniprotID, starts_with("log2")) |> 
  left_join(WP_protein_raw_sl_tmm_Jurkat, by = "UniprotID") |> 
  dplyr::select(UniprotID, log2_ratio1_OG_glycoprotein:log2_ratio3_OG_glycoprotein, ends_with("sl_tmm")) |> 
  mutate(
    log2_ratio1_WP = log2(Tuni_1_sl_tmm / Ctrl_4_sl_tmm),
    log2_ratio2_WP = log2(Tuni_2_sl_tmm / Ctrl_5_sl_tmm),
    log2_ratio3_WP = log2(Tuni_3_sl_tmm / Ctrl_6_sl_tmm)
  ) |> 
  dplyr::select(UniprotID, starts_with("log2")) |> 
  left_join(OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat, by = "UniprotID")

OG_glycoprotein_protein_supple5_Jurkat_Gene <- tibble(
  bitr(OG_glycoprotein_protein_supple5_Jurkat_1$UniprotID, fromType = 'UNIPROT', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')
) |> filter(!duplicated(UNIPROT))

OG_glycoprotein_protein_supple5_Jurkat_2 <- OG_glycoprotein_protein_supple5_Jurkat_1 |> 
  left_join(OG_glycoprotein_protein_supple5_Jurkat_Gene, by = join_by(UniprotID == UNIPROT)) |> 
  dplyr::select(
    UniProt_Accession = UniprotID, 
    Gene = SYMBOL,
    #'log2(Tuni/Ctrl)OG 1' = log2_ratio1_OG_glycoprotein, 
    #'log2(Tuni/Ctrl)OG 2' = log2_ratio2_OG_glycoprotein, 
    #'log2(Tuni/Ctrl)OG 3' = log2_ratio3_OG_glycoprotein,
    #'log2(Tuni/Ctrl)WP 1' = log2_ratio1_WP,
    #'log2(Tuni/Ctrl)WP 2' = log2_ratio2_WP,
    #'log2(Tuni/Ctrl)WP 3' = log2_ratio3_WP,
    'Δlog2(Tuni/Ctrl)' = logFC,
    'adjusted P value' = adj.P.Val
  ) |> mutate(
    category = case_when(
      `Δlog2(Tuni/Ctrl)` > 0.5 & `adjusted P value` < 0.05 ~ 'up regulated',
      `Δlog2(Tuni/Ctrl)` < -0.5 & `adjusted P value` < 0.05 ~ 'down regulated',
      .default = 'not regulated'
    )
  )

#list
supple5_list <- list(
  OG_glycoprotein_protein_supple5_HEK293T_2,
  OG_glycoprotein_protein_supple5_HepG2_2,
  OG_glycoprotein_protein_supple5_Jurkat_2
)

write_xlsx(supple5_list, path = paste0(file_path, "supplemental_table_S5.xlsx"))

site_protein_annotation_list <- list(
  OG_glycoprotein_protein_supple5_HEK293T_2 |> 
    dplyr::select(UniProt_Accession) |> 
    left_join(protein_annotation_OG_HEK293T, by = join_by('UniProt_Accession' == 'UniProt.Accession')),
  OG_glycoprotein_protein_supple5_HepG2_2 |> 
    dplyr::select(UniProt_Accession) |> 
    left_join(protein_annotation_OG_HepG2, by = join_by('UniProt_Accession' == 'UniProt.Accession')),
  OG_glycoprotein_protein_supple5_Jurkat_2 |> 
    dplyr::select(UniProt_Accession) |> 
    left_join(protein_annotation_OG_Jurkat, by = join_by('UniProt_Accession' == 'UniProt.Accession'))
)

write_xlsx(site_protein_annotation_list, path = paste0(file_path, "site_protein_annotation_list.xlsx"))
