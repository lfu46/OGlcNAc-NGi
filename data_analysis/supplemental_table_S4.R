#import packages
packages_names <- c('tidyverse', 'readxl', 'writexl')
lapply(packages_names, require, character.only = TRUE)

#import data
OG_site_Mer31_regulated_combined_sequence_parameter <- read_xlsx(
  'E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_site_Mer31_regulated_combined_sequence_parameter.xlsx'
)
OG_glycopeptide_localized_singlesite_secondary_structure_combined <- read_xlsx(
  'E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycopeptide_localized_singlesite_secondary_structure_combined.xlsx'
)
OG_site_disorder_tb <- read_xlsx(
  'E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_site_disorder_tb.xlsx'
)
OG_glycopeptide_Top_tb_singlesite_combined_adj_2_NG <- read_xlsx(
  'E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycopeptide_Top_tb_singlesite_combined_adj_2_NG.xlsx'
)

#OG site 31mer sequence feature
#singlesite combined
OG_glycopeptide_Top_tb_singlesite_combined <- bind_rows(
  OG_glycopeptide_Top_tb_HEK293T_singlesite |> mutate(cell = 'HEK293T'),
  OG_glycopeptide_Top_tb_HepG2_singlesite |> mutate(cell = 'HepG2'),
  OG_glycopeptide_Top_tb_Jurkat_singlesite |> mutate(cell = 'Jurkat')
) |> dplyr::select(Index, UniprotID, combined_site, cell, 'Avg(log2(Tuni/Ctrl))' = logFC, 'adjusted P value' = adj.P.Val)

#generate data frame
OG_glycopeptide_Top_tb_singlesite_combined_sequence_feature <- OG_glycopeptide_Top_tb_singlesite_combined |> 
  left_join(OG_site_Mer31_regulated_combined_sequence_parameter, 
            by = c('Index', 'UniprotID', 'cell')) |> 
  select(Index:`adjusted P value`, Mer31, `Isoelectric point` = Isoelectric_point, `FCR`) |> 
  left_join(OG_glycopeptide_localized_singlesite_secondary_structure_combined, 
            by = c('Index', 'UniprotID', 'combined_site', 'cell')) |> 
  select(Index:`FCR`, `Exposed/Buried` = Class.Assignment) |> 
  left_join(OG_site_disorder_tb, 
            by = c('Index', 'UniprotID', 'combined_site', 'cell')) |> 
  select(Index:`Exposed/Buried`, `Ordered/Disordered` = Disorder) |> 
  left_join(OG_glycopeptide_Top_tb_singlesite_combined_adj_2_NG,
            by = c('Index', 'UniprotID', 'combined_site', 'cell')) |> 
  select(Index, UniProt_Accession = UniprotID, combined_site, cell:`Ordered/Disordered`, `NG site/no NG site` = cluster)

#list
supple4_list <- list(
  OG_glycopeptide_Top_tb_singlesite_combined_sequence_feature |> 
    filter(cell == 'HEK293T'),
  OG_glycopeptide_Top_tb_singlesite_combined_sequence_feature |> 
    filter(cell == 'HepG2'),
  OG_glycopeptide_Top_tb_singlesite_combined_sequence_feature |> 
    filter(cell == 'Jurkat')
)

write_xlsx(supple4_list, path = paste0(file_path, "supplemental_table_S4_R.xlsx"))

