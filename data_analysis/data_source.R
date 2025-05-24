#import packages
packages_names <- c("tidyverse", "writexl", "readxl")
lapply(packages_names, require, character.only = TRUE)

#specify file path
file_path <- "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\"

#set color for plot
Color_1 = "#B5D086"
Color_2 = "#7FB2D4"
Color_3 = "#f8766d"
Color_4 = "#00bfc4"
Color_5 = "#FDE725"
Color_6 = "#440154"
Color_7 = "orange"
Color_8 = "#35b779"
Color_9 = "#31004a"

#import data source
#raw psm data from Sequest
OG_psm_raw_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_psm_raw_HepG2.xlsx"
)
OG_psm_raw_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_psm_raw_HEK293T.xlsx"
)
OG_psm_raw_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_psm_raw_Jurkat.xlsx"
)
OG_psm_2_raw_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_psm_2_raw_HepG2.xlsx"
)
OG_psm_2_raw_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_psm_2_raw_HEK293T.xlsx"
)
OG_psm_2_raw_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_psm_2_raw_Jurkat.xlsx"
)
WP_psm_raw_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\WP_psm_raw_HepG2.xlsx"
)
WP_psm_raw_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\WP_psm_raw_HEK293T.xlsx"
)
WP_psm_raw_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\WP_psm_raw_Jurkat.xlsx"
)

#OG glycopeptide data
OG_glycopeptide_raw_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycopeptide_raw_HepG2.xlsx"
)
OG_glycopeptide_raw_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycopeptide_raw_HEK293T.xlsx"
)
OG_glycopeptide_raw_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycopeptide_raw_Jurkat.xlsx"
)

#OG glycoprotein data
OG_glycoprotein_raw_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_raw_HepG2.xlsx"
)
OG_glycoprotein_raw_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_raw_HEK293T.xlsx"
)
OG_glycoprotein_raw_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_raw_Jurkat.xlsx"
)

#WP protein data
WP_protein_raw_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\WP_protein_raw_HepG2.xlsx"
)
WP_protein_raw_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\WP_protein_raw_HEK293T.xlsx"
)
WP_protein_raw_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\WP_protein_raw_Jurkat.xlsx"
)

#normalized OG glycopeptide data
OG_glycopeptide_raw_sl_tmm_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_normalization_OG_glycopeptide_HepG2\\OG_glycopeptide_raw_sl_tmm_HepG2.xlsx"
)
OG_glycopeptide_raw_sl_tmm_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_normalization_OG_glycopeptide_HEK293T\\OG_glycopeptide_raw_sl_tmm_HEK293T.xlsx"
)
OG_glycopeptide_raw_sl_tmm_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_normalization_OG_glycopeptide_Jurkat\\OG_glycopeptide_raw_sl_tmm_Jurkat.xlsx"
)

#normalized OG glycoprotein data
OG_glycoprotein_raw_sl_tmm_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_normalization_OG_glycoprotein_HepG2\\OG_glycoprotein_raw_sl_tmm_HepG2.xlsx"
)
OG_glycoprotein_raw_sl_tmm_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_normalization_OG_glycoprotein_HEK293T\\OG_glycoprotein_raw_sl_tmm_HEK293T.xlsx"
)
OG_glycoprotein_raw_sl_tmm_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_normalization_OG_glycoprotein_Jurkat\\OG_glycoprotein_raw_sl_tmm_Jurkat.xlsx"
)

#differerntial analysis OG glycopeptide data
OG_glycopeptide_Top_tb_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycopeptide_Top_tb_HepG2.xlsx"
)
OG_glycopeptide_Top_tb_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycopeptide_Top_tb_HEK293T.xlsx"
)
OG_glycopeptide_Top_tb_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycopeptide_Top_tb_Jurkat.xlsx"
)

#differerntial analysis OG glycoprotein data
OG_glycoprotein_Top_tb_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_Top_tb_HepG2.xlsx"
)
OG_glycoprotein_Top_tb_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_Top_tb_HEK293T.xlsx"
)
OG_glycoprotein_Top_tb_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_Top_tb_Jurkat.xlsx"
)

#normalized WP protein data
WP_protein_raw_sl_tmm_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_normalization_WP_protein_HepG2\\WP_protein_raw_sl_tmm_HepG2.xlsx"
)
WP_protein_raw_sl_tmm_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_normalization_WP_protein_HEK293T\\WP_protein_raw_sl_tmm_HEK293T.xlsx"
)
WP_protein_raw_sl_tmm_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_normalization_WP_protein_Jurkat\\WP_protein_raw_sl_tmm_Jurkat.xlsx"
)

#differential analysis WP protein data
WP_protein_Top_tb_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\WP_protein_Top_tb_HepG2.xlsx"
)
WP_protein_Top_tb_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\WP_protein_Top_tb_HEK293T.xlsx"
)
WP_protein_Top_tb_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\WP_protein_Top_tb_Jurkat.xlsx"
)

#OG total glycoprotein
OG_glycoprotein_total <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_total.xlsx"
)

#WP total protein
WP_protein_total <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\WP_protein_total.xlsx"
)

#OG glycoprotein cell specific/only
OG_glycoprotein_only_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_only_HepG2.xlsx"
)
OG_glycoprotein_only_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_only_HEK293T.xlsx"
)
OG_glycoprotein_only_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_only_Jurkat.xlsx"
)

#OG glycoprotein identified in all three cell lines
OG_glycoprotein_common_total <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_common_total.xlsx"
)

#OG glycoprotein up, down or middle regulated total
OG_glycoprotein_up_total <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_up_total.xlsx"
)
OG_glycoprotein_down_total <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_down_total.xlsx"
)
OG_glycoprotein_middle <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_middle.xlsx"
)

#OG glycoprotein up, down or median regulated combined
OG_glycoprotein_up_combined <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_up_combined.xlsx"
)
OG_glycoprotein_down_combined <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_down_combined.xlsx"
)
OG_glycoprotein_median_combined

#OG glycoprotein up in specific cell
OG_glycoprotein_up_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_up_HepG2.xlsx"
)
OG_glycoprotein_up_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_up_HEK293T.xlsx"
)
OG_glycoprotein_up_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_up_Jurkat.xlsx"
)

#OG glycoprotein down in specific cell
OG_glycoprotein_down_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_down_HepG2.xlsx"
)
OG_glycoprotein_down_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_down_HEK293T.xlsx"
)
OG_glycoprotein_down_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_down_Jurkat.xlsx"
)

#OG glycoprotein median in specific cell
OG_glycoprotein_median_HepG2

OG_glycoprotein_median_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_median_HEK293T.xlsx"
)
OG_glycoprotein_median_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_median_Jurkat.xlsx"
)

#single site
OG_glycopeptide_Top_tb_HepG2_singlesite <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycopeptide_Top_tb_HepG2_singlesite.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)
OG_glycopeptide_Top_tb_HEK293T_singlesite <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycopeptide_Top_tb_HEK293T_singlesite.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)
OG_glycopeptide_Top_tb_Jurkat_singlesite <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycopeptide_Top_tb_Jurkat_singlesite.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)

#OG site total
OG_site_total <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_site_total.xlsx"
)

#OG site logFC total
OG_site_logFC_total <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_site_logFC_total.xlsx"
)

#OG glycoprotein RBP
OG_glycoprotein_RBP_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_RBP_HEK293T.xlsx"
)
OG_glycoprotein_RBP_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_RBP_HepG2.xlsx"
)
OG_glycoprotein_RBP_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_RBP_Jurkat.xlsx"
)

#differential analysis OG site protein
OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T.xlsx"
)
OG_site_protein_OGlcNAc_effect_Top_tb_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_site_protein_OGlcNAc_effect_Top_tb_HepG2.xlsx"
)
OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat.xlsx"
)

#differential analysis OG glycoprotein protein
OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T.xlsx"
)
OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2.xlsx"
)
OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat.xlsx"
)

#OG glycoprotein generally up- or down-regulated 
OG_glycoprotein_generally_up <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_generally_up_across_three_cell_lines.xlsx"
)
OG_glycoprotein_generally_down <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycoprotein_generally_down_across_three_cell_lines.xlsx"
)

#OG site 31mers
OG_site_Mer31_regulated_combined <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_site_Mer31_regulated_combined.xlsx"
)
