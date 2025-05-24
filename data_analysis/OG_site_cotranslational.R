#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "rstatix", "ggpubr")
lapply(packages_names, require, character.only = TRUE)

#import data
cotranslational_site_data <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\cotranslational_site\\ac2c04779_si_005.xlsx",
  skip = 1,
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(Protein_site = paste(UniProt.ID, Site.1.Position, sep = '_'))

#generate data frame
#HEK293T
OG_glycopeptide_Top_tb_HEK293T_singlesite |> 
  mutate(Protein_site = paste(UniprotID, site_position, sep = '_')) |> 
  mutate(
    category = ifelse(Protein_site %in% cotranslational_site_data$Protein_site, 'co', 'not-co')
  ) |> 
    group_by(category) |> 
  get_summary_stats(logFC, type = 'mean')
