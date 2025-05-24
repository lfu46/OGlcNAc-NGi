#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "rstatix")
lapply(packages_names, require, character.only = TRUE)

#import data 
APR_N_site <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\APR_N_glycosylation\\adk8173_table_S2.xlsx",
  skip = 2, 
  col_names = TRUE,
  .name_repair = 'universal'
) |> 
  mutate(
    start = Position.in.protein - 15,
    end = Position.in.protein + 15
  ) |> 
  select(UniprotID = UniProt.ID, start, end)

#generate data frame
df <- OG_glycopeptide_Top_tb_Jurkat_singlesite |> 
  filter(UniprotID %in% APR_N_site$UniprotID) |> 
  left_join(APR_N_site, by = 'UniprotID', relationship = 'many-to-many') |> 
  mutate(APR = ifelse(site_position >= start & site_position <= end, 'APR', 'non-APR'))

write_xlsx(df, path = paste0(file_path, "df.xlsx"))
