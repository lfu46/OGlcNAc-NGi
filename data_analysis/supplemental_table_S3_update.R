# import packages
library(tidyverse)
library(readxl)
library(writexl)

## improt previous supplementary table S5
# HEK293T
OG_site_supple3_HEK293T_2 <- read_xlsx(
  "/Volumes/Expansion/Glycosylation_Crosstalk/Data_Analysis_Sequest/supplemental_table_S5.xlsx", 
  sheet = 1,
  skip = 1,
  col_names = TRUE
)

# HepG2
OG_site_supple3_HepG2_2 <- read_xlsx(
  "/Volumes/Expansion/Glycosylation_Crosstalk/Data_Analysis_Sequest/supplemental_table_S5.xlsx", 
  sheet = 2,
  skip = 1,
  col_names = TRUE
)

# Jurkat
OG_site_supple3_Jurkat_2 <- read_xlsx(
  "/Volumes/Expansion/Glycosylation_Crosstalk/Data_Analysis_Sequest/supplemental_table_S5.xlsx", 
  sheet = 3,
  skip = 1,
  col_names = TRUE
)

## import trypticity and miss cleavage information
# HEK293T
OG_glycopeptide_trypticity_HEK293T <- read_csv(
  "/Volumes/Expansion/Glycosylation_Crosstalk/Data/NOCT_HEK293T_08102024/OG_sequest/ronghuwulab_1744379701.csv", 
  col_names = TRUE
) |> 
  select(
    'ModScore Peptide', 'Trypticity', 'MissedCleav'
  ) |> 
  distinct()

# HepG2
OG_glycopeptide_trypticity_HepG2 <- read_csv(
  "/Volumes/Expansion/Glycosylation_Crosstalk/Data/NOCT_HepG2_07042024/OG_sequest/ronghuwulab_1744379918.csv", 
  col_names = TRUE
) |> 
  select(
    'ModScore Peptide', 'Trypticity', 'MissedCleav'
  ) |> 
  distinct()

# Jurkat
OG_glycopeptide_trypticity_Jurkat <- read_csv(
  "/Volumes/Expansion/Glycosylation_Crosstalk/Data/NOCT_Jurkat_08102024/OG_sequest/ronghuwulab_1744380017.csv", 
  col_names = TRUE
) |> 
  select(
    'ModScore Peptide', 'Trypticity', 'MissedCleav'
  ) |> 
  distinct()

## combine trypticity information with supplementary table S5
# HEK293T
OG_glycopeptide_supple5_trypticity_HEK293T <- OG_site_supple3_HEK293T_2 |> 
  left_join(
    OG_glycopeptide_trypticity_HEK293T,
    by = c("Glycopeptide" = "ModScore Peptide")
  )

# HepG2
OG_glycopeptide_supple5_trypticity_HepG2 <- OG_site_supple3_HepG2_2 |> 
  left_join(
    OG_glycopeptide_trypticity_HepG2,
    by = c("Glycopeptide" = "ModScore Peptide")
  )

# Jurkat
OG_glycopeptide_supple5_trypticity_Jurkat <- OG_site_supple3_Jurkat_2 |> 
  left_join(
    OG_glycopeptide_trypticity_Jurkat,
    by = c("Glycopeptide" = "ModScore Peptide")
  )

# output supplementary table S5
OG_glycopeptide_supple5_trypticity_list <- list(
  OG_glycopeptide_supple5_trypticity_HEK293T,
  OG_glycopeptide_supple5_trypticity_HepG2,
  OG_glycopeptide_supple5_trypticity_Jurkat
)

write_xlsx(
  OG_glycopeptide_supple5_trypticity_list, 
  path = "/Volumes/Expansion/Glycosylation_Crosstalk/Data_Analysis_Sequest/supplemental_table_S5_update.xlsx"
)
