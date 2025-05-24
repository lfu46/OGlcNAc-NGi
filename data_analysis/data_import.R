#import packages
packages_names <- c("tidyverse", "writexl")
lapply(packages_names, require, character.only = TRUE)

#specify file path
file_path <- "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\"

#import OG result
OG_psm_raw_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data\\NOCT_HepG2_07042024\\OG_sequest\\ronghuwulab_1729609351.csv",
  col_names = TRUE,
  name_repair = "universal"
)

OG_psm_raw_HEK293T <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data\\NOCT_HEK293T_08102024\\OG_sequest\\ronghuwulab_1729615074.csv",
  col_names = TRUE,
  name_repair = "universal"
)

OG_psm_raw_Jurkat <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data\\NOCT_Jurkat_08102024\\OG_sequest\\ronghuwulab_1729615323.csv",
  col_names = TRUE,
  name_repair = "universal"
)

#generate data frame for HepG2
OG_psm_raw_HepG2 <- OG_psm_raw_HepG2 |> 
  select(ScanF, SrchID, SrchName, Trimmed.Peptide, PPM, XCorr, Theo.m.z, Obs.m.z, charge = z,
         Start.Position, End.Position, Gene.Symbol, Reference, 
         ModScore.Peptide, Site.1.Position:Site.6.Score, 
         Sn.126 = ..126.Sn,
         Sn.127n = ..127n.Sn,
         Sn.128c = ..128c.Sn,
         Sn.129n = ..129n.Sn,
         Sn.130c = ..130c.Sn,
         Sn.131 = ..131.Sn,
         Sum.Sn) |> 
  filter(XCorr > 1.2, PPM > -10, PPM < 10, Sn.126 > 5, Sn.127n > 5, Sn.128c > 5, Sn.129n > 5, Sn.130c > 5, Sn.131 > 5) |> 
  filter(!str_detect(ModScore.Peptide, "C@|C#@")) |> 
  separate(Reference, into = c("sp", "UniprotID", "Gene.Symbol"), sep = "\\|") |> 
  rowwise() |> 
  mutate(
    identified = sum(c_across(c('Site.1.Position', 'Site.2.Position', 'Site.3.Position', 'Site.4.Position', 'Site.5.Position', 'Site.6.Position')) != 0),
    localized = sum(c_across(c('Site.1.Score', 'Site.2.Score', 'Site.3.Score', 'Site.4.Score', 'Site.5.Score', 'Site.6.Score')) > 13)
    ) |> 
  ungroup() |> 
  mutate(Site = str_extract_all(ModScore.Peptide, ".(?=@)"))

max_columns_HepG2 <- max(lengths(OG_psm_raw_HepG2$Site))

site_columns_HepG2 <- as.data.frame(do.call(rbind, lapply(OG_psm_raw_HepG2$Site, function(x){
  length(x) <- max_columns_HepG2
  return(x)
})))

colnames(site_columns_HepG2) <- paste0("Site", seq_len(max_columns_HepG2))

OG_psm_raw_HepG2 <- bind_cols(OG_psm_raw_HepG2, site_columns_HepG2) |> 
  mutate(Site1 = paste0(Site1, Site.1.Position),
         Site2 = paste0(Site2, Site.2.Position),
         Site3 = paste0(Site3, Site.3.Position)
  ) |> 
  mutate(
    Site1 = ifelse(str_detect(Site1, "NA"), NA, Site1),
    Site2 = ifelse(str_detect(Site2, "NA"), NA, Site2),
    Site3 = ifelse(str_detect(Site3, "NA"), NA, Site3)
  )

OG_psm_raw_HepG2 <- OG_psm_raw_HepG2 |> 
  mutate(Site1 = ifelse(Site.1.Score > 13, Site1, NA),
         Site2 = ifelse(Site.2.Score > 13, Site2, NA),
         Site3 = ifelse(Site.3.Score > 13, Site3, NA)) |> 
  unite("combined_site", Site1:Site3, remove = FALSE, na.rm = TRUE, sep = "") |> 
  mutate(
    Index = ifelse(
      combined_site == "", 
      str_c(UniprotID, Start.Position, End.Position, identified, localized, sep = "_"),
      str_c(UniprotID, Start.Position, End.Position, identified, localized, combined_site, sep = "_")
    )
  ) |> 
  select(Index, ModScore.Peptide, Trimmed.Peptide, XCorr, PPM, Theo.m.z, Obs.m.z, charge, 
         UniprotID, Site.1.Position:Site.3.Score, Sn.126:Sn.131, Sum.Sn, Start.Position, End.Position, 
         identified, localized, combined_site, Site1:Site3, ScanF:SrchName)

write_xlsx(OG_psm_raw_HepG2, path = paste0(file_path, "OG_psm_raw_HepG2.xlsx"))

#generate data frame for HEK293T
OG_psm_raw_HEK293T <- OG_psm_raw_HEK293T |> 
  select(ScanF, SrchID, SrchName, Trimmed.Peptide, PPM, XCorr, Theo.m.z, Obs.m.z, charge = z,
         Start.Position, End.Position, Gene.Symbol, Reference, 
         ModScore.Peptide, Site.1.Position:Site.6.Score, 
         Sn.126 = ..126.Sn,
         Sn.127n = ..127n.Sn,
         Sn.128c = ..128c.Sn,
         Sn.129n = ..129n.Sn,
         Sn.130c = ..130c.Sn,
         Sn.131 = ..131.Sn,
         Sum.Sn) |> 
  filter(XCorr > 1.2, PPM > -10, PPM < 10, Sn.126 > 5, Sn.127n > 5, Sn.128c > 5, Sn.129n > 5, Sn.130c > 5, Sn.131 > 5) |> 
  filter(!str_detect(ModScore.Peptide, "C@|C#@")) |> 
  separate(Reference, into = c("sp", "UniprotID", "Gene.Symbol"), sep = "\\|") |> 
  rowwise() |> 
  mutate(
    identified = sum(c_across(c('Site.1.Position', 'Site.2.Position', 'Site.3.Position', 'Site.4.Position', 'Site.5.Position', 'Site.6.Position')) != 0),
    localized = sum(c_across(c('Site.1.Score', 'Site.2.Score', 'Site.3.Score', 'Site.4.Score', 'Site.5.Score', 'Site.6.Score')) > 13)
  ) |> 
  ungroup() |> 
  mutate(Site = str_extract_all(ModScore.Peptide, ".(?=@)"))

max_columns_HEK293T <- max(lengths(OG_psm_raw_HEK293T$Site))

site_columns_HEK293T <- as.data.frame(do.call(rbind, lapply(OG_psm_raw_HEK293T$Site, function(x){
  length(x) <- max_columns_HEK293T
  return(x)
})))

colnames(site_columns_HEK293T) <- paste0("Site", seq_len(max_columns_HEK293T))

OG_psm_raw_HEK293T <- bind_cols(OG_psm_raw_HEK293T, site_columns_HEK293T) |> 
  mutate(Site1 = paste0(Site1, Site.1.Position),
         Site2 = paste0(Site2, Site.2.Position),
         Site3 = paste0(Site3, Site.3.Position)
  ) |> 
  mutate(
    Site1 = ifelse(str_detect(Site1, "NA"), NA, Site1),
    Site2 = ifelse(str_detect(Site2, "NA"), NA, Site2),
    Site3 = ifelse(str_detect(Site3, "NA"), NA, Site3)
  )

OG_psm_raw_HEK293T <- OG_psm_raw_HEK293T |> 
  mutate(Site1 = ifelse(Site.1.Score > 13, Site1, NA),
         Site2 = ifelse(Site.2.Score > 13, Site2, NA),
         Site3 = ifelse(Site.3.Score > 13, Site3, NA)) |> 
  unite("combined_site", Site1:Site3, remove = FALSE, na.rm = TRUE, sep = "") |> 
  mutate(
    Index = ifelse(
    combined_site == "", 
    str_c(UniprotID, Start.Position, End.Position, identified, localized, sep = "_"),
    str_c(UniprotID, Start.Position, End.Position, identified, localized, combined_site, sep = "_")
  )
  ) |> 
  select(Index, ModScore.Peptide, Trimmed.Peptide, XCorr, PPM, Theo.m.z, Obs.m.z, charge,
         UniprotID, Site.1.Position:Site.3.Score, Sn.126:Sn.131, Sum.Sn, Start.Position, End.Position, 
         identified, localized, combined_site, Site1:Site3, ScanF:SrchName)

OG_psm_raw_HEK293T_localized <- OG_psm_raw_HEK293T |> 
  filter(localized != 0)

OG_psm_raw_HEK293T_non_localized <- OG_psm_raw_HEK293T |> 
  filter(localized == 0) |> 
  filter(!str_detect(Trimmed.Peptide, "C[^#]"))

OG_psm_raw_HEK293T <- bind_rows(OG_psm_raw_HEK293T_localized, OG_psm_raw_HEK293T_non_localized)

write_xlsx(OG_psm_raw_HEK293T, path = paste0(file_path, "OG_psm_raw_HEK293T.xlsx"))

#generate data frame for Jurkat
OG_psm_raw_Jurkat <- OG_psm_raw_Jurkat |> 
  select(ScanF, SrchID, SrchName, Trimmed.Peptide, PPM, XCorr, Theo.m.z, Obs.m.z, charge = z,
         Start.Position, End.Position, Gene.Symbol, Reference, 
         ModScore.Peptide, Site.1.Position:Site.6.Score, 
         Sn.126 = ..126.Sn,
         Sn.127n = ..127n.Sn,
         Sn.128c = ..128c.Sn,
         Sn.129n = ..129n.Sn,
         Sn.130c = ..130c.Sn,
         Sn.131 = ..131.Sn,
         Sum.Sn) |> 
  filter(XCorr > 1.2, PPM > -10, PPM < 10, Sn.126 > 5, Sn.127n > 5, Sn.128c > 5, Sn.129n > 5, Sn.130c > 5, Sn.131 > 5) |> 
  filter(!str_detect(ModScore.Peptide, "C@|C#@")) |> 
  separate(Reference, into = c("sp", "UniprotID", "Gene.Symbol"), sep = "\\|") |> 
  rowwise() |> 
  mutate(
    identified = sum(c_across(c('Site.1.Position', 'Site.2.Position', 'Site.3.Position', 'Site.4.Position', 'Site.5.Position', 'Site.6.Position')) != 0),
    localized = sum(c_across(c('Site.1.Score', 'Site.2.Score', 'Site.3.Score', 'Site.4.Score', 'Site.5.Score', 'Site.6.Score')) > 13)
  ) |> 
  ungroup() |> 
  mutate(Site = str_extract_all(ModScore.Peptide, ".(?=@)"))

max_columns_Jurkat <- max(lengths(OG_psm_raw_Jurkat$Site))

site_columns_Jurkat <- as.data.frame(do.call(rbind, lapply(OG_psm_raw_Jurkat$Site, function(x){
  length(x) <- max_columns_Jurkat
  return(x)
})))

colnames(site_columns_Jurkat) <- paste0("Site", seq_len(max_columns_Jurkat))

OG_psm_raw_Jurkat <- bind_cols(OG_psm_raw_Jurkat, site_columns_Jurkat) |> 
  mutate(Site1 = paste0(Site1, Site.1.Position),
         Site2 = paste0(Site2, Site.2.Position),
         Site3 = paste0(Site3, Site.3.Position)
  ) |> 
  mutate(
    Site1 = ifelse(str_detect(Site1, "NA"), NA, Site1),
    Site2 = ifelse(str_detect(Site2, "NA"), NA, Site2),
    Site3 = ifelse(str_detect(Site3, "NA"), NA, Site3)
  )

OG_psm_raw_Jurkat <- OG_psm_raw_Jurkat |> 
  mutate(Site1 = ifelse(Site.1.Score > 13, Site1, NA),
         Site2 = ifelse(Site.2.Score > 13, Site2, NA),
         Site3 = ifelse(Site.3.Score > 13, Site3, NA)) |> 
  unite("combined_site", Site1:Site3, remove = FALSE, na.rm = TRUE, sep = "") |> 
  mutate(
    Index = ifelse(
      combined_site == "", 
      str_c(UniprotID, Start.Position, End.Position, identified, localized, sep = "_"),
      str_c(UniprotID, Start.Position, End.Position, identified, localized, combined_site, sep = "_")
    )
  ) |> 
  select(Index, ModScore.Peptide, Trimmed.Peptide, XCorr, PPM, Theo.m.z, Obs.m.z, charge,
         UniprotID, Site.1.Position:Site.3.Score, Sn.126:Sn.131, Sum.Sn, Start.Position, End.Position, 
         identified, localized, combined_site, Site1:Site3, ScanF:SrchName)

OG_psm_raw_Jurkat_localized <- OG_psm_raw_Jurkat |> 
  filter(localized != 0)

OG_psm_raw_Jurkat_non_localized <- OG_psm_raw_Jurkat |> 
  filter(localized == 0) |> 
  filter(!str_detect(Trimmed.Peptide, "C[^#]"))

OG_psm_raw_Jurkat <- bind_rows(OG_psm_raw_Jurkat_localized, OG_psm_raw_Jurkat_non_localized)

write_xlsx(OG_psm_raw_Jurkat, path = paste0(file_path, "OG_psm_raw_Jurkat.xlsx"))

#import WP result
WP_psm_raw_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data\\NOCT_HepG2_07042024\\WP_sequest\\ronghuwulab_1725718785.csv",
  col_names = TRUE,
  name_repair = "universal"
)

WP_psm_raw_HEK293T <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data\\NOCT_HEK293T_08102024\\WP_sequest\\ronghuwulab_1725719351.csv",
  col_names = TRUE,
  name_repair = "universal"
)

WP_psm_raw_Jurkat <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data\\NOCT_Jurkat_08102024\\WP_sequest\\ronghuwulab_1725717952.csv",
  col_names = TRUE,
  name_repair = "universal"
)

#generate data frame for HepG2
WP_psm_raw_HepG2 <- WP_psm_raw_HepG2 |> 
  select(ScanF, SrchID, SrchName, psm, Trimmed.psm, PPM, XCorr, Start.Position, End.Position, Gene.Symbol, Reference, 
         Sn.126 = ..126.Sn,
         Sn.127n = ..127n.Sn,
         Sn.128c = ..128c.Sn,
         Sn.129n = ..129n.Sn,
         Sn.130c = ..130c.Sn,
         Sn.131 = ..131.Sn,
         Sum.Sn) |> 
  filter(XCorr > 1.2, PPM > -10, PPM < 10, Sn.126 > 5, Sn.127n > 5, Sn.128c > 5, Sn.129n > 5, Sn.130c > 5, Sn.131 > 5) |> 
  separate(Reference, into = c("sp", "UniprotID", "Gene.Symbol"), sep = "\\|") |> 
  select(psm, Trimmed.psm, XCorr, PPM, UniprotID, Sn.126:Sn.131, Sum.Sn, Start.Position, End.Position, ScanF:SrchName)

write_xlsx(WP_psm_raw_HepG2, path = paste0(file_path, "WP_psm_raw_HepG2.xlsx"))

#generate data frame for HEK293T
WP_psm_raw_HEK293T <- WP_psm_raw_HEK293T |> 
  select(ScanF, SrchID, SrchName, psm, Trimmed.psm, PPM, XCorr, Start.Position, End.Position, Gene.Symbol, Reference, 
         Sn.126 = ..126.Sn,
         Sn.127n = ..127n.Sn,
         Sn.128c = ..128c.Sn,
         Sn.129n = ..129n.Sn,
         Sn.130c = ..130c.Sn,
         Sn.131 = ..131.Sn,
         Sum.Sn) |> 
  filter(XCorr > 1.2, PPM > -10, PPM < 10, Sn.126 > 5, Sn.127n > 5, Sn.128c > 5, Sn.129n > 5, Sn.130c > 5, Sn.131 > 5) |> 
  separate(Reference, into = c("sp", "UniprotID", "Gene.Symbol"), sep = "\\|") |> 
  select(psm, Trimmed.psm, XCorr, PPM, UniprotID, Sn.126:Sn.131, Sum.Sn, Start.Position, End.Position, ScanF:SrchName)

write_xlsx(WP_psm_raw_HEK293T, path = paste0(file_path, "WP_psm_raw_HEK293T.xlsx"))

#generate data frame for Jurkat
WP_psm_raw_Jurkat <- WP_psm_raw_Jurkat |> 
  select(ScanF, SrchID, SrchName, psm, Trimmed.psm, PPM, XCorr, Start.Position, End.Position, Gene.Symbol, Reference, 
         Sn.126 = ..126.Sn,
         Sn.127n = ..127n.Sn,
         Sn.128c = ..128c.Sn,
         Sn.129n = ..129n.Sn,
         Sn.130c = ..130c.Sn,
         Sn.131 = ..131.Sn,
         Sum.Sn) |> 
  filter(XCorr > 1.2, PPM > -10, PPM < 10, Sn.126 > 5, Sn.127n > 5, Sn.128c > 5, Sn.129n > 5, Sn.130c > 5, Sn.131 > 5) |> 
  separate(Reference, into = c("sp", "UniprotID", "Gene.Symbol"), sep = "\\|") |> 
  select(psm, Trimmed.psm, XCorr, PPM, UniprotID, Sn.126:Sn.131, Sum.Sn, Start.Position, End.Position, ScanF:SrchName)

write_xlsx(WP_psm_raw_Jurkat, path = paste0(file_path, "WP_psm_raw_Jurkat.xlsx"))

