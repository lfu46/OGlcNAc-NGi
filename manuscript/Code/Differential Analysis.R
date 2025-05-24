#import packages
packages_names <- c("tidyverse", "limma", "showtext", "writexl")
lapply(packages_names, require, character.only = TRUE)

#experiment model
Experiment_Model <- model.matrix(~ 0 + factor(rep(c("case1", "case2"), each = 3), levels = c("case1", "case2")))
colnames(Experiment_Model) <- c("case1", "case2")
Contrast_matrix <- makeContrasts(case1_case2 = case1 - case2, levels = Experiment_Model)

#HEK293T OG
#OG glycopeptide differential analysis
OG_glycopeptide_raw_sl_tmm_log2_HEK293T <- OG_glycopeptide_raw_sl_tmm_HEK293T |> 
  mutate(
    log2_Tuni_1_sl_tmm = log2(Tuni_1_sl_tmm),
    log2_Tuni_2_sl_tmm = log2(Tuni_2_sl_tmm),
    log2_Tuni_3_sl_tmm = log2(Tuni_3_sl_tmm),
    log2_Ctrl_4_sl_tmm = log2(Ctrl_4_sl_tmm),
    log2_Ctrl_5_sl_tmm = log2(Ctrl_5_sl_tmm),
    log2_Ctrl_6_sl_tmm = log2(Ctrl_6_sl_tmm)
  )

OG_glycopeptide_raw_sl_tmm_log2_Data_Matrix_HEK293T <- data.matrix(OG_glycopeptide_raw_sl_tmm_log2_HEK293T |> select(starts_with("log2")))
rownames(OG_glycopeptide_raw_sl_tmm_log2_Data_Matrix_HEK293T) <- OG_glycopeptide_raw_sl_tmm_log2_HEK293T$Index

OG_glycopeptide_raw_sl_tmm_log2_Fit_HEK293T <- lmFit(OG_glycopeptide_raw_sl_tmm_log2_Data_Matrix_HEK293T, Experiment_Model)
OG_glycopeptide_raw_sl_tmm_log2_Fit_Contrast_HEK293T <- contrasts.fit(OG_glycopeptide_raw_sl_tmm_log2_Fit_HEK293T, Contrast_matrix)
OG_glycopeptide_raw_sl_tmm_log2_Fit_Contrast_HEK293T <- eBayes(OG_glycopeptide_raw_sl_tmm_log2_Fit_Contrast_HEK293T)
OG_glycopeptide_Top_HEK293T <- topTable(OG_glycopeptide_raw_sl_tmm_log2_Fit_Contrast_HEK293T, number = Inf, adjust.method = "BH")

Rownames_OG_glycopeptide_Top_HEK293T <- rownames(OG_glycopeptide_Top_HEK293T)
OG_glycopeptide_Top_tb_HEK293T <- as_tibble(OG_glycopeptide_Top_HEK293T)
OG_glycopeptide_Top_tb_HEK293T$Index <- Rownames_OG_glycopeptide_Top_HEK293T

OG_glycopeptide_Top_tb_HEK293T <- OG_glycopeptide_Top_tb_HEK293T |> 
  left_join(OG_glycopeptide_raw_sl_tmm_HEK293T |> select(Index:localized), by = join_by(Index == Index)) |> 
  mutate(site_position = str_extract(combined_site, "\\d+")) %>%
  mutate(
    site_position = as.numeric(site_position)
  )

write_xlsx(OG_glycopeptide_Top_tb_HEK293T, path = paste0(file_path, "OG_glycopeptide_Top_tb_HEK293T.xlsx"))

#OG glycoprotein differential analysis
OG_glycoprotein_raw_sl_tmm_log2_HEK293T <- OG_glycoprotein_raw_sl_tmm_HEK293T |> 
  mutate(
    log2_Tuni_1_sl_tmm = log2(Tuni_1_sl_tmm),
    log2_Tuni_2_sl_tmm = log2(Tuni_2_sl_tmm),
    log2_Tuni_3_sl_tmm = log2(Tuni_3_sl_tmm),
    log2_Ctrl_4_sl_tmm = log2(Ctrl_4_sl_tmm),
    log2_Ctrl_5_sl_tmm = log2(Ctrl_5_sl_tmm),
    log2_Ctrl_6_sl_tmm = log2(Ctrl_6_sl_tmm)
  )

OG_glycoprotein_raw_sl_tmm_log2_Data_Matrix_HEK293T <- data.matrix(OG_glycoprotein_raw_sl_tmm_log2_HEK293T |> select(starts_with("log2")))
rownames(OG_glycoprotein_raw_sl_tmm_log2_Data_Matrix_HEK293T) <- OG_glycoprotein_raw_sl_tmm_log2_HEK293T$UniprotID

OG_glycoprotein_raw_sl_tmm_log2_Fit_HEK293T <- lmFit(OG_glycoprotein_raw_sl_tmm_log2_Data_Matrix_HEK293T, Experiment_Model)
OG_glycoprotein_raw_sl_tmm_log2_Fit_Contrast_HEK293T <- contrasts.fit(OG_glycoprotein_raw_sl_tmm_log2_Fit_HEK293T, Contrast_matrix)
OG_glycoprotein_raw_sl_tmm_log2_Fit_Contrast_HEK293T <- eBayes(OG_glycoprotein_raw_sl_tmm_log2_Fit_Contrast_HEK293T)
OG_glycoprotein_Top_HEK293T <- topTable(OG_glycoprotein_raw_sl_tmm_log2_Fit_Contrast_HEK293T, number = Inf, adjust.method = "BH")

Rownames_OG_glycoprotein_Top_HEK293T <- rownames(OG_glycoprotein_Top_HEK293T)
OG_glycoprotein_Top_tb_HEK293T <- as_tibble(OG_glycoprotein_Top_HEK293T)
OG_glycoprotein_Top_tb_HEK293T$UniprotID <- Rownames_OG_glycoprotein_Top_HEK293T

write_xlsx(OG_glycoprotein_Top_tb_HEK293T, path = paste0(file_path, "OG_glycoprotein_Top_tb_HEK293T.xlsx"))

#HepG2 OG
#OG glycopeptide differential analysis
OG_glycopeptide_raw_sl_tmm_log2_HepG2 <- OG_glycopeptide_raw_sl_tmm_HepG2 |> 
  mutate(
    log2_Tuni_1_sl_tmm = log2(Tuni_1_sl_tmm),
    log2_Tuni_2_sl_tmm = log2(Tuni_2_sl_tmm),
    log2_Tuni_3_sl_tmm = log2(Tuni_3_sl_tmm),
    log2_Ctrl_4_sl_tmm = log2(Ctrl_4_sl_tmm),
    log2_Ctrl_5_sl_tmm = log2(Ctrl_5_sl_tmm),
    log2_Ctrl_6_sl_tmm = log2(Ctrl_6_sl_tmm)
  )

OG_glycopeptide_raw_sl_tmm_log2_Data_Matrix_HepG2 <- data.matrix(OG_glycopeptide_raw_sl_tmm_log2_HepG2 |> select(starts_with("log2")))
rownames(OG_glycopeptide_raw_sl_tmm_log2_Data_Matrix_HepG2) <- OG_glycopeptide_raw_sl_tmm_log2_HepG2$Index

OG_glycopeptide_raw_sl_tmm_log2_Fit_HepG2 <- lmFit(OG_glycopeptide_raw_sl_tmm_log2_Data_Matrix_HepG2, Experiment_Model)
OG_glycopeptide_raw_sl_tmm_log2_Fit_Contrast_HepG2 <- contrasts.fit(OG_glycopeptide_raw_sl_tmm_log2_Fit_HepG2, Contrast_matrix)
OG_glycopeptide_raw_sl_tmm_log2_Fit_Contrast_HepG2 <- eBayes(OG_glycopeptide_raw_sl_tmm_log2_Fit_Contrast_HepG2)
OG_glycopeptide_Top_HepG2 <- topTable(OG_glycopeptide_raw_sl_tmm_log2_Fit_Contrast_HepG2, number = Inf, adjust.method = "BH")

Rownames_OG_glycopeptide_Top_HepG2 <- rownames(OG_glycopeptide_Top_HepG2)
OG_glycopeptide_Top_tb_HepG2 <- as_tibble(OG_glycopeptide_Top_HepG2)
OG_glycopeptide_Top_tb_HepG2$Index <- Rownames_OG_glycopeptide_Top_HepG2

OG_glycopeptide_Top_tb_HepG2 <- OG_glycopeptide_Top_tb_HepG2 |> 
  left_join(OG_glycopeptide_raw_sl_tmm_HepG2 |> select(Index:localized), by = join_by(Index == Index)) |> 
  mutate(site_position = str_extract(combined_site, "\\d+")) %>%
  mutate(
    site_position = as.numeric(site_position)
  )

write_xlsx(OG_glycopeptide_Top_tb_HepG2, path = paste0(file_path, "OG_glycopeptide_Top_tb_HepG2.xlsx"))

#OG glycoprotein differential analysis
OG_glycoprotein_raw_sl_tmm_log2_HepG2 <- OG_glycoprotein_raw_sl_tmm_HepG2 |> 
  mutate(
    log2_Tuni_1_sl_tmm = log2(Tuni_1_sl_tmm),
    log2_Tuni_2_sl_tmm = log2(Tuni_2_sl_tmm),
    log2_Tuni_3_sl_tmm = log2(Tuni_3_sl_tmm),
    log2_Ctrl_4_sl_tmm = log2(Ctrl_4_sl_tmm),
    log2_Ctrl_5_sl_tmm = log2(Ctrl_5_sl_tmm),
    log2_Ctrl_6_sl_tmm = log2(Ctrl_6_sl_tmm)
  )

OG_glycoprotein_raw_sl_tmm_log2_Data_Matrix_HepG2 <- data.matrix(OG_glycoprotein_raw_sl_tmm_log2_HepG2 |> select(starts_with("log2")))
rownames(OG_glycoprotein_raw_sl_tmm_log2_Data_Matrix_HepG2) <- OG_glycoprotein_raw_sl_tmm_log2_HepG2$UniprotID

OG_glycoprotein_raw_sl_tmm_log2_Fit_HepG2 <- lmFit(OG_glycoprotein_raw_sl_tmm_log2_Data_Matrix_HepG2, Experiment_Model)
OG_glycoprotein_raw_sl_tmm_log2_Fit_Contrast_HepG2 <- contrasts.fit(OG_glycoprotein_raw_sl_tmm_log2_Fit_HepG2, Contrast_matrix)
OG_glycoprotein_raw_sl_tmm_log2_Fit_Contrast_HepG2 <- eBayes(OG_glycoprotein_raw_sl_tmm_log2_Fit_Contrast_HepG2)
OG_glycoprotein_Top_HepG2 <- topTable(OG_glycoprotein_raw_sl_tmm_log2_Fit_Contrast_HepG2, number = Inf, adjust.method = "BH")

Rownames_OG_glycoprotein_Top_HepG2 <- rownames(OG_glycoprotein_Top_HepG2)
OG_glycoprotein_Top_tb_HepG2 <- as_tibble(OG_glycoprotein_Top_HepG2)
OG_glycoprotein_Top_tb_HepG2$UniprotID <- Rownames_OG_glycoprotein_Top_HepG2

write_xlsx(OG_glycoprotein_Top_tb_HepG2, path = paste0(file_path, "OG_glycoprotein_Top_tb_HepG2.xlsx"))

#Jurkat OG
#OG glycopeptide differential analysis
OG_glycopeptide_raw_sl_tmm_log2_Jurkat <- OG_glycopeptide_raw_sl_tmm_Jurkat |> 
  mutate(
    log2_Tuni_1_sl_tmm = log2(Tuni_1_sl_tmm),
    log2_Tuni_2_sl_tmm = log2(Tuni_2_sl_tmm),
    log2_Tuni_3_sl_tmm = log2(Tuni_3_sl_tmm),
    log2_Ctrl_4_sl_tmm = log2(Ctrl_4_sl_tmm),
    log2_Ctrl_5_sl_tmm = log2(Ctrl_5_sl_tmm),
    log2_Ctrl_6_sl_tmm = log2(Ctrl_6_sl_tmm)
  )

OG_glycopeptide_raw_sl_tmm_log2_Data_Matrix_Jurkat <- data.matrix(OG_glycopeptide_raw_sl_tmm_log2_Jurkat |> select(starts_with("log2")))
rownames(OG_glycopeptide_raw_sl_tmm_log2_Data_Matrix_Jurkat) <- OG_glycopeptide_raw_sl_tmm_log2_Jurkat$Index

OG_glycopeptide_raw_sl_tmm_log2_Fit_Jurkat <- lmFit(OG_glycopeptide_raw_sl_tmm_log2_Data_Matrix_Jurkat, Experiment_Model)
OG_glycopeptide_raw_sl_tmm_log2_Fit_Contrast_Jurkat <- contrasts.fit(OG_glycopeptide_raw_sl_tmm_log2_Fit_Jurkat, Contrast_matrix)
OG_glycopeptide_raw_sl_tmm_log2_Fit_Contrast_Jurkat <- eBayes(OG_glycopeptide_raw_sl_tmm_log2_Fit_Contrast_Jurkat)
OG_glycopeptide_Top_Jurkat <- topTable(OG_glycopeptide_raw_sl_tmm_log2_Fit_Contrast_Jurkat, number = Inf, adjust.method = "BH")

Rownames_OG_glycopeptide_Top_Jurkat <- rownames(OG_glycopeptide_Top_Jurkat)
OG_glycopeptide_Top_tb_Jurkat <- as_tibble(OG_glycopeptide_Top_Jurkat)
OG_glycopeptide_Top_tb_Jurkat$Index <- Rownames_OG_glycopeptide_Top_Jurkat

OG_glycopeptide_Top_tb_Jurkat <- OG_glycopeptide_Top_tb_Jurkat |> 
  left_join(OG_glycopeptide_raw_sl_tmm_Jurkat |> select(Index:localized), by = join_by(Index == Index)) |> 
  mutate(site_position = str_extract(combined_site, "\\d+")) %>%
  mutate(
    site_position = as.numeric(site_position)
  )

write_xlsx(OG_glycopeptide_Top_tb_Jurkat, path = paste0(file_path, "OG_glycopeptide_Top_tb_Jurkat.xlsx"))

#OG glycoprotein differential analysis
OG_glycoprotein_raw_sl_tmm_log2_Jurkat <- OG_glycoprotein_raw_sl_tmm_Jurkat |> 
  mutate(
    log2_Tuni_1_sl_tmm = log2(Tuni_1_sl_tmm),
    log2_Tuni_2_sl_tmm = log2(Tuni_2_sl_tmm),
    log2_Tuni_3_sl_tmm = log2(Tuni_3_sl_tmm),
    log2_Ctrl_4_sl_tmm = log2(Ctrl_4_sl_tmm),
    log2_Ctrl_5_sl_tmm = log2(Ctrl_5_sl_tmm),
    log2_Ctrl_6_sl_tmm = log2(Ctrl_6_sl_tmm)
  )

OG_glycoprotein_raw_sl_tmm_log2_Data_Matrix_Jurkat <- data.matrix(OG_glycoprotein_raw_sl_tmm_log2_Jurkat |> select(starts_with("log2")))
rownames(OG_glycoprotein_raw_sl_tmm_log2_Data_Matrix_Jurkat) <- OG_glycoprotein_raw_sl_tmm_log2_Jurkat$UniprotID

OG_glycoprotein_raw_sl_tmm_log2_Fit_Jurkat <- lmFit(OG_glycoprotein_raw_sl_tmm_log2_Data_Matrix_Jurkat, Experiment_Model)
OG_glycoprotein_raw_sl_tmm_log2_Fit_Contrast_Jurkat <- contrasts.fit(OG_glycoprotein_raw_sl_tmm_log2_Fit_Jurkat, Contrast_matrix)
OG_glycoprotein_raw_sl_tmm_log2_Fit_Contrast_Jurkat <- eBayes(OG_glycoprotein_raw_sl_tmm_log2_Fit_Contrast_Jurkat)
OG_glycoprotein_Top_Jurkat <- topTable(OG_glycoprotein_raw_sl_tmm_log2_Fit_Contrast_Jurkat, number = Inf, adjust.method = "BH")

Rownames_OG_glycoprotein_Top_Jurkat <- rownames(OG_glycoprotein_Top_Jurkat)
OG_glycoprotein_Top_tb_Jurkat <- as_tibble(OG_glycoprotein_Top_Jurkat)
OG_glycoprotein_Top_tb_Jurkat$UniprotID <- Rownames_OG_glycoprotein_Top_Jurkat

write_xlsx(OG_glycoprotein_Top_tb_Jurkat, path = paste0(file_path, "OG_glycoprotein_Top_tb_Jurkat.xlsx"))

#HEK293T WP
#WP differential analysis
WP_protein_raw_sl_tmm_log2_HEK293T <- WP_protein_raw_sl_tmm_HEK293T |> 
  mutate(
    log2_Tuni_1_sl_tmm = log2(Tuni_1_sl_tmm),
    log2_Tuni_2_sl_tmm = log2(Tuni_2_sl_tmm),
    log2_Tuni_3_sl_tmm = log2(Tuni_3_sl_tmm),
    log2_Ctrl_4_sl_tmm = log2(Ctrl_4_sl_tmm),
    log2_Ctrl_5_sl_tmm = log2(Ctrl_5_sl_tmm),
    log2_Ctrl_6_sl_tmm = log2(Ctrl_6_sl_tmm)
  )

WP_protein_raw_sl_tmm_log2_Data_Matrix_HEK293T <- data.matrix(WP_protein_raw_sl_tmm_log2_HEK293T |> select(log2_Tuni_1_sl_tmm:log2_Ctrl_6_sl_tmm))
rownames(WP_protein_raw_sl_tmm_log2_Data_Matrix_HEK293T) <- WP_protein_raw_sl_tmm_log2_HEK293T$UniprotID

WP_protein_Fit_HEK293T <- lmFit(WP_protein_raw_sl_tmm_log2_Data_Matrix_HEK293T, Experiment_Model)
WP_protein_Fit_Contrast_HEK293T <- contrasts.fit(WP_protein_Fit_HEK293T, Contrast_matrix)
WP_protein_Fit_Contrast_HEK293T <- eBayes(WP_protein_Fit_Contrast_HEK293T)

WP_protein_Top_HEK293T <- topTable(WP_protein_Fit_Contrast_HEK293T, number = Inf, adjust.method = "BH")

Rownames_WP_protein_Top_HEK293T <- rownames(WP_protein_Top_HEK293T)
WP_protein_Top_tb_HEK293T <- as_tibble(WP_protein_Top_HEK293T)
WP_protein_Top_tb_HEK293T$UniprotID <- Rownames_WP_protein_Top_HEK293T

write_xlsx(WP_protein_Top_tb_HEK293T, paste0(file_path, "WP_protein_Top_tb_HEK293T.xlsx"))

#HepG2 WP
WP_protein_raw_sl_tmm_log2_HepG2 <- WP_protein_raw_sl_tmm_HepG2 |> 
  mutate(
    log2_Tuni_1_sl_tmm = log2(Tuni_1_sl_tmm),
    log2_Tuni_2_sl_tmm = log2(Tuni_2_sl_tmm),
    log2_Tuni_3_sl_tmm = log2(Tuni_3_sl_tmm),
    log2_Ctrl_4_sl_tmm = log2(Ctrl_4_sl_tmm),
    log2_Ctrl_5_sl_tmm = log2(Ctrl_5_sl_tmm),
    log2_Ctrl_6_sl_tmm = log2(Ctrl_6_sl_tmm)
  )

WP_protein_raw_sl_tmm_log2_Data_Matrix_HepG2 <- data.matrix(WP_protein_raw_sl_tmm_log2_HepG2 |> select(log2_Tuni_1_sl_tmm:log2_Ctrl_6_sl_tmm))
rownames(WP_protein_raw_sl_tmm_log2_Data_Matrix_HepG2) <- WP_protein_raw_sl_tmm_log2_HepG2$UniprotID

WP_protein_Fit_HepG2 <- lmFit(WP_protein_raw_sl_tmm_log2_Data_Matrix_HepG2, Experiment_Model)
WP_protein_Fit_Contrast_HepG2 <- contrasts.fit(WP_protein_Fit_HepG2, Contrast_matrix)
WP_protein_Fit_Contrast_HepG2 <- eBayes(WP_protein_Fit_Contrast_HepG2)

WP_protein_Top_HepG2 <- topTable(WP_protein_Fit_Contrast_HepG2, number = Inf, adjust.method = "BH")

Rownames_WP_protein_Top_HepG2 <- rownames(WP_protein_Top_HepG2)
WP_protein_Top_tb_HepG2 <- as_tibble(WP_protein_Top_HepG2)
WP_protein_Top_tb_HepG2$UniprotID <- Rownames_WP_protein_Top_HepG2

write_xlsx(WP_protein_Top_tb_HepG2, paste0(file_path, "WP_protein_Top_tb_HepG2.xlsx"))

#Jurkat WP
#WP differential analysis
WP_protein_raw_sl_tmm_log2_Jurkat <- WP_protein_raw_sl_tmm_Jurkat |> 
  mutate(
    log2_Tuni_1_sl_tmm = log2(Tuni_1_sl_tmm),
    log2_Tuni_2_sl_tmm = log2(Tuni_2_sl_tmm),
    log2_Tuni_3_sl_tmm = log2(Tuni_3_sl_tmm),
    log2_Ctrl_4_sl_tmm = log2(Ctrl_4_sl_tmm),
    log2_Ctrl_5_sl_tmm = log2(Ctrl_5_sl_tmm),
    log2_Ctrl_6_sl_tmm = log2(Ctrl_6_sl_tmm)
  )

WP_protein_raw_sl_tmm_log2_Data_Matrix_Jurkat <- data.matrix(WP_protein_raw_sl_tmm_log2_Jurkat |> select(log2_Tuni_1_sl_tmm:log2_Ctrl_6_sl_tmm))
rownames(WP_protein_raw_sl_tmm_log2_Data_Matrix_Jurkat) <- WP_protein_raw_sl_tmm_log2_Jurkat$UniprotID

WP_protein_Fit_Jurkat <- lmFit(WP_protein_raw_sl_tmm_log2_Data_Matrix_Jurkat, Experiment_Model)
WP_protein_Fit_Contrast_Jurkat <- contrasts.fit(WP_protein_Fit_Jurkat, Contrast_matrix)
WP_protein_Fit_Contrast_Jurkat <- eBayes(WP_protein_Fit_Contrast_Jurkat)

WP_protein_Top_Jurkat <- topTable(WP_protein_Fit_Contrast_Jurkat, number = Inf, adjust.method = "BH")

Rownames_WP_protein_Top_Jurkat <- rownames(WP_protein_Top_Jurkat)
WP_protein_Top_tb_Jurkat <- as_tibble(WP_protein_Top_Jurkat)
WP_protein_Top_tb_Jurkat$UniprotID <- Rownames_WP_protein_Top_Jurkat

write_xlsx(WP_protein_Top_tb_Jurkat, paste0(file_path, "WP_protein_Top_tb_Jurkat.xlsx"))
