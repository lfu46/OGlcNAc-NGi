#import packages
packages_names <- c("tidyverse", "limma", "showtext", "writexl")
lapply(packages_names, require, character.only = TRUE)

#experiment model
Experiment_Model <- model.matrix(~ 0 + factor(rep(c("case1", "case2"), each = 3), levels = c("case1", "case2")))
colnames(Experiment_Model) <- c("case1", "case2")
Contrast_matrix <- makeContrasts(case1_case2 = case1 - case2, levels = Experiment_Model)

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
