#import packages
packages_names <- c("tidyverse", "limma", "showtext", "writexl")
lapply(packages_names, require, character.only = TRUE)

#experimental model
Experiment_Model <- model.matrix(~ 0 + factor(rep(c("case1", "case2"), each = 3), levels = c("case1", "case2")))
colnames(Experiment_Model) <- c("case1", "case2")
Contrast_matrix <- makeContrasts(case1_case2 = case1 - case2, levels = Experiment_Model)

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
