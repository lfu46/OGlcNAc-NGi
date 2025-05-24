#import packages
packages_names <- c("tidyverse", "limma", "showtext", "writexl")
lapply(packages_names, require, character.only = TRUE)

#experimental model
Experiment_Model <- model.matrix(~ 0 + factor(rep(c("case1", "case2"), each = 3), levels = c("case1", "case2")))
colnames(Experiment_Model) <- c("case1", "case2")
Contrast_matrix <- makeContrasts(case1_case2 = case1 - case2, levels = Experiment_Model)

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
