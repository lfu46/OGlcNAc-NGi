#import packages
packages_names <- c("tidyverse", "limma", "showtext", "writexl")
lapply(packages_names, require, character.only = TRUE)

#experiment model
Experiment_Model <- model.matrix(~ 0 + factor(rep(c("case1", "case2"), each = 3), levels = c("case1", "case2")))
colnames(Experiment_Model) <- c("case1", "case2")
Contrast_matrix <- makeContrasts(case1_case2 = case1 - case2, levels = Experiment_Model)

#generate data frame
OG_glycoprotein_protein_OGlcNAc_effect_HepG2 <- OG_glycoprotein_raw_sl_tmm_HepG2 |> 
  select(UniprotID, ends_with("sl_tmm")) |> 
  mutate(
    log2_ratio1_OG_glycoprotein = log2(Tuni_1_sl_tmm / Ctrl_4_sl_tmm),
    log2_ratio2_OG_glycoprotein = log2(Tuni_2_sl_tmm / Ctrl_5_sl_tmm),
    log2_ratio3_OG_glycoprotein = log2(Tuni_3_sl_tmm / Ctrl_6_sl_tmm)
  ) |> 
  select(UniprotID, starts_with("log2")) |> 
  left_join(WP_protein_raw_sl_tmm_HepG2, by = "UniprotID") |> 
  select(UniprotID, log2_ratio1_OG_glycoprotein:log2_ratio3_OG_glycoprotein, ends_with("sl_tmm")) |> 
  mutate(
    log2_ratio1_WP = log2(Tuni_1_sl_tmm / Ctrl_4_sl_tmm),
    log2_ratio2_WP = log2(Tuni_2_sl_tmm / Ctrl_5_sl_tmm),
    log2_ratio3_WP = log2(Tuni_3_sl_tmm / Ctrl_6_sl_tmm)
  ) |> 
  select(UniprotID, starts_with("log2")) |> 
  filter(!is.na(log2_ratio1_WP))

#differential analysis
data_matrix <- data.matrix(OG_glycoprotein_protein_OGlcNAc_effect_HepG2 |> select(starts_with("log2")))
rownames(data_matrix) <- OG_glycoprotein_protein_OGlcNAc_effect_HepG2$UniprotID

lm_fit <- lmFit(data_matrix, Experiment_Model)
fit_contrast <- contrasts.fit(lm_fit, Contrast_matrix)
fit_contrast <- eBayes(fit_contrast)
top_table <- topTable(fit_contrast, number = Inf, adjust.method = "BH")

rownames_HepG2 <- rownames(top_table)
top_tb <- tibble(top_table)
top_tb$UniprotID <- rownames_HepG2

OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 <- top_tb

write_xlsx(OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2, path = paste0(file_path, "OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2.xlsx"))
