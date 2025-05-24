#import packages
packages_names <- c("tidyverse", "limma", "showtext", "writexl")
lapply(packages_names, require, character.only = TRUE)

#experiment model
Experiment_Model <- model.matrix(~ 0 + factor(rep(c("case1", "case2"), each = 3), levels = c("case1", "case2")))
colnames(Experiment_Model) <- c("case1", "case2")
Contrast_matrix <- makeContrasts(case1_case2 = case1 - case2, levels = Experiment_Model)

#generate data frame
OG_site_protein_OGlcNAc_effect_HEK293T <- OG_glycopeptide_raw_sl_tmm_HEK293T |> 
  filter(!is.na(combined_site)) |> 
  select(Index, UniprotID, combined_site, ends_with("sl_tmm")) |> 
  mutate(
    log2_ratio1_OG_site = log2(Tuni_1_sl_tmm / Ctrl_4_sl_tmm),
    log2_ratio2_OG_site = log2(Tuni_2_sl_tmm / Ctrl_5_sl_tmm),
    log2_ratio3_OG_site = log2(Tuni_3_sl_tmm / Ctrl_6_sl_tmm)
  ) |> 
  select(Index, UniprotID, combined_site, starts_with("log2")) |> 
  left_join(WP_protein_raw_sl_tmm_HEK293T, by = "UniprotID") |> 
  select(Index, UniprotID, combined_site, log2_ratio1_OG_site:log2_ratio3_OG_site, ends_with("sl_tmm")) |> 
  mutate(
    log2_ratio1_WP = log2(Tuni_1_sl_tmm / Ctrl_4_sl_tmm),
    log2_ratio2_WP = log2(Tuni_2_sl_tmm / Ctrl_5_sl_tmm),
    log2_ratio3_WP = log2(Tuni_3_sl_tmm / Ctrl_6_sl_tmm)
  ) |> 
  select(Index, UniprotID, combined_site, starts_with("log2")) |> 
  filter(!is.na(log2_ratio1_WP))

#differential analysis
data_matrix <- data.matrix(OG_site_protein_OGlcNAc_effect_HEK293T |> select(starts_with("log2")))
rownames(data_matrix) <- OG_site_protein_OGlcNAc_effect_HEK293T$Index

lm_fit <- lmFit(data_matrix, Experiment_Model)
fit_contrast <- contrasts.fit(lm_fit, Contrast_matrix)
fit_contrast <- eBayes(fit_contrast)
top_table <- topTable(fit_contrast, number = Inf, adjust.method = "BH")

rownames_HEK293T <- rownames(top_table)
top_tb <- tibble(top_table)
top_tb$Index <- rownames_HEK293T

OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T <- top_tb |> 
  left_join(OG_site_protein_OGlcNAc_effect_HEK293T |> select(Index, UniprotID, combined_site), by = "Index") |> 
  mutate(
    site_position = str_extract(combined_site, "\\d+"),
    site_position = as.numeric(site_position)
  )

write_xlsx(OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T, path = paste0(file_path, "OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T.xlsx"))
