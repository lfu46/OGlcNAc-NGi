#import packages
packages_names <- c("showtext", "ggpubr", "rstatix")
lapply(packages_names, require, character.only = TRUE)

#Glycoprotein Distribution
#generate data frame
OG_glycoprotein_logFC_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  mutate(cell = "HepG2") |> 
  select(UniprotID, logFC, cell)

OG_glycoprotein_logFC_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  mutate(cell = "HEK293T") |> 
  select(UniprotID, logFC, cell)

OG_glycoprotein_logFC_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  mutate(cell = "Jurkat") |> 
  select(UniprotID, logFC, cell)

OG_glycoprotein_logFC_combined <- bind_rows(OG_glycoprotein_logFC_HepG2, OG_glycoprotein_logFC_HEK293T, OG_glycoprotein_logFC_Jurkat)

#ks.test
ks.test(OG_glycoprotein_Top_tb_HepG2$logFC, OG_glycoprotein_Top_tb_HEK293T$logFC)
ks.test(OG_glycoprotein_Top_tb_HEK293T$logFC, OG_glycoprotein_Top_tb_Jurkat$logFC)
ks.test(OG_glycoprotein_Top_tb_Jurkat$logFC, OG_glycoprotein_Top_tb_HepG2$logFC)

OG_glycoprotein_ks_test <- tribble(
  ~ .y., ~ group1, ~ group2, ~ p, ~ p.signif,
  "logFC", "HEK293T", "HepG2", 2.153e-13, "****",
  "logFC", "Jurkat", "HEK293T", 0.04457, "*",
  "logFC", "HepG2", "Jurkat", 6.29e-15, "****"
)

#wilcox test
OG_glycoprotein_logFC_wilcox_test <- OG_glycoprotein_logFC_combined |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_logFC <- OG_glycoprotein_logFC_combined |>
  ggplot() +
  geom_violin(aes(x = Cell, y = logFC, fill = Cell), color = "transparent") +
  geom_boxplot(aes(x = Cell, y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  scale_fill_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  labs(x = "", y = expression(log[2](Tuni/Ctrl))) +
  stat_pvalue_manual(data = OG_glycoprotein_ks_test, label = "p.signif", tip.length = 0, 
                     size = 40, y.position = c(2.0, 2.4, 2.8)) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_logFC.png"), plot = violin_boxplot_OG_glycoprotein_logFC, height = 3, width = 4, units = "in", dpi = 600)

#Glycoprotein vs. non-modified form
#generate data frame for HepG2
OG_WP_distribution_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  select(UniprotID, logFC_OG = logFC) |> 
  left_join(WP_protein_Top_tb_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC_OG, logFC_WP = logFC) |> 
  filter(!is.na(logFC_WP)) |> 
  pivot_longer(cols = logFC_OG:logFC_WP, names_to = "Exp", values_to = "logFC") |> 
  mutate(Cell = "HepG2")

#ks test
ks.test(
  OG_WP_distribution_HepG2 |> filter(Exp == "logFC_OG") |> pull(logFC),
  OG_WP_distribution_HepG2 |> filter(Exp == "logFC_WP") |> pull(logFC)
)

#wilcox test
OG_WP_distribution_HepG2_wilcox_test <- OG_WP_distribution_HepG2 |> 
  wilcox_test(logFC ~ Exp, p.adjust.method = "BH") |> 
  add_significance("p")

#generate data frame for HEK293T
OG_WP_distribution_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  select(UniprotID, logFC_OG = logFC) |> 
  left_join(WP_protein_Top_tb_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC_OG, logFC_WP = logFC) |> 
  filter(!is.na(logFC_WP)) |> 
  pivot_longer(cols = logFC_OG:logFC_WP, names_to = "Exp", values_to = "logFC") |> 
  mutate(Cell = "HEK293T")

#ks test
ks.test(
  OG_WP_distribution_HEK293T |> filter(Exp == "logFC_OG") |> pull(logFC),
  OG_WP_distribution_HEK293T |> filter(Exp == "logFC_WP") |> pull(logFC)
)

#wilcox test
OG_WP_distribution_HEK293T_wilcox_test <- OG_WP_distribution_HEK293T |> 
  wilcox_test(logFC ~ Exp, p.adjust.method = "BH") |> 
  add_significance("p")

#generate data frame for Jurkat
OG_WP_distribution_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  select(UniprotID, logFC_OG = logFC) |> 
  left_join(WP_protein_Top_tb_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  select(UniprotID, logFC_OG, logFC_WP = logFC) |> 
  filter(!is.na(logFC_WP)) |> 
  pivot_longer(cols = logFC_OG:logFC_WP, names_to = "Exp", values_to = "logFC") |> 
  mutate(Cell = "Jurkat")

#ks test
ks.test(
  OG_WP_distribution_Jurkat |> filter(Exp == "logFC_OG") |> pull(logFC),
  OG_WP_distribution_Jurkat |> filter(Exp == "logFC_WP") |> pull(logFC)
)

#wilcox test
OG_WP_distribution_Jurkat_wilcox_test <- OG_WP_distribution_Jurkat |> 
  wilcox_test(logFC ~ Exp, p.adjust.method = "BH") |> 
  add_significance("p")

#combine
OG_WP_distribution_combined <- bind_rows(
  OG_WP_distribution_HepG2,
  OG_WP_distribution_HEK293T,
  OG_WP_distribution_Jurkat
)

#combine ks test

#wilcox test
OG_WP_distribution_combined_wilcox_test <- OG_WP_distribution_combined |> 
  group_by(Cell) |> 
  wilcox_test(logFC ~ Exp, p.adjust.method = "BH") |> 
  add_significance("p")

#split violin plot
ks_test_label <- tibble(
  Cell = c("HEK293T", "HEK293T", "HepG2", "HepG2", "Jurkat", "Jurkat"),
  p_signif = c("****", "****", "****", "****", "****", "****"),
  logFC = 1.5,
  Name = c("glycoprotein_HEK293T", "glycoprotein_HEK293T", "glycoprotein_HepG2", "glycoprotein_HepG2", "glycoprotein_Jurkat", "glycoprotein_Jurkat"),
  Exp = c("OG_HEK293T", "WP", "OG_HepG2", "WP", "OG_Jurkat", "WP")
)

font_add(family = "Arial", regular = "arial.ttf")
showtext_auto()

split_violin_plot_OG_WP_distribution <- OG_WP_distribution_combined |> 
  mutate(Name = paste("glycoprotein", Cell, sep = "_"), Exp = ifelse(Exp == "logFC_OG", paste("OG", Cell, sep = "_"), "WP")) |> 
  ggplot(aes(x = Name, y = logFC, fill = factor(Exp, c("OG_HEK293T", "OG_HepG2", "OG_Jurkat", "WP")))) +
  geom_split_violin(color = "transparent") +
  geom_text(data = ks_test_label, aes(y = logFC, label = p_signif), size = 40, hjust = -0.2) +
  scale_fill_manual(values = c(
    "OG_HEK293T" = Color_2,
    "OG_HepG2" = Color_3,
    "OG_Jurkat" = Color_4,
    "WP" = "gray"
  )) +
  facet_grid(~ Cell, scale = "free_x") +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)"), fill = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 80, color = "black"),
    axis.line.x = element_blank(),
    axis.ticks.length.x = unit(0, "in"),
    legend.position = "none",
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "split_violin_plot_OG_WP_distribution.png"), plot = split_violin_plot_OG_WP_distribution, height = 3, width = 4, units = c("in"), dpi = 600)

#Correlation matrix
#generate OG glycopeptide logFC overlap
tb1 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_common_total$UniprotID) |> 
  select(UniprotID, HepG2_logFC = logFC)

tb2 <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_common_total$UniprotID) |> 
  select(UniprotID, HEK293T_logFC = logFC)

tb3 <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(UniprotID %in% OG_glycoprotein_common_total$UniprotID) |> 
  select(UniprotID, Jurkat_logFC = logFC)

OG_glycoprotein_logFC_common <- tb1 |> 
  left_join(tb2, by = join_by(UniprotID == UniprotID)) |> 
  left_join(tb3, by = join_by(UniprotID == UniprotID))

#calculate person's r between cell lines
correlation_HepG2_HEK293T <- cor(OG_glycoprotein_logFC_common$HepG2_logFC, OG_glycoprotein_logFC_common$HEK293T_logFC, method = "pearson")
correlation_HEK293T_Jurkat <- cor(OG_glycoprotein_logFC_common$HEK293T_logFC, OG_glycoprotein_logFC_common$Jurkat_logFC, method = "pearson")

#correlation plot between different cell lines
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

correlation_plot_OG_glycoprotein_common <- ggplot(data = OG_glycoprotein_logFC_common) +
  geom_point(aes(x = HEK293T_logFC, y = HepG2_logFC), shape = 21, color = Color_3, size = 3) +
  geom_point(aes(x = HEK293T_logFC, y = Jurkat_logFC), shape = 21, color = Color_4, size = 3) +
  labs(x = "log2(Fold change) in HEK293T", y = "log2(Fold change) in HepG2 or Jurkat") +
  coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
  annotate("point", x = -2, y = 2, shape = 21, size = 2, color = Color_3) +
  annotate("text", x = -2, y = 2, size = 25, hjust = -0.12, color = "black", label = "HepG2 (r = 0.44)") +
  annotate("point", x = -2, y = 1.8, shape = 21, size = 2, color = Color_4) +
  annotate("text", x = -2, y = 1.8, size = 25, hjust = -0.12, color = "black", label = "Jurkat (r = 0.11)") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 80),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black")
  )

ggsave(filename = paste0(file_path, "correlation_plot_OG_glycoprotein_common.png"), plot = correlation_plot_OG_glycoprotein_common, height = 4, width = 4, dpi = 600)

#correlation matrix plot
OG_glycoprotein_raw_sl_tmm_common <- OG_glycoprotein_raw_sl_tmm_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_common_total$UniprotID) |> 
  mutate(
    ratio_1_HepG2 = Tuni_1_sl_tmm/Ctrl_4_sl_tmm,
    ratio_2_HepG2 = Tuni_2_sl_tmm/Ctrl_5_sl_tmm,
    ratio_3_HepG2 = Tuni_3_sl_tmm/Ctrl_6_sl_tmm
  ) |> 
  select(UniprotID, ends_with("HepG2")) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  mutate(
    ratio_1_HEK293T = Tuni_1_sl_tmm/Ctrl_4_sl_tmm,
    ratio_2_HEK293T = Tuni_2_sl_tmm/Ctrl_5_sl_tmm,
    ratio_3_HEK293T = Tuni_3_sl_tmm/Ctrl_6_sl_tmm
  ) |> 
  select(UniprotID, starts_with("ratio")) |> 
  left_join(OG_glycoprotein_raw_sl_tmm_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  mutate(
    ratio_1_Jurkat = Tuni_1_sl_tmm/Ctrl_4_sl_tmm,
    ratio_2_Jurkat = Tuni_2_sl_tmm/Ctrl_5_sl_tmm,
    ratio_3_Jurkat = Tuni_3_sl_tmm/Ctrl_6_sl_tmm
  ) |> 
  select(UniprotID, starts_with("ratio"))

OG_glycorprotein_common_cor_matrix <- cor(OG_glycoprotein_raw_sl_tmm_common |> select(!UniprotID), method = "pearson")
rownames(OG_glycorprotein_common_cor_matrix) <- c(
  "HepG2 rep1", 
  "HepG2 rep2",
  "HepG2 rep3",
  "HEK293T rep1",
  "HEK293T rep2",
  "HEK293T rep3",
  "Jurkat rep1",
  "Jurkat rep2",
  "Jurkat rep3")

colnames(OG_glycorprotein_common_cor_matrix) <- c(
  "HepG2 rep1", 
  "HepG2 rep2",
  "HepG2 rep3",
  "HEK293T rep1",
  "HEK293T rep2",
  "HEK293T rep3",
  "Jurkat rep1",
  "Jurkat rep2",
  "Jurkat rep3")

corrplot(OG_glycorprotein_common_cor_matrix, type = "lower", method = "square",
         tl.cex = 1.2, tl.col = "black", tl.srt = 90,
         col = COL1('YlGn', n = 10),
         col.lim = c(0, 1),
         cl.cex = 0.8,
         is.corr = FALSE)

#Common terms
#RNA binding
gene_ontology_OG_glycoprotein_MF_RNA_binding_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_MF_RNA_binding.xlsx",
  sheet = 'HepG2',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UNIPROT_ACCESSION) |> 
  mutate(Cell = "HepG2") |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = join_by(UNIPROT_ACCESSION == UniprotID)) |> 
  select(UNIPROT_ACCESSION, logFC, Cell)

gene_ontology_OG_glycoprotein_MF_RNA_binding_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_MF_RNA_binding.xlsx",
  sheet = 'HEK293T',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UNIPROT_ACCESSION) |> 
  mutate(Cell = "HEK293T") |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = join_by(UNIPROT_ACCESSION == UniprotID)) |> 
  select(UNIPROT_ACCESSION, logFC, Cell)

gene_ontology_OG_glycoprotein_MF_RNA_binding_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_MF_RNA_binding.xlsx",
  sheet = 'Jurkat',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UNIPROT_ACCESSION) |> 
  mutate(Cell = "Jurkat") |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = join_by(UNIPROT_ACCESSION == UniprotID)) |> 
  select(UNIPROT_ACCESSION, logFC, Cell)

OG_glycoprotein_MF_RNA_binding_combined <- bind_rows(
  gene_ontology_OG_glycoprotein_MF_RNA_binding_HepG2,
  gene_ontology_OG_glycoprotein_MF_RNA_binding_HEK293T,
  gene_ontology_OG_glycoprotein_MF_RNA_binding_Jurkat
) |> mutate(Term = "RNA binding")

#wilcox test
OG_glycoprotein_MF_RNA_binding_wilcox_test <- OG_glycoprotein_MF_RNA_binding_combined |> 
  wilcox_test(logFC ~ Cell, p.adjust.method = "BH") |> 
  add_significance("p")

#ks test RNA binding
ks.test(gene_ontology_OG_glycoprotein_MF_RNA_binding_HepG2$logFC, gene_ontology_OG_glycoprotein_MF_RNA_binding_HEK293T$logFC)
ks.test(gene_ontology_OG_glycoprotein_MF_RNA_binding_HEK293T$logFC, gene_ontology_OG_glycoprotein_MF_RNA_binding_Jurkat$logFC)
ks.test(gene_ontology_OG_glycoprotein_MF_RNA_binding_Jurkat$logFC, gene_ontology_OG_glycoprotein_MF_RNA_binding_HepG2$logFC)

gene_ontology_OG_glycoprotein_MF_RNA_binding_ks_test <- tribble(
  ~ .y., ~ group1, ~ group2, ~ p, ~ p.signif,
  "logFC", "HEK293T", "HepG2", 0.004683, "**",
  "logFC", "Jurkat", "HEK293T", 0.04104, "*",
  "logFC", "HepG2", "Jurkat", 0.000231, "***"
)

#violin boxplot RNA binding
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_gene_ontology_OG_glycoprotein_MF_RNA_binding <- gene_ontology_OG_glycoprotein_MF_RNA_binding_combined |>
  ggplot() +
  geom_violin(aes(x = Cell, y = logFC, fill = Cell), color = "transparent") +
  geom_boxplot(aes(x = Cell, y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  scale_fill_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  labs(x = "", y = expression(log[2](Tuni/Ctrl))) +
  stat_pvalue_manual(data = gene_ontology_OG_glycoprotein_MF_RNA_binding_ks_test, label = "p.signif", tip.length = 0, size = 30,
                     y.position = c(2.5, 2.75, 3.0)) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_gene_ontology_OG_glycoprotein_MF_RNA_binding.png"), plot = violin_boxplot_gene_ontology_OG_glycoprotein_MF_RNA_binding, height = 4, width = 5, units = "in", dpi = 600)

#DNA binding
gene_ontology_OG_glycoprotein_MF_DNA_binding_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_MF_DNA_binding.xlsx",
  sheet = 'HepG2',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UNIPROT_ACCESSION) |> 
  mutate(Cell = "HepG2") |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = join_by(UNIPROT_ACCESSION == UniprotID)) |> 
  select(UNIPROT_ACCESSION, logFC, Cell)

gene_ontology_OG_glycoprotein_MF_DNA_binding_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_MF_DNA_binding.xlsx",
  sheet = 'HEK293T',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UNIPROT_ACCESSION) |> 
  mutate(Cell = "HEK293T") |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = join_by(UNIPROT_ACCESSION == UniprotID)) |> 
  select(UNIPROT_ACCESSION, logFC, Cell)

gene_ontology_OG_glycoprotein_MF_DNA_binding_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_MF_DNA_binding.xlsx",
  sheet = 'Jurkat',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UNIPROT_ACCESSION) |> 
  mutate(Cell = "Jurkat") |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = join_by(UNIPROT_ACCESSION == UniprotID)) |> 
  select(UNIPROT_ACCESSION, logFC, Cell)

OG_glycoprotein_MF_DNA_binding_combined <- bind_rows(
  gene_ontology_OG_glycoprotein_MF_DNA_binding_HepG2,
  gene_ontology_OG_glycoprotein_MF_DNA_binding_HEK293T,
  gene_ontology_OG_glycoprotein_MF_DNA_binding_Jurkat
) |> mutate(Term = "DNA binding")

#wilcox test
OG_glycoprotein_MF_DNA_binding_wilcox_test <- OG_glycoprotein_MF_DNA_binding_combined |> 
  wilcox_test(logFC ~ Cell, p.adjust.method = "BH") |> 
  add_significance("p")

#ks test DNA binding
ks.test(gene_ontology_OG_glycoprotein_MF_DNA_binding_HepG2$logFC, gene_ontology_OG_glycoprotein_MF_DNA_binding_HEK293T$logFC)
ks.test(gene_ontology_OG_glycoprotein_MF_DNA_binding_HEK293T$logFC, gene_ontology_OG_glycoprotein_MF_DNA_binding_Jurkat$logFC)
ks.test(gene_ontology_OG_glycoprotein_MF_DNA_binding_Jurkat$logFC, gene_ontology_OG_glycoprotein_MF_DNA_binding_HepG2$logFC)

gene_ontology_OG_glycoprotein_MF_DNA_binding_ks_test <- tribble(
  ~ .y., ~ group1, ~ group2, ~ p, ~ p.signif,
  "logFC", "HEK293T", "HepG2", 0.000261, "***",
  "logFC", "HepG2", "Jurkat", 0.01489, "*"
)

#violin boxplot DNA binding
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_gene_ontology_OG_glycoprotein_MF_DNA_binding <- gene_ontology_OG_glycoprotein_MF_DNA_binding_combined |>
  ggplot() +
  geom_violin(aes(x = Cell, y = logFC, fill = Cell), color = "transparent") +
  geom_boxplot(aes(x = Cell, y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  scale_fill_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  labs(x = "", y = expression(log[2](Tuni/Ctrl))) +
  stat_pvalue_manual(data = gene_ontology_OG_glycoprotein_MF_DNA_binding_ks_test, label = "p.signif", tip.length = 0, size = 30, hide.ns = "p.signif",
                     y.position = c(2.4, 2.7)) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_gene_ontology_OG_glycoprotein_MF_DNA_binding.png"), plot = violin_boxplot_gene_ontology_OG_glycoprotein_MF_DNA_binding, height = 4, width = 5, units = "in", dpi = 600)

#nucleoplasm
gene_ontology_OG_glycoprotein_CC_nucleoplasm_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_CC_nucleoplasm.xlsx",
  sheet = 'HepG2',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UNIPROT_ACCESSION) |> 
  mutate(Cell = "HepG2") |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = join_by(UNIPROT_ACCESSION == UniprotID)) |> 
  select(UNIPROT_ACCESSION, logFC, Cell)

gene_ontology_OG_glycoprotein_CC_nucleoplasm_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_CC_nucleoplasm.xlsx",
  sheet = 'HEK293T',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UNIPROT_ACCESSION) |> 
  mutate(Cell = "HEK293T") |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = join_by(UNIPROT_ACCESSION == UniprotID)) |> 
  select(UNIPROT_ACCESSION, logFC, Cell)

gene_ontology_OG_glycoprotein_CC_nucleoplasm_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_CC_nucleoplasm.xlsx",
  sheet = 'Jurkat',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UNIPROT_ACCESSION) |> 
  mutate(Cell = "Jurkat") |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = join_by(UNIPROT_ACCESSION == UniprotID)) |> 
  select(UNIPROT_ACCESSION, logFC, Cell)

OG_glycoprotein_CC_nucleoplasm_combined <- bind_rows(
  gene_ontology_OG_glycoprotein_CC_nucleoplasm_HepG2,
  gene_ontology_OG_glycoprotein_CC_nucleoplasm_HEK293T,
  gene_ontology_OG_glycoprotein_CC_nucleoplasm_Jurkat
) |> mutate(Term = "Nucleoplasm")

#wilcox test
OG_glycoprotein_CC_nucleoplasm_wilcox_test <- OG_glycoprotein_CC_nucleoplasm_combined |> 
  wilcox_test(logFC ~ Cell, p.adjust.method = "BH") |> 
  add_significance("p")

#ks test
ks.test(gene_ontology_OG_glycoprotein_CC_nucleoplasm_HepG2$logFC, gene_ontology_OG_glycoprotein_CC_nucleoplasm_HEK293T$logFC)
ks.test(gene_ontology_OG_glycoprotein_CC_nucleoplasm_HEK293T$logFC, gene_ontology_OG_glycoprotein_CC_nucleoplasm_Jurkat$logFC)
ks.test(gene_ontology_OG_glycoprotein_CC_nucleoplasm_Jurkat$logFC, gene_ontology_OG_glycoprotein_CC_nucleoplasm_HepG2$logFC)

gene_ontology_OG_glycoprotein_CC_nucleoplasm_ks_test <- tribble(
  ~ .y., ~ group1, ~ group2, ~ p, ~ p.signif,
  "logFC", "HEK293T", "HepG2", 9.636e-07, "****",
  "logFC", "Jurkat", "HEK293T", 9.279e-05, "****",
  "logFC", "HepG2", "Jurkat", 0.0002917, "***"
)

#violin boxplot nucleoplasm
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_gene_ontology_OG_glycoprotein_CC_nucleoplasm <- gene_ontology_OG_glycoprotein_CC_nucleoplasm_combined |>
  ggplot() +
  geom_violin(aes(x = Cell, y = logFC, fill = Cell), color = "transparent") +
  geom_boxplot(aes(x = Cell, y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  scale_fill_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  labs(x = "", y = expression(log[2](Tuni/Ctrl))) +
  stat_pvalue_manual(data = gene_ontology_OG_glycoprotein_CC_nucleoplasm_ks_test, label = "p.signif", tip.length = 0, size = 30, hide.ns = "p.signif",
                     y.position = c(2.4, 2.7, 3.0)) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_gene_ontology_OG_glycoprotein_CC_nucleoplasm_total.png"), plot = violin_boxplot_gene_ontology_OG_glycoprotein_CC_nucleoplasm_total, height = 4, width = 5, units = "in", dpi = 600)

#combine RNA binding, DNA binding and nucleoplasm
OG_glycoprotein_common_term_combined <- bind_rows(
  OG_glycoprotein_MF_RNA_binding_combined,
  OG_glycoprotein_MF_DNA_binding_combined,
  OG_glycoprotein_CC_nucleoplasm_combined
)

OG_glycoprotein_ks_test <- tribble(
  ~ .y., ~ group1, ~ group2, ~ p, ~ p.signif, ~ Term,
  "logFC", "HEK293T", "HepG2", 0.004683, "**", "RNA binding",
  "logFC", "Jurkat", "HEK293T", 0.04104, "*", "RNA binding",
  "logFC", "HepG2", "Jurkat", 0.000231, "***", "RNA binding",
  "logFC", "HEK293T", "HepG2", 0.000261, "***", "DNA binding",
  "logFC", "HepG2", "Jurkat", 0.01489, "*", "DNA binding",
  "logFC", "HEK293T", "HepG2", 9.636e-07, "****", "Nucleoplasm",
  "logFC", "Jurkat", "HEK293T", 9.279e-05, "****", "Nucleoplasm",
  "logFC", "HepG2", "Jurkat", 0.0002917, "***", "Nucleoplasm"
)

#violin boxplot for all three terms
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_common_term <- OG_glycoprotein_common_term_combined |>
  ggplot() +
  geom_violin(aes(x = Cell, y = logFC, fill = Cell), color = "transparent") +
  geom_boxplot(aes(x = Cell, y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  scale_fill_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  labs(x = "", y = expression(log[2](Tuni/Ctrl))) +
  stat_pvalue_manual(data = OG_glycoprotein_ks_test, label = "p.signif", tip.length = 0, 
                     size = 40, hide.ns = "p.signif", y.position = c(1.5, 1.9, 1.7, 1.5)) +
  facet_grid(~ Term, scale = "free_x") +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black"),
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_common_term.png"), plot = violin_boxplot_OG_glycoprotein_common_term, height = 4, width = 8, units = "in", dpi = 600)
