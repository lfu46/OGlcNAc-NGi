#import packages
packages_names <- c("showtext", "ggpubr", "rstatix")
lapply(packages_names, require, character.only = TRUE)

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
font_add(family = "arial", regular = "arial.ttf")
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
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  stat_pvalue_manual(data = OG_glycoprotein_ks_test, label = "p.signif", tip.length = 0, 
                     size = 10, hide.ns = "p.signif", y.position = c(1.5, 1.9, 1.7, 1.5)) +
  facet_grid(~ Term, scale = "free_x") +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(size = 15, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 15, color = "black"),
    legend.title = element_text(size = 15, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 20, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_common_term.eps"), 
       device = "eps",
       plot = violin_boxplot_OG_glycoprotein_common_term, 
       height = 4, width = 8, 
       units = "in", 
       dpi = 1200
       )
