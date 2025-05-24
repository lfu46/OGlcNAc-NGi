#import packages
packages_names <- c("tidyverse", "writexl", "readxl", "ComplexHeatmap", "circlize", 
                    "ggpubr", "rstatix", "showtext", "ggbeeswarm")
lapply(packages_names, require, character.only = TRUE)

#Nucleus, Extracellular Exosome, Mitochondrion
#impot data
Nuc_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\Nuc_Exo_Mito_HepG2\\Nuc_Exo_Mito.xlsx",
  sheet = 'Nucleus',
  col_names = TRUE,
  .name_repair = "universal"
) |> left_join(OG_glycoprotein_Top_tb_HepG2, by = join_by(UNIPROT_ACCESSION == UniprotID)) |> 
  select(UniprotID = UNIPROT_ACCESSION, logFC) |> 
  mutate(CC = 'Nucleus')

Exo_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\Nuc_Exo_Mito_HepG2\\Nuc_Exo_Mito.xlsx",
  sheet = 'Exosome',
  col_names = TRUE,
  .name_repair = "universal"
) |> left_join(OG_glycoprotein_Top_tb_HepG2, by = join_by(UNIPROT_ACCESSION == UniprotID)) |> 
  select(UniprotID = UNIPROT_ACCESSION, logFC) |> 
  mutate(CC = 'Exosome')

Mito_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\Nuc_Exo_Mito_HepG2\\Nuc_Exo_Mito.xlsx",
  sheet = 'Mitochondrion',
  col_names = TRUE,
  .name_repair = "universal"
) |> left_join(OG_glycoprotein_Top_tb_HepG2, by = join_by(UNIPROT_ACCESSION == UniprotID)) |> 
  select(UniprotID = UNIPROT_ACCESSION, logFC) |> 
  mutate(CC = 'Mitochondrion')

#combine
Nuc_Exo_Mito_combined <- bind_rows(
  Nuc_HepG2, Exo_HepG2, Mito_HepG2
)

#wilcox test
Nuc_Exo_Mito_wilcox_test <- Nuc_Exo_Mito_combined |> 
  wilcox_test(logFC ~ CC, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

violin_boxplot_nuc_exo_mito <- Nuc_Exo_Mito_combined |> 
  ggplot() +
  geom_violin(aes(x = factor(CC, levels = c("Exosome", "Nucleus", "Mitochondrion")), 
                  y = logFC, fill = factor(CC, levels = c("Exosome", "Nucleus", "Mitochondrion"))), 
              color = "transparent") +
  geom_boxplot(aes(x = factor(CC, levels = c("Exosome", "Nucleus", "Mitochondrion")), y = logFC), 
               color = "black", outliers = FALSE, width = 0.2) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(
    values = c(
      "Exosome" = Color_2,
      "Nucleus" = Color_3,
      "Mitochondrion" = Color_4
    )
  ) +
  stat_pvalue_manual(data = Nuc_Exo_Mito_wilcox_test, 
                     label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(2.0, 2.2)) +
  coord_cartesian(ylim = c(-1.8, 2.5)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 100, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_nuc_exo_mito.png"), 
       plot = violin_boxplot_nuc_exo_mito, height = 4, width = 3, units = "in", dpi = 600)

#Exosome functions
#Signaling
Signaling_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_HEK293T.xlsx",
  sheet = 'Signaling',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  mutate(cell = "HEK293T", category = "Signaling")

Signaling_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_HepG2.xlsx",
  sheet = 'Signaling',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  mutate(cell = "HepG2", category = "Signaling")

Signaling_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_Jurkat.xlsx",
  sheet = 'Signaling',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  mutate(cell = "Jurkat", category = "Signaling")

Signaling_combined <- bind_rows(
  Signaling_HEK293T,
  Signaling_HepG2,
  Signaling_Jurkat
)

#GTPase binding
GTP_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_HEK293T.xlsx",
  sheet = 'GTPase binding',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  mutate(cell = "HEK293T", category = "GTPase binding")

GTP_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_HepG2.xlsx",
  sheet = 'GTPase binding',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  mutate(cell = "HepG2", category = "GTPase binding")

GTP_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_Jurkat.xlsx",
  sheet = 'GTPase binding',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  mutate(cell = "Jurkat", category = "GTPase binding")

#ATP binding
ATP_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_HEK293T.xlsx",
  sheet = 'ATP binding',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  mutate(cell = "HEK293T", category = "ATP binding")

ATP_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_HepG2.xlsx",
  sheet = 'ATP binding',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  mutate(cell = "HepG2", category = "ATP binding")

ATP_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_Jurkat.xlsx",
  sheet = 'ATP binding',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  mutate(cell = "Jurkat", category = "ATP binding")

#protein binding
protein_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_HEK293T.xlsx",
  sheet = 'Protein binding',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  mutate(cell = "HEK293T", category = "Protein binding")

protein_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_HepG2.xlsx",
  sheet = 'Protein binding',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  mutate(cell = "HepG2", category = "Protein binding")

protein_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_Jurkat.xlsx",
  sheet = 'Protein binding',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  mutate(cell = "Jurkat", category = "Protein binding")

#combine
Exosome_cluster_combined <- bind_rows(
  Signaling_HEK293T,
  Signaling_HepG2,
  Signaling_Jurkat,
  ATP_HEK293T,
  ATP_HepG2,
  ATP_Jurkat,
  protein_HEK293T,
  protein_HepG2,
  protein_Jurkat,
  GTP_HEK293T,
  GTP_HepG2,
  GTP_Jurkat
)

write_xlsx(Exosome_cluster_combined, path = paste0(file_path, "Exosome_cluster_combined.xlsx"))

#wilcox test
Exosome_cluster_wilcox_test <- Exosome_cluster_combined |> 
  mutate(category = factor(category, 
                           levels = c("Signaling", "Protein binding", 
                                      "GTPase binding", "ATP binding"))) |> 
  group_by(category) |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_exosome_cluster <- Exosome_cluster_combined |> 
  mutate(category = factor(category, 
                           levels = c("Signaling", "Protein binding", 
                                      "GTPase binding", "ATP binding"))) |> 
  ggplot() +
  geom_violin(aes(x = cell, y = logFC, fill = cell), color = "transparent") +
  geom_boxplot(aes(x = cell, y = logFC), color = "black", outliers = FALSE, width = 0.2) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_wrap(vars(category), scales = "free_x", nrow = 1) +
  stat_pvalue_manual(data = Exosome_cluster_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.6, 1.9, 1.6, 1.9, 1.6, 1.9, 1.6)) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 90, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_exosome_cluster.png"), 
       plot = violin_boxplot_exosome_cluster, height = 4, width = 7, units = "in", dpi = 600)

#Mito. Central Dogma & Disease & Sub-mito
#import database
HumanMitoCarta_3 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\mito\\Human.MitoCarta3.xlsx",
  sheet = 'A Human MitoCarta3.0',
  col_names = TRUE,
  .name_repair = "universal"
) 

Mito_localization <- HumanMitoCarta_3 |> 
  select(Symbol, UniprotID = UniProt, MitoCarta3.0_SubMitoLocalization)

Mito_pathway <- HumanMitoCarta_3 |> 
  select(Symbol, UniprotID = UniProt, MitoCarta3.0_MitoPathways) |> 
  separate_rows(MitoCarta3.0_MitoPathways, sep = ' \\| ') |> 
  distinct()

#MitoCop
MitoCop_all_protein_group <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\mito\\MitoCoP.xlsx",
  sheet = '(A) All protein groups',
  skip = 1,
  col_names = TRUE,
  .name_repair = "universal"
)

MitoCop_all_protein_group_adj <- MitoCop_all_protein_group |> 
  select(UniprotID = Simplified.protein.IDs, MitoCoP, 
         substractive = Classification...32, spatial = Classification...41,
         importomics = Classification...66, membrane_association = Category, 
         mito_copy = Mean.mito.copies.per.cell...2.3.Reps.,
         cell_copy = Mean.copy.numbers.per.cell...2.3.Reps., half_lives = Protein.half.lives..h.) |> 
  filter(MitoCoP == 1) |> 
  separate_rows(UniprotID, sep = ";") |> 
  distinct()

MitoCop_all_protein_group_adj_2 <- MitoCop_all_protein_group |> 
  select(UniprotID = Simplified.protein.IDs, MitoCoP,
         'Mitochondrial.outer.membrane..GO.0005741...all.evidence.sources.':'Fe.S.protein.biogenesis') |> 
  filter(MitoCoP == 1) |> 
  separate_rows(UniprotID, sep = ";") |> 
  pivot_longer('Mitochondrial.outer.membrane..GO.0005741...all.evidence.sources.':'Fe.S.protein.biogenesis',
               names_to = 'submitochondrial', values_to = 'value') |> 
  filter(!is.na(value))

MitoCop_mito <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\mito\\MitoCoP.xlsx",
  sheet = '(B) MitoCoP (1,134 genes)',
  skip = 1,
  col_names = TRUE,
  .name_repair = "universal"
)

MitoCop_mito_function <- MitoCop_mito |> 
  select(UniprotID = Simplified.protein.IDs, MitoCoP, 
         Morphology..dynamics...organization:Unknown) |> 
  filter(MitoCoP == 1) |> 
  separate_rows(UniprotID, sep = ";") |> 
  pivot_longer(Morphology..dynamics...organization:Unknown, 
               names_to = 'Function', values_to = 'value') |> 
  filter(!is.na(value))

MitoCop_mito_disease <- MitoCop_mito |> 
  select(UniprotID = Simplified.protein.IDs, MitoCoP.disease.gene, 
         Central.nervous.system:Tumors) |> 
  filter(MitoCoP.disease.gene == 1) |> 
  separate_rows(UniprotID, sep = ";") |> 
  pivot_longer(Central.nervous.system:Tumors, 
               names_to = 'Function', values_to = 'value') |> 
  filter(!is.na(value)) |> 
  distinct()

#generate data frame
#HEK293T
OG_glycoprotein_mito_1_HEK293T <- MitoCop_all_protein_group_adj |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HEK293T")

OG_glycoprotein_mito_2_HEK293T <- MitoCop_all_protein_group_adj_2 |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HEK293T")

OG_glycoprotein_mito_3_HEK293T <- MitoCop_mito_disease |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HEK293T")

#HepG2
OG_glycoprotein_mito_1_HepG2 <- MitoCop_all_protein_group_adj |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HepG2")

OG_glycoprotein_mito_2_HepG2 <- MitoCop_all_protein_group_adj_2 |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HepG2")

OG_glycoprotein_mito_3_HepG2 <- MitoCop_mito_disease |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HepG2")

#Jurkat
OG_glycoprotein_mito_1_Jurkat <- MitoCop_all_protein_group_adj |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "Jurkat")

OG_glycoprotein_mito_2_Jurkat <- MitoCop_all_protein_group_adj_2 |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "Jurkat")

OG_glycoprotein_mito_3_Jurkat <- MitoCop_mito_disease |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "Jurkat")

#combined
OG_glycoprotein_mitocop_mito_1_combined <- bind_rows(
  OG_glycoprotein_mito_1_HEK293T,
  OG_glycoprotein_mito_1_HepG2,
  OG_glycoprotein_mito_1_Jurkat
)

write_xlsx(OG_glycoprotein_mitocop_mito_1_combined, 
           path = paste0(file_path, "OG_glycoprotein_mitocop_mito_1_combined.xlsx"))

OG_glycoprotein_mitocop_mito_2_combined <- bind_rows(
  OG_glycoprotein_mito_2_HEK293T,
  OG_glycoprotein_mito_2_HepG2,
  OG_glycoprotein_mito_2_Jurkat
)

write_xlsx(OG_glycoprotein_mitocop_mito_2_combined,
           path = paste0(file_path, "OG_glycoprotein_mitocop_mito_2_combined.xlsx"))

OG_glycoprotein_mitocop_mito_3_combined <- bind_rows(
  OG_glycoprotein_mito_3_HEK293T,
  OG_glycoprotein_mito_3_HepG2,
  OG_glycoprotein_mito_3_Jurkat
)

write_xlsx(OG_glycoprotein_mitocop_mito_3_combined,
           path = paste0(file_path, "OG_glycoprotein_mitocop_mito_3_combined.xlsx"))

#sub mito 
#wilcox test
OG_glycoprotein_mitocop_submito_wilcox_test <- OG_glycoprotein_mitocop_mito_2_combined |> 
  mutate(group = ifelse(submitochondrial == "Mitochondrial.inner.membrane..GO.0005743...all.evidence.sources.", 'Mito.inner', NA),
         group = ifelse(submitochondrial == "Mitochondrial.outer.membrane..GO.0005741...all.evidence.sources.", 'Mito.outer', group),
         group = ifelse(submitochondrial == "Mitochondrial.matrix..GO.0005759...all.evidence.sources.", 'Mito.matrix', group), 
         group = ifelse(is.na(group), 'Other', group)) |> 
  group_by(group) |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

point_boxplot_OG_glycoprotein_mito_sub_mito <- OG_glycoprotein_mitocop_mito_2_combined |> 
  mutate(group = ifelse(submitochondrial == "Mitochondrial.inner.membrane..GO.0005743...all.evidence.sources.", 'Mito.inner', NA),
         group = ifelse(submitochondrial == "Mitochondrial.outer.membrane..GO.0005741...all.evidence.sources.", 'Mito.outer', group),
         group = ifelse(submitochondrial == "Mitochondrial.matrix..GO.0005759...all.evidence.sources.", 'Mito.matrix', group), 
         group = ifelse(is.na(group), 'Other', group)) |> 
  ggplot(aes(x = cell, y = logFC)) +
  geom_jitter(aes(color = cell), position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", outliers = FALSE, fill = "transparent") +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_color_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_grid(~ factor(group, levels = c("Mito.outer", "Mito.inner", "Mito.matrix", "Other")), scales = "free_x") +
  stat_pvalue_manual(data = OG_glycoprotein_mitocop_submito_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(2.2, 2.1, 2.1, 1.8, 2.3, 1.7, 2.0)) +
  coord_cartesian(ylim = c(-1.2, 2.5)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent")
  )

ggsave(filename = paste0(file_path, "point_boxplot_OG_glycoprotein_mito_sub_mito.png"), 
       plot = point_boxplot_OG_glycoprotein_mito_sub_mito, height = 3, width = 6, units = "in", dpi = 600)

#disease
#wilcox test
OG_glycoprotein_mitocop_disease_wilcox_test <- OG_glycoprotein_mitocop_mito_3_combined |> 
  filter(Function == "Metabolism" | Function == "Liver") |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#point boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

point_boxplot_OG_glycoprotein_mito_disease <- OG_glycoprotein_mitocop_mito_3_combined |> 
  ggplot(aes(x = cell, y = logFC)) +
  geom_jitter(aes(color = cell), position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", fill = "transparent", outliers = FALSE) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_color_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  stat_pvalue_manual(data = OG_glycoprotein_mitocop_disease_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.4, 1.7, 2.0)) +
  coord_cartesian(ylim = c(-1.8, 2.2)) +
  ggtitle("Liver & Metabolism") +
  theme_bw() +
  theme(
    title = element_text(color = "black", size = 70, lineheight = 0.15),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent")
  )

ggsave(filename = paste0(file_path, "point_boxplot_OG_glycoprotein_mito_disease.png"), 
       plot = point_boxplot_OG_glycoprotein_mito_disease, height = 3, width = 2.5, units = "in", dpi = 600)

#generate data frame
#HEK293T
OG_glycoprotein_mito_pathway_HEK293T <- Mito_pathway |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HEK293T")

#HepG2
OG_glycoprotein_mito_pathway_HepG2 <- Mito_pathway |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "HepG2")

#Jurkat
OG_glycoprotein_mito_pathway_Jurkat <- Mito_pathway |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = 'UniprotID') |> 
  filter(!is.na(logFC)) |> 
  mutate(cell = "Jurkat")

#combined
OG_glycoprotein_mito_pathway_combined <- bind_rows(
  OG_glycoprotein_mito_pathway_HEK293T,
  OG_glycoprotein_mito_pathway_HepG2,
  OG_glycoprotein_mito_pathway_Jurkat
)

write_xlsx(OG_glycoprotein_mito_pathway_combined, path = paste0(file_path, "OG_glycoprotein_mito_pathway_combined.xlsx"))

#pathway
#wilcox test
OG_glycoprotein_mito_pathway_wilcox_test <- OG_glycoprotein_mito_pathway_combined |> 
  filter(str_detect(MitoCarta3.0_MitoPathways, "Mitochondrial central dogma")) |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

point_boxplot_OG_glycoprotein_mito_pathway <- OG_glycoprotein_mito_pathway_combined |> 
  ggplot(aes(x = cell, y = logFC)) +
  geom_jitter(aes(color = cell), position = position_jitter(width = 0.2)) +
  geom_boxplot(color = "black", fill = "transparent", outliers = FALSE) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_color_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  stat_pvalue_manual(data = OG_glycoprotein_mito_pathway_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.7, 2.2, 1.8)) +
  coord_cartesian(ylim = c(-1.8, 2.4)) +
  ggtitle("Mito. Central Dogma") +
  theme_bw() +
  theme(
    title = element_text(color = "black", size = 70, lineheight = 0.15),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent")
  )

ggsave(filename = paste0(file_path, "point_boxplot_OG_glycoprotein_mito_pathway.png"), 
       plot = point_boxplot_OG_glycoprotein_mito_pathway, height = 3, width = 2.5, units = "in", dpi = 600)
