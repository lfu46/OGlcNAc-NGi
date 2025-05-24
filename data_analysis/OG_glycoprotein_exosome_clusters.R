#import packages
packages_names <- c("tidyverse", "writexl", "readxl", "ComplexHeatmap", "circlize", 
                    "ggpubr", "rstatix", "showtext")
lapply(packages_names, require, character.only = TRUE)

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

write_xlsx(Signaling_combined, path = paste0(file_path, "Signaling_combined.xlsx"))

Rap1 <- c("P46108", "P04899", "P15153", "P63261", "P05107", "P50552")
AMPK <- c("Q53ET0", "P61019", "P13639", "P49327")
PI3K <- c("Q53ET0", "P15153", "P02751", "P14625", "P07942", "P11047", "P27348")

#Nucleus
Nucleus_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_HEK293T.xlsx",
  sheet = 'Nucleus',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  mutate(cell = "HEK293T", category = "Nucleus")

Nucleus_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_HepG2.xlsx",
  sheet = 'Nucleus',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  mutate(cell = "HepG2", category = "Nucleus")

Nucleus_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_Jurkat.xlsx",
  sheet = 'Nucleus',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  mutate(cell = "Jurkat", category = "Nucleus")

#Mitochondrion
Mitochondrion_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_HEK293T.xlsx",
  sheet = 'Mitochondrion',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  mutate(cell = "HEK293T", category = "Mitochondrion")

Mitochondrion_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_HepG2.xlsx",
  sheet = 'Mitochondrion',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  mutate(cell = "HepG2", category = "Mitochondrion")

Mitochondrion_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_Jurkat.xlsx",
  sheet = 'Signaling',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  mutate(cell = "Jurkat", category = "Mitochondrion")

#ER
ER_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_HEK293T.xlsx",
  sheet = 'ER',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  mutate(cell = "HEK293T", category = "ER")

ER_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_HepG2.xlsx",
  sheet = 'ER',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  mutate(cell = "HepG2", category = "ER")

ER_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\OG_glycoprotein_exosome_Jurkat.xlsx",
  sheet = 'ER',
  col_names = TRUE,
  .name_repair = "universal"
) |> select(UniprotID = UNIPROT_ACCESSION) |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  mutate(cell = "Jurkat", category = "ER")

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
font_add(family = "arial", regular = "arial.ttf")
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
                     tip.length = 0, size = 6,
                     hide.ns = "p", y.position = c(1.5, 1.8, 1.5, 1.8, 1.6, 1.8, 1.5)) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(size = 10, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 8, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent")
  )

ggsave(filename = paste0(file_path, "violin_boxplot_exosome_cluster.eps"), 
       device = "eps",
       plot = violin_boxplot_exosome_cluster, 
       height = 2, width = 4, 
       units = "in", 
       dpi = 1200
       )

#PI3K
#wilcox test
PI3K_wilcox_test <- Signaling_combined |> 
  filter(UniprotID %in% PI3K) |> 
  wilcox_test(logFC ~ cell, p.adjust.method = "BH") |> 
  add_significance("p")
