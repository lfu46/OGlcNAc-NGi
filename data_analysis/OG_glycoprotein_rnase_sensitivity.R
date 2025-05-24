#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "showtext", "ggpubr", "rstatix", "ComplexHeatmap", "circlize", "eulerr")
lapply(packages_names, require, character.only = TRUE)

#import database
ribo_complex_database <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ribonucleoprotein_complex_database\\TableS3_final.xlsx",
  col_name = TRUE,
  .name_repair= "universal"
) |> select(complex_id, rna_associated_accs, rnp_behavior) |> 
  separate_rows(rna_associated_accs, sep = " ") |> 
  filter(!is.na(rnp_behavior))

write_xlsx(ribo_complex_database, path = paste0(file_path, "ribo_complex_database.xlsx"))

#generate data frame
WP_protein_ribo_complex_HEK293T <- WP_protein_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_RBP_HEK293T$UniprotID) |> 
  mutate(cell = "HEK293T", group = "WP_RBP", group2 = "WP")
WP_protein_ribo_complex_HepG2 <- WP_protein_Top_tb_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_RBP_HepG2$UniprotID) |> 
  mutate(cell = "HepG2", group = "WP_RBP", group2 = "WP")
WP_protein_ribo_complex_Jurkat <- WP_protein_Top_tb_Jurkat |> 
  filter(UniprotID %in% OG_glycoprotein_RBP_Jurkat$UniprotID) |> 
  mutate(cell = "Jurkat", group = "WP_RBP", group2 = "WP")


#generate data frame
#HEK293T
OG_glycoprotein_ribo_complex_HEK293T <- ribo_complex_database |> 
  left_join(OG_glycoprotein_RBP_HEK293T, by = join_by(rna_associated_accs == UniprotID)) |> 
  filter(!is.na(logFC))

WP_protein_ribo_complex_HEK293T <- ribo_complex_database |> 
  left_join(WP_protein_ribo_complex_HEK293T, by = join_by(rna_associated_accs == UniprotID)) |> 
  filter(!is.na(logFC))

OG_glycoprotein_ribo_complex_HEK293T |> 
  wilcox_test(logFC ~ rnp_behavior, p.adjust.method = "BH") |> 
  add_significance("p")

#HepG2
OG_glycoprotein_ribo_complex_HepG2 <- ribo_complex_database |> 
  left_join(OG_glycoprotein_RBP_HepG2, by = join_by(rna_associated_accs == UniprotID)) |> 
  filter(!is.na(logFC))

WP_protein_ribo_complex_HepG2 <- ribo_complex_database |> 
  left_join(WP_protein_ribo_complex_HepG2, by = join_by(rna_associated_accs == UniprotID)) |> 
  filter(!is.na(logFC))


OG_glycoprotein_ribo_complex_HepG2 |> 
  wilcox_test(logFC ~ rnp_behavior, p.adjust.method = "BH") |> 
  add_significance("p")

#Jurkat
OG_glycoprotein_ribo_complex_Jurkat <- ribo_complex_database |> 
  left_join(OG_glycoprotein_RBP_Jurkat, by = join_by(rna_associated_accs == UniprotID)) |> 
  filter(!is.na(logFC))

WP_protein_ribo_complex_Jurkat <- ribo_complex_database |> 
  left_join(WP_protein_ribo_complex_Jurkat, by = join_by(rna_associated_accs == UniprotID)) |> 
  filter(!is.na(logFC))

OG_glycoprotein_ribo_complex_Jurkat |> 
  wilcox_test(logFC ~ rnp_behavior, p.adjust.method = "BH") |> 
  add_significance("p")

#combine
OG_WP_ribo_complex <- bind_rows(
  OG_glycoprotein_ribo_complex_HEK293T,
  WP_protein_ribo_complex_HEK293T,
  OG_glycoprotein_ribo_complex_HepG2,
  WP_protein_ribo_complex_HepG2,
  OG_glycoprotein_ribo_complex_Jurkat,
  WP_protein_ribo_complex_Jurkat
) |> mutate(group2 = ifelse(is.na(group2), "OG", group2))

write_xlsx(OG_WP_ribo_complex, path = paste0(file_path, "OG_WP_ribo_complex.xlsx"))

#wilcox test
OG_WP_complex_wilcox_test <- OG_WP_ribo_complex |> 
  group_by(cell, rnp_behavior) |> 
  wilcox_test(logFC ~ group2, p.adjust.method = "BH") |> 
  add_significance("p") |> 
  filter(rnp_behavior != "COMPOSITIONAL")

#point boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

point_boxplot_OG_WP_ribo_complex <- OG_WP_ribo_complex |> 
  filter(rnp_behavior != "COMPOSITIONAL") |> 
  ggplot(aes(x = factor(group2, levels = c("WP", "OG")), y = logFC)) +
  geom_line(aes(group = rna_associated_accs), color = "gray") +
  geom_boxplot() +
  geom_point(aes(fill = group), shape = 21,
             color = "black") +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(
    values = c(
      "OG_RBP_HEK293T" = Color_2,
      "OG_RBP_HepG2" = Color_3,
      "OG_RBP_Jurkat" = Color_4,
      "WP_RBP" = "gray"
    )
  ) +
  stat_pvalue_manual(data = OG_WP_complex_wilcox_test, label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(2.5)) +
  coord_cartesian(ylim = c(-3, 3)) +
  facet_grid(cell ~ rnp_behavior, scales = "free") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 65, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "point_boxplot_OG_WP_ribo_complex.png"), plot = point_boxplot_OG_WP_ribo_complex, height = 4, width = 3, units = "in", dpi = 600)

#heatmap
uniprot_selected_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ribonucleoprotein_complex_database\\uniprot_selected_glycoprotein_HepG2.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)

#generate data frame
OG_glycoprotein_selected_HepG2 <- uniprot_selected_HepG2 |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = join_by(From == UniprotID)) |> 
  select(Gene.Names, UniprotID = From, logFC) |> 
  mutate(group = "OG")
WP_protein_selected_HepG2 <- uniprot_selected_HepG2 |> 
  left_join(WP_protein_Top_tb_HepG2, by = join_by(From == UniprotID)) |> 
  select(Gene.Names, UniprotID = From, logFC) |> 
  mutate(group = "WP")

#combine
OG_WP_selected_HepG2 <- bind_rows(
  OG_glycoprotein_selected_HepG2,
  WP_protein_selected_HepG2
) |> 
  select(Gene.Names, logFC, group) |> 
  pivot_wider(names_from = Gene.Names, values_from = logFC)

#heatmap
mat <- t(data.matrix(OG_WP_selected_HepG2 |> select(!group)))
colnames(mat) <- c("OG", "WP")
col_mat <- colorRamp2(breaks = c(-1, 0, 1.0), colors = c("blue", "white", "red"))

h1 <- Heatmap(matrix = mat, col = col_mat, name = 'log2(Tuni/Ctrl)')
