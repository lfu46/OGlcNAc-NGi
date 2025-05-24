#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "showtext", "ggpubr", "rstatix", "ComplexHeatmap", "circlize", "eulerr")
lapply(packages_names, require, character.only = TRUE)

#import database
G3BP1_interactor <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\G3BP1_interactor\\G3BP1_interactor_database.xlsx",
  col_name = TRUE,
  .name_repair = "universal"
)
uniprot_G3BP1_interactor <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\G3BP1_interactor\\uniprot_G3BP1.xlsx",
  col_name = TRUE,
  .name_repair = "universal"
)

G3BP1_interactor_uniprot <- G3BP1_interactor |> 
  left_join(uniprot_G3BP1_interactor, by = join_by(Gene.Symbol == From)) |> 
  select(UniprotID = Entry, Gene.Symbol, Type)

#generate data frame
#HEK293T
OG_glycoprotein_G3BP1_interactor_HEK293T <- G3BP1_interactor_uniprot |> 
  left_join(OG_glycoprotein_RBP_HEK293T, by = "UniprotID") |> 
  filter(!is.na(logFC))

WP_protein_G3BP1_interactor_HEK293T <- WP_protein_Top_tb_HEK293T |> 
  filter(UniprotID %in% OG_glycoprotein_G3BP1_interactor_HEK293T$UniprotID) |> 
  left_join(G3BP1_interactor_uniprot, by = "UniprotID") |> 
  mutate(cell = "HEK293T", group = "WP_RBP", group2 = "WP")

#HepG2
OG_glycoprotein_G3BP1_interactor_HepG2 <- G3BP1_interactor_uniprot |> 
  left_join(OG_glycoprotein_RBP_HepG2, by = "UniprotID") |> 
  filter(!is.na(logFC))

WP_protein_G3BP1_interactor_HepG2 <- WP_protein_Top_tb_HepG2 |> 
  filter(UniprotID %in% OG_glycoprotein_G3BP1_interactor_HepG2$UniprotID) |> 
  left_join(G3BP1_interactor_uniprot, by = "UniprotID") |> 
  mutate(cell = "HepG2", group = "WP_RBP", group2 = "WP")

#Jurkat
OG_glycoprotein_G3BP1_interactor_Jurkat <- G3BP1_interactor_uniprot |> 
  left_join(OG_glycoprotein_RBP_Jurkat, by = "UniprotID") |> 
  filter(!is.na(logFC))

WP_protein_G3BP1_interactor_Jurkat <- WP_protein_Top_tb_Jurkat |> 
  filter(UniprotID %in% OG_glycoprotein_G3BP1_interactor_Jurkat$UniprotID) |> 
  left_join(G3BP1_interactor_uniprot, by = "UniprotID") |> 
  mutate(cell = "Jurkat", group = "WP_RBP", group2 = "WP")

#combine
OG_WP_G3BP1_interactor_combined <- bind_rows(
  OG_glycoprotein_G3BP1_interactor_HEK293T,
  WP_protein_G3BP1_interactor_HEK293T,
  OG_glycoprotein_G3BP1_interactor_HepG2,
  WP_protein_G3BP1_interactor_HepG2,
  OG_glycoprotein_G3BP1_interactor_Jurkat,
  WP_protein_G3BP1_interactor_Jurkat
) |> mutate(group2 = ifelse(is.na(group2), "OG", group2))

write_xlsx(OG_WP_G3BP1_interactor_combined, path = paste0(file_path, "OG_WP_G3BP1_interactor_combined.xlsx"))

#wilcox test
OG_WP_G3BP1_interactor_wilcox_test <- OG_WP_G3BP1_interactor_combined |> 
  group_by(cell, Type) |> 
  wilcox_test(logFC ~ group2, p.adjust.method = "BH") |> 
  add_significance("p")

#point boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

point_boxplot_OG_WP_G3BP1_interactor <- OG_WP_G3BP1_interactor_combined |> 
  ggplot(aes(x = factor(group2, levels = c("WP", "OG")), y = logFC)) +
  geom_line(aes(group = UniprotID), color = "gray") +
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
  stat_pvalue_manual(data = OG_WP_G3BP1_interactor_wilcox_test, label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.2)) +
  coord_cartesian(ylim = c(-2, 2)) +
  facet_grid(cell ~ Type, scales = "free", labeller = labeller(
    "stress-dependent" = "stress\ndependent",
    "stress-independent" = "stress\nindependent"
  )) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text.x = element_text(size = 50, color = "black"),
    strip.text.y = element_text(size = 70, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "point_boxplot_OG_WP_G3BP1_interactor.png"), plot = point_boxplot_OG_WP_G3BP1_interactor, height = 4, width = 3, units = "in", dpi = 600)

#select proteins
selected_proteins <- c("Q5T6F2", "Q8TB72", "Q14157")

#generate data frame
selected_proteins_HEK293T <- OG_glycoprotein_RBP_HEK293T |> 
  filter(UniprotID %in% selected_proteins)
selected_proteins_HepG2 <- OG_glycoprotein_RBP_HepG2 |> 
  filter(UniprotID %in% selected_proteins)
selected_proteins_Jurkat <- OG_glycoprotein_RBP_Jurkat |> 
  filter(UniprotID %in% selected_proteins)

#combine
selected_proteins_combined <- bind_rows(
  selected_proteins_HEK293T,
  selected_proteins_HepG2,
  selected_proteins_Jurkat
) |> mutate(Gene_Names = c("UBAP2", "UBAP2L", "PUM2", "UBAP2", "UBAP2L", "UBAP2", "PUM2", "UBAP2L")) |> 
  select(cell, logFC, Gene_Names) |> 
  pivot_wider(names_from = Gene_Names, values_from = logFC)

#heatmap
mat <- data.matrix(selected_proteins_combined |> select(!cell))
rownames(mat) <- selected_proteins_combined$cell

mat <- t(mat)

col_mat <- colorRamp2(breaks = c(-1, 0, 1), colors = c("blue", "white", "red"))

h1 <- Heatmap(matrix = mat, col = col_mat, name = 'log2(Tuni/Ctrl)')
