#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "showtext", "ggpubr", "rstatix", "ComplexHeatmap", "circlize", "eulerr")
lapply(packages_names, require, character.only = TRUE)

#RBP subcellular location
#import database
human_protein_atlas_subcellular <- read_tsv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_protein_atlas\\subcellular_location.tsv",
  col_names = TRUE, 
  name_repair = "universal"
) |> filter(Reliability != "Uncertain") |> 
  mutate(Main_Additional = paste(Main.location, Additional.location, sep = ";")) |> 
  select(Gene, Gene.name, Main_Additional)

#import data
OG_glycoprotein_RBP_subcellular_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_protein_atlas\\OG_glycoprotein_RBP_subcellular_HEK293T.xlsx"
) |> separate(To, into = c('Ensembl', 'v'), sep = "\\.")
OG_glycoprotein_RBP_subcellular_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_protein_atlas\\OG_glycoprotein_RBP_subcellular_HepG2.xlsx"
) |> separate(To, into = c('Ensembl', 'v'), sep = "\\.")
OG_glycoprotein_RBP_subcellular_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_protein_atlas\\OG_glycoprotein_RBP_subcellular_Jurkat.xlsx"
) |> separate(To, into = c('Ensembl', 'v'), sep = "\\.")

#generate data frame
#HEK293T
OG_glycoprotein_RBP_subcellular_logFC_HEK293T <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  left_join(OG_glycoprotein_RBP_subcellular_HEK293T, by = join_by("UniprotID" == "From")) |> 
  mutate(cell = "HEK293T") |> 
  select(UniprotID, logFC, cell, Ensembl) |> 
  left_join(human_protein_atlas_subcellular, by = join_by("Ensembl" == "Gene")) |> 
  filter(!is.na(Main_Additional)) |> 
  mutate(
    Nuc = ifelse(str_detect(Main_Additional, "Nuc"), 1, 0)
  ) |> 
  mutate(
    Cyto = ifelse(str_detect(Main_Additional, "Cyto"), 1, 0)
  ) |> 
  mutate(
    Nuc_Cyto = Nuc + Cyto
  ) |> 
  mutate(
    location = ifelse(Nuc_Cyto == 2, "Both", ifelse(Nuc == 1, "Nuc", "Cyto"))
  )

#HepG2
OG_glycoprotein_RBP_subcellular_logFC_HepG2 <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  left_join(OG_glycoprotein_RBP_subcellular_HepG2, by = join_by("UniprotID" == "From")) |> 
  mutate(cell = "HepG2") |> 
  select(UniprotID, logFC, cell, Ensembl) |> 
  left_join(human_protein_atlas_subcellular, by = join_by("Ensembl" == "Gene")) |> 
  filter(!is.na(Main_Additional)) |> 
  mutate(
    Nuc = ifelse(str_detect(Main_Additional, "Nuc"), 1, 0)
  ) |> 
  mutate(
    Cyto = ifelse(str_detect(Main_Additional, "Cyto"), 1, 0)
  ) |> 
  mutate(
    Nuc_Cyto = Nuc + Cyto
  ) |> 
  mutate(
    location = ifelse(Nuc_Cyto == 2, "Both", ifelse(Nuc == 1, "Nuc", "Cyto"))
  )

#Jurkat
OG_glycoprotein_RBP_subcellular_logFC_Jurkat <- OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  left_join(OG_glycoprotein_RBP_subcellular_Jurkat, by = join_by("UniprotID" == "From")) |> 
  mutate(cell = "Jurkat") |> 
  select(UniprotID, logFC, cell, Ensembl) |> 
  left_join(human_protein_atlas_subcellular, by = join_by("Ensembl" == "Gene")) |> 
  filter(!is.na(Main_Additional)) |> 
  mutate(
    Nuc = ifelse(str_detect(Main_Additional, "Nuc"), 1, 0)
  ) |> 
  mutate(
    Cyto = ifelse(str_detect(Main_Additional, "Cyto"), 1, 0)
  ) |> 
  mutate(
    Nuc_Cyto = Nuc + Cyto
  ) |> 
  mutate(
    location = ifelse(Nuc_Cyto == 2, "Both", ifelse(Nuc == 1, "Nuc", "Cyto"))
  )

#nuc cyto both proportion
#HEK293T
OG_glycoprotein_RBP_subcellular_logFC_HEK293T |> 
  count(location)

#HepG2
OG_glycoprotein_RBP_subcellular_logFC_HepG2 |> 
  count(location)

#Jurkat
OG_glycoprotein_RBP_subcellular_logFC_Jurkat |> 
  count(location)

#combine
OG_glycoprotein_RBP_subcellular_logFC_combined <- bind_rows(
  OG_glycoprotein_RBP_subcellular_logFC_HEK293T,
  OG_glycoprotein_RBP_subcellular_logFC_HepG2,
  OG_glycoprotein_RBP_subcellular_logFC_Jurkat
)

#bar plot
OG_glycoprotein_RBP_subcellular_logFC_combined_count <- OG_glycoprotein_RBP_subcellular_logFC_combined |> 
  count(cell, location)

font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

barplot_OG_glycoprotein_RBP_subcellular_count <- OG_glycoprotein_RBP_subcellular_logFC_combined_count |> 
  ggplot(aes(x = factor(location, levels = c("Nuc", "Cyto", "Both")), y = n)) +
  geom_bar(aes(fill = cell), color = "transparent", stat = "identity") +
  geom_text(aes(label = n), vjust = -1, size = 30, color = "black") +
  labs(x = "", y = "# of quantified RBP") +
  coord_cartesian(ylim = c(0, 180)) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_grid(~ cell, scales = "free_x") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "barplot_OG_glycoprotein_RBP_subcellular_count.png"), plot = barplot_OG_glycoprotein_RBP_subcellular_count, height = 4, width = 6, units = "in", dpi = 600)

#wilcox test
#HEK293T
OG_glycoprotein_RBP_subcellular_logFC_combined_wilcox_test <- OG_glycoprotein_RBP_subcellular_logFC_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ location, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_boxplot_OG_glycoprotein_RBP_subcellular_logFC <- OG_glycoprotein_RBP_subcellular_logFC_combined |> 
  ggplot() +
  geom_violin(aes(x = factor(location, levels = c("Nuc", "Both", "Cyto")), y = logFC, fill = cell),
              color = "transparent") +
  geom_boxplot(aes(x = factor(location, levels = c("Nuc", "Both", "Cyto")), y = logFC), outliers = FALSE,
               width = 0.2) +
  labs(x = "", y = expression(Δlog[2]*"FC")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_grid(~ cell, scales = "free_x") +
  stat_pvalue_manual(data = OG_glycoprotein_RBP_subcellular_logFC_combined_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.6, 1.9)) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "violin_boxplot_OG_glycoprotein_RBP_subcellular_logFC.png"), plot = violin_boxplot_OG_glycoprotein_RBP_subcellular_logFC, 
       height = 4, width = 6, units = "in", dpi = 600)
