#import packages
packages_names <- c("tidyverse", "writexl", "readxl", "ComplexHeatmap", "circlize", 
                    "ggpubr", "rstatix", "showtext", "ggbeeswarm")
lapply(packages_names, require, character.only = TRUE)

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
  labs(x = "", y = "") +
  scale_fill_manual(
    values = c(
      "Exosome" = Color_2,
      "Nucleus" = Color_3,
      "Mitochondrion" = Color_4
    )
  ) +
  stat_pvalue_manual(data = Nuc_Exo_Mito_wilcox_test, 
                     label = "p.signif", tip.length = 0, size = 6,
                     hide.ns = "p", y.position = c(2.0, 2.2)) +
  coord_cartesian(ylim = c(-1.8, 2.5)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "violin_boxplot_nuc_exo_mito.eps"), 
       device = "eps",
       plot = violin_boxplot_nuc_exo_mito, 
       height = 2, width = 2, 
       units = "in", 
       dpi = 1200)
