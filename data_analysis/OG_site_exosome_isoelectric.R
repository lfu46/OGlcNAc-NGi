#import packages
packages_names <- c("tidyverse", "rstatix", "showtext", "ggpubr", 
                    "readxl", "writexl", "ComplexHeatmap", "circlize")
lapply(packages_names, require, character.only = TRUE)

#import uniprot sequence data
OG_glycoprotein_uniprot_sequence <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\uniprot_sequence\\OG_site_glycoprotein_sequence.tsv",
  delim = "\t",
  col_names = TRUE,
  name_repair = "universal"
)

#site level
#generate data frame for 13mer
#exosome
exosome_13mer_HepG2 <- OG_glycopeptide_Top_tb_HepG2_singlesite |> 
  filter(UniprotID %in% (OG_glycoprotein_exosome_nucleus_HepG2_combined |> 
           filter(category == "Exosome") |> pull(UniprotID))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence) |> 
  mutate(seq = '>peptide', category = "Exosome") |> 
  pivot_longer(cols = c('seq', 'Mer13'), names_to = 'Mer13_adj_1', values_to = 'Mer13_adj_2')

write_xlsx(exosome_13mer_HepG2, path = paste0(file_path, "exosome_13mer_HepG2.xlsx"))

#nucleus
nucleus_13mer_HepG2 <- OG_glycopeptide_Top_tb_HepG2_singlesite |> 
  filter(UniprotID %in% (OG_glycoprotein_exosome_nucleus_HepG2_combined |> 
                           filter(category == "Nucleus") |> pull(UniprotID))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence) |> 
  mutate(seq = '>peptide', category = "Nucleus") |> 
  pivot_longer(cols = c('seq', 'Mer13'), names_to = 'Mer13_adj_1', values_to = 'Mer13_adj_2')

write_xlsx(nucleus_13mer_HepG2, path = paste0(file_path, "nucleus_13mer_HepG2.xlsx"))

#exosome nucleus
exosome_nucleus_13mer_HepG2 <- OG_glycopeptide_Top_tb_HepG2_singlesite |> 
  filter(UniprotID %in% (OG_glycoprotein_exosome_nucleus_HepG2_combined |> 
                           filter(category == "Exosome, Nucleus") |> pull(UniprotID))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence) |> 
  mutate(seq = '>peptide', category = "Exosome, Nucleus") |> 
  pivot_longer(cols = c('seq', 'Mer13'), names_to = 'Mer13_adj_1', values_to = 'Mer13_adj_2')

write_xlsx(exosome_nucleus_13mer_HepG2, path = paste0(file_path, "exosome_nucleus_13mer_HepG2.xlsx"))

#exosome not nucleus
exosome_not_nucleus_13mer_HepG2 <- OG_glycopeptide_Top_tb_HepG2_singlesite |> 
  filter(UniprotID %in% (OG_glycoprotein_exosome_nucleus_HepG2_combined |> 
                           filter(category == "Exosome, not Nucleus") |> pull(UniprotID))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence) |> 
  mutate(seq = '>peptide', category = "Exosome, not Nucleus") |> 
  pivot_longer(cols = c('seq', 'Mer13'), names_to = 'Mer13_adj_1', values_to = 'Mer13_adj_2')

write_xlsx(exosome_not_nucleus_13mer_HepG2, path = paste0(file_path, "exosome_not_nucleus_13mer_HepG2.xlsx"))

#nucleus not exosome
nucleus_not_exosome_13mer_HepG2 <- OG_glycopeptide_Top_tb_HepG2_singlesite |> 
  filter(UniprotID %in% (OG_glycoprotein_exosome_nucleus_HepG2_combined |> 
                           filter(category == "Nucleus, not Exosome") |> pull(UniprotID))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence) |> 
  mutate(seq = '>peptide', category = "Nucleus, not Exosome") |> 
  pivot_longer(cols = c('seq', 'Mer13'), names_to = 'Mer13_adj_1', values_to = 'Mer13_adj_2')

write_xlsx(nucleus_not_exosome_13mer_HepG2, path = paste0(file_path, "nucleus_not_exosome_13mer_HepG2.xlsx"))

#import data
#exosome
iso_exosome_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\isoelectric\\exosome.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(category = 'Exosome')

iso_nucleus_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\isoelectric\\nucleus.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(category = 'Nucleus')

iso_exosome_nucleus_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\isoelectric\\exosome_nucleus.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(category = 'Exosome, Nucleus')

iso_exosome_not_nucleus_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\isoelectric\\exosome_not_nucleus.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(category = 'Exosome, not Nucleus')

iso_nucleus_not_exosome_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\isoelectric\\nucleus_not_exosome.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(category = 'Nucleus, not Exosome')

#combine
iso_mer13_combined <- bind_rows(
  iso_exosome_HepG2,
  iso_nucleus_HepG2,
  iso_exosome_nucleus_HepG2,
  iso_exosome_not_nucleus_HepG2,
  iso_nucleus_not_exosome_HepG2
)

#wilcox test
iso_mer13_wilcox_test <- iso_mer13_combined |> 
  wilcox_test(IPC2.peptide.Conv2D ~ category, p.adjust.method = "BH") |> 
  add_significance("p") |> 
  filter(p.signif == "**")

#boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_iso_13_HepG2 <- iso_mer13_combined |> 
  ggplot() +
  geom_boxplot(aes(x = factor(category, levels = c("Exosome", "Nucleus", "Exosome, Nucleus", "Nucleus, not Exosome", "Exosome, not Nucleus")), 
                   y = IPC2.peptide.Conv2D, color = factor(category, levels = c("Exosome", "Nucleus", "Exosome, Nucleus", "Nucleus, not Exosome", "Exosome, not Nucleus"))), 
               size = 1.5,
               outlier.shape = NA,
               notch = FALSE,
               notchwidth = 0.5) +
  labs(x = "", y = "isoelectric point") +
  scale_color_manual(values = c(
    "Exosome" = Color_1,
    "Nucleus" = Color_2,
    "Exosome, Nucleus" = Color_3,
    "Nucleus, not Exosome" = Color_4,
    "Exosome, not Nucleus" = Color_7
  )) +
  stat_pvalue_manual(data = iso_mer13_wilcox_test, hide.ns = "p", label = "p.signif", 
                     tip.length = 0, size = 50,
                     y.position = c(9.3, 9.6)) +
  coord_cartesian(ylim = c(3, 10)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "boxplot_iso_13_HepG2.png"), 
       plot = boxplot_iso_13_HepG2, height = 4, width = 2.5, units = "in", dpi = 600)
