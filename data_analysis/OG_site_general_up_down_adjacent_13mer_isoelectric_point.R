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
#up HepG2
OG_site_Mer13_up_HepG2 <- OG_glycopeptide_Top_tb_HepG2_singlesite |> 
  mutate(category = ifelse(logFC > 0.5 & adj.P.Val < 0.05, "up", ifelse(logFC < -0.5 & adj.P.Val < 0.05, "down", "median"))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence, category) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence) |> 
  mutate(cell = "HepG2", seq = '>peptide') |> 
  pivot_longer(cols = c('seq', 'Mer13'), names_to = 'Mer13_adj_1', values_to = 'Mer13_adj_2')

write_xlsx(OG_site_Mer13_up_HepG2, path = paste0(file_path, "OG_site_Mer13_up_HepG2.xlsx"))

#site level
#generate data frame for 13mer 
#up HEK293T
OG_site_Mer13_up_HEK293T <- OG_glycopeptide_Top_tb_HEK293T_singlesite |> 
  filter(UniprotID %in% OG_glycoprotein_generally_up$UniprotID) |> 
  mutate(category = ifelse(logFC > 0.5 & adj.P.Val < 0.05, "up", ifelse(logFC < -0.5 & adj.P.Val < 0.05, "down", "median"))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence, category) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence) |> 
  mutate(cell = "HEK293T", seq = '>peptide') |> 
  pivot_longer(cols = c('seq', 'Mer13'), names_to = 'Mer13_adj_1', values_to = 'Mer13_adj_2')

write_xlsx(OG_site_Mer13_up_HEK293T, path = paste0(file_path, "OG_site_Mer13_up_HEK293T.xlsx"))

#site level
#generate data frame for 13mer 
#up Jurkat
OG_site_Mer13_up_Jurkat <- OG_glycopeptide_Top_tb_Jurkat_singlesite |> 
  filter(UniprotID %in% OG_glycoprotein_generally_up$UniprotID) |> 
  mutate(category = ifelse(logFC > 0.5 & adj.P.Val < 0.05, "up", ifelse(logFC < -0.5 & adj.P.Val < 0.05, "down", "median"))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence, category) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence) |> 
  mutate(cell = "Jurkat", seq = '>peptide') |> 
  pivot_longer(cols = c('seq', 'Mer13'), names_to = 'Mer13_adj_1', values_to = 'Mer13_adj_2')

write_xlsx(OG_site_Mer13_up_Jurkat, path = paste0(file_path, "OG_site_Mer13_up_Jurkat.xlsx"))

#site level
#generate data frame for 13mer
#down HepG2
OG_site_Mer13_down_HepG2 <- OG_glycopeptide_Top_tb_HepG2_singlesite |> 
  filter(UniprotID %in% OG_glycoprotein_generally_down$UniprotID) |> 
  mutate(category = ifelse(logFC > 0.5 & adj.P.Val < 0.05, "up", ifelse(logFC < -0.5 & adj.P.Val < 0.05, "down", "median"))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence, category) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence) |> 
  mutate(cell = "HepG2", seq = '>peptide') |> 
  pivot_longer(cols = c('seq', 'Mer13'), names_to = 'Mer13_adj_1', values_to = 'Mer13_adj_2')

write_xlsx(OG_site_Mer13_down_HepG2, path = paste0(file_path, "OG_site_Mer13_down_HepG2.xlsx"))

#site level
#generate data frame for 13mer 
#down HEK293T
OG_site_Mer13_down_HEK293T <- OG_glycopeptide_Top_tb_HEK293T_singlesite |> 
  filter(UniprotID %in% OG_glycoprotein_generally_down$UniprotID) |> 
  mutate(category = ifelse(logFC > 0.5 & adj.P.Val < 0.05, "up", ifelse(logFC < -0.5 & adj.P.Val < 0.05, "down", "median"))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence, category) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence) |> 
  mutate(cell = "HEK293T", seq = '>peptide') |> 
  pivot_longer(cols = c('seq', 'Mer13'), names_to = 'Mer13_adj_1', values_to = 'Mer13_adj_2')

write_xlsx(OG_site_Mer13_down_HEK293T, path = paste0(file_path, "OG_site_Mer13_down_HEK293T.xlsx"))

#site level
#generate data frame for 13mer 
#down Jurkat
OG_site_Mer13_down_Jurkat <- OG_glycopeptide_Top_tb_Jurkat_singlesite |> 
  filter(UniprotID %in% OG_glycoprotein_generally_down$UniprotID) |> 
  mutate(category = ifelse(logFC > 0.5 & adj.P.Val < 0.05, "up", ifelse(logFC < -0.5 & adj.P.Val < 0.05, "down", "median"))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence, category) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence) |> 
  mutate(cell = "Jurkat", seq = '>peptide') |> 
  pivot_longer(cols = c('seq', 'Mer13'), names_to = 'Mer13_adj_1', values_to = 'Mer13_adj_2')

write_xlsx(OG_site_Mer13_down_Jurkat, path = paste0(file_path, "OG_site_Mer13_down_Jurkat.xlsx"))

#import data
#HEK293T
iso_OG_site_Mer13_up_HEK293T <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_general_up_HEK293T.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(cell = 'HEK293T', category = 'up')
iso_OG_site_Mer13_down_HEK293T <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_general_down_HEK293T.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(cell = 'HEK293T', category = 'down')

#HepG2
iso_OG_site_Mer13_up_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_general_up_HepG2.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(cell = 'HepG2', category = 'up')
iso_OG_site_Mer13_down_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_general_down_HepG2.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(cell = 'HepG2', category = 'down')

#Jurkat
iso_OG_site_Mer13_up_Jurkat <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_general_up_Jurkat.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(cell = 'Jurkat', category = 'up')
iso_OG_site_Mer13_down_Jurkat <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_general_down_Jurkat.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(cell = 'Jurkat', category = 'down')

#combine
iso_OG_site_Mer13_combined <- bind_rows(
  iso_OG_site_Mer13_up_HEK293T,
  iso_OG_site_Mer13_down_HEK293T,
  iso_OG_site_Mer13_up_HepG2,
  iso_OG_site_Mer13_down_HepG2,
  iso_OG_site_Mer13_up_Jurkat,
  iso_OG_site_Mer13_down_Jurkat
)

#wilcox test
iso_OG_site_Mer13_wilcox_test <- iso_OG_site_Mer13_combined |> 
  group_by(cell) |> 
  wilcox_test(IPC2.peptide.Conv2D ~ category, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

violin_boxplot_iso_OG_site_Mer13 <- iso_OG_site_Mer13_combined |> 
  ggplot() +
  geom_violin(aes(x = factor(category, levels = c("up", "down")), y = IPC2.peptide.Conv2D, fill = cell), 
              color = "transparent") +
  geom_boxplot(aes(x = factor(category, levels = c("up", "down")), y = IPC2.peptide.Conv2D), color = "black",
               outliers = FALSE, width = 0.2) +
  labs(x = "", y = "isoelectric point") +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_wrap(vars(cell), scales = "free", nrow = 1) +
  stat_pvalue_manual(data = iso_OG_site_Mer13_wilcox_test, label = "p.signif", tip.length = 0, 
                     size = 10,
                     hide.ns = "p", y.position = c(8.8)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 15, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "violin_boxplot_iso_OG_site_Mer13.eps"),
       device = "eps",
       plot = violin_boxplot_iso_OG_site_Mer13, 
       height = 3.5, width = 5, units = "in", 
       dpi = 1200)
