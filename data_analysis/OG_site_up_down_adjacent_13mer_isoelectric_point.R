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
#generate data frame for 31mer
#regulated HepG2
OG_site_Mer31_regulated_HepG2 <- OG_glycopeptide_Top_tb_HepG2_singlesite |> 
  mutate(category = ifelse(logFC > 0.5 & adj.P.Val < 0.05, "up", ifelse(logFC < -0.5 & adj.P.Val < 0.05, "down", "median"))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence, category) |> 
  mutate(Mer31 = substr(Sequence, site_position-15, site_position+15)) |> 
  select(!Sequence) |> 
  mutate(cell = "HepG2", seq = '>peptide') |> 
  pivot_longer(cols = c('seq', 'Mer31'), names_to = 'Mer31_adj_1', values_to = 'Mer31_adj_2')

write_xlsx(OG_site_Mer31_regulated_HepG2, path = paste0(file_path, "OG_site_Mer31_regulated_HepG2.xlsx"))

#site level
#generate data frame for 31mer 
#regulated HEK293T
OG_site_Mer31_regulated_HEK293T <- OG_glycopeptide_Top_tb_HEK293T_singlesite |> 
  mutate(category = ifelse(logFC > 0.5 & adj.P.Val < 0.05, "up", ifelse(logFC < -0.5 & adj.P.Val < 0.05, "down", "median"))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence, category) |> 
  mutate(Mer31 = substr(Sequence, site_position-15, site_position+15)) |> 
  select(!Sequence) |> 
  mutate(cell = "HEK293T", seq = '>peptide') |> 
  pivot_longer(cols = c('seq', 'Mer31'), names_to = 'Mer31_adj_1', values_to = 'Mer31_adj_2')

write_xlsx(OG_site_Mer31_regulated_HEK293T, path = paste0(file_path, "OG_site_Mer31_regulated_HEK293T.xlsx"))

#site level
#generate data frame for 13mer 
#regulated Jurkat
OG_site_Mer13_all_Jurkat <- OG_glycopeptide_Top_tb_Jurkat_singlesite |> 
  mutate(category = ifelse(logFC > 0.5 & adj.P.Val < 0.05, "up", ifelse(logFC < -0.5 & adj.P.Val < 0.05, "down", "median"))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence) |> 
  mutate(Mer13 = substr(Sequence, site_position-6, site_position+6)) |> 
  select(!Sequence) |> 
  mutate(cell = "Jurkat", seq = '>peptide') |> 
  pivot_longer(cols = c('seq', 'Mer13'), names_to = 'Mer13_adj_1', values_to = 'Mer13_adj_2')

write_xlsx(OG_site_Mer13_all_Jurkat, path = paste0(file_path, "OG_site_Mer13_all_Jurkat.xlsx"))

#import data
#HEK293T
iso_OG_site_Mer13_up_HEK293T <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_up_regulated_HEK293T.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(cell = 'HEK293T', category = 'up')
iso_OG_site_Mer13_down_HEK293T <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_down_regulated_HEK293T.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(cell = 'HEK293T', category = 'down')

#HepG2
iso_OG_site_Mer13_up_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_up_regulated_HepG2.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(cell = 'HepG2', category = 'up')
iso_OG_site_Mer13_down_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_down_regulated_HepG2.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(cell = 'HepG2', category = 'down')

#Jurkat
iso_OG_site_Mer13_up_Jurkat <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_up_regulated_Jurkat.csv",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(cell = 'Jurkat', category = 'up')
iso_OG_site_Mer13_down_Jurkat <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_down_regulated_Jurkat.csv",
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
  facet_grid(~ cell, scales = "free_x") +
  stat_pvalue_manual(data = iso_OG_site_Mer13_wilcox_test, label = "p.signif", tip.length = 0, 
                     size = 6,
                     hide.ns = "p", y.position = c(8.8)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 10, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "violin_boxplot_iso_OG_site_Mer13.eps"),
       device = "eps",
       plot = violin_boxplot_iso_OG_site_Mer13, 
       height = 2.5, width = 2.5, units = "in"
       )

#gravy
#Jurkat
gravy_OG_site_Mer13_up_regulated_Jurkat <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gravy\\gravy_OG_site_Mer13_up_regulated_Jurkat.txt",
  delim = ";",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(Exp = "up", cell = "Jurkat")
gravy_OG_site_Mer13_down_regulated_Jurkat <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gravy\\gravy_OG_site_Mer13_down_regulated_Jurkat.txt",
  delim = ";",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(Exp = "down", cell = "Jurkat")

#HEK293T
gravy_OG_site_Mer13_up_regulated_HEK293T <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gravy\\gravy_OG_site_Mer13_up_regulated_HEK293T.txt",
  delim = ";",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(Exp = "up", cell = "HEK293T")
gravy_OG_site_Mer13_down_regulated_HEK293T <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gravy\\gravy_OG_site_Mer13_down_regulated_HEK293T.txt",
  delim = ";",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(Exp = "down", cell = "HEK293T")

#HepG2
gravy_OG_site_Mer13_up_regulated_HepG2 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gravy\\gravy_OG_site_Mer13_up_regulated_HEK293T.txt",
  delim = ";",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(Exp = "up", cell = "HepG2")
gravy_OG_site_Mer13_down_regulated_HepG2 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gravy\\gravy_OG_site_Mer13_down_regulated_HEK293T.txt",
  delim = ";",
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(Exp = "down", cell = "HepG2")

#combined
OG_site_Mer13_regulated_combined <- bind_rows(
  gravy_OG_site_Mer13_up_regulated_HEK293T,
  gravy_OG_site_Mer13_down_regulated_HEK293T,
  gravy_OG_site_Mer13_up_regulated_HepG2,
  gravy_OG_site_Mer13_down_regulated_HepG2,
  gravy_OG_site_Mer13_up_regulated_Jurkat,
  gravy_OG_site_Mer13_down_regulated_Jurkat
)

OG_site_Mer13_regulated_combined |>
  group_by(cell) |> 
  wilcox_test(GRAVY ~ Exp, p.adjust.method = "BH") |> 
  add_significance("p")

#import data
iso_OG_site_Mer13_all_Jurkat_1 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_all_Jurkat_1.csv",
  col_names = TRUE,
  name_repair = "universal"
)
iso_OG_site_Mer13_all_Jurkat_2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_all_Jurkat_2.csv",
  col_names = TRUE,
  name_repair = "universal"
)
iso_OG_site_Mer13_all_Jurkat_3 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_all_Jurkat_3.csv",
  col_names = TRUE,
  name_repair = "universal"
)
iso_OG_site_Mer13_all_Jurkat_4 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_all_Jurkat_4.csv",
  col_names = TRUE,
  name_repair = "universal"
)
iso_OG_site_Mer13_all_Jurkat_5 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_all_Jurkat_5.csv",
  col_names = TRUE,
  name_repair = "universal"
)
iso_OG_site_Mer13_all_Jurkat_6 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\ipc2\\OG_site_adjacent_mer13_all_Jurkat_6.csv",
  col_names = TRUE,
  name_repair = "universal"
)

#combined
iso_OG_site_Mer13_all_Jurkat_combined <- bind_rows(
  iso_OG_site_Mer13_all_Jurkat_1,
  iso_OG_site_Mer13_all_Jurkat_2,
  iso_OG_site_Mer13_all_Jurkat_3,
  iso_OG_site_Mer13_all_Jurkat_4,
  iso_OG_site_Mer13_all_Jurkat_5,
  iso_OG_site_Mer13_all_Jurkat_6
) |> left_join(OG_site_Mer13_all_Jurkat, by = join_by('sequence' == 'Mer13_adj_2')) |> 
  filter(!is.na(IPC2.peptide.Conv2D)) |> 
  select(sequence, IPC2.peptide.Conv2D, logFC)

cor.test(iso_OG_site_Mer13_all_Jurkat_combined$IPC2.peptide.Conv2D, iso_OG_site_Mer13_all_Jurkat_combined$logFC,
         method = "spearman")
