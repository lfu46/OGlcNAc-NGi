#import packages
packages_names <- c("tidyverse", "rstatix", "ggpubr", "showtext", "readxl")
lapply(packages_names, require, character.only = TRUE)

#Dot plot
#import data
gene_ontology_OG_glycoprotein_up_regulated_total <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_generally_up_down.xlsx",
  sheet = 'up',
  col_names = TRUE,
  .name_repair = "universal"
)
gene_ontology_OG_glycoprotein_down_regulated_total <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_generally_up_down.xlsx",
  sheet = 'down',
  col_names = TRUE,
  .name_repair = "universal"
)

#generate data frame
gene_ontology_OG_glycoprotein_up_down_regulated_total <- bind_rows(
  gene_ontology_OG_glycoprotein_up_regulated_total |> select(Term, Count, P.Value) |> mutate(Group = "up"),
  gene_ontology_OG_glycoprotein_down_regulated_total |> select(Term, Count, P.Value) |> mutate(Group = "down")
) |> mutate(Term = factor(Term, levels = Term))

#balloon plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

heatmap_gene_ontology_OG_glycoprotein_up_down_regulated_total <- gene_ontology_OG_glycoprotein_up_down_regulated_total |> 
  ggplot(aes(x = factor(Group, levels = c("up", "down")), y = Term, size = Count, fill = P.Value)) +
  geom_point(shape = 21) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_size(range = c(5, 11), breaks = c(8, 16)) +
  scale_fill_stepsn(limits = c(0, 1), 
                    breaks = c(0, 5E-3, 7.5E-3, 0.01, 0.05, 1), 
                    labels = c("0", "5E-3", "7.5E-3", "0.01", "0.05", "1"), 
                    n.breaks = 6, 
                    values = scales::rescale(c(0, 5E-3, 7.5E-3, 0.01, 0.05, 1)),
                    colours = c(Color_9, "transparent")) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.2),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black", lineheight = 0.1),
    legend.text = element_text(size = 80),
    legend.title = element_text(size = 80, color = "black"),
    legend.frame = element_rect(color = "black", linewidth = 0.2),
    legend.ticks = element_line(color = "black", linewidth = 0.2)
  )

ggsave(filename = paste0(file_path, "heatmap_gene_ontology_OG_glycoprotein_up_down_regulated_total.png"), plot = heatmap_gene_ontology_OG_glycoprotein_up_down_regulated_total, height = 4, width = 6, units = c("in"), dpi = 600)

#isoelectric point
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
  filter(UniprotID %in% OG_glycoprotein_generally_up$UniprotID) |> 
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
font_add(family = "calibri", regular = "calibri.ttf")
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
                     size = 50,
                     hide.ns = "p", y.position = c(8.8)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 100, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "violin_boxplot_iso_OG_site_Mer13.png"), plot = violin_boxplot_iso_OG_site_Mer13, 
       height = 3.5, width = 5, units = "in", dpi = 600)

#Ordered vs. Disordered
OG_glycosite_disorder_1 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_1.csv",
  col_names = TRUE,
  name_repair = "universal"
)
OG_glycosite_disorder_2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_2.csv",
  col_names = TRUE,
  name_repair = "universal"
)
OG_glycosite_disorder_3 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_3.csv",
  col_names = TRUE,
  name_repair = "universal"
)
OG_glycosite_disorder_4 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_4.csv",
  col_names = TRUE,
  name_repair = "universal"
)
OG_glycosite_disorder_5 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_5.csv",
  col_names = TRUE,
  name_repair = "universal"
)
OG_glycosite_disorder_6 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_6.csv",
  col_names = TRUE,
  name_repair = "universal"
)
OG_glycosite_disorder_7 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_7.csv",
  col_names = TRUE,
  name_repair = "universal"
)
OG_glycosite_disorder_8 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_8.csv",
  col_names = TRUE,
  name_repair = "universal"
)
OG_glycosite_disorder_9 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_9.csv",
  col_names = TRUE,
  name_repair = "universal"
)
OG_glycosite_disorder_10 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_10.csv",
  col_names = TRUE,
  name_repair = "universal"
)
OG_glycosite_disorder_11 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_11.csv",
  col_names = TRUE,
  name_repair = "universal"
)
OG_glycosite_disorder_12 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_12.csv",
  col_names = TRUE,
  name_repair = "universal"
)
OG_glycosite_disorder_13 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_13.csv",
  col_names = TRUE,
  name_repair = "universal"
)

#combine
OG_glycosite_disorder_combined <- bind_rows(
  OG_glycosite_disorder_1,
  OG_glycosite_disorder_2,
  OG_glycosite_disorder_3,
  OG_glycosite_disorder_4,
  OG_glycosite_disorder_5,
  OG_glycosite_disorder_6,
  OG_glycosite_disorder_7,
  OG_glycosite_disorder_8,
  OG_glycosite_disorder_9,
  OG_glycosite_disorder_10,
  OG_glycosite_disorder_11,
  OG_glycosite_disorder_12,
  OG_glycosite_disorder_13
)

#fix some sequence
OG_glycosite_disorder_combined <- OG_glycosite_disorder_combined |> 
  mutate(n = ifelse(str_detect(id, "Q09666_2"), n + 3000, n)) |> 
  mutate(n = ifelse(str_detect(id, "Q8IVF2_2"), n + 3000, n)) |> 
  mutate(n = ifelse(str_detect(id, "O14686_2"), n + 3000, n)) |> 
  mutate(n = ifelse(str_detect(id, "Q9Y6V0_2"), n + 3000, n))

#generate data frame
OG_glycosite_disorder_combined <- OG_glycosite_disorder_combined |> 
  mutate(UniprotID = str_extract(id, "(?<=^....).{6}"),
         Site = paste0(seq, n)) |> 
  mutate(Protein_site = paste(UniprotID, Site, sep = "_")) |> 
  select(UniprotID, Protein_site, disorder)

#generate data frame
#HepG2
OG_glycosite_disorder_HepG2 <- OG_glycopeptide_Top_tb_HepG2_singlesite |> 
  #filter(UniprotID %in% OG_glycoprotein_generally_down$UniprotID) |> 
  mutate(Protein_site = paste(UniprotID, combined_site, sep = "_"),
         Cell = "HepG2") |> 
  select(Index, Protein_site, UniprotID, combined_site, logFC, adj.P.Val, Cell) |> 
  left_join(OG_glycosite_disorder_combined, by = "Protein_site")

#HEK293T
OG_glycosite_disorder_HEK293T <- OG_glycopeptide_Top_tb_HEK293T_singlesite |> 
  #filter(UniprotID %in% OG_glycoprotein_generally_down$UniprotID) |> 
  mutate(Protein_site = paste(UniprotID, combined_site, sep = "_"),
         Cell = "HEK293T") |> 
  select(Index, Protein_site, UniprotID, combined_site, logFC, adj.P.Val, Cell) |> 
  left_join(OG_glycosite_disorder_combined, by = "Protein_site")

#Jurkat
OG_glycosite_disorder_Jurkat <- OG_glycopeptide_Top_tb_Jurkat_singlesite |> 
  #filter(UniprotID %in% OG_glycoprotein_generally_down$UniprotID) |> 
  mutate(Protein_site = paste(UniprotID, combined_site, sep = "_"),
         Cell = "Jurkat") |> 
  select(Index, Protein_site, UniprotID, combined_site, logFC, adj.P.Val, Cell) |> 
  left_join(OG_glycosite_disorder_combined, by = "Protein_site")

#combined
OG_glycosite_disorder_tb <- bind_rows(
  OG_glycosite_disorder_HepG2,
  OG_glycosite_disorder_HEK293T,
  OG_glycosite_disorder_Jurkat
) |> mutate(
  Disorder = case_when(
    disorder > 0.5 ~ "Disordered",
    disorder < 0.5 ~ "Ordered"
  )
)

write_xlsx(OG_glycosite_disorder_tb, path = paste0(file_path, "OG_glycosite_disorder_tb.xlsx"))

#wilcox test
OG_glycosite_disorder_wilcox_test <- OG_glycosite_disorder_tb |> 
  mutate(
    Disorder = case_when(
      disorder > 0.5 ~ "Disordered",
      disorder < 0.5 ~ "Ordered"
    )
  )|> 
  group_by(Cell) |> 
  wilcox_test(logFC ~ Disorder, p.adjust.method = "BH") |> 
  add_significance("p")

#violin point plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_point_OG_glycosite_disorder <- OG_glycosite_disorder_tb |> 
  filter(!is.na(disorder)) |> 
  mutate(
    Disorder = case_when(
      disorder > 0.5 ~ "Disordered",
      disorder < 0.5 ~ "Ordered"
    )
  )|> 
  ggplot() +
  geom_violin(aes(x = factor(Disorder, levels = c("Ordered", "Disordered")), y = logFC, fill = Cell), 
              color = "transparent") +
  geom_boxplot(aes(x = factor(Disorder, levels = c("Ordered", "Disordered")), y = logFC), color = "black",
               outliers = FALSE, width = 0.2) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_wrap(vars(Cell), scales = "free", nrow = 1) +
  stat_pvalue_manual(data = OG_glycosite_disorder_wilcox_test, label = "p.signif", tip.length = 0, 
                     size = 50,
                     hide.ns = "p", y.position = c(1.8)) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 100, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "violin_point_OG_glycosite_disorder.png"), plot = violin_point_OG_glycosite_disorder, height = 3.5, width = 5, units = "in", dpi = 600)
