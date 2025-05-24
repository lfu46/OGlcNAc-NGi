#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "rstatix", "ggpubr", "showtext")
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
#regulated HEK293T
OG_site_Mer31_regulated_HEK293T <- OG_glycopeptide_Top_tb_HEK293T_singlesite |> 
  mutate(category = ifelse(logFC > 0.5 & adj.P.Val < 0.05, "up", ifelse(logFC < -0.5 & adj.P.Val < 0.05, "down", "median"))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence, category) |> 
  mutate(Mer31 = substr(Sequence, site_position-15, site_position+15)) |> 
  select(!Sequence) |> 
  mutate(cell = "HEK293T")

#site level
#generate data frame for 31mer
#regulated HepG2
OG_site_Mer31_regulated_HepG2 <- OG_glycopeptide_Top_tb_HepG2_singlesite |> 
  mutate(category = ifelse(logFC > 0.5 & adj.P.Val < 0.05, "up", ifelse(logFC < -0.5 & adj.P.Val < 0.05, "down", "median"))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence, category) |> 
  mutate(Mer31 = substr(Sequence, site_position-15, site_position+15)) |> 
  select(!Sequence) |> 
  mutate(cell = "HepG2")

#site level
#generate data frame for 31mer
#regulated Jurkat
OG_site_Mer31_regulated_Jurkat <- OG_glycopeptide_Top_tb_Jurkat_singlesite |> 
  mutate(category = ifelse(logFC > 0.5 & adj.P.Val < 0.05, "up", ifelse(logFC < -0.5 & adj.P.Val < 0.05, "down", "median"))) |> 
  left_join(OG_glycoprotein_uniprot_sequence, by = join_by(UniprotID == Entry)) |> 
  select(Index, logFC, adj.P.Val, UniprotID, site_position, Sequence, category) |> 
  mutate(Mer31 = substr(Sequence, site_position-15, site_position+15)) |> 
  select(!Sequence) |> 
  mutate(cell = "Jurkat")

#combined
OG_site_Mer31_regulated_combined <- bind_rows(
  OG_site_Mer31_regulated_HEK293T,
  OG_site_Mer31_regulated_HepG2,
  OG_site_Mer31_regulated_Jurkat
)

write_xlsx(OG_site_Mer31_regulated_combined, path = paste0(file_path, "OG_site_Mer31_regulated_combined.xlsx"))

#import result from localcider
OG_site_Mer31_regulated_combined_sequence_feature <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_site_Mer31_regulated_combined_sequence_parameter.xlsx"
)

#wilcox test
#isoelectric point
isoelectric_point_wilcox_test <- OG_site_Mer31_regulated_combined_sequence_feature |> 
  filter(category %in% c('up', 'down')) |> 
  group_by(cell) |> 
  wilcox_test(Isoelectric_point ~ category, p.adjust.method = 'BH') |> 
  add_significance('p')

#violin boxplot
font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

violin_boxplot_isoelectric_point_Mer31 <- OG_site_Mer31_regulated_combined_sequence_feature |> 
  filter(category %in% c('up', 'down')) |> 
  ggplot() +
  geom_violin(aes(x = factor(category, levels = c("up", "down")), y = Isoelectric_point, fill = cell), 
              color = "transparent") +
  geom_boxplot(aes(x = factor(category, levels = c("up", "down")), y = Isoelectric_point), color = "black",
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
  stat_pvalue_manual(data = isoelectric_point_wilcox_test |> filter(), label = "p.signif", tip.length = 0, 
                     size = 6,
                     hide.ns = "p", y.position = c(14)) +
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

ggsave(filename = paste0(file_path, "violin_boxplot_isoelectric_point_Mer31.eps"),
       device = "eps",
       plot = violin_boxplot_isoelectric_point_Mer31, 
       height = 2.5, width = 2.5, units = "in"
)

#FCR
FCR_wilcox_test <- OG_site_Mer31_regulated_combined_sequence_feature |> 
  filter(category %in% c('up', 'down')) |> 
  group_by(cell) |> 
  wilcox_test(FCR ~ category, p.adjust.method = 'BH') |> 
  add_significance('p')

#violin boxplot
font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

violin_boxplot_FCR_Mer31 <- OG_site_Mer31_regulated_combined_sequence_feature |> 
  filter(category %in% c('up', 'down')) |> 
  ggplot() +
  geom_violin(aes(x = factor(category, levels = c("up", "down")), y = FCR, fill = cell), 
              color = "transparent") +
  geom_boxplot(aes(x = factor(category, levels = c("up", "down")), y = FCR), color = "black",
               outliers = FALSE, width = 0.2) +
  labs(x = "", y = "FCR") +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  coord_cartesian(ylim = c(0, 0.5)) +
  facet_grid(~ cell, scales = "free_x") +
  stat_pvalue_manual(data = FCR_wilcox_test |> filter(), label = "p.signif", tip.length = 0, 
                     size = 6,
                     hide.ns = "p", y.position = c(0.45)) +
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

ggsave(filename = paste0(file_path, "violin_boxplot_FCR_Mer31.eps"),
       device = "eps",
       plot = violin_boxplot_FCR_Mer31, 
       height = 2.5, width = 2.5, units = "in"
)
