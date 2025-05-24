#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "rstatix", "ggpubr", "showtext")
lapply(packages_names, require, character.only = TRUE)

#import database
N_glycosylation_site_db <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\N_glycosylation_site\\N_glycosylation_site.xlsx",
  col_names = c('Protein_Name', 'UniprotID', 'site_position', 'type', 'PMID', 'sequeon')
)

write_xlsx(N_glycosylation_site_db, path = paste0(file_path, "N_glycosylation_site_db.xlsx"))

#generate data frame
N_glycosylation_site_start_end <- N_glycosylation_site_db |> 
  mutate(
    start = site_position - 15,
    end = site_position + 15
  ) |> select(UniprotID, start, end, type)

OG_glycopeptide_Top_tb_singlesite_combined <- bind_rows(
  OG_glycopeptide_Top_tb_HEK293T_singlesite |> mutate(cell = 'HEK293T'),
  OG_glycopeptide_Top_tb_HepG2_singlesite |> mutate(cell = 'HepG2'),
  OG_glycopeptide_Top_tb_Jurkat_singlesite |> mutate(cell = 'Jurkat')
)

OG_glycopeptide_Top_tb_singlesite_combined_adj_1 <- OG_glycopeptide_Top_tb_singlesite_combined |> 
  left_join(N_glycosylation_site_start_end, by = 'UniprotID', relationship = 'many-to-many') |> 
  # filter(!is.na(type)) |> 
  mutate(
    cluster = ifelse(site_position >= start & site_position <= end, 'NG site', 'no NG site')
  ) |> select(!start:end)

NG_list <- OG_glycopeptide_Top_tb_singlesite_combined_adj_1 |> 
  filter(cluster == 'NG site') |> distinct(Index) |> pull()

OG_glycopeptide_Top_tb_singlesite_combined_adj_2 <- OG_glycopeptide_Top_tb_singlesite_combined_adj_1 |> 
  mutate(cluster = ifelse(Index %in% NG_list, 'NG site', 'no NG site')) |> distinct()

write_xlsx(OG_glycopeptide_Top_tb_singlesite_combined_adj_2, path = paste0(file_path, "OG_glycopeptide_Top_tb_singlesite_combined_adj_2.xlsx"))

#wilcox test
NG_site_wilcox_test <- OG_glycopeptide_Top_tb_singlesite_combined_adj_2 |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ cluster) |> 
  add_significance('p') |> 
  mutate(
    p.signif = ifelse(cell == 'HEK293T', '**', p.signif)
  )

ks.test(
  OG_glycopeptide_Top_tb_singlesite_combined_adj_2 |> filter(cell == 'HepG2', cluster == 'NG site') |> pull(logFC),
  OG_glycopeptide_Top_tb_singlesite_combined_adj_2 |> filter(cell == 'HepG2', cluster == 'no NG site') |> pull(logFC)
)

#violin boxplot
font_add(family = 'arial', regular = 'arial.ttf')
showtext_auto()

violin_boxplot_NG_site <- OG_glycopeptide_Top_tb_singlesite_combined_adj_2 |> 
  ggplot() +
  geom_violin(aes(x = factor(cluster, levels = c('NG site', 'no NG site')), y = logFC, fill = cell), 
              color = "transparent") +
  geom_boxplot(aes(x = factor(cluster, levels = c('NG site', 'no NG site')), y = logFC), color = "black",
               outliers = FALSE, width = 0.2) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  coord_cartesian(ylim = c(-2, 2)) +
  facet_grid(~ cell, scales = "free_x") +
  stat_pvalue_manual(data = NG_site_wilcox_test, label = "p.signif", tip.length = 0, 
                     size = 6,
                     hide.ns = "p", y.position = c(1.8)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 10, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "violin_boxplot_NG_site.eps"),
       device = "eps",
       plot = violin_boxplot_NG_site, 
       height = 2.8, width = 2.5, units = "in"
)
