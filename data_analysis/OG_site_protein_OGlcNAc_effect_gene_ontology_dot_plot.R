#import packages
packages_names <- c("tidyverse", "rstatix", "showtext", "ggpubr", "readxl")
lapply(packages_names, require, character.only = TRUE)

#import data
OG_site_protein_up_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_site_protein.xlsx",
  sheet = 'up_HEK293T',
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(Cell = "HEK293T")

OG_site_protein_up_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_site_protein.xlsx",
  sheet = 'up_HepG2',
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(Cell = "HepG2")

OG_site_protein_up_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_site_protein.xlsx",
  sheet = 'up_Jurkat',
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(Cell = "Jurkat")

OG_site_protein_down_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_site_protein.xlsx",
  sheet = 'down_HEK293T',
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(Cell = "HEK293T")

OG_site_protein_down_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_site_protein.xlsx",
  sheet = 'down_HepG2',
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(Cell = "HepG2")

OG_site_protein_down_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_site_protein.xlsx",
  sheet = 'down_Jurkat',
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(Cell = "Jurkat")

#generate data frame
#up
OG_site_protein_up_combined <- bind_rows(
  OG_site_protein_up_HEK293T,
  OG_site_protein_up_HepG2,
  OG_site_protein_up_Jurkat
) |> select(Term, Count, P.Value, Exp)

font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_OG_site_protein_up <- OG_site_protein_up_combined |> 
  ggplot(aes(x = Exp, y = factor(Term, levels = Term))) +
  geom_point(aes(fill = Exp, size = Count, alpha = P.Value), shape = 21) +
  scale_x_discrete(labels = c("HEK293T", "HepG2", "Jurkat")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_size(range = c(6, 10)) +
  scale_alpha_binned(range = c(1, 0.2)) +
  scale_fill_manual(
    name = "Cell",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.2),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black", lineheight = 0.1),
    legend.text = element_text(size = 80),
    legend.title = element_text(size = 80, color = "black"),
    legend.frame = element_rect(color = "black", linewidth = 0.2),
    legend.ticks = element_line(color = "black", linewidth = 0.2)
  )

ggsave(filename = paste0(file_path, "dot_plot_OG_site_protein_up.png"), 
       plot = dot_plot_OG_site_protein_up, 
       height = 7, width = 6, units = c("in"), dpi = 600)

#down
OG_site_protein_down_combined <- bind_rows(
  OG_site_protein_down_HEK293T,
  OG_site_protein_down_HepG2,
  OG_site_protein_down_Jurkat
) |> select(Term, Count, P.Value, Cell)

font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_OG_site_protein_down <- OG_site_protein_down_combined |> 
  ggplot(aes(x = Cell, y = factor(Term, levels = OG_site_protein_down_combined |> distinct(Term) |> pull()))) +
  geom_point(aes(fill = Cell, size = Count, alpha = P.Value), shape = 21) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_size(range = c(6, 10)) +
  scale_alpha_binned(range = c(1, 0.2), breaks = c(0.001, 0.005, 0.1)) +
  scale_fill_manual(
    name = "Cell",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.2),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black", lineheight = 0.1),
    legend.text = element_text(size = 80),
    legend.title = element_text(size = 80, color = "black"),
    legend.frame = element_rect(color = "black", linewidth = 0.2),
    legend.ticks = element_line(color = "black", linewidth = 0.2)
  )

ggsave(filename = paste0(file_path, "dot_plot_OG_site_protein_down.png"), 
       plot = dot_plot_OG_site_protein_down, 
       height = 7, width = 6, units = c("in"), dpi = 600)
