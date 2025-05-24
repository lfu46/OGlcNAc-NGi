#import packages
packages_names <- c("tidyverse", "rstatix", "ggpubr", "showtext", "readxl")
lapply(packages_names, require, character.only = TRUE)

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
font_add(family = "arial", regular = "arial.ttf")
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
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black", lineheight = 0.1),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, color = "black"),
    legend.frame = element_rect(color = "black", linewidth = 0.2),
    legend.ticks = element_line(color = "black", linewidth = 0.2)
  )

ggsave(filename = paste0(file_path, "heatmap_gene_ontology_OG_glycoprotein_up_down_regulated_total.eps"), 
       device = cairo_ps,
       plot = heatmap_gene_ontology_OG_glycoprotein_up_down_regulated_total, 
       height = 3, width = 5, 
       units = "in", 
       fallback_resolution = 1200
       )
