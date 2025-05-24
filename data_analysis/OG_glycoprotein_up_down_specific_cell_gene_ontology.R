#import packages
packages_names <- c("tidyverse", "rstatix", "ggpubr", "showtext", "readxl")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
#HepG2
OG_glycoprotein_up_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_up_HepG2, path = paste0(file_path, "OG_glycoprotein_up_HepG2.xlsx"))

OG_glycoprotein_down_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_down_HepG2, path = paste0(file_path, "OG_glycoprotein_down_HepG2.xlsx"))

#HEK293T
OG_glycoprotein_up_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_up_HEK293T, path = paste0(file_path, "OG_glycoprotein_up_HEK293T.xlsx"))

OG_glycoprotein_down_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_down_HEK293T, path = paste0(file_path, "OG_glycoprotein_down_HEK293T.xlsx"))

#Jurkat
OG_glycoprotein_up_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC > 0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_up_Jurkat, path = paste0(file_path, "OG_glycoprotein_up_Jurkat.xlsx"))

OG_glycoprotein_down_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC < -0.5, adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_down_Jurkat, path = paste0(file_path, "OG_glycoprotein_down_Jurkat.xlsx"))

#import data
gene_ontology_OG_glycoprotein_up_down_specific_cell_up_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_up_down_specific_cell.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "up_HepG2"
)
gene_ontology_OG_glycoprotein_up_down_specific_cell_down_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_up_down_specific_cell.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "down_Jurkat"
)
gene_ontology_OG_glycoprotein_up_down_specific_cell_up_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_up_down_specific_cell.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "up_HEK293T"
)

#bar plot
#up HepG2
font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

barplot_OG_glycoprotein_up_gene_ontology_HepG2 <- 
  gene_ontology_OG_glycoprotein_up_down_specific_cell_up_HepG2 |> 
  ggplot(aes(x = -log10(P.Value), y = fct_reorder(Term, -log10(P.Value)))) +
  geom_bar(stat = "identity", fill = Color_3, width = 0.7) +
  labs(x = expression(-log[10]*"("*paste(italic(P), " Value")*")"), y = "") +
  geom_text(aes(label = Term, x = 0), hjust = 0, size = 3) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "in")
  )

ggsave(filename = paste0(file_path, "barplot_OG_glycoprotein_up_gene_ontology_HepG2.eps"), 
       device = "eps",
       plot = barplot_OG_glycoprotein_up_gene_ontology_HepG2, 
       height = 1.5, width = 2, 
       units = "in", 
       dpi = 1200)

#down Jurkat
font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

barplot_OG_glycoprotein_down_gene_ontology_Jurkat <- 
  gene_ontology_OG_glycoprotein_up_down_specific_cell_down_Jurkat |> 
  ggplot(aes(x = -log10(P.Value), y = fct_reorder(Term, -log10(P.Value)))) +
  geom_bar(stat = "identity", fill = Color_4, width = 0.7) +
  labs(x = expression(-log[10]*"("*paste(italic(P), " Value")*")"), y = "") +
  geom_text(aes(label = Term, x = 0), hjust = 0, size = 3) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "in")
  )

ggsave(filename = paste0(file_path, "barplot_OG_glycoprotein_down_gene_ontology_Jurkat.eps"), 
       device = "eps",
       plot = barplot_OG_glycoprotein_down_gene_ontology_Jurkat, 
       height = 1.5, width = 2, 
       units = "in", 
       dpi = 1200
       )

#up HEK293T
font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

barplot_OG_glycoprotein_up_gene_ontology_HEK293T <- 
  gene_ontology_OG_glycoprotein_up_down_specific_cell_up_HEK293T |> 
  ggplot(aes(x = -log10(P.Value), y = fct_reorder(Term, -log10(P.Value)))) +
  geom_bar(stat = "identity", fill = Color_2, width = 0.7) +
  labs(x = expression(-log[10]*"("*paste(italic(P), " Value")*")"), y = "") +
  geom_text(aes(label = Term, x = 0), hjust = 0, size = 3) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "in")
  )

ggsave(filename = paste0(file_path, "barplot_OG_glycoprotein_up_gene_ontology_HEK293T.eps"), 
       device = "eps",
       plot = barplot_OG_glycoprotein_up_gene_ontology_HEK293T, 
       height = 1.5, width = 2, 
       units = "in", 
       dpi = 1200
       )
