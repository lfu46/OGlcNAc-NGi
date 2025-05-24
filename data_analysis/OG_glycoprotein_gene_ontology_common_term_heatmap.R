#import packages
packages_names <- c("tidyverse", "rstatix", "showtext", "ggpubr", "readxl")
lapply(packages_names, require, character.only = TRUE)
q
#generate dataframe for the unique glycoproteins in each cell line
OG_glycoprotein_unique_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> distinct(UniprotID)
OG_glycoprotein_unique_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> distinct(UniprotID)
OG_glycoprotein_unique_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> distinct(UniprotID)

OG_glycoprotein_only_HepG2 <- OG_glycoprotein_unique_HepG2 |> 
  filter(!(UniprotID %in% OG_glycoprotein_unique_HEK293T$UniprotID) & !(UniprotID %in% OG_glycoprotein_unique_Jurkat$UniprotID))
OG_glycoprotein_only_HEK293T <- OG_glycoprotein_unique_HEK293T |> 
  filter(!(UniprotID %in% OG_glycoprotein_unique_HepG2$UniprotID) & !(UniprotID %in% OG_glycoprotein_unique_Jurkat$UniprotID))
OG_glycoprotein_only_Jurkat <- OG_glycoprotein_unique_Jurkat |> 
  filter(!(UniprotID %in% OG_glycoprotein_unique_HepG2$UniprotID) & !(UniprotID %in% OG_glycoprotein_unique_HEK293T$UniprotID))

write_xlsx(OG_glycoprotein_only_HepG2, path = paste0(file_path, "OG_glycoprotein_only_HepG2.xlsx"))
write_xlsx(OG_glycoprotein_only_HEK293T, path = paste0(file_path, "OG_glycoprotein_only_HEK293T.xlsx"))
write_xlsx(OG_glycoprotein_only_Jurkat, path = paste0(file_path, "OG_glycoprotein_only_Jurkat.xlsx"))

#generate dataframe for the common glycoproteins in all three cell lines
OG_glycoprotein_common_total <- OG_glycoprotein_unique_HepG2 |> 
  filter((UniprotID %in% OG_glycoprotein_unique_HEK293T$UniprotID) & (UniprotID %in% OG_glycoprotein_unique_Jurkat$UniprotID))

write_xlsx(OG_glycoprotein_common_total, path = paste0(file_path, "OG_glycoprotein_common_total.xlsx"))

#import data
gene_ontology_OG_glycoprotein_common_only <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_common_only.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
) |> select(Term, P.Value, Cell) |> 
  pivot_wider(names_from = Cell, values_from = P.Value) |> 
  mutate(
    HepG2 = ifelse(is.na(HepG2), 1, HepG2),
    Jurkat = ifelse(is.na(Jurkat), 1, Jurkat),
    HEK293T = ifelse(is.na(HEK293T), 1, HEK293T),
    Common = ifelse(is.na(Common), 1, Common)
  )

colnames(gene_ontology_OG_glycoprotein_common_only) <- c("Term", "Common", "HepG2 only", "Jurkat only", "HEK293T only")

gene_ontology_OG_glycoprotein_common_only <- gene_ontology_OG_glycoprotein_common_only |> 
  pivot_longer(cols = 'Common':'HEK293T only', names_to = 'Cell', values_to = 'P.Value') |> 
  mutate(
    Term = factor(Term, levels = c(
      "RNA binding", "Nucleoplasm", "Cytoplasmic stress granule", "P-body", "DNA binding",
      "Humoral immune response", "Chromatin", "Positive regulation of canonical Wnt signaling pathway",
      "Lipid metabolic process", "Cholesterol homeostasis", "Cell-matrix adhesion", "Very-low-density lipoprotein particle assembly",
      "Cell surface receptor signaling pathway", "Positive regulation of T cell proliferation", "Plasma membrane", "Adaptive immune response"
    ))
  )

#heatmap
font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

gene_ontology_OG_glycoprotein_common_only_heatmap <- gene_ontology_OG_glycoprotein_common_only |> 
  ggplot(aes(x = Cell, y = Term, fill = P.Value)) +
  geom_tile() +
  labs(fill = "P.Value") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  scale_fill_stepsn(limits = c(0, 1), 
                    breaks = c(0, 1E-3, 0.01, 0.05, 1), 
                    labels = c("0", "1E-3", "0.01", "0.05", "1"), 
                    n.breaks = 4, 
                    values = scales::rescale(c(0, 1E-3, 0.01 ,0.05, 1)),
                    colours = c(Color_9, NA)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0),
    panel.grid.minor = element_line(color = "gray", linewidth = 0),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 12, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 10, color = "black", lineheight = 0.6),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, color = "black"),
    legend.frame = element_rect(color = "black", linewidth = 0.2),
    legend.ticks = element_line(color = "black", linewidth = 0.2),
    legend.key.size = unit(0.15, "in")
  )

ggsave(
  filename = paste0(file_path, "gene_ontology_OG_glycoprotein_common_only_heatmap.eps"), 
  device = cairo_ps,
  plot = gene_ontology_OG_glycoprotein_common_only_heatmap, 
  height = 3.5, width = 5, 
  units = "in", fallback_resolution = 1200
  )

