#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "UpSetR", "ggupset", "eulerr")
lapply(packages_names, require, character.only = TRUE)

#Upset plot
#generate total glycoprotein
OG_glycoprotein_total <- bind_rows(
  OG_glycoprotein_Top_tb_HepG2 |> distinct(UniprotID),
  OG_glycoprotein_Top_tb_HEK293T |> distinct(UniprotID),
  OG_glycoprotein_Top_tb_Jurkat |> distinct(UniprotID)
) |> distinct()

OG_glycoprotein_upset_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> select(UniprotID) |> mutate(HepG2 = 1)

OG_glycoprotein_upset_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> select(UniprotID) |> mutate(HEK293T = 1)

OG_glycoprotein_upset_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> select(UniprotID) |> mutate(Jurkat = 1)

OG_glycoprotein_total_upset <- OG_glycoprotein_total |> 
  left_join(OG_glycoprotein_upset_HepG2, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_upset_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_upset_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  mutate(
    HepG2 = ifelse(is.na(HepG2), 0, 1),
    HEK293T = ifelse(is.na(HEK293T), 0, 1),
    Jurkat = ifelse(is.na(Jurkat), 0, 1)
  )

OG_glycoprotein_total_upset_dataframe <- as.data.frame(OG_glycoprotein_total_upset)

upset(OG_glycoprotein_total_upset_dataframe, sets = c("HepG2", "HEK293T", "Jurkat"),
      main.bar.color = Color_2,
      sets.bar.color = Color_2,
      order.by = "freq", point.size = 5,
      text.scale = c(2, 1.5, 1.5, 1.5, 2, 3),
      mainbar.y.label = "# of quantified \nO-GlcNAcylated proteins")

#Euler plot
#site level
#generate data frame
#HepG2
OG_site_HepG2 <- OG_glycopeptide_Top_tb_HepG2 |> 
  filter(!is.na(combined_site)) |> 
  pull(Index)

#HEK293T
OG_site_HEK293T <- OG_glycopeptide_Top_tb_HEK293T |> 
  filter(!is.na(combined_site)) |> 
  pull(Index)

#Jurkat
OG_site_Jurkat <- OG_glycopeptide_Top_tb_Jurkat |> 
  filter(!is.na(combined_site)) |> 
  pull(Index)


#euler diagram
mat <- c(
  "HEK293T (2146)" = 1368,
  "HepG2 (439)" = 109,
  "Jurkat (1324)" = 580,
  "HEK293T (2146)&HepG2 (439)" = 70,
  "HEK293T (2146)&Jurkat (1324)" = 484,
  "HepG2 (439)&Jurkat (1324)" = 36,
  "HepG2 (439)&HEK293T (2146)&Jurkat (1324)" = 224
)

fit <- euler(mat)
plot(fit, 
     fill = list(fill = c(Color_2, Color_3, Color_4), alpha = 0.8),
     edges = list(col = "white", lwd = 3),
     labels = list(alpha = c(0)),
     legend = list(fontsize = 15))

#Heatmap
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
font_add(family = "calibri", regular = "calibri.ttf")
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
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black", lineheight = 0.1),
    legend.text = element_text(size = 80),
    legend.title = element_text(size = 80, color = "black"),
    legend.frame = element_rect(color = "black", linewidth = 0.2),
    legend.ticks = element_line(color = "black", linewidth = 0.2)
  )

ggsave(filename = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology_OG_glycoprotein_common_only_heatmap.png", plot = gene_ontology_OG_glycoprotein_common_only_heatmap, height = 5, width = 6, units = c("in"), dpi = 600)
