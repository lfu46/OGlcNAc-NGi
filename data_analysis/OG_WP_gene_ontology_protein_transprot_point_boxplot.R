#import packages
packages_names <- c("tidyverse", "showtext", "readxl", "ggpubr", "rstatix")
lapply(packages_names, require, character.only = TRUE)

#import data
#HepG2
gene_ontology_OG_glycoprotein_BP_proteintranport_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_BP_protein_transport.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "HepG2"
) |> select(UniprotID = UNIPROT_ACCESSION)

#HEK293T
gene_ontology_OG_glycoprotein_BP_proteintranport_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_BP_protein_transport.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "HEK293T"
) |> select(UniprotID = UNIPROT_ACCESSION)

#Jurkat
gene_ontology_OG_glycoprotein_BP_proteintranport_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_BP_protein_transport.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "Jurkat"
) |> select(UniprotID = UNIPROT_ACCESSION)

#generate data frame
#HepG2
gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_HepG2 <- gene_ontology_OG_glycoprotein_BP_proteintranport_HepG2 |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = "UniprotID") |> 
  select(UniprotID, OG = logFC) |> 
  left_join(WP_protein_Top_tb_HepG2, by = "UniprotID") |> 
  filter(!is.na(logFC)) |> 
  select(UniprotID, OG, WP = logFC) |> 
  pivot_longer(cols = OG:WP, names_to = 'Group', values_to = 'logFC') |> 
  mutate(Cell = "HepG2")

OG_WP_ks_test_HepG2 <- ks.test(
  gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_HepG2 |> filter(Group == "OG") |> pull(logFC),
  gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_HepG2 |> filter(Group == "WP") |> pull(logFC)
)

#HEK293T
gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_HEK293T <- gene_ontology_OG_glycoprotein_BP_proteintranport_HEK293T |> 
  left_join(OG_glycoprotein_Top_tb_HEK293T, by = "UniprotID") |> 
  select(UniprotID, OG = logFC) |> 
  left_join(WP_protein_Top_tb_HEK293T, by = "UniprotID") |> 
  filter(!is.na(logFC)) |> 
  select(UniprotID, OG, WP = logFC) |> 
  pivot_longer(cols = OG:WP, names_to = 'Group', values_to = 'logFC') |> 
  mutate(Cell = "HEK293T")

OG_WP_ks_test_HEK293T <- ks.test(
  gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_HEK293T |> filter(Group == "OG") |> pull(logFC),
  gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_HEK293T |> filter(Group == "WP") |> pull(logFC)
)

#Jurkat
gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_Jurkat <- gene_ontology_OG_glycoprotein_BP_proteintranport_Jurkat |> 
  left_join(OG_glycoprotein_Top_tb_Jurkat, by = "UniprotID") |> 
  select(UniprotID, OG = logFC) |> 
  left_join(WP_protein_Top_tb_Jurkat, by = "UniprotID") |> 
  filter(!is.na(logFC)) |> 
  select(UniprotID, OG, WP = logFC) |> 
  pivot_longer(cols = OG:WP, names_to = 'Group', values_to = 'logFC') |> 
  mutate(Cell = "Jurkat")

OG_WP_ks_test_Jurkat <- ks.test(
  gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_Jurkat |> filter(Group == "OG") |> pull(logFC),
  gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_Jurkat |> filter(Group == "WP") |> pull(logFC)
)

#combine
gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_total <- bind_rows(
  gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_HepG2,
  gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_HEK293T,
  gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_Jurkat
)

#wilcoxt test
gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_wilcox_test <- 
  gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_total |> 
  group_by(Cell) |> 
  wilcox_test(logFC ~ Group, p.adjust.method = "BH") |> 
  add_significance("p")

#ks test
gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_ks_test <- 
  gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_wilcox_test |> 
  mutate(p = ifelse(Cell == "HEK293T", 0.111, ifelse(Cell == "HepG2", 0.007656, 0.03393))) |> 
  mutate(p.signif = ifelse(Cell == "HEK293T", "ns", ifelse(Cell == "HepG2", "**", "*")))

#point plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

point_plot_gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_total <- gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_total |> 
  mutate(Group_Color = ifelse(Group == "OG", paste(Group, Cell, sep = "_"), "WP")) |> 
  ggplot(aes(x = factor(Group, levels = c("WP", "OG")), y = logFC)) +
  geom_line(aes(group = UniprotID), color = "gray") +
  geom_boxplot(outliers = FALSE) +
  geom_point(aes(fill = Group_Color), shape = 21, color = "black", size = 2) +
  scale_fill_manual(values = c(
    "OG_HEK293T" = Color_2,
    "OG_HepG2" = Color_3,
    "OG_Jurkat" = Color_4,
    "WP" = "gray"
  )) +
  stat_pvalue_manual(data = gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_ks_test, 
                     label = "p.signif", tip.length = 0, size = 50, hide.ns = "p",
                     y.position = c(1.8)) +
  facet_grid(~ Cell, scale = "free_x") +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)"), fill = "") +
  ggtitle("Protein transport") +
  theme_bw() +
  theme(
    title = element_text(size = 100, color = "black"),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    axis.ticks.length.x = unit(0, "in"),
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = NA, color = NA),
    legend.position = "none",
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "point_plot_gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_total.png"), 
       plot = point_plot_gene_ontology_OG_glycoprotein_WP_BP_proteintranport_logFC_total, 
       height = 4, width = 4, units = c("in"), dpi = 600)
