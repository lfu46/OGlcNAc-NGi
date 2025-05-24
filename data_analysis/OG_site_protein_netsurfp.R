#import packages
packages_names <- c("tidyverse", "rstatix", "showtext", "ggpubr", "readxl", "writexl")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
OG_glycopeptide_localized_glycoprotein_HepG2 <- OG_glycopeptide_Top_tb_HepG2 |> 
  filter(!is.na(combined_site)) |> 
  distinct(UniprotID)

OG_glycopeptide_localized_glycoprotein_HEK293T <- OG_glycopeptide_Top_tb_HEK293T |> 
  filter(!is.na(combined_site)) |> 
  distinct(UniprotID)

OG_glycopeptide_localized_glycoprotein_Jurkat <- OG_glycopeptide_Top_tb_Jurkat |> 
  filter(!is.na(combined_site)) |> 
  distinct(UniprotID)

OG_glycopeptide_localized_glycoprotein_total <- bind_rows(
  OG_glycopeptide_localized_glycoprotein_HepG2,
  OG_glycopeptide_localized_glycoprotein_HEK293T,
  OG_glycopeptide_localized_glycoprotein_Jurkat
) |> distinct()

write_xlsx(OG_glycopeptide_localized_glycoprotein_total, path = paste0(file_path, "OG_glycopeptide_localized_glycoprotein_total.xlsx"))

#import data
#secondary structure & buried exposed
OG_glycoprotein_localized_glycoprotein_secondary_structure_1 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_1.txt", delim = "\t",
  col_names = c("Class Assignment",
                "Amino Acid",
                "Sequence Name",
                "Amino Acid Number",
                "Relative Surface Accessibility - RSA",
                "Absolute Surface Accessibility",
                "Not Used",
                "Probability for Alpha - Helix",
                "Probability for Beta - Strand",
                "Probability for Coil"), name_repair = "universal", skip = 21)

OG_glycoprotein_localized_glycoprotein_secondary_structure_2 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_2.txt", delim = "\t",
  col_names = c("Class Assignment",
                "Amino Acid",
                "Sequence Name",
                "Amino Acid Number",
                "Relative Surface Accessibility - RSA",
                "Absolute Surface Accessibility",
                "Not Used",
                "Probability for Alpha - Helix",
                "Probability for Beta - Strand",
                "Probability for Coil"), name_repair = "universal", skip = 21)

OG_glycoprotein_localized_glycoprotein_secondary_structure_3 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_3.txt", delim = "\t",
  col_names = c("Class Assignment",
                "Amino Acid",
                "Sequence Name",
                "Amino Acid Number",
                "Relative Surface Accessibility - RSA",
                "Absolute Surface Accessibility",
                "Not Used",
                "Probability for Alpha - Helix",
                "Probability for Beta - Strand",
                "Probability for Coil"), name_repair = "universal", skip = 21)

OG_glycoprotein_localized_glycoprotein_secondary_structure_4 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_4.txt", delim = "\t",
  col_names = c("Class Assignment",
                "Amino Acid",
                "Sequence Name",
                "Amino Acid Number",
                "Relative Surface Accessibility - RSA",
                "Absolute Surface Accessibility",
                "Not Used",
                "Probability for Alpha - Helix",
                "Probability for Beta - Strand",
                "Probability for Coil"), name_repair = "universal", skip = 21)

OG_glycoprotein_localized_glycoprotein_secondary_structure_5 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_5.txt", delim = "\t",
  col_names = c("Class Assignment",
                "Amino Acid",
                "Sequence Name",
                "Amino Acid Number",
                "Relative Surface Accessibility - RSA",
                "Absolute Surface Accessibility",
                "Not Used",
                "Probability for Alpha - Helix",
                "Probability for Beta - Strand",
                "Probability for Coil"), name_repair = "universal", skip = 21)

OG_glycoprotein_localized_glycoprotein_secondary_structure_6 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_6.txt", delim = "\t",
  col_names = c("Class Assignment",
                "Amino Acid",
                "Sequence Name",
                "Amino Acid Number",
                "Relative Surface Accessibility - RSA",
                "Absolute Surface Accessibility",
                "Not Used",
                "Probability for Alpha - Helix",
                "Probability for Beta - Strand",
                "Probability for Coil"), name_repair = "universal", skip = 21)

OG_glycoprotein_localized_glycoprotein_secondary_structure_7 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_7.txt", delim = "\t",
  col_names = c("Class Assignment",
                "Amino Acid",
                "Sequence Name",
                "Amino Acid Number",
                "Relative Surface Accessibility - RSA",
                "Absolute Surface Accessibility",
                "Not Used",
                "Probability for Alpha - Helix",
                "Probability for Beta - Strand",
                "Probability for Coil"), name_repair = "universal", skip = 21)

OG_glycoprotein_localized_glycoprotein_secondary_structure_8 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_8.txt", delim = "\t",
  col_names = c("Class Assignment",
                "Amino Acid",
                "Sequence Name",
                "Amino Acid Number",
                "Relative Surface Accessibility - RSA",
                "Absolute Surface Accessibility",
                "Not Used",
                "Probability for Alpha - Helix",
                "Probability for Beta - Strand",
                "Probability for Coil"), name_repair = "universal", skip = 21)

OG_glycoprotein_localized_glycoprotein_secondary_structure_9 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_9.txt", delim = "\t",
  col_names = c("Class Assignment",
                "Amino Acid",
                "Sequence Name",
                "Amino Acid Number",
                "Relative Surface Accessibility - RSA",
                "Absolute Surface Accessibility",
                "Not Used",
                "Probability for Alpha - Helix",
                "Probability for Beta - Strand",
                "Probability for Coil"), name_repair = "universal", skip = 21)

OG_glycoprotein_localized_glycoprotein_secondary_structure_10 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_10.txt", delim = "\t",
  col_names = c("Class Assignment",
                "Amino Acid",
                "Sequence Name",
                "Amino Acid Number",
                "Relative Surface Accessibility - RSA",
                "Absolute Surface Accessibility",
                "Not Used",
                "Probability for Alpha - Helix",
                "Probability for Beta - Strand",
                "Probability for Coil"), name_repair = "universal", skip = 21)

OG_glycoprotein_localized_glycoprotein_secondary_structure_11 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_11.txt", delim = "\t",
  col_names = c("Class Assignment",
                "Amino Acid",
                "Sequence Name",
                "Amino Acid Number",
                "Relative Surface Accessibility - RSA",
                "Absolute Surface Accessibility",
                "Not Used",
                "Probability for Alpha - Helix",
                "Probability for Beta - Strand",
                "Probability for Coil"), name_repair = "universal", skip = 21)

OG_glycoprotein_localized_glycoprotein_secondary_structure_12 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_12.txt", delim = "\t",
  col_names = c("Class Assignment",
                "Amino Acid",
                "Sequence Name",
                "Amino Acid Number",
                "Relative Surface Accessibility - RSA",
                "Absolute Surface Accessibility",
                "Not Used",
                "Probability for Alpha - Helix",
                "Probability for Beta - Strand",
                "Probability for Coil"), name_repair = "universal", skip = 21)

OG_glycoprotein_localized_glycoprotein_secondary_structure_13 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\netsurfp_new\\OG_glycopeptide_localized_glycoprotein_13.txt", delim = "\t",
  col_names = c("Class Assignment",
                "Amino Acid",
                "Sequence Name",
                "Amino Acid Number",
                "Relative Surface Accessibility - RSA",
                "Absolute Surface Accessibility",
                "Not Used",
                "Probability for Alpha - Helix",
                "Probability for Beta - Strand",
                "Probability for Coil"), name_repair = "universal", skip = 21)

#combine
OG_glycoprotein_localized_glycoprotein_secondary_structure_combined <- bind_rows(
    OG_glycoprotein_localized_glycoprotein_secondary_structure_1,
    OG_glycoprotein_localized_glycoprotein_secondary_structure_2,
    OG_glycoprotein_localized_glycoprotein_secondary_structure_3,
    OG_glycoprotein_localized_glycoprotein_secondary_structure_4,
    OG_glycoprotein_localized_glycoprotein_secondary_structure_5,
    OG_glycoprotein_localized_glycoprotein_secondary_structure_6,
    OG_glycoprotein_localized_glycoprotein_secondary_structure_7,
    OG_glycoprotein_localized_glycoprotein_secondary_structure_8,
    OG_glycoprotein_localized_glycoprotein_secondary_structure_9,
    OG_glycoprotein_localized_glycoprotein_secondary_structure_10,
    OG_glycoprotein_localized_glycoprotein_secondary_structure_11,
    OG_glycoprotein_localized_glycoprotein_secondary_structure_12,
    OG_glycoprotein_localized_glycoprotein_secondary_structure_13
  )

OG_glycoprotein_localized_glycoprotein_secondary_structure_combined <- OG_glycoprotein_localized_glycoprotein_secondary_structure_combined |> 
  mutate(Amino.Acid.Number = ifelse(str_detect(Sequence.Name, "Q09666_2"), Amino.Acid.Number + 3000, Amino.Acid.Number)) |> 
  mutate(Amino.Acid.Number = ifelse(str_detect(Sequence.Name, "Q8IVF2_2"), Amino.Acid.Number + 3000, Amino.Acid.Number)) |> 
  mutate(Amino.Acid.Number = ifelse(str_detect(Sequence.Name, "O14686_2"), Amino.Acid.Number + 3000, Amino.Acid.Number)) |> 
  mutate(Amino.Acid.Number = ifelse(str_detect(Sequence.Name, "Q9Y6V0_2"), Amino.Acid.Number + 3000, Amino.Acid.Number))

OG_glycoprotein_localized_glycoprotein_secondary_structure_combined_tb <- OG_glycoprotein_localized_glycoprotein_secondary_structure_combined |> 
  mutate(ProteinID = str_extract(Sequence.Name, "(?<=^...).{6}"),
         Site = paste0(Amino.Acid, Amino.Acid.Number)) |> 
  mutate(Protein_Site = paste(ProteinID, Site, sep = "_")) |> 
  select(Protein_Site, Class.Assignment, starts_with("Probability"))

#import single site data frame
OG_glycopeptide_Top_tb_HepG2_singlesite <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycopeptide_Top_tb_HepG2_singlesite.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)
OG_glycopeptide_Top_tb_HEK293T_singlesite <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycopeptide_Top_tb_HEK293T_singlesite.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)
OG_glycopeptide_Top_tb_Jurkat_singlesite <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\OG_glycopeptide_Top_tb_Jurkat_singlesite.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)

#generate data frame
OG_glycopeptide_Top_tb_singlesite_secondary_structure_HepG2 <- OG_site_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  mutate(Protein_Site = paste(UniprotID, combined_site, sep = "_"),
         Cell = "HepG2") |> 
  select(Index, Protein_Site, UniprotID, combined_site, logFC, adj.P.Val, Cell) |> 
  left_join(OG_glycoprotein_localized_glycoprotein_secondary_structure_combined_tb, by = join_by(Protein_Site == Protein_Site))

OG_glycopeptide_Top_tb_singlesite_secondary_structure_HEK293T <- OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  mutate(Protein_Site = paste(UniprotID, combined_site, sep = "_"),
         Cell = "HEK293T") |> 
  select(Index, Protein_Site, UniprotID, combined_site, logFC, adj.P.Val, Cell) |> 
  left_join(OG_glycoprotein_localized_glycoprotein_secondary_structure_combined_tb, by = join_by(Protein_Site == Protein_Site))

OG_glycopeptide_Top_tb_singlesite_secondary_structure_Jurkat <- OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  mutate(Protein_Site = paste(UniprotID, combined_site, sep = "_"),
         Cell = "Jurkat") |> 
  select(Index, Protein_Site, UniprotID, combined_site, logFC, adj.P.Val, Cell) |> 
  left_join(OG_glycoprotein_localized_glycoprotein_secondary_structure_combined_tb, by = join_by(Protein_Site == Protein_Site))

#combine
OG_glycopeptide_localized_singlesite_secondary_structure_combined <- bind_rows(
  OG_glycopeptide_Top_tb_singlesite_secondary_structure_HepG2,
  OG_glycopeptide_Top_tb_singlesite_secondary_structure_HEK293T,
  OG_glycopeptide_Top_tb_singlesite_secondary_structure_Jurkat
)

#generate data frame
OG_glycopeptide_localized_singlesite_secondary_structure_combined <- OG_glycopeptide_localized_singlesite_secondary_structure_combined |> 
  mutate(
    Secondary = case_when(
      Probability.for.Alpha...Helix > Probability.for.Beta...Strand &  Probability.for.Alpha...Helix > Probability.for.Coil ~ "Helix",
      Probability.for.Beta...Strand > Probability.for.Alpha...Helix &  Probability.for.Beta...Strand > Probability.for.Coil ~ "Strand",
      Probability.for.Coil > Probability.for.Alpha...Helix &  Probability.for.Coil > Probability.for.Beta...Strand ~ "Coil"
    )
  ) |> 
  select(!starts_with("Probability"))

OG_glycopeptide_localized_singlesite_secondary_structure_combined <- OG_glycopeptide_localized_singlesite_secondary_structure_combined |> 
  filter(!is.na(Secondary))
  
write_xlsx(OG_glycopeptide_localized_singlesite_secondary_structure_combined, paste0(file_path, "OG_glycopeptide_localized_singlesite_secondary_structure_combined.xlsx"))

#wilcox test
OG_glycosite_secondary_structure_combined_wilcox_test <- OG_glycopeptide_localized_singlesite_secondary_structure_combined |> 
  group_by(Cell) |> 
  wilcox_test(logFC ~ Secondary, p.adjust.method = "BH") |> 
  add_significance("p")

#violin point plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_point_OG_site_secondary <- OG_glycopeptide_localized_singlesite_secondary_structure_combined |> 
  ggplot() +
  geom_violin(aes(x = factor(Secondary, levels = c("Helix", "Strand", "Coil")), y = logFC, fill = Cell), 
              color = "transparent") +
  geom_boxplot(aes(x = factor(Secondary, levels = c("Helix", "Strand", "Coil")), y = logFC), color = "black", width = 0.2, 
               outliers = FALSE) +
  labs(x = "", y = expression(Δlog[2]*"FC")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_grid(~ Cell, scales = "free_x") +
  stat_pvalue_manual(data = OG_glycosite_secondary_structure_combined_wilcox_test, 
                     label = "p.signif", tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.8)) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1, lineheight = 0.1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "violin_point_OG_site_secondary.png"), plot = violin_point_OG_site_secondary , height = 3, width = 6, units = "in", dpi = 600)

#buried & exposed
#wilcox test
OG_glycosite_BE_combined_wilcox_test <- OG_glycopeptide_localized_singlesite_secondary_structure_combined |> 
  group_by(Cell) |> 
  wilcox_test(logFC ~ Class.Assignment, p.adjust.method = "BH") |> 
  add_significance("p")

#violin point plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_point_OG_glycosite_BE <- OG_glycopeptide_localized_singlesite_secondary_structure_combined |> 
  ggplot() +
  geom_violin(aes(x = factor(Class.Assignment, levels = c("E", "B")), y = logFC, fill = Cell), 
              color = "transparent") +
  geom_boxplot(aes(x = factor(Class.Assignment, levels = c("E", "B")), y = logFC), color = "black", width = 0.2, 
               outliers = FALSE) +
  labs(x = "", y = expression(Δlog[2]*"FC")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  scale_x_discrete(labels = c("E" = "Exposed",
                              "B" = "Buried")) +
  facet_grid(~ Cell, scales = "free_x") +
  stat_pvalue_manual(data = OG_glycosite_BE_combined_wilcox_test, label = "p.signif", 
                     tip.length = 0, size = 50,
                     hide.ns = "p", y.position = c(1.8)) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1, lineheight = 0.1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "violin_point_OG_glycosite_BE.png"), plot = violin_point_OG_glycosite_BE, height = 3, width = 6, units = "in", dpi = 600)

#order & disorder
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
OG_glycosite_disorder_HepG2 <- OG_site_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  mutate(Protein_site = paste(UniprotID, combined_site, sep = "_"),
         Cell = "HepG2") |> 
  select(Index, Protein_site, UniprotID, combined_site, logFC, adj.P.Val, Cell) |> 
  left_join(OG_glycosite_disorder_combined, by = "Protein_site")

#HEK293T
OG_glycosite_disorder_HEK293T <- OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  mutate(Protein_site = paste(UniprotID, combined_site, sep = "_"),
         Cell = "HEK293T") |> 
  select(Index, Protein_site, UniprotID, combined_site, logFC, adj.P.Val, Cell) |> 
  left_join(OG_glycosite_disorder_combined, by = "Protein_site")

#Jurkat
OG_glycosite_disorder_Jurkat <- OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  mutate(Protein_site = paste(UniprotID, combined_site, sep = "_"),
         Cell = "Jurkat") |> 
  select(Index, Protein_site, UniprotID, combined_site, logFC, adj.P.Val, Cell) |> 
  left_join(OG_glycosite_disorder_combined, by = "Protein_site")

#combined
OG_glycosite_disorder_tb <- bind_rows(
  OG_glycosite_disorder_HepG2,
  OG_glycosite_disorder_HEK293T,
  OG_glycosite_disorder_Jurkat
)

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
  labs(x = "", y = expression(Δlog[2]*"FC")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_grid(~ Cell, scales = "free_x") +
  stat_pvalue_manual(data = OG_glycosite_disorder_wilcox_test, label = "p.signif", tip.length = 0, 
                     size = 50,
                     hide.ns = "p", y.position = c(1.8)) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1, lineheight = 0.1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 80, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "violin_point_OG_glycosite_disorder.png"), 
       plot = violin_point_OG_glycosite_disorder, height = 3, width = 6, units = "in", dpi = 600)
