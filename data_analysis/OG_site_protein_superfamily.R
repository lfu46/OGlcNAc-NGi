#import packages
packages_names <- c("tidyverse", "rstatix", "showtext", "ggpubr", "writexl")
lapply(packages_names, require, character.only = TRUE)

#import data
OG_glycopeptide_superfamily_1 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_1.txt",
  col_names = FALSE
)
OG_glycopeptide_superfamily_2 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_2.txt",
  col_names = FALSE
)
OG_glycopeptide_superfamily_3 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_3.txt",
  delim = " ",
  col_names = FALSE
)
OG_glycopeptide_superfamily_4 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_4.txt",
  delim = " ",
  col_names = FALSE
)
OG_glycopeptide_superfamily_5 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_5.txt",
  delim = " ",
  col_names = FALSE
)
OG_glycopeptide_superfamily_6 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_6.txt",
  delim = " ",
  col_names = FALSE
)
OG_glycopeptide_superfamily_7 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_7.txt",
  delim = " ",
  col_names = FALSE
)
OG_glycopeptide_superfamily_8 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_8.txt",
  delim = " ",
  col_names = FALSE
)
OG_glycopeptide_superfamily_9 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_9.txt",
  delim = " ",
  col_names = FALSE
)
OG_glycopeptide_superfamily_10 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_10.txt",
  delim = " ",
  col_names = FALSE
)
OG_glycopeptide_superfamily_11 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_11.txt",
  delim = " ",
  col_names = FALSE
)
OG_glycopeptide_superfamily_12 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_12.txt",
  delim = " ",
  col_names = FALSE
)
OG_glycopeptide_superfamily_13 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_13.txt",
  delim = " ",
  col_names = FALSE
)
OG_glycopeptide_superfamily_14 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_14.txt",
  col_names = FALSE
)
OG_glycopeptide_superfamily_15 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_15.txt",
  delim = "\t",
  col_names = FALSE
)
OG_glycopeptide_superfamily_16 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_16.txt",
  delim = "\t",
  col_names = FALSE
)
OG_glycopeptide_superfamily_17 <- read_delim(
  file = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\superfamily\\OG_glycopeptide_localized_glycoprotein_total_17.txt",
  delim = "\t",
  col_names = FALSE
)

#combine
OG_glycopeptide_superfamily_total <- bind_rows(
  OG_glycopeptide_superfamily_1,
  OG_glycopeptide_superfamily_2,
  OG_glycopeptide_superfamily_3,
  OG_glycopeptide_superfamily_4,
  OG_glycopeptide_superfamily_5,
  OG_glycopeptide_superfamily_6,
  OG_glycopeptide_superfamily_7,
  OG_glycopeptide_superfamily_8,
  OG_glycopeptide_superfamily_9,
  OG_glycopeptide_superfamily_10,
  OG_glycopeptide_superfamily_11,
  OG_glycopeptide_superfamily_12,
  OG_glycopeptide_superfamily_13,
  OG_glycopeptide_superfamily_14,
  OG_glycopeptide_superfamily_15,
  OG_glycopeptide_superfamily_16,
  OG_glycopeptide_superfamily_17
)

#generate data frame
OG_glycopeptide_superfamily_total_adj <- OG_glycopeptide_superfamily_total |> 
  separate_rows(X3, sep = ",") |> 
  separate(X1, sep = "\\|", into = c("sp", "UniprotID", "HUMAN")) |> 
  mutate(Region = X3) |> 
  separate(Region, sep = "-", into = c("Start", "End"), remove = FALSE) |> 
  mutate_at(vars(Start, End), as.numeric) |> 
  select(UniprotID, Region, Start, End) |> 
  mutate(Start_adj = Start - 10,
         End_adj = End + 10) |> 
  select(UniprotID, Region, Start_adj, End_adj) |> 
  mutate(Start_adj = ifelse(str_detect(UniprotID, "Q09666_2"), Start_adj + 3000, Start_adj),
         End_adj = ifelse(str_detect(UniprotID, "Q09666_2"), End_adj + 3000, End_adj)) |> 
  mutate(Start_adj = ifelse(str_detect(UniprotID, "Q8IVF2_2"), Start_adj + 3000, Start_adj),
         End_adj = ifelse(str_detect(UniprotID, "Q8IVF2_2"), End_adj + 3000, End_adj)) |> 
  mutate(Start_adj = ifelse(str_detect(UniprotID, "O14686_2"), Start_adj + 3000, Start_adj),
         End_adj = ifelse(str_detect(UniprotID, "O14686_2"), End_adj + 3000, End_adj)) |> 
  mutate(Start_adj = ifelse(str_detect(UniprotID, "Q9Y6V0_2"), Start_adj + 3000, Start_adj),
         End_adj = ifelse(str_detect(UniprotID, "Q9Y6V0_2"), End_adj + 3000, End_adj)) |> 
  mutate(UniprotID = ifelse(UniprotID %in% c("Q09666_1", "Q09666_2"), "Q09666", UniprotID)) |> 
  mutate(UniprotID = ifelse(UniprotID %in% c("Q8IVF2_1", "Q8IVF2_2"), "Q8IVF2", UniprotID)) |> 
  mutate(UniprotID = ifelse(UniprotID %in% c("O14686_1", "O14686_2"), "O14686", UniprotID)) |> 
  mutate(UniprotID = ifelse(UniprotID %in% c("Q9Y6V0_1", "Q9Y6V0_2"), "Q9Y6V0", UniprotID))

#generate data frame
OG_glycopeptide_singlesite_domain_HepG2 <- OG_site_protein_OGlcNAc_effect_Top_tb_HepG2 |> 
  mutate(site_position = str_extract(combined_site, "\\d+")) |> 
  left_join(OG_glycopeptide_superfamily_total_adj, by = "UniprotID", relationship = "many-to-many") |> 
  select(Index, UniprotID, Region, Start_adj, End_adj, combined_site, site_position, logFC, adj.P.Val) |> 
  mutate(Domain = ifelse(is.na(Region), "Not Domain", ifelse(site_position > Start_adj & site_position < End_adj, "Domain", "Not Domain"))) |> 
  mutate(Cell = "HepG2")

OG_glycopeptide_singlesite_domain_HEK293T <- OG_site_protein_OGlcNAc_effect_Top_tb_HEK293T |> 
  mutate(site_position = str_extract(combined_site, "\\d+")) |> 
  left_join(OG_glycopeptide_superfamily_total_adj, by = "UniprotID", relationship = "many-to-many") |> 
  select(Index, UniprotID, Region, Start_adj, End_adj, combined_site, site_position, logFC, adj.P.Val) |> 
  mutate(Domain = ifelse(is.na(Region), "Not Domain", ifelse(site_position > Start_adj & site_position < End_adj, "Domain", "Not Domain"))) |> 
  mutate(Cell = "HEK293T")

OG_glycopeptide_singlesite_domain_Jurkat <- OG_site_protein_OGlcNAc_effect_Top_tb_Jurkat |> 
  mutate(site_position = str_extract(combined_site, "\\d+")) |> 
  left_join(OG_glycopeptide_superfamily_total_adj, by = "UniprotID", relationship = "many-to-many") |> 
  select(Index, UniprotID, Region, Start_adj, End_adj, combined_site, site_position, logFC, adj.P.Val) |> 
  mutate(Domain = ifelse(is.na(Region), "Not Domain", ifelse(site_position > Start_adj & site_position < End_adj, "Domain", "Not Domain"))) |> 
  mutate(Cell = "Jurkat")

#combine
OG_glycopeptide_singlesite_domain_combined <- bind_rows(
  OG_glycopeptide_singlesite_domain_HepG2,
  OG_glycopeptide_singlesite_domain_HEK293T,
  OG_glycopeptide_singlesite_domain_Jurkat
)

#wilcox test
OG_glycopeptide_singlesite_domain_wilcox_test <- OG_glycopeptide_singlesite_domain_combined |> 
  group_by(Cell) |> 
  wilcox_test(logFC ~ Domain, p.adjust.method = "BH") |> 
  add_significance("p")

#violin point plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_point_OG_glycopeptide_singlesite_domain <- OG_glycopeptide_singlesite_domain_combined |> 
  ggplot() +
  geom_violin(aes(x = factor(Domain, levels = c("Domain", "Not Domain")), y = logFC, fill = Cell), 
              color = "transparent") +
  geom_boxplot(aes(x = factor(Domain, levels = c("Domain", "Not Domain")), y = logFC), color = "black",
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
  stat_pvalue_manual(data = OG_glycopeptide_singlesite_domain_wilcox_test, label = "p.signif", 
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

ggsave(filename = paste0(file_path, "violin_point_OG_glycopeptide_singlesite_domain.png"), 
       plot = violin_point_OG_glycopeptide_singlesite_domain, height = 3, width = 6, units = "in", dpi = 600)
