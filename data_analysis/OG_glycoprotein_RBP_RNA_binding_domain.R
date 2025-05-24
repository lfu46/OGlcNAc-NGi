#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "showtext", "ggpubr", "rstatix", 
                    "ComplexHeatmap", "circlize", "eulerr", "ggridges")
lapply(packages_names, require, character.only = TRUE)

#import data
#HEK293T
OG_glycoprotein_RBP_RBD_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBP_RBD_HEK293T.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(cell = "HEK293T")

#HepG2
OG_glycoprotein_RBP_RBD_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBP_RBD_HepG2.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(cell = "HepG2")

#Jurkat
OG_glycoprotein_RBP_RBD_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBP_RBD_Jurkat.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
) |> mutate(cell = "Jurkat")

#combine
OG_glycorprotein_RBP_RBD_combined <- bind_rows(
  OG_glycoprotein_RBP_RBD_HEK293T,
  OG_glycoprotein_RBP_RBD_HepG2,
  OG_glycoprotein_RBP_RBD_Jurkat
)

#dot plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

remove_domain <- c("PDZ", "PABP", "MA3", "HnRNPA1_LC", "CID", "CH")

dot_plot_OG_glycoprotein_RBP_RBD <- OG_glycorprotein_RBP_RBD_combined |> 
  filter(! Term %in% remove_domain) |> 
  ggplot() +
  geom_point(aes(x = cell, y = Term, size = Count, fill = P.Value), color = "black", shape = 21) +
  labs(x = "", y = "") +
  scale_size(range = c(4, 8)) +
  scale_fill_stepsn(limits = c(0, 1), 
                    breaks = c(0, 1E-5, 1E-4, 1E-3, 0.05, 1), 
                    labels = c("0", "1E-5", "1E-4", "1E-3", "0.05", "1"), 
                    n.breaks = 6, 
                    values = scales::rescale(c(0, 1E-5, 1E-4, 1E-3, 0.05, 1)),
                    colours = c(Color_1, NA)) +
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

ggsave(filename = paste0(file_path, "dot_plot_OG_glycoprotein_RBP_RBD.png"), plot = dot_plot_OG_glycoprotein_RBP_RBD, height = 4, width = 4, units = c("in"), dpi = 600)

#import data
#HEK293T
OG_glycoprotein_RBD_logFC_RRM_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBD_logFC_HEK293T.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = 'RRM_1'
) |> mutate(domain = "RRM_1")
OG_glycoprotein_RBD_logFC_KH_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBD_logFC_HEK293T.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = 'KH_1'
) |> mutate(domain = "KH_1")
OG_glycoprotein_RBD_logFC_SAP_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBD_logFC_HEK293T.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = 'SAP'
) |> mutate(domain = "SAP")
OG_glycoprotein_RBD_logFC_AAA_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBD_logFC_HEK293T.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = 'AAA_11'
) |> mutate(domain = "AAA_11")
OG_glycoprotein_RBD_logFC_zf_CCCH_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBD_logFC_HEK293T.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = 'zf-CCCH'
) |> mutate(domain = "zf-CCCH")

OG_glycoprotein_RBD_logFC_HEK293T_combined <- bind_rows(
  OG_glycoprotein_RBD_logFC_RRM_HEK293T,
  OG_glycoprotein_RBD_logFC_KH_HEK293T,
  OG_glycoprotein_RBD_logFC_SAP_HEK293T,
  OG_glycoprotein_RBD_logFC_AAA_HEK293T,
  OG_glycoprotein_RBD_logFC_zf_CCCH_HEK293T
) |> 
  left_join(OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HEK293T, by = join_by("UNIPROT_ACCESSION" == "UniprotID")) |> 
  select(domain, UNIPROT_ACCESSION, logFC) |> 
  mutate(cell = "HEK293T")

#wilcox test
OG_glycoprotein_RBD_logFC_HEK293T_wilcox_test <- OG_glycoprotein_RBD_logFC_HEK293T_combined |> 
  wilcox_test(logFC ~ domain, p.adjust.method = "BH") |> 
  add_significance("p")


#HepG2
OG_glycoprotein_RBD_logFC_RRM_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBD_logFC_HepG2.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = 'RRM_1'
) |> mutate(domain = "RRM_1")
OG_glycoprotein_RBD_logFC_KH_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBD_logFC_HepG2.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = 'KH_1'
) |> mutate(domain = "KH_1")
OG_glycoprotein_RBD_logFC_SAP_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBD_logFC_HepG2.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = 'SAP'
) |> mutate(domain = "SAP")
OG_glycoprotein_RBD_logFC_LIM_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBD_logFC_HepG2.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = 'LIM'
) |> mutate(domain = "LIM")
OG_glycoprotein_RBD_logFC_PDZ_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBD_logFC_HepG2.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = 'PDZ'
) |> mutate(domain = "PDZ")

OG_glycoprotein_RBD_logFC_HepG2_combined <- bind_rows(
  OG_glycoprotein_RBD_logFC_RRM_HepG2,
  OG_glycoprotein_RBD_logFC_KH_HepG2,
  OG_glycoprotein_RBD_logFC_SAP_HepG2,
  OG_glycoprotein_RBD_logFC_LIM_HepG2,
  OG_glycoprotein_RBD_logFC_PDZ_HepG2
) |> 
  left_join(OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_HepG2, by = join_by("UNIPROT_ACCESSION" == "UniprotID")) |> 
  select(domain, UNIPROT_ACCESSION, logFC) |> 
  mutate(cell = "HepG2")

#wilcox test
OG_glycoprotein_RBD_logFC_HepG2_wilcox_test <- OG_glycoprotein_RBD_logFC_HepG2_combined |> 
  wilcox_test(logFC ~ domain, p.adjust.method = "BH") |> 
  add_significance("p")

#Jurkat
OG_glycoprotein_RBD_logFC_RRM_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBD_logFC_Jurkat.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = 'RRM_1'
) |> mutate(domain = "RRM_1")
OG_glycoprotein_RBD_logFC_KH_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBD_logFC_Jurkat.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = 'KH_1'
) |> mutate(domain = "KH_1")
OG_glycoprotein_RBD_logFC_SPOC_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBD_logFC_Jurkat.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = 'SPOC'
) |> mutate(domain = "SPOC")
OG_glycoprotein_RBD_logFC_G_patch_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBD_logFC_Jurkat.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = 'G-patch'
) |> mutate(domain = "G-patch")
OG_glycoprotein_RBD_logFC_BAT2_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\human_RBP_database\\OG_glycoprotein_RBD_logFC_Jurkat.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = 'BAT2_N'
) |> mutate(domain = "BAT2_N")

OG_glycoprotein_RBD_logFC_Jurkat_combined <- bind_rows(
  OG_glycoprotein_RBD_logFC_RRM_Jurkat,
  OG_glycoprotein_RBD_logFC_KH_Jurkat,
  OG_glycoprotein_RBD_logFC_SPOC_Jurkat,
  OG_glycoprotein_RBD_logFC_G_patch_Jurkat,
  OG_glycoprotein_RBD_logFC_BAT2_Jurkat
) |> 
  left_join(OG_glycoprotein_protein_OGlcNAc_effect_Top_tb_Jurkat, by = join_by("UNIPROT_ACCESSION" == "UniprotID")) |> 
  select(domain, UNIPROT_ACCESSION, logFC) |> 
  mutate(cell = "Jurkat")

#wilcox test
OG_glycoprotein_RBD_logFC_Jurkat_wilcox_test <- OG_glycoprotein_RBD_logFC_Jurkat_combined |> 
  wilcox_test(logFC ~ domain, p.adjust.method = "BH") |> 
  add_significance("p")

#combine
OG_glycoprotein_RBD_logFC_combined <- bind_rows(
  OG_glycoprotein_RBD_logFC_HEK293T_combined,
  OG_glycoprotein_RBD_logFC_HepG2_combined,
  OG_glycoprotein_RBD_logFC_Jurkat_combined
)

#wilcox test
OG_glycoprotein_RBD_logFC_wilcox_test <- OG_glycoprotein_RBD_logFC_combined |> 
  group_by(cell) |> 
  wilcox_test(logFC ~ domain, p.adjust.method = "BH") |> 
  add_significance("p") |> 
  filter(p.signif != "ns")


#ggridges
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

ggridges_OG_glycoprotein_RBD_logFC <- OG_glycoprotein_RBD_logFC_combined |> 
  ggplot() +
  geom_density_ridges(
    aes(x = logFC, y = domain, fill = cell), color = "black",
    quantile_lines = TRUE, quantiles = 2, 
    scale = 0.55,
    jittered_points = TRUE, position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 8, point_alpha = 1
  ) +
  labs(x = expression(Δlog[2]*"FC"), y = "") +
  geom_text(data = OG_glycoprotein_RBD_logFC_wilcox_test, aes(y = group1, x = 2.0, label = p.signif), size = 70, vjust = -0.1) +
  coord_cartesian(xlim = c(-2, 2)) +
  facet_grid(cell ~ ., scale = "free_y") +
  scale_fill_manual(
    name = "",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black"),
    strip.text = element_text(size = 100, color = "black"),
    strip.background = element_rect(fill = "transparent", color = "transparent"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "ggridges_OG_glycoprotein_RBD_logFC.png"), plot = ggridges_OG_glycoprotein_RBD_logFC, height = 8, width = 4, units = "in", dpi = 600)
