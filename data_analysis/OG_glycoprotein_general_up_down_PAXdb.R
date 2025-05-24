#import packages
packages_names <- c("tidyverse", "writexl", "readxl", "ggpubr", "rstatix", "showtext")
lapply(packages_names, require, character.only = TRUE)

#import data
uniprot_OG_glycoprotein_up_string_id <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\PAXdb\\uniprot_OG_glycoprotein_up_string_id.xlsx"
)
uniprot_OG_glycoprotein_down_string_id <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\PAXdb\\uniprot_OG_glycoprotein_down_string_id.xlsx"
)

#import database
HEK293T_PAXdb <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\PAXdb\\HEK293_PAXdb.txt",
  skip = 11,
  col_names = TRUE,
  delim = "\t"
)
HepG2_PAXdb <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\PAXdb\\HepG2_PAXdb.txt",
  skip = 11,
  col_names = TRUE,
  delim = "\t"
)
Jurkat_PAXdb <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\PAXdb\\Jurkat_PAXdb.txt",
  skip = 11,
  col_names = TRUE,
  delim = "\t"
)

#generate data frame
#up HEK293T
OG_glycoprotein_PAXdb_up_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  left_join(uniprot_OG_glycoprotein_up_string_id, by = join_by(UniprotID == From)) |> 
  select(UniprotID, String_id = To) |> 
  filter(!is.na(String_id)) |> 
  left_join(HEK293T_PAXdb, by = join_by(String_id == '#string_external_id')) |> 
  filter(!is.na(abundance)) |> 
  mutate(cell = "HEK293T", category = "up")

#up HepG2
OG_glycoprotein_PAXdb_up_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  left_join(uniprot_OG_glycoprotein_up_string_id, by = join_by(UniprotID == From)) |> 
  select(UniprotID, String_id = To) |> 
  filter(!is.na(String_id)) |> 
  left_join(HepG2_PAXdb, by = join_by(String_id == '#string_external_id')) |> 
  filter(!is.na(abundance)) |> 
  mutate(cell = "HepG2", category = "up")

#up Jurkat
OG_glycoprotein_PAXdb_up_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  left_join(uniprot_OG_glycoprotein_up_string_id, by = join_by(UniprotID == From)) |> 
  select(UniprotID, String_id = To) |> 
  filter(!is.na(String_id)) |> 
  left_join(Jurkat_PAXdb, by = join_by(String_id == '#string_external_id')) |> 
  filter(!is.na(abundance)) |> 
  mutate(cell = "Jurkat", category = "up")

#down HEK293T
OG_glycoprotein_PAXdb_down_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  left_join(uniprot_OG_glycoprotein_down_string_id, by = join_by(UniprotID == From)) |> 
  select(UniprotID, String_id = To) |> 
  filter(!is.na(String_id)) |> 
  left_join(HEK293T_PAXdb, by = join_by(String_id == '#string_external_id')) |> 
  filter(!is.na(abundance)) |> 
  mutate(cell = "HEK293T", category = "down")

#down HepG2
OG_glycoprotein_PAXdb_down_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  left_join(uniprot_OG_glycoprotein_down_string_id, by = join_by(UniprotID == From)) |> 
  select(UniprotID, String_id = To) |> 
  filter(!is.na(String_id)) |> 
  left_join(HepG2_PAXdb, by = join_by(String_id == '#string_external_id')) |> 
  filter(!is.na(abundance)) |> 
  mutate(cell = "HepG2", category = "down")

#down Jurkat
OG_glycoprotein_PAXdb_down_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  left_join(uniprot_OG_glycoprotein_down_string_id, by = join_by(UniprotID == From)) |> 
  select(UniprotID, String_id = To) |> 
  filter(!is.na(String_id)) |> 
  left_join(Jurkat_PAXdb, by = join_by(String_id == '#string_external_id')) |> 
  filter(!is.na(abundance)) |> 
  mutate(cell = "Jurkat", category = "down")

#combine
OG_glycoprotein_PAXdb_combined <- bind_rows(
  OG_glycoprotein_PAXdb_up_HEK293T,
  OG_glycoprotein_PAXdb_up_HepG2,
  OG_glycoprotein_PAXdb_up_Jurkat,
  OG_glycoprotein_PAXdb_down_HEK293T,
  OG_glycoprotein_PAXdb_down_HepG2,
  OG_glycoprotein_PAXdb_down_Jurkat
) |> mutate(log2_abundance = log2(abundance))

#wilcox test
OG_glycoprotein_PAXdb_wilcox_test <- OG_glycoprotein_PAXdb_combined |> 
  group_by(cell) |> 
  wilcox_test(log2_abundance ~ category, p.adjust.method = "BH") |> 
  add_significance("p")

#violin boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

violin_point_OG_glycoprotein_PAXdb <- OG_glycoprotein_PAXdb_combined |> 
  ggplot() +
  geom_violin(aes(x = factor(category, levels = c("up", "down")), y = log2_abundance, fill =  cell), 
              color = "transparent") +
  geom_boxplot(aes(x = factor(category, levels = c("up", "down")), y = log2_abundance), color = "black",
               outliers = FALSE, width = 0.2) +
  labs(x = "", y = expression(log[2]*"(Abundance)")) +
  scale_fill_manual(
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )
  ) +
  facet_grid(~ cell, scales = "free_x") +
  stat_pvalue_manual(data = OG_glycoprotein_PAXdb_wilcox_test, label = "p.signif", tip.length = 0, 
                     size = 50,
                     hide.ns = "p", y.position = c(8.8)) +
  coord_cartesian(ylim = c(5, 9)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.position = "none",
    strip.text = element_text(size = 100, color = "black"),
    strip.background = element_rect(fill = NA, color = NA)
  )

ggsave(filename = paste0(file_path, "violin_point_OG_glycoprotein_PAXdb.png"), 
       plot = violin_point_OG_glycoprotein_PAXdb, height = 3, width = 4, units = "in", dpi = 600)
