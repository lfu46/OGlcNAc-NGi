#import packages
packages_names <- c("tidyverse", "writexl", "readxl", "ggpubr", "rstatix", "showtext")
lapply(packages_names, require, character.only = TRUE)

#import data
transcription_factor_uniprotID <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\transcription_factor\\Transcription_Factor_Uniprot.xlsx",
  col_name = TRUE,
  .name_repair = "universal"
)

#generate data frame
#HepG2
OG_glycoprotein_transcription_factor_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(UniprotID %in% transcription_factor_uniprotID$Entry) |> 
  select(UniprotID, OG = logFC) |> 
  left_join(WP_protein_Top_tb_HepG2, by = "UniprotID") |> 
  filter(!is.na(logFC)) |> 
  select(UniprotID, OG, WP = logFC) |> 
  pivot_longer(cols = OG:WP, names_to = 'Group', values_to = 'logFC') |> 
  mutate(Cell = "HepG2")

OG_WP_ks_test_HepG2 <- ks.test(
  OG_glycoprotein_transcription_factor_HepG2 |> filter(Group == "OG") |> pull(logFC),
  OG_glycoprotein_transcription_factor_HepG2 |> filter(Group == "WP") |> pull(logFC)
)

#HEK293T
OG_glycoprotein_transcription_factor_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(UniprotID %in% transcription_factor_uniprotID$Entry) |> 
  select(UniprotID, OG = logFC) |> 
  left_join(WP_protein_Top_tb_HEK293T, by = "UniprotID") |> 
  filter(!is.na(logFC)) |> 
  select(UniprotID, OG, WP = logFC) |> 
  pivot_longer(cols = OG:WP, names_to = 'Group', values_to = 'logFC') |> 
  mutate(Cell = "HEK293T")

OG_WP_ks_test_HEK293T <- ks.test(
  OG_glycoprotein_transcription_factor_HEK293T |> filter(Group == "OG") |> pull(logFC),
  OG_glycoprotein_transcription_factor_HEK293T |> filter(Group == "WP") |> pull(logFC)
)

#Jurkat
OG_glycoprotein_transcription_factor_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(UniprotID %in% transcription_factor_uniprotID$Entry) |> 
  select(UniprotID, OG = logFC) |> 
  left_join(WP_protein_Top_tb_Jurkat, by = "UniprotID") |> 
  filter(!is.na(logFC)) |> 
  select(UniprotID, OG, WP = logFC) |> 
  pivot_longer(cols = OG:WP, names_to = 'Group', values_to = 'logFC') |> 
  mutate(Cell = "Jurkat")

OG_WP_ks_test_Jurkat <- ks.test(
  OG_glycoprotein_transcription_factor_Jurkat |> filter(Group == "OG") |> pull(logFC),
  OG_glycoprotein_transcription_factor_Jurkat |> filter(Group == "WP") |> pull(logFC)
)

#combine
OG_glycoprotein_transcription_factor_logFC_total <- bind_rows(
  OG_glycoprotein_transcription_factor_HepG2,
  OG_glycoprotein_transcription_factor_HEK293T,
  OG_glycoprotein_transcription_factor_Jurkat
)

#t test
OG_glycoprotein_transcription_factor_logFC_t_test <- OG_glycoprotein_transcription_factor_logFC_total |> 
  group_by(Cell) |> 
  sign_test(logFC ~ Group, p.adjust.method = "BH") |> 
  add_significance("p")

#ks test
OG_glycoprotein_transcription_factor_logFC_ks_test <- OG_glycoprotein_transcription_factor_logFC_t_test |> 
  mutate(p = ifelse(Cell == "HEK293T", 0.7557, ifelse(Cell == "HepG2", 0.0005823, 0.4999))) |> 
  mutate(p.signif = ifelse(Cell == "HEK293T", "ns", ifelse(Cell == "HepG2", "***", "ns")))

#point plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

point_plot_OG_glycoprotein_transcription_factor_logFC_total <- OG_glycoprotein_transcription_factor_logFC_total |> 
  mutate(Group_Color = ifelse(Group == "OG", paste(Group, Cell, sep = "_"), "WP")) |> 
  ggplot(aes(x = Group, y = logFC)) +
  geom_line(aes(group = UniprotID), color = "gray") +
  geom_point(aes(fill = Group_Color), shape = 21, color = "black") +
  scale_fill_manual(values = c(
    "OG_HEK293T" = Color_2,
    "OG_HepG2" = Color_3,
    "OG_Jurkat" = Color_4,
    "WP" = "gray"
  )) +
  stat_pvalue_manual(data = OG_glycoprotein_transcription_factor_logFC_ks_test, label = "p.signif", tip.length = 0, size = 30,
                     y.position = c(2.7)) +
  facet_grid(~ Cell, scale = "free_x") +
  coord_cartesian(ylim = c(-3, 3)) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)"), fill = "") +
  ggtitle("Transcription factor") +
  theme_bw() +
  theme(
    title = element_text(size = 100, color = "black"),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    axis.ticks.length.x = unit(0, "in"),
    strip.text = element_text(size = 80, color = "black", margin = margin(b = 0.01, t = 0, unit = "line")),
    strip.background = element_rect(fill = NA, color = NA),
    legend.position = "none",
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "point_plot_OG_glycoprotein_transcription_factor_logFC_total.png"), plot = point_plot_OG_glycoprotein_transcription_factor_logFC_total, height = 4, width = 4, units = c("in"), dpi = 600)
