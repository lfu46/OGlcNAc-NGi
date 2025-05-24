#import packages
packages_names <- c("tidyverse", "rstatix", "ggpubr", "showtext")
lapply(packages_names, require, character.only = TRUE)

#gravy score
#import data
gravy_exosome_HepG2 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\gravy\\GRAVY_Exosome_HepG2.txt",
  delim = ";",
  col_names = TRUE,
  name_repair = "universal"
)

gravy_exosome_nucleus_HepG2 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\gravy\\GRAVY_Exosome_Nucleus_HepG2.txt",
  delim = ";",
  col_names = TRUE,
  name_repair = "universal"
)

gravy_nucleus_HepG2 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\gravy\\GRAVY_Nucleus_HepG2.txt",
  delim = ";",
  col_names = TRUE,
  name_repair = "universal"
)

gravy_nucleus_not_exosome_HepG2 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\gravy\\GRAVY_Nucleus_not_Exosome_HepG2.txt",
  delim = ";",
  col_names = TRUE,
  name_repair = "universal"
)

gravy_exosome_not_nucleus_HepG2 <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\gravy\\GRAVY_Exosome_not_Nucleus_HepG2.txt",
  delim = ";",
  col_names = TRUE,
  name_repair = "universal"
)

#generate data frame
gravy_exosome_HepG2_adj <- gravy_exosome_HepG2 |> 
  mutate(ProteinID = str_extract(header, "(?<=^....).{6}"),
         Exp = "Exosome")

gravy_exosome_nucleus_HepG2_adj <- gravy_exosome_nucleus_HepG2 |> 
  mutate(ProteinID = str_extract(header, "(?<=^....).{6}"),
         Exp = "Exosome, Nucleus")

gravy_nucleus_HepG2_adj <- gravy_nucleus_HepG2 |> 
  mutate(ProteinID = str_extract(header, "(?<=^....).{6}"),
         Exp = "Nucleus")

gravy_nucleus_not_exosome_HepG2_adj <- gravy_nucleus_not_exosome_HepG2 |> 
  mutate(ProteinID = str_extract(header, "(?<=^....).{6}"),
         Exp = "Nucleus, not Exosome")

gravy_exosome_not_nucleus_HepG2_adj <- gravy_exosome_not_nucleus_HepG2 |> 
  mutate(ProteinID = str_extract(header, "(?<=^....).{6}"),
         Exp = "Exosome, not Nucleus")

#combine
gravy_HepG2_combined <- bind_rows(
  gravy_exosome_HepG2_adj,
  gravy_exosome_nucleus_HepG2_adj,
  gravy_nucleus_HepG2_adj,
  gravy_nucleus_not_exosome_HepG2_adj,
  gravy_exosome_not_nucleus_HepG2_adj
)

#wilcox test
gravy_HepG2_wilcox_test <- gravy_HepG2_combined |>
  wilcox_test(GRAVY ~ Exp, p.adjust.method = "BH") |> 
  add_significance("p") |> 
  filter(p.signif != "ns") |> 
  slice(c(1, 4))

#boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_gravy_HepG2 <- gravy_HepG2_combined |> 
  ggplot() +
  geom_boxplot(aes(x = factor(Exp, levels = c("Exosome", "Nucleus", "Exosome, Nucleus", "Nucleus, not Exosome", "Exosome, not Nucleus")), 
                   y = GRAVY, color = factor(Exp, levels = c("Exosome", "Nucleus", "Exosome, Nucleus", "Nucleus, not Exosome", "Exosome, not Nucleus"))), 
               size = 1.5,
               outlier.shape = NA) +
  labs(x = "", y = "GRAVY Score") +
  scale_color_manual(values = c(
    "Exosome" = Color_1,
    "Nucleus" = Color_2,
    "Exosome, Nucleus" = Color_3,
    "Nucleus, not Exosome" = Color_4,
    "Exosome, not Nucleus" = Color_7
  )) +
  coord_cartesian(ylim = c(-1.3, 0.6)) +
  stat_pvalue_manual(data = gravy_HepG2_wilcox_test, hide.ns = "p", label = "p.signif", 
                     tip.length = 0, size = 50,
                     y.position = c(0.5, 0.1)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "boxplot_gravy_HepG2.png"), 
       plot = boxplot_gravy_HepG2, height = 4, width = 2.5, units = "in", dpi = 600)

#mass
#import data
mass_exosome_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\uniprot_mass_exosome_HepG2.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)
mass_exosome_nucleus_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\uniprot_mass_exosome_nucleus_HepG2.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)
mass_nucleus_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\uniprot_mass_nucleus_HepG2.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)
mass_nucleus_not_exosome_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\uniprot_mass_nucleus_not_exosome_HepG2.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)
mass_exosome_not_nucleus_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\uniprot_mass_exosome_not_nucleus_HepG2.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)

#generate data frame
mass_exosome_combined <- bind_rows(
  mass_exosome_HepG2 |> select(From, Mass) |> mutate(Exp = "Exosome", Mass_100000 = Mass/100000),
  mass_exosome_nucleus_HepG2 |> select(From, Mass) |> mutate(Exp = "Exosome, Nucleus", Mass_100000 = Mass/100000),
  mass_nucleus_HepG2 |> select(From, Mass) |> mutate(Exp = "Nucleus", Mass_100000 = Mass/100000),
  mass_exosome_not_nucleus_HepG2 |> select(From, Mass) |> mutate(Exp = "Exosome, not Nucleus", Mass_100000 = Mass/100000),
  mass_nucleus_not_exosome_HepG2 |> select(From, Mass) |> mutate(Exp = "Nucleus, not Exosome", Mass_100000 = Mass/100000)
)

#wilcox test
mass_exosome_wilcox_test <- mass_exosome_combined |> 
  wilcox_test(Mass ~ Exp, p.adjust.method = "BH") |> 
  add_significance("p")

#disorder percentage
#import
RAPID_exosome_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\RAPID\\exosome.csv",
  col_names = TRUE,
  name_repair = "universal"
)
RAPID_exosome_nucleus_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\RAPID\\exosome_nucleus.csv",
  col_names = TRUE,
  name_repair = "universal"
)
RAPID_nucleus_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\RAPID\\nucleus.csv",
  col_names = TRUE,
  name_repair = "universal"
)
RAPID_exosome_not_nucleus_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\RAPID\\exosome_not_nucleus.csv",
  col_names = TRUE,
  name_repair = "universal"
)
RAPID_nucleus_not_exosome_HepG2 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\exosome\\RAPID\\nucleus_not_exosome.csv",
  col_names = TRUE,
  name_repair = "universal"
)

#generate data frame
RAPID_exosome_HepG2_adj <- RAPID_exosome_HepG2 |> 
  mutate(ProteinID = str_extract(Prot..ID, "(?<=^..).{6}"),
         Disorder_Percentage = Disorder.Content..,
         .keep = "none")

RAPID_exosome_nucleus_HepG2_adj <- RAPID_exosome_nucleus_HepG2 |> 
  mutate(ProteinID = str_extract(Prot..ID, "(?<=^..).{6}"),
         Disorder_Percentage = Disorder.Content..,
         .keep = "none")

RAPID_nucleus_HepG2_adj <- RAPID_nucleus_HepG2 |> 
  mutate(ProteinID = str_extract(Prot..ID, "(?<=^..).{6}"),
         Disorder_Percentage = Disorder.Content..,
         .keep = "none")

RAPID_exosome_not_nucleus_HepG2_adj <- RAPID_exosome_not_nucleus_HepG2 |> 
  mutate(ProteinID = str_extract(Prot..ID, "(?<=^..).{6}"),
         Disorder_Percentage = Disorder.Content..,
         .keep = "none")

RAPID_nucleus_not_exosome_HepG2_adj <- RAPID_nucleus_not_exosome_HepG2 |> 
  mutate(ProteinID = str_extract(Prot..ID, "(?<=^..).{6}"),
         Disorder_Percentage = Disorder.Content..,
         .keep = "none")

#combine
RAPID_HepG2_adj_combined <- bind_rows(
  RAPID_exosome_HepG2_adj |> mutate(Exp = "Exosome"),
  RAPID_exosome_nucleus_HepG2_adj |> mutate(Exp = "Exosome, Nucleus"),
  RAPID_nucleus_HepG2_adj |> mutate(Exp = "Nucleus"),
  RAPID_exosome_not_nucleus_HepG2_adj |> mutate(Exp = "Exosome, not Nucleus"),
  RAPID_nucleus_not_exosome_HepG2_adj |> mutate(Exp = "Nucleus, not Exosome")
)

#wilcox test
RAPID_HepG2_adj_wilcox_test <- RAPID_HepG2_adj_combined |> 
  wilcox_test(Disorder_Percentage ~ Exp, p.adjust.method = "BH") |> 
  add_significance("p") |> 
  filter(p.signif != "ns") |> 
  slice(c(1, 3, 5))

#boxplot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_RPAID_HepG2 <- RAPID_HepG2_adj_combined |> 
  ggplot() +
  geom_boxplot(aes(x = factor(Exp, levels = c("Exosome", "Nucleus", "Exosome, Nucleus", "Nucleus, not Exosome", "Exosome, not Nucleus")), 
                   y = Disorder_Percentage, color = factor(Exp, levels = c("Exosome", "Nucleus", "Exosome, Nucleus", "Nucleus, not Exosome", "Exosome, not Nucleus"))), 
               size = 1.5,
               outlier.shape = NA) +
  labs(x = "", y = "Disorder Percentage") +
  scale_color_manual(values = c(
    "Exosome" = Color_1,
    "Nucleus" = Color_2,
    "Exosome, Nucleus" = Color_3,
    "Nucleus, not Exosome" = Color_4,
    "Exosome, not Nucleus" = Color_7
  )) +
  coord_cartesian(ylim = c(0, 80)) +
  stat_pvalue_manual(data = RAPID_HepG2_adj_wilcox_test, hide.ns = "p", label = "p.signif", 
                     tip.length = 0, size = 50,
                     y.position = c(65, 75, 65)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "boxplot_RPAID_HepG2.png"), 
       plot = boxplot_RPAID_HepG2, height = 4, width = 2.5, units = "in", dpi = 600)
