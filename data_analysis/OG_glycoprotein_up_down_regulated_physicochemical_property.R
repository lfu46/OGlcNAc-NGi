#import packages
packages_names <- c("tidyverse", "rstatix", "ggpubr", "showtext")
lapply(packages_names, require, character.only = TRUE)

#OG glycoprotein middle
OG_glycoprotein_median <- OG_glycoprotein_total |> 
  filter(!(UniprotID %in% OG_glycoprotein_generally_up$UniprotID)) |> 
  filter(!(UniprotID %in% OG_glycoprotein_generally_down$UniprotID))

write_xlsx(OG_glycoprotein_median, path = paste0(file_path, "OG_glycoprotein_median.xlsx"))

#disorder percentage
#import data
RAPID_OG_glycoprotein_up_regulated_total <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RAPID\\RAPID_OG_glycoprotein_generally_up.csv",
  col_names = TRUE,
  name_repair = "universal"
)
RAPID_OG_glycoprotein_down_regulated_total <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RAPID\\RAPID_OG_glycoprotein_generally_down.csv",
  col_names = TRUE,
  name_repair = "universal"
)
RAPID_OG_glycoprotein_middle <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\RAPID\\RAPID_OG_glycoprotein_median.csv",
  col_names = TRUE,
  name_repair = "universal"
)

#generate data frame
RAPID_OG_glycoprotein_up_regulated_total <- RAPID_OG_glycoprotein_up_regulated_total |> 
  mutate(ProteinID = str_extract(Prot..ID, "(?<=^..).{6}"),
         Disorder_Percentage = Disorder.Content..,
         .keep = "none")

RAPID_OG_glycoprotein_down_regulated_total <- RAPID_OG_glycoprotein_down_regulated_total |> 
  mutate(ProteinID = str_extract(Prot..ID, "(?<=^..).{6}"),
         Disorder_Percentage = Disorder.Content..,
         .keep = "none")

RAPID_OG_glycoprotein_middle <- RAPID_OG_glycoprotein_middle |> 
  mutate(ProteinID = str_extract(Prot..ID, "(?<=^..).{6}"),
         Disorder_Percentage = Disorder.Content..,
         .keep = "none")

OG_glycoprotein_RAPID_combined <- bind_rows(
  RAPID_OG_glycoprotein_up_regulated_total |> mutate(Exp = "up"),
  RAPID_OG_glycoprotein_down_regulated_total |> mutate(Exp = "down"),
  RAPID_OG_glycoprotein_middle |> mutate(Exp = "median")
)

#wilcox test
OG_glycoprotein_RAPID_wilcox_test <- OG_glycoprotein_RAPID_combined |> 
  wilcox_test(Disorder_Percentage ~ Exp, p.adjust.method = "BH") |> 
  add_significance("p")

#boxplot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_OG_glycoprotein_RAPID <- OG_glycoprotein_RAPID_combined |> 
  ggplot() +
  geom_boxplot(aes(x = factor(Exp, levels = c("up", "median", "down")), y = Disorder_Percentage, color = factor(Exp, levels = c("up", "middle", "down"))), 
               size = 1.5,
               outlier.shape = NA,
               notch = TRUE,
               notchwidth = 0.5) +
  labs(x = "", y = "Disorder Percentage") +
  scale_color_manual(values = c(
    "up" = Color_7,
    "middle" = "gray",
    "down" = Color_6
  )) +
  stat_pvalue_manual(data = OG_glycoprotein_RAPID_wilcox_test, hide.ns = "p", label = "p.signif", 
                     tip.length = 0, size = 50,
                     y.position = c(75, 79)) +
  coord_cartesian(ylim = c(0, 80)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "boxplot_OG_glycoprotein_RAPID.png"), plot = boxplot_OG_glycoprotein_RAPID, height = 3, width = 2.5, units = "in", dpi = 600)

#molecular weight
#import data
mass_OG_glycoprotein_up_regulated_total <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\uniprot_mass\\uniprot_OG_glycoprotein_generally_up.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)
mass_OG_glycoprotein_down_regulated_total <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\uniprot_mass\\uniprot_OG_glycoprotein_generally_down.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)
mass_OG_glycoprotein_middle <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\uniprot_mass\\uniprot_OG_glycoprotein_median.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)

#generate data frame
OG_glycoprotein_mass_combined <- bind_rows(
  mass_OG_glycoprotein_up_regulated_total |> select(From, Mass) |> mutate(Exp = "up", Mass_100000 = Mass/100000),
  mass_OG_glycoprotein_down_regulated_total |> select(From, Mass) |> mutate(Exp = "down", Mass_100000 = Mass/100000),
  mass_OG_glycoprotein_middle |> select(From, Mass) |> mutate(Exp = "median", Mass_100000 = Mass/100000)
)

#wilcox test
OG_glycoprotein_mass_wilcox_test <- OG_glycoprotein_mass_combined |> 
  wilcox_test(Mass ~ Exp, p.adjust.method = "BH") |> 
  add_significance("p")

#boxplot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_OG_glycoprotein_mass <- OG_glycoprotein_mass_combined |> 
  ggplot() +
  geom_boxplot(aes(x = factor(Exp, levels = c("up", "median", "down")), y = Mass_100000, color = factor(Exp, levels = c("up", "middle", "down"))), 
               size = 1.5,
               outlier.shape = NA,
               notch = TRUE,
               notchwidth = 0.5) +
  labs(x = "", y = "MW * 100000") +
  scale_color_manual(values = c(
    "up" = Color_7,
    "middle" = "gray",
    "down" = Color_6
  )) +
  coord_cartesian(ylim = c(0, 3)) +
  stat_pvalue_manual(data = OG_glycoprotein_mass_wilcox_test, hide.ns = "p", label = "p.signif", 
                     tip.length = 0, size = 40,
                     y.position = c(2.8, 3.0)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "boxplot_OG_glycoprotein_mass.png"), 
       plot = boxplot_OG_glycoprotein_mass, height = 4, width = 3, units = "in", dpi = 600)

#gravy score
#import data
gravy_OG_glycoprotein_up_regulated_total <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gravy\\gravy_OG_glycoprotein_general_up.txt",
  delim = ";",
  col_names = TRUE,
  name_repair = "universal"
)
gravy_OG_glycoprotein_down_regulated_total <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gravy\\gravy_OG_glycoprotein_general_down.txt",
  delim = ";",
  col_names = TRUE,
  name_repair = "universal"
)
gravy_OG_glycoprotein_middle <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gravy\\gravy_OG_glycoprotein_median.txt",
  delim = ";",
  col_names = TRUE,
  name_repair = "universal"
)

#generate data frame
gravy_OG_glycoprotein_up_regulated_total <- gravy_OG_glycoprotein_up_regulated_total |> 
  mutate(ProteinID = str_extract(header, "(?<=^....).{6}"),
         Exp = "up")

gravy_OG_glycoprotein_down_regulated_total <- gravy_OG_glycoprotein_down_regulated_total |> 
  mutate(ProteinID = str_extract(header, "(?<=^....).{6}"),
         Exp = "down")

gravy_OG_glycoprotein_middle <- gravy_OG_glycoprotein_middle |> 
  mutate(ProteinID = str_extract(header, "(?<=^....).{6}"),
         Exp = "median")

OG_glycoprotein_gravy_combined <- bind_rows(
  gravy_OG_glycoprotein_up_regulated_total,
  gravy_OG_glycoprotein_down_regulated_total,
  gravy_OG_glycoprotein_middle
)

#wilcox test
OG_glycoprotein_gravy_wilcox_test <- OG_glycoprotein_gravy_combined |>
  wilcox_test(GRAVY ~ Exp, p.adjust.method = "BH") |> 
  add_significance("p")

#boxplot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_OG_glycoprotein_gravy <- OG_glycoprotein_gravy_combined |> 
  ggplot() +
  geom_boxplot(aes(x = factor(Exp, levels = c("up", "median", "down")), y = GRAVY, color = factor(Exp, levels = c("up", "middle", "down"))), 
               size = 1.5,
               outlier.shape = NA,
               notch = TRUE,
               notchwidth = 0.5) +
  labs(x = "", y = "GRAVY Score") +
  scale_color_manual(values = c(
    "up" = Color_7,
    "middle" = "gray",
    "down" = Color_6
  )) +
  coord_cartesian(ylim = c(-1.3, 0.4)) +
  stat_pvalue_manual(data = OG_glycoprotein_gravy_wilcox_test, hide.ns = "p", label = "p.signif", 
                     tip.length = 0, size = 50,
                     y.position = c(0.3)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "boxplot_OG_glycoprotein_gravy.png"), plot = boxplot_OG_glycoprotein_gravy, height = 3, width = 2.5, units = "in", dpi = 600)

#isoelectric point
#import data
#up
iso_OG_glycoprotein_up_file_names <- paste0(file_path, "ipc2\\OG_glycoprotein_up_", seq(1,9), ".csv")
iso_OG_glycoprotein_up_dfs <- list()

for (i in iso_OG_glycoprotein_up_file_names) {
  df <- read_csv(i, col_names = TRUE, name_repair = "universal")
  iso_OG_glycoprotein_up_dfs[[i]] <- df
}

iso_OG_glycoprotein_up_tb <- do.call(bind_rows, iso_OG_glycoprotein_up_dfs) |> 
  mutate(Exp = "up")

#down
iso_OG_glycoprotein_down_file_names <- paste0(file_path, "ipc2\\OG_glycoprotein_down_", seq(1,9), ".csv")
iso_OG_glycoprotein_down_dfs <- list()

for (i in iso_OG_glycoprotein_down_file_names) {
  df <- read_csv(i, col_names = TRUE, name_repair = "universal")
  iso_OG_glycoprotein_down_dfs[[i]] <- df
}

iso_OG_glycoprotein_down_tb <- do.call(bind_rows, iso_OG_glycoprotein_down_dfs) |> 
  mutate(Exp = "down")

#middle
iso_OG_glycoprotein_middle_file_names <- paste0(file_path, "ipc2\\OG_glycoprotein_middle_", seq(1,21), ".csv")
iso_OG_glycoprotein_middle_dfs <- list()

for (i in iso_OG_glycoprotein_middle_file_names) {
  df <- read_csv(i, col_names = TRUE, name_repair = "universal")
  iso_OG_glycoprotein_middle_dfs[[i]] <- df
}

iso_OG_glycoprotein_middle_tb <- do.call(bind_rows, iso_OG_glycoprotein_middle_dfs) |> 
  mutate(Exp = "median")

#combine
OG_glycoprotein_iso_combined <- bind_rows(
  iso_OG_glycoprotein_up_tb,
  iso_OG_glycoprotein_down_tb,
  iso_OG_glycoprotein_middle_tb
)

#wilcox test
OG_glycoprotein_iso_wilcox_test <- OG_glycoprotein_iso_combined |> 
  wilcox_test(IPC2.protein.svr19 ~ Exp, p.adjust.method = "BH") |> 
  add_significance("p")

#boxplot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_OG_glycoprotein_iso <- OG_glycoprotein_iso_combined |> 
  ggplot() +
  geom_boxplot(aes(x = factor(Exp, levels = c("up", "median", "down")), y = IPC2.protein.svr19, color = factor(Exp, levels = c("up", "middle", "down"))), 
               size = 1.5,
               outlier.shape = NA,
               notch = TRUE,
               notchwidth = 0.5) +
  labs(x = "", y = "Isoelectric Point") +
  scale_color_manual(values = c(
    "up" = Color_7,
    "middle" = "gray",
    "down" = Color_6
  )) +
  coord_cartesian(ylim = c(4, 9.5)) +
  stat_pvalue_manual(data = OG_glycoprotein_iso_wilcox_test, hide.ns = "p", label = "p.signif", 
                     tip.length = 0, size = 40,
                     y.position = c(9.5)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "boxplot_OG_glycoprotein_iso.png"), plot = boxplot_OG_glycoprotein_iso, height = 4, width = 3, units = "in", dpi = 600)

