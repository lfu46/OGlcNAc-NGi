#import packages
packages_names <- c("tidyverse", "writexl", "readxl", "showtext")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
raw_data <- tribble(
  ~ Channel, ~ Intensity, ~ Exp,
  "Tuni 1", 313948.1417, "Exp 1",
  "Tuni 2", 315790.9971, "Exp 2",
  "Tuni 3", 302237.3513, "Exp 3",
  "Ctrl 4", 67272.67622, "Exp 4",
  "Ctrl 5", 67193.37058, "Exp 5",
  "Ctrl 6", 73063.01217, "Exp 6"
)

report_ion_plot <- raw_data |> 
  ggplot(aes(x = factor(Channel, levels = c("Tuni 1", "Tuni 2", "Tuni 3", "Ctrl 4", "Ctrl 5", "Ctrl 6")), y = Intensity, fill = Exp)) +
  geom_bar(stat = "identity", width = 0.1) +
  labs(x = "", y = "") +
  scale_fill_manual(values = c(
    "Exp 1" = "#317EC2",
    "Exp 2" = "#A83E3F",
    "Exp 3" = "#E69965",
    "Exp 4" = "#584482",
    "Exp 5" = "#6FAA54",
    "Exp 6" = "#8785BA"
  )) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

ggsave(filename =  "E:\\Glycosylation_Crosstalk\\Figures\\workflow_example.png", plot = report_ion_plot, height = 4, width = 3, units = c("in"), dpi = 600)

#generate data frame
raw_data_2 <- tribble(
  ~ Channel, ~ Intensity, ~ Exp,
  "Tuni 1", 620.345503671639, "Exp 1",
  "Tuni 2", 625.860494690048, "Exp 2",
  "Tuni 3", 541.154374298508, "Exp 3",
  "Ctrl 4", 1758.16626221336, "Exp 4",
  "Ctrl 5", 1529.11168991783, "Exp 5",
  "Ctrl 6", 1678.1482177534, "Exp 6"
)

font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

report_ion_plot_2 <- raw_data_2 |> 
  ggplot(aes(x = factor(Channel, levels = c("Tuni 1", "Tuni 2", "Tuni 3", "Ctrl 4", "Ctrl 5", "Ctrl 6")), y = Intensity, fill = Exp)) +
  geom_bar(stat = "identity", width = 0.1) +
  labs(x = "", y = "") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c(
    "Exp 1" = "#317EC2",
    "Exp 2" = "#A83E3F",
    "Exp 3" = "#E69965",
    "Exp 4" = "#584482",
    "Exp 5" = "#6FAA54",
    "Exp 6" = "#8785BA"
  )) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

ggsave(filename =  "E:\\Glycosylation_Crosstalk\\Figures\\workflow_example_2.eps", 
       plot = report_ion_plot_2,
       device = "eps",
       height = 4, width = 2, units = "in", 
       dpi = 600)

#mass spectra 
#E_LF_NOCT_HEK_OG_4_08102024_16496
#import data
spectra_16496 <- read_csv(
  "E:\\Glycosylation_Crosstalk\\Figures\\figure1\\spectrum_16496.csv",
  skip = 7,
  col_names = TRUE,
  name_repair = "universal"
) |> mutate(
  Relative_Intensity = Intensity/max(Intensity) * 100
)

#mass spectra 16496
font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

mass_spectra_16496 <- spectra_16496 |> 
  ggplot() +
  geom_col(aes(x = Mass, y = Relative_Intensity), color = "black", width = 0.00001) +
  labs(x = "", y = "Relative Abundance") +
  scale_x_continuous(breaks = c(0, 400, 800, 1200, 1600), expand = c(0.01, 0.1), limits = c(0, 1700)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.title = element_text(color = "black", size = 10),
    axis.text = element_text(color = "black", size = 8)
  )

ggsave(
  filename = "E:\\Glycosylation_Crosstalk\\Figures\\figure1\\mass_spectra_16496.eps",
  plot = mass_spectra_16496, 
  device = "eps",
  width = 3, height = 2.5, units = "in"
)
