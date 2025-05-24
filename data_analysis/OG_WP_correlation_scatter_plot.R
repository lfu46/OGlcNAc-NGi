#import packages
packages_names <- c("tidyverse", "showtext", "ggpmisc")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
#HepG2
OG_WP_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  select(UniprotID, logFC_OG = logFC) |> 
  left_join(WP_protein_Top_tb_HepG2, by = "UniprotID") |> 
  select(UniprotID, logFC_OG, logFC_WP = logFC) |> 
  filter(!is.na(logFC_WP))

#HEK293T
OG_WP_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  select(UniprotID, logFC_OG = logFC) |> 
  left_join(WP_protein_Top_tb_HEK293T, by = "UniprotID") |> 
  select(UniprotID, logFC_OG, logFC_WP = logFC) |> 
  filter(!is.na(logFC_WP))

#Jurkat
OG_WP_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  select(UniprotID, logFC_OG = logFC) |> 
  left_join(WP_protein_Top_tb_Jurkat, by = "UniprotID") |> 
  select(UniprotID, logFC_OG, logFC_WP = logFC) |> 
  filter(!is.na(logFC_WP))

#generate label dataframe
label_df <- tribble(
  ~ logFC_WP, ~ logFC_OG, ~ label,
  0.5, 1.8, "I",
  -0.9, 1.8, "II",
  -0.5, -1.8, "III",
  0.9, -1.8, "IV"
)

#correlation plot
#HepG2
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

correlation_sactter_point_plot_OG_glycoprotein_WP_HepG2 <- OG_WP_HepG2 |> 
  ggplot(aes(x = logFC_WP, y = logFC_OG)) +
  geom_point() +
  geom_smooth(method = "lm", color = Color_3) +
  stat_poly_eq(formula = y ~ x,
               use_label("eq"),
               parse = TRUE,
               size = 40,
               label.x = 0.95,
               label.y = 0.15,
               color = "red") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(slope = 2, intercept = 0) +
  geom_text(aes(x = logFC_WP, y = logFC_OG, label = label, color = label), 
            size = 80, data = label_df) +
  scale_color_manual(values = c(
    "I" = Color_7,
    "II" = "#35b779",
    "III" = "#31688e",
    "IV" = Color_6
  )) +
  labs(x = expression(log[2]*"(Tuni/Ctrl)"[WP]), y = expression(log[2]*"(Tuni/Ctrl)"[OG])) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-2, 2)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.2),
    axis.title = element_text(size = 100, color = "black"),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black", lineheight = 0.1),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "correlation_sactter_point_plot_OG_glycoprotein_WP_HepG2.png"), plot = correlation_sactter_point_plot_OG_glycoprotein_WP_HepG2, height = 4, width = 4, units = c("in"), dpi = 600)

#HEK293T
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

correlation_sactter_point_plot_OG_glycoprotein_WP_HEK293T <- OG_WP_HEK293T |> 
  ggplot(aes(x = logFC_WP, y = logFC_OG)) +
  geom_point() +
  geom_smooth(method = "lm", color = Color_2) +
  stat_poly_eq(formula = y ~ x,
               use_label("eq"),
               parse = TRUE,
               size = 40,
               label.x = 0.95,
               label.y = 0.15,
               color = "red") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(slope = 2, intercept = 0) +
  geom_text(aes(x = logFC_WP, y = logFC_OG, label = label, color = label), size = 80, 
            data = label_df) +
  scale_color_manual(values = c(
    "I" = Color_7,
    "II" = "#35b779",
    "III" = "#31688e",
    "IV" = Color_6
  )) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-2, 2)) +
  labs(x = expression(log[2]*"(Tuni/Ctrl)"[WP]), y = expression(log[2]*"(Tuni/Ctrl)"[OG])) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.2),
    axis.title = element_text(size = 100, color = "black"),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black", lineheight = 0.1),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "correlation_sactter_point_plot_OG_glycoprotein_WP_HEK293T.png"), plot = correlation_sactter_point_plot_OG_glycoprotein_WP_HEK293T, height = 4, width = 4, units = c("in"), dpi = 600)

#Jurkat
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

correlation_sactter_point_plot_OG_glycoprotein_WP_Jurkat <- OG_WP_Jurkat |> 
  ggplot(aes(x = logFC_WP, y = logFC_OG)) +
  geom_point() +
  geom_smooth(method = "lm", color = Color_4) +
  stat_poly_eq(formula = y ~ x,
               use_label("eq"),
               parse = TRUE,
               size = 40,
               label.x = 0.95,
               label.y = 0.15,
               color = "red") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(slope = 2, intercept = 0) +
  geom_text(aes(x = logFC_WP, y = logFC_OG, label = label, color = label), 
            size = 80, data = label_df) +
  scale_color_manual(values = c(
    "I" = Color_7,
    "II" = "#35b779",
    "III" = "#31688e",
    "IV" = Color_6
  )) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-2, 2)) +
  labs(x = expression(log[2]*"(Tuni/Ctrl)"[WP]), y = expression(log[2]*"(Tuni/Ctrl)"[OG])) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.2),
    axis.title = element_text(size = 100, color = "black"),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black", lineheight = 0.1),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "correlation_sactter_point_plot_OG_glycoprotein_WP_Jurkat.png"), plot = correlation_sactter_point_plot_OG_glycoprotein_WP_Jurkat, height = 4, width = 4, units = c("in"), dpi = 600)

#generate data frame for dot plot
#HepG2
OG_glycoprotein_WP_I_HepG2 <- OG_glycoprotein_WP_HepG2 |> filter(logFC_OG > 0, logFC_WP > 0, logFC_OG > 2*logFC_WP)
OG_glycoprotein_WP_II_HepG2 <- OG_glycoprotein_WP_HepG2 |> filter(logFC_OG > 0, logFC_WP < 0)
OG_glycoprotein_WP_III_HepG2 <- OG_glycoprotein_WP_HepG2 |> filter(logFC_OG < 0, logFC_WP < 0, logFC_OG < 2*logFC_WP)
OG_glycoprotein_WP_IV_HepG2 <- OG_glycoprotein_WP_HepG2 |> filter(logFC_OG < 0, logFC_WP > 0)

write_xlsx(OG_glycoprotein_WP_I_HepG2, paste0(file_path, "OG_glycoprotein_WP_I_HepG2.xlsx"))
write_xlsx(OG_glycoprotein_WP_II_HepG2, paste0(file_path, "OG_glycoprotein_WP_II_HepG2.xlsx"))
write_xlsx(OG_glycoprotein_WP_III_HepG2, paste0(file_path, "OG_glycoprotein_WP_III_HepG2.xlsx"))
write_xlsx(OG_glycoprotein_WP_IV_HepG2, paste0(file_path, "OG_glycoprotein_WP_IV_HepG2.xlsx"))

#HEK293T
OG_glycoprotein_WP_I_HEK293T <- OG_glycoprotein_WP_HEK293T |> filter(logFC_OG > 0, logFC_WP > 0, logFC_OG > 2*logFC_WP)
OG_glycoprotein_WP_II_HEK293T <- OG_glycoprotein_WP_HEK293T |> filter(logFC_OG > 0, logFC_WP < 0)
OG_glycoprotein_WP_III_HEK293T <- OG_glycoprotein_WP_HEK293T |> filter(logFC_OG < 0, logFC_WP < 0, logFC_OG < 2*logFC_WP)
OG_glycoprotein_WP_IV_HEK293T <- OG_glycoprotein_WP_HEK293T |> filter(logFC_OG < 0, logFC_WP > 0)

write_xlsx(OG_glycoprotein_WP_I_HEK293T, paste0(file_path, "OG_glycoprotein_WP_I_HEK293T.xlsx"))
write_xlsx(OG_glycoprotein_WP_II_HEK293T, paste0(file_path, "OG_glycoprotein_WP_II_HEK293T.xlsx"))
write_xlsx(OG_glycoprotein_WP_III_HEK293T, paste0(file_path, "OG_glycoprotein_WP_III_HEK293T.xlsx"))
write_xlsx(OG_glycoprotein_WP_IV_HEK293T, paste0(file_path, "OG_glycoprotein_WP_IV_HEK293T.xlsx"))

#Jurkat
OG_glycoprotein_WP_I_Jurkat <- OG_glycoprotein_WP_Jurkat |> filter(logFC_OG > 0, logFC_WP > 0, logFC_OG > 2*logFC_WP)
OG_glycoprotein_WP_II_Jurkat <- OG_glycoprotein_WP_Jurkat |> filter(logFC_OG > 0, logFC_WP < 0)
OG_glycoprotein_WP_III_Jurkat <- OG_glycoprotein_WP_Jurkat |> filter(logFC_OG < 0, logFC_WP < 0, logFC_OG < 2*logFC_WP)
OG_glycoprotein_WP_IV_Jurkat <- OG_glycoprotein_WP_Jurkat |> filter(logFC_OG < 0, logFC_WP > 0)

write_xlsx(OG_glycoprotein_WP_I_Jurkat, paste0(file_path, "OG_glycoprotein_WP_I_Jurkat.xlsx"))
write_xlsx(OG_glycoprotein_WP_II_Jurkat, paste0(file_path, "OG_glycoprotein_WP_II_Jurkat.xlsx"))
write_xlsx(OG_glycoprotein_WP_III_Jurkat, paste0(file_path, "OG_glycoprotein_WP_III_Jurkat.xlsx"))
write_xlsx(OG_glycoprotein_WP_IV_Jurkat, paste0(file_path, "OG_glycoprotein_WP_IV_Jurkat.xlsx"))

