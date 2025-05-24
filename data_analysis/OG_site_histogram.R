#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "showtext")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
OG_site_glycoprotein_HepG2 <- OG_glycopeptide_Top_tb_HepG2 |> 
  filter(!is.na(combined_site)) |> 
  select(Index, UniprotID)

OG_site_glycoprotein_HEK293T <- OG_glycopeptide_Top_tb_HEK293T |> 
  filter(!is.na(combined_site)) |> 
  select(Index, UniprotID)

OG_site_glycoprotein_Jurkat <- OG_glycopeptide_Top_tb_Jurkat |> 
  filter(!is.na(combined_site)) |> 
  select(Index, UniprotID)

#combine
OG_site_glycoprotein_total <- bind_rows(
  OG_site_glycoprotein_HepG2,
  OG_site_glycoprotein_HEK293T,
  OG_site_glycoprotein_Jurkat
) |> distinct()

#generate data frame
OG_site_glycoprotein_count <- OG_site_glycoprotein_total |> 
  group_by(UniprotID) |> 
  summarize(
    site_count = n()
  ) |> 
  arrange(desc(site_count)) |> 
  mutate(site_count = ifelse(site_count > 9, 10, site_count))

#site level
#site number histogram
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

OG_site_count_histogram <- OG_site_glycoprotein_count |> 
  ggplot(aes(x = site_count)) +
  geom_histogram(bins = 10) +
  labs(x = "Localized # of O-GlcNAcylation site\nin each glycoprotein", y = "") +
  coord_cartesian(xlim = c(1, 10)) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10), labels = c('2', '4', '6', '8', '10 or more')) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title.x = element_text(size = 80, lineheight = 0.1),
    axis.text.x = element_text(size = 80, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "OG_site_count_histogram.png"), plot = OG_site_count_histogram,
       height = 4, width = 3.5, dpi = 600)
