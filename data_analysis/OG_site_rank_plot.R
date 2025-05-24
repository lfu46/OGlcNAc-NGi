#import packages
packages_names <- c("tidyverse", "writexl", "readxl", "showtext")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
OG_glycopeptide_localized_glycoprotein_HepG2 <- OG_glycopeptide_Top_tb_HepG2 |> 
  filter(!is.na(combined_site)) |> 
  arrange(desc(logFC)) |> 
  mutate(rank = rank(logFC),
         group = ifelse((logFC > 0.5 & adj.P.Val < 0.05), "up", 
                        ifelse((logFC < -0.5 & adj.P.Val < 0.05), "down", "median"))) |> 
  select(Index, logFC, rank, adj.P.Val, group)

OG_glycopeptide_localized_glycoprotein_HEK293T <- OG_glycopeptide_Top_tb_HEK293T |> 
  filter(!is.na(combined_site)) |> 
  arrange(desc(logFC)) |> 
  mutate(rank = rank(logFC),
         group = ifelse((logFC > 0.5 & adj.P.Val < 0.05), "up", 
                        ifelse((logFC < -0.5 & adj.P.Val < 0.05), "down", "median"))) |> 
  select(Index, logFC, rank, adj.P.Val, group)

OG_glycopeptide_localized_glycoprotein_Jurkat <- OG_glycopeptide_Top_tb_Jurkat |> 
  filter(!is.na(combined_site)) |> 
  arrange(desc(logFC)) |> 
  mutate(rank = rank(logFC),
         group = ifelse((logFC > 0.5 & adj.P.Val < 0.05), "up", 
                        ifelse((logFC < -0.5 & adj.P.Val < 0.05), "down", "median"))) |> 
  select(Index, logFC, rank, adj.P.Val, group)

#site rank plot
#HepG2
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

site_rank_HepG2 <- OG_glycopeptide_localized_glycoprotein_HepG2 |> 
  ggplot(aes(x = rank, y = logFC)) +
  geom_point(aes(fill = group), shape = 21, color = "transparent") +
  scale_fill_manual(
    values = c(
      "up" = Color_7,
      "median" = "gray",
      "down" = Color_6
    )
  ) +
  ggtitle("Localized O-GlcNAc Site \nin HepG2 (439)") +
  labs(x = "Site rank", y = expression(log[2]*"(Tuni/Ctrl)")) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_bw() +
  theme(
    title = element_text(size = 100, color = "black", lineheight = 0.15),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "site_rank_HepG2.png"), plot = site_rank_HepG2, height = 4, width = 4, units = "in", dpi = 600)

#HEK293T
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

site_rank_HEK293T <- OG_glycopeptide_localized_glycoprotein_HEK293T |> 
  ggplot(aes(x = rank, y = logFC)) +
  geom_point(aes(fill = group), shape = 21, color = "transparent") +
  scale_fill_manual(
    values = c(
      "up" = Color_5,
      "median" = "gray",
      "down" = Color_6
    )
  ) +
  ggtitle("Localized O-GlcNAc Site \nin HEK293T (2146)") +
  labs(x = "Site rank", y = expression(log[2]*"(Tuni/Ctrl)")) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_bw() +
  theme(
    title = element_text(size = 100, color = "black", lineheight = 0.15),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "site_rank_HEK293T.png"), plot = site_rank_HEK293T, height = 4, width = 4, units = "in", dpi = 600)

#Jurkat
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

site_rank_Jurkat <- OG_glycopeptide_localized_glycoprotein_Jurkat |> 
  ggplot(aes(x = rank, y = logFC)) +
  geom_point(aes(fill = group), shape = 21, color = "transparent") +
  scale_fill_manual(
    values = c(
      "up" = Color_5,
      "median" = "gray",
      "down" = Color_6
    )
  ) +
  ggtitle("Localized O-GlcNAc Site \nin Jurkat (1324)") +
  labs(x = "Site rank", y = expression(log[2]*"(Tuni/Ctrl)")) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_bw() +
  theme(
    title = element_text(size = 100, color = "black", lineheight = 0.15),
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 80, color = "black"),
    axis.text.y = element_text(size = 80, color = "black"),
    legend.position = "none"
  )

ggsave(filename = paste0(file_path, "site_rank_Jurkat.png"), plot = site_rank_Jurkat, height = 4, width = 4, units = "in", dpi = 600)

