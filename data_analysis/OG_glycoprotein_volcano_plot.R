#import packages
packages_names <- c("tidyverse", "showtext", "writexl")
lapply(packages_names, require, character.only = TRUE)

#generate up or down regulated OG glycoprotein for each cell line
OG_glycoprotein_up_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC > log2(1.5), adj.P.Val < 0.05) |>   select(UniprotID)

write_xlsx(OG_glycoprotein_up_HepG2, path = paste0(file_path, "OG_glycoprotein_up_HepG2.xlsx"))

OG_glycoprotein_down_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(logFC < -log2(1.5), adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_down_HepG2, path = paste0(file_path, "OG_glycoprotein_down_HepG2.xlsx"))

OG_glycoprotein_median_HepG2 <- OG_glycoprotein_Top_tb_HepG2 |> 
  filter(! UniprotID %in% OG_glycoprotein_up_HepG2$UniprotID) |> 
  filter(! UniprotID %in% OG_glycoprotein_down_HepG2$UniprotID)

write_xlsx(OG_glycoprotein_median_HepG2, path = paste0(file_path, "OG_glycoprotein_median_HepG2.xlsx"))

OG_glycoprotein_up_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC > log2(1.5), adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_up_HEK293T, path = paste0(file_path, "OG_glycoprotein_up_HEK293T.xlsx"))

OG_glycoprotein_down_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(logFC < -log2(1.5), adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_down_HEK293T, path = paste0(file_path, "OG_glycoprotein_down_HEK293T.xlsx"))

OG_glycoprotein_median_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(! UniprotID %in% OG_glycoprotein_up_HEK293T$UniprotID) |> 
  filter(! UniprotID %in% OG_glycoprotein_down_HEK293T$UniprotID)

write_xlsx(OG_glycoprotein_median_HEK293T, path = paste0(file_path, "OG_glycoprotein_median_HEK293T.xlsx"))

OG_glycoprotein_up_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC > log2(1.5), adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_up_Jurkat, path = paste0(file_path, "OG_glycoprotein_up_Jurkat.xlsx"))

OG_glycoprotein_down_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(logFC < -log2(1.5), adj.P.Val < 0.05) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_down_Jurkat, path = paste0(file_path, "OG_glycoprotein_down_Jurkat.xlsx"))

OG_glycoprotein_median_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(! UniprotID %in% OG_glycoprotein_up_Jurkat$UniprotID) |> 
  filter(! UniprotID %in% OG_glycoprotein_down_Jurkat$UniprotID)

write_xlsx(OG_glycoprotein_median_Jurkat, path = paste0(file_path, "OG_glycoprotein_median_Jurkat.xlsx"))

#generate OG glycoprotein up or down regulated common
OG_glycoprotein_up_total <- bind_rows(
  OG_glycoprotein_up_HepG2,
  OG_glycoprotein_up_HEK293T,
  OG_glycoprotein_up_Jurkat
) |> distinct()

write_xlsx(OG_glycoprotein_up_total, path = paste0(file_path, "OG_glycoprotein_up_total.xlsx"))

OG_glycoprotein_down_total <- bind_rows(
  OG_glycoprotein_down_HepG2,
  OG_glycoprotein_down_HEK293T,
  OG_glycoprotein_down_Jurkat
) |> distinct()

write_xlsx(OG_glycoprotein_down_total, path = paste0(file_path, "OG_glycoprotein_down_total.xlsx"))

#HepG2
#volcano plot
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

volcano_plot_OG_glycoprotein_HepG2 <- ggplot() +
  geom_point(data = OG_glycoprotein_median_HepG2,
             aes(x = logFC, y = -log10(adj.P.Val)), size = 1, color = "gray", alpha = 0.5) +
  geom_point(data = OG_glycoprotein_Top_tb_HepG2 |> filter(logFC > log2(1.5), adj.P.Val < 0.05), 
             aes(x = logFC, y = -log10(adj.P.Val), size = -log10(adj.P.Val), color = -log10(adj.P.Val))) +
  geom_point(data = OG_glycoprotein_Top_tb_HepG2 |> filter(logFC < -log2(1.5), adj.P.Val < 0.05), 
             aes(x = logFC, y = -log10(adj.P.Val), size = -log10(adj.P.Val), color = -log10(adj.P.Val))) +
  geom_hline(aes(yintercept = 1.3), linetype = "dashed", color = "black") +
  geom_vline(aes(xintercept = c(-log2(1.5), log2(1.5))),  linetype= "dashed", color = "black") +
  geom_text(aes(x = -1.7, y = 1.6, label = "adj.P = 0.05"), color = "black", size = 25) +
  ggtitle("HepG2") +
  labs(x = expression(log[2]*"(Tuni/Ctrl)"), y = expression(-log[10]*"("*paste("adjusted ", italic(P), " Value")*")")) +
  coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(0, 6.5)) +
  scale_size(range = c(1, 2)) +
  scale_color_gradient(low = alpha(Color_3, 0.5), high = Color_3) +
  theme_bw() +
  theme(
    title = element_text(size = 100, color = "black"),
    panel.grid.major = element_line(color = "gray", linewidth = 0.1),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.05),
    axis.text = element_text(size = 100, color = "black"),
    axis.title.x = element_text(size = 100, color = "black"),
    axis.title.y = element_text(size = 100, color = "black"),
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.position = "none"
  )

ggsave(paste0(file_path, "volcano_plot_OG_glycoprotein_HepG2.png"), plot = volcano_plot_OG_glycoprotein_HepG2, height = 3.5, width = 3, units = c("in"), dpi = 600) 

#volcano plot HEK293T
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

volcano_plot_OG_glycoprotein_HEK293T <- ggplot() +
  geom_point(data = OG_glycoprotein_median_HEK293T,
             aes(x = logFC, y = -log10(adj.P.Val)), size = 1, color = "gray", alpha = 0.5) +
  geom_point(data = OG_glycoprotein_Top_tb_HEK293T |> filter(logFC > log2(1.5), adj.P.Val < 0.05), 
             aes(x = logFC, y = -log10(adj.P.Val), size = -log10(adj.P.Val), color = -log10(adj.P.Val))) +
  geom_point(data = OG_glycoprotein_Top_tb_HEK293T |> filter(logFC < -log2(1.5), adj.P.Val < 0.05), 
             aes(x = logFC, y = -log10(adj.P.Val), size = -log10(adj.P.Val), color = -log10(adj.P.Val))) +
  geom_hline(aes(yintercept = 1.3), linetype = "dashed", color = "black") +
  geom_vline(aes(xintercept = c(-log2(1.5), log2(1.5))),  linetype= "dashed", color = "black") +
  geom_text(aes(x = -1.7, y = 1.6, label = "adj.P = 0.05"), color = "black", size = 25) +
  ggtitle("HEK293T") +
  labs(x = expression(log[2]*"(Tuni/Ctrl)"), y = expression(-log[10]*"("*paste("adjusted ", italic(P), " Value")*")")) +
  coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(0, 6.5)) +
  scale_size(range = c(1, 2)) +
  scale_color_gradient(low = alpha(Color_2, 0.5), high = Color_2) +
  theme_bw() +
  theme(
    title = element_text(size = 100, color = "black"),
    panel.grid.major = element_line(color = "gray", linewidth = 0.1),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.05),
    legend.position = "none",
    axis.text = element_text(size = 100, color = "black"),
    axis.title.x = element_text(size = 100, color = "black"),
    axis.title.y = element_text(size = 100, color = "black"),
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black")
  )

ggsave(paste0(file_path, "volcano_plot_OG_glycoprotein_HEK293T.png"), plot = volcano_plot_OG_glycoprotein_HEK293T, height = 3.5, width = 3, units = c("in"), dpi = 600) 

#volcano plot Jurkat
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

volcano_plot_OG_glycoprotein_Jurkat <- ggplot() +
  geom_point(data = OG_glycoprotein_median_Jurkat,
             aes(x = logFC, y = -log10(adj.P.Val)), size = 1, color = "gray", alpha = 0.5) +
  geom_point(data = OG_glycoprotein_Top_tb_Jurkat |> filter(logFC > log2(1.5), adj.P.Val < 0.05), 
             aes(x = logFC, y = -log10(adj.P.Val), size = -log10(adj.P.Val), color = -log10(adj.P.Val))) +
  geom_point(data = OG_glycoprotein_Top_tb_Jurkat |> filter(logFC < -log2(1.5), adj.P.Val < 0.05), 
             aes(x = logFC, y = -log10(adj.P.Val), size = -log10(adj.P.Val), color = -log10(adj.P.Val))) +
  geom_hline(aes(yintercept = 1.3), linetype = "dashed", color = "black") +
  geom_vline(aes(xintercept = c(-log2(1.5), log2(1.5))),  linetype= "dashed", color = "black") +
  geom_text(aes(x = -1.7, y = 1.6, label = "adj.P = 0.05"), color = "black", size = 25) +
  ggtitle("Jurkat") +
  labs(x = expression(log[2]*"(Tuni/Ctrl)"), y = expression(-log[10]*"("*paste("adjusted ", italic(P), " Value")*")")) +
  coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(0, 6.5)) +
  scale_size(range = c(1, 2)) +
  scale_color_gradient(low = alpha(Color_4, 0.5), high = Color_4) +
  theme_bw() +
  theme(
    title = element_text(size = 100, color = "black"),
    panel.grid.major = element_line(color = "gray", linewidth = 0.1),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.05),
    legend.position = "none",
    axis.text = element_text(size = 100, color = "black"),
    axis.title.x = element_text(size = 100, color = "black"),
    axis.title.y = element_text(size = 100, color = "black"),
    legend.title = element_text(size = 80, color = "black"),
    legend.text = element_text(size = 80, color = "black")
  )

ggsave(paste0(file_path, "volcano_plot_OG_glycoprotein_Jurkat.png"), plot = volcano_plot_OG_glycoprotein_Jurkat, height = 3.5, width = 3, units = c("in"), dpi = 600) 
