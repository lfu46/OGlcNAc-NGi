#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "eulerr")
lapply(packages_names, require, character.only = TRUE)

#site level
#generate data frame
#HepG2
OG_site_HepG2 <- OG_glycopeptide_Top_tb_HepG2 |> 
  filter(!is.na(combined_site)) |> 
  pull(Index)

#HEK293T
OG_site_HEK293T <- OG_glycopeptide_Top_tb_HEK293T |> 
  filter(!is.na(combined_site)) |> 
  pull(Index)

#Jurkat
OG_site_Jurkat <- OG_glycopeptide_Top_tb_Jurkat |> 
  filter(!is.na(combined_site)) |> 
  pull(Index)


#euler diagram
mat <- c(
  "HEK293T (2146)" = 1368,
  "HepG2 (439)" = 109,
  "Jurkat (1324)" = 580,
  "HEK293T (2146)&HepG2 (439)" = 70,
  "HEK293T (2146)&Jurkat (1324)" = 484,
  "HepG2 (439)&Jurkat (1324)" = 36,
  "HepG2 (439)&HEK293T (2146)&Jurkat (1324)" = 224
)

fit <- euler(mat)

tiff(filename = paste0(file_path, "figure2_euler_remake.tiff"), 
     width = 2, height = 2, units = "in", compression = "lzw",
     res = 1200)

plot(fit, 
     fill = list(fill = c(Color_2, Color_3, Color_4), alpha = 0.8),
     edges = list(col = "white", lwd = 3),
     labels = list(alpha = c(0)),
     legend = list(fontsize = 0))

dev.off()
