#import packages
packages_names <- c("tidyverse", "readxl", "writexl", 
                    "showtext", "ggpubr", "rstatix", "ComplexHeatmap", "circlize", "eulerr")
lapply(packages_names, require, character.only = TRUE)

#eruler plot for OG glycoprotein up in three cell lines
mat_up_eruler <- c(
  "HEK293T (191)" = 144,
  "HepG2 (155)" = 116,
  "Jurkat (124)" = 95,
  "HEK293T (191)&HepG2 (155)" = 26,
  "HEK293T (191)&Jurkat (124)" = 16,
  "HepG2 (155)&Jurkat (124)" = 8,
  "HepG2 (155)&HEK293T (191)&Jurkat (124)" = 5
)

fit_up_eruler <- euler(mat_up_eruler)

tiff(filename = paste0(file_path, "figureS3_up_remake.tiff"), 
     width = 3, height = 2, units = "in", compression = "lzw",
     res = 1200)

plot(fit_up_eruler, 
     fill = list(fill = c(Color_2, Color_3, Color_4), alpha = 0.8),
     edges = list(col = "white", lwd = 1),
     check.overlap = TRUE,
     legend = list(fontsize = 0)
)

dev.off()

#eruler plot for OG glycoprotein down in three cell lines
mat_down_eruler <- c(
  "HEK293T (218)" = 177,
  "HepG2 (107)" = 74,
  "Jurkat (114)" = 84,
  "HEK293T (218)&HepG2 (107)" = 21,
  "HEK293T (218)&Jurkat (114)" = 18,
  "HepG2 (107)&Jurkat (114)" = 10,
  "HepG2 (107)&HEK293T (218)&Jurkat (114)" = 2
)

fit_down_eruler <- euler(mat_down_eruler)

tiff(filename = paste0(file_path, "figureS3_down_remake.tiff"), 
     width = 3, height = 2, units = "in", compression = "lzw",
     res = 1200)

plot(fit_down_eruler, 
     fill = list(fill = c(Color_2, Color_3, Color_4), alpha = 0.8),
     edges = list(col = "white", lwd = 1),
     legend = list(fontsize = 0),
     check.overlap = TRUE
)

dev.off()
