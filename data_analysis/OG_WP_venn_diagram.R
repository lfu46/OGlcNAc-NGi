#import packages
packages_names <- c("VennDiagram")
lapply(packages_names, require, character.only = TRUE)

#generate total glycoprotein
OG_glycoprotein_total <- bind_rows(
  OG_glycoprotein_Top_tb_HepG2 |> distinct(UniprotID),
  OG_glycoprotein_Top_tb_HEK293T |> distinct(UniprotID),
  OG_glycoprotein_Top_tb_Jurkat |> distinct(UniprotID)
) |> distinct()

write_xlsx(OG_glycoprotein_total, paste0(file_path, "OG_glycoprotein_total.xlsx"))

#generate total protein
WP_protein_total <- bind_rows(
  WP_protein_Top_tb_HepG2 |> distinct(UniprotID),
  WP_protein_Top_tb_HEK293T |> distinct(UniprotID),
  WP_protein_Top_tb_Jurkat |> distinct(UniprotID)
) |> distinct()

write_xlsx(WP_protein_total, paste0(file_path, "WP_protein_total.xlsx"))

#venn diagram
venn.diagram(
  x = list(OG_glycoprotein_total |> pull(), WP_protein_total |> pull()),
  category.names = c("O-GlcNAcylated proteins \n (1693)", "Total proteins \n (9728)"),
  filename = "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\OG_WP_venn_diagram.png", 
  fill = c(Color_1, Color_2),
  output = TRUE,
  imagetype = "png",
  height = 3,
  width = 3,
  units = c("in"),
  resolution = 1200,
  lty = 'blank',
  cex = 1.0,
  fontfamily = "sans",
  cat.cex = 1.0,
  cat.default.pos = "outer",
  cat.pos = c(15, 150),
  cat.dist = c(0.01, 0.04),
  cat.fontfamily = "sans",
  margin = 0.05
)

