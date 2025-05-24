#import packages
packages_names <- c("tidyverse", "writexl", "readxl", "ComplexHeatmap", "circlize")
lapply(packages_names, require, character.only = TRUE)

#import data
up_overlap <- c("Q9HCE1", "Q96P11", "Q15437", "Q12888", "P61019")
down_overlap <- c("Q8WWI1", "Q9P219")

uniprot_OG_glycoprotein_up_down_overlap <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\data_source\\uniprot_OG_glycoprotein_up_down_overlap.xlsx",
  col_names = TRUE,
  .name_repair = "universal"
)

#generate data frame
OG_glycoprotein_up_overlap <- bind_rows(
  OG_glycoprotein_Top_tb_HEK293T |> 
    filter(UniprotID %in% up_overlap) |> 
    mutate(cell = 'HEK293T'),
  OG_glycoprotein_Top_tb_HepG2 |> 
    filter(UniprotID %in% up_overlap) |> 
    mutate(cell = 'HepG2'),
  OG_glycoprotein_Top_tb_Jurkat |> 
    filter(UniprotID %in% up_overlap) |> 
    mutate(cell = 'Jurkat')
)

OG_glycoprotein_down_overlap <- bind_rows(
  OG_glycoprotein_Top_tb_HEK293T |> 
    filter(UniprotID %in% down_overlap) |> 
    mutate(cell = 'HEK293T'),
  OG_glycoprotein_Top_tb_HepG2 |> 
    filter(UniprotID %in% down_overlap) |> 
    mutate(cell = 'HepG2'),
  OG_glycoprotein_Top_tb_Jurkat |> 
    filter(UniprotID %in% down_overlap) |> 
    mutate(cell = 'Jurkat')
)

OG_glycoprotein_up_down_overlap <- bind_rows(
  OG_glycoprotein_up_overlap,
  OG_glycoprotein_down_overlap
) |> left_join(uniprot_OG_glycoprotein_up_down_overlap, by = join_by(UniprotID == From)) |> 
  separate(col = Gene.Names, into = c('Gene_Name'))

OG_glycoprotein_up_down_overlap_adj <- OG_glycoprotein_up_down_overlap |> 
  select(Gene_Name, cell, logFC) |> 
  pivot_wider(names_from = Gene_Name, values_from = logFC)

#heatmap
mat <- data.matrix(OG_glycoprotein_up_down_overlap_adj |> select(!cell))
rownames(mat) <- OG_glycoprotein_up_down_overlap_adj$cell
mat_adj <- t(mat)
col_mat <- colorRamp2(breaks = c(-1, 0 ,1), colors = c("blue", "white", "red"))

h1 <- Heatmap(matrix = mat_adj, name = 'log2(Tuni/Ctrl)', 
              row_names_side = "left", 
              col = col_mat)
