#import packages
packages_names <- c("tidyverse", "showtext", "writexl", "readxl", "networkD3", "UpSetR")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
OG_glycoprotein_median_HEK293T <- OG_glycoprotein_Top_tb_HEK293T |> 
  filter(!(UniprotID %in% OG_glycoprotein_up_HEK293T$UniprotID)) |> 
  filter(!(UniprotID %in% OG_glycoprotein_down_HEK293T$UniprotID)) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_median_HEK293T, path = paste0(file_path, "OG_glycoprotein_median_HEK293T.xlsx"))

OG_glycoprotein_median_Jurkat <- OG_glycoprotein_Top_tb_Jurkat |> 
  filter(!(UniprotID %in% OG_glycoprotein_up_Jurkat$UniprotID)) |> 
  filter(!(UniprotID %in% OG_glycoprotein_down_Jurkat$UniprotID)) |> 
  select(UniprotID)

write_xlsx(OG_glycoprotein_median_Jurkat, path = paste0(file_path, "OG_glycoprotein_median_Jurkat.xlsx"))

#upset plot for OG glycoprotein HEK293T and Jurkat
OG_glycoprotein_upset_up_HEK293T <- OG_glycoprotein_up_HEK293T |> mutate(up_HEK293T = 1)
OG_glycoprotein_upset_median_HEK293T <- OG_glycoprotein_median_HEK293T |> mutate(median_HEK293T = 1)
OG_glycoprotein_upset_down_HEK293T <- OG_glycoprotein_down_HEK293T |> mutate(down_HEK293T = 1)

OG_glycoprotein_upset_up_Jurkat <- OG_glycoprotein_up_Jurkat |> mutate(up_Jurkat = 1)
OG_glycoprotein_upset_median_Jurkat <- OG_glycoprotein_median_Jurkat |> mutate(median_Jurkat = 1)
OG_glycoprotein_upset_down_Jurkat <- OG_glycoprotein_down_Jurkat |> mutate(down_Jurkat = 1)

upset_OG_glycoprotein_HEK293T_Jurkat <- OG_glycoprotein_total |> 
  left_join(OG_glycoprotein_upset_up_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_upset_median_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_upset_down_HEK293T, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_upset_up_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_upset_median_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  left_join(OG_glycoprotein_upset_down_Jurkat, by = join_by(UniprotID == UniprotID)) |> 
  mutate(
    up_HEK293T = ifelse(is.na(up_HEK293T), 0, 1),
    median_HEK293T = ifelse(is.na(median_HEK293T), 0, 1),
    down_HEK293T = ifelse(is.na(down_HEK293T), 0, 1),
    up_Jurkat = ifelse(is.na(up_Jurkat), 0, 1),
    median_Jurkat = ifelse(is.na(median_Jurkat), 0, 1),
    down_Jurkat = ifelse(is.na(down_Jurkat), 0, 1)
  )

upset_OG_glycoprotein_HEK293T_Jurkat_dataframe <- as.data.frame(upset_OG_glycoprotein_HEK293T_Jurkat)

upset(upset_OG_glycoprotein_HEK293T_Jurkat_dataframe, 
      sets = c("up_HEK293T", "median_HEK293T", "down_HEK293T", "up_Jurkat", "median_Jurkat", "down_Jurkat"),
      order.by = "freq", point.size = 5,
      text.scale = c(1.5, 1.5, 1.5, 1.5, 2, 2))

#sankey plot
#generate data frame
link_HEK293T_Jurkat <- tribble(
  ~ source, ~ target, ~ value,
  "HEK293T up", "Jurkat median ", 73,
  "HEK293T up", "Jurkat down ", 8,
  "HEK293T median", "Jurkat down ", 51,
  "HEK293T median", "Jurkat up ", 32,
  "HEK293T down", "Jurkat median ", 78,
  "HEK293T down", "Jurkat up ", 6
)

node_HEK293T_Jurkat <- data.frame(
  name = c(as.character(link_HEK293T_Jurkat$source), as.character(link_HEK293T_Jurkat$target)) |> unique()
)

link_HEK293T_Jurkat$sourceID <- match(link_HEK293T_Jurkat$source, node_HEK293T_Jurkat$name) - 1
link_HEK293T_Jurkat$targetID <- match(link_HEK293T_Jurkat$target, node_HEK293T_Jurkat$name) - 1

link_HEK293T_Jurkat_dataframe <- as.data.frame(link_HEK293T_Jurkat)
ColourScal ='d3.scaleOrdinal() .range(["#7FB2D4","#00bfc4"])'

sankeyNetwork(Links = link_HEK293T_Jurkat_dataframe, Nodes = node_HEK293T_Jurkat,
              Source = "sourceID", Target = "targetID",
              Value = "value", NodeID = "name",
              sinksRight = FALSE, colourScale = ColourScal, 
              height = 574, width = 574, units = "in",
              nodeWidth = 40 , fontSize = 20, nodePadding = 10)

