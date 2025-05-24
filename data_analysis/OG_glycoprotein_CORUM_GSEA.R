#import packages
packages_names <- c("tidyverse", "writexl", "clusterProfiler", "org.Hs.eg.db", 
                    "showtext", "eulerr", "enrichplot")
lapply(packages_names, require, character.only = TRUE)

#CORUM database
CORUM_database <- read_delim(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\CORUM_protein_complex\\corum_allComplexes.txt",
  col_names = TRUE,
  name_repair = "universal"
) |> filter(organism == "Human") |> 
  dplyr::select(complex_name, UniprotID = subunits_uniprot_id) |> 
  separate_rows(UniprotID, sep = ";")

write_xlsx(CORUM_database, path = paste0(file_path, "CORUM_database.xlsx"))

#GSEA protein complex
#HEK293T
OG_glycoprotein_List_HEK293T <- OG_glycoprotein_Top_tb_HEK293T$logFC
names(OG_glycoprotein_List_HEK293T) <- OG_glycoprotein_Top_tb_HEK293T$UniprotID
OG_glycoprotein_List_HEK293T <- sort(OG_glycoprotein_List_HEK293T, decreasing = TRUE)

protein_complex_GSEA_HEK293T <- GSEA(OG_glycoprotein_List_HEK293T, TERM2GENE = CORUM_database)

#HepG2
OG_glycoprotein_List_HepG2 <- OG_glycoprotein_Top_tb_HepG2$logFC
names(OG_glycoprotein_List_HepG2) <- OG_glycoprotein_Top_tb_HepG2$UniprotID
OG_glycoprotein_List_HepG2 <- sort(OG_glycoprotein_List_HepG2, decreasing = TRUE)

protein_complex_GSEA_HepG2 <- GSEA(OG_glycoprotein_List_HepG2, TERM2GENE = CORUM_database)

write_xlsx(protein_complex_GSEA_HepG2@result, path = paste0(file_path, "protein_complex_GSEA_HepG2.xlsx"))

#Jurkat
OG_glycoprotein_List_Jurkat <- OG_glycoprotein_Top_tb_Jurkat$logFC
names(OG_glycoprotein_List_Jurkat) <- OG_glycoprotein_Top_tb_Jurkat$UniprotID
OG_glycoprotein_List_Jurkat <- sort(OG_glycoprotein_List_Jurkat, decreasing = TRUE)

protein_complex_GSEA_Jurkat <- GSEA(OG_glycoprotein_List_Jurkat, TERM2GENE = CORUM_database)

#HepG2 spliceosome complex
spliceosome_E_complex <- CORUM_database |> 
  filter(complex_name == "Spliceosome, E complex") |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = 'UniprotID') |> 
  select(complex_name, UniprotID, logFC) |> 
  filter(!is.na(logFC))

spliceosome_A_complex <- CORUM_database |> 
  filter(complex_name == "Spliceosome, A complex") |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = 'UniprotID') |> 
  select(complex_name, UniprotID, logFC) |> 
  filter(!is.na(logFC))

spliceosome_B_complex <- CORUM_database |> 
  filter(complex_name == "Spliceosome, B complex") |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = 'UniprotID') |> 
  select(complex_name, UniprotID, logFC) |> 
  filter(!is.na(logFC))

spliceosome_C_complex <- CORUM_database |> 
  filter(complex_name == "Spliceosome, C complex") |> 
  left_join(OG_glycoprotein_Top_tb_HepG2, by = 'UniprotID') |> 
  select(complex_name, UniprotID, logFC) |> 
  filter(!is.na(logFC))

#combined
spliceosome_complex_combined_HepG2 <- bind_rows(
  spliceosome_E_complex,
  spliceosome_A_complex
)

#wilcox test
spliceosome_complex_combined_HepG2 |> 
  wilcox_test(logFC ~ complex_name)

#GSEA enrichment plot
GSEA_enrichmentplot_spliceosomeEcomplex_HepG2 <- gseaplot(protein_complex_GSEA_HepG2, geneSetID = 1, 
         color.line = Color_3, by = "runningScore") + 
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black")
  )

ggsave(filename = paste0(file_path, "GSEA_enrichmentplot_spliceosomeEcomplex_HepG2.eps"), 
       device = "eps",
       plot = GSEA_enrichmentplot_spliceosomeEcomplex_HepG2, 
       height = 2, width = 2, 
       units = "in", 
       dpi = 1200)
