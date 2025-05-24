#import packages
packages_names <- c("tidyverse", "writexl", "clusterProfiler", "org.Hs.eg.db", "showtext", "eulerr")
lapply(packages_names, require, character.only = TRUE)

#CORUM database
Uniprot_database <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\Pfam\\uniprotkb_organism_id_9606_AND_reviewed_2024_11_25.xlsx",
  col_names = TRUE
)

Pfam_database <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\Pfam\\proteomes\\Pfam_homo_sapiens.xlsx",
  skip = 2,
  col_names = TRUE,
  .name_repair = "universal"
) |> dplyr::select(Protein_domain = .hmm.name., UniprotID = .seq.id.) |> 
  filter(UniprotID %in% Uniprot_database$Entry)

write_xlsx(Pfam_database, path = paste0(file_path, "Pfam_database.xlsx"))

#GSEA protein domain
#HEK293T
OG_glycoprotein_List_HEK293T <- OG_glycoprotein_Top_tb_HEK293T$logFC
names(OG_glycoprotein_List_HEK293T) <- OG_glycoprotein_Top_tb_HEK293T$UniprotID
OG_glycoprotein_List_HEK293T <- sort(OG_glycoprotein_List_HEK293T, decreasing = TRUE)

protein_domain_GSEA_HEK293T <- GSEA(OG_glycoprotein_List_HEK293T, TERM2GENE = Pfam_database,
                                    pvalueCutoff = 0.5)

#HepG2
OG_glycoprotein_List_HepG2 <- OG_glycoprotein_Top_tb_HepG2$logFC
names(OG_glycoprotein_List_HepG2) <- OG_glycoprotein_Top_tb_HepG2$UniprotID
OG_glycoprotein_List_HepG2 <- sort(OG_glycoprotein_List_HepG2, decreasing = TRUE)

protein_domain_GSEA_HepG2 <- GSEA(OG_glycoprotein_List_HepG2, TERM2GENE = Pfam_database,
                                  pvalueCutoff = 0.5)

write_xlsx(protein_domain_GSEA_HepG2@result, path = paste0(file_path, "protein_domain_GSEA_HepG2.xlsx"))

#Jurkat
OG_glycoprotein_List_Jurkat <- OG_glycoprotein_Top_tb_Jurkat$logFC
names(OG_glycoprotein_List_Jurkat) <- OG_glycoprotein_Top_tb_Jurkat$UniprotID
OG_glycoprotein_List_Jurkat <- sort(OG_glycoprotein_List_Jurkat, decreasing = TRUE)

protein_domain_GSEA_Jurkat <- GSEA(OG_glycoprotein_List_Jurkat, TERM2GENE = Pfam_database,
                                   pvalueCutoff = 0.5)

