#import packages
packages_names <- c("tidyverse", "readxl", "writexl", "clusterProfiler", "org.Hs.eg.db", "showtext")
lapply(packages_names, require, character.only = TRUE)

#GO analysis
#up
OG_glycoprotein_generally_up_GO <- enrichGO(
  gene = OG_glycoprotein_generally_up$UniprotID,
  OrgDb = org.Hs.eg.db,
  keyType = "UNIPROT",
  ont = "ALL",
  universe = OG_glycoprotein_total$UniprotID,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  minGSSize = 1
)

write_xlsx(OG_glycoprotein_generally_up_GO@result, path = paste0(file_path, "OG_glycoprotein_generally_up_GO.xlsx"))

#down
OG_glycoprotein_generally_down_GO <- enrichGO(
  gene = OG_glycoprotein_generally_down$UniprotID,
  OrgDb = org.Hs.eg.db,
  keyType = "UNIPROT",
  ont = "ALL",
  universe = OG_glycoprotein_total$UniprotID,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  minGSSize = 1
)

write_xlsx(OG_glycoprotein_generally_down_GO@result, path = paste0(file_path, "OG_glycoprotein_generally_down_GO.xlsx"))

#protein domain
Pfam_database <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\Pfam\\Pfam_database.xlsx"
)

OG_glycoprotein_generally_up_down <- bind_rows(
  OG_glycoprotein_generally_up,
  OG_glycoprotein_generally_down
) |> distinct() |> pull()

OG_glycoprotein_protein_domain <- enricher(
  gene = OG_glycoprotein_generally_up_down,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  minGSSize = 3,
  TERM2GENE = Pfam_database
)

write_xlsx(OG_glycoprotein_protein_domain@result, path = paste0(file_path, "OG_glycoprotein_protein_domain.xlsx"))

#protein complex
CORUM_database <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\CORUM_protein_complex\\CORUM_database.xlsx"
)

OG_glycoprotein_protein_complex <- enricher(
  gene = OG_glycoprotein_generally_up_down,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  minGSSize = 3,
  TERM2GENE = CORUM_database
)

write_xlsx(OG_glycoprotein_protein_complex@result, path = paste0(file_path, "OG_glycoprotein_protein_complex.xlsx"))

#import data
OG_glycoprotein_protein_complex <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\OG_glycoprotein_further_analysis\\OG_glycoprotein_protein_complex.xlsx",
  sheet = 'Sheet2'
) |> mutate(Category = 'Protein Complex')

OG_glycoprotein_protein_domain <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\OG_glycoprotein_further_analysis\\OG_glycoprotein_protein_domain.xlsx",
  sheet = 'Sheet2'
) |> mutate(Category = 'Protein Domain')

OG_glycoprotein_complex_domain <- bind_rows(
  OG_glycoprotein_protein_complex,
  OG_glycoprotein_protein_domain
) |> dplyr::select(Description, pvalue, Category) |> mutate(logpvalue = -log10(pvalue))

#barplot
font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

barplot_OG_glycoprotein_complex_domain <- OG_glycoprotein_complex_domain |> 
  ggplot() +
  geom_bar(aes(x = logpvalue, y = fct_reorder(Description, logpvalue), fill = Category), stat = "identity") +
  scale_fill_manual(
    values = c(
      "Protein Complex" = Color_1,
      "Protein Domain" = Color_2
    )
  ) +
  labs(x = expression(-log[10]*"(P value)")) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, color = "black")
  )

ggsave(
  filename = paste0(file_path, "barplot_OG_glycoprotein_complex_domain.eps"),
  plot = barplot_OG_glycoprotein_complex_domain,
  device = 'eps',
  height = 3, width = 6.5, unit = 'in'
)

