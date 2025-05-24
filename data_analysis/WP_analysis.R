#import packages
packages_names <- c("tidyverse", "writexl", "clusterProfiler", "org.Hs.eg.db", "showtext", "eulerr")
lapply(packages_names, require, character.only = TRUE)

#generate data frame
#HEK293T
proteinList_HEK293T <- WP_protein_Top_tb_HEK293T$logFC
names(proteinList_HEK293T) <- WP_protein_Top_tb_HEK293T$UniprotID
proteinList_HEK293T <- sort(proteinList_HEK293T, decreasing = TRUE)

#HepG2
proteinList_HepG2 <- WP_protein_Top_tb_HepG2$logFC
names(proteinList_HepG2) <- WP_protein_Top_tb_HepG2$UniprotID
proteinList_HepG2 <- sort(proteinList_HepG2, decreasing = TRUE)

#Jurkat
proteinList_Jurkat <- WP_protein_Top_tb_Jurkat$logFC
names(proteinList_Jurkat) <- WP_protein_Top_tb_Jurkat$UniprotID
proteinList_Jurkat <- sort(proteinList_Jurkat, decreasing = TRUE)

#GSEA GO
#HEK293T
GSEA_GO_WP_HEK293T <- gseGO(
  geneList = proteinList_HEK293T,
  ont = "ALL",
  OrgDb = org.Hs.eg.db,
  keyType = "UNIPROT"
)

write_xlsx(GSEA_GO_WP_HEK293T@result, path = paste0(file_path, "GSEA_GO_WP_HEK293T.xlsx"))

#HepG2
GSEA_GO_WP_HepG2 <- gseGO(
  geneList = proteinList_HepG2,
  ont = "ALL",
  OrgDb = org.Hs.eg.db,
  keyType = "UNIPROT"
)

write_xlsx(GSEA_GO_WP_HepG2@result, path = paste0(file_path, "GSEA_GO_WP_HepG2.xlsx"))

#HEK293T
GSEA_GO_WP_Jurkat <- gseGO(
  geneList = proteinList_Jurkat,
  ont = "ALL",
  OrgDb = org.Hs.eg.db,
  keyType = "UNIPROT"
)

write_xlsx(GSEA_GO_WP_Jurkat@result, path = paste0(file_path, "GSEA_GO_WP_Jurkat.xlsx"))

#conbined
GSEA_GO_WP_HEK293T_selected <- bind_rows(
  tibble(GSEA_GO_WP_HEK293T@result |> filter(enrichmentScore > 0) |> head(2)),
  tibble(GSEA_GO_WP_HEK293T@result |> filter(enrichmentScore < 0) |> head(2))
) |> mutate(cell = "HEK293T")

GSEA_GO_WP_HepG2_selected <- bind_rows(
  tibble(GSEA_GO_WP_HepG2@result |> filter(enrichmentScore > 0) |> head(2)),
  tibble(GSEA_GO_WP_HepG2@result |> filter(enrichmentScore < 0) |> head(2))
) |> mutate(cell = "HepG2")

GSEA_GO_WP_Jurkat_selected <- bind_rows(
  tibble(GSEA_GO_WP_Jurkat@result |> filter(enrichmentScore > 0) |> head(2)),
  tibble(GSEA_GO_WP_Jurkat@result |> filter(enrichmentScore < 0) |> head(2))
) |> mutate(cell = "Jurkat")

GSEA_GO_WP_combined <- bind_rows(
  GSEA_GO_WP_HEK293T_selected,
  GSEA_GO_WP_HepG2_selected,
  GSEA_GO_WP_Jurkat_selected
) |> mutate(log_adjP = -log10(p.adjust))

write_xlsx(GSEA_GO_WP_combined, path = paste0(file_path, "GSEA_GO_WP_combined.xlsx"))

#dot plot
font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

dot_plot_GSEA_GO_WP <- GSEA_GO_WP_combined |> 
  ggplot() +
  geom_point(aes(x = enrichmentScore, y = log_adjP, fill = cell, size = setSize),
             shape = 21) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed", linewidth = 0.5) +
  labs(x = "Enrichment Score", y = expression(-log[10]*"("*paste("adjusted ", italic(P), " Value")*")")) +
  scale_fill_manual(
    name = "Cell",
    values = c(
    "HEK293T" = Color_2,
    "HepG2" = Color_3,
    "Jurkat" = Color_4
  )) +
  scale_size(
    name = "Size",
    range = c(3, 6)
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.05, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.01, color = "gray"),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    legend.title = element_text(size = 10, color = "black"),
    legend.key.size = unit(0.01, "in"),
    legend.spacing = unit(0.05, "in")
  )

ggsave(
  filename = paste0(file_path, "dot_plot_GSEA_GO_WP.eps"),
  plot = dot_plot_GSEA_GO_WP,
  device = "eps",
  width = 3,
  height = 2,
  units = "in",
  dpi = 1200
)

#GSEA KEGG
#HEK293T
GSEA_KEGG_WP_HEK293T <- gseKEGG(
  geneList = proteinList_HEK293T,
  organism = 'hsa',
  keyType = 'uniprot'
)

write_xlsx(GSEA_KEGG_WP_HEK293T@result, path = paste0(file_path, "GSEA_KEGG_WP_HEK293T.xlsx"))

#HepG2
GSEA_KEGG_WP_HepG2 <- gseKEGG(
  geneList = proteinList_HepG2,
  organism = 'hsa',
  keyType = 'uniprot'
)

write_xlsx(GSEA_KEGG_WP_HepG2@result, path = paste0(file_path, "GSEA_KEGG_WP_HepG2.xlsx"))

#Jurkat
GSEA_KEGG_WP_Jurkat <- gseKEGG(
  geneList = proteinList_Jurkat,
  organism = 'hsa',
  keyType = 'uniprot'
)

write_xlsx(GSEA_KEGG_WP_Jurkat@result, path = paste0(file_path, "GSEA_KEGG_WP_Jurkat.xlsx"))

#combined
GSEA_KEGG_WP_HEK293T_selected <- tibble(GSEA_KEGG_WP_HEK293T@result) |> 
  filter(Description %in% c("Ribosome", "Protein export", "Ribosome biogenesis in eukaryotes")) |> 
  mutate(cell = "HEK293T")

GSEA_KEGG_WP_HepG2_selected <- tibble(GSEA_KEGG_WP_HepG2@result) |> 
  filter(Description %in% c("Cell cycle", "Aminoacyl-tRNA biosynthesis", "DNA replication", "Protein processing in endoplasmic reticulum")) |> 
  mutate(cell = "HepG2")

GSEA_KEGG_WP_Jurkat_selected <- tibble(GSEA_KEGG_WP_Jurkat@result) |> 
  filter(Description %in% c("Oxidative phosphorylation", "Primary immunodeficiency", "Transcriptional misregulation in cancer", "Ubiquitin mediated proteolysis"))  |> 
  mutate(cell = "Jurkat")

GSEA_KEGG_WP_combined <- bind_rows(
  GSEA_KEGG_WP_HEK293T_selected,
  GSEA_KEGG_WP_HepG2_selected,
  GSEA_KEGG_WP_Jurkat_selected
) |> mutate(log_adjP = -log10(p.adjust))

write_xlsx(GSEA_KEGG_WP_combined, path = paste0(file_path, "GSEA_KEGG_WP_combined.xlsx"))

#dot plot
font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

dot_plot_GSEA_KEGG_WP <- GSEA_KEGG_WP_combined |> 
  ggplot() +
  geom_point(aes(x = enrichmentScore, y = log_adjP, fill = cell, size = setSize),
             shape = 21) +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed", linewidth = 0.5) +
  labs(x = "Enrichment Score", y = expression(-log[10]*"("*paste("adjusted ", italic(P), " Value")*")")) +
  scale_fill_manual(
    name = "Cell",
    values = c(
      "HEK293T" = Color_2,
      "HepG2" = Color_3,
      "Jurkat" = Color_4
    )) +
  scale_size(
    name = "Size",
    range = c(3, 6)
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.05, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.01, color = "gray"),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    legend.title = element_text(size = 10, color = "black"),
    legend.key.size = unit(0.01, "in"),
    legend.spacing = unit(0.05, "in")
  )

ggsave(
  filename = paste0(file_path, "dot_plot_GSEA_KEGG_WP.eps"),
  plot = dot_plot_GSEA_KEGG_WP,
  device = "eps",
  width = 3,
  height = 2,
  units = "in",
  dpi = 1200
)

#ER stress marker
#P11021 & Q6UXH1
#HEK293T
WP_P11021_HEK293T <- WP_protein_Top_tb_HEK293T |> 
  filter(UniprotID == "P11021") |> 
  select(UniprotID, logFC) |> 
  left_join(WP_protein_raw_sl_tmm_HEK293T, by = 'UniprotID') |> 
  select(UniprotID, logFC, Tuni_1_sl_tmm:Ctrl_6_sl_tmm) |> 
  mutate(
    logFC_rep1 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    logFC_rep2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    logFC_rep3 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(UniprotID, starts_with("logFC")) |> 
  mutate(cell = "HEK293T")

WP_Q6UXH1_HEK293T <- WP_protein_Top_tb_HEK293T |> 
  filter(UniprotID == "Q6UXH1") |> 
  select(UniprotID, logFC) |> 
  left_join(WP_protein_raw_sl_tmm_HEK293T, by = 'UniprotID') |> 
  select(UniprotID, logFC, Tuni_1_sl_tmm:Ctrl_6_sl_tmm) |> 
  mutate(
    logFC_rep1 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    logFC_rep2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    logFC_rep3 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(UniprotID, starts_with("logFC")) |> 
  mutate(cell = "HEK293T")

#HepG2
WP_P11021_HepG2 <- WP_protein_Top_tb_HepG2 |> 
  filter(UniprotID == "P11021") |> 
  select(UniprotID, logFC) |> 
  left_join(WP_protein_raw_sl_tmm_HepG2, by = 'UniprotID') |> 
  select(UniprotID, logFC, Tuni_1_sl_tmm:Ctrl_6_sl_tmm) |> 
  mutate(
    logFC_rep1 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    logFC_rep2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    logFC_rep3 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(UniprotID, starts_with("logFC")) |> 
  mutate(cell = "HepG2")

WP_Q6UXH1_HepG2 <- WP_protein_Top_tb_HepG2 |> 
  filter(UniprotID == "Q6UXH1") |> 
  select(UniprotID, logFC) |> 
  left_join(WP_protein_raw_sl_tmm_HepG2, by = 'UniprotID') |> 
  select(UniprotID, logFC, Tuni_1_sl_tmm:Ctrl_6_sl_tmm) |> 
  mutate(
    logFC_rep1 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    logFC_rep2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    logFC_rep3 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(UniprotID, starts_with("logFC")) |> 
  mutate(cell = "HepG2")

#Jurkat
WP_P11021_Jurkat <- WP_protein_Top_tb_Jurkat |> 
  filter(UniprotID == "P11021") |> 
  select(UniprotID, logFC) |> 
  left_join(WP_protein_raw_sl_tmm_Jurkat, by = 'UniprotID') |> 
  select(UniprotID, logFC, Tuni_1_sl_tmm:Ctrl_6_sl_tmm) |> 
  mutate(
    logFC_rep1 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    logFC_rep2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    logFC_rep3 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(UniprotID, starts_with("logFC")) |> 
  mutate(cell = "Jurkat")

WP_Q6UXH1_Jurkat <- WP_protein_Top_tb_Jurkat |> 
  filter(UniprotID == "Q6UXH1") |> 
  select(UniprotID, logFC) |> 
  left_join(WP_protein_raw_sl_tmm_Jurkat, by = 'UniprotID') |> 
  select(UniprotID, logFC, Tuni_1_sl_tmm:Ctrl_6_sl_tmm) |> 
  mutate(
    logFC_rep1 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    logFC_rep2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    logFC_rep3 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(UniprotID, starts_with("logFC")) |> 
  mutate(cell = "Jurkat")

#combined
WP_example_combined <- bind_rows(
  WP_P11021_HEK293T,
  WP_Q6UXH1_HEK293T,
  WP_P11021_HepG2,
  WP_Q6UXH1_HepG2,
  WP_P11021_Jurkat,
  WP_Q6UXH1_Jurkat
)

write_xlsx(WP_example_combined, path = paste0(file_path, "WP_example_combined.xlsx"))

#bar plot
font_add(family = "arial", regular = "arial.ttf")
showtext_auto()

WP_example_barplot <- WP_example_combined |> 
  pivot_longer(cols = logFC_rep1:logFC_rep3, names_to = 'Exp', values_to = 'logFC_rep') |> 
  ggplot() +
  geom_bar(aes(x = cell, y = logFC/3, fill = cell), stat = "identity") +
  geom_point(aes(x = cell, y = logFC_rep)) +
  facet_wrap(vars(UniprotID), scales = "free", ncol = 1,
             labeller = labeller(UniprotID = c("P11021" = "HSPA5",
                                               "Q6UXH1" = "CRELD2"))) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  scale_fill_manual(values = c(
    "HEK293T" = Color_2,
    "HepG2" = Color_3,
    "Jurkat" = Color_4
  )) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    strip.text = element_text(size = 12, color = "black"),
    strip.background = element_rect(fill = NA, color = NA),
    legend.position = "none"
  )

ggsave(
  filename = paste0(file_path, "WP_example_barplot.eps"),
  plot = WP_example_barplot,
  device = "eps",
  width = 2,
  height = 4, 
  units = "in",
  dpi = 1200
)

#significantly upregulated proteins
#HEK293T
WP_upregulated_HEK293T <- WP_protein_Top_tb_HEK293T |> 
  filter(adj.P.Val < 0.05, logFC > 0.5) |> 
  select(UniprotID)

write_xlsx(WP_upregulated_HEK293T, path = paste0(file_path, "WP_upregulated_HEK293T.xlsx"))

#HepG2
WP_upregulated_HepG2 <- WP_protein_Top_tb_HepG2 |> 
  filter(adj.P.Val < 0.05, logFC > 0.5) |> 
  select(UniprotID)

write_xlsx(WP_upregulated_HepG2, path = paste0(file_path, "WP_upregulated_HepG2.xlsx"))

#Jurkat
WP_upregulated_Jurkat <- WP_protein_Top_tb_Jurkat |> 
  filter(adj.P.Val < 0.05, logFC > 0.5) |> 
  select(UniprotID)

write_xlsx(WP_upregulated_Jurkat, path = paste0(file_path, "WP_upregulated_Jurkat.xlsx"))

#euler plot
mat <- c(
  "HEK293T (720)" = 579,
  "HepG2 (178)" = 108,
  "Jurkat (356)" = 221,
  "HEK293T (720)&HepG2 (178)" = 29,
  "HEK293T (720)&Jurkat (356)" = 94,
  "HepG2 (178)&Jurkat (356)" = 23,
  "HepG2 (178)&HEK293T (720)&Jurkat (356)" = 18
)

fit <- euler(mat)

cairo_ps(
  filename = paste0(file_path, "figure2_euler.ps"),
  width = 2, height = 2, fallback_resolution = 1200
)

plot(fit, 
     fill = list(fill = c(Color_2, Color_3, Color_4), alpha = 0.8),
     edges = list(col = "white", lwd = 3),
     labels = list(alpha = c(0)),
     legend = list(fontsize = 0))

dev.off()

#reproducibility
#HEK293T
WP_protein_repro_HEK293T <- WP_protein_raw_sl_tmm_HEK293T |> 
  mutate(
    Rep1 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    Rep2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    Rep3 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(starts_with("rep"))

pairs.panels(WP_protein_repro_HEK293T[1:3], lm = TRUE, main = "HEK293T")

#HepG2
WP_protein_repro_HepG2 <- WP_protein_raw_sl_tmm_HepG2 |> 
  mutate(
    Rep1 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    Rep2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    Rep3 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(starts_with("Rep"))

pairs.panels(WP_protein_repro_HepG2[1:3], lm = TRUE, main = "HepG2")

#Jurkat
WP_protein_repro_Jurkat <- WP_protein_raw_sl_tmm_Jurkat |> 
  mutate(
    Rep1 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    Rep2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    Rep3 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(starts_with("Rep"))

pairs.panels(WP_protein_repro_Jurkat[1:3], lm = TRUE, main = "Jurkat")

