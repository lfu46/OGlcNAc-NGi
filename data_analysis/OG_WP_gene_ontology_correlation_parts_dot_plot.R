#import packages
packages_names <- c("tidyverse", "showtext", "readxl", "ComplexHeatmap", "circlize")
lapply(packages_names, require, character.only = TRUE)

#import data
#HepG2
gene_ontology_OG_glycoprotein_WP_scatterpoint_I_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_WP_sactterplot.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "HepG2_1"
) |> select(Term, Count, P.Value) |> mutate(cell = "HepG2", Group = "I")

gene_ontology_OG_glycoprotein_WP_scatterpoint_II_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_WP_sactterplot.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "HepG2_2"
) |> select(Term, Count, P.Value) |> mutate(cell = "HepG2", Group = "II")

gene_ontology_OG_glycoprotein_WP_scatterpoint_III_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_WP_sactterplot.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "HepG2_3"
) |> select(Term, Count, P.Value) |> mutate(cell = "HepG2", Group = "III")

gene_ontology_OG_glycoprotein_WP_scatterpoint_IV_HepG2 <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_WP_sactterplot.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "HepG2_4"
) |> select(Term, Count, P.Value) |> mutate(cell = "HepG2", Group = "IV")

#HEK293T
gene_ontology_OG_glycoprotein_WP_scatterpoint_I_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_WP_sactterplot.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "HEK293T_1"
) |> select(Term, Count, P.Value) |> mutate(cell = "HEK293T", Group = "I")

gene_ontology_OG_glycoprotein_WP_scatterpoint_II_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_WP_sactterplot.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "HEK293T_2"
) |> select(Term, Count, P.Value) |> mutate(cell = "HEK293T", Group = "II")

gene_ontology_OG_glycoprotein_WP_scatterpoint_III_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_WP_sactterplot.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "HEK293T_3"
) |> select(Term, Count, P.Value) |> mutate(cell = "HEK293T", Group = "III")

gene_ontology_OG_glycoprotein_WP_scatterpoint_IV_HEK293T <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_WP_sactterplot.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "HEK293T_4"
) |> select(Term, Count, P.Value) |> mutate(cell = "HEK293T", Group = "IV")

#Jurkat
gene_ontology_OG_glycoprotein_WP_scatterpoint_I_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_WP_sactterplot.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "Jurkat_1"
) |> select(Term, Count, P.Value) |> mutate(cell = "Jurkat", Group = "I")

gene_ontology_OG_glycoprotein_WP_scatterpoint_II_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_WP_sactterplot.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "Jurkat_2"
) |> select(Term, Count, P.Value) |> mutate(cell = "Jurkat", Group = "II")

gene_ontology_OG_glycoprotein_WP_scatterpoint_III_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_WP_sactterplot.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "Jurkat_3"
) |> select(Term, Count, P.Value) |> mutate(cell = "Jurkat", Group = "III")

gene_ontology_OG_glycoprotein_WP_scatterpoint_IV_Jurkat <- read_xlsx(
  "E:\\Glycosylation_Crosstalk\\Data_Analysis_Sequest\\gene_ontology\\gene_ontology_OG_glycoprotein_WP_sactterplot.xlsx",
  col_names = TRUE,
  .name_repair = "universal",
  sheet = "Jurkat_4"
) |> select(Term, Count, P.Value) |> mutate(cell = "Jurkat", Group = "IV")

#combine
gene_ontology_OG_glycoprotein_WP_scatterpoint_total <- bind_rows(
  gene_ontology_OG_glycoprotein_WP_scatterpoint_I_HepG2,
  gene_ontology_OG_glycoprotein_WP_scatterpoint_II_HepG2,
  gene_ontology_OG_glycoprotein_WP_scatterpoint_III_HepG2,
  gene_ontology_OG_glycoprotein_WP_scatterpoint_IV_HepG2,
  gene_ontology_OG_glycoprotein_WP_scatterpoint_I_HEK293T,
  gene_ontology_OG_glycoprotein_WP_scatterpoint_II_HEK293T,
  gene_ontology_OG_glycoprotein_WP_scatterpoint_III_HEK293T,
  gene_ontology_OG_glycoprotein_WP_scatterpoint_IV_HEK293T,
  gene_ontology_OG_glycoprotein_WP_scatterpoint_I_Jurkat,
  gene_ontology_OG_glycoprotein_WP_scatterpoint_II_Jurkat,
  gene_ontology_OG_glycoprotein_WP_scatterpoint_III_Jurkat,
  gene_ontology_OG_glycoprotein_WP_scatterpoint_IV_Jurkat
)

#generate data frame for heatmap
#HepG2
gene_ontology_OG_glycoprotein_WP_scatterpoint_total_HepG2 <- gene_ontology_OG_glycoprotein_WP_scatterpoint_total |> 
  filter(cell == "HepG2") |> 
  pivot_wider(names_from = c(cell, Group), values_from = P.Value)

gene_ontology_OG_glycoprotein_WP_scatterpoint_total_HEK293T <- gene_ontology_OG_glycoprotein_WP_scatterpoint_total |> 
  filter(cell == "HEK293T") |> 
  pivot_wider(names_from = c(cell, Group), values_from = P.Value)

gene_ontology_OG_glycoprotein_WP_scatterpoint_total_Jurkat <- gene_ontology_OG_glycoprotein_WP_scatterpoint_total |> 
  filter(cell == "Jurkat") |> 
  pivot_wider(names_from = c(cell, Group), values_from = P.Value)

#heatmap group
col_group = c("I" = Color_7, "II" = "#35b779", "III" = "#31688e", "IV" = Color_6)
mat_group <- c("I", "II", "III", "IV")

#heatmap
#HepG2
mat_HepG2 <- data.matrix(gene_ontology_OG_glycoprotein_WP_scatterpoint_total_HepG2 |> select(!Term:Count))
rownames(mat_HepG2) <- gene_ontology_OG_glycoprotein_WP_scatterpoint_total_HepG2$Term

term_number_HepG2 <- gene_ontology_OG_glycoprotein_WP_scatterpoint_total_HepG2$Count
group_number <- gene_ontology_OG_glycoprotein_WP_scatterpoint_total |> 
  filter(cell == "HepG2") |> 
  group_by(cell, Group) |> summarise(Count = sum(Count)) |> pull(Count)

ha_group_top = HeatmapAnnotation(Count = anno_barplot(group_number))
ha_group_bottom = HeatmapAnnotation(Group = mat_group,
                               col = list(Group = col_group))

Heatmap(mat_HepG2, na_col = "transparent", name = "P.Value",
        width = ncol(mat_HepG2)*unit(0.8, "in"),
        col = colorRamp2(c(0, 0.04, 0.05), c("red", "white", "blue")),
        cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE,
        row_names_side = "left", show_column_names = FALSE, 
        top_annotation = ha_group_top, bottom_annotation = ha_group_bottom)

#Dot plot
#HepG2 Part I
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_HepG2_PartI <- gene_ontology_OG_glycoprotein_WP_scatterpoint_I_HepG2 |> 
  mutate(
    Group = "PartI",
    category = paste(cell, Group, sep = "_")
    ) |> 
  ggplot(aes(x = category, y = Term, fill = P.Value, size = Count)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(4, 8)) +
  scale_fill_gradient(low = Color_7, high = "transparent", 
                      limits = c(0, 0.05), breaks = c(0.01, 0.03, 0.05)) +
  scale_x_discrete(labels = c("HepG2 PartI")) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )
  
ggsave(filename = paste0(file_path, "dot_plot_HepG2_PartI.png"), plot = dot_plot_HepG2_PartI, height = 4, width = 6, units = "in", dpi = 600)

#Dot plot
#HepG2 Part II
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_HepG2_PartII <- gene_ontology_OG_glycoprotein_WP_scatterpoint_II_HepG2 |> 
  filter(P.Value < 0.01) |> 
  mutate(
    Group = "PartI",
    category = paste(cell, Group, sep = "_")
  ) |> 
  ggplot(aes(x = category, y = Term, fill = P.Value, size = Count)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(4, 8)) +
  scale_fill_gradient(low = "#35b779", high = "transparent", 
                      limits = c(0, 0.05), breaks = c(0.01, 0.03, 0.05)) +
  scale_x_discrete(labels = c("HepG2 PartII")) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "dot_plot_HepG2_PartII.png"), plot = dot_plot_HepG2_PartII, height = 4, width = 6, units = "in", dpi = 600)

#Dot plot
#HepG2 Part III
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_HepG2_PartIII <- gene_ontology_OG_glycoprotein_WP_scatterpoint_III_HepG2 |> 
  filter(P.Value < 5E-3) |> 
  mutate(
    Group = "PartII",
    category = paste(cell, Group, sep = "_")
  ) |> 
  ggplot(aes(x = category, y = Term, fill = P.Value, size = Count)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(4, 8)) +
  scale_fill_gradient(low = "#31688e", high = "transparent", 
                      limits = c(0, 0.05), breaks = c(0.01, 0.03, 0.05)) +
  scale_x_discrete(labels = c("HepG2 PartIII")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black", lineheight = 0.12),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "dot_plot_HepG2_PartIII.png"), plot = dot_plot_HepG2_PartIII, height = 4, width = 6, units = "in", dpi = 600)

#Dot plot
#HepG2 Part IV
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_HepG2_PartIV <- gene_ontology_OG_glycoprotein_WP_scatterpoint_IV_HepG2 |> 
  filter(P.Value < 0.02) |> 
  mutate(
    Group = "PartIV",
    category = paste(cell, Group, sep = "_")
  ) |> 
  ggplot(aes(x = category, y = Term, fill = P.Value, size = Count)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(4, 8)) +
  scale_fill_gradient(low = Color_6, high = "transparent") +
  scale_x_discrete(labels = c("HepG2 PartIV")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black", lineheight = 0.12),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "dot_plot_HepG2_PartIV.png"), plot = dot_plot_HepG2_PartIV, height = 4, width = 6, units = "in", dpi = 600)

#Dot plot
#HEK293T Part I
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_HEK293T_PartI <- gene_ontology_OG_glycoprotein_WP_scatterpoint_I_HEK293T |> 
  filter(P.Value < 0.03) |> 
  mutate(
    Group = "PartI",
    category = paste(cell, Group, sep = "_")
  ) |> 
  ggplot(aes(x = category, y = Term, fill = P.Value, size = Count)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(4, 8)) +
  scale_fill_gradient(low = Color_7, high = "transparent", limits = c(0, 0.05), breaks = c(0.01, 0.03, 0.05)) +
  scale_x_discrete(labels = c("HEK293T PartI")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "dot_plot_HEK293T_PartI.png"), plot = dot_plot_HEK293T_PartI, height = 4, width = 6, units = "in", dpi = 600)

#Dot plot
#HEK293T Part II
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_HEK293T_PartII <- gene_ontology_OG_glycoprotein_WP_scatterpoint_II_HEK293T |> 
  filter(P.Value < 0.03) |> 
  mutate(
    Group = "PartII",
    category = paste(cell, Group, sep = "_")
  ) |> 
  ggplot(aes(x = category, y = Term, fill = P.Value, size = Count)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(4, 8)) +
  scale_fill_gradient(low = "#35b779", high = "transparent", limits = c(0, 0.05), breaks = c(0.01, 0.03, 0.05)) +
  scale_x_discrete(labels = c("HEK293T PartII")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black", lineheight = 0.12),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "dot_plot_HEK293T_PartII.png"), plot = dot_plot_HEK293T_PartII, height = 4, width = 6, units = "in", dpi = 600)

#Dot plot
#HEK293T Part III
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_HEK293T_PartIII <- gene_ontology_OG_glycoprotein_WP_scatterpoint_III_HEK293T |> 
  mutate(
    Group = "PartIII",
    category = paste(cell, Group, sep = "_")
  ) |> 
  ggplot(aes(x = category, y = Term, fill = P.Value, size = Count)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(4, 8)) +
  scale_fill_gradient(low = "#31688e", high = "transparent", limits = c(0, 0.05), breaks = c(0.01, 0.03, 0.05)) +
  scale_x_discrete(labels = c("HEK293T PartIII")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black", lineheight = 0.12),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "dot_plot_HEK293T_PartIII.png"), plot = dot_plot_HEK293T_PartIII, height = 4, width = 6, units = "in", dpi = 600)

#Dot plot
#HEK293T Part IV
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_HEK293T_PartIV <- gene_ontology_OG_glycoprotein_WP_scatterpoint_IV_HEK293T |> 
  filter(P.Value < 0.01) |> 
  mutate(
    Group = "PartIV",
    category = paste(cell, Group, sep = "_")
  ) |> 
  ggplot(aes(x = category, y = Term, fill = P.Value, size = Count)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(4, 8), breaks = c(20, 40, 60)) +
  scale_fill_gradient(low = Color_6, high = "transparent", limits = c(0, 0.05), breaks = c(0.01, 0.03, 0.05)) +
  scale_x_discrete(labels = c("HEK293T PartIV")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black", lineheight = 0.12),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "dot_plot_HEK293T_PartIV.png"), plot = dot_plot_HEK293T_PartIV, height = 4, width = 6, units = "in", dpi = 600)

#Dot plot
#Jurkat Part I
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_Jurkat_PartI <- gene_ontology_OG_glycoprotein_WP_scatterpoint_I_Jurkat |> 
  filter(P.Value < 0.03) |> 
  mutate(
    Group = "PartI",
    category = paste(cell, Group, sep = "_")
  ) |> 
  ggplot(aes(x = category, y = Term, fill = P.Value, size = Count)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(4, 8)) +
  scale_fill_gradient(low = Color_7, high = "transparent", limits = c(0, 0.05), breaks = c(0.01, 0.03, 0.05)) +
  scale_x_discrete(labels = c("Jurkat PartI")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black", lineheight = 0.12),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "dot_plot_Jurkat_PartI.png"), plot = dot_plot_Jurkat_PartI, height = 4, width = 6, units = "in", dpi = 600)

#Dot plot
#Jurkat Part II
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_Jurkat_PartII <- gene_ontology_OG_glycoprotein_WP_scatterpoint_II_Jurkat |> 
  filter(P.Value < 0.04) |> 
  mutate(
    Group = "PartII",
    category = paste(cell, Group, sep = "_")
  ) |> 
  ggplot(aes(x = category, y = Term, fill = P.Value, size = Count)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(4, 8)) +
  scale_fill_gradient(low = "#35b779", high = "transparent", limits = c(0, 0.05), breaks = c(0.01, 0.03, 0.05)) +
  scale_x_discrete(labels = c("Jurkat PartII")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black", lineheight = 0.12),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "dot_plot_Jurkat_PartII.png"), plot = dot_plot_Jurkat_PartII, height = 4, width = 6, units = "in", dpi = 600)

#Dot plot
#Jurkat Part III
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_Jurkat_PartIII <- gene_ontology_OG_glycoprotein_WP_scatterpoint_III_Jurkat |> 
  filter(P.Value < 0.03) |> 
  mutate(
    Group = "PartIII",
    category = paste(cell, Group, sep = "_")
  ) |> 
  ggplot(aes(x = category, y = Term, fill = P.Value, size = Count)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(4, 8)) +
  scale_fill_gradient(low = "#31688e", high = "transparent", limits = c(0, 0.05), breaks = c(0.01, 0.03, 0.05)) +
  scale_x_discrete(labels = c("Jurkat PartIII")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black", lineheight = 0.12),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "dot_plot_Jurkat_PartIII.png"), plot = dot_plot_Jurkat_PartIII, height = 4, width = 6, units = "in", dpi = 600)

#Dot plot
#Jurkat Part IV
font_add(family = "calibri", regular = "calibri.ttf")
showtext_auto()

dot_plot_Jurkat_PartIV <- gene_ontology_OG_glycoprotein_WP_scatterpoint_IV_Jurkat |> 
  filter(P.Value < 0.035) |> 
  mutate(
    Group = "PartIV",
    category = paste(cell, Group, sep = "_")
  ) |> 
  ggplot(aes(x = category, y = Term, fill = P.Value, size = Count)) +
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(4, 8)) +
  scale_fill_gradient(low = Color_6, high = "transparent", limits = c(0, 0.05), breaks = c(0.01, 0.03, 0.05)) +
  scale_x_discrete(labels = c("Jurkat PartIV")) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black", lineheight = 0.12),
    legend.text = element_text(size = 80, color = "black"),
    legend.title = element_text(size = 80, color = "black")
  )

ggsave(filename = paste0(file_path, "dot_plot_Jurkat_PartIV.png"), plot = dot_plot_Jurkat_PartIV, height = 4, width = 6, units = "in", dpi = 600)
