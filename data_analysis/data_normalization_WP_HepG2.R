#import packages
packages_names <- c("tidyverse", "limma", "edgeR", "showtext")
lapply(packages_names, require, character.only = TRUE)

#HepG2 whole proteome
#check the distribution of intensity of each channel
WP_protein_raw_log2_HepG2 <- WP_protein_raw_HepG2 |> 
  select(Tuni_1:Ctrl_6) |> 
  pivot_longer(cols = c(Tuni_1:Ctrl_6), names_to = "Exp", values_to = "Intensity") |> 
  mutate(Log2_intensity = log2(Intensity))

#boxplot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_WP_protein_raw_log2_HepG2 <- ggplot() +
  geom_boxplot(data = WP_protein_raw_log2_HepG2, 
               aes(x = factor(Exp, levels = c("Tuni_1", "Tuni_2", "Tuni_3", "Ctrl_4", "Ctrl_5", "Ctrl_6")), y = Log2_intensity)) +
  labs(x = "", y = expression(log[2](Intensity))) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 100, color = "black")
  )

ggsave(filename = paste0(file_path, "boxplot_WP_protein_raw_log2_HepG2.png"), plot = boxplot_WP_protein_raw_log2_HepG2, height = 4, width = 4, units = c("in"), dpi = 600)

#density plot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

densityplot_WP_protein_raw_log2_HepG2 <- ggplot() +
  geom_density(data = WP_protein_raw_log2_HepG2,
               aes(x = Log2_intensity, color = factor(Exp, levels = c("Tuni_1", "Tuni_2", "Tuni_3", "Ctrl_4", "Ctrl_5", "Ctrl_6")))) +
  labs(x = expression(log[2](Intensity)), y = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.text = element_text(size = 80),
    legend.title = element_blank()
  )

ggsave(filename = paste0(file_path, "densityplot_WP_protein_raw_log2_HepG2.png"), plot = densityplot_WP_protein_raw_log2_HepG2, height = 4, width = 5, units = c("in"), dpi = 600)

#check the column totals
WP_protein_raw_HepG2 |> select(Tuni_1:Ctrl_6) |> colSums()

#sample loading normalization
target_mean_WP_protein_HepG2 <- mean(colSums(WP_protein_raw_HepG2 |> select(Tuni_1:Ctrl_6)))

norm_facs_WP_protein_HepG2 <- target_mean_WP_protein_HepG2/colSums(WP_protein_raw_HepG2 |> select(Tuni_1:Ctrl_6))

WP_protein_raw_Tuni1_Ctrl6_sl_HepG2 <- sweep(WP_protein_raw_HepG2 |> select(Tuni_1:Ctrl_6), 2, norm_facs_WP_protein_HepG2, FUN = "*")

WP_protein_raw_Tuni1_Ctrl6_sl_tb_HepG2 <- as_tibble(WP_protein_raw_Tuni1_Ctrl6_sl_HepG2)
colnames(WP_protein_raw_Tuni1_Ctrl6_sl_tb_HepG2) <- c("Tuni_1_sl", "Tuni_2_sl", "Tuni_3_sl", "Ctrl_4_sl", "Ctrl_5_sl", "Ctrl_6_sl")

WP_protein_raw_sl_HepG2 <- bind_cols(WP_protein_raw_HepG2, WP_protein_raw_Tuni1_Ctrl6_sl_tb_HepG2)

write_xlsx(WP_protein_raw_sl_HepG2, path = paste0(file_path, "WP_protein_raw_sl_HepG2.xlsx"))

#check the distribution of intensity of each sl channel
WP_protein_raw_sl_log2_HepG2 <- WP_protein_raw_sl_HepG2 |> 
  select(Tuni_1_sl:Ctrl_6_sl) |> 
  pivot_longer(cols = c(Tuni_1_sl:Ctrl_6_sl), names_to = "Exp", values_to = "Intensity") |> 
  mutate(Log2_intensity = log2(Intensity))

#boxplot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_WP_protein_raw_sl_log2_HepG2 <- ggplot() +
  geom_boxplot(data = WP_protein_raw_sl_log2_HepG2, 
               aes(x = factor(Exp, levels = c("Tuni_1_sl", "Tuni_2_sl", "Tuni_3_sl", "Ctrl_4_sl", "Ctrl_5_sl", "Ctrl_6_sl")), y = Log2_intensity)) +
  labs(x = "", y = expression(log[2](Intensity))) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 100, color = "black")
  )

ggsave(filename = paste0(file_path, "boxplot_WP_protein_raw_sl_log2_HepG2.png"), plot = boxplot_WP_protein_raw_sl_log2_HepG2, height = 4, width = 4, units = c("in"), dpi = 600)

#density plot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

densityplot_WP_protein_raw_sl_log2_HepG2 <- ggplot() +
  geom_density(data = WP_protein_raw_sl_log2_HepG2,
               aes(x = Log2_intensity, color = factor(Exp, levels = c("Tuni_1_sl", "Tuni_2_sl", "Tuni_3_sl", "Ctrl_4_sl", "Ctrl_5_sl", "Ctrl_6_sl")))) +
  labs(x = expression(log[2](Intensity)), y = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.text = element_text(size = 80),
    legend.title = element_blank()
  )

ggsave(filename = paste0(file_path, "densityplot_WP_protein_raw_sl_log2_HepG2.png"), plot = densityplot_WP_protein_raw_sl_log2_HepG2, height = 4, width = 5, units = c("in"), dpi = 600)

#check the column sl totals
WP_protein_raw_sl_HepG2 |> select(Tuni_1_sl:Ctrl_6_sl) |> colSums()

#TMM normalization
norm_facs_WP_protein_raw_sl_HepG2 <- calcNormFactors(WP_protein_raw_sl_HepG2 |> select(Tuni_1_sl:Ctrl_6_sl))

WP_protein_raw_sl_Tuni1sl_Ctrl6sl_HepG2 <- sweep(WP_protein_raw_sl_HepG2 |> select(Tuni_1_sl:Ctrl_6_sl), 2, norm_facs_WP_protein_raw_sl_HepG2, FUN = "/")

WP_protein_raw_sl_Tuni1sl_Ctrl6sl_tb_HepG2 <- as_tibble(WP_protein_raw_sl_Tuni1sl_Ctrl6sl_HepG2)
colnames(WP_protein_raw_sl_Tuni1sl_Ctrl6sl_tb_HepG2) <- c("Tuni_1_sl_tmm", "Tuni_2_sl_tmm", "Tuni_3_sl_tmm", "Ctrl_4_sl_tmm", "Ctrl_5_sl_tmm", "Ctrl_6_sl_tmm")

WP_protein_raw_sl_tmm_HepG2 <- bind_cols(WP_protein_raw_sl_HepG2, WP_protein_raw_sl_Tuni1sl_Ctrl6sl_tb_HepG2)

write_xlsx(WP_protein_raw_sl_tmm_HepG2, path = paste0(file_path, "WP_protein_raw_sl_tmm_HepG2.xlsx"))

#check the distribution of intensity of each sl channel
WP_protein_raw_sl_tmm_log2_HepG2 <- WP_protein_raw_sl_tmm_HepG2 |> 
  select(Tuni_1_sl_tmm:Ctrl_6_sl_tmm) |> 
  pivot_longer(cols = c(Tuni_1_sl_tmm:Ctrl_6_sl_tmm), names_to = "Exp", values_to = "Intensity") |> 
  mutate(Log2_intensity = log2(Intensity))

#boxplot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_WP_protein_raw_sl_tmm_log2_HepG2 <- ggplot() +
  geom_boxplot(data = WP_protein_raw_sl_tmm_log2_HepG2, 
               aes(x = factor(Exp, levels = c("Tuni_1_sl_tmm", "Tuni_2_sl_tmm", "Tuni_3_sl_tmm", "Ctrl_4_sl_tmm", "Ctrl_5_sl_tmm", "Ctrl_6_sl_tmm")), y = Log2_intensity)) +
  labs(x = "", y = expression(log[2](Intensity))) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 100, color = "black")
  )

ggsave(filename = paste0(file_path, "boxplot_WP_protein_raw_sl_tmm_log2_HepG2.png"), plot = boxplot_WP_protein_raw_sl_tmm_log2_HepG2, height = 4, width = 4, units = c("in"), dpi = 600)

#density plot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

densityplot_WP_protein_raw_sl_tmm_log2_HepG2 <- ggplot() +
  geom_density(data = WP_protein_raw_sl_tmm_log2_HepG2,
               aes(x = Log2_intensity, color = factor(Exp, levels = c("Tuni_1_sl_tmm", "Tuni_2_sl_tmm", "Tuni_3_sl_tmm", "Ctrl_4_sl_tmm", "Ctrl_5_sl_tmm", "Ctrl_6_sl_tmm")))) +
  labs(x = expression(log[2](Intensity)), y = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.text = element_text(size = 80),
    legend.title = element_blank()
  )

ggsave(filename = paste0(file_path, "densityplot_WP_protein_raw_sl_tmm_log2_HepG2.png"), plot = densityplot_WP_protein_raw_sl_tmm_log2_HepG2, height = 4, width = 5, units = c("in"), dpi = 600)

#check the column sl tmm totals
WP_protein_raw_sl_tmm_HepG2 |> select(Tuni_1_sl_tmm:Ctrl_6_sl_tmm) |> colSums()

#like principle component analysis
plotMDS(WP_protein_raw_sl_tmm_HepG2 |> select(Tuni_1_sl_tmm:Ctrl_6_sl_tmm))

#check coefficients of variance
make_CVs <- function(df) {
  Tuni <- df |> select(Tuni_1_sl_tmm:Tuni_3_sl_tmm)
  Ctrl <- df |> select(Ctrl_4_sl_tmm:Ctrl_6_sl_tmm)
  
  Tuni$ave <- rowMeans(Tuni)
  Tuni$sd <- apply(Tuni[1:3], 1, sd)
  Tuni$cv <- 100 * Tuni$sd / Tuni$ave
  Ctrl$ave <- rowMeans(Ctrl)
  Ctrl$sd <- apply(Ctrl[1:3], 1, sd)
  Ctrl$cv <- 100 * Ctrl$sd / Ctrl$ave
  
  ave_df <- data.frame(Tuni$ave, Ctrl$ave)
  sd_df <- data.frame(Tuni$sd, Ctrl$sd)
  cv_df <- data.frame(Tuni$cv, Ctrl$cv)
  return(list(ave_df, sd_df, cv_df))
}

list_sl_tmm_WP_protein_raw_HepG2 <- make_CVs(WP_protein_raw_sl_tmm_HepG2)

list_sl_tmm_cv_WP_protein_raw_HepG2 <- as_tibble(list_sl_tmm_WP_protein_raw_HepG2[[3]]) |> 
  pivot_longer(ends_with("cv"), names_to = "Exp", values_to = "CV")

font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_sl_tmm_cv_WP_protein_raw_HepG2 <- ggplot() +
  geom_boxplot(data = list_sl_tmm_cv_WP_protein_raw_HepG2, 
               aes(x = factor(Exp, levels = c("Tuni.cv", "Ctrl.cv")), y = CV), outliers = FALSE) +
  labs(x = "", y = "CV (%)") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black")
  )

ggsave(filename = paste0(file_path, "boxplot_sl_tmm_cv_WP_protein_raw_HepG2.png"), plot = boxplot_sl_tmm_cv_WP_protein_raw_HepG2, height = 4, width = 3, units = c("in"), dpi = 600)

#principle component analysis
df_t_WP_protein_HepG2 <- t(WP_protein_raw_sl_tmm_HepG2 |> select(Tuni_1_sl_tmm:Ctrl_6_sl_tmm))

pca_result_WP_protein_HepG2 <- prcomp(df_t_WP_protein_HepG2, scale. = TRUE)

variance_explained_WP_protein_HepG2 <- pca_result_WP_protein_HepG2$sdev^2 / sum(pca_result_WP_protein_HepG2$sdev^2)
percent_variance_WP_protein_HepG2 <- round(variance_explained_WP_protein_HepG2 * 100, 2)

pca_data_WP_protein_HepG2 <- as.data.frame(pca_result_WP_protein_HepG2$x)

Exp_WP_protein_HepG2 <- c("Tuni", "Tuni", "Tuni", "Ctrl", "Ctrl", "Ctrl")

pca_data_WP_protein_HepG2$Exp<- factor(Exp_WP_protein_HepG2, levels = c("Tuni", "Ctrl"))

font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

PCA_WP_protein_HepG2 <- ggplot(pca_data_WP_protein_HepG2, aes(x = PC1, y = PC2, color = Exp)) +
  geom_point(shape = 21, size = 3, stroke = 1.5) +
  labs(
    x = paste0("PC1 (", percent_variance_WP_protein_HepG2[1], "%)"),
    y = paste0("PC2 (", percent_variance_WP_protein_HepG2[2], "%)")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.text = element_text(size = 80),
    legend.title = element_blank()
  )

ggsave(filename = paste0(file_path, "PCA_WP_protein_HepG2.png"), plot = PCA_WP_protein_HepG2, height = 4, width = 4, units = c("in"), dpi = 600)
