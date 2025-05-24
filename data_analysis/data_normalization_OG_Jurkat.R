#import packages
packages_names <- c("tidyverse", "limma", "edgeR", "showtext")
lapply(packages_names, require, character.only = TRUE)

#Jurkat glycopeptide
#check the distribution of intensity of each channel
boxplot_OG_glycopeptide_raw_log2_Jurkat <- OG_glycopeptide_raw_Jurkat |> 
  select(Tuni_1:Ctrl_6) |> 
  pivot_longer(cols = c(Tuni_1:Ctrl_6), names_to = "Exp", values_to = "Intensity") |> 
  mutate(Log2_intensity = log2(Intensity))

#boxplot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_OG_glycopeptide_distribution_Jurkat <- ggplot() +
  geom_boxplot(data = boxplot_OG_glycopeptide_raw_log2_Jurkat, 
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

ggsave(filename = paste0(file_path, "boxplot_OG_glycopeptide_distribution_Jurkat.png"), plot = boxplot_OG_glycopeptide_distribution_Jurkat, height = 4, width = 4, units = c("in"), dpi = 600)

#density plot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

densityplot_OG_glycopeptide_distribution_Jurkat <- ggplot() +
  geom_density(data = boxplot_OG_glycopeptide_raw_log2_Jurkat,
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

ggsave(filename = paste0(file_path, "densityplot_OG_glycopeptide_distribution_Jurkat.png"), plot = densityplot_OG_glycopeptide_distribution_Jurkat, height = 4, width = 5, units = c("in"), dpi = 600)

#check the column totals
OG_glycopeptide_raw_Jurkat |> select(Tuni_1:Ctrl_6) |> colSums()

#sample loading normalization
target_mean_glycopeptide_Jurkat <- mean(colSums(OG_glycopeptide_raw_Jurkat |> select(Tuni_1:Ctrl_6)))

norm_facs_glycopeptide_Jurkat <- target_mean_glycopeptide_Jurkat/colSums(OG_glycopeptide_raw_Jurkat |> select(Tuni_1:Ctrl_6))

OG_glycopeptide_raw_Tuni1_Ctrl6_sl_Jurkat <- sweep(OG_glycopeptide_raw_Jurkat |> select(Tuni_1:Ctrl_6), 2, norm_facs_glycopeptide_Jurkat, FUN = "*")

OG_glycopeptide_raw_Tuni1_Ctrl6_sl_tb_Jurkat <- tibble(OG_glycopeptide_raw_Tuni1_Ctrl6_sl_Jurkat)
colnames(OG_glycopeptide_raw_Tuni1_Ctrl6_sl_tb_Jurkat) <- c("Tuni_1_sl", "Tuni_2_sl", "Tuni_3_sl", "Ctrl_4_sl", "Ctrl_5_sl", "Ctrl_6_sl")

OG_glycopeptide_raw_sl_Jurkat <- bind_cols(OG_glycopeptide_raw_Jurkat, OG_glycopeptide_raw_Tuni1_Ctrl6_sl_tb_Jurkat)

write_xlsx(OG_glycopeptide_raw_sl_Jurkat, path = paste0(file_path, "OG_glycopeptide_raw_sl_Jurkat.xlsx"))

#check the distribution of intensity of each sl channel
boxplot_OG_glycopeptide_raw_sl_log2_Jurkat <- OG_glycopeptide_raw_sl_Jurkat |> 
  select(Tuni_1_sl:Ctrl_6_sl) |> 
  pivot_longer(cols = c(Tuni_1_sl:Ctrl_6_sl), names_to = "Exp", values_to = "Intensity") |> 
  mutate(Log2_intensity = log2(Intensity))

#boxplot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_OG_sl_glycopeptide_distribution_Jurkat <- ggplot() +
  geom_boxplot(data = boxplot_OG_glycopeptide_raw_sl_log2_Jurkat, 
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

ggsave(filename = paste0(file_path, "boxplot_OG_sl_glycopeptide_distribution_Jurkat.png"), plot = boxplot_OG_sl_glycopeptide_distribution_Jurkat, height = 4, width = 4, units = c("in"), dpi = 600)

#density plot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

densityplot_OG_sl_glycopeptide_distribution_Jurkat <- ggplot() +
  geom_density(data = boxplot_OG_glycopeptide_raw_sl_log2_Jurkat,
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

ggsave(filename = paste0(file_path, "densityplot_OG_sl_glycopeptide_distribution_Jurkat.png"), plot = densityplot_OG_sl_glycopeptide_distribution_Jurkat, height = 4, width = 5, units = c("in"), dpi = 600)

#check the column sl totals
OG_glycopeptide_raw_sl_Jurkat |> select(Tuni_1_sl:Ctrl_6_sl) |> colSums()

#TMM normalization
sl_tmm_OG_glycopeptide_Jurkat <- calcNormFactors(OG_glycopeptide_raw_sl_Jurkat |> select(Tuni_1_sl:Ctrl_6_sl))

OG_glycopeptide_raw_sl_Tuni1sl_Ctrl6sl_Jurkat <- sweep(OG_glycopeptide_raw_sl_Jurkat |> select(Tuni_1_sl:Ctrl_6_sl), 2, sl_tmm_OG_glycopeptide_Jurkat, FUN = "/")

OG_glycopeptide_raw_sl_Tuni1sl_Ctrl6sl_tb_Jurkat <- as_tibble(OG_glycopeptide_raw_sl_Tuni1sl_Ctrl6sl_Jurkat)
colnames(OG_glycopeptide_raw_sl_Tuni1sl_Ctrl6sl_tb_Jurkat) <- c("Tuni_1_sl_tmm", "Tuni_2_sl_tmm", "Tuni_3_sl_tmm", "Ctrl_4_sl_tmm", "Ctrl_5_sl_tmm", "Ctrl_6_sl_tmm")

OG_glycopeptide_raw_sl_tmm_Jurkat <- bind_cols(OG_glycopeptide_raw_sl_Jurkat, OG_glycopeptide_raw_sl_Tuni1sl_Ctrl6sl_tb_Jurkat)

write_xlsx(OG_glycopeptide_raw_sl_tmm_Jurkat, path = paste0(file_path, "OG_glycopeptide_raw_sl_tmm_Jurkat.xlsx"))

#check the distribution of intensity of each sl channel
boxplot_OG_glycopeptide_raw_sl_tmm_log2_Jurkat <- OG_glycopeptide_raw_sl_tmm_Jurkat |> 
  select(Tuni_1_sl_tmm:Ctrl_6_sl_tmm) |> 
  pivot_longer(cols = c(Tuni_1_sl_tmm:Ctrl_6_sl_tmm), names_to = "Exp", values_to = "Intensity") |> 
  mutate(Log2_intensity = log2(Intensity))

#boxplot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_OG_sl_tmm_glycopeptide_distribution_Jurkat <- ggplot() +
  geom_boxplot(data = boxplot_OG_glycopeptide_raw_sl_tmm_log2_Jurkat, 
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

ggsave(filename = paste0(file_path, "boxplot_OG_sl_tmm_glycopeptide_distribution_Jurkat.png"), plot = boxplot_OG_sl_tmm_glycopeptide_distribution_Jurkat, height = 4, width = 4, units = c("in"), dpi = 600)

#density plot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

densityplot_OG_sl_tmm_glycopeptide_distribution_Jurkat <- ggplot() +
  geom_density(data = boxplot_OG_glycopeptide_raw_sl_tmm_log2_Jurkat,
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

ggsave(filename = paste0(file_path, "densityplot_OG_sl_tmm_glycopeptide_distribution_Jurkat.png"), plot = densityplot_OG_sl_tmm_glycopeptide_distribution_Jurkat, height = 4, width = 5, units = c("in"), dpi = 600)

#check the column sl tmm totals
OG_glycopeptide_raw_sl_tmm_Jurkat |> select(Tuni_1_sl_tmm:Ctrl_6_sl_tmm) |> colSums()

#like principle component analysis
plotMDS(OG_glycopeptide_raw_sl_tmm_Jurkat |> select(Tuni_1_sl_tmm:Ctrl_6_sl_tmm))

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

list_sl_tmm_OG_glycopeptide_Jurkat <- make_CVs(OG_glycopeptide_raw_sl_tmm_Jurkat)

list_sl_tmm_cv_OG_glycopeptide_Jurkat <- as_tibble(list_sl_tmm_OG_glycopeptide_Jurkat[[3]]) |> 
  pivot_longer(ends_with("cv"), names_to = "Exp", values_to = "CV")

font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_sl_tmm_cv_OG_glycopeptide_Jurkat <- ggplot() +
  geom_boxplot(data = list_sl_tmm_cv_OG_glycopeptide_Jurkat, 
               aes(x = factor(Exp, levels = c("Tuni.cv", "Ctrl.cv")), y = CV), outliers = FALSE) +
  labs(x = "", y = "CV (%)") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 100, color = "black")
  )

ggsave(filename = paste0(file_path, "boxplot_sl_tmm_cv_OG_glycopeptide_Jurkat.png"), plot = boxplot_sl_tmm_cv_OG_glycopeptide_Jurkat, height = 4, width = 4, units = c("in"), dpi = 600)

#principle component analysis
df_t_OG_glycopeptide_Jurkat <- t(OG_glycopeptide_raw_sl_tmm_Jurkat |> 
                                    select(Tuni_1_sl_tmm:Ctrl_6_sl_tmm))

pca_result_glycopeptide_Jurkat <- prcomp(df_t_OG_glycopeptide_Jurkat, scale. = TRUE)

variance_explained_glycopeptide_Jurkat <- pca_result_glycopeptide_Jurkat$sdev^2 / sum(pca_result_glycopeptide_Jurkat$sdev^2)
percent_variance_glycopeptide_Jurkat <- round(variance_explained_glycopeptide_Jurkat * 100, 2)

pca_data_glycopeptide_Jurkat <- as.data.frame(pca_result_glycopeptide_Jurkat$x)

Exp_glycopeptide_Jurkat <- c("Tuni", "Tuni", "Tuni", "Ctrl", "Ctrl", "Ctrl")

pca_data_glycopeptide_Jurkat$Exp<- factor(Exp_glycopeptide_Jurkat, levels = c("Tuni", "Ctrl"))

font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

PCA_OG_glycopeptide_Jurkat <- ggplot(pca_data_glycopeptide_Jurkat, aes(x = PC1, y = PC2, color = Exp)) +
  geom_point(shape = 21, size = 3, stroke = 1.5) +
  labs(
    x = paste0("PC1 (", percent_variance_glycopeptide_Jurkat[1], "%)"),
    y = paste0("PC2 (", percent_variance_glycopeptide_Jurkat[2], "%)")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.text = element_text(size = 100),
    legend.title = element_blank()
  )

ggsave(filename = paste0(file_path, "PCA_OG_glycopeptide_Jurkat.png"), plot = PCA_OG_glycopeptide_Jurkat, height = 4, width = 4, units = c("in"), dpi = 600)

#Jurkat glycoprotein
#check the distribution of intensity of each channel
boxplot_OG_glycoprotein_raw_log2_Jurkat <- OG_glycoprotein_raw_Jurkat |> 
  select(Tuni_1:Ctrl_6) |> 
  pivot_longer(cols = c(Tuni_1:Ctrl_6), names_to = "Exp", values_to = "Intensity") |> 
  mutate(Log2_intensity = log2(Intensity))

#boxplot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_OG_glycoprotein_distribution_Jurkat <- ggplot() +
  geom_boxplot(data = boxplot_OG_glycoprotein_raw_log2_Jurkat, 
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

ggsave(filename = paste0(file_path, "boxplot_OG_glycoprotein_distribution_Jurkat.png"), plot = boxplot_OG_glycoprotein_distribution_Jurkat, height = 4, width = 4, units = c("in"), dpi = 600)

#density plot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

densityplot_OG_glycoprotein_distribution_Jurkat <- ggplot() +
  geom_density(data = boxplot_OG_glycoprotein_raw_log2_Jurkat,
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

ggsave(filename = paste0(file_path, "densityplot_OG_glycoprotein_distribution_Jurkat.png"), plot = densityplot_OG_glycoprotein_distribution_Jurkat, height = 4, width = 5, units = c("in"), dpi = 600)

#check the column totals
OG_glycoprotein_raw_Jurkat |> select(Tuni_1:Ctrl_6) |> colSums()

#sample loading normalization
target_mean_glycoprotein_Jurkat <- mean(colSums(OG_glycoprotein_raw_Jurkat |> select(Tuni_1:Ctrl_6)))

norm_facs_glycoprotein_Jurkat <- target_mean_glycoprotein_Jurkat/colSums(OG_glycoprotein_raw_Jurkat |> select(Tuni_1:Ctrl_6))

OG_glycoprotein_raw_Tuni1_Ctrl6_sl_Jurkat <- sweep(OG_glycoprotein_raw_Jurkat |> select(Tuni_1:Ctrl_6), 2, norm_facs_glycoprotein_Jurkat, FUN = "*")

OG_glycoprotein_raw_Tuni1_Ctrl6_sl_tb_Jurkat <- tibble(OG_glycoprotein_raw_Tuni1_Ctrl6_sl_Jurkat)
colnames(OG_glycoprotein_raw_Tuni1_Ctrl6_sl_tb_Jurkat) <- c("Tuni_1_sl", "Tuni_2_sl", "Tuni_3_sl", "Ctrl_4_sl", "Ctrl_5_sl", "Ctrl_6_sl")

OG_glycoprotein_raw_sl_Jurkat <- bind_cols(OG_glycoprotein_raw_Jurkat, OG_glycoprotein_raw_Tuni1_Ctrl6_sl_tb_Jurkat)

write_xlsx(OG_glycoprotein_raw_sl_Jurkat, path = paste0(file_path, "OG_glycoprotein_raw_sl_Jurkat.xlsx"))

#check the distribution of intensity of each sl channel
boxplot_OG_glycoprotein_raw_sl_log2_Jurkat <- OG_glycoprotein_raw_sl_Jurkat |> 
  select(Tuni_1_sl:Ctrl_6_sl) |> 
  pivot_longer(cols = c(Tuni_1_sl:Ctrl_6_sl), names_to = "Exp", values_to = "Intensity") |> 
  mutate(Log2_intensity = log2(Intensity))

#boxplot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_OG_sl_glycoprotein_distribution_Jurkat <- ggplot() +
  geom_boxplot(data = boxplot_OG_glycoprotein_raw_sl_log2_Jurkat, 
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

ggsave(filename = paste0(file_path, "boxplot_OG_sl_glycoprotein_distribution_Jurkat.png"), plot = boxplot_OG_sl_glycoprotein_distribution_Jurkat, height = 4, width = 4, units = c("in"), dpi = 600)

#density plot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

densityplot_OG_sl_glycoprotein_distribution_Jurkat <- ggplot() +
  geom_density(data = boxplot_OG_glycoprotein_raw_sl_log2_Jurkat,
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

ggsave(filename = paste0(file_path, "densityplot_OG_sl_glycoprotein_distribution_Jurkat.png"), plot = densityplot_OG_sl_glycoprotein_distribution_Jurkat, height = 4, width = 5, units = c("in"), dpi = 600)

#check the column sl totals
OG_glycoprotein_raw_sl_Jurkat |> select(Tuni_1_sl:Ctrl_6_sl) |> colSums()

#TMM normalization
sl_tmm_OG_glycoprotein_Jurkat <- calcNormFactors(OG_glycoprotein_raw_sl_Jurkat |> select(Tuni_1_sl:Ctrl_6_sl))

OG_glycoprotein_raw_sl_Tuni1sl_Ctrl6sl_Jurkat <- sweep(OG_glycoprotein_raw_sl_Jurkat |> select(Tuni_1_sl:Ctrl_6_sl), 2, sl_tmm_OG_glycoprotein_Jurkat, FUN = "/")

OG_glycoprotein_raw_sl_Tuni1sl_Ctrl6sl_tb_Jurkat <- as_tibble(OG_glycoprotein_raw_sl_Tuni1sl_Ctrl6sl_Jurkat)
colnames(OG_glycoprotein_raw_sl_Tuni1sl_Ctrl6sl_tb_Jurkat) <- c("Tuni_1_sl_tmm", "Tuni_2_sl_tmm", "Tuni_3_sl_tmm", "Ctrl_4_sl_tmm", "Ctrl_5_sl_tmm", "Ctrl_6_sl_tmm")

OG_glycoprotein_raw_sl_tmm_Jurkat <- bind_cols(OG_glycoprotein_raw_sl_Jurkat, OG_glycoprotein_raw_sl_Tuni1sl_Ctrl6sl_tb_Jurkat)

write_xlsx(OG_glycoprotein_raw_sl_tmm_Jurkat, path = paste0(file_path, "OG_glycoprotein_raw_sl_tmm_Jurkat.xlsx"))

#check the distribution of intensity of each sl channel
boxplot_OG_glycoprotein_raw_sl_tmm_log2_Jurkat <- OG_glycoprotein_raw_sl_tmm_Jurkat |> 
  select(Tuni_1_sl_tmm:Ctrl_6_sl_tmm) |> 
  pivot_longer(cols = c(Tuni_1_sl_tmm:Ctrl_6_sl_tmm), names_to = "Exp", values_to = "Intensity") |> 
  mutate(Log2_intensity = log2(Intensity))

#boxplot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_OG_sl_tmm_glycoprotein_distribution_Jurkat <- ggplot() +
  geom_boxplot(data = boxplot_OG_glycoprotein_raw_sl_tmm_log2_Jurkat, 
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

ggsave(filename = paste0(file_path, "boxplot_OG_sl_tmm_glycoprotein_distribution_Jurkat.png"), plot = boxplot_OG_sl_tmm_glycoprotein_distribution_Jurkat, height = 4, width = 4, units = c("in"), dpi = 600)

#density plot
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

densityplot_OG_sl_tmm_glycoprotein_distribution_Jurkat <- ggplot() +
  geom_density(data = boxplot_OG_glycoprotein_raw_sl_tmm_log2_Jurkat,
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

ggsave(filename = paste0(file_path, "densityplot_OG_sl_tmm_glycoprotein_distribution_Jurkat.png"), plot = densityplot_OG_sl_tmm_glycoprotein_distribution_Jurkat, height = 4, width = 5, units = c("in"), dpi = 600)

#check the column sl tmm totals
OG_glycoprotein_raw_sl_tmm_Jurkat |> select(Tuni_1_sl_tmm:Ctrl_6_sl_tmm) |> colSums()

#like principle component analysis
plotMDS(OG_glycoprotein_raw_sl_tmm_Jurkat |> select(Tuni_1_sl_tmm:Ctrl_6_sl_tmm))

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

list_sl_tmm_OG_glycoprotein_Jurkat <- make_CVs(OG_glycoprotein_raw_sl_tmm_Jurkat)

list_sl_tmm_cv_OG_glycoprotein_Jurkat <- as_tibble(list_sl_tmm_OG_glycoprotein_Jurkat[[3]]) |> 
  pivot_longer(ends_with("cv"), names_to = "Exp", values_to = "CV")

font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

boxplot_sl_tmm_cv_OG_glycoprotein_Jurkat <- ggplot() +
  geom_boxplot(data = list_sl_tmm_cv_OG_glycoprotein_Jurkat, 
               aes(x = factor(Exp, levels = c("Tuni.cv", "Ctrl.cv")), y = CV), outliers = FALSE) +
  labs(x = "", y = "CV (%)") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 100, color = "black")
  )

ggsave(filename = paste0(file_path, "boxplot_sl_tmm_cv_OG_glycoprotein_Jurkat.png"), plot = boxplot_sl_tmm_cv_OG_glycoprotein_Jurkat, height = 4, width = 4, units = c("in"), dpi = 600)

#principle component analysis
df_t_OG_glycoprotein_Jurkat <- t(OG_glycoprotein_raw_sl_tmm_Jurkat |> 
                                    select(Tuni_1_sl_tmm:Ctrl_6_sl_tmm))

pca_result_glycoprotein_Jurkat <- prcomp(df_t_OG_glycoprotein_Jurkat, scale. = TRUE)

variance_explained_glycoprotein_Jurkat <- pca_result_glycoprotein_Jurkat$sdev^2 / sum(pca_result_glycoprotein_Jurkat$sdev^2)
percent_variance_glycoprotein_Jurkat <- round(variance_explained_glycoprotein_Jurkat * 100, 2)

pca_data_glycoprotein_Jurkat <- as.data.frame(pca_result_glycoprotein_Jurkat$x)

Exp_glycoprotein_Jurkat <- c("Tuni", "Tuni", "Tuni", "Ctrl", "Ctrl", "Ctrl")

pca_data_glycoprotein_Jurkat$Exp<- factor(Exp_glycoprotein_Jurkat, levels = c("Tuni", "Ctrl"))

font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

PCA_OG_glycoprotein_Jurkat <- ggplot(pca_data_glycoprotein_Jurkat, aes(x = PC1, y = PC2, color = Exp)) +
  geom_point(shape = 21, size = 3, stroke = 1.5) +
  labs(
    x = paste0("PC1 (", percent_variance_glycoprotein_Jurkat[1], "%)"),
    y = paste0("PC2 (", percent_variance_glycoprotein_Jurkat[2], "%)")
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "gray", linewidth = 0.2),
    panel.grid.minor = element_line(color = "gray", linewidth = 0.1),
    axis.title = element_text(size = 100),
    axis.text.x = element_text(size = 100, color = "black"),
    axis.text.y = element_text(size = 100, color = "black"),
    legend.text = element_text(size = 100),
    legend.title = element_blank()
  )

ggsave(filename = paste0(file_path, "PCA_OG_glycoprotein_Jurkat.png"), plot = PCA_OG_glycoprotein_Jurkat, height = 4, width = 4, units = c("in"), dpi = 600)
