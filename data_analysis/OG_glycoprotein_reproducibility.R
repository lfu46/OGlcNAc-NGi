#import packages
packages_names <- c("tidyverse", "writexl", "readxl", "showtext", "psych")
lapply(packages_names, require, character.only = TRUE)

#reproducibility
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(r)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

#generate data frame for Tuni
#tuni1 vs. tuni2
df_tm1_tm2_OG_glycopeptide_HepG2 <- data.frame(
  x = OG_glycopeptide_raw_sl_tmm_HepG2$Tuni_1_sl_tmm,
  y = OG_glycopeptide_raw_sl_tmm_HepG2$Tuni_2_sl_tmm
)

#scatter point plot
#tuni1 vs. tuni2
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

repro_tm1_tm2_OG_glycopeptide_HepG2 <- OG_glycopeptide_raw_sl_tmm_HepG2 |> 
  ggplot(aes(Tuni_1_sl_tmm, Tuni_2_sl_tmm)) +
  geom_point(color = Color_1, size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "gray", linetype = "dashed", linewidth = 0.8) +
  geom_text(x = 500, y = 2500, label = lm_eqn(df_tm1_tm2_OG_glycopeptide_HepG2), parse = TRUE, size = 40) +
  labs(x = "Tuni 1", y = "Tuni 2") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 100),
    axis.title.y = element_text(size = 100)
  )

ggsave(paste0(file_path, "repro_tm1_tm2_OG_glycopeptide_HepG2.png"), plot = repro_tm1_tm2_OG_glycopeptide_HepG2, width = 4, height = 4, units = c("in"), dpi = 600)

#generate data frame for Tuni
#tuni2 vs. tuni3
df_tm2_tm3_OG_glycopeptide_HepG2 <- data.frame(
  x = OG_glycopeptide_raw_sl_tmm_HepG2$Tuni_2_sl_tmm,
  y = OG_glycopeptide_raw_sl_tmm_HepG2$Tuni_3_sl_tmm
)

#scatter point plot
#tuni2 vs. tuni3
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

repro_tm2_tm3_OG_glycopeptide_HepG2 <- OG_glycopeptide_raw_sl_tmm_HepG2 |> 
  ggplot(aes(Tuni_2_sl_tmm, Tuni_3_sl_tmm)) +
  geom_point(color = Color_1, size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "gray", linetype = "dashed", linewidth = 0.8) +
  geom_text(x = 500, y = 2500, label = lm_eqn(df_tm2_tm3_OG_glycopeptide_HepG2), parse = TRUE, size = 40) +
  labs(x = "Tuni 2", y = "Tuni 3") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 100),
    axis.title.y = element_text(size = 100)
  )

ggsave(paste0(file_path, "repro_tm2_tm3_OG_glycopeptide_HepG2.png"), plot = repro_tm2_tm3_OG_glycopeptide_HepG2, width = 4, height = 4, units = c("in"), dpi = 600)

#generate data frame for Tuni
#tuni3 vs. tuni1
df_tm3_tm1_OG_glycopeptide_HepG2 <- data.frame(
  x = OG_glycopeptide_raw_sl_tmm_HepG2$Tuni_3_sl_tmm,
  y = OG_glycopeptide_raw_sl_tmm_HepG2$Tuni_1_sl_tmm
)

#scatter point plot
#tuni3 vs. tuni1
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

repro_tm3_tm1_OG_glycopeptide_HepG2 <- OG_glycopeptide_raw_sl_tmm_HepG2 |> 
  ggplot(aes(Tuni_3_sl_tmm, Tuni_1_sl_tmm)) +
  geom_point(color = Color_1, size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "gray", linetype = "dashed", linewidth = 0.8) +
  geom_text(x = 500, y = 2500, label = lm_eqn(df_tm3_tm1_OG_glycopeptide_HepG2), parse = TRUE, size = 40) +
  labs(x = "Tuni 3", y = "Tuni 1") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 100),
    axis.title.y = element_text(size = 100)
  )

ggsave(paste0(file_path, "repro_tm3_tm1_OG_glycopeptide_HepG2.png"), plot = repro_tm3_tm1_OG_glycopeptide_HepG2, width = 4, height = 4, units = c("in"), dpi = 600)

#generate data frame for Ctrl
#ctrl4 vs. ctrl5
df_ctrl4_ctrl5_OG_glycopeptide_HepG2 <- data.frame(
  x = OG_glycopeptide_raw_sl_tmm_HepG2$Ctrl_4_sl_tmm,
  y = OG_glycopeptide_raw_sl_tmm_HepG2$Ctrl_5_sl_tmm
)

#scatter point plot
#ctrl4 vs. ctrl5
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

repro_ctrl4_ctrl5_OG_glycopeptide_HepG2 <- OG_glycopeptide_raw_sl_tmm_HepG2 |> 
  ggplot(aes(Ctrl_4_sl_tmm, Ctrl_5_sl_tmm)) +
  geom_point(color = Color_2, size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "gray", linetype = "dashed", linewidth = 0.8) +
  geom_text(x = 800, y = 3500, label = lm_eqn(df_ctrl4_ctrl5_OG_glycopeptide_HepG2), parse = TRUE, size = 40) +
  labs(x = "Ctrl 1", y = "Ctrl 2") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 100),
    axis.title.y = element_text(size = 100)
  )

ggsave(paste0(file_path, "repro_ctrl4_ctrl5_OG_glycopeptide_HepG2.png"), plot = repro_ctrl4_ctrl5_OG_glycopeptide_HepG2, width = 4, height = 4, units = c("in"), dpi = 600)

#generate data frame for Ctrl
#ctrl5 vs. ctrl6
df_ctrl5_ctrl6_OG_glycopeptide_HepG2 <- data.frame(
  x = OG_glycopeptide_raw_sl_tmm_HepG2$Ctrl_5_sl_tmm,
  y = OG_glycopeptide_raw_sl_tmm_HepG2$Ctrl_6_sl_tmm
)

#scatter point plot
#ctrl5 vs. ctrl6
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

repro_ctrl5_ctrl6_OG_glycopeptide_HepG2 <- OG_glycopeptide_raw_sl_tmm_HepG2 |> 
  ggplot(aes(Ctrl_5_sl_tmm, Ctrl_6_sl_tmm)) +
  geom_point(color = Color_2, size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "gray", linetype = "dashed", linewidth = 0.8) +
  geom_text(x = 800, y = 3500, label = lm_eqn(df_ctrl5_ctrl6_OG_glycopeptide_HepG2), parse = TRUE, size = 40) +
  labs(x = "Ctrl 2", y = "Ctrl 3") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 100),
    axis.title.y = element_text(size = 100)
  )

ggsave(paste0(file_path, "repro_ctrl5_ctrl6_OG_glycopeptide_HepG2.png"), plot = repro_ctrl5_ctrl6_OG_glycopeptide_HepG2, width = 4, height = 4, units = c("in"), dpi = 600)

#generate data frame for Ctrl
#ctrl6 vs. ctrl4
df_ctrl6_ctrl4_OG_glycopeptide_HepG2 <- data.frame(
  x = OG_glycopeptide_raw_sl_tmm_HepG2$Ctrl_6_sl_tmm,
  y = OG_glycopeptide_raw_sl_tmm_HepG2$Ctrl_4_sl_tmm
)

#scatter point plot
#ctrl6 vs. ctrl4
font_add(family = "Calibri", regular = "calibri.ttf")
showtext_auto()

repro_ctrl6_ctrl4_OG_glycopeptide_HepG2 <- OG_glycopeptide_raw_sl_tmm_HepG2 |> 
  ggplot(aes(Ctrl_6_sl_tmm, Ctrl_4_sl_tmm)) +
  geom_point(color = Color_2, size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "gray", linetype = "dashed", linewidth = 0.8) +
  geom_text(x = 800, y = 3500, label = lm_eqn(df_ctrl6_ctrl4_OG_glycopeptide_HepG2), parse = TRUE, size = 40) +
  labs(x = "Ctrl 3", y = "Ctrl 1") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 100),
    axis.title.y = element_text(size = 100)
  )

ggsave(paste0(file_path, "repro_ctrl6_ctrl4_OG_glycopeptide_HepG2.png"), plot = repro_ctrl6_ctrl4_OG_glycopeptide_HepG2, width = 4, height = 4, units = c("in"), dpi = 600)

#reproducibility for triplicates experiments
#HEK293T
OG_glycoprotein_repro_HEK293T <- OG_glycoprotein_raw_sl_tmm_HEK293T |> 
  mutate(
    Rep1 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    Rep2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    Rep3 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(starts_with("rep"))

pairs.panels(OG_glycoprotein_repro_HEK293T[1:3], lm = TRUE, main = "HEK293T")

#HepG2
OG_glycoprotein_repro_HepG2 <- OG_glycoprotein_raw_sl_tmm_HepG2 |> 
  mutate(
    Rep1 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    Rep2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    Rep3 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(starts_with("Rep"))

pairs.panels(OG_glycoprotein_repro_HepG2[1:3], lm = TRUE, main = "HepG2")

#Jurkat
OG_glycoprotein_repro_Jurkat <- OG_glycoprotein_raw_sl_tmm_Jurkat |> 
  mutate(
    Rep1 = log2(Tuni_1_sl_tmm/Ctrl_4_sl_tmm),
    Rep2 = log2(Tuni_2_sl_tmm/Ctrl_5_sl_tmm),
    Rep3 = log2(Tuni_3_sl_tmm/Ctrl_6_sl_tmm)
  ) |> 
  select(starts_with("Rep"))

pairs.panels(OG_glycoprotein_repro_Jurkat[1:3], lm = TRUE, main = "Jurkat")
