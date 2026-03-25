### TS_diagram.R
### ASGARD 2017 — T-S Diagram with Water Mass Definitions (10 clusters)
### T-Sダイアグラム（水塊定義付き、10クラスター）
###
### REQUIRES (from 00_setup.R, S01, S02):
###   meta_denovo_2  - metadata table (1163 samples × 45 cols)
###   meta_asgard    - survey metadata 181 samples (from S01)
###   clusnum10      - 10-cluster assignments (from S02)
###   cc10           - 10-colour palette (from S02)
###
### OUTPUT:
###   output/survey/TS_diagram_10clusters.pdf

library(tidyverse)
library(ggnewscale)
library(here)

# ==============================================================================
# Section 1: 水塊ボックスの定義 / Water mass definitions (Danielson et al. 2020, Table 2)
# ==============================================================================

wm_boxes <- data.frame(
  name  = c("wCW", "IMW & cCW", "cSW", "wSW", "AnW", "MWW", "WW",
            "AtlW & BBW", "AtlW & BBW"),
  s_min = c(18, 18, 30.8, 30.8, 32.5, 30.8, 30.8, 33.4, 33.8),
  s_max = c(30.8, 30.8, 32.5, 33.4, 33.8, 33.8, 33.8, 35.5, 35.5),
  t_min = c(3, -2.5, 0, 3, 0, -1, -2.5, 3, -2.5),
  t_max = c(14, 3, 3, 14, 3, 0, -1, 14, 3),
  color = c("#E53935", "#1E88E5", "#43A047", "#FB8C00", "#8E24AA",
            "#00ACC1", "#1A237E", "#6D4C41", "#6D4C41"),
  show_label = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

wm_boxes$label_x <- (wm_boxes$s_min + wm_boxes$s_max) / 2
wm_boxes$label_y <- with(wm_boxes, ifelse((t_max - t_min) > 2, t_max - 0.5, (t_min + t_max) / 2))
wm_labels <- wm_boxes[wm_boxes$show_label, ]

# ==============================================================================
# Section 2: 等密度線（σθ）の計算 / Compute sigma-theta contours
# ==============================================================================

s_grid <- seq(17, 36, length.out = 300)
t_grid <- seq(-2.5, 15, length.out = 300)
grid   <- expand.grid(S = s_grid, T = t_grid)

grid$sigma <- with(grid, {
  999.842594 + 6.793952e-2 * T - 9.095290e-3 * T^2 + 1.001685e-4 * T^3 -
    1.120083e-6 * T^4 + 6.536332e-9 * T^5 +
    (0.824493 - 4.0899e-3 * T + 7.6438e-5 * T^2 - 8.2467e-7 * T^3 + 5.3875e-9 * T^4) * S +
    (-5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * T^2) * S^1.5 +
    4.8314e-4 * S^2 - 1000
})

freeze_df <- data.frame(S = seq(17, 36, length.out = 100))
freeze_df$T <- -0.054 * freeze_df$S

# ==============================================================================
# Section 3: T-Sプロット（10クラスター） / T-S diagram with 10 clusters
# ==============================================================================

lvls <- as.character(1:10)

ts_df <- meta_asgard %>% filter(!is.na(temp), !is.na(salinity))
ts_df$cluster <- factor(clusnum10[rownames(ts_df)], levels = lvls)

n_ts <- table(ts_df$cluster)
clbl_ts <- paste0("Cluster ", names(n_ts), " (n=", n_ts, ")")
names(clbl_ts) <- names(n_ts)

dir.create(here("output", "survey"), showWarnings = FALSE, recursive = TRUE)
pdf(file = here::here("output", "survey", "TS_diagram_10clusters.pdf"), width = 14, height = 10)

ggplot() +
  geom_contour(data = grid, aes(x = S, y = T, z = sigma),
               breaks = seq(12, 28, by = 0.5), color = "gray70", linewidth = 0.3) +
  geom_line(data = freeze_df, aes(x = S, y = T),
            linetype = "dashed", color = "black", alpha = 0.5) +
  geom_rect(data = wm_boxes,
            aes(xmin = s_min, xmax = s_max, ymin = t_min, ymax = t_max,
                fill = name, color = name),
            alpha = 0.06, linewidth = 0.8, linetype = "dashed", show.legend = FALSE) +
  scale_fill_manual(values = setNames(wm_boxes$color, wm_boxes$name)) +
  scale_color_manual(values = setNames(wm_boxes$color, wm_boxes$name)) +
  geom_text(data = wm_labels,
            aes(x = label_x, y = label_y, label = name, color = name),
            fontface = "bold", size = 3.5, alpha = 0.5, show.legend = FALSE) +
  new_scale_color() +
  geom_point(data = ts_df, aes(x = salinity, y = temp, color = cluster),
             size = 3, alpha = 0.7) +
  scale_color_manual(values = cc10, labels = clbl_ts) +
  coord_cartesian(xlim = c(24, 35.5), ylim = c(-2.5, 13)) +
  labs(x = "Salinity", y = "Potential Temperature (\u00B0C)",
       title = "T-S Diagram \u2014 ASGARD 2017 Survey (10 clusters)",
       subtitle = "Water mass boxes: Danielson et al. (2020)",
       color = "Cluster") +
  theme_bw(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 13),
    legend.position      = c(0.01, 0.99),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = alpha("white", 0.9))
  )

dev.off()

message("TS_diagram.R: done. PDF: output/survey/TS_diagram_10clusters.pdf")
