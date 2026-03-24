### TS_diagram.R
### ASGARD 2017 — T-S Diagram with Water Mass Definitions
### T-Sダイアグラム（水塊定義付き）
###
### REQUIRES (from 00_setup.R):
###   meta_denovo_2  - metadata table (1163 samples × 45 cols)
###
### OUTPUT:
###   output/TS_diagram.pdf (3 pages: combined, survey, processing)

library(tidyverse)
library(ggnewscale)
library(here)

# ==============================================================================
# Section 1: ASGARD2017サンプルを抽出 / Extract ASGARD2017 samples with valid T/S
# ==============================================================================

asgard_ts <- meta_denovo_2 %>%
  filter(project == "ASGARD2017", !is.na(temp), !is.na(salinity))

n_survey     <- sum(asgard_ts$station_type == "S")
n_processing <- sum(asgard_ts$station_type == "P")
message("TS_diagram.R: Survey=", n_survey, ", Processing=", n_processing)

asgard_ts$panel <- factor(
  ifelse(asgard_ts$station_type == "S",
         paste0("Survey (0.2 um) (n=", n_survey, ")"),
         paste0("Processing (3 filters) (n=", n_processing, ")")),
  levels = c(paste0("Survey (0.2 um) (n=", n_survey, ")"),
             paste0("Processing (3 filters) (n=", n_processing, ")"))
)

# ==============================================================================
# Section 2: 水塊ボックスの定義 / Water mass definitions (Danielson et al. 2020, Table 2)
#            Mutually exclusive boxes — no overlapping regions
# ==============================================================================

wm_boxes <- data.frame(
  name  = c("wCW",     "IMW & cCW", "cSW",  "wSW",  "AnW",  "MWW",  "WW",
            "AtlW & BBW", "AtlW & BBW"),
  s_min = c(18,        18,          30.8,   30.8,   32.5,   30.8,   30.8,
            33.4,         33.8),
  s_max = c(30.8,      30.8,        32.5,   33.4,   33.8,   33.8,   33.8,
            35.5,         35.5),
  t_min = c(3,         -2.5,        0,      3,      0,      -1,     -2.5,
            3,            -2.5),
  t_max = c(14,        3,           3,      14,     3,      0,      -1,
            14,           3),
  color = c("#E53935", "#1E88E5",   "#43A047", "#FB8C00", "#8E24AA", "#00ACC1", "#1A237E",
            "#6D4C41",    "#6D4C41"),
  # show_label: only label once for AtlW & BBW
  show_label = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

# ラベル位置 / Label positions
wm_boxes$label_x <- (wm_boxes$s_min + wm_boxes$s_max) / 2
wm_boxes$label_y <- with(wm_boxes,
  ifelse((t_max - t_min) > 2, t_max - 0.5, (t_min + t_max) / 2)
)

wm_labels <- wm_boxes[wm_boxes$show_label, ]

# ==============================================================================
# Section 3: 等密度線（σθ）の計算 / Compute sigma-theta contours
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
# Section 4: 共通プロット関数 / Shared plot builder
# ==============================================================================

build_ts_base <- function() {
  ggplot() +
    geom_contour(data = grid, aes(x = S, y = T, z = sigma),
                 breaks = seq(12, 28, by = 0.5),
                 color = "gray70", linewidth = 0.3) +
    geom_line(data = freeze_df, aes(x = S, y = T),
              linetype = "dashed", color = "black", alpha = 0.5) +
    geom_rect(data = wm_boxes,
              aes(xmin = s_min, xmax = s_max, ymin = t_min, ymax = t_max,
                  fill = name, color = name),
              alpha = 0.08, linewidth = 1, linetype = "dashed", show.legend = FALSE) +
    scale_fill_manual(values = setNames(wm_boxes$color, wm_boxes$name)) +
    scale_color_manual(values = setNames(wm_boxes$color, wm_boxes$name)) +
    geom_text(data = wm_labels,
              aes(x = label_x, y = label_y, label = name, color = name),
              fontface = "bold", size = 4, vjust = 0.5, show.legend = FALSE) +
    coord_cartesian(xlim = c(24, 35.5), ylim = c(-2.5, 13)) +
    theme_bw(base_size = 14) +
    theme(
      plot.title    = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 13)
    )
}

# ==============================================================================
# Section 5: PDF出力 / Generate PDF (3 pages)
# ==============================================================================

dir.create(here("output"), showWarnings = FALSE, recursive = TRUE)
pdf(file = here::here("output", "TS_diagram.pdf"), width = 14, height = 10)

# --- Page 1: Combined (both Survey + Processing) ---
p_combined <- build_ts_base() +
  new_scale_color() +
  geom_point(data = asgard_ts,
             aes(x = salinity, y = temp, color = panel, shape = panel),
             size = 2.5, alpha = 0.7, stroke = 0.3) +
  scale_color_manual(values = setNames(c("#2C5F8A", "#E86C00"), levels(asgard_ts$panel))) +
  scale_shape_manual(values = setNames(c(17, 16), levels(asgard_ts$panel))) +
  labs(x = "Salinity", y = "Potential Temperature (\u00B0C)",
       title = "T-S Diagram \u2014 ASGARD 2017",
       subtitle = "Water mass definitions: Danielson et al. (2020)",
       color = NULL, shape = NULL) +
  theme(
    legend.position = c(0.01, 0.99),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.9)),
    legend.text = element_text(size = 12)
  )
print(p_combined)

# --- Page 2: Faceted (Survey and Processing side by side) ---
p_faceted <- build_ts_base() +
  new_scale_color() +
  geom_point(data = asgard_ts,
             aes(x = salinity, y = temp, color = panel, shape = panel),
             size = 2.5, alpha = 0.7, stroke = 0.3, show.legend = FALSE) +
  scale_color_manual(values = setNames(c("#2C5F8A", "#E86C00"), levels(asgard_ts$panel))) +
  scale_shape_manual(values = setNames(c(17, 16), levels(asgard_ts$panel))) +
  facet_wrap(~ panel) +
  labs(x = "Salinity", y = "Potential Temperature (\u00B0C)",
       title = "T-S Diagram \u2014 ASGARD 2017 (by station type)",
       subtitle = "Water mass definitions: Danielson et al. (2020)") +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    strip.background = element_rect(fill = "gray95")
  )
print(p_faceted)

dev.off()

message("TS_diagram.R: done. PDF saved to output/TS_diagram.pdf (2 pages)")
