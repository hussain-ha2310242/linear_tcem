# summary.R
# Reporting only: reads saved simulation results from output/tables and writes
# easier ARL summaries for the joint intercept x slope shift study.

library(ggplot2)
library(dplyr)
library(tidyr)

tables_dir <- if (dir.exists("results/tables")) "results/tables" else "output/tables"
plots_dir  <- if (dir.exists("results/plots"))  "results/plots"  else "output/plots"

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)

sigma_eta_levels <- c(0, 0.02, 0.05, 0.08, 0.10)
charts <- c("EWMA-R", "T2")

max_intercept <- 1.00
max_beta <- 0.15

format_arl <- function(x) {
  ifelse(
    is.na(x),
    NA,
    sprintf("%.2f", x)
  )
}

read_shift2d <- function(eta, chart) {
  path <- file.path(tables_dir, paste0("shift2d_arl_eta", eta, "_", chart, ".csv"))
  tab <- read.csv(path, check.names = FALSE)

  tab |>
    pivot_longer(
      cols = -c(chart, intercept),
      names_to = "beta",
      values_to = "ARL"
    ) |>
    mutate(
      sigma_eta = eta,
      beta = as.numeric(beta),
      chart = factor(chart, levels = charts)
    )
}

shift2d <- bind_rows(lapply(sigma_eta_levels, function(eta) {
  bind_rows(lapply(charts, function(chart) read_shift2d(eta, chart)))
}))

shift2d_small <- shift2d |>
  filter(intercept <= max_intercept, beta <= max_beta) |>
  mutate(
    arl_label = format_arl(ARL)
  )

# Quick simulation check against the intended in-control ARL0 = 200.
simulation_check <- shift2d |>
  filter(intercept == 0, beta == 0) |>
  transmute(
    chart,
    sigma_eta,
    ARL,
    target_ARL = 200,
    difference = ARL - target_ARL,
    percent_difference = 100 * difference / target_ARL,
    check = ifelse(abs(percent_difference) <= 5, "OK", "CHECK")
  )

write.csv(
  simulation_check,
  file.path(tables_dir, "summary_simulation_check.csv"),
  row.names = FALSE
)

L_tab <- read.csv(file.path(tables_dir, "calibrated_L.csv"))
UCL_tab <- read.csv(file.path(tables_dir, "calibrated_UCL_T2.csv"))

calibrated_limits <- bind_rows(
  L_tab |>
    transmute(
      sigma_eta,
      chart = "EWMA-R",
      limit_type = "L",
      limit = L_ewma_r,
      achieved_ARL0 = achieved_ARL
    ),
  UCL_tab |>
    transmute(
      sigma_eta,
      chart = "T2",
      limit_type = "UCL T2",
      limit = UCL_T2,
      achieved_ARL0 = achieved_ARL
    )
)

write.csv(
  calibrated_limits,
  file.path(tables_dir, "summary_calibrated_limits.csv"),
  row.names = FALSE
)

# In-control robustness check: keep the sigma_eta = 0 control limits fixed and
# read the saved ARL0 values as multiplicative measurement error increases.
incontrol_eta_path <- file.path(tables_dir, "summary_incontrol_eta_fixed_limits.csv")

if (!file.exists(incontrol_eta_path)) {
  stop(
    "Missing saved robustness table: ", incontrol_eta_path,
    "\nsummary.R does not run simulations. Create this CSV from the simulation script first."
  )
}
incontrol_eta <- read.csv(incontrol_eta_path)

# Long table is easiest for filtering, plotting, and paper/report import.
summary_long <- shift2d_small |>
  select(chart, sigma_eta, intercept, beta, ARL) |>
  arrange(chart, sigma_eta, intercept, beta)

write.csv(
  summary_long,
  file.path(tables_dir, "summary_slope_intercept_arl_long.csv"),
  row.names = FALSE
)

# Wide tables keep ARL values to 2 significant figures for compact reading.
for (chart_name in charts) {
  for (eta in sigma_eta_levels) {
    wide <- summary_long |>
      filter(chart == chart_name, sigma_eta == eta) |>
      mutate(ARL = format_arl(ARL)) |>
      select(intercept, beta, ARL) |>
      pivot_wider(names_from = beta, values_from = ARL)

    write.csv(
      wide,
      file.path(
        tables_dir,
        paste0("summary_slope_intercept_arl_eta", eta, "_", chart_name, ".csv")
      ),
      row.names = FALSE
    )
  }
}

axis_intercept <- c(0, 0.5, 1.0)
axis_beta <- c(0, 0.05, 0.10, 0.15)

plot_data <- shift2d_small |>
  mutate(
    sigma_eta = factor(sigma_eta, levels = sigma_eta_levels),
    ARL_capped = pmin(ARL, 200),
    beta_f = factor(beta, levels = sort(unique(beta))),
    intercept_f = factor(intercept, levels = sort(unique(intercept)))
  )

publication_theme <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey35", linewidth = 0.5),
    strip.background = element_rect(fill = "grey90", color = "grey45"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

arl_fill_scale <- scale_fill_gradientn(
  colors = c("#313695", "#4575b4", "#74add1", "#abd9e9", "#ffffbf", "#fdae61", "#d73027"),
  values = scales::rescale(c(1, 5, 10, 25, 50, 100, 200)),
  limits = c(1, 200),
  name = "ARL"
)

make_grid_plot <- function(data, title, subtitle = NULL) {
  ggplot(data, aes(x = beta, y = intercept, fill = ARL_capped)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = arl_label), size = 2.8, color = "grey10") +
    facet_wrap(~ chart, nrow = 1) +
    scale_x_continuous(breaks = axis_beta, expand = c(0, 0)) +
    scale_y_continuous(breaks = axis_intercept, expand = c(0, 0)) +
    arl_fill_scale +
    labs(
      title = title,
      subtitle = subtitle,
      x = "slope shift (beta)",
      y = "intercept shift (lambda)"
    ) +
    publication_theme
}

# Separate 1 x 2 plots: one file per measurement-error level, two panels for
# the two monitoring techniques.
for (eta in sigma_eta_levels) {
  eta_data <- plot_data |> filter(as.numeric(as.character(sigma_eta)) == eta)
  eta_plot <- make_grid_plot(
    eta_data,
    bquote("Joint slope and intercept shifts, " ~ sigma[eta] ~ "=" ~ .(eta)),
    "Raw ARL values; color scale capped at 200"
  )

  ggsave(
    file.path(plots_dir, paste0("summary_slope_intercept_arl_eta", eta, "_1x2.png")),
    eta_plot,
    width = 7.5,
    height = 3.8,
    dpi = 300
  )
}

# Combined overview retained for comparison across eta values.
summary_plot <- ggplot(plot_data, aes(x = beta, y = intercept, fill = ARL_capped)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = arl_label), size = 2.3, color = "grey10") +
  facet_grid(chart ~ sigma_eta, labeller = label_both) +
  scale_x_continuous(breaks = axis_beta, expand = c(0, 0)) +
  scale_y_continuous(breaks = axis_intercept, expand = c(0, 0)) +
  arl_fill_scale +
  labs(
    title = "Joint slope and intercept shifts",
    subtitle = "Raw ARL values; color scale capped at 200",
    x = "slope shift (beta)",
    y = "intercept shift (lambda)"
  ) +
  publication_theme +
  theme(strip.text = element_text(face = "bold", size = 8))

ggsave(
  file.path(plots_dir, "summary_slope_intercept_arl_grid.png"),
  summary_plot,
  width = 12,
  height = 5.8,
  dpi = 300
)

# Simple line view of ARL versus beta, with fewer axis values.
line_plot <- ggplot(
  shift2d_small,
  aes(x = beta, y = ARL, color = factor(intercept), group = intercept)
) +
  geom_hline(yintercept = 200, linetype = "dashed", color = "grey35") +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.6) +
  facet_grid(chart ~ sigma_eta, labeller = label_both) +
  scale_x_continuous(breaks = axis_beta, limits = c(0, max_beta)) +
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 300), limits = c(0, NA)) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "ARL by slope shift for selected intercept shifts",
    x = "slope shift (beta)",
    y = "ARL",
    color = "intercept"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "grey45"),
    strip.text = element_text(face = "bold", size = 8),
    legend.position = "right"
  )

ggsave(
  file.path(plots_dir, "summary_slope_intercept_arl_lines.png"),
  line_plot,
  width = 12,
  height = 5.8,
  dpi = 300
)

line_colors <- c("EWMA-R" = "#1b9e77", "T2" = "#d95f02")

paper_line_theme <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10, margin = margin(b = 8)),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9, color = "grey20"),
    axis.line = element_line(color = "grey25", linewidth = 0.35),
    axis.ticks = element_line(color = "grey25", linewidth = 0.35),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3),
    strip.background = element_rect(fill = "grey95", color = "grey55", linewidth = 0.35),
    strip.text = element_text(face = "bold", size = 10),
    legend.title = element_blank(),
    legend.position = "top",
    legend.margin = margin(b = 2),
    plot.margin = margin(8, 14, 8, 10)
  )

calibrated_limits <- calibrated_limits |>
  mutate(
    label = sprintf("%.3f", limit),
    label_vjust = ifelse(sigma_eta == max(sigma_eta), 1.35, -0.75),
    label_hjust = ifelse(sigma_eta == min(sigma_eta), -0.05,
                         ifelse(sigma_eta == max(sigma_eta), 1.05, 0.5))
  )

calibrated_limits_plot <- ggplot(
  calibrated_limits,
  aes(x = sigma_eta, y = limit, color = chart, shape = chart)
) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 2.2) +
  geom_text(
    aes(label = label, vjust = label_vjust, hjust = label_hjust),
    size = 2.7,
    show.legend = FALSE
  ) +
  facet_wrap(~ limit_type, scales = "free_y", nrow = 1) +
  scale_x_continuous(breaks = sigma_eta_levels, expand = expansion(mult = c(0.08, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.16))) +
  scale_color_manual(values = line_colors) +
  labs(
    title = expression("Calibrated control limits by " * sigma[eta]),
    subtitle = expression("Target in-control ARL"[0] * " = 200"),
    x = expression(sigma[eta]),
    y = "Control limit",
    color = "chart",
    shape = "chart"
  ) +
  coord_cartesian(clip = "off") +
  paper_line_theme +
  theme(legend.position = "none")

ggsave(
  file.path(plots_dir, "summary_calibrated_limits.png"),
  calibrated_limits_plot,
  width = 6.8,
  height = 3.6,
  dpi = 300
)

incontrol_eta_plot <- ggplot(
  incontrol_eta,
  aes(x = sigma_eta, y = ARL0, color = chart, shape = chart)
) +
  geom_hline(yintercept = 200, linetype = "dashed", color = "grey35", linewidth = 0.45) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  scale_x_continuous(breaks = sigma_eta_levels, expand = expansion(mult = c(0.04, 0.05))) +
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200), limits = c(0, 215))  +
  labs(
    title = expression("Effect of multiplicative measurement error on ARL"[0]),
    x = expression(sigma[eta]),
    y = expression("In-control ARL"[0]),
    shape = "chart"
  ) +
  paper_line_theme

ggsave(
  file.path(plots_dir, "summary_incontrol_eta_fixed_limits.png"),
  incontrol_eta_plot,
  width = 5.8,
  height = 3.9,
  dpi = 300
)

cat("Wrote summary tables and plots.\n")
cat(sprintf("Simulation check saved to %s/summary_simulation_check.csv\n", tables_dir))
