# ------------------------------------------------------------------------------
# Figure 3: density of Aroclor 1254-equivalent PCB concentrations by country,
# maturity class and sex.
#
# Reproduces the layout of the Inkscape-retouched version ("output/density plot.svg"):
# a small per-panel key with a coloured swatch and "N = x", no global legend,
# and the axis label "Aroclor 1254-equivalents (mg / kg l.w.)".
#
# Two errors of the previous version are corrected:
#   1. Colours were swapped. The palette used elsewhere in the paper is
#      Males = #55A868 (green), Females = #4C72B0 (blue),
#      Pregnant females = #FDD9A2; the old figure had males in blue.
#   2. The x axis was truncated with a scale limit, which DROPS the rows outside
#      the range: the four most contaminated French males (288.4, 257.7, 180.8
#      and 154.4 mg/kg l.w.) were excluded from the density estimate and from
#      the panel counts, which summed to 179 instead of 183. The truncation now
#      uses coord_cartesian(), so those animals still inform the curves.
#
# Output: GA_figure3_density.svg / .png
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(ggplot2); library(svglite) })

X_MAX <- 150
ALPHA <- 0.35

d <- read.csv("data/PCBcorrectedDB.csv", sep = ";")
d$Country <- factor(ifelse(d$Country_Code == 4, "France", "Spain"),
                    levels = c("France", "Spain"))
d$Sex_plot <- factor(c("Males", "Females", "Pregnant females")[d$Sex],
                     levels = c("Males", "Females", "Pregnant females"))
d$Maturity_plot <- factor(ifelse(d$Maturity == 0, "Immatures", "Matures"),
                          levels = c("Immatures", "Matures"))

pal      <- c("Males" = "#55A868", "Females" = "#4C72B0", "Pregnant females" = "#FDD9A2")
# the pale orange is invisible as a 0.7 pt line, so outlines use a darker tone
pal_line <- c(pal["Males"], pal["Females"], "Pregnant females" = "#D9963A")

# ---- per-panel key: counts over ALL individuals ------------------------------
key <- as.data.frame(table(d$Country, d$Maturity_plot, d$Sex_plot))
names(key) <- c("Country", "Maturity_plot", "Sex_plot", "n")
key <- key[key$n > 0, ]
key$Sex_plot <- factor(as.character(key$Sex_plot), levels = names(pal))
key <- key[order(key$Country, key$Maturity_plot, key$Sex_plot), ]
key$k <- ave(seq_len(nrow(key)), key$Country, key$Maturity_plot, FUN = seq_along)
stopifnot(sum(key$n) == nrow(d))

# highest density value across every group, to place the key near the top
ymax <- max(vapply(split(d$ICES73, list(d$Country, d$Maturity_plot, d$Sex_plot), drop = TRUE),
                   function(v) if (length(v) > 1) max(density(v, adjust = 1)$y) else 0, numeric(1)))
row_y <- function(k) ymax * (1.02 - 0.085 * (k - 1))
key$y  <- row_y(key$k)
key$xr <- X_MAX * 0.735          # swatch left edge
key$xt <- X_MAX * 0.80           # text left edge

cat("panel counts sum to", sum(key$n), "individuals\n")
out <- subset(d, ICES73 > X_MAX)
cat(sprintf("outside the plotted range but kept in the density estimate: %d - %s\n",
            nrow(out), paste(sprintf("%s %s %s %.1f", out$Country, out$Maturity_plot,
                                     out$Sex_plot, out$ICES73), collapse = "; ")))

p <- ggplot(d, aes(x = ICES73, fill = Sex_plot, colour = Sex_plot)) +
  geom_density(alpha = ALPHA, linewidth = 0.7, adjust = 1) +
  geom_rug(alpha = 0.5, linewidth = 0.3) +
  # the caption distinguishes the two thresholds by tone: 17 light, 41 dark
  geom_vline(xintercept = 17, linetype = "dashed", colour = "#E06666", linewidth = 0.6) +
  geom_vline(xintercept = 41, linetype = "dashed", colour = "#8B0000", linewidth = 0.6) +
  # per-panel key, as in the retouched figure
  geom_rect(data = key, inherit.aes = FALSE,
            aes(xmin = xr, xmax = xr + X_MAX * 0.038,
                ymin = y - ymax * 0.028, ymax = y + ymax * 0.028,
                fill = Sex_plot, colour = Sex_plot),
            alpha = ALPHA, linewidth = 0.4, show.legend = FALSE) +
  geom_text(data = key, inherit.aes = FALSE,
            aes(x = xt, y = y, label = sprintf("N = %d", n)),
            hjust = 0, vjust = 0.5, size = 3.2, colour = "grey15") +
  facet_grid(Maturity_plot ~ Country) +
  scale_fill_manual(values = pal) +
  scale_colour_manual(values = pal_line) +
  coord_cartesian(xlim = c(0, X_MAX), ylim = c(0, ymax * 1.06)) +
  labs(x = "Aroclor 1254-equivalents (mg / kg l.w.)", y = "density") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey90", colour = "grey40"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave("output/GA_figure3_density.svg", p, width = 10, height = 4.6, device = svglite::svglite)
ggsave("output/GA_figure3_density.png", p, width = 10, height = 4.6, dpi = 200)
cat("written: GA_figure3_density.svg / .png\n")
