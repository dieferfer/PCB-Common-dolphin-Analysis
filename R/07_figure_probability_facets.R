# ------------------------------------------------------------------------------
# P(mean PCB < threshold) for immature dolphins, in the same faceted format used
# for the male projections: facet_grid(Sex ~ Country), dashed line at P = 0.8,
# and the first year reaching it marked with a light green dashed line.
#
# Reads the probabilities produced by immature_projection.R, so the crossing
# years shown here are exactly those of the projection tables.
#
# Outputs: GA_prob_facets_thr17.svg / .png
#          GA_prob_facets_thr41.svg / .png
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(ggplot2); library(svglite) })

p_level    <- 0.8
col_cross  <- "#5FBF5F"          # light green, as in the male-projection figure

pb <- read.csv("output/immature_probability_below.csv")
pb$Sex     <- factor(pb$Sex, levels = c("Males", "Females"))
pb$Country <- factor(pb$Country, levels = c("France", "Spain"))

for (thr in sort(unique(pb$Threshold))) {
  dat <- pb[pb$Threshold == thr, ]

  # first year reaching P >= 0.8 within each facet
  cross <- do.call(rbind, lapply(split(dat, list(dat$Sex, dat$Country), drop = TRUE),
    function(g) {
      ok <- which(g$P_below >= p_level)
      if (!length(ok)) return(NULL)
      data.frame(Sex = g$Sex[1], Country = g$Country[1],
                 year = g$Year[order(g$Year)][min(ok)])
    }))

  p <- ggplot(dat, aes(x = Year, y = P_below)) +
    geom_hline(yintercept = p_level, linetype = "dashed", linewidth = 0.5) +
    geom_vline(data = cross, aes(xintercept = year), colour = col_cross,
               linetype = "dashed", linewidth = 0.6) +
    geom_text(data = cross, aes(x = year, y = 0.02, label = year),
              colour = col_cross, angle = 90, hjust = 0, vjust = -0.4,
              size = 3.6, fontface = "plain") +
    geom_line(linewidth = 0.7) +
    facet_grid(Sex ~ Country) +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(breaks = seq(2020, 2070, by = 10)) +
    labs(x = "Year", y = sprintf("P(mean < %d mg/kg lw)", thr)) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey90", colour = "grey40"),
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank())

  f <- sprintf("output/GA_prob_facets_thr%d", thr)
  ggsave(paste0(f, ".svg"), p, width = 10, height = 4.4, device = svglite::svglite)
  ggsave(paste0(f, ".png"), p, width = 10, height = 4.4, dpi = 200)
  cat("written:", f, ".svg / .png   crossings: ",
      paste(sprintf("%s-%s %d", cross$Country, cross$Sex, cross$year), collapse = " | "),
      "\n", sep = "")
}
