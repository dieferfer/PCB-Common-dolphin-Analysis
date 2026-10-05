# ------------------------------------------------------------------------------
# Figure 9: posterior predictive temporal trend of mean PCB concentrations in
# immature individuals, faceted by sex and coloured by country of stranding.
#
# Same visual grammar as the male projection figure (Figure 6): log10 y axis,
# Dark2 lines over Pastel2 ribbons, red dashed toxicity thresholds, theme_bw.
#
# Input : immature_projection_grid.csv  (written by immature_projection.R)
# Output: GA_figure9_projection_immatures.svg / .png
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(ggplot2); library(svglite) })

g <- read.csv("output/immature_projection_grid.csv")
g$Country_Code <- factor(g$Country_Code, levels = c("France", "Spain"))
g$Sex <- factor(g$Sex, levels = c("Males", "Females"))

p <- ggplot(g, aes(x = Year, y = mu_median, colour = Country_Code)) +
  geom_ribbon(aes(ymin = mu_q025, ymax = mu_q975, fill = Country_Code),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ Sex) +
  geom_hline(yintercept = 17, linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 41, linetype = "dashed", colour = "red") +
  scale_y_log10() +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Pastel2") +
  labs(x = "Year", y = "Predicted PCB concentration (mg/kg lw)", title = "") +
  theme_bw()

ggsave("output/GA_figure9_projection_immatures.svg", p, width = 10, height = 4.2,
       device = svglite::svglite)
ggsave("output/GA_figure9_projection_immatures.png", p, width = 10, height = 4.2, dpi = 200)

cat("written: GA_figure9_projection_immatures.svg / .png\n")
cat(sprintf("rows: %d | years %d-%d | groups: %s\n", nrow(g), min(g$Year), max(g$Year),
            paste(levels(interaction(g$Country_Code, g$Sex, sep = "-")), collapse = ", ")))
print(subset(g, Year %in% c(2020, 2045, 2075))[, c("Year", "Country_Code", "Sex", "mu_median")],
      row.names = FALSE)
