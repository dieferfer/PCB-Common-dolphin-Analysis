# ------------------------------------------------------------------------------
# Immature model (M1): partial effect of Year on the PCB scale, overall and
# split by country and by sex, plus a summary table.
#
#   log(mu) = b0 + b1*Country_Spain + b2*Sex_Females + f(Year), rw2
#
# Country and Sex are log-additive shifts, so every curve shares the same rw2
# shape and differs only by a constant factor. Ribbons come from the joint
# posterior (inla.posterior.sample), not from summing marginal intervals.
#
# Outputs: GA_year_effect_immatures.svg          (overall)
#          GA_year_effect_immatures_country.svg  (France vs Spain)
#          GA_year_effect_immatures_sex.svg      (Males vs Females)
#          GA_immature_summary_table.svg / .tex
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(INLA) })
n_draws <- 4000; set.seed(1)

d <- read.csv("data/PCBcorrectedDB.csv", sep = ";")
d$Country_Code <- factor(ifelse(d$Country_Code == 4, "France", "Spain"),
                         levels = c("France", "Spain"))
d$Sex <- factor(ifelse(d$Sex == 1, "Males", "Females"), levels = c("Males", "Females"))
d <- d[d$Maturity == 0, ]

fit <- inla(ICES73 ~ Country_Code + Sex + f(Year, model = "rw2"),
            family = "gamma", data = d,
            control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE))

# ---- joint-posterior draws of the year curve, at chosen covariate levels -----
s   <- inla.posterior.sample(n_draws, fit)
rn  <- rownames(s[[1]]$latent)
i_b0 <- grep("^\\(Intercept\\)", rn)
i_sp <- grep("^Country_CodeSpain", rn)
i_fe <- grep("^SexFemales", rn)
i_yr <- grep("^Year:", rn)

lat  <- vapply(s, function(x) x$latent[, 1], numeric(length(rn)))
b0   <- lat[i_b0, ]; b_sp <- lat[i_sp, ]; b_fe <- lat[i_fe, ]
fyr  <- lat[i_yr, , drop = FALSE]
yrs  <- fit$summary.random$Year$ID
base <- exp(sweep(fyr, 2, b0, "+"))                 # exp(b0 + f(Year))

is_sp <- as.numeric(d$Country_Code == "Spain")
is_fe <- as.numeric(d$Sex == "Females")

summ <- function(offset) {                          # offset: length n_draws
  pred <- base * rep(offset, each = length(yrs))
  data.frame(year = yrs,
             med = apply(pred, 1, median),
             lo  = apply(pred, 1, quantile, 0.025),
             hi  = apply(pred, 1, quantile, 0.975))
}

# overall: average over the observed composition
est_all    <- summ(colMeans(exp(outer(is_sp, b_sp) + outer(is_fe, b_fe))))
# by country: vary country, average sex
est_france <- summ(colMeans(exp(outer(is_fe, b_fe))))
est_spain  <- summ(colMeans(exp(outer(is_fe, b_fe))) * exp(b_sp))
# by sex: vary sex, average country
est_males  <- summ(colMeans(exp(outer(is_sp, b_sp))))
est_females<- summ(colMeans(exp(outer(is_sp, b_sp))) * exp(b_fe))

# ---- svg helpers -------------------------------------------------------------
col_base <- "#2E7B84"; col_fr <- "#1F5AA8"; col_es <- "#4E9E52"
# same two hues as the year-effect figure for males (Figure 4)
col_m <- "#4E9E52"; col_f <- "#1F5AA8"   # males green, females blue (as in the paper)
col_text <- "#1A1A1A"; col_axis <- "#4A4A4A"
mix <- function(col, p, with = "#FFFFFF") {
  a <- grDevices::col2rgb(col); b <- grDevices::col2rgb(with)
  grDevices::rgb(t(round((1 - p) * a + p * b)), maxColorValue = 255)
}
esc <- function(s) gsub("&", "&amp;", s, fixed = TRUE)
layer <- function(id, ...) c(sprintf('<g inkscape:groupmode="layer" id="%s" inkscape:label="%s">', id, id),
                             unlist(list(...)), '</g>')

plot_year <- function(series, file) {
  W <- 560; H <- 360; x0 <- 56; x1 <- W - 22; y0 <- H - 62; y1 <- 26
  xr <- range(series[[1]]$est$year)
  yl <- c(0, max(pretty(c(0, max(vapply(series, function(z) max(z$est$hi), 0))), 5)))
  sx <- function(v) round(x0 + (v - xr[1]) / diff(xr) * (x1 - x0), 1)
  sy <- function(v) round(y0 - (v - yl[1]) / diff(yl) * (y0 - y1), 1)
  txt <- function(x, y, s, size = 11, anchor = "middle", weight = "normal",
                  colour = col_text, extra = "") sprintf(
    paste0('<text x="%s" y="%s" text-anchor="%s" font-size="%s" font-weight="%s" fill="%s"%s>%s</text>'),
    x, y, anchor, size, weight, colour, extra, esc(s))
  # the label of the lower series goes under its own band, so it can never sit
  # on top of the other curve
  tops <- vapply(series, function(z) z$est$hi[1], numeric(1))
  body <- unlist(lapply(seq_along(series), function(k) {
    z <- series[[k]]; e <- z$est
    rib <- paste0("M", paste0(sx(e$year), ",", sy(e$hi), collapse = "L"),
                  "L", paste0(sx(rev(e$year)), ",", sy(rev(e$lo)), collapse = "L"), "Z")
    ln <- paste0("M", paste0(sx(e$year), ",", sy(e$med), collapse = "L"))
    below <- length(series) > 1 && tops[k] < max(tops)
    lab <- if (is.null(z$label)) NULL else
      txt(sx(e$year[1]) + 6,
          if (below) sy(e$lo[1]) + 16 else sy(e$hi[1]) - 8,
          z$label, 13, "start", "bold", z$colour)
    layer(z$id,
      sprintf('<path d="%s" fill="%s" fill-opacity="0.72" stroke="none"/>', rib, mix(z$colour, 0.72)),
      sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="2.4" stroke-linejoin="round" stroke-linecap="round"/>', ln, z$colour),
      lab)
  }))
  xt <- pretty(xr, 5); xt <- xt[xt >= xr[1] & xt <= xr[2]]
  yt <- pretty(yl, 5); yt <- yt[yt >= yl[1] & yt <= yl[2]]
  ymid <- (sy(yl[1]) + sy(yl[2])) / 2
  svg <- c('<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" width="%spx" height="%spx" viewBox="0 0 %s %s" font-family="Arial, Helvetica, sans-serif">', W, H, W, H),
    body,
    layer("axis",
      sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>', x0, sy(yl[1]), x1, sy(yl[1]), col_axis),
      sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>', x0, sy(yl[1]), x0, sy(yl[2]), col_axis),
      sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>', sx(xt), sy(yl[1]), sx(xt), sy(yl[1]) + 5, col_axis),
      txt(sx(xt), sy(yl[1]) + 18, format(xt), 10, colour = col_axis),
      sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>', x0 - 5, sy(yt), x0, sy(yt), col_axis),
      txt(x0 - 9, sy(yt) + 4, format(yt), 10, "end", colour = col_axis),
      txt((x0 + x1) / 2, H - 16, "Year", 13, weight = "bold"),
      txt(18, ymid, "PCB concentration (mg/kg lw)", 13, weight = "bold",
          extra = sprintf(' transform="rotate(-90 %s %s)"', 18, ymid))),
    local({
      key <- if (length(series) == 1) series[[1]]$colour else col_axis
      layer("legend",
        sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="2.4"/>',
                x1 - 116, y1 + 6, x1 - 94, y1 + 6, key),
        txt(x1 - 86, y1 + 10, "Posterior median", 11, "start"),
        sprintf('<rect x="%s" y="%s" width="22" height="11" fill="%s" fill-opacity="0.72"/>',
                x1 - 116, y1 + 16, mix(key, 0.72)),
        txt(x1 - 86, y1 + 25, "95% CI", 11, "start"))
    }),
    '</svg>')
  writeLines(svg, file, useBytes = TRUE)
  cat("written:", file, "\n")
}

plot_year(list(list(id = "overall", label = NULL, est = est_all, colour = col_base)),
          "output/GA_year_effect_immatures.svg")
plot_year(list(list(id = "france", label = "France", est = est_france, colour = col_fr),
               list(id = "spain",  label = "Spain",  est = est_spain,  colour = col_es)),
          "output/GA_year_effect_immatures_country.svg")
plot_year(list(list(id = "males",   label = "Males",   est = est_males,   colour = col_m),
               list(id = "females", label = "Females", est = est_females, colour = col_f)),
          "output/GA_year_effect_immatures_sex.svg")

# ---- summary table -----------------------------------------------------------
fx <- fit$summary.fixed
ef <- exp(fx[, c("mean", "0.025quant", "0.975quant")])
tab <- data.frame(
  Parameter = c("Intercept", "Country: Spain (vs France)", "Sex: Females (vs Males)"),
  Mean   = round(fx$mean, 3), SD = round(fx$sd, 3),
  CI2.5  = round(fx$`0.025quant`, 3), CI97.5 = round(fx$`0.975quant`, 3),
  ExpMean = round(ef$mean, 3), ExpCI = sprintf("%.3f - %.3f", ef$`0.025quant`, ef$`0.975quant`),
  stringsAsFactors = FALSE)
print(tab, row.names = FALSE)

if (requireNamespace("kableExtra", quietly = TRUE)) {
  k <- kableExtra::kable_styling(
    kableExtra::kbl(tab, format = "latex", booktabs = TRUE, align = "lrrrrrr",
      col.names = c("Parameter", "Mean", "SD", "IC 2.5", "IC 97.5", "exp(Mean)", "exp 95% CI"),
      caption = "Summary of the fixed effects of the selected immature model (Country + f(Year, RW2) + Sex).\\\\label{tab:immature}"),
    latex_options = "hold_position", font_size = 9)
  writeLines(as.character(k), "output/GA_immature_summary_table.tex")
  cat("written: GA_immature_summary_table.tex\n")
}

# render the same table as a small SVG for quick viewing
rows <- nrow(tab)
cw <- c(210, 70, 60, 74, 76, 84, 120); x <- cumsum(c(20, cw))
W <- tail(x, 1) + 8; rh <- 30; H <- 54 + rows * rh + 40
cell <- function(cx, cy, s, anchor = "start", weight = "normal", size = 12, col = "#1A1A1A")
  sprintf('<text x="%s" y="%s" text-anchor="%s" font-size="%s" font-weight="%s" fill="%s" font-family="Arial">%s</text>',
          cx, cy, anchor, size, weight, col, esc(s))
hdr <- c("Parameter", "Mean", "SD", "IC 2.5", "IC 97.5", "exp(Mean)", "exp 95% CI")
svg <- c('<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
  sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%spx" height="%spx" viewBox="0 0 %s %s">', W, H, W, H),
  '<text x="20" y="26" font-size="14" font-weight="bold" font-family="Arial" fill="#1A1A1A">Immature model - fixed effects (log scale, and exp = multiplicative on PCB)</text>',
  sprintf('<rect x="8" y="40" width="%s" height="%s" rx="4" fill="#2E7B84"/>', W - 16, rh),
  unlist(lapply(seq_along(hdr), function(j)
    cell(if (j == 1) x[j] else x[j] + cw[j] - 6, 60, hdr[j],
         if (j == 1) "start" else "end", "bold", 12, "#FFFFFF"))))
for (i in seq_len(rows)) {
  yb <- 40 + i * rh
  if (i %% 2 == 0) svg <- c(svg, sprintf('<rect x="8" y="%s" width="%s" height="%s" fill="#F2F5F6"/>', yb, W - 16, rh))
  vals <- c(tab$Parameter[i], sprintf("%.3f", tab$Mean[i]), sprintf("%.3f", tab$SD[i]),
            sprintf("%.3f", tab$CI2.5[i]), sprintf("%.3f", tab$CI97.5[i]),
            sprintf("%.3f", tab$ExpMean[i]), tab$ExpCI[i])
  svg <- c(svg, unlist(lapply(seq_along(vals), function(j)
    cell(if (j == 1) x[j] else x[j] + cw[j] - 6, yb + 20, vals[j],
         if (j == 1) "start" else "end", if (j == 1) "bold" else "normal"))))
}
svg <- c(svg,
  sprintf('<text x="20" y="%s" font-size="10" font-family="Arial" fill="#4A4A4A">WAIC %.2f  |  LCPO %.2f  |  DIC %.2f  |  n = %d immatures. Year modelled as RW2 (smooth).</text>',
          H - 14, fit$waic$waic, -sum(log(fit$cpo$cpo), na.rm = TRUE), fit$dic$dic, nrow(d)),
  '</svg>')
writeLines(svg, "output/GA_immature_summary_table.svg", useBytes = TRUE)
cat("written: GA_immature_summary_table.svg\n")
