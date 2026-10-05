# ------------------------------------------------------------------------------
# Graphical abstract - non-linear effect of Year in the male model
#
# Refits M10, the best-WAIC/LCPO male model of the paper
#     Y_i ~ Gamma(mu_i, phi)
#     log(mu_i) = b0 + b1*Country_Spain + b2*Maturity_Matures + f(Year), rw2
# and draws f(Year) back on the concentration scale.
#
# The ribbon comes from the JOINT posterior (inla.posterior.sample), not from
# adding up marginal intervals: exp(b0 + f(Year)) is non-linear, so summing the
# separate quantiles would give a wrong interval.
#
# Output: GA_year_effect_males.svg
# ------------------------------------------------------------------------------

suppressPackageStartupMessages(library(INLA))

out_svg <- "output/GA_year_effect_males.svg"
n_draws <- 4000
set.seed(1)

# ---- data ---------------------------------------------------------------------
d <- read.csv("data/PCBcorrectedDB.csv", sep = ";")
d <- d[d$Sex == 1, ]                                   # males
d$Country_Code <- factor(ifelse(d$Country_Code == 4, "France", "Spain"),
                         levels = c("France", "Spain"))
d$Maturity <- factor(ifelse(d$Maturity == 0, "Immatures", "Matures"),
                     levels = c("Immatures", "Matures"))
cat("males:", nrow(d), " years:", min(d$Year), "-", max(d$Year), "\n")

# ---- fit ----------------------------------------------------------------------
fit <- inla(ICES73 ~ Country_Code + Maturity + f(Year, model = "rw2"),
            family = "gamma", data = d,
            control.compute = list(config = TRUE, dic = TRUE, waic = TRUE))

# ---- partial effect on the response scale -------------------------------------
# The model is log-additive, so the shape of f(Year) is fixed; where the curve
# sits depends on the other covariates. Maturity is averaged over the observed
# composition of the 90 males in every series, so the only thing separating the
# two country curves is the country coefficient - a different intercept, and
# hence a constant ratio exp(beta_Spain) between them.
s   <- inla.posterior.sample(n_draws, fit)
rn  <- rownames(s[[1]]$latent)
i_b0 <- grep("^\\(Intercept\\)", rn)
i_sp <- grep("^Country_CodeSpain", rn)
i_mt <- grep("^MaturityMatures", rn)
i_yr <- grep("^Year:", rn)

lat  <- vapply(s, function(x) x$latent[, 1], numeric(length(rn)))
b0   <- lat[i_b0, ]; b_sp <- lat[i_sp, ]; b_mt <- lat[i_mt, ]
fyr  <- lat[i_yr, , drop = FALSE]

yrs     <- fit$summary.random$Year$ID
mat_off <- colMeans(exp(outer(as.numeric(d$Maturity == "Matures"), b_mt)))
base    <- exp(sweep(fyr, 2, b0, "+"))            # exp(b0 + f(Year)), years x draws

summarise <- function(pred) data.frame(
  year = yrs,
  med  = apply(pred, 1, median),
  lo   = apply(pred, 1, quantile, 0.025),
  hi   = apply(pred, 1, quantile, 0.975))

# pooled over the whole male sample, and one curve per country
off_pool <- colMeans(exp(outer(as.numeric(d$Country_Code == "Spain"), b_sp) +
                         outer(as.numeric(d$Maturity == "Matures"), b_mt)))

est        <- summarise(base * rep(off_pool, each = length(yrs)))
est_france <- summarise(base * rep(mat_off, each = length(yrs)))
est_spain  <- summarise(base * rep(mat_off * exp(b_sp), each = length(yrs)))

cat("\npooled over the male sample:\n"); print(round(est, 1), row.names = FALSE)
cat("\nFrance:\n"); print(round(est_france, 1), row.names = FALSE)
cat("\nSpain:\n");  print(round(est_spain, 1), row.names = FALSE)
cat(sprintf("\nSpain / France ratio: %.2f (95%% CI %.2f - %.2f)\n",
            median(exp(b_sp)), quantile(exp(b_sp), 0.025), quantile(exp(b_sp), 0.975)))

# ---- look ---------------------------------------------------------------------
col_base   <- "#2E7B84"
col_france <- "#1F5AA8"
col_spain  <- "#4E9E52"
col_text   <- "#1A1A1A"
col_axis   <- "#4A4A4A"

mix <- function(col, p, with = "#FFFFFF") {
  a <- grDevices::col2rgb(col); b <- grDevices::col2rgb(with)
  grDevices::rgb(t(round((1 - p) * a + p * b)), maxColorValue = 255)
}

esc <- function(s) gsub("&", "&amp;", s, fixed = TRUE)
layer <- function(id, ...) {
  c(sprintf('<g inkscape:groupmode="layer" id="%s" inkscape:label="%s">', id, id),
    unlist(list(...)), '</g>')
}

# series: list of list(id, label, est, colour); label = NULL draws no curve label
make_plot <- function(series, file) {
  W <- 560; H <- 360
  x0 <- 56; x1 <- W - 22; y0 <- H - 62; y1 <- 26

  xr <- range(series[[1]]$est$year)
  yr_lim <- c(0, max(pretty(c(0, max(vapply(series, function(s) max(s$est$hi), 0))), 5)))

  sx <- function(v) round(x0 + (v - xr[1]) / diff(xr) * (x1 - x0), 1)
  sy <- function(v) round(y0 - (v - yr_lim[1]) / diff(yr_lim) * (y0 - y1), 1)

  txt <- function(x, y, s, size = 11, anchor = "middle", weight = "normal",
                  colour = col_text, extra = "") {
    sprintf(paste0('<text x="%s" y="%s" text-anchor="%s" font-size="%s" ',
                   'font-weight="%s" fill="%s"%s>%s</text>'),
            x, y, anchor, size, weight, colour, extra, esc(s))
  }

  body <- unlist(lapply(series, function(s) {
    e <- s$est
    ribbon <- paste0("M", paste0(sx(e$year), ",", sy(e$hi), collapse = "L"),
                     "L", paste0(sx(rev(e$year)), ",", sy(rev(e$lo)), collapse = "L"), "Z")
    line <- paste0("M", paste0(sx(e$year), ",", sy(e$med), collapse = "L"))
    lab <- if (is.null(s$label)) NULL else
      txt(sx(e$year[1]) + 6, sy(e$hi[1]) - 8, s$label, 13, anchor = "start",
          weight = "bold", colour = s$colour)
    layer(s$id,
      sprintf('<path d="%s" fill="%s" fill-opacity="0.75" stroke="none"/>',
              ribbon, mix(s$colour, 0.72)),
      sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="2.4" stroke-linejoin="round" stroke-linecap="round"/>',
              line, s$colour),
      lab)
  }))

  xticks <- pretty(xr, 5); xticks <- xticks[xticks >= xr[1] & xticks <= xr[2]]
  yticks <- pretty(yr_lim, 5); yticks <- yticks[yticks >= yr_lim[1] & yticks <= yr_lim[2]]
  ymid <- (sy(yr_lim[1]) + sy(yr_lim[2])) / 2

  # the legend keys are neutral: the colours already say which country is which
  key <- if (length(series) == 1) series[[1]]$colour else col_axis

  svg <- c(
    '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
    sprintf(paste0('<svg xmlns="http://www.w3.org/2000/svg" ',
                   'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" ',
                   'width="%spx" height="%spx" viewBox="0 0 %s %s" ',
                   'font-family="Arial, Helvetica, sans-serif">'), W, H, W, H),
    body,
    layer("axis",
      sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>',
              x0, sy(yr_lim[1]), x1, sy(yr_lim[1]), col_axis),
      sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>',
              x0, sy(yr_lim[1]), x0, sy(yr_lim[2]), col_axis),
      sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>',
              sx(xticks), sy(yr_lim[1]), sx(xticks), sy(yr_lim[1]) + 5, col_axis),
      txt(sx(xticks), sy(yr_lim[1]) + 18, format(xticks), 10, colour = col_axis),
      sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>',
              x0 - 5, sy(yticks), x0, sy(yticks), col_axis),
      txt(x0 - 9, sy(yticks) + 4, format(yticks), 10, anchor = "end", colour = col_axis),
      txt((x0 + x1) / 2, H - 16, "Year", 13, weight = "bold"),
      txt(18, ymid, "PCB concentration (mg/kg lw)", 13, weight = "bold",
          extra = sprintf(' transform="rotate(-90 %s %s)"', 18, ymid))),
    layer("legend",
      sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="2.4"/>',
              x1 - 116, y1 + 6, x1 - 94, y1 + 6, key),
      txt(x1 - 86, y1 + 10, "Posterior median", 11, anchor = "start"),
      sprintf('<rect x="%s" y="%s" width="22" height="11" fill="%s" fill-opacity="0.75"/>',
              x1 - 116, y1 + 16, mix(key, 0.72)),
      txt(x1 - 86, y1 + 25, "95% CI", 11, anchor = "start")),
    '</svg>')

  writeLines(svg, file, useBytes = TRUE)
  cat(sprintf("written: %s (%s x %s px)\n", file, W, H))
}

make_plot(list(list(id = "trend", label = NULL, est = est, colour = col_base)),
          out_svg)

make_plot(list(
  list(id = "france", label = "France", est = est_france, colour = col_france),
  list(id = "spain",  label = "Spain",  est = est_spain,  colour = col_spain)),
  "output/GA_year_effect_country.svg")

cat(sprintf("WAIC %.1f   DIC %.1f\n", fit$waic$waic, fit$dic$dic))
