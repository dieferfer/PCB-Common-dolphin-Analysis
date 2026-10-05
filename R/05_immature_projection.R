# ------------------------------------------------------------------------------
# Projections from the immature model, following exactly the procedure used for
# the male projections (report, chunk "Predictiva a posteriori"):
#
#   horizon 2020-2075 | future rows with ICES73 = NA | refit with
#   control.predictor(compute = TRUE, link = 1) | 2000 inla.posterior.sample()
#   draws | first year with P(mu < threshold) >= 0.8
#
# Model: ICES73 ~ Country_Code + Sex + f(Year, rw2), immatures only (n = 78).
#
# Outputs: GA_projection_immatures_sex.svg   males vs females, one panel
#          immature_projection_crossings.csv threshold-crossing years
#          immature_projection_curves.csv    yearly medians and 95% CI
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(INLA) })
set.seed(1)

years_future <- 2020:2075
n_samp       <- 2000
p_level      <- 0.8
thresholds   <- c(17, 41)

# ---- data --------------------------------------------------------------------
d <- read.csv("data/PCBcorrectedDB.csv", sep = ";")
d$Country_Code <- factor(ifelse(d$Country_Code == 4, "France", "Spain"),
                         levels = c("France", "Spain"))
d$Sex <- factor(ifelse(d$Sex == 1, "Males", "Females"), levels = c("Males", "Females"))
imm <- d[d$Maturity == 0, c("ICES73", "Year", "Country_Code", "Sex")]

cat(sprintf("immatures: %d (France %d, Spain %d | Males %d, Females %d)\n",
            nrow(imm), sum(imm$Country_Code == "France"), sum(imm$Country_Code == "Spain"),
            sum(imm$Sex == "Males"), sum(imm$Sex == "Females")))

# ---- augment with the future scenarios --------------------------------------
new_future <- expand.grid(Year = years_future,
                          Country_Code = levels(imm$Country_Code),
                          Sex = levels(imm$Sex))
new_future$ICES73 <- NA_real_
data_aug <- rbind(imm, new_future[, names(imm)])
data_aug$is_pred <- is.na(data_aug$ICES73)

fit <- inla(ICES73 ~ Country_Code + Sex + f(Year, model = "rw2"),
            family = "gamma", data = data_aug,
            control.predictor = list(compute = TRUE, link = 1),
            control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE))

cat(sprintf("refitted model: WAIC %.2f | Spain %.3f | Females %.3f\n",
            fit$waic$waic, fit$summary.fixed["Country_CodeSpain", "mean"],
            fit$summary.fixed["SexFemales", "mean"]))

# ---- posterior draws of mu on the prediction rows ----------------------------
samps     <- inla.posterior.sample(n_samp, fit)
pred_rows <- which(data_aug$is_pred)
eta <- sapply(samps, function(s) as.numeric(s$latent[paste0("Predictor:", pred_rows), 1]))
if (nrow(eta) != length(pred_rows)) eta <- t(eta)
mu <- exp(eta)                                    # scenarios x draws

pred_df <- data_aug[pred_rows, c("Year", "Country_Code", "Sex")]

# ---- sex-level curves ---------------------------------------------------------
# Both sexes are averaged over the SAME country composition of the immature
# sample, so the two curves differ only by the sex coefficient.
w <- c(France = mean(imm$Country_Code == "France"),
       Spain  = mean(imm$Country_Code == "Spain"))
cat(sprintf("country weights: France %.3f, Spain %.3f\n", w["France"], w["Spain"]))

mu_sex <- function(sex) {
  iF <- which(pred_df$Sex == sex & pred_df$Country_Code == "France")
  iS <- which(pred_df$Sex == sex & pred_df$Country_Code == "Spain")
  oF <- order(pred_df$Year[iF]); oS <- order(pred_df$Year[iS])
  w["France"] * mu[iF[oF], , drop = FALSE] + w["Spain"] * mu[iS[oS], , drop = FALSE]
}
mu_m <- mu_sex("Males"); mu_f <- mu_sex("Females")

summarise_mu <- function(M) data.frame(
  year = sort(unique(pred_df$Year)),
  med  = apply(M, 1, median),
  lo   = apply(M, 1, quantile, 0.025),
  hi   = apply(M, 1, quantile, 0.975))

est_m <- summarise_mu(mu_m); est_f <- summarise_mu(mu_f)

curves <- rbind(cbind(Sex = "Males", est_m), cbind(Sex = "Females", est_f))
write.csv(curves, "output/immature_projection_curves.csv", row.names = FALSE)

# full Year x Country x Sex grid, for the manuscript figure
grid <- do.call(rbind, lapply(
  split(seq_len(nrow(pred_df)), list(pred_df$Year, pred_df$Country_Code, pred_df$Sex),
        drop = TRUE),
  function(i) {
    v <- as.numeric(mu[i, ])
    data.frame(Year = pred_df$Year[i[1]], Country_Code = pred_df$Country_Code[i[1]],
               Sex = pred_df$Sex[i[1]], mu_median = median(v),
               mu_q025 = unname(quantile(v, 0.025)), mu_q975 = unname(quantile(v, 0.975)))
  }))
write.csv(grid, "output/immature_projection_grid.csv", row.names = FALSE)

cat("\n--- projected concentration (mg/kg lw), pooled over countries ---\n")
show <- c(2020, 2030, 2040, 2050, 2060, 2075)
for (s in c("Males", "Females")) {
  e <- if (s == "Males") est_m else est_f
  cat(sprintf("%-8s ", s))
  cat(paste(sprintf("%d: %.1f (%.1f-%.1f)", show,
                    e$med[match(show, e$year)], e$lo[match(show, e$year)],
                    e$hi[match(show, e$year)]), collapse = "  "), "\n")
}

# ---- first year with P(mu < threshold) >= 0.8 --------------------------------
first_year <- function(M, yrs, thr) {
  p <- rowMeans(M < thr)
  ok <- which(p >= p_level)
  if (!length(ok)) NA_integer_ else yrs[min(ok)]
}
yrs <- sort(unique(pred_df$Year))

cross <- data.frame()
for (thr in thresholds) {
  # pooled over countries, by sex
  cross <- rbind(cross, data.frame(
    Threshold = thr, Country = "Both (pooled)",
    Males   = first_year(mu_m, yrs, thr),
    Females = first_year(mu_f, yrs, thr)))
  # by country, the direct analogue of the male projections
  for (cty in levels(imm$Country_Code)) {
    gm <- which(pred_df$Sex == "Males" & pred_df$Country_Code == cty)
    gf <- which(pred_df$Sex == "Females" & pred_df$Country_Code == cty)
    gm <- gm[order(pred_df$Year[gm])]; gf <- gf[order(pred_df$Year[gf])]
    cross <- rbind(cross, data.frame(
      Threshold = thr, Country = cty,
      Males   = first_year(mu[gm, , drop = FALSE], yrs, thr),
      Females = first_year(mu[gf, , drop = FALSE], yrs, thr)))
  }
}
cat("\n--- first year with P(mu < threshold) >= 0.8  (NA = not reached by 2075) ---\n")
print(cross, row.names = FALSE)
write.csv(cross, "output/immature_projection_crossings.csv", row.names = FALSE)

# ---- P(mu < threshold) by year, per country and sex --------------------------
group_mu <- function(cty, sex) {
  i <- which(pred_df$Sex == sex & pred_df$Country_Code == cty)
  mu[i[order(pred_df$Year[i])], , drop = FALSE]
}
groups <- list(
  list(id = "france_males",   cty = "France", sex = "Males",   lab = "France, males"),
  list(id = "france_females", cty = "France", sex = "Females", lab = "France, females"),
  list(id = "spain_males",    cty = "Spain",  sex = "Males",   lab = "Spain, males"),
  list(id = "spain_females",  cty = "Spain",  sex = "Females", lab = "Spain, females"))

prob_curves <- lapply(thresholds, function(thr)
  lapply(groups, function(g) {
    p <- rowMeans(group_mu(g$cty, g$sex) < thr)
    ok <- which(p >= p_level)
    list(id = g$id, lab = g$lab, cty = g$cty, sex = g$sex, year = yrs, p = p,
         cross = if (length(ok)) yrs[min(ok)] else NA_integer_)
  }))
names(prob_curves) <- paste0("thr", thresholds)

prob_out <- do.call(rbind, lapply(seq_along(thresholds), function(k)
  do.call(rbind, lapply(prob_curves[[k]], function(z)
    data.frame(Threshold = thresholds[k], Country = z$cty, Sex = z$sex,
               Year = z$year, P_below = round(z$p, 4))))))
write.csv(prob_out, "output/immature_probability_below.csv", row.names = FALSE)

# ---- figure ------------------------------------------------------------------
col_m <- "#B4632A"; col_f <- "#8145B5"
col_text <- "#1A1A1A"; col_axis <- "#4A4A4A"; col_thr <- "#C0392B"
mix <- function(col, p, with = "#FFFFFF") {
  a <- grDevices::col2rgb(col); b <- grDevices::col2rgb(with)
  grDevices::rgb(t(round((1 - p) * a + p * b)), maxColorValue = 255)
}
esc <- function(s) { s <- gsub("&", "&amp;", s, fixed = TRUE)
                     gsub("<", "&lt;", s, fixed = TRUE) }
layer <- function(id, ...) c(sprintf('<g inkscape:groupmode="layer" id="%s" inkscape:label="%s">', id, id),
                             unlist(list(...)), '</g>')

W <- 640; H <- 400; x0 <- 62; x1 <- W - 24; y0 <- H - 62; y1 <- 30
xr <- range(yrs)
# log10 y axis, as in the male projection figure: extrapolating the RW2 to 2075
# widens the 95% band over orders of magnitude, which a linear axis cannot show
ylo <- 0.1; yhi <- 10^ceiling(log10(max(est_m$hi)))
sx <- function(v) round(x0 + (v - xr[1]) / diff(xr) * (x1 - x0), 1)
sy <- function(v) {
  v <- pmin(pmax(v, ylo), yhi)
  round(y0 - (log10(v) - log10(ylo)) / (log10(yhi) - log10(ylo)) * (y0 - y1), 1)
}
txt <- function(x, y, s, size = 11, anchor = "middle", weight = "normal",
                colour = col_text, extra = "") sprintf(
  '<text x="%s" y="%s" text-anchor="%s" font-size="%s" font-weight="%s" fill="%s"%s>%s</text>',
  x, y, anchor, size, weight, colour, extra, esc(s))

# ribbons go underneath both lines, otherwise the second band drawn hides the
# first median line
ribbon_of <- function(e, col)
  sprintf('<path d="%s" fill="%s" fill-opacity="0.45" stroke="none"/>',
          paste0("M", paste0(sx(e$year), ",", sy(e$hi), collapse = "L"),
                 "L", paste0(sx(rev(e$year)), ",", sy(rev(e$lo)), collapse = "L"), "Z"),
          mix(col, 0.62))

line_of <- function(e, col, label)
  c(sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="2.6" stroke-linejoin="round"/>',
            paste0("M", paste0(sx(e$year), ",", sy(e$med), collapse = "L")), col),
    txt(sx(2021), sy(e$med[e$year == 2021]) - 9, label, 13, "start", "bold", col))

xt <- seq(2020, 2075, by = 10)
yt <- 10^seq(log10(ylo), log10(yhi))          # one tick per decade
ymid <- (sy(ylo) + sy(yhi)) / 2

svg <- c('<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
  sprintf('<svg xmlns="http://www.w3.org/2000/svg" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" width="%spx" height="%spx" viewBox="0 0 %s %s" font-family="Arial, Helvetica, sans-serif">', W, H, W, H),
  layer("thresholds", unlist(lapply(thresholds, function(t) c(
    sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1.2" stroke-dasharray="6,4"/>',
            x0, sy(t), x1, sy(t), col_thr),
    txt(x1 - 2, sy(t) - 5, sprintf("%d mg/kg lw", t), 10, "end", colour = col_thr))))),
  layer("bands", ribbon_of(est_m, col_m), ribbon_of(est_f, col_f)),
  layer("males",   line_of(est_m, col_m, "Immature males")),
  layer("females", line_of(est_f, col_f, "Immature females")),
  layer("axis",
    sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>', x0, sy(ylo), x1, sy(ylo), col_axis),
    sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>', x0, sy(ylo), x0, sy(yhi), col_axis),
    sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>', sx(xt), sy(ylo), sx(xt), sy(ylo) + 5, col_axis),
    txt(sx(xt), sy(ylo) + 18, as.character(xt), 10, colour = col_axis),
    sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>', x0 - 5, sy(yt), x0, sy(yt), col_axis),
    txt(x0 - 9, sy(yt) + 4, formatC(yt, format = "fg", drop0trailing = TRUE),
        10, "end", colour = col_axis),
    txt((x0 + x1) / 2, H - 16, "Year", 13, weight = "bold"),
    txt(18, ymid, "Predicted PCB concentration (mg/kg lw, log scale)", 12, weight = "bold",
        extra = sprintf(' transform="rotate(-90 %s %s)"', 18, ymid))),
  layer("title",
    txt(x0, 22, "Projected PCB concentrations in immature common dolphins", 13,
        "start", "bold")),
  '</svg>')
writeLines(svg, "output/GA_projection_immatures_sex.svg", useBytes = TRUE)

# ==============================================================================
# P(mu < threshold) vs year: one panel per threshold, four curves each
# colour = country (as elsewhere in the paper), dashes = females
# ==============================================================================
col_fr <- "#1F5AA8"; col_es <- "#4E9E52"
PW <- 900; PH <- 430
pan_w <- 366; pan_h <- 268
pan_x <- c(74, 500); pan_top <- 74
py0 <- pan_top + pan_h                      # p = 0
py1 <- pan_top                              # p = 1

px <- function(k, v) round(pan_x[k] + (v - xr[1]) / diff(xr) * pan_w, 1)
py <- function(v) round(py0 - v * (py0 - py1), 1)

style_of <- function(z) list(
  col = if (z$cty == "France") col_fr else col_es,
  dash = if (z$sex == "Females") ' stroke-dasharray="7,4"' else "")

panel <- function(k, thr) {
  cur <- prob_curves[[k]]
  out <- c(
    # frame and grid
    sprintf('<rect x="%s" y="%s" width="%s" height="%s" fill="#FFFFFF" stroke="%s" stroke-width="1"/>',
            pan_x[k], py1, pan_w, pan_h, "#9AA3AA"),
    unlist(lapply(seq(0, 1, by = 0.25), function(v)
      sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="#E8EBED" stroke-width="1"/>',
              pan_x[k], py(v), pan_x[k] + pan_w, py(v)))),
    # the 0.8 decision level
    sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1.4" stroke-dasharray="6,4"/>',
            pan_x[k], py(p_level), pan_x[k] + pan_w, py(p_level), col_thr),
    txt(pan_x[k] + pan_w - 4, py(p_level) - 6, sprintf("P = %.1f", p_level), 10, "end",
        colour = col_thr),
    txt(pan_x[k] + pan_w / 2, py1 - 12, sprintf("Threshold %d mg/kg lw", thr), 13,
        weight = "bold"))
  # curves, then the crossing markers on top
  for (z in cur) {
    s <- style_of(z)
    out <- c(out, sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="2.2"%s stroke-linejoin="round"/>',
      paste0("M", paste0(px(k, z$year), ",", py(z$p), collapse = "L")), s$col, s$dash))
  }
  for (z in cur) if (!is.na(z$cross)) {
    s <- style_of(z)
    out <- c(out,
      sprintf('<circle cx="%s" cy="%s" r="4" fill="%s" stroke="#FFFFFF" stroke-width="1.2"/>',
              px(k, z$cross), py(p_level), s$col))
  }
  # axes
  xtk <- seq(2020, 2075, by = 10)
  out <- c(out,
    sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>',
            px(k, xtk), py0, px(k, xtk), py0 + 5, col_axis),
    txt(px(k, xtk), py0 + 18, as.character(xtk), 10, colour = col_axis),
    txt(pan_x[k] + pan_w / 2, py0 + 40, "Year", 12, weight = "bold"))
  if (k == 1) {
    ytk <- seq(0, 1, by = 0.25)
    out <- c(out,
      sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>',
              pan_x[1] - 5, py(ytk), pan_x[1], py(ytk), col_axis),
      txt(pan_x[1] - 9, py(ytk) + 4, sprintf("%.2f", ytk), 10, "end", colour = col_axis),
      txt(22, (py0 + py1) / 2, "P(mean PCB < threshold)", 12, weight = "bold",
          extra = sprintf(' transform="rotate(-90 %s %s)"', 22, (py0 + py1) / 2)))
  }
  layer(sprintf("panel-%d", thr), out)
}

leg <- function(x, y, col, dash, lab) c(
  sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="2.2"%s/>',
          x, y, x + 26, y, col, dash),
  txt(x + 32, y + 4, lab, 11, "start"))

svg2 <- c('<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
  sprintf('<svg xmlns="http://www.w3.org/2000/svg" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" width="%spx" height="%spx" viewBox="0 0 %s %s" font-family="Arial, Helvetica, sans-serif">', PW, PH, PW, PH),
  sprintf('<rect width="%s" height="%s" fill="#FFFFFF"/>', PW, PH),
  layer("title",
    txt(74, 26, "Probability that immature dolphins fall below each toxicity threshold", 14,
        "start", "bold"),
    txt(74, 44, "Posterior predictive probability by year; the marker shows the first year reaching P = 0.8",
        11, "start", colour = "#4A4A4A")),
  panel(1, thresholds[1]),
  panel(2, thresholds[2]),
  layer("legend",
    leg(pan_x[1],       PH - 14, col_fr, "",                     "France, males"),
    leg(pan_x[1] + 150, PH - 14, col_fr, ' stroke-dasharray="7,4"', "France, females"),
    leg(pan_x[2],       PH - 14, col_es, "",                     "Spain, males"),
    leg(pan_x[2] + 150, PH - 14, col_es, ' stroke-dasharray="7,4"', "Spain, females")),
  '</svg>')
writeLines(svg2, "output/GA_probability_below_threshold.svg", useBytes = TRUE)

cat("\nwritten: GA_projection_immatures_sex.svg, GA_probability_below_threshold.svg,\n")
cat("         immature_projection_crossings.csv, immature_projection_curves.csv,\n")
cat("         immature_probability_below.csv\n")
