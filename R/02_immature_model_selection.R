# ------------------------------------------------------------------------------
# Model selection for IMMATURE individuals only (Maturity == 0), pooling sexes.
#
# Same procedure as the males / females selection in the report: enumerate every
# subset of the candidate terms across the four candidate sets, fit each Gamma
# GLM in INLA, rank by WAIC and LCPO. Among immatures, Maturity is constant, so
# it is dropped; Sex now varies and takes its place as a candidate covariate.
#
# The declared model space treats Year, Month and Length as smooth (RW2) effects
# (methods section), so models entering any of them linearly are excluded from
# the reported top-10 - exactly as for the males / females tables.
#
# Outputs: top10_immatures.csv        best-10 models, WAIC / LCPO / relevance
#          immature_best_model.txt    full numeric summary of the best model
#          immature_model_selection.rds
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(INLA); library(stringr) })

# ---- data --------------------------------------------------------------------
d <- read.csv("data/PCBcorrectedDB.csv", sep = ";")
d$Country_Code <- factor(ifelse(d$Country_Code == 4, "France", "Spain"),
                         levels = c("France", "Spain"))
d$Sex <- factor(ifelse(d$Sex == 1, "Males", "Females"), levels = c("Males", "Females"))

imm <- d[d$Maturity == 0, ]
cat(sprintf("immatures: %d  (Males %d, Females %d | France %d, Spain %d)\n",
            nrow(imm), sum(imm$Sex == "Males"), sum(imm$Sex == "Females"),
            sum(imm$Country_Code == "France"), sum(imm$Country_Code == "Spain")))

# ---- candidate term sets (Sex replaces Maturity) -----------------------------
sets <- list(
  c("Country_Code", "Year", "Length", "Sex"),
  c("Country_Code", "f(Year, model = 'rw2')", "f(Month, model = 'rw2', cyclic = TRUE)", "Length", "Sex"),
  c("Country_Code", "f(Year, model = 'rw2')", "f(Month, model = 'rw2', cyclic = TRUE)", "f(Length, model = 'rw2')", "Sex"),
  c("Country_Code", "f(Year, model = 'iid')", "f(Month, model = 'rw2', cyclic = TRUE)", "f(Length, model = 'rw2')", "Sex")
)

# all subsets of one term set, as "resp ~ 1 + ..." formulas
all_formulas <- function(v) {
  out <- "resp ~ 1"
  for (m in seq_along(v))
    out <- c(out, apply(combn(v, m), 2, function(z)
      paste("resp ~ 1", paste(z, collapse = " + "), sep = " + ")))
  out
}
forms <- unique(unlist(lapply(sets, all_formulas)))
cat("unique candidate models:", length(forms), "\n")

# ---- fit every model ---------------------------------------------------------
imm$resp <- imm$ICES73
fit_one <- function(f)
  try(inla(as.formula(f), family = "gamma", data = imm,
           control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
           control.inla = list(strategy = "simplified.laplace")), silent = TRUE)

res <- vector("list", length(forms))
WAIC <- LCPO <- DIC <- rep(NA_real_, length(forms))
for (i in seq_along(forms)) {
  m <- fit_one(forms[i])
  if (!inherits(m, "try-error")) {
    res[[i]] <- m
    WAIC[i] <- m$waic$waic
    DIC[i]  <- m$dic$dic
    if (!is.null(m$cpo$cpo)) LCPO[i] <- -sum(log(m$cpo$cpo), na.rm = TRUE)
  }
  cat(sprintf("\r  fitted %d / %d", i, length(forms)))
}
cat("\n")

# ---- relevance of covariates: 95% CI excluding 0 -----------------------------
# fixed effects (Country, Sex, linear Year/Length): CI of the coefficient
# smooth effects (rw2 Year/Month/Length): the effect's 95% band leaves 0 for
#   at least one level (i.e. the smooth is not a flat zero line)
relevance <- function(m) {
  bad <- character(0)
  fx <- m$summary.fixed
  fx <- fx[rownames(fx) != "(Intercept)", , drop = FALSE]
  if (nrow(fx))
    for (r in rownames(fx))
      if (fx[r, "0.025quant"] <= 0 && fx[r, "0.975quant"] >= 0)
        bad <- c(bad, r)
  for (nm in names(m$summary.random)) {
    s <- m$summary.random[[nm]]
    if (!any(s[["0.025quant"]] > 0 | s[["0.975quant"]] < 0))
      bad <- c(bad, nm)                       # band never departs from 0
  }
  list(all_relevant = length(bad) == 0, not_relevant = bad)
}

# ---- covariate label ---------------------------------------------------------
linear_part <- function(x) str_remove_all(x, "f\\s*\\([^)]*\\)")
norm_covariates <- function(x) {
  lin <- linear_part(x); parts <- c()
  if (str_detect(x, "\\bCountry(?:_Code)?\\b"))    parts <- c(parts, "Country")
  if (str_detect(x, "f\\s*\\(\\s*Year[^)]*iid"))   parts <- c(parts, "Year (iid)")
  if (str_detect(x, "f\\s*\\(\\s*Year[^)]*rw2"))   parts <- c(parts, "Year (rw2)")
  if (str_detect(lin, "\\bYear\\b"))               parts <- c(parts, "Year (linear)")
  if (str_detect(x, "f\\s*\\(\\s*Month[^)]*rw2") &&
      str_detect(x, "f\\s*\\(\\s*Month[^)]*cyclic")) parts <- c(parts, "Month (rw2, cc)")
  if (str_detect(x, "f\\s*\\(\\s*Length[^)]*rw2"))  parts <- c(parts, "Length (rw2)")
  if (str_detect(lin, "\\bLength\\b"))              parts <- c(parts, "Length (linear)")
  if (str_detect(x, "\\bSex\\b"))                   parts <- c(parts, "Sex")
  if (length(parts) == 0) "Intercept only" else paste(parts, collapse = ", ")
}

# ---- assemble, restrict to the declared (smooth) model space, rank -----------
ok <- !is.na(WAIC)
tab <- data.frame(idx = which(ok), formula = forms[ok],
                  WAIC = WAIC[ok], LCPO = LCPO[ok], DIC = DIC[ok],
                  stringsAsFactors = FALSE)

linear_smooth <- str_detect(linear_part(tab$formula), "\\bYear\\b|\\bMonth\\b|\\bLength\\b")
cat(sprintf("fitted %d models; %d use a linear Year/Month/Length -> excluded, %d left\n",
            nrow(tab), sum(linear_smooth), sum(!linear_smooth)))
tab <- tab[!linear_smooth, ]
tab <- tab[order(tab$WAIC), ]

rel <- lapply(tab$idx, function(i) relevance(res[[i]]))
tab$Covariates    <- vapply(tab$formula, norm_covariates, character(1))
tab$AllRelevant   <- vapply(rel, `[[`, logical(1), "all_relevant")
tab$NotRelevant   <- vapply(rel, function(z) paste(z$not_relevant, collapse = "; "),
                            character(1))

top <- head(tab, 10)
top$Model <- paste0("M", seq_len(nrow(top)))

out <- top[, c("Model", "Covariates", "WAIC", "LCPO", "AllRelevant", "NotRelevant")]
out$WAIC <- round(out$WAIC, 2); out$LCPO <- round(out$LCPO, 2)

cat("\n================ TOP 10 IMMATURE MODELS (by WAIC) ================\n")
print(out, row.names = FALSE, right = FALSE)

write.csv(out, "output/top10_immatures.csv", row.names = FALSE)
saveRDS(tab, "output/immature_model_selection.rds")

# ---- best model = best WAIC among the fully-relevant models ------------------
best_i <- tab$idx[which(tab$AllRelevant)[1]]
if (is.na(best_i)) best_i <- tab$idx[1]
best <- res[[best_i]]

cat("\n================ BEST MODEL (best WAIC with all covariates relevant) ====\n")
cat("formula:", sub("resp", "ICES73", forms[best_i]), "\n")
cat(sprintf("WAIC %.2f | LCPO %.2f | DIC %.2f\n\n",
            best$waic$waic, -sum(log(best$cpo$cpo), na.rm = TRUE), best$dic$dic))

sink("immature_best_model.txt")
cat("Best immature model\n")
cat("formula:", sub("resp", "ICES73", forms[best_i]), "\n")
cat(sprintf("WAIC %.2f | LCPO %.2f | DIC %.2f\n\n",
            best$waic$waic, -sum(log(best$cpo$cpo), na.rm = TRUE), best$dic$dic))
cat("--- Fixed effects (log scale) ---\n")
print(round(best$summary.fixed, 4))
cat("\n--- Fixed effects (multiplicative, exp) ---\n")
ef <- exp(best$summary.fixed[, c("mean", "0.025quant", "0.975quant")])
print(round(ef, 3))
cat("\n--- Hyperparameters ---\n")
print(round(best$summary.hyperpar, 4))
for (nm in names(best$summary.random)) {
  cat(sprintf("\n--- Smooth effect: %s ---\n", nm))
  print(round(best$summary.random[[nm]][, c("ID","mean","0.025quant","0.975quant")], 4))
}
sink()

cat("--- Fixed effects (log scale) ---\n")
print(round(best$summary.fixed, 4))
cat("\n--- Fixed effects, exp(coef) = multiplicative effect on PCB ---\n")
print(round(exp(best$summary.fixed[, c("mean", "0.025quant", "0.975quant")]), 3))
cat("\n--- Hyperparameters ---\n")
print(round(best$summary.hyperpar, 4))

cat("\nwritten: top10_immatures.csv, immature_best_model.txt, immature_model_selection.rds\n")
