# ------------------------------------------------------------------------------
# Top-10 model selection tables, restricted to the model space the paper's
# methods section actually declares. Two exclusion criteria:
#
#   1. Length AND Maturity in the same linear predictor. Maturity is derived
#      from Length, so they must not enter together.
#   2. Year, Month or Length entering as a plain linear effect. The methods
#      state these three carry RW2 priors, i.e. they are smooth effects.
#
# Dropping rows and re-ranking is equivalent to having run the selection over
# the restricted space: WAIC/LCPO are per-model and do not depend on which
# other models were fitted.
#
# Outputs: top10_males_filtered.csv / .tex
#          top10_females_filtered.csv / .tex
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(dplyr); library(stringr) })

# whatever sits outside the f(...) terms is entering the model linearly
linear_part <- function(x) str_remove_all(x, "f\\s*\\([^)]*\\)")

norm_covariates <- function(x) {
  lin <- linear_part(x)
  parts <- c()
  if (str_detect(x, "\\bCountry(?:_Code)?\\b"))    parts <- c(parts, "Country")
  if (str_detect(x, "f\\s*\\(\\s*Year[^)]*iid"))   parts <- c(parts, "Year (iid)")
  if (str_detect(x, "f\\s*\\(\\s*Year[^)]*rw2"))   parts <- c(parts, "Year (rw2)")
  if (str_detect(lin, "\\bYear\\b"))               parts <- c(parts, "Year (linear)")
  if (str_detect(x, "f\\s*\\(\\s*Month[^)]*rw2") &&
      str_detect(x, "f\\s*\\(\\s*Month[^)]*cyclic")) parts <- c(parts, "Month (rw2, cc)")
  if (str_detect(x, "f\\s*\\(\\s*Length[^)]*rw2"))  parts <- c(parts, "Length (rw2)")
  if (str_detect(lin, "\\bLength\\b"))              parts <- c(parts, "Length (linear)")
  if (str_detect(x, "\\bMaturity\\b"))              parts <- c(parts, "Maturity")
  paste(parts, collapse = ", ")
}

drop_len_mat <- function(x) str_detect(x, "\\bLength\\b") & str_detect(x, "\\bMaturity\\b")
drop_linear  <- function(x) str_detect(linear_part(x), "\\bYear\\b|\\bMonth\\b|\\bLength\\b")

# The male table describes models in `Model` ("Int + Country_Code + Year + ...");
# the female one carries the INLA formula in `Formula` (`Model` is just a label).
build <- function(rds, text_col, label, out_stem) {
  x <- readRDS(rds)
  n_all <- nrow(x)

  ok <- !is.na(x[[text_col]]) & !is.na(x$WAIC) & !is.na(x$LCPO)
  if (any(!ok))
    cat(sprintf("\n== %s: %d row(s) with no WAIC/LCPO (failed fits) removed\n",
                label, sum(!ok)))
  x   <- x[ok, , drop = FALSE]
  txt <- x[[text_col]]

  d1 <- drop_len_mat(txt)
  d2 <- drop_linear(txt)
  cat(sprintf("\n== %s (%s): %d models fitted\n", label, rds, n_all))
  cat(sprintf("   Length + Maturity together : %d\n", sum(d1)))
  cat(sprintf("   Year/Month/Length linear   : %d\n", sum(d2)))
  cat(sprintf("   (both reasons at once      : %d)\n", sum(d1 & d2)))
  cat(sprintf("   -> %d models left\n", sum(!(d1 | d2))))

  keep <- !(d1 | d2)
  kept <- x[keep, , drop = FALSE]
  txt  <- txt[keep]

  ord  <- order(kept$WAIC)
  kept <- kept[ord, , drop = FALSE]
  txt  <- txt[ord]

  top <- head(data.frame(
    Model      = paste0("M", seq_len(nrow(kept))),
    Covariates = vapply(txt, norm_covariates, character(1), USE.NAMES = FALSE),
    WAIC       = round(kept$WAIC, 2),
    LCPO       = round(kept$LCPO, 2),
    stringsAsFactors = FALSE), 10)

  print(top, row.names = FALSE)
  write.csv(top, paste0(out_stem, ".csv"), row.names = FALSE)

  if (requireNamespace("kableExtra", quietly = TRUE)) {
    tex <- kableExtra::kable_styling(
      kableExtra::kbl(top, format = "latex", booktabs = TRUE, longtable = TRUE,
                      align = "clcc",
                      caption = sprintf(paste0("Top 10 Models for %s Based on WAIC ",
                                               "and LCPO.\\\\label{tab:modelsel%s}"),
                                        label, substr(tolower(label), 1, 1))),
      font_size = 8, latex_options = c("hold_position"))
    writeLines(as.character(tex), paste0(out_stem, ".tex"))
  }
  invisible(top)
}

males   <- build("data/unique_modelsm.rds", "Model",   "Males",   "output/top10_males_filtered")
females <- build("data/unique_modelsf.rds", "Formula", "Females", "output/top10_females_filtered")

# ---- where the models the paper selects end up -------------------------------
rank_of <- function(rds, col, pattern, label) {
  x <- readRDS(rds)
  x <- x[!is.na(x$WAIC) & !is.na(x[[col]]), ]
  x <- x[order(x$WAIC), ]
  before <- grep(pattern, x[[col]])
  y <- x[!(drop_len_mat(x[[col]]) | drop_linear(x[[col]])), ]
  after <- grep(pattern, y[[col]])
  cat(sprintf("%-28s rank %s of %d  ->  rank %s of %d\n", label,
              ifelse(length(before), before, "-"), nrow(x),
              ifelse(length(after), after, "-"), nrow(y)))
}

cat("\n---- selected models, before -> after filtering ----\n")
rank_of("data/unique_modelsm.rds", "Model",
        "^Int \\+ Country_Code \\+ f\\(Year, model = 'rw2'\\) \\+ Maturity$",
        "Males  Country+Yr(rw2)+Mat")
rank_of("data/unique_modelsf.rds", "Formula",
        "^y ~ 1 \\+ Country_Code \\+ Maturity$", "Females Country+Maturity")
