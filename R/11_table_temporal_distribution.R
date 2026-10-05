# ------------------------------------------------------------------------------
# Temporal distribution of the samples by country, sex and maturity.
# Checks whether sampling year is confounded with the grouping variables.
#
# Outputs: table_temporal_distribution.csv / .tex   (periods x group)
#          table_samples_per_year.csv                (year x country, supplement)
# ------------------------------------------------------------------------------

setwd("D:/Usuarios/dfernandez/Dropbox/Diego/PCB DDE paper")
d <- read.csv("data/PCBcorrectedDB.csv", sep = ";")
d$Country <- factor(ifelse(d$Country_Code == 4, "France", "Spain"), levels = c("France", "Spain"))
d$Sex_l <- factor(c("Males", "Females", "Pregnant females")[d$Sex],
                  levels = c("Males", "Females", "Pregnant females"))
d$Mat <- factor(ifelse(d$Maturity == 0, "Immature", "Mature"), levels = c("Immature", "Mature"))

# BIOCET period (2001-2004), then five-year blocks
d$Period <- cut(d$Year, breaks = c(2000, 2004, 2009, 2014, 2019),
                labels = c("2001-2004", "2005-2009", "2010-2014", "2015-2019"))

# ---- main table: group x period, with n, median and range of year ----------
grp <- interaction(d$Country, d$Sex_l, d$Mat, sep = "|", drop = TRUE)
tab <- as.data.frame.matrix(table(grp, d$Period))
key <- do.call(rbind, strsplit(rownames(tab), "|", fixed = TRUE))
out <- data.frame(Country = key[, 1], Sex = key[, 2], Maturity = key[, 3], tab,
                  check.names = FALSE, stringsAsFactors = FALSE)
out$n <- rowSums(tab)
yrs <- split(d$Year, grp)[rownames(tab)]
out$`Median year` <- vapply(yrs, function(y) sprintf("%.0f", median(y)), character(1))
out$`Year range` <- vapply(yrs, function(y) sprintf("%d-%d", min(y), max(y)), character(1))
out <- out[order(match(out$Country, levels(d$Country)), match(out$Sex, levels(d$Sex_l)),
                 match(out$Maturity, levels(d$Mat))), ]
rownames(out) <- NULL

# country subtotals and grand total
sub <- function(sel, label) {
  y <- d$Year[sel]; p <- table(d$Period[sel])
  data.frame(Country = label, Sex = "", Maturity = "", t(as.matrix(p)), n = sum(p),
             `Median year` = sprintf("%.0f", median(y)),
             `Year range` = sprintf("%d-%d", min(y), max(y)),
             check.names = FALSE, stringsAsFactors = FALSE)
}
names_fix <- function(x) { names(x) <- names(out); x }
final <- rbind(out[out$Country == "France", ], names_fix(sub(d$Country == "France", "France (all)")),
               out[out$Country == "Spain", ],  names_fix(sub(d$Country == "Spain", "Spain (all)")),
               names_fix(sub(rep(TRUE, nrow(d)), "Total")))
stopifnot(final$n[nrow(final)] == nrow(d))

print(final, row.names = FALSE)
write.csv(final, "output/table_temporal_distribution.csv", row.names = FALSE)
write.csv(addmargins(table(Year = d$Year, Country = d$Country)),
          "output/table_samples_per_year.csv")

if (requireNamespace("kableExtra", quietly = TRUE)) {
  k <- kableExtra::kable_styling(
    kableExtra::kbl(final, format = "latex", booktabs = TRUE, align = "lllrrrrrcc",
      caption = paste("Temporal distribution of the sampled common dolphins by country, sex",
                      "and inferred maturity status. Periods: BIOCET project (2001-2004) and",
                      "subsequent five-year blocks.\\\\label{tab:temporal}")),
    latex_options = "hold_position", font_size = 8)
  writeLines(as.character(k), "output/table_temporal_distribution.tex")
}

# ---- is sampling year confounded with country, within sex? ------------------
cat("\n=== Year by country, within sex (Wilcoxon rank-sum) ===\n")
for (s in levels(d$Sex_l)) {
  a <- d$Year[d$Sex_l == s & d$Country == "France"]
  b <- d$Year[d$Sex_l == s & d$Country == "Spain"]
  if (length(a) < 2 || length(b) < 2) { cat(sprintf("%-17s not testable (n = %d / %d)\n", s, length(a), length(b))); next }
  w <- suppressWarnings(wilcox.test(a, b))
  cat(sprintf("%-17s France median %.0f (n=%d) vs Spain median %.0f (n=%d)  p = %.2g\n",
              s, median(a), length(a), median(b), length(b), w$p.value))
}
cat("\n=== share of each country's sample taken in 2001-2002 ===\n")
for (c_ in levels(d$Country)) {
  y <- d$Year[d$Country == c_]
  cat(sprintf("%-7s %d of %d (%.0f%%)\n", c_, sum(y <= 2002), length(y), 100 * mean(y <= 2002)))
}
cat("\nwritten: table_temporal_distribution.csv/.tex, table_samples_per_year.csv\n")
