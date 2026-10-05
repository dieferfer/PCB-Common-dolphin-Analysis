# ------------------------------------------------------------------------------
# Is the lower PCB burden of immature females an artefact of classifying
# maturity by length? (Graham's comments 108 / 112)
#
# If some "immature" females were in fact mature (and had offloaded), the
# largest immature females would carry the LOWEST burdens. Under genuine
# accumulation, burden should rise, not fall, with length.
#
# Checks: 1. length distributions of immatures by sex
#         2. PCB vs length within immature females and males
#         3. immature model refitted with Length as a covariate
#         4. sensitivity: drop the largest immature females and refit
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(INLA))

d <- read.csv("data/PCBcorrectedDB.csv", sep = ";")
d$Country_Code <- factor(ifelse(d$Country_Code == 4, "France", "Spain"), levels = c("France", "Spain"))
d$Sex <- factor(ifelse(d$Sex == 1, "Males", "Females"), levels = c("Males", "Females"))
imm <- d[d$Maturity == 0, ]
mat <- d[d$Maturity != 0, ]

cat("=== 1. Length (cm) of immatures vs matures, by sex ===\n")
for (s in levels(d$Sex)) {
  a <- imm$Length[imm$Sex == s]; b <- mat$Length[mat$Sex == s]
  cat(sprintf("%-8s immature n=%2d  median %3.0f  range %3.0f-%3.0f  |  mature n=%2d  min %3.0f\n",
              s, length(a), median(a), min(a), max(a), length(b), min(b)))
}
cat("-> the largest immature of each sex sits right at the length cut-off used for maturity\n")

cat("\n=== 2. PCB vs length WITHIN immatures (Spearman) ===\n")
for (s in levels(d$Sex)) {
  x <- imm[imm$Sex == s, ]
  ct <- suppressWarnings(cor.test(x$Length, x$ICES73, method = "spearman"))
  cat(sprintf("%-8s n=%2d  rho = %+.2f  p = %.3f\n", s, nrow(x), ct$estimate, ct$p.value))
}
cat("\n   median PCB by length tercile (immature females):\n")
f <- imm[imm$Sex == "Females", ]
f$ter <- cut(f$Length, quantile(f$Length, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
             labels = c("smallest third", "middle third", "largest third"))
print(aggregate(ICES73 ~ ter, data = f, FUN = function(v) c(n = length(v), median = round(median(v), 1))))
cat("\n   same for immature males:\n")
m <- imm[imm$Sex == "Males", ]
m$ter <- cut(m$Length, quantile(m$Length, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
             labels = c("smallest third", "middle third", "largest third"))
print(aggregate(ICES73 ~ ter, data = m, FUN = function(v) c(n = length(v), median = round(median(v), 1))))

fit <- function(dat, form) inla(form, family = "gamma", data = dat,
                                control.compute = list(waic = TRUE))
fx <- function(m, row) { s <- m$summary.fixed[row, ]
  sprintf("%+.3f (%+.3f, %+.3f)  => x%.2f", s$mean, s$`0.025quant`, s$`0.975quant`, exp(s$mean)) }

cat("\n=== 3. Sex effect among immatures, with and without Length in the model ===\n")
m0 <- fit(imm, ICES73 ~ Country_Code + Sex + f(Year, model = "rw2"))
m1 <- fit(imm, ICES73 ~ Country_Code + Sex + Length + f(Year, model = "rw2"))
cat("without Length : SexFemales", fx(m0, "SexFemales"), " | WAIC", round(m0$waic$waic, 1), "\n")
cat("with Length    : SexFemales", fx(m1, "SexFemales"), " | WAIC", round(m1$waic$waic, 1), "\n")
cat("                 Length     ", fx(m1, "Length"), " per cm\n")

cat("\n=== 4. Sensitivity: drop the largest immature females (those closest to LSM) ===\n")
cut_f <- quantile(imm$Length[imm$Sex == "Females"], 0.75)
keep <- !(imm$Sex == "Females" & imm$Length > cut_f)
m2 <- fit(imm[keep, ], ICES73 ~ Country_Code + Sex + f(Year, model = "rw2"))
cat(sprintf("dropping immature females > %.0f cm (n dropped = %d): SexFemales %s\n",
            cut_f, sum(!keep), fx(m2, "SexFemales")))
keep2 <- !(imm$Sex == "Females" & imm$Length > median(imm$Length[imm$Sex == "Females"]))
m3 <- fit(imm[keep2, ], ICES73 ~ Country_Code + Sex + f(Year, model = "rw2"))
cat(sprintf("dropping immature females > median (%.0f cm, n dropped = %d): SexFemales %s\n",
            median(imm$Length[imm$Sex == "Females"]), sum(!keep2), fx(m3, "SexFemales")))
