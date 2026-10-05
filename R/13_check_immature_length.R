# ------------------------------------------------------------------------------
# Does body length interact with sex among immatures? (follow-up to the size
# check requested by Graham). Compares the selected immature model with
# versions adding Length as a main effect and as a Sex x Length interaction.
# LCPO = -sum(log CPO), as in the model-selection tables of the paper.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages(library(INLA))

d <- read.csv("data/PCBcorrectedDB.csv", sep = ";")
d$Country_Code <- factor(ifelse(d$Country_Code == 4, "France", "Spain"), levels = c("France", "Spain"))
d$Sex <- factor(ifelse(d$Sex == 1, "Males", "Females"), levels = c("Males", "Females"))
imm <- d[d$Maturity == 0, ]
imm$Length_c <- imm$Length - mean(imm$Length)     # centred, so Sex = effect at mean length

fit <- function(form) inla(form, family = "gamma", data = imm,
                           control.compute = list(waic = TRUE, dic = TRUE, cpo = TRUE))
models <- list(
  "M1  Country + Sex + f(Year)"                     = fit(ICES73 ~ Country_Code + Sex + f(Year, model = "rw2")),
  "M1L Country + Sex + Length + f(Year)"            = fit(ICES73 ~ Country_Code + Sex + Length_c + f(Year, model = "rw2")),
  "M1X Country + Sex * Length + f(Year)"            = fit(ICES73 ~ Country_Code + Sex * Length_c + f(Year, model = "rw2")),
  "M1R Country + Sex + f(Length, rw2) + f(Year)"    = fit(ICES73 ~ Country_Code + Sex + f(Length, model = "rw2") + f(Year, model = "rw2")))

cat(sprintf("%-46s %8s %8s %8s %6s\n", "model", "WAIC", "DIC", "LCPO", "p.eff"))
for (nm in names(models)) {
  m <- models[[nm]]
  cat(sprintf("%-46s %8.2f %8.2f %8.2f %6.2f\n", nm, m$waic$waic, m$dic$dic,
              -sum(log(m$cpo$cpo), na.rm = TRUE), m$waic$p.eff))
}

show <- function(m, rows) {
  s <- m$summary.fixed[rows, c("mean", "0.025quant", "0.975quant"), drop = FALSE]
  out <- round(s, 4)
  out$relevant <- ifelse(s[, 2] > 0 | s[, 3] < 0, "yes", "no")
  print(out)
}
cat("\n--- M1L: Length as main effect ---\n"); show(models[[2]], c("SexFemales", "Length_c"))
cat("\n--- M1X: Sex x Length interaction ---\n")
show(models[[3]], c("SexFemales", "Length_c", "SexFemales:Length_c"))

# implied length slopes by sex (log scale per cm), from the interaction model
b <- models[[3]]$summary.fixed
cat(sprintf("\nimplied slope per cm (log scale):  males %+.4f   females %+.4f\n",
            b["Length_c", "mean"], b["Length_c", "mean"] + b["SexFemales:Length_c", "mean"]))
cat(sprintf("over the immature length range (%d-%d cm, %d cm): males x%.2f, females x%.2f\n",
            min(imm$Length), max(imm$Length), diff(range(imm$Length)),
            exp(b["Length_c", "mean"] * diff(range(imm$Length))),
            exp((b["Length_c", "mean"] + b["SexFemales:Length_c", "mean"]) * diff(range(imm$Length)))))
