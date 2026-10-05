# ------------------------------------------------------------------------------
# Builds the two composite figures needed to fit the display-item limit:
#
#   Figure 3 = (a) year effect in males by country   [was Figure 4]
#              (b) year effect in immatures by sex   [was Figure 5, recoloured]
#   Figure 4 = (a) posterior predictive trend, male model      [was Figure 6]
#              (b) posterior predictive trend, immature model  [was Figure 9]
#
# Panel (a) of Figure 4 is taken from the raster in the manuscript (862 px wide),
# so it is the resolution-limiting element; regenerate it from the .qmd at
# 300 dpi before submission if the journal objects.
#
# Outputs: Figure3_year_effects.svg/.png, Figure4_projections.png
# ------------------------------------------------------------------------------

suppressPackageStartupMessages(library(png))
inkscape <- "C:/Program Files/Inkscape/bin/inkscape.exe"

# ---- Figure 3: two vector panels stacked ------------------------------------
panel_svg <- function(f) {
  s <- paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  s <- sub("^.*?<svg[^>]*>", "", s)          # drop declaration and opening tag
  sub("</svg>\\s*$", "", s)
}
PW <- 560; PH <- 360; GAP <- 14; LAB <- 20   # panel size, gap, label strip

# Nature style: the panel marker is a bold lower-case letter on its own; the
# description of each panel lives in the legend, not in the figure.
panels <- list(
  list(f = "output/GA_year_effect_country.svg",         lab = "a"),
  list(f = "output/GA_year_effect_immatures_sex.svg",   lab = "b"))

H <- length(panels) * (PH + LAB) + GAP
out <- c('<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
         sprintf(paste0('<svg xmlns="http://www.w3.org/2000/svg" ',
                        'xmlns:xlink="http://www.w3.org/1999/xlink" ',
                        'width="%spx" height="%spx" viewBox="0 0 %s %s" ',
                        'font-family="Arial, Helvetica, sans-serif">'), PW, H, PW, H),
         sprintf('<rect width="%s" height="%s" fill="#FFFFFF"/>', PW, H))
y <- 0
for (p in panels) {
  out <- c(out,
    sprintf('<text x="6" y="%s" font-size="17" font-weight="bold" fill="#1A1A1A">%s</text>',
            y + 16, p$lab),
    sprintf('<g transform="translate(0,%s)">', y + LAB),
    panel_svg(p$f), '</g>')
  y <- y + PH + LAB + GAP
}
out <- c(out, '</svg>')
writeLines(out, "output/Figure3_year_effects.svg", useBytes = TRUE)

system2(inkscape, c(shQuote("output/Figure3_year_effects.svg"), "--export-type=png",
                    "--export-width=1950", "--export-background=white",
                    paste0("--export-filename=", shQuote("output/Figure3_year_effects.png"))),
        stdout = FALSE, stderr = FALSE)
cat("Figure3_year_effects.svg/.png written\n")

# ---- Figure 4: two rasters stacked, scaled to a common width ----------------
src <- c("a" = "output/projection_males.png",                 # from analysis.qmd
         "b" = "output/GA_figure9_projection_immatures.png")  # from R/06
labs <- c(a = "a", b = "b")

imgs <- lapply(src, readPNG)
W    <- max(vapply(imgs, function(x) dim(x)[2], numeric(1)))
scl  <- vapply(imgs, function(x) W / dim(x)[2], numeric(1))
hs   <- vapply(seq_along(imgs), function(i) dim(imgs[[i]])[1] * scl[i], numeric(1))
lab_h <- 40; gap <- 16
H    <- sum(hs) + length(imgs) * lab_h + gap

# y is measured from the top of the canvas; the viewport is bottom-up, so each
# element is placed at H - y (grid.raster itself must not be flipped)
png("output/Figure4_projections.png", width = W, height = H, res = 300, bg = "white")
grid::grid.newpage()
grid::pushViewport(grid::viewport(xscale = c(0, W), yscale = c(0, H)))
y <- 0
for (i in seq_along(imgs)) {
  grid::grid.text(labs[i], x = grid::unit(10, "native"),
                  y = grid::unit(H - (y + lab_h * 0.55), "native"),
                  just = c("left", "centre"), gp = grid::gpar(fontsize = 13, fontface = "bold"))
  grid::grid.raster(imgs[[i]], x = grid::unit(W / 2, "native"),
                    y = grid::unit(H - (y + lab_h + hs[i] / 2), "native"),
                    width = grid::unit(W, "native"), height = grid::unit(hs[i], "native"))
  y <- y + lab_h + hs[i] + gap
}
invisible(dev.off())
cat(sprintf("Figure4_projections.png written (%d x %d px)\n", W, round(H)))
