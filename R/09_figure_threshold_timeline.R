# ------------------------------------------------------------------------------
# Graphical abstract - years in which each group drops below the PCB thresholds,
# now including immature females from the immature-only model.
#
# Same design as threshold_timeline.R (two progress bars, France above and
# Spain below), with three groups per country instead of two:
#   immature females (immature model), immature males and mature males
#   (male model). A sex glyph next to the dolphin tells the groups apart.
# Years of one country that would collide are merged into one label
# ("2020-21") that lists all the symbols of the cluster.
#
# Output: GA_timeline_thresholds_immF.svg
# ------------------------------------------------------------------------------

out_svg <- "output/GA_timeline_thresholds_immF.svg"

# ---- results: first year with P(mean PCB < threshold) >= 0.8 -----------------
# males:   male model (Figure 8);  immature females: immature model (Figure 9)
events <- rbind(
  data.frame(thr = 17, country = "Spain",  group = "ImmF", year = 2020),
  data.frame(thr = 17, country = "Spain",  group = "ImmM", year = 2020),
  data.frame(thr = 17, country = "Spain",  group = "MatM", year = 2028),
  data.frame(thr = 17, country = "France", group = "ImmF", year = 2032),
  data.frame(thr = 17, country = "France", group = "ImmM", year = 2046),
  data.frame(thr = 17, country = "France", group = "MatM", year = 2063),
  data.frame(thr = 41, country = "Spain",  group = "ImmF", year = 2020),
  data.frame(thr = 41, country = "Spain",  group = "ImmM", year = 2020),
  data.frame(thr = 41, country = "Spain",  group = "MatM", year = 2020),
  data.frame(thr = 41, country = "France", group = "ImmF", year = 2020),
  data.frame(thr = 41, country = "France", group = "ImmM", year = 2021),
  data.frame(thr = 41, country = "France", group = "MatM", year = 2030)
)
grp_order <- c("ImmF", "ImmM", "MatM")

# ---- look --------------------------------------------------------------------
col_france <- "#1F5AA8"
col_spain  <- "#4E9E52"
col_thr    <- "#2E7B84"
col_track  <- "#E6E8EA"
col_text   <- "#1A1A1A"
col_axis   <- "#6A6A6A"

W <- 820; H <- 350
x0 <- 196; x1 <- W - 60            # bar starts after the threshold label box
bar_y <- c(138, 264)               # one bar per threshold
bar_h <- 13

year_size <- 21                    # the years carry the message, so they are big
dig_w     <- 0.556 * year_size     # Arial digit advance, to centre each label
len_imm   <- 26                    # dolphin length: immature
len_mat   <- 42                    #                 mature
glyph_w   <- 10                   # sex glyph width next to the dolphin

# The year axis is not drawn, so the left end only needs to leave enough room
# for the first label (2020) to clear the country name and the bar start.
yr_lim <- c(2013, 2066)
sx <- function(v) round(x0 + (v - yr_lim[1]) / diff(yr_lim) * (x1 - x0), 1)

esc <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;",  s, fixed = TRUE)
  gsub(">", "&gt;", s, fixed = TRUE)
}
layer <- function(id, ...) {
  c(sprintf('<g inkscape:groupmode="layer" id="%s" inkscape:label="%s">', id, id),
    unlist(list(...)), '</g>')
}
txt <- function(x, y, s, size = 11, anchor = "middle", weight = "normal",
                colour = col_text, raw = FALSE) {
  sprintf(paste0('<text x="%s" y="%s" text-anchor="%s" font-size="%s" ',
                 'font-weight="%s" fill="%s">%s</text>'),
          x, y, anchor, size, weight, colour, if (raw) s else esc(s))
}
marker <- function(x, y, col) {
  sprintf('<circle cx="%s" cy="%s" r="8.5" fill="%s"/>', x, y, col)
}

# ---- dolphin symbol: small = immature, large = mature -------------------------
sil <- paste(readLines("R/dolphin_silhouette.svg", warn = FALSE), collapse = " ")
sil_d <- gsub("\\s+", " ", sub('.*<path d="([^"]+)".*', "\\1", sil))
sil_box <- as.numeric(strsplit(
  trimws(sub('.*viewBox="([^"]+)".*', "\\1", sil)), "\\s+")[[1]])
sil_w <- sil_box[3]; sil_h <- sil_box[4]
sil_ar <- sil_h / sil_w

dolphin <- function(x, y, len, col) {
  s <- len / sil_w
  sprintf(paste0('<g transform="translate(%s,%s) scale(%s,%s)">',
                 '<path d="%s" fill="%s" fill-rule="evenodd"/></g>'),
          round(x + len, 2), round(y - len * sil_ar / 2, 2),
          round(-s, 6), round(s, 6), sil_d, col)
}

# sex glyphs drawn as clean vector strokes (font glyphs render ragged in
# PowerPoint). Both are 10 px wide; (x, y) is the left edge / vertical centre.
glyph <- function(sex, x, y, col) {
  r <- 3.1; cx <- x + 4.4; sw <- 1.5
  if (sex == "F") {
    cy <- y - 2.2
    sprintf(paste0('<g fill="none" stroke="%s" stroke-width="%s" stroke-linecap="round">',
                   '<circle cx="%s" cy="%s" r="%s"/>',
                   '<line x1="%s" y1="%s" x2="%s" y2="%s"/>',       # stem
                   '<line x1="%s" y1="%s" x2="%s" y2="%s"/></g>'),  # crossbar
            col, sw, cx, cy, r,
            cx, cy + r, cx, cy + r + 4.6,
            cx - 2.6, cy + r + 2.6, cx + 2.6, cy + r + 2.6)
  } else {
    cy <- y + 1.4; cx <- x + 3.6
    ax <- cx + r * 0.707; ay <- cy - r * 0.707             # circle edge, 45 deg
    tx <- ax + 4.2; ty <- ay - 4.2                         # arrow tip
    sprintf(paste0('<g fill="none" stroke="%s" stroke-width="%s" stroke-linecap="round" stroke-linejoin="round">',
                   '<circle cx="%s" cy="%s" r="%s"/>',
                   '<line x1="%s" y1="%s" x2="%s" y2="%s"/>',       # shaft
                   '<polyline points="%s,%s %s,%s %s,%s"/></g>'),  # head
            col, sw, cx, cy, r,
            ax, ay, tx, ty,
            tx - 3.4, ty, tx, ty, tx, ty + 3.4)
  }
}

# one group symbol = dolphin + sex glyph at its tail; returns svg and its width
sym_len <- c(ImmF = len_imm, ImmM = len_imm, MatM = len_mat)
sym_sex <- c(ImmF = "F", ImmM = "M", MatM = "M")
sym_w <- function(g) sym_len[g] + 3 + glyph_w
symbol <- function(g, x, y, col) {
  L <- sym_len[g]
  c(dolphin(x, y, L, col), glyph(sym_sex[g], x + L + 3, y, col))
}

# ---- cluster the years of one country whose labels would collide -------------
label_w <- function(grp) max(4 * dig_w, sum(sym_w(grp)) + 6 * (length(grp) - 1))
clusters <- function(ec) {
  ys <- sort(unique(ec$year)); cl <- list(); cur <- ys[1]
  for (yy in ys[-1]) {
    if (sx(yy) - sx(max(cur)) < 4 * dig_w + 6) cur <- c(cur, yy) else { cl <- c(cl, list(cur)); cur <- yy }
  }
  c(cl, list(cur))
}

bar_layer <- function(thr, y) {
  e <- events[events$thr == thr, ]
  last_es <- max(e$year[e$country == "Spain"])
  last_fr <- max(e$year[e$country == "France"])

  out <- c(
    sprintf('<rect x="%s" y="%s" width="%s" height="%s" rx="%s" fill="%s"/>',
            x0, y - bar_h / 2, round(x1 - x0, 1), bar_h, bar_h / 2, col_track),
    sprintf('<rect x="%s" y="%s" width="%s" height="%s" rx="%s" fill="%s"/>',
            x0, y - bar_h / 2, round(sx(last_fr) - x0, 1), bar_h, bar_h / 2, col_france),
    sprintf('<rect x="%s" y="%s" width="%s" height="%s" rx="%s" fill="%s"/>',
            x0, y - bar_h / 2, round(sx(last_es) - x0, 1), bar_h, bar_h / 2, col_spain),
    sprintf('<rect x="%s" y="%s" width="150" height="46" rx="8" fill="%s"/>',
            26, y - 23, col_thr),
    txt(101, y - 4, sprintf("Threshold %d", thr), 13, weight = "bold", colour = "#FFFFFF"),
    txt(101, y + 13, "mg/kg lw", 12, colour = "#FFFFFF"))

  for (cty in c("France", "Spain")) {
    ec  <- e[e$country == cty, ]
    col <- if (cty == "France") col_france else col_spain
    up  <- cty == "France"

    for (cl in clusters(ec)) {
      grp <- grp_order[grp_order %in% ec$group[ec$year %in% cl]]
      xc  <- mean(sx(cl))                          # label centre
      for (yy in cl) out <- c(out, marker(sx(yy), y, col))

      lab <- if (length(cl) == 1) as.character(cl) else
        sprintf("%d–%02d", min(cl), max(cl) %% 100)
      ty   <- if (up) y - 28 else y + 50           # baseline of the year
      ysym <- if (up) ty - 23 else y + 26

      wsym <- sum(sym_w(grp)) + 6 * (length(grp) - 1)
      xs <- xc - wsym / 2
      for (g in grp) { out <- c(out, symbol(g, xs, ysym, col)); xs <- xs + sym_w(g) + 6 }
      for (yy in cl) out <- c(out,
        sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1"/>',
                sx(yy), if (up) y - 11 else y + 11, sx(yy), if (up) y - 21 else y + 17, col))
      out <- c(out, txt(round(xc, 1), ty, lab, year_size, anchor = "middle",
                        weight = "bold", colour = col))
    }
    out <- c(out, txt(x0 - 10, if (up) y - 28 else y + 50, cty, 12,
                      anchor = "end", weight = "bold", colour = col))
  }
  layer(sprintf("threshold-%d", thr), out)
}

# ---- layout check ------------------------------------------------------------
for (t in unique(events$thr)) for (cty in c("France", "Spain")) {
  ec <- events[events$thr == t & events$country == cty, ]
  spans <- lapply(clusters(ec), function(cl) {
    grp <- ec$group[ec$year %in% cl]; w <- label_w(grp)
    if (length(cl) > 1) w <- max(w, 7 * dig_w)
    c(mean(sx(cl)) - w / 2, mean(sx(cl)) + w / 2)
  })
  for (i in seq_along(spans)) {
    if (spans[[i]][1] < x0 - 10) warning(sprintf("thr %d %s: label starts left of the bar", t, cty))
    if (spans[[i]][2] > W - 6)   warning(sprintf("thr %d %s: label runs off the canvas", t, cty))
    if (i > 1 && spans[[i]][1] < spans[[i - 1]][2] + 6)
      warning(sprintf("thr %d %s: labels %d and %d overlap", t, cty, i - 1, i))
  }
}

key_x <- c(26, 196, 356)
svg <- c(
  '<?xml version="1.0" encoding="UTF-8" standalone="no"?>',
  sprintf(paste0('<svg xmlns="http://www.w3.org/2000/svg" ',
                 'xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" ',
                 'width="%spx" height="%spx" viewBox="0 0 %s %s" ',
                 'font-family="Arial, Helvetica, sans-serif">'), W, H, W, H),
  layer("title",
    txt(26, 34, "Year at which each group falls below the threshold", 15,
        anchor = "start", weight = "bold")),
  layer("symbol-key",
    symbol("ImmF", key_x[1], 60, "#7A7A7A"),
    txt(key_x[1] + sym_w("ImmF") + 5, 65, "immature females", 13, anchor = "start", colour = "#4A4A4A"),
    symbol("ImmM", key_x[2], 60, "#7A7A7A"),
    txt(key_x[2] + sym_w("ImmM") + 5, 65, "immature males", 13, anchor = "start", colour = "#4A4A4A"),
    symbol("MatM", key_x[3], 60, "#7A7A7A"),
    txt(key_x[3] + sym_w("MatM") + 5, 65, "mature males", 13, anchor = "start", colour = "#4A4A4A")),
  bar_layer(17, bar_y[1]),
  bar_layer(41, bar_y[2]),
  '</svg>')

writeLines(svg, out_svg, useBytes = TRUE)
cat(sprintf("written: %s (%s x %s px)\n", out_svg, W, H))
