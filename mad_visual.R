# Οπτικοποίηση της μέσης απόλυτης απόκλισης και της διακύμανσης, βήμα-βήμα:
#   1. η πυκνότητα f(x) και ο μέσος όρος μ
#   2. η «απόσταση» κάθε x από τον μέσο: |x − μ|  ή  (x − μ)²
#   3. το γινόμενο απόσταση·f(x) — το εμβαδόν του είναι το ζητούμενο μέτρο
#   4. το αποτέλεσμα πάνω στην κατανομή

# Συνάρτηση πυκνότητας πιθανότητας της N(0, 1)
f <- function(x) (1 / sqrt(2 * pi)) * exp(-x^2 / 2)

# Χρώματα
col_ink   <- "#0b0b0b"
col_muted <- "#898781"
col_grid  <- "#e1e0d9"
col_dens  <- "#2a78d6"   # πυκνότητα
col_dist  <- "#eb6834"   # απόσταση
col_prod  <- "#1baf7a"   # γινόμενο απόσταση·f(x)
fade <- function(col, a) adjustcolor(col, alpha.f = a)

fmt_bound <- function(v) if (is.infinite(v)) ifelse(v < 0, "−∞", "∞") else format(v)

# Γεμίζει το εμβαδόν κάτω από την fun στο [from, to]
shade <- function(fun, from, to, col) {
  xs <- seq(from, to, length.out = 300)
  polygon(c(from, xs, to), c(0, fun(xs), 0), col = col, border = NA)
}

# Τα δύο μέτρα διασποράς: τι «απόσταση» χρησιμοποιούν και πώς περιγράφονται
measures <- list(
  mad = list(
    name     = "Μέση απόλυτη απόκλιση",
    dist     = function(d) abs(d),
    dist_lab = "|x − μ|",
    ex_lab   = function(d) sprintf("%.1f", d),
    note2    = "πάντα θετική, και από τις δύο πλευρές",
    note3    = "τα δύο μισά είναι πάντα ίσα",
    width    = function(v) v,          # πόσο απλώνεται η ζώνη γύρω από το μ
    title4   = function(v) sprintf("4. Μέση απόλυτη απόκλιση = %.4f", v),
    band_lab = "μ ± απόκλιση"
  ),
  var = list(
    name     = "Διακύμανση",
    dist     = function(d) d^2,
    dist_lab = "(x − μ)²",
    ex_lab   = function(d) sprintf("%.1f² = %.2f", d, d^2),
    note2    = "τετράγωνο: οι μακρινές τιμές μετράνε πολύ περισσότερο",
    note3    = "τα δύο μισά δεν είναι απαραίτητα ίσα",
    width    = sqrt,                   # τυπική απόκλιση σ = √διακύμανση
    title4   = function(v) sprintf("4. Διακύμανση = %.4f   →   σ = √%.4f = %.4f", v, v, sqrt(v)),
    band_lab = "μ ± σ"
  )
)

# dens:      η συνάρτηση πυκνότητας (αν το εμβαδόν της στο [a, b] δεν είναι 1, κανονικοποιείται)
# dens_name: πώς αναφέρεται στον τίτλο
# examples:  αποστάσεις από το μ για τα δύο σημεία-παραδείγματα του βήματος 2
# dens_lab:  πώς λέγεται η πυκνότητα στις ετικέτες: f(x), ή g(x) όταν κανονικοποιείται σε διάστημα
visualize_spread <- function(kind, a, b, xlim, file, dens = f, dens_name = "της N(0, 1)",
                             examples = c(-0.8, 1.2), dens_lab = "f(x)") {
  m <- measures[[kind]]

  # Υπολογισμοί με ολοκληρώματα (όπως στα προηγούμενα scripts)
  Z   <- integrate(dens, a, b)$value
  g   <- function(x) ifelse(x >= a & x <= b, dens(x) / Z, 0)
  mu  <- integrate(function(x) x * g(x), a, b)$value
  h   <- function(x) m$dist(x - mu) * g(x)
  left  <- integrate(h, a, mu)$value   # συνεισφορά των x < μ
  right <- integrate(h, mu, b)$value   # συνεισφορά των x > μ
  value <- left + right
  w     <- m$width(value)

  lo <- max(a, xlim[1]); hi <- min(b, xlim[2])   # ορατό κομμάτι του [a, b]
  x  <- sort(c(seq(xlim[1], xlim[2], length.out = 800), lo, hi))

  # Κοινό πλαίσιο για κάθε πάνελ
  panel <- function(y_max, title, ylab, note = NULL) {
    plot(NA, xlim = xlim, ylim = c(0, y_max * 1.12), axes = FALSE,
         xlab = "x", ylab = ylab, main = "")
    title(main = title, adj = 0, line = 1.8, font.main = 1, cex.main = 1.25, col.main = col_ink)
    if (!is.null(note)) mtext(note, side = 3, line = 0.2, adj = 0, cex = 0.9, col = col_muted)
    abline(h = pretty(c(0, y_max)), col = col_grid)
    axis(1, col = col_muted, col.axis = col_muted)
    axis(2, las = 1, col = col_muted, col.axis = col_muted)
    abline(v = mu, col = col_ink, lty = 2)
  }

  png(file, width = 1500, height = 1050, res = 110)
  par(mfrow = c(2, 2), mar = c(4.2, 4.5, 3.8, 1.5), oma = c(0, 0, 3, 0),
      col.lab = col_muted)

  # 1. Πυκνότητα και μέσος όρος
  y1 <- max(g(x))
  panel(y1, sprintf("1. Η πυκνότητα %s και ο μέσος όρος μ", dens_lab), dens_lab)
  shade(g, lo, hi, fade(col_dens, 0.15))
  lines(x, g(x), col = col_dens, lwd = 2)
  text(mu, y1 * 1.1, sprintf("μ = %.3f", mu), pos = 4, col = col_ink)

  # 2. Απόσταση από τον μέσο όρο
  y2 <- max(m$dist(xlim - mu))
  panel(y2, sprintf("2. Η απόσταση κάθε x από τον μέσο: %s", m$dist_lab), m$dist_lab, m$note2)
  # Αχνά από πίσω το σχήμα της g(x), για να φαίνεται ποιες αποστάσεις θα «ζυγίσουν»
  # στο βήμα 3. Μόνο το σχήμα: ο άξονας y εδώ μετράει απόσταση, όχι πυκνότητα.
  s  <- 0.7 * y2 / y1
  polygon(c(x[1], x, x[length(x)]), c(0, g(x) * s, 0), col = fade(col_dens, 0.1), border = NA)
  lines(x, g(x) * s, col = fade(col_dens, 0.4), lwd = 1.5)
  legend("top", inset = 0.01, cex = 0.85, bg = fade("white", 0.9), box.col = col_grid,
         legend = c(sprintf("%s: η γραφική παράσταση", m$dist_lab),
                    sprintf("αχνό: σχήμα της %s, μόνο για σύγκριση —", dens_lab),
                    "δεν είναι μέρος της γραφικής παράστασης"),
         col = c(col_dist, fade(col_dens, 0.35), NA),
         lty = c(1, NA, NA), lwd = c(2, NA, NA), pch = c(NA, 15, NA), pt.cex = 2)
  xin <- x[x >= lo & x <= hi]
  lines(x, m$dist(x - mu), col = fade(col_dist, 0.3), lwd = 2)   # εκτός [a, b]: δεν μετράει
  lines(xin, m$dist(xin - mu), col = col_dist, lwd = 2)
  for (x0 in mu + examples) {
    d <- abs(x0 - mu)
    segments(x0, 0, x0, m$dist(d), col = col_dist, lty = 3)
    points(x0, m$dist(d), pch = 19, col = col_dist, cex = 1.2)
    arrows(mu, 0.06 * y2, x0, 0.06 * y2, code = 3, length = 0.07, col = col_ink)
    text(x0, m$dist(d), m$ex_lab(d), pos = if (x0 < mu) 2 else 4, col = col_ink)
  }

  # 3. Γινόμενο: κάθε απόσταση ζυγισμένη με την πυκνότητα
  y3 <- max(h(x))
  panel(y3, sprintf("3. Γινόμενο %s·%s   →   εμβαδόν = %.4f", m$dist_lab, dens_lab, value),
        sprintf("%s·%s", m$dist_lab, dens_lab), m$note3)
  shade(h, lo, mu, fade(col_prod, 0.25))
  shade(h, mu, hi, fade(col_prod, 0.5))
  lines(x, h(x), col = col_prod, lwd = 2)
  cx_left  <- integrate(function(x) x * h(x), a, mu)$value / left
  cx_right <- integrate(function(x) x * h(x), mu, b)$value / right
  text(cx_left,  0.3 * h(cx_left),  sprintf("αριστερά\n%.4f", left),  col = col_ink)
  text(cx_right, 0.3 * h(cx_right), sprintf("δεξιά\n%.4f",    right), col = col_ink)

  # 4. Το αποτέλεσμα πάνω στην κατανομή
  panel(y1, m$title4(value), dens_lab)
  shade(g, max(mu - w, lo), min(mu + w, hi), fade(col_dens, 0.3))
  lines(x, g(x), col = col_dens, lwd = 2)
  arrows(mu, y1 * 0.45, mu - w, y1 * 0.45, length = 0.08, col = col_ink, lwd = 1.5)
  arrows(mu, y1 * 0.45, mu + w, y1 * 0.45, length = 0.08, col = col_ink, lwd = 1.5)
  text(mu - w / 2, y1 * 0.45, sprintf("%.3f", w), pos = 3, col = col_ink)
  text(mu + w / 2, y1 * 0.45, sprintf("%.3f", w), pos = 3, col = col_ink)
  text(mu, y1 * 1.1, sprintf("%s = [%.3f, %.3f]", m$band_lab, mu - w, mu + w),
       pos = 4, col = col_ink)

  interval <- sprintf("%s%s, %s%s", if (is.infinite(a)) "(" else "[", fmt_bound(a),
                      fmt_bound(b), if (is.infinite(b)) ")" else "]")
  mtext(sprintf("%s %s στο διάστημα %s", m$name, dens_name, interval),
        outer = TRUE, line = 1, cex = 1.3, font = 2, col = col_ink)
  invisible(dev.off())
  cat(sprintf("%s %s: μ = %.6f, αριστερά = %.6f, δεξιά = %.6f, σύνολο = %.6f → %s\n",
              m$name, interval, mu, left, right, value, file))
}

# Ολόκληρη η κατανομή και το διάστημα [-0.5, 2] του δεύτερου παραδείγματος
visualize_spread("mad", -Inf, Inf, xlim = c(-4, 4),   file = "mad_steps_full.png")
visualize_spread("mad", -0.5,  2,  xlim = c(-1.5, 3), file = "mad_steps_interval.png", dens_lab = "g(x)")
visualize_spread("var", -Inf, Inf, xlim = c(-4, 4),   file = "var_steps_full.png")
visualize_spread("var", -0.5,  2,  xlim = c(-1.5, 3), file = "var_steps_interval.png", dens_lab = "g(x)")

# Τριγωνική κατανομή: ισοσκελές τρίγωνο με κορυφή στο 0, από το -1 έως το +1.
# f(x) = 1 − |x|  — εμβαδόν = βάση·ύψος/2 = 2·1/2 = 1, άρα είναι ήδη πυκνότητα.
triangle <- function(x) pmax(1 - abs(x), 0)
visualize_spread("mad", -1, 1, xlim = c(-1.5, 1.5), file = "mad_steps_triangle.png",
                 dens = triangle, dens_name = "της τριγωνικής κατανομής",
                 examples = c(-0.4, 0.6))
