# Μέσος όρος, διάμεσος και μέση απόλυτη απόκλιση των ΤΙΜΩΝ μιας συνάρτησης h(x)
# σε ένα διάστημα [a, b], με ολοκληρώματα.
#
# Διαφορά από τα προηγούμενα scripts: εδώ η καμπύλη δεν είναι πυκνότητα (βάρος).
# Η καμπύλη δίνει τις ίδιες τις τιμές y = h(x), και κάθε x του [a, b] έχει
# το ίδιο βάρος 1/(b − a).
#
#   μέσος όρος             ȳ = 1/(b − a) · ∫_a^b h(x) dx
#   διάμεσος               M : για το μισό μήκος του [a, b] ισχύει h(x) ≤ M
#                          για αύξουσα h αυτό σημαίνει  M = h((a + b) / 2)
#   μέση απόλυτη απόκλιση  1/(b − a) · ∫_a^b |h(x) − ȳ| dx

# Χρώματα
col_ink    <- "#0b0b0b"
col_muted  <- "#898781"
col_grid   <- "#e1e0d9"
col_fn     <- "#2a78d6"   # η συνάρτηση h(x)
col_dist   <- "#eb6834"   # απόσταση |h(x) − ȳ|
col_area   <- "#1baf7a"   # εμβαδόν ανάμεσα σε h(x) και ȳ
col_mean   <- "#e34948"   # μέσος όρος
col_median <- "#008300"   # διάμεσος
fade <- function(col, a) adjustcolor(col, alpha.f = a)

# Όλοι οι υπολογισμοί
describe_values <- function(h, a, b) {
  x <- seq(a, b, length.out = 1000)
  stopifnot(all(diff(h(x)) >= -1e-12))   # η h πρέπει να είναι αύξουσα

  L    <- b - a
  ybar <- integrate(h, a, b)$value / L
  mid  <- (a + b) / 2
  M    <- h(mid)
  # Το σημείο όπου η h περνάει το ȳ: αριστερά του h < ȳ, δεξιά του h > ȳ
  xs    <- uniroot(function(x) h(x) - ybar, c(a, b), tol = 1e-10)$root
  below <- integrate(function(x) ybar - h(x), a, xs)$value   # εμβαδόν κάτω από το ȳ
  above <- integrate(function(x) h(x) - ybar, xs, b)$value   # εμβαδόν πάνω από το ȳ
  mad   <- (below + above) / L
  # Έλεγχος: το ίδιο ολοκλήρωμα απευθείας με απόλυτη τιμή
  mad_direct <- integrate(function(x) abs(h(x) - ybar), a, b)$value / L

  list(L = L, ybar = ybar, mid = mid, M = M, xs = xs,
       below = below, above = above, mad = mad, mad_direct = mad_direct)
}

# Κοινό πλαίσιο για κάθε πάνελ
panel <- function(xlim, ylim, title, ylab, note = NULL) {
  plot(NA, xlim = xlim, ylim = ylim, axes = FALSE, xlab = "x", ylab = ylab, main = "")
  title(main = title, adj = 0, line = 1.8, font.main = 1, cex.main = 1.25, col.main = col_ink)
  if (!is.null(note)) mtext(note, side = 3, line = 0.2, adj = 0, cex = 0.9, col = col_muted)
  abline(h = axTicks(2), col = col_grid)
  axis(1, col = col_muted, col.axis = col_muted)
  axis(2, las = 1, col = col_muted, col.axis = col_muted)
}

# Γεμίζει την περιοχή ανάμεσα σε δύο καμπύλες y1(x), y2(x) στο [from, to]
fill_between <- function(from, to, y1, y2, col) {
  xs <- seq(from, to, length.out = 300)
  polygon(c(xs, rev(xs)), c(y1(xs), rev(y2(xs))), col = col, border = NA)
}

# Πάνελ μέσου όρου: η περιοχή κάτω από την h και το ορθογώνιο ύψους ȳ με το ίδιο εμβαδόν
draw_mean_panel <- function(h, a, b, title, note) {
  s  <- describe_values(h, a, b)
  x  <- seq(a, b, length.out = 500)
  panel(c(a, b), c(min(0, h(a)), h(b) * 1.12), title, "h(x)", note)
  fill_between(a, b, h, function(x) rep(0, length(x)), fade(col_fn, 0.18))
  rect(a, 0, b, s$ybar, border = col_mean, lty = 2, lwd = 1.5)
  lines(x, h(x), col = col_fn, lwd = 2)
  segments(a, s$ybar, b, s$ybar, col = col_mean, lwd = 2)
  text(a, s$ybar, sprintf("  ȳ = %.4f", s$ybar), adj = c(0, -0.6), col = col_ink)
}

# ---------------------------------------------------------------------------
# Εικόνα Riemann: ο μέσος όρος n ισαπεχουσών τιμών → ολοκλήρωμα καθώς n → ∞
# ---------------------------------------------------------------------------
plot_riemann <- function(h, a, b, name, interval_lab, ns, file) {
  s  <- describe_values(h, a, b)
  x  <- seq(a, b, length.out = 500)
  yl <- c(min(0, h(a)), h(b) * 1.12)
  png(file, width = 1800, height = 620, res = 110)
  par(mfrow = c(1, length(ns)), mar = c(4.2, 4.5, 3.8, 1.5), oma = c(0, 0, 3, 0), col.lab = col_muted)
  par(cex = 1)   # το mfrow με 3 πάνελ μικραίνει αυτόματα τα γράμματα
  for (n in ns) {
    w   <- (b - a) / n
    xi  <- a + (seq_len(n) - 0.5) * w        # τα μέσα των n κομματιών
    avg <- mean(h(xi))                        # Σ h(xᵢ) / n
    panel(c(a, b), yl, sprintf("n = %d τιμές", n), "h(x)",
          sprintf("πλάτος ορθογωνίου (b − a)/n = %g", w))
    rect(xi - w / 2, 0, xi + w / 2, h(xi), col = fade(col_fn, 0.25), border = "white", lwd = 1)
    lines(x, h(x), col = col_fn, lwd = 2)
    if (n <= 12) points(xi, h(xi), pch = 19, col = col_fn, cex = 1.1)
    segments(a, avg, b, avg, col = col_mean, lwd = 2, lty = 2)
    text(a, avg, sprintf("  Σ h(xᵢ)/n = %.4f", avg), adj = c(0, -0.6), col = col_ink)
  }
  mtext(sprintf("Ο μέσος όρος n τιμών της %s πλησιάζει το 1/(b − a)·∫ h(x) dx = %.4f",
                name, s$ybar),
        outer = TRUE, line = 1, cex = 1.3, font = 2, col = col_ink)
  invisible(dev.off())
}

# ---------------------------------------------------------------------------
# Εικόνα 0: μόνο ο μέσος όρος, πολλές συναρτήσεις δίπλα-δίπλα
# ---------------------------------------------------------------------------
plot_means_row <- function(fns, a, b, interval_lab, file) {
  png(file, width = 1800, height = 620, res = 110)
  par(mfrow = c(1, length(fns)), mar = c(4.2, 4.5, 3.8, 1.5), oma = c(0, 0, 3, 0), col.lab = col_muted)
  par(cex = 1)   # το mfrow με 3 πάνελ μικραίνει αυτόματα τα γράμματα
  for (fn in fns) {
    s <- describe_values(fn$h, a, b)
    draw_mean_panel(fn$h, a, b, fn$name,
                    sprintf("εμβαδόν = %.4f = ȳ · (b − a)", s$ybar * s$L))
  }
  mtext(sprintf("Μέσος όρος των τιμών: ορθογώνιο ύψους ȳ με το ίδιο εμβαδόν, στο %s", interval_lab),
        outer = TRUE, line = 1, cex = 1.3, font = 2, col = col_ink)
  invisible(dev.off())
}

# ---------------------------------------------------------------------------
# Εικόνα 1: μέσος όρος και διάμεσος
# ---------------------------------------------------------------------------
plot_mean_median <- function(h, a, b, name, interval_lab, file) {
  s  <- describe_values(h, a, b)
  x  <- seq(a, b, length.out = 500)
  yl <- c(min(0, h(a)), h(b) * 1.12)
  const <- function(v) function(x) rep(v, length(x))

  png(file, width = 1500, height = 620, res = 110)
  par(mfrow = c(1, 2), mar = c(4.2, 4.5, 3.8, 1.5), oma = c(0, 0, 3, 0), col.lab = col_muted)

  # Αριστερά: μέσος όρος = ύψος ορθογωνίου με το ίδιο εμβαδόν
  draw_mean_panel(h, a, b, "Μέσος όρος: ορθογώνιο με το ίδιο εμβαδόν",
                  sprintf("εμβαδόν κάτω από την h = %.4f = ȳ · (b − a) = %.4f · %.4f",
                          s$ybar * s$L, s$ybar, s$L))

  # Δεξιά: διάμεσος = η τιμή στο μέσο του διαστήματος
  # (λίγος χώρος κάτω από το 0 για τις ετικέτες των μισών του διαστήματος)
  panel(c(a, b), c(yl[1] - 0.1 * diff(yl), yl[2]), "Διάμεσος: η τιμή στο μέσο του διαστήματος", "h(x)",
        "οι μισές x δίνουν τιμές ≤ M, οι άλλες μισές ≥ M")
  segments(a, 0, s$mid, 0, lwd = 7, col = fade(col_fn, 0.6), lend = 1)
  segments(s$mid, 0, b, 0, lwd = 7, col = fade(col_dist, 0.6), lend = 1)
  segments(a, h(a), a, s$M, lwd = 7, col = fade(col_fn, 0.6), lend = 1)
  segments(a, s$M, a, h(b), lwd = 7, col = fade(col_dist, 0.6), lend = 1)
  lines(x, h(x), col = col_fn, lwd = 2)
  segments(s$mid, 0, s$mid, s$M, lty = 3, col = col_ink)
  segments(a, s$M, b, s$M, lty = 2, lwd = 2, col = col_median)
  segments(a, s$ybar, b, s$ybar, lwd = 1, col = col_mean)
  points(s$mid, s$M, pch = 19, col = col_median, cex = 1.3)
  # Ετικέτες απευθείας πάνω στα στοιχεία (η h ανεβαίνει, οπότε κάτω δεξιά από το M είναι άδειο)
  text(s$mid, s$M, sprintf("Διάμεσος M = h(%.3f) = %.4f", s$mid, s$M),
       adj = c(-0.05, 1.6), col = col_ink)
  text(a, s$ybar, sprintf("   Μέσος όρος ȳ = %.4f", s$ybar), adj = c(0, -0.6), col = col_ink, cex = 0.9)
  text((a + s$mid) / 2, 0, "[a, μέσο]", pos = 1, col = col_fn, cex = 0.9)
  text((s$mid + b) / 2, 0, "[μέσο, b]", pos = 1, col = col_dist, cex = 0.9)
  text(a, (h(a) + s$M) / 2, "τιμές ≤ M", pos = 4, col = col_fn, cex = 0.9)
  text(a, (s$M + h(b)) / 2, "τιμές ≥ M", pos = 4, col = col_dist, cex = 0.9)

  mtext(sprintf("Μέσος όρος και διάμεσος των τιμών της %s στο %s", name, interval_lab),
        outer = TRUE, line = 1, cex = 1.3, font = 2, col = col_ink)
  invisible(dev.off())
}

# ---------------------------------------------------------------------------
# Εικόνα 2: μέση απόλυτη απόκλιση σε 4 βήματα
# ---------------------------------------------------------------------------
plot_mad_steps <- function(h, a, b, name, interval_lab, file, examples = c(0.2, 0.9)) {
  s  <- describe_values(h, a, b)
  x  <- seq(a, b, length.out = 500)
  yl <- c(min(0, h(a)), h(b) * 1.12)
  d  <- function(x) h(x) - s$ybar            # προσημασμένη απόσταση
  const <- function(v) function(x) rep(v, length(x))
  xe <- a + examples * s$L                   # δύο σημεία-παραδείγματα

  png(file, width = 1500, height = 1050, res = 110)
  par(mfrow = c(2, 2), mar = c(4.2, 4.5, 3.8, 1.5), oma = c(0, 0, 3, 0), col.lab = col_muted)

  # 1. Η συνάρτηση και ο μέσος όρος
  panel(c(a, b), yl, "1. Η συνάρτηση h(x) και ο μέσος όρος ȳ", "h(x)",
        "τα βελάκια: πόσο απέχει η τιμή h(x) από το ȳ")
  lines(x, h(x), col = col_fn, lwd = 2)
  abline(h = s$ybar, col = col_ink, lty = 2)
  text(a, s$ybar, sprintf("  ȳ = %.3f", s$ybar), adj = c(0, -0.6), col = col_ink)
  arrows(xe, s$ybar, xe, h(xe), length = 0.07, col = col_dist, lwd = 1.5)
  points(xe, h(xe), pch = 19, col = col_fn, cex = 1.1)

  # 2. Η απόσταση |h(x) − ȳ|
  y2 <- max(abs(d(x)))
  panel(c(a, b), c(min(d(x)), y2 * 1.12), "2. Η απόσταση κάθε τιμής από τον μέσο: |h(x) − ȳ|",
        "|h(x) − ȳ|", "πάντα θετική: το κομμάτι κάτω από το ȳ «διπλώνει» προς τα πάνω")
  abline(h = 0, col = col_muted)
  lines(x, d(x), col = fade(col_dist, 0.35), lwd = 1.5, lty = 2)
  lines(x, abs(d(x)), col = col_dist, lwd = 2)
  segments(xe, 0, xe, abs(d(xe)), col = col_dist, lty = 3)
  points(xe, abs(d(xe)), pch = 19, col = col_dist, cex = 1.2)
  for (x0 in xe)   # ετικέτα στην πλευρά όπου δεν περνάει η καμπύλη
    text(x0, abs(d(x0)), sprintf("%.2f", abs(d(x0))),
         adj = if (x0 < s$xs) c(-0.15, -0.6) else c(1.15, -0.6), col = col_ink)
  legend("top", inset = 0.01, cex = 0.85, bg = fade("white", 0.9), box.col = col_grid,
         legend = c("|h(x) − ȳ|: η γραφική παράσταση",
                    "αχνό: h(x) − ȳ πριν την απόλυτη τιμή —",
                    "μόνο για σύγκριση"),
         col = c(col_dist, fade(col_dist, 0.35), NA),
         lty = c(1, 2, NA), lwd = c(2, 1.5, NA))

  # 3. Εμβαδόν ανάμεσα στην h(x) και στο ȳ
  panel(c(a, b), yl, sprintf("3. Εμβαδόν ανάμεσα σε h(x) και ȳ = %.4f", s$below + s$above),
        "h(x)", "τα δύο κομμάτια είναι πάντα ίσα")
  fill_between(a, s$xs, h, const(s$ybar), fade(col_area, 0.25))
  fill_between(s$xs, b, h, const(s$ybar), fade(col_area, 0.5))
  lines(x, h(x), col = col_fn, lwd = 2)
  abline(h = s$ybar, col = col_ink, lty = 2)
  # ετικέτες στο κέντρο βάρους κάθε κομματιού
  cx_b <- integrate(function(x) x * (s$ybar - h(x)), a, s$xs)$value / s$below
  cy_b <- s$ybar - integrate(function(x) (s$ybar - h(x))^2 / 2, a, s$xs)$value / s$below
  cx_a <- integrate(function(x) x * (h(x) - s$ybar), s$xs, b)$value / s$above
  cy_a <- s$ybar + integrate(function(x) (h(x) - s$ybar)^2 / 2, s$xs, b)$value / s$above
  text(cx_b, cy_b, sprintf("κάτω από ȳ\n%.4f", s$below), col = col_ink, cex = 0.9)
  text(cx_a, cy_a, sprintf("πάνω από ȳ\n%.4f", s$above), col = col_ink, cex = 0.9)

  # 4. Το αποτέλεσμα: διαιρούμε με το μήκος του διαστήματος
  panel(c(a, b), yl,
        sprintf("4. Απόκλιση = εμβαδόν ÷ (b − a) = %.4f ÷ %.4f = %.4f",
                s$below + s$above, s$L, s$mad),
        "h(x)", sprintf("ȳ ± απόκλιση = [%.3f, %.3f]", s$ybar - s$mad, s$ybar + s$mad))
  rect(a, s$ybar - s$mad, b, s$ybar + s$mad, col = fade(col_fn, 0.18), border = NA)
  lines(x, h(x), col = col_fn, lwd = 2)
  abline(h = s$ybar, col = col_ink, lty = 2)
  xa <- a + 0.06 * s$L
  arrows(xa, s$ybar, xa, s$ybar + s$mad, length = 0.08, col = col_ink, lwd = 1.5)
  arrows(xa, s$ybar, xa, s$ybar - s$mad, length = 0.08, col = col_ink, lwd = 1.5)
  text(xa, s$ybar + s$mad / 2, sprintf("%.3f", s$mad), pos = 4, col = col_ink)
  text(xa, s$ybar - s$mad / 2, sprintf("%.3f", s$mad), pos = 4, col = col_ink)

  mtext(sprintf("Μέση απόλυτη απόκλιση των τιμών της %s στο %s", name, interval_lab),
        outer = TRUE, line = 1, cex = 1.3, font = 2, col = col_ink)
  invisible(dev.off())
}

# ---------------------------------------------------------------------------
# Εκτύπωση αποτελεσμάτων και εικόνες
# ---------------------------------------------------------------------------
run <- function(h, a, b, name, interval_lab, tag, examples = c(0.2, 0.9)) {
  s <- describe_values(h, a, b)
  cat(sprintf("\n%s στο %s\n", name, interval_lab))
  cat(sprintf("  Μέσος όρος             ȳ = 1/(b−a)·∫ h(x) dx        = %.6f\n", s$ybar))
  cat(sprintf("  Διάμεσος               M = h((a+b)/2) = h(%.4f)  = %.6f\n", s$mid, s$M))
  cat(sprintf("  Εμβαδόν κάτω / πάνω από ȳ                        = %.6f / %.6f\n", s$below, s$above))
  cat(sprintf("  Μέση απόλυτη απόκλιση  1/(b−a)·∫ |h(x) − ȳ| dx    = %.6f  (απευθείας: %.6f)\n",
              s$mad, s$mad_direct))
  plot_mean_median(h, a, b, name, interval_lab, sprintf("values_mean_median_%s.png", tag))
  plot_mad_steps(h, a, b, name, interval_lab, sprintf("values_mad_steps_%s.png", tag), examples)
}

# 1. Απλή αύξουσα ευθεία: h(x) = 2x στο [0, 1]
#    Με το χέρι: ȳ = ∫_0^1 2x dx = 1,  M = h(0.5) = 1,  απόκλιση = ∫_0^1 |2x − 1| dx = 0.5
run(function(x) 2 * x, 0, 1, "h(x) = 2x", "[0, 1]", "linear")

# 2. Κυρτή (λυγίζει προς τα πάνω): h(x) = x² στο [0, 1]
#    Με το χέρι: ȳ = ∫_0^1 x² dx = 1/3,  M = h(0.5) = 1/4  →  μέσος > διάμεσος
#    Η h περνάει το ȳ στο x = 1/√3, απόκλιση = 4/(9√3)
run(function(x) x^2, 0, 1, "h(x) = x²", "[0, 1]", "convex")

# 3. Κοίλη (λυγίζει προς τα κάτω): h(x) = 2x − x² στο [0, 1]
#    Με το χέρι: ȳ = 1 − 1/3 = 2/3,  M = h(0.5) = 3/4  →  μέσος < διάμεσος
#    Είναι το x² «αναποδογυρισμένο» (h(x) = 1 − (1 − x)²), άρα ίδια απόκλιση 4/(9√3)
run(function(x) 2 * x - x^2, 0, 1, "h(x) = 2x − x²", "[0, 1]", "concave")

# Από το άθροισμα στο ολοκλήρωμα, με την x²
plot_riemann(function(x) x^2, 0, 1, "h(x) = x²", "[0, 1]", c(4, 10, 50), "values_riemann.png")

# Οι τρεις μαζί, μόνο με τον μέσο όρο
plot_means_row(list(list(h = function(x) 2 * x,       name = "h(x) = 2x"),
                    list(h = function(x) x^2,         name = "h(x) = x²"),
                    list(h = function(x) 2 * x - x^2, name = "h(x) = 2x − x²")),
               0, 1, "[0, 1]", "values_mean_all.png")

cat(sprintf("\nΈλεγχος με το χέρι:  1/3 = %.6f,  2/3 = %.6f,  4/(9√3) = %.6f\n",
            1 / 3, 2 / 3, 4 / (9 * sqrt(3))))
