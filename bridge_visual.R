# Η μετάβαση από τις απλές συναρτήσεις στις πυκνότητες, σε τρία πάνελ στο [0, 6]:
#   1. ίσα βάρη w(x) = 1/6              (ομοιόμορφη κατανομή: όπως στις απλές συναρτήσεις)
#   2. άνισα βάρη, διακριτά wᵢ          (σταθμισμένος μέσος όρος, μπάρες)
#   3. άνισα βάρη, συνεχή f(x)          (πυκνότητα f(x) = (6 − x)/18)
# Σε όλα: η τιμή είναι η θέση στον άξονα x και το ύψος είναι το βάρος της.
# Ο μέσος όρος είναι το σημείο ισορροπίας (▲). Με περισσότερο βάρος αριστερά,
# το ▲ μετακινείται από το 3 προς το 2.

# Χρώματα
col_ink   <- "#0b0b0b"
col_muted <- "#898781"
col_grid  <- "#e1e0d9"
col_w     <- "#2a78d6"   # βάρος
col_mean  <- "#e34948"   # μέσος όρος
fade <- function(col, a) adjustcolor(col, alpha.f = a)

a <- 0; b <- 6
xlim  <- c(-0.4, 6.4)
y_max <- 0.35   # ίδια κλίμακα βάρους και στα τρία πάνελ

# 1. Ίσα βάρη (ομοιόμορφη κατανομή)
w_flat  <- function(x) ifelse(x >= a & x <= b, 1 / (b - a), 0)
mu_flat <- integrate(function(x) x * w_flat(x), a, b)$value

# 3. Πυκνότητα που φθίνει: f(x) = (6 − x)/18, εμβαδόν = 6 · (1/3) / 2 = 1
f_dec  <- function(x) ifelse(x >= a & x <= b, (b - x) / 18, 0)
mu_dec <- integrate(function(x) x * f_dec(x), a, b)$value

# 2. Διακριτά βάρη: τα μέσα των διαστημάτων πλάτους 1, με βάρος f(xᵢ) · 1
#    (δηλαδή 11/36, 9/36, ..., 1/36 — άθροισμα 1)
xi <- seq(0.5, 5.5, by = 1)
ni <- 12 - 2 * xi            # αριθμητές: 11, 9, 7, 5, 3, 1
wi <- ni / 36
mu_bars <- sum(xi * wi)

cat(sprintf("Συνολικό βάρος:  %.4f  %.4f  %.4f\n",
            integrate(w_flat, a, b)$value, sum(wi), integrate(f_dec, a, b)$value))
cat(sprintf("Μέσος όρος:      %.4f  %.4f  %.4f\n", mu_flat, mu_bars, mu_dec))

# Κοινό πλαίσιο
panel <- function(title, note, formula, mu, ghost = TRUE) {
  plot(NA, xlim = xlim, ylim = c(-0.12 * y_max, y_max * 1.2), axes = FALSE,
       xlab = "x  (η τιμή)", ylab = "βάρος", main = "")
  title(main = title, adj = 0, line = 1.8, font.main = 1, cex.main = 1.25, col.main = col_ink)
  mtext(note, side = 3, line = 0.2, adj = 0, cex = 0.9, col = col_muted)
  ticks <- seq(0, 0.3, by = 0.1)
  abline(h = ticks, col = col_grid)
  axis(1, at = 0:6, col = col_muted, col.axis = col_muted)
  axis(2, at = ticks, las = 1, col = col_muted, col.axis = col_muted)
  text(xlim[2], y_max * 1.12, formula, adj = 1, col = col_ink)
  if (ghost) {   # πού ήταν ο μέσος όρος με ίσα βάρη
    points(mu_flat, -0.06 * y_max, pch = 2, cex = 2, col = col_muted)
    text(mu_flat, -0.06 * y_max, "   ίσα βάρη: 3", adj = 0, col = col_muted, cex = 0.85)
  }
  points(mu, -0.06 * y_max, pch = 17, cex = 2.2, col = col_mean)   # σημείο ισορροπίας
  text(mu, -0.06 * y_max, sprintf("μ = %.2f   ", mu), adj = 1, col = col_ink)
}

png("bridge_weights.png", width = 1800, height = 680, res = 110)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.8, 1.5), oma = c(2, 0, 3, 0), col.lab = col_muted)
par(cex = 1)   # το mfrow με 3 πάνελ μικραίνει αυτόματα τα γράμματα

# 1. Ίσα βάρη
y1 <- 1 / (b - a)
panel("1. Ίσα βάρη (ομοιόμορφη)", "όπως στις απλές συναρτήσεις · συνολικό βάρος = 1",
      "μ = ∫ x · 1/(b − a) dx", mu_flat, ghost = FALSE)
rect(a, 0, b, y1, col = fade(col_w, 0.18), border = NA)
segments(c(xlim[1], a, b), c(0, y1, 0), c(a, b, xlim[2]), c(0, y1, 0), col = col_w, lwd = 2)
segments(c(a, b), 0, c(a, b), y1, col = col_w, lwd = 2)
text(1.5, y1, "κάθε x: βάρος 1/6", pos = 3, col = col_ink)

# 2. Διακριτά άνισα βάρη
panel("2. Άνισα βάρη, διακριτά", "σταθμισμένος μέσος όρος · Σ wᵢ = 1",
      "μ = Σ xᵢ · wᵢ", mu_bars)
rect(xi - 0.4, 0, xi + 0.4, wi, col = fade(col_w, 0.6), border = NA)
text(xi, wi, sprintf("%d/36", ni), pos = 3, cex = 0.85, col = col_ink)

# 3. Συνεχή άνισα βάρη: πυκνότητα
x <- sort(c(seq(xlim[1], xlim[2], length.out = 500), a, b))
panel("3. Άνισα βάρη, συνεχή: πυκνότητα", "f(x) = (6 − x)/18 · συνολικό εμβαδόν = 1",
      "μ = ∫ x · f(x) dx", mu_dec)
polygon(c(x[1], x, x[length(x)]), c(0, f_dec(x), 0), col = fade(col_w, 0.18), border = NA)
lines(x, f_dec(x), col = col_w, lwd = 2)

mtext("Από τα ίσα βάρη στην πυκνότητα: το βάρος μετακινεί τον μέσο όρο",
      outer = TRUE, line = 1, cex = 1.3, font = 2, col = col_ink)
mtext("Σε όλα τα πάνελ: η τιμή είναι η θέση στον άξονα x, το ύψος είναι το βάρος της.  ▲ = μέσος όρος (σημείο ισορροπίας)",
      side = 1, outer = TRUE, line = 0.5, cex = 0.95, col = col_muted)
invisible(dev.off())
cat("Το γράφημα αποθηκεύτηκε στο bridge_weights.png\n")
