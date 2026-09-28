# Τυπική κανονική κατανομή N(0, 1):
# υπολογισμός μέσου όρου και διαμέσου με ολοκληρώματα.

# Συνάρτηση πυκνότητας πιθανότητας (γραμμένη ρητά αντί για dnorm)
f <- function(x) (1 / sqrt(2 * pi)) * exp(-x^2 / 2)

# 1. Έλεγχος ότι είναι έγκυρη πυκνότητα: ∫ f(x) dx = 1
total <- integrate(f, -Inf, Inf)$value
cat(sprintf("Ολικό εμβαδόν  ∫ f(x) dx      = %.6f\n", total))

# 2. Μέσος όρος: μ = ∫ x·f(x) dx
mu <- integrate(function(x) x * f(x), -Inf, Inf)$value
cat(sprintf("Μέσος όρος     μ = ∫ x·f(x) dx = %.6f\n", mu))

# 3. Διάμεσος: η τιμή m για την οποία F(m) = ∫_{-∞}^{m} f(x) dx = 0.5
F <- function(m) integrate(f, -Inf, m)$value
m <- uniroot(function(m) F(m) - 0.5, interval = c(-5, 5), tol = 1e-10)$root
cat(sprintf("Διάμεσος       m : F(m) = 0.5  → m = %.6f\n", m))

# 4. Μέση απόλυτη απόκλιση: ο μέσος όρος της νέας ποσότητας |x − μ|,
#    με την ίδια πυκνότητα f ως βάρος: E|X − μ| = ∫ |x − μ|·f(x) dx
mad <- integrate(function(x) abs(x - mu) * f(x), -Inf, Inf, rel.tol = 1e-10)$value
cat(sprintf("Μέση απόλ. απόκλιση ∫ |x − μ|·f(x) dx = %.6f  (τύπος √(2/π) = %.6f)\n",
            mad, sqrt(2 / pi)))

# Χρώματα (ίδια με τα άλλα scripts)
col_ink    <- "#0b0b0b"
col_muted  <- "#898781"
col_grid   <- "#e1e0d9"
col_fn     <- "#2a78d6"   # πυκνότητα, μισό εμβαδόν αριστερά του διαμέσου
col_dist   <- "#eb6834"   # μισό εμβαδόν δεξιά του διαμέσου
col_mean   <- "#e34948"
col_median <- "#008300"
fade <- function(col, a) adjustcolor(col, alpha.f = a)

# Γράφημα: τα δύο μισά του εμβαδού, ο διάμεσος που τα χωρίζει και ο μέσος όρος
png("normal_mean_median.png", width = 1200, height = 650, res = 110)
par(mar = c(4.2, 4.5, 3.8, 1.5), col.lab = col_muted)
x    <- seq(-4, 4, length.out = 500)
ymax <- max(f(x))
plot(NA, xlim = c(-4, 4), ylim = c(0, ymax * 1.15), axes = FALSE, xlab = "x", ylab = "f(x)", main = "")
title(main = "Τυπική κανονική κατανομή N(0, 1): μέσος όρος και διάμεσος",
      adj = 0, line = 1.8, font.main = 1, cex.main = 1.25, col.main = col_ink)
mtext("διάμεσος: χωρίζει το εμβαδόν σε δύο μισά · μέσος όρος: το σημείο ισορροπίας",
      side = 3, line = 0.2, adj = 0, cex = 0.9, col = col_muted)
abline(h = axTicks(2), col = col_grid)
axis(1, col = col_muted, col.axis = col_muted)
axis(2, las = 1, col = col_muted, col.axis = col_muted)
xl <- seq(-4, m, length.out = 300)
xr <- seq(m, 4, length.out = 300)
polygon(c(-4, xl, m), c(0, f(xl), 0), col = fade(col_fn, 0.3), border = NA)
polygon(c(m, xr, 4), c(0, f(xr), 0), col = fade(col_dist, 0.3), border = NA)
lines(x, f(x), col = col_fn, lwd = 2)
text(-1.2, 0.3 * f(-1.2), sprintf("εμβαδόν\n%.2f", F(m)), col = col_ink)
text(1.2, 0.3 * f(1.2), sprintf("εμβαδόν\n%.2f", 1 - F(m)), col = col_ink)
# Η κόκκινη πρώτα, ώστε η πράσινη διακεκομμένη να φαίνεται πάνω της όταν συμπίπτουν
segments(mu, 0, mu, ymax * 1.05, col = col_mean, lwd = 2)
segments(m, 0, m, ymax * 1.05, col = col_median, lwd = 2, lty = 2)
if (abs(mu - m) < 1e-6) {
  text(mu, ymax * 1.1, sprintf("μέσος μ = διάμεσος m = %.4f", abs(mu)), col = col_ink)
} else {
  # ο μεγαλύτερος από τους δύο δεξιά της γραμμής του, ο μικρότερος αριστερά
  text(m,  ymax * 1.1, sprintf("διάμεσος m = %.4f ", m), adj = if (m <= mu) 1 else 0, col = col_median)
  text(mu, ymax * 1.1, sprintf(" μέσος μ = %.4f", mu),   adj = if (m <= mu) 0 else 1, col = col_mean)
}
invisible(dev.off())
cat("Το γράφημα αποθηκεύτηκε στο normal_mean_median.png\n")
