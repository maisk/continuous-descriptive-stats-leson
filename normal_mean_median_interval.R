# Τυπική κανονική κατανομή N(0, 1) περιορισμένη σε διάστημα [a, b]:
# υπολογισμός μέσου όρου και διαμέσου με ολοκληρώματα.

# Το διάστημα (άλλαξέ το όπως θέλεις)
a <- -0.5
b <- 2

# Συνάρτηση πυκνότητας πιθανότητας της N(0, 1)
f <- function(x) (1 / sqrt(2 * pi)) * exp(-x^2 / 2)

# 1. Εμβαδόν κάτω από την καμπύλη στο [a, b]: Z = ∫_a^b f(x) dx
#    Διαιρώντας με το Z, η f γίνεται πυκνότητα με ολικό εμβαδόν 1 στο [a, b].
Z <- integrate(f, a, b)$value
g <- function(x) f(x) / Z
cat(sprintf("Διάστημα       [a, b] = [%g, %g]\n", a, b))
cat(sprintf("Εμβαδόν        Z = ∫_a^b f(x) dx   = %.6f\n", Z))
cat(sprintf("Έλεγχος        ∫_a^b g(x) dx       = %.6f\n", integrate(g, a, b)$value))

# 2. Μέσος όρος: μ = ∫_a^b x·g(x) dx
mu <- integrate(function(x) x * g(x), a, b)$value
cat(sprintf("Μέσος όρος     μ = ∫_a^b x·g(x) dx = %.6f\n", mu))

# 3. Διάμεσος: η τιμή m στο [a, b] για την οποία ∫_a^m g(x) dx = 0.5
G <- function(m) integrate(g, a, m)$value
m <- uniroot(function(m) G(m) - 0.5, interval = c(a, b), tol = 1e-10)$root
cat(sprintf("Διάμεσος       m : ∫_a^m g = 0.5  → m = %.6f\n", m))

# 4. Μέση απόλυτη απόκλιση: ο μέσος όρος της νέας ποσότητας |x − μ|,
#    με την ίδια πυκνότητα g ως βάρος: E|X − μ| = ∫_a^b |x − μ|·g(x) dx
mad <- integrate(function(x) abs(x - mu) * g(x), a, b, rel.tol = 1e-10)$value
cat(sprintf("Μέση απόλ. απόκλιση ∫_a^b |x − μ|·g(x) dx = %.6f\n", mad))

# Επαλήθευση με τους γνωστούς κλειστούς τύπους (dnorm/pnorm/qnorm)
mu_exact  <- (dnorm(a) - dnorm(b)) / (pnorm(b) - pnorm(a))
m_exact   <- qnorm((pnorm(a) + pnorm(b)) / 2)
mad_exact <- 2 * (mu * (pnorm(mu) - pnorm(a)) - (dnorm(a) - dnorm(mu))) / Z
cat(sprintf("\nΕπαλήθευση:    μ (τύπος) = %.6f,  m (τύπος) = %.6f,  απόλ. απόκλιση (τύπος) = %.6f\n",
            mu_exact, m_exact, mad_exact))

# Χρώματα (ίδια με τα άλλα scripts)
col_ink    <- "#0b0b0b"
col_muted  <- "#898781"
col_grid   <- "#e1e0d9"
col_fn     <- "#2a78d6"   # πυκνότητα, μισό εμβαδόν αριστερά του διαμέσου
col_dist   <- "#eb6834"   # μισό εμβαδόν δεξιά του διαμέσου
col_mean   <- "#e34948"
col_median <- "#008300"
fade <- function(col, a) adjustcolor(col, alpha.f = a)

# Γράφημα: η κανονικοποιημένη g(x) στο [a, b], τα δύο μισά του εμβαδού,
# ο διάμεσος που τα χωρίζει και ο μέσος όρος
png("normal_mean_median_interval.png", width = 1200, height = 650, res = 110)
par(mar = c(4.2, 4.5, 3.8, 1.5), col.lab = col_muted)
xlim <- c(-1.5, 3)
g0   <- function(x) ifelse(x >= a & x <= b, g(x), 0)   # 0 έξω από το [a, b]
x    <- sort(c(seq(xlim[1], xlim[2], length.out = 500), a, b))
ymax <- max(g0(x))
plot(NA, xlim = xlim, ylim = c(0, ymax * 1.15), axes = FALSE, xlab = "x", ylab = "g(x)", main = "")
title(main = sprintf("Η N(0, 1) στο [%g, %g]: μέσος όρος και διάμεσος", a, b),
      adj = 0, line = 1.8, font.main = 1, cex.main = 1.25, col.main = col_ink)
mtext("g(x) = f(x)/Z · τα δύο χρωματιστά εμβαδά είναι ίσα, όμως ο μέσος όρος είναι δεξιότερα",
      side = 3, line = 0.2, adj = 0, cex = 0.9, col = col_muted)
abline(h = axTicks(2), col = col_grid)
axis(1, col = col_muted, col.axis = col_muted)
axis(2, las = 1, col = col_muted, col.axis = col_muted)
xl <- seq(a, m, length.out = 300)
xr <- seq(m, b, length.out = 300)
polygon(c(a, xl, m), c(0, g(xl), 0), col = fade(col_fn, 0.3), border = NA)
polygon(c(m, xr, b), c(0, g(xr), 0), col = fade(col_dist, 0.3), border = NA)
lines(x, g0(x), col = col_fn, lwd = 2)
text((a + m) / 2, 0.3 * g((a + m) / 2), sprintf("εμβαδόν\n%.2f", G(m)), col = col_ink)
text((m + b) / 2, 0.3 * g((m + b) / 2), sprintf("εμβαδόν\n%.2f", 1 - G(m)), col = col_ink)
segments(m, 0, m, ymax * 1.05, col = col_median, lwd = 2, lty = 2)
segments(mu, 0, mu, ymax * 1.05, col = col_mean, lwd = 2)
# Ετικέτες: ο μεγαλύτερος από τους δύο δεξιά της γραμμής του, ο μικρότερος αριστερά
text(m,  ymax * 1.1, sprintf("διάμεσος m = %.4f ", m), adj = if (m <= mu) 1 else 0, col = col_median)
text(mu, ymax * 1.1, sprintf(" μέσος μ = %.4f", mu),   adj = if (m <= mu) 0 else 1, col = col_mean)
invisible(dev.off())
cat("Το γράφημα αποθηκεύτηκε στο normal_mean_median_interval.png\n")
