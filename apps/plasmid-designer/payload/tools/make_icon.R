#!/usr/bin/env Rscript
# make_icon.R — plasmid/vector-map icon for the "Plasmid Designer" .app bundle.
NCOL <- 1024L
SKY <- '#0284c7'
EMERA <- '#047857'
AMBER <- '#f59e0b'
VIOLA <- '#7c3aed'
INK <- '#1e293b'
SLATE <- '#94a3b8'
BG <- '#ffffff'

dir.create('iconset', showWarnings = FALSE)
master <- 'iconset/plasmid_master.png'

draw_feature_arc <- function(x, y, r, a1, a2, col, lw = 5) {
  pts <- seq(a1, a2, length.out = 40L)
  os <- sapply(pts, function(a) c(x + r * cos(a), y + r * sin(a)))
  lines(os[1, ], os[2, ], col = col, lwd = lw)
  ea <- a2
  ax <- x + r * cos(ea)
  ay <- y + r * sin(ea)
  head_ang <- seq(ea - 0.15, ea + 0.15, length.out = 10L)
  hx <- sapply(head_ang, function(a) x + r * cos(a))
  hy <- sapply(head_ang, function(a) y + r * sin(a))
  lines(c(ax, hx[2], hx[10]), c(ay, hy[2], hy[10]), col = col, lwd = lw)
  lines(c(ax, hx[1], hx[9]), c(ay, hy[3], hy[9]), col = col, lwd = lw)
}

p2 <- 0.96
polygon(c(1 - p2, 1, p2, 0, 0, p2, 1, 1 - p2), c(0, p2, p2, 0, 0, p2, p2, 1), col = BG, border = SLATE, lwd = 1.2)

ring <- function(cx, cy, R, n, col = SLATE) {
  for (i in seq_len(n)) {
    i1 <- 2 * pi * (i - 1) / n
    i2 <- 2 * pi * i / n
    draw_feature_arc(cx, cy, R, i1, i2, col, lw = 4)
  }
}
ring(cx, cy, R_out <- NCOL / 2 - 14, 24L, SLATE)
ring(cx, cy, R_out * 0.90, 24L, SLATE)

for (i in seq_len(4L)) {
  j <- c(SKY, EMERA, AMBER, VIOLA)[i]
  draw_feature_arc(cx, cy, R_out * 0.68, (i - 1) * pi / 2, i * pi / 2, rep(j, 1), lw = 5)
}

draw_feature_arc(cx, cy, R_out * 0.64, 5 * pi / 6, pi / 6, AMBER, lw = 4)
ea <- pi / 6
ax <- cx + R_out * 0.58 * cos(ea)
ay <- cy + R_out * 0.58 * sin(ea)
head_ang <- seq(ea - 0.14, ea + 0.14, length.out = 10L)
hx <- sapply(head_ang, function(a) cx + R_out * 0.58 * cos(a))
hy <- sapply(head_ang, function(a) cy + R_out * 0.58 * sin(a))
lines(c(ax, hx[2], hx[10]), c(ay, hy[2], hy[10]), col = AMBER, lwd = 4)
lines(c(ax, hx[1], hx[9]), c(ay, hy[3], hy[9]), col = AMBER, lwd = 4)

for (i in seq_len(12L)) {
  a <- 2 * pi * (i - 1) / 12
  lines(cx + R_out * 0.74 * cos(a), cy + R_out * 0.74 * sin(a), cx + R_out * 0.83 * cos(a), cy + R_out * 0.83 * sin(a), col = SLATE, lwd = 1.1)
}
for (i in seq_len(4L)) {
  a <- 2 * pi * (i - 1) / 4
  lines(cx + R_out * 0.74 * cos(a), cy + R_out * 0.74 * sin(a), cx + R_out * 0.88 * cos(a), cy + R_out * 0.88 * sin(a), col = SLATE, lwd = 1.1)
}

text(cx, cy - R_out * 0.10, 'pJAT', cex = 1.6, col = INK, font = 2)
text(cx, cy - R_out * 0.26, 'donor payload', cex = 0.75, col = INK)
dev.off()

for (s in c('16 16', '32 32', '64 64', '128 128', '256 256', '512 512')) system2('sips', c('-z', s, master, '-s', 'format', 'png'), stdout = TRUE, stderr = TRUE)
file.rename(master, 'iconset/icon_512x512@2x.png')
icns <- 'iconset/icon.icns'
if (file.exists(icns)) file.remove(icns)
system2('iconutil', c('-c', 'icns', 'iconset'), stdout = TRUE, stderr = TRUE)
if (file.exists(icns)) message('Wrote iconset/icon.icns (', file.info(icns)$size, ' bytes)')

par(mar = c(0, 0, 0, 0), xaxs = 'i', yaxs = 'i', fin = c(1, 1))
plot(0, 0, type = 'n', asp = 1, axes = FALSE, xlab = '', ylab = '',
     bty = 'n', col = BG, fg = BG, xlim = c(0, 1), ylim = c(0, 1))
