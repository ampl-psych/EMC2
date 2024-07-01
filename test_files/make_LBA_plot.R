rm(list = ls())

get_drift_angle <- function(x0, x1, y0, y1){
  width <- x1 - x0
  height <- y1 - y0
  w_plot <- par("pin")[1] / diff(par("usr")[1:2])
  h_plot <- par("pin")[2] / diff(par("usr")[3:4])
  return(atan((h_plot/w_plot) * (height/width)) * 180/pi)
}

par(mfrow =c(2,2))

B <- 2
v <- 3
max <- 2
rt <- 1.5
t0 <- .3
x <- c(0, t0)
y <- c(0, 0)
ylim <- c(-.1, B+.15)
plot(x, y, type = "l", xlim = c(0, max), xlab = "", ylab = "", axes = F, ylim = ylim, lwd = 3)
arrows(x0 = 0, x1 = max, y0 = -.1, y1 = -.1, lwd = 3, length = .1)
arrows(x0 = t0, x1 = rt, y0 = 0, y1 = B-.05, lwd = 3, length = .1, col = "purple")
arrows(x0 = t0, x1 = rt, y0 = 0, y1 = B - B/3, lwd = 3, length = .1, col = "darkgreen")
lines(c(0, max), c(B, B), lwd = 3, lty = 2)
title(xlab="RT", line=0, cex.lab = 1.5)
text(t0/2, .15, "t0", cex = 1.5)
text(t0/2, B + .15, "B", cex = 1.5)
# debug(get_drift_angle)
text((t0 + rt)/1.5 - .1, B/1.5, "test_word", cex = 1.5, srt = get_drift_angle(t0, rt, 0, B-.05))
text((t0 + rt)/1.5 - .1, (B - B/3)/1.5, "test_word", cex = 1.5, srt = get_drift_angle(t0, rt, 0, B-B/3))

