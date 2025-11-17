plot(c(-1,0,1,-1,0,1),c(0,0,0,1,1,1), lwd = 2, pch = 19, 
     xlim = c(-1.5,1.5), ylim = c(-0.3, 1.7), asp = 1, ylab = "", xlab = "", axes = FALSE)
axis(1, at = -2:2)
axis(2, at = -1:1)
box()

lines(c(-3,3), c(0,0), lty = 2, lwd = 0.2)
lines(c(-3,3), c(1,1), lty = 2, lwd = 0.2)
lines(c(0,0), c(-3,3), lty = 2, lwd = 0.2)
lines(c(1,1), c(-3,3), lty = 2, lwd = 0.2)
lines(c(-1,-1), c(-3,3), lty = 2, lwd = 0.2)

lines(c(-1,1), c(0,0), lwd = 3, lty = 1)
lines(c(0,1), c(0,1), lwd = 3, lty = 1)
lines(c(0,-1), c(0,1), lwd = 3, lty = 1)
lines(c(0,0), c(0,1), lwd = 3, lty = 1)

lines(c(0,4), c(0,2), lwd = 1, lty = 2)
lines(c(0,2), c(0,4), lwd = 1, lty = 2)
lines(c(0,-2), c(0,4), lwd = 1, lty = 2)
lines(c(0,-4), c(0,2), lwd = 1, lty = 2)
lines(c(-4,4), c(0,0), lwd = 1, lty = 2)

text(0, 1.3,"|m| > 2")
text(-1.6,1.25,"-2 < m < -0.5")
text(-2,0.4,"-0.5 < m < 0")
text(1.6,1.25,"0.5 < m < 2")
text(2,0.4,"0 < m < 0.5")

