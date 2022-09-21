library(data.table)
library(rgl)


s <- seq(1, 100, 0.01)
s = sqrt(s^3)
s <- sin(s)
plot(s)

start_point <- 1
offset <- 33

initial_point <- start_point + (2 * offset)
end_padding <- length(s) - initial_point

x <- s[start_point:end_padding]
y <- s[(start_point + offset):(end_padding+offset)]
z <- s[(start_point + (offset * 2)):(end_padding+(offset * 2))]


plot(x, y, cex = 0.2)
## dot product a.b = a %*%

a = c(2, 3, 4)
b = c(5, 6, 7)

c = a %*% b

library(rgl)
plot3d(x,y,z)


##

s <- read.csv('~/projects11/3d2d/samples.csv')
colnames(s) <- c('t', 'mv')
s = s[2:nrow(s), ]


s <- read.csv('~/Documents/data/phaseEncodingData/export_1.csv')
#s <- read.csv('~/projects11/3d2d/gluej_1.csv')
#s <- read.csv('~/projects11/3d2d/gluek.csv')
s <- s[, 2][1:1000]
# for (j in c(1:length(s))) {
#   
#   if (s[j] < 4) {
#     s[j]=s[j]^-2
#   } else {
#     s[j] = s[j]^2
#   }
#   
# }

plot(s, cex = 0)
lines(s)


start_point <- 1
offset <- 32

initial_point <- start_point + (2 * offset)
end_padding <- length(s) - initial_point

x <- s[start_point:end_padding]
y <- s[(start_point + offset):(end_padding+offset)]
z <- s[(start_point + (offset * 2)):(end_padding+(offset * 2))]

scatterplot3d(x,y,z)

plot3d(x,y,z,
       type = 'p',
       col = ifelse(x < 3, 'red', 'black'),
       size = 2
)

plot3d(x,y,z,
       type = 'l',
       lwd = 0.8,
       col = rgb(0,0,0,0.1, maxColorValue = 1),
       size=2,
       decorate = T
       )

## glu

s <- c(rep(5.2, 10), rep(8.5, 12), rep(10.4, 6), rep(11.5, 20), rep(5.3, 10))
s <- c(rep(10.9, 10), rep(8.4, 12), rep(2.6, 6), rep(5.2, 20), rep(11.1, 10))
plot(s)

start_point <- 1
offset <- 20


initial_point <- start_point + (2 * offset)
end_padding <- length(s) - initial_point

x <- s[start_point:end_padding]
y <- s[(start_point + offset):(end_padding+offset)]
z <- s[(start_point + (offset * 2)):(end_padding+(offset * 2))]


plot3d(x,y,z,
       type = 'p',
       col = ifelse(x < 3, 'red', 'black'),
       size = 5
)

plot3d(x,y,z,
       type = 'l',
       lwd = 0.8,
       col = rgb(0,0,0,0.1, maxColorValue = 1),
       size=5,
       decorate = F
)

##

for (i in c(0:1000)) {
  
  Xc <- 1000
  Yc <- 1000
  Zc <- 1000
  
  X = x
  Y = y
  Z = z
  
  F = 50
  
  x_dev <- ((X - Xc) * (F/Z)) + Xc
  y_dev <-  ((Y - Yc) * (F/Z)) + Yc
  
  png(paste0('./Documents/data/plots/camera_', i,'.png' ))
  plot(x_dev, y_dev, main = i)
  dev.off()
  
}

