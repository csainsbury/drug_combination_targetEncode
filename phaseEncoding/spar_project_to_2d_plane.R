library(data.table)
library(rgl)
library(mosaic)
library(scatterplot3d)

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
# steps to generate orthogonal projection onto plane:
# take the normal of the projection plan (1,1,1)
# project the vector (origin -> s, where s is each point) onto the normal vector
# subtract this 

vector_length <- max(x) * 10
normal_vector <- c(vector_length, vector_length, vector_length)
 
# points in 3d space
points <- data.frame(x, y, z)

#projected vectors
proj_vec <- data.frame(rep(0, length(x)), rep(0, length(x)), rep(0, length(x)))
colnames(proj_vec) <- c('proj_x', 'proj_y', 'proj_z')

for (i in c(1:length(x))) {
  # project returns the projection of x onto u
  # project(x, u)
  
  proj_vec[i, ] <- project(c(x[i], y[i], z[i]), normal_vector)
  
}

plane_projection <- points - proj_vec
plot3d(plane_projection)

# rotate plane to 2D plot
cross_prod <- function(a, b) {
  
  result = c(a[2]*b[3] - a[3]*b[2],
             a[3]*b[1] - a[1]*b[3],
             a[1]*b[2] - a[2]*b[1])
  return(result)
  
}

M = normal_vector # normal to current plane
M = c(1,1,1)
N = c(0, 0, 1) # normal to the plane to rotate into

costheta <-dot(M, N) / (norm(matrix(M))*norm(matrix(N)))

axis <- cross_prod(M, N) / norm(matrix(cross_prod(M, N)))

c = costheta
s = sqrt(1 - c*c)
C = 1-c

rotated_points <- matrix(0, nrow = nrow(plane_projection), ncol = 3)
  
for (j in c(1:nrow(plane_projection))) {
  
  x = plane_projection$x[j]
  y = plane_projection$y[j]
  z = plane_projection$z[j]
  
  matrix_components <- c(x*x*C+c, x*y*C-z*s, x*z*C+y*s, y*x*C+z*s,  y*y*C+c, y*z*C-x*s, z*x*C-y*s,  z*y*C+x*s,  z*z*C+c)
  rmat <- matrix(matrix_components, nrow = 3, ncol = 3, byrow = TRUE)
  
  v <- c(x,y,z)
  
  newpoint  = rmat %*% v
  
  rotated_points[j, ] <- newpoint
  
}

plot3d(rotated_points)

## 2D rotation
# functions in degrees
# R operates in radians by default
coss <- function(a){
  a=cos(a*pi/180)
  return(a)
}
sinn <- function(a){
  a= sin(a*pi/180)
  return(a)
}
tann <- function(a){
  a= tan(a*pi/180)
  return(a)
}

rows = 1000
p = data.frame(matrix(nrow = rows, ncol = 0))
p$x = rnorm(rows, 4, 1)
p$y = rnorm(rows, 4, 4)

plot(p$x, p$y, xlim = c(-20,20), ylim = c(-20,20))

#p = c(4, 5)


for (j in c(1:nrow(p))) {
  
  ang = 120
  
  R <- c(coss(ang), -sinn(ang), sinn(ang), coss(ang))
  rmat <- matrix(R, nrow = 2, ncol = 2, byrow = TRUE)
  
  t <- as.matrix(p[j,])
  pr <- rmat %*% as.vector(t)
  
  # pr1 <- (x * cos(j)) - (y * sin(j))
  # pr2 <- (x * sin(j)) + (y * cos(j))
  # 
  # pr <- c(pr1, pr2)
  # 
  if (j == 1) {
    plot(p$x, p$y, xlim = c(-20,20), ylim = c(-20,20))
  }
    points(pr[1], pr[2], col = 5)
}

# 
# for (j in seq(1, 360, 1)) {
#   
#   ang = j
#   
#   R <- c(coss(ang), -sinn(ang), sinn(ang), coss(ang))
#   rmat <- matrix(R, nrow = 2, ncol = 2, byrow = TRUE)
#   
#   pr <- rmat %*% matrix(p)
#   
#   # pr1 <- (x * cos(j)) - (y * sin(j))
#   # pr2 <- (x * sin(j)) + (y * cos(j))
#   # 
#   # pr <- c(pr1, pr2)
#   # 
#   if (j == 1) {
#     plot(p$x, p$y, xlim = c(-5,10), ylim = c(-15,15))
#   } else
#     {
#     points(pr[1], pr[2], col = 5)
#   }
# }
# 
# 
# 
