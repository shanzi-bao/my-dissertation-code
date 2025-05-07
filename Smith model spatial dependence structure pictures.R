library(SpatialExtremes)  

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))  # 增加标题边距


set.seed(8)
x <- y <- seq(0, 10, length = 100)
coord <- cbind(x, y)

data1 <- rmaxstab(1, coord, "gauss", cov11 = 9/8, cov12 = 0, cov22 = 9/8, grid = TRUE)
filled.contour(x, y, log(data1), color.palette = terrain.colors, 
               main = "Isotropic\ncov11=9/8, cov12=0, cov22=9/8")

data2 <- rmaxstab(1, coord, "gauss", cov11 = 2, cov12 = 0, cov22 = 0.5, grid = TRUE)
filled.contour(x, y, log(data2), color.palette = terrain.colors,
               main = "Stronger x-dependence\ncov11=2, cov12=0, cov22=0.5")

data3 <- rmaxstab(1, coord, "gauss", cov11 = 0.5, cov12 = 0, cov22 = 2, grid = TRUE)
filled.contour(x, y, log(data3), color.palette = terrain.colors,
               main = "Stronger y-dependence\ncov11=0.5, cov12=0, cov22=2")

data4 <- rmaxstab(1, coord, "gauss", cov11 = 1, cov12 = 0.5, cov22 = 1, grid = TRUE)
filled.contour(x, y, log(data4), color.palette = terrain.colors,
               main = "With correlation\ncov11=1, cov12=0.5, cov22=1")