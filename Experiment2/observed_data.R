library(SpatialExtremes)


data(rainfall)
#rain
coord <- coord[,-3]
sim_coords <- coord
### fit Smith's with SpatialExtremes package 
loc.form <- loc ~ lon + lat 
scale.form <- scale ~ lon + lat 
shape.form <- shape ~ 1
mod.spe <- fitmaxstab(rain,coord=coord, cov.mod="gauss",loc.form, scale.form,shape.form, fit.marge=TRUE)
hat.param = mod.spe$fitted.values
hat.param
# set the correlation true value as zero
hat.param[2] = 0
hat.param
set.seed(6)


real_cov11 <- hat.param[1]
real_cov12 <- hat.param[2]
real_cov22 <- hat.param[3]
real_beta_loc <- hat.param[4:6]
real_beta_scale <- hat.param[7:9]
real_shape <- hat.param[10]

real_beta_loc <- as.numeric(real_beta_loc)
real_beta_scale <- as.numeric(real_beta_scale) 
real_cov11 <- as.numeric(real_cov11)
real_cov12 <- as.numeric(real_cov12)
real_cov22 <- as.numeric(real_cov22)


simData_spExt <- function(n, K, coord, cov11, cov12, cov22, beta.loc, beta.scale, shape)
{
  locs <- as.double(cbind(rep(1, K), coord)%*%beta.loc)
  scales <- as.double(cbind(rep(1, K), coord)%*%beta.scale)
  sim.data <- rmaxstab(n, coord, "gauss", grid = F, cov11=cov11, cov12=cov12, cov22=cov22)
  for(j in 1:n){
    sim.data[j,] = frech2gev(sim.data[j,], loc=locs, scale=scales, shape=shape) 
  }
  return(sim.data)
}
n_observations <- 47
observed_data <- simData_spExt(
  n = n_observations, 
  K = nrow(sim_coords), 
  coord = sim_coords,
  cov11 = real_cov11, 
  cov12 = real_cov12, 
  cov22 = real_cov22,
  beta.loc = real_beta_loc, 
  beta.scale = real_beta_scale, 
  shape = real_shape
)


print(head(observed_data))
dim(observed_data)


colnames(sim_coords) <- c("lon", "lat")

loc.form1 <- loc ~ lon + lat 
scale.form1 <- scale ~ lon + lat 
shape.form1 <- shape ~ 1
mod.spe1 <- fitmaxstab(observed_data, coord = sim_coords, cov.mod = "gauss", loc.form1, scale.form1, shape.form1, fit.marge = TRUE, locCoeff1 = 20.65249877, locCoeff2 = 0.06445933, locCoeff3 = -0.15529697, scaleCoeff1 = 3.54050776, scaleCoeff2 = 0.02308246, scaleCoeff3 = -0.03934086, shapeCoeff1 = 0.19174323)
#mod.spe1 <- fitmaxstab(observed_data, coord = sim_coords, cov.mod = "gauss", loc.form1, scale.form1, shape.form1, fit.marge = TRUE)
#mod.spe1 <- fitmaxstab(observed_data, coord = sim_coords, cov.mod = "gauss", loc.form1, scale.form1, shape.form1, fit.marge = TRUE,cov22 = 184.62, locCoeff1 = 20.65249877, locCoeff2 = 0.06445933, locCoeff3 = -0.15529697, scaleCoeff1 = 3.54050776, scaleCoeff2 = 0.02308246, scaleCoeff3 = -0.03934086, shapeCoeff1 = 0.19174323)
#mod.spe1 <- fitmaxstab(observed_data, coord = sim_coords, cov.mod = "gauss", loc.form1, scale.form1, shape.form1, fit.marge = TRUE,  scaleCoeff1 = 3.54050776, scaleCoeff2 = 0.02308246, scaleCoeff3 = -0.03934086, shapeCoeff1 = 0.19174323)




full_mean <- numeric(10)


full_mean[1:3] <- mod.spe1$fitted.values  # [304.25028, 84.00783, 154.94213]


full_mean[4] <- 20.65   
full_mean[5] <- 0.064
full_mean[6] <- -0.155
full_mean[7] <- 3.54
full_mean[8] <- 0.023
full_mean[9] <- -0.039
full_mean[10] <- 0.192


full_vcov <- matrix(0, 10, 10)


full_vcov[1:3, 1:3] <- 2.5 * mod.spe1$var.cov  


is.param <- list(mean = full_mean, vcov = full_vcov)



