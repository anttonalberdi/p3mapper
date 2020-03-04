#' @title Predator population density estimation
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords XXXXX
#' @description XXXXX
#' @param density Matrix or data frame with average, minimum and maximum population density estimations, and the corresponding geographical coordinates.
#' @param enm Distribution map shapefile of the predator.
#' @param distribution Geographical extension of the analysis.
#' @param iterations Geographical resolution if the analysis.
#' @usage predator_density(density,enm,distribution,iterations)
#' @return A raster (if one iteration) or rasterstack (if multiple iterations) object.
#' @importFrom XXXXX
#' @examples
#' predator_density()
#' predator_density()
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' XXXXXX
#' @export

predator_density <- function(density,enm,distribution,iterations){

if(missing(iterations)){iterations=100}

#Sample from density table and create density matrix
density <- as.data.frame(density)
rownames(density) <- paste("S",c(1:nrow(density)),sep="")

density_matrix <- c()
for(r in c(1:nrow(density))){
  row <- density[r,]
  density_iter <- rnorm(iterations, row$Average, row$SD)
  density_matrix <- rbind(density_matrix,density_iter)
}
rownames(density_matrix) <- rownames(density)
colnames(density_matrix) <- paste("I",c(1:ncol(density_matrix)),sep="")

#Convert data to Spatial Points
coordinates(density) <- ~x + y

#Create reference Spatial Points
reference_coord <- as(enm, "SpatialPixelsDataFrame")

#Iterate calculations
for(i in c(1:iterations)){

  cat("Iteration",i,"\n")

  #Interpolate densities data
  density_interpol <- raster(idw(formula = density_matrix[,i] ~ 1, locations = density, newdata = reference_coord, idp = 1, debug.level = 0))

  #Crop by distribution
  density_interpol <- density_interpol * distribution
  enm <- enm * distribution

  #Refine density distribution using enm
  enm_density_temp <- enm * density_interpol
  totalpop <- cellStats(density_interpol, stat = "sum", na.rm = TRUE)
  total_enm_density <- cellStats(enm_density_temp, stat = "sum", na.rm = TRUE)
  enm_density <- enm_density_temp * totalpop / total_enm_density
  enm_density[enm_density < 0] <- 0 

  names(enm_density) <- paste("Iter",i,sep="")

  if(i == 1){
  enm_density_stack <- enm_density
  }

  if(i > 1){
  enm_density_stack <- addLayer(enm_density_stack, enm_density)
  }

}

return(enm_density_stack)

}
