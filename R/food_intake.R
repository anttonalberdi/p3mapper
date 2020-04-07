#' @title Food intake estimation
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords mass  energy budget
#' @description Computes pest consumption estimations based on body mass of predators and average energy content of predators.
#' @param avgmass Average body mass of the predator (grams).
#' @param sdmass Standard deviation of the body mass of the predator (grams).
#' @param avgenergy Average energy content of the prey (kJ/g).
#' @param sdenergy Standard deviation of the energy content of the prey (kJ/g).
#' @param constanta Constant a from the allometric equations by Nagy et al. 1999.
#' @param constantb Constant b from the allometric equations by Nagy et al. 1999.
#' @param iterations Number of iterations (default 100).
#' @usage food_intake(avgmass,minmass,maxmass,avgenergy,minenergy,maxenergy,iterations)
#' @return A value (if one iteration) or vector (if multiple iterations) of prey consumption estimations.
#' @import raster gstat rgdal
#' @examples
#' food_intake()
#' food_intake()
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' XXXXXX
#' @export

food_intake <- function(avgmass,sdmass,avgenergy,sdenergy,constanta,constantb,iterations){

if(missing(iterations)){iterations=100}

#Generate mass Gaussian distribution
mass_dist <- rnorm(iterations, avgmass, sdmass)
mass_dist[mass_dist < 0.0001] <- 0
#Generate energy Gaussian distribution
energy_dist <- rnorm(iterations, avgenergy, sdenergy)

#Calculate wet arthropod mass consumption
intake_vector <- c()
for (i in c(1:iterations)){
  mass <- mass_dist[i]
  energy <- energy_dist[i]

  #Apply energy budget equation (Nagy 1999)
  budget <- constanta * mass^constantb

  #Obtain arthropod mass consumption
  intake <- budget / energy

  #Add value to vector
  intake_vector <- c(intake_vector,intake)
}

return(intake_vector)
}
