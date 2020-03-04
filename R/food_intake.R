#' @title Food intake estimation
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords mass efficiency energy budget
#' @description Computes pest consumption estimations based on body mass of predators, digestion efficincy and average energy content of predators.
#' @param avgmass Average body mass of the predator (grams).
#' @param sdmass Standard deviation of the body mass of the predator (grams).
#' @param maxmass Maximum body mass (top 5% percentile) of the predator (grams).
#' @param efficiency Digestion efficiency of the predator (percentage).
#' @param avgenergy Average energy content of the prey (kJ/g).
#' @param sdenergy Standard deviation of the energy content of the prey (kJ/g).
#' @param iterations Number of iterations (default 100).
#' @usage prey_consumption(avgmass,minmass,maxmass,efficiency,avgenergy,minenergy,maxenergy,iterations)
#' @return A value (if one iteration) or vector (if multiple iterations) of prey consumption estimations.
#' @import raster gstat rgdal
#' @examples
#' food_intake()
#' food_intake()
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' XXXXXX
#' @export

food_intake <- function(avgmass,sdmass,avgenergy,sdenergy,efficiency,constanta,constantb,iterations){

if(missing(iterations)){iterations=100}

#Generate mass Gaussian distribution
mass_dist <- rnorm(iterations, avgmass, sdmass)

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
