# p3mapper: Pest Predation Pressure Mapper

p3mapper is still a tool under development.

## Workflow

### Declare basic information
````R
#Install p3mapper
library(devtools)
install_github('anttonalberdi/p3mapper')

#Set working directory
#workdir = 'path/to/directory'
workdir = '/Users/anttonalberdi/Google\ Drive/Projects/2020_Pest_Bats/p3mapper/'

#Specify global parameters
iterations = 100 #Number of iterations

#Load libraries
library(p3mapper)
library(raster)
library(gstat)
library(rgdal)
library(viridis)
color <- rev(magma(55))
````

### Estimate predator densities

````R
setwd(workdir)

species="Rhinolophus_ferrumequinum"

#Specify function-specific parameters
density=read.csv(paste("density/",species,".csv",sep=""))
enm=raster(paste("enm/",species,".asc",sep=""))
distribution=raster(paste("distribution/",species,".asc",sep=""))

#Run function and save object
predator_density(base,density,enm,distribution,iterations)
save(sp_density, file=paste("density/",species,".Rdata",sep=""))
#sp_density <- readRDS(paste("density/",species,".Rdata",sep=""))

#Get average
sp_density_avg <- calc(sp_density, fun = mean, na.rm=TRUE)
#pdf(paste("density/",species,"_avg.pdf",sep=""),height=8,width=8)
#plot(sp_density_avg, col = color,zlim=c(0,12))
#dev.off()
writeRaster(sp_density_avg, paste("density/",species,"_avg.asc",sep=""), "ascii", overwrite=TRUE)

#Get SD
sp_density_sd <- calc(sp_density, fun=sd, na.rm=TRUE)
writeRaster(sp_density_sd, paste("density/",species,"_sd.asc",sep=""), "ascii", overwrite=TRUE)

#Get SE
sp_density_se <- sp_density_sd/sqrt(nlayers(sp_density_sd))
writeRaster(sp_density_se, paste("density/",species,"_se.asc",sep=""), "ascii", overwrite=TRUE)

````


### Estimate overall prey consumption

````R
#Specify function-specific parameters
species="Myotis emarginatus"
species2 <- gsub(" ","_",species)
bodymass <- read.csv("intake/species_bodymass.csv")

avgmass=bodymass[bodymass[,1]== species,2]
sdmass=bodymass[bodymass[,1]== species,3]
efficiency=0.88
avgenergy=5.20
sdenergy=0.73
constanta=6.49
constantb=0.681
iternumber=100

#Run function
food_consumption <- food_intake(avgmass,sdmass,avgenergy,sdenergy,efficiency,constanta,constantb,iterations)
saveRDS(food_consumption, paste("intake/",species2,".Rdata",sep=""))

````

### Estimate pest proportion

````R
#Declare species list
species.list <- c("Miniopterus_schreibersii","Myotis_capaccinii","Myotis_daubentonii","Myotis_emarginatus","Myotis_myotis","Rhinolophus_euryale","Rhinolophus_ferrumequinum")

for(species in species.list){

  cat(species,"\n")
  species2 <- gsub("_"," ",species)
  otutable <- read.table(paste("diet/",species,".tsv",sep=""),header=TRUE,sep="\t")
  sampleinfo <- read.table("diet/sample.info.csv",sep=",",header=TRUE)
  sampleinfo <- sampleinfo[sampleinfo$Species == species2,c(1,3)]
  siteinfo <- read.table("diet/site_coordinate.csv",sep=",",header=TRUE)
  OTU_taxonomy_pest <- read.table("diet/OTU_taxonomy_pest.csv",sep=",",header=TRUE)
  OTU_taxonomy_pest <- OTU_taxonomy_pest[,c(1,6)]
  OTU_taxonomy_pest[OTU_taxonomy_pest$Pest > 0,2] <- 1

  counttable=otutable
  sampleinfo=sampleinfo
  siteinfo=siteinfo
  pestinfo=OTU_taxonomy_pest
  distribution=raster(paste("distribution/",species,".asc",sep=""))
  iterations=100

  #Run function
  pest_proportion <- pest_proportion(counttable,sampleinfo,siteinfo,pestinfo,distribution,iterations)

  saveRDS(pest_proportion, paste("diet/",species,".Rdata",sep=""))
  #sp_density <- readRDS(paste("density/",species,".Rdata",sep=""))

  #Get average
  pest_proportion_avg <- calc(pest_proportion, fun = mean, na.rm=TRUE)
  #pdf(paste("density/",species,"_avg.pdf",sep=""),height=8,width=8)
  #plot(sp_density_avg, col = color,zlim=c(0,12))
  #dev.off()
  writeRaster(pest_proportion_avg, paste("diet/",species,"_avg.asc",sep=""), "ascii", overwrite=TRUE)

  #Get SD
  pest_proportion_sd <- calc(pest_proportion, fun=sd, na.rm=TRUE)
  writeRaster(pest_proportion_sd, paste("diet/",species,"_sd.asc",sep=""), "ascii", overwrite=TRUE)

  #Get SE
  pest_proportion_se <- pest_proportion_sd/sqrt(nlayers(pest_proportion))
  writeRaster(pest_proportion_se, paste("diet/",species,"_se.asc",sep=""), "ascii", overwrite=TRUE)
}

````
### Merge all data
````R

species="Rhinolophus_ferrumequinum"

#Memory issues: https://stackoverflow.com/questions/51248293/error-vector-memory-exhausted-limit-reached-r-3-5-0-macos
density <- readRDS(paste("density/",species,".Rdata",sep=""))
food_intake <- readRDS(paste("intake/",species,".Rdata",sep=""))
pest_proportion <- readRDS(paste("diet/",species,".Rdata",sep=""))

for (i in c(1:iterations)){
predator_pressure[[i]] <- density[[i]] * food_intake[i] * pest_proportion[[i]]
}

i=1
predator_pressure <- density[[i]] * food_intake[i] * pest_proportion[[i]]
````



### Plot overall statistics
````R

density_avg <- raster(paste("density/",species,"_avg.asc",sep=""))
plot(density(density_avg[density_avg > 0], adjust = 3))

````

### Food intake distribution per species

````R
species.list <- c("Miniopterus_schreibersii","Myotis_capaccinii","Myotis_daubentonii","Myotis_emarginatus","Myotis_myotis","Rhinolophus_euryale","Rhinolophus_ferrumequinum")

for(species in species.list){
  food_intake <- readRDS(paste("intake/",species,".Rdata",sep=""))
  pdf(paste("results/intake_",species,".pdf",sep=""),width=10,height=6)
  plot(density(food_intake,adjust = 3),xlim=c(0,20))
  dev.off()
}
````


### Overlay with agricultural intensity
````R

#http://www.earthstat.org/cropland-pasture-area-2000/
agriintensity <- raster("agriculture/Cropland2000_5m.tif")
agriintensity_europe <- crop(agriintensity,distribution)
agriintensity_europe <- disaggregate(agriintensity_europe, fact=10)
agriintensity_europe <- crop(agriintensity_europe,distribution)
extent(agriintensity_europe) <- extent(distribution)

agriintensity_europe_MSc <- agriintensity_europe * distribution
predator_pressure_MSc <- (predator_pressure - cellStats(predator_pressure,stat=min))/(cellStats(predator_pressure,stat=max)-cellStats(predator_pressure,stat=min))

agriintensity_europe_MSC_service <- agriintensity_europe_MSc - predator_pressure_MSc

density_avg <- raster()
hist(density_avg[density_avg > 0])
````
