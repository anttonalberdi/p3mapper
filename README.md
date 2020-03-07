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
setwd(workdir)

#Specify global parameters
iterations = 100 #Number of iterations

#Load libraries
library(p3mapper)
library(raster)
library(gstat)
library(rgdal)
library(viridis)
color <- rev(magma(55))

#Declare species list
species.list <- c("Miniopterus_schreibersii","Myotis_capaccinii","Myotis_daubentonii","Myotis_emarginatus","Myotis_myotis","Rhinolophus_euryale","Rhinolophus_ferrumequinum")

````

### Estimate predator densities

````R

species="Miniopterus_schreibersii"

#Specify function-specific parameters
base=paste(workdir,"density/",species,sep="")
density=read.csv(paste("density/",species,".csv",sep=""))
enm=raster(paste("enm/",species,".asc",sep=""))
distribution=raster(paste("distribution/",species,".asc",sep=""))

#Run function
predator_density(base,density,enm,distribution,iterations)
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

species="Miniopterus_schreibersii"

species2 <- gsub("_"," ",species)
otutable <- read.table(paste("diet/",species,".tsv",sep=""),header=TRUE,sep="\t")
sampleinfo <- read.table("diet/sample.info.csv",sep=",",header=TRUE)
sampleinfo <- sampleinfo[sampleinfo$Species == species2,c(1,3)]
siteinfo <- read.table("diet/site_coordinate.csv",sep=",",header=TRUE)
OTU_taxonomy_pest <- read.table("diet/OTU_taxonomy_pest.csv",sep=",",header=TRUE)
OTU_taxonomy_pest <- OTU_taxonomy_pest[,c(1,6)]
OTU_taxonomy_pest[OTU_taxonomy_pest$Pest > 0,2] <- 1

base=paste(workdir,"diet/",species,sep="")
counttable=otutable
sampleinfo=sampleinfo
siteinfo=siteinfo
pestinfo=OTU_taxonomy_pest  
distribution=raster(paste("distribution/",species,".asc",sep=""))
iterations=100

#Run function
pest_proportion(base,counttable,sampleinfo,siteinfo,pestinfo,distribution,iterations)
````


### Merge all data
````R

species="Miniopterus_schreibersii"

base=paste(workdir,"index/",species,sep="")

for (i in c(1:iterations)){
density <- raster(paste("density/",species,"_",i,".asc",sep=""))
food_intake <- readRDS(paste("intake/",species,".Rdata",sep=""))
pest_proportion <- raster(paste("diet/",species,"_",i,".asc",sep=""))

#Calculate index
predator_pressure <- density * food_intake[i] * pest_proportion

#Sum to previous data to calculate average
if(i == 1){
  predator_pressure_avg <- predator_pressure
  }
if(i > 1){
  predator_pressure_avg <- predator_pressure_avg + predator_pressure
  }
}

#Divide by iteration number to obtain average
predator_pressure_avg <- predator_pressure_avg / iterations
writeRaster(predator_pressure_avg, paste(base,"_avg.asc",sep=""), "ascii", overwrite=TRUE)

#Generate SD
cat("Generating Standard Deviation",i,"\n")
for(i in c(1:iterations)){
  cat("   Iteration",i,"\n")
raster_iter <- raster(paste(base,"_",i,".asc",sep=""))
  if(i == 1){
    raster_sum <- (raster_iter - predator_pressure_avg)^2
  }
  if(i > 1){
    raster_sum <- raster_sum + (raster_iter - predator_pressure_avg)^2
  }
}

predator_pressure_sd <- sqrt(raster_sum/iterations)
writeRaster(predator_pressure_sd, paste(base,"_sd.asc",sep=""), "ascii", overwrite=TRUE)

#Generate SE
predator_pressure_se <- predator_pressure_sd/sqrt(iterations)
writeRaster(predator_pressure_se, paste(base,"_se.asc",sep=""), "ascii", overwrite=TRUE)




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
