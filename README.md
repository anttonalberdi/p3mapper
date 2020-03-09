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
library(RColorBrewer)
color <- rev(magma(55))

#Declare species list
species.list <- c("Miniopterus_schreibersii","Myotis_capaccinii","Myotis_daubentonii","Myotis_emarginatus","Myotis_myotis","Rhinolophus_euryale","Rhinolophus_ferrumequinum")

````

### Estimate predator densities

````R

species="Rhinolophus_euryale"

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

species="Myotis_daubentonii"

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


### Compute models
````R

species="Myotis_daubentonii"

base=paste(workdir,"index/",species,sep="")

cat("Generating iterative maps \n")
for (i in c(1:iterations)){
  cat("   Iteration",i,"\n")
  density <- raster(paste("density/",species,"_",i,".asc",sep=""))
  food_intake <- readRDS(paste("intake/",species,".Rdata",sep=""))
  pest_proportion <- raster(paste("diet/",species,"_",i,".asc",sep=""))

  #Calculate index
  predator_pressure <- density * food_intake[i] * pest_proportion

  #Name raster
  names(predator_pressure) <- paste("Iter",i,sep="")

  #Save to file
  writeRaster(predator_pressure, paste(base,"_",i,".asc",sep=""), "ascii", overwrite=TRUE)

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
cat("Generating Standard Deviation \n")
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

### Pest proportion statistics of predator species (based on empirical data)

````R
species="Rhinolophus_ferrumequinum"

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

#Intersect sample with counttable information
sampleinfo <- sampleinfo[sampleinfo[,1] %in% colnames(counttable),]
sampleinfo[,1] <- as.character(sampleinfo[,1])
sampleinfo[,2] <- as.character(sampleinfo[,2])
site.list <- unique(sampleinfo[,2])
siteinfo <- siteinfo[siteinfo[,1] %in% site.list,]

site.list <- unique(sampleinfo[,2])

  site.vector <- c()
  for (site in site.list){
    #Subset
    samples <- sampleinfo[sampleinfo$Site == site,1]
    counttable.sub <- counttable[,samples]

    #Compute abundance
    abundance <- rowMeans(counttable.sub)
    abundance <- abundance[abundance > 0]

    #Compute relative pest abundance
    pestinfo.subset <- pestinfo[pestinfo$OTU %in% names(abundance),]
    pestinfo.subset.merged <- merge(t(t(abundance)),pestinfo.subset,by.x="row.names",by.y="OTU")
    pestinfo.subset.merged.filtered <- pestinfo.subset.merged[pestinfo.subset.merged$Pest > 0,]
    pest_rel_abun <- sum(pestinfo.subset.merged.filtered[,2])

    #
    site.vector <- c(site.vector,pest_rel_abun)
  }


site.table <- cbind(site.list,site.vector)
colnames(site.table)[1] <- "Site"
site.table <- as.data.frame(site.table)
site.table[,2] <- as.numeric(as.character(site.table[,2]))

#Proportion of pests
mean(site.table[,2])
sd(site.table[,2])

#Add intake data
food_intake <- readRDS(paste("intake/",species,".Rdata",sep=""))

site.table.intake <- c()
for(r in c(1:nrow(site.table))){
allvalues <- site.table[r,2] * food_intake
row <- c(mean(allvalues),var(allvalues))
site.table.intake <- rbind(site.table.intake,row)
}

#Consumption of pests per animal
mean(site.table.intake[,1])
sqrt(mean(site.table[,2]))

````

### Pest proportion statistics of predator species (based on the models)

````R
#Declare predator species list
species.list <- c("Miniopterus_schreibersii","Myotis_capaccinii","Myotis_daubentonii","Myotis_emarginatus","Myotis_myotis","Rhinolophus_euryale","Rhinolophus_ferrumequinum")

#Declare stats table
stat.table <- c()

#Iterate across predator species
for(species in species.list){
  avg <- raster(paste("diet/",species,"_avg.asc",sep=""))
  sd <- raster(paste("diet/",species,"_sd.asc",sep=""))
  distribution <- raster(paste("distribution/",species,".asc",sep=""))

  #Calculate distribution area
  dis_area <- cellStats(distribution, stat=sum)

  #Calculate variance
  variance <- sd^2

  #Average pest proportion
  cell_avg <- cellStats(avg, stat=sum)/dis_area*100
  cell_sd <- sqrt(cellStats(variance, stat=sum)/dis_area)*100

  row <- c(cell_avg, cell_sd)
  stat.table <- rbind(stat.table,row)
}
rownames(stat.table) <- species.list
colnames(stat.table) <- c("km2_avg","km2_sd","total_avg","total_sd","total_area")

````

### Predation pressure statistics of predator species
The following script calculates the overall pest predation pressure statistics (average and standard deviation) per species, both at cell (km2) and total scales.

````R
#Declare predator species list
species.list <- c("Miniopterus_schreibersii","Myotis_capaccinii","Myotis_daubentonii","Myotis_emarginatus","Myotis_myotis","Rhinolophus_euryale","Rhinolophus_ferrumequinum")

#Load study area
studyarea <- raster("studyarea.asc")

#Declare stats table
stat.table <- c()
#Iterate across predator species
for(species in species.list){
  #Load species-specific rasters
  avg <- raster(paste("index/",species,"_avg.asc",sep=""))
  sd <- raster(paste("index/",species,"_sd.asc",sep=""))
  distribution <- raster(paste("distribution/",species,".asc",sep=""))

  #Crop distribution by study area
  distribution <- distribution * studyarea

  #Calculate distribution area
  dis_area <- cellStats(distribution, stat=sum)

  #Calculate variance
  variance <- sd^2

  #Average pressure
  cell_avg <- cellStats(avg, stat=sum)/dis_area
  cell_sd <- sqrt(cellStats(variance, stat=sum)/dis_area)

  #Total Pressure
  total_avg <- cellStats(avg, stat=sum)
  total_sd <- sqrt(cellStats(variance, stat=sum))

  row <- c(cell_avg, cell_sd, total_avg, total_sd,dis_area)
  stat.table <- rbind(stat.table,row)
}
rownames(stat.table) <- species.list
colnames(stat.table) <- c("km2_avg","km2_sd","total_avg","total_sd","total_area")

````



### Merge species averages and SDs
````R
species.list <- c("Miniopterus_schreibersii","Myotis_capaccinii","Myotis_daubentonii","Myotis_emarginatus","Myotis_myotis","Rhinolophus_euryale","Rhinolophus_ferrumequinum")

rm(index_avg_all)
rm(index_sd_all)
for(species in species.list){
index_avg <- raster(paste("index/",species,"_avg.asc",sep=""))
index_sd <- raster(paste("index/",species,"_sd.asc",sep=""))
index_var <- index_sd^2
if(!exists("index_avg_all")){index_avg_all <- index_avg}
if(exists("index_avg_all")){index_avg_all <- index_avg_all + index_avg}
if(!exists("index_var_all")){index_var_all <- index_var}
if(exists("index_var_all")){index_var_all <- index_var_all + index_var}
}

index_avg_all <- index_avg_all
index_avg_all[index_avg_all > 40] <- 40
index_sd_all <- sqrt(index_var_all)
index_sd_all[index_sd_all > 40] <- 40

#Plot maps
pdf(paste("results/all_index_avg.pdf",sep=""),height=8,width=8)
plot(index_avg_all, col = color,zlim=c(0,40))
dev.off()

pdf(paste("results/all_index_sd.pdf",sep=""),height=8,width=8)
plot(index_sd_all, col = color,zlim=c(0,40))
dev.off()
````


### Plot average maps
````R
for(species in species.list){
index_avg <- raster(paste("index/",species,"_avg.asc",sep=""))
index_avg[index_avg > 20] <- 20
index_sd <- raster(paste("index/",species,"_sd.asc",sep=""))
index_sd[index_sd > 20] <- 20

#Define color pallette
color <- rev(magma(55))
color[1] <- "#f4f4f4ff"

pdf(paste("results/",species,"_index_avg.pdf",sep=""),height=8,width=8)
plot(index_avg, col = color,zlim=c(0,20))
dev.off()

pdf(paste("results/",species,"_index_sd.pdf",sep=""),height=8,width=8)
plot(index_sd, col = color,zlim=c(0,20))
dev.off()
}
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

### Modify agricultural intensity
````R
#
#http://www.earthstat.org/cropland-pasture-area-2000/
agriintensity <- raster("agriculture/Cropland2000_5m.tif")
agriintensity_europe <- crop(agriintensity,distribution)
#agriintensity_europe <- disaggregate(agriintensity_europe, fact=10)
studyarea <- raster("studyarea.asc")
studyarea <- aggregate(studyarea, fact=10)
agriintensity_europe <- crop(agriintensity_europe,studyarea)
extent(agriintensity_europe) <- extent(studyarea)
agriintensity_europe <- agriintensity_europe * studyarea
writeRaster(agriintensity_europe,"agriculture/cropland_10km.asc",overwrite=TRUE)
````

### Plot cropland intensity
````R
agriculture <- raster("agriculture/cropland_10km.asc")

#Plot agricultural intensity
color <- colorRampPalette(brewer.pal(n = 9, name = 'Greens'))(55)
color[1] <- "#f4f4f4ff"
pdf(paste("results/agriculture_cropland_10km.pdf",sep=""),height=8,width=8)
plot(agriculture, col = color,zlim=c(0,1))
dev.off()
````

### Correlate predation pressure with cropland intensity (10km resolution)

````R
agriculture <- raster("agriculture/cropland_10km.asc")
#Change to total index (all species)
pressure <- raster(paste("index/",species,"_avg.asc",sep=""))
pressure <- aggregate(pressure, fact=10, fun=mean)

#Set analysis area
area <- pressure
area[area > 0] <- 1

#Crop agriculture by area
agriculture_cropped <- agriculture * area

#Calculate correlation (cell-wise)
corlocal <- corLocal(agriculture_cropped, pressure, method="pearson", na.rm=TRUE)

#Calculate correlation (global)
stack <- stack(agriculture_cropped,pressure)
subsample <- sampleRegular(stack, 10000000, na.rm=TRUE)
subsample <- subsample[rowSums(subsample) > 0,]
subsample <- subsample[complete.cases(subsample),]
cor.test(subsample[,1], subsample[,2])

#Sliding (not working with agri intensity sliding)
subsample <- subsample[subsample[,1] > 0.3 & subsample[,1] <= 0.6,]
cor.test(subsample[,1], subsample[,2])

model = lm(layer ~ Miniopterus_schreibersii_avg, data = as.data.frame(subsample))
summary(model)
int =  model$coefficient["(Intercept)"]
slope =model$coefficient["Miniopterus_schreibersii_avg"]

plot(subsample[,1] ~ subsample[,2],
     data = as.data.frame(subsample),
     pch=16,
     xlab = "Pressure",
     ylab = "Intensity")

abline(int, slope, lty=1, lwd=2, col="blue")

````
