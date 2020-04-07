# p3mapper: Pest Predation Pressure Mapper

Pest Predation Pressure Mapper (p3mapper) is an R package that incorporates a set of tools to estimate the predation-pressure of predators on crop-damaging arthropods by combining dietary information generated through DNA metabarcoding, energy budget of predators and density estimations refined with spatial distribution models. p3mapper is still a tool under development, and the article is under review.

This readme page describes the entire process to 1) Generate and 2) Analyse the data. The workflow described in the following enables replicating all the results shown in the original article. 

### Installation
p3mapper is available at github, and can be therefore installed in any R environment using devtools.
````R
#Install p3mapper
library(devtools)
#Install
install_github('anttonalberdi/p3mapper')
#Load
library(p3mapper)
````

### Dependencies
p3mapper has a number of dependencies that need to be installed and loaded in order to perform all operations.
````R
library(raster)
library(gstat)
library(rgdal)
library(viridis)
library(RColorBrewer)
library(hilldiv)
````

### Declare basic information
Before starting with any operation, it is necessary to declare some basic information
````R
#Working directory
#workdir = 'path/to/directory'
workdir = '/Users/anttonalberdi/Google\ Drive/Projects/2020_Pest_Bats/p3mapper/'
setwd(workdir)

#Predator species list
species.list <- c("Miniopterus_schreibersii","Myotis_capaccinii","Myotis_daubentonii","Myotis_emarginatus","Myotis_myotis","Rhinolophus_euryale","Rhinolophus_ferrumequinum")

#Default number of iterations
iterations = 100 #Number of iterations

#Default color pallette
color <- rev(magma(55))
````

## 1) Data generation

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

### Sum all overall indices and obtain average stats
````R
setwd("/Volumes/Haundi/p3mapper/index/")
sumvector <- c()
for(i in c(1:iterations)){
cat("Iteration",i,"\n")
rasterfiles <- list.files(pattern = paste("_",i,".asc",sep=""))
rasterstack <- stack(rasterfiles)
rastersum <-  calc(rasterstack, sum)
writeRaster(rastersum, paste("all_",i,".asc",sep=""), "ascii", overwrite=TRUE)
sumvector <- c(sumvector,cellStats(rastersum,stat=sum))
}
saveRDS(sumvector,"/Users/anttonalberdi/Google\ Drive/Projects/2020_Pest_Bats/p3mapper/index/all_index_avgs.Rdata")
````

### Obtain species index stats
````R
setwd("/Volumes/Haundi/p3mapper/index/")
species.list <- c("Miniopterus_schreibersii","Myotis_capaccinii","Myotis_daubentonii","Myotis_emarginatus","Myotis_myotis","Rhinolophus_euryale","Rhinolophus_ferrumequinum")

#Declare stats table
sum.table <- c()
#Iterate across predator species
for(species in species.list){
cat("     Species:",species,"\n")
sumvector <- c()
for(i in c(1:iterations)){
cat("     Iteration",i,"\n")
raster <- raster(paste(species,"_",i,".asc",sep=""))
sumvector <- c(sumvector,cellStats(raster,stat=sum))
}
sum.table <- rbind(sum.table,sumvector)
}
rownames(sum.table) <- species.list
write.table(sum.table,"/Users/anttonalberdi/Google\ Drive/Projects/2020_Pest_Bats/p3mapper/index/species_index_avgs.txt",row.names=TRUE,quote=FALSE,col.names=FALSE)
saveRDS(sum.table,"/Users/anttonalberdi/Google\ Drive/Projects/2020_Pest_Bats/p3mapper/index/species_index_avgs.Rdata")

apply(sum.table, 1, function(x) mean(x))
apply(sum.table, 1, function(x) sd(x))

#> apply(sum.table, 1, function(x) mean(x))
 #Miniopterus_schreibersii         Myotis_capaccinii        Myotis_daubentonii
#               14095814.1                  592123.3                31691767.8
#       Myotis_emarginatus             Myotis_myotis       Rhinolophus_euryale
#                1530369.8                 8279959.2                 3095295.2
#Rhinolophus_ferrumequinum
#                4011746.4
#> apply(sum.table, 1, function(x) sd(x))
 #Miniopterus_schreibersii         Myotis_capaccinii        Myotis_daubentonii
#                3005362.7                  213841.6                12796548.4
#       Myotis_emarginatus             Myotis_myotis       Rhinolophus_euryale
#                 315180.3                 2318694.5                  580221.0
#Rhinolophus_ferrumequinum
#                 858716.4

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

#Visualise and save table
stat.table
write.table(stat.table,"results/p3_species_statistics.csv",col.names=TRUE,row.names=TRUE,quote=FALSE)

#Obtain total values
stat.table <- read.table("results/p3_species_statistics.csv")
colSums(stat.table)

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
index_sd_all <- sqrt(index_var_all)

#Write rasters
writeRaster(index_avg_all,"index/all_avg.asc", "ascii", overwrite=TRUE)
writeRaster(index_sd_all,"index/all_sd.asc", "ascii", overwrite=TRUE)

#Prepare for plotting
index_avg_all[index_avg_all > 40] <- 40
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

### Species richness vs. predation pressure
````R
#Richness calculation
rm('richness')
for(species in species.list){
  distribution <- raster(paste("distribution/",species,".asc",sep=""))
  if(exists('richness')){richness <- richness + distribution}
  if(!exists('richness')){richness <- distribution}
}
richness <- richness * studyarea
writeRaster(richness,"richness.asc")
pressure <- raster(paste("index/all_avg.asc",sep=""))

stack <- stack(richness,pressure)
subsample <- sampleRegular(stack, 5000000, na.rm=TRUE)
colnames(subsample) <- c("Richness","Pressure")
subsample <- subsample[rowSums(subsample) > 0,]
subsample <- subsample[complete.cases(subsample),]
subsample <- subsample[subsample[,1] != 0,]
subsample <- subsample[subsample[,2] != 0,]
subsample <- subsample[sample(c(1:nrow(subsample)),10000),]
cor.test(subsample[,1], subsample[,2])

model = lm(Richness ~ Pressure, data = as.data.frame(subsample))
summary(model)
int =  model$coefficient["(Intercept)"]
slope = model$coefficient["Pressure"]

plot(subsample[,1] ~ subsample[,2],
     data = as.data.frame(subsample),
     pch=16,
     xlab = "Pressure",
     ylab = "Richness")

abline(int, slope, lty=1, lwd=2, col="blue")

````

### Overall ESE analysis normalised by richness
````R
pressure_all <- raster(paste("index/all_avg.asc",sep=""))

#Create index stack
rm(pressure_stack)
for(species in species.list){
  pressure_sp <- raster(paste("index/",species,"_avg.asc",sep=""))
  if(exists('pressure_stack')){pressure_stack <- stack(pressure_stack,pressure_sp)}
  if(!exists('pressure_stack')){pressure_stack <- pressure_sp}
}


#Create fine distribution
for(species in species.list){
  density <- raster(paste("density/",species,"_avg.asc",sep=""))
  density[density > 0] <- 1
  writeRaster(density,paste("distribution_fine/",species,".asc",sep=""),overwrite=TRUE)
}

#Calculate fine richness
rm(distribution_stack)
for(species in species.list){
  distribution_sp <- raster(paste("distribution_fine/",species,".asc",sep=""))
  if(exists('distribution_stack')){distribution_stack <- stack(distribution_stack,distribution_sp)}
  if(!exists('distribution_stack')){distribution_stack <- distribution_sp}
}
richness_fine <- calc(distribution_stack, fun=sum)
writeRaster(richness_fine,"richness_fine.asc",overwrite=TRUE)

#Obtain relative contribution
pressure_stack_rel <- pressure_stack / pressure_all

#Calculate ESE
fun <- function(x) {hill_div(x,qvalue=1)}

#Shannon diversity
D1 <- calc(pressure_stack_rel, fun)
richness <- raster("richness_fine.asc")
ESE <- (D1 - 1) / (richness - 1)
writeRaster(ESE,"ESE/all_ESE_norm.asc",overwrite=TRUE)

#Statsitics
cellStats(ESE,stat=mean)
cellStats(ESE,stat=sd)

#Plot
color <- rev(viridis(25))
pdf(paste("results/ESE_norm.pdf",sep=""),width=8,height=8)
plot(ESE, col = color,zlim=c(0,1))
dev.off()

````


### Overall ESE analysis
````R
pressure_all <- raster(paste("index/all_avg.asc",sep=""))

#Create index stack
rm(pressure_stack)
for(species in species.list){
  pressure_sp <- raster(paste("index/",species,"_avg.asc",sep=""))
  if(exists('pressure_stack')){pressure_stack <- stack(pressure_stack,pressure_sp)}
  if(!exists('pressure_stack')){pressure_stack <- pressure_sp}
}

#Obtain relative contribution
pressure_stack_rel <- pressure_stack / pressure_all

#Calculate ESE
fun <- function(x) {hill_div(x,qvalue=1)}
#Shannon diversity
D1 <- calc(pressure_stack_rel, fun)

#Normalised
D1max <- hill_div(c(0.142,0.142,0.142,0.142,0.142,0.142,0.142),qvalue=1)
D1min <- hill_div(c(1,0,0,0,0,0,0),qvalue=1)
ESE <- (D1 - D1min) / (D1max - D1min)
writeRaster(ESE,"ESE/all_ESE.asc",overwrite=TRUE)

#Plot
color <- rev(viridis(25))
pdf(paste("results/ESE.pdf",sep=""),width=8,height=8)
plot(ESE, col = color,zlim=c(0,1))
dev.off()

studyarea <- raster("studyarea.asc")
pdf(paste("results/studyarea.pdf",sep=""),width=8,height=8)
plot(studyarea, col = '#f4f4f4ff',zlim=c(1))
dev.off()

#ESE vs. richness

ESE <- raster("ESE/all_ESE.asc")
richness <- raster("richness.asc")
stack <- stack(ESE,richness)
subsample <- sampleRegular(stack, 1000000, na.rm=TRUE)
subsample <- subsample[rowSums(subsample) > 0,]
subsample <- subsample[complete.cases(subsample),]
subsample <- subsample[sample(c(1:nrow(subsample)),10000),]
cor.test(subsample[,1], subsample[,2])
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

### Ecosystem service potential - ESP (10km resolution)

````R
agriculture <- raster("agriculture/cropland_10km.asc")
#Change to total index (all species)
pressure <- raster(paste("index/all_avg.asc",sep=""))
pressure <- aggregate(pressure, fact=10, fun=mean)

#Set analysis area
area <- pressure
area[area > 0] <- 1

#Crop agriculture by area
agriculture_cropped <- agriculture * area

#Transfrom pressure into relative allvalues
pressure_rel <- pressure / cellStats(pressure,stat=max)

#Ecosystem service potential
ESP <- agriculture_cropped * pressure_rel

color <- rev(colorRampPalette(brewer.pal(n = 9, name = 'RdYlBu'))(60))
color <- color[c(20:60)]
color[1] <- "#f4f4f4ff"
pdf(paste("results/ESP.pdf",sep=""),height=8,width=8)
plot(ESP, col = color, zlim=c(0,0.8))
dev.off()

#80%
ESP80 <- ESP
ESP80[ESP80 > 0.80] <- 0
ESP80[ESP80 <= 0.70] <- 0
ESP80[ESP80 > 0] <- 1
agriculture_cropped2 <- agriculture_cropped * ESP80
agriculture_cropped2[agriculture_cropped2 == 0] <- NA
cellStats(agriculture_cropped2,stat=mean)
cellStats(agriculture_cropped2,stat=sd)
#Calculate area
pressure_perc[pressure_perc > 0] <- 1
cellStats(pressure_perc,stat=sum)

#70%
ESP70 <- ESP
ESP70[ESP70 > 0.70] <- 0
ESP70[ESP70 <= 0.60] <- 0
ESP70[ESP70 > 0] <- 1
agriculture_cropped2 <- agriculture_cropped * ESP70
agriculture_cropped2[agriculture_cropped2 == 0] <- NA
cellStats(agriculture_cropped2,stat=mean)
cellStats(agriculture_cropped2,stat=sd)
#Calculate area
agriculture_cropped2[agriculture_cropped2 > 0] <- 1
cellStats(agriculture_cropped2,stat=sum)

#60%
ESP60 <- ESP
ESP60[ESP60 > 0.60] <- 0
ESP60[ESP60 <= 0.50] <- 0
ESP60[ESP60 > 0] <- 1
agriculture_cropped2 <- agriculture_cropped * ESP60
agriculture_cropped2[agriculture_cropped2 == 0] <- NA
cellStats(agriculture_cropped2,stat=mean)
cellStats(agriculture_cropped2,stat=sd)
#Calculate area
agriculture_cropped2[agriculture_cropped2 > 0] <- 1
cellStats(agriculture_cropped2,stat=sum)


#50%
ESP50 <- ESP
ESP50[ESP50 > 0.50] <- 0
ESP50[ESP50 <= 0.40] <- 0
ESP50[ESP50 > 0] <- 1
agriculture_cropped2 <- agriculture_cropped * ESP50
agriculture_cropped2[agriculture_cropped2 == 0] <- NA
cellStats(agriculture_cropped2,stat=mean)
cellStats(agriculture_cropped2,stat=sd)
#Calculate area
agriculture_cropped2[agriculture_cropped2 > 0] <- 1
cellStats(agriculture_cropped2,stat=sum)





stack <- stack(agri_pressure,agriculture_cropped)
subsample <- sampleRegular(stack, 10000000, na.rm=TRUE)
colnames(subsample) <- c("Crops","Pressure")
subsample <- subsample[rowSums(subsample) > 0,]
subsample <- subsample[complete.cases(subsample),]
subsample <- subsample[sample(c(1:nrow(subsample)),10000),]
plot(subsample[,1], subsample[,2])



#Calculate correlation (global)
stack <- stack(agriculture_cropped,pressure)
subsample <- sampleRegular(stack, 10000000, na.rm=TRUE)
colnames(subsample) <- c("Crops","Pressure")
subsample <- subsample[rowSums(subsample) > 0,]
subsample <- subsample[complete.cases(subsample),]
subsample <- subsample[sample(c(1:nrow(subsample)),10000),]
cor.test(subsample[,1], subsample[,2])

#Sliding (not working with agri intensity sliding)
#subsample <- subsample[subsample[,1] > 0.3 & subsample[,1] <= 0.6,]
#cor.test(subsample[,1], subsample[,2])

model = lm(Crops ~ Pressure, data = as.data.frame(subsample))
summary(model)
int =  model$coefficient["(Intercept)"]
slope = model$coefficient["Pressure"]

plot(subsample[,1] ~ subsample[,2],
     data = as.data.frame(subsample),
     pch=16,
     xlab = "Pressure",
     ylab = "Intensity")

abline(int, slope, lty=1, lwd=2, col="blue")

````

### Average pressure on different crop intensities

````R
pressure <- raster(paste("index/Myotis_daubentonii_avg.asc",sep=""))
pressure <- aggregate(pressure, fact=10, fun=mean)

agriculture <- raster("agriculture/cropland_10km.asc")

#10%
agri10 <- agriculture
agri10[agri10 > 0.10] <- 0
agri10[agri10 > 0] <- 1
pressure_perc <- pressure * agri10
pressure_perc[pressure_perc == 0] <- NA
cellStats(pressure_perc,stat=mean)
cellStats(pressure_perc,stat=sd)
#Calculate area
pressure_perc[pressure_perc > 0] <- 1
cellStats(pressure_perc,stat=sum)

#20%
agri20 <- agriculture
agri20[agri20 > 0.20] <- 0
agri20[agri20 <= 0.10] <- 0
agri20[agri20 > 0] <- 1
pressure_perc <- pressure * agri20
pressure_perc[pressure_perc == 0] <- NA
cellStats(pressure_perc,stat=mean)
cellStats(pressure_perc,stat=sd)
#Calculate area
pressure_perc[pressure_perc > 0] <- 1
cellStats(pressure_perc,stat=sum)


#30%
agri30 <- agriculture
agri30[agri30 > 0.30] <- 0
agri30[agri30 <= 0.20] <- 0
agri30[agri30 > 0] <- 1
pressure_perc <- pressure * agri30
pressure_perc[pressure_perc == 0] <- NA
cellStats(pressure_perc,stat=mean)
cellStats(pressure_perc,stat=sd)
#Calculate area
pressure_perc[pressure_perc > 0] <- 1
cellStats(pressure_perc,stat=sum)

#40%
agri40 <- agriculture
agri40[agri40 > 0.40] <- 0
agri40[agri40 <= 0.30] <- 0
agri40[agri40 > 0] <- 1
pressure_perc <- pressure * agri40
pressure_perc[pressure_perc == 0] <- NA
cellStats(pressure_perc,stat=mean)
cellStats(pressure_perc,stat=sd)
#Calculate area
pressure_perc[pressure_perc > 0] <- 1
cellStats(pressure_perc,stat=sum)

#50%
agri50 <- agriculture
agri50[agri50 > 0.50] <- 0
agri50[agri50 <= 0.40] <- 0
agri50[agri50 > 0] <- 1
pressure_perc <- pressure * agri50
pressure_perc[pressure_perc == 0] <- NA
cellStats(pressure_perc,stat=mean)
cellStats(pressure_perc,stat=sd)
#Calculate area
pressure_perc[pressure_perc > 0] <- 1
cellStats(pressure_perc,stat=sum)

#60%
agri60 <- agriculture
agri60[agri60 > 0.60] <- 0
agri60[agri60 <= 0.50] <- 0
agri60[agri60 > 0] <- 1
pressure_perc <- pressure * agri60
pressure_perc[pressure_perc == 0] <- NA
cellStats(pressure_perc,stat=mean)
cellStats(pressure_perc,stat=sd)
#Calculate area
pressure_perc[pressure_perc > 0] <- 1
cellStats(pressure_perc,stat=sum)

#70%
agri70 <- agriculture
agri70[agri70 > 0.70] <- 0
agri70[agri70 <= 0.60] <- 0
agri70[agri70 > 0] <- 1
pressure_perc <- pressure * agri70
pressure_perc[pressure_perc == 0] <- NA
cellStats(pressure_perc,stat=mean)
cellStats(pressure_perc,stat=sd)
#Calculate area
pressure_perc[pressure_perc > 0] <- 1
cellStats(pressure_perc,stat=sum)

#80%
agri80 <- agriculture
agri80[agri80 > 0.80] <- 0
agri80[agri80 <= 0.70] <- 0
agri80[agri80 > 0] <- 1
pressure_perc <- pressure * agri80
pressure_perc[pressure_perc == 0] <- NA
cellStats(pressure_perc,stat=mean)
cellStats(pressure_perc,stat=sd)
#Calculate area
pressure_perc[pressure_perc > 0] <- 1
cellStats(pressure_perc,stat=sum)

#90%
agri90 <- agriculture
agri90[agri90 > 0.90] <- 0
agri90[agri90 <= 0.80] <- 0
agri90[agri90 > 0] <- 1
pressure_perc <- pressure * agri90
pressure_perc[pressure_perc == 0] <- NA
cellStats(pressure_perc,stat=mean)
cellStats(pressure_perc,stat=sd)
#Calculate area
pressure_perc[pressure_perc > 0] <- 1
cellStats(pressure_perc,stat=sum)

#100%
agri100 <- agriculture
agri100[agri100 <= 0.90] <- 0
agri100[agri100 > 0] <- 1
pressure_perc <- pressure * agri100
pressure_perc[pressure_perc == 0] <- NA
cellStats(pressure_perc,stat=mean)
cellStats(pressure_perc,stat=sd)
#Calculate area
pressure_perc[pressure_perc > 0] <- 1
cellStats(pressure_perc,stat=sum)

#Change to total index (all species)


````



### Local analyses

````R
color <- rev(magma(55))
color[1] <- "#f4f4f4ff"

pressure <- raster(paste("index/all_avg.asc",sep=""))

#balkans <- raster(xmn=19, xmx=25, ymn=41, ymx=43,resolution=0.008333334)
#wales <- raster(xmn=-5, xmx=-2, ymn=51.3, ymx=53.3, resolution=0.008333334)

#Basque
basque <- raster(xmn=-3.3, xmx=-0.6, ymn=41.8, ymx=43.5, resolution=0.008333334)
italy <- raster(xmn=10.5, xmx=14.5, ymn=42, ymx=44, resolution=0.008333334)
england <- raster(xmn=-3, xmx=-0.21, ymn=53, ymx=54.5, resolution=0.008333334)
provence <- raster(xmn=4.2, xmx=7.3, ymn=43, ymx=44.5, resolution=0.008333334)
britany <- raster(xmn=-4.9, xmx=-0.9, ymn=47.1, ymx=48.9, resolution=0.008333334)

region=britany
regionname="britany"

pressure_all <- raster(paste("index/all_avg.asc",sep=""))
pressure_all_region <- crop(pressure_all,region)
cellStats(pressure_all_region,stat='mean')
cellStats(pressure_all_region,stat='sd')

pdf(paste("results/pressure_",regionname,".pdf",sep=""),height=8,width=8)
plot(pressure_all_region, col = color,zlim=c(0,40))
dev.off()

richness <- raster("richness_fine.asc")
richness_crop <- crop(richness,region)
cellStats(richness_crop,stat=mean)
cellStats(richness_crop,stat=sd)

rm(pressure_sp_region_stack)
pressure.vector <- c()
for(species in species.list){
  pressure_sp <- raster(paste("index/",species,"_avg.asc",sep=""))
  pressure_sp_region <- crop(pressure_sp,region)
  pressure.vector <- c(pressure.vector,cellStats(pressure_sp_region,stat=mean))

  if(exists('pressure_sp_region_stack')){pressure_sp_region_stack <- stack(pressure_sp_region_stack,pressure_sp_region)}
  if(!exists('pressure_sp_region_stack')){pressure_sp_region_stack <- pressure_sp_region}

  #Plot
  pressure_sp_region[pressure_sp_region > 40] <- 40
  pdf(paste("results/pressure_",regionname,"_",species,".pdf",sep=""),width=10,height=6)
  plot(pressure_sp_region, col = color,zlim=c(0,40))
  dev.off()
}
pressure.vector.rel <- round(pressure.vector/sum(pressure.vector),3)
names(pressure.vector.rel) <- species.list
pressure.vector.rel

#Calculate ESE
fun <- function(x) {hill_div(x,qvalue=1)}
#Shannon diversity
D1_region <- calc(pressure_sp_region_stack, fun)
ESE_region <- (D1_region - 1) / (richness_crop - 1)
cellStats(ESE_region,stat=mean)
cellStats(ESE_region,stat=sd)
writeRaster(ESE_region,paste("ESE/",regionname,".asc",sep=""),overwrite=TRUE)

````


### Plots for Figure 1
````R
#Myotis capaccinii distribution
distribution <- raster("distribution/Myotis_capaccinii.asc")
studyarea <- raster("studyarea.asc")
distribution2 <- distribution * studyarea
pdf("results/distribution_Mca.pdf",width=8,height=8)
plot(distribution2, col = c("#f4f4f4","#a57c39"))
dev.off()

#Myotis capaccinii interpolated density
iterations=100
density=read.csv(paste("density/Myotis_capaccinii.csv",sep=""))
rownames(density) <- paste("S",c(1:nrow(density)),sep="")
density_matrix <- c()
for(r in c(1:nrow(density))){
  row <- density[r,]
  density_iter <- rnorm(iterations, row$Average, row$SD)
  density_matrix <- rbind(density_matrix,density_iter)
}
rownames(density_matrix) <- rownames(density)
colnames(density_matrix) <- paste("I",c(1:ncol(density_matrix)),sep="")
coordinates(density) <- ~x + y
reference_coord <- as(raster("enm/Myotis_capaccinii.asc"), "SpatialPixelsDataFrame")
i=1
density_interpol <- raster(idw(formula = density_matrix[,i] ~ 1, locations = density, newdata = reference_coord, idp = 1, debug.level = 0))
density_interpol <- density_interpol * distribution

color <- colorRampPalette(c("#e5ca9a","#cc9f45", "#a57c39","#7a5723","#492f0d"))(50)
color[1] <- "#f4f4f4ff"
pdf("results/density_raw_Mca.pdf",width=8,height=8)
plot(density_interpol, col = color)
dev.off()

#Myotis capaccinii refined density
density <- raster("density/Myotis_capaccinii_avg.asc")
studyarea <- raster("studyarea.asc")
density2 <- density * studyarea
color <- colorRampPalette(c("#ede491","#7c7322"))(50)
color[1] <- "#f4f4f4ff"
pdf("results/density_Mca.pdf",width=8,height=8)
plot(density2, col = color)
dev.off()


foreste <- SpatialPoints(matrix(c(11.6629,43.9719),ncol=2), proj4string=CRS(as.character(NA)), bbox = NULL)

extract(density, foreste, method='simple')
````
