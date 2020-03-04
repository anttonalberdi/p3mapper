#' @title Pest proportion estimation
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords mass efficiency energy budget
#' @description Computes pest consumption estimations.
#' @param counttable Average body mass of the predator (grams).
#' @param sampleinfo Standard deviation of the body mass of the predator (grams).
#' @param siteinfo Maximum body mass (top 5% percentile) of the predator (grams).
#' @param pestinfo Digestion efficiency of the predator (percentage).
#' @param distribution Average energy content of the prey (kJ/g).
#' @param iterations Number of iterations (default 100).
#' @usage pest_proportion(counttable,sampleinfo,siteinfo,pestinfo,distribution,iterations)
#' @return A value (if one iteration) or vector (if multiple iterations) of prey consumption estimations.
#' @import raster gstat rgdal
#' @examples
#' pest_proportion()
#' pest_proportion()
#' @references
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.\cr\cr
#' XXXXXX
#' @export

pest_proportion <- function(counttable,sampleinfo,siteinfo,pestinfo,distribution,iterations){

#Intersect sample with counttable information
sampleinfo <- sampleinfo[sampleinfo[,1] %in% colnames(counttable),]
sampleinfo[,1] <- as.character(sampleinfo[,1])
sampleinfo[,2] <- as.character(sampleinfo[,2])
site.list <- unique(sampleinfo[,2])
siteinfo <- siteinfo[siteinfo[,1] %in% site.list,]
sampleinfo.new <- sampleinfo

#Create reference coordinates
reference_coord <- as(distribution, "SpatialPixelsDataFrame")
siteinfo.sd <- siteinfo
coordinates(siteinfo.sd) <- ~x + y

#Iterate measurements
site.table <- c()
for(i in c(1:iterations)){

###
# Pest consumption bootstrapping
###

cat("Iteration",i,"\n")

#Resample (randomly replace 1 sample )
for(site in site.list){
  in.samples <- sampleinfo[sampleinfo[,2] == site,1]
  out.samples <- sampleinfo[sampleinfo[,2] != site,1]
  in.samples[sample(c(1:length(in.samples)),1)] <- out.samples[sample(c(1:length(out.samples)),1)]
  sampleinfo.new[sampleinfo.new[,2] == site,1] <- in.samples
}

#Create site-specific data
site.list <- unique(sampleinfo[,2])

  site.vector <- c()
  for (site in site.list){
    #Subset
    samples <- sampleinfo.new[sampleinfo.new$Site == site,1]
    counttable.sub <- counttable[,samples]

    #Compute abundance
    abundance <- rowMeans(counttable.sub)
    abundance <- abundance[abundance > 0]

    #Compute relativa pest abundance
    pestinfo.subset <- pestinfo[pestinfo$OTU %in% names(abundance),]
    pestinfo.subset.merged <- merge(t(t(abundance)),pestinfo.subset,by.x="row.names",by.y="OTU")
    pestinfo.subset.merged.filtered <- pestinfo.subset.merged[pestinfo.subset.merged$Pest > 0,]
    pest_rel_abun <- sum(pestinfo.subset.merged.filtered[,2])

    #
    site.vector <- c(site.vector,pest_rel_abun)
  }


site.table <- cbind(site.list,site.vector)
colnames(site.table)[1] <- "Site"

###
# Pest consumption interpolating
###

pest_location <- merge(site.table,siteinfo,by="Site")
pest_location[,2] <- as.numeric(as.character(pest_location[,2]))
pest_location[,3] <- as.numeric(as.character(pest_location[,3]))
pest_location[,4] <- as.numeric(as.character(pest_location[,4]))

pest_interpol <- raster(idw(formula = pest_location[,2] ~ 1, locations = siteinfo.sd, newdata = reference_coord, idp = 1, debug.level = 0))

#Crop distribution
pest_interpol <- pest_interpol * distribution

names(pest_interpol) <- paste("Iter",i,sep="")

if(i == 1){
pest_interpol_stack <- pest_interpol
}

if(i > 1){
pest_interpol_stack <- addLayer(pest_interpol_stack, pest_interpol)

}

}

return(pest_interpol_stack)

}
