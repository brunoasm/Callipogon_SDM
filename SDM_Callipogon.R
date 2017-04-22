# This script takes as input a table with occurrences of Callipogon relictus
# It creates distribution maps and uses R package biomod2 to do species distribution modelling

###############################################################################
# 1 - load packages
library(biomod2)
library(rworldmap)
library(rgdal)
library(rgeos)
 
###############################################################################
# 2 - read occurence data
locs <- read.table("relictus_distribution_final.txt", header = T, sep = "\t")
#uncomment below to remove suspect record
#locs <- locs[-55,] #removing suspect record
#rownames(locs)<-as.character(seq(1,dim(locs)[1]))

###############################################################################
# 3 - make an occurence map with utm projection
# this way, coordinates are in meters
# map will be cropped to a bounding box 200 km around coordinate limits
lonlat <- SpatialPoints(as.matrix(locs[c("longitude", "latitude")]), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "))
utm_locs <- spTransform(lonlat,CRS(" +proj=utm +zone=52N +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0"))

worldmap <- getMap(resolution = 'high')
mapborders <- extent(lonlat) + c(-50,50,-50,50)
ll_map <- crop(worldmap,mapborders) #has to do a preliminary crop in other to project to utm
utm_map <- spTransform(ll_map,CRS(" +proj=utm +zone=52N +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0"))

utm_borders <- extent(utm_locs) + c(-200000,200000,-200000,200000)
#utm_map <- crop(utm_map,utm_borders)

dir.create('graphs')
pdf(file = 'graphs/occurences.pdf', width = 7*(extent(utm_map)@xmax-extent(utm_map)@xmin)/(extent(utm_map)@ymax-extent(utm_map)@ymin),
                             height = 7)


mapCountryData(utm_map,
               xlim = c(utm_borders@xmin, utm_borders@xmax),
               ylim = c(utm_borders@ymin,utm_borders@ymax),
               colourPalette=c('cornsilk','cornsilk'),
               addLegend = F,
               borderCol=gray(level = 0.3),
               lwd =1,
               mapTitle = "")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col="lightblue")
mapCountryData(utm_map,
               xlim = c(utm_borders@xmin, utm_borders@xmax),
               ylim = c(utm_borders@ymin,utm_borders@ymax),
               colourPalette=c('cornsilk','cornsilk'),
               addLegend = F,
               borderCol=gray(level = 0.3),
               lwd =1,
               mapTitle = "",
               add=T)




box()

points(utm_locs, col = 'red')
# add some nice state labels ...
text(x=-1000000, y=4300000, "China", col="cornsilk4", cex=1)
text(x=900000, y=5700000, "Russia", col="cornsilk4", cex=1)
text(x=850000, y=4500000, "N. Korea", col="cornsilk4", cex=0.8, adj = 1)
text(x=800000, y=4200000, "S. Korea", col="cornsilk4", cex=0.8, adj = 1)

dev.off()

###############################################################################
# 4 - subsample occurences to reduce spatial autocorrelation
# We will remove closest samples until no two samples are closest than mindist (in meters)
mindist = 20000
utm_subs = utm_locs

distmatrix <- as.matrix(dist(utm_locs@coords,method = 'euclidean', diag = FALSE, upper = TRUE))
diag(distmatrix) <- NA

while(min(distmatrix,na.rm = T) < mindist){
  closest <- which(distmatrix == min(distmatrix, na.rm = T), arr.ind=T)[1:2]
  remove_first = sort(distmatrix[closest[1],])[2] < sort(distmatrix[closest[2],])[2]
  distmatrix <- distmatrix[-closest[remove_first+1],-closest[remove_first+1]]
}

utm_subs@coords <- utm_subs@coords[as.integer(rownames(distmatrix)),]

pdf(file = 'graphs/occurences_subsampled.pdf', width = 7*(extent(utm_map)@xmax-extent(utm_map)@xmin)/(extent(utm_map)@ymax-extent(utm_map)@ymin),
    height = 7)
mapCountryData(utm_map,
               xlim = c(utm_borders@xmin, utm_borders@xmax),
               ylim = c(utm_borders@ymin,utm_borders@ymax),
               colourPalette=c('cornsilk','cornsilk'),
               addLegend = F,
               borderCol=gray(level = 0.3),
               lwd =1,
               mapTitle = "")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col="lightblue")
mapCountryData(utm_map,
               xlim = c(utm_borders@xmin, utm_borders@xmax),
               ylim = c(utm_borders@ymin,utm_borders@ymax),
               colourPalette=c('cornsilk','cornsilk'),
               addLegend = F,
               borderCol=gray(level = 0.3),
               lwd =1,
               mapTitle = "",
               add=T)
box()
points(utm_locs, col = 'goldenrod2')
points(utm_subs, col = 'red')
# add some nice state labels ...
text(x=-1000000, y=4300000, "China", col="cornsilk4", cex=1)
text(x=900000, y=5700000, "Russia", col="cornsilk4", cex=1)
text(x=850000, y=4500000, "N. Korea", col="cornsilk4", cex=0.8, adj = 1)
text(x=800000, y=4200000, "S. Korea", col="cornsilk4", cex=0.8, adj = 1)

dev.off()


###############################################################################
# 5 - get environmental data and project to utm
# citation for worldclim: Hijmans, R.J., S.E. Cameron, J.L. Parra, P.G. Jones and A. Jarvis, 2005. Very high resolution interpolated climate surfaces for global land areas. International Journal of Climatology 25: 1965-1978.

#comment lines below if running for second time. BLICM already saved
#dir.create('data')
#BClim = getData("worldclim", var="bio", res=5, path="data/")
#BClim = crop(BClim, mapborders)
#BClim = projectRaster(BClim, res = c(8000,8000), crs="+proj=utm +zone=52N +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0") #value of res is sort of arbitrary: projection was creating stack with resolution c(7360, 9270), and having different resolutions in x and y was resulting in trouble
#BClim = crop(BClim,extent(utm_locs) + c(-1,1,-1,1)*501000) #crop data to 500km beyond ocurrence localities
#BClim = stack(BClim)
#area_to_mask = readOGR('./islands_exclude/islands.kml','islands') #mask islands from environmental data
#area_to_mask = spTransform(area_to_mask,"+proj=utm +zone=52N +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0")
#BClim = stack(mask(BClim,area_to_mask,inverse=T))
#writeRaster(BClim,'data/bclim_utm52_cropped.grd', overwrite=T)

#to avoid computations above, uncomment below to read saved raster
BClim = stack('data/bclim_utm52_cropped.grd')


###############################################################################
# 6 - run models
# 6A - format data
biomod_data = BIOMOD_FormatingData(resp.var = utm_subs, #object with coordinates
                     expl.var = BClim, #object with explanatory variables
                     resp.name = "C_relictus", #name of species
                     PA.nb.rep = 1, #number of pseudo-absence replicates
                     PA.nb.absences = 5000,#number of pseudo-absence
                     PA.strategy = 'disk', #will select pseudo-absences from disks around occurences
                     PA.dist.min = 80000, #minimal distance for pseudoabsence is 80 km
                     PA.dist.max = 400000, #maximal distance for pseudoabsence is 400 km
                     na.rm = T) #remove data points with no climatic data
#check data
biomod_data
pdf('graphs/presences_and_absences.pdf')
plot(biomod_data)
dev.off()

#6B - set modeling options
#will leave all at the default, except for the path to maxent executable
biomod_options = BIOMOD_ModelingOptions(MAXENT.Phillips = list( path_to_maxent.jar = './maxent/',
                                                                memory_allocated = 512,
                                                                background_data_dir = 'default',
                                                                maximumbackground = 'default',
                                                                maximumiterations = 200,
                                                                visible = FALSE,
                                                                linear = TRUE,
                                                                quadratic = TRUE,
                                                                product = TRUE,
                                                                threshold = TRUE,
                                                                hinge = TRUE,
                                                                lq2lqptthreshold = 80,
                                                                l2lqthreshold = 10,
                                                                hingethreshold = 15,
                                                                beta_threshold = -1,
                                                                beta_categorical = -1,
                                                                beta_lqp = -1,
                                                                beta_hinge = -1,
                                                                betamultiplier = 1,
                                                                defaultprevalence = 0.5))
  

#6C - Compute models
# default values for models and evaluation metrics
biomod_modelsout <- BIOMOD_Modeling(data = biomod_data,
                                    models.options = biomod_options,
                                    Prevalence = 0.5,
                                    NbRunEval = 5,
                                    DataSplit = 80,
                                    VarImport = 5,
                                    SaveObj = T,
                                    rescal.all.models = F,
                                    do.full.models = F,
                                    models.eval.meth=c('ROC','TSS'),
                                    models = c("GLM","GAM","GBM","CTA","ANN","SRE","FDA","MARS","RF","MAXENT.Phillips"))
# tried to run all models, but Maxent.Tsuruoka fails to project. Ended up removing it from this step

###############################################################################
#7 - Evaluate and combine models

biomod_evaluations <- get_evaluations(biomod_modelsout)
biomod_evaluations
# It seems none of the models is particularly good.
# Selecting models with  TSS > 0.8

biomod_EM <- BIOMOD_EnsembleModeling(modeling.output = biomod_modelsout,
                                     chosen.models = 'all',
                                     em.by='all',
                                     eval.metric = c('TSS'),
                                     eval.metric.quality.threshold = c(0.75),
                                     prob.mean = T,
                                     prob.cv = T,
                                     prob.ci = T,
                                     prob.ci.alpha = 0.05,
                                     prob.median = T,
                                     committee.averaging = T,
                                     prob.mean.weight = T,
                                     prob.mean.weight.decay = 'proportional',
                                     models.eval.meth = c('ROC','TSS'))
                                    

biomod_EM_evaluations <- get_evaluations(biomod_EM)
biomod_EM_evaluations

##############################################################################
#8 - Project models to the present


biomod_present_projection <- BIOMOD_Projection(modeling.output = biomod_modelsout,
                                               new.env=BClim,
                                               proj.name = 'C_relictus_current',
                                               do.stack = T,
                                               keep.in.memory = T,
                                               output.format = '.grd',
                                               on_0_1000 = F)

predictions <- get_predictions(biomod_present_projection)

for (i in 1:dim(predictions)[3]){
  pdf(paste('graphs/projection_present_model_', i, '.pdf',sep = ""),
      width = 7*(extent(utm_map)@xmax-extent(utm_map)@xmin)/(extent(utm_map)@ymax-extent(utm_map)@ymin),
      height = 7)
  plot(predictions@layers[[i]], main = names(predictions)[i])
  points(utm_subs, col = 'red', cex = 0.5)
  plot(utm_map, add = T)
  dev.off()
}

##############################################################################
#8 - Project ensemble models to the present
biomod_present_EM_projection <- BIOMOD_EnsembleForecasting(EM.output = biomod_EM,
                                                           projection.output = biomod_present_projection,
                                                           proj.name = 'C_relictus_current',
                                                           on_0_1000 = F)

predictions_EM <- get_predictions(biomod_present_EM_projection)

for (i in 1:dim(predictions_EM)[3]){
  pdf(paste('graphs/projection_EM_present', i, '.pdf',sep = ""),
      width = 7*(extent(utm_map)@xmax-extent(utm_map)@xmin)/(extent(utm_map)@ymax-extent(utm_map)@ymin),
      height = 7)
  plot(predictions_EM@layers[[i]], main = names(predictions_EM)[i])
  points(utm_subs, col = 'red', cex = 0.5)
  plot(utm_map, add = T)
  dev.off()
}

##############################################################################
#8 - Get data for future

#comment below if running for second time, area already downloaded
#future_K_best = getData("CMIP5", var="bio", res=5, rcp=26, model='HD', year=50, path="data/")
#future_K_best = crop(future_K_best, mapborders)
#future_K_best = projectRaster(future_K_best, res = c(8000,8000), crs="+proj=utm +zone=52N +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0") 
#future_K_best = crop(future_K_best,extent(utm_locs) + c(-1,1,-1,1)*501000) #crop data to 500km beyond ocurrence localities
#future_K_best = stack(future_K_best)
#future_K_best = stack(mask(future_K_best,area_to_mask,inverse=T))
#names(future_K_best) = names(BClim)
#writeRaster(future_K_best,'data/hd26_2050_utm52_cropped.grd', overwrite=T)

#future_K_worst = getData("CMIP5", var="bio", res=5, rcp=85, model='HD', year=50, path="data/")
#future_K_worst = crop(future_K_worst, mapborders)
#future_K_worst = projectRaster(future_K_worst, res = c(8000,8000), crs="+proj=utm +zone=52N +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0") 
#future_K_worst = crop(future_K_worst,extent(utm_locs) + c(-1,1,-1,1)*501000) #crop data to 500km beyond ocurrence localities
#future_K_worst = stack(future_K_worst)
#future_K_worst = stack(mask(future_K_worst,area_to_mask,inverse=T))
#names(future_K_worst) = names(BClim)
#writeRaster(future_K_worst,'data/hd85_2050_utm52_cropped.grd', overwrite=T)

#future_J_best = getData("CMIP5", var="bio", res=5, rcp=26, model='MC', year=50, path="data/")
#future_J_best = crop(future_J_best, mapborders)
#future_J_best = projectRaster(future_J_best, res = c(8000,8000), crs="+proj=utm +zone=52N +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0") 
#future_J_best = crop(future_J_best,extent(utm_locs) + c(-1,1,-1,1)*501000) #crop data to 500Jm beyond ocurrence localities
#future_J_best = stack(future_J_best)
#future_J_best = stack(mask(future_J_best,area_to_mask,inverse=T))
#names(future_J_best) = names(BClim)
#writeRaster(future_J_best,'data/mc26_2050_utm52_cropped.grd', overwrite=T)

#future_J_worst = getData("CMIP5", var="bio", res=5, rcp=85, model='MC', year=50, path="data/")
#future_J_worst = crop(future_J_worst, mapborders)
#future_J_worst = projectRaster(future_J_worst, res = c(8000,8000), crs="+proj=utm +zone=52N +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0") 
#future_J_worst = crop(future_J_worst,extent(utm_locs) + c(-1,1,-1,1)*501000) #crop data to 500Jm beyond ocurrence localities
#future_J_worst = stack(future_J_worst)
#future_J_worst = stack(mask(future_J_worst,area_to_mask,inverse=T))
#names(future_J_worst) = names(BClim)
#writeRaster(future_J_worst,'data/mc85_2050_utm52_cropped.grd', overwrite=T)

#future_C_best = getData("CMIP5", var="bio", res=5, rcp=26, model='BC', year=50, path="data/")
#future_C_best = crop(future_C_best, mapborders)
#future_C_best = projectRaster(future_C_best, res = c(8000,8000), crs="+proj=utm +zone=52N +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0") 
#future_C_best = crop(future_C_best,extent(utm_locs) + c(-1,1,-1,1)*501000) #crop data to 500Cm beyond ocurrence localities
#future_C_best = stack(future_C_best)
#future_C_best = stack(mask(future_C_best,area_to_mask,inverse=T))
#names(future_C_best) = names(BClim)
#writeRaster(future_C_best,'data/bc26_2050_utm52_cropped.grd', overwrite=T)

#future_C_worst = getData("CMIP5", var="bio", res=5, rcp=85, model='BC', year=50, path="data/")
#future_C_worst = crop(future_C_worst, mapborders)
#future_C_worst = projectRaster(future_C_worst, res = c(8000,8000), crs="+proj=utm +zone=52N +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0") 
#future_C_worst = crop(future_C_worst,extent(utm_locs) + c(-1,1,-1,1)*501000) #crop data to 500Cm beyond ocurrence localities
#future_C_worst = stack(future_C_worst)
#future_C_worst = stack(mask(future_C_worst,area_to_mask,inverse=T))
#names(future_C_worst) = names(BClim)
#writeRaster(future_C_worst,'data/bc85_2050_utm52_cropped.grd', overwrite=T)

#uncomment below to read downloaded and formatted data:
future_K_best = stack('data/hd26_2050_utm52_cropped.grd')
future_J_best = stack('data/mc26_2050_utm52_cropped.grd')
future_C_best = stack('data/bc26_2050_utm52_cropped.grd')

future_K_worst = stack('data/hd85_2050_utm52_cropped.grd')
future_J_worst = stack('data/mc85_2050_utm52_cropped.grd')
future_C_worst = stack('data/bc85_2050_utm52_cropped.grd')
           

##############################################################################
#8 - Project models to future
#save here first, because there is some error going on
save.image("SDM_Callipogon.RData")
#load("SDM_Callipogon.RData")
for (scenario in list(list("future_K_best",future_K_best),
                      list("future_J_best",future_J_best),
                      list("future_C_best",future_C_best),
                      list("future_K_worst",future_K_worst),
                      list("future_J_worst",future_J_worst),
                      list("future_C_worst",future_C_worst))){
  scenario_name = scenario[[1]]
  scenario_data = scenario[[2]]
  
  biomod_proj <- BIOMOD_Projection(modeling.output = biomod_modelsout,
                                   new.env=scenario_data,
                                   proj.name = paste('C_relictus', scenario_name,sep="_"),
                                   do.stack = T,
                                   keep.in.memory = T,
                                   output.format = '.grd',
                                   on_0_1000 = F)
  
  predictions <- get_predictions(biomod_proj)
  
  for (i in 1:dim(predictions)[3]){
    pdf(paste('graphs/projection_',scenario_name,'_model_', i, '.pdf',sep = ""),
        width = 7*(extent(utm_map)@xmax-extent(utm_map)@xmin)/(extent(utm_map)@ymax-extent(utm_map)@ymin),
        height = 7)
    plot(predictions@layers[[i]], main = names(predictions)[i])
    points(utm_subs, col = 'red', cex = 0.5)
    plot(utm_map, add = T)
    dev.off()
  }
  
  ##############################################################################
  #8 - Project ensemble models to the future
  biomod_EM_projection <- BIOMOD_EnsembleForecasting(EM.output = biomod_EM,
                                                     projection.output = biomod_proj,
                                                     proj.name = paste('C_relictus', scenario_name,sep="_"),
                                                     on_0_1000 = F)
  
  predictions_EM <- get_predictions(biomod_EM_projection)
  
  for (i in 1:dim(predictions_EM)[3]){
    pdf(paste('graphs/projection_EM_',scenario_name,'_' , i, '.pdf',sep = ""),
        width = 7*(extent(utm_map)@xmax-extent(utm_map)@xmin)/(extent(utm_map)@ymax-extent(utm_map)@ymin),
        height = 7)
    plot(predictions_EM@layers[[i]], main = names(predictions_EM)[i])
    points(utm_subs, col = 'red', cex = 0.5)
    plot(utm_map, add = T)
    dev.off()
  }
}
save.image("SDM_Callipogon.RData")
