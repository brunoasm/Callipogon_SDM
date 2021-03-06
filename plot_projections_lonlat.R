# This script takes models generated by SDM_Callipogon.R and projects them to latlon projection to make figures
###############################################################################

# 1 - load packages
library(biomod2)
library(rworldmap)
library(rgdal)
library(rgeos)

load("SDM_Callipogon.RData")

predictions_present = stack('C.relictus/proj_C_relictus_current/proj_C_relictus_current_C.relictus.grd')
ext_latlon = extent(projectRaster(predictions_present,crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "))
ll_subs = spTransform(utm_subs,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  
sm_ll_map = crop(ll_map, ext_latlon + c(-10,10,-10,10))


ll_BClim = stack(crop(getData("worldclim", var="bio", res=5, path="data/"),ext_latlon))
ll_future_K_best = stack(crop(getData("CMIP5", var="bio", res=5, rcp=26, model='HD', year=50, path="data/"), ext_latlon))
ll_future_J_best = stack(crop(getData("CMIP5", var="bio", res=5, rcp=26, model='MC', year=50, path="data/"), ext_latlon))
ll_future_C_best = stack(crop(getData("CMIP5", var="bio", res=5, rcp=26, model='BC', year=50, path="data/"), ext_latlon))
ll_future_K_worst = stack(crop(getData("CMIP5", var="bio", res=5, rcp=85, model='HD', year=50, path="data/"), ext_latlon))
ll_future_J_worst = stack(crop(getData("CMIP5", var="bio", res=5, rcp=85, model='MC', year=50, path="data/"), ext_latlon))
ll_future_C_worst = stack(crop(getData("CMIP5", var="bio", res=5, rcp=85, model='BC', year=50, path="data/"), ext_latlon))

#now redo projection of ensemble models
for (scenario in list(list("present",ll_BClim),
                      list("future_K_best",ll_future_K_best),
                      list("future_J_best",ll_future_J_best),
                      list("future_C_best",ll_future_C_best),
                      list("future_K_worst",ll_future_K_worst),
                      list("future_J_worst",ll_future_J_worst),
                      list("future_C_worst",ll_future_C_worst))){
  scenario_name = scenario[[1]]
  scenario_data = scenario[[2]]
  names(scenario_data) = names(ll_BClim) #cmip5 models use different names for bioclim variables, this causes an error
  
  biomod_proj <- BIOMOD_Projection(modeling.output = biomod_modelsout,
                                   new.env=scenario_data,
                                   proj.name = paste('C_relictus_ll_', scenario_name,sep="_"),
                                   do.stack = T,
                                   keep.in.memory = T,
                                   output.format = '.grd',
                                   on_0_1000 = F)
  
  biomod_EM_projection <- BIOMOD_EnsembleForecasting(EM.output = biomod_EM,
                                                     projection.output = biomod_proj,
                                                     proj.name = paste('C_relictus_ll_', scenario_name,sep="_"),
                                                     on_0_1000 = F)
  
  predictions_EM <- get_predictions(biomod_EM_projection)
  
  for (i in 1:dim(predictions_EM)[3]){
    pdf(paste('graphs/projection_EM_ll_',scenario_name,'_' , i, '.pdf',sep = ""),
        width = 7*(ext_latlon@xmax-ext_latlon@xmin)/(ext_latlon@ymax-ext_latlon@ymin),
        height = 7)
    plot(predictions_EM@layers[[i]], main = names(predictions_EM)[i])
    points(ll_subs, col = 'red', cex = 0.5)
    plot(ll_map, add = T)
    dev.off()
  }
}
