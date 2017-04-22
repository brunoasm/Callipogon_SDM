# This script makes supplementary figures
###############################################################################

library(raster)
library(ggplot2)
library(rgeos)
library(rgdal)

load('SDM_Callipogon.RData')

#1 - plot presence and pseudo-absence
map.points = fortify(utm_map)
pres_pa = cbind(biomod_data@coord, list('type' = c(rep('occurrence',52),rep('pseudo-absence',5000))))

pdf('sup_data_and_pa.pdf')
ggplot(pres_pa) +
  geom_polygon(aes(long/1000,lat/1000,group=group), data = map.points, fill = 'cornsilk2', colour = 'black') +
  geom_point(aes(x = longitude/1000, y = latitude/1000, color = type), size = 0.5) +
  coord_fixed(ratio = 1, xlim =c(-1600.197, 1528.176), ylim = c(3563.031,6386.691)) +
  theme(legend.position = 'bottom', panel.grid = element_blank(), panel.background = element_rect(fill='white')) +
  labs(title = 'Occurrence and Pseudo-Absence data', x = 'Longitude (km)', y = 'Longitude (km)')
dev.off()



#2 - plot ensemble projections for all scenarios
#preparing map polygons to ggplot
ext_latlon = extent(stack('C.relictus/proj_C_relictus_ll__present/proj_C_relictus_ll__present_C.relictus_ensemble.grd'))
sm_ll_map = crop(ll_map, ext_latlon + c(-10,10,-10,10))
map.points = fortify(sm_ll_map)


#read weighted average ensemble projection for present, removing islands
#area_to_mask = readOGR('./islands_exclude/islands.kml','islands') #mask islands from environmental data
proj_EM_ll = stack('C.relictus/proj_C_relictus_ll__present/proj_C_relictus_ll__present_C.relictus_ensemble.grd')
#convert raster to points for plotting weighted mean
plot_data = data.frame(rasterToPoints(proj_EM_ll[[c(6,7)]]))
names(plot_data) = c('Longitude','Latitude','Committee_Averaging','Weighted_Mean')
plot_data[[c('RCP')]] = 'current'
plot_data[[c('Model')]] = 'current'
plot_data = rbind(cbind(plot_data[c('Longitude','Latitude','RCP','Model')], data.frame(Ensemble_model = 'Committee Averaging', Probability = plot_data$Committee_Averaging)),
                  cbind(plot_data[c('Longitude','Latitude','RCP','Model')], data.frame(Ensemble_model = 'Weighted Mean', Probability = plot_data$Weighted_Mean)))

#make a graph showing committee averaging, weighted mea, present and future
for (scenario in paste(c('C','J','K','C','J','K'), c('best','worst'), sep = '_')){
  filename = paste('C.relictus/proj_C_relictus_ll__future_',scenario, '/proj_C_relictus_ll__future_', scenario, '_C.relictus_ensemble.grd', sep = '')
  rcp = c('2.6','8.5')[grepl('worst',scenario) + 1]
  if (grepl('K_', scenario)){mod_group = 'HadGEM2-AO'} else if(grepl('J_', scenario)) {mod_group = 'MIROC5'}  else if(grepl('C_', scenario)) {mod_group = 'BCC-CSM1-1'}
  proj_EM_ll_fut = stack(filename)
  plot_data_future = data.frame(rasterToPoints(proj_EM_ll_fut[[c(6,7)]]))
  names(plot_data_future) = c('Longitude','Latitude','Committee_Averaging','Weighted_Mean')
  plot_data_future[[c('RCP')]] = rcp
  plot_data_future[[c('Model')]] = mod_group
  plot_data_future = rbind(cbind(plot_data_future[c('Longitude','Latitude','RCP','Model')], data.frame(Ensemble_model = 'Committee Averaging', Probability = plot_data_future$Committee_Averaging)),
                           cbind(plot_data_future[c('Longitude','Latitude','RCP','Model')], data.frame(Ensemble_model = 'Weighted Mean', Probability = plot_data_future$Weighted_Mean)))
  plot_data = rbind(plot_data, plot_data_future)
}

plot_data$RCP = factor(plot_data$RCP, levels = c('current','2.6','8.5'), ordered = T)
plot_data$Model = factor(plot_data$Model, levels = c("current","MIROC5","BCC-CSM1-1","HadGEM2-AO"), ordered = T)


pdf('sup_map_projections.pdf',width = 7, height = 10)
base_plot = ggplot(plot_data) + 
  geom_raster(aes(Longitude, Latitude, fill=Probability)) +
  geom_path(aes(long,lat,group=group), data = map.points) +
  coord_fixed(ratio = 1.367336*24/30, xlim =c(110, 140), ylim = c(31,55)) +
  scale_fill_gradientn(colours = c('#4575b4','#ffffbf', '#d73027')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "ghostwhite")) +
  geom_point(aes(x = longitude, y = latitude), data = as.data.frame(spTransform(utm_subs, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "))), pch = 21, colour = 'black', fill = "cornsilk2", size = 0.8)
print(base_plot + facet_grid(Model + RCP ~ Ensemble_model, labeller = labeller(RCP = label_both, Ensemble_model = label_value, Model = label_value)))
dev.off()


##function that I saw on stackoverflow to calculate aspect ratio - not sure if I have to multiply result by 24/30, but looks good
map_aspect = function(x, y) {
  x.center <- sum(range(x)) / 2
  y.center <- sum(range(y)) / 2
  x.dist <- ggplot2:::dist_central_angle(x.center + c(-0.5, 0.5), rep(y.center, 2))
  y.dist <- ggplot2:::dist_central_angle(rep(x.center, 2), y.center + c(-0.5, 0.5))
  y.dist / x.dist
}

map_aspect(c(110, 140),c(31,55))
