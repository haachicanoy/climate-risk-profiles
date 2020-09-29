# Obtain elevation map
# H. Achicanoy, A. Ghosh, & A. Esquivel
# Alliance Bioversity-CIAT, 2020

# R options
options(warn = -1, scipen = 999)

# Load libraries
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, tmap, fst))

# Paths
OSys <- Sys.info()[1]
root <<- switch(OSys,
                'Linux'   = '/home/jovyan/work/cglabs',
                'Windows' = '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles')

getAltitude <- function(iso3 = 'KEN', country = 'Kenya', county = 'Vihiga')
{
  
  alt <- raster::getData('alt', country = iso3, path = paste0(root,'/data/shps/',country))
  adm <- raster::getData('GADM', country = iso3, level = 2, path = paste0(root,'/data/shps/',country))
  if(county == 'Muranga'){
    adm$NAME_1[grep(pattern = 'Murang', x = adm$NAME_1)] <- county
  }
  if(county == 'Trans-Nzoia'){
    adm$NAME_1[grep(pattern = 'Trans Nzoia', x = adm$NAME_1)] <- county
  }
  adm_c <- adm[adm$NAME_1 == county,]
  
  alt_c <- alt %>% raster::crop(., adm_c) %>% raster::mask(., adm_c)
  
  pp <- tm_shape(alt_c) +
    tm_raster(title = "Elevation(m)", palette = terrain.colors(50), style="cont",
              legend.is.portrait = TRUE) +
    tm_shape(adm_c) + 
    tm_borders(lwd = 0.8, col = "gray50") +
    tm_text("NAME_2", size = 2, col = "gray30") +
    tm_layout(frame = FALSE,
              legend.title.size = 2.5,
              legend.text.size = 2,
              legend.position = c("left","top"),
              legend.bg.color = "white")
  
  tmap_save(pp, paste0(root,'/results/',country,'/graphs/',tolower(county),'/',tolower(county),'_elevation.png'), width = 10, height = 10, units = "in", dpi = 300)
  
  return('Map saved successfully\n')
  
}

cnty_list <- c('Bungoma','Kiambu','Kirinyaga','Kisii','Kitui','Migori','Muranga','Nandi','Narok','Nyamira','Samburu','Trans-Nzoia','Turkana','Vihiga')
for(i in 1:length(cnty_list))
{
  getAltitude(iso3 = 'KEN', country = 'Kenya', county = cnty_list[i])
}
