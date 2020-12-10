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
  adm@data$NAME_1 <- iconv(adm@data$NAME_1,from="UTF-8",to="ASCII//TRANSLIT")
  adm@data$NAME_1 <- case_when(adm@data$NAME_1 == 'Trans Nzoia' ~ 'Trans-Nzoia',
                               adm@data$NAME_1 == "Murang'a" ~ 'Muranga', 
                               TRUE ~ adm@data$NAME_1)
  if(county == 'Muranga'){
    adm$NAME_1[grep(pattern = 'Murang', x = adm$NAME_1)] <- county
  }
  if(county == 'Trans-Nzoia'){
    adm$NAME_1[grep(pattern = 'Trans Nzoia', x = adm$NAME_1)] <- county
  }
  if(country == 'Malawi'){
    
    mlw_t <- c('Dedza','Dowa','Kasungu','Lilongwe','Mchinji','Nkhotakota','Ntcheu','Ntchisi','Salima')
    ply_t <- which(adm$NAME_1 %in% mlw_t)
    
    adm$NAME_2 =  adm$NAME_1
    adm$NAME_1[ply_t] <- "Central Region"
    trat <- adm %>% sf::st_as_sf()
    trat_1 <- trat %>% group_by(NAME_1, NAME_2) %>%
      summarise(geometry = sf::st_union(geometry)) %>% ungroup()
    
    adm <- trat_1 %>% as_Spatial()
    
  }
  if(country %in% c('Kenya', 'Ghana')){
    
    adm <- adm %>% sf::st_as_sf()
    adm$NAME_2 =  adm$NAME_1
    
    adm <-  adm %>% group_by(NAME_1, NAME_2) %>%
      summarise(geometry = sf::st_union(geometry)) %>% ungroup() %>% as_Spatial()
    
  }
  if(country == 'Tunisia'){
    
    NW <- c('Beja', 'Le Kef', 'Siliana', 'Jendouba')
    Pos_NW <- which(adm$NAME_1 %in% NW)
    CW <- c('Sidi Bou Zid', 'Kairouan', 'Kasserine')
    Pos_CW <- which(adm$NAME_1 %in% CW)
    
    adm$NAME_2 =  adm$NAME_1
    adm$NAME_1[Pos_NW] <- "North West"
    adm$NAME_1[Pos_CW] <- "Central West"
    
    trat <- adm %>% sf::st_as_sf()
    trat_1 <- trat %>% group_by(NAME_1, NAME_2) %>%
      summarise(geometry = sf::st_union(geometry)) %>% ungroup()
    
    adm <- trat_1 %>% as_Spatial()
  }
  if(country == 'Vietnam'){
    adm <- readRDS("//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/data/shps/shps_from_R/Vietnam/gadm36_VNM_1_sp.rds")
    
    adm$NAME_2 =  adm$NAME_1
  }
  if(country == 'Ethiopia'){
    adm$NAME_1 =  adm$NAME_2
  }
  
  
  adm_c <- adm[adm$NAME_1 == county,]
  
  alt_c <- alt %>% raster::crop(., adm_c) %>% raster::mask(., adm_c)
  
  # pp <- tm_shape(alt_c) + tm_raster(title = "Elevation(m)", palette = terrain.colors(50), style="cont",
  #             legend.is.portrait = TRUE) + tm_shape(adm_c) + tm_borders(lwd = 0.8, col = "gray50") +
  #   tm_text("NAME_2", size = 2, col = "gray30") + tm_layout(frame = FALSE, legend.title.size = 2.5,  
  #   legend.text.size = 2, legend.position = c("left","top"), legend.bg.color = "white")
  
  # tmap_save(pp, paste0(root,'/results/',country,'/graphs/',tolower(county),'/',tolower(county),'_elevation.png'), width = 10, height = 10, units = "in", dpi = 300)
  
  alt_c <- alt_c %>% rasterToPoints() %>% as_tibble() 
  xlims <- sf::st_bbox(adm_c)[c(1, 3)]
  ylims <- sf::st_bbox(adm_c)[c(2, 4)]
  adm_c <- adm_c %>% sf::st_as_sf()
  test <- as_tibble(st_centroid(adm_c) %>% st_coordinates())  %>%
    mutate(name =  iconv(adm_c$NAME_2,from="UTF-8",to="ASCII//TRANSLIT")) 
  
  glwd1 <- raster::shapefile('//dapadfs/workspace_cluster_8/climateriskprofiles/data/shps/GLWD/glwd_1.shp' ) 
  crs(glwd1) <- crs(adm_c)
  glwd1 <- glwd1 %>% sf::st_as_sf()
  
  
  glwd2 <- raster::shapefile('//dapadfs/workspace_cluster_8/climateriskprofiles/data/shps/GLWD/glwd_2.shp' ) 
  crs(glwd2) <- crs(adm_c)
  glwd2 <- glwd2 %>% sf::st_as_sf()
  
  
  pp <-  ggplot() +
    geom_tile(data = alt_c %>% setNames(c('x', 'y', 'alt')), aes(x = x, y = y, fill =  alt)) +
    geom_sf(data = adm_c, fill = NA, color = gray(.2)) +
    geom_sf(data = glwd1, fill = 'lightblue', color = 'lightblue') +
    geom_sf(data = glwd2, fill = 'lightblue', color = 'lightblue') +
    coord_sf(xlim = round(xlims, 2), ylim = round(ylims, 2)) +
    scale_fill_gradientn(colours = terrain.colors(20),  
                         guide = guide_colourbar(barheight = 12 ,
                                                 barwidth = 2, label.theme = element_text(size = 35))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    labs(fill = glue::glue('Elevation (m)'), x = 'Longitude', y = 'Latitude') + # title = county ,
    ggrepel::geom_label_repel(data=test, aes(x=X, y=Y, label=name), 
                              # arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "first"),
                              force = 10, 
                              size = 12) +
    theme_bw() + theme(text = element_text(size=35), 
                       legend.title=element_text(size=35), 
                       legend.spacing = unit(5, units = 'cm'),
                       legend.spacing.x = unit(1.0, 'cm'), plot.title = element_text(hjust = 0.5)) 
  
  ggsave(paste0(root,'/results/',country,'/graphs/',tolower(county),'/',tolower(county),'_elevation.png'), width = 12, height = 12,  dpi = 300)
  
  
  return('Map saved successfully\n')
  
}

cnty_list <- c('Brong Ahafo',
               'Ashanti',
               'Eastern',
               'Volta')
for(i in 1:length(cnty_list))
{
  getAltitude(iso3 = 'GHA', country = 'Ghana', county = cnty_list[i])
}