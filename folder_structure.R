# Download shp and save. 
path_shp <- glue::glue('{root}/Project/{iso3c}/adm_lvl_{adm_lvl}/Study_region')
if(!dir.exists(path_shp)){dir.create(path_shp, recursive = T)}

shp_data <- raster::getData('GADM', country= iso3c, level=adm_lvl, path =  glue::glue('{path_shp}/'))

# in this part, the names of the countries of interest are arranged.
shp_data@data[, glue::glue('NAME_{adm_lvl}')] <-
  dplyr::pull(shp_data@data, glue::glue('NAME_{adm_lvl}')) %>%
  stringr::str_replace_all(., ' ', '-') %>%
  stringr::str_remove('[("!\"#$%()*,.:;<=>@[\\]^`{|}~.\']]') %>%
  stringi::stri_trans_general(., "Latin-ASCII")

# shp_data@data  <- shp_data@data %>% dplyr::select(GID_0, NAME_0, glue::glue('NAME_{1:adm_lvl}'))
rgdal::writeOGR(shp_data,dsn = path_shp,  glue::glue('{iso3c}_{adm_lvl}'), driver="ESRI Shapefile") 

# Erasing .rds format. 
rm(shp_data)
file.remove(glue::glue('{root}/Project/{iso3c}/adm_lvl_{adm_lvl}/Study_region/gadm36_BFA_{adm_lvl}_sp.rds'))

# create all folders structure. 
folder_names <- c('/Raw_data/Fut_climate/', '/Raw_data/His_climate/', 
                  '/Raw_data/Soil/', '/Procesed_data/', '/Results/Tables/', 
                  '/Results/Future/', '/Results/His_indices/', 
                  '/Graphs/Maps/', '/Graphs/Time_series/')

for(i in 1:length(folder_names)){
 dir.create(glue::glue('{root}Project/{iso3c}/adm_lvl_{adm_lvl}{folder_names[i]}'), recursive = TRUE)
}
rm(i, folder_names)




