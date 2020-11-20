# Hazards. 
# H. Achicanoy & A. Esquivel
# Alliance Bioversity-CIAT
# Nov - 2020. 

rm(list = ls()); gc(reset = TRUE)

# =--------------------
# Packages 
options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyr, dplyr, tibble, ggplot2, raster, ncdf4, sf, lubridate, glue, cowsay, fst, ggspatial, vroom, sp, compiler))
# =--------------------
# =----------------------------------
# Identificacion de pixel para ETH
# =----------------------------------
country <- co <-  'Burkina_Faso'
count_i  <- C_shp <- county <-   c('Sud-Ouest')
iso3c <- 'BFA'
Big <- 'B'
adm_lvl <- 1



# Ruta Principal para guardados: 
root <- '//dapadfs/workspace_cluster_8/climateriskprofiles/'

# Load county shapefile
if(country == 'India'){
  # India 
  country1 <- raster::shapefile(glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/results/India/states/Admin2.shp'))
  shp <- raster::shapefile(glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/results/India/states/Admin2.shp'))
  shp <- shp[shp@data$ST_NM %in% C_shp,]
  plot(shp)
  shp@data$ISO <- iso3c
}else{
  country1 <- readRDS(glue::glue('{root}data/shps/shps_from_R/{country}/gadm36_{iso3c}_{adm_lvl}_sp.rds')) 
  shp <- readRDS(glue::glue('{root}data/shps/shps_from_R/{country}/gadm36_{iso3c}_{adm_lvl}_sp.rds')) 
  shp@data$NAME_1 <- iconv(shp@data$NAME_1,from="UTF-8",to="ASCII//TRANSLIT")
  shp@data$NAME_1 <- case_when(shp@data$NAME_1 == 'Trans Nzoia' ~ 'Trans-Nzoia',
                               shp@data$NAME_1 == "Murang'a" ~ 'Muranga', 
                               TRUE ~ shp@data$NAME_1)
  shp <- shp[shp@data$NAME_1 %in% C_shp,]
  plot(shp)
  shp@data$ISO <- iso3c
}

# Load id coords
crd <- vroom('//dapadfs/workspace_cluster_8/climateriskprofiles/data/id_all_country.csv')
crd <- crd %>%
  dplyr::filter(Country == country)
pnt <- crd %>% dplyr::select(x,y) %>% sp::SpatialPoints(coords = .)
crs(pnt) <- crs(shp)
# Filter coordinates that are present in the county
pnt <- sp::over(pnt, shp) %>% data.frame %>% dplyr::select(ISO) %>% complete.cases() %>% which()
crd <- crd[pnt,]
crd <<- crd


# =--------------------

tictoc::tic()
all_climate <-  fst::fst(glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/data/observational_data/{co}/{county}.fst')) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(-id1)  %>% 
  tidyr::nest(-id, -x, -y, -ISO3) %>% 
  dplyr::mutate(Country = country) %>% 
  dplyr::rename(climate = 'data') %>% 
  dplyr::select(id, x, y, ISO3, Country, climate)
tictoc::toc()


historic <-  all_climate %>%
  dplyr::mutate(summary =  purrr::map(.x = climate, .f = function(z){
    z <- z %>% dplyr::mutate(year = lubridate::year(Date), month = lubridate::month(Date) ,
                             tmean = (tmax + tmin)/2 )  %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(prec = sum(prec), tmean = mean(tmean)) %>%
      dplyr::ungroup() %>% dplyr::select(-year) %>%
      dplyr::summarise_all(.funs = function(x){round( mean(x, na.rm = TRUE), 1)}) %>%
      dplyr::ungroup()})) %>% dplyr::select(-climate) %>% tidyr::unnest()  %>%
  dplyr::group_by( id, x, y, ISO3, Country) %>%
  dplyr::summarise_all(~mean(.))

historic <- historic %>% filter(id %in% crd$id)



# =----------------------------------------------------
path <- '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/'

# Prueba


reading_data <- function(list_p){
  
  time <- list_p$time
  gcm     <- list_p$gcm
  period  <- list_p$period
  
  if(time == 'past'){
    # Load calculated indices for 30% of pixels
    indices_all <- fst::read_fst(paste0(root,'/results/',country,'/',time,'/',county,'_',period,'_corrected_idw.fst')) %>% 
      as_tibble()
  } else {
    if(time == 'future'){
      # Load calculated indices for 30% of pixels
      indices_all <- fst::read_fst(paste0(root,'/results/',country,'/',time,'/',gcm,'/',period,'/',county,'_',period,'_corrected_idw.fst'))
      as_tibble()
    }
  }
  
  data_index <- indices_all %>% 
    mutate(time = time, gcm = ifelse(time == 'future', gcm, time), 
           period = period) 
  return(data_index)}

# reading_data(list_p = tibble(time = 'past', gcm = 'ipsl_cm5a_mr', period = '1985_2015'))

parameters <- tibble(time = c('past', 'future', 'future', 'future', 'future', 'future', 'future'), 
                     gcm = c('---', 'ipsl_cm5a_mr','miroc_esm_chem','ncc_noresm1_m', 'ipsl_cm5a_mr','miroc_esm_chem','ncc_noresm1_m'), 
                     period = c('1985_2015', '2021_2045', '2021_2045', '2021_2045', '2041_2065', '2041_2065', '2041_2065'))

data_lect <- parameters %>% group_split(row_number()) %>%  
  purrr::map(.f = reading_data) %>% purrr::map(.f = as_tibble) %>% 
  bind_rows()



data_lect_basic <-  data_lect %>% 
  dplyr::select(id,x,y,ISO3,Country,county,season,CDD,P5D,P95,NT35, ndws, time, gcm, period) %>% 
  group_by(id,x,y,ISO3,Country,county, season, time, gcm, period) %>% 
  summarise_all(.funs = funs(mean, sd))


data_lect_basic %>% ungroup() %>% 
  dplyr::select(-gcm) %>% 
  group_by(id,x,y,ISO3,Country,county,season, time, period) %>% 
  summarise_all(.funs = funs(mean)) %>% 
  ungroup() %>% 
  dplyr::select(-period) %>% 
  group_by(id,x,y,ISO3,Country,county,season, time) %>% 
  summarise_all(.funs = funs(mean)) %>% 
  arrange(time)
