# Load climate files per county
# By: H. Achicanoy & A. Esquivel
# CIAT, 2020

options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, vroom, sp, fst))

# Please check the correct county name within: Country_Counts.xlsx
get_observational_data <- function(country = 'Pakistan',
                                   county  = 'Kashmore',
                                   iso3    = 'PAK',
                                   adm_lvl = 3){
  
  country <<- country
  county  <<- county
  # Paths
  OSys <- Sys.info()[1]
  root <<- switch(OSys,
                  'Linux'   = '/home/jovyan/work/cglabs',
                  'Windows' = '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles')
  
  # Scripts
  source(paste0(root,'/scripts/win_parallelization.R'))
  
  # Load county shapefile
  shp <<- raster::shapefile(paste0(root,'/data/shps/',country,'/',iso3,'_adm',adm_lvl,'.shp'))
  glue::glue('shp <- shp[shp@data$NAME_{adm_lvl} == county,]') %>%
    as.character %>%
    parse(text = .) %>%
    eval(expr = ., envir = .GlobalEnv)
  
  # Load id coords
  crd <<- vroom::vroom(paste0(root,'/data/id_all_country.csv'), delim = ',')
  crd <<- crd %>%
    dplyr::filter(Country == country)
  pnt <<- crd %>% dplyr::select('x','y') %>% sp::SpatialPoints(coords = .)
  raster::crs(pnt) <<- raster::crs(shp)
  # Filter coordinates that are present in the county
  pnt <<- sp::over(pnt, shp) %>% data.frame %>% dplyr::select('ISO') %>% complete.cases() %>% which()
  crd <<- crd[pnt,]
  crd <<- crd
  
  # His obs
  out  <- paste0(root,'/data/observational_data/',tolower(country),'/',tolower(county),'_prec_temp.fst')
  out2 <- paste0(root,'/data/observational_data/',tolower(country),'/',tolower(county),'.fst')
  if(!dir.exists(dirname(out))){dir.create(dirname(out), recursive = T)}
  if(!file.exists(out)){
    
    cat(paste0('Precipitation and temperature data file does not exists for county: ',county,'. Creating it ...\n'))
    
    if(country %in% c('Ivory_Coast','Ghana')){
      ts <- list.files(path=paste0(root,'/data/Climate_other_C'),full.names=F,pattern='*.fst$') %>% sort()
    } else {
      ts <- list.files(path=paste0(root,'/data/Chirps_Chirts'),full.names=F,pattern='*.fst$') %>% sort()
    }
    cl <- createCluster(30, export = list("ts","root","crd"), lib = list("tidyverse","fst"))
    
    cat('>>> Loading and organizing time series\n')
    
    temp_prec <- ts %>% parallel::parLapply(cl, ., function(i){
      df <- fst::read_fst(paste0(root,'/data/Chirps_Chirts/',i))
      df <- df[df$id %in% crd$id,]
      return(df)
    }) %>%
      do.call(rbind, .)
    parallel::stopCluster(cl)
    temp_prec2 <- temp_prec %>% dplyr::group_by(id) %>% dplyr::arrange(Date) %>% dplyr::group_split(id)
    
    his_obs <- tibble::tibble(id      = crd$id,
                              x       = crd$x,
                              y       = crd$y,
                              ISO3    = crd$ISO3,
                              Country = crd$Country,
                              Climate = temp_prec2 %>% purrr::map(function(tb){tbl <- tb %>% dplyr::select('id','Date','prec','tmax','tmin'); return(tbl)}))
    his_obs <- his_obs %>% tidyr::unnest(.)
    outDir <- paste0(root,'/data/observational_data/',tolower(country))
    if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}
    out <- paste0(outDir,'/',tolower(county),'_prec_temp.fst')
    
    cat('>>> Saving climate data file\n')
    fst::write_fst(his_obs, out)
    
  } else {
    
    cat(paste0('Precipitation and temperature data file exists for county: ',county,'. Reading it ...\n'))
    his_obs <- fst::read_fst(out)
    his_obs <- his_obs %>%
      tidyr::nest(Climate = c('id','Date','prec','tmax','tmin')) %>%
      dplyr::rename(id = 'id1') %>%
      dplyr::select(id, everything(.))
    
    if(!file.exists(out2)){
      
      cat(paste0('Precipitation, temperature, and solar radiation data file does not exists for county: ',tolower(county),'. Creating it ...\n'))
      
      srad_fls <- paste0(root,"/data/NASA/",crd$id,'.fst')
      if(sum(file.exists(srad_fls)) == nrow(crd)){
        
        cat('>>> Adding solar radiation time series\n')
        srad <- srad_fls %>% purrr::map(.f = function(x){
          df <- fst::read_fst(x)
        })
        all_clim <- purrr::map2(.x = his_obs$Climate, .y = srad, .f = function(x, y){
          z <- dplyr::left_join(x = x, y = y %>% dplyr::select('id', 'Date', 'srad'), by = c('id','Date')) # Previous: srad
          return(z)
        })
        his_obs$Climate <- all_clim
        his_obs <- his_obs %>% tidyr::unnest(.)
        
        cat('>>> Saving climate data file\n')
        fst::write_fst(his_obs, out2)
        
      } else {
        cat('Not all SRAD files are available. Skipping reading ...\n')
      }
    } else {
      cat('Final file already exists\n')
    }
    
  }
  
  return(cat('Process successfully finished\n'))
  
}
# Run twice
for(i in 1:2){get_observational_data(country='India',county='Andhra Pradesh',iso3='IND',adm_lvl=1)}
