# Interpolation strategy for big counties - GIZ risk profiles
# H. Achicanoy & A. Esquivel
# Alliance Bioversity-CIAT

# R options
options(warn = -1, scipen = 999)

# Load libraries
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, gstat, rgdal, fst, glue))

# Paths
OSys <- Sys.info()[1]
root <<- switch(OSys,
                'Linux'   = '/home/jovyan/work/cglabs',
                'Windows' = '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles')

interpolation_strtg <- function(country = 'Zambia',
                                county  = 'Eastern',
                                iso3c   = 'ZMB',
                                gcm     = 'ipsl_cm5a_mr',
                                period  = '1985_2015',
                                time    = 'past'){
  
  # Load observational data
  clim_data <- fst::read_fst(paste0(root,'/data/observational_data/',tolower(country),'/',tolower(county),'.fst'))
  clim_data <- clim_data %>%
    tidyr::nest(Climate = c('id','Date','prec','tmax','tmin','srad')) %>%
    dplyr::rename(id = 'id1') %>%
    dplyr::select(id, everything(.))
  
  if(time == 'past'){
    # Load calculated indices for 30% of pixels
    indices_all <- fst::read_fst(paste0(root,'/results/',country,'/',time,'/',county,'_',period,'_corrected_idw.fst'))
    # Load calculated indices for 100% of pixels just for step 05
    #indices_prt <- fst::read_fst(paste0(root,'/results/',country,'/',time,'/',county,'_',period,'_prec_temp.fst'))
  } else {
    if(time == 'future'){
      # Load calculated indices for 30% of pixels
      indices_all <- fst::read_fst(paste0(root,'/results/',country,'/',time,'/',gcm,'/',period,'/',county,'_',period,'_corrected_idw.fst'))
      # Load calculated indices for 100% of pixels just for step 05
      # indices_prt <- fst::read_fst(paste0(root,'/results/',country,'/',time,'/',gcm,'/',period,'/',county,'_',period,'_prec_temp.fst'))
    }
  }
  
  cat('>>> Obtain raster for all coordinates of big county\n')
  r <- raster::rasterFromXYZ(xyz = clim_data %>% dplyr::select(x, y, id))
  raster::crs(r) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  r.empty <- r
  r.empty[] <- NA
  
  cat('>>> Obtain averages per growing season\n')
  tbl_ref <- indices_all  %>% 
    dplyr::select(id, x, y, ISO3, Country, county, year, season, CDD, P5D, P95, NT35, ndws, gSeason, SLGP, LGP)
  tbl_gSeasons <- tbl_ref %>%
    dplyr::group_by(id, gSeason) %>%
    dplyr::summarise(SLGP = mean(SLGP, na.rm = T),
                     LGP  = mean(LGP, na.rm = T))
  tbl_cSeasons <- tbl_ref %>%
    dplyr::mutate(season = factor(season)) %>%
    dplyr::group_by(id, season) %>%
    dplyr::summarise(CDD  = mean(CDD, na.rm = T),
                     P5D  = mean(P5D, na.rm = T),
                     P95  = mean(P95, na.rm = T),
                     NT35 = mean(NT35, na.rm = T),
                     ndws = mean(ndws, na.rm = T))
  tbl_cSeasons <- tbl_cSeasons %>%
    dplyr::filter(!is.na(season))
  tbl <- dplyr::full_join(x = tbl_cSeasons, y = tbl_gSeasons, by = 'id')
  rm(tbl_cSeasons, tbl_gSeasons)
  if(length(tbl$season %>% unique) == 2 & length(tbl$gSeason %>% unique) == 2){
    d1 <- tbl %>%
      dplyr::filter(season == 's1' & gSeason == 1)
    d2 <- tbl %>%
      dplyr::filter(season == 's2' & gSeason == 2)
    tbl <- dplyr::bind_rows(d1, d2)
  }
  
  tbl <- dplyr::left_join(x = tbl, y = tbl_ref %>% dplyr::select(id, x, y) %>% unique, by = 'id')
  tbl <- tbl %>% dplyr::select(id, x, y, dplyr::everything(.))
  if(sum(tbl$id %>% duplicated & tbl$gSeason %>% is.na) > 0){
    tbl <- tbl[-which(tbl$id %>% duplicated & tbl$gSeason %>% is.na),]
  }
  
  cSeasons_idcs <- c('CDD','P5D','P95','NT35','ndws')
  gSeasons_idcs <- c('gSeason','SLGP','LGP')
  
  cat('>>> Interpolate 70% of pixels ...\n')
  
  cseasons <- tbl$season %>% unique
  gseasons <- tbl$gSeason %>% unique
  if(length(gseasons) >= 3){
    gseasons <- gseasons[1:2]
  }
  
  cat('>>> Calculating interpolated surfaces for agro-climatic indices for climatology season indices\n')
  c_idx <- cseasons %>%
    purrr::map(.f = function(i){
      
      cat(paste0(' --- Filter indices per growing season: ',i,'\n'))
      tbl2 <<- tbl[which(tbl$season == i),] %>%
        tidyr::drop_na()
      
      cat(paste0(' --- Create SpatialDataFrame object\n'))
      spdf <<- sp::SpatialPointsDataFrame(coords      = tbl2[,c('x','y')] %>% data.frame,
                                          data        = tbl2 %>% data.frame,
                                          proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
      
      surfaces <- 1:length(cSeasons_idcs) %>%
        purrr::map(.f = function(j){
          
          cat(paste0(' --- Fit inverse distance weighted interpolation for: ',cSeasons_idcs[j],'\n'))
          glue::glue('idw_fit <- gstat::gstat(formula = {cSeasons_idcs[j]} ~ 1, locations = spdf)') %>%
            as.character %>%
            parse(text = .) %>%
            eval(expr = ., envir = .GlobalEnv)
          idw_int <- raster::interpolate(r.empty, idw_fit)
          idw_msk <- raster::mask(idw_int, r)
          return(idw_msk)
          
        })
      surfaces <- raster::stack(surfaces)
      names(surfaces) <- cSeasons_idcs; rm(spdf)
      
      idx_df          <- clim_data %>% dplyr::select(x, y, id)
      idx_df          <- idx_df[idx_df$id %in% base::setdiff(clim_data$id, indices_all$id),]
      idx_df$ISO3     <- iso3c
      idx_df$Country  <- country
      idx_df$county   <- county
      idx_df          <- cbind(idx_df, raster::extract(surfaces, idx_df %>% dplyr::select(x, y) %>% data.frame))
      idx_df$season   <- i
      return(idx_df)
      cat('\n')
      cat('\n')
      
    })
  c_idx       <- dplyr::bind_rows(c_idx)
  
  cat('>>> Calculating interpolated surfaces for agro-climatic indices for growing season indices\n')
  g_idx <- gseasons %>%
    purrr::map(.f = function(i){
      
      cat(paste0(' --- Filter indices per growing season: ',i,'\n'))
      tbl2 <<- tbl[which(tbl$gSeason == i),] %>%
        tidyr::drop_na()
      
      cat(paste0(' --- Create SpatialDataFrame object\n'))
      spdf <<- sp::SpatialPointsDataFrame(coords      = tbl2[,c('x','y')] %>% data.frame,
                                          data        = tbl2 %>% data.frame,
                                          proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
      
      surfaces <- 1:length(gSeasons_idcs) %>%
        purrr::map(.f = function(j){
          
          cat(paste0(' --- Fit inverse distance weighted interpolation for: ',gSeasons_idcs[j],'\n'))
          glue::glue('idw_fit <- gstat::gstat(formula = {gSeasons_idcs[j]} ~ 1, locations = spdf)') %>%
            as.character %>%
            parse(text = .) %>%
            eval(expr = ., envir = .GlobalEnv)
          idw_int <- raster::interpolate(r.empty, idw_fit)
          idw_msk <- raster::mask(idw_int, r)
          return(idw_msk)
          
        })
      surfaces <- raster::stack(surfaces)
      names(surfaces) <- gSeasons_idcs; rm(spdf)
      
      idx_df          <- clim_data %>% dplyr::select(x, y, id)
      idx_df          <- idx_df[idx_df$id %in% base::setdiff(clim_data$id, indices_all$id),]
      idx_df$ISO3     <- iso3c
      idx_df$Country  <- country
      idx_df$county   <- county
      idx_df          <- cbind(idx_df, raster::extract(surfaces, idx_df %>% dplyr::select(x, y) %>% data.frame))
      idx_df$gSeason  <- i
      return(idx_df)
      cat('\n')
      cat('\n')
      
    })
  g_idx       <- dplyr::bind_rows(g_idx)
  
  interpolated_indices <- dplyr::left_join(x = c_idx, y = g_idx %>% dplyr::select(id, gSeason, SLGP, LGP), by = 'id')
  if(length(interpolated_indices$season %>% unique) == 2 & length(interpolated_indices$gSeason %>% unique) == 2){
    d1 <- interpolated_indices %>%
      dplyr::filter(season == 's1' & gSeason == 1)
    d2 <- interpolated_indices %>%
      dplyr::filter(season == 's2' & gSeason == 2)
    interpolated_indices <- dplyr::bind_rows(d1, d2)
  }
  
  tbl$ISO3    <- iso3c
  tbl$Country <- country
  tbl$county  <- county
  
  tbl <- dplyr::bind_rows(interpolated_indices, tbl)
  tbl$CDD  <- round(tbl$CDD)
  tbl$NT35 <- round(tbl$NT35)
  tbl$ndws <- round(tbl$ndws)
  tbl$SLGP <- round(tbl$SLGP)
  tbl$LGP  <- round(tbl$LGP)
  
  tbl <- tbl %>% as_tibble() %>%
    dplyr::select(  x, y, id, ISO3,  Country, county, season, everything()) %>% 
    dplyr::arrange(id) %>% 
    mutate_at(vars(CDD:LGP), ~ifelse(. == 'NaN', mean(., na.rm = TRUE), .) %>% round()) %>% 
    mutate(ndws = ifelse(ndws == 0, median(ndws), ndws))
  
  
  
  
  if(time == 'past'){
    outFut <- paste0(root,'/results/',country,'/',time,'/',county,'_',period,'_idw.fst')
    if(!file.exists(outFut)){
      fst::write_fst(tbl, outFut)
    }
  } else {
    if(time == 'future'){
      outFut <- paste0(root,'/results/',country,'/',time,'/',gcm,'/',period,'/',county,'_',period,'_idw.fst')
      if(!file.exists(outFut)){
        fst::write_fst(tbl, outFut)
      }
    }
  }
  
  return(cat('Process finished successfully ...\n'))
  
}

counties <- c('Eastern', 
              'Southern')
for(cnty in counties){
  interpolation_strtg(country = 'Zambia',
                      county  = cnty,
                      iso3c   = 'ZMB',
                      gcm     = NULL,
                      period  = '1985_2015',
                      time    = 'past')
}
periods <- c('2021_2045','2041_2065')
gcmList <- c('ipsl_cm5a_mr','miroc_esm_chem','ncc_noresm1_m')
for(cnty in counties){
  for(p in periods){
    for(g in gcmList){
      interpolation_strtg(country = 'Zambia',
                          county  = cnty,
                          iso3c   = 'ZMB',
                          gcm     = g,
                          period  = p,
                          time    = 'future')
    }
  }
}