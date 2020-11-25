# Hazards. 
# H. Achicanoy & A. Esquivel
# Alliance Bioversity-CIAT
# Nov - 2020. 

rm(list = ls()); gc(reset = TRUE)

# =--------------------
# Packages 
options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyr, dplyr, tibble, ggplot2, raster, ncdf4, sf, lubridate, glue, cowsay, fst, ggspatial, vroom, sp, compiler, FactoMineR, factoextra, furrr, future))

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

# =----------------------------------------------------
path <- '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/'

# Prueba
reading_data <- function(list_p){
  time <- list_p$time ; gcm <- list_p$gcm ; period  <- list_p$period
  
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

# =------------------------------------------------------------------------------------
parameters <- tibble(time = 'past', gcm = '---', period = '1985_2015')

data_lect <- parameters %>% group_split(row_number()) %>%  
  purrr::map(.f = reading_data) %>% purrr::map(.f = as_tibble) %>% bind_rows()  %>% 
  dplyr::select(-gcm, -period, -time)

library(trend)

slope_f <- function(x){
  slope <- trend::sens.slope(x = data_lect$P5D, conf.level = 0.95)
  #p_value = as.numeric(slope$p.value)
  # sl <- tibble(slope = as.numeric(slope$estimates),  p_value)
  slope = as.numeric(slope$estimates)
  return(slope)}

lect1 <- data_lect %>% 
  dplyr::select(id, x, y, ISO3, Country, county, season, CDD:ndws) %>% 
  group_by(id, x, y, ISO3, Country, county, season) %>% 
  summarise_all(.funs = funs(mean, sd, slope_f)) 

tictoc::tic()
prueba_1 <- data_lect %>% 
  dplyr::select(id, x, y, ISO3, Country, county,  gSeason, SLGP, LGP) %>% 
  group_by(id, x, y, ISO3, Country, county, gSeason) %>% 
  summarise_all(.funs = funs(mean, sd))
tictoc::toc()



# =- 

ncores <- 12 # Modify this part if you want run in a server.  
plan(cluster, workers = ncores)

prueba_1a <- data_lect %>% 
  dplyr::select(id, gSeason, SLGP, LGP) %>% # filter(id == 223024, gSeason == 1) %>% 
  group_split(id, gSeason) %>% 
  furrr::future_map(.f = function(x){group_by(x, id, gSeason) %>% 
      summarise(SLGP_slope_f = slope_f(SLGP), LGP_slope_f = slope_f(LGP) )})


lect2 <- prueba_1a %>% bind_rows() %>% full_join(prueba_1 , .)


mutate(lect1 , season = stringr::str_remove(season, 's') %>% as.numeric())

lect1 %>% mutate(season = stringr::str_remove(season, 's'))


full_test <- full_join(mutate(lect1 , season = stringr::str_remove(season, 's') %>% as.numeric()), 
                       rename(lect2, season = 'gSeason') %>% dplyr::select(-x, -y) ) %>% ungroup()

# full_test %>% dplyr::select(ends_with('_slope_f'))


# =----------------------------------------------------------------------- PCA

numeric <- full_test %>% dplyr::select(-id, -x, -y, -ISO3, -Country, -county, -season)
corrplot::corrplot.mixed(corr = numeric%>%
                           cor(use = 'pairwise.complete.obs'))

normalization <- function(x){
  y <- (x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T))
  return(y)
}


numeric_norm <- numeric %>% apply(X = ., MARGIN = 2, FUN = normalization) %>% as_tibble()

corrplot::corrplot.mixed(corr = numeric_norm %>%
                           cor(use = 'pairwise.complete.obs'))

clm_pca <- numeric %>% FactoMineR::PCA(X = ., scale.unit = T, ncp = ncol(.), graph = F)

fviz_screeplot(clm_pca, ncp = 4)
fviz_pca_var(clm_pca, col.var = "contrib")

coord <- clm_pca$ind$coord[, 1:2] %>% as_tibble() %>% 
  bind_cols(dplyr::select(full_test, id, x, y, ISO3, Country, county, season), .)

tbl <- coord 
cSeasons_idcs <- c('Dim.1','Dim.2')

cat('>>> Obtain raster for all coordinates of big county\n')
r <- raster::rasterFromXYZ(xyz = all_climate %>% dplyr::select(x, y, id))
raster::crs(r) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
r.empty <- r
r.empty[] <- NA


# =------
c_idx <- unique(tbl$season) %>%
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
    
    idx_df          <- all_climate %>% dplyr::select(x, y, id)
    idx_df$ISO3     <- iso3c
    idx_df$Country  <- country
    idx_df$county   <- county
    idx_df          <- cbind(idx_df, raster::extract(surfaces, idx_df %>% dplyr::select(x, y) %>% data.frame)) %>% 
      as_tibble()
    idx_df$season   <- i
    return(idx_df) }) %>% bind_rows


ggplot() + geom_tile(data = c_idx, aes(x = x, y =  y, fill = Dim.1)) + theme_bw()

ggplot() + geom_tile(data = c_idx, aes(x = x, y =  y, fill = Dim.2)) + theme_bw()

