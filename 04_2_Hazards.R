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
  slope <- trend::sens.slope(x = x, conf.level = 0.95)
  #p_value = as.numeric(slope$p.value)
  # sl <- tibble(slope = as.numeric(slope$estimates),  p_value)
  slope = as.numeric(slope$estimates)
  return(slope)}

lect1 <- data_lect %>% 
  dplyr::select(id, x, y, ISO3, Country, county, season, CDD:ndws) %>% 
  group_by(id, x, y, ISO3, Country, county, season) %>% 
  summarise_all(.funs = funs(mean, sd, slope_f)) 

tictoc::tic()
lect2 <- data_lect %>% 
  dplyr::select(id, x, y, ISO3, Country, county,  gSeason, SLGP, LGP) %>% 
  group_by(id, x, y, ISO3, Country, county, gSeason) %>% 
  summarise_all(.funs = funs(mean, sd, slope_f))
tictoc::toc()



# =- 

full_test <- full_join(mutate(lect1 , season = stringr::str_remove(season, 's') %>% as.numeric()), 
                       rename(lect2, season = 'gSeason') %>% dplyr::select(-x, -y) ) %>% ungroup()

# full_test %>% dplyr::select(ends_with('_slope_f'))
write.csv(x = full_test, file = glue::glue('{path}{country}/graphs/{county}/full_data.csv'))


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
ggsave(glue::glue('{path}{country}/graphs/{county}/bar.png'))
fviz_pca_var(clm_pca, col.var = "contrib")
ggsave(glue::glue('{path}{country}/graphs/{county}/contrib.png'))



# =--- Pruebas de indices

pca_index <- function(pca = natural_pca, percent = 80)
{
  ncomp    <- which(pca$eig[,3] >= percent)[1]
  variance <- pca$eig[1:ncomp,2] %>% matrix(data = ., ncol = 1, byrow = F)
  raw_indx <- pca$ind$coord[,1:ncomp] %*% variance
  fnl_indx <- 0.1 + ((raw_indx-min(raw_indx))/(max(raw_indx)-min(raw_indx))*0.9)
  fnl_indx <- fnl_indx %>% as.numeric()
  return(fnl_indx)
}



# =-- Pruebas de interpolacion 
tbl <- bind_cols(dplyr::select(full_test, id, x, y, ISO3, Country, county, season), 
                 ind = pca_index(clm_pca, percent = 70)) 
cSeasons_idcs <- c('ind')

cat('>>> Obtain raster for all coordinates of big county\n')
r <- raster::rasterFromXYZ(xyz = all_climate %>% dplyr::select(x, y, id))
raster::crs(r) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
r.empty <- r
r.empty[] <- NA


# =------ Interpolation
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


country1 <- country1 %>% sf::st_as_sf()
shp_sf <- shp  %>% sf::st_as_sf()
xlims <- sf::st_bbox(shp_sf)[c(1, 3)]
ylims <- sf::st_bbox(shp_sf)[c(2, 4)]
map_world <- raster::shapefile(glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/data/shps/all_country/all_countries.shp')) %>% 
  sf::st_as_sf()



write.csv(x = c_idx, file = glue::glue('{path}{country}/graphs/{county}/Hazards.csv'))

ggplot() + 
  # geom_sf(data = map_world, fill = gray(.9), color = gray(.9)) +
  geom_sf(data = country1, fill = 'white', color = gray(.5)) +
  geom_tile(data = c_idx, aes(x = x, y =  y, fill = ind )) + 
  geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
  coord_sf(xlim = xlims, ylim = ylims) +
  scale_fill_gradient2(low = "#E8E8E8", mid = "#64ACBE", high = "#985356", 
                       guide = guide_colourbar(barwidth = 12, 
                                               label.theme = element_text(angle = 25, size = 15))) +
  labs(fill = glue::glue('Hazard'), title = 'Historic', x = 'Longitude', y = 'Latitude') +
  theme_bw() + theme(legend.position = 'bottom', text = element_text(size=15), 
                     legend.title=element_text(size=15), 
                     legend.spacing = unit(5, units = 'cm'),
                     legend.spacing.x = unit(1.0, 'cm'), plot.title = element_text(hjust = 0.5))

# #574249 --- # #64ACBE  ---- # #E8E8E8


ggsave(glue::glue('{path}{country}/graphs/{county}/Hazards.png') , width = 8, height = 5.5, dpi = 300)



# =----------------------------------------------------------

comp <- clm_pca$ind$coord[,1:2] %>% as_tibble() %>% bind_cols(dplyr::select(full_test, id, x, y), .)


Dq1 <- quantile(comp$Dim.1, c(0.33, 0.66))
Dq2 <- quantile(comp$Dim.2, c(0.33, 0.66))


Julian_graph <- comp %>% 
  mutate(D1_m =  dplyr::case_when(Dim.1 < Dq1[1] ~ 1, 
                                  Dim.1 > Dq1[2] ~ 3, 
                                  TRUE ~ 2), 
         D2_m =  dplyr::case_when(Dim.2 < Dq2[1] ~ 1, 
                                  Dim.2 > Dq2[2] ~ 3, 
                                  TRUE ~ 2)) %>% 
  mutate(all_two = dplyr::case_when(D1_m == 1 & D2_m == 1 ~ 1, 
                                    D1_m == 1 & D2_m == 2 ~ 2,
                                    D1_m == 1 & D2_m == 3 ~ 3,
                                    D1_m == 2 & D2_m == 1 ~ 4, 
                                    D1_m == 2 & D2_m == 2 ~ 5,
                                    D1_m == 2 & D2_m == 3 ~ 6,
                                    D1_m == 3 & D2_m == 1 ~ 7, 
                                    D1_m == 3 & D2_m == 2 ~ 8,
                                    D1_m == 3 & D2_m == 3 ~ 9,
                                    TRUE ~ NA_real_)) # %>% count(all_two)


# =--- Interpolar. 
tbl2 <- Julian_graph
cSeasons_idcs <- 'all_two' 

tbl_f <- 1 %>% purrr::map(.f = function(i){
  tbl2 <<- tbl2 %>% tidyr::drop_na()
  
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

tbl_f <- tbl_f %>% mutate(all_two = round(all_two, 0))



# =--- 


ggplot() + 
  geom_sf(data = country1, fill = 'white', color = gray(.5)) +
  geom_tile(data = tbl_f, aes(x = x, y =  y, fill = as.factor(all_two ))) + 
  geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
  scale_fill_manual(values = c("#E8E8E8", "#E4ACAC", "#C85A5A", "#B0D5DF", "#AD9EA5", "#985356",
                               "#64ACBE", "#627F8C", "#574249"), guide =) +
  coord_sf(xlim = xlims, ylim = ylims) +
  scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
  scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
  labs(fill = NULL, title = 'Hazards', x = 'Longitude', y = 'Latitude') +
  theme_bw() + theme(legend.position = 'none', text = element_text(size=15), 
                     legend.title=element_text(size=15), 
                     legend.spacing = unit(5, units = 'cm'),
                     legend.spacing.x = unit(1.0, 'cm'), plot.title = element_text(hjust = 0.5))
ggsave(glue::glue('{path}{country}/graphs/{county}/Hazards_biv.png') , width = 8, height = 5.5, dpi = 300)
