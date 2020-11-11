# Maps for 
# H. Achicanoy & A. Esquivel
# Alliance Bioversity-CIAT
# June - 2020. 

rm(list = ls())
gc(reset = TRUE)

# =--------------------
# Packages 
options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, tibble, raster, ncdf4, sf, lubridate, glue, cowsay, fst, ggspatial, vroom, sp, compiler))
# =--------------------


# =----------------------------------
# Identificacion de pixel para ETH
# =----------------------------------
country <- 'Nigeria'
county <-   c('Abia','Adamawa','Akwa Ibom','Anambra','Bauchi','Benue','Borno',
              'Cross River','Delta','Ebonyi','Edo','Ekiti','Enugu','FCT Abuja','Gombe','Imo',
              'Jigawa','Kaduna','Kano','Katsina','Kebbi','Kogi','Kwara','Lagos','Nassarawa','Niger',
              'Ogun','Ondo','Osun','Oyo','Plateau','Rivers','Sokoto','Taraba','Yobe','Zamfara', 'Bayelsa')
adm_lvl <- 1
iso3c <- 'NGA'
Big <- 'B'

chain <- TRUE # Este
# Una cadena a la vez por ahora. 
value_chain <- 'Wheat' #  Cassava Soybean Maize Cotton Wheat Sesame


# =---------------------------------------------------
# Ruta Principal para guardados: 
root <- '//dapadfs/workspace_cluster_8/climateriskprofiles/'

# Incertando todo el tema de cadena de valor. 
if(isTRUE(chain)){
  vc_tibble <- read.csv(paste0(root, '/NIRSAL/chain.csv')) %>% 
    as_tibble() %>% 
    dplyr::select(County, value_chain) 
  
  county <- vc_tibble[which(vc_tibble[,2] == 'X'), ] %>% pull(County) %>% as.character()
}else{
  print('All Country')
}


m_coords <- tibble(county = county) %>% 
  mutate(coords = purrr::map(.x = county, 
                             .f = function(x){
                               fst::fst(glue::glue('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/NIRSAL/{country}_{value_chain}/past/{x}_1985_2015_idw.fst')) %>% 
                                 as.tibble() %>% dplyr::select(id, x, y)})) %>% 
  unnest(coords)



# Load county shapefile
country1 <- raster::shapefile(paste0(root,'/data/shps/',country,'/',iso3c,'_adm',adm_lvl,'.shp'))
country1@data$Regions_GP <- case_when(country1@data$NAME_1 %in% c("Edo", "Delta", "Bayelsa", "Rivers","Akwa Ibom", "Cross River" ) ~ 'SS',
                                      country1@data$NAME_1 %in% c("Oyo", "Ogun", "Lagos", "Osun", "Ondo", "Ekiti") ~ 'SW',
                                      country1@data$NAME_1 %in% c("Anambra", "Enugu", "Ebonyi", "Imo", "Abia") ~ 'SE',
                                      country1@data$NAME_1 %in% c("Kwara", "Niger", "FCT Abuja", "Kogi", "Nassarawa", "Plateau", "Benue") ~ 'NC',
                                      country1@data$NAME_1 %in% c("Bauchi", "Gombe", "Yobe", "Borno", "Adamawa", "Taraba") ~ 'NE',
                                      country1@data$NAME_1 %in% c("Kebbi", "Sokoto", "Zamfara" , "Katsina", "Kano", "Jigawa", "Kaduna") ~ 'NW',
                                      TRUE ~ country1@data$NAME_1) 
shp <- raster::shapefile(paste0(root,'/data/shps/',country,'/',iso3c,'_adm',adm_lvl,'.shp'))
shp <- raster::shapefile(paste0(root,'/data/shps/',country,'/',iso3c,'_adm',adm_lvl,'.shp'))
glue::glue('shp <- shp[shp@data$NAME_{adm_lvl} %in% county,]') %>%
  as.character %>%
  parse(text = .) %>%
  eval(expr = ., envir = .GlobalEnv)
plot(shp)


# Load id coords
crd <- vroom('//dapadfs/workspace_cluster_8/climateriskprofiles/data/id_all_country.csv')
crd <- crd %>%
  dplyr::filter(Country %in% country)
pnt <- crd %>% dplyr::select(x,y) %>% sp::SpatialPoints(coords = .)
crs(pnt) <- crs(shp)
# Filter coordinates that are present in the county
pnt <- sp::over(pnt, shp) %>% data.frame %>% dplyr::select(ISO) %>% complete.cases() %>% which()
crd <- crd[pnt,]
crd <<- crd



# Lectura de datos...

# Observed data for each country... 
tictoc::tic()
all_climate <- tibble(county) %>% 
  mutate(data = purrr::map(.x = county, .f =  function(z){
    data <- fst::fst(glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/data/observational_data/{country}/{z}.fst')) %>% 
      as_tibble() %>% dplyr::select(-id1)  %>% nest(-id, -x, -y, -ISO3) %>% mutate(Country = country) %>% 
      rename(climate = 'data') %>% 
      dplyr::select(id, x, y, ISO3, Country, climate)}))  %>% unnest()
tictoc::toc()


all_climate <-  all_climate %>% dplyr::select(-x, -y) %>% inner_join(m_coords, . , by = c('county', 'id')) 


# Datos historicos. 
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

##  =-----------------------------------------
Clim_graph <- function(historic, value_chain = NULL){
  ISO3 <- unique(historic$ISO3); 
  Country <- unique(historic$Country)
  
  shp_sf <- shp  %>% sf::st_as_sf()
  country <- country1 %>% sf::st_as_sf()
  country_2 <- country %>% group_by(Regions_GP) %>% summarize()
  xlims <- sf::st_bbox(shp_sf)[c(1, 3)]
  ylims <- sf::st_bbox(shp_sf)[c(2, 4)]
  
  #### Voy por aqui. 
  Country <- ifelse(!is.null(value_chain), glue::glue('{Country}_{value_chain}'))
  path <- glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/NIRSAL/{Country}/graphs/maps/')
  if(!dir.exists(path)){dir.create(path, recursive = TRUE)}else{print('ok')}
  # =--------------------------------------------------------
  historic <- historic %>% replace_na(list(z = mean)) %>% 
    dplyr::group_by( id, x, y, ISO3, Country) %>% # dplyr::select_if(is.numeric) %>% 
    dplyr::summarise_all(~round(. , 1))
  # =--------------------------------------------------------
  # Texto...
  af <- as_tibble(st_centroid(shp_sf) %>% st_coordinates()) %>%
    mutate( name = shp_sf$NAME_1) %>%
    mutate(Initals = substr(name, start = 1, stop = 3))
  # geom_text(data = af, aes(X, Y, label = Initals), colour ='black')
  
  af_1 <- as_tibble(st_centroid(country_2) %>% st_coordinates()) %>%
    mutate(Initals = country_2$Regions_GP) %>% 
    filter(Initals != 'Water body')
  
  # =-------------------------------------------------------------------------------
  pais <- ggplot() +
    geom_sf(data = shp_sf, aes(fill = NAME_1), color = gray(.8), alpha = 0.8) +
    geom_sf(data = country_2, fill = NA, color = gray(.1)) +
    theme_bw() +
    # labs(title = country) +
    labs(x = NULL, y = NULL, fill = 'County') +
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black') +
    geom_text(data = af_1, aes(X, Y, label = Initals), colour ='lightgray') +
    # coord_sf(xlim = xlims, ylim = ylims) +
    # scale_fill_brewer('County',palette="Spectral") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_blank(), axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 7, face = "bold"))
  
  # pais + labs(subtitle = glue::glue('{NULL}'))
  
  ggsave(glue::glue('{path}/Country.png') , width = 8, height = 5.5)
  
  
  
  # =--------------------------------------------------------
  historic <- historic %>% dplyr::select(-x, -y, -county ) %>% inner_join(m_coords)
  # =--------------------------------------------------------
  
  graph <-  historic[historic$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = prec ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = historic[historic$county%in% county[i],] , aes(x = x, y = y, fill = prec )) 
  }
  
  prec <- graph   +
    geom_sf(data = country_2, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.5)) +
    coord_sf(xlim = xlims, ylim = ylims) +
    geom_text(data = af_1, aes(X, Y, label = Initals), colour ='black') +
    labs(fill = glue::glue('(mm)  '), title = 'Historical Annual Mean Precipitation (mm/year)',x = 'Longitude', y = 'Latitude') +
    scale_fill_gradientn(colours = blues9, 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() + theme(legend.position = 'bottom')
  
  
  ggsave(glue::glue('{path}/H_prec.png') , width = 8, height = 5.5, dpi = 300)
  
  ##########################################################
  
  
  graph <-  historic[historic$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = tmean ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = historic[historic$county%in% county[i],] , aes(x = x, y = y, fill = tmean )) 
  }
  
  
  tmn <- graph +
    geom_sf(data = country_2, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.5)) +
    coord_sf(xlim = xlims, ylim = ylims) +
    geom_text(data = af_1, aes(X, Y, label = Initals), colour ='black') +
    labs(fill = expression('('*~degree*C*')'), title = expression('Historical Annual Mean Temperature ('*~degree*C*')'),x = 'Longitude', y = 'Latitude') +
    scale_fill_gradient(low = "yellow", high = "red",
                        guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0)))+
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() + theme(legend.position = 'bottom')
  
  ggsave(glue::glue('{path}/H_tmn.png') , width = 8, height = 5.5, dpi = 300)
  
  
  png(filename = glue::glue('{path}/A_Multi_Anual.png'), width = 1580, height = 720)
  print(gridExtra::grid.arrange(prec, tmn, ncol=2,
                                top = glue::glue('{Country}',
                                                 bottom =   "Data source: Alliance Bioversity-CIAT")))
  dev.off()
  
}


# =---------------------------------
# =---------------------------------
# graph de clima. 
Clim_graph(historic, value_chain= value_chain)


# =-------------------------------------

do_clim_country <- function(data_split){
  path <- glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/NIRSAL/{country}_{value_chain}')
  ISO3 <- unique(data_split$ISO3);# county <- unique(data_split$county)
  # time <- unique(data_split$time);
  semester <-  unique(data_split$semester) ;  Country <- unique(data_split$Country)
  
  #### Voy por aqui. 
  Country <- ifelse(!is.null(value_chain), glue::glue('{country}_{value_chain}'))
  
  if(dir.exists(glue::glue('{path}/graphs/maps/'))==FALSE){dir.create(glue::glue('{path}/graphs/maps/'))}else{print('ok')}
  
  # Conditions...
  if(Big =='N'){
    
    median_data <- data_split %>%
      dplyr::group_by(id, county, Country, x, y, ISO3, semester, time) %>%
      dplyr::select(-year) %>%
      dplyr::summarise_all(mean) %>%  dplyr::ungroup() %>%
      dplyr::mutate(NT35 = round(x = NT35, digits = 0), 
                    ndws = round(ndws,0))
    
  }else if(Big == 'B'){
    # =---------------------------------------------------------
    median_data <- data_split 
    
  } else{print('Change big argument... >.<')}
  
  limits <- dplyr::select(median_data, CDD, P5D, P95, NT35, ndws) %>%
    dplyr::summarise_all(.funs = c('min', 'max'))
  
  median_data <- median_data %>% dplyr::filter(time == 'future') %>%
    dplyr::rename('CDD_f' = 'CDD'   , 'P5D_f' = 'P5D'  , 'P95_f'= 'P95' , 'NT35_f' = 'NT35', 
                  'ndws_f' = 'ndws' ) %>%
    dplyr::select(-time) %>%
    dplyr::inner_join(dplyr::filter(median_data , time == 'past') %>% dplyr::select(-time), . ) %>%
    dplyr::mutate(CDD_c = CDD_f - CDD, P5D_c = P5D_f - P5D, P95_c = P95_f - P95, NT35_c = NT35_f -NT35, ndws_c = ndws_f - ndws)
  
  
  # Aqui se hace solo la figura base...
  shp_sf <- shp  %>% sf::st_as_sf()
  country <- country1 %>% sf::st_as_sf() %>% group_by(Regions_GP) %>% summarize()
  xlims <- sf::st_bbox(shp_sf)[c(1, 3)]
  ylims <- sf::st_bbox(shp_sf)[c(2, 4)]
  
  
  #===---------------------------------------------------------
  #=------------------------------------------------------------
  af <- as_tibble(st_centroid(country) %>% st_coordinates()) %>%
    mutate(Initals = country$Regions_GP) %>% 
    filter(Initals != 'Water body')
  
  #===---------------------------------------------------------
  #=------------------------------------------------------------
  
  # Primero dejaré hechos los de presente... luego repito los de futuro...
  # Esta función va a quedar super manual.
  index_a <- 'CDD'
  
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = CDD ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = CDD )) 
  }
  a <- graph +
    geom_sf(data = country, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) + #  aes(colour = NAME_1)     geom_text(data = af, aes(X, Y, label = Initals), colour ='black')+  #  aes(colour = NAME_1)
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black')+ 
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('{index_a}\n(days)'), title = 'Historic', x = 'Longitude', y = 'Latitude') +
    scale_fill_viridis_c(limits = c(round(limits$CDD_min,2)-0.1, round(limits$CDD_max,2)+0.1), 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  
  ggsave(glue::glue('{path}/graphs/maps/{index_a}_past_S{semester}.png') , width = 8, height = 5.5)
  
  
  # =- Lo mismo para futuro...
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = CDD_f ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = CDD_f )) 
  }
  a1 <- graph +
    geom_sf(data = country, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) + #  aes(colour = NAME_1)     
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black')+ 
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('{index_a}\n(days)'), title = 'Future',x = 'Longitude', y = 'Latitude') +
    scale_fill_viridis_c(limits = c(round(limits$CDD_min,2)-0.1, round(limits$CDD_max,2)+0.1), 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  ggsave(glue::glue('{path}/graphs/maps/{index_a}_future_S{semester}.png') , width = 8, height = 5.5)
  
  
  # =- 
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = CDD_c ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = CDD_c )) 
  }
  a_d <- graph+
    geom_sf(data = country, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) + #  aes(colour = NAME_1)     
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black')+ 
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('{index_a}\n(days) '), title = 'Change', x = 'Longitude', y = 'Latitude') +
    scale_fill_gradient2(low = '#000099', mid = 'white', high = '#A50026', 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  ggsave(glue::glue('{path}/graphs/maps/Dif_{index_a}_S{semester}.png') , width = 8, height = 5.5)
  
  
  
  png(filename = glue::glue('{path}/graphs/maps/Dif_{index_a}_{semester}.png') , width = 1580, height = 720)
  print(gridExtra::grid.arrange(a, a1, a_d, ncol=3,  
                                top = glue::glue('{Country}\nS:{semester}',
                                                 bottom =   "Data source: Alliance Bioversity-CIAT")))
  dev.off()
  
  
  #===---------------------------------------------------------
  #=------------------------------------------------------------
  
  #   # =---------------------------------------
  #   # Siguiente indice...
  index_c <- 'P5D'
  
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = P5D ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = P5D )) 
  }
  c <- graph +
    geom_sf(data = country, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) + #  aes(colour = NAME_1)     
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black')+ 
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('{index_c}\n(mm) '), title = 'Historic',x = 'Longitude', y = 'Latitude') +
    scale_fill_gradientn(colours = blues9, limits = c(round(limits$P5D_min, 2)- 0.1, 
                                                      round(limits$P5D_max, 2)+0.1), 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  ggsave(glue::glue('{path}/graphs/maps//{index_c}_past_S{semester}.png') , width = 8, height = 5.5)
  
  # =- Futuro.
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = P5D_f ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = P5D_f )) 
  }
  c1 <- graph +
    geom_sf(data = country, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) + #  aes(colour = NAME_1)     
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black')+ 
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('{index_c}\n(mm)'), title = 'Future', x = 'Longitude', y = 'Latitude') +
    scale_fill_gradientn(colours = blues9, limits = c(round(limits$P5D_min, 2)- 0.1, 
                                                      round(limits$P5D_max, 2)+0.1), 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  ggsave(glue::glue('{path}/graphs/maps//{index_c}_future_S{semester}.png') , width = 8, height = 5.5)
  
  # =----
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = P5D_c ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = P5D_c )) 
  }
  c_d <- graph +
    geom_sf(data = country, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) + #  aes(colour = NAME_1)     
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black')+ 
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('{index_c}\n(mm) '), title = 'Change', x = 'Longitude', y = 'Latitude') +
    scale_fill_gradient2(low = '#A50026', mid = 'white', high = '#000099', 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  ggsave(glue::glue('{path}/graphs/maps//Dif_{index_c}_S{semester}.png') , width = 8, height = 5.5)
  
  
  png(filename = glue::glue('{path}/graphs/maps//Dif_{index_c}_{semester}.png') , width = 1580, height = 720)
  print(gridExtra::grid.arrange(c, c1, c_d, ncol=3,  
                                top = glue::glue('{Country}\nS:{semester}',
                                                 bottom =   "Data source: Alliance Bioversity-CIAT")))
  dev.off()
  
  
  
  
  #===---------------------------------------------------------
  # =------------
  index_d <- 'P95'
  
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = P95 ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = P95 )) 
  }
  
  d <- graph +
    geom_sf(data = country, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) + #  aes(colour = NAME_1)     
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black')+ 
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('{index_d}\n(mm)'), title = 'Historic', x = 'Longitude', y = 'Latitude') +
    scale_fill_gradientn(colours = blues9, limits = c(round(limits$P95_min,2)-0.1, 
                                                      round(limits$P95_max, 2)+0.1), 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  ggsave(glue::glue('{path}/graphs/maps//{index_d}_past_S{semester}.png') , width = 8, height = 5.5)
  
  # =----
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = P95_f ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = P95_f )) 
  }
  d1 <- graph + geom_sf(data = country, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) + #  aes(colour = NAME_1)     
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black')+ 
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('{index_d}\n(mm) '), title = 'Future', x = 'Longitude', y = 'Latitude') +
    scale_fill_gradientn(colours = blues9, limits = c(round(limits$P95_min,2)-0.1, 
                                                      round(limits$P95_max, 2)+0.1), 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  
  ggsave(glue::glue('{path}/graphs/maps//{index_d}_future_S{semester}.png') , width = 8, height = 5.5)
  
  
  # =----
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = P95_c ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = P95_c )) 
  }
  d_d <- graph + geom_sf(data = country, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) + #  aes(colour = NAME_1)     
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black')+ 
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('{index_d}\n(mm) '), title = 'Change', x = 'Longitude', y = 'Latitude') +
    scale_fill_gradient2(low = '#A50026', mid = 'white', high = '#000099', 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  ggsave(glue::glue('{path}/graphs/maps//Dif_{index_d}_S{semester}.png') , width = 8, height = 5.5)
  
  
  
  png(filename = glue::glue('{path}/graphs/maps//Dif_{index_d}_{semester}.png') , width = 1580, height = 720)
  print(gridExtra::grid.arrange(d, d1, d_d, ncol=3,  
                                top = glue::glue('{Country}\nS:{semester}',
                                                 bottom =   "Data source: Alliance Bioversity-CIAT")))
  dev.off()
  
  
  
  #===---------------------------------------------------------
  #=------------------------------------------------------------
  
  # =------------
  index_e <- 'NT35'
  
  # =----
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = NT35 ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = NT35 )) 
  }
  e <- graph + geom_sf(data = country, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) + #  aes(colour = NAME_1)     
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black')+ 
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('{index_e}\n(days)'), title = 'Historic', x = 'Longitude', y = 'Latitude') +
    scale_fill_viridis_c(limits = c(round(limits$NT35_min, 2) - 0.1, round(limits$NT35_max, 2)+0.1), 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  ggsave(glue::glue('{path}/graphs/maps//{index_e}_past_S{semester}.png') , width = 8, height = 5.5)
  
  # =------------
  # =----
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = NT35_f ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = NT35_f )) 
  }
  e1 <- graph + geom_sf(data = country, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) + #  aes(colour = NAME_1)     
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black')+ 
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('{index_e}\n(days)'), title = 'Future', x = 'Longitude', y = 'Latitude') +
    scale_fill_viridis_c(limits = c(round(limits$NT35_min, 2) - 0.1, round(limits$NT35_max, 2)+0.1), 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  ggsave(glue::glue('{path}/graphs/maps//{index_e}_future_S{semester}.png') , width = 8, height = 5.5)
  
  
  # =----
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = NT35_c ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = NT35_c )) 
  }
  e_d <- graph + geom_sf(data = country, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) + #  aes(colour = NAME_1)     
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black')+ 
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('{index_e}\n(days)'), title = 'Change', x = 'Longitude', y = 'Latitude') +
    scale_fill_gradient2(low = '#000099', mid = 'white', high = '#A50026', 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  ggsave(glue::glue('{path}/graphs/maps//Dif_{index_e}_S{semester}.png') , width = 8, height = 5.5)
  
  
  png(filename = glue::glue('{path}/graphs/maps//Dif_{index_e}_{semester}.png') , width = 1580, height = 720)
  print(gridExtra::grid.arrange(e, e1, e_d, ncol=3,  
                                top = glue::glue('{Country}\nS:{semester}',
                                                 bottom =   "Data source: Alliance Bioversity-CIAT")))
  dev.off()
  
  
  
  
  
  #===---------------------------------------------------------
  #=------------------------------------------------------------
  cowsay::say('here')
  # =------------
  index_f <- 'ndws'
  
  # =----
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = ndws ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = ndws )) 
  }
  
  f <- graph + geom_sf(data = country, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black')+ 
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('{index_f}\n(days)  '), title = 'Historic', x = 'Longitude', y = 'Latitude') +
    scale_fill_viridis_c(limits = c(round(limits$ndws_min, 2), round(limits$ndws_max, 2)), 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() + theme(legend.position = 'bottom')
  
  ggsave(glue::glue('{path}/graphs/maps//{index_f}_past_S{semester}.png') , width = 8, height = 5.5)
  
  # =------------
  # =----
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = ndws_f ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = ndws_f )) 
  }
  
  f1 <- graph + geom_sf(data = country, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black')+ 
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('{index_f}\n(days)  '), title = 'Future', x = 'Longitude', y = 'Latitude') +
    scale_fill_viridis_c(limits = c(round(limits$ndws_min, 2), round(limits$ndws_max, 2)), 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  ggsave(glue::glue('{path}/graphs/maps//{index_e}_future_S{semester}.png') , width = 8, height = 5.5)
  
  
  # =----
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = ndws_c ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = ndws_c )) 
  }
  f_d <- graph + geom_sf(data = country, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black')+ 
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('{index_f}\n(days)  '), title = 'Change', x = 'Longitude', y = 'Latitude') +
    scale_fill_gradient2(low = '#000099', mid = 'white', high = '#A50026', 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  ggsave(glue::glue('{path}/graphs/maps//Dif_{index_f}_S{semester}.png') , width = 8, height = 5.5)
  
  
  png(filename = glue::glue('{path}/graphs/maps//Dif_{index_f}_{semester}.png') , width=12.5,height=4.5,units="in", res = 300)
  print(gridExtra::grid.arrange(f, f1, f_d, ncol=3,  
                                top = glue::glue('{Country}\nS:{semester}',
                                                 bottom =   "Data source: Alliance Bioversity-CIAT")))
  dev.off()
  
  
  return(median_data)}



# =------------
# Funcion for run basic maps. 

do_srad_country <- function(data_split){
  path <- glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/NIRSAL/{country}_{value_chain}')
  ISO3 <- unique(data_split$ISO3); #county <- unique(data_split$county) 
  Country <- ifelse(!is.null(value_chain), glue::glue('{country}_{value_chain}'))
  Big <- unique(data_split$Big)
  
  # Aqui se hace solo la figura base...
  shp_sf <- shp  %>% sf::st_as_sf()
  pais <- country1 %>% sf::st_as_sf() %>% group_by(Regions_GP) %>% summarize()
  xlims <- sf::st_bbox(shp_sf)[c(1, 3)]
  ylims <- sf::st_bbox(shp_sf)[c(2, 4)]
  
  #### Voy por aqui. 
  Country <- ifelse(!is.null(value_chain), glue::glue('{Country}_{value_chain}'))
  
  # Texto... 
  af <- as_tibble(st_centroid(pais) %>% st_coordinates()) %>%
    mutate(Initals = pais$Regions_GP) %>% 
    filter(Initals != 'Water body')
  
  if(Big =='N'){
    # =----------------------------------------------------------------------
    # gSeason...
    gSeason_i <- data_split %>% dplyr::select(time, x, y,gSeason) %>%
      group_by( time, x, y) %>% summarise(gSeason = max(gSeason)) %>%
      ungroup() %>% group_by(time, x, y) %>%
      summarise(gSeason = round(mean(gSeason, na.rm = TRUE), 0)) %>% ungroup()
    
    # =----------------------------------------------------------------------
    # SLGP -- LGP 
    two_index <- data_split %>% dplyr::select( time, x, y, gSeason, SLGP, LGP) %>% 
      group_by(time, x, y, gSeason) %>% 
      summarise_all(.f = function(x){round(mean(x, na.rm = TRUE), 0)}) %>% 
      ungroup(year) %>%  arrange(gSeason) %>% 
      nest(data = c('x', 'y','SLGP', 'LGP')) %>% 
      drop_na() %>% filter(gSeason < 3) %>% unnest() %>% ungroup()
    
  }else if(Big == 'B'){
    # =----------------------------------------------------------------------
    # gSeason...
    gSeason_i <- data_split %>% dplyr::select( time, x, y, gSeason) %>%
      group_by( time, x, y) %>% summarise(gSeason = max(gSeason, na.rm =  TRUE)) %>%
      ungroup() %>% group_by(time, x, y) %>%
      summarise(gSeason = round(mean(gSeason, na.rm = TRUE), 0)) %>% ungroup()
    
    # =----------------------------------------------------------------------
    # gSeason...
    two_index <- data_split %>% dplyr::select( time, x, y, gSeason, SLGP, LGP) %>% 
      group_by(time, x, y, gSeason) %>% 
      summarise_all(.f = function(x){round(mean(x, na.rm = TRUE), 0)}) %>% 
      # ungroup(year) %>%  
      arrange(gSeason) %>% 
      nest(data = c('x', 'y','SLGP', 'LGP')) %>% 
      drop_na() %>% filter(gSeason < 3) %>% unnest() %>% ungroup()
  } else{print('Change big argument... >.<')}
  
  # =---------------------------------------------------------------------
  # =---------------------------------------------------------------------
  # gSeason
  limits_gs <- dplyr::select(gSeason_i, gSeason) %>% summarise_all(.funs = c('min', 'max'))
  
  # # Primero dejaré hechos los de presente... luego repito los de futuro...
  
  median_data <- inner_join(filter(gSeason_i, time == 'past'), m_coords)
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = gSeason ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = gSeason )) 
  }
  
  
  gs <- graph + geom_sf(data = pais, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black') +
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('gSeason\n(day)  '), title = 'Historic', x = 'Longitude', y = 'Latitude') +
    scale_fill_gradientn(colours = blues9, limits = as.integer(limits_gs), 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  
  ggsave(glue::glue('{path}/graphs/maps//gSeason_past.png') , width = 8, height = 5.5)
  
  # # =- Futuro.
  median_data <- inner_join(filter(gSeason_i, time == 'future'), m_coords)
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = gSeason ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = gSeason )) 
  }
  
  gs_f <- graph + geom_sf(data = pais, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black') +
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('gSeason\n(day)  '), title = 'Future', x = 'Longitude', y = 'Latitude') +
    scale_fill_gradientn(colours = blues9, limits = as.integer(limits_gs), 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  ggsave(glue::glue('{path}/graphs/maps//gSeason_future.png') , width = 8, height = 5.5)
  
  # =--------
  gSeason_dif <- gSeason_i %>%
    pivot_wider(names_from = time, values_from = gSeason) %>%
    mutate(gSeason = future - past)
  
  median_data <- inner_join(gSeason_dif, m_coords)
  graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = gSeason ))
  for(i in 2:length(county) ){
    graph <- graph + geom_tile(data = median_data[median_data$county%in% county[i],] , aes(x = x, y = y, fill = gSeason )) 
  }
  
  c <- graph + geom_sf(data = pais, fill = NA, color = gray(.1)) +
    # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
    geom_text(data = af, aes(X, Y, label = Initals), colour ='black') +
    coord_sf(xlim = xlims, ylim = ylims) +
    labs(fill = glue::glue('gSeason\n(days)'), title = 'Change', x = 'Longitude', y = 'Latitude') +
    scale_fill_gradient2(low = '#A50026', mid = 'white', high = '#000099', 
                         guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
    scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
    scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
    theme_bw() +
    theme(legend.position = 'bottom')
  
  ggsave(glue::glue('{path}/graphs/maps//Dif_gSeason.png') , width = 8, height = 5.5)
  
  
  png(filename = glue::glue('{path}/graphs/maps//all_gSeason.png') , width = 1580, height = 720)
  print(gridExtra::grid.arrange(gs, gs_f, c, ncol=3,  
                                top = glue::glue('{Country}',
                                                 bottom =   "Data source: Alliance Bioversity-CIAT")))
  dev.off()
  
  
  # =---------------------------------------------------------------------
  # =---------------------------------------------------------------------
  # Aqui hay que hacer los filtros de acuerdo al grupo ...
  
  # rows_n <- data_split %>% filter(time == 'past', year == '1985') %>% dplyr::select(x, y, id) %>% unique() %>% nrow()
  
  SLGP_dif <- two_index %>% dplyr::select(-LGP)  %>%
    pivot_wider(names_from = time, values_from = SLGP) %>%
    mutate(SLGP = future - past, gSeason = glue::glue('gSeason = {gSeason}'))
  
  limits_two <-   two_index %>% group_by(gSeason) %>%
    dplyr::select(SLGP, LGP) %>% summarise_all(.funs = c('min', 'max')) %>% 
    ungroup()
  
  gS <- unique(two_index$gSeason)
  
  for(i in gS){
    
    median_data <- inner_join(filter(two_index %>% mutate(gSeason = glue::glue('gSeason = {gSeason}')),
                                     time == 'past', gSeason == glue::glue('gSeason = {i}')), m_coords)
    graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = SLGP))
    for(k in 2:length(county) ){
      graph <- graph + geom_tile(data = median_data[median_data$county%in% county[k],] , aes(x = x, y = y, fill = SLGP)) 
    }
    
    SLGP_p_1 <-graph + geom_sf(data = pais, fill = NA, color = gray(.1)) +
      # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      geom_text(data = af, aes(X, Y, label = Initals), colour ='black') +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('SLGP\n(Day of\nthe year)  '), 
           title = glue::glue('gSeason = {i}; Historic'),
           x = 'Longitude', y = 'Latitude') +
      scale_fill_gradientn(colours = blues9, limits = c(limits_two$SLGP_min[i], limits_two$SLGP_max[i]), 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}/graphs/maps//SLGP_past_{i}.png') , width = 8, height = 4)
    
    # =----
    median_data <- inner_join(filter(two_index %>% mutate(gSeason = glue::glue('gSeason = {gSeason}')),
                                     time == 'future', gSeason == glue::glue('gSeason = {i}')), m_coords)
    graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = SLGP))
    for(k in 2:length(county) ){
      graph <- graph + geom_tile(data = median_data[median_data$county%in% county[k],] , aes(x = x, y = y, fill = SLGP)) 
    }
    
    SLGP_f_1 <-graph + geom_sf(data = pais, fill = NA, color = gray(.1)) +
      # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      geom_text(data = af, aes(X, Y, label = Initals), colour ='black') +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('SLGP\n(Day of\nthe year)  '), 
           title = glue::glue('gSeason = {i}; Future'), x = 'Longitude', y = 'Latitude') +
      scale_fill_gradientn(colours = blues9, limits = c(limits_two$SLGP_min[i], limits_two$SLGP_max[i]), 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    
    ggsave(glue::glue('{path}/graphs/maps//SLGP_future_{i}.png') , width = 8, height = 4)
    
    # =----
    median_data <- inner_join(filter(SLGP_dif, gSeason == glue::glue('gSeason = {i}')), m_coords)
    graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = SLGP))
    for(k in 2:length(county) ){
      graph <- graph + geom_tile(data = median_data[median_data$county%in% county[k],] , aes(x = x, y = y, fill = SLGP)) 
    }
    
    d <- graph + geom_sf(data = pais, fill = NA, color = gray(.1)) +
      # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      geom_text(data = af, aes(X, Y, label = Initals), colour ='black') +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('SLGP\n(days)  '), 
           title = glue::glue('gSeason = {i}; Change'),x = 'Longitude', y = 'Latitude') +
      scale_fill_gradient2(low = '#A50026', mid = 'white', high = '#000099', 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}/graphs/maps//Dif_SLGP_{i}.png') , width = 8, height = 4)
    
    
    png(filename = glue::glue('{path}/graphs/maps//all_SLGP_{i}.png') , width = 1580, height = 720)
    print(gridExtra::grid.arrange(SLGP_p_1, SLGP_f_1, d, ncol=3,  
                                  top = glue::glue('{Country}',
                                                   bottom =   "Data source: Alliance Bioversity-CIAT")))
    dev.off()
  }
  
  
  #### =---------------------------------------------------------------------------
  
  # =--------------------------------------------------------------------
  # =--------------------------------------------------------------------
  LGP_dif <- two_index %>% dplyr::select(-SLGP)  %>%
    pivot_wider(names_from = time, values_from = LGP) %>%
    mutate(LGP = future - past, gSeason = glue::glue('gSeason = {gSeason}'))
  
  for(i in gS){
    median_data <- inner_join(filter(two_index %>% mutate(gSeason = glue::glue('gSeason = {gSeason}')),
                                     time == 'past', gSeason == glue::glue('gSeason = {i}')), m_coords)
    graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = LGP))
    for(k in 2:length(county) ){
      graph <- graph + geom_tile(data = median_data[median_data$county%in% county[k],] , aes(x = x, y = y, fill = LGP)) 
    }
    
    LGP_p <- graph + geom_sf(data = pais, fill = NA, color = gray(.1)) +
      # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      geom_text(data = af, aes(X, Y, label = Initals), colour ='black') +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('LGP\n(days)  '), 
           title = glue::glue('gSeason = {i}; Historic'),x = 'Longitude', y = 'Latitude') +
      scale_fill_gradientn(colours = blues9, limits = c(limits_two$LGP_min[i], limits_two$LGP_max[i]), 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}/graphs/maps//LGP_past_{i}.png') , width = 8, height = 4)
    
    # =------
    median_data <- inner_join(filter(two_index %>% mutate(gSeason = glue::glue('gSeason = {gSeason}')),
                                     time == 'future', gSeason == glue::glue('gSeason = {i}')), m_coords)
    graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = LGP))
    for(k in 2:length(county) ){
      graph <- graph + geom_tile(data = median_data[median_data$county%in% county[k],] , aes(x = x, y = y, fill = LGP)) 
    }
    
    LGP_f <- graph + geom_sf(data = pais, fill = NA, color = gray(.1)) +
      # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      geom_text(data = af, aes(X, Y, label = Initals), colour ='black') +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('LGP\n(days)  '), 
           title = glue::glue('gSeason = {i}; Future'), x = 'Longitude', y = 'Latitude') +
      scale_fill_gradientn(colours = blues9, limits = c(limits_two$LGP_min[i], limits_two$LGP_max[i]), 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}/graphs/maps//LGP_future_{i}.png') , width = 8, height = 4)
    
    # =------
    median_data <- inner_join(filter(LGP_dif, gSeason == glue::glue('gSeason = {i}')) , m_coords)
    graph <-  median_data[median_data$county%in% county[1],]  %>% ggplot(.) + geom_tile(aes(x = x, y = y, fill = LGP))
    for(k in 2:length(county) ){
      graph <- graph + geom_tile(data = median_data[median_data$county%in% county[k],] , aes(x = x, y = y, fill = LGP)) 
    }
    
    e <- graph + geom_sf(data = pais, fill = NA, color = gray(.1)) +
      # geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      geom_text(data = af, aes(X, Y, label = Initals), colour ='black') +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('LGP\n(days)  '), 
           title = glue::glue('gSeason = {i}; Change'),x = 'Longitude', y = 'Latitude') +
      scale_fill_gradient2(low = '#A50026', mid = 'white', high = '#000099', 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}/graphs/maps//Dif_LGP_{i}.png') , width = 8, height = 4)
    
    png(filename = glue::glue('{path}/graphs/maps//all_LGP_{i}.png') , width = 1580, height = 720)
    print(gridExtra::grid.arrange(LGP_p, LGP_f, e, ncol=3,  
                                  top = glue::glue('{Country}',
                                                   bottom =   "Data source: Alliance Bioversity-CIAT")))
    dev.off()
    
  }
  
  
}







# =---------------------------------------------------------------------------------------
# =--------------------------------------
# =-------------------------------------

# Optimization. 
ag <- all_climate %>% dplyr::select(-climate) %>% 
  nest(-county) %>% 
  group_split(county)

# county <- ag$county

tag <- ag %>% 
  purrr::map(.f = function(z){
    
    # glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/NIRSAL/{Country}_{value_chain}')
    path <- '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/NIRSAL/'
    # Si... se tienen los idw. 
    past_c <- fst::fst(glue::glue('{path}/{country}_{value_chain}/past/{z$county}_1985_2015_idw.fst')) %>%
      as_tibble() %>% 
      dplyr::mutate(time = 'past') 
    
    futDir  <- paste0(path,country, '_', value_chain,'/future')
    fut_fls <- list.files(futDir, pattern = paste0('^',z$county,'_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_idw.fst'), recursive = T)
    fut_fls <- paste0(futDir,'/',fut_fls)  
    
    future_c  <- fut_fls %>%
      purrr::map(.f = function(x){df <- fst(x) %>% as_tibble() %>% dplyr::mutate(time = 'future'); return(df)}) %>%
      dplyr::bind_rows()
    
    
    data_c <- bind_rows(past_c, future_c) %>% dplyr::select(-x, -y)
    return(data_c)  }) %>% 
  purrr::map(.f = nest) 

data_for_graphs <- list()
for(i in 1:length(ag)){
  data_for_graphs[[i]] <- bind_cols(ag[[i]] , tag[[i]]) %>% setNames(c('county', 'data', 'data1'))
}

data_for_graphs <- bind_rows(data_for_graphs) %>% rename(data_graph = 'data1') %>% unnest(data) %>% 
  nest(id, x, y, ISO3, Country) 


# =----- 
library(future)
library(furrr)

# Puedo leer esto en paralelo, para que no haya problema... 
cores <- 5
plan(cluster, workers = cores)

probando <- data_for_graphs %>% 
  mutate(data_graph = furrr::future_map2(.x = data_graph, .y = data, .f = function(z,r){
    data_n <- inner_join(z, r)  %>% 
      separate(season, c('ok','semester'), sep = 's')})) %>% 
  dplyr::select(-data) %>% 
  unnest %>% dplyr::select(-county1, -ok)

data_all <- probando %>%  
  dplyr::select(id, ISO3, county, Country, x, y, time ,semester, CDD, P5D, P95, NT35, ndws) %>% 
  unique() %>% group_split(semester)

# data_split <- data_all[[1]]

# =--------------------------------------

data_all %>% purrr::walk(.f = do_clim_country)
# =--------------------------------------
# =-------------------------------------


# Srad index....
index_complete <- probando %>% 
  dplyr::select(id, ISO3, county, Country, x, y, time, gSeason, SLGP,LGP) %>% 
  mutate(Big = 'B')


# =-----------------------

do_srad_country(data_split = index_complete)
