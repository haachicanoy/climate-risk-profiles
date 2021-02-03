# -------------------------------------------------- #
# Climate Risk Profiles -- Master code
# A. Esquivel, H. Achicanoy, and J. Ramirez-Villegas
# Alliance Bioversity-CIAT, 2021
# -------------------------------------------------- #

# R options and load packages
options(warn = -1, scipen = 999)
suppressMessages(if(!require(pacman)){install.packages('pacman'); library(pacman)})
suppressMessages(pacman::p_load(tidyverse, raster,
                                terra, sp, compiler))

root <- 'D:/OneDrive - CGIAR/PARM toolkit/notebook'

# Establish working directory
setwd(root)

# Import functions
if(!dir.exists('./scripts')){
  # Download updated repository
  download.file(url = 'https://github.com/haachicanoy/climate-risk-profiles/archive/master.zip', destfile = "./crp.zip")
  # Unzip the .zip file
  unzip(zipfile = "crp.zip")
  dir.create('./scripts', recursive = T)
  file.copy2 <- Vectorize(FUN = file.copy, vectorize.args = c('from','to'))
  file.copy2(from = list.files(path = './climate-risk-profiles-master/', pattern = '*.R$', full.names = T),
             to   = paste0('./scripts/',list.files(path = './climate-risk-profiles-master/', pattern = '*.R$', full.names = F)))
  unlink('./climate-risk-profiles-master', recursive = T)
  file.remove('./crp.zip')
  source('./scripts/00_functions.R')
  source('./scripts/_get_soil_data.R')
} else {
  source('./scripts/00_functions.R')
  source('./scripts/_get_soil_data.R')
}

# Object created from shiny app
inputs <- list(country   = 'Burkina Faso',
               county    = 'Sud-Ouest',
               iso3c     = 'BFA',
               adm_lvl   = 1,
               seasons   = 'All year',
               m_seasons = NULL,
               n_wtts    = NULL,
               big_cnt   = T,
               ncores    = 2)

## Step 1. Setting up the study region
# Load country shapefile
shp <- raster::shapefile(paste0(root,'/data/shps/',inputs$iso3c,'.shp'))
# Load a reference raster to obtain geographical coordinates
ref <- raster::raster(paste0(root,'/data/rasters/tmplt.tif'))
# Crop the reference raster and fit it to the shapefile region
ref <- ref %>% raster::crop(x = ., y = raster::extent(shp)) %>% raster::mask(x = ., mask = shp)
# Get the geographical coordinates and their id
crd <- ref %>% raster::rasterToPoints() %>% base::as.data.frame() %>% dplyr::select(x, y)
crd$id <- raster::cellFromXY(object = ref, xy = crd[,c('x','y')])
crd <- crd %>% dplyr::select(id, x, y)

## Step 2. Get soil data
# This function requires to be connected to the dapadfs storage cluster
get_soil(crd = crd, root_depth = 60, outfile = './soilcp_data.fst')

## Step 3. Get daily climate data. FIX
# This function requires to be connected to the dapadfs storage cluster
get_observational_data(country = inputs$country, county  = inputs$county, iso3 = inputs$iso3c, adm_lvl = inputs$iso3c)

## Step 4. Calculate indices. FIX
calc_indices(country = inputs$country, county = inputs$county, iso3c = inputs$iso3c, adm_lvl = inputs$adm_lvl,
             seasons = NULL, # Seasons manually defined
             n_ssns  = 2,    # 2-seasons automatically defined
             n_wtts  = 100,  # 100-wettest days
             big_cnt = TRUE,
             ncores  = 10)

## Step 5. Do graphs
## Step 5.1. Climatology
do_climatology(country = 'Kenya',
               county  = 'Vihiga',
               seasons = TRUE, # Climatology without any season (manual or automatic)
               manual  = NULL, # Seasons defined manually e.g. list(s1 = c(11:12,1:4)
               auto    = list(n = 2))

## Step 5.2. Time series plots
time_series_plot(country = inputs$country, county = inputs$county)

## Step 5.3. Maps
do_maps(country = inputs$country, county = inputs$county)

## Step 5.4. Elevation map
do_alt_map(country = inputs$country, county = inputs$county)

## Step 5.6. Country maps
do_country_maps(country = inputs$country)
