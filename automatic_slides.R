# Create slides automatically
# A. Esquivel and H. Achicanoy
# Alliance Bioversity-CIAT, 2020

options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, rmarkdown))

sld_dir <- 'D:/Slides'

countiesList <- c('Bungoma','Kiambu','Kirinyaga','Kisii','Kitui','Migori','Muranga','Nandi','Narok','Nyamira','Samburu','Trans-Nzoia','Turkana','Vihiga')

create_slides <- function(country = 'Kenya', iso3 = 'KEN', countiesList = countiesList){
  
  for(county in countiesList){
    ppt <- readLines('D:/ppt_template_kcp.Rmd')
    ppt <- ppt %>%
      purrr::map(.f = function(l){
        l <- gsub(pattern = 'COUNTY', replacement = county, x = l, fixed = T)
        l <- gsub(pattern = 'COUNTRY', replacement = country, x = l, fixed = T)
        return(l)
      }) %>%
      unlist()
    writeLines(ppt, paste0(sld_dir,'/',country,'/',county,'-',iso3,'.Rmd'))
    rmarkdown::render(paste0(sld_dir,'/',country,'/',county,'-',iso3,'.Rmd'))
  }
  
  return('Process done!\n')
  
}
create_slides(country = 'Kenya', iso3 = 'KEN', countiesList = countiesList)
