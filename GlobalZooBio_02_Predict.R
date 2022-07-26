
source("utils.R") # Harmonic function

library(raster)
library(ggplot2)
library(ncdf4)
library(dplyr)
library(tidyverse)
library(splines)
library(lmerTest)

# IMPORT COPEPOD DATABASE BIOMASS GLM FROM HENEGHAN ET AL. 2020
m7 <- readRDS(file.path("Output","m7.rds"))

## IMPORT BATHYMETRY DATA AND ORIENT LATITUDES TO BE SOUTH TO NORTH
bathy_data <- readRDS(file.path("one","Bathy_raster_oneDeg.rds"))
bathy_matrix <- t(as.matrix(bathy_data$Bathy))
bathy_matrix <- bathy_matrix[,180:1]

## Calculate areas of grid cells
glob_area <- t(as.matrix(raster::area(raster())))
vars <- c("SST", "Chl")

## Set up array to make predictions and save outputs
save_array <- array(NA, dim = c(12,13,64800)) # Lon, lat, month, variable
dimnames(save_array)[[1]] <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
dimnames(save_array)[[2]] <- c("Lon", "Lat", "BiomassMethod", "Mesh", "Depth", "Gear", "Institution", "HarmTOD", "Bathy", "HarmDOY", "SST", "Chl", "GLM_Mesozoo")

## Set mesh, start depth and time of day
save_array[,"Mesh",] <- 100
save_array[,"Depth",] <- 1
save_array[,"HarmTOD",] <- 0

## Final output matrix, with non-important factors removed
save_array2 <- save_array[,-c(3,4,5,6,7,8,10),]

lons <- -179.5:179.5
lats <- -89.5:89.5
lonlat <- as.matrix(expand.grid("lons" = lons, "lats" = lats))

## Harmonic day of year for each month
days_of_year_harmonic <- seq(15,365,30)/365*2*pi

## Depths to integrate over to cover top 200m
depths <- 0.5:199.5

for(k in 1:12){ # Loop over months

      print(paste("Now working on month", k,  sep = " "))
      mons_nh <- k # Northern hemisphere month
      mons_sh <- mons_nh + 6 # Southern hemisphere month, corrected to equivalent northern hemisphere month

      if(mons_sh[length(mons_sh)] > 12){ # If after June in SH, equivalent northern hemisphere month in first half of year
        mons_sh <- mons_sh - 12
      }

      ## Import current month sst and chlo climatology
      curr_sst <- t(as.matrix(readRDS(list.files(path='./one/',
                                     pattern=glob2rx(paste("SST*", dimnames(save_array)[[1]][mons_nh], "*", sep = "")), full.names = TRUE))))[,180:1]

      curr_chl <- t(as.matrix(readRDS(list.files(path='./one/',
                                                    pattern=glob2rx(paste("Chl*", dimnames(save_array)[[1]][mons_nh], "*", sep = "")), full.names = TRUE))))[,180:1]

      ## Fill in this month's slice of save_array
      save_array[k,"Lon",] <- lonlat[,1]
      save_array[k,"Lat",] <- lonlat[,2]
      save_array[k,"Bathy",] <- as.vector(bathy_matrix)
      save_array[k,"HarmDOY",c(1:32400)] <- days_of_year_harmonic[mons_sh] # Southern hemisphere months
      save_array[k,"HarmDOY",c(32401:64800)] <- days_of_year_harmonic[mons_nh] # Northern hemisphere months
      save_array[k,"SST",] <- as.vector(curr_sst)
      save_array[k,"Chl",] <- as.vector(curr_chl)

      # Align SST and chlo maps with bathy (where bathy is land, mask sst and chlo)
      save_array[k,"SST", which(is.na(save_array[k,"Bathy",] == TRUE))] <- NA
      save_array[k,"Chl", which(is.na(save_array[k,"Bathy",] == TRUE))] <- NA

      # Convert to dataframe and add random effects factors and 0.5m depth for glmm prediction
      kk <- as.data.frame(t(save_array[k,,]))
      kk$BiomassMethod <- as.factor("Carbon")
      kk$Gear <- as.factor("116")
      kk$Institution <- as.factor("103")
      kk$Depth <- depths[1]

      ####### SURFACE LAYER PREDICTIONS ########
      # Get surface layer estimate
      kk$GLM_Mesozoo <- (10^(predict(m7, newdata = kk, re.form = NA, se.fit = FALSE))) # (re.form = NA sets random effects to zero), use se.fit = TRUE to calculate confidence interval, use bootnsim = 100 for number of bootstrap samples

      ### BELOW SURFACE TO 200M PREDICTIONS
      # Which cells are less than or more than 200m deep?
      kk$bathy_check <- NA
      kk[which(kk$Bathy < 200), "bathy_check"] <- 1
      kk[which(kk$Bathy > 200), "bathy_check"] <- 2

      # Split into shallow and deep data sets
      kk_shallow <- kk %>%
        dplyr::filter(bathy_check == 1)
      kk_deep <- kk %>%
        dplyr::filter(bathy_check == 2)

      ## Do easy deep cells first
      print("Now working on deep cells")
      pb <- txtProgressBar(min = 0, max = length(depths), style = 3)
      for(m in 2:length(depths)){ # Loop over depths
        setTxtProgressBar(pb, m)
        kk_deep$Depth <- depths[m]
        kk_deep$GLM_Mesozoo <- kk_deep$GLM_Mesozoo + (10^(predict(m7, newdata = kk_deep, re.form = NA, se.fit = FALSE))) # re.form = NA sets random effects to zero
      }

      ## Do shallow cells. This is the slow bit. Probably the code is lazily written, but I didn't think I'd need to run it very much so it seemed good enough.
      print("Now working on shallow cells")
      pb <- txtProgressBar(min = 0, max = dim(kk_shallow)[1], style = 3)
      for(n in 1:dim(kk_shallow)[1]){ # Loop over shallow grid cells

        setTxtProgressBar(pb, n)

        curr_kk <- kk_shallow[n,] # Select current grid cell

        curr_depths <- 0.5:(floor(curr_kk$Bathy)+0.5) # Bespoke depth intervals for this cell

        if(is.na(curr_kk$Chl) == FALSE & is.na(curr_kk$SST) == FALSE){ # If you have both sst and chl (some cells don't have both due to seasonality, daylight etc)
          for(m in 2:length(curr_depths)){ # Loop over depths
            curr_kk$Depth <- curr_depths[m]
            kk_shallow[n,"GLM_Mesozoo"] <- kk_shallow[n,"GLM_Mesozoo"] + (10^(predict(m7, newdata = curr_kk, re.form = NA, se.fit = FALSE))) # re.form = NA sets random effects to zero
          }
        }
      }

      ## Join shallow and deep back together
      kk <- left_join(kk, kk_deep %>% dplyr::select(Lon,Lat,GLM_Mesozoo), by = c("Lon","Lat"), suffix = c("","_deep"))
      kk <- left_join(kk, kk_shallow %>% dplyr::select(Lon,Lat,GLM_Mesozoo), by = c("Lon","Lat"), suffix = c("","_shallow"))

      kk <- kk %>% mutate(GLM_Mesozoo = coalesce(GLM_Mesozoo_deep, GLM_Mesozoo_shallow))
      kk <- dplyr::select(kk, -c(GLM_Mesozoo_deep, GLM_Mesozoo_shallow))

      #kk$GLM_Mesozoo <- depth_int_pred

      ## Dump output into this month's time slice of save array2
      save_array2[k,"Lon",] <- as.vector(kk$Lon)
      save_array2[k,"Lat",] <- as.vector(kk$Lat)
      save_array2[k,"Bathy",] <- as.vector(kk$Bathy)
      save_array2[k,"SST",] <- as.vector(kk$SST)
      save_array2[k,"Chl",] <- as.vector(kk$Chl)
      save_array2[k,"GLM_Mesozoo",] <- as.vector(kk$GLM_Mesozoo)
    } # End month loop


# Save the output array
saveRDS(save_array2, file = file.path("Output", "Outputglm_mesozoo_obs_100um.RDS", version = 3))


#
# #### PLOT OUTPUT
# library(rasterImage)
#
# col <- rev(rasterImage::colorPalette(n = 14, type = "jet.colors"))
#
# theme_opts <- list(theme(panel.grid.major = element_line(colour = "transparent"),
#                          panel.background = element_blank(),
#                          plot.background = element_rect(fill="white"),
#                          plot.title = element_text(hjust = 0.5, size = rel(0.5)),
#                          panel.border = element_blank(),
#                          axis.line = element_blank(),
#                          axis.text.x = element_blank(),
#                          axis.text.y = element_blank(),
#                          axis.ticks = element_blank(),
#                          axis.title.x = element_blank(),
#                          axis.title.y = element_blank(),
#                          legend.position = "bottom",
#                          legend.text = element_text(size = rel(1))
# ))
#
# mons <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
# plot_list1 <- plot_list2 <- plot_list3 <- list()
# for(i in 1:12){
#   kk <- data.frame(t(save_array2[i,,]))
#
# plot_list1[[i]] <- ggplot(kk) +  geom_raster(aes(x = Lon, y = Lat, fill = log10(Chl))) +
#   scale_fill_gradientn(name = expression(paste("mg m"^-3)),
#                        colours = rev(col),
#                        position = "bottom",
#                        na.value = "grey80")+ theme_opts + ggtitle(mons[i]) +
#   theme(plot.title = element_text(size = rel(1.5)))
#
# plot_list2[[i]] <- ggplot(kk) +  geom_raster(aes(x = Lon, y = Lat, fill = SST)) +
#   scale_fill_gradientn(name = expression(paste("Celsius")),
#                        colours = rev(col),
#                        position = "bottom",
#                        na.value = "grey80") + theme_opts + ggtitle(mons[i]) +
#   theme(plot.title = element_text(size = rel(1.5)))
#
# plot_list3[[i]] <- ggplot(kk) +  geom_raster(aes(x = Lon, y = Lat, fill = log10(GLM_Mesozoo))) +
#   scale_fill_gradientn(name = expression(paste("log10 Mesozoo Biomass mg m"^-2)),
#                        colours = rev(col),
#                        position = "bottom",
#                        na.value = "grey80")+ theme_opts + ggtitle(mons[i]) +
#   theme(plot.title = element_text(size = rel(1.5)))
#
# }
#
# figure <- ggpubr::ggarrange(plotlist = plot_list1,  common.legend = TRUE, legend = "bottom", nrow = 4, ncol = 3)
# ggpubr::ggexport(figure, filename = paste("./modis_aqua_chlo.jpeg", sep = ""), width = 1200, height =1600)
#
# figure <- ggpubr::ggarrange(plotlist = plot_list2,  common.legend = TRUE, legend = "bottom", nrow = 4, ncol = 3)
# ggpubr::ggexport(figure, filename = paste("./modis_aqua_sst.jpeg", sep = ""), width = 1200, height =1600)
#
# figure <- ggpubr::ggarrange(plotlist = plot_list3,  common.legend = TRUE, legend = "bottom", nrow = 4, ncol = 3)
# ggpubr::ggexport(figure, filename = paste("./glm_mesozoo2.jpeg", sep = ""), width = 1200, height =1600)

