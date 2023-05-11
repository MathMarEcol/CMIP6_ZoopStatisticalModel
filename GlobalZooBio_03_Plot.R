library(tidyverse)
library(sf)
library(patchwork)
library(cmocean)


dat <- readRDS(file = file.path("Output", "glm_mesozoo_obs_100um.RDS")) %>% # Load array
  as.data.frame.table() %>% # Convert to dataframe
  pivot_wider(names_from = "Var2", values_from = "Freq") %>% # Rearrange
  dplyr::rename(Month = Var1) %>% # Rename variables
  # dplyr::filter(!is.na(GLM_Mesozoo)) %>% 
  dplyr::mutate(GLM_Mesozoo_log10 = log10(GLM_Mesozoo)) %>% 
  dplyr::select(Lon, Lat, GLM_Mesozoo_log10, Month) %>% 
  pivot_wider(names_from = Month, values_from = GLM_Mesozoo_log10) %>% 
  terra::rast(crs = "EPSG:4326") %>% 
  terra::as.polygons(dissolve = FALSE, na.rm = FALSE) %>% 
  st_as_sf() %>% 
  pivot_longer(cols = !tidyselect::starts_with("geometry"), names_to = "Month", values_to = "GLM_Mesozoo_log10") %>% 
  mutate(Month = factor(Month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")))

land <- rnaturalearth::ne_countries(returnclass = "sf")


gg <- ggplot() +
  geom_sf(data = dat, aes(fill = GLM_Mesozoo_log10, colour = GLM_Mesozoo_log10)) + 
  # labs(title = Month) +
  facet_wrap("Month", nrow = 4) + 
  scale_fill_cmocean(name = "amp", 
                     aesthetics = c("colour", "fill"), 
                     limits = c(2.2, 3.6), 
                     breaks = seq(2.2, 3.6, 0.2),
                     oob = scales::squish,
                     guide = guide_colourbar(title.position = "right", 
                                             title.theme = element_text(angle = 270), 
                                             title.vjust = 0.5,
                                             title = "Zooplankton Biomass (log\u2081\u2080 mg m\u00B3)")) + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) + 
  theme_bw() + 
  theme(legend.key.height = unit(0.06, "npc"), legend.key.width = unit(0.01, "npc")) + 
  geom_sf(data = land, colour = "grey20", fill = "grey20")


library(Cairo)
ggsave(file.path("Figures", "GlobalMap.pdf"), width = 420, height = 297, units = "mm", device = cairo_pdf)



# Some additional code for another project looking at zooplankton biomass in
# Western Boundary Currents
# 
# SpatPlan_Create_Polygon <- function(Limits, cCRS = "EPSG:4326", res = 0.1){
#   
#   x <- dplyr::tibble(x = seq(as.numeric(Limits["xmin"]), as.numeric(Limits["xmax"]), by = res), y = as.numeric(Limits["ymin"])) %>%
#     dplyr::bind_rows(dplyr::tibble(x = as.numeric(Limits["xmax"]), y = seq(as.numeric(Limits["ymin"]), as.numeric(Limits["ymax"]), by = res))) %>%
#     dplyr::bind_rows(dplyr::tibble(x = seq(as.numeric(Limits["xmax"]), as.numeric(Limits["xmin"]), by = -res), y = as.numeric(Limits["ymax"]))) %>%
#     dplyr::bind_rows(dplyr::tibble(x = as.numeric(Limits["xmin"]), y = seq(as.numeric(Limits["ymax"]), as.numeric(Limits["ymin"]), by = -res)))
#   
#   x <- x %>%
#     as.matrix() %>%
#     list() %>%
#     sf::st_polygon() %>%
#     sf::st_sfc(crs = "EPSG:4326") %>%
#     sf::st_transform(crs = cCRS)
# }
# 
# 
# 
# # Analyse the 5 WBC for Iain
# 
# WBC <- data.frame(WBC = NULL, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL)
# WBC <- bind_rows(WBC, 
#                  data.frame(WBC = "Agulhas", xmin = 20, xmax = 40, ymin = -38, ymax = -22), # Agulhas region 22-38*S, 20-40*E
#                  data.frame(WBC = "EAC", xmin = 140, xmax = 160, ymin = -40, ymax = -20), # EAC region 24-40*S, 140-160*E
#                  data.frame(WBC = "Loop-Florida", xmin = -90, xmax = -70, ymin = 24, ymax = 40), # Loop-Florida 24-40*N, 90-70*W
#                  data.frame(WBC = "Kuroshio", xmin = 125, xmax = 145, ymin = 24, ymax = 40), # Kuroshio 24-40*N, 125-145*E
#                  data.frame(WBC = "Brazil", xmin = -65, xmax = -45, ymin = -40, ymax = -20)) # Brazil Current 24-40*S, 65-45*W
# 
# mon <- levels(dat$Month)
# 
# 
# WBC_Bio <- bind_rows(
#   WBC[1,] %>% SpatPlan_Create_Polygon() %>% st_crop(dat, .) %>% st_drop_geometry() %>% filter(!is.na(GLM_Mesozoo_log10)) %>% mutate(WBC = WBC$WBC[1]),
#   WBC[2,] %>% SpatPlan_Create_Polygon() %>% st_crop(dat, .) %>% st_drop_geometry() %>% filter(!is.na(GLM_Mesozoo_log10)) %>% mutate(WBC = WBC$WBC[2]),
#   WBC[3,] %>% SpatPlan_Create_Polygon() %>% st_crop(dat, .) %>% st_drop_geometry() %>% filter(!is.na(GLM_Mesozoo_log10)) %>% mutate(WBC = WBC$WBC[3]),
#   WBC[4,] %>% SpatPlan_Create_Polygon() %>% st_crop(dat, .) %>% st_drop_geometry() %>% filter(!is.na(GLM_Mesozoo_log10)) %>% mutate(WBC = WBC$WBC[4]),
#   WBC[5,] %>% SpatPlan_Create_Polygon() %>% st_crop(dat, .) %>% st_drop_geometry() %>% filter(!is.na(GLM_Mesozoo_log10)) %>% mutate(WBC = WBC$WBC[5])
# ) %>% 
#   # mutate(Area = st_area(.)) %>% 
#   summarise(GLM_Mesozoo = mean(10^GLM_Mesozoo_log10, na.rm = TRUE), .by = c(Month, WBC))
# 
# gg2 <- ggplot(data = WBC_Bio, aes(x = Month, y = GLM_Mesozoo_log10, fill = WBC), ) + 
#   geom_col(show.legend = FALSE) +
#   facet_wrap(facets = "WBC") + 
#   theme_bw() + 
#   scale_x_discrete(expand = c(0,0)) +
#   scale_y_continuous(expand = c(0,0)) +
#   ylab("Zooplankton Biomass (mg m\u00B3)") + 
#   theme(axis.title.x = element_blank())
# 
# 
# library(Cairo)
# ggsave(file.path("Figures", "WBC_Seasonal.pdf"), plot = gg2, width = 297, height = 210, units = "mm", device = cairo_pdf)

