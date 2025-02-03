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



