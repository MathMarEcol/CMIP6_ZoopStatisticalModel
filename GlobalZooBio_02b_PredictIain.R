
# Some additional code for another project looking at zooplankton biomass in
# Western Boundary Currents

library(tidyverse)
library(sf)

SpatPlan_Create_Polygon <- function(Limits, cCRS = "EPSG:4326", res = 0.1){
  
  x <- dplyr::tibble(x = seq(as.numeric(Limits["xmin"]), as.numeric(Limits["xmax"]), by = res), y = as.numeric(Limits["ymin"])) %>%
    dplyr::bind_rows(dplyr::tibble(x = as.numeric(Limits["xmax"]), y = seq(as.numeric(Limits["ymin"]), as.numeric(Limits["ymax"]), by = res))) %>%
    dplyr::bind_rows(dplyr::tibble(x = seq(as.numeric(Limits["xmax"]), as.numeric(Limits["xmin"]), by = -res), y = as.numeric(Limits["ymax"]))) %>%
    dplyr::bind_rows(dplyr::tibble(x = as.numeric(Limits["xmin"]), y = seq(as.numeric(Limits["ymax"]), as.numeric(Limits["ymin"]), by = -res)))
  
  x <- x %>%
    as.matrix() %>%
    list() %>%
    sf::st_polygon() %>%
    sf::st_sfc(crs = "EPSG:4326") %>%
    sf::st_transform(crs = cCRS)
}


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


###############################################

# Analyse MEOW 


MEOW <- rbind(st_read("Data/MEOW-subset-Jase.shp") %>% select("PROVINCE") %>% 
                 filter(PROVINCE == "Agulhas") %>% st_union() %>% sf::st_sf() %>% mutate(PROVINCE = "Agulhas"),
               st_read("Data/MEOW-subset-Jase.shp") %>% select("PROVINCE") %>% 
                 filter(PROVINCE %in% c("Warm Temperate Northwest Pacific", "Cold Temperate Northwest Pacific")) %>% 
                 st_union() %>% sf::st_sf() %>% mutate(PROVINCE = "Kuroshio"))


fn <- function(x){
  st_crop(dat, x) %>% st_drop_geometry() %>% filter(!is.na(GLM_Mesozoo_log10)) %>% mutate(MEOW = x$PROVINCE)}

MEOW_Bio <- bind_rows(
  fn(MEOW[1,]),
  fn(MEOW[2,]),
)

out1 <- MEOW_Bio %>% 
  summarise(Mth_Mesozoo_mgC_m3 = mean(10^GLM_Mesozoo_log10, na.rm = TRUE),
            Mth_Mesozoo_sd = sd(10^GLM_Mesozoo_log10, na.rm = TRUE), 
            .by = c(Month, MEOW))

out2 <- MEOW_Bio %>% 
  summarise(Mesozoo_mgC_m3 = mean(10^GLM_Mesozoo_log10, na.rm = TRUE),
            Mesozoo_sd = sd(10^GLM_Mesozoo_log10, na.rm = TRUE), 
            .by = c(MEOW))

gg2 <- ggplot(data = out1, aes(x = Month, y = Mth_Mesozoo_mgC_m3, fill = MEOW), ) +
  geom_col(show.legend = FALSE) +
  geom_errorbar(aes(ymin = Mth_Mesozoo_mgC_m3 - Mth_Mesozoo_sd, ymax = Mth_Mesozoo_mgC_m3 + Mth_Mesozoo_sd), width=.2,
                position = position_dodge(.9)) +
  facet_wrap(facets = "MEOW") +
  theme_bw() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Zooplankton Biomass (mg C m\u00B3)") +
  theme(axis.title.x = element_blank())

ggsave(file.path("Figures", "MEOW_MonthlyBiomass.png"), plot = gg2, dpi = 600)


gg3 <- ggplot(data = out2, aes(x = MEOW, y = Mesozoo_mgC_m3, fill = MEOW), ) +
  geom_col(show.legend = FALSE) +
  geom_errorbar(aes(ymin = Mesozoo_mgC_m3 - Mesozoo_sd, ymax = Mesozoo_mgC_m3 + Mesozoo_sd), width=.2,
                position = position_dodge(.9)) +
  theme_bw() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Zooplankton Biomass (mg C m\u00B3)") +
  theme(axis.title.x = element_blank())

ggsave(file.path("Figures", "MEOW_TotalBiomass.png"), plot = gg3, dpi = 600)

out <- left_join(out1, out2, by = "MEOW")
write_csv(out, file.path("Output","MEOW_ZooBiomass.csv"))

###############################################

# Analyse LME 

LME <- st_read("Data/LME-subset-Jase.shp")

fn <- function(x){
  st_crop(dat, x) %>% st_drop_geometry() %>% filter(!is.na(GLM_Mesozoo_log10)) %>% mutate(LME = x$LME_NAME)}

LME_Bio <- bind_rows(
  fn(LME[1,]),
  fn(LME[2,]),
  fn(LME[3,])
)

out1 <- LME_Bio %>% 
  summarise(Mth_Mesozoo_mgC_m3 = mean(10^GLM_Mesozoo_log10, na.rm = TRUE),
            Mth_Mesozoo_sd = sd(10^GLM_Mesozoo_log10, na.rm = TRUE), 
            .by = c(Month, LME))

out2 <- LME_Bio %>% 
  summarise(Mesozoo_mgC_m3 = mean(10^GLM_Mesozoo_log10, na.rm = TRUE),
            Mesozoo_sd = sd(10^GLM_Mesozoo_log10, na.rm = TRUE), 
            .by = c(LME))

gg2 <- ggplot(data = out1, aes(x = Month, y = Mth_Mesozoo_mgC_m3, fill = LME), ) +
  geom_col(show.legend = FALSE) +
  geom_errorbar(aes(ymin = Mth_Mesozoo_mgC_m3 - Mth_Mesozoo_sd, ymax = Mth_Mesozoo_mgC_m3 + Mth_Mesozoo_sd), width=.2,
                position = position_dodge(.9)) +
  facet_wrap(facets = "LME", ncol = 2) +
  theme_bw() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Zooplankton Biomass (mg C m\u00B3)") +
  theme(axis.title.x = element_blank())

ggsave(file.path("Figures", "LME_MonthlyBiomass.png"), plot = gg2, dpi = 600)


gg3 <- ggplot(data = out2, aes(x = LME, y = Mesozoo_mgC_m3, fill = LME), ) +
  geom_col(show.legend = FALSE) +
  geom_errorbar(aes(ymin = Mesozoo_mgC_m3 - Mesozoo_sd, ymax = Mesozoo_mgC_m3 + Mesozoo_sd), width=.2,
                position = position_dodge(.9)) +
  theme_bw() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Zooplankton Biomass (mg C m\u00B3)") +
  theme(axis.title.x = element_blank())

ggsave(file.path("Figures", "LME_TotalBiomass.png"), plot = gg3, dpi = 600)

out <- left_join(out1, out2, by = "LME")
write_csv(out, file.path("Output","LME_ZooBiomass.csv"))




###############################################

## ORIGINAL ANALYSIS

# Analyse the 5 WBC
WBC <- data.frame(WBC = NULL, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL)
WBC <- bind_rows(WBC,
                 data.frame(WBC = "Agulhas", xmin = 20, xmax = 40, ymin = -38, ymax = -22), # Agulhas region 22-38*S, 20-40*E
                 data.frame(WBC = "EAC", xmin = 140, xmax = 160, ymin = -40, ymax = -20), # EAC region 24-40*S, 140-160*E
                 data.frame(WBC = "Loop-Florida", xmin = -90, xmax = -70, ymin = 24, ymax = 40), # Loop-Florida 24-40*N, 90-70*W
                 data.frame(WBC = "Kuroshio", xmin = 125, xmax = 145, ymin = 24, ymax = 40), # Kuroshio 24-40*N, 125-145*E
                 data.frame(WBC = "Brazil", xmin = -65, xmax = -45, ymin = -40, ymax = -20)) # Brazil Current 24-40*S, 65-45*W

mon <- levels(dat$Month)


WBC_Bio <- bind_rows(
  WBC[1,] %>% SpatPlan_Create_Polygon() %>% st_crop(dat, .) %>% 
    st_drop_geometry() %>% filter(!is.na(GLM_Mesozoo_log10)) %>% mutate(WBC = WBC$WBC[1]),
  WBC[2,] %>% SpatPlan_Create_Polygon() %>% st_crop(dat, .) %>% 
    st_drop_geometry() %>% filter(!is.na(GLM_Mesozoo_log10)) %>% mutate(WBC = WBC$WBC[2]),
  WBC[3,] %>% SpatPlan_Create_Polygon() %>% st_crop(dat, .) %>% 
    st_drop_geometry() %>% filter(!is.na(GLM_Mesozoo_log10)) %>% mutate(WBC = WBC$WBC[3]),
  WBC[4,] %>% SpatPlan_Create_Polygon() %>% st_crop(dat, .) %>% 
    st_drop_geometry() %>% filter(!is.na(GLM_Mesozoo_log10)) %>% mutate(WBC = WBC$WBC[4]),
  WBC[5,] %>% SpatPlan_Create_Polygon() %>% st_crop(dat, .) %>% 
    st_drop_geometry() %>% filter(!is.na(GLM_Mesozoo_log10)) %>% mutate(WBC = WBC$WBC[5])
)

out1 <- WBC_Bio %>% 
  summarise(Mth_Mesozoo_mgC_m3 = mean(10^GLM_Mesozoo_log10, na.rm = TRUE),
            Mth_Mesozoo_sd = sd(10^GLM_Mesozoo_log10, na.rm = TRUE), 
            .by = c(Month, WBC))
out2 <- WBC_Bio %>% 
  summarise(Mesozoo_mgC_m3 = mean(10^GLM_Mesozoo_log10, na.rm = TRUE),
            Mesozoo_sd = sd(10^GLM_Mesozoo_log10, na.rm = TRUE), 
            .by = c(WBC))

gg2 <- ggplot(data = out1, aes(x = Month, y = Mth_Mesozoo_mgC_m3, fill = WBC), ) +
  geom_col(show.legend = FALSE) +
  geom_errorbar(aes(ymin=Mth_Mesozoo_mgC_m3 - Mth_Mesozoo_sd, ymax=Mth_Mesozoo_mgC_m3+Mth_Mesozoo_sd), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(facets = "WBC") +
  theme_bw() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Zooplankton Biomass (mg C m\u00B3)") +
  theme(axis.title.x = element_blank())

ggsave(file.path("Figures", "MonthlyBiomass.png"), plot = gg2, dpi = 600)


gg3 <- ggplot(data = out2, aes(x = WBC, y = Mesozoo_mgC_m3, fill = WBC), ) +
  geom_col(show.legend = FALSE) +
  geom_errorbar(aes(ymin = Mesozoo_mgC_m3 - Mesozoo_sd, ymax = Mesozoo_mgC_m3 + Mesozoo_sd), width=.2,
                position=position_dodge(.9)) +
  theme_bw() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Zooplankton Biomass (mg C m\u00B3)") +
  theme(axis.title.x = element_blank())

ggsave(file.path("Figures", "TotalBiomass.png"), plot = gg3, dpi = 600)


out <- left_join(out1, out2, by = "WBC")

write_csv(out, file.path("Output","WBC_ZooBiomass.csv"))


library(Cairo)
ggsave(file.path("Figures", "WBC_Seasonal.pdf"), plot = gg2, width = 297, height = 210, units = "mm", device = cairo_pdf)
