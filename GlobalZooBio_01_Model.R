# Zooplankton Biomass Models from COPEPOD
# Jason Everett and Anthony Richardson (UQ)
# Created: 26th September 2019 
# Last Updated: 26th July 2022 

# Model Used for Heneghan et al 2020 and Petrik et al (submitted) is m7

############################  Preliminaries ############################  
library(tidyverse)
library(visreg)
library(splines)
library(lme4)
library(lattice)

source("uils.R") # Functions for Harmonic and plotting Biomass LMs

dat <- readRDS("Data/GlobalBiomassData.rds") # n = 197,413

min_val <- min(dat$Biomass[dat$Biomass>0])/2

## Reduce (but don't remove) some extreme values
dat <- dat %>% 
  mutate(
    HarmTOD = (TimeLocal/24)*2*pi, # Convert to radians
    HarmDOY = (DOY2/365)*2*pi, # Convert to radians
    Latitude2 = abs(Latitude),
    Mesh = replace(Mesh, Mesh > 1000, 1000),
    Depth = replace(Depth, Depth > 1500, 1500),
    Bathy = replace(Bathy, Bathy > 7000, 7000),
    SST = replace(SST, SST > 31, 31),
    Biomass = replace(Biomass, Biomass > 10000, 10000),
    Biomass = Biomass + min_val) %>% 
  unite(Gear_Mesh, Gear, Mesh, remove = FALSE) %>%  # Calculate Gear_Mesh
  mutate(Gear_Mesh = as.factor(Gear_Mesh)) %>% 
  droplevels()

# Investigate factors converning institutions/projects/gear
nlevels(dat$ShpCruise)
nlevels(dat$Gear)
nlevels(dat$Project)
nlevels(dat$Institution)
nlevels(dat$DatasetID)

table(dat$Gear, dat$Project)
table(dat$Gear, dat$Institution)
table(dat$Gear, dat$DatasetID)

sum(table(dat$Gear, dat$Project)>0) / (nlevels(dat$Gear) * nlevels(dat$Project)) * 100
sum(table(dat$Gear, dat$Institution)>0) / (nlevels(dat$Gear) * nlevels(dat$Institution)) * 100 # Best data coverage, but still only 2.5%
sum(table(dat$Gear, dat$DatasetID)>0) / (nlevels(dat$Gear) * nlevels(dat$DatasetID)) * 100

# How many Projects use Multiple Gears?
# Most Projects have one Gear, but several Projects use multiple Gears
tempory <- dat %>% 
  group_by(Project) %>% 
  summarise(N = length(unique(Gear)))

tempory <- dat %>% 
  group_by(Institution) %>% 
  summarise(N = length(unique(Gear)))
write_csv(tempory, "InstitutionxGear.csv")
  
####################
##### MODELS
# Just Gear
m1 <- lmer(log10(Biomass) ~ BiomassMethod + Mesh + log10(Chl) + exp(-Depth/1000) + 
             fHarmonic(HarmTOD, k = 1) + ns(Bathy, df = 3) + 
             fHarmonic(HarmDOY, k = 1) * ns(SST, 3) +
             (1|Gear), 
           data = dat)
# r.squaredGLMM(m1) # Fixed r2 = 83.8%, Total r2 = 90.5%
anova(m1)
summary(m1)
fPlotBiomassLM(m1, "Biomass_lm1") # Bathy goes up offshore (does not asymptote)
RE <- ranef(m1, whichel = "Gear", condVar = TRUE)
x11(height = 8, width = 5)
dotplot(RE, lattice.options = list(layout = c(1,2,3))) # Looks reasonable
saveRDS(m1, "Output/m1.rds")
plot(m1) # Residuals reasonable
qqnorm(residuals(m1))
qqline(residuals(m1)) # Normality not that good

####################
# Would a family = Gamma(link = log) be better? 
# Seems to be worse for fitting high biomass values
m2 <- glmer(Biomass ~ BiomassMethod + Mesh + log10(Chl) + exp(-Depth/1000) + 
              fHarmonic(HarmTOD, k = 1) + ns(Bathy, df = 3) + 
              fHarmonic(HarmDOY, k = 1) * ns(SST, 3) +
              (1|Gear), family = Gamma(link = log), nAGQ = 0, data = dat)
# r.squaredGLMM(m2) # lognormal: Fixed r2 = 65.2%, Total r2 = 87.2%; trigamma (recommended in help): Fixed r2 = 22.8%, Total r2 = 30.5%
anova(m2)
summary(m2)
fPlotBiomassLM(m2, "Biomass_lm2") # Bathy goes up offshore (does not asymptote)
RE <- ranef(m2, whichel = "Gear", condVar = TRUE)
x11(height = 8, width = 5)
dotplot(RE, lattice.options = list(layout = c(1,2,3))) # Looks reasonable
saveRDS(m2, "Output/m2.rds")
plot(m2) # Residuals?
qqnorm(residuals(m2))
qqline(residuals(m2)) # Normality not good for higher biomass values

# Go with Normal error structure with log10 response
# Which random effects make most sense?


####################
## Try Gear and Institution
m3 <- lmer(log10(Biomass) ~ BiomassMethod + Mesh + log10(Chl) + exp(-Depth/1000) +
             fHarmonic(HarmTOD, k = 1) + ns(Bathy, df = 3) + 
             fHarmonic(HarmDOY, k = 1) * ns(SST, 3) +
             (1|Gear) + (1|Institution), 
           data = dat)
summary(m3) # Both Gear and Institution important (Gear slightly more imp.) - Best model so far
# r.squaredGLMM(m3) # Fixed r2 = 82.0%, Total r2 = 91.4%
anova(m3)
fPlotBiomassLM(m3, "Biomass_lm3") # Looks good and Bathy is better (asymptotes offshore)
saveRDS(m3, "Output/m3.rds")
## *** m3 is best model so far

####################
## Try Gear and Project
m4 <- lmer(log10(Biomass) ~ BiomassMethod + Mesh + log10(Chl) + exp(-Depth/1000) +
             fHarmonic(HarmTOD, k = 1) + ns(Bathy, df = 3) + 
             fHarmonic(HarmDOY, k = 1) * ns(SST, 3) +
             (1|Gear) + (1|Project), 
           data = dat)
summary(m4) # Project more important than Gear
# r.squaredGLMM(m4) # Fixed r2 = 82.0%, Total r2 = 91.4%
anova(m4)
fPlotBiomassLM(m4, "Biomass_lm4") # Bathy becomes somewhat unimodal peaking at 3000 m
saveRDS(m4, "Output/m4.rds")

####################
# Try Gear and ShpCruise
m5 <- lmer(log10(Biomass) ~ BiomassMethod + Mesh + log10(Chl) + exp(-Depth/1000) +
             fHarmonic(HarmTOD, k = 1) + ns(Bathy, df = 3) + 
             fHarmonic(HarmDOY, k = 1) * ns(SST, 3) +
             (1|Gear) + (1|ShpCruise), 
           data = dat)
summary(m5) # Gear more important than ShpCruise (drop ShpCruise)
# r.squaredGLMM(m5) # Fixed r2 = 82.0%, Total r2 = 91.4%
anova(m5)
fPlotBiomassLM(m5, "Biomass_lm5") # Bathy becomes very unimodal peaking at 3000 m and seasonality increases ijn the Tropics
saveRDS(m5, "Output/m5.rds")

####################
# Try Gear and DatasetId
m6 <- lmer(log10(Biomass) ~ BiomassMethod + Mesh + log10(Chl) + exp(-Depth/1000) +
             fHarmonic(HarmTOD, k = 1) + ns(Bathy, df = 3) + 
             fHarmonic(HarmDOY, k = 1) * ns(SST, 3) +
             (1|Gear) + (1|DatasetID), 
           data = dat)
summary(m6) # DatasetID more important than Gear
# r.squaredGLMM(m6) # Fixed r2 = 82.0%, Total r2 = 91.4%
anova(m6)
fPlotBiomassLM(m6, "Biomass_lm6") # Bathy becomes very unimodal peaking at 3000 m and seasonality increases ijn the Tropics
# x11(height = 8, width = 12)
# RE1 <- ranef(m6, whichel = "DatasetId", condVar = TRUE)
# dotplot(RE1, lattice.options = list(layout = c(1,2,3)))
# x11(height = 8, width = 6)
# RE2 <- ranef(m6, whichel = "Gear", condVar = TRUE)
# dotplot(RE2, lattice.options = list(layout = c(1,2,3)))
saveRDS(m6, "Output/m6.rds")

####################
# So use Gear and Institution because has highest amount of joint data and reasonable model
## *** m3 best so far
# Add interaction between Depth and TOD
m7 <- lmer(log10(Biomass) ~ BiomassMethod + Mesh + log10(Chl) + 
             exp(-Depth/1000) * fHarmonic(HarmTOD, k = 1) + 
             ns(Bathy, df = 3) +
             fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             (1|Gear) + 
             (1|Institution), data = dat)

summary(m7)
# r.squaredGLMM(m7) # Fixed r2 = 82%, Total r2 = 91%
anova(m7)
fPlotBiomassLM(m7, "Biomass_lm7") # Looks good with Depth*TOD
saveRDS(m7, "Output/m7.rds")

plot(m7) # Residuals?
qqnorm(residuals(m7))
qqline(residuals(m7))

# RE1 <- ranef(m7, whichel = "Institution", condVar = TRUE)
# RE2 <- ranef(m7, whichel = "Gear", condVar = TRUE)
# x11(height = 8, width = 12)
# dotplot(RE1, lattice.options = list(layout = c(1,2,3)))
# x11(height = 8, width = 6)
# dotplot(RE2, lattice.options = list(layout = c(1,2,3)))

## *** m7 is new best model



####################
# Try k = 2 for TOD and DOY
m8 <- lmer(log10(Biomass) ~ BiomassMethod + Mesh + log10(Chl) +
             exp(-Depth/1000) * fHarmonic(HarmTOD, k = 2) + 
             ns(Bathy, df = 3) + 
             fHarmonic(HarmDOY, k = 2) * ns(SST, 3) + 
             (1|Gear) + (1|Institution), 
           data = dat)
summary(m8)
# r.squaredGLMM(m8) # Fixed r2 = 82%, Total r2 = 91%
anova(m8)
fPlotBiomassLM(m8, "Biomass_lm8") # With k = 2 for DOY, makes surface go wavy at Poles
# With k = 2 for TOD, makes deeper Z go up during the afternoon
saveRDS(m8, "Output/m8.rds")



####################
# Try k = 1 for DOY and k = 3 for TOD
m9 <- lmer(log10(Biomass) ~ BiomassMethod + Mesh + log10(Chl) +
             exp(-Depth/1000) * fHarmonic(HarmTOD, k = 3) + 
             ns(Bathy, df = 3) + 
             fHarmonic(HarmDOY, k = 1) * ns(SST, 3) + 
             (1|Gear) + (1|Institution), 
           data = dat)
summary(m9)
# r.squaredGLMM(m9) # Fixed r2 = 82%, Total r2 = 91%
anova(m9)
fPlotBiomassLM(m9, "Biomass_lm9") # Shows midnight sinking, but probably best to stay simple
saveRDS(m9, "Output/m9.rds")
