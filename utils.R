# Harmonic is to fit Hour and DOY

fHarmonic <- function (theta, k = 4) {
  X <- matrix(0, length(theta), 2 * k)
  nam <- as.vector(outer(c("c", "s"), 1:k, paste, sep = ""))
  dimnames(X) <- list(names(theta), nam)
  m <- 0
  for (j in 1:k) {
    X[, (m <- m + 1)] <- cos(j * theta)
    X[, (m <- m + 1)] <- sin(j * theta)
  }
  X
}


# Function to plot the linear models using visreg
fPlotBiomassLM <- function (mdl, Name) {
  
  # If lm use this
  # Terms <- as.character(m1$terms)[3] # Terms from the model so we can print blank if n.s.
  
  # If lmer
  Terms <- as.character(mdl@call)[2]
  
  x11(width = 12, height = 6)
  if (length(mdl@frame) <= 8){r <- 2} else{r <- 3}
  par(mfrow = c(r,4), mar = c(4,4,2,2))
  
  if(grepl("BiomassMethod", Terms, fixed = TRUE)) {
    visreg(mdl, "BiomassMethod", rug = FALSE, scale = "response", xlab = "Method", ylab = expression("log"[10]*"(Biomass)"))}
  
  
  if(grepl("Tow", Terms, fixed = TRUE)) {
    visreg(mdl, "Tow", rug = FALSE, scale = "response", xlab = "Tow", ylab = expression("log"[10]*"(Biomass)"))}
  
  
  if(grepl("Mesh", Terms, fixed = TRUE)) { 
    visreg(mdl, "Mesh", scale = "response", xlab = "Mesh  (microns)", ylab = expression("log"[10]*"(Biomass)"))}
  
  
  if(grepl("Chl", Terms, fixed = TRUE)) {
    visreg(mdl, "Chl", scale = "response", xlab = expression("Chl-a (mg m"^-3*")"), ylab = expression("log"[10]*"(Biomass)"))
    } 
  
  
  if(grepl("Bathy", Terms, fixed = TRUE)) {
    visreg(mdl, "Bathy", scale = "response", xlab = "Bathy (m)", ylab = expression("log"[10]*"(Biomass)"))}
  
  if(grepl("HarmTOD", Terms, fixed = TRUE)) {
    visreg(mdl, "HarmTOD", rug = FALSE, scale = "response", xlab = "Time of Day", xaxt = 'n', ylab = expression("log"[10]*"(Biomass)"))
    axis(side=1, at=c(0, pi/2 , pi, pi + pi/2, pi*2), labels=c("00:00","06:00","12:00","18:00","00:00"))
  }
  
  if(grepl("Depth", Terms, fixed = TRUE)) {
    visreg(mdl, "Depth", scale = "response", xlab = "Depth", ylab = expression("log"[10]*"(Biomass)"))}
  
  
  if(grepl("HarmDOY", Terms, fixed = TRUE)) {
    visreg(mdl, "HarmDOY", scale = "response", xlab = "Day of Year", xaxt = 'n', ylab = expression("log"[10]*"(Biomass)"))
    axis(side=1, at=c(0, pi/2 , pi, pi + pi/2, pi*2), labels=c("1","91","182","273","365"))
  }
  
  if(grepl("SST", Terms, fixed = TRUE)) {
    visreg(mdl, "SST", scale = "response", xlab = "SST (ºC)", ylab = expression("log"[10]*"(Biomass)"))}
  
  
  # if(grepl("fHarmonic\\(HarmDOY, k = \\d\\) \\* ns\\(SST, \\d\\)", Terms)){
  #   visreg2d(mdl, yvar = "HarmDOY", xvar = "SST", scale = "response", 
  #          plot.type = "persp", theta = 45, phi = 10, r = 100, 
  #          type = "conditional", 
  #          ticktype = "detailed", xlab = "\nSST (ºC)", 
  #          ylab = "\nDay of Year", zlab = "\nlog10(Biomass)", 
  #          color = "deepskyblue2")}
  
  if(grepl("fHarmonic\\(HarmDOY, k = \\d\\) \\* ns\\(SST, \\d\\)", Terms)){
  visreg(mdl, "HarmDOY", by = "SST",
         type = "conditional",
         scale = "response",
         overlay = TRUE, rug = 0, 
         breaks = c(2, 15, 30), 
         xlab = "Day of Year", 
         ylab = expression("log"[10]*"(Biomass)"),
         strip.names = c("2 ºC", "15 ºC", "30 ºC"),
         xaxt = 'n')
    axis(side=1, at=c(0, pi/2 , pi, pi + pi/2, pi*2), labels=c("1","91","182","273","365"))
  }
  
  
  if(grepl('exp\\(-Depth\\/1000\\) \\* fHarmonic\\(HarmTOD, k = \\d\\)', Terms) |
     grepl('fHarmonic\\(HarmTOD, k = \\d\\) \\* exp\\(-Depth\\/1000', Terms)) {
    visreg(mdl, "HarmTOD", by = "Depth", breaks = c(0, 100, 500), 
           xlab = "Time of Day", ylab = expression("log"[10]*"(Biomass)"),
           type = "conditional", scale = "response", overlay = TRUE, rug = 0, 
           strip.names = c("Depth=0","Depth=100", "Depth=500"), xaxt = 'n')
    axis(side=1, at=c(0, pi/2 , pi, pi + pi/2, pi*2), labels=c("00:00","06:00","12:00","18:00","00:00"))
    }
  
  
  ### DO RANDOM EFFECTS
  if(grepl('(1 | DatasetId)', Terms, fixed=TRUE )){
    RE <- ranef(mdl, whichel = "DatasetId", condVar = TRUE)
    dotchart(RE$DatasetId$`(Intercept)`, ylab = "DatasetId")
  }
  
  if(grepl('(1 | Gear)', Terms, fixed=TRUE )){
    RE <- ranef(mdl, whichel = "Gear", condVar = TRUE)
    dotchart(RE$Gear$`(Intercept)`, ylab = "Gear")
    # labels = sort(unique(dat$Gear)))
  }
  
  if(grepl('(1 | Institution)', Terms, fixed=TRUE )){
    RE <- ranef(mdl, whichel = "Institution", condVar = TRUE)
    dotchart(RE$Institution$`(Intercept)`, ylab = "Institution")
    # labels = sort(unique(dat$Institution))
    
  }
  
  if(grepl('(1 | ShpCruise)', Terms, fixed=TRUE )){
    RE <- ranef(mdl, whichel = "ShpCruise", condVar = TRUE)
    dotchart(RE$ShpCruise$`(Intercept)`, ylab = "ShpCruise")
  }
  
  dev.print(pdf, paste0("Figures/", Name, ".pdf"))
}


