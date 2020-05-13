## PML-NEDI model -----
# Penman-Menteith-Leuning model was intergrated with NEDI function

# main equations are after Leuning et al. (2008) WRR
# energy allocation is different from Leuning 2008
# soil water content limiting function of soil evaporation is after Morilas et al. 2013
# NEDI function is modified from Change et al. (2018) AFM

## parameters
# parameter need optimizatoin: gmax
# site specific: hc, zr
# parameter fixed: q50= 30 W/m2, d50=0.7,  # Zhang et al. (2010) WRR

# input:daily flux data and LAI
# NEDI should be given separately from savannna and grassland
# output: daily estimated ET, T, Esoil
  
et_model <- function (driver, pars, out='et_total')  
  {
  # parameters needing optimization
  gmax= pars[1] # max stomatal conductance (m/s)
  m   = pars[2] # multifier for NEDI
  
  # driver: fsd, rn, rg, lai, swc, ta, u, vpd
  day_length = driver$day_length   # number of hours in a whole day
  q   = driver$fsd * 0.45  # visible energy unit: W/m2
  lai = driver$lai         # leaf area index, m2/m2
  rg  = driver$rg          # ground heat flux, W/m2    
  rn  = driver$rn          # net radiation, W/m2
  swc = driver$swc         # soil water content, cm3/cm3
  ta  = driver$ta          # air temperature
  u   = driver$u           # wind speed, m/s
  vpd = driver$vpd         # vapor pressure deficit, kPa
  NEDI = driver$NEDI
  a_wavail =  ifelse(!is.na(driver$NEDI) & driver$NEDI > 0, 1, driver$NEDI * m + 1)
  
  ## slope of the relationship between saturated vapor pressure and temperature
  delta   = 17.081 * 234.175 * (6.1078 * exp(17.081 * ta / (ta + 234.175))) / (ta + 234.175)^2 / 10 ## unit kPA / K
  
  ## canopy aerodynamic conductance
  hc  = max(hc_tree, hc_grass) 
  d0  = 2/3 * hc 
  z0  = 0.123 * hc
  z0h = 0.1 * z0
  ra  = (log((zr - d0) / z0)) * (log((zr - d0) / z0h)) / (0.41^2 * u) #  m/s type value allen  1/30 m/s
  
  ## canopy conductance
  gc  =  gmax / kq * log((q + q50) / (q * exp(-kq * lai) + q50)) / (1 + vpd / d50)
  gc  =  gc * a_wavail
  rsc = ifelse(gc ==0, 10^6, 1/gc)
  
  ## total available energy
  aa = rn - rg 
  ## available energy for soil
  as = rn*exp(-kn * lai) - rg
  ## available energy for canopy
  ac = aa - as  
  
  ## tc: canopy transpiration
  tc  = (delta * ac + rho * cp * vpd / ra ) / (delta + gamma*(1 + rsc / ra))
  
  # es: soil evaporation
  epsilon = delta/gamma
  ### soil water constraint function is at the depth of 10 cm
  fswc= ifelse(swc < swc_min, 0, ifelse(swc >swc_max, 1, (swc - swc_min) / (swc_max - swc_min))) 
  es  = fswc * epsilon * as / (epsilon + 1)  

  # total ET
  et_total   = tc + es # unit W/m2
  
  scale_day  = day_length * 3600 / 2.45 * 10^-6 # transform to unit mm/d
  
  ## return results
  et_pre <- switch(out,
               "trans" = tc * scale_day,
               "esoil" = es * scale_day,
               "et_total" = et_total * scale_day) # mm/d
  return(et_pre)
}
