## SW-NEDI model
# Shuttleworth-Wallace 1985 model was intergrated with NEDI function
# NEDI function is modified from Change et al. (2018) AFM
# Add water availability limitation on canopy conductance by Wei 2019/10/25

## Parameters
# parameter need optimizatoin: gmax
# paramter site specific: hc, zr
# parameter fixed: q50= 50, d50=0.7,  # Zhang et al. (2010) WRR

# input:daily flux data and LAI
# output: daily estimated ET, T, Esoil

et_model <-function(driver, pars, out='et_total'){
  # parameters needing optimization
    gmax = pars[1] # max stomatal conductance (m/s)
    m    = pars[2] # multifier for NEDI   
      
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
	
  # available energy partitioning
	aa   = rn - rg
	asubs= rn*exp(-kq * lai) - rg

  ## slope of the relationship between saturated vapor pressure and temperature
    delta   = 17.081 * 234.175 * (6.1078 * exp(17.081 * ta / (ta + 234.175))) / (ta + 234.175)^2 / 10 ## unit kPA / K

  # 5 resistances
	# r1: raa, aerodynamic resistance
	hc   = max(hc_tree, hc_grass)
	disp =  hc * 2/3  # zero-plane displacement (m)
	z0   = 0.123 * hc # roughness parameter (m)  
	ustar= kr * u / (log((zr - disp) / z0))   ## kr is Karman constance
	kh   = kr * ustar * (hc - disp)
	raa  = log((zr - disp) / (hc - disp)) / (kr * ustar) + (hc / (km * kh)) * (-1 + exp(km * (hc - disp - z0) / hc))	

	# r2: ras, in-canopy aerodynamic resistance
	ras  = (hc * exp(km) / (km * kh)) * (exp(-km * z0gs / hc) - exp(-km * (z0 + disp) / hc))

	# r3: rbc, boundary-layer resistance
	uh   = (ustar / kr) * log((hc - disp) / z0)
	rbc  = 70 / lai * (lwidth / uh)^0.5

	# r4: rsc, surface canopy resistance 
	gc   = gmax / kq * log((q + q50) / (q * exp(-kq * lai) + q50)) / (1 + vpd / d50)
	gc   = gc * a_wavail
	gc   = ifelse(gc == 0, 0.000001, gc)
	rsc  = 1/gc

	# r5: rss, surface soil resistance
	## soil water constraint function:
	fswc = ifelse(swc <= swc_min, 0.000001, ifelse(swc >swc_max, 1, (swc - swc_min) / (swc_max - swc_min)))
	rss  = b1*(1/fswc)^b2 + b3

  # Two source SW model calculation
	rs   = (delta + gamma) * ras + gamma * rss
	rc   = (delta + gamma) * rbc + gamma * rsc
	ra   = (delta + gamma) * raa
	ccs  = 1 / (1 + rs * ra / (rc * (rs + ra)))
	ccc  = 1 / (1 + rc * ra / (rs * (rc + ra)))
  
    pmc  = (delta * aa +(rho * cp * vpd - delta * rbc * asubs)/(raa+rbc)) /(delta + gamma * (1 + rsc/(raa+rbc)))  
    pms  = (delta * aa +(rho * cp * vpd - delta * ras *(aa-asubs))/(raa+ras))/(delta+gamma*(1+rss/(raa+ras)))
   
  # total ET
    et_total   = ccc * pmc + ccs * pms # unit W/m2
    
  # unit scale for ET 
    scale_day  =  day_length * 3600 / 2.45 * 10^-6 # transform to unit mm/d
  
  ## return results
  et_pre <- switch(out,
               "tc" = tc * scale_day,
               "es" = es * scale_day,
               "et_total" = et_total * scale_day) # mm/d
  
  return(et_pre)
}
