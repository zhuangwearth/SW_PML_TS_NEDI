## TS-NEDI model
# The three-source model (TS) with drought indicator NEDI
# based on the three-source (TS) model, drought indicators were added to TS to constrain transpiratoin
# The TS model is after Lhomme et al. (2013).

## parameter
# set the same gmax for grass (herb) and tree (wood) 
# two parameters need optimizatoin: gmax, m
# paramter site specific: hc, zr
# parameter fixed: q50= 50, d50=0.7,  # Zhang 2010 WRR

# input:daily flux data and LAI
# output: daily estimated t_tree, t_grass, e_soil and et_total

## model equations:
# PM - Penman-Monteith latent heat flux density (W/m2)
PM <- function(AA, VPD, Delta, ra, rs)
             {
              Pm <- ( Delta * AA + cp * rho * VPD / ra) / (Delta + gamma * (1 + rs/ra))
              return(Pm)
}

et_model <-function(driver, pars, out='et_total'){
  # parameters needing optimization
    gmax  = pars[1]    # max stomatal conductance (m/s)
	gmax_tree = gmax
	gmax_grass = gmax
	m   = pars[2]      # multifier for NEDI
  
  # driver: 
	day_length = driver$day_length
	q   = driver$fsd *0.45   # visible radiation
	lai = driver$lai
	lai_tree = driver$lai_tree
	lai_grass = driver$lai_grass
	rg  = driver$rg
	rn  = driver$rn
	swc = driver$swc
	ta  = driver$ta
	u   = driver$u
	vpd = driver$vpd 
	NEDI_tree  = driver$NEDI_tree
	NEDI_grass = driver$NEDI_grass
    a_wavail_tree  = ifelse(!is.na(driver$NEDI_tree) & driver$NEDI_tree > 0, 1, driver$NEDI_tree * m + 1)
    a_wavail_grass = ifelse(!is.na(driver$NEDI_grass) & driver$NEDI_grass > 0, 1, driver$NEDI_grass * m + 1)

######## BOTH tree and grass
if (hc_tree != 0) 
        {
      # available energy partitioning
        ## total available energy
        AA = rn - rg 
        ## available energy for soil
        A3 = rn*exp(-kn * lai) - rg
        ## available energy for overstory tree layer
        A1 = rn*(1- exp(-kn * lai_tree))
        ## available energy for understory grass layer  
        A2 = AA - A1 - A3

      ## slope of the relationship between saturated vapor pressure and temperature
        delta = 17.081 * 234.175 * (6.1078 * exp(17.081 * ta / (ta + 234.175))) / (ta + 234.175)^2 / 10 ## unit kPA / K

      # 8 resistances
        ## r1: ra0, aerodynamic resistance
        hc   = hc_tree
        disp = 2/3  * hc # zero-plane displacement (m)
        z0   = 0.123* hc # roughness parameter (m)  
        ustar= kr * u / (log((zr - disp) / z0)) #  kr is Karman constance
        kh   = kr * ustar * (hc - disp) 
        ra0  = log((zr - disp) / (hc - disp)) / (kr * ustar) + (hc / (km * kh)) * (-1 + exp(km * (hc - disp - z0) / hc))	

        ## r2: rau after Dolman et al 1993
        # aerodynamic resistance between grass layer and mean source height
        z0u  = 0.123  * hc_grass + 2/3  * hc_grass # roughness parameter (m)   
        ra2  = (hc * exp(km) / (km * kh)) * (exp(-km * z0u / hc) - exp(-km * (z0 + disp) / hc))

        ## r3: rb1, boundary-layer resistance of trees
        u_tree   = (ustar / kr) * log((hc - disp) / z0) # wind speed at the wood canopy
        rb1   = 70 * (lwidth / u_tree) ^ 0.5 / lai_tree 

        ## r4: rb2, boundary-layer resistance of grasses
        u_grass   = u_tree * exp(km * (hc_grass / hc - 1))
        rb2   = ifelse(lai_grass ==0, 10^10, 70 * (lwidth / u_grass) ^ 0.5 / lai_grass)

        ## r5: rs1, tree canopy surface resistance 
        gc_tree  = gmax_tree / kq * log((q + q50) / (q * exp(-kq * lai_tree) + q50)) / (1 + vpd / d50)
        gc_tree  = gc_tree * a_wavail_tree # constrain tree Gc with NEDI
        rs1  = ifelse(gc_tree==0, 10^6, 1/gc_tree)

        ## r6: rs2, rs_grass, grass canopy surface resistance 
        q    = q * exp(-kq * lai_tree)
        gc_grass   = gmax_grass / kq * log((q + q50) / (q * exp(-kq * lai_grass) + q50)) / (1 + vpd / d50)
        gc_grass  = gc_grass * a_wavail_grass # constrain grass Gc with NEDI
        gc_grass  = ifelse(gc_grass <= 0, 0.000000001, gc_grass)
        rs2  = 1/gc_grass

        # r7: ra3, ras, aerodynamic resistance between soil surface and mean source height
        ra3  = (hc * exp(km) / (km * kh)) * (exp(-km * z0gs / hc) - exp(-km * (z0 + disp) / hc))

        # r8: rs3, rss, soil surface resistance
        # soil water constraint function:
        fswc = ifelse(swc < swc_min, 0.000001, ifelse(swc >swc_max, 1, (swc - swc_min) / (swc_max - swc_min)))
        rs3  = b1*(1/fswc)^b2 + b3

      # potential ET
        Ep = (delta * AA + rho * cp * vpd/ra0)/(delta + gamma)
        # epsilon = delta/gamma
        # change to Shuttleworth 1993
        # Ep  = epsilon * (rn - rg) / (epsilon + 1)  +  6.43 * (1 + 0.536*u) *vpd  / (epsilon + 1)   # unit W/m2 
        
      # Three source model calculation	
      # Calculate total ET
        ra1 = rb1
        ra2 = rb2 + ra2 
        ra3 = ra3

        R0 = (1 + delta/gamma) * ra0
        R1 = rs1 + (1 + delta/gamma) * ra1
        R2 = rs2 + (1 + delta/gamma) * ra2
        R3 = rs3 + (1 + delta/gamma) * ra3

        P1 = 1/R1/(1 + R0 * (1/R1 + 1/R2 + 1/R3))  
        P2 = 1/R2/(1 + R0 * (1/R1 + 1/R2 + 1/R3))  
        P3 = 1/R3/(1 + R0 * (1/R1 + 1/R2 + 1/R3)) 

        # total ET
        ET = R0 * (P1 + P2 + P3) * Ep + delta * (P1 * A1 * ra1 + P2 * A2 * ra2 + P3 * A3 * ra3)/gamma

      # separate components	
        # tree transpiration
        T1 = (R0 * (Ep - ET) + delta/gamma*ra1*A1)/R1 
        # grass transpiration
        T2  = ifelse(gc_grass > 0, (R0 * (Ep - ET) + delta/gamma*ra2*A2)/R2, 0)
        # soil evaporation 
        ES = (R0 * (Ep - ET) + delta/gamma*ra3*A3)/R3

      # unit scale for ET 	
        scale_day =  day_length * 3600 / 2.45 * 10^-6  # transform to unit mm/d

      ## return results
      et_pre <- switch(out,
                   "t_tree"  = T1 * scale_day,
                   "t_grass" = T2 * scale_day,
                   "e_soil"  = ES * scale_day,
                   "et_total"= ET * scale_day,
                   "Ep"      = Ep * scale_day,
                   "EpAA"    = EpAA * scale_day)
    }
    

######## ONLY grass
 else #if (hc_tree == 0 )
        {
      # available energy partitioning
        aa   = rn - rg
        asubs= rn*exp(-kq * lai) - rg

      ## slope of the relationship between saturated vapor pressure and temperature
        delta   = 17.081 * 234.175 * (6.1078 * exp(17.081 * ta / (ta + 234.175))) / (ta + 234.175)^2 / 10 ## unit kPA / K

      # 5 resistances
        # r1: raa, aerodynamic resistance
        hc   =  max(hc_tree, hc_grass)
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
        gc   = gmax_grass / kq * log((q + q50) / (q * exp(-kq * lai) + q50)) / (1 + vpd / d50)
        gc   = gc * a_wavail_grass # constrain grass Gc with NEDI
        rsc  = 1/gc

        # r5: rss, surface soil resistance
        ## soil water constraint function:
        fswc = ifelse(swc < swc_min, 0.000001, ifelse(swc >swc_max, 1, (swc - swc_min) / (swc_max - swc_min)))
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
        
        Dm = vpd + (delta * aa - (delta + gamma)*et_total)*raa/(rho*cp)
      # transpiration 
        tc = PM(AA= aa-asubs, VPD=Dm, delta, ra=rbc, rs=rsc)
        
      # soil evaporation 
        es = PM(AA= asubs, VPD=Dm, delta, ra=ras, rs=rss)               
        
      # unit scale for ET 
        scale_day  =  day_length * 3600 / 2.45 * 10^-6 # transform to unit mm/d

      ## return results
      et_pre <- switch(out,
                   "t_tree"  = 0,
                   "t_grass" = tc * scale_day,
                   "e_soil"  = es * scale_day,
                   "et_total" = et_total * scale_day) # mm/d
  
    } 
           
    return(et_pre)
}




