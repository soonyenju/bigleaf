# %%R
# ==============================================================================
# WUE_metrics.r
WUE.metrics <- function(data,GPP="GPP",NEE="NEE",LE="LE",VPD="VPD",Tair="Tair",
                        constants=bigleaf.constants()){

  check.input(data,list(GPP,NEE,LE,VPD,Tair))

  ET  <- LE.to.ET(LE,Tair)                 # kg H2O m-2 s-1
  GPP <- (GPP * constants$umol2mol * constants$Cmol) * constants$kg2g  # gC m-2 s-1
  NEE <- (NEE * constants$umol2mol * constants$Cmol) * constants$kg2g  # gC m-2 s-1

  WUE     <- median(GPP/ET,na.rm=TRUE)
  WUE_NEE <- median(abs(NEE)/ET,na.rm=TRUE)
  IWUE    <- median((GPP*VPD)/ET,na.rm=TRUE)
  uWUE    <- median((GPP*sqrt(VPD))/ET,na.rm=TRUE)

  return(c(WUE=WUE,WUE_NEE=WUE_NEE,IWUE=IWUE,uWUE=uWUE))
}
# ==============================================================================
# aerodynamic_conductance.r
aerodynamic.conductance <- function(data,Tair="Tair",pressure="pressure",wind="wind",ustar="ustar",H="H",
                                    zr,zh,d,z0m=NULL,Dl,N=2,fc=NULL,LAI,Cd=0.2,hs=0.01,wind_profile=FALSE,
                                    stab_correction=TRUE,stab_formulation=c("Dyer_1970","Businger_1971"),
                                    Rb_model=c("Thom_1972","Choudhury_1988","Su_2001","constant_kB-1"),
                                    kB_h=NULL,Sc=NULL,Sc_name=NULL,constants=bigleaf.constants()){

  pv <- packageVersion("bigleaf")
  if (pv > "0.7.5"){
    cat("Note new column 'z0h' in the function output for 'bigleaf' versions > 0.7.5.",fill=TRUE)
  }

  Rb_model         <- match.arg(Rb_model)
  stab_formulation <- match.arg(stab_formulation)

  check.input(data,list(Tair,pressure,wind,ustar,H))

  ## calculate canopy boundary layer conductance (Gb)
  if (Rb_model %in% c("Thom_1972","Choudhury_1988","Su_2001")){

    if (Rb_model == "Thom_1972"){

      Gb_mod <- Gb.Thom(ustar=ustar,Sc=Sc,Sc_name=Sc_name,constants=constants)

    } else if (Rb_model == "Choudhury_1988"){

      Gb_mod <- Gb.Choudhury(data,Tair=Tair,pressure=pressure,wind=wind,ustar=ustar,
                             H=H,leafwidth=Dl,LAI=LAI,zh=zh,zr=zr,d=d,z0m=z0m,
                             stab_formulation=stab_formulation,Sc=Sc,Sc_name=Sc_name,
                             constants=constants)

    } else if (Rb_model == "Su_2001"){

      Gb_mod <- Gb.Su(data=data,Tair=Tair,pressure=pressure,ustar=ustar,wind=wind,
                      H=H,zh=zh,zr=zr,d=d,z0m=z0m,Dl=Dl,N=N,fc=fc,LAI=LAI,Cd=Cd,hs=hs,
                      stab_formulation=stab_formulation,Sc=Sc,Sc_name=Sc_name,
                      constants=constants)

    }

    kB_h <- Gb_mod[,"kB_h"]
    Rb_h <- Gb_mod[,"Rb_h"]
    Gb_h <- Gb_mod[,"Gb_h"]
    Gb_x <- data.frame(Gb_mod[,grep(colnames(Gb_mod),pattern="Gb_")[-1]])
    colnames(Gb_x) <- grep(colnames(Gb_mod),pattern="Gb_",value=TRUE)[-1]


  } else if (Rb_model == "constant_kB-1"){

    if(is.null(kB_h)){
      stop("value of kB-1 has to be specified if Rb_model is set to 'constant_kB-1'!")
    } else {
      Rb_h <- kB_h/(constants$k * ustar)
      Gb_h <- 1/Rb_h

      if (!is.null(Sc) | !is.null(Sc_name)){
        if (length(Sc) != length(Sc_name)){
          stop("arguments 'Sc' and 'Sc_name' must have the same length")
        }
        if (!is.numeric(Sc)){
          stop("argument 'Sc' must be numeric")
        }
      }

      Sc   <- c(constants$Sc_CO2,Sc)
      Gb_x <- data.frame(lapply(Sc,function(x) Gb_h / (x/constants$Pr)^0.67))
      colnames(Gb_x) <- paste0("Gb_",c("CO2",Sc_name))

    }

  }

  ## calculate aerodynamic conductance for momentum (Ga_m)
  if (wind_profile){

    if (is.null(z0m) & Rb_model %in% c("constant_kB-1","Thom_1972")){
      stop("z0m must be provided if wind_profile=TRUE!")
    } else if (is.null(z0m) & Rb_model %in% c("Choudhury_1988","Su_2001")){
      # z0m estimated as in Choudhury_1988 or Su_2001
      z0m <- roughness.parameters(method="wind_profile",zh=zh,zr=zr,d=d,data=data,
                                  Tair=Tair,pressure=pressure,wind=wind,ustar=ustar,H=H,
                                  stab_roughness=TRUE,stab_formulation=stab_formulation,
                                  constants=constants)[,"z0m"]
    }

    if (stab_correction){

      zeta  <-  stability.parameter(data=data,Tair=Tair,pressure=pressure,ustar=ustar,
                                    H=H,zr=zr,d=d,constants=constants)

      if (stab_formulation %in% c("Dyer_1970","Businger_1971")){

        psi_h <- stability.correction(zeta,formulation=stab_formulation)[,"psi_h"]

      } else {
        stop("'stab_formulation' has to be one of 'Dyer_1970' or 'Businger_1971'.
             Choose 'stab_correction = FALSE' if no stability correction should be applied.")
      }

      Ra_m  <- pmax((log((zr - d)/z0m) - psi_h),0) / (constants$k*ustar)

      } else {

        Ra_m  <- pmax((log((zr - d)/z0m)),0) / (constants$k*ustar)
        zeta = psi_h <- rep(NA_integer_,length=length(Ra_m))

      }

  } else {

    if ((!missing(zr) | !missing(d) | !missing(z0m)) & Rb_model %in% c("constant_kB-1","Thom_1972")){
      warning("Provided roughness length parameters (zr,d,z0m) are not used if 'wind_profile = FALSE' (the default). Ra_m is calculated as wind / ustar^2")
    }

    Ra_m <- wind / ustar^2
    zeta = psi_h <- rep(NA_integer_,length=length(Ra_m))

  }

  Ga_m   <- 1/Ra_m
  Ra_h   <- Ra_m + Rb_h
  Ga_h   <- 1/Ra_h
  Ga_x   <- 1/(Ra_m + 1/Gb_x)
  Ra_CO2 <- 1/Ga_x[,1]
  colnames(Ga_x) <- paste0("Ga_",c("CO2",Sc_name))

  if(!is.null(z0m)){
    z0h <- roughness.length.heat(z0m,kB_h)
  } else {
    z0h <- rep(NA_integer_,length=length(Ra_m))
  }

  Gab_x <- cbind(Ga_x,Gb_x)
  Gab_x <- Gab_x[rep(c(1,ncol(Gab_x)-(ncol(Gab_x)/2-1)),ncol(Gab_x)/2) + sort(rep(0:(ncol(Gab_x)/2-1),2))] # reorder columns

  return(data.frame(Ga_m,Ra_m,Ga_h,Ra_h,Gb_h,Rb_h,kB_h,z0h,zeta,psi_h,Ra_CO2,Gab_x))

}
# ==============================================================================
# bigleaf_constants.r
bigleaf.constants <- function(
  ## Physical constants
  cp         = 1004.834,        # specific heat of air for constant pressure (J K-1 kg-1)
  Rgas       = 8.31451,         # universal gas constant (J mol-1 K-1)
  Rv         = 461.5,           # gas constant of water vapor (J kg-1 K-1) (Stull 1988 p.641)
  Rd         = 287.0586,        # gas constant of dry air (J kg-1 K-1) (Foken 2008 p. 245)
  Md         = 0.0289645,       # molar mass of dry air (kg mol-1)
  Mw         = 0.0180153,       # molar mass of water vapor (kg mol-1)
  eps        = 0.622,           # ratio of the molecular weight of water vapor to dry air (=Mw/Md)
  g          = 9.81,            # gravitational acceleration (m s-2)
  solar_constant = 1366.1,      # solar constant, i.e. solar radation at earth distance from the sun (W m-2)
  pressure0  = 101325,          # reference atmospheric pressure at sea level (Pa)
  Tair0      = 273.15,          # reference air temperature (K)
  k          = 0.41,            # von Karman constant
  Cmol       = 0.012011,        # molar mass of carbon (kg mol-1)
  Omol       = 0.0159994,       # molar mass of oxygen (kg mol-1)
  H2Omol     = 0.01801528,      # molar mass of water (kg mol-1)
  sigma      = 5.670367e-08,    # Stefan-Boltzmann constant (W m-2 K-4)
  Pr         = 0.71,            # Prandtl number
  Sc_CO2     = 1.07,            # Schmidt number for CO2 (Hicks et al. 1987)
  Le067      = 0.93,            # Lewis number for water vapor to the power of 0.67

  ## Conversion constants
  Kelvin       = 273.15,         # conversion degree Celsius to Kelvin
  DwDc         = 1.6,            # Ratio of the molecular diffusivities for water vapor and CO2
  days2seconds = 86400,          # seconds per day
  kPa2Pa       = 1000,           # conversion kilopascal (kPa) to pascal (Pa)
  Pa2kPa       = 0.001,          # conversion pascal (Pa) to kilopascal (kPa)
  umol2mol     = 1e-06,          # conversion micromole (umol) to mole (mol)
  mol2umol     = 1e06,           # conversion mole (mol) to micromole (umol)
  kg2g         = 1000,           # conversion kilogram (kg) to gram (g)
  g2kg         = 0.001,          # conversion gram (g) to kilogram (kg)
  kJ2J         = 1000,           # conversion kilojoule (kJ) to joule (J)
  J2kJ         = 0.001,          # conversion joule (J) to kilojoule (kJ)
  se_median    = 1.253,          # conversion standard error (SE) of the mean to SE of the median (http://influentialpoints.com/Training/standard_error_of_median.htm)
  frac2percent = 100             # conversion between fraction and percent
){

  list(
    cp = cp, Rgas = Rgas, Rv = Rv, Rd = Rd, Md = Md, Mw = Mw, eps = eps, g = g,
    solar_constant = solar_constant,
    pressure0 = pressure0, Tair0 = Tair0, k = k, Cmol = Cmol, Omol = Omol,
    H2Omol = H2Omol,
    sigma = sigma, Pr = Pr, Sc_CO2 = Sc_CO2, Le067 = Le067, Kelvin = Kelvin, DwDc = DwDc,
    days2seconds = days2seconds, kPa2Pa = kPa2Pa, Pa2kPa = Pa2kPa, umol2mol = umol2mol,
    mol2umol = mol2umol, kg2g = kg2g, g2kg = g2kg, kJ2J = kJ2J, J2kJ = J2kJ,
    se_median = se_median, frac2percent = frac2percent
  )

}
# ==============================================================================
# bigleaf_physiology.r
intercellular.CO2 <- function(data,Ca="Ca",GPP="GPP",Gs="Gs_mol",Rleaf=NULL,
                              missing.Rleaf.as.NA=FALSE,constants=bigleaf.constants()){

  check.input(data,list(Ca,GPP,Gs))

  if(!is.null(Rleaf)){
    if(!missing.Rleaf.as.NA){Rleaf[is.na(Rleaf)] <- 0 }
  } else {
    cat("Respiration from the leaves is ignored and set to 0.",fill=TRUE)
    Rleaf <- 0
  }

  Ci <- Ca - (GPP - Rleaf)/(Gs/constants$DwDc)

  return(Ci)

}
photosynthetic.capacity <- function(data,C3=TRUE,Temp,GPP="GPP",Ci,PPFD="PPFD",PPFD_j=c(200,500),PPFD_c=1000,
                                    Rleaf=NULL,Oi=0.21,Kc25=404.9,Ko25=278.4,Gam25=42.75,
                                    Kc_Ha=79.43,Ko_Ha=36.38,Gam_Ha=37.83,Vcmax_Ha=65.33,Vcmax_Hd=200,
                                    Vcmax_dS=0.635,Jmax_Ha=43.9,Jmax_Hd=200,Jmax_dS=0.640,
                                    Theta=0.7,alpha_canopy=0.8,missing.Rleaf.as.NA=FALSE,Ci_C4=100,
                                    constants=bigleaf.constants()){

  check.input(data,list(Temp,GPP,Ci,PPFD))

  Temp <- Temp + constants$Kelvin
  Tref <- 25.0 + constants$Kelvin

  if (C3){  # C3 vegetation
    Kc_Ha    <- Kc_Ha * constants$kJ2J
    Ko_Ha    <- Ko_Ha * constants$kJ2J
    Gam_Ha   <- Gam_Ha * constants$kJ2J

    # Temperature dependencies of photosynthetic parameters
    Kc  <- Kc25 * exp(Kc_Ha * (Temp - Tref) / (Tref*constants$Rgas*Temp))
    Ko  <- Ko25 * exp(Ko_Ha * (Temp - Tref) / (Tref*constants$Rgas*Temp))
    Gam <- Gam25 * exp(Gam_Ha * (Temp - Tref) / (Tref*constants$Rgas*Temp))
    Ko  <- Ko * constants$J2kJ

    # basic filtering on Ci
    Ci[Ci < 80 | is.na(Ci)] <- NA

    # Presumed limitation states
    GPPc = GPPj <- GPP
    GPPj[PPFD < PPFD_j[1] | PPFD > PPFD_j[2] | is.na(PPFD)] <- NA
    GPPc[PPFD < PPFD_c | is.na(PPFD)] <- NA

    if(!is.null(Rleaf)){
      if(!missing.Rleaf.as.NA){Rleaf[is.na(Rleaf)] <- 0 }
    } else {
      cat("Respiration from the leaves is ignored and set to 0.",fill=TRUE)
      Rleaf <- 0
    }

    # calculate Vcmax and J (electron transport rate)
    Vcmax <- (GPPc-Rleaf) * (Ci + Kc*(1.0 + Oi/Ko)) / (Ci - Gam)
    J     <- (GPPj-Rleaf) * (4.0 * Ci + 8.0 * Gam) / (Ci - Gam)


  } else {  # C4 vegetation

    # Presumed limitation states (C4)
    GPPc = GPPj <- GPP
    GPPj[PPFD < PPFD_j[1] | PPFD > PPFD_j[2] | is.na(PPFD) | Ci < 0] <- NA
    GPPc[PPFD < PPFD_c | Ci < Ci_C4 | is.na(PPFD)] <- NA

    Vcmax <- GPPc
    J     <- 3 * GPPj / (1 - 0.5)

  }


  # calculate Jmax from J
  APPFD_PSII <- PPFD * alpha_canopy * 0.85 * 0.5

  calcJmax  <- which(complete.cases(J,APPFD_PSII))
  if (length(calcJmax) > 0){
    Jmax <- sapply(calcJmax, function(i) tryCatch(optimize(function(Jmax){abs(J[i] - c((APPFD_PSII[i] + Jmax -
                                                                                          sqrt((APPFD_PSII[i] + Jmax)^2 -
                                                                                                 4.0 * Theta * APPFD_PSII[i] * Jmax)) /
                                                                                         (2.0 * Theta)))},
                                                           interval=c(0,1000),tol=1e-02)$minimum,
                                                  error=function(err){NA}
    )
    )
  } else {
    warning("Not enough observations to calculate Jmax!")
    Jmax <- NA
  }

  # calculate Vcmax25 and Jmax25
  Vcmax25 <- Arrhenius.temp.response(Vcmax,Temp-constants$Kelvin,Ha=Vcmax_Ha,
                                     Hd=Vcmax_Hd,dS=Vcmax_dS,constants=constants)

  Jmax25 <- Arrhenius.temp.response(Jmax,Temp[calcJmax]-constants$Kelvin,Ha=Jmax_Ha,
                                    Hd=Jmax_Hd,dS=Jmax_dS,constants=constants)


  # calculate medians and standard errors of the median
  Vcmax25_Median <- median(Vcmax25,na.rm=TRUE)
  Vcmax25_SE     <- constants$se_median * sd(Vcmax25,na.rm=TRUE)/sqrt((sum(!is.na(Vcmax25))))
  Jmax25_Median  <- median(Jmax25,na.rm=TRUE)
  Jmax25_SE      <- constants$se_median * sd(Jmax25,na.rm=TRUE)/sqrt((sum(!is.na(Jmax25))))

  return(c("Vcmax25"=round(Vcmax25_Median,2),"Vcmax25_SE"=round(Vcmax25_SE,2),
           "Jmax25"=round(Jmax25_Median,2),"Jmax25_SE"=round(Jmax25_SE,2)))

}

Arrhenius.temp.response <- function(param,Temp,Ha,Hd,dS,constants=bigleaf.constants()){

  Temp <- Temp + constants$Kelvin
  Tref <- 25.0 + constants$Kelvin

  Ha <- ifelse(missing(Ha),NA,Ha*constants$kJ2J)
  Hd <- ifelse(missing(Hd),NA,Hd*constants$kJ2J)
  dS <- ifelse(missing(dS),NA,dS*constants$kJ2J)

  if (is.na(Ha)){

    stop("Activation energy (Ha) has to be provided!")

  }

  if (is.na(Hd) & is.na(dS)){

    param25 <- param / exp(Ha * (Temp - Tref) / (Tref*constants$Rgas*Temp))

  } else if (!is.na(Hd) & !is.na(dS)){

    param25 <- param /
      ( exp(Ha * (Temp - Tref) / (Tref*constants$Rgas*Temp)) *
        (1 + exp((Tref*dS - Hd) / (Tref * constants$Rgas))) /
        (1 + exp((Temp*dS - Hd) / (Temp * constants$Rgas)))
      )

  } else if ((!is.na(Hd) & is.na(dS)) | (is.na(Hd) & !is.na(dS)) ){

    warning("Both Hd and dS have to be provided for a temperature response
             that considers a temperature optimum and a deactivation term!
             Continue considering activation energy (Ha) only...")

    param25 <- param / exp(Ha * (Temp - Tref) / (Tref*constants$Rgas*Temp))

  }

  return(param25)

}

stomatal.slope <- function(data,Tair="Tair",pressure="pressure",GPP="GPP",Gs="Gs_mol",
                           VPD="VPD",Ca="Ca",Rleaf=NULL,model=c("USO","Ball&Berry","Leuning"),
                           robust.nls=FALSE,nmin=40,fitg0=FALSE,g0=0,fitD0=FALSE,
                           D0=1.5,Gamma=50,missing.Rleaf.as.NA=FALSE,
                           constants=bigleaf.constants(),...){

  model <- match.arg(model)

  check.input(data,list(Tair,pressure,GPP,Gs,VPD,Ca))

  df   <- data.frame(Tair,pressure,GPP,Gs,VPD,Ca)
  DwDc <- constants$DwDc  # ...to work within nls()


  if (model == "Leuning"){
    check.input(data,Gamma)
    df$Gamma <- Gamma
  }



  if(!is.null(Rleaf)){
    if(!missing.Rleaf.as.NA){Rleaf[is.na(Rleaf)] <- 0 }
  } else {
    cat("Respiration from the leaves is ignored and set to 0.",fill=TRUE)
    Rleaf <- 0
  }

  GPP <- (GPP - Rleaf)


  if (model == "Leuning"){
    nr_data <- sum(!is.na(GPP) & !is.na(Gs) & !is.na(VPD) & !is.na(Ca) & !is.na(Gamma))
  } else {
    nr_data <- sum(!is.na(GPP) & !is.na(Gs) & !is.na(VPD) & !is.na(Ca))
  }


  if (nr_data < nmin){
    stop("number of data is less than 'nmin'. g1 is not fitted to the data.")
  } else {

    if (model == "USO"){

      if (fitg0){
        if (robust.nls){
          df$DwDc <- rep(DwDc,nrow(df))
          mod_weights <- nlrob(Gs ~ g0 + DwDc*(1.0 + g1/sqrt(VPD))*GPP/Ca,data=df,start=list(g0=0,g1=3),
                               na.action=na.exclude,...)$w
          mod <- nls(Gs ~ g0 + DwDc*(1.0 + g1/sqrt(VPD))*GPP/Ca,start=list(g0=0,g1=3),weights=mod_weights,...)
        } else {
          mod <- nls(Gs ~ g0 + DwDc*(1.0 + g1/sqrt(VPD))*GPP/Ca,start=list(g0=0,g1=3),...)
        }
      } else {
        if (robust.nls){
          df$g0   <- rep(g0,nrow(df))    # g0 as constant does not work in the nlrob function...
          df$DwDc <- rep(DwDc,nrow(df))  # same with constants$DwDc
          mod_weights <- nlrob(Gs ~ g0 + DwDc*(1.0 + g1/sqrt(VPD))*GPP/Ca,data=df,start=list(g1=3),
                               na.action=na.exclude,...)$w
          mod <- nls(Gs ~ g0 + DwDc*(1.0 + g1/sqrt(VPD))*GPP/Ca,start=list(g1=3),weights=mod_weights,...)
        } else {
          mod <- nls(Gs ~ g0 + DwDc*(1.0 + g1/sqrt(VPD))*GPP/Ca,start=list(g1=3),...)
        }
      }

    } else if (model == "Leuning"){

      if (fitg0){
        if (fitD0){
          if (robust.nls){
            mod_weights <- nlrob(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),data=df,
                                 start=list(g0=0,g1=9,D0=1.5),na.action=na.exclude,...)$w
            mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g0=0,g1=9,D0=1.5),
                       weights=mod_weights,...)
          } else {
            mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g0=0,g1=9,D0=1.5),...)
          }
        } else {
          if (robust.nls){
            df$D0  <- rep(D0,nrow(df))
            mod_weights <- nlrob(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),data=df,
                                 start=list(g0=0,g1=9),na.action=na.exclude,...)$w
            mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g0=0,g1=9),
                       weights=mod_weights,...)
          } else {
            mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g0=0,g1=9),...)
          }
        }
      } else {
        if (fitD0){
          if (robust.nls){
            df$g0    <- rep(g0,nrow(df))
            mod_weights <- nlrob(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),data=df,
                                 start=list(g1=9,D0=1.5),na.action=na.exclude,...)$w
            mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g1=9,D0=1.5),
                       weights=mod_weights,...)
          } else {
            mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g1=9,D0=1.5),...)
          }
        } else {
          if (robust.nls){
            df$g0  <- rep(g0,nrow(df))
            df$D0  <- rep(D0,nrow(df))
            mod_weights <- nlrob(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),data=df,
                                 start=list(g1=9),na.action=na.exclude,...)$w
            mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g1=9),
                       weights=mod_weights,...)
          } else {
            mod <- nls(Gs ~ g0 + g1*GPP / ((Ca - Gamma) * (1 + VPD/D0)),start=list(g1=9),...)
          }
        }
      }

    } else if (model == "Ball&Berry"){

      rH <- VPD.to.rH(VPD,Tair)
      df$rH <- rH

      if (fitg0){
        if (robust.nls){
          mod_weights <- nlrob(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g0=0,g1=9),data=df,
                               na.action=na.exclude,...)$w
          mod <- nls(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g0=0,g1=9),weights=mod_weights,...)
        } else {
          mod <- nls(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g0=0,g1=9),...)
        }
      } else {
        if (robust.nls){
          df$g0   <- rep(g0,nrow(df))
          mod_weights <- nlrob(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g1=9),data=df,
                               na.action=na.exclude,...)$w
          mod <- nls(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g1=9),weights=mod_weights,...)
        } else {
          mod <- nls(Gs ~ g0 + g1 * (GPP * rH) / Ca,start=list(g1=9),...)
        }
      }

    }

  }

  return(mod)

}
light.response <- function(data,NEE="NEE",Reco="Reco",PPFD="PPFD",PPFD_ref=2000,...){

  check.input(data,list(NEE,Reco,PPFD))

  mod <- nls(-NEE ~ alpha * PPFD / (1 - (PPFD / PPFD_ref) + (alpha * PPFD / GPP_ref)) - Reco,
             start=list(alpha=0.05,GPP_ref=30),...)

  return(mod)
}

light.use.efficiency <- function(GPP,PPFD){

  comp <- complete.cases(GPP,PPFD)

  LUE <- sum(GPP[comp],na.rm=TRUE)/sum(PPFD[comp],na.rm=TRUE)

  return(c("LUE"=LUE))
}

stomatal.sensitivity <- function(data,Gs="Gs_mol",VPD="VPD",...){

  check.input(data,list(Gs,VPD))

  mod <- nls(Gs ~ -m * log(VPD) + b,start=list(m=0.05,b=0.2),...)

  return(mod)
}
# ==============================================================================
# boundary_layer_conductance.r
Gb.Thom <- function(ustar,Sc=NULL,Sc_name=NULL,constants=bigleaf.constants()){

  check.input(NULL,ustar)

  Rb_h <- 6.2*ustar^-0.667
  Gb_h <- 1/Rb_h
  kB_h <- Rb_h*constants$k*ustar

  if (!is.null(Sc) | !is.null(Sc_name)){
    if (length(Sc) != length(Sc_name)){
      stop("arguments 'Sc' and 'Sc_name' must have the same length")
    }
    if (!is.numeric(Sc)){
      stop("argument 'Sc' must be numeric")
    }
  }

  Sc   <- c(constants$Sc_CO2,Sc)
  Gb_x <- data.frame(lapply(Sc,function(x) Gb_h / (x/constants$Pr)^0.67))
  colnames(Gb_x) <- paste0("Gb_",c("CO2",Sc_name))

  return(data.frame(Gb_h,Rb_h,kB_h,Gb_x))
}

Gb.Choudhury <- function(data,Tair="Tair",pressure="pressure",wind="wind",ustar="ustar",H="H",
                         leafwidth,LAI,zh,zr,d,z0m=NULL,stab_formulation=c("Dyer_1970","Businger_1971"),
                         Sc=NULL,Sc_name=NULL,constants=bigleaf.constants()){

  stab_formulation <- match.arg(stab_formulation)

  check.input(data,list(Tair,pressure,wind,ustar,H))

  alpha   <- 4.39 - 3.97*exp(-0.258*LAI)

  if (is.null(z0m)){
    estimate_z0m <- TRUE
    z0m <- NULL
  } else {
    estimate_z0m <- FALSE
  }

  wind_zh <- wind.profile(data=data,z=zh,Tair=Tair,pressure=pressure,ustar=ustar,H=H,
                          zr=zr,estimate_z0m=estimate_z0m,zh=zh,d=d,z0m=z0m,frac_z0m=NULL,
                          stab_correction=TRUE,stab_formulation=stab_formulation)

  ## avoid zero windspeed
  wind_zh <- pmax(0.01,wind_zh)

  if (!is.null(Sc) | !is.null(Sc_name)){
    if (length(Sc) != length(Sc_name)){
      stop("arguments 'Sc' and 'Sc_name' must have the same length")
    }
    if (!is.numeric(Sc)){
      stop("argument 'Sc' must be numeric")
    }
  }

  Gb_h <- LAI*((0.02/alpha)*sqrt(wind_zh/leafwidth)*(1-exp(-alpha/2)))
  Rb_h <- 1/Gb_h
  kB_h <- Rb_h*constants$k*ustar

  Sc   <- c(constants$Sc_CO2,Sc)
  Gb_x <- data.frame(lapply(Sc,function(x) Gb_h / (x/constants$Pr)^0.67))
  colnames(Gb_x) <- paste0("Gb_",c("CO2",Sc_name))


  return(data.frame(Gb_h,Rb_h,kB_h,Gb_x))
}

Gb.Su <- function(data,Tair="Tair",pressure="pressure",ustar="ustar",wind="wind",
                  H="H",zh,zr,d,z0m=NULL,Dl,fc=NULL,LAI=NULL,N=2,Cd=0.2,hs=0.01,
                  stab_formulation=c("Dyer_1970","Businger_1971"),
                  Sc=NULL,Sc_name=NULL,constants=bigleaf.constants()){

  stab_formulation <- match.arg(stab_formulation)

  check.input(data,list(Tair,pressure,ustar,wind,H))

  if (is.null(fc)){
    if (is.null(LAI)){
      stop("one of 'fc' or 'LAI' must be provided",call.=FALSE)
    } else {
      fc <- (1-exp(-LAI/2))
    }
  }

  if (is.null(z0m)){
    estimate_z0m <- TRUE
    z0m <- NULL
  } else {
    estimate_z0m <- FALSE
  }

  wind_zh <- wind.profile(data=data,z=zh,Tair=Tair,pressure=pressure,ustar=ustar,H=H,
                          zr=zr,estimate_z0m=estimate_z0m,zh=zh,d=d,z0m=z0m,frac_z0m=NULL,
                          stab_correction=TRUE,stab_formulation=stab_formulation)

  v   <- kinematic.viscosity(Tair,pressure,constants)
  Re  <- Reynolds.Number(Tair,pressure,ustar,hs,constants)
  kBs <- 2.46 * (Re)^0.25 - log(7.4)
  Reh <- Dl * wind_zh / v
  Ct  <- 1*constants$Pr^-0.6667*Reh^-0.5*N

  kB_h <- (constants$k*Cd)/(4*Ct*ustar/wind_zh)*fc^2 + kBs*(1 - fc)^2
  Rb_h <- kB_h/(constants$k*ustar)
  Gb_h <- 1/Rb_h

  if (!is.null(Sc) | !is.null(Sc_name)){
    if (length(Sc) != length(Sc_name)){
      stop("arguments 'Sc' and 'Sc_name' must have the same length")
    }
    if (!is.numeric(Sc)){
      stop("argument 'Sc' must be numeric")
    }
  }

  Sc   <- c(constants$Sc_CO2,Sc)
  Gb_x <- data.frame(lapply(Sc,function(x) Gb_h / (x/constants$Pr)^0.67))
  colnames(Gb_x) <- paste0("Gb_",c("CO2",Sc_name))

  return(data.frame(Gb_h,Rb_h,kB_h,Gb_x))
}

roughness.length.heat <- function(z0m,kB_h){

  z0h <- z0m / exp(kB_h)

  return(z0h)
}
# ==============================================================================
# check_input.r
check.input <- function(data,...){

  vars <- check.length(list(...))

  if (missing(data)){
    data <- NULL
  }

  varlist  <- match.call()[-c(1:2)]
  varnames <- c(unlist(sapply(varlist,as.character)))
  varnames <- varnames[!varnames %in% c("c","list")]

  for (v in seq_along(vars)){

    var     <- vars[[v]]
    varname <- ifelse(varnames[v] %in% c("var","var_qc"),gsub("\"","",deparse(substitute(var))),varnames[v])

    if (is.character(var)){
      if (!missing(data) & !is.null(data)){
        if (length(var) == 1){
          if (var %in% colnames(data)){
            var <- data[,var]
            if (is.numeric(var)){
              assign(varname,var,pos=sys.frame(-1))
            } else {
              stop("column representing '",varname,"' in the input matrix/data.frame must be numeric",call.=FALSE)
            }
          } else {
            stop ("there is no column named '",var,"' in the input matrix/data.frame. Indicate the name of the column representing variable '",varname,"', or alternatively, provide a numeric vector of the same length as the input matrix/data.frame or of length 1.",call.=FALSE)
          }
        } else {
          stop("name of variable '",varname,"' must have length 1",call.=FALSE)
        }
      } else {
        if ("data" %in% names(formals(sys.function(which=-1)))){
          if (var %in% as.character(unlist(match.call(definition=sys.function(-1),call=sys.call(-1))[-1]))){
            stop("variable '",var,"' is of type character and interpreted as a column name, but no input matrix/data.frame is provided. Provide '",var,"' as a numeric vector, or an input matrix/data.frame with a column named '",var,"'",call.=FALSE)
          } else {
            stop("variable '",var,"' is not provided",call.=FALSE)
          }
        } else {
          stop("variable '",var,"' must be numeric",call.=FALSE)
        }
      }
    } else {
      if (length(var) < 2){
        if (is.null(var)){
          assign(varname,var,pos=sys.frame(-1))
          next
        } else if (is.na(var)){
          assign(varname,var,pos=sys.frame(-1))
          next
        }
      }
      if (!missing(data) & !is.null(data)){
        if (is.numeric(var) & length(var) == nrow(data)){
          assign(varname,var,envir=sys.frame(-1))
        } else if (is.numeric(var) & length(var) != nrow(data)) {
          if (length(var) == 1){
            var <- rep(var,length=nrow(data))
            assign(varname,var,envir=sys.frame(-1))
          } else {
            stop("variable '",varname,"' must have the same length as the input matrix/data.frame or length 1. Do NOT provide an input matrix/data.frame if none of its variables are used!",call.=FALSE)
          }
        } else if (!is.numeric(var)){
          stop("variable '",varname,"' must be numeric",call.=FALSE)
        }
      } else {
        if (is.numeric(var)){
          assign(varname,var,envir=sys.frame(-1))
        } else {
          stop("variable '",varname,"' must be numeric",call.=FALSE)
        }
      }
    }
  }
}
check.length <- function(varlist){

  if (is.list(unlist(varlist,recursive=FALSE))){
    varlist <- unlist(varlist,recursive=FALSE)
  }

  length.vars <- sapply(varlist,length)
  length.vars <- length.vars[length.vars > 0]

  if (length(unique(length.vars)) >= 2){
    if (sort(unique(length.vars))[1] != 1 | length(unique(length.vars)) > 2){
      stop("All input variables must have the same length or a length of 1!",call.=FALSE)
    }
  }
  return(varlist)
}
# ==============================================================================
# decoupling.r
decoupling <- function(data,Tair="Tair",pressure="pressure",Ga="Ga_h",Gs="Gs_ms",
                       approach=c("Jarvis&McNaughton_1986","Martin_1989"),
                       LAI,Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                       constants=bigleaf.constants()){

  approach    <- match.arg(approach)

  check.input(data,list(Tair,pressure,Ga,Gs))

  Delta   <- Esat.slope(Tair,Esat.formula,constants)[,"Delta"]
  gamma   <- psychrometric.constant(Tair,pressure,constants)
  epsilon <- Delta/gamma

  if (approach == "Jarvis&McNaughton_1986"){

    Omega <- (epsilon + 1) / (epsilon + 1 + Ga/Gs)

  } else if (approach == "Martin_1989") {

    if (is.null(LAI)){

      stop("LAI is not provided!")

    } else {

      Gr    <- longwave.conductance(Tair,LAI,constants)
      Omega <- (epsilon + 1 + Gr/Ga) / (epsilon + 1 + Ga/Gs + Gr/Gs + Gr/Ga)

    }
  }

  return(Omega)

}
longwave.conductance <- function(Tair,LAI,constants=bigleaf.constants()){

  Tair <- Tair + constants$Kelvin

  Gr <- 4 * constants$sigma * Tair^3 * LAI / constants$cp

  return(Gr)
}
# ==============================================================================
# energy_balance.r
biochemical.energy <- function(NEE,alpha=0.422){
  Sp <- alpha*-NEE
  return(Sp)
}

energy.use.efficiency <- function(GPP,alpha=0.422,Rn){

  Sp <- biochemical.energy(-GPP,alpha)

  comp  <- complete.cases(Sp,Rn)

  Sp_sum <- sum(Sp[comp],na.rm=T)
  Rn_sum <- sum(Rn[comp],na.rm=T)

  EUE <- Sp_sum/Rn_sum

  return(c("EUE"=EUE))
}

energy.closure <- function(data,Rn="Rn",G=NULL,S=NULL,LE="LE",H="H",instantaneous=FALSE,
                           missing.G.as.NA=FALSE,missing.S.as.NA=FALSE){

  check.input(data,list(Rn,LE,H,G,S))

  if(!is.null(G)){
    if (!missing.G.as.NA){G[is.na(G)] <- 0}
  } else {
    cat("Ground heat flux G is not provided and set to 0.",fill=TRUE)
    G <- rep(0,nrow(data))
  }

  if(!is.null(S)){
    if(!missing.S.as.NA){S[is.na(S)] <- 0 }
  } else {
    cat("Energy storage fluxes S are not provided and set to 0.",fill=TRUE)
    S <- rep(0,nrow(data))
  }

  if (!instantaneous){
    comp <- complete.cases(Rn,G,S,LE,H)
    n    <- sum(comp)

    EBR <- sum(LE[comp] + H[comp]) / sum(Rn[comp] - G[comp] - S[comp])

    emod <- lm(c(LE + H) ~ c(Rn - G - S))
    intercept <- summary(emod)$coef[1,1]
    slope     <- summary(emod)$coef[2,1]
    r_squared <- summary(emod)$r.squared

    return(c("n"=n,"intercept"=round(intercept,3),"slope"=round(slope,3),"r^2"=round(r_squared,3),"EBR"=round(EBR,3)))

  } else {

    EBR <- (LE + H) /(Rn - G - S)

    return(EBR)
  }
}

isothermal.Rn <- function(data,Rn="Rn",Tair="Tair",Tsurf="Tsurf",emissivity,
                          constants=bigleaf.constants()){

  check.input(data,list(Rn,Tair,Tsurf))

  Tair  <- Tair + constants$Kelvin
  Tsurf <- Tsurf + constants$Kelvin

  Rni <- Rn + emissivity * constants$sigma * (Tsurf^4 - Tair^4)

  return(Rni)

}
# ==============================================================================
# evapotranspiration.r
potential.ET <- function(data,Tair="Tair",pressure="pressure",Rn="Rn",G=NULL,S=NULL,
                         VPD="VPD",Ga="Ga_h",approach=c("Priestley-Taylor","Penman-Monteith"),
                         alpha=1.26,Gs_pot=0.6,missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,
                         Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                         constants=bigleaf.constants()){

  approach <- match.arg(approach)

  check.input(data,list(Tair,pressure,Rn,G,S))

  if(!is.null(G)){
    if (!missing.G.as.NA){G[is.na(G)] <- 0}
  } else {
    cat("Ground heat flux G is not provided and set to 0.",fill=TRUE)
    G <- 0
  }

  if(!is.null(S)){
    if(!missing.S.as.NA){S[is.na(S)] <- 0 }
  } else {
    cat("Energy storage fluxes S are not provided and set to 0.",fill=TRUE)
    S <- 0
  }

  gamma  <- psychrometric.constant(Tair,pressure,constants)
  Delta  <- Esat.slope(Tair,Esat.formula,constants)[,"Delta"]


  if (approach == "Priestley-Taylor"){

    LE_pot <- (alpha * Delta * (Rn - G - S)) / (Delta + gamma)
    ET_pot <- LE.to.ET(LE_pot,Tair)

  } else if (approach == "Penman-Monteith"){

    check.input(data,list(Gs_pot,VPD,Ga))

    Gs_pot <- mol.to.ms(Gs_pot,Tair=Tair,pressure=pressure,constants=constants)
    rho    <- air.density(Tair,pressure,constants)

    LE_pot <- (Delta * (Rn - G - S) + rho * constants$cp * VPD * Ga) /
      (Delta + gamma * (1 + Ga / Gs_pot))
    ET_pot <- LE.to.ET(LE_pot,Tair)
  }

  return(data.frame(ET_pot,LE_pot))

}

reference.ET <- function(data,Gs_ref=0.0143,Tair="Tair",pressure="pressure",VPD="VPD",Rn="Rn",Ga="Ga_h",
                         G=NULL,S=NULL,missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,
                         Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                         constants=bigleaf.constants()){

  stop("this function is deprecated (since bigleaf version 0.6.0). For the calculation of potential ET from the Penman-Monteith equation (as formerly calculated with this function), use function potential.ET() with the argument approach='Penman-Monteith'. Note that the default value for argument 'Gs_pot' is expressed now in mol m-2 s-1 for simplicity (0.6 mol m-2 s-1).")

}

equilibrium.imposed.ET <- function(data,Tair="Tair",pressure="pressure",VPD="VPD",Gs="Gs_ms",
                                   Rn="Rn",G=NULL,S=NULL,missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,
                                   Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                                   constants=bigleaf.constants()){

  check.input(data,list(Tair,pressure,VPD,Rn,Gs,G,S))

  if(!is.null(G)){
    if (!missing.G.as.NA){G[is.na(G)] <- 0}
  } else {
    cat("Ground heat flux G is not provided and set to 0.",fill=TRUE)
    G <- 0
  }

  if(!is.null(S)){
    if(!missing.S.as.NA){S[is.na(S)] <- 0 }
  } else {
    cat("Energy storage fluxes S are not provided and set to 0.",fill=TRUE)
    S <- 0
  }

  rho    <- air.density(Tair,pressure,constants)
  gamma  <- psychrometric.constant(Tair,pressure,constants)
  Delta  <- Esat.slope(Tair,Esat.formula,constants)[,"Delta"]

  LE_eq  <- (Delta * (Rn - G - S)) / (gamma + Delta)
  LE_imp <- (rho * constants$cp * Gs * VPD) / gamma

  ET_imp <- LE.to.ET(LE_imp,Tair)
  ET_eq  <- LE.to.ET(LE_eq,Tair)

  return(data.frame(ET_eq,ET_imp,LE_eq,LE_imp))
}
# ==============================================================================
# filter_data.r
filter.data <- function(data,quality.control=TRUE,filter.growseas=FALSE,
                        filter.precip=FALSE,filter.vars=NULL,
                        filter.vals.min,filter.vals.max,NA.as.invalid=TRUE,
                        vars.qc=NULL,quality.ext="_qc",good.quality=c(0,1),
                        missing.qc.as.bad=TRUE,GPP="GPP",doy="doy",
                        year="year",tGPP=0.5,ws=15,min.int=5,precip="precip",
                        tprecip=0.01,precip.hours=24,records.per.hour=2,
                        filtered.data.to.NA=TRUE,constants=bigleaf.constants()){


  ### I) Quality control filter
  if (quality.control){

    if (is.null(vars.qc)){
      stop("quality.control (qc) is TRUE, but no qc variables are provided!")
    }

    if (any(!vars.qc %in% colnames(data))){

      missing_vars <- vars.qc[which(!vars.qc %in% colnames(data))]
      stop(paste("Variable ",missing_vars," is included in 'vars.qc', but does not exist in the input data!"))

    }

    vars.qc_qc <- paste0(vars.qc,quality.ext)
    if (any(!vars.qc_qc %in% colnames(data))){

      missing_vars_qc <- vars.qc_qc[which(!vars.qc_qc %in% colnames(data))]
      missing_vars2   <- substr(missing_vars_qc,1,nchar(missing_vars_qc) - nchar(quality.ext))
      stop(paste("Quality control for variable ",missing_vars2,"(",missing_vars_qc,") does not exist in the input data!"))
    }

    ## data quality
    cat("Quality control:",fill=TRUE)
    for (var in vars.qc){
      var_qc <- paste0(var,quality.ext)
      check.input(data,var)
      check.input(data,var_qc)

      if (missing.qc.as.bad){
        data[get(paste0(var,quality.ext)) > max(good.quality) | is.na(get(paste0(var,quality.ext))),var] <- NA   # exclude bad quality data or those where qc flag is not available
        qc_invalid      <- sum(get(paste0(var,quality.ext)) > max(good.quality) | is.na(get(paste0(var,quality.ext)))) # count & report
      } else { # same, but consider missing quality flag variables as good
        data[get(paste0(var,quality.ext)) > max(good.quality) & !is.na(get(paste0(var,quality.ext))),var] <- NA
        qc_invalid      <- sum(get(paste0(var,quality.ext)) > max(good.quality) & !is.na(get(paste0(var,quality.ext))))
      }

      qc_invalid_perc <- round((qc_invalid/nrow(data))*constants$frac2percent,2)

      cat(var,": ",qc_invalid," data points (",qc_invalid_perc,"%) set to NA",fill=TRUE,sep="")
    }
  }


  #### II) Data filter
  valid <- rep(1L,nrow(data))

  # 1) GPP
  growseas_invalid <- numeric()
  if(filter.growseas){
    check.input(data,doy,year,GPP)
    date             <- strptime(paste0(year,"-",doy),format="%Y-%j")
    GPP_daily        <- aggregate(GPP,by=list(strftime(date)),function(x){if (sum(!is.na(x)) < 3){NA}else{mean(x,na.rm=TRUE)}})[,2]
    growing_season   <- filter.growing.season(GPP_daily,tGPP=tGPP,ws=ws,min.int=min.int)
    tsperday         <- table(as.character(date))
    growseas_invalid <- which(rep(growing_season,tsperday) == 0)
  }

  # 2) precipitation
  precip_invalid <- numeric()
  if (filter.precip){
    check.input(data,precip)
    if (NA.as.invalid){
      precip_events <- which(precip > tprecip | is.na(precip))
    } else {
      precip_events <- which(precip > tprecip)
    }
    precip_invalid <- unique(as.numeric(unlist(sapply(precip_events, function(x) x:(min(x+precip.hours*records.per.hour,nrow(data),na.rm=TRUE))))))
  }

  # 3) all other filter variables (as defined in filter.vars)
  invalids <- list(growseas_invalid,precip_invalid)

  if (!is.null(filter.vars)){
    for (var in filter.vars){
      v  <- which(filter.vars == var)
      vf <- v + 2
      check.input(data,var)
      if (NA.as.invalid){
        invalids[[vf]] <- which(get(var) < filter.vals.min[v] | get(var) > filter.vals.max[v] | is.na(get(var)))
      } else {
        invalids[[vf]] <- which(get(var) < filter.vals.min[v] | get(var) > filter.vals.max[v] & !is.na(get(var)))
      }
    }
  }

  # 4) calculate number and percentage of filtered values
  invalids_perc <- sapply(invalids, function(x) round((length(x)/nrow(data))*constants$frac2percent,2))

  additional_invalids <- sapply(2:length(invalids), function(x)
    length(setdiff(invalids[[x]],unique(unlist(invalids[1:(x-1)])))))

  additional_invalids_perc <- round(additional_invalids/nrow(data)*constants$frac2percent,2)


  # 5) write to output
  if (filter.growseas | filter.precip | length(filter.vars) > 0){

    var.names <- c("growing season","precipitation",filter.vars)

    if (quality.control){
      cat("-------------------------------------------------------------------",fill=TRUE)
    }

    cat("Data filtering:",fill=TRUE)

    cat(length(growseas_invalid)," data points (",invalids_perc[1],"%) excluded by growing season filter",fill=TRUE,sep="")

    invisible(sapply(c(1:(length(invalids)-1)), function(x) cat(additional_invalids[x]," additional data points (",
                                                                additional_invalids_perc[x],"%) excluded by ",var.names[x+1],
                                                                " filter (",length(unlist(invalids[x+1]))," data points = ",
                                                                invalids_perc[x+1]," % in total)",fill=TRUE,sep="")))


    invalid        <- unique(unlist(invalids))
    valid[invalid] <- 0

    excl_perc <- round((length(invalid)/nrow(data))*constants$frac2percent,2)

    cat(length(invalid)," data points (",excl_perc,"%) excluded in total",fill=TRUE,sep="")
    cat(nrow(data) - length(invalid)," valid data points (",constants$frac2percent-excl_perc,"%) remaining.",fill=TRUE,sep="")


    # 6) return input data frame with filtered time steps set to NA or an additional 'valid' column
    if (filtered.data.to.NA){
      data_filtered <- data
      data_filtered[valid < 1,] <- NA
    } else {
      data_filtered <- data.frame(data,valid)
    }

  } else {

    data_filtered <- data

  }

  return(data_filtered)
}


filter.growing.season <- function(GPPd,tGPP,ws=15,min.int=5){

  if(sum(is.na(GPPd)) < 0.5*length(GPPd)){

    growseas      <- rep(1,length(GPPd))
    GPP_threshold <- quantile(GPPd,probs=0.95,na.rm=TRUE)*tGPP

    ## smooth GPP
    GPPd_smoothed <- filter(GPPd,method="convolution",filter=rep(1/ws,ws))

    ## set values at the beginning and end of the time series to the mean of the original values
    wsd <- floor(ws/2)
    GPPd_smoothed[1:wsd] <- mean(GPPd[1:(2*wsd)],na.rm=TRUE)
    GPPd_smoothed[(length(GPPd)-(wsd-1)):length(GPPd)] <- mean(GPPd[(length(GPPd)-(2*wsd-1)):length(GPPd)],na.rm=TRUE)

    # check for occurrence of missing values and set them to mean of the values surrounding them
    missing <- which(is.na(GPPd_smoothed))
    if (length(missing) > 0){
      if (length(missing) > 10){warning("Attention, there is a gap in 'GPPd' of length n = ",length(missing))}
      replace_val <- mean(GPPd_smoothed[max(1,missing[1] - 4):min((missing[length(missing)] + 4),length(GPPd_smoothed))],na.rm=TRUE)
      GPPd_smoothed[missing] <- replace_val
    }

    # filter daily GPP
    growseas[GPPd_smoothed < GPP_threshold] <- 0

    ## change short intervals to the surrounding values to avoid 'wrong' fluctuations
    intervals <- rle(growseas)
    short_int <- which(intervals$lengths <= min.int)

    if (length(short_int) > 0){
      start <- numeric()
      end   <- numeric()

      for (i in 1:length(short_int)){

        start[i] <- sum(intervals$lengths[1:short_int[i]-1]) + 1
        end[i]   <- start[i]+intervals$lengths[short_int[i]] - 1

        val <- unique(growseas[start[i]:end[i]])

        if (val == 0 & growseas[start[i]-1] == 1){
          growseas[start[i]:end[i]] <- 1
        } else if (val == 1 & growseas[start[i]-1] == 0){
          growseas[start[i]:end[i]] <- 0
        }
      }
    }

    growseas <- as.integer(growseas)

  } else {

    warning("number of available GPPd data is less than half the total number of days per year. Filter is not applied!")
    growseas <- as.integer(rep(1,length(GPPd)))

  }

  return(growseas)
}
# ==============================================================================
# meteorological_variables.r
air.density <- function(Tair,pressure,constants=bigleaf.constants()){

  Tair     <- Tair + constants$Kelvin
  pressure <- pressure * constants$kPa2Pa

  rho <- pressure / (constants$Rd * Tair)

  return(rho)
}

pressure.from.elevation <- function(elev,Tair,VPD=NULL,constants=bigleaf.constants()){

  Tair     <- Tair + constants$Kelvin

  if(is.null(VPD)){

    pressure <- constants$pressure0 / exp(constants$g * elev / (constants$Rd*Tair))

  } else {

    pressure1   <- constants$pressure0 / exp(constants$g * elev / (constants$Rd*Tair))
    Tv          <- virtual.temp(Tair - constants$Kelvin,pressure1 * constants$Pa2kPa,
                                VPD,Esat.formula="Sonntag_1990",constants) + constants$Kelvin

    pressure    <- constants$pressure0 / exp(constants$g * elev / (constants$Rd*Tv))
  }

  pressure <- pressure * constants$Pa2kPa
  return(pressure)
}

Esat.slope <- function(Tair,formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                       constants=bigleaf.constants()){

  formula <- match.arg(formula)

  if (formula == "Sonntag_1990"){
    a <- 611.2
    b <- 17.62
    c <- 243.12
  } else if (formula == "Alduchov_1996"){
    a <- 610.94
    b <- 17.625
    c <- 243.04
  } else if (formula == "Allen_1998"){
    a <- 610.8
    b <- 17.27
    c <- 237.3
  }

  # saturation vapor pressure
  Esat <- a * exp((b * Tair) / (c + Tair))
  Esat <- Esat * constants$Pa2kPa

  # slope of the saturation vapor pressure curve
  Delta <- a * (exp((b * Tair)/(c + Tair)) * (b/(c + Tair) - (b * Tair)/(c + Tair)^2))
  Delta <- Delta * constants$Pa2kPa

  return(data.frame(Esat,Delta))
}

psychrometric.constant <- function(Tair,pressure,constants=bigleaf.constants()){

  lambda <- latent.heat.vaporization(Tair)
  gamma  <- (constants$cp * pressure) / (constants$eps * lambda)

  return(gamma)
}

latent.heat.vaporization <- function(Tair) {

  k1 <- 2.501
  k2 <- 0.00237
  lambda <- ( k1 - k2 * Tair ) * 1e+06

  return(lambda)
}

wetbulb.solver <- function(ea,Tair,gamma,accuracy,Esat.formula,constants=bigleaf.constants()){
  wetbulb.optim <- optimize(function(Tw){abs(ea - c((Esat.slope(Tw,Esat.formula,constants)[,"Esat"] - constants$Le067*gamma*(Tair - Tw))))},
                            interval=c(-100,100),tol=accuracy)
  return(wetbulb.optim)
}

wetbulb.temp <- function(Tair,pressure,VPD,accuracy=1e-03,Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                         constants=bigleaf.constants()){

  if (!is.numeric(accuracy)){
    stop("'accuracy' must be numeric!")
  }

  if (accuracy > 1){
    print("'accuracy' is set to 1 degC")
    accuracy <- 1
  }

  # determine number of digits to print
  ndigits <- as.numeric(strsplit(format(accuracy,scientific = TRUE),"-")[[1]][2])
  ndigits <- ifelse(is.na(ndigits),0,ndigits)


  gamma  <- psychrometric.constant(Tair,pressure)
  ea     <- VPD.to.e(VPD,Tair,Esat.formula)

  Tw <- sapply(seq_along(ea),function(i){ if (any(c(ea[i],Tair[i],gamma[i]) %in% c(NA,NaN,Inf))){
                                              NA
                                          } else {
                                             round(wetbulb.solver(ea[i],Tair[i],gamma[i],
                                                   accuracy=accuracy,Esat.formula,constants)$minimum,ndigits)
                                          }
                                        }
               )

  return(Tw)

}

dew.point.solver <- function(ea,accuracy,Esat.formula,constants=bigleaf.constants()){

  Td.optim <- optimize(function(Td){abs(ea - Esat.slope(Td,Esat.formula,constants)[,"Esat"])},
                       interval=c(-100,100),tol=accuracy)
  return(Td.optim)
}

dew.point <- function(Tair,VPD,accuracy=1e-03,Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                      constants=bigleaf.constants()){

  if (!is.numeric(accuracy)){
    stop("'accuracy' must be numeric!")
  }

  if (accuracy > 1){
    print("'accuracy' is set to 1 degC")
    accuracy <- 1
  }

  # determine number of digits to print
  ndigits <- as.numeric(strsplit(format(accuracy,scientific = TRUE),"-")[[1]][2])
  ndigits <- ifelse(is.na(ndigits),0,ndigits)

  ea <- VPD.to.e(VPD,Tair,Esat.formula)
  Td <- sapply(seq_along(ea),function(i){ if (ea[i] %in% c(NA,NaN,Inf)){
                                             NA
                                          } else {
                                            round(dew.point.solver(ea[i],accuracy=accuracy,
                                                                   Esat.formula,constants)$minimum,ndigits)
                                          }
                                        }
               )

  return(Td)
}

virtual.temp <- function(Tair,pressure,VPD,Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                         constants=bigleaf.constants()){

  e    <- VPD.to.e(VPD,Tair,Esat.formula)
  Tair <- Tair + constants$Kelvin

  Tv <- Tair / (1 - (1 - constants$eps) * e/pressure)
  Tv <- Tv - constants$Kelvin

  return(Tv)
}

kinematic.viscosity <- function(Tair,pressure,constants=bigleaf.constants()){

  Tair     <- Tair + constants$Kelvin
  pressure <- pressure * constants$kPa2Pa

  v  <- 1.327e-05*(constants$pressure0/pressure)*(Tair/constants$Tair0)^1.81
  return(v)
}
# ==============================================================================
# optimum_temperature.r
optimum.temperature <- function(data, GPP="GPP", Tair="Tair", BLine=0.9, Obs_filter=30){

  check.input(data, list(GPP, Tair))

  #round to 1degC temp bins
  Tair_bin <- trunc(Tair+sign(Tair)*0.5)

  #get boundary line
  df.bl <- aggregate(x=GPP,
                     by = list(Tair_bin = Tair_bin),
                     data = data,
                     FUN = quantile,
                     probs = BLine,
                     type = 8)


  n_obs <-aggregate(GPP ~ Tair_bin,
                    FUN = length)

  df.bl <-merge(df.bl, n_obs, by= c("Tair_bin"))

  colnames(df.bl) <- c("Tair_bin", "GPP_Bline", "n_obs")

  #Remove Tair bins with n_obs below filter
  df.bl <- subset(df.bl, n_obs >= Obs_filter)

  #get the smoothed boundary line
  df.bl$smooth_bl <- predict(loess(GPP_Bline~Tair_bin, df.bl), df.bl$Tair_bin)
  colnames(df.bl) <- c("Tair_bin", "GPP_Bline", "n_obs", "GPP_Bline_smooth")

  #get topt
  Topt.df <- df.bl[order(df.bl$GPP_Bline_smooth, decreasing = TRUE), ]
  Topt <- Topt.df[1,1]
  GPP_bl <- Topt.df[1,2]
  opt.temp <- c("Topt" = Topt, "GPP_bl" = GPP_bl)

  #output
  optimum.temp <- list(df.bl, opt.temp)
  names(optimum.temp) <- c("df.bl", "opt.temp")
 return(optimum.temp)

}
# ==============================================================================
# potential_radiation.r

extraterrestrial.radiation <- function(doy,constants = bigleaf.constants()){

  # Fractional year in radians
  FracYearRad <- 2 * pi * (doy - 1) / 365.24

  #Eccentricity correction
  ExtRadiation <- constants$solar_constant * (
    1.00011 + 0.034221 * cos(FracYearRad) + 0.00128 * sin(FracYearRad)
     + 0.000719 * cos(2 * FracYearRad) + 0.000077 * sin(2 * FracYearRad)
     )

  return(ExtRadiation)
}

potential.radiation <- function(doy, hour, latDeg, longDeg, timezone, useSolartime = TRUE){

  # Calculate potential radiation from solar elevation and extraterrestrial solar radiation
  solElevRad <- computeSunPositionDoyHour(
    doy, hour, latDeg, longDeg, timezone
    , isCorrectSolartime = useSolartime)[,"elevation"]

  extRadiation <- extraterrestrial.radiation(doy)

  potRad <- ifelse(
    solElevRad <= 0, 0, extRadiation * sin(solElevRad) )

  return(potRad)
}
# ==============================================================================
# stability_correction.r
Monin.Obukhov.length <- function(data,Tair="Tair",pressure="pressure",ustar="ustar",
                                 H="H",constants=bigleaf.constants()){

  check.input(data,list(Tair,pressure,ustar,H))

  rho  <- air.density(Tair,pressure,constants)
  Tair <- Tair + constants$Kelvin
  MOL  <- (-rho*constants$cp*ustar^3*Tair) / (constants$k*constants$g*H)

  return(MOL)
}

stability.parameter <- function(data,Tair="Tair",pressure="pressure",ustar="ustar",
                                H="H",zr,d,constants=bigleaf.constants()){

  check.input(data,list(Tair,pressure,ustar,H))

  MOL  <- Monin.Obukhov.length(data,Tair,pressure,ustar,H,constants)
  zeta <- (zr - d) / MOL

  return(zeta)

}

stability.correction <- function(zeta,formulation=c("Dyer_1970","Businger_1971")){

  formulation  <- match.arg(formulation)

  check.input(NULL,zeta)

  psi_h = psi_m <- rep(NA_real_,length(zeta))

  # universal functions
  if (formulation == "Businger_1971"){
    x_h <- -7.8
    x_m <- -6
    y_h <- 0.95 * ( 1 - 11.6 * zeta)^0.5
    y_m <- (1 - 19.3*zeta)^0.25
  } else if (formulation == "Dyer_1970"){
    x_h = x_m <- -5
    y_h       <- (1 - 16 * zeta)^0.5
    y_m       <- (1 - 16 * zeta)^0.25
  }

  # integration of universal functions (after Paulson_1970 and Foken 2008)
  # stable
  stable <- zeta >= 0 & !is.na(zeta)
  psi_h[stable] <- x_h * zeta[stable]
  psi_m[stable] <- x_m * zeta[stable]
  # unstable
  unstable <- zeta < 0 & !is.na(zeta)
  psi_h[unstable] <- 2 * log( (1 + y_h[unstable] ) / 2)
  psi_m[unstable] <- 2 * log( (1 + y_m[unstable] ) / 2) +
                     log( ( 1 + y_m[unstable]^2 ) / 2)
                     -2 * atan(y_m[unstable]) + pi/2

  return(data.frame(psi_h,psi_m))

}
# ==============================================================================
# surface_conditions.r
surface.conditions <- function(data,Tair="Tair",pressure="pressure",LE="LE",H="H",
                               VPD="VPD",Ga="Ga_h",calc.surface.CO2=FALSE,Ca="Ca",Ga_CO2="Ga_CO2",
                               NEE="NEE",Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                               constants=bigleaf.constants()){

  check.input(data,list(Tair,pressure,LE,H,VPD,Ga))

  rho   <- air.density(Tair,pressure,constants)
  gamma <- psychrometric.constant(Tair,pressure,constants)

  # 1) Temperature
  Tsurf <- Tair + H / (rho * constants$cp * Ga)

  # 2) Humidity
  esat      <- Esat.slope(Tair,Esat.formula,constants)[,"Esat"]
  e         <- esat - VPD
  esat_surf <- Esat.slope(Tsurf,Esat.formula,constants)[,"Esat"]
  esurf     <- e + (LE * gamma)/(Ga * rho * constants$cp)
  VPD_surf  <- pmax(esat_surf - esurf,0)
  qsurf     <- VPD.to.q(VPD_surf,Tsurf,pressure,Esat.formula,constants)
  rH_surf   <- VPD.to.rH(VPD_surf,Tsurf,Esat.formula)

  # 3) CO2 concentration
  if (calc.surface.CO2){
    check.input(data,Ca,NEE,Ga_CO2)
    Ca_surf <- surface.CO2(Ca,NEE,Ga_CO2,Tair,pressure)
  } else {
    Ca_surf <- rep(NA_integer_,length(Tair))
  }

  return(data.frame(Tsurf,esat_surf,esurf,VPD_surf,qsurf,rH_surf,Ca_surf))
}

surface.CO2 <- function(Ca,NEE,Ga_CO2,Tair,pressure){

  Ga_CO2 <- ms.to.mol(Ga_CO2,Tair,pressure)

  Ca_surf <- Ca + NEE/Ga_CO2

  return(Ca_surf)

}

radiometric.surface.temp <- function(data,LW_up="LW_up",LW_down="LW_down",
                                     emissivity,constants=bigleaf.constants()){

  check.input(data,list(LW_up,LW_down))

  Trad_K    <- ((LW_up - (1 - emissivity)*LW_down) / (constants$sigma * emissivity))^(1/4)
  Trad_degC <- Trad_K - constants$Kelvin

  return(data.frame(Trad_K,Trad_degC))
}

# ==============================================================================
# surface_conductance.r
surface.conductance <- function(data,Tair="Tair",pressure="pressure",Rn="Rn",G=NULL,S=NULL,
                                VPD="VPD",LE="LE",Ga="Ga_h",missing.G.as.NA=FALSE,missing.S.as.NA=FALSE,
                                formulation=c("Penman-Monteith","Flux-Gradient"),
                                Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                                constants=bigleaf.constants()){

  formulation <- match.arg(formulation)

  if (formulation == "Flux-Gradient"){

    check.input(data,list(Tair,pressure,VPD,LE))

    Gs_mol <- (LE.to.ET(LE,Tair)/constants$Mw) * pressure / VPD
    Gs_ms  <- mol.to.ms(Gs_mol,Tair,pressure)

  } else if (formulation == "Penman-Monteith"){

    check.input(data,list(Tair,pressure,VPD,LE,Rn,Ga,G,S))

    if(!is.null(G)){
      if (!missing.G.as.NA){G[is.na(G)] <- 0}
    } else {
      cat("Ground heat flux G is not provided and set to 0.",fill=TRUE)
      G <- 0
    }

    if(!is.null(S)){
      if(!missing.S.as.NA){S[is.na(S)] <- 0 }
    } else {
      cat("Energy storage fluxes S are not provided and set to 0.",fill=TRUE)
      S <- 0
    }

    Delta <- Esat.slope(Tair,Esat.formula,constants)[,"Delta"]
    gamma <- psychrometric.constant(Tair,pressure,constants)
    rho   <- air.density(Tair,pressure,constants)

    Gs_ms  <- ( LE * Ga * gamma ) / ( Delta * (Rn-G-S) + rho * constants$cp * Ga * VPD - LE * ( Delta + gamma ) )
    Gs_mol <- ms.to.mol(Gs_ms,Tair,pressure)

  }

  return(data.frame(Gs_ms,Gs_mol))

}
# ==============================================================================
# surface_roughness.r
Reynolds.Number <- function(Tair,pressure,ustar,z0m,constants=bigleaf.constants()){

  v  <- kinematic.viscosity(Tair,pressure,constants)
  Re <- z0m*ustar/v

  return(Re)
}

roughness.parameters <- function(method=c("canopy_height","canopy_height&LAI","wind_profile"),zh,
                                 frac_d=0.7,frac_z0m=0.1,LAI,zr,cd=0.2,hs=0.01,data,Tair="Tair",pressure="pressure",
                                 wind="wind",ustar="ustar",H="H",d=NULL,z0m=NULL,
                                 stab_roughness=TRUE,stab_formulation=c("Dyer_1970","Businger_1971"),
                                 constants=bigleaf.constants()){

  method           <- match.arg(method)
  stab_formulation <- match.arg(stab_formulation)

  if (method == "canopy_height"){

    d      <- frac_d*zh
    z0m    <- frac_z0m*zh
    z0m_se <- NA

  } else if (method == "canopy_height&LAI"){

    X <- cd * LAI
    d <- 1.1 * zh * log(1 + X^(1/4))

    if (X >= 0 & X <= 0.2){
      z0m <- hs + 0.3 * X^(1/2)
    } else {
      z0m <- 0.3 * zh * (1 - d/zh)
    }
    z0m_se <- NA

  } else if (method == "wind_profile"){

    check.input(data,Tair,pressure,wind,ustar,H)

    if (is.null(d)){

      d <- frac_d * zh

    }

    if (stab_roughness){

      zeta  <- stability.parameter(data=data,Tair=Tair,pressure=pressure,ustar=ustar,H=H,
                                   zr=zr,d=d,constants=constants)
      psi_m <- stability.correction(zeta,formulation=stab_formulation)[,"psi_m"]
      z0m_all <- (zr - d) * exp(-constants$k*wind / ustar - psi_m)

    } else {

      z0m_all <- (zr - d) * exp(-constants$k*wind / ustar)

    }

    z0m_all[z0m_all > zh] <- NA

    z0m    <- median(z0m_all,na.rm=TRUE)
    z0m_se <- constants$se_median * (sd(z0m_all,na.rm=TRUE) / sqrt(length(z0m_all[complete.cases(z0m_all)])))

  }

  return(data.frame(d,z0m,z0m_se))
}

wind.profile <- function(data,z,Tair="Tair",pressure="pressure",ustar="ustar",H="H",wind="wind",
                         zr,zh,d=NULL,frac_d=0.7,z0m=NULL,frac_z0m=NULL,estimate_z0m=TRUE,
                         stab_correction=TRUE,stab_formulation=c("Dyer_1970","Businger_1971"),
                         constants=bigleaf.constants()){

  stab_formulation <- match.arg(stab_formulation)

  check.input(data,ustar)

  ## determine roughness parameters
  if (is.null(d)){
    if (is.null(frac_d)){
      stop("Either 'd' or 'frac_d' must be specified")
    }
    d <- frac_d * zh
  }

  if (is.null(z0m) & !estimate_z0m){
    if (is.null(frac_z0m)){
      stop("Either 'z0m' or 'frac_z0m' must be specified if 'estimate_z0m' = FALSE")
    }
    z0m <- frac_z0m * zh
  }


  if (estimate_z0m){

    if (!is.null(z0m) | !is.null(frac_z0m)){
      cat("Note that arguments 'z0m' and 'frac_z0m' are ignored if 'estimate_z0m' = TRUE. z0m is
           calculated from the logarithmic wind_profile equation.",fill=TRUE)
    }

    check.input(data,Tair,pressure,wind,ustar,H)

    z0m <- roughness.parameters(method="wind_profile",zh=zh,zr=zr,d=d,data=data,
                                Tair=Tair,pressure=pressure,wind=wind,ustar=ustar,H=H,
                                stab_roughness=TRUE,stab_formulation=stab_formulation,
                                constants=constants)[,"z0m"]
  }

  if ( any(z < (d + z0m) & !is.na(d + z0m)) ){
    warning("function is only valid for heights above d + z0m! Wind speed for heights below d + z0m will return 0!")
  }

  ## calculate wind speeds at given heights z
  if (stab_correction){

    zeta  <- stability.parameter(data=data,Tair=Tair,pressure=pressure,ustar=ustar,H=H,
                                 zr=z,d=d,constants=constants)
    psi_m <- stability.correction(zeta,formulation=stab_formulation)[,"psi_m"]
    wind_heights <- pmax(0,(ustar / constants$k) * (log(pmax(0,(z - d)) / z0m) - psi_m))

  } else {

    wind_heights <- pmax(0,(ustar / constants$k) * (log(pmax(0,(z - d)) / z0m)))

  }

  return(wind_heights)
}

# ==============================================================================
# unit_conversions.r
LE.to.ET <- function(LE,Tair){

  lambda <- latent.heat.vaporization(Tair)
  ET     <- LE/lambda

  return(ET)
}

ET.to.LE <- function(ET,Tair){

  lambda <- latent.heat.vaporization(Tair)
  LE     <- ET*lambda

  return(LE)
}

ms.to.mol <- function(G_ms,Tair,pressure,constants=bigleaf.constants()){

  Tair     <- Tair + constants$Kelvin
  pressure <- pressure * constants$kPa2Pa

  G_mol  <- G_ms * pressure / (constants$Rgas * Tair)

  return(G_mol)
}

mol.to.ms <- function(G_mol,Tair,pressure,constants=bigleaf.constants()){

  Tair     <- Tair + constants$Kelvin
  pressure <- pressure * constants$kPa2Pa

  G_ms  <- G_mol * (constants$Rgas * Tair) / (pressure)

  return(G_ms)
}

VPD.to.rH <- function(VPD,Tair,Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                      constants=bigleaf.constants()){

  esat <- Esat.slope(Tair,Esat.formula,constants)[,"Esat"]
  rH   <- 1 - VPD/esat
  return(rH)
}

rH.to.VPD <- function(rH,Tair,Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                      constants=bigleaf.constants()){

  if(any(rH > 1 & !is.na(rH))){
    warning("relative humidity (rH) has to be between 0 and 1.")
  }

  esat <- Esat.slope(Tair,Esat.formula,constants)[,"Esat"]
  VPD  <- esat - rH*esat
  return(VPD)
}

e.to.rH <- function(e,Tair,Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                    constants=bigleaf.constants()){

  esat <- Esat.slope(Tair,Esat.formula,constants)[,"Esat"]

  if (any(e > esat + .Machine$double.eps^0.5 & !is.na(e))){
    warning("Provided vapour pressure that was higher than saturation.
             Returning rH=1 for those cases.")
  }

  rH  <- pmin(1, e/esat)
  return(rH)
}

VPD.to.e <- function(VPD,Tair,Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                     constants=bigleaf.constants()){

  esat <- Esat.slope(Tair,Esat.formula,constants)[,"Esat"]
  e    <- esat - VPD
  return(e)
}

e.to.VPD <- function(e,Tair,Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                     constants=bigleaf.constants()){

  esat <- Esat.slope(Tair,Esat.formula,constants)[,"Esat"]
  VPD  <- esat - e
  return(VPD)
}

e.to.q <- function(e,pressure,constants=bigleaf.constants()){
  q <- constants$eps * e / (pressure - (1-constants$eps) * e)
  return(q)
}

q.to.e <- function(q,pressure,constants=bigleaf.constants()){
  e <- q * pressure / ((1-constants$eps) * q + constants$eps)
  return(e)
}

q.to.VPD <- function(q,Tair,pressure,Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                     constants=bigleaf.constants()){

  esat <- Esat.slope(Tair,Esat.formula,constants)[,"Esat"]
  e    <- q.to.e(q,pressure,constants)
  VPD  <- esat - e
  return(VPD)
}

VPD.to.q <- function(VPD,Tair,pressure,Esat.formula=c("Sonntag_1990","Alduchov_1996","Allen_1998"),
                     constants=bigleaf.constants()){

  esat <- Esat.slope(Tair,Esat.formula,constants)[,"Esat"]
  e    <- esat - VPD
  q    <- e.to.q(e,pressure,constants)
  return(q)
}

Rg.to.PPFD <- function(Rg,J_to_mol=4.6,frac_PAR=0.5){
  PPFD <- Rg * frac_PAR * J_to_mol
  return(PPFD)
}

PPFD.to.Rg <- function(PPFD,J_to_mol=4.6,frac_PAR=0.5){
  Rg <- PPFD / frac_PAR / J_to_mol
  return(Rg)
}

kg.to.mol <- function(mass, molarMass=bigleaf.constants()$H2Omol){

  moles <- mass / molarMass

  return(moles)

}
umolCO2.to.gC <- function(CO2_flux,constants=bigleaf.constants()){

  C_flux <- CO2_flux * constants$umol2mol * constants$Cmol * constants$kg2g * constants$days2seconds

  return(C_flux)
}

gC.to.umolCO2 <- function(C_flux,constants=bigleaf.constants()){

  CO2_flux <- (C_flux * constants$g2kg / constants$days2seconds) / constants$Cmol * constants$mol2umol

  return(CO2_flux)
}
# ==============================================================================

