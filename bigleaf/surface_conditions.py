import numpy as np
import pandas as pd

def surface_conditions(Tair, pressure, H, LE, VPD, Ga,
                       calc_surface_CO2=False, Ca=None, NEE=None, Ga_CO2=None,
                       Esat_formula="Sonntag_1990", constants=None):
    """
    Calculate big-leaf surface conditions (temperature, vapor pressure, humidity, CO2).
    """
    if constants is None:
        raise ValueError("Constants must be provided.")

    rho = air_density(Tair, pressure, constants)
    gamma = psychrometric_constant(Tair, pressure, constants)

    Tsurf = Tair + H / (rho * constants["cp"] * Ga)

    esat = esat_slope(Tair, Esat_formula, constants)["Esat"]
    e = esat - VPD
    esat_surf = esat_slope(Tsurf, Esat_formula, constants)["Esat"]
    esurf = e + (LE * gamma) / (Ga * rho * constants["cp"])
    VPD_surf = np.maximum(esat_surf - esurf, 0)
    qsurf = VPD_to_q(VPD_surf, Tsurf, pressure, Esat_formula, constants)
    rH_surf = VPD_to_rH(VPD_surf, Tsurf, Esat_formula)

    if calc_surface_CO2:
        if Ca is None or NEE is None or Ga_CO2 is None:
            raise ValueError("Ca, NEE, and Ga_CO2 must be provided if calc_surface_CO2 is True.")
        Ca_surf = surface_CO2(Ca, NEE, Ga_CO2, Tair, pressure)
    else:
        Ca_surf = np.full_like(np.array(Tair), np.nan)

    return pd.DataFrame({
        "Tsurf": Tsurf,
        "esat_surf": esat_surf,
        "esurf": esurf,
        "VPD_surf": VPD_surf,
        "qsurf": qsurf,
        "rH_surf": rH_surf,
        "Ca_surf": Ca_surf
    })

def surface_CO2(Ca, NEE, Ga_CO2, Tair, pressure):
    """
    Calculate CO2 concentration at the canopy surface.
    """
    Ga_CO2_mol = ms_to_mol(Ga_CO2, Tair, pressure)
    return Ca + NEE / Ga_CO2_mol

def radiometric_surface_temp(LW_up, LW_down, emissivity, constants):
    """
    Calculate radiometric surface temperature from longwave radiation.
    """
    Trad_K = ((LW_up - (1 - emissivity) * LW_down) / (constants["sigma"] * emissivity)) ** 0.25
    Trad_degC = Trad_K - constants["Kelvin"]
    return pd.DataFrame({"Trad_K": Trad_K, "Trad_degC": Trad_degC})