import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, minimize_scalar

# Constants
# def bigleaf_constants():
#     return {
#         'DwDc': 1.6,
#         'Kelvin': 273.15,
#         'Rgas': 8.314,
#         'kJ2J': 1000,
#         'J2kJ': 0.001,
#         'se_median': 1.253
#     }

# 1. Intercellular CO2
def intercellular_CO2(data, Ca='Ca', GPP='GPP', Gs='Gs_mol', Rleaf=None, missing_Rleaf_as_NA=False, constants=None):
    constants = constants or bigleaf_constants()
    Rleaf = 0 if Rleaf is None else data[Rleaf] if missing_Rleaf_as_NA else data[Rleaf].fillna(0)
    Ci = data[Ca] - (data[GPP] - Rleaf) / (data[Gs] / constants['DwDc'])
    return Ci

# 2. Arrhenius Temperature Response
def arrhenius_temp_response(param, Temp, Ha, Hd=None, dS=None, constants=None):
    constants = constants or bigleaf_constants()
    Temp += constants['Kelvin']
    Tref = 25.0 + constants['Kelvin']
    Ha *= constants['kJ2J']
    if Hd is not None: Hd *= constants['kJ2J']
    if dS is not None: dS *= constants['kJ2J']
    if Hd is None or dS is None:
        return param / np.exp(Ha * (Temp - Tref) / (Tref * constants['Rgas'] * Temp))
    term1 = np.exp(Ha * (Temp - Tref) / (Tref * constants['Rgas'] * Temp))
    term2 = (1 + np.exp((Tref * dS - Hd) / (Tref * constants['Rgas'])))
    term3 = (1 + np.exp((Temp * dS - Hd) / (Temp * constants['Rgas'])))
    return param / (term1 * term2 / term3)

# 3. Photosynthetic Capacity
def photosynthetic_capacity(data, Temp, GPP='GPP', Ci=None, PPFD='PPFD', Rleaf=None,
                            Oi=0.21, Kc25=404.9, Ko25=278.4, Gam25=42.75,
                            Kc_Ha=79.43, Ko_Ha=36.38, Gam_Ha=37.83,
                            Vcmax_Ha=65.33, Vcmax_Hd=200, Vcmax_dS=0.635,
                            Jmax_Ha=43.9, Jmax_Hd=200, Jmax_dS=0.640,
                            Theta=0.7, alpha_canopy=0.8, missing_Rleaf_as_NA=False,
                            Ci_C4=100, C3=True, PPFD_j=(200, 500), PPFD_c=1000,
                            constants=None):
    constants = constants or bigleaf_constants()
    TempK = data[Temp] + constants['Kelvin']
    Tref = 25.0 + constants['Kelvin']
    Rleaf = 0 if Rleaf is None else data[Rleaf] if missing_Rleaf_as_NA else data[Rleaf].fillna(0)
    Kc_Ha *= constants['kJ2J']
    Ko_Ha *= constants['kJ2J']
    Gam_Ha *= constants['kJ2J']
    Kc = Kc25 * np.exp(Kc_Ha * (TempK - Tref) / (Tref * constants['Rgas'] * TempK))
    Ko = Ko25 * np.exp(Ko_Ha * (TempK - Tref) / (Tref * constants['Rgas'] * TempK))
    Gam = Gam25 * np.exp(Gam_Ha * (TempK - Tref) / (Tref * constants['Rgas'] * TempK))
    Ko *= constants['J2kJ']
    Ci = Ci.copy()
    Ci[Ci < 80] = np.nan
    GPP_data = data[GPP]
    PPFD_data = data[PPFD]
    GPPc = GPP_data.where(PPFD_data >= PPFD_c)
    GPPj = GPP_data.where((PPFD_data >= PPFD_j[0]) & (PPFD_data <= PPFD_j[1]))
    Vcmax = (GPPc - Rleaf) * (Ci + Kc * (1 + Oi / Ko)) / (Ci - Gam)
    J = (GPPj - Rleaf) * (4 * Ci + 8 * Gam) / (Ci - Gam)
    APPFD_PSII = PPFD_data * alpha_canopy * 0.85 * 0.5
    Jmax = []
    for i in J.dropna().index:
        result = minimize_scalar(
            lambda jmax: abs(J[i] - ((APPFD_PSII[i] + jmax - np.sqrt((APPFD_PSII[i] + jmax)**2 - 4 * Theta * APPFD_PSII[i] * jmax)) / (2 * Theta))),
            bounds=(0, 1000), method='bounded'
        )
        Jmax.append(result.x if result.success else np.nan)
    Jmax = np.array(Jmax)
    Vcmax25 = arrhenius_temp_response(Vcmax, TempK - constants['Kelvin'], Vcmax_Ha, Vcmax_Hd, Vcmax_dS, constants)
    Jmax25 = arrhenius_temp_response(Jmax, TempK.loc[J.dropna().index] - constants['Kelvin'], Jmax_Ha, Jmax_Hd, Jmax_dS, constants)
    return {
        "Vcmax25": round(np.nanmedian(Vcmax25), 2),
        "Vcmax25_SE": round(constants['se_median'] * np.nanstd(Vcmax25) / np.sqrt(np.sum(~np.isnan(Vcmax25))), 2),
        "Jmax25": round(np.nanmedian(Jmax25), 2),
        "Jmax25_SE": round(constants['se_median'] * np.nanstd(Jmax25) / np.sqrt(np.sum(~np.isnan(Jmax25))), 2)
    }

# 4. Stomatal Slope (Medlyn model)
def stomatal_slope(data, GPP='GPP', Ca='Ca', VPD='VPD', Ci=None, constants=None):
    constants = constants or bigleaf_constants()
    Ci = data[Ci] if isinstance(Ci, str) else Ci
    VPD = data[VPD]
    Ca = data[Ca]
    GPP = data[GPP]
    m = GPP * Ca / (Ci - constants['DwDc'] * Ci * VPD)
    return m

# 5. Light Response (Rectangular hyperbola)
def light_response(PPFD, alpha, Amax, theta):
    return ((alpha * PPFD + Amax) - np.sqrt((alpha * PPFD + Amax)**2 - 4 * alpha * PPFD * Amax * theta)) / (2 * theta)

# 6. Light Use Efficiency
def light_use_efficiency(data, PPFD='PPFD', GPP='GPP', model=False, start_vals=(0.01, 20, 0.85), bounds=None):
    PPFD = data[PPFD].values
    GPP = data[GPP].values
    LUE = GPP / PPFD
    if not model:
        return np.nanmedian(LUE)
    bounds = bounds or ([0, 0, 0], [np.inf, np.inf, 1])
    popt, _ = curve_fit(light_response, PPFD, GPP, p0=start_vals, bounds=bounds)
    return popt[0]  # alpha

# 7. Stomatal Sensitivity (Slope of ln(GPP) vs ln(VPD))
def stomatal_sensitivity(data, GPP='GPP', VPD='VPD', constants=None):
    constants = constants or bigleaf_constants()
    ln_GPP = np.log(data[GPP])
    ln_VPD = np.log(data[VPD])
    mask = (~np.isnan(ln_GPP)) & (~np.isnan(ln_VPD)) & (ln_VPD != -np.inf)
    if mask.sum() < 2:
        return np.nan
    slope, _ = np.polyfit(ln_VPD[mask], ln_GPP[mask], 1)
    return slope