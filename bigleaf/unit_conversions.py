# Unit conversions - Python version of R's bigleaf conversion functions

def LE_to_ET(LE, Tair):
    lambda_ = latent_heat_vaporization(Tair)
    return LE / lambda_

def ET_to_LE(ET, Tair):
    lambda_ = latent_heat_vaporization(Tair)
    return ET * lambda_

def ms_to_mol(G_ms, Tair, pressure, constants):
    Tair_K = Tair + constants['Kelvin']
    pressure_Pa = pressure * constants['kPa2Pa']
    return G_ms * pressure_Pa / (constants['Rgas'] * Tair_K)

def mol_to_ms(G_mol, Tair, pressure, constants):
    Tair_K = Tair + constants['Kelvin']
    pressure_Pa = pressure * constants['kPa2Pa']
    return G_mol * (constants['Rgas'] * Tair_K) / pressure_Pa

def VPD_to_rH(VPD, Tair, Esat_slope, Esat_formula, constants):
    esat = Esat_slope(Tair, Esat_formula, constants)["Esat"]
    return 1 - VPD / esat

def rH_to_VPD(rH, Tair, Esat_slope, Esat_formula, constants):
    if any(r > 1 for r in rH if r is not None):
        print("Warning: relative humidity (rH) has to be between 0 and 1.")
    esat = Esat_slope(Tair, Esat_formula, constants)["Esat"]
    return esat - rH * esat

def e_to_rH(e, Tair, Esat_slope, Esat_formula, constants):
    esat = Esat_slope(Tair, Esat_formula, constants)["Esat"]
    if any(val > es + 1e-15 for val, es in zip(e, esat) if val is not None):
        print("Warning: Provided vapour pressure was higher than saturation. Returning rH=1 for those cases.")
    return [min(1, val/es) for val, es in zip(e, esat)]

def VPD_to_e(VPD, Tair, Esat_slope, Esat_formula, constants):
    esat = Esat_slope(Tair, Esat_formula, constants)["Esat"]
    return esat - VPD

def e_to_VPD(e, Tair, Esat_slope, Esat_formula, constants):
    esat = Esat_slope(Tair, Esat_formula, constants)["Esat"]
    return esat - e

def e_to_q(e, pressure, constants):
    return constants['eps'] * e / (pressure - (1 - constants['eps']) * e)

def q_to_e(q, pressure, constants):
    return q * pressure / ((1 - constants['eps']) * q + constants['eps'])

def q_to_VPD(q, Tair, pressure, Esat_slope, Esat_formula, constants):
    esat = Esat_slope(Tair, Esat_formula, constants)["Esat"]
    e = q_to_e(q, pressure, constants)
    return esat - e

def VPD_to_q(VPD, Tair, pressure, Esat_slope, Esat_formula, constants):
    esat = Esat_slope(Tair, Esat_formula, constants)["Esat"]
    e = esat - VPD
    return e_to_q(e, pressure, constants)

def Rg_to_PPFD(Rg, J_to_mol=4.6, frac_PAR=0.5):
    return Rg * frac_PAR * J_to_mol

def PPFD_to_Rg(PPFD, J_to_mol=4.6, frac_PAR=0.5):
    return PPFD / (frac_PAR * J_to_mol)

def kg_to_mol(mass, molar_mass):
    return mass / molar_mass

def umolCO2_to_gC(CO2_flux, constants):
    return CO2_flux * constants['umol2mol'] * constants['Cmol'] * constants['kg2g'] * constants['days2seconds']

def gC_to_umolCO2(C_flux, constants):
    return (C_flux * constants['g2kg'] / constants['days2seconds']) / constants['Cmol'] * constants['mol2umol']