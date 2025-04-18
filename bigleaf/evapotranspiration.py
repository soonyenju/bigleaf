import numpy as np
import pandas as pd

def potential_ET(data=None, Tair="Tair", pressure="pressure", Rn="Rn", G=None, S=None,
                 VPD="VPD", Ga="Ga", approach="Priestley-Taylor", alpha=1.26,
                 Gs_pot=0.6, missing_G_as_NA=False, missing_S_as_NA=False,
                 Esat_formula="Sonntag_1990", constants=None):

    if constants is None:
        constants = bigleaf_constants()

    if approach not in ["Priestley-Taylor", "Penman-Monteith"]:
        raise ValueError("approach must be either 'Priestley-Taylor' or 'Penman-Monteith'")

    # Extract data if necessary
    def extract(var):
        return data[var] if isinstance(var, str) else var

    Tair = extract(Tair)
    pressure = extract(pressure)
    Rn = extract(Rn)
    if G is not None:
        G = extract(G)
        if not missing_G_as_NA:
            G = np.nan_to_num(G)
    else:
        print("Ground heat flux G is not provided and set to 0.")
        G = 0

    if S is not None:
        S = extract(S)
        if not missing_S_as_NA:
            S = np.nan_to_num(S)
    else:
        print("Energy storage fluxes S are not provided and set to 0.")
        S = 0

    gamma = psychrometric_constant(Tair, pressure, constants)
    Delta = Esat_slope(Tair, Esat_formula, constants)["Delta"]

    if approach == "Priestley-Taylor":
        LE_pot = (alpha * Delta * (Rn - G - S)) / (Delta + gamma)
        ET_pot = LE_to_ET(LE_pot, Tair)

    elif approach == "Penman-Monteith":
        VPD = extract(VPD)
        Ga = extract(Ga)
        Gs_pot = mol_to_ms(Gs_pot, Tair=Tair, pressure=pressure, constants=constants)
        rho = air_density(Tair, pressure, constants)
        numerator = Delta * (Rn - G - S) + rho * constants['cp'] * VPD * Ga
        denominator = Delta + gamma * (1 + Ga / Gs_pot)
        LE_pot = numerator / denominator
        ET_pot = LE_to_ET(LE_pot, Tair)

    return pd.DataFrame({'ET_pot': ET_pot, 'LE_pot': LE_pot})