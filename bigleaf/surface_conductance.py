import numpy as np
import pandas as pd

def surface_conductance(data, Tair="Tair", pressure="pressure", Rn="Rn", G=None, S=None,
                        VPD="VPD", LE="LE", Ga="Ga_h", missing_G_as_NA=False, missing_S_as_NA=False,
                        formulation="Penman-Monteith", Esat_formula="Sonntag_1990",
                        constants=bigleaf_constants()):

    formulation = formulation if formulation in ["Penman-Monteith", "Flux-Gradient"] else "Penman-Monteith"

    if formulation == "Flux-Gradient":
        check_input(data, [Tair, pressure, VPD, LE])

        Gs_mol = (LE_to_ET(data[LE], data[Tair], constants) / constants["Mw"]) * data[pressure] / data[VPD]
        Gs_ms = mol_to_ms(Gs_mol, data[Tair], data[pressure], constants)

    elif formulation == "Penman-Monteith":
        check_input(data, [Tair, pressure, VPD, LE, Rn, Ga])

        G_col = data[G] if G is not None else 0
        if G is not None and not missing_G_as_NA:
            G_col = G_col.fillna(0)

        S_col = data[S] if S is not None else 0
        if S is not None and not missing_S_as_NA:
            S_col = S_col.fillna(0)

        Delta = Esat_slope(data[Tair], Esat_formula, constants)["Delta"]
        gamma = psychrometric_constant(data[Tair], data[pressure], constants)
        rho = air_density(data[Tair], data[pressure], constants)

        Rn_net = data[Rn] - G_col - S_col
        Ga_col = data[Ga]
        VPD_col = data[VPD]
        LE_col = data[LE]

        numerator = LE_col * Ga_col * gamma
        denominator = Delta * Rn_net + rho * constants["cp"] * Ga_col * VPD_col - LE_col * (Delta + gamma)

        Gs_ms = numerator / denominator
        Gs_mol = ms_to_mol(Gs_ms, data[Tair], data[pressure], constants)

    return pd.DataFrame({"Gs_ms": Gs_ms, "Gs_mol": Gs_mol})