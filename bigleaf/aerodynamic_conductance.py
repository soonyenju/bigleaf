import numpy as np
import pandas as pd
from .bigleaf_constants import *
from .boundary_layer_conductance import Gb_Thom, Gb_Choudhury, Gb_Su
from .surface_roughness import roughness_parameters
from .boundary_layer_conductance import roughness_length_heat
from .stability_correction import stability_parameter, stability_correction

def aerodynamic_conductance(data, Tair="Tair", pressure="pressure", wind="wind", ustar="ustar", H="H",
                            zr=None, zh=None, d=None, z0m=None, Dl=None, N=2, fc=None, LAI=None,
                            Cd=0.2, hs=0.01, wind_profile=False, stab_correction=True,
                            stab_formulation="Dyer_1970", Rb_model="Thom_1972", kB_h=None,
                            Sc=None, Sc_name=None, constants=None):
    """
    Python version of bigleaf::aerodynamic.conductance
    """

    if constants is None:
        constants = bigleaf_constants()

    if Rb_model not in ["Thom_1972", "Choudhury_1988", "Su_2001", "constant_kB-1"]:
        raise ValueError("Invalid Rb_model")

    if stab_formulation not in ["Dyer_1970", "Businger_1971"]:
        raise ValueError("Invalid stab_formulation")

    # Convert columns to arrays
    T = data[Tair].values
    P = data[pressure].values
    U = data[wind].values
    ustar_vals = data[ustar].values
    H_vals = data[H].values

    if Rb_model in ["Thom_1972", "Choudhury_1988", "Su_2001"]:
        if Rb_model == "Thom_1972":
            Gb_mod = Gb_Thom(data[ustar], Sc=Sc, Sc_name=Sc_name, constants=constants)
        elif Rb_model == "Choudhury_1988":
            Gb_mod = Gb_Choudhury(data, Dl, LAI, zh, zr, d, Tair, pressure, wind, ustar, H,
                                  z0m, stab_formulation, Sc, Sc_name, constants)
        elif Rb_model == "Su_2001":
            Gb_mod = Gb_Su(data, T, P, ustar_vals, U, H_vals, zh, zr, d, z0m, Dl, N, fc,
                           LAI, Cd, hs, stab_formulation, Sc, Sc_name, constants)

        kB_h = Gb_mod["kB_h"]
        Rb_h = Gb_mod["Rb_h"]
        Gb_h = Gb_mod["Gb_h"]
        Gb_x = Gb_mod.filter(regex="^Gb_").drop(columns=["Gb_h"])

    elif Rb_model == "constant_kB-1":
        if kB_h is None:
            raise ValueError("kB-1 must be provided for constant_kB-1 model")

        Rb_h = kB_h / (constants["k"] * ustar_vals)
        Gb_h = 1 / Rb_h

        if Sc is not None or Sc_name is not None:
            if len(Sc) != len(Sc_name):
                raise ValueError("Sc and Sc_name must have same length")
            if not all(isinstance(s, (int, float)) for s in Sc):
                raise TypeError("Sc must be numeric")

        Sc_full = [constants["Sc_CO2"]] + list(Sc or [])
        names = ["CO2"] + list(Sc_name or [])
        Gb_x = pd.DataFrame({
            f"Gb_{name}": Gb_h / (sc / constants["Pr"])**0.67
            for name, sc in zip(names, Sc_full)
        })

    # Calculate aerodynamic resistance
    if wind_profile:
        if z0m is None and Rb_model in ["constant_kB-1", "Thom_1972"]:
            raise ValueError("z0m must be provided when wind_profile=True")

        if z0m is None and Rb_model in ["Choudhury_1988", "Su_2001"]:
            z0m = roughness_parameters("wind_profile", zh, zr, d, data, Tair, pressure,
                                       wind, ustar, H, True, stab_formulation, constants)["z0m"]

        if stab_correction:
            zeta = stability_parameter(data, Tair, pressure, ustar, H, zr, d, constants)
            psi_h = stability_correction(zeta, stab_formulation)["psi_h"]
            Ra_m = np.maximum((np.log((zr - d) / z0m) - psi_h), 0) / (constants["k"] * ustar_vals)
        else:
            Ra_m = np.maximum(np.log((zr - d) / z0m), 0) / (constants["k"] * ustar_vals)
            zeta = psi_h = np.full_like(Ra_m, np.nan)
    else:
        Ra_m = U / ustar_vals**2
        zeta = psi_h = np.full_like(Ra_m, np.nan)

    # Calculate conductances
    Ga_m = 1 / Ra_m
    Ra_h = Ra_m + Rb_h
    Ga_h = 1 / Ra_h
    Ga_x = 1 / (Ra_m[:, np.newaxis] + 1 / Gb_x.values)
    Ra_CO2 = 1 / Ga_x[:, 0]
    Ga_x_df = pd.DataFrame(Ga_x, columns=[f"Ga_{col[3:]}" for col in Gb_x.columns])

    if z0m is not None:
        z0h = roughness_length_heat(z0m, kB_h)
    else:
        z0h = np.full_like(Ra_m, np.nan)

    # Interleave Ga and Gb columns
    Gab_x = pd.concat([Ga_x_df, Gb_x], axis=1)
    Ga_cols = Ga_x_df.columns.tolist()
    Gb_cols = Gb_x.columns.tolist()
    interleaved_cols = [val for pair in zip(Ga_cols, Gb_cols) for val in pair]
    Gab_x = Gab_x[interleaved_cols]

    # Assemble final output
    result = pd.DataFrame({
        "Ga_m": Ga_m,
        "Ra_m": Ra_m,
        "Ga_h": Ga_h,
        "Ra_h": Ra_h,
        "Gb_h": Gb_h,
        "Rb_h": Rb_h,
        "kB_h": kB_h,
        "z0h": z0h,
        "zeta": zeta,
        "psi_h": psi_h,
        "Ra_CO2": Ra_CO2
    })
    result = pd.concat([result, Gab_x], axis=1)

    return result
