import numpy as np
import pandas as pd
from .bigleaf_constants import *
from .stability_correction import stability_parameter, stability_correction
from .meteorological_variables import kinematic_viscosity

def reynolds_number(Tair, pressure, ustar, z0m, constants):
    v = kinematic_viscosity(Tair, pressure, constants)
    Re = z0m * ustar / v
    return Re

def roughness_parameters(method, zh, frac_d=0.7, frac_z0m=0.1, LAI=None, zr=None, cd=0.2, hs=0.01,
                         data=None, Tair="Tair", pressure="pressure", wind="wind",
                         ustar="ustar", H="H", d=None, z0m=None,
                         stab_roughness=True, stab_formulation="Dyer_1970", constants=None):

    if constants is None:
        constants = bigleaf_constants()

    if method == "canopy_height":
        d = frac_d * zh
        z0m = frac_z0m * zh
        z0m_se = np.nan

    elif method == "canopy_height&LAI":
        X = cd * LAI
        d = 1.1 * zh * np.log(1 + X**0.25)
        if 0 <= X <= 0.2:
            z0m = hs + 0.3 * X**0.5
        else:
            z0m = 0.3 * zh * (1 - d / zh)
        z0m_se = np.nan

    elif method == "wind_profile":
        # check_input(data, [Tair, pressure, wind, ustar, H])

        if d is None:
            d = frac_d * zh

        wind = data[wind]
        ustar = data[ustar]

        if stab_roughness:
            zeta = stability_parameter(data, Tair, pressure, ustar, H, zr, d, constants)
            psi_m = stability_correction(zeta, formulation=stab_formulation)["psi_m"]
            z0m_all = (zr - d) * np.exp(-constants['k'] * wind / ustar - psi_m)
        else:
            z0m_all = (zr - d) * np.exp(-constants['k'] * wind / ustar)

        z0m_all[z0m_all > zh] = np.nan
        z0m = np.nanmedian(z0m_all)
        valid = z0m_all.dropna()
        z0m_se = constants['se_median'] * (np.std(valid) / np.sqrt(len(valid)))

    return pd.DataFrame({"d": [d], "z0m": [z0m], "z0m_se": [z0m_se]})

def wind_profile(data, z, Tair="Tair", pressure="pressure", ustar="ustar", H="H",
                 wind="wind", zr=None, zh=None, d=None, frac_d=0.7, z0m=None, frac_z0m=0.1,
                 estimate_z0m=True, stab_correction=True, stab_formulation="Dyer_1970",
                 constants=None):

    if constants is None:
        constants = bigleaf_constants()

    # check_input(data, [ustar])

    if d is None:
        d = frac_d * zh

    if estimate_z0m:
        rough = roughness_parameters("wind_profile", zh=zh, zr=zr, frac_d=frac_d, data=data,
                                     Tair=Tair, pressure=pressure, wind=wind,
                                     ustar=ustar, H=H, stab_roughness=stab_correction,
                                     stab_formulation=stab_formulation, constants=constants)
        z0m = rough["z0m"].values[0]

    z_diff = z - d
    z_diff[z_diff <= z0m] = np.nan

    if stab_correction:
        zeta = stability_parameter(data, Tair, pressure, ustar, H, zr, d, constants)
        psi_m = stability_correction(zeta, formulation=stab_formulation)["psi_m"]
    else:
        psi_m = 0

    ustar = data[ustar]
    uz = (ustar / constants['k']) * (np.log(z_diff / z0m) - psi_m)
    uz[z_diff <= z0m] = 0

    return uz