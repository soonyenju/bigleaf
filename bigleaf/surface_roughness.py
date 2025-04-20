import numpy as np
import pandas as pd
from .bigleaf_constants import *
from .stability_correction import stability_parameter, stability_correction
from .meteorological_variables import kinematic_viscosity

def reynolds_number(Tair, pressure, ustar, z0m, constants=None):
    if constants is None:
        constants = bigleaf_constants()
    v = kinematic_viscosity(Tair, pressure, constants)
    Re = z0m * ustar / v
    return Re

def roughness_parameters(method, zh, frac_d=0.7, frac_z0m=0.1, LAI=None, zr=None, cd=0.2, hs=0.01,
                         data=None, Tair="Tair", pressure="pressure", wind="wind",
                         ustar="ustar", H="H", d=None, z0m=None,
                         stab_roughness=True, stab_formulation="Dyer_1970", constants=None):
    if zr is None:
        raise ValueError("zr must be specified")
    if zh is None:
        raise ValueError("zh must be specified")

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

        # wind = data[wind]
        # ustar = data[ustar]

        if stab_roughness:
            zeta = stability_parameter(data, Tair, pressure, ustar, H, zr, d, constants)
            psi_m = stability_correction(zeta, formulation=stab_formulation)["psi_m"]
            z0m_all = (zr - d) * np.exp(-constants['k'] * data[wind] / data[ustar] - psi_m)
        else:
            z0m_all = (zr - d) * np.exp(-constants['k'] * data[wind] / data[ustar])

        print(z0m_all)
        z0m_all[z0m_all > zh] = np.nan
        z0m = np.nanmedian(z0m_all)
        valid = z0m_all.dropna()
        z0m_se = constants['se_median'] * (np.std(valid) / np.sqrt(len(valid)))

    return pd.DataFrame({"d": [d], "z0m": [z0m], "z0m_se": [z0m_se]})

# def wind_profile(data, z, Tair="Tair", pressure="pressure", ustar="ustar", H="H",
#                  wind="wind", zr=None, zh=None, d=None, frac_d=0.7, z0m=None, frac_z0m=0.1,
#                  estimate_z0m=True, stab_correction=True, stab_formulation="Dyer_1970",
#                  constants=None):

#     if constants is None:
#         constants = bigleaf_constants()

#     # check_input(data, [ustar])

#     if d is None:
#         d = frac_d * zh

#     if estimate_z0m:
#         rough = roughness_parameters("wind_profile", zh=zh, zr=zr, frac_d=frac_d, data=data,
#                                      Tair=Tair, pressure=pressure, wind=wind,
#                                      ustar=ustar, H=H, stab_roughness=stab_correction,
#                                      stab_formulation=stab_formulation, constants=constants)
#         z0m = rough["z0m"].values[0]

#     z_diff = z - d
#     z_diff[z_diff <= z0m] = np.nan

#     if stab_correction:
#         zeta = stability_parameter(data, Tair, pressure, ustar, H, zr, d, constants)
#         psi_m = stability_correction(zeta, formulation=stab_formulation)["psi_m"]
#     else:
#         psi_m = 0

#     ustar = data[ustar]
#     uz = (ustar / constants['k']) * (np.log(z_diff / z0m) - psi_m)
#     uz[z_diff <= z0m] = 0

#     return uz

# import numpy as np
# import pandas as pd

def wind_profile(data, z, Tair="Tair", pressure="pressure", ustar="ustar", H="H", wind="wind",
                 zr=None, zh=None, d=None, frac_d=0.7, z0m=None, frac_z0m=None, estimate_z0m=True,
                 stab_correction=True, stab_formulation="Dyer_1970", constants=None):

    if constants is None:
        constants = bigleaf_constants()

    if stab_formulation not in ["Dyer_1970", "Businger_1971"]:
        raise ValueError("stab_formulation must be 'Dyer_1970' or 'Businger_1971'")

    # check_input(data, [ustar])

    # Determine roughness parameters
    if d is None:
        if frac_d is None:
            raise ValueError("Either 'd' or 'frac_d' must be specified")
        d = frac_d * zh

    if z0m is None and not estimate_z0m:
        if frac_z0m is None:
            raise ValueError("Either 'z0m' or 'frac_z0m' must be specified if 'estimate_z0m' = False")
        z0m = frac_z0m * zh

    if estimate_z0m:
        if z0m is not None or frac_z0m is not None:
            print("Note that arguments 'z0m' and 'frac_z0m' are ignored if 'estimate_z0m' = True. "
                  "z0m is calculated from the logarithmic wind profile equation.")

        # check_input(data, [Tair, pressure, wind, ustar, H])

        z0m = roughness_parameters(
            method="wind_profile", zh=zh, zr=zr, d=d, data=data,
            Tair=Tair, pressure=pressure, wind=wind, ustar=ustar, H=H,
            stab_roughness=True, stab_formulation=stab_formulation,
            constants=constants
        )["z0m"]

    # Check for heights below d + z0m
    if np.any(z < (d + z0m)):
        print("Warning: function is only valid for heights above d + z0m! "
              "Wind speed for heights below d + z0m will return 0!")

    # Calculate wind speeds at given heights z
    if stab_correction:
        zeta = stability_parameter(
            data=data, Tair=Tair, pressure=pressure,
            ustar=ustar, H=H, zr=z, d=d, constants=constants
        )

        psi_m = stability_correction(zeta, formulation=stab_formulation)["psi_m"]
        wind_heights = np.maximum(0, (data[ustar] / constants["k"]) *
                                  (np.log(np.maximum(0.001, (z - d) / z0m)) - psi_m))
    else:
        wind_heights = np.maximum(0, (data[ustar] / constants["k"]) *
                                  np.log(np.maximum(0.001, (z - d) / z0m)))

    return wind_heights
