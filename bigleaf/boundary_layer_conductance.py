import numpy as np
import pandas as pd
from .surface_roughness import wind_profile, reynolds_number
from .meteorological_variables import kinematic_viscosity

def Gb_Thom(ustar, Sc=None, Sc_name=None, constants=None):
    Rb_h = 6.2 * ustar ** -0.667
    Gb_h = 1 / Rb_h
    kB_h = Rb_h * constants['k'] * ustar

    if Sc is not None or Sc_name is not None:
        if len(Sc) != len(Sc_name):
            raise ValueError("Arguments 'Sc' and 'Sc_name' must have the same length")
        if not isinstance(Sc, (list, np.ndarray)):
            raise TypeError("Argument 'Sc' must be numeric")

    Sc_full = [constants['Sc_CO2']] + (Sc if Sc else [])
    Gb_x = {f"Gb_{name}": Gb_h / ((s / constants['Pr']) ** 0.67)
            for s, name in zip(Sc_full, ['CO2'] + (Sc_name if Sc_name else []))}

    return pd.DataFrame({'Gb_h': Gb_h, 'Rb_h': Rb_h, 'kB_h': kB_h, **Gb_x})


def Gb_Choudhury(data, leafwidth, LAI, zh, zr, d,
                 Tair='Tair', pressure='pressure', wind='wind', ustar='ustar', H='H',
                 z0m=None, stab_formulation='Dyer_1970',
                 Sc=None, Sc_name=None, constants=None):

    alpha = 4.39 - 3.97 * np.exp(-0.258 * LAI)
    estimate_z0m = z0m is None

    wind_zh = wind_profile(data, z=zh, Tair=Tair, pressure=pressure, ustar=ustar, H=H,
                           zr=zr, estimate_z0m=estimate_z0m, zh=zh, d=d, z0m=z0m,
                           frac_z0m=None, stab_correction=True, stab_formulation=stab_formulation)

    wind_zh = np.maximum(0.01, wind_zh)

    if Sc is not None or Sc_name is not None:
        if len(Sc) != len(Sc_name):
            raise ValueError("Arguments 'Sc' and 'Sc_name' must have the same length")
        if not isinstance(Sc, (list, np.ndarray)):
            raise TypeError("Argument 'Sc' must be numeric")

    Gb_h = LAI * ((0.02 / alpha) * np.sqrt(wind_zh / leafwidth) * (1 - np.exp(-alpha / 2)))
    Rb_h = 1 / Gb_h
    kB_h = Rb_h * constants['k'] * data[ustar]

    Sc_full = [constants['Sc_CO2']] + (Sc if Sc else [])
    Gb_x = {f"Gb_{name}": Gb_h / ((s / constants['Pr']) ** 0.67)
            for s, name in zip(Sc_full, ['CO2'] + (Sc_name if Sc_name else []))}

    return pd.DataFrame({'Gb_h': Gb_h, 'Rb_h': Rb_h, 'kB_h': kB_h, **Gb_x})


def Gb_Su(data, zh, zr, d, Dl,
          Tair='Tair', pressure='pressure', ustar='ustar', wind='wind', H='H',
          z0m=None, fc=None, LAI=None, N=2, Cd=0.2, hs=0.01,
          stab_formulation='Dyer_1970', Sc=None, Sc_name=None, constants=None):

    if fc is None:
        if LAI is None:
            raise ValueError("One of 'fc' or 'LAI' must be provided")
        fc = 1 - np.exp(-LAI / 2)

    estimate_z0m = z0m is None

    wind_zh = wind_profile(data, z=zh, Tair=Tair, pressure=pressure, ustar=ustar, H=H,
                           zr=zr, estimate_z0m=estimate_z0m, zh=zh, d=d, z0m=z0m,
                           frac_z0m=None, stab_correction=True, stab_formulation=stab_formulation)

    v = kinematic_viscosity(data[Tair], data[pressure], constants)
    Re = reynolds_number(data[Tair], data[pressure], data[ustar], hs, constants)
    kBs = 2.46 * (Re ** 0.25) - np.log(7.4)
    Reh = Dl * wind_zh / v
    Ct = constants['Pr'] ** -0.6667 * Reh ** -0.5 * N

    kB_h = (constants['k'] * Cd) / (4 * Ct * data[ustar] / wind_zh) * fc ** 2 + kBs * (1 - fc) ** 2
    Rb_h = kB_h / (constants['k'] * data[ustar])
    Gb_h = 1 / Rb_h

    if Sc is not None or Sc_name is not None:
        if len(Sc) != len(Sc_name):
            raise ValueError("Arguments 'Sc' and 'Sc_name' must have the same length")
        if not isinstance(Sc, (list, np.ndarray)):
            raise TypeError("Argument 'Sc' must be numeric")

    Sc_full = [constants['Sc_CO2']] + (Sc if Sc else [])
    Gb_x = {f"Gb_{name}": Gb_h / ((s / constants['Pr']) ** 0.67)
            for s, name in zip(Sc_full, ['CO2'] + (Sc_name if Sc_name else []))}

    return pd.DataFrame({'Gb_h': Gb_h, 'Rb_h': Rb_h, 'kB_h': kB_h, **Gb_x})


def roughness_length_heat(z0m, kB_h):
    return z0m / np.exp(kB_h)