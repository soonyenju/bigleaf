import numpy as np
import pandas as pd
from .bigleaf_constants import *

def air_density(Tair, pressure, constants = None):
    """Calculate air density."""
    if constants is None:
        constants = bigleaf_constants()
    Tair_K = Tair + constants['Kelvin']
    pressure_Pa = pressure * constants['kPa2Pa']
    return pressure_Pa / (constants['Rd'] * Tair_K)

def pressure_from_elevation(elev, Tair, VPD=None, constants=None, virtual_temp_func=None):
    """Calculate pressure from elevation."""
    if constants is None:
        constants = bigleaf_constants()
    Tair_K = Tair + constants['Kelvin']

    if VPD is None:
        pressure = constants['pressure0'] / np.exp(constants['g'] * elev / (constants['Rd'] * Tair_K))
    else:
        pressure1 = constants['pressure0'] / np.exp(constants['g'] * elev / (constants['Rd'] * Tair_K))
        Tv = virtual_temp_func(Tair, pressure1 * constants['Pa2kPa'], VPD, "Sonntag_1990", constants)
        Tv_K = Tv + constants['Kelvin']
        pressure = constants['pressure0'] / np.exp(constants['g'] * elev / (constants['Rd'] * Tv_K))

    return pressure * constants['Pa2kPa']

def esat_slope(Tair, formula="Sonntag_1990", constants=None):
    """Calculate saturation vapor pressure and its slope."""
    if constants is None:
        constants = bigleaf_constants()
    if formula == "Sonntag_1990":
        a, b, c = 611.2, 17.62, 243.12
    elif formula == "Alduchov_1996":
        a, b, c = 610.94, 17.625, 243.04
    elif formula == "Allen_1998":
        a, b, c = 610.8, 17.27, 237.3
    else:
        raise ValueError("Invalid formula name.")

    exp_term = np.exp((b * Tair) / (c + Tair))
    Esat = a * exp_term * constants['Pa2kPa']
    Delta = a * exp_term * (b / (c + Tair) - (b * Tair) / (c + Tair)**2) * constants['Pa2kPa']

    return pd.DataFrame({'Esat': Esat, 'Delta': Delta})

def latent_heat_vaporization(Tair):
    """Calculate latent heat of vaporization."""
    k1 = 2.501
    k2 = 0.00237
    return (k1 - k2 * Tair) * 1e6  # J/kg

def psychrometric_constant(Tair, pressure, constants = None):
    """Calculate psychrometric constant."""
    if constants is None:
        constants = bigleaf_constants()
    lambda_val = latent_heat_vaporization(Tair)
    return (constants['cp'] * pressure) / (constants['eps'] * lambda_val)

def wetbulb_solver(ea, Tair, gamma, accuracy, formula, constants = None):
    """Solve for wet bulb temperature."""
    if constants is None:
        constants = bigleaf_constants()
    def objective(Tw):
        Esat_val = esat_slope(Tw, formula, constants)['Esat']
        return np.abs(ea - (Esat_val - constants['Le067'] * gamma * (Tair - Tw)))

    result = minimize_scalar(objective, bounds=(-100, 100), method='bounded', options={'xatol': accuracy})
    return result

def wetbulb_temp(Tair, pressure, VPD, accuracy=1e-3, formula="Sonntag_1990", constants=None, vpd_to_e=None):
    """Calculate wet-bulb temperature."""
    if constants is None:
        constants = bigleaf_constants()
    if not isinstance(accuracy, (float, int)):
        raise ValueError("'accuracy' must be numeric")
    if accuracy > 1:
        print("'accuracy' is set to 1 degC")
        accuracy = 1

    ndigits = int(np.ceil(-np.log10(accuracy))) if accuracy < 1 else 0

    gamma = psychrometric_constant(Tair, pressure, constants)
    ea = vpd_to_e(VPD, Tair, formula, constants)

    Tw = []
    for i in range(len(Tair)):
        if any(np.isnan([ea[i], Tair[i], gamma[i]])) or any(np.isinf([ea[i], Tair[i], gamma[i]])):
            Tw.append(np.nan)
        else:
            res = wetbulb_solver(ea[i], Tair[i], gamma[i], accuracy, formula, constants)
            Tw.append(round(res.x, ndigits))

    return np.array(Tw)


def dew_point_solver(ea, accuracy, Esat_formula, constants = None):
    """
    Solves for the dew point temperature (Td) given actual vapor pressure (ea) by 
    minimizing the difference between ea and the saturation vapor pressure (Esat).

    Parameters:
    ----------
    ea : float
        Actual vapor pressure (Pa or hPa, depending on the Esat_slope output).
    accuracy : float
        Desired accuracy (tolerance) for the optimization (in °C).
    Esat_formula : str or callable
        The formula or method used to compute saturation vapor pressure.
    constants : dict
        A dictionary of physical constants used in the Esat_slope calculation.

    Returns:
    -------
    result : OptimizeResult
        The result object from scipy.optimize.minimize_scalar, containing the 
        optimized dew point temperature and additional metadata.
    """
    def objective(Td):
        Esat = esat_slope(Td, Esat_formula, constants)["Esat"]
        return abs(ea - Esat)
    
    result = minimize_scalar(objective, bounds=(-100, 100), method='bounded', options={'xatol': accuracy})
    return result


def dew_point(Tair, VPD, accuracy=1e-3, Esat_formula="Sonntag_1990", constants=None):
    """
    Calculates the dew point temperature from air temperature and vapor pressure deficit (VPD)
    using numerical optimization to solve for the temperature at which saturation vapor pressure 
    equals the actual vapor pressure.

    Parameters:
    ----------
    Tair : array-like
        Air temperature in degrees Celsius.
    VPD : array-like
        Vapor Pressure Deficit (same units as returned by `VPD_to_e`, typically Pa or hPa).
    accuracy : float, optional (default=1e-3)
        Desired accuracy (tolerance in °C) for dew point calculation.
        If greater than 1, it will be reset to 1 and a message will be printed.
    Esat_formula : str, optional (default="Sonntag_1990")
        The formula or method used to compute saturation vapor pressure (e.g., "Sonntag_1990").
    constants : dict, optional
        Dictionary of physical constants used in saturation vapor pressure calculations.
        If None, defaults will be used via `bigleaf_constants()`.

    Returns:
    -------
    Td : numpy.ndarray
        Dew point temperature array (°C) with same length as input `Tair`/`VPD`.
        If `ea` is missing or invalid, the corresponding Td value is NaN.
    """
    if constants is None:
        constants = bigleaf_constants()

    if not isinstance(accuracy, (int, float)):
        raise ValueError("'accuracy' must be numeric!")

    if accuracy > 1:
        print("'accuracy' is set to 1 degC")
        accuracy = 1

    # determine number of digits to print
    try:
        ndigits = int(str(accuracy).split('e-')[1])
    except IndexError:
        ndigits = 0

    ea = VPD_to_e(VPD, Tair, Esat_formula)

    Td = []
    for i in range(len(ea)):
        if ea[i] is None or isinstance(ea[i], float) and (np.isnan(ea[i]) or np.isinf(ea[i])):
            Td.append(np.nan)
        else:
            result = dew_point_solver(ea[i], accuracy=accuracy, Esat_formula=Esat_formula, constants=constants)
            Td.append(round(result['minimum'], ndigits))

    return np.array(Td)


def virtual_temp(Tair, pressure, VPD, formula, constants=None):
    """Calculate virtual temperature."""
    if constants is None:
        constants = bigleaf_constants()
    esat = esat_slope(Tair, formula, constants)['Esat']
    VPD_Pa = VPD * constants['kPa2Pa']
    return Tair * (1 + constants['eps'] * VPD_Pa / esat)

# def vpd_to_e(VPD, Tair, formula, constants=None):
#     """Convert VPD to actual vapor pressure (Pa)."""
#     if constants is None:
#         constants = bigleaf_constants()
#     esat = esat_slope(Tair, formula, constants)['Esat']
#     return VPD * esat * constants['kPa2Pa']

def kinematic_viscosity(Tair, pressure, constants=None):
    """
    Calculate kinematic viscosity of air.

    Parameters:
    - Tair: Air temperature in °C
    - pressure: Air pressure in kPa
    - constants: Dictionary of physical constants

    Returns:
    - Kinematic viscosity (m² s⁻¹)
    """
    if constants is None:
        constants = bigleaf_constants()
    Tair_K = Tair + constants["Kelvin"]
    pressure_Pa = pressure * constants["kPa2Pa"]

    v = 1.327e-05 * (constants["pressure0"] / pressure_Pa) * (Tair_K / constants["Tair0"]) ** 1.81
    return v