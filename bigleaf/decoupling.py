import numpy as np
import pandas as pd


def longwave_conductance(Tair, LAI, constants):
    """
    Longwave Radiative Transfer Conductance of the Canopy (Martin, 1989)

    Parameters
    ----------
    Tair : array-like or float
        Air temperature (°C)
    LAI : array-like or float
        Leaf area index (m2 m-2)
    constants : dict
        Must include:
            - 'Kelvin': conversion from °C to K
            - 'sigma': Stefan-Boltzmann constant (W m-2 K-4)
            - 'cp': specific heat of air for constant pressure (J K-1 kg-1)

    Returns
    -------
    Gr : ndarray or float
        Longwave radiative transfer conductance of the canopy (m s-1)
    """
    Tair_K = np.asarray(Tair) + constants['Kelvin']
    Gr = 4 * constants['sigma'] * Tair_K**3 * LAI / constants['cp']
    return Gr


def decoupling(data, Tair_col="Tair", pressure_col="pressure", Ga_col="Ga_h", Gs_col="Gs_ms",
               approach="Jarvis&McNaughton_1986", LAI=None, Esat_formula="Sonntag_1990",
               constants=None):
    """
    Canopy-Atmosphere Decoupling Coefficient 'Omega'

    Parameters
    ----------
    data : pd.DataFrame
        Input data containing columns for Tair, pressure, Ga, Gs
    Tair_col, pressure_col, Ga_col, Gs_col : str
        Column names for respective variables
    approach : str
        Either 'Jarvis&McNaughton_1986' or 'Martin_1989'
    LAI : float or array-like, optional
        Leaf area index (m2 m-2), required if using 'Martin_1989'
    Esat_formula : str
        Method to compute esat slope; passed to Esat_slope
    constants : dict
        Dictionary with keys: Kelvin, cp, eps, sigma, Pa2kPa

    Returns
    -------
    Omega : np.ndarray
        Decoupling coefficient
    """
    # Input checks
    if constants is None:
        raise ValueError("constants dictionary must be provided.")

    required_cols = [Tair_col, pressure_col, Ga_col, Gs_col]
    if not all(col in data.columns for col in required_cols):
        raise ValueError("Missing required columns in input data.")

    # Retrieve input data
    Tair = data[Tair_col].values
    pressure = data[pressure_col].values
    Ga = data[Ga_col].values
    Gs = data[Gs_col].values

    # Calculate slope of saturation vapor pressure (Pa/K)
    Delta = Esat_slope(Tair, formula=Esat_formula, constants=constants)  # returns a vector
    gamma = psychrometric_constant(Tair, pressure, constants)            # returns a vector
    epsilon = Delta / gamma

    if approach == "Jarvis&McNaughton_1986":
        Omega = (epsilon + 1) / (epsilon + 1 + Ga / Gs)

    elif approach == "Martin_1989":
        if LAI is None:
            raise ValueError("LAI must be provided for Martin_1989 approach.")

        Gr = longwave_conductance(Tair, LAI, constants)
        Omega = (epsilon + 1 + Gr / Ga) / (epsilon + (1 + Ga / Gs) * (1 + Gr / Ga))

    else:
        raise ValueError("Invalid approach selected.")

    return Omega