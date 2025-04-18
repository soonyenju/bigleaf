import numpy as np
import pandas as pd

# Constants (example values, adjust according to your constants)
def bigleaf_constants():
    return {
        'Kelvin': 273.15,  # Celsius to Kelvin conversion
        'cp': 1005,        # Specific heat of air at constant pressure (J/kg·K)
        'k': 0.4,          # von Karman constant
        'g': 9.81          # Gravitational acceleration (m/s^2)
    }

# Helper function to check if the necessary input columns are in the data
def check_input(data, required_columns):
    missing_columns = [col for col in required_columns if col not in data.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

# Air density function
def air_density(Tair, pressure, constants):
    # Ideal gas law: rho = p / (R * T)
    R = 287.05  # Specific gas constant for dry air in J/(kg·K)
    return pressure * 1e3 / (R * Tair)  # converting pressure from kPa to Pa

# Monin-Obukhov Length function
def monin_obukhov_length(data, Tair="Tair", pressure="pressure", ustar="ustar", H="H", constants=None):
    if constants is None:
        constants = bigleaf_constants()

    check_input(data, [Tair, pressure, ustar, H])

    rho = air_density(data[Tair], data[pressure], constants)
    Tair_K = data[Tair] + constants['Kelvin']  # Celsius to Kelvin
    MOL = - (rho * constants['cp'] * data[ustar]**3 * Tair_K) / (constants['k'] * constants['g'] * data[H])

    return MOL

# Stability Parameter "zeta" function
def stability_parameter(data, Tair="Tair", pressure="pressure", ustar="ustar", H="H", zr="zr", d="d", constants=None):
    if constants is None:
        constants = bigleaf_constants()

    check_input(data, [Tair, pressure, ustar, H, zr, d])

    MOL = monin_obukhov_length(data, Tair, pressure, ustar, H, constants)
    zeta = (data[zr] - data[d]) / MOL

    return zeta

# Stability Correction function for Heat and Momentum
def stability_correction(zeta, formulation="Dyer_1970"):
    if formulation not in ["Dyer_1970", "Businger_1971"]:
        raise ValueError("Formulation must be either 'Dyer_1970' or 'Businger_1971'")

    psi_h = np.full_like(zeta, np.nan)
    psi_m = np.full_like(zeta, np.nan)

    if formulation == "Businger_1971":
        x_h = -7.8
        x_m = -6
        y_h = 0.95 * (1 - 11.6 * zeta) ** 0.5
        y_m = (1 - 19.3 * zeta) ** 0.25
    elif formulation == "Dyer_1970":
        x_h = x_m = -5
        y_h = (1 - 16 * zeta) ** 0.5
        y_m = (1 - 16 * zeta) ** 0.25

    # Stable conditions
    stable = zeta >= 0
    psi_h[stable] = x_h * zeta[stable]
    psi_m[stable] = x_m * zeta[stable]

    # Unstable conditions
    unstable = zeta < 0
    psi_h[unstable] = 2 * np.log((1 + y_h[unstable]) / 2)
    psi_m[unstable] = (2 * np.log((1 + y_m[unstable]) / 2) +
                       np.log((1 + y_m[unstable]**2) / 2) -
                       2 * np.arctan(y_m[unstable]) + np.pi / 2)

    return pd.DataFrame({'psi_h': psi_h, 'psi_m': psi_m})

# # Example usage:
# data = pd.DataFrame({
#     'Tair': [25],
#     'pressure': [100],
#     'ustar': [0.3],
#     'H': [100],
#     'zr': [40],
#     'd': [15]
# })

# zeta_values = stability_parameter(data, zr="zr", d="d")
# corrections = stability_correction(zeta_values)
# print(corrections)