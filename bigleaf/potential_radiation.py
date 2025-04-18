import numpy as np

def extraterrestrial_radiation(doy, solar_constant=1361):
    """
    Compute the extraterrestrial solar radiation with eccentricity correction.

    Parameters:
    doy (array-like): Day of year (DoY)
    solar_constant (float): Solar constant (W m^-2)

    Returns:
    np.ndarray: Extraterrestrial radiation (W m^-2)
    """
    # Fractional year in radians
    FracYearRad = 2 * np.pi * (doy - 1) / 365.24

    # Eccentricity correction
    ExtRadiation = solar_constant * (
        1.00011 + 0.034221 * np.cos(FracYearRad) + 0.00128 * np.sin(FracYearRad)
        + 0.000719 * np.cos(2 * FracYearRad) + 0.000077 * np.sin(2 * FracYearRad)
    )

    return ExtRadiation


def potential_radiation(doy, hour, latDeg, longDeg, timezone, useSolartime=True):
    """
    Compute potential radiation for a given geolocation and day of year.

    Parameters:
    doy (array-like): Day of year (start at 1)
    hour (array-like): Daytime as decimal hour of local time zone
    latDeg (float): Latitude (decimal degrees)
    longDeg (float): Longitude (decimal degrees)
    timezone (int): Time zone (hours)
    useSolartime (bool): If True, correct for solar time. Otherwise, use local winter time

    Returns:
    np.ndarray: Potential radiation (W m^-2)
    """
    # Assuming you have a function `compute_sun_position_doy_hour` to get solar elevation
    solElevRad = compute_sun_position_doy_hour(doy, hour, latDeg, longDeg, timezone, useSolartime)["elevation"]

    # Compute extraterrestrial radiation
    extRadiation = extraterrestrial_radiation(doy)

    # Calculate potential radiation
    potRad = np.where(solElevRad <= 0, 0, extRadiation * np.sin(solElevRad))

    return potRad

# # Helper function: this function should be defined based on your solar time computation
# def compute_sun_position_doy_hour(doy, hour, latDeg, longDeg, timezone, isCorrectSolartime=True):
#     # Placeholder for actual implementation of sun position computation
#     # This should return a DataFrame or structured array with 'elevation' data
#     pass
