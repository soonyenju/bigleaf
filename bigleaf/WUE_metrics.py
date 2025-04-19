import numpy as np
import pandas as pd
from .bigleaf_constants import *
from .unit_conversions import LE_to_ET

def WUE_metrics(data, GPP_col="GPP", NEE_col="NEE", LE_col="LE", VPD_col="VPD", Tair_col="Tair",
                constants=None):
    """
    Calculate WUE, WUE_NEE, IWUE, and uWUE metrics from flux data.

    Parameters:
    - data: pandas DataFrame with flux data
    - GPP_col, NEE_col, LE_col, VPD_col, Tair_col: column names in data
    - constants: dictionary of unit conversion constants

    Returns:
    - Dictionary of median WUE metrics
    """
    if constants is None:
        constants = bigleaf_constants()

    # Extract relevant columns
    GPP = data[GPP_col]
    NEE = data[NEE_col]
    LE = data[LE_col]
    VPD = data[VPD_col]
    Tair = data[Tair_col]

    # Convert LE to ET
    ET = LE_to_ET(LE, Tair)

    # Convert GPP and NEE to gC m-2 s-1
    GPP_gC = GPP * constants["umol2mol"] * constants["Cmol"] * constants["kg2g"]
    NEE_gC = NEE * constants["umol2mol"] * constants["Cmol"] * constants["kg2g"]

    # Avoid divide-by-zero or NaN by masking invalid ET values
    valid_mask = (ET > 0) & GPP_gC.notna() & NEE_gC.notna() & VPD.notna()

    # Compute metrics
    WUE = np.median((GPP_gC[valid_mask] / ET[valid_mask]).dropna())
    WUE_NEE = np.median((np.abs(NEE_gC[valid_mask]) / ET[valid_mask]).dropna())
    IWUE = np.median(((GPP_gC[valid_mask] * VPD[valid_mask]) / ET[valid_mask]).dropna())
    uWUE = np.median(((GPP_gC[valid_mask] * np.sqrt(VPD[valid_mask])) / ET[valid_mask]).dropna())

    return {
        "WUE": WUE,
        "WUE_NEE": WUE_NEE,
        "IWUE": IWUE,
        "uWUE": uWUE
    }