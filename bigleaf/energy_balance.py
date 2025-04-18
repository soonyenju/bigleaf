import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

def biochemical_energy(NEE, alpha=0.422):
    """
    Radiant energy absorbed in photosynthesis or heat release by respiration from NEE.
    """
    return alpha * -np.array(NEE)

def energy_use_efficiency(GPP, Rn, alpha=0.422):
    """
    Fraction of net radiation fixed by primary productivity.
    """
    GPP = np.array(GPP)
    Rn = np.array(Rn)
    Sp = biochemical_energy(-GPP, alpha)
    mask = np.isfinite(Sp) & np.isfinite(Rn)
    return np.sum(Sp[mask]) / np.sum(Rn[mask])

def energy_closure(data, Rn='Rn', G=None, S=None, LE='LE', H='H',
                   instantaneous=False, missing_G_as_NA=False, missing_S_as_NA=False):
    """
    Calculates the energy balance closure from energy flux terms.
    """
    df = data.copy()
    df['Rn'] = df[Rn]
    df['LE'] = df[LE]
    df['H'] = df[H]
    df['G'] = df[G] if G else 0
    df['S'] = df[S] if S else 0

    if not missing_G_as_NA and G:
        df['G'] = df['G'].fillna(0)
    elif G is None:
        print("Ground heat flux G is not provided and set to 0.")

    if not missing_S_as_NA and S:
        df['S'] = df['S'].fillna(0)
    elif S is None:
        print("Energy storage fluxes S are not provided and set to 0.")

    if instantaneous:
        EBR = (df['LE'] + df['H']) / (df['Rn'] - df['G'] - df['S'])
        return EBR.to_numpy()

    mask = df[['Rn', 'LE', 'H']].notna().all(axis=1)
    if G: mask &= df[[G]].notna().all(axis=1)
    if S: mask &= df[[S]].notna().all(axis=1)
    df = df[mask]
    n = len(df)

    sum_LE_H = np.sum(df['LE'] + df['H'])
    sum_Rn_G_S = np.sum(df['Rn'] - df['G'] - df['S'])
    EBR = sum_LE_H / sum_Rn_G_S

    X = (df['Rn'] - df['G'] - df['S']).values.reshape(-1, 1)
    y = (df['LE'] + df['H']).values
    model = LinearRegression().fit(X, y)

    return {
        "n": n,
        "intercept": round(model.intercept_, 3),
        "slope": round(model.coef_[0], 3),
        "r^2": round(model.score(X, y), 3),
        "EBR": round(EBR, 3)
    }

def isothermal_Rn(Rn, Tair, Tsurf, emissivity, constants=None):
    """
    Calculates the isothermal net radiation assuming surface and air at same temperature.
    """
    if constants is None:
        constants = {"sigma": 5.670374419e-8, "Kelvin": 273.15}
    Tair_K = np.array(Tair) + constants["Kelvin"]
    Tsurf_K = np.array(Tsurf) + constants["Kelvin"]
    return np.array(Rn) + emissivity * constants["sigma"] * (Tsurf_K**4 - Tair_K**4)