import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import Binarizer

def optimum_temperature(data, GPP="GPP", Tair="Tair", BLine=0.9, Obs_filter=30):
    # Check input
    if GPP not in data.columns or Tair not in data.columns:
        raise ValueError(f"Columns {GPP} and {Tair} must be present in the data")

    # Round to 1°C temperature bins
    Tair_bin = np.floor(data[Tair] + np.sign(data[Tair]) * 0.5).astype(int)

    # Get boundary line using quantiles
    df_bl = data.groupby(Tair_bin)[GPP].quantile(BLine).reset_index(name='GPP_Bline')

    # Get the number of observations in each temperature bin
    n_obs = data.groupby(Tair_bin).size().reset_index(name='n_obs')

    # Merge the data frames
    df_bl = pd.merge(df_bl, n_obs, on=Tair_bin)

    # Remove temperature bins with n_obs below the threshold
    df_bl = df_bl[df_bl['n_obs'] >= Obs_filter]

    # Get the smoothed boundary line using loess (locally weighted scatterplot smoothing)
    from sklearn.preprocessing import PolynomialFeatures
    from sklearn.linear_model import LinearRegression

    poly = PolynomialFeatures(degree=2)
    X_poly = poly.fit_transform(df_bl[[Tair_bin]])
    model = LinearRegression().fit(X_poly, df_bl['GPP_Bline'])

    # Predict smoothed values
    df_bl['GPP_Bline_smooth'] = model.predict(X_poly)

    # Find the thermal optimum (Topt)
    df_bl_sorted = df_bl.sort_values(by='GPP_Bline_smooth', ascending=False)
    Topt = df_bl_sorted.iloc[0]['Tair_bin']
    GPP_bl = df_bl_sorted.iloc[0]['GPP_Bline']

    opt_temp = {"Topt": Topt, "GPP_bl": GPP_bl}

    # Return results
    optimum_temp = {"df_bl": df_bl, "opt_temp": opt_temp}
    return optimum_temp
