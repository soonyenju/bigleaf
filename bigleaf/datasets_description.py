
"""
AT_Neu_Jul_2010: Half-hourly eddy covariance data from Neustift (Austria)

Description
-----------
Half-hourly eddy covariance measurements from the FLUXNET site AT-Neu, a mountain meadow in Austria.
Data are from July 2010.

Source: https://sites.fluxdata.org/AT-Neu/

Format
------
DataFrame with 1488 rows and 31 columns:

- year : int
    Year of measurement
- month : int
    Month of measurement
- doy : int
    Day of year
- hour : float
    Hour (0 to 23.5)
- Tair : float
    Air temperature (°C) [TA_F]
- Tair_qc : int
    Quality control flag for Tair [TA_F_QC]
- PPFD : float
    Photosynthetic photon flux density (μmol m⁻² s⁻¹) [PPFD_IN]
- PPFD_qc : int
    QC flag for PPFD [PPFD_IN_QC]
- VPD : float
    Vapor pressure deficit (kPa) [VPD_F]
- VPD_qc : int
    QC flag for VPD [VPD_F_QC]
- pressure : float
    Atmospheric pressure (kPa) [PA_F]
- precip : float
    Precipitation (mm) [P_F]
- precip_qc : int
    QC flag for precip [P_F_QC]
- ustar : float
    Friction velocity (m s⁻¹) [USTAR]
- wind : float
    Horizontal wind speed (m s⁻¹) [WS_F]
- wind_qc : int
    QC flag for wind [WS_F_QC]
- Ca : float
    CO₂ concentration (ppm) [CO2_F_MDS]
- Ca_qc : int
    QC flag for Ca [CO2_F_MDS_QC]
- LW_up : float
    Upward longwave radiation (W m⁻²) [LW_OUT]
- Rn : float
    Net radiation (W m⁻²) [NETRAD]
- LE : float
    Latent heat flux (W m⁻²) [LE_F_MDS]
- LE_qc : int
    QC flag for LE [LE_F_MDS_QC]
- H : float
    Sensible heat flux (W m⁻²) [H_F_MDS]
- H_qc : int
    QC flag for H [H_F_MDS_QC]
- G : float
    Ground heat flux (W m⁻²) [G_F_MDS]
- G_qc : int
    QC flag for G [G_F_MDS_QC]
- NEE : float
    Net ecosystem exchange (μmol m⁻² s⁻¹) [NEE_VUT_USTAR50]
- NEE_qc : int
    QC flag for NEE [NEE_VUT_USTAR50_QC]
- GPP : float
    Gross primary productivity (μmol m⁻² s⁻¹) [GPP_NT_VUT_USTAR50]
- GPP_qc : int
    QC flag for GPP [NEE_VUT_USTAR50_QC]
- Reco : float
    Ecosystem respiration (μmol m⁻² s⁻¹) [RECO_NT_VUT_USTAR50]

Notes
-----
Some units have been converted (e.g., VPD from hPa to kPa).
Original variable names from the FLUXNET2015 dataset are provided in square brackets.

Source
------
https://fluxnet.org/ (accessed 09 November 2016)
"""

"""
DE_Tha_Jun_2014: Half-hourly eddy covariance data from Tharandt (Germany)

Description
-----------
Measurements from DE-Tha, a spruce forest site in Eastern Germany. Data collected in June 2014.

Source: https://sites.fluxdata.org/DE-Tha/

Format
------
DataFrame with 1440 rows and 32 columns, similar to AT_Neu_Jul_2010, with the addition of:

- LW_down : float
    Downward longwave radiation (W m⁻²) [LW_IN_F]

Notes
-----
Units have been converted where necessary. Original FLUXNET2015 variable names are included in square brackets.

Source
------
https://fluxnet.org/ (accessed 09 November 2016)
"""

"""
FR_Pue_May_2012: Half-hourly eddy covariance data from Puechabon (France)

Description
-----------
Data from FR-Pue, a Mediterranean evergreen oak forest in Southern France. Measurements are from May 2012.

Source: https://sites.fluxdata.org/FR-Pue/

Format
------
DataFrame with 1488 rows and 29 columns.

Structure is mostly identical to AT_Neu_Jul_2010, but this dataset lacks `G`, `G_qc`, and `LW_down`.

Notes
-----
Units have been converted from original FLUXNET2015 values where necessary.

Source
------
https://fluxnet.org/ (accessed 09 November 2016)
"""