def filter_data(df,
                vars=None,
                filter_vars=True,
                quality_control=None,
                quality_threshold=0,
                precip=None,
                exclude_precip=True,
                u_star=None,
                u_star_threshold=None,
                growing_season=None,
                include_growing_season_only=True):
    """
    Filter eddy covariance data based on quality control, precipitation, u* threshold, and growing season.

    Parameters
    ----------
    df : pd.DataFrame
        Input data with all variables as columns.
    vars : list of str, optional
        Variable names to include in the filtered output. If None, all columns are returned.
    filter_vars : bool
        If True, return only the filtered variables in `vars`; else, return all variables.
    quality_control : str or list of str, optional
        Column(s) with quality control flags (0, 1, 2). Can be one or multiple columns.
    quality_threshold : int
        Maximum quality control value to allow (0 = best).
    precip : str, optional
        Name of the precipitation column.
    exclude_precip : bool
        If True, exclude time steps with precipitation.
    u_star : str, optional
        Column name of friction velocity (u*).
    u_star_threshold : float or None
        Minimum acceptable u* value.
    growing_season : str, optional
        Name of column (bool or 0/1) indicating growing season.
    include_growing_season_only : bool
        If True, retain only growing season records.

    Returns
    -------
    pd.DataFrame
        Filtered data.
    """

    data = df.copy()
    mask = pd.Series(True, index=data.index)

    # Quality control
    if quality_control is not None:
        if isinstance(quality_control, str):
            quality_control = [quality_control]
        for qc_col in quality_control:
            mask &= data[qc_col] <= quality_threshold

    # Precipitation
    if exclude_precip and precip is not None:
        mask &= data[precip] == 0

    # u* threshold
    if u_star is not None and u_star_threshold is not None:
        mask &= data[u_star] >= u_star_threshold

    # Growing season
    if growing_season is not None and include_growing_season_only:
        mask &= data[growing_season].astype(bool)

    filtered_data = data[mask]

    if vars is not None and filter_vars:
        filtered_data = filtered_data[vars]

    return filtered_data