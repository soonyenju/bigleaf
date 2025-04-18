import pandas as pd
import numpy as np

def check_length(varlist):
    """
    Test Variables for Equal Length

    Parameters
    ----------
    varlist : list
        List of variables for which the length has to be compared.

    Returns
    -------
    list
        Validated list of variables
    """
    flat_list = []
    for v in varlist:
        if isinstance(v, list):
            flat_list.extend(v)
        else:
            flat_list.append(v)

    lengths = [len(v) for v in flat_list if hasattr(v, '__len__') and not isinstance(v, str)]

    unique_lengths = list(set(lengths))

    if len(unique_lengths) >= 2:
        if sorted(unique_lengths)[0] != 1 or len(unique_lengths) > 2:
            raise ValueError("All input variables must have the same length or a length of 1!")

    return flat_list

def check_input(data=None, **kwargs):
    """
    Check Input for Functions in the bigleaf Package (Python version)

    Parameters
    ----------
    data : pd.DataFrame or np.ndarray or None
        Optional input dataset
    **kwargs : dict
        Input variables as keyword arguments

    Returns
    -------
    dict
        Dictionary of validated variables
    """
    vars_list = check_length(list(kwargs.values()))
    validated_vars = {}

    for varname, var in kwargs.items():
        if isinstance(var, str):
            if data is not None:
                if len(var) == 1:
                    if var in data.columns:
                        column = data[var]
                        if pd.api.types.is_numeric_dtype(column):
                            validated_vars[varname] = column.values
                        else:
                            raise TypeError(f"Column '{var}' representing '{varname}' in the input must be numeric.")
                    else:
                        raise ValueError(f"There is no column named '{var}' in the input data. "
                                         f"Provide a valid column name or a numeric array of correct length.")
                else:
                    raise ValueError(f"Variable name '{varname}' must have length 1.")
            else:
                raise ValueError(f"Variable '{var}' is a string and interpreted as a column name, "
                                 "but no input DataFrame was provided.")
        else:
            # Handle numeric or other directly provided input
            if var is None or (isinstance(var, float) and np.isnan(var)):
                validated_vars[varname] = var
                continue

            if data is not None:
                if isinstance(var, (list, np.ndarray, pd.Series)):
                    if len(var) == len(data):
                        validated_vars[varname] = np.array(var)
                    elif len(var) == 1:
                        validated_vars[varname] = np.repeat(var, len(data))
                    else:
                        raise ValueError(f"Variable '{varname}' must have same length as data or be of length 1.")
                else:
                    raise TypeError(f"Variable '{varname}' must be numeric.")
            else:
                if isinstance(var, (int, float, np.ndarray, list, pd.Series)):
                    validated_vars[varname] = np.array(var)
                else:
                    raise TypeError(f"Variable '{varname}' must be numeric.")

    return validated_vars