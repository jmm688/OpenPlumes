def create_well_points(df,X ,Y ,Z):
    X, Y, Z = df['X'],df['Y'],df['depth']

    well_points = np.column_stack(
        (
            X.to_numpy(),
            Y.to_numpy(),
            -Z.to_numpy())
    )
    return well_points

def add_contaminants(wells, data, list_of_contaminants):
    if len(data) != wells.n_points:
        raise ValueError(
            "The number of DataFrame rows does not match "
            "the number of PyVista points."
        )
    for contaminant in list_of_contaminants:
        if contaminant not in data.columns:
            raise ValueError(
                f"Contaminant '{contaminant}' is not in the original data."
            )
        wells[contaminant] = data[contaminant].to_numpy()
    return wells