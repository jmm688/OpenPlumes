import numpy as np
import pandas as pd

from scipy.spatial.distance import pdist, squareform

#Here but not very usefull
#X, Y, Z = df['X'],df['Y'],df['depth']
def coordinate_diagnostics(X,Y,Z):   
    # calculating minimum, maximum, and span for each coordinate
    diagnostics = {
    #"X_min": X.min(),
    #"X_max": X.max(),
    "X_span": X.max() - X.min(),
    #"Y_min": Y.min(),
    #"Y_max": Y.max(),
    "Y_span":Y.max() - Y.min(),
    "Z_min": Z.min(),
    "Z_max": Z.max(),
    "Z_span": Z.max() - Z.min()}
    diagnostics= {k: float(v) for k, v in diagnostics.items()}
    return diagnostics

# Very important for profiles!!! makes the distance matrix
def pairwise_horizontal_distances(df_, X, Y):
    coordinates = df_[[X, Y]].to_numpy(dtype=float)
    condensed_distance = pdist(coordinates, metric="euclidean",)
    distance_matrix = squareform(condensed_distance)

    return condensed_distance, distance_matrix

# Is this really a diagnostics? because its used later to generate outputs like the profiles...
def nearest_neighbor_diagnostics(distance_matrix,well_ids,):
    distances = distance_matrix.copy()
    # Exclude each well's distance to itself.
    np.fill_diagonal(distances, np.inf)
    # Find the column containing each row's minimum.
    nearest_indices = np.argmin(distances, axis=1)
    # Retrieve those minimum distances.
    nearest_distances = np.min(distances, axis=1)
    # Making sure positional NumPy indexing works.
    well_ids = np.asarray(well_ids)
    results = pd.DataFrame(
        {
            "Well_ID": well_ids,
            "Nearest_Well": well_ids[nearest_indices],
            "Distance_m": nearest_distances,
        }
    )

    return results