import numpy as np

def get_k_profile_points(well_a, nearest_k, nearest_dist, points):
    t = np.linspace(0, 1, points)
    export = {}
    well_data_a = df.loc[df["Well_ID"] == well_a]
    
    xa = well_data_a["X"].to_numpy()[0]
    ya = well_data_a["Y"].to_numpy()[0]
    za = -well_data_a["depth"].to_numpy()[0]

    # getting starting point (well of choice)
    start_point = np.array([xa,ya,za,])
    
    for well, dist in zip(nearest_k, nearest_dist):
        # getting end point (nearest neighbor)
        well_data_b = df.loc[df["Well_ID"] == well]
    
        xb = well_data_b["X"].to_numpy()[0]
        yb = well_data_b["Y"].to_numpy()[0]
        zb = -well_data_b["depth"].to_numpy()[0]
    
        end_point = np.array([xb,yb,zb,])
    
        profile_points = (
            start_point
            + t[:, None] * (end_point - start_point)
        )
        export[well] = profile_points
    return export

def predict_profile_points(model, profile_dict, neighbors):
    scalar_dict = {}
    for neighbor in neighbors:
        print(neighbor)
        predict = model(profile_dict[neighbor])
        scalar_dict[neighbor] = predict
    
    return scalar_dict

def generate_nearest_neighbor_profiles(well, k_neighbors=1, points_per_profile=50):
    well_of_choice = well
    ordered_distances = ( 
        distance_df.loc[well_of_choice].drop(well_of_choice) # getting distances to well of choice and removing itself
        .sort_values() # sorting by numerical values so we have order of nearest neighbors
                    )
    # extracting the names of those wells
    well_names = ordered_distances.index.to_numpy()
    distances = ordered_distances.to_numpy()

    #retaining only the selected number of nearest neabohrs
    nearest = ordered_distances.head(k_neighbors)
    neighbor_names = nearest.index.to_numpy()
    neighbor_distances = nearest.to_numpy()

    profile_points =  get_k_profile_points(well_of_choice,
                                           neighbor_names,
                                           neighbor_distances,
                                           points_per_profile)

    profile_predict = predict_profile_points(model_fit, profile_points, neighbor_names)

    return profile_points, profile_predict