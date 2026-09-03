import sys
import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import geopandas as gpd

from pyproj import CRS, Transformer
from scipy.spatial.distance import pdist, squareform
import pyvista as pv
from pathlib import Path

from scipy.spatial.distance import pdist, squareform
from scipy.interpolate import RBFInterpolator
from engine.coordinates import project_coordinates
from engine.wells import create_well_points
from engine.wells import add_contaminants
from engine.diagnostics import pairwise_horizontal_distances
from engine.diagnostics import nearest_neighbor_diagnostics
from engine.boundaries import get_domain_bounds
from engine.interpolation import interpolateRBF
from engine.grids import generate_grid
from engine.profiles import generate_nearest_neighbor_profiles
from engine.plot import generate_scene

def main():
    project_root = Path(__file__).resolve().parent
    examples_directory = project_root / "Examples"
    #Defining path to file
    path = examples_directory / "new_dummy_data3.csv"
    df = pd.read_csv(path)

    # File path for the model domain which will be used later
    boundary_path = (examples_directory / "model_domain_2" / "model_domain_2.shp")

    source_crs = "EPSG:4326"
    model_crs = "EPSG:32619"
    Longitude = 'long'
    Latitude = 'lat'

    df = project_coordinates(df,Longitude,Latitude,source_crs,model_crs)

    # Defining x,y,z points using our df data and projected coordinates
    well_points = create_well_points(df, X='X', Y='Y', Z='depth',)

    # Transforming the x,y,z, array into a pvista array and assigning a COC as a scalar
    pdata = pv.PolyData(well_points)
    COC_list = ["PCB"] # COC means contaminant/s of concern
    pdata = add_contaminants(
        wells=pdata,
        data=df,
        list_of_contaminants=COC_list,
    )


    # Really just calling the distance matrix generator
    condensed_distance, distance_matrix = (pairwise_horizontal_distances(df,X="X",Y="Y",))
    distance_df = pd.DataFrame(distance_matrix,index=df["Well_ID"],columns=df["Well_ID"],) ### hmmm i can see some redundancy here. especially because distance is the only one we really need..

    # creating a df of the nearest neighbors and their distances... Pretty ysefull TBH!!
    nearest_neighbors = nearest_neighbor_diagnostics(distance_matrix, df["Well_ID"],)
    nearest_neighbors.sort_values("Distance_m") # sorting by increasing distance

    # Using a shapefile to determine the model boundary
    boundary_path = ("/mnt/c/Users/josem/OneDrive/Documents/Dummy_data/model_domain/shape_file/model_domain_2/model_domain_2.shp")
    bounds = get_domain_bounds(boundaries=boundary_path,z_bounds=(-110, -20),model_crs=model_crs,) # make sure you define mode_crs...

    # Generating an RBF model and fitting the wells to the model
    neighbors_num = None
    epsilon_num = .001
    contaminant = 'PCB'

    model_fit = interpolateRBF(xyz_points = pdata.points,
                            contaminant = pdata.point_data[contaminant],
                            kernel='thin_plate_spline',
                            neighbors_num = neighbors_num, 
                            epsilon_num = epsilon_num,
                            )

    # Generating a model grid
    requested_spacing = np.array([
                                    5.0, #X meters
                                    5.0, #Y meters
                                    5.0  #Z meters
                                ])
    grid = generate_grid(
                        boundary = bounds,
                        xyz_resolution = requested_spacing
                        )

    # making predictions at each grid point
    prediction_points = grid.points
    predicted_values = model_fit(prediction_points)

    grid.point_data[contaminant] = predicted_values # assigning scalar name to the predictions should be COC list but we shall see....

    # defining isofurfaces or "making plumes"
    #isosurfaces = grid.contour(isosurfaces=[0, 500, 1000], scalars=contaminant, # again calls contaminant but should be one global var maybe COC list.. but we shall see...
    #)

    # Evaluating and making profiles for a given well
    profile_points, profile_predictions = generate_nearest_neighbor_profiles(well = 'MW-10',
                                                                            dist_matrix = distance_df,
                                                                            data_ = df,
                                                                            k_neighbors = 5,
                                                                            points_per_profile=10,
                                                                            model_fit = model_fit,
                                                                            )

    # Use one color range for both predicted and observed PCB.
    color_limits = (
        min(grid[contaminant].min(), pdata[contaminant].min()),
        max(grid[contaminant].max(), pdata[contaminant].max()),
    )

    filename = "/mnt/c/Users/josem/OneDrive/Documents/Dummy_data/DEM.tif"


    a = generate_scene(
        wells=pdata,
        grid=grid,
        contaminant_of_concern=contaminant,
        isosurface_list=[0,1000],
        cmap="turbo",
        color_limits=color_limits,
        surface_map=filename,
    )

    a.show()


if __name__ == "__main__":
    main()