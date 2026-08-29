import numpy as np
import pyvista as pv

def generate_grid(boundary, xyz_resolution):
    # making a numpy array using min and max from model boundaries
    #mins
    domain_min = np.array([
    round(boundary["x_min"]),
    round(boundary["y_min"]),
    round(boundary["z_min"]),
    ])
    # maxs
    domain_max = np.array([
        round(boundary["x_max"]),
        round(boundary["y_max"]),
        round(boundary["z_max"]),
    ])
    
    domain_span = domain_max - domain_min
    
    number_of_cells = np.ceil(domain_span / xyz_resolution).astype(int)
    
    # Include both endpoints.
    number_of_points = number_of_cells + 1
    
    # Adjust spacing so the grid ends exactly at domain_max.
    actual_spacing = domain_span / number_of_cells

    grid = pv.ImageData()

    ### diagnostics!!!!
    grid.dimensions = tuple(number_of_points)
    grid.origin = tuple(domain_min)
    grid.spacing = tuple(actual_spacing)

    return grid