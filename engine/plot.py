import pyvista as pv

def generate_well_plot(plotter, pdata, cmap,color_limits,grid,contaminant_of_concern):
    # Add observed monitoring wells as spheres.
    sphere = pv.Sphere(radius=2, phi_resolution=10, theta_resolution=10)
    pc = pdata.glyph(scale=False, geom=sphere, orient=False)
    plotter.add_mesh(pc, scalars=contaminant_of_concern, cmap="turbo",
                          clim=color_limits, show_scalar_bar=False,
)
    return plotter

def generate_isosurfaces_plot(plotter,isosurface_list,color_limits,cmap,contaminant_of_concern,grid):
    iso = grid.contour(isosurfaces=isosurface_list, scalars=contaminant_of_concern)
    plotter.add_mesh(iso, scalars=contaminant_of_concern, cmap=cmap,
                         clim=color_limits, show_edges=True,
                         opacity=0.70, scalar_bar_args={"title": contaminant_of_concern, "vertical": True, "title_font_size": 12, "label_font_size": 10},
)
    return plotter

def add_profiles_to_plot(plot, dict_profile_coordinates, dict_profile_predictions,linewidth=0.05):
    
    for neighbor in dict_profile_predictions.keys():
        profile_plotter = pv.Plotter()
        
        points = dict_profile_coordinates[neighbor]
        predictions = dict_profile_predictions[neighbor]

        line = pv.MultipleLines(points)
        line['PCB'] = predictions # ***NEEDS TO BE REPLACED WITH A MORE GLOBAL VARIABLE! ***
        
        tube = line.tube(
            radius=linewidth,
            radius_factor=1.0,  # Prevent scalar-dependent thickness
        )
        
        plot.add_mesh(tube, scalars='PCB', cmap='turbo')

    return plot

def generate_scene(
    wells=None,
    grid=None,
    contaminant_of_concern=None,
    isosurface_list=None,
    cmap="turbo",
    color_limits=None,
    surface_map=None,
):
    # At least one layer must be requested.
    if wells is None and isosurface_list is None:
        raise ValueError(
            "Nothing was provided to visualize.\n"
            "Provide wells, isosurface_list, or both."
        )
    if contaminant_of_concern is None:
        raise ValueError(
            "A contaminant name must be provided.\n"
            "Example: contaminant_of_concern='PCB'"
        )
    # Isosurfaces require a prediction grid.
    if isosurface_list is not None and grid is None:
        raise ValueError(
            "A prediction grid is required to generate isosurfaces."
        )
    if (isosurface_list is not None and not isinstance(isosurface_list, list)):
        raise TypeError(
            "isosurface_list must be a list.\n"
            "Examples: [100] or [100, 500, 1000]"
        )
    if isosurface_list == []:
        raise ValueError(
            "isosurface_list cannot be empty."
        )
    plotter = pv.Plotter()

    if wells is not None:
        plotter = generate_well_plot(
            plotter=plotter,
            pdata=wells,
            grid=grid,
            contaminant_of_concern=contaminant_of_concern,
            cmap=cmap,
            color_limits=color_limits,
        )

    if isosurface_list is not None:
        plotter = generate_isosurfaces_plot(
            plotter=plotter,
            grid=grid,
            isosurface_list=isosurface_list,
            contaminant_of_concern=contaminant_of_concern,
            cmap=cmap,
            color_limits=color_limits,
        )

    if surface_map is not None:
        # TODO after v0.1: add optional surface-map workflow.
        pass

    return plotter

   