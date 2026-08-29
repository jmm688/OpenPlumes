from pathlib import Path

import geopandas as gpd
import numpy as np


# validating to make sure bounds are in the correct format
def validate_bounds(bounds):
    """Validate and return six numeric 3D domain bounds."""
    expected_format = ("[x_min, x_max, y_min, y_max, z_min, z_max]")

    if not isinstance(bounds, (list, tuple, np.ndarray)):
        raise TypeError(
            "Manual domain bounds must be supplied as a list, "
            "tuple, or NumPy array.\n"
            f"Received type: {type(bounds).__name__}\n"
            f"Expected format: {expected_format}\n"
            "Example: [445020, 445180, 4046750, 4046910, -110, -20]"
        )
    if len(bounds) != 6:
        raise ValueError(
            "Manual domain bounds require exactly six values.\n"
            f"Received {len(bounds)} values: {bounds}\n"
            f"Expected format: {expected_format}"
        )
    invalid_values = [value for value in bounds if not isinstance(value, (int, float, np.number))]

    if invalid_values:
        raise TypeError(
            "Every domain-bound value must be numeric.\n"
            f"Invalid values: {invalid_values}\n"
            f"Received bounds: {bounds}\n"
            f"Expected format: {expected_format}"
        )

    bounds = np.asarray(bounds, dtype=float)

    if not np.isfinite(bounds).all():
        invalid_positions = np.flatnonzero(
            ~np.isfinite(bounds)
        ).tolist()
        raise ValueError(
            "Domain bounds contain NaN or infinite values.\n"
            f"Invalid positions: {invalid_positions}\n"
            f"Received bounds: {bounds.tolist()}"
        )
    x_min, x_max, y_min, y_max, z_min, z_max = bounds

    if x_min >= x_max:
        raise ValueError(
            "Invalid X bounds.\n"
            f"Received x_min={x_min:g} and x_max={x_max:g}.\n"
            "x_min must be smaller than x_max.\n"
            f"Expected format: {expected_format}"
        )
    if y_min >= y_max:
        raise ValueError(
            "Invalid Y bounds.\n"
            f"Received y_min={y_min:g} and y_max={y_max:g}.\n"
            "y_min must be smaller than y_max.\n"
            f"Expected format: {expected_format}"
        )
    if z_min >= z_max:
        raise ValueError(
            "Invalid Z bounds.\n"
            f"Received z_min={z_min:g} and z_max={z_max:g}.\n"
            "z_min must be smaller than z_max.\n"
            "If depth is represented as negative Z, a valid example "
            "is z_min=-110 and z_max=-20."
        )

    return {
        "x_min": float(x_min),
        "x_max": float(x_max),
        "y_min": float(y_min),
        "y_max": float(y_max),
        "z_min": float(z_min),
        "z_max": float(z_max),
    }


def domain_bounds_from_values(boundaries):
    """Create domain bounds from six manually supplied values."""
    return validate_bounds(boundaries)

def domain_bounds_from_file(boundaries_path,z_bounds,model_crs,):
    """Read horizontal bounds from a polygon boundary file."""

    path = Path(boundaries_path)
    accepted_formats = {".shp"}

    if not path.is_file():
        raise FileNotFoundError(
            "The model-boundary file could not be found.\n"
            f"Received path: {path}\n"
            "Check the spelling, directory, and file extension."
        )
    if path.suffix.lower() not in accepted_formats:
        raise ValueError(
            "The model-boundary file format is not currently supported.\n"
            f"Received extension: {path.suffix or '[no extension]'}\n"
            f"Supported extensions: {sorted(accepted_formats)}\n"
            "Convert the boundary to a supported format and try again."
        )
    if z_bounds is None:
        raise ValueError(
            "Vertical bounds are required when loading a horizontal "
            "boundary file.\n"
            "Expected format: z_bounds=(z_min, z_max)\n"
            "Example: z_bounds=(-110, -20)"
        )
    try:
        boundary_data = gpd.read_file(path)
        
    except Exception as error:
        raise RuntimeError(
            "The model-boundary file exists but could not be read.\n"
            f"File: {path}\n"
            "For a shapefile, confirm that its companion files such as "
            ".dbf and .shx are present in the same directory."
        ) from error

    if boundary_data.empty:
        raise ValueError(
            "The model-boundary file contains no features.\n"
            f"File: {path}\n"
            "Add at least one polygon feature before using this file."
        )

    if boundary_data.crs is None:
        raise ValueError(
            "The model-boundary file does not declare a CRS.\n"
            f"File: {path}\n"
            "OpenPlumes cannot safely combine this boundary with projected "
            "well coordinates until its CRS is known."
        )

    valid_geometry_types = {"Polygon", "MultiPolygon"}

    invalid_geometry_types = sorted(
        set(boundary_data.geom_type)
        - valid_geometry_types
    )
    if invalid_geometry_types:
        raise ValueError(
            "The model-boundary file must contain polygon geometry.\n"
            f"Invalid geometry types found: {invalid_geometry_types}\n"
            "Use a Polygon or MultiPolygon boundary."
        )
    if not boundary_data.geometry.is_valid.all():
        invalid_count = int((~boundary_data.geometry.is_valid).sum())

        raise ValueError(
            "The model-boundary file contains invalid geometry.\n"
            f"Number of invalid features: {invalid_count}\n"
            "Repair the polygon geometry before creating the model grid."
        )
    if model_crs is None:
        raise ValueError(
            "A projected modeling CRS is required when loading a "
            "boundary file.\n"
            "Example: model_crs='EPSG:32613'"
        )
    # Transform the boundary into the same CRS as the modeling coordinates.
    boundary_data = boundary_data.to_crs(model_crs)

    x_min, y_min, x_max, y_max = boundary_data.total_bounds
    z_min, z_max = z_bounds

    bounds = [
        x_min,
        x_max,
        y_min,
        y_max,
        z_min,
        z_max,
    ]

    return validate_bounds(bounds)

def get_domain_bounds(boundaries,z_bounds=None,model_crs=None,):
    """
    Create 3D domain bounds from either a boundary file or manual values.

    File input:
        boundaries="model_boundary.shp"
        z_bounds=(-110, -20)
        model_crs="EPSG:32613"

    Manual input:
        boundaries=[
            x_min,
            x_max,
            y_min,
            y_max,
            z_min,
            z_max,
        ]
    """

    if isinstance(boundaries, (str, Path)):
        return domain_bounds_from_file(boundaries_path=boundaries,z_bounds=z_bounds,model_crs=model_crs,)
    
    if isinstance(boundaries, (list, tuple, np.ndarray)):
        return domain_bounds_from_values(boundaries)
    
    raise TypeError(
        "OpenPlumes could not interpret the domain-boundary input.\n"
        f"Received type: {type(boundaries).__name__}\n"
        "Supply either:\n"
        "1. A boundary-file path, such as 'model_boundary.shp'; or\n"
        "2. Six numeric values in this order:\n"
        "   [x_min, x_max, y_min, y_max, z_min, z_max]"
    )