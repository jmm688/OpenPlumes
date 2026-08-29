# CRS
import numpy as np
from pyproj import Transformer

def validate_projection(df_, X1, Y1, source_crs, model_crs,):
    
    """Verify projected X/Y by transforming them back to source coordinates."""
    reverse_transformer = Transformer.from_crs(model_crs, source_crs, always_xy=True,)
    X_test, Y_test = reverse_transformer.transform(df_["X"].to_numpy(), df_["Y"].to_numpy(),)
    coordinates_match = (np.allclose(X_test, df_[X1].to_numpy(), atol=1e-8,)
        and np.allclose(Y_test,df_[Y1].to_numpy(),atol=1e-8,))
    if not coordinates_match:
        raise ValueError('Projected coordinates failed the round-trip check.')
    
    return True

def project_coordinates(df, Longitude, Latitude, source_crs, model_crs,):
    if Longitude not in df.columns:
        raise ValueError(f"Missing longitude column: {Longitude}")

    if Latitude not in df.columns:
        raise ValueError(f"Missing latitude column: {Latitude}")  
    projected_df = df.copy()

    transformer = Transformer.from_crs(source_crs, model_crs, always_xy = True)
    projected_df['X'], projected_df['Y'] = transformer.transform(projected_df[Longitude].to_numpy(), projected_df[Latitude].to_numpy(),)
    
    # checking if projection transformation worked out
    validate_projection(projected_df, X1=Longitude, Y1=Latitude, source_crs=source_crs, model_crs=model_crs,)
    
    return projected_df