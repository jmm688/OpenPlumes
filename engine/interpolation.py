from scipy.interpolate import RBFInterpolator

def interpolateRBF(xyz_points, contaminant, kernel, neighbors_num, epsilon_num):
    
    rbf_fit = RBFInterpolator(
        xyz_points, 
        contaminant,
        kernel = kernel, 
        neighbors = neighbors_num,
        epsilon=epsilon_num,
        )

    return rbf_fit