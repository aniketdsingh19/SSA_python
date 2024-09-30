# import numpy as np
# def multidimgrid(ndim_interp, unit_grid, axes, vector):

#     nodes = 2 ** ndim_interp
#     axes_offset = axes[0, :]
#     axes_scale = axes[1] - axes[0]
#     unit_axes = np.zeros((2,ndim_interp))
#     unit_axes[0, :] = 0.0 #initialize it
#     unit_axes[1, :] = 1.0

#     unit_vector = (vector-axes_offset)/axes_scale
#     #print("unit vector:", vector-axes_offset, axes_scale)
#     reform_grid = np.zeros(nodes)
#     ii = 0

#     for i in range(2):
#         for j in range(2):
#             reform_grid[ii]  =unit_grid[i][j]
#             ii += 1


#     multi_interp = 0.0
#     for i in range(0, nodes):
#         index_grid_value = reform_grid[i]
#         interpolate_product=1.0
#         index_product=1.0

#         for j in range(ndim_interp):
#             i1 = int((((i-1)/2**(j-1)) % 2))
#             i2 = int((((i-1)*ndim_interp + (j-1)) % ndim_interp ))
#             index_product = ( 1.0 - abs(unit_vector[j] - unit_axes[i1,i2]) )
#             interpolate_product = interpolate_product * index_product

#         index_sum = index_grid_value * interpolate_product
#         multi_interp = multi_interp + index_sum
#         #print("multi:", multi_interp)

#     return multi_interp

import numpy as np

def multidimgrid(ndim_interp, unit_grid, axes, vector):
    nodes = 2 ** ndim_interp
    axes_offset = axes[0, :]
    axes_scale = axes[1] - axes[0]
    
    unit_vector = (vector - axes_offset) / axes_scale

    # Flatten the grid
    reform_grid = unit_grid.flatten()

    # Precompute unit axes once
    unit_axes = np.array([0.0, 1.0])

    multi_interp = 0.0

    # Loop over all grid points
    for i in range(nodes):
        index_grid_value = reform_grid[i]
        interpolate_product = 1.0

        for j in range(ndim_interp):
            i1 = (i // (2 ** j)) % 2  # Simplified index
            index_product = 1.0 - abs(unit_vector[j] - unit_axes[i1])
            interpolate_product *= index_product

        # Accumulate interpolated value
        multi_interp += index_grid_value * interpolate_product

    return multi_interp
