import numpy as np
import os

def penrose_coords(r, t):
    R = 1./np.pi*(np.arctan(t+r)-np.arctan(t-r))
    T = 1./np.pi*(np.arctan(t+r)+np.arctan(t-r))
    return R, T

# Set height function
def set_height(surface_type, r):
    if surface_type == "standard":
        height = 0
    elif surface_type == "outgoing_null":
        height = -r # Outgoing null surfaces u = t - r
    elif surface_type == "ingoing_null":
        height = r # Ingoing null surfaces v = t + r
    elif surface_type == "future_hyperboloidal":
        # Future hyperboloidal tau = t - \sqrt{L^2+r^2} + L
        # Change value of L_CMC to change K = 3/L_CMC
        L_CMC = 1.0
        height = -(np.sqrt(L_CMC**2 + r**2) - L_CMC)
    elif surface_type == "past_hyperboloidal":
        # Past hyperboloidal tau = t + \sqrt{1+r^2}
        height = np.sqrt(1 + r**2)
    elif surface_type == "future_asymptotically_null":
        height = -np.sqrt(1 / (1 + r) + r**2)
    else:
        raise ValueError("Invalid surface_type")
    return height


t = np.tan(np.linspace(-np.pi/2., np.pi/2, 8, endpoint=False)[1:])
r = np.tan(np.linspace(0, np.pi/2, 300))
height = set_height("future_hyperboloidal", r)

# Create an array to store the data
all_data = np.empty((len(r), 2 * len(t)))

# Iterate over each value of t
for i, t_val in enumerate(t):
    tau = t_val - height # Note the opposite sign due to definition of t vs tau
    R, T = penrose_coords(r, tau)
    
    # Store the data in the array
    all_data[:, 2 * i] = R
    all_data[:, 2 * i + 1] = T

# Save the data to a CSV file
# Create data directory if it doesn't exist
if not os.path.exists('data'):
    os.makedirs('data')
column_headers = [f'{axis}_{line}' for line in range(len(t)) for axis in ['R', 'T']]
np.savetxt('data/minkowski.csv', all_data, delimiter=',', header=','.join(column_headers), comments='', fmt='%f')
