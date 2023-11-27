from scipy.special import lambertw      # For the tortoise coordinate
import numpy as np
import os

# Suppress all warnings to avoid the conditional rho and tau functions 
# giving annoting warnings for the np.sqrt(1-r)
import warnings
warnings.filterwarnings("ignore")

def penrose_coords(r, t):
    # https://en.wikipedia.org/wiki/Kruskal%E2%80%93Szekeres_coordinates#Definition
    # With G=1, M=1/2 (dimensionless coordinates)
    condition = r > 1
    rho = np.where(condition, np.sqrt(r - 1) * np.exp(r / 2) * np.cosh(t / 2),
                   np.sqrt(1 - r) * np.exp(r / 2) * np.sinh(t / 2))
    tau = np.where(condition, np.sqrt(r - 1) * np.exp(r / 2) * np.sinh(t / 2),
                   np.sqrt(1 - r) * np.exp(r / 2) * np.cosh(t / 2))
    # rho = np.sqrt(r - 1) * np.exp(r / 2) * np.cosh(t / 2)
    # tau = np.sqrt(r - 1) * np.exp(r / 2) * np.sinh(t / 2)
    R = 2./np.pi*(np.arctan(tau+rho)-np.arctan(tau-rho))
    T = 2./np.pi*(np.arctan(tau+rho)+np.arctan(tau-rho))
    return R, T

# Set height function
def set_height(surface_type, r):
    r_tort =  r + np.log(np.abs(r-1)); 
    if surface_type == "standard":
        height = 0
    elif surface_type == "outgoing_null":
        height = -r_tort # Outgoing null surfaces u = t - r
    elif surface_type == "ingoing_null":
        height = r_tort # Ingoing null surfaces v = t + r
    elif surface_type == "future_hyperboloidal":
        # Future hyperboloidal tau = t - \sqrt{L^2+r^2} + L
        # Change value of L_CMC to change K = 3/L_CMC
        L_CMC = 1.0
        height = -(np.sqrt(L_CMC**2 + r_tort**2) - L_CMC)
    elif surface_type == "past_hyperboloidal":
        # Past hyperboloidal tau = t + \sqrt{1+r^2}
        height = np.sqrt(1 + r_tort**2)
    elif surface_type == "gullstrand_painleve":
        height = 2*np.sqrt(r) - np.log( np.abs((np.sqrt(r)+1)/(np.sqrt(r)-1)) )
    elif surface_type == "ingoing_EF":
        height = np.log(np.abs(r-1))
    elif surface_type == "minimal_gauge":
        height = - r - 2*np.log(r) + np.log(np.abs(r - 1)) + 2
    else:
        raise ValueError("Invalid surface_type")
    return height

t = np.tan(np.linspace(-np.pi/2., np.pi/2, 8, endpoint=False)[1:])
# Tortoise and areal coordinates
nr_points = 100
r_tort = np.tan(np.linspace(-np.pi/2, np.pi/2, nr_points, endpoint=False)[1:])
r_in = 0.5 * (-np.cos((np.arange(nr_points) * np.pi) / (nr_points - 1)) + 1)[1:-1]
r_out = np.real(lambertw(np.exp(r_tort-1))+1)[1:]
r = np.concatenate((r_in, r_out))

height =  set_height("ingoing_EF", r)

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
np.savetxt('data/schwarzschild.csv', all_data, delimiter=',', header=','.join(column_headers), comments='', fmt='%f')
