from scipy.special import lambertw      # For the tortoise coordinate
import numpy as np
import os

# Suppress all warnings to avoid the conditional rho and tau functions 
# giving annoting warnings for the np.sqrt(1-r)
import warnings
warnings.filterwarnings("ignore")

# def penrose_coords(r, t):
#     condition = r <= 1
#     T = np.where(condition, np.sqrt(1-r**2)*np.sinh(t),
#                    np.sqrt(r**2-1)*np.cosh(t))
#     R = np.where(condition, r/np.sqrt(1-r**2)/np.cosh(t),
#                    -r/np.sqrt(r**2-1)/np.sinh(t))
#     # Compactification with arctan and 
#     # normalization to avoid using pi in the TikZ code
#     R=2./np.pi*np.arctan(R)
#     T=2./np.pi*np.arctan(T)
#     return R, T

def penrose_coords(r, t):
    condition = r <= 1
    # tau = np.where(condition, np.sqrt((1-r)/(1+r))*np.sinh(t), np.sqrt(-(1-r)/(1+r))*np.cosh(t))
    # rho = np.where(condition, -np.sqrt((1-r)/(1+r))*np.cosh(t), -np.sqrt(-(1-r)/(1+r))*np.sinh(t))

    # T = 2./np.pi*(np.arctan(tau+rho)+np.arctan(tau-rho))
    # R = 2./np.pi*(np.arctan(tau+rho)-np.arctan(tau-rho)) + 1

    # r_tort = 0.5*(np.log(1+r)-np.log(np.abs(1-r)))
    # ub = np.where(condition, t-r_tort, t+r_tort)
    # vb = np.where(condition, t+r_tort, t-r_tort)

    ub = np.where(condition, np.sqrt((1-r)/(1+r))*np.exp(t), np.sqrt(-(1-r)/(1+r))*np.exp(-t))
    vb = np.where(condition, -np.sqrt((1-r)/(1+r))*np.exp(-t), np.sqrt(-(1-r)/(1+r))*np.exp(t))

    T = 2./np.pi*(np.arctan(vb)+np.arctan(ub))
    R = 2./np.pi*(np.arctan(vb)-np.arctan(ub)) + 1

    return R, T

# Set height function
def set_height(surface_type, r):
    condition = r <= 1
    r_tort = np.where(condition, 0.5*(np.log(1+r)-np.log(np.abs(1-r))), 
                      -0.5*(np.log(1+r)-np.log(np.abs(1-r))))
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
    elif surface_type == "kerr_schild":
        height = np.where(condition, r - r_tort, -r - r_tort)
    # elif surface_type == "minimal_gauge":
    #     height = - r - 2*np.log(r) + np.log(np.abs(r - 1)) + 2
    elif surface_type == "misner":
        height = np.where(condition, r - r_tort - np.sqrt(1+r**2), -r - r_tort+np.sqrt(1+r**2))
    else:
        raise ValueError("Invalid surface_type")
    return height

t = np.tan(np.linspace(-1, 1, 7, endpoint=True))
# Tortoise and areal coordinates
nr_points = 100
r_in = 0.5 * (-np.cos((np.arange(nr_points) * np.pi) / (nr_points - 1)) + 1)[:-1]
r_out = 2 * (-np.cos((np.arange(nr_points) * np.pi) / (nr_points - 1)) + 1)+1
r = np.concatenate((r_in, r_out))

height =  set_height("misner", r)

# Create an array to store the data
all_data = np.empty((len(r), 2 * len(t)))

# Iterate over each value of t
for i, t_val in enumerate(t):
    tau = t_val - height # Note the opposite sign due to definition of t vs tau
    R, T = penrose_coords(r, tau)
    
    # Store the data in the array
    all_data[:, 2 * i] = R
    all_data[:, 2 * i + 1] = T

# Iterate over each value of t
# for i, t_val in enumerate(t):
#     i += len(t)
#     tau = t_val - height # Note the opposite sign due to definition of t vs tau
#     R, T = penrose_coords(r, tau)
#     R = -R
#     R += 2 # Flip the sign of T to get the other side of the diamond
#     # Store the data in the array
#     all_data[:, 2 * i] = R
#     all_data[:, 2 * i + 1] = T


# Save the data to a CSV file
# Create data directory if it doesn't exist
if not os.path.exists('data'):
    os.makedirs('data')
# column_headers = [f'{axis}_{line}' for line in range(2*len(t)) for axis in ['R', 'T']]
column_headers = [f'{axis}_{line}' for line in range(2*len(t)) for axis in ['R', 'T']]
np.savetxt('data/deSitter.csv', all_data, delimiter=',', header=','.join(column_headers), comments='', fmt='%f')
