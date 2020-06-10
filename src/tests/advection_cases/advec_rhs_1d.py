import numpy as np 
from dg import Line

def advec_rhs_1d(u,time,speed):
    """
    Evaluate RHS flux in 1D advection
    """
    alpha = 1
    du = np.zeros([n_fp*n_faces,k_elem])
    du = (u[vmap_m]-u[vmap_p])*(speed*nx-(1-alpha)*np.abs(speed*nx))/2

    uin =  -np.sin(speed*time)
    du[map_i] = (u[vmap_i]-uin)*(a*nx[map_i]-(1-alpha)*np.abs(a*nx[map_i]))/2
    du[map_o] = 0

    rhs_u = -a*rx*np.matmul(diff_matrix,u)+lift*(Fscale*(du))
    return rhs_u
