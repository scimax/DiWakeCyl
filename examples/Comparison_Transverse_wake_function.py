import logging

log_format = "%(asctime)s [%(threadName)-10.10s] [%(name)-14.14s] [%(levelname)-5.5s]  %(message)s"
logging.basicConfig( #filename='{:%Y-%m-%dT%H%M%S}.log'.format(datetime.now()),
                    level=logging.INFO,
                    format=log_format,)

# child logger for local use
logger = logging.getLogger(__name__ + '.nb')
logger.setLevel(logging.DEBUG)

def printl(*args, sep=' '):
    logging.info(sep.join([str(val) for val in args]) )

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.constants import e as qelec, epsilon_0
import math


sys.path.append('..') 
from DWFA_cyl_func_Ng import FindMode, Long_GreenFunction, Trans_GreenFunction, WakePotential, BunchDistG, BunchDistrib, BunchDistU, multimode_trans_greens_function

a = 0.69e-3
b = 1.29e-3
eps_r = 9.6
l_dlw = 0.09

m_range = (1,8)
n_max = 8
# large_al2o3_streaker.find_modes((1,8), 10)
# N_z_wake = 1001


from scipy.special import jnp_zeros
u_11 = jnp_zeros(1,1)[0]
assert m_range[1] > m_range[0] and len(m_range) == 2

m_range = tuple(m_range)
num_m = m_range[1] - m_range[0]
num_n = n_max

mode_ampl = np.empty((num_m, num_n))
mode_wave_vec = np.empty((num_m, num_n))

for m in range(*m_range):
    logger.info(f"Computing {n_max} amplitudes and wavevectors for angular mode number m={m}...")
    mode_ampl[m - m_range[0]], mode_wave_vec[m - m_range[0]] =\
        FindMode(b, a, n=m, epsilon=eps_r, Nmode = n_max, num_k_sampling = 500) 

max_wvl = 2*np.pi/np.min(mode_wave_vec)
logger.debug(f"Max. wavelength: {max_wvl:1.2e} m")

##----------
r_offset = 0.45e-3
N_z = 1001

kwargs= {
        "r0" : r_offset, "r": r_offset, 
        "b": a, "a" : b,
        "zmin": 0, "zmax" : max_wvl,
        "Nz": N_z , "epsilon": eps_r
    }
m = 1
# for m in range(*m_range):
kwargs["RootAmplit"] = mode_ampl[m - m_range[0]]
kwargs["RootWavVec"] = mode_wave_vec[m - m_range[0]]
kwargs["n"] = m

if __name__ == '__main__':
    import timeit

    t1 = timeit.timeit("Trans_GreenFunction(**kwargs, mode='legacy')", globals=locals(), number=10000)
    print(f'{t1*1000:2.3f} ms')
    
    
    t2 = timeit.timeit("Trans_GreenFunction(**kwargs)", globals=locals(), number=10000)
    print(f'{t2*1000:2.3f} ms')

    print(f'The second method requires { (1 - t2/t1) * 100 : 2.0f}% less time.')
    print('Note that the legacy version requires even more time if the WARN-log messages are printed. Set the log level accordingly for reasonable results.')