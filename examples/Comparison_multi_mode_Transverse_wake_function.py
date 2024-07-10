# author: M. Kellermeier
# date: 2024.07.10
# Desc: Comparing the legacy version of the full transverse wake function `multimode_trans_greens_function` 
#       and the new loop-less numpy-style version, both in terms of data (set boolean flag accordingly)
#       and in computation time


import logging

log_format = "%(asctime)s [%(threadName)-10.10s] [%(name)-14.14s] [%(levelname)-5.5s]  %(message)s"
logging.basicConfig( #filename='{:%Y-%m-%dT%H%M%S}.log'.format(datetime.now()),
                    level=logging.INFO,
                    format=log_format,)

# child logger for local use
logger = logging.getLogger(__name__ + '.nb')
logger.setLevel(logging.DEBUG)
logging.getLogger('DWFA_cyl_func_Ng').setLevel(logging.ERROR)

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

sys.path.append('..') 
from DWFA_cyl_func_Ng import FindMode,  multimode_trans_greens_function

a = 0.69e-3
b = 1.29e-3
eps_r = 9.6
l_dlw = 0.09

m_range = (1,8)
n_max = 8

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
Δr_trail = 50e-6
theta = np.deg2rad(0)
N_z = 1001

kwargs= {
        "r0" : r_offset, "r": r_offset + Δr_trail, 'theta': theta,
        "b": a, "a" : b,
        "zmin": 0, "zmax" : max_wvl,
        "Nz": N_z , "epsilon": eps_r
    }
m = 1
# for m in range(*m_range):
kwargs["RootAmplit"] = mode_ampl
kwargs["RootWavVec"] = mode_wave_vec
kwargs["m_modes"] = np.arange(*m_range)

if __name__ == '__main__':
    # If `True` the functions are compared in terms of computational time aka benchmarked.
    # If `False` the returned data is compared in a plot.
    bool_benchmark = True
    if bool_benchmark:
        import timeit

        t1 = timeit.timeit("multimode_trans_greens_function(**kwargs, mode='legacy')", globals=locals(), number=10000)
        print(f'legacy version:    {t1*1000:2.3f} ms')
        
        t2 = timeit.timeit("multimode_trans_greens_function(**kwargs)", globals=locals(), number=10000)
        print(f'loop-less version: {t2*1000:2.3f} ms')

        print(f'The second method requires { (1 - t2/t1) * 100 : 2.0f}% less time.')
        print('Note that the legacy version requires even more time if the WARN-log messages are printed. Set the log level accordingly for reasonable results.')

        print('NOTE: If the legacy-version of `multimode_trans_greens_function` uses the legacy version of `Trans_GreenFunction`, the speed improvment is about 60%.')
        print('      If the legacy-version of `multimode_trans_greens_function` uses the numpy-style version of `Trans_GreenFunction`, the speed improvment is about 30%.')
    else:
        z_trans_wake, trans_wake = multimode_trans_greens_function(**kwargs)
        z_trans_wake_legacy, trans_wake_legacy = multimode_trans_greens_function(**kwargs, mode='legacy')
        fig, ax = plt.subplots()
        ax.plot(z_trans_wake_legacy, trans_wake_legacy, label='legacy')
        ax.plot(z_trans_wake, trans_wake, ls='--', label = 'numpy-style')
        ax.xaxis.set_major_formatter(ticker.EngFormatter(unit='m'))
        ax.set_ylabel('wake field function / V/(m C)')
        ax.set_title('Multi polar modes, multi radial modes')
        ax.legend()
        ax.grid()
        plt.show()