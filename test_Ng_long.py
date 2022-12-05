from DWFA_cyl_func_Ng import *
import matplotlib.pyplot as plt
from scipy.special import jnp_zeros

#from DWFA_cyl_func_08032016 import *
# structure parameters:
# becareful Ng switch convention wrt Gai (a>b in Ng!!!)
b = 2e-3
a = 5e-3
epsilon = 3.0#dielectric constant of the medium

# bunch length (rms, fwhm, etc.. depends on selected distribution)
sigmaz=1.e-3
Q=100.e-9 #? source particle carrying charge Q :100nC

# Start parameter for root finding method
u_11 = jnp_zeros(1,1)[0]

# number of mode we want to get
Nmode = 4

# compute the Green's function
Nz=10000
zmin=0
zmax=160*sigmaz

r=b
r0=1.e-3#b
azimuthal_mode=0  #m=0

RootAmplit, RootWavVec= FindMode(a,b,azimuthal_mode,epsilon,Nmode,k_0=u_11/(2*np.pi * a))
n=azimuthal_mode
zGreen, WlGreen = Long_GreenFunction(RootAmplit, RootWavVec, r0,r, b, a, n, zmin, zmax,Nz,epsilon)
zGreen, WtGreen = Trans_GreenFunction(RootAmplit, RootWavVec, r0,r, b, a, n, zmin, zmax,Nz,epsilon)#(RootAmplit, RootWavVec, r, r0, a, azimuthal_mode, zmin, zmax,Nz)
# case of a Gaussian
zWakePotG, WakePotG_l= WakePotential (BunchDistG, WlGreen, zGreen, sigmaz)
zWakePotG_t, WakePotG_t= WakePotential (BunchDistG, WtGreen, zGreen, sigmaz)
WakePotG_l=WakePotG_l*Q #why do we need this?
WakePotG_t=WakePotG_t*Q
zBunchG, BunchG= BunchDistrib (BunchDistG, zGreen, sigmaz)
BunchG=BunchG*Q




dataE0 = np.loadtxt("wake_commercial/wake_gai_long_m=0_1C.txt")
dataE1 = np.loadtxt("wake_commercial/wake_gai_long_m=1_1C.txt")
dataE2 = np.loadtxt("wake_commercial/wake_gai_long_m=2_1C.txt")
plt.style.use("tex.mplstyle")
plt.subplot(111)
plt.plot (zGreen/sigmaz, -WlGreen,'-',label='Ng')
plt.plot(dataE0[:,0]*10.0,dataE0[:,1],'--', label='commercial 0')
plt.plot(dataE1[:,0]*10.0,dataE1[:,1],'--', label='commercial 1')
plt.plot(dataE2[:,0]*10.0,dataE2[:,1],'--', label='commercial 2')
plt.ylabel(r'$w_l(r,\zeta)$ (V/m/C)')
plt.xlabel(r'$\zeta$ (mm)')
plt.legend(loc='upper right')
plt.tight_layout()


# this should be the LAST line
plt.show()
