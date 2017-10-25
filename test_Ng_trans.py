from pydefault import *
from DWFA_cyl_func_Ng import *
#from DWFA_cyl_func_08032016 import *
# structure parameters:
# be careful Ng switch convention wrt Gai (a>b in Ng!!!)
b = 2e-3
a = 5e-3
epsilon = 3.0#dielectric constant of the medium
#mu=200000.0# relative magnetic permeability,Aluminum	1.000022 ;Iron (99.95% pure Fe annealed in H) 200000,Iron (99.8% pure) 5000 from Wikipedia
mu=1.0000

# bunch length (rms, fwhm, etc.. depends on selected distribution)
sigmaz=1.e-3
Q=100.e-9 #? source particle carrying charge Q :100nC


# number of mode we want to get
Nmode = 4

# compute the Green's function
Nz=10000
zmin=0
zmax=160*sigmaz

r=b
r0=1.e-3#b
azimuthal_mode=1  #m=0

RootAmplit, RootWavVec= FindMode(a,b,azimuthal_mode,mu,epsilon,Nmode,0.1)
n=azimuthal_mode
zGreen, WlGreen = Long_GreenFunction(RootAmplit, RootWavVec, r0,r, b, a, n, zmin, zmax,Nz,mu,epsilon)
zGreen, WtGreen = Trans_GreenFunction(RootAmplit, RootWavVec, r0,r, b, a, n, zmin, zmax,Nz,mu,epsilon)#(RootAmplit, RootWavVec, r, r0, a, azimuthal_mode, zmin, zmax,Nz)
# case of a Gaussian
zWakePotG, WakePotG_l= WakePotential (BunchDistG, WlGreen, zGreen, sigmaz)
zWakePotG_t, WakePotG_t= WakePotential (BunchDistG, WtGreen, zGreen, sigmaz)
WakePotG_l=WakePotG_l*Q #why do we need this?
WakePotG_t=WakePotG_t*Q
zBunchG, BunchG= BunchDistrib (BunchDistG, zGreen, sigmaz)
BunchG=BunchG*Q




dataE1 = np.loadtxt("wake_commercial/wake_gai_trans_m=1_1C.txt")
dataE2 = np.loadtxt("wake_commercial/wake_gai_trans_m=2_1C.txt")

plt.subplot(111)
plt.plot (zGreen/sigmaz, WtGreen,'-',label='Ng')
plt.plot(dataE1[:,0]*10.0,dataE1[:,1],'--', label='commercial')
plt.ylabel(r'$w_t (r,\zeta)$ (V/m/C)', fontsize=18)
plt.xlabel(r'$\zeta$ (mm)', fontsize=18)
plt.legend(loc='upper right')

plt.tight_layout()

# this should be the LAST line
plt.show()
