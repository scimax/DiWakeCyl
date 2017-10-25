from pydefault import *
from DWFA_cyl_func_Ng import *
#from DWFA_cyl_func_08032016 import *
# structure parameters:
# becareful Ng switch convention wrt Gai (a>b in Ng!!!)
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
azimuthal_mode=0  #m=0

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




dataE0 = np.loadtxt("wake_gai_long_m=0_1C.txt")
#dataFr0 = np.loadtxt("Fr0.txt")
#dataE1 = np.loadtxt("E1.txt")
#dataFr1 = np.loadtxt("Fr1.txt")
# plots 
plt.subplot(211)
plt.plot (zGreen/sigmaz, -WlGreen,'-',label='Ng')
plt.plot(dataE0[:,0]*10.0,dataE0[:,1],'--', label='commercial 0')
plt.ylabel(r'$w_l(r,\zeta)$ (V/m/C)', fontsize=18)
plt.xlabel(r'$\zeta$ (mm)', fontsize=18)
#plt.legend(loc='upper center')
plt.subplot(212)
plt.plot (zGreen/sigmaz, WtGreen,'-',label='Ng')
#plt.plot(dataFr1[:,0]*10.0,dataFr1[:,1]*1.e7,label='commercial')
plt.ylabel(r'$w_t(r,\zeta)$ (V/m/C)', fontsize=18)
plt.xlabel(r'$\zeta$ (mm)', fontsize=18)
#plt.legend(loc='upper center')

plt.figure()
plt.subplot(211)
plt.plot (zBunchG/sigmaz, BunchG/max(BunchG)*max(WakePotG_l),'-')
plt.plot (zWakePotG/sigmaz, -WakePotG_l,'-')
plt.ylabel(r'$E(z)$ (V/m)', fontsize=18)
plt.xlabel(r'$z/\sigma_z$', fontsize=18)

plt.xlim(-4, zmax/sigmaz/2.)
plt.subplot(212)
plt.plot (zBunchG/sigmaz, BunchG/max(BunchG)*max(WakePotG_t),'-')
plt.plot (zWakePotG_t/sigmaz, -WakePotG_t,'-')
plt.ylabel(r'$F_r(z)$ (eV/m)', fontsize=18)
plt.xlabel(r'$z/\sigma_z$', fontsize=18)
plt.xlim(-4, zmax/sigmaz/2.)


# this should be the LAST line
plt.show()
