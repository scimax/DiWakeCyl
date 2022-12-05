from DWFA_cyl_func_Ng import *
import matplotlib.pyplot as plt
#from DWFA_cyl_func_08032016 import *
# structure parameters:
# becareful Ng switch convention wrt Gai (a>b in Ng!!!)
b = 450e-6
a = 550e-6
epsilon = 4.41 #dielectric constant of the medium
#mu=200000.0# relative magnetic permeability,Aluminum	1.000022 ;Iron (99.95% pure Fe annealed in H) 200000,Iron (99.8% pure) 5000 from Wikipedia
mu=1.0000

# bunch length (rms, fwhm, etc.. depends on selected distribution)
sigmaz=3e-3
Q=100.e-9 #? source particle carrying charge Q :100nC


# number of mode we want to get
Nmode = 4

# compute the Green's function
Nz=10000
zmin=0
zmax=16*sigmaz

r=b
r0=1.e-3#b
azimuthal_mode=0  #m=0

RootAmplit, RootWavVec= FindMode(a,b,azimuthal_mode,epsilon,Nmode,k_step=.3)
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




#dataE0 = np.loadtxt("wake_gai_long_m=0_1C.txt")
#dataFr0 = np.loadtxt("Fr0.txt")
#dataE1 = np.loadtxt("E1.txt")
#dataFr1 = np.loadtxt("Fr1.txt")
# plots 
plt.style.use("tex.mplstyle")
plt.style.use("default")
plt.subplot(211)
plt.plot (zGreen/sigmaz, -WlGreen,'-',label='Ng')
plt.ylabel(r'$w_l(r,\zeta)$ (V/m/C)')
plt.xlabel(r'$\zeta/\sigma_z$')
#plt.legend(loc='upper center')
plt.subplot(212)
plt.plot (zGreen/sigmaz, WtGreen,'-',label='Ng')
#plt.plot(dataFr1[:,0]*10.0,dataFr1[:,1]*1.e7,label='commercial')
plt.ylabel(r'$w_t(r,\zeta)$ (V/m/C)')
plt.xlabel(r'$\zeta/\sigma_z$')
#plt.legend(loc='upper center')
plt.tight_layout()

plt.figure()
plt.subplot(211)
plt.plot (zBunchG/sigmaz, BunchG/max(BunchG)*max(WakePotG_l),'-')
plt.plot (zWakePotG/sigmaz, -WakePotG_l,'-')
plt.ylabel(r'$E(z)$ (V/m)')
plt.xlabel(r'$z/\sigma_z$')

plt.xlim(-4, zmax/sigmaz/2.)
plt.subplot(212)
plt.plot (zBunchG/sigmaz, BunchG/max(BunchG)*max(WakePotG_t),'-')
plt.plot (zWakePotG_t/sigmaz, -WakePotG_t,'-')
plt.ylabel(r'$F_r(z)$ (eV/m)')
plt.xlabel(r'$z/\sigma_z$')
plt.xlim(-4, zmax/sigmaz/2.)
plt.tight_layout()


f = open('wake_tmp.dat', 'w')

f.write (str(Nz)+"\t"+str(0)+"\n")
f.write (str(0)+"\t"+str(0)+"\n")
f.write (str(0)+"\t"+str(0)+"\n")
for i in range(Nz):
   f.write(str(zGreen[i])+"    "+str(WlGreen[i])+"\n")

f.close()


# this should be the LAST line
plt.show()
