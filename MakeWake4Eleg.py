import sys
# argument are inner radius, outer radius, epsilon, Nmodes, zmax
from pydefault import *
from DWFA_cyl_func_Ng import *
import os
print(sys.argv)

# structure parameters:
# becareful Ng switch convention wrt Gai (a>b in Ng!!!)
b      = float(sys.argv[1])
a      = float(sys.argv[2])
epsilon= float(sys.argv[3])
Nmode  = int(float(sys.argv[4]))
zmax   = float(sys.argv[5])
fin    = sys.argv[6]
fout   = sys.argv[7]

command='cp '+fin+' '+fout
os.system(command)

print('argument:')
print('inner_radius =', b)
print('outer_radius =', a)
print('expsilon_rel =', epsilon)
print('NumberofModes=', Nmode)
print('Maximum z pos=', zmax)

### initial set that works
#b = 450e-6
#a = 550e-6
#epsilon = 4.41 #dielectric constant of the medium
mu=1.0000
# number of mode we want to get
#Nmode = 1


# if fin sdds file is given zmax is OVERWRTTTEN ** need a condition here eventually **
command = 'sddsprocess '+fin+' -process=t,largest,max_%s  ' \
                            +' -process=t,smallest,min_%s ' 
os.system(command)
min_t=float(os.popen('sdds2stream '+fin+' -par=min_t').read())
max_t=float(os.popen('sdds2stream '+fin+' -par=max_t').read())

print('minimum', min_t)
print('maximum', max_t)

cms=299792458.0
zmax= cms*np.abs((max_t-min_t))

print('Maximum z pos=', zmax)

# compute the Green's function
Nz=10000
zmin=0
# zmax is changed to have conservative headroom
zmax=2.*zmax

r=b
r0=1.e-3  
azimuthal_mode=1  #m=0

RootAmplit, RootWavVec= FindMode(a,b,azimuthal_mode,mu,epsilon,Nmode,10)
n=azimuthal_mode
zGreen, WlGreen = Long_GreenFunction (RootAmplit, RootWavVec, r0, r, b, a, n, zmin, zmax, Nz, mu, epsilon)
zGreen, WtGreen = Trans_GreenFunction(RootAmplit, RootWavVec, r0, r, b, a, n, zmin, zmax, Nz, mu, epsilon)


f = open('tmp_wakeT.sdds', 'w')

f.write ("SDDS1\n")
f.write ("&column name=z, units=m,   type=double, &end  \n")
f.write ("&column name=t, units=s,    type=double, &end  \n")
f.write ("&column name=W, units=V/C/m,  type=double, &end  \n")
f.write ("&data mode=ascii, &end          \n")              
f.write ("! page number 1                      \n")         
f.write (str(Nz)+"\n")

for i in range(Nz):
   f.write(str(zGreen[i])+"\t"+str(zGreen[i]/cms)+"\t"+str(WtGreen[i])+"\n")
   
f.close()

f = open('tmp_wakeL.sdds', 'w')

cms=299792458.0

f.write ("SDDS1\n")
f.write ("&column name=z, units=m,   type=double, &end  \n")
f.write ("&column name=t, units=s,    type=double, &end  \n")
f.write ("&column name=W, units=V/C,  type=double, &end  \n")
f.write ("&data mode=ascii, &end          \n")              
f.write ("! page number 1                      \n")         
f.write (str(Nz)+"\n")

for i in range(Nz):
   f.write(str(zGreen[i])+"\t"+str(zGreen[i]/cms)+"\t"+str(WlGreen[i])+"\n")
   
f.close()


