import numpy as np
import scipy.special as spe
import scipy.optimize as opt
import math
import random
from scipy.constants import c as cms, e as qelec, epsilon_0 as epsilon0

# from pydefault import *
'''
set of functions to compute wakes in cylindrical-symmetric DLW
04-May-2017
'''
#---------- function that solves the dispersion equation

def DispersionEqn(s,a,b,n, mu_r,epsilon):
# in this equation x corresponds to s in Ng's paper
    x = s*a
    xi = b/a
    if n==0:
     P = P_rg  (s,a,b,n)
     Pp= Pp_rg (s,a,b,n)
     return (x*Pp+x*x*xi*P/(2.0*epsilon)) #return equation (3.10) in Ng's paper
    if n>0:
     P=  P_rg (s,a,b,n)
     Pp= Pp_rg(s,a,b,n)
     R=  R_rg (s,a,b,n)
     Rp= Rp_rg(s,a,b,n)
     return ((x*x*xi*xi/(n+1.0)-n*(mu_r*epsilon+1.0))*P*R+x*xi*(epsilon*Pp*R+mu_r*Rp*P)) #return equation (4.12) in  Ng paper

# r function defined in Ng paper:equation (2.17)
def R_rg(x,a,b,n):
  return (spe.jvp(n,x*a)*spe.yn(n,x*b) - spe.yvp(n,x*a)*spe.jn(n,x*b))

# r prime function defined in Ng paper:equation (2.20)
def Rp_rg(x,a,b,n):
  return (spe.jvp(n,x*a)*spe.yvp(n,x*b) - spe.yvp(n,x*a)*spe.jvp(n,x*b))

# P function defined in Ng paper:equation (2.17)
def P_rg(x,a,b,n):
  #return (spe.jn(n,x)*spe.yn(n,x*b/a) - spe.yn(n,x)*spe.jn(n,x*b/a))
  return (spe.jn(n,x*a)*spe.yn(n,x*b) - spe.yn(n,x*a)*spe.jn(n,x*b))

# P prime function defined in Ng paper:equation (2.20)
def Pp_rg(x,a,b,n):
  #return (spe.jn(n,x)*spe.yvp(n,x*b/a) - spe.yn(n,x)*spe.jvp(n,x*b/a))
  return (spe.jn(n,x*a)*spe.yvp(n,x*b) - spe.yn(n,x*a)*spe.jvp(n,x*b))


def heaviside(x):
    x = np.array(x)
    if x.shape != ():
        y = np.zeros(x.shape)
        y[x > 0.0] = 1
        y[x == 0.0] = 0.5
    else: # special case for 0d array (a number)
        if x > 0: y = 1
        elif x == 0: y = 0.5
        else: y = 0
    return y
    

def FindMode(a,b,n,epsilon,Nmode, k_0=1, k_step=1, mu_r=1.):
   '''
   a: float
      outer radius of the DLW
   b: float
      inner radius of the DLW
   n: integer
      angular mode index (sometimes called *polar* or *azimuthal*)
   epsilon: float
      permittivity of the dielectric
   Nmode: int
      number of modes considered in the wakefield calculation
   k_0: float
      Start parameter for the root finding method
   k_step: float
      k is a steping parameter that scans through the dispersion 
      equation (smaller k is more precise but lengthier)
   mu_r: float, optional
      relative permeability of the dielectric
   
   '''

# here k is a steping parameter that scans through the dispersion 
# equation (smaller k is more precise but lengthier...

   # define internal variable for taking numerical derivative

   CurrentNmode=0
   k = k_0
   RootEqn=np.zeros(Nmode)
   
   while (CurrentNmode<Nmode):  
     Dkmin=DispersionEqn(k*1.,a,b,n,mu_r,epsilon)             #DispersionEqn(x,a,b,n, mu_r,epsilon)
     Dkmax=DispersionEqn(k+k_step,a,b,n,mu_r,epsilon)
     
     if (Dkmax>0.) and (Dkmin<0.):
        RootEqn[CurrentNmode]=opt.fsolve(DispersionEqn,k,args=(a,b,n,mu_r,epsilon),xtol=1e-6)
        print(k, 'pos. sl', CurrentNmode, Dkmax, Dkmin, RootEqn[CurrentNmode])
        CurrentNmode=CurrentNmode+1
     
     if (Dkmax<0.) and (Dkmin>0.):
        RootEqn[CurrentNmode]=opt.fsolve(DispersionEqn,k,args=(a,b,n,mu_r,epsilon), xtol=1e-6)
        print(k, 'neg. sl', CurrentNmode, Dkmax, Dkmin, RootEqn[CurrentNmode])
        CurrentNmode=CurrentNmode+1
     
     k=k+k_step
  
# get the amplitude of each mode and wavevector
    
   if n==0:
      NormalizationCGS2SI=4.0*qelec/(a*b*epsilon)/(4*math.pi*epsilon0)  # Eq 39
      delta=(a-b) #relative to delta?

      while (np.any(np.abs(DispersionEqn(RootEqn+delta,a,b,n,mu_r,epsilon)-DispersionEqn(RootEqn,a,b,n,mu_r,epsilon)))>1.):#if the difference between two DispersionEqns is not small enough, then make delta smaller
          delta=delta/10.0
      print('step for derivative', delta)
# renormalize the field to the charge so units are now V/m/nC
      D_s=(DispersionEqn(RootEqn+delta,a,b,n,mu_r,epsilon)- \
           DispersionEqn(RootEqn,a,b,n,mu_r,epsilon))/delta/a # d/dx D(x) in Equation (4.11) in SI units
      Field2wake=1.0/qelec 
      RootAmplit=a*RootEqn*P_rg(RootEqn,a,b,n)/D_s*NormalizationCGS2SI*Field2wake
      RootWavVec=RootEqn/(np.sqrt(epsilon*mu_r-1.0))
      RootWavLen=2.0*math.pi/RootWavVec

 
   if n>0:
      delta=(a-b) #relative to delta?
      while (np.any(np.abs(DispersionEqn(RootEqn+delta,a,b,n,mu_r,epsilon)-DispersionEqn(RootEqn,a,b,n,mu_r,epsilon)))>1.):
          delta=delta/10.0

      print('step for derivative', delta)
      #NormalizationCGS2SI=qelec*qelec/(a*a)/(4*math.pi*epsilon0)  # Eq 5.3
      D_s=(DispersionEqn(RootEqn+delta,a,b,n,mu_r,epsilon)- \
           DispersionEqn(RootEqn,a,b,n,mu_r,epsilon))/delta/a # d/dx D(x) in Equation (4.12) in CGS
      Field2wake=1.0/qelec
      RootAmplit=a*RootEqn*P_rg(RootEqn,a,b,n)*R_rg(RootEqn,a,b,n)/D_s*Field2wake
      RootWavVec=RootEqn/(np.sqrt(epsilon*mu_r-1.0))
      RootWavLen=2.0*math.pi/RootWavVec
      

# renormalize the field to the charge so units are now V/m/NC
      

   print("----- Summary ------")
   print("mode order =",n)
   print("Roots:", RootEqn)
   print("dDisp:", D_s)
   print("Mode Amplitudes:",  RootAmplit)
   print("Mode WaveVectors / m^-1 :", RootWavVec)
   print("Mode Wavelengths /m :", 2.*np.pi/RootWavVec)
   print("--------------------")
   
   return(RootAmplit, RootWavVec)

# ---------------- definition of distributions ----------------------------------------------
#            the distribution are normalized to 1 C
#            so they need to be multiplied by the charge in the main program

def BunchDistG(zz, sigmaz):
# Gaussian bunch with rms value sigmaz
  global cms
  sigmat=sigmaz/cms
  return(1./np.sqrt(2.*np.pi*sigmat**2)*
         np.exp(-zz**2/(2.*(cms*sigmat)**2)))

def BunchDistU(zz,sigmaz):
# Uniform bunch with width 2*sigmaz
  global cms
  sigmat=sigmaz/cms
  zzz=zz-0.*sigmaz
  return (1./(2.*sigmat)*
         (-heaviside(zzz-sigmaz)+heaviside(zzz+sigmaz)))
    
def BunchDistL(zz,sigmaz):
# linearly-ramped  bunch with width 2*sigmaz
  global cms
  sigmat=sigmaz/cms
  zzz=zz-0.*sigmaz
  return (cms/(2.*sigmaz**2)*(zzz+sigmaz)*
         (-heaviside(zzz-sigmaz)+heaviside(zzz+sigmaz)))
    
def BunchDistP(zz, sigmaz):
# parabolic bunch with width 2*sigmaz
  global cms
  sigmat=sigmaz/cms
  zzz=zz-0.*sigmaz
  return (1.*cms/(4./3.*sigmaz**3)*(sigmaz**2-zzz**2)*
         (-heaviside(zzz-sigmaz)+heaviside(zzz+sigmaz)))


def BunchDistrib (Distribution, zz, sigmaz):
# evaluation of function "Distribution" with length sigmaz on array zz
#     mainly for plotting prupose 
   zzmean=np.mean(zz)
   zeval=zz-zzmean;
   dz=np.abs(zz[1]-zz[0])
   print("dz:", dz)
   MyBunch=Distribution(zeval,sigmaz)
   return(zeval,MyBunch)


#---------- for tracking purpose

def MonteCarlo (Distribution, sigmaz, Nsample):
#   print Distribution(0) 
# we assume the function is max at 0 for now // improve this later  
   MyMax=Distribution(0,sigmaz)
   index=0
   Zarray=np.zeros(Nsample)
   while (index<Nsample):
      zt=random.uniform(-4*sigmaz,4*sigmaz)
      yt=random.uniform(0, MyMax)
      if yt < Distribution(zt,sigmaz):
         Zarray[index]=zt
         index=index+1
   return Zarray	 


#---------- Green's function and Wake potential calculation 

def GreenFunction(RootAmplit, RootWavVec, zmin, zmax,Nz):

   zz=np.linspace(zmin,zmax,Nz)
   WakeGreen=0.0*zz
#   dz=np.abs(zz[1]-zz[0])
#   print "dz:", dz
   Nmode=len(RootAmplit)
   for i in range(Nmode):
      print(i, RootAmplit[i], RootWavVec[i])
      WakeGreen=WakeGreen+RootAmplit[i]*(r0*r/(a*a))**n*np.cos(RootWavVec[i]*zz)
   return(zz,WakeGreen)

def Long_GreenFunction(RootAmplit, RootWavVec, r0,r, b, a, n, zmin, zmax,Nz,epsilon, mu_r=1):
   zz=np.linspace(zmin,zmax,Nz)
   WakeGreen=0.0*zz
   Nmode=len(RootAmplit)
   for i in range(Nmode):
      print('LGreen:', i, RootAmplit, RootWavVec)
      if n>0:
        NormalizationCGS2SI=8.0*qelec/(a*a)*(r0*r/(b*b))**n/(4*math.pi*epsilon0)
        WakeGreen=WakeGreen+NormalizationCGS2SI*RootAmplit[i]*np.cos(RootWavVec[i]*zz)#equation (4.11) in Ng paper
      if n==0:
        WakeGreen= WakeGreen + RootAmplit[i]*np.cos(RootWavVec[i]*zz)  #equation (3.9) in Ng paper
   return(zz,WakeGreen)

def Trans_GreenFunction(RootAmplit, RootWavVec, r0,r, b, a, n, zmin, zmax,Nz,epsilon, mu_r=1):
   '''
   mu_r: float
      relative permeability of the dielectric
   '''

   zz=np.linspace(zmin,zmax,Nz)
   WakeGreen=0.0*zz
   Nmode=len(RootAmplit)
   F_=np.zeros(Nmode)
   for i in range(Nmode):
      print(i, RootAmplit[i], RootWavVec[i])
      #P=P_rg(RootEqn[i],a,b,n)
      #R=R_rg(RootEqn[i],a,b,n)
      NormalizationCGS2SI=qelec*qelec/(a*a)*  (r0/a)**n  *(r/a)**(n-1)  *8.0*n*np.sqrt(epsilon-1.0)/(b/a)**(2.0*n) /(4*math.pi*epsilon0)
      #F_[i]=8.0*n*np.sqrt(epsilon-1.)/((b/a)**(2*n))*P*R/D_s[i]*a        #equation (5.4)
      WakeGreen=WakeGreen+RootAmplit[i]/(RootWavVec[i]*np.sqrt(epsilon*mu_r-1.0)*a) * np.sin(RootWavVec[i]*zz)*NormalizationCGS2SI #equation (5.3), RootWavVec=RootEqn/(np.sqrt(epsilon*mu-1.0)) RootAmplit=a*RootEqn*P_rg(RootEqn,a,b,n)*R_rg(RootEqn,a,b,n)/D_s*Field2wake
   WakeGreen=WakeGreen/qelec#missing mu here ,different from equation (4.12)?
   return(zz,WakeGreen)

def WakePotential (Distribution, WakeGreen, zz, sigmaz):
   zzmean=np.mean(zz)
   zeval=zz-zzmean;
   dz=np.abs(zz[1]-zz[0])
   print("dz:", dz)
   MyBunch=Distribution(zeval,sigmaz)
   WakePot=np.convolve(MyBunch,WakeGreen)*dz/cms
   WakePot=WakePot[0:len(zeval)]
   zWakePot=zeval;
   return(zeval,WakePot)
    

def MonteCarlo (Distribution, Interval, Nsample):
#  Montecarlo generator 
# need to estimate the max value of the distribution
   dummy=Distribution(np.linspace(np.min(Interval),np.max(Interval),10000))
   MyMax=np.max(dummy)
   index=0
   Zarray=np.zeros(Nsample)
   while (index<Nsample):
      zt=random.uniform(np.min(Interval),np.max(Interval))
      yt=random.uniform(0, MyMax)
      if yt < Distribution(zt):
         Zarray[index]=zt
         index=index+1
   return Zarray	 


