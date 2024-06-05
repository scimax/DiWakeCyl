import logging
import random
import math
import numpy as np
import scipy.special as spe
import scipy.optimize as opt
from scipy.constants import c as cms, e as qelec, epsilon_0 as epsilon0

logger = logging.getLogger(__name__)

# from pydefault import *
'''
set of functions to compute wakes in cylindrical-symmetric DLW
04-May-2017
'''
#---------- function that solves the dispersion equation

def DispersionEqn(s,a,b,n, epsilon, mu_r=1):
    """
    Characteristic equation for the dispersion relation of the wake function 
    of the dielectric lined waveguide (DLW). See equations (3.10) and (4.12) in Ng's paper.
    To find the modes, the roots of the characteristic equation must be found.
    
    Parameters
    ----------
    s : float
        Scanning parameter for the wave number, at which the characteristic equation is zero.
    a: float
        outer radius of the DLW.
    b: float
        inner radius of the DLW.
    n: float
        angular mode index of interest.
    epsilon: float
        relative permittivity of the dielectric. Must be larger than 1.
    mu_r: float, optional
        relative permeability of the dielectric. Must be >= 1. By default
        it is set to 1.

    Returns
    -------
    float:
        value of the characteristic equation evaluated at s.

    """
    # in this equation x corresponds to s in Ng's paper
    x = s*a
    xi = b/a
    if n==0:
        P = P_rg  (s,a,b,n)
        Pp= Pp_rg (s,a,b,n)
        return x*Pp+x*x*xi*P/(2.0*epsilon) #return equation (3.10) in Ng's paper
    if n>0:
        P=  P_rg (s,a,b,n)
        Pp= Pp_rg(s,a,b,n)
        R=  R_rg (s,a,b,n)
        Rp= Rp_rg(s,a,b,n)
        #return equation (4.12) in  Ng paper
        return (x*x*xi*xi/(n+1.0)-n*(mu_r*epsilon+1.0))*P*R+x*xi*(epsilon*Pp*R+mu_r*Rp*P)
    
def dDispersionEqn_dx_n0(s,a,b, epsilon, mu_r=1):
    '''
    Derivative of the characteristic equation of the dispersion function w.r.t. x, as required in Eq. (3.9).
    Only valid for n=0.
    '''
    x = s*a
    xi = b/a
    xi_x = s*b

    P_0 = P_rg  (s,a,b,0)
    Pp_0= Pp_rg (s,a,b,0)

    dP_dx = (-spe.j1(x) * spe.y0( xi_x) + spe.y1( x) * spe.j0( xi_x)) +\
        xi * Pp_0
    dPp_dx = spe.j1(x) * spe.y1(xi_x) - xi * spe.j0(x) * spe.yvp(1, xi_x) +\
          -spe.y1( x) * spe.j1(xi_x) + xi * spe.y0(x) * spe.jvp(1, xi_x) 
    return Pp_0 + x * dPp_dx + xi/(2.0*epsilon) * ( 2 * x * P_0 + x**2 * dP_dx )

def R_rg(x,a,b,n):
    '''
    r_m(s) function defined in Ng paper: equation (2.17)
    '''
    return (spe.jvp(n,x*a)*spe.yn(n,x*b) - spe.yvp(n,x*a)*spe.jn(n,x*b))

def Rp_rg(x,a,b,n):
    '''
    r_m'(s) prime function defined in Ng paper: equation (2.20)
    '''
    return (spe.jvp(n,x*a)*spe.yvp(n,x*b) - spe.yvp(n,x*a)*spe.jvp(n,x*b))

def P_rg(x,a,b,n):
    '''
    p_m(s) function defined in Ng paper:equation (2.17)
    '''
    #return (spe.jn(n,x)*spe.yn(n,x*b/a) - spe.yn(n,x)*spe.jn(n,x*b/a))
    return (spe.jn(n,x*a)*spe.yn(n,x*b) - spe.yn(n,x*a)*spe.jn(n,x*b))

def Pp_rg(x,a,b,n):
    '''
    p_m'(s) prime function defined in Ng paper: equation (2.20)
    '''
    #return (spe.jn(n,x)*spe.yvp(n,x*b/a) - spe.yn(n,x)*spe.jvp(n,x*b/a))
    return (spe.jn(n,x*a)*spe.yvp(n,x*b) - spe.yn(n,x*a)*spe.jvp(n,x*b))


def heaviside(x):
    '''
    Heaviside function. See also `np.heaviside`.
    https://numpy.org/doc/stable/reference/generated/numpy.heaviside.html
    '''
    return np.heaviside(x, 0.5)

def FindMode(a,b,n,epsilon,Nmode, k_lower= None, k_upper=None, num_k_sampling=50, mu_r=1.):
    '''
    Parameters
    ----------
    a: float
        outer radius of the DLW
    b: float
        inner radius of the DLW
    n: integer
        angular mode index (sometimes called *polar* or *azimuthal*)
    epsilon: float
        permittivity of the dielectric
    Nmode: int
        number of radial modes considered in the wakefield calculation
    k_lower, k_upper: float, optional
        lower and upper boundaries of the initial k-space interval for coarse
        root-finding of the dispersion characteric equation. Default is None.
        In this case, the initial interval is estimated from the waveguide dimensions.
    num_k_sampling: int, optional
        number of sampling points in the initial k-space interval for coarse 
        root-finding
    mu_r: float, optional
        relative permeability of the dielectric

    Returns
    -------
    RootAmplit: ndarray
        Array of shape (Nmode,), containing the amplitudes for all modes
        with angular mode index n, and radial mode index from 0 to Nmode (exclusively)
    RootWavVec:
        Array of shape (Nmode,), containing the wavevectors for all modes
        with angular mode index n, and radial mode index from 0 to Nmode (exclusively)

    Notes
    -----
    The amplitudes and wavevectors refer to the terms of the sum in Eq. (3.9) or Eq. (5.3),
    depending on the angular mode index.
    
    **Important**: Check that the modes found are basically all by gradually increasing the 
    number of sampling points in the k-space, `num_k_sampling`.

    '''
    # define internal variable for taking numerical derivative
    assert b < a   # Ensure that the order of inner and outer radii is correct

    #### Goal: Find intervals with sign changes similar to bisection method, but use faster  
    ####       root-finding method in each interval.
    from scipy.special import jnp_zeros, jn_zeros
    u_11 = jnp_zeros(1,1)[0]    # first root of J_1'(x)
    u_mn = jn_zeros(n, Nmode)[-1]  # N'th root of J_n(x). Mind that this is not the derivative

    RootEqn=np.zeros(Nmode)
    if k_lower is None:
        # The wavevector must be larger than the one of the fundamental mode of 
        # a completely filled waveguid
        k_lower = u_11/a
    if k_upper is None:
        # The wavevector of the TM-contribution is always larger than the one of the corresponding
        # TE-contribution for a given pair (m,n)
        # multiplied by sqrt(epsilon*mu_r-1.0) to find an upper bound, close to a unfilled waveguide
        k_upper = math.sqrt((epsilon*mu_r-1.0)) * u_mn/a

    k_samples = np.linspace(k_lower, k_upper, num_k_sampling)
    num_sign_changes = 0
    Dk_values = DispersionEqn(k_samples, a, b , n, epsilon, mu_r)

    toogle_doubling_samples = True

    while num_sign_changes < Nmode:
        k_samples = np.linspace(k_lower, k_upper, num_k_sampling)
        Dk_values = DispersionEqn(k_samples, a, b , n, epsilon, mu_r)
        
        # determine indices at which the sign changes and count them
        asign = np.sign(Dk_values)
        idx_sign_change = np.nonzero((np.roll(asign, 1) - asign) != 0)[0]
        if idx_sign_change[0] == 0:
            idx_sign_change = idx_sign_change[1:]
        
        num_sign_changes = idx_sign_change.shape[0]
        # signchange[0] = 0
        # n_sign_changes = np.sum(signchange)
        # print(signchange)
        logger.debug(f'number of sign changes in k-space interval: {num_sign_changes}')

        # Switch between refining sampling of the interval and resizing the interval
        # Only matters if the 'while'-criterion is not met
        if toogle_doubling_samples:
            logger.debug(f'doubling sampling points of k-space interval from {num_k_sampling} to {num_k_sampling * 2}')
            num_k_sampling *= 2
            toogle_doubling_samples = False
        else:
            logger.debug(f'Resizing k-space interval')
            k_upper *=2
            toogle_doubling_samples = True
        
        # termination criterion to prevent from endless loop
        if num_k_sampling > 5e4:
            logger.error(f'Termination criterion reached. Sampling size of k-space for root finding too large!')
            return None
    
    idx_sign_change = idx_sign_change[:Nmode]
    for i, (k_l, k_u) in enumerate(zip(k_samples[idx_sign_change - 1 ], k_samples[idx_sign_change - 1 ] )):
        RootEqn[i]=opt.fsolve(DispersionEqn, 0.5*(k_l + k_u) , args=(a,b,n,epsilon, mu_r), xtol=1e-6)
        if logger.level <= logging.DEBUG:
            D_k_l = DispersionEqn(k_l, a, b , n, epsilon, mu_r)
            D_k_u = DispersionEqn(k_u, a, b , n, epsilon, mu_r)
            slope_str = 'pos. slope' if D_k_l < D_k_u else 'neg. slope'
            logger.debug(((k_l, k_u), slope_str, i, D_k_l, D_k_u, RootEqn[i]))
    logger.debug('All roots: '+repr(RootEqn))

# get the amplitude of each mode and wavevector

    if n==0:
        NormalizationCGS2SI=4.0*qelec/(a*b*epsilon)/(4*math.pi*epsilon0)  # Eq 39
        D_s = dDispersionEqn_dx_n0(RootEqn, a, b, epsilon, mu_r)
        logger.debug(f'D_s = {D_s}')
        Field2wake=1.0/qelec
        RootAmplit=a*RootEqn*P_rg(RootEqn,a,b,n)/D_s*NormalizationCGS2SI*Field2wake
        RootWavVec=RootEqn/(math.sqrt(epsilon*mu_r-1.0))

    if n>0:
        delta= a-b #relative to delta?
        while (np.any(np.abs(DispersionEqn(RootEqn+delta,a,b,n, epsilon, mu_r)-DispersionEqn(RootEqn,a,b,n, epsilon, mu_r)))>1.):
            delta=delta/10.0

        logger.debug(f'step for derivative = {delta}')
        #NormalizationCGS2SI=qelec*qelec/(a*a)/(4*math.pi*epsilon0)  # Eq 5.3
        D_s=(DispersionEqn(RootEqn+delta,a,b,n, epsilon, mu_r)- \
             DispersionEqn(RootEqn,a,b,n, epsilon, mu_r))/delta/a # d/dx D(x) in Equation (4.12) in CGS
        Field2wake=1.0/qelec
        RootAmplit=a*RootEqn*P_rg(RootEqn,a,b,n)*R_rg(RootEqn,a,b,n)/D_s*Field2wake
        RootWavVec=RootEqn/(np.sqrt(epsilon*mu_r-1.0))

# renormalize the field to the charge so units are now V/m/NC

    summary_str = \
 f"""
    ----- Summary ------
    mode order ={n}
    Roots:{RootEqn}
    dDisp:{D_s}
    Mode Amplitudes:{RootAmplit}
    Mode WaveVectors / m^-1 :{RootWavVec}
    Mode Wavelengths /m :{2.*np.pi/RootWavVec}
    --------------------"""
    logger.info(summary_str)
    return(RootAmplit, RootWavVec)


# ---------------- definition of distributions ----------------------------------------------
#            the distribution are normalized to 1 C
#            so they need to be multiplied by the charge in the main program

def BunchDistG(zz, sigmaz):
    """
    Gaussian Bunch Distribution with RMS width `sigmaz` centered around 0, evaluated at `zz`.

    Parameters
    ----------
    zz : array_like
        Points in meters at which to evaluate the distribution.
    sigmaz : float
        RMS width of the bunch distribution in meters. This is the scale of the Gaussian.

    Returns
    -------
    array_like:
        normalized weights of the distribution at the points `zz`.
    """
    sigmat=sigmaz/cms
    return 1./np.sqrt(2.*np.pi*sigmat**2)  *  np.exp(-zz**2/(2.*(cms*sigmat)**2))

def BunchDistU(zz,sigmaz):
    """
    Uniform Bunch Distribution with total hard-edge width 2*`sigmaz` centered around 0, evaluated at `zz`.

    Parameters
    ----------
    zz : array_like
        Points in meters at which to evaluate the distribution.
    sigmaz : float
        half hard-edge width of the bunch distribution in meters. This is the scale of distribution.

    Returns
    -------
    array_like:
        normalized weights of the distribution at the points `zz`.
    """
    sigmat=sigmaz/cms
    zzz=zz-0.*sigmaz
    return (1./(2.*sigmat)*
           (-heaviside(zzz-sigmaz)+heaviside(zzz+sigmaz)))

def BunchDistL(zz,sigmaz):
    """
    Linearly-ramped Bunch Distribution with total hard-edge width 2*`sigmaz` centered around 0, evaluated at `zz`.

    Parameters
    ----------
    zz : array_like
        Points in meters at which to evaluate the distribution.
    sigmaz : float
        half hard-edge width of the bunch distribution in meters. This is the scale of distribution.

    Returns
    -------
    array_like:
        normalized weights of the distribution at the points `zz`.
    """
   #  sigmat=sigmaz/cms
    zzz=zz-0.*sigmaz
    return (cms/(2.*sigmaz**2)*(zzz+sigmaz)*
           (-heaviside(zzz-sigmaz)+heaviside(zzz+sigmaz)))

def BunchDistP(zz, sigmaz):
    """
    Parabolic Bunch Distribution with total hard-edge width 2*`sigmaz` centered around 0, evaluated at `zz`.

    Parameters
    ----------
    zz : array_like
        Points in meters at which to evaluate the distribution.
    sigmaz : float
        half hard-edge width of the bunch distribution in meters. This is the scale of distribution.

    Returns
    -------
    array_like:
        normalized weights of the distribution at the points `zz`.
    """
   #  sigmat=sigmaz/cms
    zzz=zz-0.*sigmaz
    return (1.*cms/(4./3.*sigmaz**3)*(sigmaz**2-zzz**2)*
           (-heaviside(zzz-sigmaz)+heaviside(zzz+sigmaz)))


def BunchDistrib (Distribution, zz, sigmaz):
    '''
    Evaluate the distribution with characteristic size `sigmaz` at the positions `zz`, centered around zero.
    The meaning of `sigmaz` depends on the chosen distribution. `zz` and `sigmaz` must have the same dimensions.
    The unit of the distribution is 1/second, such that it is normalized to
    ```
    np.sum(MyBunch *dz/c) == 1
    ```

    Parameters
    ----------
    Distribution: callable
        distribution function accepting two arguments, where the first one refers to the positions and
        the second one refers to characteristic size. Currently, the module supports `BunchDistG`,
        `BunchDistU`, `BunchDistL` and `BunchDistP`.
    zz: array_like
        array of shape (N,) of positions.
    sigmaz: float
        parameter characterizing the size of the distribution.

    Returns
    -------
    '''
# evaluation of function "Distribution" with length sigmaz on array zz
#     mainly for plotting prupose
    zzmean=np.mean(zz)
    zeval=zz-zzmean
    dz=np.abs(zz[1]-zz[0])
    logger.debug(("dz:", dz))
    MyBunch=Distribution(zeval,sigmaz)
    return zeval, MyBunch


#---------- for tracking purpose

def Long_GreenFunction(RootAmplit, RootWavVec, r0,r, b, a, n, zmin, zmax,Nz,epsilon, mu_r=1):
    """
    Single-mode Green's function of the longitudinal wake field, i.e. the wake field function. For the wake
    potential it most be convolved with the longitudinal charge distribution and integrated along the z-coordinate.

    Parameters
    ----------
    RootAmplit: ndarray of float,
        One-dimensional array of shape (Nmode,) containing the amplitudes of the radial mode expansion for modes 
        with a polar mode index `n`. Those are the amplitudes in the sums of Eq. (4.11) or (3.9), resp.
        Retrieved from `FindModes`.
    RootWavVec:
        One-dimensional array of shape (Nmode,) containing the wavevectors of the radial mode expansion for modes 
        with a polar mode index `n`. See sums of Eq. (4.11) or (3.9), resp.
        Retrieved from `FindModes`.
    r0: float
        Position of the leading source particle (Dirac delta function) exciting the wake. In meter.
    r: float
        Position of the trailing particle (Dirac delta function) on which the wake acts. In meter.
    b: float
        Inner radius of the dielectric lined waveguide. In meter
    a: float
        Outer radius of the dielectric lined waveguide. In meter
    n: int
        polar mode index
    zmin: float
        lower boundary of the z-interval over which the wake is computed
    zmax: float
        upper boundary of the z-interval over which the wake is computed
    Nz: integer
        number of sampling points of the z-interval over which the wake is computed 
    epsilon: float
        dielectric permittivity of the dielectric lined waveguide 
    mu_r: float, optional
        permeability of the dielectric lined waveguide.
    
    """
    assert b < a   # Ensure that the order of inner and outer radii is correct
    assert r0 < b and r < b
    zz=np.linspace(zmin,zmax,Nz)
    WakeGreen=0.0*zz
    Nmode=len(RootAmplit)
    for i in range(Nmode):
        logger.debug(('LGreen:', i, RootAmplit, RootWavVec))
        if n>0:
            NormalizationCGS2SI=8.0*qelec/(a*a)*(r0*r/(b*b))**n/(4*math.pi*epsilon0)
            WakeGreen=WakeGreen+NormalizationCGS2SI*RootAmplit[i]*np.cos(RootWavVec[i]*zz)#equation (4.11) in Ng paper
        if n==0:
            WakeGreen= WakeGreen + RootAmplit[i]*np.cos(RootWavVec[i]*zz)  #equation (3.9) in Ng paper
    return zz, WakeGreen

def multimode_long_greens_function(RootAmplit, RootWavVec, m_modes, r0, r, theta, zmin, zmax, Nz, b, a, epsilon, mu_r=1):
    """
    Longitudinal multimode Green's function, i. e. the wake field function (in V/ ( m C )  ).

    Parameters
    ----------
    RootAmplit: ndarray of float,
        Two-dimensional array of shape (Mmode, Nmode) containing the amplitudes of the mode expansion. The row index 
        iterates to the polar mode index, the column index iterates the radial mode index. For each polar mode index 
        the amplitudes and wavevectors are computed with `FindModes`. The radial mode expansion in the radial index is given 
        in Eq. (4.11) or (3.9), resp. The polar mode expansion is given in Eq. (2.9).
    RootWavVec:
        Two-dimensional array of shape (Mmode, Nmode) containing the wave vectors of the mode expansion
    m_modes: ndarray of type int
        One-dimensional array of shape (Mmode,), containing the polar mode indices taken into account in the
        polar mode expansion, Eq. (2.9).
    r0: float
        Position of the leading source particle (Dirac delta function) exciting the wake. In meter.
    r: float
        Position of the trailing particle (Dirac delta function) on which the wake acts. In meter.
    theta: float
        Polar angle between leading and trailing particle in the cross-sectional plane. In radian.
    zmin: float
        lower boundary of the z-interval over which the wake is computed
    zmax: float
        upper boundary of the z-interval over which the wake is computed
    Nz: integer
        number of sampling points of the z-interval over which the wake is computed 
    b: float
        Inner radius of the dielectric lined waveguide. In meter.
    a: float
        Outer radius of the dielectric lined waveguide. In meter.
    epsilon: float
        dielectric permittivity of the dielectric lined waveguide 
    mu_r: float, optional
        permeability of the dielectric lined waveguide.
    """
    z_wake_func = np.zeros(Nz)
    long_wake_func = np.zeros(Nz)

    kwargs= {
                    "b": b, "a" : a, "zmin": zmin, "zmax" : zmax,
                    "Nz": Nz, "epsilon": epsilon
                } 
    for i, m in enumerate(m_modes):        
        kwargs["RootAmplit"] = RootAmplit[i]
        kwargs["RootWavVec"] = RootWavVec[i]
        kwargs["n"] = m
        kwargs["r0"] = r0
        kwargs["r"] = r
        logger.debug(f"Computing longitudinal wake potential function between z=0 and z={zmax}\n"+
                     f"with {m_modes.shape[0]} angular modes.")
        z_wake_func_temp, long_wake_func_temp = Long_GreenFunction(**kwargs)
        long_wake_func_temp *= np.cos(m * theta) 
        z_wake_func[:] = z_wake_func_temp
        long_wake_func[:] += long_wake_func_temp

    return z_wake_func, long_wake_func
    # wake potential function from wake field function
    pot_long_func = long_wake_func*(l_dlw - z_wake_func )
    # wake potential from wake potential function and bunch profile
    z_wake, long_wake_pot = WakePotential(BunchDistG, pot_long_func, z_wake_func, σ_z)
    return np.vstack((z_wake, long_wake_pot*1e-12))


def Trans_GreenFunction(RootAmplit, RootWavVec, r0,r, b, a, n, zmin, zmax,Nz,epsilon, mu_r=1):
    '''
    mu_r: float
       relative permeability of the dielectric
    '''
    assert b < a   # Ensure that the order of inner and outer radii is correct
    assert r0 < b and r < b
    zz=np.linspace(zmin,zmax,Nz)
    WakeGreen=0.0*zz
    Nmode=len(RootAmplit)
    F_=np.zeros(Nmode)
    for i in range(Nmode):
        logger.debug((i, RootAmplit[i], RootWavVec[i]))
        #P=P_rg(RootEqn[i],a,b,n)
        #R=R_rg(RootEqn[i],a,b,n)
        NormalizationCGS2SI=qelec*qelec/(a*a)*  (r0/a)**n  *(r/a)**(n-1)  *8.0*n*np.sqrt(epsilon-1.0)/(b/a)**(2.0*n) /(4*math.pi*epsilon0)
        #F_[i]=8.0*n*np.sqrt(epsilon-1.)/((b/a)**(2*n))*P*R/D_s[i]*a        #equation (5.4)
        logger.debug(f"NormalizationCGS2SI = {NormalizationCGS2SI}")
        WakeGreen=WakeGreen+RootAmplit[i]/(RootWavVec[i]*np.sqrt(epsilon*mu_r-1.0)*a) * np.sin(RootWavVec[i]*zz)*NormalizationCGS2SI #equation (5.3), RootWavVec=RootEqn/(np.sqrt(epsilon*mu-1.0)) RootAmplit=a*RootEqn*P_rg(RootEqn,a,b,n)*R_rg(RootEqn,a,b,n)/D_s*Field2wake
    WakeGreen=WakeGreen/qelec#missing mu here ,different from equation (4.12)?
    return zz, WakeGreen

def multimode_trans_greens_function(RootAmplit, RootWavVec, m_modes, r0, r, theta, zmin, zmax, Nz, b, a, epsilon, mu_r=1):
    """
    Transverse multimode Green's function, i. e. the wake field function (in V/ ( m C )  ), depending on the longitudinal position. Note that the transverse wake function is a vector in cylindrical coordinates. The returned array is two-dimensional, where the first row is the radial component, and the second row is the polar component.

    Parameters
    ----------
    RootAmplit: ndarray of float,
        Two-dimensional array of shape (Mmode, Nmode) containing the amplitudes of the mode expansion. The row index 
        iterates to the polar mode index, the column index iterates the radial mode index. For each polar mode index 
        the amplitudes and wavevectors are computed with `FindModes`. The radial mode expansion in the radial index is given 
        in Eq. (4.11) or (3.9), resp. The polar mode expansion is given in Eq. (2.9).
    RootWavVec:
        Two-dimensional array of shape (Mmode, Nmode) containing the wave vectors of the mode expansion
    m_modes: ndarray of type int
        One-dimensional array of shape (Mmode,), containing the polar mode indices taken into account in the
        polar mode expansion, Eq. (2.9).
    r0: float
        Position of the leading source particle (Dirac delta function) exciting the wake. In meter.
    r: float
        Position of the trailing particle (Dirac delta function) on which the wake acts. In meter.
    theta: float
        Polar angle between leading and trailing particle in the cross-sectional plane. In radian.
    zmin: float
        lower boundary of the z-interval over which the wake is computed
    zmax: float
        upper boundary of the z-interval over which the wake is computed
    Nz: integer
        number of sampling points of the z-interval over which the wake is computed 
    b: float
        Inner radius of the dielectric lined waveguide. In meter.
    a: float
        Outer radius of the dielectric lined waveguide. In meter.
    epsilon: float
        dielectric permittivity of the dielectric lined waveguide 
    mu_r: float, optional
        permeability of the dielectric lined waveguide.
    """
    z_wake_func = np.zeros(Nz)
    trans_wake_func = np.zeros(Nz)

    kwargs= {
                    "b": b, "a" : a, "zmin": zmin, "zmax" : zmax,
                    "Nz": Nz, "epsilon": epsilon, "mu_r": mu_r 
                } 
    for i, m in enumerate(m_modes):        
        kwargs["RootAmplit"] = RootAmplit[i]
        kwargs["RootWavVec"] = RootWavVec[i]
        kwargs["n"] = m
        kwargs["r0"] = r0
        kwargs["r"] = r
        logger.debug(f"Computing longitudinal wake potential function between z=0 and z={zmax}\n"+
                     f"with {m_modes.shape[0]} angular modes.")
        z_wake_func_temp, trans_wake_func_temp = Trans_GreenFunction(**kwargs)
        trans_wake_func_temp *= np.cos(m * theta) 
        z_wake_func[:] = z_wake_func_temp
        trans_wake_func[:] += trans_wake_func_temp

    return z_wake_func, trans_wake_func

def WakePotential (Distribution, WakeGreen, zz, sigmaz):
    zzmean=np.mean(zz)
    zeval=zz-zzmean
    dz=np.abs(zz[1]-zz[0])
    logger.debug(("dz:", dz))
    MyBunch=Distribution(zeval,sigmaz)
    WakePot=np.convolve(MyBunch,WakeGreen)*dz/cms
    WakePot=WakePot[0:len(zeval)]
    return zeval, WakePot

def MonteCarlo(Distribution, sigmaz, Nsample):
#   print Distribution(0)
# we assume the function is max at 0 for now // improve this later
    MyMax=Distribution(0,sigmaz)
    index=0
    Zarray=np.zeros(Nsample)
    while index<Nsample:
        zt=random.uniform(-4*sigmaz,4*sigmaz)
        yt=random.uniform(0, MyMax)
        if yt < Distribution(zt,sigmaz):
            Zarray[index]=zt
            index=index+1
    return Zarray
