"""splitStepSolver.py

Module file implementing split step solver for the simplified
forward model for the analytic signal (sFMAS) discussed in Ref. [1] (see Eq.
(18)). Following a divide-and-conquer strategy, the solver work in the 
Fourier domain where it is possible to solve the linear part exactly. The 
nonlinear part is extrapolated using a Runge-Kutta extrapolation scheme.

Refs:
    [1] Supercontinuum generation by multiple scatterings at a group 
        velocity horizon
        Demircan, A. and Amiranashvili, S. and Bree, C. and Morgner, U. and
        Steinmeyer, G.
        Optics Express, 22 (2014) 3866

"""
import numpy as np
import numpy.fft as nfft
import nonLin as ioni
from nonLin import *
#import math

def sFMASDQRK4_2(t,w,Eps0_w,dz,dt,Nz, wMin, wMax, mProp,ioni,ionization,vBoost,callbackFunc,n2,p_he):
    """Divide-and-conquer based solver for the simple field-based FMAS 

    Implements divide-and-conquer strategy to solve simplified field-based
    forward model for the analytic signal (FMAS) discussed in Ref. [1] (see Eq.
    (18)).  Signal propagation is performed in Fourier domain. While the linear
    part ist solved exactly, the nonlinear part is solved via a fourth-order
    Runge-Kutta method of fixed step-size, see chapter 16.1 of Ref. [2] (or
    Chapter II.1 of Ref. [3]).

    Args:
        t (numpy array, ndim=1): time axis
        w (numpy array, ndim=2): angular frequency axis
        Eps0_w (numpy array, ndim=2): analytic signal in Fourier domain
        dz (float): increment for z-axis
        Nz (int): number of propagation steps
        mProp (object): data structure handling refractive index and 
            derived quantities
        vBoost (float): velocity of moving reference frame
        callbackFunc (object): callback function handling measurements
        n2 (float): n2-value characterizing the fiber in units of microns^2/W                    

        callback function (callee callbackFunc) takes 5 parameters in the form
        callbackFunc(n, z, t, w, u), where:

        n (int): current iteration step
        z (float): current z-position 
        t (numpy array, ndim=1): t-mesh
        w (numpy array, ndim=1): frequency-mesh
        u (numpy array, ndim=1): signal in Fourier domain 
    
    Refs:
        [1] Supercontinuum generation by multiple scatterings at a group 
            velocity horizon
            Demircan, A. and Amiranashvili, S. and Bree, C. and Morgner, U. and
            Steinmeyer, G.
            Optics Express, 22 (2014) 3866

        [2] Numerical Recipes in Fortran 77
            Press, WH and Teukolsky, SA and Vetterling, WT and Flannery, BP
            Cambridge University Press (2nd Edition, 2002)

        [3] Solving Ordinary Differential Equations I: Nonstiff Problems
            Hairer, E. and Norsett, S. P. and Wanner, G.
            Springer Series in Computational Mathematics, 3rd Edition, 2008

	[4] Babushkin, Opt. Express 18, 9660 (2010)
    """

    # CONVENIENT INITIALIZATION AND DECLARATION OF AUXILIARY OBJECTS ##########
    betaDS =  mProp.beta(w)- w/vBoost#+ 1.0j*mProp.loss(w)  # profile of propagation const.
                                                            # including loss
    c0 = mProp.c0                                           # free-space speed of light 
    wMask = np.logical_and(w>0,w<w.max()*0.65)              # mask for dealiasing
    c1 = np.where(wMask, np.exp(0.5j*betaDS*dz), 0j)        # coeff-array of lin. op.
    c2 = np.where(wMask, 1.0j*w*dz*n2/c0, 0j)               # coeff-array of nonlin. op.
    u = np.where(wMask, Eps0_w, 0j)                         # dealiased initial condition
    gamma=n2/(c0)

    def _RK4Extrap(func,u):
        """helper function implementing custom 4th order Runge-Kutta solver

        Implements fourth-order Runge-Kutta scheme using the increments 
        presented in 16.1 of Ref. [1]

        Args:
            func (object): function implementing nonlinear operator
            u (numpy array): signal in Fourier-domain 

        Returns:
            uNew (numpy array): extrapolated Fourier-domain signal  

        Refs:
            [1] Numerical Recipes in Fortran 77
                Press, WH and Teukolsky, SA and Vetterling, WT and Flannery, BP
                Cambridge University Press (2nd Edition, 2002)
        """
        k1 = func(u) 
        k2 = func(u+k1/2)
        k3 = func(u+k2/2)
        k4 = func(u+k3)
        return u + k1/6 + k2/3 + k3/3 + k4/6

    def rk4_step(f,s):
        
        ds1=f(s)
        ds2=f(s+0.5*ds1)
        ds3=f(s+0.5*ds2)
        ds4=f(s+ds3)
        return s + (1./6.)*(ds1+2.*ds2+2.*ds3+ds4)


    

    # LINEAR PART (HALF-STEP) - 2ND TERM OF LHS IN EQ. (33) OF REF. [1]
    _linPart = lambda u: c1*u    
    
    # TIP: Code runs slow with "fdz" Module. 
    # If just Kerr effect is required, commented line runs much faster.
  
   # field = lambda E: nfft.fft(fdz(E,gamma, t, w, ionization))
    field = lambda E: c2*nfft.fft(E*np.abs(E)**2)
    
 
    def _nlinPart(u):
        return field(nfft.ifft(u))
    

    callbackFunc(0,0,t,w,u)
    for zi in range(1,Nz):
        u = _linPart(u)               # LINEAR PART - HALF STEP 
        u = _RK4Extrap(_nlinPart, u)  # NONLINEAR PART - FULL STEP
        u = _linPart(u)               # LINEAR PART - HALF STEP 
        callbackFunc(zi,zi*dz,t,w,u)   
	
    
    return u




# EOF: splitStepSolverFMAS.py 
