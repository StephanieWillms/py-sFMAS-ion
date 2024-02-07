""" mediumDispersion.py

Module implementing four different refractive index profiles
and class data structure providing convenient access to derived properties        
"""
import numpy as np
import scipy.optimize as so
from ini import *
import matplotlib.pyplot as plt
from numpy.polynomial import polynomial as P
from scipy import special
import csv

def _derivative(f, x):
    """helper function  
    
    Implements first derivative via finite differences

    Args:
        f (object) callable function
        x (float) argument at which first derivative should be calculated

    Returns:
        dfdx (float) first derivative via finite differences
    """
    h = 4.0e-3 
    xp = x + h
    xm = x - h
    dx = xp - xm
    return (f(xp)-f(xm))/dx
    

def refractiveIndex_sinusiod_square():
    """refractive index profile

    Implements polynomial function exhibiting several regions of anomalous dispersion, due to
    sinus-like function around expansion frequency w0=3 rad/fs.
    It is manually adapted to have similar values for beta2 as refractive index profile
    " WD2017 ". 

    Returns:
        ratFunc (function) rational function representing polynom

    Refs:
        [1] PhD Thesis, S.Willms Ch. "Periodic dispersion characteristics "
    """    
    return lambda x: (0.3/x)*(0.1*(1827414433229931786240000000*x**10-54822432996897953587200000000*x**9+733524153498494618996736000000*x**8-
5762934156633912881086464000000*x**7+29441225405404511737389711360000*x**6-102221820850123129221967380480000*x**5+
244429405096684843070613946368000*x**4-397809187176764794603251695616000*x**3+422208859892562485554534105743360*x**2-
263760816265074052677065796747264*x+2384185791015625*(x-3)**30+82969665527343750*(x-3)**28+2509002685546875000*
(x-3)**26-65234069824218750000*(x-3)**24+1440368261718750000000*(x-3)**22-26618005476562500000000*(x-3)**20+
404593683243750000000000*(x-3)**18-4952226682903500000000000*(x-3)**16+47541376155873600000000000*(x-3)**14-
346101218414759808000000000*(x-3)**12)/33952366055960455505447485440)+5 


def refractiveIndex_WD2017():
    """refractive index profile

    Implements rational Pade-approximant of order [3/3] for the medium 
    refractive index as used by Ref. [1] 

    Returns:
        ratFunc (function) rational function representing Pade-approximant

    Refs:
        [1] Soliton-Soliton Kollision am optischen Ereignishorizont
            Willms, S.
            Bachelor-Thesis, LUH (IQO), 2016 
            
    Note:
        - beta_2 exhibits three roots 
    """
    pCoeffs = (9.653881, -39.738626, 16.8848987, -2.745456)
    qCoeffs = (1.000000,  -9.496406,  4.2206250, -0.703437)
    b2Roots = (1.51174072266, 2.511328125, 5.46130371094)
    return coeffs2RatFunc(pCoeffs,qCoeffs)


def benabid_multires(pre,R,st):
    """refractive index profile

    Implements refractive index as used by Ref. [1] 

    Returns:
        (function) rational function 

    Refs:
        [1] Analytic model for the complex effective index of the leaky modes of 
            tube-type anti-resonant hollow core fibers
            Matthias Zeisberger, Markus A. Schmidt
            Scientific Reports 7, 11761 (2017)            
    """
    _safe_inv = lambda x: np.sign(x)/(np.abs(x)+1e-13)        # Prevents division by zero

    def f_wavelength(omega):
        """ 
        Args.:
            angular frequency in rad/fs
            
        Returns:
            vacuum wavelength in um
        """
        return 2.*np.pi*0.3*_safe_inv(omega)

    def _nhe(omega):
        '''dispersion of Helium
        '''
        lambda_um = 2.*np.pi*0.3*_safe_inv(omega)        # wavelength in um
        lm2 = _safe_inv(lambda_um**2)
        return 1.0+pre*0.014700910*_safe_inv(423.98-lm2) # note here 1+delta n

    na = lambda omega: _nhe(omega)
    ng = 1.45                    # refractive index of glass
    br=special.jn_zeros(0,1)[-1] # Bessel function for single mode
    
    Phi =lambda omega: st * (2 *np.pi)*_safe_inv(f_wavelength(omega))*np.sqrt( np.abs(ng**2 - na(omega)**2))                   # Definition as in Ref. [1]
    
    nw_re =lambda omega: - (f_wavelength(omega)/(2*np.pi))**2 * ( (br**2)/(2*na(omega)*R**2) ) - (f_wavelength(omega)/(2*np.pi) )**3 * (br**2/(na(omega)**2*R**3)) * ((_safe_inv(np.tan(Phi(omega))))/1) * ( ((ng/na(omega))**2 + 1)/(2*np.sqrt(np.abs(((ng/na(omega))**2 - 1))) ))      # Definition real refractive index as in Ref. [1]


    alpha =lambda omega: (1.+_safe_inv(np.tan(Phi(omega))**2))/(np.abs(((ng/na(omega))**2 - 1)))*br**3*(f_wavelength(omega)/(2*np.pi))**3/na(omega)**3/R**4*((ng/na(omega))**2 + 1)/2.0
                                 # Definition loss as in Ref. [1]
    nw_im = lambda omega:alpha(omega)*0.3*_safe_inv(omega) # Resulting imaginary refractive index
    return lambda x: nw_re(x) + 1.0j*nw_im(x) + _nhe(x) 


def Kagome_PCF(pre, Re, t):
    """refractive index profile

    Implements refractive index as used by Ref. [1] 

    Returns:
        (function) rational function 

    Refs:
        [1] Influence of ionization on ultrafast gas-based nonlinear fiber optics
            W. Chang, A. Nazarkin, J. C. Travers, J. Nold, P. Hölzer,N. Y. Joly, 
            and P. St.J. Russell 
            Optics Express 19, 21018 (2011)  
            
        [2] The journal of physical chemistry 79 (1975)
    """
    _safe_inv = lambda x: np.sign(x)/((np.abs(x)+1e-13))
    
    def _nhe(omega):
        '''dispersion of Helium
        '''
       # p_he=5.0
        lambda_um = 2.*np.pi*0.3*_safe_inv(omega)        # wavelength in um
        lm2 = _safe_inv(lambda_um**2)
        return 1.0+pre*0.014700910*_safe_inv(423.98-lm2) # note here 1+delta n

    def _nar(omega):
        '''dispersion of Argon as defined in Ref. [2]    
        '''
        #pre = 4.0
        lambda_um = 2.*np.pi*0.3*_safe_inv(omega)
        lm2 = _safe_inv(lambda_um**2)
        return 1.0+pre*(2.e-6)*120.625*(14.0*_safe_inv(87.892 - lm2) + 14.00*_safe_inv(91.011 - lm2) + 430.630*_safe_inv(217.6 - lm2))  # note here 1+delta n

    br=special.jn_zeros(0,1)[-1] # Bessel function for single mode
    diam=2*Re                      # fiber diameter in um
    Aeff=np.pi*(diam/2.0)**2     # effective area in um^2
    c0 =0.3                       # speed of light
    nw_re =lambda omega: - 0.5*((br**2*c0**2)/(omega**2*(0.5*diam)**2)) # refractive index 
                                                                        # of waveguide
    return lambda x: nw_re(x) + _nar(x)
    

def coeffs2RatFunc(pCoeffs,qCoeffs):
    """coefficients-to-rational function

    Helper function converting numerator and denominator coefficients to 
    Pade-approximation of order [m/n]

    Args:
        pCoeffs (tuple) length m list containing numerator coefficients 
        qCoeffs (tuple) length n list containing denominator coefficients

    Returns:
        ratFunc (function) rational function representing Pade-approximation of          
                 order [m/n]

    Note: 
        - uses numpys polynomial class poly1d to construct numerator and 
          denominator polynoms
        - numpys poly1d expects coefficients to be provided in reverse order
    """
    p = np.poly1d(pCoeffs[::-1])
    q = np.poly1d(qCoeffs[::-1])
    return lambda x: p(x)/q(x)


def fetchProfile(mode='Zhang', pre=3, Re=20, t=0.8):
    """refractive index profile handler

    Hack emulating switch statement of other languages for easy access to 
    implemented refractive index profiles

    Args:
        mode (str) type of refractive index profile 
                   (choices: WD2017, NLPM750, ESM, AD2010; default: WD2017)

    Returns:
        n (object) Pade approximant of refractive index profile 
    """
    return {'WD2017': refractiveIndex_WD2017(), 'kagome':Kagome_PCF(pre, Re, t),'sinusoid':refractiveIndex_sinusiod_square(),
       'Benabid': benabid_multires(pre, Re, t)
              }.get(mode)


def boundedIndex(nRefIdx, xMin, xMax):
    """bounding for refractive index

    Args:
        nRefIdx (object) function returning refractive index values
        xMin (float) lower bound for refractive index profile
        xMax (float) upper bound for refractive index profile
    
    Returns:
        nBound (object) vectorized numpy function returning bounded refractive
            index values

    Notes:
        - refractive index profile n(x) is set to n(xMin) for all x<xMin and to
          n(xMax) for all x>xMax
    """

    def _nBounded(x):
        """helper function """
      #  np.where(x < xMin, x, nRefIdx(xMin))
        if x < xMin :
           val = nRefIdx(xMin)
        elif x > xMax:
           val = nRefIdx(xMax)
        else: 
           val = nRefIdx(x)
        return val
    return np.vectorize(_nBounded) 


class MediumDispersion(object):
    """Medium dispersion 

    Implements convenient access to refractive index and derived properties
    """
    def __init__(self,profType="Benabid",wMin=0.0,wMax=40, c0=0.29979,pre=4.0, Re=1, t=0.1):
        """initializes instance of medium dispersion class 

        Args:
            profType (str) refractive index profile type (default: WD2017)
            wMin (float) lower bound for refractive index profile
            wMax (float) upper bound for refractive index profile

        Attributes:
            c0 (float) free space velocity 
            n (func) refractive index profile
            R (int) fiber core radius
            st (float) fiber bridge diameter
            pre (int) pressure of gas in fiber in atm
        """
        self.c0 = c0
        self.nUnbound = fetchProfile(profType, pre,Re, t)
        self.n = boundedIndex(self.nUnbound,wMin,wMax)#
        
    def loss(self,w):
        """frequency dependent loss in 1/um
        Args:
            w (numpy array, float): angular frequency value/array
        Returns:
            im(beta) (numpy array, float): dispersion profile

        Remarks: 
            loss is calculated for abs(omega); That is,
            loss(omega)=loss(-omega) ensured.

        """
        return np.abs(w)*np.imag(self.n(np.abs(w)))/self.c0


    def beta(self,w):
        """frequency dependent propagation profile
        
        Args:
            w (numpy array/float) angular frequency value/array 

        Returns:
            beta (numpy array/float) dispersion profile
        """
        return w*np.real(self.n(np.abs(w)))/self.c0


    def beta1(self,w):
        """dispersion profile (order 1)

        Args:
            w (numpy array/float) angular frequency value/array 

        Returns:
            beta1 (numpy array/float) dispersion profile (order 1)
        """
        return _derivative(self.beta,w)
       

    def beta2(self,w):
        """dispersion profile (order 2)

        Args:
            w (numpy array/float) angular frequency value/array 

        Returns:
            beta1 (numpy array/float) dispersion profile (order 2)
        """
        return _derivative(self.beta1,w)
 

    def beta2Root(self,Min,wMax):
        """find bracketed root of 2nd order dispersion profile

        Note:
            - uses scipy.optimize.bisect for bracketed root finding

        Args:
            wMin (float) lower bound for root finding procedure
            wMax (float) upper bound for root finding procedure

        Returns:
            wc (float) root of 2nd order dispersion profile in bracketed
                 intervall 
        """
        return so.bisect(self.beta2, wMin, wMax)
        

    def vRef(self, p01_w0):
        return 1./self.beta1(p01_w0)

    def resonanceFrequency(self,w0,wMin,wMax):
        return so.minimize_scalar(lambda w: np.abs(self.beta1(w)-self.beta1(w0)), bounds=(wMin, wMax), method='bounded')



# EOF: mediumDispersion.py
