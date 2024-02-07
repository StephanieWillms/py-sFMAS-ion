import sys; sys.path.append('../')
import numpy as np
import numpy.fft as nfft
from ini import *
from mediumDispersion import MediumDispersion 
import math


def f_rho(rho,W,t,dt):
    ''' 
    Function returning electron density.
    Function describes static model for tunneling ionization. Normalized!
    nitially, rho[0]=0. If everything is ionized: rho=1.0.

    Ref.: Landau, Lifshitz, Quantum Mechanics, 2nd ed. p. 276

    Args.: 
        rho (array): initial ionization units of full density
        W (array): ionization rate
        t (array): time array
        dt (float): time step

    Returns: (array): electron density rho        
    '''
    for i in range(1,len(t)):     
        im1 = i-1
        rho[i] = rho[im1] + dt*W[im1]*(1.0-rho[im1])
        if rho[i]>=1.0:
            rho[i] = 1.0
    return rho

def fi_st(x,U):
    ''' Function returning ionization rate.

    ef.: Babushkin, Opt. Express 18, 9660 (2010)
    Args.:
    eps: arbitrary constant, field in a.u.
    x (array): field in a.u
    U (float): ionization potential, eV

    Returns: ionization rate in fs**-1
    '''    
    # selecting strong field
    eps = 1e-5 #(small enough for non-ionization)
    ax = abs(x)
    f = np.zeros(ax.shape)
    ind = ax<eps
    f[ind] = 0.0
    ind = ax>=eps # where field is strong
    # calculating coeffs
    alpha_st,beta_st=tunnel_ion_coefs(U)
    # ionization rate main formula
    f[ind] =  (alpha_st/ax[ind])*(np.exp(-beta_st/(ax[ind])))
    return f
    

def tunnel_ion_coefs(U_ion):
    """ Defines tunnel ionization coefficients alpha, beta
    
        Args.: 
            U_ion (int): ionization energy of gas ion
            
        Returns.:
            beta_st (int): E[au]
            alpha_st (int): in fs^-1
    """
    U_H = 13.6                  # ionization potential of hydrogen in eV
    omega_a = 4.13e+1           # fs-1
    beta_st =  (2.0/3.0)*(U_ion/U_H)**(3./2) 
    alpha_st =  4.0*omega_a*(U_ion/U_H)**(5./2)  
    return (alpha_st,beta_st)
    

def convert_int_tw_to_field_au(inten):
    """converts field [W/um^2] in a.u.
        3.5*10^16 W/cm^2 =1 a.u. -> 3.4*10^8 W/um^2
        
        Args.:
        	inten: field in intensity
        	
        Ref.:  Taken from I.Babushkin code
        
        RETURNS:
        	returns field in a.u.        
    """
    return np.sqrt(inten/100.)*0.053
    

def time_deriv5(s,dt):
    """ Numerical differentiation in time domain
        
        Args.: 
        	s (array): field^3 - Kerr term
        	dt (int): time increment 
        	
        Ref.: https://en.wikipedia.org/wiki/Numerical_differentiation
    """
    s1 = np.roll(s,-1)
    s2 = np.roll(s1,-1)
    sm1 = np.roll(s,1)
    sm2 = np.roll(sm1,1)
    return 1./(12.*dt)*(-s2+8.*s1-8.*sm1+sm2)

def fdz(Et, gamma, t, w, ioni):
    """ Todo 
    """
    U_ion = 24.4
    j_norm = 0.279			#/np.sqrt(refractive_index_0) 
    					# refractive_index_0 to be seto above to !=1 for solids
    j_norm = j_norm*pre 		# pressure correction
    j_loss_norm = 0.21			#/np.sqrt(refractive_index_0)
    j_loss_norm = j_loss_norm*(2.0*U_ion) #(2.0*U_ion/refractive_index_0) 
    					# taking another constant coefficient to be a fixed multiplier
    j_loss_norm = j_loss_norm*5 	# pressure correction
    gamma_j = 1./190.*pre 		# 1/fs -- plasma current decay
    field = (np.real(Et)).reshape(t.shape)  #  E=real(analyt signal)
    field2 = field**2
    afield_au = convert_int_tw_to_field_au(field2)   
    W = fi_st(afield_au*1.5,U_ion).reshape(t.shape)
    rho =np.zeros(t.shape,dtype=float) # ionization array
    rho[0] = 0.0 			# initial condition for integration
    rho = f_rho(rho=rho,W=W,t=t,dt=dt)
    s2 = Et*np.conj(Et)
    s3 = Et*s2
    ds3 = time_deriv5(s3,dt)
    Kerr= (-1.0j*1j*gamma*dz)*ds3 	#nfft.ifft(c2)*Et*np.abs(Et)**2  
    dj = (rho*field).reshape(t.shape)      
    jrf = np.fft.fft(dj)
    jr =np.fft.ifft(jrf/(1.0j*w + gamma_j))
    safe_inv = lambda x: np.sign(x)/((np.abs(x)+1e-13))   # safe version of division 
    jloss = field*W*(1.0-rho)*safe_inv(abs(field2)) # all const factors are in _loss_norm          
    rj = -(j_norm*jr + j_loss_norm*jloss)          
    ## puttig all together
    r_bound =  Kerr.copy()
    if ioni in 'on':
        r =  r_bound+rj
    elif ioni in 'off':
        r =  r_bound	
    else:
        print ('IONIZATION: Just the values ''on'' and ''off'' are allowed for ionization')
    return r 


    return fdz(Et)
     
