"""inputPulse.py

Module file implementing functions to determine input pulse properties and
routines to convert input pulse to initial condition for use with sFMAS 
propagation algorithms.

Refs:
    [1] Hamiltonian structure of propagation equations for ultrashort
        optical pulses
        Amirasnashvili, Sh. and Demircan, A.
        Phys. Rev. A, 82 (2010) 013812

    [2] Supercontinuum generation by multiple scatterings at a group 
        velocity horizon
        Demircan, A. and Amiranashvili, S. and Bree, C. and Morgner, U. and
        Steinmeyer, G.
        Optics Express, 22 (2014) 3866

"""
import numpy as np
import numpy.fft as nfft
from ini import pulseshape
from mediumDispersion import MediumDispersion 


def tREF_to_wAS(w,E):
    """real electric field (t-domain) to analytic signal (w-domain)

    Implements complexification of real electric field in time-domain to
    analytic signal in Fourier-domain following section II.B of Ref. [1]

    Args:
        w (numpy-array, ndim=1): angular frequency axis
        E (numpy-array, ndim=1): real electric field in time domain 

    Returns:
        Ew (numpy-array, ndim=1): Analytic signal in Fourier domain 

    Note:
        - computes standard discrete-time 'analytic' signal following 
        Eq. (8) in sect. IV of Ref. [2] via the two-step procedure:
        (1) Compute N-point DFT Ew[:] using FFT of original N=w.size samples:
            Ew = fft(E)
        (2) Construct N-point one-sided discrete-time 'analytic' singal trafo:
            Eps_w[0]           =   Ew[0]
            Eps_w[w.size/2]    =   Ew[w.size/2]
            Eps_w[1:w.size/2]  = 2*Ew[1:w.size/2]
            Eps_w[w.size/2+1:] = 0j
        - the above procedure ensures the following two properties:
        (1) real part of reverse transform yields the original (real) signal:
            E = np.real(np.nfft.ifft(Eps_w))
        (2) real and imaginary component are orthogonal over interval:
            E = np.real(np.nfft.ifft(Eps_w))
            EIm = np.imag(np.nfft.ifft(Eps_w))
            np.sum(E*EIm) < 1e-8
        see unit test module test_fieldRepresentations.py

    Refs:
        [1] Hamiltonian structure of propagation equations for ultrashort
            optical pulses
            Amirasnashvili, Sh. and Demircan, A.
            Phys. Rev. A, 82 (2010) 013812

        [2] Computing the Discrete-Time 'Analytic' Signal via FFT
            Marple, S. L. Jr.
            IEEE Trans. Sig. Proc., 47 (1999) 2600
    """
    # OBTAIN INITIAL CONDITION FROM INPUT PULSE, SEE REF [1], SECT. V.B
    # (1) DIRECTLY OBTAIN INITIAL CONDITION IN FOURIER-DOMAIN
    Eps_w = nfft.fft(E)
    # (2) CORRECT NEGATIVE- AND ZERO-FREQUENCY COMPONENTS 
   
    Eps_w[1:int(0.5*w.size)] = 2*Eps_w[1:int(0.5*w.size)]
    Eps_w[int(w.size*0.5)+1:] = 0j
    return Eps_w
    
def inputPulseAmplitude_sFMAS2(mProp,n2, t0,w0):
    """amplitude of initial pulse

    Determine amplitude of the initial pulse envelope from pulse width and
    central angular frequency for simplified forward model.

    Args:
        mProp (object): instance of MediumDispersion class
        n2 (float): nonlinear refraction  
        t0 (float): width of the initial pulse
        w0 (float): central anglular frequency of the initial pulse

    Returns:
        Psi0 (float): amplitude of the initial pulse envelope
    """
    return np.sqrt( np.abs(mProp.beta2(w0)))/t0/np.sqrt(n2*w0/(mProp.c0))

def inputPulseAmplitude_sFMAS(factor,Ep, t0,w0, aeff,mProp,n2,energy_stat):
    """amplitude of initial pulse
    # amplitude in sqrt(TW/cm^2)

    Determine amplitude of the initial pulse envelope from pulse width and
    central angular frequency for simplified forward model.

    Args:
        mProp (object): instance of MediumDispersion class
        n2 (float): nonlinear refraction  
        t0 (float): width of the initial pulse
        w0 (float): central anglular frequency of the initial pulse

    Returns:
        Psi0 (float): amplitude of the initial pulse envelope
    """
    if energy_stat in 'True':
        amp=np.sqrt((factor*(Ep/(t0*10**-15)))/(aeff*10**4)) # peak power[W]->[TW] # um->cm 
    else:
        amp=np.sqrt( np.abs(mProp.beta2(w0)))/t0/np.sqrt(n2*w0/(mProp.c0))
    return amp

#### SELECTION OF DIFFERENT PREPARED ANALYTICAL SIGNALS ############################################

def initialElectricField(Psi01, p01_t0, p01_w0):
    """initial electric field

    Implements initial electric field pulse of hyperbolic-secant shape 
    following section V.A "Simulations of Pulse Propagation - Numerical 
    Procedure" of Ref. [1].

    Args: 
        Psi0 (float): amplitude of the initial pulse envelope
        t0 (float): width of the initial pulse
        w0 (float): central anglular frequency of the initial pulse

    Returns:
        E0 (float): input pulse electric field with hyperbolic-secant shape

    Refs:
        [1] Hamiltonian structure of propagation equations for ultrashort
            optical pulses
            Amirasnashvili, Sh. and Demircan, A.
            Phys. Rev. A, 82 (2010) 013812

    """
    sech = lambda x: 1./np.cosh(x)
    tau=p01_t0*1.7627
    for pulseshape in 'sech':
        return lambda t: np.real(Psi01*sech(t/p01_t0)*np.exp(-1j*p01_w0*t))
    else:
        return lambda t: np.real(Psi01*np.exp(-t**2/(2*tau**2))*np.exp(-1j*p01_w0*t))   
 
    



def initialCondition(t,E0):
    """get initial condition from input pulse

    Obtain initial conditions associated with provided input pulse as
    discussed by Ref [1]

    Args: 
        t (numpy array) array representing the time domain
        E0 (numpy array) array containing input pulse

    Returns:
        (t, Eps0_t) (tuple, numpy arrays) time-axis and initial condition 
            in time domain
        (w, Eps0_w) (tuple, numpy arrays) frequency-axis and initial 
            condition in Fourier domain

    Refs:
        [1] Hamiltonian structure of propagation equations for ultrashort
            optical pulses
            Amirasnashvili, Sh. and Demircan, A.
            Phys. Rev. A, 82 (2010) 013812

        [2] Supercontinuum generation by multiple scatterings at a group
            velocity horizon
            Demircan, A. and Amiranashvili, S. and Bree, C. and Morgner, U.
            and Steinmeyer, G.
            Optics Express 22 (2014) 3866

    Note: 
        - forward/reverse FFT procedure is needed to filter out all negative
          frequency coefficients and to obtain proper initial conditions 
          for later propagation with FMAS (forward model for analytic signal) 
          propagator
        - albeit initial condition is composed of positive frequency terms
          only, see Eq. 16 of Ref. [2], the full Fourier-domain, i.e. including
          zero- and negative-frequency components (all set to zero) are used to
          represent the initial condition 
    """
    # Note: albeit the initial condition exhibits only positive frequencies
    # the full Fourier domain is used to represent the initial condition
    w = nfft.fftfreq(t.size,d=t[1]-t[0])*2*np.pi

    # obtain initial condition from input pulse, see Ref [1], sect. V.B
    # (1) directly obtain initial condition in Fourier-domain
    #Eps0_w = 2*nfft.fft(E0)
    # (2) correct negative- and zero-frequency components 
    #Eps0_w[w<=0]=0
    Eps0_w = tREF_to_wAS(w,E0)
    # (3) yield initial condition for the analytic signal in time domain.
    # Note: signal is just a sum of positive-frequency components, cf. Eq (16)
    # in Ref. [2]
    Eps0_t = nfft.ifft(Eps0_w)
    return (t,Eps0_t, w, Eps0_w)


def estimateReferenceGroupVelocity(t, mProp, Psi01, p01_t0, p01_w0, p02_tOff):
    """determine reference group velocity 

    Computes NUMERICAL estimate of the reference group velocity at which the
    reference system is supposed to move so as to match the speed of the 
    reference pulse

    Args:
        t (numpy array) array representing the time domain
        mProp (object): instance of MediumDispersion class
        Psi0 (float): amplitude of the initial pulse envelope
        t0 (float): width of the initial pulse
        w0 (float): central anglular frequency of the initial pulse
        tOff (float): time offset for preparation of the inital pulse

    Returns:
        vRef (float): numerical estimate of reference group velocity 

    Note:
        - numerical estimate of reference group velocity might differ
          slightly from the supposed reference velocity as directly computed
          from the central angular frequency of the reference pulse. However,
          for an optimal velocity-matching of the reference frame to the 
          reference pulse, the numerical estimate is most adequate.
    """
    E0 = initialElectricField(Psi01, p01_t0, p01_w0)
    (t, Eps0_t, w, Eps0_w) = initialCondition(t, E0(t+p02_tOff))
    # ESTIMATE VIA INTENSITY AVERAGED FREQUENCY 
    # OBSERVATION: UNDERESTIMATES REFERENCE VELOCITY 
    # I0 = np.abs(Eps0_w)**2
    # w0Ref = np.trapz(w*I0,x=w)/np.trapz(I0,x=w)
    # ESTIMATE VIA SIMPLE PEAK FREQUENCY OF INTENSITY PROFILE
    # OBSERVATION: MOST ADEQUATE ESTIMATE OF REFERENCE VELOCITY
    w0Ref = w[np.argmax(np.abs(Eps0_w)**2)]
    return mProp.vRef(w0Ref)

# EOF: inputPulse.py
