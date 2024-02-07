"""splitStepSolverFMAS.py

Minimalistic example to describe propagation in nonlinear medium including effect of ionization.
Example:
1	 can be used to reprodiuce example in PhD thesis of S. willms [CITE]
2	Chang paper

Refs:
    [1] Diss Stephanie

    [2] W. Chang, A. Nazarkin, J. C. Travers, J. Nold, P. Hölzer, N. Y. Joly, and P. S. J. Russell,
       “Influence of ionization on ultrafast gas-based nonlinear fiber optics,” Opt. Exp., vol. 19,
        p. 21018, 2011.

AUTHOR: S. Willms
"""

import sys; sys.path.append('./src')
from ini import *
from splitStepSolver import *
from inputPulse import *
from mediumDispersion import * 
from observables import *
from nonLin import *

def main():

    #' INITILAZE COMPUTATIONAL DOMAIN AND SYSTEM PARAMETRS ####################
    print ('initializing system ... ')
    print ('Pulsetype:' + pulseshape)
    print ('Pulse energy:' + str(Energ) + 'J')
    print ('Ionization:' + str(ionization) )
    print ('------ Computational Domain -------')
    print ('tMax:' + str(tMax))
    print ('Nt:' + str(Nt))
    print ('zMax:' + str(zMax))
    print ('Nz:' + str(Nz))
    print ('------ Fiber Properties -------')
    print ('Fiber Radius:' + str(R) + 'um')
    print ('strut:' + str(st) + 'um')
    print ('Gas pressure:' + str(pre) + 'bar')
    print ('... finished')
   
   
    ## PULSE CHARACTERIZATION ################################################    
    
    """ Pulse 1 """
    p01_t0 = float(sys.argv[1])                    # pulsewidth
    p01_w0   = float(sys.argv[2])                  # centralfrequency
    p01_N    = float(sys.argv[3])                  # soliton order                 
    p01_tOff = 0.00                                # temporal offset between p01 and p02
   
    """ Pulse 2 """
    p02_t0 = float(sys.argv[4])                    # pulsewidth
    p02_w0   = float(sys.argv[5])                  # centralfrequency    
    p02_N    = float(sys.argv[6])                  # order
    p02_tOff = float(sys.argv[7])                  # temporal offset between p01 and p02
    
    p_n = 'one'                                    # number of pulses to propagate

    print ('------ Parameters -------')
    print ('Nonlinearity:' + str(n2) + 'um^2/W')
    print ('Involved Pulses:' + str(p_n) )
    print ('--- Pulse #1 ---')
    print ('duration:' + str(p01_t0) + 'fs')
    print ('frequency:' + str(p01_w0) + 'rad/fs')
    print ('order:' + str(p01_N) )
    print ('--- Pulse #2 ---')
    print ('duration:' + str(p02_t0) + 'fs')
    print ('frequency:' + str(p02_w0) + 'rad/fs')
    print ('order:' + str(p02_N) )
    ## GIVE NAME TO FILE #####################################################
    
    """ Stores important parameters in file-output Name """
    fdBase = '_%s_Nt%d_Nz%d'%('sFMAS-RK4-da',Nt,Nz) # custom file descriptor    
   
    ## SET AUXILIARY DATA STRUCTURES ##########################################  

    t = np.linspace( -tMax, tMax, Nt, endpoint=False) # SET TIME-AXIS 
    mProp = MediumDispersion( "kagome", wMinFiber,wMaxFiber, c0,pre=4.0, Re=R, t=st) # FIBER PROPS
    

    ## SET INITIAL CONDITION ##################################################   
    
    E0 = np.zeros(t.size, dtype=complex)   
    
    ## CALCULATE AMPLITUDE #####################################################    
    """ Calculate amplitude from pulse energy (if known: True), or calculate pulse 
        energy from the parameters.
    """    
    
 #   if Pulseenergy_known in 'True' :
 #       Ep=Energ
 #   else:
 #       """https://www.rp-photonics.com/solitons.html
 #       """
 #       n22=n2#*10**-4       # nonlinear refractive index [cm^2/TW -> um^2/W]
 #       Ep=((2*np.abs(mProp.beta2(p01_w0))/((n22*p01_w0)/(0.3*np.pi*20**2)*p01_t0)))*10**-15#
        
    if pulseshape in 'sech':
        """https://www.rp-photonics.com/peak_power.html
        """
        factor=0.88
    else:
        factor=0.94
    
    Amp1 =p01_N*inputPulseAmplitude_sFMAS(factor,Energ, p01_t0,p01_w0, aeff,mProp,n2,Pulseenergy_known) # amplitude in sqrt(TW/cm^2)
    Amp2 =p02_N*inputPulseAmplitude_sFMAS(factor,Energ, p02_t0,p02_w0,aeff, mProp,n2,Pulseenergy_known) # amplitude in sqrt(TW/cm^2)

    ## PROPAGATION TYPE #######################################################
    """ Reads in the given pulse configuration, calculates an analytical signal and
        runs the propagation routine
    """

    if p_n in 'one':
        E0Curr1 = initialElectricField(Amp1, p01_t0, p01_w0)
        E0 = E0Curr1(t + p01_tOff) 
        fdBase += '_%s_t0%4.3lf_w0%4.3lf_N%4.3lf'%('p01',p01_t0, p01_w0, p01_N)

    elif p_n in 'two':
        E0Curr1 = initialElectricField(Amp1, p01_t0, p01_w0)
        E0Curr2 = initialElectricField(Amp2, p02_t0, p02_w0)
        E0 += E0Curr1(t) + E0Curr2(t + p02_tOff)
        fdBase += '_%s_t0%4.3lf_w0%4.3lf_N%4.3lf_%s_t0%4.3lf_w0%4.3lf_N%4.3lf_tOff%4.3lf'%('p01',p01_t0, p01_w0, p01_N, 'p02', p02_t0, p02_w0, p02_N, p02_tOff)

   
    # CONVERT INITIAL PULSES TO PROPER INITIAL CONDITION
    (t, Eps0_t, w, Eps0_w )= initialCondition(t, E0)
   
    ## VELOCITY FOR RETARDED FRAME OF REFERENCE ###############################        
    vRef = estimateReferenceGroupVelocity(t, mProp, Amp1, p01_t0, p01_w0, p01_tOff)
    
    ## IONIZATION #############################################################
    gamma=n2/(c0)
    ioni = fdz(Eps0_t,gamma, t, w, ionization)
    
    ## SET UP MEASUREMENT DATA STRUCTURE ###################################### 
    zSkip = 100  
    monitor = Observable(t, dz,Nz,zSkip, mProp)

    ## PROPAGATE INITIAL CONDITION ############################################   
    sFMASDQRK4_2(t,w,Eps0_w,dz,dt,Nz, wMinFiber, wMaxFiber,mProp,ioni,ionization, vRef,monitor.measure,n2, pressure_he)  	
    fdBase += '_%s_%s'%('_ioni', ionization) 
   
        
    ## WRITE DATA TO OUTFILES ################################################# 
    monitor.saveNpz('./Data/', fdBase, Nt)


main()
# EOF: main_FMAS1D.py
