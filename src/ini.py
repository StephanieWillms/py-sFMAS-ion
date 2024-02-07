import numpy as np
from mediumDispersion import MediumDispersion 



## COMPUTATIONAL DOMAIN #################################################  
'''Here we set up the computational domain
'''
pulseshape='sech'                    # sech or gauss possible
Pulseenergy_known='True'

ionization='on' 		      # Two choices: 'on' or 'off'

Energ=2.5*10**-6                     # pulseenergy in Joule
tMax = 2000.0                        # bound for time-axis: (-tMax, tMax) 
Nt =  2**16                          # number of mesh points for time-axis
zMax = 10000                        # maximal propagation distance
Nz =   2000                          # number of mesh points for z-axis
dz =   zMax/Nz  
dt = 2*tMax/Nt
c0 =   0.29979                       # speed of light    
# ------- Fiber Props. ---------------      
wMinFiber = 0.2                      # fiber boundaries
wMaxFiber = 60.0                     # fiber boundaries
R = 8                               # Fiber core radius in um
st = 0.8                              # Fiber bridge diameter in um
pre= 3
aeff=np.pi*(R)**2                    # approximated effective fiber area 
# ------- Gas Props. ---------------
pressure_he = 4.0                    # gas pressure in atm 
#insetad of 0.4 we use 10.4 here, because was found in paper
n2 = pressure_he*10.4*1.0e-20*1e+12   # scaled fiber nonlinearity 0.4*1.0e-20*1e+10 # cm^2/W -> um^2/W
n2 = n2#*0.566                        # effective area factor
mProp = MediumDispersion( "kagome", (wMinFiber,wMaxFiber), c0, pre=4) 
 

 
  
