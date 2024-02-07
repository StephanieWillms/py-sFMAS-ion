import sys
import os
import numpy as np
import numpy.fft as nfft


class Observable(object):
    """data-structure facilitating measurement
    """
    def __init__(self, t, dz,Nz,zSkip, mProp):
        """generates instance of Observable data structure

        Note: 
            - fiber properties, available through argument mProp, are 
              required in order to compute photon number and energy later on

        Args:
            t  (numpy-array, ndim=1): time axis
            dz (float): z-axis, i.e. propagation direction, increment
            Nz (int): number of propagation steps
            zSkip (int): store only configurations in steps modulo zSkip
            mProp (object): data structure holding fiber properties

        Attrib:
            w (numpy-array, ndim=1): anglular frequency axis
            nw (numpy-array, ndim=1): refractive index profile
            c0 (float): free-space wave speed
            z (numpy-array, ndim=1): z-axis, i.e. propagation direction axis
            t  (numpy-array, ndim=1): time axis
            Eps_w (numpy-array, ndim=2): frequency components of field 

        """
        self.w = nfft.fftfreq(t.size,d=t[1]-t[0])*2*np.pi
        self.nw = mProp.n(self.w)
        self.c0 = mProp.c0
        self.dnz = zSkip
        self.z = np.arange(0,Nz,self.dnz)*dz
        self.t = t
        self.Eps_w = np.zeros((self.z.size,t.size),dtype=complex)

    def measure(self, n, z, t, w, Aw): 
        """measure

        Callback function facilitating measurement

        Args:
            n (int): current propagation step
            z (flozt): current z position
            t (numpy-array, ndim=1): time-axis
            w (numpy-array, ndim=1): angular frequency axis
            Aw (numpy-array, ndim=1): frequency components of field
        """
        if n%self.dnz==0:
            myId = n/self.dnz 
            self.Eps_w[int(myId),:] = Aw[:]


    def saveNpz(self, oPath, fdBase, auxStr = ''):
        """save data in numpy format 

        Saves data in numpy native compressed npz format

        Args:
            oPath (str): path to folder where output will be stored
            fdBase (str): base name for output files
            auxStr (str): auxiliary information to store with simulated data
        """
        try:
            os.makedirs(oPath)
        except OSError:
            pass

        np.savez_compressed(
            oPath + 'obs' + fdBase + '.npz', z = self.z, t = self.t, w = self.w, 
            Eps_w = self.Eps_w, nw = self.nw, c0 = self.c0, info=auxStr)

    def saveGpi(self, oPath, fdBase, auxStr = '', xSamp=100, ySamp=100):
        """save data in gnuplot format

        Saves data in gnuplot readable nonuniform matrix format

        Note:
            - in gnuplot, data can be viewed with command:
              p 'dataset' nonuniform matrix with image
        
        Args:
            oPath (str): path to folder where output will be stored
            fdBase (str): base name for output files
            auxStr (str): auxiliary information to store with simulated data
            xSamp (int): sample points in x direction
            ySamp (int): sample points in y direction
        """
        # BEGIN HELPER FUNCTION ###############################################
        def _dumpGpi(fName, y, x, xSamp, ySamp, func, auxStr=''):
             """helper function converting data to gnuplot format

             Helper function storing possibly reduced dataset in gnuplot 
             readable nonuniform matrix format

             Args:
                 fName (str): full path to output filename 
                 y (numpy-array, ndim=1): full y-axis
                 x (numpy-array, ndim=1): full x-axis
                 xSamp (int): sample points in x direction
                 ySamp (int): sample points in y direction
                 func (object): function to apply to data row-wise
             """
             Nx = int(float(x.size)/min(x.size,xSamp))
             Ny = int(float(y.size)/min(y.size,ySamp))
             xIdx = range(0,x.size,Nx) 
             yIdx = range(0,y.size,Ny) 
             x = x[xIdx]
             y = y[yIdx]

             fStream = open(fName,'w')  

             fStream.write("%s\n"%(auxStr))

             fStream.write("%d "%(x.size))
             for i in range(x.size):
                 fStream.write("%g "%(x[i]))
             fStream.write("\n")

             for j in range(y.size):
                 fStream.write("%g "%(y[j]))
                 E = func(self.Eps_w[yIdx[j],:])[xIdx]
                 for i in range(x.size):
                     fStream.write("%g "%(np.abs(E[i])**2))
                 fStream.write("\n")
        # END HELPER FUNCTION #################################################

        try:
            os.makedirs(oPath)
        except OSError:
            pass

        _wDomFunc = nfft.fftshift
        _tDomFunc = nfft.ifft

        _dumpGpi(oPath+'wDomain'+fdBase+'.dat', self.z, nfft.fftshift(self.w), (xSamp,ySamp),_wDomFunc, auxStr)
        _dumpGpi(oPath+'tDomain'+fdBase+'.dat', self.z, self.t, (xSamp,ySamp),_tDomFunc, auxStr)



