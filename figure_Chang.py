import sys
import os
import matplotlib as mpl
import numpy as np
import numpy.fft as nfft
import matplotlib.pyplot as plt
import matplotlib.colors as col
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.ticker import AutoMinorLocator
#import sys; sys.path.append('../src')
#from inputPulseTRY import *
#from mediumDispersion import *

def fetchNpz(iPath):
    data = np.load(iPath)
    return data['z'], data['t'], data['w'], data['Eps_w']

def getClosestIndex(z,z0):
    z0 = z.min() if z0<z.min() else z0
    z0 = z.max() if z0>z.max() else z0
    dz = z[1]-z[0]
    zIdx = int((z0-z.min())/dz)
    return zIdx

def subfigHelper(ax, x, z, u, cmap):
   # u = u/u.max()
    im = ax.imshow(
            u,
            origin = 'lower',
            interpolation = 'bilinear',
            aspect = 'auto',
            extent = [x.min(), x.max(), z.min(), z.max()],
	    vmin=-50,
	    vmax=0,
           # norm = col.LogNorm(vmin=u.max()*1e-5,vmax=u.max()),
            cmap=cmap
        )
    return im

def subfigHelper_lin(ax, x, z, u, cmap):
    u = u/u.max()
    im = ax.imshow(
            u,
            origin = 'lower',
            interpolation = 'bilinear',
            aspect = 'auto',
            extent = [x.min(), x.max(), z.min(), z.max()],
	    vmin=0,
	    vmax=1,
           # norm = col.LogNorm(vmin=u.max()*1e-5,vmax=u.max()),
            cmap=cmap
        )
    return im




def save_figure(fig_format = None, fig_name = 'test'):
    """ save figure

    Function that saves figure or shows interactive plot

    Note:
    - if no valid option is provided, an interactive plot is shown

    Args:
      fig_format (str): format to save figure in (options: png, pdf, svg)
      fig_name (str): name for figure (default: 'test')
    """
    if fig_format == 'png':
        plt.savefig(fig_name+'.png', format='png', dpi=600)
    elif fig_format == 'pdf':
        plt.savefig(fig_name+'.pdf', format='pdf', dpi=600)
    elif fig_format == 'svg':
        plt.savefig(fig_name+'.svg', format='svg')
    else:
        plt.show()



def set_style():
    """set figure style

    Function that customizes the default style to be conform with the Physical
    Review style and notation guide [1]. For instructions on how to set the
    default style using style sheets see [2].

    Notes:
    - main font size is chosen as 8pt, matching the fontsize of figure captions
    - fontsize of legends and auxiliary text labels are set to 6pt, exceeding
      the minimally required pointsize of 4.25 pt. (1.5 mm)
    - default rc (rc = "run commands", i.e. startup information) settings are
      changed dynamically
    - the custom font-scheme 'type2' depends on your latex installation and
      is not guaranteed to run on your specific system

    Refs:
      [1] https://journals.aps.org/prl/authors
      [2] https://matplotlib.org/3.3.1/tutorials/introductory/customizing.html
    """

    fig_width_1col = 3.4        # figure width in inch
    fig_width_2col = 7.0        # figure width in inch
    fig_aspect_ratio = 0.8       # width to height aspect ratio
    font_size = 8               # font size in pt
    font_size_small = 6         # font size in pt
    font_scheme = 'type2'          # options: 
                                #   None    - default matplotlib fonts
                                #   'type1' - text: Helvetica, math: Computer modern
                                #   'type2' - text: Helvetica, math: Helvetica 

    mpl.rcParams['figure.figsize'] = fig_width_2col, fig_aspect_ratio*fig_width_2col
    mpl.rcParams['axes.labelsize'] = font_size
    mpl.rcParams['font.size'] = font_size
    mpl.rcParams['legend.fontsize'] = font_size_small
    mpl.rcParams['xtick.labelsize'] = font_size
    mpl.rcParams['ytick.labelsize'] = font_size
    mpl.rcParams['xtick.direction'] = 'out'
    mpl.rcParams['ytick.direction'] = 'out'
    mpl.rcParams['lines.linewidth'] = 1.1
    mpl.rcParams['axes.linewidth'] =  0.5

    if font_scheme == 'type1':
        mpl.rcParams['text.usetex'] = True
        mpl.rcParams['font.family'] = 'sans-serif'
        mpl.rcParams['font.sans-serif'] = 'Helvetica'
        mpl.rcParams['mathtext.fontset'] = 'cm'

    if font_scheme == 'type2':
        mpl.rcParams['text.usetex'] = True
        mpl.rcParams['text.latex.preamble'] = r'\usepackage{siunitx}'  r'\sisetup{detect-all}'          r'\usepackage{helvet}'   r'\usepackage{sansmath}'  r'\sansmath'



def set_colorbar(fig, img, ax, ax2):
    """set colorbar

    Function that generates a custom colorbar. For more options on
    colorbars and tick parameters see Refs. [1,2]

    Refs:
      [1] https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.colorbar.html
      [2] https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.tick_params.html

    Args:
      fig (object): main figure object
      img (object): image that should be described
      ax (object): axes object with position info for colorbar placement
    """
    # -- extract position information for colorbar placement
    refPos = ax.get_position()
    x0, y0, w, h = refPos.x0, refPos.y0, refPos.width, refPos.height
    refPos2 = ax2.get_position()
    w2, h2 = refPos2.width, refPos2.height
    # -- set new axes as reference for colorbar
    colorbar_axis = fig.add_axes([x0+1.05*w, y0, 0.07*w, h2])

    # -- set custom colorbar
    colorbar = fig.colorbar(img,        # image described by colorbar
            cax = colorbar_axis,        # reference axex
            orientation = 'vertical' # colorbar orientation
            #extend = 'both'             # ends with out-of range values
            )
    colorbar.ax.tick_params(
            color = 'k',                # tick color 
            labelcolor = 'k',           # label color
            bottom = False,             # no ticks at bottom
            labelbottom = False,        # no labels at bottom
            labeltop = True,            # labels on top
            top = True,                 # ticks on top
            direction = 'out',          # place ticks outside
            length = 2,                 # tick length in pts. 
            labelsize = 6.,             # tick font in pts.
            pad = 1.5                    # tick-to-label distance in pts.
            )
    return colorbar

def set_colorbar_right(fig, img, ax, ax2):
    """set colorbar

    Function that generates a custom colorbar. For more options on
    colorbars and tick parameters see Refs. [1,2]

    Refs:
      [1] https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.colorbar.html
      [2] https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.tick_params.html

    Args:
      fig (object): main figure object
      img (object): image that should be described
      ax (object): axes object with position info for colorbar placement
    """
    # -- extract position information for colorbar placement
    refPos = ax.get_position()
    x0, y0, w, h = refPos.x0, refPos.y0, refPos.width, refPos.height
    refPos2 = ax2.get_position()
    w2, h2 = refPos2.width, refPos2.height
    # -- set new axes as reference for colorbar
    colorbar_axis = fig.add_axes([x0, y0+1.03*h, w, 0.05*h2])

    # -- set custom colorbar
    colorbar = fig.colorbar(img,        # image described by colorbar
            cax = colorbar_axis,        # reference axex
            orientation = 'horizontal' # colorbar orientation
            #extend = 'both'             # ends with out-of range values
            )
    colorbar.ax.tick_params(
            color = 'k',                # tick color 
            labelcolor = 'k',           # label color
            bottom = False,             # no ticks at bottom
            labelbottom = False,        # no labels at bottom
            labeltop = True,            # labels on top
            top = True,                 # ticks on top
            direction = 'out',          # place ticks outside
            length = 2,                 # tick length in pts. 
            labelsize = 6.,             # tick font in pts.
            pad = 1.5                    # tick-to-label distance in pts.
            )
    return colorbar
    


def set_legend(ax, lines):
    """set legend

    Function generating a custom legend, see [1] for more options

    Refs:
      [1] https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.legend.html

    Args:
      ax (object): figure part for which the legend is intended
      lines (list): list of  Line2D objects
    """
    # -- extract labels from lines
    labels = [x.get_label() for x in lines]
    # -- customize legend 
    ax.legend(lines,                # list of Line2D objects
              labels,               # labels 
              title = '',           # title shown on top of legend 
              loc = (0.03,0.1),    # location of the legend
	      ncol = 1,             # number of columns
              labelspacing = 0.01,  # vertical space between handles in font-size units
              borderpad = 0.3,      # distance to legend border in font-size units
              handletextpad = 0.6,  # distance between handle and label in font-size units
              handlelength = 1.0,   # length of handle in font-size units
              frameon = True,     # remove background patch
	      edgecolor='none',
	      facecolor='white',
	      framealpha=0.7,
              fontsize=5,
	      )

def getXfrogTrace(t,w,Et, s0=None, tauMin=-np.inf, tauMax=np.inf, Ntau=2000, Nw=1500):
      """Crosscorrelation XFROG trace for time-domain analytic signal

      Computes crosscorrelation XFROG trace for the amplitude and phase
      analysis of an input test signal (in time-domain) applying the XFROG
      method discussed in Ref. [1].  Technically, the crosscorrelation XFROG
      trace is computed via a short time Fourier transform [2] of the test
      signal using a Gaussian window function.

      Args:
          t (numpy-array, ndim=1): t-axis
          w (numpy-array, ndim=1): angular frequency axis
          Et (numpy-array, ndim=1): analytic signal in time-domain
          s0 (float): width of Gaussian gate function (default: 1/sqrt(2*pi),
              i.e. special parameterless STFT referred to as Gabor transform)
          tauMin (float): lower bound for delay time (default: -inf)
          tauMax (float): upper bound for delay time (default: inf)
          Ntau (int): number of delay times in [tMin, tMax) (default: 500)
          Nw (int): number of stored samples in w-domain (default: 500)

      Returns:
          delayTimes (numpy-array, ndim=1): delay times used in XFROG analysis
          w (numpy-array, ndim=1): angular frequencies kept for output
          IXFROG (numpy-array, ndim=2): XFROG trace in (tau, omega) plane

      Note:
          - XFROG method implemented via short time Fourier transform (STFT)
            using a Gaussian window function
          - STFT resolution issues [2]: narrow window-function yields
            good time resolution but poor frequency resolution; wide window
            function yields good frequency resolution but bad time resolution.
          - Parameterless Gabor transform obtained by setting s0=1/sqrt(2*pi)

      Refs:
          [1] XFROG - A New Method for Amplitude and Phase Characterization of
              Weak Ultrashort Pulses
              Linden, S. and Giessen, H. and Kuhl, J.
              Phys. Stat. Sol., 206 (1998) 119

          [2] https://en.wikipedia.org/wiki/Short-time_Fourier_transform
      """
      # INITIALIZE ANALYSIS PARAMETERS ########################################
      tauMin = np.max((tauMin, np.min(t)))
      tauMax = np.min((tauMax, np.max(t)))
      Ntau = np.min((Ntau, t.size))
      Nw = 2**(int(np.min((Nw, w.size/2)))-1).bit_length()/2
      
      
      s0 = s0 if s0 else 1./np.sqrt(2*np.pi) # 10*(t[1]-t[0])

      # SET DELAY TIMES IN DESIRED RANGE ######################################
      delayTimes = np.linspace(tauMin, tauMax, Ntau)

      # SET ONE-SIDED INDEX MASK FOR OUTPUT ###################################
      wIdx = slice(0,int(w.size/2), int(w.size/2/Nw))
      

      # SET EMPTY SPECTROGRAM #################################################
      IXFROG = np.zeros((int(Ntau),int(Nw)))

      # PHASE UNWRAPPING FOR TESTPULSE ########################################
      absEt = np.abs(Et)
      phaseEt = np.unwrap(np.angle(Et))

      # SET PROPERLY NORMALIZED GAUSSIAN WINDOW FUNCTION ######################
      gateFunc = lambda t: np.exp(-t**2/2/s0/s0)/np.sqrt(2.*np.pi*s0*s0)

      # YIELD CROSSCORRELATION XFROG TRACE, SEE EQ. (2) OF REF. [1] ########### 
      for tauId, tau in enumerate(delayTimes):
          Ecc = absEt*gateFunc(t-tau)*np.exp(1j*phaseEt)
          I = np.abs(nfft.ifft(Ecc))**2
          IXFROG[tauId,:] = I[wIdx]

      return delayTimes, w[wIdx],  np.add(IXFROG,10**-10)

def custom_colormap():
    # -- CREATE CUSTOM COLORMAP
    from matplotlib.colors import ListedColormap
    cmap_base = mpl.cm.jet(np.arange(256))
    blank = np.ones((25,4))
    for i in range(3):
        blank[:,i] = np.linspace(1,cmap_base[0,i], blank.shape[0])
    my_cmap = ListedColormap(np.vstack((blank, cmap_base )))
    return my_cmap
    


def generate_figure( fig_format=None, fig_name='fig01'):
    """generate figure

    Function generating a figure reproducing FIG. 3 of [1].

    Refs:
      [1] Experimental Observation of Picosecond Pulse Narrowing and Solitons in Optical Fibers
          L. F. Mollenauer, R. H. Stolen, and J. P. Gordon
          Phys. Rev. Lett. 45, 1095 (1980)

    Args:
      data_set_01 (tuple): data set of the form (z, t, Azt, N), where
                              z (1D array): z samples
                              t (1D array): t samples
                              Azt (2D array): complex-valued field A(z,t)
                              N (float): soliton order
      data_set_02 (tuple): data set of the form (z, t, Azt, N)
      fig_format (str): format for output figure
                        (choices: png, pdf, svg; default: interactive figure)
      fig_name (str): name for output figure wihtout suffix (default='fig01')
    """

    # SET A STYLE THAT FITS THE TARGET JOURNAL
    set_style()

   
    # SET FIGURE LAYOUT
    fig = plt.figure()
    plt.subplots_adjust(left = 0.08, bottom = 0.07,
                        right =0.98, top = 0.91,
                        wspace = 0.2, hspace = 0.3)

    gs = GridSpec(12,11)

    
   
    axa = fig.add_subplot(gs[0:5,0:3])
    axb = fig.add_subplot(gs[0:5,3:6], sharey=axa)  

    axd= fig.add_subplot(gs[7:,3:6], sharex=axb,sharey=axa)
    axc = fig.add_subplot(gs[7:,0:3],sharex=axa, sharey=axa)

    ax1 = fig.add_subplot(gs[0:3,8:])
    ax2 = fig.add_subplot(gs[3:6,8:],sharex=ax1,sharey=ax1)
    ax3 = fig.add_subplot(gs[6:9,8:],sharex=ax1,sharey=ax1)
    ax4 = fig.add_subplot(gs[9:12,8:],sharex=ax1,sharey=ax1)
    



 # SUBPLOT 1 - SET AXIS DETAILS ---------------------------------------------------

 
 


     # SUBPLOT 1 - SET AXIS DETAILS ---------------------------------------------------

    axa.text(0.02,0.92,r"(a)",color='k', weight='bold', fontsize=8,
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=axa.transAxes)

    # -- customize x-axis
    x_lima = (-0.03,0.03)
    x_ticksa = (-0.02,0,0.02)
    axa.set_xlim(x_lima)
    axa.set_xticks(x_ticksa)
    axa.tick_params(axis='x', length=3, pad=2)   
    axa.set_xlabel(r"Time $t$ (ps)")

    # -- customize y-axis
    y_lima = (0,0.1)
    y_ticksa = (0,0.02,0.04,0.06,0.08,0.1)
    axa.set_ylim(y_lima)
    axa.set_yticks(y_ticksa) 
    axa.tick_params(axis='y', length=3., pad=2, top=False)
    axa.set_ylabel(r"Propagation distance $z$ (m)")  
    

    # SUBPLOT 2 - SET AXIS DETAILS --------------------------------------------------

    axb.text(0.02,0.92,r"(b)",color='k', weight='bold', fontsize=8,
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=axb.transAxes)

    # -- customize x-axis
    x_limb = (0,12)
    x_ticksb = (3,6,9,12)
    axb.tick_params(axis='x', length=3., pad=2, right=False)
    axb.set_xlim(x_limb)
    axb.set_xticks(x_ticksb)
    axb.set_xlabel(r"Angular frequency $\omega$ (rad/fs)")

    axbt=axb.twiny()
    new_tick_locations = np.array([3,6,9,12])

    def tick_function(X):
        lam = 2*np.pi*0.3/X
        return ["%.3f" % z for z in lam]

    axbt.set_xlim(axb.get_xlim())
    axbt.set_xticks(new_tick_locations)
    axbt.set_xticklabels(tick_function(new_tick_locations))
    axbt.set_xlabel(r"Wavelength $\lambda$ ($\mu$m)")

    # -- customize y-axis
    axb.tick_params(axis='y', length=3., pad=2, top=False, labelleft=False)
  
 
    # SUBPLOT 3 - SET AXIS DETAILS ---------------------------------------------------

    axc.text(0.02,0.92,r"(c)",color='k', weight='bold', fontsize=8,
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=axc.transAxes)

    # -- customize x-axis
    axc.tick_params(axis='x', length=3, pad=2)
    axc.set_xlabel(r"Time $t$ (ps)")

    # -- customize y-axis
   
    axc.tick_params(axis='y', length=3., pad=2, top=False)
    axc.set_ylabel(r"Propagation distance $z$ (m)") 

     # SUBPLOT 4 - SET AXIS DETAILS ---------------------------------------------------

    axd.text(0.02,0.92,r"(d)",color='k', weight='bold', fontsize=8,
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=axd.transAxes)

    # -- customize x-axis
    axd.tick_params(axis='x', length=3, pad=2)

    # -- customize y-axis    
    axd.tick_params(axis='y', length=3., pad=2, top=True, labelleft=False)
    axd.set_xlabel(r"Angular frequency $\omega$ (rad/fs)")

    axdt=axd.twiny()
    new_tick_locations = np.array([3,6,9,12])

    def tick_function(X):
        lam = 2*np.pi*0.3/X
        return ["%.3f" % z for z in lam]

    axdt.set_xlim(axb.get_xlim())
    axdt.set_xticks(new_tick_locations)
    axdt.set_xticklabels(tick_function(new_tick_locations))
    axdt.set_xlabel(r"Wavelength $\lambda$ ($\mu$m)")

    # SUBPLOT 5 - SET AXIS DETAILS ---------------------------------------------------

    ax1.text(0.02,0.85,r"(e)",color='k', weight='bold', fontsize=8,
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax1.transAxes)

    # -- customize x-axis
    x_lim1 = (-0.03,0.03)
    x_ticks1 = (-0.02,0,0.02)
    ax1.tick_params(axis='x', length=3., pad=2, right=False, labelbottom=False)
    ax1.set_xlim(x_lim1)
    ax1.set_xticks(x_ticks1)

    # -- customize y-axis    
    y_lim1 = (0,220)
    y_ticks1 = (0,50,100,150,200)
    ax1.set_ylim(y_lim1)
    ax1.set_yticks(y_ticks1)
    ax1.tick_params(axis='y', length=3., pad=2, top=True)
    ax1.set_ylabel(r"$|\mathcal{E}|^2$")


    # SUBPLOT 6 - SET AXIS DETAILS ---------------------------------------------------

    ax2.text(0.02,0.85,r"(f)",color='k', weight='bold', fontsize=8,
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax2.transAxes)

    # -- customize x-axis
    ax2.tick_params(axis='x', length=3., pad=2, right=False, labelbottom=False)
  

    # -- customize y-axis    
    ax2.tick_params(axis='y', length=3., pad=2, top=True)
    ax2.set_ylabel(r"$|\mathcal{E}|^2$")

    # SUBPLOT 6 - SET AXIS DETAILS ---------------------------------------------------

    ax3.text(0.02,0.85,r"(g)",color='k', weight='bold', fontsize=8,
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax3.transAxes)

    # -- customize x-axis
    ax3.tick_params(axis='x', length=3., pad=2, right=False, labelbottom=False)
  

    # -- customize y-axis    
    ax3.tick_params(axis='y', length=3., pad=2, top=True)
    ax3.set_ylabel(r"$|\mathcal{E}|^2$")

    # SUBPLOT 6 - SET AXIS DETAILS ---------------------------------------------------

    ax4.text(0.02,0.85,r"(h)",color='k', weight='bold', fontsize=8,
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax4.transAxes)

    # -- customize x-axis
    ax4.tick_params(axis='x', length=3., pad=2, right=False)
    ax4.set_xlabel(r"Time $t$ (ps)")
  

    # -- customize y-axis    
    ax4.tick_params(axis='y', length=3., pad=2, top=True)
    ax4.set_ylabel(r"$|\mathcal{E}|^2 (\mathrm{TW}/\mathrm{cm^2})$")
   

  
 
    def phasematch_sub(we, w0,width,w02,width2):
        """ PM of individual subpulse (w0) corrected for secondary subpulse (w02)
        """
        gamma= ((3*w0)/(8*0.3*mProp.n(w0)))
        gamma2= ((3*w02)/(8*0.3*mProp.n(w02)))
        a1=np.abs(mProp.beta2(w0))/(width**2*gamma)
        a2=np.abs(mProp.beta2(w02))/(width2**2*gamma2)	
        nw=mProp.beta(w0)+mProp.beta1(w0)*(we-w0)
        k1=nw+0.5*gamma*(a1**2+2*a2**2)
        k2=nw+0.5*gamma2*(a2**2+2*a1**2)	
        return  k1 -mProp.beta(we)#=0 -> chrenkov
    

    def phasematch_soliton(w, w0,width):
        """ PM of individual soliton (w0) 
        """
        gamma= ((3*w0)/(8*0.3*mProp.n(w0)))	
        a1=np.abs(mProp.beta2(w0))/(width**2*gamma)	
        nw=mProp.beta(w0)+mProp.beta1(w0)*(w-w0)
        kn=0.5*gamma*a1**2+nw
        return mProp.beta(w)-kn #=0 -> chrenkov

    def phasematch_sub_FWM(w, w0,width,w02,width2):
        """ PM of individual subpulse (w0) corrected for secondary subpulse (w02)
        """
        gamma= ((3*w0)/(8*0.3*mProp.n(w0)))
        gamma2= ((3*w02)/(8*0.3*mProp.n(w02)))
        a1=np.abs(mProp.beta2(w0))/(width**2*gamma)
        a2=np.abs(mProp.beta2(w02))/(width2**2*gamma2)	
        nw=mProp.beta(w0)+mProp.beta1(w0)*(w-w0)
        k1=0.5*gamma*(a1**2+2*a2**2)+nw
        k2=0.5*gamma2*(a2**2+2*a1**2)+nw	
        return mProp.beta(w)+k1-k2-mProp.beta(2.886) #=0 -> chrenkov
    



    ## SET RANGE PARAMETERS #################################################
    fileLoc = str(sys.argv[1]) 
    
    # File without ionization effect
    fileLoc2 = 'obs_OHNE-sFMAS-RK4-da_Nt65536_Nz20000_p01_t030.000_w02.356_N1.000__ioni_on.npz'#str(sys.argv[1]) 
  
    ## FETCH DATA FROM FILE #################################################
    (z, t, w,  Eps_w) = fetchNpz(fileLoc)
    (z2, t2, w2,  Eps_w2) = fetchNpz(fileLoc2)
   
    z = z*1e-6 # change scale to meters (m)
    t = t*1e-3 # change scale to picoseconds (ps)
    z2 = z2*1e-6 # change scale to meters (m)
    t2 = t2*1e-3 # change scale to picoseconds (ps)
    

    
    wMinFiber = 0.5                   
    wMaxFiber = 50.0
    c0 = 0.29979
  #  mProp = MediumDispersion( "WD2017", wMinFiber, wMaxFiber)
  #  we=np.linspace(0.5,5,10000) 
  #  b1 =mProp.beta1(we)
  #  b2= mProp.beta2(we)

    ZDF1=1.511
    ZDF2=2.51
  


    # CONVENIENT FUNCTIONS ##################################################################
    _norm = lambda x: x/x[0].max()
    def t_pos(Eps):
        Et=nfft.fft(Eps)      
        It = np.abs(Et)**2        
        return np.where(It>=0.99*max(It))[0]
    factor = t_pos(Eps_w[0])[0]*0.0152# Multiplication factor x maps t to tau (for XFROG) with x =Ntau/Nt (e.g. 200/131072)
    _normXFROG = lambda x: x/200*max(x[int(factor)])
    _truncate = lambda x: np.where(x>1.e-20,x,1.e-20)
    _dB = lambda x: np.where(x>1e-20,10.*np.log10(x),10*np.log10(1e-20))
    cmap=custom_colormap()#mpl.cm.get_cmap('jet') 

    # CONVENIENT FIELD DESCRIPTIONS ############################################

    ## Mask in frequency domain #########################################################    
  
    A1Mask = np.logical_and(w < ZDF1, w > wMinFiber)
    A2Mask = np.logical_and(w < wMaxFiber, w > ZDF2)
    NormMask = np.logical_and(w < ZDF2, w > ZDF1)    

   ## TRAJECTORIES ###################################################################
    def Max_intens(E):
        I=[]
        for i in range(0,np.size(z)):
            Imax=max(np.abs(nfft.fft(E[i]))**2)
            x=np.where(np.abs(nfft.fft(E[i]))**2>=0.99*Imax)[0][0]
            I.append(t[x])
        return I
    
    def t_centroid(s,Eps):
        E=np.copy(Eps)	
        I=np.abs(nfft.fft(E[s]))**2
        Imax=max(I)
        x=np.where(I>=0.99*Imax)[0]
        C=sum(t[x]*I[x])/sum(I[x])	
        return C
  
    def w_centroid(s,Eps,wmask):
        E=np.copy(Eps)
        E[:,~wmask]=0	
        Iw=np.abs(E[s])**2
        Imax=max(Iw)
        x=np.where(Iw>=0.02*Imax)[0]
        C=sum(w[x]*Iw[x])/sum(Iw[x])	
        return C

    def energy_part(Eps_part,wmask):
        Eps_ges = Eps_w[0]  
       # we = np.linspace(wMinFiber, wMaxFiber,65536)        
        mu = 1.0
        betaDS = mProp.beta(wR4) #- w/vBoost='--')
        A = betaDS/(2*mu*wR4)
        return np.trapz(np.abs(A[wmask])*np.abs(Eps_part[wmask])**2, wR4[wmask])

   
   
  
   


    # SUBFIGS 1 ###################################################  
    im1 = subfigHelper_lin(axa,t, z, np.abs(nfft.fft(Eps_w))**2, cmap)
    

    # SUBFIGS 2 ################################################### 

    im2 = subfigHelper(axb,nfft.fftshift(w[::1]), z, nfft.fftshift(_dB(_truncate(_norm(np.abs(Eps_w[::1])**2))), axes=-1), cmap)

    axb.axvline(x=4.47, color='k', linestyle='--') #lambda =421 nm  ZDW

    # SUBFIGS 3 ###################################################  
    im3 = subfigHelper_lin(axc,t, z, np.abs(nfft.fft(Eps_w2))**2, cmap)
    

    # SUBFIGS 4 ################################################### 

    im4 = subfigHelper(axd,nfft.fftshift(w[::1]), z, nfft.fftshift(_dB(_truncate(_norm(np.abs(Eps_w2[::1])**2))), axes=-1), cmap)

    axd.axvline(x=4.47, color='k', linestyle='--') #lambda =421 nm  ZDW

    # SUBFIGS 5 ################################################### 

     
    ax1.plot(t[::-1],(nfft.ifft(Eps_w[62]))**2, color='blue', label='$\mathcal{E}^2$')   
    ax1.plot(t[::-1],np.abs(nfft.ifft(Eps_w[62]))**2, color='red') 
    axa.axhline(y=0.0311, color='k', linestyle='--')   


    # SUBFIGS 6 ################################################### 

    
    ax2.plot(t[::-1],(nfft.ifft(Eps_w[82]))**2, color='blue', label='$\mathcal{E}^2$')  
    ax2.plot(t[::-1],np.abs(nfft.ifft(Eps_w[82]))**2, color='red')   
    axa.axhline(y=0.0413, color='k', linestyle='--')   


    # SUBFIGS 7 ################################################### 

    
    ax3.plot(t[::-1],(nfft.ifft(Eps_w[89]))**2, color='blue', label='$\mathcal{E}^2$')    
    ax3.plot(t[::-1],np.abs(nfft.ifft(Eps_w[89]))**2, color='red') 
    axa.axhline(y=0.0447, color='k', linestyle='--')   


    # SUBFIGS 8 ################################################### 

    
    l2=ax4.plot(t[::-1],(nfft.ifft(Eps_w[96]))**2, color='blue', label='$\mathcal{E}^2$')   
    l1=ax4.plot(t[::-1],np.abs(nfft.ifft(Eps_w[96]))**2, color='red', label='$|\mathcal{E}|^2$') 
    axa.axhline(y=0.0483, color='k', linestyle='--')   


 #   set_legend(ax1, l1+l2)




    # COLORBARS AND LABELS #############################################
    
    cb = set_colorbar(fig, im2, axb,axb)
    cb.set_ticks((-40,-20,0))
    cb.ax.set_title('$|\mathcal{E}_{\omega}|^2\,\mathrm{(dB)}$', fontsize=6., y=0.5,x=3.5, rotation=90)

    cb2	 = set_colorbar_right(fig, im1, axa,axa)
    cb2.set_ticks((0,0.2,0.4,0.6,0.8,1))
    cb2.ax.set_title('$|\mathcal{E}|^2\,(\mathrm{norm.})$', fontsize=6., y=2)


   
    # SAVE FIGURE ########################################################
    # - Set output file name    
    path, inFileName = os.path.split(fileLoc)
    inFileBasename = os.path.splitext(inFileName)[0]
    oName = '../FIGS/'+inFileBasename+'_T'
    # - save
    save_figure(fig_format, oName)


def main():
    
    generate_figure( fig_format='pdf', fig_name='bla')
    plt.show()

if __name__ == '__main__':
    main()

