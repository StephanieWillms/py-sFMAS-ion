import sys
import os
import matplotlib as mpl
import numpy as np
import numpy.fft as nfft
import matplotlib.pyplot as plt
import matplotlib.colors as col
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.colors import ListedColormap, BoundaryNorm


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
  #  u = u/u.max()
    im = ax.imshow(
            u,
            origin = 'lower',
            interpolation = 'bilinear',
            aspect = 'auto',
            extent = [x.min(), x.max(), z.min(), z.max()],
	    vmin=-35,
	    vmax=0,
          #  norm = col.LogNorm(vmin=u.max()*1e-5,vmax=u.max()),
            cmap=cmap
        )
    return im

def subfigHelper_lin(ax, x, z, u, cmap):
   # u = u/u.max()
    im = ax.imshow(
            u,
            origin = 'lower',
            interpolation = 'bilinear',
            aspect = 'auto',
            extent = [x.min(), x.max(), z.min(), z.max()],
           # norm = col.LogNorm(vmin=u.max()*1e-5,vmax=u.max()),
            cmap=cmap
        )
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.tick_params(axis='both',length=2.)
    return im

def getXfrogTrace(t,w,Et, s0=None, tauMin=-np.inf, tauMax=np.inf, Ntau=2000, Nw=2500):
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
    fig_aspect_ratio = 0.65      # width to height aspect ratio
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
    mpl.rcParams['lines.linewidth'] = 0.8
    mpl.rcParams['axes.linewidth'] =  0.5

    if font_scheme == 'type1':
        mpl.rcParams['text.usetex'] = True
        mpl.rcParams['font.family'] = 'sans-serif'
        mpl.rcParams['font.sans-serif'] = 'Helvetica'
        mpl.rcParams['mathtext.fontset'] = 'cm'

    if font_scheme == 'type2':
        mpl.rcParams['text.usetex'] = True
        mpl.rcParams['text.latex.preamble'] =  r'\usepackage{siunitx}'  r'\sisetup{detect-all}'          r'\usepackage{helvet}'   r'\usepackage{sansmath}'  r'\sansmath'
        
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
    colorbar_axis = fig.add_axes([x0, y0 + 1.*h +0.05*h2, w, 0.1*h2])

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
            pad = 1.                    # tick-to-label distance in pts.
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
              loc = (0.73,0.75),    # location of the legend
	      ncol = 1,             # number of columns
              labelspacing = 0.01,  # vertical space between handles in font-size units
              borderpad = 0.3,      # distance to legend border in font-size units
              handletextpad = 0.6,  # distance between handle and label in font-size units
              handlelength = 1.0,   # length of handle in font-size units
              frameon = True,     # remove background patch
	      edgecolor='none',
	      facecolor='none',
	      framealpha=0.7,
              fontsize=6 ,
	      )
    


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
    plt.subplots_adjust(left = 0.07, bottom = 0.08,
                        right =0.92, top = 1.02,
                        wspace = 0.3, hspace = 0.2)

    gs = GridSpec(7,21)
    ax2 =fig.add_subplot(gs[3:,0:6])
    ax2w =fig.add_subplot(gs[3:,6:14],sharey=ax2)
    ax3w =fig.add_subplot(gs[4:,14:21],sharex=ax2w)
    ax2tw= fig.add_subplot(gs[1:3,6:14],sharex=ax2w)
    ax22tw= fig.add_subplot(gs[1:3,0:6],sharex=ax2,sharey=ax2tw)
    ax3tw= fig.add_subplot(gs[1:4,14:21],sharex=ax2w,sharey=ax3w)
    

    # SUBPLOT 1 - SET AXIS DETAILS --------------------------------------------------



    ax2.text(0.85,0.02,r"$|\mathcal{E}|^2$",color='w', weight='bold', fontsize=8,
        horizontalalignment='left',
        verticalalignment='bottom',
	#bbox=dict(facecolor='white', alpha=1,edgecolor='none',boxstyle='square,pad=0.0'),
        transform=ax2.transAxes)

    ax2w.text(0.5,0.02,r"$|\mathcal{E}_{\omega}|^2$",color='w', weight='bold', fontsize=8,
        horizontalalignment='left',
        verticalalignment='bottom',
	#bbox=dict(facecolor='white', alpha=1,edgecolor='none',boxstyle='square,pad=0.0'),
        transform=ax2w.transAxes)


    # -- customize x-axis
    x_lim1 = (-0.6,0.6)
    x_ticks1 = (-0.4,0,0.4)
    ax2.tick_params(axis='x', length=3., pad=2, right=False)
    ax2.set_xlim(x_lim1)
    ax2.set_xticks(x_ticks1)
    ax2.set_xlabel(r"Time $t$ (ps)")

    # -- customize x-axis
    x_lim1w = (0.8,3.5)
    x_ticks1w = (1,2,3)
    ax2w.tick_params(axis='x', length=3., pad=2, right=False)
    ax2w.set_xlim(x_lim1w)
    ax2w.set_xticks(x_ticks1w)
    ax2w.set_xlabel(r"Angular frequency $\omega$ (rad/fs)")

    # -- customize y-axis
    y_lim1 = (0.0,0.1)
    y_ticks1 = (0,0.02,0.04,0.06,0.08,0.1)
    ax2.tick_params(axis='y', length=3., pad=2, top=False)
    ax2.set_ylim(y_lim1)
    ax2.set_yticks(y_ticks1)
    ax2.set_ylabel("Propagation distance $z\,\mathrm{(m)}$")
  
    # -- customize y-axis
    ax2w.tick_params(axis='y', length=3, pad=2, top=False, labelleft=False)

    # -- customize y-axis
    ax2tw.tick_params(axis='y', length=3, pad=2, top=False, labelbottom=False)

    # SUBPLOT 2 - SET AXIS DETAILS ----------------------------------------------------

    ax2tw.text(0.02,0.88,r"(c)",color='k', weight='bold', fontsize=8,
        horizontalalignment='left',
        verticalalignment='bottom',
        bbox=dict(facecolor='white', alpha=1,edgecolor='k',boxstyle='square,pad=0.2'),
        transform=ax2tw.transAxes)
        
    ax22tw.text(0.02,0.88,r"(a)",color='k', weight='bold', fontsize=8,
        horizontalalignment='left',
        verticalalignment='bottom',
        bbox=dict(facecolor='white', alpha=1,edgecolor='k',boxstyle='square,pad=0.2'),
        transform=ax22tw.transAxes)

    ax2.text(0.02,0.94,r"(b)",color='k', weight='bold', fontsize=8,
        horizontalalignment='left',
        verticalalignment='bottom',
        bbox=dict(facecolor='white', alpha=1,edgecolor='k',boxstyle='square,pad=0.2'),
        transform=ax2.transAxes)

    ax2w.text(0.02,0.94,r"(d)",color='k', weight='bold', fontsize=8,
        horizontalalignment='left',
        verticalalignment='bottom',
        bbox=dict(facecolor='white', alpha=1,edgecolor='k',boxstyle='square,pad=0.2'),
        transform=ax2w.transAxes)

 

  
    # -- customize x-axis
    ax3w.tick_params(axis='x', length=3., pad=2, right=False)


    ax3tw.text(0.02,0.93,r"(e)",color='k', weight='bold', fontsize=8,
        horizontalalignment='left',
        verticalalignment='bottom',
        bbox=dict(facecolor='white', alpha=1,edgecolor='k',boxstyle='square,pad=0.2'),
        transform=ax3tw.transAxes)


    ax3w.text(0.02,0.93,r"(f)",color='k', weight='bold', fontsize=8,
        horizontalalignment='left',
        verticalalignment='bottom',
        bbox=dict(facecolor='white', alpha=1,edgecolor='k',boxstyle='square,pad=0.2'),
        transform=ax3w.transAxes)

    # -- customize x-axis
    y_lim3w = (-0.6,0.6)
    y_ticks3w = (-0.4,0,0.4)
    ax3w.tick_params(axis='y', length=3., pad=2, right=True, labelleft=False, labelright=True, left=False)
    ax3w.set_ylim(y_lim3w)
    ax3w.set_yticks(y_ticks3w)
    ax3w.set_ylabel(r"Delay $\tau$ (ps)")
    ax3w.yaxis.set_label_position("right")
   
    # -- customize y-axis
    ax3w.tick_params(axis='x', length=3, pad=2, top=False, labelleft=False)
    ax3w.set_xlabel(r"Angular frequency $\omega$ (rad/fs)")

    # -- customize x-axis
    ax2tw.tick_params(axis='x', length=3., pad=2, right=False, labelbottom=False)

    # -- customize y-axis
    y_lim2tw = (-35,1.1)
    y_ticks2tw = (-30,-15,0)
    ax2tw.tick_params(axis='y', length=3., pad=2, top=False, labelleft=False)
    ax2tw.set_ylim(y_lim2tw)
    ax2tw.set_yticks(y_ticks2tw)
    
    # -- customize x-axis
    ax22tw.tick_params(axis='x', length=3., pad=2, right=False, labelbottom=False)

    # -- customize y-axis
    ax22tw.tick_params(axis='y', length=3., pad=2, top=False)
    ax22tw.set_ylabel(r"$I$ (dB)")

   
    # -- customize x-axis
    ax3tw.tick_params(axis='x', length=3, pad=2, top=False, labelbottom=False)

    # -- customize y-axis
    ax3tw.tick_params(axis='y', length=3., pad=2, right=True, labelleft=False, labelright=True, left=False)
    ax3tw.set_ylabel(r"Delay $\tau$ (ps)")
    ax3tw.yaxis.set_label_position("right")



    ## SET RANGE PARAMETERS #################################################
    fileLoc = str(sys.argv[1]) 
    
    
    ## FETCH DATA FROM FILE #################################################
    (z, t, w,  Eps_w) = fetchNpz(fileLoc)    
   
    z = z*1e-6 # change scale to meters (cm)
    t = t*1e-3 # change scale to picoseconds (ps)   

    # CONVENIENT FUNCTIONS ##################################################################
    _norm = lambda x: x/x[0].max()
    _normX = lambda x: x/np.amax(x)
    _truncate = lambda x: np.where(x>1.e-20,x,1.e-20)
    _dB = lambda x: np.where(x>1e-20,10.*np.log10(x),10*np.log10(1e-20))


   
    # CONVENIENT FIELD DESCRIPTIONS ############################################
    

    utz = nfft.fft(Eps_w,axis=-1) #*np.exp(-1j*w*z[:,np.newaxis]*b10)
    utz2 = nfft.fft(Eps_w,axis=-1)
   

   
     # SUBFIGS ###################################################     
   
    # - Plot subfigure (b) - t-domain propagation
    im1 = subfigHelper(ax2,t, z, _dB(_truncate(_norm(np.abs(utz2)**2))), cmap='turbo')
   
    # - Plot subfigure (c) - spectra
    l1=ax2tw.plot(w,_dB(_truncate(np.abs(Eps_w[-1])**2/max(np.abs(Eps_w[0])**2))),  color='gray', linewidth=1.1, label='$|\mathcal{E}_{\omega}|^2_{out}$')
    l2=ax2tw.plot(w,_dB(_truncate(np.abs(Eps_w[0])**2/max(np.abs(Eps_w[0])**2))),  color='k', linewidth=1.1, label='$|\mathcal{E}_{\omega}|^2_{in}$')
    
   # I_w=_dB(_truncate(np.abs(Eps_w)**2/max(np.abs(Eps_w[0])**2)))
 
    # - Plot subfigure (d) - w-domain propagation
    im2=subfigHelper(ax2w,w, z, _dB(_truncate(_norm(np.abs(nfft.fftshift(Eps_w, axes=-1))**2))), cmap='turbo')

    # - Plot subfigure (e,f) - spectrograms
    (tau, wu, IXFROG) = getXfrogTrace(t,w,nfft.fft(Eps_w[0]),s0=300*(t[1]-t[0]),tauMin=-1, tauMax=1)
    (tau2, wu2, IXFROG2) = getXfrogTrace(t,w,nfft.fft(Eps_w[-1])/max(nfft.fft(Eps_w[0])),s0=300*(t[1]-t[0]),tauMin=-1, tauMax=1.5)
   

    imXFROG=ax3w.pcolorfast(wu, tau, _dB(_truncate(_normX(IXFROG[:-1,:-1]))),vmin=-55, cmap='turbo')
    imXFROG2=ax3tw.pcolorfast(wu, tau2, _dB(_truncate(_normX(IXFROG2[:-1,:-1]))),vmin=-55, cmap='turbo')
  
    # - Plot subfigure (a) - t-domain profiles
    l3=ax22tw.plot(t,_dB(_truncate(np.abs(nfft.fft(Eps_w[-1]))**2/max(np.abs(nfft.fft(Eps_w[0]))**2))),  color='gray', linewidth=1.1, label='$|\mathcal{E}|^2_{out}$')
    l4=ax22tw.plot(t,_dB(_truncate(np.abs(nfft.fft(Eps_w[0]))**2/max(np.abs(nfft.fft(Eps_w[0]))**2))),  color='k', linewidth=1.1, label='$|\mathcal{E}|^2_{in}$')
    colors=  plt.cm.Blues(np.linspace(0.5,1,3))   


    ### ZERO DISPERSION LINES #######################################
    ax2tw.axvline(x=1.51,linestyle='--',color='k') 
    ax2tw.axvline(x=2.55,linestyle='--',color='k')
  
    ax2w.axvline(x=1.51,linestyle='--',color='w') 
    ax2w.axvline(x=2.55,linestyle='--',color='w')
  
    ax3tw.axvline(x=1.51,linestyle='--',color='w') 
    ax3tw.axvline(x=2.55,linestyle='--',color='w')  

    ax3w.axvline(x=1.51,linestyle='--',color='w') 
    ax3w.axvline(x=2.55,linestyle='--',color='w')
    

    


   
    # COLORBARS AND LABELS #############################################

    cb2= set_colorbar(fig, im2, ax2tw,ax2tw)
    cb2.set_ticks((-30,-15,0))
    cb2.ax.set_title('$|\mathcal{E}_\omega|^2\,\mathrm{(dB)}$',fontsize=6., y=2)
    
    cb3= set_colorbar(fig, im2, ax22tw,ax22tw)
    cb3.set_ticks((-30,-15,0))
    cb3.ax.set_title('$|\mathcal{E}|^2\,\mathrm{(dB)}$',fontsize=6., y=2)

    cb = set_colorbar(fig, imXFROG, ax3tw,ax2tw)
    cb.set_ticks((-30,-15,0))
    cb.ax.set_title('$P_S^2\,\mathrm{(dB)}$', fontsize=6., y=2)   

    set_legend(ax2tw, l1+l2)
    set_legend(ax22tw, l3+l4)

  
 
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


