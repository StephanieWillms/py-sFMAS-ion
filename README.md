# py-sFMAS-ion
The software package `py-sFMAS-ion.py` solves a z-propagation equation for an analytic signal in single mode fibers. The code represents a basic version of the one that was used for the preparation of the dissertation of the author. Single mode propagation for higher order dispersion, and advanced nonlinear contribution such as the impact of ionization can be modeled.

In this repository the source code to reproduce many of the simulations that can be found in the dissertation []. An user example for a simulation from [] is included here.

An extended user guide for a similar software package containing e.g., the nonlinear Raman contribution, can be found: ..
and is suggested by the author as a comprehensive introduction to the topic of soliton and molecule state propagation in nonlinear media.

# Userguide

The top-level code `py-sFMAS-ion.py` imports all the modules  to execute the simulation. The results are stored in a python native compressed `npz` file in the directory labeled “\Data”.

- inputPulse.py: Module generates the initial condition from values set to “ini.py”
- splitStepSolver.py: Implements a simple split-step scheme with Runge-Kutta 4th order for the nonlinear part. Nonlinear coefficient n2 is assumed to be a constant.
- observables : stores generated data in `npz`-file format
- mediumDispersion: Implements 4 different propagation constants as a refractive index. Three of those (WD2017,Benabid,sinusoid) are found in the main text of [], while kagome is used for the test-calculation of the paper [].
- ini.py: Defines all necessary initial parameters for setting up the calculation domain and the system
- ioni.py: Implements the nonlinear part, responsible for ionization


The folder structure is given as:

```bash
├── py-sFMAS-ion.py
├── src
│   ├── ini.py
│   ├── mediumDispersion.py
│   ├── splitStepSolver.py
│   ├── ioni.py
│   ├── observables.py
│   └── inputPulse.py
├── Data
│   ├── x.npz
│   └── figure.py
├── LICENSE
└── README.md

```

# Example

The following exmple describes the direct superposition of two fundamental solitons, generating a multi-frequency state, as discussed in Ch. 6.3 of Ref [].

------ Computational Domain -------

tMax:5000.0

Nt:131072

zMax:100000

Nz:10000

------ Fiber Properties -------

Fiber Radius:20um

strut:0.8um

Gas pressure:3bar


""" Pulse 1 """
    p01_t0   = 40                      # pulsewidth
    p01_w0   = 3.44                    # centralfrequency
    p01_N    = 1.0                     # soliton order                 
    p01_tOff = 0.00                    # temporal offset between p01 and p02
   
""" Pulse 2 """
    p02_t0   = 40                      # pulsewidth
    p02_w0   = 4.54                    # centralfrequency    
    p02_N    = 1.0                     # order
    p02_tOff = 0.0                     # temporal offset between p01 and p02


# Generating figure

In addition to the source code for modelling the z-propagation, a script `figure.py` is included allowing for visualizing the results.

This script, located in the folder \Data, can be excecuted in the same folder by including a string system value, containing the calculated data that shall be presented. The generated figure results in Fig.1:
