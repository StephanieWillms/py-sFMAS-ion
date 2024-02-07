# py-sFMAS-ion
The software package `py-sFMAS-ion` solves a z-propagation equation for an analytic signal in single mode fibers. The code represents a basic version of the one that was used for the preparation of the dissertation of the author. Single mode propagation for higher order dispersion, and advanced nonlinear contribution such as the impact of ionization can be modeled.

In this repository the source code to reproduce many of the simulations that can be found in the dissertation [1]. An user example for a simulation from [2] is included here.

An extended user guide for a similar software package containing e.g., the nonlinear Raman contribution, can be found as methods manuscript in Ref.[3] and the respective software in Ref.[4]
and is suggested by the author as a comprehensive introduction to the topic of soliton and molecule state propagation in nonlinear media.

# Userguide

The top-level code `main.py` imports all the modules to execute the simulation. The results are stored in a python native compressed `npz` file in the directory labeled “\Data”.

- inputPulse.py: Module generates the initial condition from values set to “ini.py”
- splitStepSolver.py: Implements a simple split-step scheme with Runge-Kutta 4th order for the nonlinear part. Nonlinear coefficient n2 is assumed to be a constant.
- observables : stores generated data in `npz`-file format
- mediumDispersion: Implements 4 different propagation constants as a refractive index. Three of those (WD2017,Benabid,sinusoid) are found in the main text of [1], while kagome is used for the test-calculation of the paper [2].
- ini.py: Defines all necessary initial parameters for setting up the calculation domain and the system
- ioni.py: Implements the nonlinear part, responsible for ionization


The folder structure is given as:

```bash
├── py-sFMAS-ion.py
├── src
│   ├── ini.py
│   ├── mediumDispersion.py
│   ├── splitStepSolver.py
│   ├── nonLin.py
│   ├── observables.py
│   └── inputPulse.py
├── Data
│   ├── x.npz
│   ├── figure_Chang.py
│   └── figure.py
├── LICENSE
└── README.md

```

# Examples

## Example #1

The following exmple describes the propagation of a direct superposition of two fundamental solitons, generating a multi-frequency state, as discussed in Ch. 2.2 (see also Fig.2.2 therein) of Ref [1]. The used fiber is 'WD2017'.

```bash
Pulsetype:sech
Pulse energy:2.4999999999999998e-06J
Ionization:off

------ Computational Domain -------
tMax:2000.0
Nt:65536
zMax:150000
Nz:15000

------ Fiber Properties -------
Fiber Radius:8um
strut:0.8um
Gas pressure:3bar

------ Parameters -------
Nonlinearity:2.3545599999999996e-07um^2/W
Involved Pulses:two
```

```bash
""" Pulse 1 """
    p01_t0   = 30                      # pulsewidth
    p01_w0   = 1.2                     # centralfrequency
    p01_N    = 1.0                     # soliton order         
    p01_tOff = 0.00                    # temporal offset between p01 and p02
   
""" Pulse 2 """
    p02_t0   = 30                      # pulsewidth
    p02_w0   = 2.939                   # centralfrequency    
    p02_N    = 1.0                     # order
    p02_tOff = 0.0                     # temporal offset between p01 and p02
```

Running the software with the above shown parameters results in the generation of a file: `/Data/obs_sFMAS-RK4-da_Nt65536_Nz15000_p01_t030.000_w01.200_N1.000_p02_t030.000_w02.939_N1.000_tOff1.000__ioni_off.npz`


## Example #2

The following exmple describes propagation of a soliton pulse under the impact of ionization. The example is taken from Ref. [2] in order to test the implemented ionization code. A detailed description can be found in the appendix D of Ref. [1] and Fig. D.3 therein. The used fiber is 'kagome'.

```bash
Pulsetype:sech
Pulse energy:2.4999999999999998e-06J
Ionization:on
------ Computational Domain -------
tMax:2000.0
Nt:65536
zMax:100000
Nz:20000
------ Fiber Properties -------
Fiber Radius:8um
strut:0.8um
Gas pressure:3bar
... finished
------ Parameters -------
Nonlinearity:4.1599999999999997e-07um^2/W
Involved Pulses:one

```

```bash
""" Pulse 1 """
    p01_t0   = 30                      # pulsewidth
    p01_w0   = 2.356                   # centralfrequency
    p01_N    = 1.0                     # soliton order         
    p01_tOff = 0.00                    # temporal offset between p01 and p02
   
```
Running the software with the above shown parameters results in the generation of a file: `/Data/obs_sFMAS-RK4-da_Nt65536_Nz20000_p01_t030.000_w02.356_N1.000__ioni_on.npz`

# Generating figure

In addition to the source code for modelling the z-propagation, the scripts `figure.py` and `figure_Chang.py` are included allowing for visualizing the results.

This script, located in the folder \Data, can be excecuted in the same folder by including a string system value, containing the calculated data that shall be presented. 

The generated result from Example #1 evaluated with `figure.py` results in ![Fig.1](FIGS/obs_sFMAS-RK4-da_Nt65536_Nz15000_p01_t030.000_w01.200_N1.000_p02_t030.000_w02.939_N1.000_tOff1.000__ioni_off_T.pdf)

The generated result from Example #2 evaluated with `figure_Chang.py` results in ![Fig.2](FIGS/obs_sFMAS-RK4-da_Nt65536_Nz20000_p01_t030.000_w02.356_N1.000__ioni_on_T.pdf)

# References

[1] S. Willms, Dissertation, 2024

[2] W. Chang, A. Nazarkin, J. C. Travers, J. Nold, P. Hölzer, N. Y. Joly, and P. S. J. Russell, “Influence of ionization on ultrafast gas-based nonlinear fiber optics,” Opt. Exp., vol. 19, p. 21018, 2011.

[3] O. Melchert, A. Demircan, A python package for ultrashort optical pulse propagation in terms of forward models for the analytic signal, arXiv:2104.11649, https://arxiv.org/abs/2104.11649.

[4] O. Melchert, py-fmas,  GitHub repository, [https://github.com/charlespwd/project-title](https://github.com/omelchert/py-fmas)https://github.com/omelchert/py-fmas


