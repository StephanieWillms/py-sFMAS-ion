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
