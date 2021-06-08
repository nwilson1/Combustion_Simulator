# Combustion_Simulator

Created - 12/10/2020

4/8/2021 - Added description

## DESCRIPTION ## 

This code is meant to simulate the combustion process with multicomponent mass diffusion.
It was originally written as part of the completion of a graduate course on combustion. 
While the code executes without a problem given appropriate inputs, the physical simulation
is incomplete and contains some as-yet unidentified computational errors.  

##  HOW TO RUN  ##

There are three files required to run the program: Comb_sim.py, Comb_aux.py, and therm.dat.  
All three files must be in the same directory from which the code is to be run.  

There is currently no GUI implemented and the file is not an executable, so it must be run 
from a command line interface such as a Unix terminal or a Python IDE.  Python 3.0 or later is required.
Apart from the included files, the only modules being loaded by the script are available with every 
Python distribution.

## FILES INCLUDED ##

Comb_aux.py    —   This file contains a variety of auxiliary functions called by the main script.
			This file calls no functions and has no input. attempting to run this file directly 
			will not result in any output.

therm.dat        —      This file contains all of the coefficients for the NASA 7-polynomial fit for 
			calculation of heat capacities and enthalpies as a function of temperature.
			The data was taken from the Extended Third Millennium
    			Thermodynamic Database for Combustion with updates from  Active
    			Thermochemical Tables by Elke Goos, Alexander Burcat and Branko Ruscic.

Comb_sim.py    —   This is the primary script that must be called by a Python interpreter.   
			Upon running the code the user will be presented with eight sequential input
			choices.  Hitting “Enter” chooses the default input for all prompts.  
      
These inputs are:
			
1. Fuel.  At the moment there is only methane and propane.  A choice of 0 for CH4 or 1 for C3H8.

2. Oxidizer.  A choice of 0 for Air or 1 for oxygen.
	
3. Equivalence ratio.  Any number between 0 and 1 may be given.

4. Pressure.  Any value in pascals may be given.

5. Initial temperature.  Any value in kelvin may be given.

6. Whether the multicomponent or mixture-averaged diffusion coefficients are used for computing diffusion velocity.  A choice of 0 for multicomponent or 1 for mixture-averaged.
		
7. Initial guess of the total mass flow rate m’’.  Any positive value 
				     may be given.

8. Number of Runge-Kutta steps to take for the solution.  Any
				     positive integer may be given.  The system is unstable past 
				     even a few R-K steps, so a choice of only 1 is strongly recommended.

## SAMPLE INPUT-OUTPUT ##

The code will output a number of different system and species properties as of the end of the last R-K iteration.  
Here is a sample output from a command line:

	\directory usr$ python Comb_sim.py

	Choose Fuel (CH4:0 C3H8:1) Enter for default (CH4) : 0

	Choose Oxidizer (Air:0 O2:1) Enter for default (Air) : 0

	Input equivalence ratio between 0 and 1.  Enter for default (1) : 1

	Input pressure in Pa. Enter for default (101 kPa) : 101325

	Input initial temperature in  K. Enter for default (298 K) : 300

	Multicompoent (0) or mixture averaged (1) diffusion? Enter for default (Averaged) : 1

	Input inital guess mass flux rate. Enter for default (0.001 kg/(s m^2)) : .001

	Warning -- Becomes unstable after multiple steps, only one step is recommended.
	Enter an integer for number of ODE steps. Enter for default (1) : 1



 ---- SYSTEM PARAMETERS FOLLOWING FINAL ITERATION ---- 



	Temperature in K :  300.0

	Total mass flux rate in kg/(s m^2) :  0.001

	Mass fracitions (F is fuel): 

				YCO  :  0.0
				YCO2  :  0.0
				YF  :  0.0552
				YH2O  :  0.0
				YN2  :  0.7247
				YO2  :  0.2201

	Specific heat capacities in J/(mol K):

				CH4  :  35.68
				O2  :  29.39
				N2  :  29.13
				CO2  :  37.22
				H2O  :  33.6
				CO  :  37.22

	Thermal conductivities in J/(m K):
	Note: This is only calculated if the species is present in the system.

				CH4  :  0.93
				O2  :  0.66
				N2  :  0.67
				CO2  :  0.0
				H2O  :  0.0
				CO  :  0.0
	\Mass production/use rates in kg/(s m^3) : 

				CH4  :  -1.5957032740346134e-34
				O2  :  -7.1609708324976e-34
				N2  :  0.0
				CO2  :  0.0
				H2O  :  8.959232837291517e-35
				CO  :  1.1143962776465631e-33

	Diffusion coefficients in m^2/s :
	Note: The diffusion coefficient is set to zero for any species not present in the system.
	Species order is [Fuel O2 N2 CO2 H2O CO]

	[0.02310664 0.02002711 0.02035977 0.         0.         0.        ]
