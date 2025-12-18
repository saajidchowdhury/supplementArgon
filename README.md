Code and data used for CSA and DMC in our paper "Coordination-driven magic numbers in protonated argon clusters" 
by Saajid Chowdhury, Maria Judit Montes de Oca-Estevez, Florian Foitzik, Elisabeth Gruber, Paul Scheier, Pablo Villarreal, Rita Prosmiti, Tomas Gonzalez-Lezana, and Jesus Perez-Rios.

CSA and DMC implementations in MATLAB (with parallel computing toolbox) are used to find global minimum structures and ground state energies of Ar_nH^+ structures for n = 3,4,...,57.

fortranLinear.m evaluates the 4-body potential term for the Ar_2H^+-Ar interaction at 100,000 points: 1,000 different distances and 100 different angles.
It stores this data as interpolation.mat. In the paper, our interpolation uses 10,000 different distances and 1,000 different angles. 

interpolation.txt contains a Google drive link to the interpolation used in the paper, interpolation.mat (https://drive.google.com/file/d/1RTpr9AMMzGvffZoOfy31XLkDKe_VYC_h/view?usp=drive_link). 

fortranLinearTest.m tests the linear interpolation against the original Fortran implementation of the 4-body Ar_2H^+-Ar interaction term by Maria Judit Montes de Oca-Estevez, Rita Prosmiti, and Pablo Villarreal.

csa.m finds the global minimum of Ar_nH^+ for n = 3,4,5,...,12, using the many-body potential.
In the code, the variable "n" represents the number of particles, which is the number of argon atoms plus one for the proton. 

csa2B.m does the same, using the 2-body potential. 

dmcread.m runs Diffusion Monte Carlo to find the ground state energies of Ar_nH^+ for n = 3,4,5,...,12, starting with the structures found by csa.m stored in JuditRita4B and using the many-body potential.
In the code, the variable "n" represents the number of particles, which is the number of argon atoms plus one for the proton. 

dmcread2B.m does the same, starting with the structures found by csa2B.m stored in nArNH+2B and using the 2-body potential.

dmcPlots.m contains the experimental data collected by Florian Foitzik, Elisabeth Gruber, Paul Scheier, 
the potential energy global minima and ground state energies computed by CSA and DMC for the 2-body and many-body potentials, and the PIMC data by Tomas Gonzalez-Lezana. 
It plots the experimental data and evaporation energies for Figures 1 and 2 in the paper. 

The directory JuditRita4B contains the global minimum geometries of Ar_nH^+ structures for n = 3,...,57, with distances measured in Bohr radii, using the many-body potential.

The directory nArNH+2B contains the global minimum geometries of Ar_nH^+ structures for n = 1,...,59, with distances measured in Bohr radii, using the 2-body potential.
