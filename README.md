# Wakefield in Axi-Symmetric Dielectric-Lined Waveguides
## Synopsis

This set of PYTHON scripts provide a toolkit to compute the Green's function
associated to electromagnetic wakefield produced in an axi-symmetric dielectric
lined waveguide.The present version of the code used is based on the formalism
described in [K.-Y. Ng, Phys. Rev. D42, 1819 (1990)]. 

## Code Example

The simple implementation in test_Ng_long.py and test_Ng_trans.py demonstrate
the use of the code to reproduce the wakefield computed in the paper [M Rosing 
and W. Gai, et al., Phys. Rev. D42, 1829 (1990)].  In addition these two code 
compare the computed Green's function withthe one obtained with a trusted commercial 
program. 

The script MakeWake4Eleg.py generate a Green's function formatted to be compatible 
with ELEGANT longitudinal wake element.  

test_Ng.py provide an example of longitudinal wake (monopole) computation along with its
convoultion with a given charge distribution (in the present case a Gaussian distribution)         

pitz.py produced the Green's function and saved is as an ASTRA-compatible file for use 
in the &wakefield element.This program was employed to perform simulations reported in https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.044801 

## Tests

executing the command line 
python test_Ng_long.py should produce a plot displaying the Green's function 
associated to the sum of the first four monopole modes. 

## Contributors 

P. Piot (NIU and Fermilab), J. Wang (NIU), F. Lemery (NIU, now DESY), 

## Citing this work

Please cite this code as [use the tab on the right to generate bibtex format]:
P. Piot, F. Lemery, and J. Wang,  DiWakeCyl: a computer program to compute Green's functions for dielectric-lined cylindrical-symmetric waveguides, PYTHON program available from GitHub at  https://github.com/NIUaard/DiWakeCyl 
