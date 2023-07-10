# potential-subroutine
This is the subroutine that reads the coefficients for the short and long range terms obtained with the fitting procedure programmed in fortran and included in my repository. This subroutine 
reads the coefficients and provides with the fitting that can be used in the dynamical calculations. 
Read carefullly the commented lines at the beginning of the file so that the values of the expansions correspond with the fitting of your system.
The fitting is based on the equations included in JCP, 139, 144305 (2013) used to fit the typical approach of a third atom to a rigid rotor molecule. The fitting is a multipolar expansion 
in Legendre polynomials. 

Do not hesitate to write to me if need some help: cristina.sanz@uam.es

To help understanding the use of the subroutine an example is given. main.f90 is the main program that reads fit.inp, where coordinates and ab initio points are read and used for comparison 
with the fitted coefficients. 

This subroutine so far is only working when the internuclear distance of the diatomic fragment is kept constant. The slight opening of the diatomic molecules is under construction.
