# Ti_model
A global model of the ion temperature (FORTRAN 77). 
Please note that this version is before integration into IRI code, so it does not include either Booker profiles (only a linear interpolation in altitude) or connections to IRI Te and Tn. 
iontif.for - main subroutine
iontif_test.for - a short test program: TEST OUTPUT:    1528.286       84.09644 
additional *.for files are also needed for compilation
The code was tested on Windows 10 64 bit and Intel Fortran compiler.
