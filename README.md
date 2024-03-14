**Fortran Code for Willard and Chandler (2010) J. Phys. Chem.**

Written by EagleLeslie for VASP post-processing. Only computes equations 1 - 5.

**Input files:** 
  - VASP XDATCAR
  - INTFC (optional)

**Output files:**
  - intf1
  - intf2
  - ADATCAR
  - DENFIELD
  - NPHI

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

**File information**

The XDATCAR file should only have the cell information printed at the header (i.e. ISIF = 2). This program does not support a changing cell shape (i.e. NPT; ISIF = 3 or ISIF = 4), and only supports an orthogonal coordinate system. The program assumes that the long axis of the system is in the z-axis. For example: a == b & c > a.

The simulation setup should be done in a way where there is a clear interface between the two phases. This is due to the program only finding two Gibbs Dividing Surfaces in the system. A sample system is visualized below.

If the file INTFC exists in the directory, it should have the desired interface constant value. The default is roughly 1/2 the bulk density of the system and does not require the INTFC file. 

Calculating the density field is the most expensive portion of the simulation, so the density field and density field gradient are output to files DENFIELD and NPHI, respectively. This enables the program to read in the files if they exist and saving time on the computation if the user would like to find the Gibbs Dividing Surface at different interface constant values. 

The user should keep this in mind when calculating the Gibbs Dividing Surface for another system and make sure that any exisitng INTFC, DENFIELD, and/or NPHI files are removed from the directory.

The file ADATCAR contains the same information as XDATCAR but with an additional column printed after the atomic positions. The last column represents the distance of the atom with respect to the nearest interface.

**Comile**
To compile the code, enter the command "make". 
Compiler flags in makefile use gfortran. MacOSX library required and uses the flag: GLIB = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib

Once compiled, enter the command "main" to run. Example command line script to run:

      python clean_wc_files.py
      make clean
      make
      ./main

Enter "make clean" to remove any *.o object files.

Output after entering the command "make" should look similar to:

      gfortran -ffixed-line-length-136 -fdefault-real-8 -Wall -fno-automatic -Ofast -g -fbounds-check -fbacktrace -finit-real=nan -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -c main.f90
      gfortran -ffixed-line-length-136 -fdefault-real-8 -Wall -fno-automatic -Ofast -g -fbounds-check -fbacktrace -finit-real=nan -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -o main leibniz.o
      find_params.o read_xdat.o common.o calc_denfield.o brakzero.o interface.o proximity.o main.o
      Model compiled


**Additional files**

The main.py file contains an example python script for automating a chosen interface value with respect to the atoms in the system. clean_wc_files.py contains a python script the cleans up output files.

**Example interface figures for visualization purposes**

_Fe-MgSiO3 system at 5500 K and 60 GPa after 15 ps_
  - Fe == White
  - Mg == Green
  - Si == Blue
  - O == Red
    
![tf_5000K](https://github.com/EagleLeslie/willchand/assets/120432106/06ed030b-93b2-40eb-ad0a-1e26146de2e7)

