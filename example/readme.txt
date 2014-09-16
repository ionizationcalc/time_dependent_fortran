Example:

Purpose:
This example show how to run this code to get the non-equilibrium ionization states.

Update:
2014-09-14

-------------------------------------------------------------------------------
1. Prepared files
-------------------------------------------------------------------------------
(a) 'streamline.dat' is unformatted f77 format, which including four lines:
The first line: n, one integer variable, defines the size of the following one dimentional arraies;
The second line: doulb real array, te(n) means the temperature history (unit: K);
The third line: double real array, ne(n) means the electron history (unit: cm^-3);
The 4th line: double real array, time(n) (unit: s).

(b) 'inicondition.dat'
If set starting from ionization equilibrium, then this file is ignored.

(c) 'input.txt'
Define all parameters required in calculations;

(d) all files including eigen matirx.
Such as the folder '../chianti_7_te501/'

-------------------------------------------------------------------------------
2. Compare and Run
-------------------------------------------------------------------------------
mpif90 -mkl -O3 time_depen_ionization.f90 -o nei_ionic.out

mpiexe -np 4 ./nei_ionic.out

-------------------------------------------------------------------------------
3. Output
-------------------------------------------------------------------------------
Double real, four dimensional array: conce(1:30, 1:30, 1:2, 1:number of records)
The first element is for charge states, the second element is for atomic number, the 3rd element defines the non-equilibrium or equilibrium case, and the 4th element is the index of records. For example, the non-equilibrium ion fraction of FeXXI for the second record is conce(21, 26,  1, 2).
