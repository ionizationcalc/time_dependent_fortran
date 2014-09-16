Example:

Purpose:
This example show how to run this code to get the non-equilibrium ionization states.

Update:
2014-09-14

-------------------------------------------------------------------------------
1. Prepared files
-------------------------------------------------------------------------------
(a)‘te_ne_history.dat’ is used to define the temperature and electron density histories. The data structure including at least five lines:
The first line define how many records contained in this file.
Then each four lines define one record, which includes: size, Te, ne, and time.
size: n, integer, defines the size of the following one diminutional arrays;
Temperature: te(n), one dimensional double real array, defines the temperature history (unit: K);
Electron density: ne(n), one dimensional double real array, defines the electron history (unit: cm^-3);
Time: time(n), double real array, define time series (unit: s)

(b) 'inicondition.dat'
If set starting from ionization equilibrium, then this file is ignored.

(c) The folder containing all eigen matirx.
Such as the folder ‘../chianti_7_te501/‘, which is calculated according to chianti 7.

(d) 'input.txt'
Define all parameters required in calculations.
-------------------------------------------------------------------------------
2. Compile and Run
-------------------------------------------------------------------------------
mpif90 -mkl -O3 time_depen_ionization.f90 -o nei_ionic.out

mpiexe -np 4 ./nei_ionic.out

-------------------------------------------------------------------------------
3. Output
-------------------------------------------------------------------------------
Double real, four dimensional array: conce(1:30, 1:30, 1:2, 1:number of records)
The first element is for charge states, the second element is for atomic number, the 3rd element defines the non-equilibrium or equilibrium case, and the 4th element is the index of records. For example, the non-equilibrium ion fraction of FeXXI for the second record is conce(21, 26,  1, 2).
