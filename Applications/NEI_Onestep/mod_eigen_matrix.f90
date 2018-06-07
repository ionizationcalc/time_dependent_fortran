!     *****************************************************************
!     Modules
!     2017-10-05:
!       Add new parameter: tetable_start.
!     2017-10-06
!       Read temperature grid parameters by subroutine
!       'su_rad_eigen_matrix.f90'.
!       Including: tetable_start, nte, dlogte
!     *****************************************************************
!     -----------------------------------------------------------------
      module mod_input_parameter
!     -----------------------------------------------------------------
!           atom list for calculation
      integer nelem
      integer,allocatable:: natom_list(:)
!           data path
      character(len=150) path_eigen
      end module mod_input_parameter

!     -----------------------------------------------------------------
      module mod_eigen_matrix
!     -----------------------------------------------------------------
!         err_ratetable is the accracy of temperature in ionization
!         and recombination rate tables. The unit is log10(Te).
      parameter(err_ratetable = 1.0d-03)
!         n_element is the number of chemical elements computed in this
!         program.
      parameter(n_element = 16)
!         char_element and index_element are name and atom index for
!         these chemical elements.
      character(len=2):: char_element(n_element)
      integer index_element(n_element)
      data index_element /1,2,6,7,8,10,11,12,13,14,16,18,20,24,26,28/
      data char_element /'H','He','C','N','O','Ne','Na','Mg','Al',&
      'Si','S','Ar','Ca','Cr','Fe','Ni'/
!         nte is the total number of temperature sample points in the
!         ionization rate table and eigenmatrix tables. e.g., Te ranges
!         from 10^3 K to 10^9 K.
      real*8:: tetable_start
      integer:: nte
      real*8:: dlogte
      real*8, allocatable:: te_grid_table(:)
!         eigen_type is used to save eigen values and eigen vectors
!         for each chemical element. c is ionization rate, and r is
!         recombination rate. eqis saved ion fractions for equilibrium
!         assumption.
      type eigen_type
      real*8,pointer:: eqis(:,:),&
                       evalues(:,:),&
                       evector(:,:,:),&
                       evector_invers(:,:,:),&
                       c(:,:),&
                       r(:,:)
      end type eigen_type
!         eigen includs all chemical elements.
      type(eigen_type),dimension(n_element):: eigen
      save eigen
      end module mod_eigen_matrix
