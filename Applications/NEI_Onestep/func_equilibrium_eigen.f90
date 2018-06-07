!     *****************************************************************
!     subroutine: func_equilibrium_eigen
!     *****************************************************************
!     Input: ichemi, natom,te
!     Output: ion charge state in equilibrium case
!     Use: mod_eigen_matrix
!     Write: 2013-09-14
      subroutine func_equilibrium_eigen(ichemi,natom,te,ft)
      use mod_eigen_matrix
      implicit none
      integer natom,ichemi
      real*8 te,ft(natom+1)

!     Local varialbes
      integer:: i
      integer:: index_te, index_2(1)
      real*8, dimension(3):: terange, dterange
!     -----------------------------------------------------------------
!           Get eigen matrix according to the inputted te
!     -----------------------------------------------------------------
!     Pre-compute:
      index_te = int(anint((dlog10(te)-tetable_start)*(1.0/dlogte))+1)

!     Confirm
!      if ((index_te .ge. 2) .or. (index_te .le. nte-1)) then
!        terange(1:3) = te_grid_table(index_te-1:index_te+1)
!        dterange(:) = DABS(terange(:) - te)
!        index_2 = MINLOC(dterange)

        ! Reset the index_te to the nearest Temperature node
!        index_te = index_te + index_2(1) - 2
!      endif

      ft(1:natom+1) = eigen(ichemi)%eqis(1:natom+1,index_te)
      return
      end subroutine func_equilibrium_eigen
