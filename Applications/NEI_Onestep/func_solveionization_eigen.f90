!     *****************************************************************
!     subroutine: func_solveionization_eigen
!     *****************************************************************
!     Input:
!           natom, atom index of this chemical element.
!           ichemi, index of element
!           rho, electron density.
!           te, temperature.
!           f0, ion fraction at initial time for this element.
!           dt, advanced time step.
!     Output:
!           ft, ion fraction at finnal time for this element.
!     Use:
!           mod_eigen_matrix, the eigenvalues and eigenvectors.
!     Call subroutines:
!           matrix_mul, matrix operation function.
!     Update:
!           2013-09-13
!     -----------------------------------------------------------------
      subroutine func_solveionization_eigen(ichemi,natom,te,rho,&
      f0,dt,ft)
      use mod_eigen_matrix
      implicit double precision(a-h,o-z)
      integer ichemi,natom
      real*8 dt,ene
      real*8 f0(natom+1),ft(natom+1)
      real*8 evals(natom+1),evect(natom+1,natom+1),&
             evect_invers(natom+1,natom+1)

!           Temprorary arrays
      real*8 diagona_evals(natom+1,natom+1)
      real*8 matrix_1(natom+1,natom+1),matrix_2(natom+1,natom+1)
      real*8 temparr(natom+1)
!     -----------------------------------------------------------------
!           Get eigen matrix according to the inputted te
!     -----------------------------------------------------------------
      index_te = int(anint((dlog10(te)-tetable_start)*(1.0/dlogte))+1)

!           approximate eigen matric
      evals(:) = eigen(ichemi)%evalues(:,index_te)
      evect(:,:) = eigen(ichemi)%evector(:,:,index_te)
      evect_invers(:,:) = eigen(ichemi)%evector_invers(:,:,index_te)

!     -----------------------------------------------------------------
!           Define a diagonal matrix including eigenvalues
!     -----------------------------------------------------------------
      diagona_evals = 0.0d0
      do ii = 1,natom+1
      diagona_evals(ii,ii) = dexp(evals(ii)*dt*rho)
      enddo
!     -----------------------------------------------------------------
!           matric operation
!     -----------------------------------------------------------------
      call matrix_mul(natom+1,evect,diagona_evals,matrix_1)

      call matrix_mul(natom+1,matrix_1,evect_invers,matrix_2)

      do j = 1,natom+1
      temparr(:) = matrix_2(:,j)*f0(:)
      ft(j) = sum(temparr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(ft(j) .le. 0.0d0)ft(j) = 0.0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      enddo
      return
      end subroutine func_solveionization_eigen

!     *****************************************************************
!     subroutine matrix_mul
!     *****************************************************************
      subroutine matrix_mul(n,a,b,c)
      implicit double precision(a-h,o-z)
      integer n
      real*8 a(n,n),b(n,n),c(n,n)
      real*8 temparr(n)
      do 10 i = 1,n
      do 10 j = 1,n
      temparr(:) = a(:,j)*b(i,:)
      c(i,j) = sum(temparr)
   10 continue
      return
      end subroutine matrix_mul
