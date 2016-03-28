!   Name: 
!       creat_eigen_matrix
!   Purpose:
!       write eigen matrix into files
!   Use:
!       Lapack mathlib
!       ionization and recombination rate table
!   Output:
!       unformatted data files
!   Update:
!       2013-11-19
!       2013-11-20
!       2014-05-14
!       Add hydrogen in elements list.
!       2014-05-21
!       Bugs modified: func_eqi did not return fraction for natom = 1.
!       2016-03-20
!       Read te_arr form data files.
!-------------------------------------------------------------------------------

program creat_eigen_matrix
    implicit none
!       file path and name
    character(len=350):: cpath,cfile
!       ioinization and recombination rate table
    integer nte
    real*8 dte
    real*8,allocatable:: te_arr(:),c_ori(:,:,:),r_ori(:,:,:)

!       atom name list
    integer:: nelems
    parameter(nelems = 15)
    integer,dimension(nelems):: arr_elemt
    character(len=2),dimension(nelems):: char_elemt

!       eigen matrix
    real*8,allocatable:: eqistate(:,:),&
                         eigenvalues(:,:),&
                         eigenvector(:,:,:),&
                         eigenvector_invers(:,:,:),&
                         c_rate(:,:),&
                         r_rate(:,:)
!       local arrays
    real*8,allocatable,dimension(:,:):: matrix_a,matrix_rate
    real*8,allocatable,dimension(:):: ionizat_rate,recombi_rate
    real*8,allocatable,dimension(:,:):: vr,vr_v
    real*8,allocatable,dimension(:):: wr,wi,work,ipiv
    real*8 dumm(1,1)
    integer lwork,info
    integer i,j,k,ite,icharg,ichemi,natom

!----------------------------------------------------------------------
!       (0) setting paramteres
!----------------------------------------------------------------------
    open(11,file='input.txt')
    read(11,*)cpath
    close(11)

    cfile=trim(cpath)//'ionrecomb_rate.dat'
    char_elemt(1:nelems)= (/'H','He','C','N','O','Ne','Na','Mg','Al','Si','S','Ar','Ca','Fe','Ni'/)
    arr_elemt(1:nelems) = (/1,2,6,7,8,10,11,12,13,14,16,18,20,26,28/)

!----------------------------------------------------------------------
!       (1) read ionization rate 'c' and recombination rate 'r'
!----------------------------------------------------------------------


!       open file and read    
    open(15,file=cfile,form='unformatted')
    read(15)nte
!       set variables
    allocate(te_arr(nte),c_ori(30,30,nte),r_ori(30,30,nte))

    read(15)te_arr
    read(15)c_ori
    read(15)r_ori
    close(15)
  
!--------------------------------------------------------------------
!       (2) computer eigen values
!--------------------------------------------------------------------
!       enter a cycyle for nelems chemical elements
!       the label is ichemi as following:
    do ichemi = 1,nelems
!           atom number is 'natom' and ionizatio states is up to 'natom+1'
        natom = arr_elemt(ichemi)
!           elemt = func_elemt_iv(natom)
        write(*,*),natom
    
!           define the size of local variables according to the atom 
!           number 'natom+1'
        allocate(eqistate(natom+1,nte),&
                 eigenvalues(natom+1,nte),&
                 eigenvector(natom+1,natom+1,nte),&
                 eigenvector_invers(natom+1,natom+1,nte),&
                 c_rate(natom+1,nte),&
                 r_rate(natom+1,nte))

!           local array for eigen calculations
        allocate(ionizat_rate(natom+1),&
                 recombi_rate(natom+1),&
                 matrix_a(natom+1,natom+1),&
                 matrix_rate(natom+1,natom+1),&
                 wr(natom+1),wi(natom+1),&
                 vr(natom+1,natom+1),vr_v(natom+1,natom+1),&
                 ipiv(natom+1))
        lwork = (natom+1)*(natom+1)*20
        allocate(work(lwork))

!      ----------------------------------------------------------------
!           inner cycle for different temperature
!      ----------------------------------------------------------------
       do ite = 1,nte
!          ------------------------------------------------------------
!              set a matrix including ionization and recombination rate
           do icharg = 1,natom
               ionizat_rate(icharg)  = c_ori(icharg,natom,ite)
               recombi_rate(icharg+1)= r_ori(icharg,natom,ite)
           enddo
           c_rate(:,ite) = ionizat_rate(:)
           r_rate(:,ite) = recombi_rate(:)

!          ------------------------------------------------------------
!              get ion fraction in equilibrium case
           call func_eqi(ionizat_rate,recombi_rate,natom,eqistate(1:natom+1,ite))

!          ------------------------------------------------------------
!              matrix operation
           call func_matrix_ionrecom(ionizat_rate,recombi_rate,natom,matrix_rate)

!              reset the order of matrx_rate to call dgeev
           do k = 1,natom+1
               matrix_a(1:natom+1,k) = matrix_rate(k,1:natom+1)
           enddo

           call dgeev('N','V',natom+1,matrix_a,natom+1,&
               wr,wi,dumm,1,vr,natom+1,work,lwork,info)
 
           if(info .eq. 0)then
!                   order eigenvectors side by side
               do j = 1,natom+1
                  i = j
                  vr_v(i,1:natom+1) = vr(1:natom+1,j)
               enddo
!                   save eigenvalues and eigenvectors
               eigenvalues(1:natom+1,ite)   = wr(1:natom+1)
               eigenvector(1:natom+1,1:natom+1,ite) = vr_v(1:natom+1,1:natom+1)

!                   calculate inverse eigenvectors
               call dgetrf(natom+1,natom+1,vr_v,natom+1,ipiv,info)
               call dgetri(natom+1,vr_v,natom+1,ipiv,work,natom+1,info)

               eigenvector_invers(1:natom+1,1:natom+1,ite) = vr_v(1:natom+1,1:natom+1)
           else
               print*,'Failure in DGEEV.  INFO = ', info
               stop
           endif
       enddo
    
!      ----------------------------------------------------------------
!          save results
!      ----------------------------------------------------------------
       open(15,file=trim(cpath)//trim(char_elemt(ichemi))//'eigen.dat',form='unformatted')
       write(15)eqistate
       write(15)eigenvalues
       write(15)eigenvector
       write(15)eigenvector_invers
       write(15)c_rate
       write(15)r_rate
       close(15)

!      ----------------------------------------------------------------      
!          release variables
       deallocate(eqistate,&
                  eigenvalues,&
                  eigenvector,&
                  eigenvector_invers,&
                  c_rate,&
                  r_rate)
       deallocate(matrix_a,matrix_rate,ionizat_rate,recombi_rate,&
                  wr,wi,vr,vr_v,work,ipiv)

    enddo
    
    deallocate(te_arr,c_ori,r_ori) 
!   -------------------------------------------------------------------
    print*,'normal stop.'
    stop
end 

subroutine func_matrix_ionrecom(ionizat_rate,recombi_rate,natom,matrix_rate)
    implicit none
    integer,intent(in):: natom
    real*8,dimension(natom+1),intent(in):: ionizat_rate,recombi_rate
    real*8,dimension(natom+1,natom+1),intent(out):: matrix_rate

!       local variables
    integer i_eq,iterm

!       initial matrix_rate
    matrix_rate = 0.0d0 

    do i_eq = 2,natom
    iterm = i_eq
    matrix_rate(iterm,  i_eq) = -(ionizat_rate(iterm) + recombi_rate(iterm))
    matrix_rate(iterm-1,i_eq) = ionizat_rate(iterm-1)
    matrix_rate(iterm+1,i_eq) = recombi_rate(iterm+1)
    enddo
!       boundary terms
    i_eq = 1
    iterm = i_eq
    matrix_rate(iterm,  i_eq) = -ionizat_rate(iterm)
    matrix_rate(iterm+1,i_eq) = recombi_rate(iterm+1)
    i_eq = natom + 1
    iterm = i_eq
    matrix_rate(iterm,  i_eq) = -recombi_rate(iterm)
    matrix_rate(iterm-1,i_eq) = ionizat_rate(iterm-1)
    return
end subroutine func_matrix_ionrecom

subroutine func_eqi(c,r,natom,conce)
    implicit none
   integer,intent(in):: natom
   real*8,dimension(natom+1),intent(in):: c,r
   real*8,dimension(natom+1),intent(out):: conce
!      local variables
   integer k
   real*8,dimension(natom+1):: f
!

!      set f1
   f(1) = 1.0d0
!      f2 = c1*f1/r2
   f(2) = c(1)*f(1)/r(2)
!
!      For Hydrogen, the following loop is not necessary.    
!
   if(natom .le. 1)then
      f(1) = 1.0d0/(1.0d0 + c(1)/r(2)) 
      f(2) = c(1)*f(1)/r(2)
      conce(1:2) = f(1:2)
      return
   endif
!         
!      For other elements:
!
!      f(i+1) = -(c(i-1)*f(i-1) - (c(i)+r(i)*f(i)))/r(i+1)
   do k = 2,natom-1
       f(k+1) = (-c(k-1)*f(k-1) + (c(k)+r(k))*f(k))/r(k+1)
   enddo
!      f(natom+1) = c(natom)*f(natom)/r(natom+1)
   f(natom+1) = c(natom)*f(natom)/r(natom+1)
!      f1 = 1/sum(f(*))
   f(1) = 1.0d0/sum(f)
!      f2 = c1*f1/r2
   f(2) = c(1)*f(1)/r(2)
!      f(i+1) = -(c(i-1)*f(i-1) - (c(i)+r(i)*f(i)))/r(i+1)
   do k = 2,natom-1
       f(k+1) = (-c(k-1)*f(k-1) + (c(k)+r(k))*f(k))/r(k+1)
   enddo
!      f(natom+1) = c(natom)*f(natom)/r(natom+1)
   f(natom+1) = c(natom)*f(natom)/r(natom+1)
!
   conce = f 
   return
end
