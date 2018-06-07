!*******************************************************************************
! Name:
!   main
! Purpose:
!   Test onestep ionization calculations routines.
! Update:
!   2017-02-28
!   initialized.
!   2017-03-03
!   Add two examples.
!   2017-10-05
!   Test for Jin-Yi's cases.
!   2018-04-19
!   Test element: Cr.
!   2018-04-20
!   Test EI.
!******************************************************************************
program main
  use mod_eigen_matrix
  use mod_input_parameter
  implicit none
  ! Ion charge states
  real*8, dimension(30,30):: conce_ini,conce_nei,conce_ei,conce_temp_ini
  integer, allocatable:: i_chemi_eigen(:)

  ! Local variables
  real*8:: dt
  integer:: ii, ichemi, natom, jelem

  ! Jin-Yi's cases
  real*8:: te_start, te_end, rhone
  real*8:: te_arr(2), ne_arr(2)
  integer:: itime

  ! Test EI output
  integer:: i
  integer, parameter:: nsamp=401
  real*8:: tesample(nsamp), fe8sample(nsamp)
  real*8:: fesample(27, nsamp)
  real*8:: conce_feei(27)

  !
  ! (1) Read input parameters
  !
  call sub_read_prameter
  allocate(i_chemi_eigen(nelem))
  i_chemi_eigen = 0
  do jelem = 1,nelem
    natom = natom_list(jelem)
  ! find the index 'ichemi' of atom in eigen matrix
    do ii = 1,16
      if(natom .eq. index_element(ii)) i_chemi_eigen(jelem) = ii
    end do
  ! check the atom index
    if(i_chemi_eigen(jelem) .eq. 0)then
      print*,'Error: unmatched atom number "',natom,'"'
      stop
    endif
  end do

  !
  ! (2) Load eigen tables
  !
  call sub_read_eigen_matrix(path_eigen)

  !----------------------------------------------------------------------------
  ! Example (1): Define temperature, density and dt.
  !----------------------------------------------------------------------------
  te_start = 101165.0    ! (K)
  te_end   = 10.0**6.50  ! (K)
  rhone = 1.0e8          ! (cm^-3)
  dt = 1000.0            ! (s)

  !
  ! (1-1) Set initial (and final) condition: equilibrium ionization states
  !
  te_arr(1) = te_end   ! (K)
  te_arr(2) = te_end     ! (K)
  ne_arr(1) = rhone      ! (cm^-3)
  ne_arr(2) = rhone      ! (cm^-3)

  do jelem = 1,nelem
    natom = natom_list(jelem)
    print*, 'jelem=',jelem, 'natom=', natom

    ! Initial fraction in EI cases
    call func_equilibrium_eigen(i_chemi_eigen(jelem),&
    natom,te_start,conce_ini(1:natom+1,natom))
    print*, conce_ini(1:natom+1, natom)

    ! Final fraction in EI cases
    call func_equilibrium_eigen(i_chemi_eigen(jelem),&
    natom,te_end,conce_ei(1:natom+1,natom))
  end do

  !
  ! (1-2) Enter the main loop
  !
  ! The first step: Te = Te_end(e.g.,10^6)
  call sub_solve_ionic_onestep(nelem, natom_list, i_chemi_eigen, &
            te_arr, ne_arr, dt, &
            conce_ini, conce_nei)
  conce_temp_ini = conce_nei

  !
  ! (1-3) Print results
  !
  do jelem = 1, nelem
    natom = natom_list(jelem)
    print*, '---------------------------------'
    print*, 'One step, atomic index = ',natom
    print*, '---------------------------------'
    print*, 'dt = ', dt
    print*, 'Initial time:'
    print*, conce_ini(1:natom+1, natom)
    print*, 'Final time: NEI'
    print*, conce_nei(1:natom+1, natom)
    print*, 'Final time: EI'
    print*, conce_ei(1:natom+1, natom)
  end do

  !
  ! (1-4) Save results
  !
  open(12, file='test_onestep_ei.dat', form='unformatted')
  write(12)nelem
  write(12)dt
  write(12)te_start
  write(12)te_end
  write(12)conce_ini
  write(12)conce_nei
  write(12)conce_ei
  close(12)

  !----------------------------------------------------------------------------
  ! Free memory
  !----------------------------------------------------------------------------
  deallocate(i_chemi_eigen)
end
