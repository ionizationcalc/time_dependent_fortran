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
  real*8, dimension(30,30):: conce_ini,conce_nei,conce_ei
  integer, allocatable:: i_chemi_eigen(:)

  ! Local variables
  real*8:: time
  integer:: ii, ichemi, natom, jelem, ion
  real*8:: oversum

  ! Te cases
  real*8:: te_start, te_end, rhone
  real*8:: te_arr(2), ne_arr(2)
  integer:: itime

  ! Test output
  integer:: i
  integer, parameter:: ntime=10
  real*8, parameter:: dt0 = 0.75

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
  ! Example (1): Define temperature, density.
  !----------------------------------------------------------------------------
  te_start = 1.0e6       ! (K) used to set the initial charge state 
  te_end   = 10.0**6.8   ! (K)
  rhone = 1.0e7          ! (cm^-3)

  !
  ! (1-1) Set initial (and final) condition: equilibrium ionization states
  !
  te_arr(1) = te_end   ! (K)
  te_arr(2) = te_end     ! (K)
  ne_arr(1) = rhone      ! (cm^-3)
  ne_arr(2) = rhone      ! (cm^-3)

  do jelem = 1,nelem
    natom = natom_list(jelem)
    print*, "Initial State (natom=", natom, ")"

    ! Set initial charge states: case 1 or case 2
    ! case 1: Initial fraction in EI cases
    call func_equilibrium_eigen(i_chemi_eigen(jelem),&
    natom,te_start,conce_ini(1:natom+1,natom))

    ! case 2: Initial fraction using random number
    !do ion = 1, natom+1
    !  call RANDOM_NUMBER(conce_ini(ion, natom))
    !enddo
    !oversum = sum(conce_ini(1:natom+1, natom))
    !do ion = 1, natom+1
    !  conce_ini(ion, natom) = conce_ini(ion, natom)/oversum
    !enddo

    print*, 'sta: ', conce_ini(1:natom+1, natom)

    ! Get final charge states in EI cases
    call func_equilibrium_eigen(i_chemi_eigen(jelem),&
    natom,te_end,conce_ei(1:natom+1,natom))
    print*, 'end_ei: ',conce_ei(1:natom+1, natom)
  end do

  !
  ! (1-2) Enter the main loop for different time
  !
  open(12, file='test_onestep_ei.dat', form='unformatted')
  write(12)te_start, te_end, rhone
  write(12)conce_ini
  write(12)conce_ei
  write(12)ntime
  do itime = 1, ntime
    time = 10.0**(itime*dt0)
    call sub_solve_ionic_onestep(nelem, natom_list, i_chemi_eigen, &
            te_arr, ne_arr, time, &
            conce_ini, conce_nei)
    write(12)time
    write(12)conce_nei
  ! print results
    do jelem = 1, nelem
      natom = natom_list(jelem)
      print*, "itime=", itime, "time=", time
      print*, conce_nei(1:natom+1, natom)
    end do
  enddo
  close(12)

  !----------------------------------------------------------------------------
  ! Free memory
  !----------------------------------------------------------------------------
  deallocate(i_chemi_eigen)
end
