!*******************************************************************************
! subroutine: sub_solve_ionic_onestep
!*******************************************************************************
! Name:
!   sub_solve_ionic_onestep
! Purpose:
!   Performing one-step time-dependent ionization calculations. The time-step is
!   from the input argument 'dt_input'. If 'dt_input' is larger than the
!   requirement for time-dependent ionization progress, the sub-timestep will
!   be actived based on the Temperature tables. The time-step estimation using
!   density then can be removed in this routine.
! Authors:
!   Chengcai Shen
! Update:
!   2017-02-28
!   initialize: set the inputs using two time nodes.
!   2017-05-19
!   Bug fixed: Including the final time node, and using two-sub-steps at least.
!-------------------------------------------------------------------------------
subroutine sub_solve_ionic_onestep(nelem, natom_array, i_chemi_eigen, &
          te_arr, ne_arr, dt_input, &
          conce_ini, conce_nei)

use mod_eigen_matrix
implicit none
integer, parameter:: ntime = 2
integer:: nelem
integer:: natom_array(nelem),i_chemi_eigen(nelem)
real*8:: te_arr(ntime), ne_arr(ntime), dt_input
real*8:: conce_ini(30,30), conce_nei(30,30)
real*8:: conce_temp_ini(30,30), conce_temp(30)

!real*8:: dt_est_ne(ntime), dt_est_ne_min
real*8:: eval_arr(30), dt_est_elem(nelem)

! Local variables
integer:: i, n_inter, jelem
integer:: natom
real*8:: te_sta, te_end, te_now, dd_te
real*8:: ne_now, ne_s, ne_e, step_ne
real*8:: dt_sub, dt_now

!-------------------------------------------------------------------------------
! (1) estimate time-step using eigenvalues and density
!-------------------------------------------------------------------------------
! check ne*t & dne
!do i = 1,ntime
!  index_te = int(anint((dlog10(te_arr(i))-4.0d0)/dlogte)+1)
!  rho  = ne_arr(i)
!  do jelem = 1,nelem
!    natom = natom_array(jelem)
!    eval_arr(1:natom+1) = eigen(jelem)%evalues(1:natom+1,index_te)
!    eval_max = maxval(dabs(eval_arr(1:natom+1)))
!    dt_est_elem(jelem)  = chang_perct/(eval_max*rho)
!  end do
!         set saft factor 0.5
!  dt_est_ne(i) = 0.5d0*minval(dt_est_elem)
!end do
! dt_est_ne_min = minval(dt_est_ne)

!-------------------------------------------------------------------------------
! (2) enter main cycle
!-------------------------------------------------------------------------------
! Check if the dt_input is suitable?
te_sta = dlog10(te_arr(1))
te_end = dlog10(te_arr(ntime))
dd_te  = te_end - te_sta
ne_now = 0.5*(ne_arr(1)+ne_arr(ntime))

n_inter = int(dabs(dd_te)/dlogte) + 2
conce_temp_ini = conce_ini

if (n_inter .le. 2) then
  ! single timestep
  do jelem = 1,nelem
    natom = natom_array(jelem)

    ! half-step for te_arr(1)
    te_now = te_arr(1)
    dt_now = 0.5*dt_input
    call func_solveionization_eigen(i_chemi_eigen(jelem),&
    natom,te_now,ne_now,&
    conce_temp_ini(1:natom+1,natom),dt_now,conce_nei(1:natom+1,natom))

    ! half-step for te_arr(2)
    conce_temp_ini(1:natom+1,natom) = conce_nei(1:natom+1,natom)
    te_now = te_arr(2)
    dt_now = 0.5*dt_input
    call func_solveionization_eigen(i_chemi_eigen(jelem),&
    natom,te_now,ne_now,&
    conce_temp_ini(1:natom+1,natom),dt_now,conce_nei(1:natom+1,natom))

  end do
else
  ! multi-timestep
  dt_sub = dt_input/dfloat(n_inter-1)
  step_ne = (ne_arr(2)-ne_arr(1))/dfloat(n_inter-1)

  do i = 1, n_inter

    te_now = dfloat(i-1)*(te_arr(2)-te_arr(1))/dfloat(n_inter-1) + te_arr(1)

    if (i .eq. 1) then
      ne_now = ne_arr(1) + 0.5*step_ne
      dt_now = 0.5*dt_sub
    else if (i .eq. n_inter) then
      ne_now = ne_arr(2) - 0.5*step_ne
      dt_now = 0.5*dt_sub
    else
      ne_now = ne_arr(1) + (i-1)*step_ne
      dt_now = dt_sub
    endif

    do jelem = 1,nelem
      natom = natom_array(jelem)
      call func_solveionization_eigen(i_chemi_eigen(jelem),&
      natom,te_now,ne_now,&
      conce_temp_ini(1:natom+1,natom),dt_now,conce_nei(1:natom+1,natom))
    end do
    conce_temp_ini = conce_nei
  end do
endif
return
end subroutine sub_solve_ionic_onestep
