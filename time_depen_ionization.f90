!     *****************************************************************
!     Modules
!     *****************************************************************
!     -----------------------------------------------------------------
      module mod_input_parameter
!     -----------------------------------------------------------------
!           switch_adap_step, adaptive time step
      logical switch_eqi
      logical switch_log
!           parameter mode:
!           = 1, use original time array
!           = 2, according to accuracy value, comparing time-forward and
!           time-backward scheme
!           = 3, according sample points of eigen matric, check Te and
!           ne history and define time-step.
      integer mode
!           error, time-dependent ionization step
      real*8 accuracy
      real*8 chang_perct
!           atom list for calculation
      integer nelem
      integer,allocatable:: natom_list(:)
!           data path
      character(len=150) path_data,path_eigen,jobname
      character(len=150) file_streamline,file_ionfraction,&
          file_inicondition,file_log
      end module mod_input_parameter

!     -----------------------------------------------------------------
      module mod_eigen_matrix
!     -----------------------------------------------------------------
!         err_ratetable is the accracy of temperature in ionization
!         and recombination rate tables. The unit is log10(Te).
      parameter(err_ratetable = 1.0d-03)
!         n_element is the number of chemical elements computed in this
!         program.
      parameter(n_element = 15)
!         char_element and index_element are name and atom index for 
!         these chemical elements.
      character(len=2):: char_element(n_element)
      integer index_element(n_element)
      data index_element /1,2,6,7,8,10,11,12,13,14,16,18,20,26,28/
      data char_element /'He','He','C','N','O','Ne','Na','Mg','Al',&
      'Si','S','Ar','Ca','Fe','Ni'/
!         nte is the total number of temperature sample points in the 
!         ionization rate table and eigenmatrix tables. The Te ranges 
!         from 10^4 K to 10^9 K.
      parameter(nte=501)
      parameter(dlogte=5.0d0/dfloat(nte-1))
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

!     -----------------------------------------------------------------
      module mod_mpi_global
!     -----------------------------------------------------------------
      integer n_nodes,myid,mpi_err
      integer, parameter :: my_root=0
      end module mod_mpi_global
      
!     *****************************************************************
!     subroutine: sub_mpiinitial
!     *****************************************************************
      subroutine sub_mpiinitial
!  
      use mpi
      use mod_mpi_global
      implicit none
      integer status(MPI_STATUS_SIZE),ierr
!         do the mpi init stuff
      call MPI_INIT( ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, n_nodes, ierr )
      call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
      end subroutine sub_mpiinitial

!     *****************************************************************
!     Main program: time_depen_ionization
!     *****************************************************************
      program time_depen_ionization
      use mod_eigen_matrix
      use mod_input_parameter
      use mod_mpi_global
      use mpi
      implicit double precision(a-h,o-z)

      character(len=10)systime_sta,systime_end
!           History of temperature and density
      real*8,allocatable,save:: te(:),rho(:),time(:)
!           Ion charge states
      real*8 conce_ini(30,30),&
             conce_nei(30,30),conce_eqi(30,30)
      real*8,allocatable:: conce_arr(:,:,:,:),conce_root(:,:,:,:)
      integer,allocatable:: i_chemi_eigen(:)
!           MPI variable
      character(len=4)char_myid

!     -----------------------------------------------------------------
!           MPI initial
!     -----------------------------------------------------------------
      call sub_mpiinitial

!     -----------------------------------------------------------------
!           input parameters:
!     -----------------------------------------------------------------
      call sub_read_prameter

      allocate(i_chemi_eigen(nelem))
      i_chemi_eigen = 0
      do jelem = 1,nelem
      natom = natom_list(jelem)
!           find the index 'ichemi' of atom in eigen matrix
      do ii = 1,15
      if(natom .eq. index_element(ii)) i_chemi_eigen(jelem) = ii
      enddo
!         check the atom index
      if(i_chemi_eigen(jelem) .eq. 0)then
      print*,'Error: unmatched atom number "',natom,'"'
      stop
      endif
      enddo

!     -----------------------------------------------------------------
!         split datafile 
!     -----------------------------------------------------------------
      if(myid .eq. 0)then
      call sub_split_data(n_nodes,nlines_total,nlines_pnode,&
          file_streamline,file_inicondition,switch_eqi)
      endif
      call mpi_barrier(mpi_comm_world,ierr)

      call mpi_bcast(nlines_pnode, 1, mpi_integer,my_root,&
           mpi_comm_world,ierr)
!      if(switch_log)open(51,file=file_log)

!     -----------------------------------------------------------------
!         print system time
!     -----------------------------------------------------------------
      call date_and_time(TIME=systime_sta)

!     -----------------------------------------------------------------
!         Read eigen_matrix 
!     -----------------------------------------------------------------
      call sub_read_eigen_matrix(path_eigen)

!     -----------------------------------------------------------------
!           Read the history of Te and rho
!     -----------------------------------------------------------------
      write(char_myid,'(i04)')myid
      my_file = myid + 21
      open(my_file, file=trim(file_streamline)//char_myid, &
          form='unformatted',action='read')

!           the number of streamlines
      read(my_file)npoints
!           screen print message:
      write(*,'(a5,i03,a14,i05)'),'node:',myid,'; streamlines:',npoints

!           datafile: save results
      allocate(conce_arr(30,30,2,nlines_pnode))
      call mpi_barrier(mpi_comm_world,ierr)

!     -----------------------------------------------------------------
!           Enter cycle label 10, each streamline
!     -----------------------------------------------------------------
      do 10 ip = 1,npoints
      read(my_file)ntime
      allocate(te(ntime),rho(ntime),time(ntime))
!           read Te, rho hsitory      
      read(my_file)te
      read(my_file)rho
      read(my_file)time
      if(.not. switch_eqi)then
          read(my_file)conce_ini
      endif

!           chose a numerical method
      select case(mode)
          case (1)
          call sub_solve_ioni_mode1(nelem,natom_list,i_chemi_eigen,&
          ntime,te,rho,time,&
          conce_ini,conce_nei,conce_eqi,&
          switch_eqi)
          case (2)
          call sub_solve_ioni_mode2(nelem,natom_list,i_chemi_eigen,&
          ntime,te,rho,time,&
          conce_ini,conce_nei,conce_eqi,&
          switch_eqi,&
          accuracy)
          case (3)
          call sub_solve_ioni_mode3(nelem,natom_list,i_chemi_eigen,&
          ntime,te,rho,time,&
          conce_ini,conce_nei,conce_eqi,&
          switch_eqi,&
          chang_perct)
      end select

!           results of one points, including all elements
      conce_arr(1:30,1:30,1,ip) = conce_nei(1:30,1:30)
      conce_arr(1:30,1:30,2,ip) = conce_eqi(1:30,1:30)
!           release te, rho and time for this streamline
      deallocate(te,rho,time)
   10 continue

!     -----------------------------------------------------------------
!           gather from each node
!     -----------------------------------------------------------------
      call mpi_barrier(mpi_comm_world,ierr)
      if(myid .eq. my_root)then
          allocate(conce_root(30,30,2,nlines_pnode*n_nodes))
      endif
      call mpi_gather(conce_arr, 30*30*2*nlines_pnode, mpi_real8, &
          conce_root, 30*30*2*nlines_pnode, mpi_real8, my_root, &
          mpi_comm_world, ierr)

      deallocate(conce_arr)

!     -----------------------------------------------------------------
!           close my datafile for each node
!     -----------------------------------------------------------------
      close(my_file,status='delete')

!     -----------------------------------------------------------------
!           my_root write results into file
!     -----------------------------------------------------------------
      if(myid .eq. my_root)then
          open(7,file=trim(file_ionfraction),form='unformatted')
          write(7)nlines_total
          do iline = 1,nlines_total
              write(7)conce_root(1:30,1:30,1,iline)
              write(7)conce_root(1:30,1:30,2,iline)
          enddo
          close(7)
          deallocate(conce_root)
      endif

!     -----------------------------------------------------------------
!	    print system time
!     -----------------------------------------------------------------
      call date_and_time(TIME=systime_end)
      if(myid .eq. my_root)then
          print*,'Compt time:',systime_sta,'-',systime_end
      endif

!     -----------------------------------------------------------------
!           MPI close
!     -----------------------------------------------------------------
      call MPI_Finalize(ierr)

      end

!     *****************************************************************
!     subroutine: sub_read_prameter
!     *****************************************************************
      subroutine sub_read_prameter
      use mod_input_parameter
      implicit none
      character (len=150)charhead
      open(15,file="input.txt")

      read(15,fmt="(a150)")charhead
      read(15,*)switch_eqi
      write(*,fmt="(a11,l2)")'Switch_eqi=',switch_eqi

      read(15,fmt="(a150)")charhead
      read(15,*)switch_log
      write(*,fmt="(a11,l2)")'Switch_log=',switch_log

      read(15,fmt="(a150)")charhead
      read(15,*)mode
      write(*,fmt="(a5,i2)")'Mode=',mode

      read(15,fmt="(a150)")charhead
      read(15,*)accuracy
      write(*,fmt="(a18,e10.3)")"Accuracy in mode2=",accuracy

      read(15,fmt="(a150)")charhead
      read(15,*)chang_perct
      write(*,fmt="(a18,e10.3)")"Accuracy in mode3=",chang_perct

      read(15,fmt="(a150)")charhead
      read(15,*)nelem
      allocate(natom_list(nelem))
      read(15,*)natom_list
      write(*,fmt="(a14,30i02)")"Elements list:",natom_list

      read(15,fmt="(a150)")charhead
      read(15,*)path_eigen
      read(15,*)file_streamline
      read(15,*)file_inicondition 
      read(15,*)path_data
      read(15,fmt="(a150)")charhead
      read(15,*)jobname
      write(*,fmt="(a8,a30)")"Jobname:",trim(jobname)

      close(15)

      chang_perct = -log(chang_perct)
      file_ionfraction = trim(path_data)//trim(jobname)//".dat"
      file_log = trim(path_data)//trim(jobname)//".log"
      return
      end subroutine sub_read_prameter


!     *****************************************************************
!     subroutine sub_split_data
!     *****************************************************************
      subroutine sub_split_data(n_nodes,nlines_total,nlines_pnode,&
          file_streamline,file_inicondition,switch_eqi)
      implicit double precision(a-h,o-z)
      integer n_nodes,nlines_total,nlines_pnode
      character(len = 150)file_streamline,file_inicondition
      character(len = 4)char_myid
      real*8,allocatable:: te(:),rho(:),time(:)
      real*8 conce_ini(30,30)
      logical switch_eqi

      open(11,file=trim(file_streamline),action="read")
!         the total number of streamlines
      read(11,*)nlines_total

!         load initical fraction if switch_eqi = .f.
      if(.not. switch_eqi)then
          open(13,file=trim(file_inicondition),form="unformatted",&
              action="read")
          read(13)n_ini_record
          if(n_ini_record .ne. nlines_total)stop
      endif

!         the number of streamlines in each process
      nlines_pnode = nlines_total/n_nodes
      fnlines_pnode = dfloat(nlines_total)/dfloat(n_nodes)
      if(fnlines_pnode .gt. dfloat(nlines_pnode))then
          nlines_pnode = nlines_pnode + 1
      endif
!     -----------------------------------------------------------------
!         (1) normal nodes
!     -----------------------------------------------------------------    
      do 20 id_now = 1,n_nodes-1
!         openfile for the current process
      myid = id_now-1
      write(char_myid,'(i04)')myid
      open(12,file=trim(file_streamline)//char_myid,action="read")
      write(12,*)nlines_pnode

      do 21 ip = 1,nlines_pnode    
!         read one record from file '11'            
      read(11,*)ntime
      allocate(te(ntime),rho(ntime),time(ntime))      
      read(11,*)te
      read(11,*)rho
      read(11,*)time
!         write one record into file '12'
      write(12)ntime
      write(12)te
      write(12)rho
      write(12)time
      deallocate(te,rho,time)
!         write inicondition
      if(.not. switch_eqi)then
          read(13)conce_ini
          write(12)conce_ini
      endif
   21 continue
      close(12)
   20 continue

!     -----------------------------------------------------------------
!         (2) the last node
!     -----------------------------------------------------------------
      myid = n_nodes-1
      write(char_myid,'(i04)')myid
      open(12,file=trim(file_streamline)//char_myid,form='unformatted')
      ip = (n_nodes-1)*nlines_pnode + 1

      write(12) nlines_total-ip+1

      do 30 icycle = ip, nlines_total
!         read one record from file '11' 
      read(11,*)ntime
      allocate(te(ntime),rho(ntime),time(ntime))      
      read(11,*)te
      read(11,*)rho
      read(11,*)time
!         write one record into file '12'
      write(12)ntime
      write(12)te
      write(12)rho
      write(12)time
      deallocate(te,rho,time)
!         write inicondition
      if(.not. switch_eqi)then
          read(13)conce_ini
          write(12)conce_ini
      endif
   30 continue
      close(12)

      close(11)
      if(.not. switch_eqi)close(13)
      return 
      end subroutine sub_split_data

!     *****************************************************************
!     subroutine: sub_read_eigen_matrix
!     *****************************************************************
      subroutine sub_read_eigen_matrix(path_eigen)
      use mod_eigen_matrix
      implicit double precision(a-h,o-z)
      character(len=150)datafile,path_eigen
!     -----------------------------------------------------------------
!
      do 10 ichemi = 1,n_element
      natom = index_element(ichemi)
!
      allocate(eigen(ichemi)%eqis(natom+1,nte),&
               eigen(ichemi)%evalues(natom+1,nte),&
               eigen(ichemi)%evector(natom+1,natom+1,nte),&
               eigen(ichemi)%evector_invers(natom+1,natom+1,nte),&
               eigen(ichemi)%c(natom+1,nte),&
               eigen(ichemi)%r(natom+1,nte))
!
      datafile = trim(path_eigen)//trim(char_element(ichemi))//&
      'eigen.dat'
      open(11,file=trim(datafile),form='unformatted',action='read')
      read(11)eigen(ichemi)%eqis
      read(11)eigen(ichemi)%evalues
      read(11)eigen(ichemi)%evector
      read(11)eigen(ichemi)%evector_invers
      read(11)eigen(ichemi)%c
      read(11)eigen(ichemi)%r
      close(11)
   10 continue
!
      return
      end subroutine sub_read_eigen_matrix

!     *****************************************************************
!     subroutine: sub_solve_ioni_mode1
!     *****************************************************************
!     Input:
!           nelem,
!           natom_array,
!           i_chemi_eigen,
!           ntime,
!           te_array,
!           rho_array,
!           time_array,
!           conce_ini,
!           switch_eqi
!     Output:
!           conce_nei,
!           conce_eqi
!     Use:
!           mod_eigen_matrix
!           mod_input_parameter
!     Call:
!           func_equilibrium_eigen
!           func_solveionization_eigen
!     Update:
!           2013-09-14
!           2013-09-24
!           add adaptive time-step.
!     Update: 
!           2013-10-22
!           remove use mode_input_parameter
!     -----------------------------------------------------------------
      subroutine sub_solve_ioni_mode1(nelem,natom_array,&
      i_chemi_eigen,&
      ntime,te_array,rho_array,time_array,&
      conce_ini,conce_nei,conce_eqi,&
      switch_eqi)
      use mod_eigen_matrix
      implicit double precision(a-h,o-z)
      integer nelem,natom_array(nelem),i_chemi_eigen(nelem)
      integer ntime
      real*8 te_array(ntime),&
             rho_array(ntime),&
             time_array(ntime)
      logical switch_eqi
      real*8 conce_ini(30,30),&
             conce_nei(30,30),conce_eqi(30,30)
!     -----------------------------------------------------------------
!           (1) equilibrium ionization at initial and final time
!     -----------------------------------------------------------------
      do 10 jelem = 1,nelem
      natom = natom_array(jelem)

!           initial time:
      if(switch_eqi)then
      call func_equilibrium_eigen(i_chemi_eigen(jelem),&
      natom,te_array(1),conce_ini(1:natom+1,natom))
      endif

!           final time:
      call func_equilibrium_eigen(i_chemi_eigen(jelem),&
      natom,te_array(ntime),conce_eqi(1:natom+1,natom))

   10 continue

!     -----------------------------------------------------------------
!           (2) mode 1: Origin time-steps
!     -----------------------------------------------------------------
      do 20 itime = 2,ntime
      dt = time_array(itime) - time_array(itime-1)
      rho_prev = rho_array(itime-1)
      te_prev  = te_array(itime-1)

!           solve the ionization equation using eigenvector. advance
!           one time-step from 't0' to 't1'.
      do 201 jelem = 1,nelem
      natom = natom_array(jelem)

      call func_solveionization_eigen(i_chemi_eigen(jelem),&
      natom,te_prev,rho_prev,&
      conce_ini(1:natom+1,natom),dt,conce_nei(1:natom+1,natom))

  201 continue
!           prepare for next step
      conce_ini = conce_nei

!           write to logfile
      write(51,*)itime,dt
   20 continue
!
      return
      end subroutine sub_solve_ioni_mode1

     
!     *****************************************************************
!     subroutine: sub_solve_ioni_mode2
!     *****************************************************************
!     Input:
!           nelem,
!           natom_array,
!           ntime,
!           te_array,
!           rho_array,
!           time_array,
!           conce_ini,
!           switch_eqi
!     Output:
!           conce_nei,
!           conce_eqi
!     Use:
!           mod_eigen_matrix
!           mod_input_parameter
!     Call:
!           func_equilibrium_eigen
!           func_solveionization_eigen
!     Update:
!           2013-09-14
!           2013-09-24
!           add mod_input_parameter,
!           add adaptive time-step.
!     Update: 
!           2013-10-22
!           remove use mode_input_parameter
!           transfer accuracy
!     -----------------------------------------------------------------
      subroutine sub_solve_ioni_mode2(nelem,natom_array,&
          i_chemi_eigen,&
          ntime,te_array,rho_array,time_array,&
          conce_ini,conce_nei,conce_eqi,&
          switch_eqi,&
          accuracy)
      use mod_eigen_matrix
      implicit double precision(a-h,o-z)
      integer nelem,ntime
      integer natom_array(nelem),i_chemi_eigen(nelem)
      real*8 te_array(ntime),&
             rho_array(ntime),&
             time_array(ntime)
      logical switch_eqi
      real*8 accuracy
      real*8 conce_ini(30,30),conce_out(30,30),&
             conce_nei(30,30),conce_eqi(30,30)

      real*8 conce_bak(30,30),conce_diff(30,30)

!     -----------------------------------------------------------------
!           (1) equilibrium ionization at initial and final time
!     -----------------------------------------------------------------
      do 10 jelem = 1,nelem
      natom = natom_array(jelem)

!           initial time:
      if(switch_eqi)then
      call func_equilibrium_eigen(i_chemi_eigen(jelem),&
      natom,te_array(1),conce_ini(1:natom+1,natom))
      endif

!           final time:
      call func_equilibrium_eigen(i_chemi_eigen(jelem),&
      natom,te_array(ntime),conce_eqi(1:natom+1,natom))

   10 continue

!     -----------------------------------------------------------------
!           (2) mode 2: Adaptive time-step
!     -----------------------------------------------------------------
!           2-2.1 enter main time-cycle: 'i30'
      time = time_array(1)
      te_prev = te_array(1)
      rho_prev = rho_array(1)
      do 30 itime = 2,ntime
          time_next = time_array(itime)
!             prepare for linear interpolation
          te_1 = te_array(itime-1)
          rho_1 = rho_array(itime-1)
          time_1 = time_array(itime-1)
          te_2 = te_array(itime)
          rho_2 = rho_array(itime)
          time_2 = time_next
          
!             set a initial time-step for inner cycles
          dt = time_next - time

!             inner cycle '31' performs between two neighbouring samples points 
          do 31 while (time .lt. time_next)
              diff_max = accuracy
              icyc = 0

              do 32 while (diff_max .ge. accuracy)
                  dt = dt/(2**icyc)
                  conce_diff = 0.0d0
                  te_next = func_interpol_linear(time+dt,time_1,time_2,te_1,te_2)
                  rho_next = func_interpol_linear(time+dt,time_1,time_2,rho_1,rho_2)

                  do 33 jelem = 1,nelem
                      natom = natom_array(jelem)
!                         forward, using Te_prev, ne_prev
                      call func_solveionization_eigen(i_chemi_eigen(jelem),&
                           natom,te_prev,rho_prev,&
                           conce_ini(1:natom+1,natom),dt,conce_nei(1:natom+1,natom))
!                         backward, using Te_next, ne_next
                      call func_solveionization_eigen(i_chemi_eigen(jelem),&
                           natom,te_next,rho_next,&
                           conce_ini(1:natom+1,natom),dt,conce_bak(1:natom+1,natom))
!                         difference between forward method and backward method
                      conce_diff(1:natom+1,natom) = dabs(conce_nei(1:natom+1,natom)&
                          -conce_bak(1:natom+1,natom))
               33 continue
!           find the maximum difference, and icyc + 1 for next
!           inner cycle
                  diff_max = maxval(conce_diff)
                  icyc = icyc + 1
           32 continue
!         now the time advance to time + dt
      time = time + dt

!           write logfile
      write(51,*)itime,time,dt

!           prepare for next cycle of '30'
      dt = time_next - time
      te_prev = te_next
      rho_prev = rho_next
      conce_ini = conce_nei
   31 continue
   30 continue
      
      return
      end subroutine sub_solve_ioni_mode2

!     *****************************************************************
!     subroutine: sub_solve_ioni_mode3
!     *****************************************************************
!     update:
!           2013-10-03
!     update:
!           2013-10-17
!           2013-10-22
!           transfer chang_perct   
      subroutine sub_solve_ioni_mode3(nelem,natom_array,&
          i_chemi_eigen,&
          ntime,te_array,rho_array,time_array,&
          conce_ini,conce_nei,conce_eqi,&
          switch_eqi,&
          chang_perct)
      use mod_eigen_matrix
      implicit double precision(a-h,o-z)
      integer nelem,ntime
      integer natom_array(nelem),i_chemi_eigen(nelem)
      real*8 te_array(ntime),&
             rho_array(ntime),&
             time_array(ntime)
      logical switch_eqi
      real*8 chang_perct
      real*8 conce_ini(30,30),&
             conce_nei(30,30),conce_eqi(30,30)

      real*8 dt_est_ne(ntime)
      real*8 eval_arr(30),dt_est_elem(nelem)

!     -----------------------------------------------------------------
!           (1) equilibrium ionization at initial and final time
!     -----------------------------------------------------------------
      do 10 jelem = 1,nelem
      natom = natom_array(jelem)

!           initial time:
      if(switch_eqi)then
      call func_equilibrium_eigen(i_chemi_eigen(jelem),&
      natom,te_array(1),conce_ini(1:natom+1,natom))
      endif

!           final time:
      call func_equilibrium_eigen(i_chemi_eigen(jelem),&
      natom,te_array(ntime),conce_eqi(1:natom+1,natom))

   10 continue

!     -----------------------------------------------------------------
!         (2) estimate time-step using eigenvalues and density
!     -----------------------------------------------------------------
!         check ne*t & dne
      do 20 i = 1,ntime
      index_te = int(anint((dlog10(te_array(i))-4.0d0)/dlogte)+1)
      rho  = rho_array(i)
      do 21 jelem = 1,nelem
      natom = natom_array(jelem)
      eval_arr(1:natom+1) = eigen(jelem)%evalues(1:natom+1,index_te)
      eval_max = maxval(dabs(eval_arr(1:natom+1)))
      dt_est_elem(jelem)  = chang_perct/(eval_max*rho)
   21 continue
!         set saft factor 0.5
      dt_est_ne(i) = 0.5d0*minval(dt_est_elem)
   20 continue

!     ----------------------------------------------------------------- 
!         (3) enter main cycle
!     -----------------------------------------------------------------
      i_loc = 1
      time = time_array(i_loc)
      do 40 while (i_loc .le. ntime-1)
          te_prev = te_array(i_loc)
!             for density: how many points can be skipped?
          i_plus = 1
          dt_est = dmin1(dt_est_ne(i_loc),dt_est_ne(i_loc + i_plus))
          dt     = time_array(i_loc + i_plus) - time
          do 41 while ((dt_est .ge. dt) .and. (i_loc + i_plus .le. ntime-1))
              i_plus = i_plus + 1
              dt = time_array(i_loc + i_plus) - time
              dt_est = dt_est_ne(i_loc + i_plus)
       41 continue
          i_plus_ne = i_plus

!             for temperature: how many points can be skipped?
          i_plus = 1
          st = anint(dlog10(te_prev)/dlogte)
          ed = anint(dlog10(te_array(i_loc+i_plus))/dlogte)
          do 42 while ((dabs(ed-st) .le. 0.) .and. (i_loc+iplus .le. ntime-1))
              i_plus = i_plus + 1
              ed = anint(dlog10(te_array(i_loc+i_plus))/dlogte)
       42 continue
          i_plus_te = i_plus

!             the number of skipped points
          i_plus = max(i_plus_ne,i_plus_te) - 1
          if(i_plus .le. 0)i_plus = 1
          if(i_loc + i_plus .ge. ntime)i_plus = 1
          dt = time_array(i_loc + i_plus) - time

!         -------------------------------------------------------------
!             case 1: i_plus > 1, skip several points
          if(i_plus .ge. 2)then
          rhot = 0.0d0
!             estimate avarage of density
          do i = i_loc,i_loc+i_plus-1
             rhot_temp = (time_array(i+1)-time_array(i))*0.5d0*(rho_array(i+1) + rho_array(i))
             rhot = rhot + rhot_temp 
          enddo
          rho_now = rhot/dt
          te_now  = te_prev
          do 50 jelem = 1,nelem
              natom = natom_array(jelem)
              call func_solveionization_eigen(i_chemi_eigen(jelem),&
              natom,te_now,rho_now,&
              conce_ini(1:natom+1,natom),dt,conce_nei(1:natom+1,natom))
       50 continue
          write(51,*)i_loc,time+dt,dt,i_plus
          conce_ini = conce_nei
          endif
!         -------------------------------------------------------------              
!             case 2: i_plus = 1, need interpolation
          if(i_plus .eq. 1)then
!!                 a. simple interpolation
              te_sta = dlog10(te_prev)
              te_end = dlog10(te_array(i_loc+1))
              dd_te  = te_end - te_sta
              n_inter = int(dabs(dd_te)/dlogte) + 2
              write(51,*)te_sta,te_end,n_inter
!                 enter cycle for interpolations
              dd_time = dt/dfloat(n_inter-1)
              do 60 i = 1,n_inter - 1
              time_inner = dfloat(i-1)*dd_time

              te_now = func_interpol_linear(time_inner,0.0d0,dt,&
                  te_array(i_loc),te_array(i_loc+1))
              rho_now=func_interpol_linear(time_inner+0.5d0*dd_time,0.0d0,dt,&
                  rho_array(i_loc),rho_array(i_loc+1))

!             perform nei-ionization calculations according to the 
!             ionization rate at te_now and rho_now
              do 61 jelem = 1,nelem
              natom = natom_array(jelem)
              call func_solveionization_eigen(i_chemi_eigen(jelem),&
              natom,te_now,rho_now,&
              conce_ini(1:natom+1,natom),dd_time,conce_nei(1:natom+1,natom))
           61 continue
              write(51,*)i_loc,time+time_inner,dd_time,'inter',n_inter
              conce_ini = conce_nei

           60 continue
          endif

!             time add 'one/several step' when cycle 50/60 has been 
!             finished. It does not depend on the number of 'n_inter'.
          i_loc = i_loc + i_plus
          time = time_array(i_loc)
   40 continue

      return
      end subroutine sub_solve_ioni_mode3

!     *****************************************************************
!     subroutine: func_equilibrium_eigen
!     *****************************************************************
!     Input: ichemi, natom,te
!     Output: ion charge state in equilibrium case
!     Use: mod_eigen_matrix
!     Write: 2013-09-14
      subroutine func_equilibrium_eigen(ichemi,natom,te,ft)
      use mod_eigen_matrix
      implicit double precision(a-h,o-z)
      integer natom,ichemi
      real*8 te,ft(natom+1)
!     -----------------------------------------------------------------
!           Get eigen matrix according to the inputted te
!     -----------------------------------------------------------------
      index_te = int(anint((dlog10(te)-4.0d0)*(1.0/dlogte))+1)

      ft(1:natom+1) = eigen(ichemi)%eqis(1:natom+1,index_te)
      return
      end subroutine func_equilibrium_eigen

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
      index_te = int(anint((dlog10(te)-4.0d0)*(1.0/dlogte))+1)
      
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

!     *****************************************************************
!     function_func_interpol_atk
!     *****************************************************************
      function func_interpol_atk(x,xarr,yarr,n)
      implicit double precision(a-h,o-z)
      integer n
      real*8 x,xarr(n),yarr(n)
      real*8 z(5),f(5)
      call atk(xarr,yarr,x,n,5,z,f,y_out)
      func_interpol_atk = y_out
      return
      end
      
      subroutine atk(x,y,t,n,m,z,f,a)
      implicit double precision (a-h,o-z)
      dimension f(m),x(n),y(n),z(m)
      if(m.gt.n) m=n
      do 1 i=1,n
      if(t.le.x(i)) go to 2
1     continue
      i=n
2     if(t.ne.x(i)) go to 3
      f(m)=y(i)
      go to 9
3     if(mod(m,2).eq.0) go to 4
      if(i.eq.1) go to 4
      if(t-x(i-1).ge.x(i)-t) go to 4
      i=i-1
4     i=i-m/2
      if(i.gt.0) go to 5
      i=1
      go to 6
5     if(i+m.gt.n) i=n-m+1
6     do 7 j=1,m
      z(j)=t-x(i)
      f(j)=y(i)
7     i=i+1
      m1=m-1
      do 8 i=1,m1
      fi=f(i)
      zi=z(i)
      i1=i+1
      do 8 j=i1,m
8     f(j)=fi+zi*(f(j)-fi)/(zi-z(j))
9     a=f(m)
      return
      end

!     *****************************************************************
!     function_func_interpol_linear
!     *****************************************************************
      function func_interpol_linear(x,x1,x2,y1,y2)
      implicit double precision(a-h,o-z)
      real*8 x,x1,x2,y1,y2
      func_interpol_linear = y1+(x-x1)*(y2-y1)/(x2-x1)
      return
      end

!     *****************************************************************
!     subroutine: sub_solve_ioni_develop
!     *****************************************************************
!     Input:
!           nelem,
!           natom_array,
!           ntime,
!           te_array,
!           rho_array,
!           time_array,
!           conce_ini,
!           switch_eqi
!     Output:
!           conce_nei,
!           conce_eqi
!     Use:
!           mod_eigen_matrix
!           mod_input_parameter
!     Call:
!           func_equilibrium_eigen
!           func_solveionization_eigen
!     Update:
!           2013-09-14
!           2013-09-24
!           add mod_input_parameter,
!           add adaptive time-step.
!     Update: 
!           2013-10-22
!           remove use mode_input_parameter
!     -----------------------------------------------------------------
      subroutine sub_solve_ioni_develop(nelem,natom_array,&
          i_chemi_eigen,&
          ntime,te_array,rho_array,time_array,&
          conce_ini,conce_nei,conce_eqi,&
          switch_eqi)
      use mod_eigen_matrix
      implicit double precision(a-h,o-z)
      integer nelem,ntime
      integer natom_array(nelem),i_chemi_eigen(nelem)
      real*8 te_array(ntime),&
             rho_array(ntime),&
             time_array(ntime)
      logical switch_eqi
      real*8 conce_ini(30,30),conce_out(30,30),&
             conce_nei(30,30),conce_eqi(30,30)

      real*8 conce_bak(30,30),conce_diff(30,30)
      real*8 te_deriv(ntime),rho_deriv(ntime)

      real*8,allocatable:: te_arr(:),rho_arr(:),time_arr(:)

      real*8 dt_est_rhot(ntime),&
             dt_est_rho(ntime),&
             dt_est_te(ntime),&
             dt_est_array(ntime)
      real*8,allocatable:: dt_rhotime(:),dt_drho(:)
      real*8 eval_arr(30),dt_est_elem(nelem),dt_est_elem_2(nelem)
      real*8,allocatable:: te_interpol(:),&
                           rho_interpol(:),&
                           time_interpol(:)
      real*8,allocatable:: te_arr_te(:),&
                           rho_arr_te(:),&
                           time_arr_te(:),&
                           te_arr_rho(:),&
                           rho_arr_rho(:),&
                           time_arr_rho(:)
!     -----------------------------------------------------------------
!           (1) equilibrium ionization at initial and final time
!     -----------------------------------------------------------------
      do 10 jelem = 1,nelem
      natom = natom_array(jelem)

!           initial time:
      if(switch_eqi)then
      call func_equilibrium_eigen(i_chemi_eigen(jelem),&
      natom,te_array(1),conce_ini(1:natom+1,natom))
      endif

!           final time:
      call func_equilibrium_eigen(i_chemi_eigen(jelem),&
      natom,te_array(ntime),conce_eqi(1:natom+1,natom))

   10 continue

!     ----------------------------------------------------------------- 
!           (2) check dte
!     -----------------------------------------------------------------
      do itime = 1,ntime-1
      dt_ori = time_array(itime+1)-time_array(itime)
      d_te = dabs(dlog10(te_array(itime+1)) - dlog10(te_array(itime)))
!      dt_est_te(itime) = dt_ori/(d_te/1.0d-03 + 1.0d0)
      dt_est_te(itime) = dt_ori*0.001d0/d_te
      enddo
      dt_est_te(ntime) = dt_est_te(ntime-1)

!           interpolation Te
      n_maxstep = 2*int((time_array(ntime)-time_array(1))/minval(dt_est_te))
      allocate(te_interpol(n_maxstep),&
               rho_interpol(n_maxstep),&
               time_interpol(n_maxstep))

      te_interpol(1) = te_array(1)
      rho_interpol(1) = rho_array(1)
      time_interpol(1)= time_array(1)

      i_loc = 1
      itime = 1
      time = time_array(1)
      do 40 while (i_loc .le. ntime-1)
          dt_ori = time_array(i_loc+1) - time
          dt_est = dmin1(dt_est_te(i_loc),dt_est_te(i_loc+1))
          ii = i_loc
          if(dt_est .ge. dt_ori)then
              do while ((dt_est .ge. dt_ori) .and. (ii .le. ntime-1))
                  ii = ii + 1
                  dt_ori = time_array(ii) - time
                  dt_est   = dt_est_te(ii)
              enddo
              i_loc = ii
              itime = itime + 1
              time = time_array(i_loc)
              time_interpol(itime) = time_array(i_loc)
              te_interpol(itime) = te_array(i_loc)
              rho_interpol(itime) = rho_array(i_loc)
          else
              n_inter = int(dt_ori/dt_est + 2.0d0)
              dt = dt_ori/dfloat(n_inter-1)
              do i = 1,n_inter-1
                  time_now = time + dt*dfloat(i)
                  itime = itime + 1
                  time_interpol(itime) = time_now
                  te_interpol(itime) = func_interpol_atk(time_now,time_array,te_array,ntime)
                  rho_interpol(itime) = func_interpol_atk(time_now,time_array,rho_array,ntime)
              enddo
              i_loc = i_loc+1
              time = time_now
          endif
   40 continue
!
      n_step = itime
!           define new history for interpolated Te
      allocate(te_arr_te(n_step),&
               rho_arr_te(n_step),&
               time_arr_te(n_step),&
               dt_rhotime(n_step),dt_drho(n_step))
      te_arr_te(1:n_step) = te_interpol(1:n_step)
      rho_arr_te(1:n_step) = rho_interpol(1:n_step)
      time_arr_te(1:n_step) = time_interpol(1:n_step)
      deallocate(te_interpol,rho_interpol,time_interpol)

!           check ne*t & dne
      do i = 1,n_step-1
      index_te = int(anint((dlog10(te_arr_te(i))-4.0d0)*(1.0/dlogte))+1)
      rho  = rho_arr_te(i)
      drho = dabs(rho_arr_te(i+1) - rho_arr_te(i))
      do jelem = 1,nelem
      natom = natom_array(jelem)
      eval_arr(1:natom+1) = eigen(jelem)%evalues(1:natom+1,index_te)
      eval_max = maxval(dabs(eval_arr(1:natom+1)))
      dt_est_elem(jelem)  = chang_perct/(eval_max*rho)
      dt_est_elem_2(jelem) = chang_perct/(eval_max*drho)
      enddo
      dt_rhotime(i) = minval(dt_est_elem)
      dt_drho(i)    = minval(dt_est_elem_2)
      enddo
      dt_rhotime(n_step) = dt_rhotime(n_step-1)
      dt_drho(n_step)    = dt_drho(n_step-1)

!           interpolation using ne*t and ne
      m_maxstep = 2*int((time_array(ntime)-time_array(1))/minval(dt_rhotime))
      allocate(te_arr_rho(m_maxstep),&
               rho_arr_rho(m_maxstep),&
               time_arr_rho(m_maxstep))

      te_arr_rho(1) = te_arr_te(1)
      rho_arr_rho(1) = rho_arr_te(1)
      time_arr_rho(1)= time_arr_te(1)

      i_loc = 1
      itime = 1
      time = time_arr_te(1)
      do 41 while (i_loc .le. n_step-1)
          dt_ori = time_arr_te(i_loc+1) - time
          dt_est = dmin1(dt_rhotime(i_loc),dt_rhotime(i_loc+1))
          ii = i_loc
          if(dt_est .ge. dt_ori)then
              do while ((dt_est .ge. dt_ori) .and. (ii .le. n_step-1))
                  ii = ii + 1
                  dt_ori = time_arr_te(ii) - time
                  dt_est = dt_rhotime(ii)
              enddo
              i_loc = ii
              itime = itime + 1
              time = time_arr_te(i_loc)
              time_arr_rho(itime) = time_arr_te(i_loc)
              te_arr_rho(itime) = te_arr_te(i_loc)
              rho_arr_rho(itime) = rho_arr_te(i_loc)
          else
              dt_est_2 = dt_drho(i_loc)
              n_inter = int(dt_ori/dt_est_2 + 2.0d0)
              dt = dt_ori/dfloat(n_inter-1)
              do i = 1,n_inter-1
                  time_now = time + dt*dfloat(i)
                  itime = itime + 1
                  time_arr_rho(itime) = time_now
                  te_arr_rho(itime) = func_interpol_atk(time_now,time_array,te_array,ntime)
                  rho_arr_rho(itime) = func_interpol_atk(time_now,time_array,rho_array,ntime)
              enddo
              i_loc = i_loc+1
              time = time_now
          endif
   41 continue
      m_step = itime

      deallocate(te_arr_te,rho_arr_te,time_arr_te,dt_rhotime,dt_drho)

!           ion fractions advance with time
      do 50 itime = 2,m_step
      te = te_arr_rho(itime-1)
      rho = rho_arr_rho(itime-1)
      dt = time_arr_rho(itime) - time_arr_rho(itime-1)

      do 51 jelem = 1,nelem
      natom = natom_array(jelem)
      call func_solveionization_eigen(i_chemi_eigen(jelem),&
      natom,te,rho,&
      conce_ini(1:natom+1,natom),dt,conce_nei(1:natom+1,natom))
   51 continue

      conce_ini = conce_nei
      write(51,*)itime,dt
   50 continue
      deallocate(te_arr_rho,rho_arr_rho,time_arr_rho)

      return
      end subroutine sub_solve_ioni_develop

