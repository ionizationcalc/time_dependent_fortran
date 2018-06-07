!     *****************************************************************
!     subroutine: sub_read_eigen_matrix
!     Update:
!       2017-05-23
!       Read Chianti table 8_601
!       2017-10-05
!       Set the tetable_start in log value.
!       2017-10-06
!       Read temperature grid information: nte, dlogte, tetable_start
!       from eigentable.
!       2018-04-19
!       Bug fixed for reading Hydrogen table.
!       Read te_grid_table
!     *****************************************************************
      subroutine sub_read_eigen_matrix(path_eigen)
      use mod_eigen_matrix
      implicit double precision(a-h,o-z)
      character(len=150)datafile,path_eigen
!      integer:: nn_atom
      real*8:: tetable_end
!     -----------------------------------------------------------------
!
      do 10 ichemi = 1,n_element
!      natom = index_element(ichemi)
!
      datafile = trim(path_eigen)//trim(char_element(ichemi))//&
      'eigen.dat'
      open(11,file=trim(datafile),form='unformatted',action='read')
      read(11)nte, natom
!       Define array according to natom and nte
      if (.not. allocated(te_grid_table))allocate(te_grid_table(nte))
      allocate(eigen(ichemi)%eqis(natom+1,nte),&
               eigen(ichemi)%evalues(natom+1,nte),&
               eigen(ichemi)%evector(natom+1,natom+1,nte),&
               eigen(ichemi)%evector_invers(natom+1,natom+1,nte),&
               eigen(ichemi)%c(natom+1,nte),&
               eigen(ichemi)%r(natom+1,nte))
      read(11)te_grid_table
      read(11)eigen(ichemi)%eqis
      read(11)eigen(ichemi)%evalues
      read(11)eigen(ichemi)%evector
      read(11)eigen(ichemi)%evector_invers
      read(11)eigen(ichemi)%c
      read(11)eigen(ichemi)%r
      close(11)

   10 continue

!     Get table infors
      tetable_start = dlog10(te_grid_table(1))
      tetable_end   = dlog10(te_grid_table(nte))
      dlogte = (tetable_end-tetable_start)/dfloat(nte-1)
!
      return
      end subroutine sub_read_eigen_matrix
