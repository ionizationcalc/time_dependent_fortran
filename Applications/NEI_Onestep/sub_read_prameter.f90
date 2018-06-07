!     *****************************************************************
!     subroutine: sub_read_prameter
!     *****************************************************************
      subroutine sub_read_prameter
      use mod_input_parameter
      implicit none
      character (len=150)charhead
      open(15,file="input.txt")
      read(15,fmt="(a150)")charhead
      read(15,*)nelem
      allocate(natom_list(nelem))
      read(15,*)natom_list
      write(*,fmt="(a14,30i02)")"Elements list:",natom_list
      read(15,fmt="(a150)")charhead
      read(15,*)path_eigen
      close(15)
      return
      end subroutine sub_read_prameter
