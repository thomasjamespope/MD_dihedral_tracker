!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      module mytypes
      type coordinate
       character(3)     :: l
       double precision :: x,y,z
      contains
       procedure :: coor_minus
       procedure :: coor_cross
       procedure :: norm => coor_norm
       generic   :: operator(-) => coor_minus
       generic   :: operator(*) => coor_cross
      endtype coordinate
      contains
!---------------------------------------------------------------------!
      function coor_minus(a,b) result(c)
      implicit none
      class(coordinate), intent(in) :: a,b
      type(coordinate)              :: c
      c%x = a%x - b%x
      c%y = a%y - b%y
      c%z = a%z - b%z
      c%l = "VEC"
      endfunction coor_minus
!---------------------------------------------------------------------!
      function coor_cross(a,b) result(c)
      implicit none
      class(coordinate), intent(in) :: a,b
      type(coordinate)              :: c
      c%x = a%y * b%z - b%y * a%z
      c%y = a%z * b%x - b%z * a%x
      c%z = a%x * b%y - b%x * a%y
      c%l = "VEC"
      endfunction coor_cross
!---------------------------------------------------------------------!
      function coor_norm(a) result(b)
      implicit none
      class(coordinate), intent(in) :: a
      type(coordinate)              :: b
      double precision              :: n
      n = 1.d0 / sqrt(a%x**2 + a%y**2 + a%z**2) 
      b%x = n * a%x
      b%y = n * a%y
      b%z = n * a%z
      b%l = a%l
      endfunction coor_norm
      endmodule mytypes
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      module functions
      use mytypes
      double precision, parameter :: deg = 180.d0 / (4.d0*atan(1.d0))
      contains
      function read_atoms_ijkl(u,i,j,k,l,n) result(c)
      integer , intent(in)           :: i,j,k,l
      type(coordinate), dimension(4) :: c
      integer                        :: u,ii,n
      do ii=1,n
       if(ii.eq.i)     then; read(u,*) c(1)
       elseif(ii.eq.j) then; read(u,*) c(2)
       elseif(ii.eq.k) then; read(u,*) c(3)
       elseif(ii.eq.l) then; read(u,*) c(4)
       else; read(u,*)
       endif
      enddo
      endfunction read_atoms_ijkl
!---------------------------------------------------------------------!
      function dihedral(c) result(d)
      implicit none
      type(coordinate), dimension(4) :: c
      type(coordinate)               :: b_ij, b_jk, b_kl, n_a, n_b, n_c
      double precision               :: x,y,d
      b_ij = c(2) - c(1)
      b_jk = c(3) - c(2)
      b_kl = c(4) - c(3)
      n_a  = b_ij * b_jk; n_a = n_a%norm()
      n_b  = b_jk * b_kl; n_b = n_b%norm()
      n_c  = n_a * b_jk%norm()
      x    = n_a%x * n_b%x + n_a%y * n_b%y + n_a%z * n_b%z
      y    = n_c%x * n_b%x + n_c%y * n_b%y + n_c%z * n_b%z
      d    = atan2(y,x) * deg
      endfunction dihedral
      endmodule functions
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
!---------------------------------------------------------------------!
      program MD_dihedral
      use functions
      implicit none
      integer, parameter             :: u = 60
      integer                        :: i,j,k,l
      type(coordinate), dimension(4) :: c
      integer                        :: ii,io,n
      character(100)                 :: file_name
      if (iargc().ne.5) stop "Incorrect Number of Arguments"
      call getarg(1,file_name); read(file_name,*,iostat=io) i    
      call getarg(2,file_name); read(file_name,*,iostat=io) j
      call getarg(3,file_name); read(file_name,*,iostat=io) k
      call getarg(4,file_name); read(file_name,*,iostat=io) l
      call getarg(5,file_name); open(u,file=trim(adjustl(file_name)))
      read(u,*,iostat=io) n
      read(u,*,iostat=io)
      do while (io.eq.0)
       write(*,*) dihedral(read_atoms_ijkl(u,i,j,k,l,n))
       read(u,*,iostat=io) n
       read(u,*,iostat=io)
      enddo
      endprogram MD_dihedral
!---------------------------------------------------------------------!

