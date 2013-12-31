subroutine lamballais(u, uk, p, vor, nlk)
  use share_vars
  use FieldExport
  implicit none
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent (inout) :: u, uk, nlk
  real(kind=pr), dimension(0:nx-1,0:ny-1), intent (inout) :: vor, p
  real(kind=pr), dimension(0:nx-1,0:ny-1) :: vortk
  real(kind=pr) :: r0,we,d,r1,r2
  integer :: ix,iy
  
  
end subroutine lamballais