subroutine lamballais(u, uk, p, vor, nlk)
  use share_vars
  use FieldExport
  implicit none
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent (inout) :: u, uk, nlk
  real(kind=pr), dimension(0:nx-1,0:ny-1), intent (inout) :: vor, p
  real(kind=pr), dimension(0:nx-1,0:ny-1) :: vortk
  integer :: ix,iy
  character(len=10) :: string
  
  write(string,'(i3,"x",i3)') nx, ny
  
  !-- load exact solution
  call LoadField ( uex(:,:,1), 'us_'//trim(string)//'.ux')
  call LoadField ( uex(:,:,2), 'us_'//trim(string)//'.uy')

  write(*,*) 'Lamballais: loaded exact solution.'

  u = 0.d0
  uk = 0.d0
  vor = 0.d0
  
end subroutine lamballais