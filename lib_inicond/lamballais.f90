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



subroutine lamballais_error(u)
  use share_vars
  use fieldexport
  implicit none
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent (in) :: u
  real(kind=pr) :: R,R1,R2,R3, e1,e2, diffu, uabs
  integer :: ix, iy

  R1=0.50d0
  R2=1.75d0
  R3=2.25d0
  x0=xl/2.d0
  y0=yl/2.d0
  e1=0.d0
  e2=0.d0
  
  do ix=0,nx-1
    do iy=0,ny-1
      R = dsqrt( (dble(ix)*dx-x0)**2 + (dble(iy)*dy-y0)**2 )
      if ((R>=R1).and.(R<=R2-dx)) then
        diffu = dsqrt(  u(ix,iy,1)**2+  u(ix,iy,2)**2) - &
                dsqrt(uex(ix,iy,1)**2+uex(ix,iy,2)**2)
        uabs  = dsqrt(uex(ix,iy,1)**2+uex(ix,iy,2)**2)
     
        e1 = e1 + diffu**2
        e2 = e2 + uabs**2
      endif
    enddo
  enddo
  
  e1 = dsqrt(e1) / dsqrt(e2)
  
  write(*,*) "rel error=", e1
  
  
end subroutine lamballais_error