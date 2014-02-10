! ==========================================================================================================================
function timestep(time,it, u)
  use share_vars
  implicit none
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: u
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: work2
  real (kind=pr) :: dt1, u_max, timestep
  real (kind=pr), intent (in) :: time
  integer, intent (in) :: it
  integer :: iy
  
  !-- CFL condition for the fluid  
  !$omp parallel do private(iy)
  do iy=0,ny-1
    work2(:,iy) = u(:,iy,1)**2 + u(:,iy,2)**2
  enddo
  !$omp end parallel do
  
  u_max = dsqrt(maxval(work2))  
  dt1 = cfl*min(dx,dy)/u_max

  if (isnan(dt1)) then
    write(*,*) "Warning: dt_CFL = NaN"
    write(*,*) u_max, cfl, dt1
    stop
  endif
  
  !-- u_max is very small
  if (u_max <= 1.0d-10) then 
    dt1 = 1.0d-3 
  endif

  !-- Time stepping control for volume penalization
  if (( dt1 >= 0.9d0*eps ).and.(iMethod.ne."RK2_implicit")) then
    dt1 = min(0.9d0*eps,dt1)
  endif
  
  !-- Don't jump past final time
  if ((Tmax - time) < dt1) then
    dt1 = Tmax - time
  endif  
  
  !-- use fixed time step is set in params
  if (dt_fixed>0.d0) then
    dt1 = dt_fixed
  endif

  if (modulo(it,10)==0) then
    !-- save max and mean velocities
    open (14, file = trim(name)//'u_max', status = 'unknown', access = 'append')
    write (14,'(4(es11.4,1x))') time, u_max, sum(u(:,:,1))*dx*dy, sum(u(:,:,2))*dx*dy
    close (14)   
    !-- save time step
    open (14, file = trim(name)//'dt', status = 'unknown', access = 'append')
    write (14,'(es11.4,1x,i6,1x,es11.4)') time, it, dt1
    close (14)
  endif
  timestep = dt1
end function  


subroutine checknan(field)
use share_vars
implicit none
real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) :: field
integer :: ix,iy
do ix=0,nx-1
do iy=0,ny-1
if (isnan(field(ix,iy))) then
write (*,'("NaN! ix=",i4," iy=",i4)') ix,iy
stop
endif
enddo
enddo
end subroutine
