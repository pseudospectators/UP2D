! ==========================================================================================================================
function timestep(time, u)
  use share_vars
  implicit none
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: u
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: work2
  real (kind=pr) :: dt1, u_max,timestep
  real (kind=pr), intent (in) :: time
  integer :: iy
  
  !========================================================================================
  !=              Calculate time step                                                     =
  !========================================================================================
  !-- CFL condition for the fluid
  
  !$omp parallel do private(iy)
  do iy=0,ny-1
    work2(:,iy) = u(:,iy,1)**2 + u(:,iy,2)**2
  enddo
  !$omp end parallel do
  
  u_max = sqrt(maxval(work2))  
  dt1 = cfl*min(xl/real(nx),yl/real(ny))/u_max

  !-- u_max is very small
  if (u_max <= 1.0e-10) then 
    dt1 = 1.e-3 
  endif

  !-- Time stepping control for volume penalization
  if ( dt1 >= 0.9*eps ) then
    dt1 = min(0.9*eps,dt1) ! time step is smaller than eps (Dmitry, 26 feb 08)
  endif
  
  !-- don't jump past final time
  if ((Tmax - time) < dt1) then
    dt1 = Tmax - time
  endif  

  ! save max and mean velocities
  open (14, file = trim(simulation_name)//'u_max', status = 'unknown', access = 'append')
  write (14,'(4(es11.4,1x))') time, u_max, sum(u(:,:,1))*dx*dy, sum(u(:,:,2))*dx*dy
  close (14)
  
  ! save time step
  open (14, file = trim(simulation_name)//'dt', status = 'unknown', access = 'append')
  write (14,'(2(es11.4,1x))') time, dt1
  close (14)
  
  timestep = dt1
end function  
