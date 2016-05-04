!-------------------------------------------------------------------------------
! Set the time step based on the CFL condition and penalization
! stability contidion. The following limitations exist:
! 1 - CFL condition
! 2 - fixed time step dt_fixed, ignoring all other constraints, if set in params
! 3 - penalization restriction dt<eps
! 4 - maximum time step dt_max, if set in params
! 5 - dt is smaller than tsave and tintegral
!-------------------------------------------------------------------------------
function adjust_dt(time,it, u)
  use share_vars
  implicit none
  real(kind=pr), dimension(0:nx-1, 0:ny-1,1:2), intent (in) :: u
  real(kind=pr), dimension(0:nx-1, 0:ny-1) :: work2
  real(kind=pr) :: dt1, u_max, adjust_dt, t
  real(kind=pr), intent (in) :: time
  integer, intent (in) :: it
  integer :: iy

  if (dt_fixed>0.0d0) then
    !-- fix the time step no matter what. the result may be unstable.
    dt1 = dt_fixed
  else
    ! CFL condition for the fluid
    ! compute mag(u)**2
    !$omp parallel do private(iy)
    do iy=0,ny-1
      work2(:,iy) = u(:,iy,1)**2 + u(:,iy,2)**2
    enddo
    !$omp end parallel do

    ! find maximum velocity
    u_max = dsqrt(maxval(work2))

    if(u_max>=1.0d3) then
      write(*,*) "Umax is very big, surely this is an error, ", u_max
      stop
    endif

    !-- Impose the CFL condition.
    if (u_max >= 1.0d-8) then
      dt1 = cfl*min(dx,dy)/u_max
    else
      !-- u_max is very very small
      dt1 = 1.0d-3
    endif

    !-- impose max dt, if specified
    if (dt_max>0.d0) dt1 = min(dt1,dt_max)
    !-- Impose penalty stability condition: dt cannot be larger than eps
    dt1 = min(0.99d0*eps,dt1)

    !*************************************************************************
    ! we respect all necessary restrictions, the following ones are optional
    !*************************************************************************
    if (intelligent_dt=="yes") then
      ! intelligent dt means we make sure not to jump past multiples of tsave
      ! tend tintegral tslice.
      if (time>=tsave_first) then
        ! Don't jump past save-points: if the time-step is larger than
        ! the time interval between outputs, decrease the time-step.
        dt1 = min(dt1,0.98d0*tsave)
        t = dble(ceiling(time/tsave))*tsave
        if ((time+dt1>t).and.(abs(time-t)>=1.0d-6)) then
          dt1=t-time
        endif
      endif

      ! Don't jump past save-points: if the time-step is larger than
      ! the time interval between outputs, decrease the time-step.
      dt1 = min(dt1,0.98d0*tdrag)
      t = dble(ceiling(time/tdrag))*tdrag
      if ((time+dt1>t).and.(abs(time-t)>=1.0d-6)) then
        dt1=t-time
      endif

      t = tmax
      if ((time+dt1>t).and.(abs(time-t)>=1.0d-6)) then
        dt1 = tmax-time
      endif
    endif

    !-- save max and mean velocities
    open (14, file = 'u_max.t', status = 'unknown', access = 'append')
    write (14,'(4(es11.4,1x))') time, u_max, sum(u(:,:,1))*dx*dy, sum(u(:,:,2))*dx*dy
    close (14)
    !-- save time step
    open (14, file = 'dt.t', status = 'unknown', access = 'append')
    write (14,'(es11.4,1x,i6,1x,es11.4)') time, it, dt1
    close (14)
  endif


  adjust_dt = dt1
end function
