subroutine time_step (u, uk, nlk, pk, vort, mask, us, mask_sponge)
  use vars
  use hdf5_wrapper
  use RK2_module
  use timing_module
  implicit none
  real(kind=pr) :: time=0.0d0, dt1=0.0d0, max_divergence
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2), intent(inout) :: u, uk, us, nlk
  real(kind=pr),dimension(0:nx-1,0:ny-1), intent(inout) :: pk, vort, mask, mask_sponge
  real(kind=pr) :: T_lastdrag=0.0d0, T_lastsave=0.0d0, t1
  integer :: it=0

  character(len=5) :: timestring

  if (FD_2nd) write (*,*) "!!! ATTENTION; RUNNING IN REDUCED ACCURACY MODE"
  t1 = MPI_wtime()

  !-- Initialize fields or read values from a backup file
  call init_fields(time, u, uk, pk, vort, nlk, mask, us, mask_sponge)
  call sponge_mask(time, mask_sponge)

  write (*,'("Initialization done, looping now.")')
  write (*,'("time=",es12.4," Tmax=",es12.4," it=",i2," nt=",i9)') time, Tmax, it, nt

  !----------------------------------------------------------------
  ! loop over time steps
  !----------------------------------------------------------------
  do while ((time<Tmax) .and. (it<=nt))
    !----------------------------------------------------------------
    !-- Actual time step
    !----------------------------------------------------------------
    call RK2 (time, dt1, it, u, uk, pk, vort, nlk, mask, us, mask_sponge)
    ! Advance in time
    time = time + dt1
    it = it + 1

    ! compute lift/drag
    open (14, file = 'ekin.t', status = 'unknown', access = 'append')
    write (14,'(2(es15.8,1x))') time, sum( (u(:,:,1)**2 + u(:,:,2)**2)/2.d0 )*dx*dy
    close (14)

    if ( time_for_output( time, dt1, it, tsave, itsave, Tmax, 0.d0) ) then
      ! save output fields to disk
      call save_fields(time, it, u, uk, vort, mask, us, mask_sponge)
      T_lastsave=time
    endif

    ! remaining time
    if (modulo(it,5)==0) then
      call are_we_there_yet( time, t1, dt1, it )
    endif
  enddo

  write (*,*) "Loop done."
end subroutine time_step


!-------------------------------------------------------------------------------
! Output how much time remains in the simulation.
!-------------------------------------------------------------------------------
subroutine are_we_there_yet(time, wtime_tstart, dt1, it)
  use vars
  use timing_module
  implicit none

  real(kind=pr),intent(inout) :: time,wtime_tstart,dt1
  integer,intent(inout) :: it
  real(kind=pr):: time_left, t2, time_nt

  ! elapsed time since time stepping started
  t2 = MPI_wtime() - wtime_tstart
  ! estimate remaining time until we reach tmax
  time_left = (tmax-time) * (t2/(time-tstart))
  ! estimate remaining time until we real nt time steps
  time_nt = 9e9*dble(nt-it) * (t2/dble(it))
  ! remaining time is minimum of both
  time_left = min( time_left, time_nt )

  write(*,'("time left: ",i2,"d ",i2,"h ",i2,"m ",i2,"s dt=",es10.2,"s t=",g10.2)') &
  floor(time_left/(24.d0*3600.d0))   ,&
  floor(mod(time_left,24.d0*3600.d0)/3600.d0),&
  floor(mod(time_left,3600.d0)/60.d0),&
  floor(mod(mod(time_left,3600.d0),60.d0)),&
  dt1,time
end subroutine are_we_there_yet
