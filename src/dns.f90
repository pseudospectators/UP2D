program dns
  use share_vars
  use FieldExport
  implicit none
  character(len=40) :: infile

  ! read parameters from ini file
  call get_command_argument(1,infile)
  call get_params(infile)
  ! and run the simulation
  call StartSimulation()

end program dns



subroutine StartSimulation()
  use share_vars
  implicit none
  write (*,*) "*** information: entering StartSimulation"

  allocate ( dealiase(0:nx-1,0:ny-1) )
  allocate ( mask(0:nx-1,0:ny-1) )
  allocate ( us(0:nx-1,0:ny-1,1:2) )

  if (use_sponge==1) allocate ( mask_sponge(0:nx-1,0:ny-1) )
  ! for safety:
  mask = 0.d0
  us = 0.d0

  write (*,*) "*** information: allocated memory"

  if (FD_2nd) write (*,*) "!!! ATTENTION; RUNNING IN REDUCED ACCURACY MODE"

! Initialize fft
  call fft_initialize
  write (*,*) "*** information: did fft_initialize"
! Set up mask for dealiasing
  dealiase = 1.0d0
  call dealiase_mask
! Step forward in time
  write (*,*) "*** information: entering time_step"
  call time_step  !after time_step, the last vort field is saved in mask_sponge
  deallocate (dealiase, us, mask)
  call fft_free

  if (use_sponge==1) deallocate ( mask_sponge )
end subroutine StartSimulation
