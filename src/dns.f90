program dns
  use share_vars
  use FieldExport
  use parameters
  implicit none
  character(len=40) :: infile

  call get_command_argument(1,infile)
  call get_params(infile)
  call StartSimulation()

end program dns



subroutine StartSimulation()
  use share_vars
  implicit none
  write (*,*) "*** information: entering StartSimulation"

  allocate ( dealiase(0:nx-1,0:ny-1) )
  allocate ( mask(0:nx-1,0:ny-1) )
  allocate ( phi(0:nx-1,0:ny-1) )
  allocate ( normals(0:nx-1,0:ny-1,1:2) )
  allocate ( us(0:nx-1,0:ny-1,1:2) )
  allocate ( uex(0:nx-1,0:ny-1,1:2) )
  allocate ( u_BC(0:nx-1,0:ny-1,1:2) )

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
  deallocate (dealiase,us, mask, uex, phi, normals, u_BC)
  call fft_free

end subroutine StartSimulation
