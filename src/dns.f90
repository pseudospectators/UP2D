program dns
  use share_vars
  implicit none
  character(len=strlen) :: infile

  ! get the first command line argument
  call get_command_argument(1,infile)

  ! check what to do
  if ( index(infile,'.ini') .ne. 0) then
      !-------------------------------------------------------------------------
      ! the file is an *.ini file -> we run a normal simulation
      !-------------------------------------------------------------------------
      ! read parameters from ini file
      call get_params(infile)
      ! and run the simulation
      call Start_Simulation()

    elseif ((infile=="--postprocess").or.(infile=="-p").or.(infile=="-h").or.(infile=="--help")) then
      !-------------------------------------------------------------------------
      ! the first argument tells us that we're postprocessing
      !-------------------------------------------------------------------------
      call postprocessing()

    elseif ( infile == "--dry-run" ) then
      !-------------------------------------------------------------------------
      ! dry run that only generates and dumps the mask function, without
      ! allocating or computing the fluid.
      !-------------------------------------------------------------------------
      ! call dry_run()
  endif

end program dns



subroutine Start_Simulation()
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
end subroutine Start_Simulation
