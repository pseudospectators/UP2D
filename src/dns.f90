program dns
  use vars
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
  use vars
  implicit none
  real(kind=pr), dimension(:,:,:), allocatable :: u, uk, nlk, us
  real(kind=pr), dimension(:,:), allocatable :: pk, vort,mask
  real(kind=pr), dimension(:,:), allocatable, save :: mask_sponge

  write (*,*) "*** information: entering StartSimulation"

  ! allocate memory
  allocate( dealiase(0:nx-1,0:ny-1) )
  allocate( mask(0:nx-1,0:ny-1) )
  allocate( us(0:nx-1,0:ny-1,1:2) )
  allocate( pk(0:nx-1,0:ny-1) )
  allocate( vort(0:nx-1,0:ny-1) )
  allocate( u(0:nx-1,0:ny-1,1:2) )
  allocate( uk(0:nx-1,0:ny-1,1:2) )
  allocate( nlk(0:nx-1,0:ny-1,1:2) )

  if (use_sponge==1) allocate ( mask_sponge(0:nx-1,0:ny-1) )
  ! for safety:
  mask = 0.d0
  us = 0.d0

  write (*,*) "*** information: allocated memory"

  if (FD_2nd) write (*,*) "!!! ATTENTION; RUNNING IN REDUCED ACCURACY MODE"

  ! Initialize fft
  call fft_initialize()
  write (*,*) "*** information: did fft_initialize"
  ! Set up mask for dealiasing
  dealiase = 1.0d0
  call dealiase_mask()

  ! Step forward in time
  write (*,*) "*** information: entering time_step"
  call time_step (u, uk, nlk, pk, vort, mask, us, mask_sponge)

  ! free memory
  deallocate (dealiase, us, mask,pk,vort,u,uk,nlk)
  ! free memory allocated by fft wrapper
  call fft_free

  if (allocated(mask_sponge)) deallocate ( mask_sponge )
end subroutine Start_Simulation
