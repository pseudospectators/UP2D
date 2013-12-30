program dns
  use share_vars
  use FieldExport
  use parameters
  implicit none
  
  call get_params("PARAMS.ini")
  dir_name = '.'
  simulation_name='test'
  call StartSimulation()   

end program dns



subroutine StartSimulation()
  use share_vars
  implicit none
  write (*,*) "*** information: entering StartSimulation"

  call system('mkdir '//trim(dir_name) )
  call system('mkdir '//trim(dir_name)//'/fields' )
  call system('mkdir '//trim(dir_name)//'/vor' )
  write (*,*) "*** information: created subdirectories"

  allocate ( dealiase(0:nx-1,0:ny-1) )
  allocate ( mask(0:nx-1,0:ny-1) )    
  allocate ( us(0:nx-1,0:ny-1,1:2) )

  write (*,*) "*** information: allocated memory"

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

end subroutine StartSimulation
