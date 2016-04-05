subroutine save_fields(time, it, u, uk, vort)
  use share_vars
  implicit none

  real(kind=pr), intent (in) :: time
  integer, intent(in) :: it
  real(kind=pr), dimension (0:nx-1, 0:ny-1, 1:2), intent (in) :: u,uk
  real(kind=pr), dimension(0:nx-1,0:ny-1), intent(in) :: vort

  real(kind=pr),dimension (:,:), allocatable :: work, pk
  character(len=strlen) :: timestring

  allocate(work(0:nx-1, 0:ny-1), pk(0:nx-1, 0:ny-1))

  write(timestring,'(i6.6)') nint(time*100.d0)
  write(timestring,'(i6.6)') it
  write(*,'("Saving. time=",es12.4," vormax=",es12.4," fname=",A)') time, maxval(vort), timestring


  if ( iSaveVorticity == 1) then
    call SaveField( time, "vor_"//trim(timestring), vort)
  endif

  if ( iSaveVelocity == 1) then
    call SaveField( time, "ux_"//trim(timestring), u(:,:,1))
    call SaveField( time, "uy_"//trim(timestring), u(:,:,2))
  endif

  if ( iSaveMask == 1) then
    call SaveField( time, "mask_"//trim(timestring), mask)
  endif

  if ( iSavePressure == 1 ) then
    call cal_pressure (time,u,uk,pk)
    call cofitxy(pk, work)
    call SaveField( time, "p_"//trim(timestring), work)
  endif

  deallocate(work, pk)
end subroutine



! save a single field to file (this is just a wrapper)
subroutine SaveField( time, filename, field_out)
  use share_vars
  use hdf5_wrapper
  implicit none

  ! The field to be written to disk:
  real(kind=pr),intent(in) :: field_out(0:nx-1,0:ny-1)
  real(kind=pr),intent(in) :: time
  character(len=*), intent (in) :: filename

  call write_flusi_hdf5_2d_openmp( time, filename, field_out)
end subroutine



! checks if a given file ("fname") exists. if not, code is stopped brutally
subroutine check_file_exists(fname)
  implicit none

  character (len=*), intent(in) :: fname
  logical :: exist1

  inquire ( file=fname, exist=exist1 )

  if ( exist1 .eqv. .false.) then
    write (*,'("ERROR! file: ",A," not found")') trim(adjustl(fname))
    stop
  endif

end subroutine check_file_exists



! overwrite and initialize file
subroutine init_empty_file( fname )
  implicit none
  character (len=*), intent(in) :: fname

  open (15, file=fname,status='replace')
  close(15)
end subroutine
