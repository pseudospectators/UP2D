!-------------------------------------------------------------------------------
! Wrapper for different postprocessing tools
!-------------------------------------------------------------------------------
subroutine postprocessing()
  use share_vars
  implicit none
  character(len=strlen) :: postprocessing_mode, mode, filename, key1, key2
  logical :: help
  help = .false.

  ! the first argument is either "-p" or "-h"
  call get_command_argument(1,mode)
  ! the second argument tells us what to do
  call get_command_argument(2,postprocessing_mode)

  !-----------------------------------------------------------------------------
  ! if mode is help, then we'll call the routines with help=.true. and skip
  ! all other output
  !-----------------------------------------------------------------------------
  if ((mode=="-h").or.(mode=="--help")) then
    ! we'll just show the help for the routine (we even skip the header)
    help=.true.
  else
    ! we'll actually do something, postprocessing
    write(*,'(80("~"))')
    write(*,'(A)') "~~~ UP2D is running in postprocessing mode"
    write(*,'(80("~"))')
    ! show the call from the command line in output
    call postprocessing_ascii_header( 6 )
  endif

  !-----------------------------------------------------------------------------
  ! check what to do
  !-----------------------------------------------------------------------------
  select case (postprocessing_mode)
  case ("--keyvalues")
    call get_command_argument(3,filename)
    call keyvalues(filename)
  case ("--compare-keys")
    call get_command_argument(3,key1)
    call get_command_argument(4,key2)
    call compare_key (key1,key2)
  case default
    write(*,*) "Available Postprocessing tools are:"
    write(*,*) "--keyvalues"
    write(*,*) "--compare-keys"
    write(*,*) "Postprocessing option is "// trim(adjustl(postprocessing_mode))
    write(*,*) "But I don't know what to do with that"
  end select
end subroutine postprocessing


!-------------------------------------------------------------------------------
! ./up2d --postprocess --keyvalues mask_00000.h5
!-------------------------------------------------------------------------------
! load the specified *.h5 file and creates a *.key file that contains
! min / max / mean / L2 norm of the field data. This is used for unit testing
! so that we don't need to store entire fields but rather the *.key only
subroutine keyvalues(filename)
  use hdf5_wrapper
  use share_vars
  implicit none
  character(len=*), intent(in) :: filename
  real(kind=pr) :: time, npoints, x,y
  real(kind=pr) :: maxl,minl,squarl,meanl,ql
  real(kind=pr), dimension(:,:), allocatable :: field
  integer :: ix,iy,iz

  call check_file_exists( filename )

  write (*,*) "analyzing file "//trim(adjustl(filename))//" for keyvalues"

  !---------------------------------------------------------
  ! in a first step, we fetch the attributes from the dataset
  ! namely the resolution is whats important
  ! this routine was created in the mpi2vis repo -> convert_hdf2xmf.f90
  !---------------------------------------------------------
  call Fetch_attributes_2d_openmp( filename, nx, ny, xl, yl, time , nu )
  write(*,*) nx,ny
  allocate( field(0:nx-1,0:ny-1) )

  call read_flusi_hdf5_2d_openmp (filename, field)

  npoints = dble(nx)*dble(ny)

  ! compute an additional quantity that depends also on the position
  ! (the others are translation invariant)
  ql = 0.d0
  do iy = 0, ny-1
    do ix = 0, nx-1
      x = dble(ix)*xl/dble(nx)
      y = dble(iy)*yl/dble(ny)
      ql = ql + x*field(ix,iy) + y*field(ix,iy)
    enddo
  enddo

  maxl = maxval(field)
  minl = minval(field)
  squarl = sum(field**2)
  meanl  = sum(field)

  ql = ql / npoints
  squarl = squarl / npoints
  meanl = meanl / npoints

  open  (14, file = filename(1:index(filename,'.'))//'key', status = 'replace')
  write (14,'(6(es15.8,1x))') time, maxl, minl, meanl, squarl, ql
  write (* ,'(6(es15.8,1x))') time, maxl, minl, meanl, squarl, ql
  close (14)

  deallocate (field)
end subroutine keyvalues



!-------------------------------------------------------------------------------
! ./flusi --postprocess --compare-keys mask_00000.key saved.key
!-------------------------------------------------------------------------------
! compares to *.key files if they're equal
subroutine compare_key(key1,key2)
  use share_vars
  implicit none
  character(len=*), intent(in) :: key1,key2
  real(kind=pr) :: a1,a2,b1,b2,c1,c2,d1,d2,t1,t2,q1,q2
  real(kind=pr) :: e1,e2,e3,e4,e0,e5

  open  (14, file = key1, status = 'unknown', action='read')
  read (14,'(6(es15.8,1x))') t1,a1,b1,c1,d1,q1
  close (14)

  open  (14, file = key2, status = 'unknown', action='read')
  read (14,'(6(es15.8,1x))') t2,a2,b2,c2,d2,q2
  close (14)

  write (*,'("present  : time=",es15.8," max=",es15.8," min=",es15.8," sum=",es15.8," sum**2=",es15.8," q=",es15.8)') &
  t1,a1,b1,c1,d1,q1

  write (*,'("reference: time=",es15.8," max=",es15.8," min=",es15.8," sum=",es15.8," sum**2=",es15.8," q=",es15.8)') &
  t2,a2,b2,c2,d2,q2

  ! errors:
  if (dabs(t2)>=1.0d-7) then
    e0 = dabs( (t2-t1) / t2 )
  else
    e0 = dabs( (t2-t1) )
  endif

  if (dabs(a2)>=1.0d-7) then
    e1 = dabs( (a2-a1) / a2 )
  else
    e1 = dabs( (a2-a1) )
  endif

  if (dabs(b2)>=1.0d-7) then
    e2 = dabs( (b2-b1) / b2 )
  else
    e2 = dabs( (b2-b1) )
  endif

  if (dabs(c2)>=1.0d-7) then
    e3 = dabs( (c2-c1) / c2 )
  else
    e3 = dabs( (c2-c1) )
  endif

  if (dabs(d2)>=1.0d-7) then
    e4 = dabs( (d2-d1) / d2 )
  else
    e4 = dabs( (d2-d1) )
  endif

  if (dabs(q2)>=1.0d-7) then
    e5 = dabs( (q2-q1) / q2 )
  else
    e5 = dabs( (q2-q1) )
  endif

  write (*,'("err(rel) : time=",es15.8," max=",es15.8," min=",es15.8," sum=",es15.8," sum**2=",es15.8," q=",es15.8)') &
  e0,e1,e2,e3,e4,e5

  if ((e1<1.d-4) .and. (e2<1.d-4) .and. (e3<1.d-4) .and. (e4<1.d-4) .and. (e0<1.d-4) .and. (e5<1.d-4)) then
    ! all cool
    write (*,*) "OKAY..."
    call exit(0)
  else
    ! very bad
    write (*,*) "ERROR"
    call exit(1)
  endif
end subroutine compare_key




!-------------------------------------------------------------------------------
! write the call to flusi (i.e. the command line arguments) to an ascii file
! use this in conjunction with small ascii output files to document a little
! bit what you have been doing.
!-------------------------------------------------------------------------------
subroutine postprocessing_ascii_header( io_stream )
  use share_vars
  implicit none
  integer, intent(in) :: io_stream
  integer :: i
  character(len=strlen) :: arg

  ! MATLAB comment character:
  write(io_stream,'(A,1x)',advance='no') "% CALL: ./UP2D "

  arg = "-p"
  i=1
  ! loop over command line arguments and dump them to the file:
  do while ( arg /= "" )
    call get_command_argument(i,arg)
    write(io_stream,'(A,1x)',advance='no') trim(adjustl(arg))
    i=i+1
  end do

  write(io_stream,'(A,1x)',advance='yes') "%"
end subroutine postprocessing_ascii_header
