module utils

  implicit none

contains
  
  function make_timestring(time) result(res)

  use vars
  implicit none
  
  character(:), allocatable  :: res
  real(kind=pr), intent (in) :: time

  character(len=strlen)      :: tmp
  write(tmp,'(i6.6)') nint(time*100.d0)

  res = trim(tmp)
  
end function make_timestring

end module utils

subroutine write_xdmf(time, it)
  use vars
  use utils
  implicit none

  real(kind=pr), intent (in) :: time
  integer,       intent (in) :: it

  character(LEN=80) :: filename
  character(LEN=80) :: filebase
  character(LEN=80) :: charbuf, char_nx_ny
  character(1), parameter :: endl  = char(10)  ! end of line
  
  integer           :: error

  integer, parameter    :: IO_NODES_PER_CELL=4
  character(len=80), parameter :: IO_TOPOLOGY_TYPE = 'Quadrilateral'
  
  filebase = 'UP2D_' // make_timestring(time)
  filename = 'UP2D_' // make_timestring(time) // '.xmf'

  open(10,file=filename,iostat=error)

  write(10,'(a)') '<?xml version="1.0" ?>'
  write(10,'(a)') '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">'
  write(10,'(a)') '  <Domain>'
  write(10,'(a)') '    <Grid Name="' // trim(filebase) // '" GridType="Uniform">'

  write(charbuf,fmt='(1(f7.2))',iostat=error) time
  write(10,'(a)') '      <Time TimeType="Single" Value="' // trim(charbuf) // '" />'

  ! Topology
  write(char_nx_ny,fmt='((i7,a,i7))',iostat=error) nx,' ',ny
  write(10,'(a)') '      <Topology TopologyType="2DCoRectMesh" NumberOfElements="' // trim(char_nx_ny) // '"  />'
  
  write(10,'(a)') '    <Geometry Type="ORIGIN_DXDY">'
  write(10,'(a)') '    <DataStructure'
  write(10,'(a)') '       Name="Origin"'
  write(10,'(a)') '       DataType="' // 'Double'  // '"'
  write(10,'(a)') '       Dimensions="2"'
  write(10,'(a)') '       Format="XML">'
  write(10,'(a)') '       0 0'
  write(10,'(a)') '    </DataStructure>'
  write(10,'(a)') '    <DataStructure'
  write(10,'(a)') '       Name="Spacing"'
  write(10,'(a)') '       DataType="' // 'Double' // '"'
  write(10,'(a)') '       Dimensions="2"'
  write(10,'(a)') '       Format="XML">'
  write(10,'(a)') '       1 1'
  write(10,'(a)') '    </DataStructure>'
  write(10,'(a)') '    </Geometry>' 

  ! write data attributes
  if ( iSaveVorticity == 1) then
     write(10,'(a)') '      <Attribute Center="Node" Name="' // 'vor' // '">'
     write(10,'(a)') '        <DataStructure'
     write(10,'(a)') '           DataType="' // 'Double' // '"'
     write(10,'(a)') '           Dimensions="' // trim(char_nx_ny) // '"'
     write(10,'(a)') '           Format="HDF">'
     write(10,'(a)') '           ' // 'vor_' // make_timestring(time) // '.h5' // ':/vor'
     write(10,'(a)') '        </DataStructure>'
     write(10,'(a)') '      </Attribute>'
  end if

  if ( iSaveVelocity == 1) then
     write(10,'(a)') '      <Attribute Center="Node" Name="' // 'ux' // '">'
     write(10,'(a)') '        <DataStructure'
     write(10,'(a)') '           DataType="' // 'Double' // '"'
     write(10,'(a)') '           Dimensions="' // trim(char_nx_ny) // '"'
     write(10,'(a)') '           Format="HDF">'
     write(10,'(a)') '           ' // 'ux_' // make_timestring(time) // '.h5' // ':/ux'
     write(10,'(a)') '        </DataStructure>'
     write(10,'(a)') '      </Attribute>'

     write(10,'(a)') '      <Attribute Center="Node" Name="' // 'uy' // '">'
     write(10,'(a)') '        <DataStructure'
     write(10,'(a)') '           DataType="' // 'Double' // '"'
     write(10,'(a)') '           Dimensions="' // trim(char_nx_ny) // '"'
     write(10,'(a)') '           Format="HDF">'
     write(10,'(a)') '           ' // 'uy_' // make_timestring(time) // '.h5' // ':/uy'
     write(10,'(a)') '        </DataStructure>'
     write(10,'(a)') '      </Attribute>'
  end if

  if ( iSaveMask == 1) then
     write(10,'(a)') '      <Attribute Center="Node" Name="' // 'mask' // '">'
     write(10,'(a)') '        <DataStructure'
     write(10,'(a)') '           DataType="' // 'Double' // '"'
     write(10,'(a)') '           Dimensions="' // trim(char_nx_ny) // '"'
     write(10,'(a)') '           Format="HDF">'
     write(10,'(a)') '           ' // 'mask_' // make_timestring(time) // '.h5' // ':/mask'
     write(10,'(a)') '        </DataStructure>'
     write(10,'(a)') '      </Attribute>'
  end if

  if ( iSavePressure == 1 ) then
     write(10,'(a)') '      <Attribute Center="Node" Name="' // 'p' // '">'
     write(10,'(a)') '        <DataStructure'
     write(10,'(a)') '           DataType="' // 'Double' // '"'
     write(10,'(a)') '           Dimensions="' // trim(char_nx_ny) // '"'
     write(10,'(a)') '           Format="HDF">'
     write(10,'(a)') '           ' // 'p_' // make_timestring(time) // '.h5' // ':/p'
     write(10,'(a)') '        </DataStructure>'
     write(10,'(a)') '      </Attribute>'
  end if
     
  ! write footer
  write(10,'(a)') '    </Grid>'
  write(10,'(a)') '  </Domain>'
  write(10,'(a)') '</Xdmf>'

  
  close(10)
  
end subroutine write_xdmf

subroutine save_fields(time, it, u, uk, vort, mask, us, mask_sponge)
  use vars
  use utils
  implicit none

  real(kind=pr), intent (in) :: time
  integer, intent(in) :: it
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent(inout) :: vort
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent(inout) :: mask, mask_sponge
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(inout) :: u,uk
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2),intent(inout) :: us

  real(kind=pr),dimension(:,:), allocatable :: work, pk
  character(:),allocatable :: timestring

  allocate(work(0:nx-1, 0:ny-1), pk(0:nx-1, 0:ny-1))

  !write(timestring,'(i6.6)') nint(time*100.d0)
  !write(timestring,'(i6.6)') it
  timestring=make_timestring(time)
  
  write(*,'("Saving. time=",es12.4," vormax=",es12.4," fname=",A)') time, maxval(vort), timestring

  if ( iSaveVorticity == 1) then
    call curl (uk, work)
    call ifft (work, vort)
    call SaveField (time, "vor_"//trim(timestring), vort)
  endif

  if ( iSaveVelocity == 1) then
    call SaveField (time, "ux_"//trim(timestring), u(:,:,1))
    call SaveField (time, "uy_"//trim(timestring), u(:,:,2))
  endif

  if ( iSaveMask == 1) then
    call SaveField (time, "mask_"//trim(timestring), mask)
  endif

  if ( iSavePressure == 1 ) then
    call cal_pressure (time, u, uk, pk, mask, us, mask_sponge)
    call ifft(pk, work)
    call SaveField( time, "p_"//trim(timestring), work)
  endif

  deallocate(work, pk)
end subroutine



! save a single field to file (this is just a wrapper)
subroutine SaveField( time, filename, field_out)
  use vars
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
