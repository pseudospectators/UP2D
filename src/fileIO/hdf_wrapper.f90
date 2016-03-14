module hdf5_wrapper

  use hdf5
  use share_vars, only : nx,ny,xl,yl,nu,pr,strlen,eps

contains

!-------------------------------------------------------------------------------
! this subroutine writes an array (and just one) to a HDF5 file according to FLUSI
! conventions. Thus, e.g., mask_0000000.h5 is an array at time t=000.000 and con-
! tains a datset "mask" with attributes "nxyz", "domain_size", "epsi", "time" and
! "viscosity"
!-------------------------------------------------------------------------------
subroutine write_flusi_hdf5_2d_openmp( time, filename, field_out)
  implicit none

  ! The field to be written to disk:
  real(kind=pr),intent(in) :: field_out(0:nx-1,0:ny-1)
  real(kind=pr),intent(in) :: time
  character(len=*), intent (in) :: filename

  integer, parameter :: rank = 2 ! data dimensionality (2D or 3D)
  integer(hid_t) :: file_id   ! file identifier
  integer(hid_t) :: dset_id   ! dataset identifier
  integer(hid_t) :: filespace ! dataspace identifier in file
  integer(hid_t) :: memspace  ! dataspace identifier in memory
  integer(hid_t) :: plist_id  ! property list identifier

  ! dataset dimensions in the file.
  integer(hsize_t), dimension(rank) :: dimensions_file
  ! hyperslab dimensions
  integer(hsize_t), dimension(rank) :: dimensions_local
  ! chunk dimensions
  integer(hsize_t), dimension(rank) :: chunking_dims
  ! how many blocks to select from dataspace
  integer(hsize_t),  dimension(rank) :: count  = 1
  integer(hssize_t), dimension(rank) :: offset
  ! stride is spacing between elements, this is one here. striding is done in the
  ! caller; here, we just write the entire (possibly downsampled) field to disk.
  integer(hsize_t),  dimension(rank) :: stride = 1
  integer :: error  ! error flags

  ! HDF attribute variables
  integer, parameter :: arank = 1
  integer(hsize_t), DIMENSION(1) :: adims  ! Attribute dimension

  integer :: mpierror, i, mindim, maxdim, nxred_file, nyred_file, nzred_file, mpicode

  ! ----------------------------------------------------------------------------
  ! Compute the dimension of the complete field (i.e. the union of all CPU's)
  ! which we will write to file.
  ! ----------------------------------------------------------------------------

  ! Tell HDF5 how our  data is organized:
  dimensions_file = (/nx, ny/)
  offset = 0
  dimensions_local = (/nx, ny/)
  chunking_dims = (/nx, ny/)

  !!! Set up the HDF data structures:

  ! Initialize HDF5 library and Fortran interfaces.
  call h5open_f(error)
  ! Setup file access property list with parallel I/O access.
  ! this sets up a property list (plist_id) with standard values for
  ! FILE_ACCESS
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! setup serial driver (no MPI)
  call H5Pset_fapl_stdio_f(plist_id, error)


  if ( index(filename,'.h5')==0 ) then
    ! file does not contain *.h5 ending -> add suffix
    call h5fcreate_f(trim(adjustl(filename))//'.h5', H5F_ACC_TRUNC_F, &
    file_id, error, access_prp = plist_id)
  else
    ! field DOES contain .h5 ending -> just write
    call h5fcreate_f(trim(adjustl(filename)), H5F_ACC_TRUNC_F, &
    file_id, error, access_prp = plist_id)
  endif

  ! this closes the property list plist_id (we'll re-use it)
  call h5pclose_f(plist_id, error)

  !-----------------------------------------------------------------------------
  ! create dataspace "filespace" to write to
  !-----------------------------------------------------------------------------
  ! Create the data space for the  dataset.
  ! Dataspace in the file: contains all data from all procs
  call h5screate_simple_f(rank, dimensions_file, filespace, error)

  ! Create chunked dataset.
  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call h5pset_chunk_f(plist_id, rank, chunking_dims, error)
  ! Output files in single precision
  call h5dcreate_f(file_id, get_dsetname(filename), H5T_NATIVE_REAL, filespace, &
       dset_id, error, plist_id)
  call h5sclose_f(filespace, error)

  ! Select hyperslab in the file.
  call h5dget_space_f(dset_id, filespace, error)
  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
  error, stride, dimensions_local)

  ! Create property list for collective dataset write
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)

  !-----------------------------------------------------------------------------
  ! create dataspace "memspace" to be written
  !-----------------------------------------------------------------------------
  ! dataspace in memory: contains only local data
  call h5screate_simple_f(rank, dimensions_local, memspace, error)

  !-----------------------------------------------------------------------------
  ! actual writing of heavy data
  !-----------------------------------------------------------------------------
  ! Write the dataset collectively, double precision in memory
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field_out, dimensions_file, &
  error, file_space_id = filespace, mem_space_id = memspace,&
  xfer_prp = plist_id)

  !!! Write the attributes to the HDF files.
  ! The attributes written are time, penalisation parameter,
  ! computational resolution, and physical domain size.
  adims = (/1/)
  call write_attribute_dble(adims,"time",(/time/),1,dset_id)
  call write_attribute_dble(adims,"viscosity",(/nu/),1,dset_id)
  call write_attribute_dble(adims,"epsi",(/eps/),1,dset_id)
  adims = (/2/)
  call write_attribute_dble(adims,"domain_size",(/xl,yl/),2,dset_id)
  call write_attribute_int(adims,"nxyz",(/nx,ny/),2,dset_id)

  !!! Close dataspaces:
  call h5sclose_f(filespace, error)
  call h5sclose_f(memspace, error)
  call h5dclose_f(dset_id, error) ! Close the dataset.
  call h5pclose_f(plist_id, error) ! Close the property list.
  call h5fclose_f(file_id, error) ! Close the file.
  call h5close_f(error) ! Close Fortran interfaces and HDF5 library.


end subroutine write_flusi_hdf5_2d_openmp



! Write a given attribute with attribute name aname and dimensions
! adims/dims to a given dataset identifier dset_id. Double version.
subroutine write_attribute_dble(adims,aname,attribute,dim,dset_id)
  implicit none

  integer, intent(in) :: dim
  integer(hsize_t), DIMENSION(dim), intent(in) :: adims  ! Attribute dimension
  real (kind=pr), DIMENSION(dim), intent (in) :: attribute
  character(len=*), intent(in) :: aname ! attribute name
  integer(hid_t),intent(in) :: dset_id  ! dataset identifier
  integer, parameter :: arank = 1

  integer :: error  ! error flags
  integer(hid_t) :: aspace_id ! Attribute Dataspace identifier
  integer(hid_t) :: attr_id   ! Attribute identifier

  ! Determine the dataspace identifier aspace_id
  call h5screate_simple_f(arank,adims,aspace_id,error)

  ! set attr_id, ie create an attribute attached to the object dset_id
  call h5acreate_f(dset_id,aname,H5T_NATIVE_DOUBLE,aspace_id,attr_id,error)

  ! Write the attribute data attribute to the attribute identifierattr_id
  call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attribute, adims, error)

  call h5aclose_f(attr_id, error) ! Close the attribute.
  call h5sclose_f(aspace_id, error) ! Terminate access to the data space.
end subroutine write_attribute_dble


! Write a given attribute with attribute name aname and dimensions
! adims/dims to a given dataset identifier dset_id. Integer version.
subroutine write_attribute_int(adims,aname,attribute,dim,dset_id)
  implicit none

  integer, intent(in) :: dim
  integer(hsize_t), DIMENSION(dim), intent(in) :: adims  ! Attribute dimension
  integer, DIMENSION(dim), intent (in) :: attribute
  character(len=*), intent(in) :: aname ! attribute name
  integer(hid_t),intent(in) :: dset_id  ! dataset identifier
  integer, parameter :: arank = 1

  integer :: error  ! error flags
  integer(hid_t) :: aspace_id ! Attribute Dataspace identifier
  integer(hid_t) :: attr_id   ! Attribute identifier

  ! Determine the dataspace identifier aspace_id
  call h5screate_simple_f(arank,adims,aspace_id,error)

  ! set attr_id, ie create an attribute attached to the object dset_id
  call h5acreate_f(dset_id,aname,H5T_NATIVE_INTEGER,aspace_id,attr_id,error)

  ! Write the attribute data attribute to the attribute identifier attr_id.
  call h5awrite_f(attr_id,H5T_NATIVE_INTEGER,attribute,adims,error)

  call h5aclose_f(attr_id,error) ! Close the attribute.
  call h5sclose_f(aspace_id,error) ! Terminate access to the data space.
end subroutine write_attribute_int


! Read in a single file that follows the naming convention
! note you need to know the dimensions and domain decomposition before
! calling it.
subroutine read_flusi_hdf5_2d_openmp( filename, field )
  implicit none

  character(len=*),intent(in) :: filename
  real(kind=pr),dimension(0:nx-1,0:ny-1),intent (out) :: field

  integer, parameter        :: rank = 2 ! data dimensionality (2D or 3D)
  real (kind=pr)            :: time, xl_file, yl_file
  real (kind=pr)            :: fmax,fmin,favg,viscosity_dummy
  character(len=80)         :: dsetname
  integer                   :: nx_file, ny_file, i

  integer(hid_t) :: file_id       ! file identifier
  integer(hid_t) :: dset_id       ! dataset identifier
  integer(hid_t) :: filespace     ! dataspace identifier in file
  integer(hid_t) :: memspace      ! dataspace identifier in memory
  integer(hid_t) :: plist_id      ! property list identifier

  ! dataset dimensions in the file.
  integer(hsize_t), dimension(rank) :: dimensions_file
  integer(hsize_t), dimension(rank) :: dimensions_local  ! chunks dimensions
  integer(hsize_t), dimension(rank) :: chunking_dims  ! chunks dimensions

  integer(hsize_t),  dimension(rank) :: count  = 1
  integer(hssize_t), dimension(rank) :: offset
  integer(hsize_t),  dimension(rank) :: stride = 1
  integer :: error, mpicode  ! error flags

  ! what follows is for the attribute "time"
  integer, parameter :: arank = 1

  ! the dataset is named the same way as the file: (this is convention)
  dsetname = get_dsetname(filename)

  write(*,'(40("~"))')
  write(*,'("Reading from file ",A)') trim(adjustl(filename))
  write(*,'("dsetname=",A)') trim(adjustl(dsetname))

  !-----------------------------------------------------------------------------
  ! perform tests
  !-----------------------------------------------------------------------------
  call check_file_exists ( filename )

  ! fetch attributes from file to see if it is a good idea to load it
  call Fetch_attributes_2d_openmp( filename,nx_file,ny_file,&
  xl_file,yl_file,time, viscosity_dummy )

  write(*,'("nx=",i4," ny=",i4," time=",g12.4," viscosity=",g16.4)')&
   nx_file,ny_file,time,viscosity_dummy
  write(*,'("xl=",g12.4," yl=",g12.4)') xl_file,yl_file

  ! if the domain size doesn't match, proceed, but yell.
  if ((xl.ne.xl_file).or.(yl.ne.yl_file)) then
    write (*,'(A)') " WARNING! Domain size mismatch."
    write (*,'("in memory:   xl=",es12.4,"yl=",es12.4)') xl,yl
    write (*,'("but in file: xl=",es12.4,"yl=",es12.4)') xl_file,yl_file
    write (*,'(A)') "proceed, with fingers crossed."
  endif

  !-----------------------------------------------------------------------------
  ! load the file
  !-----------------------------------------------------------------------------
  ! Initialize HDF5 library and Fortran interfaces.
  call h5open_f(error)

  ! Setup file access property list with SERIAL I/O access.  this
  ! sets up a property list ("plist_id") with standard values for
  ! FILE_ACCESS
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! setup stdio (not mpi)
  call H5Pset_fapl_stdio_f(plist_id, error)
  ! open the file in parallel
  call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error, plist_id)
  ! this closes the property list (we'll re-use it)
  call h5pclose_f(plist_id, error)

  ! Definition of memory distribution
  dimensions_file = (/nx,ny/)
  dimensions_local = (/nx,ny/)
  chunking_dims = (/nx,ny/)
  offset = 0

  ! if the resolutions do not match, yell and hang yourself
  if ((nx/=nx_file).or.(ny/=ny_file)) then
    write (*,'(A)') "ERROR! Resolution mismatch"
    write (*,'("in memory:   nx=",i4," ny=",i4)') nx,ny
    write (*,'("but in file: nx=",i4," ny=",i4)') nx_file,ny_file
    stop
  endif

  !----------------------------------------------------------------------------
  ! Read actual field from file (dataset)
  !----------------------------------------------------------------------------
  ! dataspace in the file: contains all data from all procs
  call h5screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  call h5screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset
  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call h5pset_chunk_f(plist_id, rank, chunking_dims, error)

  ! Open an existing dataset.
  call h5dopen_f(file_id, dsetname, dset_id, error)

  ! Select hyperslab in the file.
  call h5dget_space_f(dset_id, filespace, error)
  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
  error, stride, dimensions_local)

  ! Create property list for collective dataset read
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)

  call h5dread_f( dset_id, H5T_NATIVE_DOUBLE, field, dimensions_local, error, &
  mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

  call h5sclose_f(filespace, error)
  call h5sclose_f(memspace, error)
  call h5pclose_f(plist_id, error) ! note the dataset remains opened

  ! Close dataset
  call h5dclose_f(dset_id, error)
  call h5fclose_f(file_id,error)
  call H5close_f(error)

  fmax=maxval(field)
  fmin=minval(field)
  favg=sum(field)/dble(nx*ny)

  write (*,'("max=",g12.4," min=",g12.4," mean=",g12.4)') fmax,fmin,favg
  write (*,'("Done reading file")')
  write(*,'(40("~"))')

end subroutine read_flusi_hdf5_2d_openmp


!----------------------------------------------------
! This routine fetches the resolution, the domain size and the time
! form a *.h5 file
!----------------------------------------------------
! filename: a *.h5 file to read from.
! dsetname: a dataset inside the file. in our system, we only have one per file
!           and this matches the prefix: mask_00010.h5  --> dsetname = "mask"
! note:
!           the file must contain the dataset
!           but especially the attributes "nxyz", "time", "domain_size"
!----------------------------------------------------
subroutine Fetch_attributes_2d_openmp( filename, nx, ny, xl, yl, time, viscosity )
  implicit none

  integer, parameter :: pr = 8
  integer, intent (out) :: nx, ny
  real (kind=pr), intent(out) :: xl,yl,time,viscosity

  character(len=*) :: filename  ! file name
  character(len=80) :: dsetname  ! dataset name

  integer(hid_t) :: file_id       ! file identifier
  integer(hid_t) :: dset_id       ! dataset identifier
  integer(hid_t) :: attr_id       ! attribute identifier
  integer(hid_t) :: aspace_id     ! attribute dataspace identifier

  real (kind=pr) ::  attr_data  ! attribute data
  real (kind=pr), dimension (1:2) :: attr_data2
  integer, dimension (1:2) :: attr_data3

  integer     ::   error ! error flag
  integer(hsize_t), dimension(1) :: data_dims
  logical :: exists


  call check_file_exists ( filename )
  ! get dataset name from file name
  dsetname = get_dsetname( filename )

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)

  ! Open an existing file.
  CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
  ! Open an existing dataset.
  CALL h5dopen_f(file_id, dsetname, dset_id, error)


  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! open attribute (time)
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL h5aopen_f(dset_id, "time", attr_id, error)

  ! Get dataspace and read
  CALL h5aget_space_f(attr_id, aspace_id, error)
  data_dims(1) = 1
  CALL h5aread_f( attr_id, H5T_NATIVE_DOUBLE, attr_data, data_dims, error)

  time = attr_data
  CALL h5aclose_f(attr_id, error) ! Close the attribute.
  CALL h5sclose_f(aspace_id, error) ! Terminate access to the data space.

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! open attribute (domain_length)
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL h5aopen_f(dset_id, "domain_size", attr_id, error)

  ! Get dataspace and read
  CALL h5aget_space_f(attr_id, aspace_id, error)
  data_dims(1) = 2
  CALL h5aread_f( attr_id, H5T_NATIVE_DOUBLE, attr_data2, data_dims, error)

  xl = attr_data2(1)
  yl = attr_data2(2)

  CALL h5aclose_f(attr_id, error) ! Close the attribute.
  CALL h5sclose_f(aspace_id, error) ! Terminate access to the data space.

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! open attribute (viscosity)
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! only recent files contain the attribute viscosity, so we need to check if it
  ! is actually there
  call h5aexists_f(dset_id, "viscosity", exists, error)

  if (exists) then
    CALL h5aopen_f(dset_id, "viscosity", attr_id, error)
    ! Get dataspace and read
    CALL h5aget_space_f(attr_id, aspace_id, error)
    data_dims(1) = 1
    CALL h5aread_f( attr_id, H5T_NATIVE_DOUBLE, attr_data, data_dims, error)

    viscosity = attr_data
    CALL h5aclose_f(attr_id, error) ! Close the attribute.
    CALL h5sclose_f(aspace_id, error) ! Terminate access to the data space.
  else
    write (*,*) "the file did not contain the viscosity!"
    viscosity = 0.d0
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! open attribute (sizes)
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL h5aopen_f(dset_id, "nxyz", attr_id, error)

  ! Get dataspace and read
  CALL h5aget_space_f(attr_id, aspace_id, error)
  data_dims(1) = 2
  CALL h5aread_f( attr_id, H5T_NATIVE_INTEGER, attr_data3, data_dims, error)

  nx = attr_data3(1)
  ny = attr_data3(2)

  CALL h5aclose_f(attr_id, error) ! Close the attribute.
  CALL h5sclose_f(aspace_id, error) ! Terminate access to the data space.

  CALL h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
  CALL h5fclose_f(file_id, error) ! Close the file.
  CALL h5close_f(error)  ! Close FORTRAN interface.
end subroutine Fetch_attributes_2d_openmp


!-----------------------------------------------------------------------------
! This function returns, to a given filename, the corresponding dataset name
! in the hdf5 file, following flusi conventions (folder/ux_0000.h5 -> "ux")
!-----------------------------------------------------------------------------
character(len=strlen)  function get_dsetname(fname)
  implicit none
  character(len=*), intent(in) :: fname
  ! extract dsetname (from "/" until "_", excluding both)
  get_dsetname  = fname  ( index(fname,'/',.true.)+1:index( fname, '_',.true. )-1 )
  return
end function get_dsetname

end module
