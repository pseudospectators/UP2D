module share_vars
  implicit none
  integer, parameter :: pr = kind (0.d0)
  integer, parameter :: strlen=80

  integer, save :: nx, ny
  integer, save :: nt
  integer, save :: iDealias, nPalettes=14

  real(kind=pr), save :: colorscale = 0.d0 ! scaling for farge palette (vorticity)
  real(kind=pr), save :: xl,yl,dx,dy,x0,y0, delta
  real(kind=pr), save :: Tmax, CFL, tsave, tdrag, dt_fixed
  real(kind=pr), save :: nu, eps, pi
  real(kind=pr), save :: ux_mean, uy_mean

  character (len=strlen), save :: inicond, iMask, iMeanFlow, iMethod
  character (len=strlen), save :: name, ipressure, iActive

  ! deliberately reduce code to second order FD?
  logical, save :: FD_2nd = .false.

  integer,parameter :: nlines=2048 ! maximum number of lines in PARAMS-file

  ! memory
  real (kind=pr), dimension (:,:), allocatable, save :: dealiase, mask
  real (kind=pr), dimension (:,:,:), allocatable, save :: us, uex, normals, u_BC
  real (kind=pr), dimension (:,:), allocatable, save :: phi

contains
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
end module share_vars
