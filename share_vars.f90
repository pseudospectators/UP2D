module share_vars
  implicit none  
  integer, parameter :: pr = kind (0.d0) 
  
  integer, save :: nx, ny, nt
  integer, save :: iDealias, nPalettes=14

  real(kind=pr), save :: colorscale = 0.d0 ! scaling for farge palette (vorticity)
  real(kind=pr), save :: xl,yl,dx,dy,x0,y0
  real(kind=pr), save :: Tmax, CFL, tsave, tdrag, dt_fixed
  real(kind=pr), save :: nu, eps, pi
  real(kind=pr), save :: ux_mean, uy_mean
  
  character (len=40), save :: inicond, iMask, iMeanFlow, iMethod
  character (len=40), save :: name
  
  integer,parameter :: nlines=2048 ! maximum number of lines in PARAMS-file

  ! memory
  real (kind=pr), dimension (:,:), allocatable, save :: dealiase, mask
  real (kind=pr), dimension (:,:,:), allocatable, save :: us, uex, normals
  real (kind=pr), dimension (:,:), allocatable, save :: dist
end module share_vars

