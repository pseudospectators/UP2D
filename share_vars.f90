module share_vars
  implicit none
  
  integer, parameter :: pr = kind (0.d0) 
  

  character (len=40), save :: simulation_name
  integer, save :: nx, ny, nt
  integer, save :: iDealias, iPenalization, iMoving, nPalettes=14

  real(kind=pr), save :: colorscale = 0.0 ! scaling for farge palette (vorticity)
  real(kind=pr), save :: xl,yl,dx,dy,x0,y0
  real(kind=pr), save :: Tmax, CFL, tsave, tdrag
  real(kind=pr), save :: nu, eps, pi, scalex, scaley
  real(kind=pr), save :: ux_mean, uy_mean
  
  character (len=40), save :: inicond, iMask, iMeanFlow, iMethod
  
  integer,parameter :: nlines=2048 ! maximum number of lines in PARAMS-file

  ! memory
  real (kind=pr), dimension (:,:), allocatable, save :: dealiase, mask
  real (kind=pr), dimension (:,:,:), allocatable, save :: us
  
end module share_vars

