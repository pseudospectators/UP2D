subroutine dealiase_mask
!  Elliptical dealiasing (truncation 2/3 rule)
  use vars
  implicit none
  integer :: kx, ky
  integer :: kx_max, ky_max 
  real(kind=pr)    :: kx_trunc, ky_trunc

  kx_max = (nx/2-1) 
  ky_max = (ny/2-1)
  kx_trunc = 2.d0 / 3.d0 * dble (kx_max)
  ky_trunc = 2.d0 / 3.d0 * dble (ky_max)
  dealiase = 1.d0

  !$omp parallel do private(kx,ky)
  do  ky = 0, ky_max
     do  kx = 0, kx_max
        if (dble(kx)**2/kx_trunc**2 + dble(ky)**2/ky_trunc**2 >= 1.d0) then
           dealiase (2*kx, 2*ky)     = 0.d0
           dealiase (2*kx+1, 2*ky)   = 0.d0
           dealiase (2*kx, 2*ky+1)   = 0.d0
           dealiase (2*kx+1, 2*ky+1) = 0.d0
        end if
     end do
  end do
  !$omp end parallel do

  write (*,*) "dealiase_mask is set..."
end subroutine