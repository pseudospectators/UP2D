subroutine cofdx (fk, fk_dx)
!---------------------------------------------------------------
!     calculation of d/dx in the fourier-space ==> *ik
!     for the first index
!     scaling included
!     ck = ak + i bk with ak= fk(2k,l) and bk= fk(2k+1,l)
!                         for k=0, kx-1 and all l=0,ny-1
!---------------------------------------------------------------
  use share_vars
  implicit none
  integer kx, ky, k
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) :: fk   
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: fk_dx   
  real (kind=pr) :: scale1, keff

  scale1 = 2.d0*pi/xl
  
  !$omp parallel do private(kx,ky,keff,k)
  do ky = 0, ny-1
     k = 0
     do kx = 0, nx-2, 2
        if (FD_2nd) then
          ! second order
          keff = dsin( dx * dble(k) ) / dx
          fk_dx (kx, ky) = -fk (kx+1, ky) * keff * scale1
          k = k+1
        else
          ! spectral order        
          fk_dx (kx, ky) = -fk (kx+1, ky) * dble (kx/2) * scale1
        endif
     end do  
     k=0
     do kx = 1, nx-1, 2
        if (FD_2nd) then
          ! second order
          keff = dsin( dx * dble(k) ) / dx
          fk_dx (kx, ky) = fk (kx-1, ky) * keff * scale1
          k = k+1
        else
          ! spectral order        
          fk_dx (kx, ky) = fk (kx-1, ky) * dble ((kx-1)/2) * scale1
        endif
     end do
  end do
  !$omp end parallel do

end subroutine cofdx


subroutine cofdy (fk, fk_dy)
!---------------------------------------------------------------
!     calculation of d/dy in the fourier-space ==> *ik
!     for the second index
!     scaling included
!     ck = ak + i bk with ak= fk(k,2l) and bk= fk(k,2l+1)
!                         for l=0, ky-1 and all k=0,nx-1
!---------------------------------------------------------------
  use share_vars
  implicit none
  integer :: kx, ky, k
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) :: fk   
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: fk_dy   
  real (kind=pr) :: scale1, keff

  scale1 = 2.d0*pi/yl

  !$omp parallel do private(kx,ky,keff,k)
  do kx = 0, nx-1
     k = 0
     do ky = 0, ny-2, 2 ! loop over real parts of the FFT signal
        if (FD_2nd) then
          keff = dsin( dy * dble(k) ) / dy
          fk_dy (kx, ky) = -fk (kx, ky+1) * keff * scale1
          k = k +1
        else        
          fk_dy (kx, ky) = -fk (kx, ky+1) * dble (ky/2) * scale1
        endif
     end do
     
     k = 0
     do ky = 1, ny-1, 2 ! loop over imaginary parts of the FFT
        if (FD_2nd) then
          keff = dsin( dy * dble(k) ) / dy
          fk_dy (kx, ky) = fk (kx, ky-1) * keff * scale1
          k = k +1
        else        
          ! actually, both (real/imag) have of course the same wavenumber. 
          fk_dy (kx, ky) = fk (kx, ky-1) * dble ((ky-1)/2) * scale1
        endif
     end do
  end do
  !$omp end parallel do

end subroutine


subroutine cofdxdx (fk, fk_dxdx)
!---------------------------------------------------------------
!     Calculation of d/dx in the Fourier-space ==> -kx^2
!     for the first index
!     Scaling included
!     Ck = Ak + i Bk with Ak= Fk(2k,l) and Bk= Fk(2k+1,l)
!                         for k=0, KX-1 and all l=0,NY-1
!---------------------------------------------------------------
  use share_vars
  implicit none
  integer kx, ky
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) :: fk
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: fk_dxdx 
  real (kind=pr) :: scale

  scale = (2.d0*pi/xl)**2

  !$omp parallel do private(kx,ky)
  do ky = 0, ny-1
     do kx = 0, nx-2, 2
        fk_dxdx (kx, ky)   = - fk (kx, ky)   * dble (kx/2)**2 * scale
        fk_dxdx (kx+1, ky) = - fk (kx+1, ky) * dble (kx/2)**2 * scale
     end do
  end do
  !$omp end parallel do

end subroutine


subroutine cofdxdy (fk, fk_dxdy)
!---------------------------------------------------------------
!     Calculation of d/dx in the Fourier-space ==> -kx^2
!     for the first index
!     Scaling included
!     Ck = Ak + i Bk with Ak= Fk(2k,l) and Bk= Fk(2k+1,l)
!                         for k=0, KX-1 and all l=0,NY-1
!---------------------------------------------------------------
  use share_vars
  implicit none
  integer kx, kx_max, ky, ky_max
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) :: fk
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: fk_dxdy 
  real (kind=pr) :: fac, scale

  scale = (2.d0*pi/xl) * (2.d0*pi/yl)
  kx_max = (nx/2-1) 
  ky_max = (ny/2-1)

  !$omp parallel do private(kx,ky,fac)
  do ky = 0, ky_max
     do kx = 0, kx_max
        fac = - dble(kx) * dble (ky) * scale
        fk_dxdy (2*kx, 2*ky)     =  fac * fk (2*kx+1, 2*ky+1)
        fk_dxdy (2*kx+1, 2*ky)   = -fac * fk (2*kx, 2*ky+1)
        fk_dxdy (2*kx, 2*ky+1)   = -fac * fk (2*kx+1, 2*ky)
        fk_dxdy (2*kx+1, 2*ky+1) =  fac * fk (2*kx, 2*ky)
     end do
  end do
  !$omp end parallel do

end subroutine


subroutine cofdydy (fk, fk_dydy)
!---------------------------------------------------------------
!     Calculation of d/dx in the Fourier-space ==> -kx^2
!     for the first index
!     Scaling included
!     Ck = Ak + i Bk with Ak= Fk(2k,l) and Bk= Fk(2k+1,l)
!                         for k=0, KX-1 and all l=0,NY-1
!---------------------------------------------------------------
  use share_vars
  implicit none
  integer kx, ky
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in)  :: fk
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: fk_dydy 
  real (kind=pr) :: scale

  scale = (2.d0*pi/yl)**2

  !$omp parallel do private(kx,ky)
  do kx = 0, nx-1
     do ky = 0, ny-2, 2
        fk_dydy (kx, ky)   = - fk (kx, ky)   * dble (ky/2)**2 * scale
        fk_dydy (kx, ky+1) = - fk (kx, ky+1) * dble (ky/2)**2 * scale
     end do
  end do
  !$omp end parallel do

end subroutine


subroutine poisson (f, ans)
! Calculate solution to the Poisson equation f = -grad^2 ans in Fourier space
  use share_vars
  implicit none
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) :: f
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: ans
  integer :: kx, kx_max, ky, ky_max
  real (kind=pr) :: quot, scalex,scaley, kx2eff, ky2eff
  kx_max = (nx/2-1) 
  ky_max = (ny/2-1)
  
  scalex = (2.d0*pi/xl)**2
  scaley = (2.d0*pi/yl)**2

  !$omp parallel do private(kx,ky,quot,kx2eff,ky2eff) 
  do ky = 0, ky_max
     do kx = 0, kx_max
        if ( (kx == 0) .and. (ky == 0) ) then
           ans (2*kx, 2*ky)     = 0.d0
           ans (2*kx+1, 2*ky)   = 0.d0
           ans (2*kx, 2*ky+1)   = 0.d0
           ans (2*kx+1, 2*ky+1) = 0.d0
        else
          if (FD_2nd) then
            ! second order 
            kx2eff = ( 2.0 - 2.0*cos( dx * real(kx) ) ) / dx**2
            ky2eff = ( 2.0 - 2.0*cos( dy * real(ky) ) ) / dy**2
            quot = (kx2eff*scalex + ky2eff*scaley)
            ans (2*kx, 2*ky)     =  f (2*kx, 2*ky) / quot
            ans (2*kx+1, 2*ky)   =  f (2*kx+1, 2*ky) / quot
            ans (2*kx, 2*ky+1)   =  f (2*kx, 2*ky+1) / quot
            ans (2*kx+1, 2*ky+1) =  f (2*kx+1, 2*ky+1) / quot            
          else
            ! spectral accuracy
            quot = (dble(kx**2)*scalex + dble(ky**2)*scaley)
            ans (2*kx, 2*ky)     =  f (2*kx, 2*ky) / quot
            ans (2*kx+1, 2*ky)   =  f (2*kx+1, 2*ky) / quot
            ans (2*kx, 2*ky+1)   =  f (2*kx, 2*ky+1) / quot
            ans (2*kx+1, 2*ky+1) =  f (2*kx+1, 2*ky+1) / quot
          end if
        end if
     end do
  end do
  !$omp end parallel do
end subroutine poisson




subroutine curl( uk, vortk )
  ! this routine computes the streamfunction and the velocity field given a vorticity in fourier space
  use share_vars
  implicit none
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: uk
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out)  :: vortk  
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: work1, work2
  integer :: iy
  
  call cofdx(uk(:,:,2), work1)
  call cofdy(uk(:,:,1), work2)
  
  !$omp parallel do private(iy)
  do iy=0,ny-1
      vortk(:,iy) = work1(:,iy) - work2(:,iy)
  enddo
  !$omp end parallel do    
end subroutine



subroutine vorticity2velocity( vortk, uk )
  ! this routine computes the streamfunction and the velocity field given a vorticity in fourier space
  use share_vars
  implicit none
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (out) :: uk
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in)  :: vortk  
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: stream
  integer :: iy

  !--Solve poisson equation for stream function
  call poisson (vortk, stream)
  !--Calculate x-derivative of stream function
  call cofdx (-stream, uk(:,:,2))
  !--Calculate y-derivative of stream function
  call cofdy ( stream, uk(:,:,1)) 
end subroutine vorticity2velocity


subroutine cal_vis (dt, vis)
!---------------------------------------------------------------
!  Calculate viscous term for time advancement
!  exp (-nu*k^2*dt)
!  Optimized for vectorization (Dmitry, Feb 1, 2008)
!---------------------------------------------------------------
  use share_vars
  implicit none
  integer :: kx, kx_max, ky, ky_max
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent(out) :: vis
  real (kind=pr), intent (in) :: dt
  real (kind=pr) :: coefx, coefy,scalex,scaley, kx2eff, ky2eff

  scalex = (2.d0*pi/xl)**2
  scaley = (2.d0*pi/yl)**2
  
  kx_max = (nx/2-1) 
  ky_max = (ny/2-1)
  coefx = - dt * nu * scalex
  coefy = - dt * nu * scaley

  !$omp parallel do private(kx,ky)
  do ky = 0, ky_max
    do kx = 0, kx_max
        if (FD_2nd) then
          ! second order accuracy
          kx2eff = ( 2.0 - 2.0*cos( dx * real(kx) ) ) / dx**2
          ky2eff = ( 2.0 - 2.0*cos( dy * real(ky) ) ) / dy**2
          vis (2*kx, 2*ky) = exp( coefx * kx2eff + coefy * ky2eff )
        else
          ! spectral accuracy
          vis (2*kx, 2*ky) = dexp( coefx * dble(kx**2) + coefy * dble(ky**2) )
        endif
    enddo
  enddo
  !$omp end parallel do 


  !$omp parallel do private(kx,ky)
  do ky=0,ny/2-1
     vis (1:nx-1:2,2*ky) = vis (0:nx-2:2,2*ky)
     vis (0:nx-2:2,2*ky+1) = vis (0:nx-2:2,2*ky)
     vis (1:nx-1:2,2*ky+1) = vis (0:nx-2:2,2*ky)
  end do
  !$omp end parallel do 

end subroutine cal_vis


function max_divergence(uk)
  use share_vars
  implicit none
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: uk
  real(kind=pr):: max_divergence
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: work1, work2, work3
  call cofdx( uk(:,:,1),work1 )
  call cofdy( uk(:,:,2),work2 )
  call cofitxy(work1+work2,work3)
  max_divergence = maxval( dabs((1.0d0-mask*eps)*work3) )
end function