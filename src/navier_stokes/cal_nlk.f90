module rhs

  implicit none
  contains

!-------------------------------------------------------------------------------
! Compute right hand side of penalized Navier-Stokes, without the pressure
! the resulting field is NOT divergence free, add the pressure later!
!-------------------------------------------------------------------------------
subroutine cal_nlk (time, u, uk, vor, nlk, mask, us, mask_sponge)
  use share_vars
  implicit none
  real(kind=pr),intent (in) :: 					time
  real(kind=pr),dimension(0:nx-1, 0:ny-1),intent(inout) :: vor, mask, mask_sponge
  real(kind=pr),dimension(0:nx-1, 0:ny-1,1:2),intent(inout) :: u, us
  real(kind=pr),dimension(0:nx-1, 0:ny-1,1:2),intent(inout) :: uk
  real(kind=pr),dimension(0:nx-1, 0:ny-1,1:2),intent(inout) :: nlk

  real(kind=pr),dimension(:,:), allocatable :: work1, work2
  integer :: ix,iy

  allocate(work1(0:nx-1, 0:ny-1), work2(0:nx-1, 0:ny-1))

  ! compute vorticity. the curl is computed in F-space, then the result is
  ! transformed back to x-space (for the non-linear term)
  call curl(uk, work1)
  call ifft(work1, vor)

  ! Compute non-linear term and penalization term. Note products are evaluated
  ! in x-space (pseudospectral code). this is the classical term, vor x u - chi/eta (u-us)
  !$omp parallel do private(iy)
  do iy=0,ny-1
    work1(:,iy) = +vor(:,iy)*u(:,iy,2) -mask(:,iy)*(u(:,iy,1)-us(:,iy,1))
    work2(:,iy) = -vor(:,iy)*u(:,iy,1) -mask(:,iy)*(u(:,iy,2)-us(:,iy,2))
  enddo
  !$omp end parallel do

  !-- non-linear terms to fourier space
  call fft(work1, nlk(:,:,1))
  call fft(work2, nlk(:,:,2))

  !-- sponge term
  call add_sponge_term(time, nlk, vor, mask_sponge, work1, work2 )

  !$omp parallel do private(iy)
  do iy=0,ny-1
    nlk(:,iy,1) = nlk(:,iy,1)*dealiase(:,iy)
    nlk(:,iy,2) = nlk(:,iy,2)*dealiase(:,iy)
  enddo
  !$omp end parallel do

  deallocate(work1, work2)
end subroutine cal_nlk

end module rhs
