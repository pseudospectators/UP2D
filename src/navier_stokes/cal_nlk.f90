module rhs

  implicit none
  contains

!-------------------------------------------------------------------------------
! Compute right hand side of penalized Navier-Stokes, without the pressure
! the resulting field is NOT divergence free, add the pressure later!
!-------------------------------------------------------------------------------
subroutine cal_nlk (time, u, uk, vor, nlk)
  use share_vars
  implicit none
  real (kind=pr), intent (in) :: 					time
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: 		vor
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: 	u
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: 	uk
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (out) :: 	nlk
  real (kind=pr), dimension (:,:), allocatable :: work1, work2
  real (kind=pr), dimension (:,:,:), allocatable :: sp_tmp
  integer :: ix,iy
  real(kind=pr) :: theta, u_parallel

  allocate(work1(0:nx-1, 0:ny-1), work2(0:nx-1, 0:ny-1))
  allocate( sp_tmp(0:nx-1, 0:ny-1,1:2) )

  ! compute vorticity. the curl is computed in F-space, then the result is
  ! transformed back to x-space (for the non-linear term)
  call curl (uk, work1)
  call ifft( work1, vor)

  ! Compute non-linear term and penalization term. Note products are evaluated
  ! in x-space (pseudospectral code). this is the classical term, vor x u - chi/eta (u-us)
  !$omp parallel do private(iy)
  do iy=0,ny-1
    work1(:,iy) = +vor(:,iy)*u(:,iy,2) -mask(:,iy)*(u(:,iy,1)-us(:,iy,1))
    work2(:,iy) = -vor(:,iy)*u(:,iy,1) -mask(:,iy)*(u(:,iy,2)-us(:,iy,2))
  enddo
  !$omp end parallel do

  !-- non-linear terms to fourier space
  call fft ( work1, nlk(:,:,1) )
  call fft ( work2, nlk(:,:,2) )

  !-- sponge term
  if (use_sponge == 1) then
    ! apply sponge penalization to vorticity
    !$omp parallel do private(iy)
    do iy=0,ny-1
      work1(:,iy) = -mask_sponge(:,iy)*vor(:,iy)/eps_sponge
    enddo
    !$omp end parallel do
    call fft(work1,work2)
    ! obtain the velocity
    call vorticity2velocity( work2, sp_tmp )
    ! add sponge term to NL terms (in F-space)
    nlk = nlk + sp_tmp
  endif

  !$omp parallel do private(iy)
  do iy=0,ny-1
    nlk(:,iy,1) = nlk(:,iy,1)*dealiase(:,iy)
    nlk(:,iy,2) = nlk(:,iy,2)*dealiase(:,iy)
  enddo
  !$omp end parallel do

  deallocate(work1, work2)
  deallocate(sp_tmp)
end subroutine cal_nlk

end module rhs
