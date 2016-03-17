module rhs

  implicit none
  contains

!-------------------------------------------------------------------------------
! Compute right hand side of penalized Navier-Stokes, without the pressure
! the resulting field is NOT divergence free, add the pressure later!
!-------------------------------------------------------------------------------
subroutine cal_nlk (time, u, uk, vor, nlk)
  use share_vars
  use FieldExport
  implicit none
  real (kind=pr), intent (in) :: 					time
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: 		vor
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: 	u
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: 	uk
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (out) :: 	nlk
  real (kind=pr), dimension (:,:), allocatable :: work1, work2
  integer :: ix,iy
  real(kind=pr) :: theta, u_parallel

  allocate(work1(0:nx-1, 0:ny-1), work2(0:nx-1, 0:ny-1))

  !-- compute vorticity
  call curl (uk, work1)
  call cofitxy(work1, vor)

  !-- NL term + penal
  ! !$omp parallel do private(iy,ix,theta,u_parallel)
  ! do iy=0,ny-1
  !   do ix=0,nx-1
  !     theta = atan2(dble(iy)*dy-y0,dble(ix)*dx-x0)
  !     u_parallel = -u(ix,iy,1)*sin(theta)   + u(ix,iy,2)*cos(theta)
  !     work1(ix,iy) = +vor(ix,iy)*u(ix,iy,2) -mask(ix,iy)*(u(ix,iy,1)+u_parallel*sin(theta))
  !     work2(ix,iy) = -vor(ix,iy)*u(ix,iy,1) -mask(ix,iy)*(u(ix,iy,2)-u_parallel*cos(theta))
  !   enddo
  ! enddo
  ! !$omp end parallel do


  !$omp parallel do private(iy)
  do iy=0,ny-1
    work1(:,iy) = +vor(:,iy)*u(:,iy,2) -mask(:,iy)*(u(:,iy,1)-us(:,iy,1))
    work2(:,iy) = -vor(:,iy)*u(:,iy,1) -mask(:,iy)*(u(:,iy,2)-us(:,iy,2))
  enddo
  !$omp end parallel do

  !-- NLK to fourier space
  call coftxy ( work1, nlk(:,:,1) )
  call coftxy ( work2, nlk(:,:,2) )

  !$omp parallel do private(iy)
  do iy=0,ny-1
    nlk(:,iy,1) = nlk(:,iy,1)*dealiase(:,iy)
    nlk(:,iy,2) = nlk(:,iy,2)*dealiase(:,iy)
  enddo
  !$omp end parallel do

  deallocate(work1, work2)
end subroutine cal_nlk

end module rhs
