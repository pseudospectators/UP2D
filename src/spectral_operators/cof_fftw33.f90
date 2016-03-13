!====================================================================
!====================================================================
!
!     Fourier transform subroutines using MKL library
!
!====================================================================
!====================================================================


module fftw33_descriptors
  use share_vars  
  implicit none
  integer, save :: nthreads ! this holds (globally) the number of threads
  integer*8, save :: Desc_Handle_xy_f, Desc_Handle_xy_b
end module fftw33_descriptors


subroutine fft_initialize
!====================================================================
!     Allocate memory and initialize FFT
!====================================================================
  use share_vars
  use fftw33_descriptors
  use omp_lib
  implicit none
  include 'fftw3.f'

  integer :: ierr, init_flag = FFTW_MEASURE, i,id,ix,iy
  real (kind=pr), dimension (0:nx-1, 0:ny-1) ::  fftwdata_r
  complex (kind=pr), dimension(0:nx/2, 0:ny-1) :: fftwdata_c
  real(kind=pr), dimension(0:nx-1,0:ny-1) :: work1, work2, work3
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2) :: u
  real (kind=pr)::e

  !----------------------------------------------------
  ! here, we determine how many thready we have available
  !----------------------------------------------------  
  !$omp parallel private(id)
  id = omp_get_thread_num()
  !$omp barrier
  if ( id == 0 ) then
    nthreads = omp_get_num_threads()
  end if
  !$omp end parallel  
  write (*,*) "fftw wrapper: nthreads=", nthreads
  
  
  ! Initialize two-dimensional FFT
    ! set up threads
    call dfftw_init_threads(ierr)
    call dfftw_plan_with_nthreads(nthreads)
    ! forward
    call dfftw_plan_dft_r2c_2d(Desc_Handle_xy_f, nx, ny, fftwdata_r, fftwdata_c, init_flag)
    ! backward
    call dfftw_plan_dft_c2r_2d(Desc_Handle_xy_b, nx, ny, fftwdata_c, fftwdata_r, init_flag)

  !-----------------------------------------
  ! UNIT TESTING
  !-----------------------------------------  
  do ix=0,nx-1
    work1(ix,:) = dsin( dble(ix)*dx*2.d0*pi/xl )
    work3(ix,:) = dcos( dble(ix)*dx*2.d0*pi/xl )*2.d0*pi/xl
  enddo
  
  call coftxy(work1,work2)
  call cofdx(work2,work1)
  call cofitxy(work1,work2)
  e = sum((work2-work3)**2)*dx*dy
  write(*,*) "FFT testing error=", e

  if (e >= 1e-10) then
    write(*,*) "FFT unit test failed..."
    if (FD_2nd.eqv..false.) then
    stop
    endif
  endif

  !-----------------------------------------
  ! UNIT TESTING (y)
  !-----------------------------------------  
  do iy=0,ny-1
    work1(:,iy) = dsin( dble(iy)*dy*2.d0*pi/yl )
    work3(:,iy) = dcos( dble(iy)*dy*2.d0*pi/yl )*2.d0*pi/yl
  enddo
  
  call coftxy(work1,work2)
  call cofdy(work2,work1)
  call cofitxy(work1,work2)  
  e = sum((work2-work3)**2)*dx*dy
  write(*,*) "FFT testing error=", e

  if (e >= 1e-10) then
    write(*,*) "FFT unit test failed..."
    if (FD_2nd.eqv..false.) then
    stop
    endif
  endif
    
  
end subroutine fft_initialize



subroutine fft_initialize_best
!====================================================================
!     Allocate memory and initialize FFT
!====================================================================
  use share_vars
  use fftw33_descriptors
  use omp_lib
  implicit none
  include 'fftw3.f'

  integer :: ierr, init_flag = FFTW_PATIENT, i,id
  real (kind=pr), dimension (0:nx-1, 0:ny-1) ::  fftwdata_r
  complex (kind=pr), dimension(0:nx/2, 0:ny-1) :: fftwdata_c

  !----------------------------------------------------
  ! here, we determine how many thready we have available
  !----------------------------------------------------  
  !$omp parallel private(id)
  id = omp_get_thread_num()
  !$omp barrier
  if ( id == 0 ) then
    nthreads = omp_get_num_threads()
  end if
  !$omp end parallel  
  write (*,*) "fftw wrapper: nthreads=", nthreads
  
  
  ! Initialize two-dimensional FFT
    write (*,*) 'double precision'
    ! set up threads
    call dfftw_init_threads(ierr)
    call dfftw_plan_with_nthreads(nthreads)
    ! forward
    call dfftw_plan_dft_r2c_2d(Desc_Handle_xy_f, nx, ny, fftwdata_r, fftwdata_c, init_flag)
    ! backward
    call dfftw_plan_dft_c2r_2d(Desc_Handle_xy_b, nx, ny, fftwdata_c, fftwdata_r, init_flag)

end subroutine



subroutine fft_free
!====================================================================
!     Free memory allocated for FFT
!====================================================================
  use share_vars
  use fftw33_descriptors
  implicit none

    call dfftw_destroy_plan(Desc_Handle_xy_f)
    call dfftw_destroy_plan(Desc_Handle_xy_b)
    call dfftw_cleanup_threads


end subroutine fft_free




subroutine coftxy (f, fk)
!====================================================================
!     Calculation of the Fourier-coefficients of a real function
!     along x (1st index) and y (2nd index)
!     FILTERING
!====================================================================
  use share_vars
  use fftw33_descriptors
  implicit none

  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  f
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  fk
  complex (kind=pr), dimension(0:nx/2, 0:ny-1) :: ft
  real (kind=pr) :: norm
  integer :: ix, iy

  ! Compute forward FFT
     call dfftw_execute_dft_r2c(Desc_Handle_xy_f,f,ft)


  ! normalization
  norm = 1.d0 / ( real(nx) * real(ny) )

  ! complex 2d Fourier coefficients are used to
  ! compute coefficients before sin and cos (Temperton style)
  !$omp parallel do private(ix)
  do ix=0,nx/2-1
     fk(2*ix,0) = real( ft(ix,0) ) * norm
     fk(2*ix+1,0) = imag( ft(ix,0) ) * norm
     fk(2*ix,1) = 0.d0
     fk(2*ix+1,1) = 0.d0
  end do
  !$omp end parallel do

  !$omp parallel do private(iy)
  do iy=1,ny/2-1
     fk(0:nx-2:2,2*iy) = 0.5d0 * ( real( ft(0:nx/2-1,iy)) &
                            + real( ft(0:nx/2-1,ny-iy)) ) * norm
     fk(0:nx-2:2,2*iy+1) = 0.5d0 * ( imag( ft(0:nx/2-1,iy)) &
                            - imag( ft(0:nx/2-1,ny-iy)) ) * norm
     fk(1:nx-1:2,2*iy) = 0.5d0 * ( imag( ft(0:nx/2-1,iy)) &
                            + imag( ft(0:nx/2-1,ny-iy)) ) * norm
     fk(1:nx-1:2,2*iy+1) = 0.5d0 * ( - real( ft(0:nx/2-1,iy)) &
                            + real( ft(0:nx/2-1,ny-iy)) ) * norm
  end do
  !$omp end parallel do

  ! mode KF left unconsidered => filtering
end subroutine coftxy



subroutine cofitxy (fk, f)
!====================================================================
!     Calculation of a real function from its Fourier coefficients
!     along x (1st index) and y (2nd index)
!     FILTERING
!====================================================================
  use share_vars
  use fftw33_descriptors
  implicit none

  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  fk
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  f
  complex (kind=pr), dimension(0:nx/2, 0:ny-1) :: ft
  integer :: ix, iy

  ! Arrange the data in a Temperton-like ordering

  !$omp parallel do private(ix)
  do ix=0,nx/2-1
     ft(ix,0) = cmplx( fk(2*ix,0), fk(2*ix+1,0), pr )
     ft(ix,ny/2) = 0.0d0
  end do
  !$omp end parallel do

  ft(nx/2,0) = 0.0d0
  ft(nx/2,ny/2) = 0.0d0

  !$omp parallel do private(iy)
  do iy=1,ny/2-1
     ft(0:nx/2-1,iy) = cmplx( fk(0:nx-2:2,2*iy) - fk(1:nx-1:2,2*iy+1), &
                          fk(0:nx-2:2,2*iy+1) + fk(1:nx-1:2,2*iy) )
     ft(0:nx/2-1,ny-iy) = cmplx( fk(0:nx-2:2,2*iy) + fk(1:nx-1:2,2*iy+1), &
                         -fk(0:nx-2:2,2*iy+1) + fk(1:nx-1:2,2*iy) )
     ft(nx/2,iy) = 0.0d0
     ft(nx/2,ny-iy) = 0.0d0
  end do
  !$omp end parallel do

  ! Compute inverse FFT
     call dfftw_execute_dft_c2r(Desc_Handle_xy_b,ft,f)


  ! mode KF left unconsidered => filtering
end subroutine cofitxy




subroutine cofts (f, fk, L, n)
!====================================================================
!     Calculation of the Fourier coefficients of n real functions
!     aligned in columns
!     FILTERING
!     MK library is used
!====================================================================
  use share_vars
  use fftw33_descriptors
  implicit none
  include 'fftw3.f'

  integer, intent (in) :: L, n
  integer :: j
  real (kind=pr), dimension (0:L+1, 0:n-1) :: ft
  real (kind=pr), dimension (0:L-1, 0:n-1), intent (in) ::  f
  real (kind=pr), dimension (0:L-1, 0:n-1), intent (out) ::  fk
  real (kind=pr) :: norm
  integer*8 :: Desc_Handle_1D_f


     call dfftw_plan_many_dft_r2c(Desc_Handle_1D_f, 1, L, n, f, 0, 1, L, ft, 0, 1, L/2+1, FFTW_ESTIMATE)
     call dfftw_execute(Desc_Handle_1D_f)
     call dfftw_destroy_plan(Desc_Handle_1D_f)


  ! Output
  norm = 1.d0 / real(L)
  !$omp parallel do private(j)
  do j=0,n-1
     fk(:,j) = ft(0:L-1,j) * norm
  end do
  !$omp end parallel do

!      last mode M=L, L+1; mode KF left unconsidered => filtering
end subroutine cofts



subroutine cofits (fk, f, L, n)
!====================================================================
!     Calculation of n real functions from their Fourier coefficients
!     aligned in columns
!     FILTERING
!     MK library is used
!====================================================================
  use share_vars
  use fftw33_descriptors
  implicit none
  include 'fftw3.f'

  integer, intent (in) :: L, n
  integer :: j
  real (kind=pr), dimension (0:L+1, 0:n-1) :: ft
  real (kind=pr), dimension (0:L-1, 0:n-1), intent (in) ::  fk
  real (kind=pr), dimension (0:L-1, 0:n-1), intent (out) ::  f
  integer*8 :: Desc_Handle_1D_b

  ! Input
  !$omp parallel do private(j)
  do j=0,n-1
     ft(0:L-1,j) = fk(:,j)
     ft(L:L+1,j) = 0.d0
  end do
  !$omp end parallel do

     call dfftw_plan_many_dft_c2r(Desc_Handle_1D_b, 1, L, n, ft, 0, 1, L/2+1, f, 0, 1, L, FFTW_ESTIMATE)
     call dfftw_execute(Desc_Handle_1D_b)
     call dfftw_destroy_plan(Desc_Handle_1D_b)


  ! last mode M=L, L+1; mode KF left unconsidered => filtering
end subroutine cofits

