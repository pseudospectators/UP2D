
!====================================================================
!====================================================================
!
!     Fourier transform subroutines using MKL library
!
!====================================================================
!====================================================================


module mkl_fft_descriptors
  use mkl_dfti
  implicit none
  type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_xy, Desc_Handle_x, Desc_Handle_y, Desc_Handle_xy_small, Desc_Handle_xy_big
end module mkl_fft_descriptors


subroutine fft_initialize
!====================================================================
!     Allocate memory and initialize FFT
!====================================================================
  use share_vars
  use mkl_dfti
  use mkl_fft_descriptors
  implicit none
  integer :: Status, L(2), strides_in(3)
  ! ------------------------------------------------------------------------------
  ! ------ Two-dimensional FFT ------
  ! ------------------------------------------------------------------------------
  ! Prepare FFT data
  ! Dataset dimensions
  L(1) = nx
  L(2) = ny

  ! Strides for a (nx+2)*(ny+2) array
  strides_in(1) = 0
  strides_in(2) = 1
  strides_in(3) = nx+2

  ! Initialize FFT
  if ( precision(0.0) .eq. precision(0d0)) then
!    write (*,*) 'double precision'
    Status = DftiCreateDescriptor( Desc_Handle_xy, DFTI_DOUBLE, DFTI_REAL, 2, L )
    Status = DftiSetValue( Desc_Handle_xy, DFTI_PRECISION, DFTI_DOUBLE) ! Set double precision
  else
!    write (*,*) 'single precision'
    Status = DftiCreateDescriptor( Desc_Handle_xy, DFTI_SINGLE, DFTI_REAL, 2, L )
    Status = DftiSetValue( Desc_Handle_xy, DFTI_PRECISION, DFTI_SINGLE) ! Set single precision
  endif

  Status = DftiSetValue( Desc_Handle_xy, DFTI_FORWARD_SCALE, 1.0 / ( real(nx) * real(ny) ) ) ! Set forward scaling
  Status = DftiSetValue( Desc_Handle_xy, DFTI_BACKWARD_SCALE, 1.0 ) ! Set backward scaling
  Status = DftiSetValue( Desc_Handle_xy, DFTI_INPUT_STRIDES, strides_in) ! Set strides

  Status = DftiCommitDescriptor( Desc_Handle_xy )


  ! ------------------------------------------------------------------------------
  ! ------ One-dimensional FFT in x ------
  ! ------------------------------------------------------------------------------
  
  ! Strides for a (nx+2)*ny array
  strides_in(1) = 0
  strides_in(2) = 1

  ! Initialize FFT
  if ( precision(0.0) .eq. precision(0d0)) then
!    write (*,*) 'double precision'
    Status = DftiCreateDescriptor( Desc_Handle_x, DFTI_DOUBLE, DFTI_REAL, 1, nx )
    Status = DftiSetValue( Desc_Handle_x, DFTI_PRECISION, DFTI_DOUBLE) ! Set double precision
  else
!    write (*,*) 'single precision'
    Status = DftiCreateDescriptor( Desc_Handle_x, DFTI_SINGLE, DFTI_REAL, 1, nx )
    Status = DftiSetValue( Desc_Handle_x, DFTI_PRECISION, DFTI_SINGLE) ! Set single precision
  endif

  Status = DftiSetValue( Desc_Handle_x, DFTI_FORWARD_SCALE, 1.0 / real(nx)  ) ! Set forward scaling
  Status = DftiSetValue( Desc_Handle_x, DFTI_BACKWARD_SCALE, 1.0 ) ! Set backward scaling
  Status = DftiSetValue( Desc_Handle_x, DFTI_INPUT_STRIDES, strides_in) ! Set strides
  Status = DftiSetValue( Desc_Handle_x, DFTI_NUMBER_OF_TRANSFORMS, ny)
  Status = DftiSetValue( Desc_Handle_x, DFTI_INPUT_DISTANCE, nx+2)

  Status = DftiCommitDescriptor( Desc_Handle_x )



  ! ------ One-dimensional FFT in y ------

  ! Strides for a nx*(ny+2) array
  strides_in(1) = 0
  strides_in(2) = nx

  ! Initialize FFT
  if ( precision(0.0) .eq. precision(0d0)) then
    !write (*,*) 'double precision'
    Status = DftiCreateDescriptor( Desc_Handle_y, DFTI_DOUBLE, DFTI_REAL, 1, ny )
    Status = DftiSetValue( Desc_Handle_y, DFTI_PRECISION, DFTI_DOUBLE) ! Set double precision
  else
    !write (*,*) 'single precision'
    Status = DftiCreateDescriptor( Desc_Handle_y, DFTI_SINGLE, DFTI_REAL, 1, ny )
    Status = DftiSetValue( Desc_Handle_y, DFTI_PRECISION, DFTI_SINGLE) ! Set single precision
  endif

  Status = DftiSetValue( Desc_Handle_y, DFTI_FORWARD_SCALE, 1.0 / real(ny)  ) ! Set forward scaling
  Status = DftiSetValue( Desc_Handle_y, DFTI_BACKWARD_SCALE, 1.0 ) ! Set backward scaling
  Status = DftiSetValue( Desc_Handle_y, DFTI_INPUT_STRIDES, strides_in) ! Set strides
  Status = DftiSetValue( Desc_Handle_y, DFTI_NUMBER_OF_TRANSFORMS, nx)
  Status = DftiSetValue( Desc_Handle_y, DFTI_INPUT_DISTANCE, 1)

  Status = DftiCommitDescriptor( Desc_Handle_y )


end subroutine fft_initialize




subroutine fft_free
!====================================================================
!     Free memory allocated for FFT
!====================================================================
  use share_vars
  use mkl_dfti
  use mkl_fft_descriptors
  implicit none
  integer :: Status

  Status = DftiFreeDescriptor(Desc_Handle_xy)
!    Status = DftiFreeDescriptor(Desc_Handle_xy_small)
!      Status = DftiFreeDescriptor(Desc_Handle_xy_big)
  Status = DftiFreeDescriptor(Desc_Handle_x)
  Status = DftiFreeDescriptor(Desc_Handle_y)

end subroutine fft_free



subroutine coftx (f, fk)
!====================================================================
!     Calculation of the Fourier coefficients of a real function
!     along x (1st index)
!     FILTERING
!====================================================================
  use share_vars
  use mkl_dfti
  use mkl_fft_descriptors
  implicit none
  integer :: j
  real (kind=pr), dimension ( (nx+2) * ny ) :: ft
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  f
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  fk
  integer :: Status

  ! Arrange the data in a 1D array
  !$omp parallel do private(j)
  do j = 0, ny-1
     ft( (nx+2)*j+1:(nx+2)*j+nx ) = f(:,j)
     ft( (nx+2)*j+nx+1:(nx+2)*j+nx+2 ) = 0.0
  enddo
  !$omp end parallel do

  ! Compute forward FFT
  Status = DftiComputeForward( Desc_Handle_x, ft)

  ! Arrange the data in a 2D array
  !$omp parallel do private(j)
  do j = 0, ny-1
     fk(:,j) = ft( (nx+2)*j+1:(nx+2)*j+nx )
  enddo
  !$omp end parallel do

  ! last mode M=NX, NX+1; mode KF left unconsidered => filtering
end subroutine coftx



subroutine cofitx (fk, f)
!====================================================================
!     Calculation of a real function from its Fourier coefficients
!     along x (1st index)
!     FILTERING
!====================================================================
  use share_vars
  use mkl_dfti
  use mkl_fft_descriptors
  implicit none
  integer :: j
  real (kind=pr), dimension ( (nx+2) * ny ) :: ft
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  fk
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  f
  integer :: Status

  ! Arrange the data in a 1D array
  !$omp parallel do private(j)
  do j = 0, ny-1
     ft( (nx+2)*j+1:(nx+2)*j+nx ) = fk(:,j)
     ft( (nx+2)*j+nx+1:(nx+2)*j+nx+2 ) = 0.0
  enddo
  !$omp end parallel do

  ! Compute forward FFT
  Status = DftiComputeBackward( Desc_Handle_x, ft)

  ! Arrange the data in a 2D array
  !$omp parallel do private(j)
  do j = 0, ny-1
     f(:,j) = ft( (nx+2)*j+1:(nx+2)*j+nx )
  enddo
  !$omp end parallel do

!      last mode M=NX, NX+1; mode KF left unconsidered => filtering
end subroutine cofitx



subroutine cofty (f, fk)
!====================================================================
!     Calculation of the Fourier-coefficients of a real function
!     along y (2nd index)
!     FILTERING
!====================================================================
  use share_vars
  use mkl_dfti
  use mkl_fft_descriptors
  implicit none
  integer :: j
  real (kind=pr), dimension ( nx * (ny+2) ) :: ft
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  f
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  fk
  integer :: Status

  ! Arrange the data in a 1D array
  !$omp parallel do private(j)
  do j = 0, ny-1
     ft( nx*j+1:nx*j+nx ) = f(:,j)
  enddo
  !$omp end parallel do
  ft( nx*ny+1:nx*(ny+2) ) = 0.0

  ! Compute forward FFT
  Status = DftiComputeForward( Desc_Handle_y, ft)

  ! Arrange the data in a 2D array
  !$omp parallel do private(j)
  do j = 0, ny-1
     fk(:,j) = ft( nx*j+1:nx*j+nx )
  enddo
  !$omp end parallel do

!      last mode M=NY, NY+1; mode KF left unconsidered => filtering
end subroutine cofty



subroutine cofity (fk, f)
!====================================================================
!     Calculation of a real function from its Fourier coefficients
!     along y (2nd index)
!     FILTERING
!====================================================================
  use share_vars
  use mkl_dfti
  use mkl_fft_descriptors
  implicit none
  integer :: j
  real (kind=pr), dimension ( nx * (ny+2) ) :: ft
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  fk
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  f
  integer :: Status

  ! Arrange the data in a 1D array
  !$omp parallel do private(j)
  do j = 0, ny-1
     ft( nx*j+1:nx*j+nx ) = fk(:,j)
  enddo
  !$omp end parallel do
  ft( nx*ny+1:nx*(ny+2) ) = 0.0

  ! Compute forward FFT
  Status = DftiComputeBackward( Desc_Handle_y, ft)

  ! Arrange the data in a 2D array
  !$omp parallel do private(j)
  do j = 0, ny-1
     f(:,j) = ft( nx*j+1:nx*j+nx )
  enddo
  !$omp end parallel do

!      last mode M=NY, NY+1; mode KF left unconsidered => filtering
end subroutine cofity





subroutine coftxy (f, fk)
!====================================================================
!     Calculation of the Fourier-coefficients of a real function
!     along x (1st index) and y (2nd index)
!     FILTERING
!====================================================================
  use share_vars
  use mkl_dfti
  use mkl_fft_descriptors
  implicit none
  integer :: j
  real (kind=pr), dimension ( (nx+2) * (ny+2) ) :: ft
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  f
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  fk
  integer :: Status

  ! Arrange the data in a 1D array
  !$omp parallel do private(j)
  do j = 0, ny-1
     ft( (nx+2)*j+1:(nx+2)*j+nx ) = f(:,j)
     ft( (nx+2)*j+nx+1:(nx+2)*j+nx+2 ) = 0.0
  enddo
  !$omp end parallel do
  ft( (nx+2)*ny+1:(nx+2)*(ny+2) ) = 0.0

  ! Compute forward FFT
  Status = DftiComputeForward( Desc_Handle_xy, ft)

  ! Arrange the result in a Temperton-like ordering
  fk(0:nx-1,0) = ft(1:nx)
  fk(:,1) = 0.0

  fk(0,0:ny-1) = ft(1:ny*(nx+2):(nx+2))
  fk(1,:) = 0.0

  !$omp parallel do private(j)
  do j = 1, ny/2-1
    fk(2:nx-2:2,2*j) = 0.5 * ( ft( (nx+2)*j+3:(nx+2)*j+nx-1:2 ) &
                           + ft( (nx+2)*(ny-j)+3:(nx+2)*(ny-j)+nx-1:2 ) )
    fk(2:nx-2:2,2*j+1) = 0.5 * ( ft( (nx+2)*j+4:(nx+2)*j+nx:2 ) &
                           - ft( (nx+2)*(ny-j)+4:(nx+2)*(ny-j)+nx:2 ) )
    fk(3:nx-1:2,2*j) = 0.5 * ( ft( (nx+2)*j+4:(nx+2)*j+nx:2 ) &
                           + ft( (nx+2)*(ny-j)+4:(nx+2)*(ny-j)+nx:2 ) )
    fk(3:nx-1:2,2*j+1) = 0.5 * ( - ft( (nx+2)*j+3:(nx+2)*j+nx-1:2 ) &
                           + ft( (nx+2)*(ny-j)+3:(nx+2)*(ny-j)+nx-1:2 ) )
  enddo
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
  use mkl_dfti
  use mkl_fft_descriptors
  implicit none
  integer :: j
  real (kind=pr), dimension ( (nx+2) * (ny+2) ) :: ft
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  fk
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  f
  integer :: Status

  ! Arrange the result in a Temperton-like ordering
  ft = 0.0
  ft(1:nx) = fk(0:nx-1,0)
  ft(2:ny*(nx+2)+1:(nx+2)) = 0.0
  ft(1:ny*(nx+2):(nx+2)) = fk(0,0:ny-1)

  !$omp parallel do private(j)
  do j = 1, ny/2-1
    ft( (nx+2)*j+3:(nx+2)*j+nx-1:2 ) = fk(2:nx-2:2,2*j) - fk(3:nx-1:2,2*j+1) ! real part, first columns (0..ny/2-1)
    ft( (nx+2)*j+4:(nx+2)*j+nx:2 ) =  fk(2:nx-2:2,2*j+1) + fk(3:nx-1:2,2*j) ! imaginary part, first columns (0..ny/2-1)
    ft( (nx+2)*(ny-j)+3:(nx+2)*(ny-j)+nx-1:2 ) = fk(2:nx-2:2,2*j) + fk(3:nx-1:2,2*j+1) ! real part, last columns
    ft( (nx+2)*(ny-j)+4:(nx+2)*(ny-j)+nx:2 ) = - fk(2:nx-2:2,2*j+1) + fk(3:nx-1:2,2*j) ! imaginary part, last columns
  enddo
  !$omp end parallel do

  ! Compute backward FFT
  Status = DftiComputeBackward( Desc_Handle_xy, ft)

  ! Arrange the data in a 2D array
  !$omp parallel do private(j)
  do j = 0, ny-1
     f(:,j) = ft( (nx+2)*j+1:(nx+2)*j+nx )
  enddo
  !$omp end parallel do

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
  use mkl_dfti
  use mkl_fft_descriptors
  implicit none
  integer, intent (in) :: L, n
  integer :: j
  real (kind=pr), dimension ( (L+2) * n ) :: ft
  real (kind=pr), dimension (0:L-1, 0:n-1), intent (in) ::  f
  real (kind=pr), dimension (0:L-1, 0:n-1), intent (out) ::  fk
  integer :: Status, strides_in(3)
  type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_1D

  ! Arrange the data in a 1D array
  !$omp parallel do private(j)
  do j = 0, n-1
     ft( (L+2)*j+1:(L+2)*j+L ) = f(:,j)
     ft( (L+2)*j+L+1:(L+2)*j+L+2 ) = 0.0
  enddo
  !$omp end parallel do

  ! Initialize FFT
  strides_in(1) = 0   ! Strides for a (L+2)*n array
  strides_in(2) = 1

  if ( precision(0.0) .eq. precision(0d0)) then
    Status = DftiCreateDescriptor( Desc_Handle_1D, DFTI_DOUBLE, DFTI_REAL, 1, L )
    Status = DftiSetValue( Desc_Handle_1D, DFTI_PRECISION, DFTI_DOUBLE) ! Set double precision
  else
    Status = DftiCreateDescriptor( Desc_Handle_1D, DFTI_SINGLE, DFTI_REAL, 1, L )
    Status = DftiSetValue( Desc_Handle_1D, DFTI_PRECISION, DFTI_SINGLE) ! Set single precision
  endif

  Status = DftiSetValue( Desc_Handle_1D, DFTI_FORWARD_SCALE, 1.0 / real(L)  ) ! Set forward scaling
  Status = DftiSetValue( Desc_Handle_1D, DFTI_BACKWARD_SCALE, 1.0 ) ! Set backward scaling
  Status = DftiSetValue( Desc_Handle_1D, DFTI_INPUT_STRIDES, strides_in) ! Set strides
  Status = DftiSetValue( Desc_Handle_1D, DFTI_NUMBER_OF_TRANSFORMS, n)
  Status = DftiSetValue( Desc_Handle_1D, DFTI_INPUT_DISTANCE, L+2)

  Status = DftiCommitDescriptor( Desc_Handle_1D )

  ! Compute forward FFT
  Status = DftiComputeForward( Desc_Handle_1D, ft)

  ! Arrange the data in a 2D array
  !$omp parallel do private(j)
  do j = 0, n-1
     fk(:,j) = ft( (L+2)*j+1:(L+2)*j+L )
  enddo
  !$omp end parallel do

  ! Deallocate 
  Status = DftiFreeDescriptor(Desc_Handle_1D)

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
  use mkl_dfti
  use mkl_fft_descriptors
  implicit none
  integer, intent (in) :: L, n
  integer :: j
  real (kind=pr), dimension ( (L+2) * n ) :: ft
  real (kind=pr), dimension (0:L-1, 0:n-1), intent (in) ::  fk
  real (kind=pr), dimension (0:L-1, 0:n-1), intent (out) ::  f
  integer :: Status, strides_in(3)
  type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_1D

  ! Arrange the data in a 1D array
  !$omp parallel do private(j)
  do j = 0, n-1
     ft( (L+2)*j+1:(L+2)*j+L ) = fk(:,j)
     ft( (L+2)*j+L+1:(L+2)*j+L+2 ) = 0.0
  enddo
  !$omp end parallel do

  ! Initialize FFT
  strides_in(1) = 0   ! Strides for a (L+2)*n array
  strides_in(2) = 1

  if ( precision(0.0) .eq. precision(0d0)) then
    Status = DftiCreateDescriptor( Desc_Handle_1D, DFTI_DOUBLE, DFTI_REAL, 1, L )
    Status = DftiSetValue( Desc_Handle_1D, DFTI_PRECISION, DFTI_DOUBLE) ! Set double precision
  else
    Status = DftiCreateDescriptor( Desc_Handle_1D, DFTI_SINGLE, DFTI_REAL, 1, L )
    Status = DftiSetValue( Desc_Handle_1D, DFTI_PRECISION, DFTI_SINGLE) ! Set single precision
  endif

  Status = DftiSetValue( Desc_Handle_1D, DFTI_FORWARD_SCALE, 1.0 / real(L)  ) ! Set forward scaling
  Status = DftiSetValue( Desc_Handle_1D, DFTI_BACKWARD_SCALE, 1.0 ) ! Set backward scaling
  Status = DftiSetValue( Desc_Handle_1D, DFTI_INPUT_STRIDES, strides_in) ! Set strides
  Status = DftiSetValue( Desc_Handle_1D, DFTI_NUMBER_OF_TRANSFORMS, n)
  Status = DftiSetValue( Desc_Handle_1D, DFTI_INPUT_DISTANCE, L+2)

  Status = DftiCommitDescriptor( Desc_Handle_1D )

  ! Compute forward FFT
  Status = DftiComputeBackward( Desc_Handle_1D, ft)

  ! Arrange the data in a 2D array
  !$omp parallel do private(j)
  do j = 0, n-1
     f(:,j) = ft( (L+2)*j+1:(L+2)*j+L )
  enddo
  !$omp end parallel do

  ! Deallocate 
  Status = DftiFreeDescriptor(Desc_Handle_1D)

  ! last mode M=L, L+1; mode KF left unconsidered => filtering
end subroutine cofits





