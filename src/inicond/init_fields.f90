subroutine init_fields (time, u, uk, pk, vor, nlk)
  use share_vars
  use FieldExport
  implicit none
  real(kind=pr), intent(out) :: time
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent (inout) :: u, uk, nlk
  real(kind=pr), dimension(0:nx-1,0:ny-1), intent (inout) :: vor, pk
  real(kind=pr), dimension(0:nx-1,0:ny-1) :: vortk
  real(kind=pr) :: r0,we,d,r1,r2, max_divergence, dummy1,dummy2,dummy3
  integer :: ix,iy
  time = 0.d0
  u   = 0.d0
  uk  = 0.d0
  vor = 0.d0
  pk  = 0.d0

  select case (inicond)
  !*****************
  case ('quiescent')
  !*****************
    u = 0.d0
    uk= 0.d0
    vor = 0.d0
    pk = 0.d0

  !*******************
  case ('lamballais')
  !*******************
    call lamballais(u,uk,pk,vor,nlk)

  !*****************
  case ('turbulent')
  !*****************
    call random_seed()
    do ix=0,nx-1
    do iy=0,ny-1
       call RANDOM_NUMBER(d)
       vor(ix,iy) = 200.d0*(2.0d0*d - 1.d0)
    enddo
    enddo

    call coftxy ( vor, vortk )
    call vorticity2velocity ( vortk, uk )
    call cofitxy( uk(:,:,1), u(:,:,1))
    call cofitxy( uk(:,:,2), u(:,:,2))

  !*****************
  case ('bckp')
  !*****************
    call Fetch_attributes_2d_openmp( 'vor_00071.h5', ix, iy, dummy1, dummy2, time, dummy3 )
    call read_flusi_hdf5_2d_openmp( 'vor_00071.h5', vor )
    call coftxy ( vor, vortk )
    call vorticity2velocity ( vortk, uk )
    call cofitxy( uk(:,:,1), u(:,:,1))
    call cofitxy( uk(:,:,2), u(:,:,2))

  !*****************
  case ('dipole')
  !*****************
    x0 = 0.5d0*xl
    y0 = 0.5d0*yl
    r0 = 0.1d0
    we = 299.528385375226d0
    d  = 0.1d0

    do ix=0,nx-1
    do iy=0,ny-1
       r1 = sqrt( (real(ix)*dx-x0)**2 + (real(iy)*dy-y0-d)**2 ) / r0
       r2 = sqrt( (real(ix)*dx-x0)**2 + (real(iy)*dy-y0+d)**2 ) / r0
       vor(ix,iy) = we * (1.d0-r1**2)*exp(-r1**2) - we * (1.d0-r2**2)*exp(-r2**2)
    enddo
    enddo

    call coftxy ( vor, vortk )
    call vorticity2velocity ( vortk, uk )
    call cofitxy( uk(:,:,1), u(:,:,1))
    call cofitxy( uk(:,:,2), u(:,:,2))
  end select


  !-----------------------------------------------------------------------------
  ! ensure the initial field is divergence-free
  !-----------------------------------------------------------------------------
  write(*,'("inicond=",A," max field divergence ",es12.4)') trim(inicond), max_divergence(uk)
  if ( ipressure == "classic") then
    ! there's a problem doing this with "modified"
    call add_pressure( uk, uk, u, vor )
  endif
  write(*,'("inicond=",A," max field divergence ",es12.4)') trim(inicond), max_divergence(uk)

end subroutine
