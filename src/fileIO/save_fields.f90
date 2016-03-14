subroutine SaveField( time, filename, field_out)
  use share_vars
  use hdf5
  implicit none

  ! The field to be written to disk:
  real(kind=pr),intent(in) :: field_out(0:nx-1,0:ny-1)
  real(kind=pr),intent(in) :: time
  character(len=*), intent (in) :: filename

  call write_flusi_hdf5_2d_openmp( time, filename, field_out)
end subroutine

! subroutine save_fields (n1, time, dt1, vortk, nlk, workvis, nbackup, beam, bpressure, ivideo, u, press, tau_beam_old)
! !======================================================================
! !     fields are stored in one file per time step
! !======================================================================
!   use share_vars
!   use fieldexport
!   use spectral_solver
!   use spectral_essentials
!   implicit none
!   integer :: ix, iy
!   integer, intent (inout) :: nbackup
!   integer, intent (in) :: n1,ivideo
!   real (kind=pr), intent (in) :: time, dt1
!   real (kind=pr), dimension (0:ns-1), intent(in) :: bpressure, tau_beam_old
!   real (kind=pr), dimension (0:nx-1, 0:ny-1, 0:1), intent (in)  :: nlk, vortk, workvis
!   real (kind=pr), dimension (0:nx-1, 0:ny-1, 1:2), intent (in) :: u
!   real (kind=pr), dimension (0:nx-1,0:ny-1) ::  vort, work1, work2, work3
!   real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) :: press
!   real (kind=pr), dimension (0:ns-1, 1:6), intent (in) :: beam
!   character (len=11) :: name
!   character (len=1) :: name1
!   real (kind=pr) :: Mean_ux, Mean_uy
! !--Set up file name base
!   write (name, '(es10.4)') time
!
!
!   !=================================================================================
!   !--Save pressure
!   if (iSavePress > 0) then
!   call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'p_'//name, work3, iSavePress, xl,yl, "pressure")
!   endif
!   !=================================================================================
!   !--Save velocity and streamfunction
!   if (iSaveVel>0) then
!   call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'ux_'//name, u(:,:,1), iSaveVel, xl,yl, "x-velocity")
!   call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'uy_'//name, u(:,:,2), iSaveVel, xl,yl, "y-velocity")
!   endif
!
!   if (iSaveSTR>0) then
!   call poisson (vortk, vort) ! vort = streamfunction for the moment
!   call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'str_'//name, vort, iSaveSTR, xl,yl, "streamfunction")
!   endif
!   !=================================================================================
!   !--Save vorticity
!   if (iSaveVort > 0) then
!   call cofitxy (vortk, vort)
!   call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'vor_'//name, vort, iSaveVort, xl,yl,"vorticity")
!   endif
!   !=================================================================================
!   !save mask
!   if (iSaveMask > 0) then
!   call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'mask_'//name, mask, iSaveMask, xl,yl,"mask")
!   endif
!
!   !=================================================================================
!   if (iSaveMaskVel >0) then
!   call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'mask_vy_'//name, maskvy, iSaveMaskVel, xl,yl, "mask_vx")
!   call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'mask_vx_'//name, maskvx, iSaveMaskVel, xl,yl, "mask_vy")
!   endif
!
! !=================================================================================
!    if (iSaveStress>0) then
!        call coftxy  ( u(:,:,1), vort  )
!        call cofdx   ( vort    , work1 )
!        call cofitxy ( work1   , work3 )
!
!        work3 = 2.0 * nu * work3
!        call SaveField( './fields/'//trim(simulation_name)//'stress_a_'//name, work3, iSaveStress, xl,yl,"viscous-stress-a")
!        !------------------------
!        call coftxy  ( u(:,:,1), vort  )
!        call cofdy   ( vort    , work1 )
!        call cofitxy ( work1   , work2 )
!
!        call coftxy  ( u(:,:,2), vort  )
!        call cofdx   ( vort    , work1 )
!        call cofitxy ( work1   , work3 )
!
!        work3 = nu * ( work2 + work3 )
!        call SaveField( './fields/'//trim(simulation_name)//'stress_b_'//name, work3, iSaveStress, xl,yl,"viscous-stress-b")
!        !----------finite differences
!        do ix=1,nx-2
!        do iy=1,ny-2
! 	work3(ix,iy)=(u(ix+1,iy,1)-u(ix-1,iy,1))/2.0/dx
!        enddo
!        enddo
!        work3 = 2.0 * nu * work3
!        call SaveField( './fields/'//trim(simulation_name)//'stress_a_FD_'//name, work3, iSaveStress, xl,yl,"viscous-stress-a")
!
!        do ix=1,nx-2
!        do iy=1,ny-2
! 	work2(ix,iy)=(u(ix,iy+1,1)-u(ix,iy-1,1))/2.0/dy
! 	work3(ix,iy)=(u(ix+1,iy,2)-u(ix-1,iy,2))/2.0/dx
!        enddo
!        enddo
!        work3 = nu * ( work2 + work3 )
!        call SaveField( './fields/'//trim(simulation_name)//'stress_b_FD_'//name, work3, iSaveStress, xl,yl,"viscous-stress-b")
!
!    endif
!
!
! !=================================================================================
! !--Backup data
!   if (nbackup == 2) then
!      nbackup = 0 ! do not backup this time
!   else
!     call MakeRuntimeBackup(n1, time, dt1, vortk, nlk, workvis, nbackup, beam, bpressure, ivideo, u, tau_beam_old)
!   endif
!
! end subroutine save_fields
!
! !=====================================================================================================================================
!
! subroutine MakeRuntimeBackup(n1, time, dt1, vortk, nlk, workvis, nbackup, beam,bpressure,ivideo, u, tau_beam_old)
!   use share_vars
!   implicit none
!   integer, 					   intent (inout) :: nbackup
!   integer, 					   intent (in) :: n1,ivideo
!   real (kind=pr),                                  intent (in) :: time, dt1
!   real (kind=pr), dimension (0:ns-1),              intent (in) :: bpressure, tau_beam_old
!   real (kind=pr), dimension (0:nx-1, 0:ny-1, 0:1), intent (in) :: nlk, vortk, workvis, u
!   real (kind=pr), dimension (0:ns-1, 1:6),         intent (in) :: beam
!   character (len=1) :: name1
!
!   write (33,*) '*** Making a backup, time=',time
!
!   write (name1, '(I1)') nbackup
!   open (15, file = trim(dir_name)//'/runtime_backup'//name1//'.in', form='unformatted', status='replace')
!   write (15) time
!   write (15) n1, dt1, vortk, nlk, workvis, mask, maskvx, maskvy, beam, bpressure, ivideo, u, tau_beam_old
!   close (15)
!   nbackup = 1 - nbackup
!
!   write (33,*) '*** excecuting: tar czf backup_'//name1//'.tar.gz '//trim(simulation_name)//'*'
!   call system('tar czf backup_'//name1//'.tar.gz '//trim(simulation_name)//'*')
!
!
!   write (33,*) '*** information: backup done'
!
! end subroutine MakeRuntimeBackup
!
