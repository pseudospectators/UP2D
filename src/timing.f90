module timing_module
  !-------------------------------------------------------
  implicit none
  !-------------------------------------------------------

  contains

  ! returns the time, in seconds, since an arbitrary time in the past.
  ! NOTE: compatible with the MPI routine MPI_wtime()
  function MPI_wtime()
    implicit none
    double precision :: MPI_wtime
    integer*8 :: time, rate, dummy

    call system_clock( time, rate, dummy )
    MPI_wtime = dble(time) / rate
  end function

  !----------------------------------------------------------------------------
end module timing_module
