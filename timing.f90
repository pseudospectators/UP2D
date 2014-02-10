module timing
  !-------------------------------------------------------
  implicit none
  save 
  integer, parameter                        :: N_avg =100 ! each measurement is average of N_avg individual measurements
  integer, parameter                        :: N_chans = 20 ! how many channels?
  real, dimension (1:N_chans, 1:N_avg)      :: time_measurements = 0.0 !contains N_avg values for each channel (for averaging)
  integer*8, dimension (1:N_chans)          :: start_times = 0  !contains the time when the measurement started, for each channel
  integer, dimension (1:N_chans)            :: counters = 1  !counters for the measurements

  !-------------------------------------------------------
  
  contains

  function performance(Start, channel)
    implicit none
    integer, intent (in) :: channel
    integer*8 :: rate, dummy, time_current
    real :: performance
    character (len=*), intent(in) :: Start
    performance = 0.0
    if (Start == "start") then
      !-- save start time to the "start_times"-array (there's only one value per channel)
      call system_clock(start_times(channel), rate, dummy) 
      !-- return dummy value
      performance = 0.0
      
    elseif (Start == "stop") then
      !-- get current time    
      call system_clock(time_current, rate, dummy) 
      
      !-- difference between stop and start in seconds
      time_measurements(channel,counters(channel)) =  real(time_current-start_times(channel)) / real(rate)
      
      !-- now return the avg value
      performance = sum(time_measurements(channel,:)) / real(N_avg)
      
      !-- increase counter
      counters(channel) = counters(channel) + 1
      !-- when called N_avg start_times, start over and overwrite the first
      if (counters(channel)>N_avg) counters(channel) = 1
      
    endif
    
  end function performance

  
  !----------------------------------------------------------------------------
end module timing

