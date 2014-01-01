module parameters

 implicit none
 contains

! Wrapper to read parameters from an ini file for fsi.  Also reads
! parameters which are set in the vars module.
subroutine get_params(paramsfile) 
  use share_vars

  character (len=*) :: paramsfile  ! The file we read the PARAMS from
  integer :: i  
  character PARAMS(nlines)*256 ! this array will contain the ascii-params file
  logical :: exist1
  
  ! check if the specified file exists
  inquire ( file=paramsfile, exist=exist1 )
  
  if ((paramsfile=="").or.(exist1.eqv..false.)) then
    write(*,*) "PARAMS file not found!"
    stop
  endif
  
  ! Read the paramsfile and put the length i and the text in PARAMS
  call read_params_file(PARAMS,i,paramsfile,.true.)

  ! Get parameter values from PARAMS
  call get_params_common(PARAMS,i)
end subroutine get_params


! Read the file paramsfile, count the lines (output in i) and put the
! text in PARAMS.
subroutine read_params_file(PARAMS,i,paramsfile, verbose)
  use share_vars
  implicit none
  
  integer,intent(out) :: i
  integer :: io_error
  ! This array will contain the ascii-params file
  character,intent(inout) :: PARAMS(nlines)*256
  character(len=*) :: paramsfile ! this is the file we read the PARAMS from
  logical, intent(in) :: verbose

  ! Read in the params file (root only)
  io_error=0
  if (verbose) then
    write (*,*) "*************************************************"
    write (*,*) "*** info: reading params from ",trim(paramsfile)
    write (*,*) "*************************************************"
  endif
  i = 1
  open(unit=14,file=trim(adjustl(paramsfile)),action='read',status='old')    
  do while ((io_error==0).and.(i<=nlines))
      read (14,'(A)',iostat=io_error) PARAMS(i)  
      i = i+1
  enddo
  close (14)
  i = i-1 ! counted one too far
  
end subroutine read_params_file


! Read individual parameter values from the PARAMS string for the vars
! module.
subroutine get_params_common(PARAMS,i)
  use share_vars
  implicit none

  integer,intent(in) :: i
  character,intent(in) :: PARAMS(nlines)*256 ! Contains the ascii-params file

  ! Resolution section
  call GetValue_Int(PARAMS,i,"Resolution","nx",nx, 4)
  call GetValue_Int(PARAMS,i,"Resolution","ny",ny, 4)
  call GetValue_Int(PARAMS,i,"Time","nt",nt, 9999999)
  call GetValue_Real(PARAMS,i,"Time","Tmax",Tmax,1.d9)
  call GetValue_Real(PARAMS,i,"Time","CFL",cfl,0.1d0)
  call GetValue_String(PARAMS,i,"Time","iMethod",iMethod,"RK2")
  call GetValue_Real(PARAMS,i,"ReynoldsNumber","nu",nu,1.d-2)  

  ! Initial conditions section
  call GetValue_String(PARAMS,i,"InitialCondition","inicond",inicond, "none")

  ! Dealasing section
  call GetValue_Int(PARAMS,i,"Dealiasing","iDealias",iDealias, 1)

  ! Penalization section
  call GetValue_String(PARAMS,i,"Penalization","iMask",iMask, "none")
  call GetValue_Real(PARAMS,i,"Penalization","eps",eps, 1.d-2)

  ! Geometry section
  call GetValue_Real(PARAMS,i,"Geometry","xl",xl, 1.d0)
  call GetValue_Real(PARAMS,i,"Geometry","yl",yl, 1.d0)

  ! saving section
  call GetValue_Real(PARAMS,i,"Saving","tsave",tsave, 1.d0)
  call GetValue_Real(PARAMS,i,"Saving","tdrag",tdrag, 1.d0)

  ! mean flow
  call GetValue_String(PARAMS,i,"MeanFlow","iMeanFlow",iMeanFlow, "none")
  call GetValue_Real(PARAMS,i,"MeanFlow","ux_mean",ux_mean, 0.d0)
  call GetValue_Real(PARAMS,i,"MeanFlow","uy_mean",uy_mean, 0.d0)
  
  ! Set other parameters (all procs)
  pi=4.d0 *datan(1.d0)
  ! lattice spacing is global
  dx=xl/dble(nx)
  dy=yl/dble(ny)
end subroutine get_params_common


!-------------------------------------------------------------------------------
! Fetches a REAL VALUED parameter from the PARAMS.ini file.
! Displays what it does on stdout (so you can see whats going on)
! Input:
!       PARAMS: the complete *.ini file as array of characters*256
!               with maximum nlines=2048 lines (this value is set in share_vars)
!       actual_lines: number of lines in the params file
!                     because we don't need to loop over empty lines
!       section: the section we're looking for
!       keyword: the keyword we're looking for
!       defaultvalue: if the we can't find the parameter, we return this and warn
! Output:
!       params_real: this is the parameter you were looking for
subroutine GetValue_real (PARAMS, actual_lines, section, keyword, params_real,defaultvalue)
  use share_vars
  implicit none
  character section*(*) ! What section do you look for? for example [Resolution]
  character keyword*(*)   ! what keyword do you look for? for example nx=128
  character (len=80)  value    ! returns the value
  character PARAMS(nlines)*256  ! this is the complete PARAMS.ini file
  real (kind=pr) :: params_real, defaultvalue 
  integer actual_lines

  ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
  call GetValue(PARAMS, actual_lines, section, keyword, value)
  if (value .ne. '') then
    read (value, *) params_real    
    write (value,'(g10.3)') params_real
    write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
  else
    write (value,'(g10.3)') defaultvalue
    write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))//&
          " (THIS IS THE DEFAULT VALUE!)"
    params_real = defaultvalue
  endif    
end subroutine GetValue_real




!-------------------------------------------------------------------------------
! Fetches a STRING VALUED parameter from the PARAMS.ini file.
! Displays what it does on stdout (so you can see whats going on)
! Input:
!       PARAMS: the complete *.ini file as array of characters*256
!               with maximum nlines=2048 lines (this value is set in share_vars)
!       actual_lines: number of lines in the params file
!                     because we don't need to loop over empty lines
!       section: the section we're looking for
!       keyword: the keyword we're looking for
!       defaultvalue: if the we can't find the parameter, we return this and warn
! Output:
!       params_string: this is the parameter you were looking for
subroutine GetValue_string (PARAMS, actual_lines, section, keyword, params_string, defaultvalue)
  use share_vars
  implicit none
  
  character section*(*) ! what section do you look for? for example [Resolution]
  character keyword*(*)   ! what keyword do you look for? for example nx=128
  character (len=80)  value    ! returns the value
  character PARAMS(nlines)*256  ! this is the complete PARAMS.ini file
  character (len=*), intent (inout) :: params_string
  character (len=*), intent (in) :: defaultvalue 
  integer actual_lines

  ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
  call GetValue(PARAMS, actual_lines, section, keyword, value)
  if (value .ne. '') then
    params_string = value
    ! its a bit dirty but it avoids filling the screen with "nothing" anytime we check
    ! the runtime control file
    if (keyword.ne."runtime_control") then 
    write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
    endif
  else
    write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))//&
          " (THIS IS THE DEFAULT VALUE!)"
    params_string = defaultvalue
  endif
end subroutine GetValue_string



!-------------------------------------------------------------------------------
! Fetches a VECTOR VALUED parameter from the PARAMS.ini file.
! Displays what it does on stdout (so you can see whats going on)
! Input:
!       PARAMS: the complete *.ini file as array of characters*256
!               with maximum nlines=2048 lines (this value is set in share_vars)
!       actual_lines: number of lines in the params file
!                     because we don't need to loop over empty lines
!       section: the section we're looking for
!       keyword: the keyword we're looking for
!       defaultvalue: if the we can't find a vector, we return this and warn
! Output:
!       params_vector: this is the parameter you were looking for
subroutine GetValue_vector (PARAMS, actual_lines, section, keyword, params_vector, &
     defaultvalue)
  use share_vars
  implicit none
  character section*(*) ! What section do you look for? for example [Resolution]
  character keyword*(*)   ! what keyword do you look for? for example nx=128
  character (len=80)  value    ! returns the value
  character PARAMS(nlines)*256  ! this is the complete PARAMS.ini file
  real (kind=pr) :: params_vector(1:3), defaultvalue(1:3)
  integer actual_lines

  ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
  call GetValue(PARAMS, actual_lines, section, keyword, value)
  if (value .ne. '') then
    ! read the three values from the vector string
    read (value, *) params_vector
    write (value,'(3(g10.3,1x))') params_vector
    write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
  else
    write (value,'(3(g10.3,1x))') defaultvalue
    write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))//&
          " (THIS IS THE DEFAULT VALUE!)"
    params_vector = defaultvalue
  endif
end subroutine GetValue_vector



!-------------------------------------------------------------------------------
! Fetches a INTEGER VALUED parameter from the PARAMS.ini file.
! Displays what it does on stdout (so you can see whats going on)
! Input:
!       PARAMS: the complete *.ini file as array of characters*256
!               with maximum nlines=2048 lines (this value is set in share_vars)
!       actual_lines: number of lines in the params file
!                     because we don't need to loop over empty lines
!       section: the section we're looking for
!       keyword: the keyword we're looking for
!       defaultvalue: if the we can't find the parameter, we return this and warn
! Output:
!       params_int: this is the parameter you were looking for
subroutine GetValue_Int(PARAMS, actual_lines, section, keyword, params_int,&
     defaultvalue)
  use share_vars
  implicit none

  character section*(*) ! What section do you look for? for example [Resolution]
  character keyword*(*)   ! what keyword do you look for? for example nx=128
  character (len=80)  value    ! returns the value
  character PARAMS(nlines)*256  ! this is the complete PARAMS.ini file
  integer params_int, actual_lines, defaultvalue

  !------------------
  ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
  !------------------
  call GetValue(PARAMS, actual_lines, section, keyword, value)
  if (value .ne. '') then
    read (value, *) params_int
    write (value,'(i7)') params_int
    write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
  else
    write (value,'(i7)') defaultvalue
    write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))//&
          " (THIS IS THE DEFAULT VALUE!)"
    params_int = defaultvalue
  endif
end subroutine GetValue_Int





!-------------------------------------------------------------------------------
! Extracts a value from the PARAMS.ini file, which is in "section" and
! which is named "keyword"
! Input:
!       PARAMS: the complete *.ini file as array of characters*256
!               with maximum nlines=2048 lines (this value is set in share_vars)
!       actual_lines: number of lines in the params file
!                     because we don't need to loop over empty lines
!       section: the section we're looking for
!       keyword: the keyword we're looking for
! Output:
!       value: is a string contaiing everything between '=' and ';'
!              to be processed further, depending on the expected type
!              of variable (e.g. you read an integer from this string)
subroutine GetValue (PARAMS, actual_lines, section, keyword, value)
  use share_vars
  implicit none

  character section*(*) ! What section do you look for? for example [Resolution]
  character keyword*(*)   ! what keyword do you look for? for example nx=128
  character value*(*)   ! returns the value
  character PARAMS(nlines)*256  ! this is the complete PARAMS.ini file
  integer actual_lines   ! how many lines did you actually read?  
  integer :: maxline = 256  ! how many characters per line?
  integer i, j,k
  logical foundsection

  foundsection = .false.
  value = ''

  !------------------------------------------------------------------
  do i=1, actual_lines     ! loop over the lines of PARAMS.ini file
     if ((PARAMS(i)(1:1).ne.'#').and.&
          (PARAMS(i)(1:1).ne.';').and.&
          (PARAMS(i)(1:1).ne.'!')) then   ! ignore commented lines compleetly


        if (PARAMS(i)(1:1) == '[') then   ! the first char would have to be '['
           do j = 2, maxline    ! then we look fot the corrresponding ']'
              if (PARAMS(i)(j:j) == ']') then  ! we found it
                 if (section == PARAMS(i)(2:j-1)) then ! is this the section we"re looking for?
                    foundsection = .true.   ! yes, it is
                    exit
                 endif
              endif
           enddo
        else      
           if (foundsection .eqv. .true.) then  ! yes we found the section, now we're looking for the keyword
              do j=1, maxline    ! scan the line
                 if (PARAMS(i)(j:j) == '=') then  ! found the '='
                    if (keyword == PARAMS(i)(1:j-1)) then ! is this the keyword you're looking for?
                       do k = j+1, maxline   ! everything behind the '=' and before ';' is the value
                          if (PARAMS(i)(k:k) == ';') then  ! found the delimiter
                             value = PARAMS(i)(j+1:k-1)  ! value is between '=', and ';'   
                             exit
                          endif
                       enddo
                       if ((value == '')) then
                          write (*,'(A)') "??? Though found the keyword, I'm unable to find value for variable --> "& 
                               //trim(keyword)//" <-- maybe missing delimiter (;)?"    
                       endif
                    endif
                 endif
              enddo
           endif
        endif

     endif
  enddo
end subroutine getvalue

end module parameters