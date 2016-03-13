module fieldexport
  implicit none
  contains

! -------------------------------------------------------------------------------------------------------------------------------------------------
! 		Field Operations Module
! 
!  ******************************************
!  * Import / Export routines for 2D fields *
!  ******************************************
!   
!   Export Routines:
!   ----------------
!     - SavePPM: Saves a field with dimensions nx,ny to a PPM image file by using a specified palette and scaling
!     - SavePNG: Same as SavePPM, but converts the image file to PNG
!     - SaveGIF: Directly writes a compressed GIF-imagefile. to be preferred!!!!!!!!!!!!!!
!     - ConvertFieldToRGB: Generates a RGB array with the color values 0..255 out of the data array by using a specified palette and scaling.
!     - FargePalette: Fills an array with the color tables used in ConvertFieldToRGB, where:
! 	01: Farge Table blue-yellow-red
! 	02: Farge Table gold-red-blue (pressure)
! 	03: Farge Table gold-turquoise-purple
! 	04: Farge Table gold-turquoise-red (sharp)
! 	05: Farge Table blue-yellow-red
! 	06: Farge Table gold-purple-turquoise
! 	07: Farge Table gold-green-blue
! 	08: Farge Table pink-yellow-light green
! 	09: Farge Table red-b/w-green
! 	10: Farge Table blue-yellow-red
! 	11: Matlab Spring Table
! 	12: Flameball Table
! 	13: Matlab Standard (Jet) Table
! 	14: Matlab Standard (Jet) Table, contours for pressure
! 	
!   Import Routines:
!   ----------------
!     - LoadField: Loads a field from a specified file. There are two file formats: one contains only the data, the second one also the resultion.
!       The latter can be obtained using GetSize. The array must be allocated, LoadField gets an assumed shaped array.
!     - GetSize: extracts the resolution information from a file, if present. If not, it returns a flag indicating that the file has an old format without
!       resolution in the file.
!     - ResizeField: Sometimes, it may be nessesairy to load a field from a file with a different resolution, for example when computing a steady state 
!       flow. In this case, it is very convenient to take the final result from a lower resolution as starting point for higher resolutions. This routine
!       resizes a field to an arbitrary new resolution, defined by the target array itself.
!       

  
subroutine ResizeAvg(field_source, field_target)
  ! this routine is the simplest way to upsample a field to a higher resolution.
  ! attention, it only works for doubling the resolution (e.g. from 100x100 to 200x200)
  ! we copy all values from the old field to the coinciding points on the new one
  ! missing points are initialized as the average of their neighbors.
  ! start with "diagonal" points, where the neighbors on the four corners are the values of the old field
  ! then go row-wise in x-direction (note everything has to be periodic)
  ! then go column-wise in y-direction.
  !
  ! It has been devolped to solve LAPLACE on multigrids. 
  ! NOTE: its essentially different from PolynomUpsampling!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! the latter must use ODD grids, 
  use share_vars
  implicit none
  integer :: ix,iy,nx1,nx2,ny1,ny2
  real (kind=pr), intent (in )    :: field_source (0:,0:)
  real (kind=pr), intent (out)    :: field_target (0:,0:)
  
  nx1=size(field_source,1)
  ny1=size(field_source,2)  
  nx2=size(field_target,1)
  ny2=size(field_target,2)   
  
  field_target = 0.0 
  
  ! fill matching points, the points where the finer grid matches the coarser one
  field_target(0:nx2-1:2,0:ny2-1:2) = field_source  
  
  !$omp parallel do private(iy, ix)
  do ix=1,nx2-1,2
  do iy=1,ny2-1,2
  field_target(ix,iy) = ( field_target( getindex(ix+1,nx2),getindex(iy+1,ny2)) &
    + field_target(getindex(ix-1,nx2),getindex(iy-1,ny2)) &
    + field_target(getindex(ix-1,nx2),getindex(iy+1,ny2)) &
    + field_target(getindex(ix+1,nx2),getindex(iy-1,ny2)) )/4.0
  enddo
  enddo
  !$omp end parallel do
  
  !$omp parallel do private(iy, ix)
  do ix=0,nx2-1,2
  do iy=1,ny2-1,2
  field_target(ix,iy) = ( field_target( getindex(ix+1,nx2),getindex(iy,ny2)) &
   + field_target(getindex(ix-1,nx2),getindex(iy,ny2)) &
   + field_target(getindex(ix,nx2),getindex(iy+1,ny2)) &
   + field_target(getindex(ix,nx2),getindex(iy-1,ny2)) )/4.0
  enddo
  enddo
  !$omp end parallel do

  !$omp parallel do private(iy, ix)
  do ix=1,nx2-1,2
  do iy=0,ny2-1,2
  field_target(ix,iy) = ( field_target( getindex(ix+1,nx2),getindex(iy,ny2)) &
  + field_target(getindex(ix-1,nx2),getindex(iy,ny2)) &
  + field_target(getindex(ix,nx2),getindex(iy+1,ny2)) &
  + field_target(getindex(ix,nx2),getindex(iy-1,ny2)) )/4.0
  enddo
  enddo
  !$omp end parallel do
  
  
end subroutine


integer pure function GetIndex(ix,nx)
  implicit none
  integer, intent (in) :: ix,nx
  integer :: tmp
  tmp=ix
  if (tmp<0) tmp = tmp+nx
  if (tmp>nx-1) tmp = tmp-nx
  GetIndex=tmp
  return
end function GetIndex




  
  subroutine PolynomUpsampling(field_source,field_target)
    use share_vars
    implicit none
    real (kind=pr), intent (in )    :: field_source (0:,0:)
    real (kind=pr), intent (out)    :: field_target (0:,0:)
    real (kind=pr), dimension (:,:), allocatable :: field_source_k, field_target_k
    integer :: nx1,nx2,ny1,ny2,ix,iy,i,j, nx_tmp,ny_tmp, ix2,iy2,R,RR,L,LL,K
    
    write (*,*) "source:", shape(field_source)
    write (*,*) "target:", shape(field_target)
  
    nx1=size(field_source,1)
    ny1=size(field_source,2)    
    nx2=size(field_target,1)
    ny2=size(field_target,2)   
    
    field_target = 0.0
    ! fill matching points, the points where the finer grid matches the coarser one
    field_target(0:nx2-1:2,0:ny2-1:2) = field_source
    
    !+===========================================================================================
    ! PolynomUpsampling
    !===========================================================================================
    ! then we do th elines
    do iy=0,ny2-1, 2 ! gerade zahlen entlang y( hier gibts schon wert))
    do ix=1,nx2-2, 2 ! ungrade zahlen entlang x (zwischen den org gitter))
      if ((ix>=3).and.(ix<=nx2-1-3)) then
      K   = int( 0.5*real(ix) ) ! k
      field_target(ix,iy)=(9.d0*(field_source(K,iy/2)+field_source(K+1,iy/2)) & 
                          -field_source(K-1,iy/2)-field_source(K+2,iy/2)   )/16.d0
      elseif (ix==1) then 
      field_target(ix,iy)=(5.d0*field_source(0,iy/2)+15.d0*field_source(1,iy/2)&
                          -5.d0*field_source(2,iy/2)+ 1.d0*field_source(3,iy/2) )/16.d0
      elseif (ix==nx2-2) then 
       field_target(ix,iy)=(5.d0*field_source(nx1-4,iy/2)+15.d0*field_source(nx1-3,iy/2)&
                          - 5.d0*field_source(nx1-2,iy/2)+ 1.d0*field_source(nx1-1,iy/2) )/16.d0                           
      endif
    enddo
    enddo

    ! then we do th elines
    do ix=0,nx2-1, 2 ! gerade zahlen entlang y( hier gibts schon wert))
    do iy=1,ny2-2, 2 ! ungrade zahlen entlang x (zwischen den org gitter))
      if ((iy>=3).and.(iy<+ny2-1-3)) then
      K   = int( 0.5*real(iy) ) ! k
      field_target(ix,iy)=(9.d0*(field_source(ix/2,K)+field_source(ix/2,K+1)) & 
                          -field_source(ix/2,K-1)-field_source(ix/2,K+2)   )/16.d0
      elseif (iy==1) then 
      field_target(ix,iy)=(5.d0*field_source(ix/2,0)+15.d0*field_source(ix/2,1)&
                          -5.d0*field_source(ix/2,2)+1.d0*field_source(ix/2,3))/16.d0
      elseif (iy==ny2-2) then 
       field_target(ix,iy)=(+5.d0*field_source(ix/2,ny1-1)+15.d0*field_source(ix/2,ny1-2)&
                            -5.d0*field_source(ix/2,ny1-3)+ 1.d0*field_source(ix/2,ny1-4))/16.d0                             
      endif
    enddo
    enddo
    
    
    do iy=1,ny2-2, 2 
    do ix=1,nx2-2, 2 
      if ((ix>=3).and.(ix<+nx2-1-3)) then
      field_target(ix,iy)=(9.d0*(field_target(ix-1,iy)+field_target(ix+1,iy)) & 
                          -field_target(ix-3,iy)-field_target(ix+3,iy)   )/16.d0
      endif
    enddo
    enddo
    
      

    !=========================================================================================
    ! linear
    !=========================================================================================
    
!    ! then we do th elines
!    do iy=0,ny2-1, 2 ! gerade zahlen entlang y( hier gibts schon wert))
!    do ix=1,nx2-2, 2 ! ungrade zahlen entlang x (zwischen den org gitter))
!      L=ix/2
!      LL = L-1 !left left
!      R  = L+1 ! right
!      RR = L+2 !right right
!      field_target(ix,iy)=0.5d0*(field_source( L, iy/2 )+field_source( R, iy/2 ) )
!    enddo
!    enddo

!    do ix=0,nx2-1, 2 ! gerade zahlen entlang x( hier gibts schon wert))
!    do iy=1,ny2-2, 2 ! ungrade zahlen entlang y (zwischen den org gitter))
!      L=iy/2
!      LL = L-1 !left left
!      R  = L+1 ! right
!      RR = L+2 !right right
!      field_target(ix,iy)=0.5d0*(field_source( ix/2,L )+field_source( ix/2,R ) )
!    enddo
!    enddo
    
!    ! points that are in the center of the original grid cubes
!    do ix=1,nx2-2,2
!    do iy=1,ny2-2,2
!      field_target(ix,iy) = 0.5d0*( field_target(ix,iy-1) + field_target(ix,iy+1) )
!    enddo
!    enddo

    
    
    
end subroutine PolynomUpsampling

!====================================================================================================

subroutine SaveGIF(field, filename, iColorTable, color_min, color_max)
    !-------------------------------------------------------------
    ! Save gif subroutine. What is the essential difference? 
    ! the write_gif module works with a array for the color table, ColorMap.
    ! this contains the 256 possible colors in a gif. In the palettes (created in FargePalette)
    ! we create exactly this amount of different colors - perfect! now in the gif, we have a color table
    ! with the RGB values of all the 256 colors, and in the field we are left with an array of indices
    ! that indicate which color the pixel gets. 
    !
    ! IMPROVED 08 apr 2012: now indepedent of nx,ny, can save fields of any size
    ! IMPROVED: arguments iColorTable, color_min,color_max are now OPTIONAL
    !
    !-------------------------------------------------------------
    
    use share_vars
    use gif_util
    implicit none

    character*(*) 				:: filename
    real (kind=pr), intent (in) 		:: field (:,:)
    real (kind=pr), intent (in), optional 	:: color_min,color_max
    integer, intent (in), optional 		:: iColorTable
    real (kind=pr) 				:: minf,maxf
    integer, dimension(1:nPalettes,1:256*3) 	:: palette  !array for color values for FARGE palettes
    integer, dimension(1:3,0:255) 		:: ColorMap
    integer, allocatable 			:: Pixel(:,:) ! RGB data array    
    integer i, j, iPalette, intensity,nx1,ny1
    
    nx1=size(field,1)
    ny1=size(field,2)
    
    allocate (Pixel(1:nx1,1:ny1))
     
    !--------------------------------------------
    ! optional parameters
    !--------------------------------------------
    if ( present(color_min) ) then    
      minf=color_min
    else
      minf=minval(field)
    end if
    
    if ( present(color_max) ) then    
      maxf=color_max
    else
      maxf=maxval(field)
    end if    
    
    if ( present(iColorTable) ) then    
      iPalette=iColorTable
    else
      iPalette=13
    end if
    
    !create palettes
    call FargePalette(palette)
    ! fill colormap, the 256 colors for the gif.
    do i=1,256
      ColorMap(1,i-1) = palette(iPalette,i) !Red
      ColorMap(2,i-1) = palette(iPalette,i+256)!green
      ColorMap(3,i-1) = palette(iPalette,i+512)!blue
    enddo
    !loop over the field, intensity is the index of the color
    do i = 1, nx1
    do j = 1, ny1
      !intensity is the index of the color place in the palette array
      intensity = 1+nint(256.0*(field(i,1+ny1-j)-minf)/(maxf-minf)) !remember that we use inverse ordering in the code
      intensity = 1+nint(255.0*(field(i,1+ny1-j)-minf)/(maxf-minf)) !remember that we use inverse ordering in the code
      if (intensity<1) intensity=1 
      if (intensity>256) intensity=256
      Pixel(i,j)=intensity-1
    enddo
    enddo

    call writegif (trim(filename)//".gif", Pixel, ColorMap)
end subroutine SaveGIF

!==================================================================================================================================================
!==================================================================================================================================================

subroutine LoadField(field, filename)
!----------------------------------------------------------------
! Load a raw data field. in the new code versions, such a file contains a massive amount of comments (the PARAMS file used
! in the simulation). The beginning of the file is indicated by the string "%format1"
! then follows the resolution nx, ny, the domain size xl, yl and then three comments. the one in the middle
! is a string indicating the type of data, e.g. "% vorticity"
!
!	summary file format:
! ==========================================================
!     .
!     .
!  COMMENTS
!     .
!     .
!  %format1
!  nx
!  ny
!  xl
!  yl
!  COMMENT
!  % vorticity  (for example, there may be others)
!     .
!     .
!   DATA
!     .
!     .
! ==========================================================
!----------------------------------------------------------------
    use share_vars
    implicit none
    integer i,j,nx1,ny1,nx1_file,ny1_file
    character*(*), intent(in) :: filename
    character(len=256) :: dummy !to read from the file
    real (kind=pr), intent (inout) :: field (:,:)
    real(kind=pr)::xl_file, yl_file
    
    nx1=size(field,1)
    ny1=size(field,2)
    
    i=1
    open(unit=2,file=filename,status='old')
    do while ((dummy.ne."%format1").and.(i<150))
      read (2,'(A)') dummy
      i=i+1
    end do
    if (dummy.ne."%format1") then
      close(2)
      Write (*,*) "Unfortunately I did not find the string 'format1' in the file. This means it does NOT contain the resolution."
      stop
    endif

    read (2,*) nx1_file
    read (2,*) ny1_file
    if ((nx1_file.ne.nx1).or.(ny1_file.ne.ny1)) then
    write (*,*) "??? file resolution does not match the field you gave my to read in..."
    write (*,'("file: ",i4,"x",i4," field: ",i4,"x",i4)') nx1_file,ny1_file, nx1,ny1
    stop
    endif
    read (2,*) xl_file
    read (2,*) yl_file
    
    if ((xl_file.ne.xl).or.(yl_file.ne.yl)) then
      write (*,*) "domain sizes in file and PARAMS.ini differ"
      stop
    endif
    
    read (2,'(A)') dummy
    read (2,'(A)') dummy
    do j = ny1, 1, -1
      do i = 1, nx1
	read (2,*) field(i,j)
      enddo
    enddo

    close(2)  
end subroutine LoadField

!==================================================================================================================================================

subroutine GetSize(filename, FileContainsResolution)
    use share_vars
    implicit none
    character*(*) :: filename
    character(len=256) :: dummy !to read from the file
    integer :: i=1
    logical 	     :: FileContainsResolution
    
    open(unit=2,file=filename,status='old')
    FileContainsResolution=.true.
    
    do while ((dummy.ne."%format1").and.(i<150))
      read (2,'(A)') dummy
      i=i+1
    end do
    
    if (dummy.ne."%format1") then
      close(2)
      Write (*,*) "Unfortunately I did not find the string 'format1' in the file. This means it does NOT contain the resolution."
      FileContainsResolution=.false. 
    else
      read (2,*) nx
      read (2,*) ny
      read (2,*) xl
      read (2,*) yl
    endif
    close(2)
    
end subroutine GetSize

!==================================================================================================================================================
!==================================================================================================================================================

subroutine ConvertFieldToRGB(field,rgb,n,minf,maxf)
  use share_vars
    implicit none
    integer 					:: i,j, intensity,n
    real (kind=pr) 				:: minf,maxf
    real (kind=pr), dimension(1:nx,1:ny) 	::  field
    character*1, dimension(1:3,1:nx,1:ny) 	:: rgb
    integer, dimension(1:nPalettes,1:256*3) 		:: palette  !array for color values for FARGE palettes
    integer :: ixmin,ixmax,iymin,iymax

    if ( (n<0).or.(n>nPalettes) ) then
	write(*,*) n, "ConvertFieldToRGB:: Unknown Palette. Stop."
	stop
    endif

    !create palettes
    call FargePalette(palette)
    !loop over the field
    do i = 1, nx
    do j = 1, ny
      !intensity is the index of the color place in the palette array
      intensity = 1+nint(256.0*(field(i,j)-minf)/(maxf-minf)) !indexing starts with 1 so add 1
      if (intensity<1) intensity=1 
      if (intensity>255) intensity=255
      rgb(3,i,j) = char(Palette(n,intensity+512) ) !blue
      rgb(2,i,j) = char(Palette(n,intensity+256) ) !green
      rgb(1,i,j) = char(Palette(n,intensity)     ) !red
  
    enddo
    enddo


end subroutine ConvertFieldToRGB

!=============================================================================================================================================
subroutine ResizeField(field_source, field_target,lx,ly)
! this subroutine resizes a field to another resolution.
! for simplicity, it uses only bi-linear interpolation, even if the grids match at most of the points.
! as it is usually called only once, its performance is not important

  use share_vars
  real (kind=pr), intent (inout)    :: field_source (:,:)
  real (kind=pr), intent (out)   :: field_target (:,:)
  real (kind=pr), intent (in) :: lx,ly
  real (kind=pr) :: dx1,dx2,dy1,dy2,x,y,x1,y1,x2,y2,R1,R2, summe
  integer :: nx1,nx2,ny1,ny2,ix,iy,i,j
  write (*,*) "----------------- resize field --------------------"
  write (*,*) "source:", shape(field_source)
  write (*,*) "target:", shape(field_target)
  
  ! the size of assumed-shaped arrays is dynamic. remember that indices always start with 1!
  nx1=size(field_source,1)
  ny1=size(field_source,2)
  nx2=size(field_target,1)
  ny2=size(field_target,2)
  write(*,'(4(i5,1x))') nx1,ny1,nx2,ny2

  ! determine lattice spacing. here, the domain is not periodic so its lx/nx-1
  dx1=lx/real(nx1-1,kind=pr)
  dy1=ly/real(ny1-1,kind=pr)
  dx2=lx/real(nx2-1,kind=pr)
  dy2=ly/real(ny2-1,kind=pr)    
  field_target = 9.0
  write (*,'(4(F7.4,1x))') dx1,dy1,lx,ly

  ! fill target field
  do ix=1,nx2
    do iy=1,ny2
      ! coordinate in target field
      x=real(ix,kind=pr)*dx2 !as indices start with one, subtract one
      y=real(iy,kind=pr)*dy2

      ! coordinate in source field
      i=int(x/dx1)
      if (i==nx1) i=nx1-1 ! if the last column and row is reached, change interpolation so that it is not out of bounds
      if (i==0) i=1
      j=int(y/dy1)
      if (j==ny1) j=ny1-1
      if (j==0) j=1

      ! coordinates of interpolation points 1&2
      x1=real(i,kind=pr)  *dx1
      y1=real(j,kind=pr)  *dy1
      x2=real(i+1,kind=pr)*dx1
      y2=real(j+1,kind=pr)*dy1

      ! interpolation in x-direction
      R1 = (x2-x)*field_source(i,j  )/dx1 + (x-x1)*field_source(i+1,j  )/dx1
      R2 = (x2-x)*field_source(i,j+1)/dx1 + (x-x1)*field_source(i+1,j+1)/dx1      
      ! interpolation in y-direction
      field_target(ix,iy) = (y2-y)*R1/dy1 + (y-y1)*R2/dy1

    enddo
  enddo
end subroutine ResizeField

!=============================================================================================================================================

subroutine SaveField(field,filename, NSave, lx,ly)
! SEE LOAD_FIELD for details on the format 
  use share_vars
  implicit none
  real (kind=pr), intent (in)   :: field(:,:)
  real (kind=pr), intent (in), optional   :: lx,ly
  real (kind=pr) :: llx,lly
  integer, intent(in), optional          :: NSave
  integer	                :: ix,iy,nx1,ny1, iSave
  character*(*), intent(in)     :: filename
  

  nx1=size(field,1)
  ny1=size(field,2)  

  if (present(NSave) ) then
    iSave=NSave
  else
    iSave =1
  endif
  
  if (present(lx) ) then
    llx=lx
  else
    llx=xl
  endif
  
  if (present(ly) ) then
    lly=ly
  else
    lly=yl
  endif

  open  (10, file=filename, form='formatted', status='replace')
  write (10,'(A)') "%format1"
  write (10,'(i5.5,3x,A)') nx1, "% nx"
  write (10,'(i5.5,3x,A)') ny1, "% ny"
  write (10,'(es11.4,3x,A)') llx, "% xl, domain size in x-direction"
  write (10,'(es11.4,3x,A)') lly, "% yl, domain size in y-direction"
  write (10,'(A)') "% The next entry can be used to identify the type of the field: "
  write (10,'(A)') "% OBSOLETE"


  do iy = 0, ny1-1, iSave  
  write (10, '(es16.9)') (field(ix, ny1-iy) , ix=1,nx1, iSave)
  end do
  
  close (10)

end subroutine SaveField


!=============================================================================================================================================

subroutine FargePalette(palette)
      use share_vars
	implicit none
       integer palette(nPalettes,256*3)
       integer i,sector, n,j
       reaL :: intensity
       Palette(1,:) =  (/ Z'00', Z'00', Z'00', Z'00', Z'01', Z'01', Z'01', Z'01',&   ! Rot */
			  Z'02', Z'02', Z'02', Z'02', Z'03', Z'03', Z'03', Z'03',&
			  Z'03', Z'04', Z'04', Z'04', Z'04', Z'04', Z'04', Z'05',&
			  Z'05', Z'05', Z'05', Z'05', Z'05', Z'06', Z'06', Z'06',&
			  Z'06', Z'06', Z'06', Z'06', Z'07', Z'07', Z'07', Z'07',&
			  Z'07', Z'07', Z'07', Z'08', Z'08', Z'08', Z'08', Z'08',&
			  Z'08', Z'08', Z'09', Z'09', Z'09', Z'09', Z'09', Z'09',&
			  Z'09', Z'0a', Z'0a', Z'0a', Z'0a', Z'0a', Z'0a', Z'0a',&
			  Z'0b', Z'0b', Z'0b', Z'0b', Z'0b', Z'0c', Z'0c', Z'0c',&
			  Z'0c', Z'0c', Z'0c', Z'0d', Z'0d', Z'0d', Z'0d', Z'0e',&
			  Z'0e', Z'0e', Z'0e', Z'0f', Z'0f', Z'0f', Z'0f', Z'10',&
			  Z'00', Z'04', Z'09', Z'0e', Z'12', Z'17', Z'1b', Z'1f',&
			  Z'23', Z'27', Z'2b', Z'2e', Z'32', Z'36', Z'39', Z'3c',&
			  Z'40', Z'43', Z'46', Z'49', Z'4c', Z'4f', Z'52', Z'55',&
			  Z'58', Z'5b', Z'5f', Z'62', Z'65', Z'68', Z'6b', Z'6d',&
			  Z'70', Z'72', Z'74', Z'76', Z'78', Z'7a', Z'fd', Z'fd',&
			  Z'fd', Z'fd', Z'84', Z'85', Z'87', Z'89', Z'8b', Z'8e',&
			  Z'90', Z'92', Z'95', Z'98', Z'9a', Z'9e', Z'a1', Z'a4',&
			  Z'a8', Z'ab', Z'ae', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd',&
			  Z'c1', Z'c4', Z'c7', Z'cb', Z'cf', Z'd2', Z'd6', Z'da',&
			  Z'de', Z'e2', Z'e6', Z'eb', Z'ef', Z'f4', Z'f9', Z'fd',&
			  Z'40', Z'44', Z'47', Z'4a', Z'4d', Z'50', Z'53', Z'56',&
			  Z'58', Z'5b', Z'5e', Z'60', Z'63', Z'65', Z'68', Z'6a',&
			  Z'6d', Z'6f', Z'71', Z'73', Z'75', Z'77', Z'79', Z'7b',&
			  Z'7d', Z'7f', Z'81', Z'83', Z'85', Z'87', Z'89', Z'8a',&
			  Z'8c', Z'8e', Z'90', Z'91', Z'93', Z'95', Z'96', Z'98',&
			  Z'99', Z'9b', Z'9d', Z'9e', Z'a0', Z'a2', Z'a3', Z'a5',&
			  Z'a6', Z'a8', Z'aa', Z'ab', Z'ad', Z'af', Z'b0', Z'b2',&
			  Z'b4', Z'b6', Z'b7', Z'b9', Z'bb', Z'bd', Z'bf', Z'c1',&
			  Z'c3', Z'c5', Z'c7', Z'c9', Z'cb', Z'cd', Z'cf', Z'd2',&
			  Z'd4', Z'd6', Z'd9', Z'db', Z'de', Z'e0', Z'e3', Z'e6',&
			  Z'e9', Z'eb', Z'ee', Z'f1', Z'f4', Z'f7', Z'fb', Z'fe',&

			  Z'3f', Z'41', Z'42', Z'44', Z'45', Z'47', Z'48', Z'4a',&   ! Gruen */
			  Z'4b', Z'4c', Z'4e', Z'4f', Z'50', Z'51', Z'53', Z'54',&
			  Z'55', Z'56', Z'57', Z'58', Z'59', Z'5a', Z'5b', Z'5c',&
			  Z'5d', Z'5e', Z'5f', Z'60', Z'61', Z'62', Z'63', Z'64',&
			  Z'65', Z'66', Z'66', Z'67', Z'68', Z'69', Z'6a', Z'6a',&
			  Z'6b', Z'6c', Z'6d', Z'6e', Z'6f', Z'6f', Z'70', Z'71',&
			  Z'72', Z'73', Z'73', Z'74', Z'75', Z'76', Z'77', Z'78',&
			  Z'78', Z'79', Z'7a', Z'7b', Z'7c', Z'7d', Z'7e', Z'7f',&
			  Z'80', Z'81', Z'82', Z'83', Z'84', Z'85', Z'86', Z'87',&
			  Z'89', Z'8a', Z'8b', Z'8c', Z'8d', Z'8f', Z'90', Z'91',&
			  Z'93', Z'94', Z'96', Z'97', Z'99', Z'9a', Z'9c', Z'9d',&
			  Z'00', Z'04', Z'09', Z'0e', Z'12', Z'17', Z'1b', Z'1f',&
			  Z'23', Z'27', Z'2b', Z'2e', Z'32', Z'36', Z'39', Z'3c',&
			  Z'40', Z'43', Z'46', Z'49', Z'4c', Z'4f', Z'52', Z'55',&
			  Z'58', Z'5b', Z'5f', Z'62', Z'65', Z'68', Z'6b', Z'6d',&
			  Z'70', Z'72', Z'74', Z'76', Z'78', Z'7a', Z'd8', Z'd8',&
			  Z'd8', Z'd8', Z'84', Z'85', Z'87', Z'89', Z'8b', Z'8e',&
			  Z'90', Z'92', Z'95', Z'98', Z'9a', Z'9e', Z'a1', Z'a4',&
			  Z'a8', Z'ab', Z'ae', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd',&
			  Z'c1', Z'c4', Z'c7', Z'cb', Z'cf', Z'd2', Z'd6', Z'da',&
			  Z'de', Z'e2', Z'e6', Z'eb', Z'ef', Z'f4', Z'f9', Z'fd',&
			  Z'00', Z'02', Z'04', Z'06', Z'08', Z'0a', Z'0c', Z'0e',&
			  Z'10', Z'12', Z'14', Z'15', Z'17', Z'19', Z'1a', Z'1c',&
			  Z'1e', Z'1f', Z'21', Z'22', Z'23', Z'25', Z'26', Z'28',&
			  Z'29', Z'2a', Z'2c', Z'2d', Z'2e', Z'2f', Z'31', Z'32',&
			  Z'33', Z'34', Z'35', Z'36', Z'37', Z'39', Z'3a', Z'3b',&
			  Z'3c', Z'3d', Z'3e', Z'3f', Z'40', Z'41', Z'42', Z'44',&
			  Z'45', Z'46', Z'47', Z'48', Z'49', Z'4a', Z'4b', Z'4d',&
			  Z'4e', Z'4f', Z'50', Z'51', Z'53', Z'54', Z'55', Z'57',&
			  Z'58', Z'59', Z'5b', Z'5c', Z'5d', Z'5f', Z'60', Z'62',&
			  Z'64', Z'65', Z'67', Z'68', Z'6a', Z'6c', Z'6e', Z'70',&
			  Z'71', Z'73', Z'75', Z'77', Z'79', Z'7c', Z'7e', Z'80',&

			  Z'53', Z'56', Z'59', Z'5c', Z'5f', Z'62', Z'64', Z'67',&   ! Blau */
			  Z'69', Z'6c', Z'6e', Z'70', Z'73', Z'75', Z'77', Z'79',&
			  Z'7b', Z'7e', Z'80', Z'81', Z'83', Z'85', Z'87', Z'89',&
			  Z'8b', Z'8d', Z'8e', Z'90', Z'92', Z'93', Z'95', Z'97',&
			  Z'98', Z'9a', Z'9b', Z'9d', Z'9e', Z'a0', Z'a1', Z'a3',&
			  Z'a4', Z'a6', Z'a7', Z'a9', Z'aa', Z'ac', Z'ad', Z'af',&
			  Z'b0', Z'b1', Z'b3', Z'b4', Z'b6', Z'b8', Z'b9', Z'bb',&
			  Z'bc', Z'be', Z'bf', Z'c1', Z'c3', Z'c4', Z'c6', Z'c8',&
			  Z'ca', Z'cc', Z'cd', Z'cf', Z'd1', Z'd3', Z'd5', Z'd7',&
			  Z'd9', Z'dc', Z'de', Z'e0', Z'e2', Z'e5', Z'e7', Z'e9',&
			  Z'ec', Z'ef', Z'f1', Z'f4', Z'f7', Z'fa', Z'fc', Z'ff',&
			  Z'00', Z'04', Z'09', Z'0e', Z'12', Z'17', Z'1b', Z'1f',&
			  Z'23', Z'27', Z'2b', Z'2e', Z'32', Z'36', Z'39', Z'3c',&
			  Z'40', Z'43', Z'46', Z'49', Z'4c', Z'4f', Z'52', Z'55',&
			  Z'58', Z'5b', Z'5f', Z'62', Z'65', Z'68', Z'6b', Z'6d',&
			  Z'70', Z'72', Z'74', Z'76', Z'78', Z'7a', Z'00', Z'00',&
			  Z'00', Z'00', Z'84', Z'85', Z'87', Z'89', Z'8b', Z'8e',&
			  Z'90', Z'92', Z'95', Z'98', Z'9a', Z'9e', Z'a1', Z'a4',&
			  Z'a8', Z'ab', Z'ae', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd',&
			  Z'c1', Z'c4', Z'c7', Z'cb', Z'cf', Z'd2', Z'd6', Z'da',&
			  Z'de', Z'e2', Z'e6', Z'eb', Z'ef', Z'f4', Z'f9', Z'fd',&
			  Z'00', Z'02', Z'04', Z'06', Z'08', Z'0a', Z'0c', Z'0e',&
			  Z'10', Z'12', Z'14', Z'15', Z'17', Z'19', Z'1a', Z'1c',&
			  Z'1e', Z'1f', Z'21', Z'22', Z'23', Z'25', Z'26', Z'28',&
			  Z'29', Z'2a', Z'2c', Z'2d', Z'2e', Z'2f', Z'31', Z'32',&
			  Z'33', Z'34', Z'35', Z'36', Z'37', Z'39', Z'3a', Z'3b',&
			  Z'3c', Z'3d', Z'3e', Z'3f', Z'40', Z'41', Z'42', Z'44',&
			  Z'45', Z'46', Z'47', Z'48', Z'49', Z'4a', Z'4b', Z'4d',&
			  Z'4e', Z'4f', Z'50', Z'51', Z'53', Z'54', Z'55', Z'57',&
			  Z'58', Z'59', Z'5b', Z'5c', Z'5d', Z'5f', Z'60', Z'62',&
			  Z'64', Z'65', Z'67', Z'68', Z'6a', Z'6c', Z'6e', Z'70',&
			  Z'71', Z'73', Z'75', Z'77', Z'79', Z'7c', Z'7e', Z'80' /)
!                                                                ! Palette 2 */
Palette(2,:)=(/ Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&   ! Rot */
		Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
		Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
		Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
		Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
		Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
		Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
		Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
		Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
		Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
		Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
		Z'00', Z'04', Z'09', Z'0e', Z'12', Z'17', Z'1b', Z'1f',&
		Z'23', Z'27', Z'2b', Z'2e', Z'32', Z'36', Z'39', Z'3c',&
		Z'40', Z'43', Z'46', Z'49', Z'4c', Z'4f', Z'52', Z'55',&
		Z'58', Z'5a', Z'5d', Z'60', Z'62', Z'65', Z'67', Z'6a',&
		Z'6c', Z'6f', Z'71', Z'74', Z'76', Z'79', Z'fe', Z'fe',&
		Z'fe', Z'82', Z'84', Z'87', Z'89', Z'8c', Z'8e', Z'91',&
		Z'93', Z'96', Z'98', Z'9b', Z'9d', Z'a0', Z'a3', Z'a5',&
		Z'a8', Z'ab', Z'ae', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd',&
		Z'c1', Z'c4', Z'c7', Z'cb', Z'cf', Z'd2', Z'd6', Z'da',&
		Z'de', Z'e2', Z'e6', Z'eb', Z'ef', Z'f4', Z'f9', Z'fd',&
		Z'53', Z'56', Z'59', Z'5c', Z'5f', Z'62', Z'64', Z'67',&
		Z'69', Z'6c', Z'6e', Z'70', Z'73', Z'75', Z'77', Z'79',&
		Z'7b', Z'7e', Z'80', Z'81', Z'83', Z'85', Z'87', Z'89',&
		Z'8b', Z'8d', Z'8e', Z'90', Z'92', Z'93', Z'95', Z'97',&
		Z'98', Z'9a', Z'9b', Z'9d', Z'9e', Z'a0', Z'a1', Z'a3',&
		Z'a4', Z'a6', Z'a7', Z'a9', Z'aa', Z'ac', Z'ad', Z'af',&
		Z'b0', Z'b1', Z'b3', Z'b4', Z'b6', Z'b8', Z'b9', Z'bb',&
		Z'bc', Z'be', Z'bf', Z'c1', Z'c3', Z'c4', Z'c6', Z'c8',&
		Z'ca', Z'cc', Z'cd', Z'cf', Z'd1', Z'd3', Z'd5', Z'd7',&
		Z'd9', Z'dc', Z'de', Z'e0', Z'e2', Z'e5', Z'e7', Z'e9',&
		Z'ec', Z'ef', Z'f1', Z'f4', Z'f7', Z'fa', Z'fc', Z'ff',&
		
		Z'3f', Z'40', Z'42', Z'43', Z'44', Z'45', Z'47', Z'48',&   ! Gruen */
		Z'49', Z'4a', Z'4b', Z'4c', Z'4d', Z'4e', Z'4f', Z'50',&
		Z'51', Z'52', Z'53', Z'54', Z'55', Z'56', Z'57', Z'57',&
		Z'58', Z'59', Z'5a', Z'5b', Z'5b', Z'5c', Z'5d', Z'5e',&
		Z'5e', Z'5f', Z'60', Z'60', Z'61', Z'62', Z'62', Z'63',&
		Z'64', Z'64', Z'65', Z'66', Z'66', Z'67', Z'68', Z'68',&
		Z'69', Z'6a', Z'6b', Z'6b', Z'6c', Z'6d', Z'6d', Z'6e',&
		Z'6f', Z'6f', Z'70', Z'71', Z'72', Z'72', Z'73', Z'74',&
		Z'75', Z'76', Z'77', Z'77', Z'78', Z'79', Z'7a', Z'7b',&
		Z'7c', Z'7d', Z'7e', Z'7f', Z'80', Z'81', Z'82', Z'83',&
		Z'85', Z'86', Z'87', Z'88', Z'89', Z'8b', Z'8c', Z'8d',&
		Z'00', Z'04', Z'09', Z'0e', Z'12', Z'17', Z'1b', Z'1f',&
		Z'23', Z'27', Z'2b', Z'2e', Z'32', Z'36', Z'39', Z'3c',&
		Z'40', Z'43', Z'46', Z'49', Z'4c', Z'4f', Z'52', Z'55',&
		Z'58', Z'5a', Z'5d', Z'60', Z'62', Z'65', Z'67', Z'6a',&
		Z'6c', Z'6f', Z'71', Z'74', Z'76', Z'79', Z'00', Z'00',&
		Z'00', Z'82', Z'84', Z'87', Z'89', Z'8c', Z'8e', Z'91',&
		Z'93', Z'96', Z'98', Z'9b', Z'9d', Z'a0', Z'a3', Z'a5',&
		Z'a8', Z'ab', Z'ae', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd',&
		Z'c1', Z'c4', Z'c7', Z'cb', Z'cf', Z'd2', Z'd6', Z'da',&
		Z'de', Z'e2', Z'e6', Z'eb', Z'ef', Z'f4', Z'f9', Z'fd',&
		Z'3c', Z'3f', Z'42', Z'46', Z'49', Z'4c', Z'4f', Z'51',&
		Z'54', Z'57', Z'5a', Z'5c', Z'5f', Z'61', Z'64', Z'66',&
		Z'68', Z'6b', Z'6d', Z'6f', Z'71', Z'73', Z'75', Z'77',&
		Z'79', Z'7b', Z'7d', Z'7f', Z'81', Z'83', Z'85', Z'86',&
		Z'88', Z'8a', Z'8c', Z'8d', Z'8f', Z'91', Z'92', Z'94',&
		Z'96', Z'97', Z'99', Z'9b', Z'9c', Z'9e', Z'9f', Z'a1',&
		Z'a3', Z'a4', Z'a6', Z'a8', Z'a9', Z'ab', Z'ad', Z'af',&
		Z'b0', Z'b2', Z'b4', Z'b6', Z'b8', Z'b9', Z'bb', Z'bd',&
		Z'bf', Z'c1', Z'c3', Z'c6', Z'c8', Z'ca', Z'cc', Z'ce',&
		Z'd1', Z'd3', Z'd6', Z'd8', Z'db', Z'dd', Z'e0', Z'e3',&
		Z'e5', Z'e8', Z'eb', Z'ee', Z'f1', Z'f4', Z'f8', Z'fb',&
		
		Z'53', Z'56', Z'59', Z'5c', Z'5f', Z'61', Z'64', Z'67',&   ! Blau */
		Z'69', Z'6c', Z'6e', Z'70', Z'72', Z'75', Z'77', Z'79',&
		Z'7b', Z'7d', Z'7f', Z'81', Z'83', Z'85', Z'87', Z'89',&
		Z'8a', Z'8c', Z'8e', Z'8f', Z'91', Z'93', Z'94', Z'96',&
		Z'98', Z'99', Z'9b', Z'9c', Z'9e', Z'9f', Z'a1', Z'a2',&
		Z'a4', Z'a5', Z'a6', Z'a8', Z'a9', Z'ab', Z'ac', Z'ae',&
		Z'af', Z'b1', Z'b2', Z'b4', Z'b5', Z'b7', Z'b8', Z'ba',&
		Z'bb', Z'bd', Z'bf', Z'c0', Z'c2', Z'c3', Z'c5', Z'c7',&
		Z'c9', Z'cb', Z'cc', Z'ce', Z'd0', Z'd2', Z'd4', Z'd6',&
		Z'd8', Z'da', Z'dd', Z'df', Z'e1', Z'e3', Z'e6', Z'e8',&
		Z'eb', Z'ed', Z'f0', Z'f3', Z'f5', Z'f8', Z'fb', Z'fe',&
		Z'00', Z'04', Z'09', Z'0e', Z'12', Z'17', Z'1b', Z'1f',&
		Z'23', Z'27', Z'2b', Z'2e', Z'32', Z'36', Z'39', Z'3c',&
		Z'40', Z'43', Z'46', Z'49', Z'4c', Z'4f', Z'52', Z'55',&
		Z'58', Z'5a', Z'5d', Z'60', Z'62', Z'65', Z'67', Z'6a',&
		Z'6c', Z'6f', Z'71', Z'74', Z'76', Z'79', Z'00', Z'00',&
		Z'00', Z'82', Z'84', Z'87', Z'89', Z'8c', Z'8e', Z'91',&
		Z'93', Z'96', Z'98', Z'9b', Z'9d', Z'a0', Z'a3', Z'a5',&
		Z'a8', Z'ab', Z'ae', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd',&
		Z'c1', Z'c4', Z'c7', Z'cb', Z'cf', Z'd2', Z'd6', Z'da',&
		Z'de', Z'e2', Z'e6', Z'eb', Z'ef', Z'f4', Z'f9', Z'fd',&
		Z'00', Z'00', Z'01', Z'01', Z'02', Z'02', Z'03', Z'03',&
		Z'04', Z'04', Z'05', Z'05', Z'05', Z'06', Z'06', Z'07',&
		Z'07', Z'07', Z'08', Z'08', Z'09', Z'09', Z'09', Z'0a',&
		Z'0a', Z'0a', Z'0b', Z'0b', Z'0b', Z'0c', Z'0c', Z'0c',&
		Z'0c', Z'0d', Z'0d', Z'0d', Z'0e', Z'0e', Z'0e', Z'0e',&
		Z'0f', Z'0f', Z'0f', Z'0f', Z'10', Z'10', Z'10', Z'11',&
		Z'11', Z'11', Z'11', Z'12', Z'12', Z'12', Z'13', Z'13',&
		Z'13', Z'13', Z'14', Z'14', Z'14', Z'15', Z'15', Z'15',&
		Z'16', Z'16', Z'16', Z'17', Z'17', Z'17', Z'18', Z'18',&
		Z'19', Z'19', Z'19', Z'1a', Z'1a', Z'1b', Z'1b', Z'1c',&
		Z'1c', Z'1d', Z'1d', Z'1d', Z'1e', Z'1f', Z'1f', Z'20' /)
!                                                                ! Palette 3 */
Palette(3,:)=(/  Z'39', Z'3b', Z'3d', Z'3e', Z'40', Z'41', Z'43', Z'44',&   ! Rot */
 Z'46', Z'47', Z'48', Z'4a', Z'4b', Z'4c', Z'4d', Z'4f',&
 Z'50', Z'51', Z'52', Z'53', Z'54', Z'55', Z'56', Z'57',&
 Z'58', Z'59', Z'5a', Z'5b', Z'5c', Z'5d', Z'5e', Z'5f',&
 Z'60', Z'61', Z'61', Z'62', Z'63', Z'64', Z'65', Z'66',&
 Z'67', Z'67', Z'68', Z'69', Z'6a', Z'6b', Z'6b', Z'6c',&
 Z'6d', Z'6e', Z'6f', Z'70', Z'70', Z'71', Z'72', Z'73',&
 Z'74', Z'75', Z'76', Z'77', Z'78', Z'78', Z'79', Z'7a',&
 Z'7b', Z'7c', Z'7d', Z'7f', Z'80', Z'81', Z'82', Z'83',&
 Z'84', Z'85', Z'87', Z'88', Z'89', Z'8a', Z'8c', Z'8d',&
 Z'8f', Z'90', Z'91', Z'93', Z'95', Z'96', Z'98', Z'99',&
 Z'00', Z'04', Z'09', Z'0e', Z'12', Z'17', Z'1b', Z'1f',&
 Z'23', Z'27', Z'2b', Z'2e', Z'32', Z'36', Z'39', Z'3c',&
 Z'40', Z'43', Z'46', Z'49', Z'4c', Z'4f', Z'52', Z'55',&
 Z'58', Z'5a', Z'5d', Z'60', Z'62', Z'65', Z'67', Z'6a',&
 Z'6c', Z'6f', Z'71', Z'74', Z'76', Z'79', Z'00', Z'00',&
 Z'00', Z'82', Z'84', Z'87', Z'89', Z'8c', Z'8e', Z'91',&
 Z'93', Z'96', Z'98', Z'9b', Z'9d', Z'a0', Z'a3', Z'a5',&
 Z'a8', Z'ab', Z'ae', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd',&
 Z'c1', Z'c4', Z'c7', Z'cb', Z'cf', Z'd2', Z'd6', Z'da',&
 Z'de', Z'e2', Z'e6', Z'eb', Z'ef', Z'f4', Z'f9', Z'fd',&
 Z'33', Z'37', Z'3a', Z'3e', Z'41', Z'44', Z'47', Z'4a',&
 Z'4d', Z'50', Z'53', Z'56', Z'59', Z'5b', Z'5e', Z'60',&
 Z'63', Z'65', Z'68', Z'6a', Z'6c', Z'6f', Z'71', Z'73',&
 Z'75', Z'77', Z'79', Z'7b', Z'7d', Z'7f', Z'81', Z'83',&
 Z'85', Z'87', Z'89', Z'8a', Z'8c', Z'8e', Z'90', Z'92',&
 Z'93', Z'95', Z'97', Z'99', Z'9a', Z'9c', Z'9e', Z'9f',&
 Z'a1', Z'a3', Z'a5', Z'a7', Z'a8', Z'aa', Z'ac', Z'ae',&
 Z'b0', Z'b2', Z'b4', Z'b6', Z'b7', Z'ba', Z'bc', Z'be',&
 Z'c0', Z'c2', Z'c4', Z'c6', Z'c9', Z'cb', Z'cd', Z'd0',&
 Z'd2', Z'd5', Z'd7', Z'da', Z'dd', Z'e0', Z'e3', Z'e5',&
 Z'e8', Z'eb', Z'ef', Z'f2', Z'f5', Z'f8', Z'fc', Z'ff',&
 
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&   ! Gruen */
 Z'00', Z'00', Z'01', Z'01', Z'01', Z'01', Z'01', Z'01',&
 Z'01', Z'01', Z'01', Z'01', Z'01', Z'01', Z'01', Z'01',&
 Z'01', Z'02', Z'02', Z'02', Z'02', Z'02', Z'02', Z'02',&
 Z'02', Z'02', Z'02', Z'02', Z'02', Z'02', Z'02', Z'02',&
 Z'02', Z'02', Z'02', Z'02', Z'03', Z'03', Z'03', Z'03',&
 Z'03', Z'03', Z'03', Z'03', Z'03', Z'03', Z'03', Z'03',&
 Z'03', Z'03', Z'03', Z'03', Z'03', Z'03', Z'03', Z'04',&
 Z'04', Z'04', Z'04', Z'04', Z'04', Z'04', Z'04', Z'04',&
 Z'04', Z'04', Z'04', Z'04', Z'04', Z'04', Z'05', Z'05',&
 Z'05', Z'05', Z'05', Z'05', Z'05', Z'05', Z'05', Z'05',&
 Z'00', Z'04', Z'09', Z'0e', Z'12', Z'17', Z'1b', Z'1f',&
 Z'23', Z'27', Z'2b', Z'2e', Z'32', Z'36', Z'39', Z'3c',&
 Z'40', Z'43', Z'46', Z'49', Z'4c', Z'4f', Z'52', Z'55',&
 Z'58', Z'5a', Z'5d', Z'60', Z'62', Z'65', Z'67', Z'6a',&
 Z'6c', Z'6f', Z'71', Z'74', Z'76', Z'79', Z'fc', Z'fc',&
 Z'fc', Z'82', Z'84', Z'87', Z'89', Z'8c', Z'8e', Z'91',&
 Z'93', Z'96', Z'98', Z'9b', Z'9d', Z'a0', Z'a3', Z'a5',&
 Z'a8', Z'ab', Z'ae', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd',&
 Z'c1', Z'c4', Z'c7', Z'cb', Z'cf', Z'd2', Z'd6', Z'da',&
 Z'de', Z'e2', Z'e6', Z'eb', Z'ef', Z'f4', Z'f9', Z'fd',&
 Z'0c', Z'0f', Z'12', Z'16', Z'19', Z'1c', Z'1e', Z'21',&
 Z'24', Z'27', Z'29', Z'2c', Z'2f', Z'31', Z'34', Z'36',&
 Z'38', Z'3b', Z'3d', Z'3f', Z'41', Z'43', Z'45', Z'47',&
 Z'49', Z'4b', Z'4d', Z'4f', Z'51', Z'53', Z'54', Z'56',&
 Z'58', Z'5a', Z'5b', Z'5d', Z'5f', Z'60', Z'62', Z'64',&
 Z'65', Z'67', Z'69', Z'6a', Z'6c', Z'6e', Z'6f', Z'71',&
 Z'72', Z'74', Z'76', Z'77', Z'79', Z'7b', Z'7d', Z'7e',&
 Z'80', Z'82', Z'84', Z'85', Z'87', Z'89', Z'8b', Z'8d',&
 Z'8f', Z'91', Z'93', Z'95', Z'97', Z'99', Z'9c', Z'9e',&
 Z'a0', Z'a3', Z'a5', Z'a8', Z'aa', Z'ad', Z'af', Z'b2',&
 Z'b5', Z'b8', Z'bb', Z'be', Z'c1', Z'c4', Z'c7', Z'ca',&
 
 Z'7d', Z'7f', Z'81', Z'83', Z'85', Z'87', Z'89', Z'8b',&   ! Blau */
 Z'8d', Z'8f', Z'91', Z'93', Z'94', Z'96', Z'98', Z'99',&
 Z'9b', Z'9c', Z'9e', Z'9f', Z'a1', Z'a2', Z'a3', Z'a5',&
 Z'a6', Z'a7', Z'a9', Z'aa', Z'ab', Z'ac', Z'ae', Z'af',&
 Z'b0', Z'b1', Z'b2', Z'b3', Z'b5', Z'b6', Z'b7', Z'b8',&
 Z'b9', Z'ba', Z'bb', Z'bc', Z'bd', Z'bf', Z'c0', Z'c1',&
 Z'c2', Z'c3', Z'c4', Z'c5', Z'c6', Z'c7', Z'c9', Z'ca',&
 Z'cb', Z'cc', Z'cd', Z'cf', Z'd0', Z'd1', Z'd2', Z'd4',&
 Z'd5', Z'd6', Z'd8', Z'd9', Z'db', Z'dc', Z'de', Z'df',&
 Z'e1', Z'e2', Z'e4', Z'e6', Z'e7', Z'e9', Z'eb', Z'ed',&
 Z'ee', Z'f0', Z'f2', Z'f4', Z'f6', Z'f9', Z'fb', Z'fd',&
 Z'00', Z'04', Z'09', Z'0e', Z'12', Z'17', Z'1b', Z'1f',&
 Z'23', Z'27', Z'2b', Z'2e', Z'32', Z'36', Z'39', Z'3c',&
 Z'40', Z'43', Z'46', Z'49', Z'4c', Z'4f', Z'52', Z'55',&
 Z'58', Z'5a', Z'5d', Z'60', Z'62', Z'65', Z'67', Z'6a',&
 Z'6c', Z'6f', Z'71', Z'74', Z'76', Z'79', Z'8a', Z'8a',&
 Z'8a', Z'82', Z'84', Z'87', Z'89', Z'8c', Z'8e', Z'91',&
 Z'93', Z'96', Z'98', Z'9b', Z'9d', Z'a0', Z'a3', Z'a5',&
 Z'a8', Z'ab', Z'ae', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd',&
 Z'c1', Z'c4', Z'c7', Z'cb', Z'cf', Z'd2', Z'd6', Z'da',&
 Z'de', Z'e2', Z'e6', Z'eb', Z'ef', Z'f4', Z'f9', Z'fd',&
 Z'00', Z'00', Z'01', Z'01', Z'02', Z'02', Z'03', Z'03',&
 Z'04', Z'04', Z'05', Z'05', Z'05', Z'06', Z'06', Z'07',&
 Z'07', Z'07', Z'08', Z'08', Z'09', Z'09', Z'09', Z'0a',&
 Z'0a', Z'0a', Z'0b', Z'0b', Z'0b', Z'0c', Z'0c', Z'0c',&
 Z'0c', Z'0d', Z'0d', Z'0d', Z'0e', Z'0e', Z'0e', Z'0e',&
 Z'0f', Z'0f', Z'0f', Z'0f', Z'10', Z'10', Z'10', Z'11',&
 Z'11', Z'11', Z'11', Z'12', Z'12', Z'12', Z'13', Z'13',&
 Z'13', Z'13', Z'14', Z'14', Z'14', Z'15', Z'15', Z'15',&
 Z'16', Z'16', Z'16', Z'17', Z'17', Z'17', Z'18', Z'18',&
 Z'19', Z'19', Z'19', Z'1a', Z'1a', Z'1b', Z'1b', Z'1c',&
 Z'1c', Z'1d', Z'1d', Z'1d', Z'1e', Z'1f', Z'1f', Z'20'/)
!                                                                ! Palette 4 */
Palette(4,:)=(/  Z'40', Z'44', Z'47', Z'4a', Z'4d', Z'50', Z'53', Z'56',&   ! Rot */
 Z'58', Z'5b', Z'5e', Z'60', Z'63', Z'65', Z'68', Z'6a',&
 Z'6d', Z'6f', Z'71', Z'73', Z'75', Z'77', Z'79', Z'7b',&
 Z'7d', Z'7f', Z'81', Z'83', Z'85', Z'87', Z'89', Z'8a',&
 Z'8c', Z'8e', Z'90', Z'91', Z'93', Z'95', Z'96', Z'98',&
 Z'99', Z'9b', Z'9d', Z'9e', Z'a0', Z'a2', Z'a3', Z'a5',&
 Z'a6', Z'a8', Z'aa', Z'ab', Z'ad', Z'af', Z'b0', Z'b2',&
 Z'b4', Z'b6', Z'b7', Z'b9', Z'bb', Z'bd', Z'bf', Z'c1',&
 Z'c3', Z'c5', Z'c7', Z'c9', Z'cb', Z'cd', Z'cf', Z'd2',&
 Z'd4', Z'd6', Z'd9', Z'db', Z'de', Z'e0', Z'e3', Z'e6',&
 Z'e9', Z'eb', Z'ee', Z'f1', Z'f4', Z'f7', Z'fb', Z'fe',&
 Z'00', Z'04', Z'09', Z'0e', Z'12', Z'17', Z'1b', Z'1f',&
 Z'23', Z'27', Z'2b', Z'2e', Z'32', Z'36', Z'39', Z'3c',&
 Z'40', Z'43', Z'46', Z'49', Z'4c', Z'4f', Z'52', Z'55',&
 Z'58', Z'5a', Z'5d', Z'60', Z'62', Z'65', Z'67', Z'6a',&
 Z'6c', Z'6f', Z'71', Z'74', Z'76', Z'79', Z'7b', Z'00',&
 Z'80', Z'82', Z'84', Z'87', Z'89', Z'8c', Z'8e', Z'91',&
 Z'93', Z'96', Z'98', Z'9b', Z'9d', Z'a0', Z'a3', Z'a5',&
 Z'a8', Z'ab', Z'ae', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd',&
 Z'c1', Z'c4', Z'c7', Z'cb', Z'cf', Z'd2', Z'd6', Z'da',&
 Z'de', Z'e2', Z'e6', Z'eb', Z'ef', Z'f4', Z'f9', Z'fd',&
 Z'53', Z'56', Z'59', Z'5c', Z'5f', Z'61', Z'64', Z'67',&
 Z'69', Z'6c', Z'6e', Z'70', Z'72', Z'75', Z'77', Z'79',&
 Z'7b', Z'7d', Z'7f', Z'81', Z'83', Z'85', Z'87', Z'89',&
 Z'8a', Z'8c', Z'8e', Z'8f', Z'91', Z'93', Z'94', Z'96',&
 Z'98', Z'99', Z'9b', Z'9c', Z'9e', Z'9f', Z'a1', Z'a2',&
 Z'a4', Z'a5', Z'a6', Z'a8', Z'a9', Z'ab', Z'ac', Z'ae',&
 Z'af', Z'b1', Z'b2', Z'b4', Z'b5', Z'b7', Z'b8', Z'ba',&
 Z'bb', Z'bd', Z'bf', Z'c0', Z'c2', Z'c3', Z'c5', Z'c7',&
 Z'c9', Z'cb', Z'cc', Z'ce', Z'd0', Z'd2', Z'd4', Z'd6',&
 Z'd8', Z'da', Z'dd', Z'df', Z'e1', Z'e3', Z'e6', Z'e8',&
 Z'eb', Z'ed', Z'f0', Z'f3', Z'f5', Z'f8', Z'fb', Z'fe',&
 
 Z'00', Z'02', Z'04', Z'06', Z'08', Z'0a', Z'0c', Z'0e',&   ! Gruen */
 Z'10', Z'12', Z'14', Z'15', Z'17', Z'19', Z'1a', Z'1c',&
 Z'1e', Z'1f', Z'21', Z'22', Z'23', Z'25', Z'26', Z'28',&
 Z'29', Z'2a', Z'2c', Z'2d', Z'2e', Z'2f', Z'31', Z'32',&
 Z'33', Z'34', Z'35', Z'36', Z'37', Z'39', Z'3a', Z'3b',&
 Z'3c', Z'3d', Z'3e', Z'3f', Z'40', Z'41', Z'42', Z'44',&
 Z'45', Z'46', Z'47', Z'48', Z'49', Z'4a', Z'4b', Z'4d',&
 Z'4e', Z'4f', Z'50', Z'51', Z'53', Z'54', Z'55', Z'57',&
 Z'58', Z'59', Z'5b', Z'5c', Z'5d', Z'5f', Z'60', Z'62',&
 Z'64', Z'65', Z'67', Z'68', Z'6a', Z'6c', Z'6e', Z'70',&
 Z'71', Z'73', Z'75', Z'77', Z'79', Z'7c', Z'7e', Z'80',&
 Z'00', Z'04', Z'09', Z'0e', Z'12', Z'17', Z'1b', Z'1f',&
 Z'23', Z'27', Z'2b', Z'2e', Z'32', Z'36', Z'39', Z'3c',&
 Z'40', Z'43', Z'46', Z'49', Z'4c', Z'4f', Z'52', Z'55',&
 Z'58', Z'5a', Z'5d', Z'60', Z'62', Z'65', Z'67', Z'6a',&
 Z'6c', Z'6f', Z'71', Z'74', Z'76', Z'79', Z'7b', Z'fb',&
 Z'80', Z'82', Z'84', Z'87', Z'89', Z'8c', Z'8e', Z'91',&
 Z'93', Z'96', Z'98', Z'9b', Z'9d', Z'a0', Z'a3', Z'a5',&
 Z'a8', Z'ab', Z'ae', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd',&
 Z'c1', Z'c4', Z'c7', Z'cb', Z'cf', Z'd2', Z'd6', Z'da',&
 Z'de', Z'e2', Z'e6', Z'eb', Z'ef', Z'f4', Z'f9', Z'fd',&
 Z'3c', Z'3f', Z'41', Z'44', Z'46', Z'49', Z'4b', Z'4e',&
 Z'50', Z'52', Z'55', Z'57', Z'59', Z'5b', Z'5d', Z'5f',&
 Z'61', Z'63', Z'65', Z'67', Z'68', Z'6a', Z'6c', Z'6d',&
 Z'6f', Z'71', Z'72', Z'74', Z'75', Z'77', Z'78', Z'7a',&
 Z'7b', Z'7d', Z'7e', Z'80', Z'81', Z'83', Z'84', Z'85',&
 Z'87', Z'88', Z'89', Z'8b', Z'8c', Z'8d', Z'8f', Z'90',&
 Z'91', Z'93', Z'94', Z'96', Z'97', Z'98', Z'9a', Z'9b',&
 Z'9d', Z'9e', Z'a0', Z'a1', Z'a3', Z'a4', Z'a6', Z'a8',&
 Z'a9', Z'ab', Z'ad', Z'ae', Z'b0', Z'b2', Z'b4', Z'b6',&
 Z'b8', Z'ba', Z'bc', Z'be', Z'c0', Z'c2', Z'c4', Z'c7',&
 Z'c9', Z'cb', Z'ce', Z'd0', Z'd3', Z'd5', Z'd8', Z'db',&
 
 Z'00', Z'02', Z'04', Z'06', Z'08', Z'0a', Z'0c', Z'0e',&   ! Blau */
 Z'10', Z'12', Z'14', Z'15', Z'17', Z'19', Z'1a', Z'1c',&
 Z'1e', Z'1f', Z'21', Z'22', Z'23', Z'25', Z'26', Z'28',&
 Z'29', Z'2a', Z'2c', Z'2d', Z'2e', Z'2f', Z'31', Z'32',&
 Z'33', Z'34', Z'35', Z'36', Z'37', Z'39', Z'3a', Z'3b',&
 Z'3c', Z'3d', Z'3e', Z'3f', Z'40', Z'41', Z'42', Z'44',&
 Z'45', Z'46', Z'47', Z'48', Z'49', Z'4a', Z'4b', Z'4d',&
 Z'4e', Z'4f', Z'50', Z'51', Z'53', Z'54', Z'55', Z'57',&
 Z'58', Z'59', Z'5b', Z'5c', Z'5d', Z'5f', Z'60', Z'62',&
 Z'64', Z'65', Z'67', Z'68', Z'6a', Z'6c', Z'6e', Z'70',&
 Z'71', Z'73', Z'75', Z'77', Z'79', Z'7c', Z'7e', Z'80',&
 Z'00', Z'04', Z'09', Z'0e', Z'12', Z'17', Z'1b', Z'1f',&
 Z'23', Z'27', Z'2b', Z'2e', Z'32', Z'36', Z'39', Z'3c',&
 Z'40', Z'43', Z'46', Z'49', Z'4c', Z'4f', Z'52', Z'55',&
 Z'58', Z'5a', Z'5d', Z'60', Z'62', Z'65', Z'67', Z'6a',&
 Z'6c', Z'6f', Z'71', Z'74', Z'76', Z'79', Z'7b', Z'fb',&
 Z'80', Z'82', Z'84', Z'87', Z'89', Z'8c', Z'8e', Z'91',&
 Z'93', Z'96', Z'98', Z'9b', Z'9d', Z'a0', Z'a3', Z'a5',&
 Z'a8', Z'ab', Z'ae', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd',&
 Z'c1', Z'c4', Z'c7', Z'cb', Z'cf', Z'd2', Z'd6', Z'da',&
 Z'de', Z'e2', Z'e6', Z'eb', Z'ef', Z'f4', Z'f9', Z'fd',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00'/)
!                                                                ! Palette 5 */
Palette(5,:)=(/  Z'00', Z'00', Z'00', Z'00', Z'01', Z'01', Z'01', Z'01',&   ! Rot */
 Z'02', Z'02', Z'02', Z'02', Z'03', Z'03', Z'03', Z'03',&
 Z'03', Z'04', Z'04', Z'04', Z'04', Z'04', Z'04', Z'05',&
 Z'05', Z'05', Z'05', Z'05', Z'05', Z'06', Z'06', Z'06',&
 Z'06', Z'06', Z'06', Z'06', Z'07', Z'07', Z'07', Z'07',&
 Z'07', Z'07', Z'07', Z'08', Z'08', Z'08', Z'08', Z'08',&
 Z'08', Z'08', Z'09', Z'09', Z'09', Z'09', Z'09', Z'09',&
 Z'09', Z'0a', Z'0a', Z'0a', Z'0a', Z'0a', Z'0a', Z'0a',&
 Z'0b', Z'0b', Z'0b', Z'0b', Z'0b', Z'0c', Z'0c', Z'0c',&
 Z'0c', Z'0c', Z'0c', Z'0d', Z'0d', Z'0d', Z'0d', Z'0e',&
 Z'0e', Z'0e', Z'0e', Z'0f', Z'0f', Z'0f', Z'0f', Z'10',&
 Z'00', Z'04', Z'09', Z'0e', Z'12', Z'17', Z'1b', Z'1f',&
 Z'23', Z'27', Z'2b', Z'2e', Z'32', Z'36', Z'39', Z'3c',&
 Z'40', Z'43', Z'46', Z'49', Z'4c', Z'4f', Z'52', Z'55',&
 Z'58', Z'5b', Z'5f', Z'62', Z'65', Z'68', Z'6b', Z'6d',&
 Z'70', Z'72', Z'74', Z'76', Z'78', Z'7a', Z'fd', Z'fd',&
 Z'fd', Z'fd', Z'84', Z'85', Z'87', Z'89', Z'8b', Z'8e',&
 Z'90', Z'92', Z'95', Z'98', Z'9a', Z'9e', Z'a1', Z'a4',&
 Z'a8', Z'ab', Z'ae', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd',&
 Z'c1', Z'c4', Z'c7', Z'cb', Z'cf', Z'd2', Z'd6', Z'da',&
 Z'de', Z'e2', Z'e6', Z'eb', Z'ef', Z'f4', Z'f9', Z'fd',&
 Z'40', Z'44', Z'47', Z'4a', Z'4d', Z'50', Z'53', Z'56',&
 Z'58', Z'5b', Z'5e', Z'60', Z'63', Z'65', Z'68', Z'6a',&
 Z'6d', Z'6f', Z'71', Z'73', Z'75', Z'77', Z'79', Z'7b',&
 Z'7d', Z'7f', Z'81', Z'83', Z'85', Z'87', Z'89', Z'8a',&
 Z'8c', Z'8e', Z'90', Z'91', Z'93', Z'95', Z'96', Z'98',&
 Z'99', Z'9b', Z'9d', Z'9e', Z'a0', Z'a2', Z'a3', Z'a5',&
 Z'a6', Z'a8', Z'aa', Z'ab', Z'ad', Z'af', Z'b0', Z'b2',&
 Z'b4', Z'b6', Z'b7', Z'b9', Z'bb', Z'bd', Z'bf', Z'c1',&
 Z'c3', Z'c5', Z'c7', Z'c9', Z'cb', Z'cd', Z'cf', Z'd2',&
 Z'd4', Z'd6', Z'd9', Z'db', Z'de', Z'e0', Z'e3', Z'e6',&
 Z'e9', Z'eb', Z'ee', Z'f1', Z'f4', Z'f7', Z'fb', Z'fe',&
 
 Z'3f', Z'41', Z'42', Z'44', Z'45', Z'47', Z'48', Z'4a',&   ! Gruen */
 Z'4b', Z'4c', Z'4e', Z'4f', Z'50', Z'51', Z'53', Z'54',&
 Z'55', Z'56', Z'57', Z'58', Z'59', Z'5a', Z'5b', Z'5c',&
 Z'5d', Z'5e', Z'5f', Z'60', Z'61', Z'62', Z'63', Z'64',&
 Z'65', Z'66', Z'66', Z'67', Z'68', Z'69', Z'6a', Z'6a',&
 Z'6b', Z'6c', Z'6d', Z'6e', Z'6f', Z'6f', Z'70', Z'71',&
 Z'72', Z'73', Z'73', Z'74', Z'75', Z'76', Z'77', Z'78',&
 Z'78', Z'79', Z'7a', Z'7b', Z'7c', Z'7d', Z'7e', Z'7f',&
 Z'80', Z'81', Z'82', Z'83', Z'84', Z'85', Z'86', Z'87',&
 Z'89', Z'8a', Z'8b', Z'8c', Z'8d', Z'8f', Z'90', Z'91',&
 Z'93', Z'94', Z'96', Z'97', Z'99', Z'9a', Z'9c', Z'9d',&
 Z'00', Z'04', Z'09', Z'0e', Z'12', Z'17', Z'1b', Z'1f',&
 Z'23', Z'27', Z'2b', Z'2e', Z'32', Z'36', Z'39', Z'3c',&
 Z'40', Z'43', Z'46', Z'49', Z'4c', Z'4f', Z'52', Z'55',&
 Z'58', Z'5b', Z'5f', Z'62', Z'65', Z'68', Z'6b', Z'6d',&
 Z'70', Z'72', Z'74', Z'76', Z'78', Z'7a', Z'd8', Z'd8',&
 Z'd8', Z'd8', Z'84', Z'85', Z'87', Z'89', Z'8b', Z'8e',&
 Z'90', Z'92', Z'95', Z'98', Z'9a', Z'9e', Z'a1', Z'a4',&
 Z'a8', Z'ab', Z'ae', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd',&
 Z'c1', Z'c4', Z'c7', Z'cb', Z'cf', Z'd2', Z'd6', Z'da',&
 Z'de', Z'e2', Z'e6', Z'eb', Z'ef', Z'f4', Z'f9', Z'fd',&
 Z'00', Z'02', Z'04', Z'06', Z'08', Z'0a', Z'0c', Z'0e',&
 Z'10', Z'12', Z'14', Z'15', Z'17', Z'19', Z'1a', Z'1c',&
 Z'1e', Z'1f', Z'21', Z'22', Z'23', Z'25', Z'26', Z'28',&
 Z'29', Z'2a', Z'2c', Z'2d', Z'2e', Z'2f', Z'31', Z'32',&
 Z'33', Z'34', Z'35', Z'36', Z'37', Z'39', Z'3a', Z'3b',&
 Z'3c', Z'3d', Z'3e', Z'3f', Z'40', Z'41', Z'42', Z'44',&
 Z'45', Z'46', Z'47', Z'48', Z'49', Z'4a', Z'4b', Z'4d',&
 Z'4e', Z'4f', Z'50', Z'51', Z'53', Z'54', Z'55', Z'57',&
 Z'58', Z'59', Z'5b', Z'5c', Z'5d', Z'5f', Z'60', Z'62',&
 Z'64', Z'65', Z'67', Z'68', Z'6a', Z'6c', Z'6e', Z'70',&
 Z'71', Z'73', Z'75', Z'77', Z'79', Z'7c', Z'7e', Z'80',&
 
 Z'53', Z'56', Z'59', Z'5c', Z'5f', Z'62', Z'64', Z'67',&   ! Blau */
 Z'69', Z'6c', Z'6e', Z'70', Z'73', Z'75', Z'77', Z'79',&
 Z'7b', Z'7e', Z'80', Z'81', Z'83', Z'85', Z'87', Z'89',&
 Z'8b', Z'8d', Z'8e', Z'90', Z'92', Z'93', Z'95', Z'97',&
 Z'98', Z'9a', Z'9b', Z'9d', Z'9e', Z'a0', Z'a1', Z'a3',&
 Z'a4', Z'a6', Z'a7', Z'a9', Z'aa', Z'ac', Z'ad', Z'af',&
 Z'b0', Z'b1', Z'b3', Z'b4', Z'b6', Z'b8', Z'b9', Z'bb',&
 Z'bc', Z'be', Z'bf', Z'c1', Z'c3', Z'c4', Z'c6', Z'c8',&
 Z'ca', Z'cc', Z'cd', Z'cf', Z'd1', Z'd3', Z'd5', Z'd7',&
 Z'd9', Z'dc', Z'de', Z'e0', Z'e2', Z'e5', Z'e7', Z'e9',&
 Z'ec', Z'ef', Z'f1', Z'f4', Z'f7', Z'fa', Z'fc', Z'ff',&
 Z'00', Z'04', Z'09', Z'0e', Z'12', Z'17', Z'1b', Z'1f',&
 Z'23', Z'27', Z'2b', Z'2e', Z'32', Z'36', Z'39', Z'3c',&
 Z'40', Z'43', Z'46', Z'49', Z'4c', Z'4f', Z'52', Z'55',&
 Z'58', Z'5b', Z'5f', Z'62', Z'65', Z'68', Z'6b', Z'6d',&
 Z'70', Z'72', Z'74', Z'76', Z'78', Z'7a', Z'00', Z'00',&
 Z'00', Z'00', Z'84', Z'85', Z'87', Z'89', Z'8b', Z'8e',&
 Z'90', Z'92', Z'95', Z'98', Z'9a', Z'9e', Z'a1', Z'a4',&
 Z'a8', Z'ab', Z'ae', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd',&
 Z'c1', Z'c4', Z'c7', Z'cb', Z'cf', Z'd2', Z'd6', Z'da',&
 Z'de', Z'e2', Z'e6', Z'eb', Z'ef', Z'f4', Z'f9', Z'fd',&
 Z'00', Z'02', Z'04', Z'06', Z'08', Z'0a', Z'0c', Z'0e',&
 Z'10', Z'12', Z'14', Z'15', Z'17', Z'19', Z'1a', Z'1c',&
 Z'1e', Z'1f', Z'21', Z'22', Z'23', Z'25', Z'26', Z'28',&
 Z'29', Z'2a', Z'2c', Z'2d', Z'2e', Z'2f', Z'31', Z'32',&
 Z'33', Z'34', Z'35', Z'36', Z'37', Z'39', Z'3a', Z'3b',&
 Z'3c', Z'3d', Z'3e', Z'3f', Z'40', Z'41', Z'42', Z'44',&
 Z'45', Z'46', Z'47', Z'48', Z'49', Z'4a', Z'4b', Z'4d',&
 Z'4e', Z'4f', Z'50', Z'51', Z'53', Z'54', Z'55', Z'57',&
 Z'58', Z'59', Z'5b', Z'5c', Z'5d', Z'5f', Z'60', Z'62',&
 Z'64', Z'65', Z'67', Z'68', Z'6a', Z'6c', Z'6e', Z'70',&
 Z'71', Z'73', Z'75', Z'77', Z'79', Z'7c', Z'7e', Z'80'/)
!                                                                ! Palette 6 */
Palette(6,:)=(/  Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&   ! Rot */
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'04', Z'09', Z'0d', Z'11', Z'15', Z'19', Z'1d',&
 Z'21', Z'24', Z'28', Z'2b', Z'2f', Z'32', Z'35', Z'39',&
 Z'3c', Z'3f', Z'42', Z'45', Z'47', Z'4a', Z'4d', Z'4f',&
 Z'52', Z'55', Z'57', Z'5a', Z'5c', Z'5f', Z'61', Z'63',&
 Z'66', Z'68', Z'6a', Z'6c', Z'6f', Z'71', Z'b9', Z'b9',&
 Z'b9', Z'b9', Z'7c', Z'7e', Z'81', Z'83', Z'85', Z'87',&
 Z'8a', Z'8c', Z'8e', Z'91', Z'93', Z'96', Z'98', Z'9b',&
 Z'9e', Z'a0', Z'a3', Z'a6', Z'a8', Z'ab', Z'ae', Z'b1',&
 Z'b4', Z'b8', Z'bb', Z'be', Z'c2', Z'c5', Z'c9', Z'cc',&
 Z'd0', Z'd4', Z'd8', Z'dc', Z'e0', Z'e4', Z'e9', Z'ed',&
 Z'33', Z'37', Z'3a', Z'3e', Z'41', Z'44', Z'47', Z'4a',&
 Z'4d', Z'50', Z'53', Z'55', Z'58', Z'5b', Z'5d', Z'60',&
 Z'62', Z'65', Z'67', Z'69', Z'6c', Z'6e', Z'70', Z'72',&
 Z'74', Z'76', Z'78', Z'7a', Z'7c', Z'7e', Z'80', Z'82',&
 Z'84', Z'86', Z'88', Z'89', Z'8b', Z'8d', Z'8f', Z'90',&
 Z'92', Z'94', Z'96', Z'97', Z'99', Z'9b', Z'9c', Z'9e',&
 Z'a0', Z'a2', Z'a3', Z'a5', Z'a7', Z'a9', Z'ab', Z'ac',&
 Z'ae', Z'b0', Z'b2', Z'b4', Z'b6', Z'b8', Z'ba', Z'bc',&
 Z'be', Z'c0', Z'c2', Z'c5', Z'c7', Z'c9', Z'cc', Z'ce',&
 Z'd0', Z'd3', Z'd5', Z'd8', Z'db', Z'de', Z'e0', Z'e3',&
 Z'e6', Z'e9', Z'ec', Z'ef', Z'f3', Z'f6', Z'f9', Z'fd',&
 
 Z'35', Z'38', Z'3c', Z'3f', Z'42', Z'45', Z'48', Z'4c',&   ! Gruen */
 Z'4e', Z'51', Z'54', Z'57', Z'5a', Z'5c', Z'5f', Z'61',&
 Z'64', Z'66', Z'68', Z'6b', Z'6d', Z'6f', Z'71', Z'74',&
 Z'76', Z'78', Z'7a', Z'7c', Z'7e', Z'80', Z'81', Z'83',&
 Z'85', Z'87', Z'89', Z'8b', Z'8c', Z'8e', Z'90', Z'92',&
 Z'93', Z'95', Z'97', Z'99', Z'9a', Z'9c', Z'9e', Z'9f',&
 Z'a1', Z'a3', Z'a5', Z'a6', Z'a8', Z'aa', Z'ac', Z'ae',&
 Z'af', Z'b1', Z'b3', Z'b5', Z'b7', Z'b9', Z'bb', Z'bd',&
 Z'bf', Z'c1', Z'c4', Z'c6', Z'c8', Z'ca', Z'cd', Z'cf',&
 Z'd1', Z'd4', Z'd7', Z'd9', Z'dc', Z'df', Z'e1', Z'e4',&
 Z'e7', Z'ea', Z'ed', Z'f0', Z'f4', Z'f7', Z'fa', Z'fe',&
 Z'00', Z'04', Z'09', Z'0d', Z'11', Z'15', Z'19', Z'1d',&
 Z'21', Z'24', Z'28', Z'2b', Z'2f', Z'32', Z'35', Z'39',&
 Z'3c', Z'3f', Z'42', Z'45', Z'47', Z'4a', Z'4d', Z'4f',&
 Z'52', Z'55', Z'57', Z'5a', Z'5c', Z'5f', Z'61', Z'63',&
 Z'66', Z'68', Z'6a', Z'6c', Z'6f', Z'71', Z'05', Z'05',&
 Z'05', Z'05', Z'7c', Z'7e', Z'81', Z'83', Z'85', Z'87',&
 Z'8a', Z'8c', Z'8e', Z'91', Z'93', Z'96', Z'98', Z'9b',&
 Z'9e', Z'a0', Z'a3', Z'a6', Z'a8', Z'ab', Z'ae', Z'b1',&
 Z'b4', Z'b8', Z'bb', Z'be', Z'c2', Z'c5', Z'c9', Z'cc',&
 Z'd0', Z'd4', Z'd8', Z'dc', Z'e0', Z'e4', Z'e9', Z'ed',&
 Z'1c', Z'1e', Z'21', Z'23', Z'25', Z'28', Z'2a', Z'2c',&
 Z'2e', Z'30', Z'32', Z'34', Z'36', Z'38', Z'3a', Z'3b',&
 Z'3d', Z'3f', Z'40', Z'42', Z'44', Z'45', Z'47', Z'48',&
 Z'4a', Z'4b', Z'4d', Z'4e', Z'4f', Z'51', Z'52', Z'54',&
 Z'55', Z'56', Z'57', Z'59', Z'5a', Z'5b', Z'5c', Z'5e',&
 Z'5f', Z'60', Z'61', Z'63', Z'64', Z'65', Z'66', Z'67',&
 Z'69', Z'6a', Z'6b', Z'6c', Z'6e', Z'6f', Z'70', Z'71',&
 Z'73', Z'74', Z'75', Z'77', Z'78', Z'7a', Z'7b', Z'7c',&
 Z'7e', Z'7f', Z'81', Z'83', Z'84', Z'86', Z'87', Z'89',&
 Z'8b', Z'8d', Z'8f', Z'90', Z'92', Z'94', Z'96', Z'98',&
 Z'9a', Z'9c', Z'9f', Z'a1', Z'a3', Z'a6', Z'a8', Z'aa',&
 
 Z'00', Z'03', Z'07', Z'0a', Z'0e', Z'11', Z'14', Z'17',&   ! Blau */
 Z'1b', Z'1e', Z'20', Z'23', Z'26', Z'29', Z'2c', Z'2e',&
 Z'31', Z'33', Z'36', Z'38', Z'3b', Z'3d', Z'3f', Z'41',&
 Z'44', Z'46', Z'48', Z'4a', Z'4c', Z'4e', Z'50', Z'52',&
 Z'54', Z'56', Z'58', Z'5a', Z'5c', Z'5d', Z'5f', Z'61',&
 Z'63', Z'65', Z'67', Z'68', Z'6a', Z'6c', Z'6e', Z'70',&
 Z'71', Z'73', Z'75', Z'77', Z'79', Z'7b', Z'7d', Z'7e',&
 Z'80', Z'82', Z'84', Z'86', Z'88', Z'8b', Z'8d', Z'8f',&
 Z'91', Z'93', Z'96', Z'98', Z'9a', Z'9d', Z'9f', Z'a2',&
 Z'a4', Z'a7', Z'aa', Z'ac', Z'af', Z'b2', Z'b5', Z'b8',&
 Z'bb', Z'be', Z'c1', Z'c5', Z'c8', Z'cc', Z'cf', Z'd3',&
 Z'00', Z'04', Z'09', Z'0d', Z'11', Z'15', Z'19', Z'1d',&
 Z'21', Z'24', Z'28', Z'2b', Z'2f', Z'32', Z'35', Z'39',&
 Z'3c', Z'3f', Z'42', Z'45', Z'47', Z'4a', Z'4d', Z'4f',&
 Z'52', Z'55', Z'57', Z'5a', Z'5c', Z'5f', Z'61', Z'63',&
 Z'66', Z'68', Z'6a', Z'6c', Z'6f', Z'71', Z'fd', Z'fd',&
 Z'fd', Z'fd', Z'7c', Z'7e', Z'81', Z'83', Z'85', Z'87',&
 Z'8a', Z'8c', Z'8e', Z'91', Z'93', Z'96', Z'98', Z'9b',&
 Z'9e', Z'a0', Z'a3', Z'a6', Z'a8', Z'ab', Z'ae', Z'b1',&
 Z'b4', Z'b8', Z'bb', Z'be', Z'c2', Z'c5', Z'c9', Z'cc',&
 Z'd0', Z'd4', Z'd8', Z'dc', Z'e0', Z'e4', Z'e9', Z'ed',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00'/)
!                                                                ! Palette 7 */
Palette(7,:)=(/  Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&   ! Rot */
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'04', Z'09', Z'0d', Z'11', Z'15', Z'19', Z'1d',&
 Z'21', Z'24', Z'28', Z'2b', Z'2f', Z'32', Z'35', Z'39',&
 Z'3c', Z'3f', Z'42', Z'45', Z'47', Z'4a', Z'4d', Z'4f',&
 Z'52', Z'55', Z'57', Z'5a', Z'5c', Z'5f', Z'61', Z'63',&
 Z'66', Z'68', Z'6a', Z'6c', Z'6f', Z'71', Z'73', Z'00',&
 Z'00', Z'7a', Z'7c', Z'7e', Z'81', Z'83', Z'85', Z'87',&
 Z'8a', Z'8c', Z'8e', Z'91', Z'93', Z'96', Z'98', Z'9b',&
 Z'9e', Z'a0', Z'a3', Z'a6', Z'a8', Z'ab', Z'ae', Z'b1',&
 Z'b4', Z'b8', Z'bb', Z'be', Z'c2', Z'c5', Z'c9', Z'cc',&
 Z'd0', Z'd4', Z'd8', Z'dc', Z'e0', Z'e4', Z'e9', Z'ed',&
 Z'33', Z'37', Z'3a', Z'3e', Z'41', Z'44', Z'47', Z'4a',&
 Z'4d', Z'50', Z'53', Z'55', Z'58', Z'5b', Z'5d', Z'60',&
 Z'62', Z'65', Z'67', Z'69', Z'6c', Z'6e', Z'70', Z'72',&
 Z'74', Z'76', Z'78', Z'7a', Z'7c', Z'7e', Z'80', Z'82',&
 Z'84', Z'86', Z'88', Z'89', Z'8b', Z'8d', Z'8f', Z'90',&
 Z'92', Z'94', Z'96', Z'97', Z'99', Z'9b', Z'9c', Z'9e',&
 Z'a0', Z'a2', Z'a3', Z'a5', Z'a7', Z'a9', Z'ab', Z'ac',&
 Z'ae', Z'b0', Z'b2', Z'b4', Z'b6', Z'b8', Z'ba', Z'bc',&
 Z'be', Z'c0', Z'c2', Z'c5', Z'c7', Z'c9', Z'cc', Z'ce',&
 Z'd0', Z'd3', Z'd5', Z'd8', Z'db', Z'de', Z'e0', Z'e3',&
 Z'e6', Z'e9', Z'ec', Z'ef', Z'f3', Z'f6', Z'f9', Z'fd',&
 
 Z'35', Z'34', Z'33', Z'32', Z'31', Z'30', Z'30', Z'2f',&   ! Gruen */
 Z'2e', Z'2d', Z'2d', Z'2c', Z'2b', Z'2a', Z'2a', Z'29',&
 Z'28', Z'28', Z'27', Z'27', Z'26', Z'25', Z'25', Z'24',&
 Z'24', Z'23', Z'23', Z'22', Z'22', Z'21', Z'21', Z'20',&
 Z'20', Z'1f', Z'1f', Z'1e', Z'1e', Z'1d', Z'1d', Z'1c',&
 Z'1c', Z'1b', Z'1b', Z'1a', Z'1a', Z'1a', Z'19', Z'19',&
 Z'18', Z'18', Z'17', Z'17', Z'16', Z'16', Z'15', Z'15',&
 Z'14', Z'14', Z'13', Z'13', Z'12', Z'12', Z'11', Z'11',&
 Z'10', Z'10', Z'0f', Z'0f', Z'0e', Z'0d', Z'0d', Z'0c',&
 Z'0b', Z'0b', Z'0a', Z'09', Z'09', Z'08', Z'07', Z'06',&
 Z'06', Z'05', Z'04', Z'03', Z'02', Z'01', Z'01', Z'00',&
 Z'00', Z'04', Z'09', Z'0d', Z'11', Z'15', Z'19', Z'1d',&
 Z'21', Z'24', Z'28', Z'2b', Z'2f', Z'32', Z'35', Z'39',&
 Z'3c', Z'3f', Z'42', Z'45', Z'47', Z'4a', Z'4d', Z'4f',&
 Z'52', Z'55', Z'57', Z'5a', Z'5c', Z'5f', Z'61', Z'63',&
 Z'66', Z'68', Z'6a', Z'6c', Z'6f', Z'71', Z'73', Z'fe',&
 Z'fe', Z'7a', Z'7c', Z'7e', Z'81', Z'83', Z'85', Z'87',&
 Z'8a', Z'8c', Z'8e', Z'91', Z'93', Z'96', Z'98', Z'9b',&
 Z'9e', Z'a0', Z'a3', Z'a6', Z'a8', Z'ab', Z'ae', Z'b1',&
 Z'b4', Z'b8', Z'bb', Z'be', Z'c2', Z'c5', Z'c9', Z'cc',&
 Z'd0', Z'd4', Z'd8', Z'dc', Z'e0', Z'e4', Z'e9', Z'ed',&
 Z'1c', Z'1e', Z'21', Z'23', Z'25', Z'28', Z'2a', Z'2c',&
 Z'2e', Z'30', Z'32', Z'34', Z'36', Z'38', Z'3a', Z'3b',&
 Z'3d', Z'3f', Z'40', Z'42', Z'44', Z'45', Z'47', Z'48',&
 Z'4a', Z'4b', Z'4d', Z'4e', Z'4f', Z'51', Z'52', Z'54',&
 Z'55', Z'56', Z'57', Z'59', Z'5a', Z'5b', Z'5c', Z'5e',&
 Z'5f', Z'60', Z'61', Z'63', Z'64', Z'65', Z'66', Z'67',&
 Z'69', Z'6a', Z'6b', Z'6c', Z'6e', Z'6f', Z'70', Z'71',&
 Z'73', Z'74', Z'75', Z'77', Z'78', Z'7a', Z'7b', Z'7c',&
 Z'7e', Z'7f', Z'81', Z'83', Z'84', Z'86', Z'87', Z'89',&
 Z'8b', Z'8d', Z'8f', Z'90', Z'92', Z'94', Z'96', Z'98',&
 Z'9a', Z'9c', Z'9f', Z'a1', Z'a3', Z'a6', Z'a8', Z'aa',&
 
 Z'00', Z'04', Z'08', Z'0c', Z'11', Z'15', Z'18', Z'1c',&   ! Blau */
 Z'20', Z'24', Z'27', Z'2b', Z'2e', Z'31', Z'35', Z'38',&
 Z'3b', Z'3e', Z'41', Z'44', Z'47', Z'49', Z'4c', Z'4f',&
 Z'52', Z'54', Z'57', Z'59', Z'5c', Z'5e', Z'61', Z'63',&
 Z'65', Z'68', Z'6a', Z'6c', Z'6e', Z'71', Z'73', Z'75',&
 Z'77', Z'79', Z'7b', Z'7e', Z'80', Z'82', Z'84', Z'86',&
 Z'89', Z'8b', Z'8d', Z'8f', Z'91', Z'94', Z'96', Z'98',&
 Z'9b', Z'9d', Z'9f', Z'a2', Z'a4', Z'a7', Z'a9', Z'ac',&
 Z'af', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd', Z'c0', Z'c3',&
 Z'c6', Z'c9', Z'cc', Z'cf', Z'd3', Z'd6', Z'da', Z'de',&
 Z'e1', Z'e5', Z'e9', Z'ed', Z'f1', Z'f5', Z'f9', Z'fe',&
 Z'00', Z'04', Z'09', Z'0d', Z'11', Z'15', Z'19', Z'1d',&
 Z'21', Z'24', Z'28', Z'2b', Z'2f', Z'32', Z'35', Z'39',&
 Z'3c', Z'3f', Z'42', Z'45', Z'47', Z'4a', Z'4d', Z'4f',&
 Z'52', Z'55', Z'57', Z'5a', Z'5c', Z'5f', Z'61', Z'63',&
 Z'66', Z'68', Z'6a', Z'6c', Z'6f', Z'71', Z'73', Z'00',&
 Z'00', Z'7a', Z'7c', Z'7e', Z'81', Z'83', Z'85', Z'87',&
 Z'8a', Z'8c', Z'8e', Z'91', Z'93', Z'96', Z'98', Z'9b',&
 Z'9e', Z'a0', Z'a3', Z'a6', Z'a8', Z'ab', Z'ae', Z'b1',&
 Z'b4', Z'b8', Z'bb', Z'be', Z'c2', Z'c5', Z'c9', Z'cc',&
 Z'd0', Z'd4', Z'd8', Z'dc', Z'e0', Z'e4', Z'e9', Z'ed',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00'/)
!                                                                ! Palette 8 */
Palette(8,:)=(/  Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&   ! Rot */
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'1a', Z'1d', Z'21', Z'24', Z'28', Z'2b', Z'2e', Z'31',&
 Z'34', Z'37', Z'3a', Z'3d', Z'3f', Z'42', Z'44', Z'47',&
 Z'49', Z'4c', Z'4e', Z'50', Z'53', Z'55', Z'57', Z'59',&
 Z'5b', Z'5d', Z'5f', Z'61', Z'63', Z'65', Z'67', Z'69',&
 Z'6b', Z'6d', Z'6e', Z'70', Z'72', Z'74', Z'f6', Z'f6',&
 Z'f6', Z'f6', Z'7d', Z'7e', Z'80', Z'82', Z'84', Z'86',&
 Z'87', Z'89', Z'8b', Z'8d', Z'8f', Z'91', Z'93', Z'95',&
 Z'97', Z'99', Z'9b', Z'9e', Z'a0', Z'a2', Z'a4', Z'a7',&
 Z'a9', Z'ac', Z'ae', Z'b1', Z'b4', Z'b6', Z'b9', Z'bc',&
 Z'bf', Z'c2', Z'c5', Z'c9', Z'cc', Z'cf', Z'd3', Z'd6',&
 Z'33', Z'37', Z'3a', Z'3e', Z'41', Z'44', Z'47', Z'4a',&
 Z'4d', Z'50', Z'53', Z'56', Z'58', Z'5b', Z'5e', Z'60',&
 Z'63', Z'65', Z'67', Z'6a', Z'6c', Z'6e', Z'70', Z'73',&
 Z'75', Z'77', Z'79', Z'7b', Z'7d', Z'7f', Z'81', Z'82',&
 Z'84', Z'86', Z'88', Z'8a', Z'8c', Z'8d', Z'8f', Z'91',&
 Z'93', Z'94', Z'96', Z'98', Z'99', Z'9b', Z'9d', Z'9f',&
 Z'a0', Z'a2', Z'a4', Z'a6', Z'a7', Z'a9', Z'ab', Z'ad',&
 Z'af', Z'b1', Z'b3', Z'b5', Z'b7', Z'b9', Z'bb', Z'bd',&
 Z'bf', Z'c1', Z'c3', Z'c5', Z'c8', Z'ca', Z'cc', Z'cf',&
 Z'd1', Z'd4', Z'd6', Z'd9', Z'dc', Z'de', Z'e1', Z'e4',&
 Z'e7', Z'ea', Z'ed', Z'f0', Z'f4', Z'f7', Z'fa', Z'fe',&
 
 Z'4e', Z'51', Z'54', Z'57', Z'5a', Z'5c', Z'5f', Z'62',&   ! Gruen */
 Z'64', Z'67', Z'69', Z'6c', Z'6e', Z'70', Z'72', Z'75',&
 Z'77', Z'79', Z'7b', Z'7d', Z'7f', Z'81', Z'83', Z'85',&
 Z'86', Z'88', Z'8a', Z'8c', Z'8d', Z'8f', Z'91', Z'92',&
 Z'94', Z'96', Z'97', Z'99', Z'9a', Z'9c', Z'9d', Z'9f',&
 Z'a0', Z'a2', Z'a3', Z'a5', Z'a6', Z'a8', Z'a9', Z'ab',&
 Z'ac', Z'ae', Z'af', Z'b1', Z'b3', Z'b4', Z'b6', Z'b7',&
 Z'b9', Z'ba', Z'bc', Z'be', Z'c0', Z'c1', Z'c3', Z'c5',&
 Z'c7', Z'c9', Z'ca', Z'cc', Z'ce', Z'd0', Z'd2', Z'd4',&
 Z'd7', Z'd9', Z'db', Z'dd', Z'e0', Z'e2', Z'e5', Z'e7',&
 Z'ea', Z'ec', Z'ef', Z'f2', Z'f5', Z'f7', Z'fa', Z'fd',&
 Z'22', Z'25', Z'29', Z'2c', Z'2f', Z'32', Z'35', Z'38',&
 Z'3b', Z'3e', Z'41', Z'43', Z'46', Z'48', Z'4b', Z'4d',&
 Z'50', Z'52', Z'54', Z'56', Z'58', Z'5b', Z'5d', Z'5f',&
 Z'61', Z'63', Z'65', Z'66', Z'68', Z'6a', Z'6c', Z'6e',&
 Z'6f', Z'71', Z'73', Z'75', Z'76', Z'78', Z'da', Z'da',&
 Z'da', Z'da', Z'81', Z'82', Z'84', Z'86', Z'88', Z'89',&
 Z'8b', Z'8d', Z'8f', Z'90', Z'92', Z'94', Z'96', Z'98',&
 Z'9a', Z'9c', Z'9e', Z'a0', Z'a2', Z'a5', Z'a7', Z'a9',&
 Z'ac', Z'ae', Z'b0', Z'b3', Z'b6', Z'b8', Z'bb', Z'be',&
 Z'c1', Z'c4', Z'c7', Z'ca', Z'cd', Z'd0', Z'd3', Z'd7',&
 Z'1c', Z'1b', Z'1b', Z'1a', Z'1a', Z'1a', Z'19', Z'19',&
 Z'18', Z'18', Z'17', Z'17', Z'17', Z'16', Z'16', Z'16',&
 Z'15', Z'15', Z'15', Z'14', Z'14', Z'14', Z'13', Z'13',&
 Z'13', Z'12', Z'12', Z'12', Z'12', Z'11', Z'11', Z'11',&
 Z'11', Z'10', Z'10', Z'10', Z'10', Z'0f', Z'0f', Z'0f',&
 Z'0f', Z'0e', Z'0e', Z'0e', Z'0e', Z'0d', Z'0d', Z'0d',&
 Z'0d', Z'0c', Z'0c', Z'0c', Z'0c', Z'0b', Z'0b', Z'0b',&
 Z'0b', Z'0a', Z'0a', Z'0a', Z'0a', Z'09', Z'09', Z'09',&
 Z'08', Z'08', Z'08', Z'08', Z'07', Z'07', Z'07', Z'06',&
 Z'06', Z'06', Z'05', Z'05', Z'04', Z'04', Z'04', Z'03',&
 Z'03', Z'02', Z'02', Z'02', Z'01', Z'01', Z'00', Z'00',&
 
 Z'00', Z'03', Z'07', Z'0b', Z'0e', Z'12', Z'15', Z'18',&   ! Blau */
 Z'1c', Z'1f', Z'22', Z'25', Z'28', Z'2b', Z'2e', Z'30',&
 Z'33', Z'36', Z'38', Z'3b', Z'3d', Z'40', Z'42', Z'44',&
 Z'47', Z'49', Z'4b', Z'4d', Z'50', Z'52', Z'54', Z'56',&
 Z'58', Z'5a', Z'5c', Z'5e', Z'60', Z'62', Z'64', Z'66',&
 Z'67', Z'69', Z'6b', Z'6d', Z'6f', Z'71', Z'73', Z'75',&
 Z'76', Z'78', Z'7a', Z'7c', Z'7e', Z'80', Z'82', Z'84',&
 Z'86', Z'88', Z'8a', Z'8c', Z'8f', Z'91', Z'93', Z'95',&
 Z'98', Z'9a', Z'9c', Z'9f', Z'a1', Z'a4', Z'a6', Z'a9',&
 Z'ac', Z'ae', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd', Z'c0',&
 Z'c4', Z'c7', Z'ca', Z'ce', Z'd1', Z'd5', Z'd9', Z'dc',&
 Z'1a', Z'1d', Z'21', Z'24', Z'28', Z'2b', Z'2e', Z'31',&
 Z'34', Z'37', Z'3a', Z'3d', Z'3f', Z'42', Z'45', Z'47',&
 Z'4a', Z'4c', Z'4e', Z'51', Z'53', Z'55', Z'57', Z'59',&
 Z'5b', Z'5d', Z'5f', Z'61', Z'63', Z'65', Z'67', Z'69',&
 Z'6b', Z'6d', Z'6f', Z'70', Z'72', Z'74', Z'00', Z'00',&
 Z'00', Z'00', Z'7d', Z'7f', Z'80', Z'82', Z'84', Z'86',&
 Z'88', Z'8a', Z'8b', Z'8d', Z'8f', Z'91', Z'93', Z'95',&
 Z'97', Z'9a', Z'9c', Z'9e', Z'a0', Z'a2', Z'a5', Z'a7',&
 Z'aa', Z'ac', Z'af', Z'b1', Z'b4', Z'b7', Z'ba', Z'bd',&
 Z'c0', Z'c3', Z'c6', Z'c9', Z'cc', Z'd0', Z'd3', Z'd7',&
 Z'00', Z'03', Z'06', Z'0a', Z'0d', Z'10', Z'13', Z'16',&
 Z'19', Z'1b', Z'1e', Z'21', Z'23', Z'26', Z'28', Z'2b',&
 Z'2d', Z'30', Z'32', Z'34', Z'36', Z'38', Z'3b', Z'3d',&
 Z'3f', Z'41', Z'43', Z'45', Z'47', Z'48', Z'4a', Z'4c',&
 Z'4e', Z'50', Z'51', Z'53', Z'55', Z'57', Z'58', Z'5a',&
 Z'5c', Z'5d', Z'5f', Z'61', Z'62', Z'64', Z'66', Z'67',&
 Z'69', Z'6b', Z'6c', Z'6e', Z'70', Z'72', Z'73', Z'75',&
 Z'77', Z'79', Z'7b', Z'7c', Z'7e', Z'80', Z'82', Z'84',&
 Z'86', Z'88', Z'8b', Z'8d', Z'8f', Z'91', Z'93', Z'96',&
 Z'98', Z'9b', Z'9d', Z'a0', Z'a2', Z'a5', Z'a8', Z'aa',&
 Z'ad', Z'b0', Z'b3', Z'b6', Z'b9', Z'bd', Z'c0', Z'c3'/)
!                                                                ! Palette 9 */
palette(9,:)=(/  Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&   ! Rot */
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'04', Z'09', Z'0d', Z'11', Z'15', Z'19', Z'1d',&
 Z'21', Z'24', Z'28', Z'2b', Z'2f', Z'32', Z'35', Z'39',&
 Z'3c', Z'3f', Z'42', Z'45', Z'47', Z'4a', Z'4d', Z'4f',&
 Z'52', Z'55', Z'57', Z'5a', Z'5c', Z'5f', Z'61', Z'63',&
 Z'66', Z'68', Z'6a', Z'6c', Z'6f', Z'71', Z'fd', Z'fd',&
 Z'fd', Z'fd', Z'7c', Z'7e', Z'81', Z'83', Z'85', Z'87',&
 Z'8a', Z'8c', Z'8e', Z'91', Z'93', Z'96', Z'98', Z'9b',&
 Z'9e', Z'a0', Z'a3', Z'a6', Z'a8', Z'ab', Z'ae', Z'b1',&
 Z'b4', Z'b8', Z'bb', Z'be', Z'c2', Z'c5', Z'c9', Z'cc',&
 Z'd0', Z'd4', Z'd8', Z'dc', Z'e0', Z'e4', Z'e9', Z'ed',&
 Z'59', Z'5c', Z'5e', Z'60', Z'62', Z'64', Z'66', Z'68',&
 Z'6a', Z'6b', Z'6d', Z'6f', Z'71', Z'72', Z'74', Z'76',&
 Z'77', Z'79', Z'7a', Z'7c', Z'7d', Z'7f', Z'80', Z'81',&
 Z'83', Z'84', Z'85', Z'86', Z'88', Z'89', Z'8a', Z'8b',&
 Z'8d', Z'8e', Z'8f', Z'90', Z'91', Z'92', Z'93', Z'94',&
 Z'96', Z'97', Z'98', Z'99', Z'9a', Z'9b', Z'9c', Z'9d',&
 Z'9e', Z'9f', Z'a1', Z'a2', Z'a3', Z'a4', Z'a5', Z'a6',&
 Z'a7', Z'a9', Z'aa', Z'ab', Z'ac', Z'ae', Z'af', Z'b0',&
 Z'b1', Z'b3', Z'b4', Z'b6', Z'b7', Z'b9', Z'ba', Z'bc',&
 Z'bd', Z'bf', Z'c0', Z'c2', Z'c4', Z'c5', Z'c7', Z'c9',&
 Z'cb', Z'cd', Z'cf', Z'd1', Z'd3', Z'd5', Z'd7', Z'd9',&
 
 Z'35', Z'38', Z'3c', Z'3f', Z'42', Z'45', Z'48', Z'4c',&   ! Gruen */
 Z'4e', Z'51', Z'54', Z'57', Z'5a', Z'5c', Z'5f', Z'61',&
 Z'64', Z'66', Z'68', Z'6b', Z'6d', Z'6f', Z'71', Z'74',&
 Z'76', Z'78', Z'7a', Z'7c', Z'7e', Z'80', Z'81', Z'83',&
 Z'85', Z'87', Z'89', Z'8b', Z'8c', Z'8e', Z'90', Z'92',&
 Z'93', Z'95', Z'97', Z'99', Z'9a', Z'9c', Z'9e', Z'9f',&
 Z'a1', Z'a3', Z'a5', Z'a6', Z'a8', Z'aa', Z'ac', Z'ae',&
 Z'af', Z'b1', Z'b3', Z'b5', Z'b7', Z'b9', Z'bb', Z'bd',&
 Z'bf', Z'c1', Z'c4', Z'c6', Z'c8', Z'ca', Z'cd', Z'cf',&
 Z'd1', Z'd4', Z'd7', Z'd9', Z'dc', Z'df', Z'e1', Z'e4',&
 Z'e7', Z'ea', Z'ed', Z'f0', Z'f4', Z'f7', Z'fa', Z'fe',&
 Z'00', Z'04', Z'09', Z'0d', Z'11', Z'15', Z'19', Z'1d',&
 Z'21', Z'24', Z'28', Z'2b', Z'2f', Z'32', Z'35', Z'39',&
 Z'3c', Z'3f', Z'42', Z'45', Z'47', Z'4a', Z'4d', Z'4f',&
 Z'52', Z'55', Z'57', Z'5a', Z'5c', Z'5f', Z'61', Z'63',&
 Z'66', Z'68', Z'6a', Z'6c', Z'6f', Z'71', Z'aa', Z'aa',&
 Z'aa', Z'aa', Z'7c', Z'7e', Z'81', Z'83', Z'85', Z'87',&
 Z'8a', Z'8c', Z'8e', Z'91', Z'93', Z'96', Z'98', Z'9b',&
 Z'9e', Z'a0', Z'a3', Z'a6', Z'a8', Z'ab', Z'ae', Z'b1',&
 Z'b4', Z'b8', Z'bb', Z'be', Z'c2', Z'c5', Z'c9', Z'cc',&
 Z'd0', Z'd4', Z'd8', Z'dc', Z'e0', Z'e4', Z'e9', Z'ed',&
 Z'00', Z'00', Z'01', Z'02', Z'02', Z'03', Z'03', Z'04',&
 Z'04', Z'05', Z'05', Z'06', Z'07', Z'07', Z'07', Z'08',&
 Z'08', Z'09', Z'09', Z'0a', Z'0a', Z'0b', Z'0b', Z'0b',&
 Z'0c', Z'0c', Z'0d', Z'0d', Z'0d', Z'0e', Z'0e', Z'0e',&
 Z'0f', Z'0f', Z'0f', Z'10', Z'10', Z'10', Z'11', Z'11',&
 Z'11', Z'12', Z'12', Z'12', Z'13', Z'13', Z'13', Z'14',&
 Z'14', Z'14', Z'15', Z'15', Z'15', Z'16', Z'16', Z'16',&
 Z'17', Z'17', Z'17', Z'18', Z'18', Z'18', Z'19', Z'19',&
 Z'1a', Z'1a', Z'1a', Z'1b', Z'1b', Z'1c', Z'1c', Z'1d',&
 Z'1d', Z'1e', Z'1e', Z'1e', Z'1f', Z'20', Z'20', Z'21',&
 Z'21', Z'22', Z'22', Z'23', Z'23', Z'24', Z'25', Z'25',&
 
 Z'00', Z'03', Z'07', Z'0a', Z'0e', Z'11', Z'14', Z'17',&   ! Blau */
 Z'1b', Z'1e', Z'20', Z'23', Z'26', Z'29', Z'2c', Z'2e',&
 Z'31', Z'33', Z'36', Z'38', Z'3b', Z'3d', Z'3f', Z'41',&
 Z'44', Z'46', Z'48', Z'4a', Z'4c', Z'4e', Z'50', Z'52',&
 Z'54', Z'56', Z'58', Z'5a', Z'5c', Z'5d', Z'5f', Z'61',&
 Z'63', Z'65', Z'67', Z'68', Z'6a', Z'6c', Z'6e', Z'70',&
 Z'71', Z'73', Z'75', Z'77', Z'79', Z'7b', Z'7d', Z'7e',&
 Z'80', Z'82', Z'84', Z'86', Z'88', Z'8b', Z'8d', Z'8f',&
 Z'91', Z'93', Z'96', Z'98', Z'9a', Z'9d', Z'9f', Z'a2',&
 Z'a4', Z'a7', Z'aa', Z'ac', Z'af', Z'b2', Z'b5', Z'b8',&
 Z'bb', Z'be', Z'c1', Z'c5', Z'c8', Z'cc', Z'cf', Z'd3',&
 Z'00', Z'04', Z'09', Z'0d', Z'11', Z'15', Z'19', Z'1d',&
 Z'21', Z'24', Z'28', Z'2b', Z'2f', Z'32', Z'35', Z'39',&
 Z'3c', Z'3f', Z'42', Z'45', Z'47', Z'4a', Z'4d', Z'4f',&
 Z'52', Z'55', Z'57', Z'5a', Z'5c', Z'5f', Z'61', Z'63',&
 Z'66', Z'68', Z'6a', Z'6c', Z'6f', Z'71', Z'00', Z'00',&
 Z'00', Z'00', Z'7c', Z'7e', Z'81', Z'83', Z'85', Z'87',&
 Z'8a', Z'8c', Z'8e', Z'91', Z'93', Z'96', Z'98', Z'9b',&
 Z'9e', Z'a0', Z'a3', Z'a6', Z'a8', Z'ab', Z'ae', Z'b1',&
 Z'b4', Z'b8', Z'bb', Z'be', Z'c2', Z'c5', Z'c9', Z'cc',&
 Z'd0', Z'd4', Z'd8', Z'dc', Z'e0', Z'e4', Z'e9', Z'ed',&
 Z'3f', Z'43', Z'46', Z'49', Z'4c', Z'4f', Z'52', Z'55',&
 Z'58', Z'5b', Z'5d', Z'60', Z'62', Z'65', Z'67', Z'6a',&
 Z'6c', Z'6e', Z'71', Z'73', Z'75', Z'77', Z'79', Z'7b',&
 Z'7d', Z'7f', Z'81', Z'83', Z'85', Z'87', Z'89', Z'8a',&
 Z'8c', Z'8e', Z'90', Z'91', Z'93', Z'95', Z'96', Z'98',&
 Z'9a', Z'9b', Z'9d', Z'9f', Z'a0', Z'a2', Z'a3', Z'a5',&
 Z'a7', Z'a8', Z'aa', Z'ac', Z'ad', Z'af', Z'b1', Z'b3',&
 Z'b4', Z'b6', Z'b8', Z'ba', Z'bc', Z'be', Z'c0', Z'c2',&
 Z'c4', Z'c6', Z'c8', Z'ca', Z'cc', Z'ce', Z'd0', Z'd3',&
 Z'd5', Z'd7', Z'da', Z'dc', Z'df', Z'e2', Z'e4', Z'e7',&
 Z'ea', Z'ed', Z'f0', Z'f3', Z'f6', Z'f9', Z'fc', Z'ff'/)
!                                                                ! Palette 10 */
palette(10,:)=(/  Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&   ! Rot */
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'02', Z'04', Z'06', Z'08', Z'0a', Z'0c', Z'0e',&
 Z'10', Z'11', Z'13', Z'15', Z'17', Z'19', Z'1a', Z'1c',&
 Z'1e', Z'1f', Z'21', Z'23', Z'24', Z'26', Z'27', Z'29',&
 Z'2a', Z'2c', Z'2d', Z'2e', Z'30', Z'31', Z'33', Z'34',&
 Z'35', Z'37', Z'38', Z'39', Z'3a', Z'3c', Z'3d', Z'3e',&
 Z'3f', Z'40', Z'42', Z'43', Z'44', Z'45', Z'46', Z'47',&
 Z'48', Z'4a', Z'4b', Z'4c', Z'4d', Z'4e', Z'4f', Z'50',&
 Z'51', Z'52', Z'53', Z'54', Z'55', Z'56', Z'57', Z'58',&
 Z'5a', Z'5b', Z'5c', Z'5d', Z'5e', Z'5f', Z'60', Z'61',&
 Z'62', Z'63', Z'64', Z'65', Z'66', Z'67', Z'68', Z'6a',&
 Z'6b', Z'6c', Z'6d', Z'6e', Z'6f', Z'70', Z'72', Z'73',&
 Z'74', Z'75', Z'76', Z'78', Z'79', Z'7a', Z'7b', Z'7d',&
 Z'7e', Z'7f', Z'81', Z'82', Z'84', Z'85', Z'86', Z'88',&
 Z'89', Z'8b', Z'8c', Z'8e', Z'8f', Z'91', Z'93', Z'94',&
 Z'96', Z'98', Z'99', Z'9b', Z'9d', Z'9f', Z'a1', Z'a2',&
 Z'a4', Z'a6', Z'a8', Z'aa', Z'ac', Z'ae', Z'b0', Z'b2',&
 Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff',&
 Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff',&
 Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff',&
 Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff',&
 Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff',&
 Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff',&
 Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff',&
 Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff',&
 
 Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff',&   ! Gruen */
 Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff',&
 Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff',&
 Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff',&
 Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff',&
 Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff',&
 Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff',&
 Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff', Z'ff',&
 Z'00', Z'02', Z'04', Z'06', Z'08', Z'0a', Z'0c', Z'0e',&
 Z'10', Z'11', Z'13', Z'15', Z'17', Z'19', Z'1a', Z'1c',&
 Z'1e', Z'1f', Z'21', Z'23', Z'24', Z'26', Z'27', Z'29',&
 Z'2a', Z'2c', Z'2d', Z'2e', Z'30', Z'31', Z'33', Z'34',&
 Z'35', Z'37', Z'38', Z'39', Z'3a', Z'3c', Z'3d', Z'3e',&
 Z'3f', Z'40', Z'42', Z'43', Z'44', Z'45', Z'46', Z'47',&
 Z'48', Z'4a', Z'4b', Z'4c', Z'4d', Z'4e', Z'4f', Z'50',&
 Z'51', Z'52', Z'53', Z'54', Z'55', Z'56', Z'57', Z'58',&
 Z'5a', Z'5b', Z'5c', Z'5d', Z'5e', Z'5f', Z'60', Z'61',&
 Z'62', Z'63', Z'64', Z'65', Z'66', Z'67', Z'68', Z'6a',&
 Z'6b', Z'6c', Z'6d', Z'6e', Z'6f', Z'70', Z'72', Z'73',&
 Z'74', Z'75', Z'76', Z'78', Z'79', Z'7a', Z'7b', Z'7d',&
 Z'7e', Z'7f', Z'81', Z'82', Z'84', Z'85', Z'86', Z'88',&
 Z'89', Z'8b', Z'8c', Z'8e', Z'8f', Z'91', Z'93', Z'94',&
 Z'96', Z'98', Z'99', Z'9b', Z'9d', Z'9f', Z'a1', Z'a2',&
 Z'a4', Z'a6', Z'a8', Z'aa', Z'ac', Z'ae', Z'b0', Z'b2',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&   ! Blau */
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'02', Z'04', Z'06', Z'08', Z'0a', Z'0c', Z'0e',&
 Z'10', Z'11', Z'13', Z'15', Z'17', Z'19', Z'1a', Z'1c',&
 Z'1e', Z'1f', Z'21', Z'23', Z'24', Z'26', Z'27', Z'29',&
 Z'2a', Z'2c', Z'2d', Z'2e', Z'30', Z'31', Z'33', Z'34',&
 Z'35', Z'37', Z'38', Z'39', Z'3a', Z'3c', Z'3d', Z'3e',&
 Z'3f', Z'40', Z'42', Z'43', Z'44', Z'45', Z'46', Z'47',&
 Z'48', Z'4a', Z'4b', Z'4c', Z'4d', Z'4e', Z'4f', Z'50',&
 Z'51', Z'52', Z'53', Z'54', Z'55', Z'56', Z'57', Z'58',&
 Z'5a', Z'5b', Z'5c', Z'5d', Z'5e', Z'5f', Z'60', Z'61',&
 Z'62', Z'63', Z'64', Z'65', Z'66', Z'67', Z'68', Z'6a',&
 Z'6b', Z'6c', Z'6d', Z'6e', Z'6f', Z'70', Z'72', Z'73',&
 Z'74', Z'75', Z'76', Z'78', Z'79', Z'7a', Z'7b', Z'7d',&
 Z'7e', Z'7f', Z'81', Z'82', Z'84', Z'85', Z'86', Z'88',&
 Z'89', Z'8b', Z'8c', Z'8e', Z'8f', Z'91', Z'93', Z'94',&
 Z'96', Z'98', Z'99', Z'9b', Z'9d', Z'9f', Z'a1', Z'a2',&
 Z'a4', Z'a6', Z'a8', Z'aa', Z'ac', Z'ae', Z'b0', Z'b2',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00',&
 Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00', Z'00' /)
 
 !=================================================== Palette 11: linear
 
  !this palette divides the range in 12 linear parts
  
  do i=1, 256
    !-----------------_Flameball palette
    intensity = real(i-1)/255.
    if ((intensity>0.0).and.(intensity<0.25)) then
      palette(12, i)     = 0 !red  
      palette(12, i+256) = 0 !green
      palette(12, i+512) = nint(intensity*255.0*4.0) !blue
    elseif ((intensity>0.25).and.(intensity<0.50)) then
      palette(12, i)     = nint((intensity-0.25)*4.0*255.0)!red  
      palette(12, i+256) = 0!green
      palette(12, i+512) = 255-nint((intensity-0.25)*4.0*255.0)!blue
    elseif ((intensity>0.50).and.(intensity<0.75)) then
      palette(12, i)     = 255 !red  
      palette(12, i+256) = nint ((intensity-0.5)*4.0*255.0)!green
      palette(12, i+512) = 0 !blue
    elseif ((intensity>0.75).and.(intensity<1.00)) then
      palette(12, i)     = 255!red  
      palette(12, i+256) = 255!green
      palette(12, i+512) = nint((intensity-0.75)*4.0*255.0)!blue   
  end if
  
  !------------------ Matlab "spring" palette
  palette(11, i)     = nint(intensity*255.) !red  
  palette(11, i+256) = nint(255.-intensity*255.) !green
  palette(11, i+512) = 255 !blue
  
  !------------------- Matlab standard palette
  if ((intensity<=0.20)) then
      palette(13, i)     = 0 !red  
      palette(13, i+256) = 0 !green
      palette(13, i+512) = 143 + nint (intensity*(255.-143.)*5.)
    elseif ((intensity>0.20).and.(intensity<=0.40)) then
      palette(13, i)     = 0!red  
      palette(13, i+256) = nint ((intensity-0.20)*255.*5. )!green
      palette(13, i+512) = 255!blue
    elseif ((intensity>0.40).and.(intensity<=0.60)) then
      palette(13, i)     = nint ((intensity-0.40)*255.*5. ) !red  
      palette(13, i+256) = 255!green
      palette(13, i+512) = 255-nint ( (intensity-0.40 )*255.*5. ) !blue
    elseif ((intensity>0.60).and.(intensity<=0.80)) then
      palette(13, i)     = 255!red  
      palette(13, i+256) = 255-nint ( (intensity-0.60 )*255.*5. )!green
      palette(13, i+512) = 0!blue   
    elseif (intensity>=0.80) then
      palette(13, i)     = 255 - nint ((intensity-0.80)*(255.-128.)*5. )!red  
      palette(13, i+256) = 0!green
      palette(13, i+512) = 0!blue        
  end if 
  enddo
  
  n=16  
  sector = nint(255./real(n))
  do j=0,n-1 !8 sections
    do i=j*sector+1, (j+1)*sector !width of a section
    palette(14, i)     = palette(13,j*sector+1) !red  
    palette(14, i+256) = palette(13,j*sector+1+256)!green
    palette(14, i+512) = palette(13,j*sector+1+512)!blue       
    enddo
  enddo

  

end subroutine FargePalette


end module FieldExport
