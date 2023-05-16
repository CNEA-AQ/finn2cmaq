program finn2cmaq

  use netcdf

  implicit none

  INTEGER, PARAMETER :: ascii = selected_char_KIND ("ascii")
  INTEGER, PARAMETER :: ucs4  = selected_char_KIND ('ISO_10646') 

  real, parameter :: R_EARTH = 6370000.
  real, parameter :: PI = 3.141592653589793
  real, parameter :: RAD2DEG = 180./PI
  real, parameter :: DEG2RAD = PI/180.

  type proj_type
     character(16)    :: pName     ! Projection name
     integer          :: typ       ! Integer code for projection TYPE (2=lcc, 6=stere, 7=merc)
     real             :: alp,bet,gam,xcent,ycent !proj parameters.
     real             :: p1,p2,p3,p4 !extra parameters to speed up calculation once p%typ is defined.
  end type proj_type

  type grid_type
      character(12)   :: gName        !grid-name
      integer         :: nx,ny,nz     !number of cells in x-y direction (ncols, nrows, nlevs)
      real            :: dx,dy        !x-y cell dimension (x_cell, y_cell)
      real            :: xmin,ymin,xmax,ymax,xc,yc
      real            :: lonmin,latmin,lonmax,latmax
   end type grid_type

  type(proj_type) :: proj
  type(grid_type) :: grid

  integer :: status,iostat
  integer :: ncid,tstep_dim_id,date_time_dim_id,col_dim_id,row_dim_id,lay_dim_id,var_dim_id,pollut_var_id
  logical :: file_exists
  
  !character(256) :: command
  character(256) :: griddesc_file,finn_data_directory,finnFile,outFile
  character(16)  :: gridname, chemistry, pollut
  character(16), allocatable :: var_list(:), var_units(:) !lista de polluts
  character(800) :: var_list_string
  character(80) :: att_var_desc
  integer :: nvars
  
  !finnFile:
  character(len=512) :: header!, line
  character(10)      :: colnames(6)           !columnames de finnFile
  integer            :: day,time,genveg       !vars de finnFile
  real               :: lati,longi,area       !vars de finnFile
  real, allocatable  :: emis(:)               !emision de cada fila de finnFile.
 
  integer :: i,j,k,h,ii,ij                    !indices.
  real    :: xi,yi
  real   , allocatable  :: data(:,:,:,:,:)    !buffer donde meter la grilla con valores de emision [nt,nz,nx,ny,nvars]
  integer, allocatable  :: tflag(:,:,:)       !buffer donde meter los valores de TFLAG [nt,nvars,2]

  real, dimension(24) :: diurnal_cycle

  character(len=17) :: start_date, end_date
  integer :: end_date_s, current_date_s 
  character(4) :: YYYY
  character(3) :: DDD
  character(2) :: MM,HH!,DD 
  integer      :: todays_date(8)
  
  namelist /control/ chemistry,start_date,end_date,finn_data_directory,griddesc_file,gridname,diurnal_cycle
  
  call date_and_time(values=todays_date)       !fecha de hoy.

  !Leo namelist:
  read(*,nml=control, iostat=iostat)
  if( iostat /= 0 ) then
    write(*,*) 'finn2cmaq: failed to read namelist; error = ',iostat
    stop
  end if

  !Leo GRIDDESC:
  call read_GRIDDESC(griddesc_file,gridname, proj, grid)    !(!) TO-DO: mejorar esta funcion basado en lo que haga IOAPI

  !Loop over each day
   current_date_s = atoi( date(start_date, "%s") )
       end_date_s = atoi( date(  end_date, "%s") )

  do while (current_date_s <= end_date_s)
 
    YYYY=date("@"//itoa(current_date_s), "%Y")   !año
      MM=date("@"//itoa(current_date_s), "%m")   !mes
     DDD=date("@"//itoa(current_date_s), "%j")   !dia juliano

    write(*,'(A,A4,A,A3,A,A2,A)') "Day: ",YYYY,"-",DDD," (month: ", MM,")"
    
    finnFile=trim(finn_data_directory)//"GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt"
    inquire(file=finnFile, exist=file_exists)
    
    !Descargo finn file:
    !if ( .not. file_exists) then
    !  if ( atoi(YYYY) < todays_date(1)-2 ) then
    !     command="wget https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/"//YYYY//"/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt.gz -P finn_data/"
    !     call system(command)
    !  else
    !     command="wget https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt.gz -P finn_data/"
    !     call system(command)
    !  end if
    !  call system ("gzip -d finn_data/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt.gz")
    !  inquire(file=finnFile, exist=file_exists)
    !end if

    if ( file_exists ) then
       
       print*," Reading Finn file: ",trim(finnFile)
       !Abro finnFile:
       open(1,file=finnFile, status='old', action='read')
          read(1,'(A)') header                                                  !Aca asumo que FinnFile header es siempre:
          nvars = COUNT((/ (header(i:i) == ',', i=1,len(header)) /)) - 6 + 1    !DAY,TIME,GENVEG,LATI,LONGI,AREA,CO2,CO,...,PM25

          allocate(var_list(nvars))                         !array con nombre de polluts
          allocate(emis(nvars))                             !array con emisiones de los polluts
          allocate(var_units(nvars))                        !array con unidades de emisiones
          allocate(data(grid%nx,grid%ny,1,24,nvars))        !grilla de emisiones
          allocate(tflag(2,nvars,24))                       !variable tflag
          data=0.0

          read(header,*) colnames,var_list
          write(var_list_string,*) var_list                !este es un global attr importante.

          iostat=0
          do while(iostat == 0)  !loop por cada fila de finnFile:
              read(1,*,iostat=iostat) day,time,genveg,lati,longi,area,emis
           
              if ( lati > grid%latmax .or. lati < grid%latmin .or. longi > grid%lonmax .or. longi < grid%lonmin ) then
                continue
              else    
                call ll2xy(proj,longi,lati,xi,yi)  !transformo lati y longi a proyectada xi, yi
                !print*,"longi,lati:", longi,lati;print*,"xi,y:", xi,yi
                ii=floor((xi-grid%xmin)/(grid%xmax-grid%xmin)*grid%nx) !calculo posición-X en la grilla
                ij=floor((yi-grid%ymin)/(grid%ymax-grid%ymin)*grid%ny) !calculo posición-Y en la grilla
                
                do k=1,nvars
                        pollut=var_list(k)
                        if ( trim(pollut) == "OC" .or. trim(pollut)  == "BC" .or.  trim(pollut) == "PM25" .or. trim(pollut) == "PM10" ) then
                                emis(k) = emis(k) / 3600000.0      !  kg/day -> g/s    (estrictamente es cierto cuando aplico el ciclo diurno)
                                var_units(k)="g/s"             
                        else
                                emis(k) = emis(k) / 3600.0         !mole/day -> mole/s (estrictamente es cierto cuando aplico el ciclo diurno)
                                var_units(k) = "mole/s"           
                        endif

                        do h=1,24
                            write(HH, '(I0.2)') h-1
                                       
                            tflag(1,:,h) = atoi(YYYY//DDD) 
                            tflag(2,:,h) = atoi(HH//"0000")
                        
                            data(ii,ij,1,h,k) = data(ii,ij,1,h,k) + emis(k)*diurnal_cycle(h)
                         enddo
                enddo
              endif
          enddo
       close(1) !Cierro finnFile
       
       outFile="./emis_fires_"//YYYY//DDD//"_d01.nc"

       print*," Creating NetCDF file: ",trim(outFile)
       ! Create the NetCDF file                                                                                          
       call check(nf90_create(outFile, NF90_CLOBBER, ncid))
           !! Defino dimensiones
           call check(nf90_def_dim(ncid, "TSTEP"    ,   24   , tstep_dim_id    )) 
           call check(nf90_def_dim(ncid, "DATE_TIME",   2    , date_time_dim_id)) 
           call check(nf90_def_dim(ncid, "COL"      , grid%nx, col_dim_id      )) 
           call check(nf90_def_dim(ncid, "ROW"      , grid%ny, row_dim_id      )) 
           call check(nf90_def_dim(ncid, "LAY"      ,   1    , lay_dim_id      )) 
           call check(nf90_def_dim(ncid, "VAR"      , nvars  , var_dim_id      )) 
           !! Defino attributos
           call check(nf90_put_att(ncid, nf90_global,"IOAPI_VERSION", "ioapi-3.2: \$Id: init3" ))
           call check(nf90_put_att(ncid, nf90_global,"EXEC_ID"  , "????????????????"   ))
           call check(nf90_put_att(ncid, nf90_global,"FTYPE"    , 1                    ))
           call check(nf90_put_att(ncid, nf90_global,"SDATE"    , atoi(YYYY//DDD)      ))   !int
           call check(nf90_put_att(ncid, nf90_global,"STIME"    , 000000               ))
           call check(nf90_put_att(ncid, nf90_global,"WDATE"    , 2023001              ))
           call check(nf90_put_att(ncid, nf90_global,"WTIME"    , 000000               ))
           call check(nf90_put_att(ncid, nf90_global,"CDATE"    , 2023001              ))
           call check(nf90_put_att(ncid, nf90_global,"CTIME"    , 000000               ))
           call check(nf90_put_att(ncid, nf90_global,"TSTEP"    , 10000                ))
           call check(nf90_put_att(ncid, nf90_global,"NTHIK"    , 1                    ))   
           call check(nf90_put_att(ncid, nf90_global,"NCOLS"    , grid%nx              ))
           call check(nf90_put_att(ncid, nf90_global,"NROWS"    , grid%ny              ))
           call check(nf90_put_att(ncid, nf90_global,"NLAYS"    , 1                    ))!grid%nz              
           call check(nf90_put_att(ncid, nf90_global,"NVARS"    , nvars                ))
           call check(nf90_put_att(ncid, nf90_global,"GDTYP"    , proj%typ             ))
           call check(nf90_put_att(ncid, nf90_global,"P_ALP"    , proj%alp             ))                    
           call check(nf90_put_att(ncid, nf90_global,"P_BET"    , proj%bet             ))                    
           call check(nf90_put_att(ncid, nf90_global,"P_GAM"    , proj%gam             ))                    
           call check(nf90_put_att(ncid, nf90_global,"XCENT"    , proj%xcent           ))
           call check(nf90_put_att(ncid, nf90_global,"YCENT"    , proj%ycent           ))
           call check(nf90_put_att(ncid, nf90_global,"XORIG"    , grid%xmin            ))
           call check(nf90_put_att(ncid, nf90_global,"YORIG"    , grid%ymin            ))
           call check(nf90_put_att(ncid, nf90_global,"XCELL"    , grid%dx              ))
           call check(nf90_put_att(ncid, nf90_global,"YCELL"    , grid%dy              ))
           call check(nf90_put_att(ncid, nf90_global,"VGTYP"    , -9999                ))!averiguar que poner
           call check(nf90_put_att(ncid, nf90_global,"VGTOP"    , 5000.                ))!averiguar que poner
           call check(nf90_put_att(ncid, nf90_global,"VGLVLS"   , [1., 0.9938147 ]     ))!averiguar que poner
           call check(nf90_put_att(ncid, nf90_global,"GDNAM"    , grid%gName           ))
           call check(nf90_put_att(ncid, nf90_global,"UPNAM"    , "OUTCM3IO"           ))!averiguar que poner
           call check(nf90_put_att_any(ncid, nf90_global,"VAR-LIST",nf90_char, nvars*16, adjustl(var_list_string))) 
           call check(nf90_put_att(ncid, nf90_global,"FILEDESC" , "Fire emission file" ))
           call check(nf90_put_att(ncid, nf90_global,"HISTORY"  , ""                   ))
           !!Defino variables
           call check(nf90_def_var(ncid,"TFLAG"       ,NF90_FLOAT    , [date_time_dim_id,var_dim_id,tstep_dim_id], pollut_var_id))
           call check(nf90_put_att(ncid, pollut_var_id, "units"      , "<YYYYDDD,HHMMSS>" ))
           call check(nf90_put_att(ncid, pollut_var_id, "long_name"  , "TFLAG           " ))
           call check(nf90_put_att(ncid, pollut_var_id, "var_desc"   , "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                "))
           do k=1, nvars             
             pollut=var_list(k) 
             att_var_desc=trim(pollut)//"[1]"
             call check(nf90_def_var(ncid, pollut       ,NF90_FLOAT , [col_dim_id,row_dim_id,lay_dim_id,tstep_dim_id], pollut_var_id)) !
             call check(nf90_put_att(ncid, pollut_var_id,"units"    , trim(var_units(k)) ))
             call check(nf90_put_att(ncid, pollut_var_id,"long_name", pollut             ))
             call check(nf90_put_att(ncid, pollut_var_id,"var_desc" , att_var_desc       ))
           end do

       call check(nf90_enddef(ncid))
       !End NetCDF define mode   
        
       !Abro NetCDF outFile
       call check(nf90_open(outFile, nf90_write, ncid))
       do k=1,nvars
           pollut=var_list(k)
           call check(nf90_inq_varid(ncid, trim(pollut), pollut_var_id))     !Obtengo id de variable
           call check(nf90_put_var(ncid, pollut_var_id, data(:,:,:,:,k)  ))  !Escribo valores en NetCDF
       enddo
       call check(nf90_inq_varid(ncid, "TFLAG"    , pollut_var_id))                                
       call check(nf90_put_var(ncid, pollut_var_id, tflag(:,:,:) ))

       call check(nf90_close(ncid))
       !Cierro NetCDF outFile

    deallocate(data)       !Libero memoria
    deallocate(tflag)      !Libero memoria
    deallocate(var_list)   !Libero memoria
    deallocate(var_units)  !Libero memoria
    deallocate(emis)       !Libero memoria
   
    endif
    
    current_date_s=current_date_s + 86400  !siguiente día!
  end do

print*, "==================================="
print*, " finn2cmaq: Completed successfully "
print*, "==================================="
contains

 subroutine check(status)
   integer, intent(in) :: status
   if (status /= nf90_noerr) then
     write(*,*) nf90_strerror(status)
     stop 'netcdf error'
   end if
 end subroutine check

 !Interfaz a "date"
 function date(date_str, fmt_str) result(output)
   implicit none
   character(*), intent(in) :: date_str, fmt_str
   character(256)           :: command
   character(20)            :: output
   command="date -d "//trim(date_str)//" '+"//trim(fmt_str)//"'  > tmp_date.txt"
   call system( trim(command) )
   !print*,trim(command)
   open(9, file='tmp_date.txt', status='old',action='read'); read(9, '(A)', iostat=status) output;  close(9)
   call system('rm tmp_date.txt')
 end function

 function atoi(str)     !string -> int
   implicit none
   character(len=*), intent(in) :: str
   integer :: atoi
   read(str,*) atoi
 end function
 function itoa(i)       !int -> string
    implicit none
    integer, intent(in) :: i
    character(len=20) :: itoa
    write(itoa, '(i0)') i
    itoa = adjustl(itoa)
 end function
 function rtoa(r)       !real -> string
    implicit none
    real, intent(in) :: r
    character(len=16) :: rtoa
    write(rtoa, '(F16.3)') r
    rtoa = adjustl(rtoa)
 end function

 subroutine read_GRIDDESC(griddescFile,gridName, p, g)
    implicit none
    character(200),intent(in) :: griddescFile
    character(*) ,intent(in)  :: gridName
    type(proj_type), intent(inout) :: p
    type(grid_type), intent(inout) :: g
    character(20) :: row
    iostat=0
    open(unit=2,file=griddescFile,status='old',action='read',access='sequential')
    do while(iostat == 0)  !loop por cada fila
       read(2,*,iostat=iostat) row
       if ( trim(row) == trim(gridname)) then
         g%gName=row
         read(2,*) p%pName,g%xmin,g%ymin,g%dx,g%dy,g%nx,g%ny !projName xorig yorig xcell ycell nrows ncols
         rewind(2)
       endif
       if (trim(row) == trim(p%pName)) then
         read(2,*) p%typ,p%alp,p%bet,p%gam,p%xcent,p%ycent   !map_proj truelat1 truelat2 stand_lon ref_lon ref_lat
         iostat=1
       endif
    enddo
    close(2)

    !Calculate proj parameters used then for coordinate transformations:
    call set_additional_proj_params(p)

    !Obtener coordenadas del centro de la grilla, min y max:
    g%xc=0.0;g%yc=0.0;g%xmax=g%xmin+g%dx*g%nx; g%ymax=g%ymin+g%dy*g%ny

    !transformo boundaries a latlon
    call xy2ll(p,g%xmin,g%ymin,g%lonmin,g%latmin)
    call xy2ll(p,g%xmax,g%ymax,g%lonmax,g%latmax)
    
    !!debug:-------------------------------
    !print*,"Test: xy -> ll ..."
    !print*,"g%xmin, g%xmax, g%ymin, g%ymax      ",g%xmin, g%xmax, g%ymin, g%ymax
    !print*,"g%lonmin,g%lonmax,g%latmin,g%latmax ",g%lonmin,g%lonmax,g%latmin,g%latmax
    !print*,"Test: ll -> xy ..."
    !call ll2xy(p,g%lonmin,g%latmin,g%xmin,g%ymin)
    !call ll2xy(p,g%lonmax,g%latmax,g%xmax,g%ymax)
    !print*,"g%lonmin,g%lonmax,g%latmin,g%latmax ",g%lonmin,g%lonmax,g%latmin,g%latmax
    !print*,"g%xmin, g%xmax, g%ymin, g%ymax      ",g%xmin, g%xmax, g%ymin, g%ymax
    !!------------------------------------- 
 end subroutine

 !COORDINATE TRANSFORMATION FUNCTIONS:======================================
 subroutine xy2ll(p,x,y,lon,lat)
     implicit none                            
     type(proj_type) ,intent(in) :: p
     real, intent(in)   :: x,y
     real, intent(inout):: lon,lat
 
     if      ( p%typ == 2 ) then  !Lambert Conformal Conic:
        call xy2ll_lcc(p,x,y,lon,lat)
     else if ( p%typ == 6 ) then  !polar secant stereographic
        call xy2ll_stere(p,x,y,lat,lon)
     else if ( p%typ == 7 ) then  !equatorial mercator
        call xy2ll_merc(p,x,y,lon,lat)
     else
        print*, "codigo de proyección invalido:",p%typ,"."; stop
     end if
 end subroutine
 subroutine ll2xy(p,lon,lat,x,y)
       implicit none                            
       type(proj_type) ,intent(in) :: p
       real, intent(in):: lon,lat
       real, intent(inout)   :: x,y
 
       if      ( p%typ == 2 ) then  !Lambert Conformal Conic:
          call ll2xy_lcc(p,lon,lat,x,y)
       else if ( p%typ == 6 ) then  !Polar Secant Stereographic
          call ll2xy_stere(p,lon,lat,x,y)
       else if ( p%typ == 7 ) then  !Equatorial Mercator
          call ll2xy_merc(p,lon,lat,x,y)
       else
          print*, "codigo de proyección invalido:",p%typ,"."; stop
       end if                             
 end subroutine

 subroutine set_additional_proj_params(p)
    implicit none                            
    type(proj_type) ,intent(inout) :: p
 
    if ( p%typ == 2 ) then       !lambert conformal conic:        
       if ( ABS(p%alp - p%bet) > 0.1 ) then  !secant proj case
          p%p2=     LOG( COS(p%alp           *deg2rad )/ COS(p%bet           *deg2rad)   )                                    
          p%p2=p%p2/LOG( TAN((45.0+0.5*p%bet)*deg2rad )/ TAN((45.0+0.5*p%alp)*deg2rad)   ) !n
        else                                 !tangent proj case
          p%p2=SIN(p%alp*deg2rad) !n
       endif
       p%p3=R_EARTH*(COS(p%alp*deg2rad)*TAN((45+0.5*p%alp)*deg2rad)**p%p2)*(1/p%p2)  !F
       p%p1=p%p3/(TAN((45 + 0.5*p%ycent)*deg2rad)**p%p2)                             !rho0 
                                                                                                                               
    else if ( p%typ == 6 ) then  !polar secant stereographic
       print*, "Todavia no desarrollado soporte para proyeccion polar stereografica (ptype=",p%typ,")."; stop

    else if ( p%typ == 7 ) then  !equatorial mercator
       p%p1=COS(p%alp*deg2rad)   !k0

    else
        print*, "codigo de proyección invalido:",p%typ,"."; stop
    end if
 end subroutine
 !--------------------------------------------------------------------------
 !LAMBERT CONFORMAL CONIC:
 subroutine xy2ll_lcc(p,x,y,lon,lat)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)   :: x,y
   real, intent(inout):: lon,lat
   real :: n,F,rho0,rho,theta
   
   rho0=p%p1
   n=p%p2
   F=p%p3
   
   theta=ATAN(x/(rho0-y))*rad2deg
   rho=SIGN(1.0,n) * SQRT( x*x + (rho0-y)*(rho0-y))
   
   lon=p%gam+theta/n
   lat=2.0 * ATAN( (F/rho)**(1/n) )*rad2deg - 90.0 
 end subroutine
 subroutine ll2xy_lcc(p,lon,lat,x,y)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)      :: lon,lat
   real, intent(inout)   :: x,y
   real :: n,F,rho0,rho,dlon

   !interm params:
   rho0=p%p1
   n=p%p2
   F=p%p3

   rho=F/(TAN((45.0 + 0.5*lat)*deg2rad)**n)
   dlon=lon-p%gam
   !
   x=     rho*SIN(n*dlon*deg2rad )
   y=rho0-rho*COS(n*dlon*deg2rad )
 end subroutine
 !--------------------------------------------------------------------------
 !MERCATOR                
 subroutine xy2ll_merc(p,x,y,lon,lat)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)   :: x,y
   real, intent(inout):: lon,lat
   real :: k0R,phi
 
   k0R=p%p1*R_EARTH     !es una longitud (k0*R_EARTH)
   phi=y/k0R            !es un angulo
   
   lon=p%gam + x/k0R*rad2deg
   lat=90.0-2*ATAN( EXP(-phi) )*rad2deg
 end subroutine
 subroutine ll2xy_merc(p,lon,lat,x,y)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)      :: lon,lat
   real, intent(inout)   :: x,y
   real :: k0,lon0

   k0=p%p1              !adminesional
   lon0=p%gam           !es un angulo

   x=k0*R_EARTH*(lon-lon0)*deg2rad
   y=k0*R_EARTH*LOG(TAN((45.0+0.5*lat)*deg2rad))
 end subroutine
!--------------------------------------------------------------------------
 !POLAR STEREOGRAPHIC     
 subroutine xy2ll_stere(p,x,y,lon,lat)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)   :: x,y
   real, intent(inout):: lon,lat
   real :: k,rho

   rho = sqrt(x*x+y*y)
   k = 2.0*ATAN( rho/2.0/R_EARTH )

   lat =         ASIN(   COS(k)*SIN(p%alp*deg2rad) + y*SIN(k)*COS(p%alp*deg2rad)/rho )               * rad2deg
   lon = p%gam + ATAN( x*SIN(k)  / ( rho*COS(p%alp*deg2rad)*COS(k) - y*SIN(p%alp*deg2rad)*SIN(k) ) ) * rad2deg

 end subroutine
 subroutine ll2xy_stere(p,lon,lat,x,y)
   implicit none                            
   type(proj_type) ,intent(in) :: p
   real, intent(in)      :: lon,lat
   real, intent(inout)   :: x,y
   real :: k!,hemi

   !hemi=SIGN(1.0,p%alp)
   k = 2.0*R_EARTH / (1 + SIN(p%alp*deg2rad)*SIN(lat*deg2rad) + COS(p%alp*deg2rad)*COS(lat*deg2rad)*COS( (lon-p%gam)*deg2rad ))

   x = k *   COS( lat *deg2rad) * SIN( (lon - p%gam)*deg2rad )
   y = k * ( COS(p%alp*deg2rad) * SIN(  lat         *deg2rad ) - SIN(p%alp*deg2rad)*COS(lat*deg2rad)*COS((lon-p%gam)*deg2rad) )
 end subroutine
!!END COORDINATE TRANFORMATION FUNCTIONS====================================

end program finn2cmaq
