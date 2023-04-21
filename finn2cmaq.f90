program finn2cmaq

  use netcdf

  implicit none

  INTEGER, PARAMETER :: ascii = selected_char_KIND ("ascii")
  INTEGER, PARAMETER :: ucs4  = selected_char_KIND ('ISO_10646') 

  type proj_type
     character(16)    :: pName     ! nombre de la proyeccion
     character(7)     :: typ_str   ! String  code for projection TYPE
     integer          :: typ       ! Integer code for projection TYPE
     real             :: ref_lat,ref_lon,truelat1,truelat2,stand_lon,pole_lat,pole_lon
     character(125)   :: proj4     ! PROJ4 srs definition.
  end type proj_type

  type grid_type
      character(12)    :: gName        !grid-name
      integer          :: nx,ny,nz     !number of cells in x-y direction (ncols, nrows, nlevs)
      real             :: dx,dy        !x-y cell dimension (x_cell, y_cell)
      real             :: xmin,ymin,xmax,ymax,xc,yc
      real             :: lonmin,latmin,lonmax,latmax
   end type grid_type

  type(proj_type) :: proj
  type(grid_type) :: grid

  integer :: status,iostat
  integer :: ncid,tstep_dim_id,date_time_dim_id,col_dim_id,row_dim_id,lay_dim_id,var_dim_id,pollut_var_id
  logical :: file_exists
  
  !character(256) :: command
  character(256) :: griddesc_file,finn_data_directory,finnFile,outFile
  character(16)  :: chemistry, pollut
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
  
  namelist /control/ chemistry,start_date,end_date,finn_data_directory,griddesc_file,diurnal_cycle
  
  call date_and_time(values=todays_date)       !fecha de hoy.

  !Leo namelist:
  !open(7, file='parameters'); read(7,parameters); close(7) !leo namelist
  read(*,nml=control, iostat=iostat)
  if( iostat /= 0 ) then
    write(*,*) 'finn2cmaq: failed to read namelist; error = ',iostat
    stop
  end if

  !Leo GRIDDESC:
  call read_GRIDDESC(griddesc_file, proj, grid)    !(!) TO-DO: mejorar esta funcion basado en lo que haga IOAPI

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
    !          if ( atoi(YYYY) < todays_date(1)-2 ) then
    !                    command="wget https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/"//YYYY//"/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt.gz -P finn_data/"
    !                   call system(command)
    !           else
    !                   command="wget https://www.acom.ucar.edu/acresp/MODELING/finn_emis_txt/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt.gz -P finn_data/"
    !                   call system(command)
    !           end if
    !           call system ("gzip -d finn_data/GLOB_"//trim(chemistry)//"_"//YYYY//DDD//".txt.gz")
    !end if

    if ( file_exists ) then
       
       print*," Reading Finn file: ",trim(finnFile)
       !Abro finnFile:
       open(1,file=finnFile, status='old', action='read')
          read(1,'(A)') header                                                  !Aca asumo que FinnFile header es siempre:
          nvars = COUNT((/ (header(i:i) == ',', i=1,len(header)) /)) - 6 + 1    !DAY,TIME,GENVEG,LATI,LONGI,AREA,CO2,CO,...,PM25
          allocate(var_list(nvars))
          allocate(var_units(nvars))
          allocate(emis(nvars))
          read(header,*) colnames,var_list
          
          write(var_list_string,*) var_list
   
          allocate(data(grid%nx,grid%ny,1,24,nvars))        !Asigno memoria a grilla con emisiones
          allocate(tflag(2,nvars,24))                       !Asigno memoria a tflag
          data=0.0
   
          iostat=0
          do while(iostat == 0)  !loop por cada fila de finnFile:
              read(1,*,iostat=iostat) day,time,genveg,lati,longi,area,emis
           
              if ( lati > grid%latmax .or. lati < grid%latmin .or. longi > grid%lonmax .or. longi < grid%lonmin ) then
                continue
              else    
                call gdalTransform(longi,lati,xi,yi,'epsg:4326',proj%proj4) !transformo lati y longi a proyectada xi, yi
                
                ii=floor((xi-grid%xmin)/(2*grid%xmax)*grid%nx) !calculo posición-X en la grilla
                ij=floor((yi-grid%ymin)/(2*grid%ymax)*grid%ny) !calculo posición-Y en la grilla
                
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
           ! Defino dimensiones
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
           call check(nf90_put_att(ncid, nf90_global,"GDTYP"    , 1                    ))
           call check(nf90_put_att(ncid, nf90_global,"P_ALP"    , -50.                 ))
           call check(nf90_put_att(ncid, nf90_global,"P_BET"    , -20.                 ))
           call check(nf90_put_att(ncid, nf90_global,"P_GAM"    , -65.                 ))
           call check(nf90_put_att(ncid, nf90_global,"XCENT"    , proj%ref_lon         ))
           call check(nf90_put_att(ncid, nf90_global,"YCENT"    , proj%ref_lat         ))
           call check(nf90_put_att(ncid, nf90_global,"XORIG"    , grid%xmin            ))
           call check(nf90_put_att(ncid, nf90_global,"YORIG"    , grid%ymin            ))
           call check(nf90_put_att(ncid, nf90_global,"XCELL"    , grid%dx              ))
           call check(nf90_put_att(ncid, nf90_global,"YCELL"    , grid%dy              ))
           call check(nf90_put_att(ncid, nf90_global,"VGTYP"    , -9999                ))
           call check(nf90_put_att(ncid, nf90_global,"VGTOP"    , 5000.                ))
           call check(nf90_put_att(ncid, nf90_global,"VGLVLS"   , [1., 0.9938147 ]     ))
           call check(nf90_put_att(ncid, nf90_global,"GDNAM"    , grid%gName           ))
           call check(nf90_put_att(ncid, nf90_global,"UPNAM"    , "OUTCM3IO"           ))
           call check(nf90_put_att_any(ncid, nf90_global,"VAR-LIST",nf90_char, nvars*16, adjustl(var_list_string))) 
           call check(nf90_put_att(ncid, nf90_global,"FILEDESC" , "Fire emission file" ))
           call check(nf90_put_att(ncid, nf90_global,"HISTORY"  , ""                   ))
       
           !Defino variables
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

print*, "================================="
print*, "finn2cmaq: Completed successfully"
print*, "================================="

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

 subroutine read_GRIDDESC(griddescFile, p, g)                                             
        implicit none                                                                          
        character(256),intent(in) :: griddescFile                                                 
        type(proj_type) ,intent(inout) :: p                                                      
        type(grid_type) ,intent(inout) :: g
        character(10) :: dummyvar
        open(2,file=griddescFile, status='old', action='read')                                  !GRIDDESC:
           read(2,*) dummyvar;                                                   !' '
           read(2,*) p%pName;                                                    !projName
           read(2,*) p%typ,p%truelat1,p%truelat2,p%stand_lon,p%ref_lon,p%ref_lat !map_proj truelat1 truelat2 stand_lon ref_lon ref_lat
           read(2,*) dummyvar;                                                   !' '
           read(2,*) g%gName;                                                    !gridName
           read(2,*) p%pName,g%xmin,g%ymin,g%dx,g%dy,g%nx,g%ny                   !projName xorig yorig xcell ycell nrows ncols
        close(2)
        
        !Calcular otros parametros:
        if (p%typ == 1 ) then           !Geographic:
                p%typ_str='ll';   p%proj4="+proj=latlong +a=6370000.0 +b=6370000.0"  
        else if ( p%typ == 2 ) then     !Lambert Conformal Conic:
                p%typ_str='lcc';  p%proj4="+proj=lcc +lat_1="//trim(rtoa(p%truelat1))//" +lat_2="//trim(rtoa(p%truelat2))//" +lon_0="//trim(rtoa(p%stand_lon))//" +lat_0="//trim(rtoa(p%ref_lat))//" +a=6370000.0 +b=6370000.0 +units=m"
        else if ( p%typ == 3 ) then     !General Mercator
                p%typ_str="merc"; p%proj4="+proj=merc +lat_ts=${truelat1} +a=6370000.0 +b=6370000.0"
        else if ( p%typ == 4 ) then     !General tangent Stereografic
                p%typ_str='stere';p%proj4="+proj=stere +lat_0=${ref_lat} +lon_0=${stand_lon} +lat_ts=lat_ts +a=6370000.0 +b=6370000.0 +k_0=1.0"
        else if ( p%typ == 5 ) then     !UTM
                p%typ_str='utm';  p%proj4="+proj=utm +zone="  
        else if ( p%typ == 6 ) then     !Polar Secant Stereographic
                p%typ_str='stere';p%proj4="+proj=stere +lat_0=${ref_lat} +lon_0=${stand_lon} +lat_ts=lat_ts +a=6370000.0 +b=6370000.0 +k_0=1.0"
        else if ( p%typ == 7 ) then     !Equatorial Mercator
                p%typ_str="merc"; p%proj4="+proj=merc +lat_ts=${truelat1} +a=6370000.0 +b=6370000.0"
        else if ( p%typ == 8 ) then     !Transverse Mercator
                p%typ_str='ll';   p%proj4="+proj=latlong +a=6370000.0 +b=6370000.0"  
        else if ( p%typ == 9 ) then     !Lambert Azimuthal Equal-Area
                print*, "proyección: 9 (Lambert Azimuthal Equal-Area) no soportada por esta aplicación."; stop
        else
                print*, "codigo de proyección invalido.", p%typ; stop
        end if
        
        !Obtener coordenadas del centro de la grilla, min y max:
        g%xc=0.0;g%yc=0.0; g%xmax=(g%xmin)*(-1); g%ymax=(g%ymin)*(-1)

        !transformo boundaries a latlon
        call gdalTransform(g%xmin,g%ymin,g%lonmin,g%latmin,p%proj4,'epsg:4326')
        call gdalTransform(g%xmax,g%ymax,g%lonmax,g%latmax,p%proj4,'epsg:4326')

 end subroutine

 subroutine gdalTransform(x1,y1,x2,y2,srs1,srs2)
        implicit none
        real, intent(in)   :: x1,y1
        real, intent(inout):: x2,y2
        character(*)   :: srs1,srs2
        character(10)  :: ellipsoidh
        character(256) :: command
        command="echo "//rtoa(x1)//" "//rtoa(y1)//" | gdaltransform -s_srs '"//trim(srs1)//"' -t_srs '"//trim(srs2)//"'  > tmp_gdal.txt";
        call system(trim(command))
        !print*,trim(command)
        open(9, file='tmp_gdal.txt', status='old',action='read'); read(9,*, iostat=status) x2, y2, ellipsoidh;  close(9)
        call system('rm tmp_gdal.txt')
 end subroutine

end program finn2cmaq
