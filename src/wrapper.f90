!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program is part of the T-RECS code
! Author: A. Bonaldi
! see Bonaldi et al. (2018) for more details
! reads the results from sampler and collate them 
! it should be run separately for SFGs and AGNs, as they have different catalogue columns
! projects the Euclidian coordinates for the objects 
! on the sphere and rotates to specified centre of FoV
! generates spherical random coordinates for big fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



program wrapper

  use healpix_types
  USE fitstools
  USE utilities
  USE pix_tools
!  USE random_tools,only: ran_mwc,randgauss_boxmuller
  use sampler_io, ONLY: rows_catalogue, read_catalogue, columns_catalogue
  USE paramfile_io, ONLY : paramfile_handle, parse_init, parse_int, &
       parse_string, parse_double, parse_lgt, concatnl
  USE extension, ONLY : getEnvironment, getArgument, nArguments
  implicit none

  real, parameter::sterad_deg=3283. 
  character(LEN=filenamelen)::paramfile,description,dummy,slice_file,arg,outdir,indir,tag
  character(LEN=filenamelen)::filename,outcat,ctmp
  character(LEN=filenamelen),allocatable::catfiles(:),zstr(:)
  character*16,allocatable::tagnames(:),tform(:),tunit(:)
  character*16::tagname
  character*16::lat_tag,lon_tag,x_tag,y_tag,x_tag_1,y_tag_1
  real(sp),allocatable::data(:,:),column(:)
  real(sp)::sim_area
  real(dp)::v(3),v0(3),v_diff(3),v_new(3)
  real(dp)::lat_astro,lon_astro,lat_temp,lon_temp,lat,lon,lat_target,lon_target
  real(dp)::thetam,phim,vnorm,dx,dy,dz,sim_side,sim_radius,pix_side_rad,npixels,npixtot
  CHARACTER(LEN=5) :: output,output2
  CHARACTER(LEN=10) ::output3
  CHARACTER(LEN=80), DIMENSION(1:120) :: header
  CHARACTER(LEN=80)::line
  integer(I4B)::Nargs,l_tag,l_indir,l_outdir,l_str
  integer(I4B),allocatable::l_catname(:)
  integer(I4B) :: nest,test(1)
  integer(I8B):: Nrows
  integer::Nfiles,i,Ncols,iostat,l,ii,l_tagname,i_lat,i_lon,i_x,i_y,i_x1,i_y1
  integer status,unit,rows,nrows_i,readwrite,blocksize,hdutype,tfields,IERR
  integer::bitpix,naxis,naxes ,jj,do_flatuniform
  integer(I4B), dimension(:),  allocatable :: listpix,listpix_copy(:)
  integer(I4B)::nside,nlist,ipix,ipix2
  INTEGER, DIMENSION(8,2) :: values_time
  REAL(SP) :: clock_time
  !REAL(SP) :: lat_target_rad,lon_target_rad,lat_rad,lon_rad
  integer(4) :: ic4, crate4, cmax4,ni,iseed
  CHARACTER(LEN=*), PARAMETER :: code ="Wrapper"
  TYPE(paramfile_handle) :: handle
  character extname*16
  integer varidat, colnum,frow,felem

  !double precision ispline
  save iseed

  !getting seed from the clock
  call system_clock(count=ic4, count_rate=crate4, count_max=cmax4)
  iseed=ic4

  ! setting the random seed
  call srand(iseed)

  call date_and_time(values = values_time(:,1))

  !1)input file:
  !  -bins in redshift
  !  -number of objects
  
  ! parse command line arguments
  Nargs = nArguments()
  if (mod(Nargs,2) /= 0) then
     print '(2a, /)', 'probable wrong command-line options'
     call print_help()
     stop 1
  endif
  i = 1
  tag = 'continuum'
  paramfile = ''
  do while ( i <= Nargs )

     call getArgument(i, arg)
     select case (arg)
     case ('-p', '--params')
        i = i+1
        call getArgument(i, paramfile)

     case ('-t', '--tag')
        i = i+1
        call getArgument(i, tag)

     case ('-h', '--help')
        call print_help()
        stop

     case default
        print '(2a, /)', 'unrecognised command-line option: ', arg
        call print_help()
        stop 1
     end select
     i = i+1

  end do

  handle = parse_init(paramfile)
  
  description = concatnl( &
       & " Enter the tag of input catalogues")
  tag=parse_string(handle,'tag', default=tag, descr=description)

  description = concatnl( &
       & " Enter the output directory path")
  outdir=parse_string(handle,'outdir', default='.', descr=description)

  description = concatnl( &
       & " Enter the input directory base path")
  indir=parse_string(handle,'wrap_indir', default=outdir, descr=description)

  description = concatnl( &
       & " Enter the central latitude coordinate ")
  lat_target = parse_double(handle, 'lat_target', default=0.d0, vmin=-90.d0, vmax=90.d0,descr=description)

  description = concatnl( &
       & " Enter the central longitude coordinate ")
  lon_target = parse_double(handle, 'lon_target', default=0.d0, vmin=0.d0, vmax=360.d0,descr=description)

  description = concatnl( &
       & " Enter the simulation size: ")
  sim_side = parse_double(handle, 'sim_side', default=5.d0, vmin=0.d0, descr=description)
  ! to compute the sky area. for the random coordinates, the FoV is taken as a circle

  description = concatnl( &
       & " Do you want to simulate coordinates with flat approximation (0=no, 1=yes)?")
  do_flatuniform = parse_int(handle, 'do_flatuniform', default=0, vmax=1, descr=description)

  description = concatnl( &
       & " Enter the slices file listing all the redshifts to wrap")
  slice_file=parse_string(handle,'slice_file', default='none', descr=description)

  !string formatting: eliminate spaces between path and file name
  tag=ADJUSTL(tag)
  l_tag=LEN_TRIM(tag)
  
  indir=ADJUSTL(indir)
  l_indir=LEN_TRIM(indir)

  outdir=ADJUSTL(outdir)
  l_outdir=LEN_TRIM(outdir)

  

  ! cont the number of files to wrap from the summary slice file
  if (slice_file=='none') then
     slice_file = indir(1:l_indir)//'/slices_'//tag(1:l_tag)//'.dat'
  endif
  open( UNIT=42, file=slice_file, status='old', form='formatted' )
  Nfiles = 0
  do while (ierr == 0)
     Nfiles = Nfiles + 1
     read(42, *, iostat=IERR) ctmp
  end do
  Nfiles = Nfiles - 1
  print*, 'Number of files = ', Nfiles
  allocate(catfiles(Nfiles))
  allocate(l_catname(Nfiles))
  allocate(zstr(Nfiles))

  if (iostat /=0) then
     print*,'Error allocating filenames'
     stop
  endif

  !here read zstr and write the input file names
  rewind ( 42 )
  do i=1,NFiles
     read( 42, '(a)') zstr(i)
     zstr(i) = adjustl(zstr(i))
     l_str = len_trim(zstr(i))
     catfiles(i)=indir(1:l_indir)//'/raw_'//tag(1:l_tag)//'/catalogue_'//tag(1:l_tag)//'_z'//zstr(i)(1:l_str)//'.fits'
     catfiles(i)=adjustl(catfiles(i))
     l_catname(i) = len_trim(catfiles(i))
  enddo
  close( 42 )
  deallocate(zstr)

  outcat=outdir(1:l_outdir)//'/catalogue_'//tag(1:l_tag)//'_wrapped.fits'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call columns_catalogue(catfiles(1)(1:l_catname(1)),Ncols,status)
  
  if (status ==0) then 
     print*,'number of columns in file',Ncols
  endif

  allocate(tagnames(Ncols),tunit(Ncols),tform(Ncols),stat=iostat)
  tunit(:)=''
  tform(:)='1E'

  ! generate position vector for the chosen central coordinate
  thetam=lat_target*pi/180.d0 !in radians
  thetam=-thetam+pi/2.d0  ! in healpix convention
  phim=lon_target*pi/180.d0  !in radians and healpix convention

  ! coo vector pointing target coordinates
  v0(1)=dsin(thetam)*dcos(phim)
  v0(2)=dsin(thetam)*dsin(phim)
  v0(3)=dcos(thetam)
  vnorm=dsqrt(v0(1)**2.d0+v0(2)**2.d0+v0(3)**2.d0)
  v0(1)=v0(1)/vnorm
  v0(2)=v0(2)/vnorm
  v0(3)=v0(3)/vnorm

  ! preparing to read the catalogues
  !count total number of rows
  Nrows=0

  do i=1,Nfiles
     call rows_catalogue(catfiles(i)(1:l_catname(i)),Ncols,nrows_i,tagnames,tunit)
     nrows=nrows+nrows_i
  enddo
  print*, 'total number of rows',Nrows



  if (do_flatuniform==0) then
     !initialization necessary for generating random coordinates
     !getting the list of healpix pixels in the area

     sim_area=sim_side*sim_side/sterad_deg  ! simulation area steradians

     !select appropriate Nside depending on number of galaxies to place

   

     do i=8,12
        nside=2.**i
        npixtot=12.*(2.**i)**2.
        npixels=npixtot/4./pi*sim_area  ! number of pixels in the sky area for t

        if (npixels > Nrows) goto 100 ! choosing so that roughly 1 galaxy per pixel

       enddo
100     continue

       print*,'nside=',nside
        sim_radius=sqrt(sim_area/pi)
        pix_side_rad=(0.1145/512*nside)*pi/180. ! mean spacing between pixels

        ! allocate pixel list (bigger than necessary make sure query_disc does not exceed it)
        allocate(listpix_copy(nint(npixels*1.1)-1))
        call query_disc(nside, v0, sim_radius, listpix_copy, nlist)
        !resize list of pixels to the actual number found by query_disc
        allocate(listpix(nlist))
        listpix=listpix_copy(1:nlist)
        deallocate(listpix_copy)

        print*,'simulated area [sterad]=',sim_area
        print*,'selected area [sterad]=',nlist*4.*pi/npixtot

     endif

     ! identify columns with the lat and lon coordinates from the catalogues
     lat_tag='latitude'
     lon_tag='longitude'
     x_tag='x_coord'
     y_tag='y_coord'
     x_tag_1='x_coord_1'
     y_tag_1='y_coord_1'
     i_x1=-1
     do ii=1,Ncols
        tagname=tagnames(ii)
        if (tagname ==lat_tag) i_lat=ii
        if (tagname ==lon_tag) i_lon=ii
        if (tagname ==x_tag) i_x=ii
        if (tagname ==y_tag) i_y=ii
        if (tagname ==x_tag_1) i_x1=ii
        if (tagname ==y_tag_1) i_y1=ii
     enddo
     if (i_x1 /=-1) print*,'Second set of coordinates found'



     tform(i_lat)='1D'
     tform(i_lon)='1D'
     !open first catalogue and read header
     status=0
     filename=outcat

     !opening new file
     unit = 66
     call ftgiou(unit,status)
     if (status .gt. 0) call printerror(status)
     !C     open the FITS file, with write access
     readwrite=1

     blocksize=1
     call ftinit(unit,filename,blocksize,status)
     if (status .gt. 0) call printerror(status)

     ! write primary header
     bitpix=8
     naxis=0
     naxes=0
     call FTPHPS(unit,bitpix,naxis,naxes,status)

     !C     move to the last (2nd) HDU in the file
     call ftmahd(unit,1,hdutype,status)
     if (status .gt. 0) call printerror(status)

     !C     append/create a new empty HDU onto the end of the file and move to it
     call ftcrhd(unit,status)
     if (status .gt. 0) call printerror(status)

     tfields=Ncols
     extname='Catalogue'
     varidat=0

     !write table header
     call ftphbn(unit,nrows,tfields,tagnames,tform,tunit,extname,varidat,status)
     if (status .gt. 0) call printerror(status)



     select case(do_flatuniform)

     case(0)
        ! generate random coordinates for a disc of specified area 

        ! write rows one by one
        frow=1
        felem=1
        do ii=1,nfiles
           print*,'copying catalogue',catfiles(ii)(1:l_catname(ii))
           call rows_catalogue(catfiles(ii)(1:l_catname(ii)),Ncols,nrows_i,tagnames,tunit)
           allocate(data(nrows_i,Ncols),column(nrows_i))
           call read_catalogue(catfiles(ii)(1:l_catname(ii)),Ncols,Nrows_i,data)
           print*,'number of rows',nrows_i
           ! select a random pixel of the list
           do jj=1,nrows_i

              ipix=int(rand()*(nlist-1))+1  !extraction of one random pixel position
              call pix2ang_ring(nside,listpix(ipix),lat_astro,lon_astro)
111           continue

              !add a random shift to coordinates to avoid pixellization effects. Uniform distribution
              lat_temp=lat_astro+(rand()-0.5)*pix_side_rad
              lon_temp=lon_astro+(rand()-0.5)*pix_side_rad/cos(lat_temp)
 
              ! ensure new coordinates are within the healpix range 
              if (lat_temp <0.) lat_temp=lat_temp+Pi
              if (lat_temp >pi) lat_temp=lat_temp-Pi

              if (lon_temp <0.) lon_temp=lon_temp+2.*Pi
              if (lon_temp >2*pi) lon_temp=lon_temp-2.*Pi

              !ensure the coordinates are still within the same pixel. If not, new random shift
              call ang2pix_ring(nside,lat_temp,lon_temp,ipix2)
              test=minval(abs(ipix2-listpix))
              if (test(1) /=0) then 
                 goto 111
              endif


              !conversion to degrees and astronomical convention on latitude range
              lon_astro=lon_temp/PI*180.d0
              lat_astro=(-lat_temp+pi/2.d0)/pi*180.d0
              data(jj,i_lat)=real(lat_astro)
              data(jj,i_lon)=real(lon_astro)
           enddo
           print*,'writing'
           print*,'file ',ii,' of',nfiles
           ! write the content on the new file
           felem=1
           i=1
           do while (i <= Ncols)
              column=data(:,i)
              call ftpcle(unit,i,frow,felem,nrows_i,column,status)  
              if (status .gt. 0) call printerror(status)
              i=i+1
           end do
           i=i-2

           frow=frow+nrows_i

           deallocate(data,column)
        enddo
     case(1)
        ! Project the Euclidian coordinates on the sphere with given central coordinates
        ! this works well for small FoVs and is the option to use for simulation with clustering

        ! write rows one by one
        frow=1
        felem=1
        do ii=1,nfiles
           print*,'copying catalogue',catfiles(ii)(1:l_catname(ii))

           call rows_catalogue(catfiles(ii)(1:l_catname(ii)),Ncols,nrows_i,tagnames,tunit)

           allocate(data(nrows_i,Ncols),column(nrows_i))

           call read_catalogue(catfiles(ii)(1:l_catname(ii)),Ncols,Nrows_i,data)

           ! project coordinates on sphere and rotate to the chosen position
           do jj=1,nrows_i
              if (data(jj,i_x)==-100.) then ! this is a crossmatched catalogue, Use the recond set of coordinates
                 data(jj,i_x)=data(jj,i_x1) 
                 data(jj,i_y)=data(jj,i_y1) 
              endif
              lat=dble(data(jj,i_y))
              lon=dble(data(jj,i_x))  
              lat_astro=lat+lat_target
              thetam=lat_astro*pi/180.d0
              lon_astro=lon/dcos(thetam)+lon_target
              data(jj,i_lat)=real(lat_astro)
              data(jj,i_lon)=real(lon_astro)
           enddo


           ! write the content on the new file
           felem=1

           i=1
           do while (i <= Ncols)

              column=data(:,i)
              call ftpcle(unit,i,frow,felem,nrows_i,column,status)  
              if (status .gt. 0) call printerror(status)
              i=i+1
           end do
           i=i-2

           frow=frow+nrows_i

           deallocate(data,column)

        enddo

     end select

     !C     close the FITS file and free the unit number
     call ftclos(unit, status)
     call ftfiou(unit, status)

     ! deallocate the raw catalogue names array
     deallocate(catfiles)
     deallocate(l_catname)

     !C     check for any error, and if so print out error messages
     if (status .gt. 0)call printerror(status)


     call date_and_time(values = values_time(:,2))

     values_time(:,1) = values_time(:,2) - values_time(:,1)
     clock_time =  (  (values_time(3,1)*24 &
          &           + values_time(5,1))*60. &
          &           + values_time(6,1))*60. &
          &           + values_time(7,1) &
          &           + values_time(8,1)/1000.
     PRINT*,"Total clock time [m]:",clock_time/60.
     PRINT*,"               "//code//" > Normal completion"

   contains
     subroutine print_help()
       print '(a, /)', 'command-line options:'
       print '(a)',    '  -p, --params /path/to/parameter_file.ini path to the input parameter file'
       print '(a)',    '  -t, --tag tagname                        tag name'
       print '(a, /)', '  -h, --help        print usage information and exit'
     end subroutine print_help
   end program wrapper
