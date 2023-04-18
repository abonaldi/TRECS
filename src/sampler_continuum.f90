
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program is part of the T-RECS code
! Author: A. Bonaldi
! see Bonaldi & Hartley (2023) MNRAS for more details
! generate samples of radio galaxies from models
! save outputs to files per redshift bin 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program sampler

  use healpix_types
  USE fitstools
  USE utilities
  USE pix_tools
  USE random_tools
  use random,only:random_poisson,random_normal
  use sampler_io
  use trecs_mods
  use interpolations
  use cosmology
  USE paramfile_io, ONLY : paramfile_handle, parse_init, parse_int, &
       parse_string, parse_double, parse_lgt, concatnl, parse_real
  USE extension, ONLY : getEnvironment, getArgument, nArguments
  implicit none

  !parameters
  real(sp),parameter::zmax=8.
  CHARACTER(LEN=*), PARAMETER :: code ="T-RECS"
  real(dp),parameter:: fullsky_area=41253._dp
  REAL, PARAMETER :: pc=3.0856776E18 !Pc in cm
  real, parameter::sterad_deg=3283. 
  REAL, PARAMETER :: Mpc=10.**6*pc
  integer,parameter::Nmbins=10
  integer,parameter::mass_scatter=1
  real(sp),parameter::flagvalue=-100.

  !character variables
  character(LEN=filenamelen)::paramfile,description,dummy,outdir,indir,rawdir
  real(dp),allocatable::frequencies(:),frequencies_rest(:),Lsyn(:),Lfree(:),Ld(:)
  real(sp),allocatable::spec(:)
  real(dp)::nu,deltanu,sfr,mn,mx,volume,integ,volumetot,fom,fom_old
  character(LEN=filenamelen)::SFR_filename,SFR2Mh_filename,filename,freq_filename,SFR2Mstar_filename
  character(LEN=filenamelen)::filename1,filename2,filename3,filename4,filename5,filename6
  character(LEN=filenamelen)::filename7,filename8,filename9,filename10
  character(LEN=filenamelen)::AGN_filename,LERG_filename,HERG_filename
  character(LEN=filenamelen)::chline,filestat,cat_filename
  character(LEN=16),allocatable::tagnames(:),tunit(:),tform(:)
  character(LEN=10)::names(3)
  CHARACTER(LEN=5) :: output,output2,tag
  CHARACTER(LEN=10) ::output3
  CHARACTER(LEN=80), DIMENSION(1:120) :: header
  CHARACTER(LEN=7)::line,col_str,zstr_long
  CHARACTER(LEN=4),allocatable::redshift_names(:)
  CHARACTER(LEN=7),allocatable::redshift_names_lum(:)
  CHARACTER(LEN=4) ::zstr

  !single precision variables
  real(sp)::z_min,z_max,draw
  real(sp)::squarealphas(3)
  real(sp),allocatable::delta(:),work(:,:),diff(:)
  real(sp),allocatable::Radioflux(:,:),catout(:,:),catout_copy(:,:),samplex(:),lums(:),z_gals(:)
  real(sp),allocatable::Polaflux(:,:),inclinations(:),qrat(:),ellipticity1(:)
  real(sp),allocatable::ellipticity2(:),polafracs(:),lums_save(:,:),thetas(:),freq_rest(:),lums_i(:)!
  real(sp)::sample100(100),rn,l_iii,nu_iii
  real(sp),allocatable::Radioflux_copy(:,:),samplex_copy(:),sizes(:),radioflux_slice(:,:),lums_slice(:),mstar_slice(:)
  real(sp),allocatable::angles(:),sizes_3d(:),lums_copy(:),spot_dist(:),himass(:),mstar_copy(:)
  real(sp),allocatable::darkmass(:),latitudes(:),longitudes(:),mstar(:)
  real(sp),allocatable::samplex_run(:),redshifts(:),L14(:),SFRtab(:),mstartab(:)
  real(sp),allocatable::masses_L14(:),redshifts_lum(:),samplex_slice(:),ids(:)
  real(sp),allocatable::masses_lerg(:),masses_herg(:),lerg_p(:),herg_p(:)
  REAL(SP) :: clock_time,coo_max,coo_min,dlat,dlon,dz,delta_parall
  REAL(SP) ::deltaz,mu,flux14,dmtol,mgal,ncone_binned,m_star
  real(sp)::dm_model,dm_model2,dm_best,lat,lon,maxflux(1),redshift,minmass,minmass_cone,q
  real(sp)::dim,lum_threshold,theta_high,theta_low,sin_i,cos_i,logs,alpha_low
  real(sp)::alpha_high,minfreq,norm_f,dlogl,col,b,a,arg1,arg2,th,freq_norm
  real(sp)::minm(Nmbins),maxm(Nmbins),mc(1),halfmax(1),mw(1),sg_m(1)
  real(sp)::a_sfr,b_sfr
  real(sp),parameter::sigma_alpha=2.35660,mean_alpha=0.609847 

  !double precision variables
  real(dp)::sim_area,skyfrac
  real(dp)::sim_side
  real(dp)::fluxlim,fluxlim_freq,current_count
  real(dp),allocatable::data(:,:),x(:),px(:),dustsed(:,:),dif(:)
  real(dp),allocatable::data2(:,:),x2(:),px2(:),x3(:),px3(:)
  real(dp),allocatable::alphascaling_14_48(:,:),alphascaling_48_20(:,:),alphascaling_015_14(:,:)
  real(dp),allocatable::poladistr(:,:),alphascaling_14_48_pola(:,:),alphascaling_48_20_pola(:,:)
  real(dp)::alphamed(3),z,alpha,conv_flux,zlow,zhigh,nfullsky,it,Ngen_db,z_i
  real(dp)::dx,xmin,xmax,val,test,thr,norm,over,dx2,lumlim,deltax
  real(dp)::zlow_lum,zhigh_lum,zlum,par_2,sigma_lin,L_stop_2

  !integer variables
  integer::no_AGN,no_SFG,save_lums,Nsample_binned,istart(1),iend(1),istart_i,iend_i,iii,jstart
  integer*8::Nsample,Nsample_lerg,Nsample_herg,join,k,n_hist,nside,ntiles,t,Nsample100
  integer::iseed,Nfreq,Nfreq_out,i,ii,nrows,nrows_sfrtab,Ncolumns,Ncolumns_lf,nskip,Ngen,Ngen2,jj,kk,nrows_sf,seed_fix
  integer::buffer_size,buffer_free,jbuf,buffer_free_old,buffer_size_old,l_outdir,l_indir,l_rawdir
  integer::count,Ncat_common,nrows_lerg,nrows_herg,seed(34),nrows_old,nrows_lf
  integer::Nsample_old,Nfunction,Nsample_surv,Nhaloes,Nsample_mass,nreds_out
  integer::l,j,nreds,nreds2,nreds_cone,zi,nloop,ist,iostat,nreds_lum,islice
  integer(4) :: ic4, crate4, cmax4,ni,sum_plus,p14(1),i14,ilim,i48,p(1),p_i,try
  INTEGER, DIMENSION(8,2) :: values_time
  integer::N,N3,Nran,reason,iunit
  integer,allocatable::indx(:,:),indx2(:,:),poisson_numbers(:),optclass(:)
  logical::first=.true.
  integer::system_status
  !types variables
  TYPE(paramfile_handle) :: handle


  save iseed

  call date_and_time(values = values_time(:,1))
  call system_clock(count=ic4, count_rate=crate4, count_max=cmax4)
  iseed=ic4
  test=rand(time())


  ! check calling sequence
  if (nArguments() == 0) then
     paramfile=''
  else if (nArguments() == 1) then

     call getArgument(1,paramfile)
     filestat='old'
  else
     print '("Usage: sampler: [parameter file name]")'
     stop 1
  endif


  ! --- declaration of intent
  PRINT*," "
  PRINT*,"               "//code
  PRINT*," "

  !*************************************
  !read input parameters from file or interactively
  !*************************************

  handle = parse_init(paramfile)

  description = concatnl( &
       & " Enter the name of the file with frequency list: ")
  freq_filename=  parse_string(handle, 'freqs_file', default=chline, descr=description)

  description = concatnl( &
       & " Do you want to skip the AGN simulation (0=no, 1=yes)")
  no_AGN = parse_int(handle, 'no_AGN', default=0, vmax=1, descr=description)

  description = concatnl( &
       & " Do you want to skip the SFG simulation (0=no, 1=yes)")
  no_SFG = parse_int(handle, 'no_SFG', default=0, vmax=1, descr=description)

  description = concatnl( &
       & " Enter the simulation size: ")
  sim_side = parse_double(handle, 'sim_side', default=5.d0, vmin=0.d0, descr=description)

  description = concatnl( &
       & " Enter the flux limit [Jy]: ")
  fluxlim = parse_double(handle, 'fluxlim_cont', default=10.d-9, descr=description)

  description = concatnl( &
       & " Enter the frequency at which the fluxlimit is imposed [MHz]: ")
  fluxlim_freq = parse_double(handle, 'fluxlim_cont_freq', default=1400.d0, descr=description)

  description = concatnl( &
       & " Do you want to output the luminosities (0=no, 1=yes)? ")
  save_lums = parse_int(handle, 'save_lums', default=0, vmax=1, descr=description)

  description = concatnl( &
       & " Enter the minimum redhift to consider: ")
  z_min = parse_real(handle, 'z_min', default=0., vmin=0., descr=description)

  description = concatnl( &
       & " Enter the maximum redshift to consider: ")
  z_max = parse_real(handle, 'z_max', default=8., vmin=0., descr=description)

  description = concatnl( &
       & " Enter the name of the the output file directory:")
  outdir=  parse_string(handle, 'outdir', default='.', descr=description)

  description = concatnl( &
       & " Enter the name of the the input files directory:")
  indir=  parse_string(handle, 'TRECSinputs', default='.', descr=description)


  description = concatnl( &
       & " Set the seed for random number generation (default=automatic)" )
  seed_fix = parse_int(handle, 'seed', default=-1, descr=description)

  ! end reading input parameters

  !string formatting: eliminate spaces between path and file name
  indir=ADJUSTL(indir)
  l_indir=LEN_TRIM(indir)

  !string formatting: eliminate spaces between path and file name
  outdir=ADJUSTL(outdir)
  l_outdir=LEN_TRIM(outdir)

  ! rename output directory and create it if does not exist
  rawdir=outdir(1:l_outdir)//'/raw_continuum'
  rawdir=ADJUSTL(rawdir)
  l_rawdir=LEN_TRIM(rawdir)
  system_status = SYSTEM( 'mkdir -p '//rawdir(1:l_rawdir) )
  if ( system_status /= 0 ) then
     stop system_status
  endif

  !***************************************************
  !Global filenames for the input data to be read:
  !
  !    AGNs
  !
  !effective spectral indices for AGNs
  filename1=indir(1:l_indir)//'/TRECS_Inputs/alphaeff/alphascaling_1.4_4.8.txt'
  filename2=indir(1:l_indir)//'/TRECS_Inputs/alphaeff/alphascaling_4.8_20.txt'
  filename3=indir(1:l_indir)//'/TRECS_Inputs/alphaeff/alphascaling_1.4_0.150.txt'



  !polarization fractions, Galluzzi et l. and Hales et al. 
  filename7=indir(1:l_indir)//'/TRECS_Inputs/AGN_polafraction/Polfrac_AGNs_1e4_hales.dat' 

  !characteristic luminosity: redshift bins
  filename8=indir(1:l_indir)//'/TRECS_Inputs/LF/CharLum/Characteristic_Luminosity/z_bins.dat'

  !Intrinsic AGN size distribution from DiPompeo et al.
  filename9=indir(1:l_indir)//'/TRECS_Inputs/LF/AGN_sizes_new.dat'  

  !   SFGs
  !
  !dust SED
  filename10=indir(1:l_indir)//'/TRECS_Inputs/LF/SEDdust.dat'!nu, UVgal, spheroids, lens_spheroids

  ! end filenames
  !**********************************************************

  !Seeding random number generators
  !Starts from a user-specified initial seed or from the execution clock
  
  if (seed_fix ==-1) then
     call system_clock(count=ic4, count_rate=crate4, count_max=cmax4)
     do i=1,34
        seed(i)=ic4*i/12.+iseed  ! combining two clock readings to fill 12 seeds
     enddo
  else
     seed_fix = 555 * seed_fix
     do i=1,34
        seed(i)=1284350*i+seed_fix  ! combining two clock readings to fill 12 seeds
     enddo
     iseed=seed_fix
  endif


  call random_seed(PUT=seed)  ! feeding to Poisson and Gaussian generator 
  rn=rand(iseed+807820) ! initialise random uniform generator

  !simulation area 
  sim_area=sim_side*sim_side
  skyfrac=sim_area/fullsky_area

  print*,'simulation area [square degrees]=',sim_area


  !read input files with effective spectral indices for the AGN frequency behaviour
  !1.4-4.8 GHz range

  Ncolumns=5
  n=rows_number(filename1,1,nskip)
  allocate(data(n,ncolumns),alphascaling_14_48(n,5)) !flux, flux, alphaflat,alphaflat,alphasteep
  call read_columns(filename1,2,n,Ncolumns,nskip,data)


  do i=1,n
     alphascaling_14_48(i,1)=data(i,1)           !flux flat
     alphascaling_14_48(i,2)=data(i,2)-data(i,3)           !
     alphascaling_14_48(i,3)=data(i,2)-data(i,3)           !alpha_flat
     alphascaling_14_48(i,4)=data(i,4)-data(i,5)           !alpha_steep
  enddo


  deallocate(data)
  !4.8-20 GHz range
  Ncolumns=5
  n=rows_number(filename2,1,nskip)
  allocate(data(n,ncolumns),alphascaling_48_20(n,5)) !flux, flux, alphaflat,alphaflat,alphasteep
  call read_columns(filename2,2,n,Ncolumns,nskip,data)

  do i=1,n
     alphascaling_48_20(i,1)=data(i,1)           !flux flat
     alphascaling_48_20(i,2)=data(i,2)-data(i,3)           !
     alphascaling_48_20(i,3)=data(i,2)-data(i,3)           !alpha_flat
     alphascaling_48_20(i,4)=data(i,4)-data(i,5)           !alpha_steep
  enddo

  deallocate(data)

  !150MHz-1.4 GHz range
  Ncolumns=5
  n=rows_number(filename3,1,nskip)
  allocate(data(n,ncolumns),alphascaling_015_14(n,5)) !flux, flux, alphaflat,alphaflat,alphasteep
  call read_columns(filename3,2,n,Ncolumns,nskip,data)


  do i=1,n
     alphascaling_015_14(i,1)=data(i,1)           !flux flat
     alphascaling_015_14(i,2)=data(i,2)+data(i,3) !different sign because 1.4->150           !
     alphascaling_015_14(i,3)=data(i,2)+data(i,3)           !alpha_flat
     alphascaling_015_14(i,4)=data(i,4)+data(i,5)           !alpha_steep
  enddo

  deallocate(data)


  ! reading file with polarization fraction information
  ! using Hales et al. eq 5 for steep spectrum sources and Galluzzi for flat

  Ncolumns=4
  n=rows_number(filename7,1,nskip)
  N3=n

  allocate(poladistr(n,ncolumns),x3(n),px3(n)) 
  !polarization fraction, N_FSRQ,N_BLLac, N_steep
  call read_columns(filename7,2,n,Ncolumns,nskip,poladistr)

  ! finished reading file with polarization fraction information


  !read input file with list of frequencies:
  Ncolumns=1
  nfreq=rows_number(freq_filename,1,nskip)

  if (nfreq==0) then !stopping if file emply or not provided
     print*,'ERROR: no frequencies specified!!'
     stop
  endif


  allocate(data(nfreq,ncolumns))
  call read_columns(freq_filename,2,nfreq,Ncolumns,nskip,data)
  Nfreq_out=Nfreq
  ! add to the frequency vector the frequencies: 
  ! 1400 and 4800 MHz and the frequency at which the flux cut is done, whether they are already 
  !included or not, because they are needed for the modelling. those added frequencies
  ! will not be included in the outputs, but just used internally by the code.

  Nfreq=Nfreq+3 ! add 1400 MHz, 4800 MHz and selection frequency
  allocate(frequencies(Nfreq),spec(Nfreq))
  frequencies(1)=1400.
  frequencies(2)=4800.
  frequencies(3)=fluxlim_freq
  frequencies(4:Nfreq)=data(:,1)

  deallocate(data)

  !position of the reference frequencies in the frequency vector
  i14=1
  i48=2
  ilim=3

  allocate(frequencies_rest(Nfreq),Lsyn(Nfreq),Lfree(Nfreq),delta(Nfreq),Ld(Nfreq))


  ! debug message on the size of the field
  coo_max=sim_side/2.
  coo_min=-1.*coo_max
  print*,'coordinates',coo_max,coo_min



  !structure of the catalogue
  Ncat_common=(2+save_lums)*Nfreq_out+20 !number of catalogue columns
  !  commmon format for SFGs and AGNs

  !creating tag names for the catalogue
  allocate(tagnames(Ncat_common),tunit(Ncat_common),tform(Ncat_common))
  j=1
  ! Unique source ID
  tagnames(j)='ID_cont'
  tunit(j)='none'
  j=2
  ! Luminosity at 1.4 GHz
  tagnames(j)='Lum1400'
  tunit(j)='log(erg/s/Hz)'
  j=3
  ! SFR 
  tagnames(j)='logSFR'
  tunit(j)='log(Msun/yr)'
  j=4
  ! total intensity flux density for all selected frequencies 
  do i=1,Nfreq_out
     write(output,"(i5)")int(frequencies(i+3))
     output=ADJUSTL(output)
     l=LEN_TRIM(output)
     dummy = 'I'//output(:l)
     tagnames(j)=dummy
     tunit(j)='mJy'
     j=j+1
  enddo

  !polarised intensity for all selected frequencies
  do i=1,Nfreq_out
     write(output,"(i5)")int(frequencies(i+3))
     output=ADJUSTL(output)
     l=LEN_TRIM(output)
     dummy = 'P'//output(:l)
     tagnames(j)=dummy
     tunit(j)='mJy'
     j=j+1
  enddo

  ! dark mass
  tagnames(j)='Mh'
  tunit(j)='log(Msun)'
  j=j+1
  ! stellar mass
  tagnames(j)='Mstar'
  tunit(j)='log(Msun)'
  j=j+1
  ! HI mass proxy
  tagnames(j)='MHI_pred'
  tunit(j)='log(Msun)'
  j=j+1
  ! coordinate shift with respect to field center in the x direction
  tagnames(j)='x_coord'
  tunit(j)='degs'
  j=j+1
    ! coordinate shift with respect to field center in the y direction
  tagnames(j)='y_coord'
  tunit(j)='degs'
  j=j+1
  ! latitude
  tagnames(j)='latitude'
  tunit(j)='degs'
  j=j+1
  !longitude
  tagnames(j)='longitude'
  tunit(j)='degs'
  j=j+1
  !redshift
  tagnames(j)='redshift'
  tunit(j)='none'
  j=j+1
  ! apparent size 
  tagnames(j)='size'
  tunit(j)='arcsec'
  j=j+1
  ! inclination 
  tagnames(j)='inclination'
  tunit(j)='degs'
  j=j+1
  ! axis ratio
  tagnames(j)='axis ratio'
  tunit(j)='none'
  j=j+1
  ! major axis
  tagnames(j)='bmaj'
  tunit(j)='arcsec'
  j=j+1
  ! monor axis
  tagnames(j)='bmin'
  tunit(j)='arcsec'
  j=j+1
  ! position angle
  tagnames(j)='PA'
  tunit(j)='degs'
  j=j+1
  ! Rs FRI/FRII classification parameter
  tagnames(j)='Rs'
  tunit(j)='none'
  j=j+1
  ! radio classification
  tagnames(j)='RadioClass'
  tunit(j)='none'
  j=j+1
  ! optical classification
  tagnames(j)='OptClass'
  tunit(j)='none'
  j=j+1

  if (save_lums ==1) then
     !optional tag names for the luminosities
     !Luminosity at all considered frequencies
     do i=1,Nfreq_out
        write(output,"(i5)")int(frequencies(i+3))
        output=ADJUSTL(output)
        l=LEN_TRIM(output)
        dummy = 'L'//output(:l)
        tagnames(j)=dummy
        tunit(j)='log(erg/s/Hz)'
        j=j+1
     enddo
  endif
  tform(:)='1E'


  !*************************
  !AGN simulation starts
  !*************************


  current_count=0. !number of sources generated for unique ID

  ! reading bin information for the AGN characteristic luminosity model
  nrows=rows_number(filename8,1,nskip)

  close(iunit)
  open (UNIT=iunit,file=filename8,status='unknown',form="formatted")

  do i=1,nskip
     READ(iunit,*,IOSTAT=Reason)  col
  enddo

  nreds_lum=nrows
  allocate(redshifts_lum(nreds_lum),redshift_names_lum(nreds_lum))


  DO i=1,nrows
     READ(iunit,*,IOSTAT=Reason)  col_str,col
     redshift_names_lum(i)=col_str
     redshifts_lum(i)=col
  ENDDO
  close(iunit)

  !reading the first characteristic luminosity file to get dimension of files for later
  AGN_filename=indir(1:l_indir)//'/TRECS_Inputs/LF/CharLum/Characteristic_Luminosity/luminosity_0.01000.dat' 
  Ncolumns=6
  Ncolumns_lf=Ncolumns
  nrows=rows_number(AGN_filename,1,nskip)
  Nrows_lf=nrows
  n_hist=nrows

  !endif

  !read dust SEDs for SFGs
  Ncolumns=4
  nrows=rows_number(filename10,1,nskip)
  allocate(dustsed(nrows,ncolumns),dif(nrows))
  call read_columns(filename10,2,nrows,Ncolumns,nskip,dustsed)


  dustsed(:,2)=10.**dustsed(:,2)
  dustsed(:,3)=10.**dustsed(:,3)
  dustsed(:,4)=10.**dustsed(:,4)

  ! definition of redshoft intervals
  nreds=57 ! gals until z=8 (where I have cone)
  nreds_out=nreds

  allocate(redshifts(nreds),redshift_names(nreds))

  ! redshift slices 
  redshifts=(/0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,&
       0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00,1.20,1.40,1.60,1.80,2.00,2.20,2.40,&
       2.60,2.80,3.00,3.20,3.40,3.60,3.80,4.00,4.20,4.40,4.60,4.80,5.00,5.20,5.40,5.60,5.80,&
       6.00,6.20,6.40,6.60,6.80,7.00,7.20,7.40,7.60,7.80,8.00,8.20,8.40,8.60,8.80,9.00,9.20,9.40,9.60,9.80,10.00/)

  redshift_names=(/'0.01','0.02','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50',&
       '0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00','1.20','1.40','1.60','1.80','2.00'&
       ,'2.20','2.40','2.60','2.80','3.00','3.20','3.40','3.60','3.80','4.00','4.20','4.40','4.60','4.80','5.00'&
       ,'5.20','5.40','5.60','5.80','6.00','6.20','6.40','6.60','6.80','7.00','7.20','7.40','7.60','7.80','8.00',&
       '8.20','8.40','8.60','8.80','9.00','9.20','9.40','9.60','9.80','10.0'/)

  if (zmax < maxval(redshifts)) then
     nreds_out=0
     i=1
     do  while (redshifts(i)<=zmax)
        nreds_out=nreds_out+1
        i=i+1
     enddo
     nreds_out=nreds_out+1
  endif

  ! names of AGN populations
  !  names(1)='FSRQ'
  !  names(2)='BLLac'
  !  names(3)='SS_AGN'

  alphamed=(/-0.1,-0.1,-0.72/) ! median spectral index for the three populations. SS spectral index Bonato et al. 2021 

  ! open summary file for wrapper
  open(42, file = outdir(1:l_outdir)//'/slices_continuum.dat', status = 'new')
  do zi=1,nreds_out-1   ! Main redshift loop 

     z=redshifts(zi)

     if ((z_min <= z) .and. (z_max >= z)) then    ! redshift slice with center z is processed

        !getting the size of the redshift slice
        select case (zi)
        case (1)
           zlow=0.
           zhigh=(redshifts(zi+1)+z)/2.
        case default
           zlow=(redshifts(zi-1)+z)/2.
           zhigh=(redshifts(zi+1)+z)/2.
        end select

        print*,'************************'
        print*,'AGNs: Processing redshift',z,' Range',zlow,zhigh
        print*,'************************'

        zstr=redshift_names(zi) ! redshift tag for files

        !files for AGN mass modelling
        LERG_filename=indir(1:l_indir)//'/TRECS_Inputs/AbundanceMatch/results_AGNs/AGNprob_LERG_z'//zstr//'.txt' 
        HERG_filename=indir(1:l_indir)//'/TRECS_Inputs/AbundanceMatch/results_AGNs/AGNprob_HERG_z'//zstr//'.txt' 

        !reading files for AGN mass modelling 
        Ncolumns=2
        nrows=rows_number(LERG_filename,1,nskip)
        nrows_lerg=nrows
        allocate(data(nrows,ncolumns),lerg_p(nrows),masses_lerg(nrows))
        call read_columns(LERG_filename,2,nrows,Ncolumns,nskip,data)
        lerg_p=data(:,2)
        masses_lerg=data(:,1)
        
        deallocate(data)

        Ncolumns=2
        nrows=rows_number(HERG_filename,1,nskip)
        nrows_herg=nrows
        allocate(data(nrows,ncolumns),herg_p(nrows),masses_herg(nrows))
        call read_columns(HERG_filename,2,nrows,Ncolumns,nskip,data)
        herg_p=data(:,2)
        masses_herg=data(:,1)
        deallocate(data)
        !reading done

        Nsample_old=0
        do ii=1,3 !loop on AGN populations  FSRQ, BL-LAC, SteepS
           ! information for limiting ram usage
           ! processing long files in chuncks of lenght buffer_size
           buffer_size=1000
           buffer_free=buffer_size
           jbuf=0 ! index to fill buffer 

           ! here allocate Radioflux to store total intensity results
           allocate(Radioflux(Nfreq,buffer_size),samplex(buffer_size))

           allocate(data(n_hist,ncolumns_lf),x(n_hist),px(n_hist),poisson_numbers(n_hist))  ! allocate arrays to read in luminosity files


           ! view angle for the sources depending on the sub-population
           select case(ii)
           case(3)

              theta_low=cos(90.*pi/180.)
              theta_high=cos(5.*pi/180.)
              
           case default

              theta_high=cos(0.)
              theta_low=cos(5.*pi/180.)
              
           end select

           !redshift bin is split into thin redshift slices
           do islice=1,nreds_lum-1
              !loading files for thin redshift slices within the redshift bin
              zlum=redshifts_lum(islice)
              select case (islice)
              case (1)
                 zlow_lum=0.
                 zhigh_lum=(redshifts_lum(islice+1)+zlum)/2.
              case default
                 zlow_lum=(redshifts_lum(islice-1)+zlum)/2.
                 zhigh_lum=(redshifts_lum(islice+1)+zlum)/2.
              end select

              if ((zlow_lum >= zlow) .and. (zlum <= zhigh)) then ! here zlum and not zhigh_lum because i dont want to lose any small slice 

                 zstr_long=redshift_names_lum(islice) !string tag 

                 conv_flux=1./(4.*pi*(lumr(zlum)*Mpc)**2)*1.e26*(1+zlum)!L [erg/s/Hz] ->S [mJy]

                 frequencies_rest=frequencies*(1.+zlum) !freq in the rest frame (for K correction)

                 !read luminosities 
                 AGN_filename=indir(1:l_indir)//'/TRECS_Inputs/LF/CharLum/Characteristic_Luminosity/luminosity_'//zstr_long//'.dat'  
                 call read_columns(AGN_filename,2,nrows_lf,Ncolumns_lf,nskip,data)

                 x=data(:,2) ! log of characteristic luminosity
                 deltax=x(2)-x(1) ! luminosity bin size    
                 px=10.**data(:,ii+3) !N(logl) per sterad 
                 norm=sim_area/sterad_deg  

                 !poission realization to give the number of objects givin the distribution
                 do i=1,Nrows_lf
                    mu=real(Px(i)*norm)
                    Nran=random_Poisson(mu, first)! integer from poisson distribution with mean mu
                    poisson_numbers(i)=Nran
                 enddo

                 Nsample=sum(poisson_numbers)  ! this is the number of galaxies to generate

                 if (Nsample==0) goto 100

                 allocate(samplex_slice(Nsample),radioflux_slice(Nfreq,Nsample))! arrays where to store objects

                 
                 !assigning luminosity to the objects
                 iii=1
                 do i=1,Nrows_lf
                    Nran=poisson_numbers(i)
                    do j=1,Nran
                       samplex_slice(iii)=x(i)+rand()*dx 
                       iii=iii+1
                    enddo

                 enddo


                 alpha=alphamed(ii) ! starting point for the spectral index 
                 ! compute flux on just one spectral index and count number of objects above flux limit
                 ! SED is refined later on

                 Nsample_surv=0
                 do i=1,Nsample
                    !converting luminosity to flux, scale flux in frequency
                    freq_norm=1400.*(1.+random_normal()/1000.) ! to introduce some scatter in flux14 from characteristic luminosity

                    flux14=(10.**samplex_slice(i)*(frequencies_rest(i14)/freq_norm)**alpha)*conv_flux !no scatter


                    Radioflux_slice(:,i)=flux14*(frequencies/1400.)**alpha
                    spec=radioflux_slice(:,i)

                    !refine frequency scaling with the "effective indics" and add scatter
                    call effective_index(i14,i48,ii,alpha,alphascaling_14_48&
                         &,alphascaling_48_20,alphascaling_015_14,spec,frequencies)

                    radioflux_slice(:,i)=spec


                    !implement flux threshold
                    if (Radioflux_slice(ilim,i)>= fluxlim*1000.) Nsample_surv= Nsample_surv+1

                 enddo


                 if (buffer_free < Nsample_surv) then 
                    ! expand buffer
                    Print*,'resize buffer'
                    buffer_size_old=buffer_size
                    buffer_free_old=buffer_free
                    buffer_size=buffer_size_old+(buffer_size_old+Nsample_surv)*10
                    buffer_free=buffer_free+(buffer_size_old+Nsample_surv)*10

                    allocate(Radioflux_copy(Nfreq,buffer_size_old),samplex_copy(buffer_size_old))
                    radioflux_copy=radioflux
                    samplex_copy=samplex
                    deallocate(radioflux,samplex)
                    allocate(Radioflux(Nfreq,buffer_size),samplex(buffer_size))
                    Radioflux(:,1:buffer_size_old)=Radioflux_copy(:,:)
                    samplex(1:buffer_size_old)=samplex_copy(:)
                    deallocate(radioflux_copy,samplex_copy)
                 endif

                 ! fill buffer
                 do i=1,Nsample
                    if (Radioflux_slice(ilim,i)>= fluxlim*1000.) then
                       jbuf=jbuf+1
                       Radioflux(:,jbuf)=Radioflux_slice(:,i)
                       samplex(jbuf)=samplex_slice(i)
                    endif
                 enddo
                 buffer_free=buffer_size-jbuf

                 deallocate(samplex_slice,radioflux_slice)

100              continue
              endif  ! small slice in big slice

           enddo ! loop on small slices

           deallocate(data,x,px,poisson_numbers) 



           Nsample=buffer_size-buffer_free
           if (Nsample==0) then
              !skip resize and output catalogue if no object is found
              deallocate(radioflux,samplex)
              goto 400 
           endif
           ! resize Radioflux to the final sample size
           if (buffer_free /=0) then
              print*,'resize final flux ' 
              allocate(radioflux_copy(Nfreq,buffer_size),samplex_copy(buffer_size),stat=iostat)
              if (iostat /=0) then
                 print*,'Allocation error'
                 stop
              endif
              radioflux_copy=radioflux
              samplex_copy=samplex
              deallocate(radioflux,samplex,stat=iostat)
              if (iostat /=0) then
                 print*,'Deallocation error'
                 stop
              endif

              Nsample=buffer_size-buffer_free
              allocate(radioflux(Nfreq,Nsample),samplex(nsample),stat=iostat)
              if (iostat /=0) then
                 print*,'Allocation error'
                 stop
              endif
              radioflux(:,:)=radioflux_copy(:,1:Nsample)
              samplex(:)=samplex_copy(1:Nsample)
              deallocate(radioflux_copy,samplex_copy,stat=iostat)
              if (iostat /=0) then
                 print*,'Deallocation error'
                 stop
              endif
           endif
           ! total intensity flux completed, flux selection completed


           !From here everything is done on the big slice

           !begin polarised intensity model
           !polarization fraction PDF 
           x3=poladistr(:,1)
           px3=poladistr(:,ii+1)
           px3=px3/sum(px3)*Nsample  !normalization of pdf

           ! allocate work arrays and array to store the results
           allocate(polaflux(Nfreq,Nsample),polafracs(Nsample),stat=iostat) 
           if (iostat /=0) then
              print*,'Allocation error'
              stop
           endif

           ! poisson-sample from the polarization PDF but on the same Nsample objects already extracted
           ! not very accurate for very small samples: I use a minimum sample of 100 end extract randomly from that sample if Nsample<100 
           if (Nsample >=100) then 
              call poisson_constrained(iseed,x3,Px3,0.,100.,polafracs,Nsample) 
           else
              Nsample100=100
              sample100(:)=0.
              call poisson_constrained(iseed,x3,Px3,0.,100.,sample100,Nsample100)
              call reordersample(iseed,sample100)
              do i=1,Nsample
                 polafracs(i)=sample100(i)
              enddo
           endif

           do i=1,Nfreq
              polaflux(i,:)=polafracs/100.*radioflux(i,:)
           enddo

           alpha=alphamed(ii)
           !frequency scaling in polarization
           do i=1,Nsample
              spec=polaflux(:,i)

              call effective_index(i14,i48,ii,alpha,alphascaling_14_48&
                   &,alphascaling_48_20,alphascaling_015_14,spec,frequencies)
              polaflux(:,i)=spec
           enddo

           deallocate(polafracs,stat=iostat)
           if (iostat /=0) then
              print*,'Dellocation error'
              stop
           endif
           !end polarization model

           !allocate arrays to store the other galaxy attributes
           allocate(Darkmass(Nsample),himass(Nsample),latitudes(Nsample),longitudes(Nsample),&
                &z_gals(Nsample),sizes_3d(Nsample),sizes(Nsample),&
                &angles(Nsample),spot_dist(Nsample),optclass(Nsample),mstar(Nsample),ids(Nsample),stat=iostat)
           if (iostat /=0) then
              print*,'Allocation error'
              stop
           endif

           !todo:change
           himass(:)=0.
           optclass(:)=1 ! Starting point for AGN is to have elliptical optical class. This is refiled ;ater

           ! size and view angle from precomputed distributions
           spot_dist(:)=0. !Distance between two bright spots (FR classification) 
           !initialised to 0 because BLLAC and FSRQ will not have two lobes

           !AGN size model
           ! load file with intrinsic size distribution
           ! from DiPompeo et al. 2013 Table 2 
           ! for flat-spectrum AGN Di Pompeo et al. 2013 Table 2 "narrow"
           ! for steep-spectrum AGN Di Pompeo et al. 2013 Table 2 "wide" 
           Ncolumns=4
           nrows=rows_number(filename9,1,nskip)


           allocate(data2(nrows,ncolumns),x2(nrows),px2(nrows),stat=iostat)
           if (iostat /=0) then
              print*,'sizes allocation failed'
              stop
           endif
           call read_columns(filename9,2,nrows,Ncolumns,nskip,data2)

           x2=data2(:,1) !size Kpc
           px2=data2(:,ii+1) !N(size) (different model for flat and steep-spectrum)
           ! in different columns of the file

           dx2=abs(x2(2)-x2(1))
           N=Nrows



           ! sample from size distribution

           if (Nsample >=100) then 
              call poisson_constrained(iseed,x2,Px2,0.,2010.,sizes_3d,Nsample)
           else
              Nsample100=100
              sample100(:)=0.
              call poisson_constrained(iseed,x2,Px2,0.,2010.,sample100,Nsample100)
              call reordersample(iseed,sample100)
              do i=1,Nsample
                 sizes_3d(i)=sample100(i)
              enddo
           endif

           print*,'sizes minmax',minval(sizes_3d),maxval(sizes_3d)
           !stop
           deallocate(data2,x2,px2,stat=iostat)
           if (iostat /= 0) then
              print*,'size dealloc failed'
              stop
           endif

           print*,'sizes ends'
           ! end AGN size model

           ! Dark mass and clustering model
           print*,'dark mass starts'

           ! generate masses from haloes for the distributions dervied for HERG/LERG
           select case(ii) !associate HERG/LERG to our 3 populations 

           case(1)
              !FSRQ -> HERG
              iii=1
              darkmass(:)=0.
              dx=masses_herg(2)-masses_herg(1)
              norm=dble(Nsample)

              do i=1,Nrows_herg
                 mu=real(herg_p(i)*norm)
                 Nran=random_Poisson(mu, first)! integer from poisson distribution with mean mu
                 xmax=masses_herg(i)+dx/2
                 xmin=masses_herg(i)-dx/2

                 do j=1,Nran
                    Darkmass(iii)=rand()*(xmax-xmin)+xmin
                    iii=iii+1
                    if (iii > Nsample ) then
                       goto 120 
                    endif
                 enddo
              enddo

              ! fill remaining masses (if any) with uniform distributions
              xmax=masses_herg(Nrows_herg)
              xmin=masses_herg(1)

              do iii=1,Nsample
                 if (Darkmass(iii)==0.) Darkmass(iii)=rand()*(xmax-xmin)+xmin
              enddo

120           continue

           case (2)

              !BLLAc ->   !LERG
              iii=1
              darkmass(:)=0.
              dx=masses_lerg(2)-masses_lerg(1)
              norm=dble(Nsample)

              do i=1,Nrows_lerg
                 mu=real(lerg_p(i)*norm)
                 Nran=random_Poisson(mu, first)
                 xmax=masses_lerg(i)+dx/2
                 xmin=masses_lerg(i)-dx/2

                 do j=1,Nran
                    Darkmass(iii)=rand()*(xmax-xmin)+xmin
                    iii=iii+1
                    if (iii > Nsample ) then
                       goto 220 
                    endif
                 enddo

              enddo

              ! fill remaining masses if any - unif distr.

              xmax=masses_lerg(Nrows_lerg)
              xmin=masses_lerg(1)

              do iii=1,Nsample
                 if (Darkmass(iii)==0.) Darkmass(iii)=rand()*(xmax-xmin)+xmin
              enddo

220           continue

           case(3)
              !SS AGN -> LERG, HERG

              ! SS herg/lerg depending on 1.4 GHz luminosity
              !§ for this herg/lerg also generate rs (distance between lobes) according to Lin et al:
              !HERG 0.62 +-0.18 -> FRII
              !LERG: 0.17 +-0.11 ->FRI


              ! count how many lergs and hergs based on luminority threshold
              lum_threshold=24.6+7 
              ! 1.3e26 at 180 MHz (Lin et al.), scaled with -0.8 to 1.4 GHz
              ! between power and luminosity there is also the 4pi sterad factor
              !the +7 is to convert from W to erg/s
              Nsample_herg=0

             
              do i=1,Nsample
                 if (samplex(i) > lum_threshold) Nsample_herg=Nsample_herg+1
              enddo
              Nsample_lerg=Nsample-Nsample_herg



              print*,'herg/lerg samples:',Nsample,Nsample_herg,Nsample_lerg


              if (Nsample_herg > 0) then

                 ! herg 
                 iii=1
                 darkmass(:)=0.
                 dx=masses_herg(2)-masses_herg(1) ! this is same for herg and lerg
                 norm=dble(Nsample_herg)


                 do i=1,Nrows_herg
                    mu=real(herg_p(i)*norm)
                    Nran=random_Poisson(mu, first)! integer from poisson distribution with mean mu
                    xmax=masses_herg(i)+dx/2
                    xmin=masses_herg(i)-dx/2

                    do j=1,Nran

                       ! jump to next element with flux over threshold
                       do while (samplex(iii) <= lum_threshold)
                          iii=iii+1  
                       enddo
                       if (iii > Nsample )                       goto 320 
                       ! these sources are morphologically described as FRII
                       ! generate distance between the bright spots
                       spot_dist(iii)=random_normal()*0.18+0.62 !FRII
                       Darkmass(iii)=rand()*(xmax-xmin)+xmin
                       iii=iii+1
                       if (iii > Nsample )                       goto 320 

                    enddo
                 enddo
320              continue
              endif

              if (Nsample_lerg > 0) then

                 ! lerg 
                 iii=1
                 darkmass(:)=0.
                 dx=masses_lerg(2)-masses_lerg(1) ! this is same for herg and lerg
                 norm=dble(Nsample_lerg)

                 do i=1,Nrows_lerg
                    mu=real(lerg_p(i)*norm)
                    Nran=random_Poisson(mu, first)! integer from poisson distribution with mean mu
                    xmax=masses_lerg(i)+dx/2
                    xmin=masses_lerg(i)-dx/2
                    do j=1,Nran
                       ! jump to next element with flux under threshold
                       do while (samplex(iii) > lum_threshold)
                          iii=iii+1  
                       enddo
                       if (iii > Nsample )                       goto 330 
                       spot_dist(iii)=random_normal()*0.11+0.17 !FRI
                       Darkmass(iii)=rand()*(xmax-xmin)+xmin
                       iii=iii+1
                       if (iii > Nsample )                       goto 330 
                    enddo
                 enddo
330              continue
              endif

              !filling any remaining masses
              xmax=masses_lerg(Nrows_lerg)
              xmin=masses_lerg(1)

              do iii=1,Nsample
                 if (Darkmass(iii)==0.) then    
                    Darkmass(iii)=rand()*(xmax-xmin)+xmin
                    ! these sources are morphologically described as FRI
                    ! generate distance between the bright spots
                    spot_dist(iii)=random_normal()*0.11+0.17 !FRI
                 endif

                 if ((spot_dist(iii) <=0.) .or. (spot_dist(iii)>1.)) spot_dist(iii)=rand()  ! correct for out of bounds values - unif distribution

              enddo

           end select
           print*,'dark mass ends'

           !start coordinates and redhifts
           !initialize coordinate and redshifts 
           do i=1,Nsample
              latitudes(i)=(rand()-0.5)*sim_side
              longitudes(i)=(rand()-0.5)*sim_side
              z_gals(i)=rand()*(zhigh-zlow)+zlow

              mstar(i)=starmass(darkmass(i),z_gals(i)) !log10 stellar mass, Aversa et al. (2015)

              !himass with a stellar mass-hi mass relation 
              if (z_gals(i) <=0.525) then
                 ! HI mass modeled only up to z=0.5 bin

                 himass(i)=(0.71 +random_normal()*0.01)*mstar(i)+2.22 !Naluminsa et al. 2021

                 if (mstar(i) > 10.) himass(i)=(0.02+random_normal()*0.002)*mstar(i)+9.52 !Catinella et al. 2010


                 himass(i)=himass(i)+random_normal()*0.21 
                 
              endif


              ! refine the optical class following Koziel-Wierzbowska et al. 2020
              ! 98% of AGN are ellipticals, 2% are spiral

              draw=rand()

              if (draw >= 0.98) optclass(i)=2 !spiral  , Koziel-Wierzbowska et al. 2020

           enddo


           ! Convert intrinsic size to projected angular size
           ! taking into account view angle and redshift
           do i=1,Nsample
              cos_i=rand()*(theta_high-theta_low)+theta_low !uniform distribution sin(i)di
              angles(i)=acos(cos_i)
              sin_i=sin(angles(i))
              z_i=dble(z_gals(i))
              dim=sizes_3d(i)*sin_i/1000. ! size in Mpc corrected for view angle. 
              sizes(i)=theta_p(dim,z_i)
              current_count=current_count+1.
              ids(i)=current_count
           enddo

           print*,'coordinates'
           print*,minval(latitudes),maxval(latitudes)
           print*,minval(longitudes),maxval(longitudes)


           !preparing to output the data in catalogue format
           jstart=Nsample_old+1


           if (.not. allocated(catout)) then
              !start filling the output catalogue array

              allocate(catout(Ncat_common,Nsample),stat=iostat)
              if (iostat/=0) then
                 print*,'sampler: allocation error'
                 stop
              endif

           else
              !copy the old part of the output catalogue and add the new part

              allocate(catout_copy(Ncat_common,Nsample_old))

              catout_copy=catout
              deallocate(catout)
              allocate(catout(Ncat_common,Nsample_old+Nsample))
              catout(:,1:Nsample_old)=catout_copy
              deallocate(catout_copy)


           endif
           catout(1,jstart:jstart+Nsample-1)=ids(1:Nsample)        !
           catout(2,jstart:jstart+Nsample-1)=samplex(1:Nsample)        !lum_1.4 GHz
           catout(3,jstart:jstart+Nsample-1)=-100. !logSFR
           catout(4:nfreq_out+2,jstart:jstart+Nsample-1)=radioflux(4:Nfreq,:)  ! total intensity
           catout(nfreq_out+4:2*nfreq_out+3,jstart:jstart+Nsample-1)=polaflux(4:Nfreq,:)  !polarization
           catout(2*nfreq_out+4,jstart:jstart+Nsample-1)=darkmass      ! dark mass
           catout(2*nfreq_out+5,jstart:jstart+Nsample-1)=mstar         ! stellar mass
           catout(2*nfreq_out+6,jstart:jstart+Nsample-1)=himass        ! hi mass
           catout(2*nfreq_out+7,jstart:jstart+Nsample-1)=latitudes     !cartesian coordinates - to be projected on the sphere by wrapper
           catout(2*nfreq_out+8,jstart:jstart+Nsample-1)=longitudes
           catout(2*nfreq_out+9,jstart:jstart+Nsample-1)=0. ! spherical coordinates - to be filled by wrapper
           catout(2*nfreq_out+10,jstart:jstart+Nsample-1)=0. 
           catout(2*nfreq_out+11,jstart:jstart+Nsample-1)=z_gals       ! redshift
           catout(2*nfreq_out+12,jstart:jstart+Nsample-1)=sizes       ! projected angular size
           catout(2*nfreq_out+13,jstart:jstart+Nsample-1)=-100.     !inclinations
           catout(2*nfreq_out+14,jstart:jstart+Nsample-1)=-100.     !axis ratio
           catout(2*nfreq_out+15,jstart:jstart+Nsample-1)=-100.      !bmaj (host)
           catout(2*nfreq_out+16,jstart:jstart+Nsample-1)=-100.       !bmin (host)
           catout(2*nfreq_out+17,jstart:jstart+Nsample-1)=-100.       !PA (of the host)
           catout(2*nfreq_out+18,jstart:jstart+Nsample-1)=spot_dist   ! distance between bright spots (Rs)
           catout(2*nfreq_out+19,jstart:jstart+Nsample-1)=ii+3        ! flag to identify population (FSRQ, BL-LAc, SS)
           catout(2*nfreq_out+20,jstart:jstart+Nsample-1)=optclass        ! flag to identify optical 1 elliptical, 2 spiral


           ! compute intrinsic luminosities and store them in the catalogue if requested
           if (save_lums ==1) then

              allocate(lums_save(Nfreq,Nsample),freq_rest(Nfreq),lums_i(Nfreq))

              do i=1,Nsample
                 zlum=z_gals(i)

                 conv_flux=1./(4.*pi*(lumr(zlum)*Mpc)**2)*1.e26*(1.+zlum) !L [erg/s/Hz] ->S [mJy]

                 lums_i=log10(radioflux(:,i)/conv_flux) !luminosities at the rest frequencies 



                 freq_rest=log10(frequencies*(1.+z_gals(i))) ! rest frequencies


                 !interpolate the luminosity as a function of frequency from rest frequencies to frequencies


                 do iii=4,Nfreq
                    nu_iii=log10(frequencies(iii))
                    l_iii=interpol_nosort(nu_iii,freq_rest,lums_i,Nfreq)
                    lums_save(iii,i)=l_iii

                 enddo


              enddo

              catout(2*nfreq_out+21:3*nfreq_out+20,jstart:jstart+Nsample-1)=lums_save(4:Nfreq,:)

              deallocate(lums_save,lums_i,freq_rest)

           endif

           Nsample_old=Nsample_old+Nsample ! size of the catalogue filled so far. It needs to be completed with SFGs

           !free memory
           deallocate(radioflux,polaflux,samplex,darkmass,mstar,himass,latitudes,&
                &longitudes,z_gals,sizes,sizes_3d,angles,spot_dist,optclass,ids,stat=iostat)
           if (iostat/=0) then
              print*,'sampler: deallocation error'
              stop
           endif


400        continue

        enddo ! end loop AGN populations

        deallocate(masses_lerg,masses_herg,lerg_p,herg_p)



        ! SFGs modelling starts here
        print*,'************************'
        print*,'SFGs: Processing redshift',z,' Range',zlow,zhigh
        print*,'************************'

        zstr=redshift_names(zi) 


        conv_flux=1./(4.*pi*(lumr(z)*Mpc)**2)*1.e26*(1+z) !L [erg/s/Hz] ->S [mJy]
        frequencies_rest=frequencies*(1.+z) !freq in the rest frame (for K correction)

        !z-distance conversion at this redshift
        deltaz=zhigh-zlow
        delta_parall=r(zhigh)-r(zlow)

        !z-angular distance conversion at this redshift

        volumetot=4.*pi/3.*(r(zhigh)**3-r(zlow)**3)
        volume=volumetot*skyfrac  !volume corresponding to FoV
        print*,'volume=',volume

        !relation between L14 and mass of dark halo, from abundance matching
        ! reading from a file
        SFR2Mh_filename=indir(1:l_indir)//'/TRECS_Inputs/AbundanceMatch/results_LSFR/L2mh_z'//zstr//'.txt' 
        Ncolumns=2
        nrows=rows_number(SFR2Mh_filename,1,nskip)
        Nfunction=nrows

        allocate(data(nrows,ncolumns),L14(nrows),masses_L14(nrows))

        call read_columns(SFR2Mh_filename,2,nrows,Ncolumns,nskip,data)
        L14=data(:,1)
        masses_L14=data(:,2)
        deallocate(data)
        minmass_cone=9.2 ! this is the minimum halo mass in the lightcone




        !relation between SFR and stellar mass, from Aversa et al. 2015
        ! reading from a file
        SFR2Mstar_filename=indir(1:l_indir)//'/TRECS_Inputs/results_SFR_Mstar/SFR2mstar_z'//zstr//'.txt' 
        Ncolumns=2
        nrows=rows_number(SFR2Mstar_filename,1,nskip)


        if(allocated(SFRtab)) deallocate(SFRTab,mstartab)
        allocate(data(nrows,ncolumns),SFRtab(nrows),mstartab(nrows),stat=iostat)
        if (iostat /=0) then
           print*,'Allocation error'
           stop
        endif
        nrows_sfrtab=nrows
        call read_columns(SFR2Mstar_filename,2,nrows,Ncolumns,nskip,data)
        SFRtab=data(:,1)
        mstartab=data(:,2)
        deallocate(data)




        !Reading SFR rate functins, from Mancuso et al. 
        SFR_filename=indir(1:l_indir)//'/TRECS_Inputs/SFRF/SFRF_z'//zstr//'.dat'  

        Ncolumns=5
        nrows=rows_number(SFR_filename,1,nskip)
        Nrows_sf=nrows

        allocate(data(nrows,ncolumns),x(nrows),px(nrows))
        call read_columns(SFR_filename,2,nrows,Ncolumns,nskip,data)
        !data(:,1) =log(SFR)
        !data(:,2) =log(phitot)
        !data(:,3) =log(phi UVgal)
        !data(:,4) =log(phi spheroids)
        !data(:,5) =log(phi lensed spheroids)


        x=data(:,1) !log(SFR)
        names(1)='UVgal'
        names(2)='spheroids'
        names(3)='lens_spheroids'

        ! morphological coefficients for UVgal, spheroids, lensed spheroids 
        squarealphas(1)=0.2**2. ! spiral
        squarealphas(2)=0.5**2. ! elliptical
        squarealphas(3)=0.5**2. ! elliptical

        do ii=1,3 !loop on SFR populations  
           ! information for limiting ram usage
           ! processing long files in chuncks of lenght buffer_size
           buffer_size=1000
           buffer_free=buffer_size
           jbuf=0 ! index to fill buffer 
           px=(10.**data(:,ii+2)) !phi(logSFR)

           call Ldust(frequencies_rest,dustsed,dif,ii,Ld) !dust luminosity is from a template. needs to be multiplied by sfr


           ! select the portion of SFR function that's relevant given the flux limit
           ! this avoids generating too many galaxies that are subsequently discarded.

           do iii=1,nrows
              sfr=10.**x(iii)

              ! use the tabulated values to find Mstar corresponding to given SFR following Aversa et al. 2015
              p=minloc(abs(x(iii)-SFRtab))
              m_star=mstartab(p(1))

              call Lsynch(frequencies_rest,sfr,m_star,Lsyn) !Lsyn with mass dependence
              call Lff(frequencies_rest,sfr,Lfree)
              test=(Lsyn(ilim)+Lfree(ilim)+Ld(ilim)*sfr)*conv_flux ! add dust SED
              if (test < fluxlim*1000.*0.95) px(iii)=0.  ! threshold below the flux limit to allow for scatter
           enddo

           integ=trapint(x,px) ! integral with trapezium rule for the PDF, to give the number of galaxies 

           print*,'number of '//names(ii)//'galaxies for Mpc**3',integ

           norm=volume*integ/sum(px) ! normalisation for histogram 
           Ngen_db=volume*integ
           N=Nrows_sf 

           allocate(poisson_numbers(N),stat=iostat)

           if (iostat /= 0.) then 
              print*,'deallocation error poisson numbers!'
           endif

           ! Poisson sampling to get number of object per SFR bin from the PDF
           do i=1,N
              mu=real(Px(i)*norm)
              Nran=random_Poisson(mu, first)! integer from poisson distribution with mean mu
              poisson_numbers(i)=Nran
           enddo
           Nsample=sum(poisson_numbers) ! this is the number of galaxies to be generated


           if (Nsample==0) then
              deallocate(poisson_numbers)
              goto 500
           endif
           if (Nsample < buffer_size) buffer_size=Nsample

           allocate(samplex(buffer_size),radioflux(Nfreq,buffer_size),lums(buffer_size),mstar(buffer_size),stat=iostat)
           if (iostat/=0) then
              print*,'sampler: allocation error'
              stop
           endif

           samplex(:)=0.

           dx=abs(x(2)-x(1))

           Nsample_surv=0


           !loop over SFR function bins
           do i=1,N
              Nran=poisson_numbers(i) ! number of objects in this bin 
              xmax=x(i)+dx/2
              xmin=x(i)-dx/2

              allocate(radioflux_slice(Nfreq,Nran),samplex_slice(Nran),lums_slice(Nran),mstar_slice(Nran))

              do j=1,Nran
                 samplex_slice(j)=rand()*(xmax-xmin)+xmin  ! random generation of SRF
                 sfr=10.**samplex_slice(j)


                 !associate stellar mass to SFR with the lookup table
                 !an intepolation is used between tabulated values and the galaxy's actual SFR
                 mstar_slice(j)=interpol(samplex_slice(j),SFRtab,mstartab,nrows_sfrtab)




                 call Lsynch(frequencies_rest,sfr,mstar_slice(j),Lsyn) !Smith et al. 2021 eq 2 relation, with scatter
                 call Lff(frequencies_rest,sfr,Lfree) 

                 Lums_slice(j)=dlog10(Lsyn(i14)+Lfree(i14)+Ld(i14)*sfr) ! 1.4 GHz luminosity. Used later for size modelling
                 Radioflux_slice(:,j)=(Lsyn+Lfree+Ld*sfr)*conv_flux !flux= synchrotron + free-free+dust
                 !scatter is included in Lsynch 

                 if (Radioflux_slice(ilim,j)>= fluxlim*1000.) Nsample_surv= Nsample_surv+1  ! implement flux selection
              enddo




              if (buffer_free < Nsample_surv) then 
                 ! expand buffer
                 Print*,'resize buffer'
                 buffer_size_old=buffer_size
                 buffer_free_old=buffer_free
                 buffer_size=buffer_size_old+(buffer_size_old+Nsample_surv)*10
                 buffer_free=buffer_free+(buffer_size_old+Nsample_surv)*10

                 allocate(Radioflux_copy(Nfreq,buffer_size_old),samplex_copy(buffer_size_old))
                 allocate(lums_copy(buffer_size_old),mstar_copy(buffer_size_old))
                 radioflux_copy=radioflux
                 samplex_copy=samplex
                 lums_copy=lums
                 mstar_copy=mstar
                 deallocate(radioflux,samplex,lums,mstar)
                 allocate(Radioflux(Nfreq,buffer_size),samplex(buffer_size),lums(buffer_size),mstar(buffer_size))
                 Radioflux(:,1:buffer_size_old)=Radioflux_copy(:,:)
                 samplex(1:buffer_size_old)=samplex_copy(:)
                 lums(1:buffer_size_old)=lums_copy(:)
                 mstar(1:buffer_size_old)=mstar_copy(:)
                 deallocate(lums_copy,radioflux_copy,samplex_copy,mstar_copy)
              endif

              ! fill buffer
              do j=1,Nran
                 if (Radioflux_slice(ilim,j)>= fluxlim*1000.) then
                    jbuf=jbuf+1
                    Radioflux(:,jbuf)=Radioflux_slice(:,j)
                    samplex(jbuf)=samplex_slice(j)
                    lums(jbuf)=lums_slice(j)
                    mstar(jbuf)=mstar_slice(j)
                 endif
              enddo
              buffer_free=buffer_size-jbuf

              deallocate(samplex_slice,radioflux_slice,lums_slice,mstar_slice)

           enddo



           deallocate(poisson_numbers)



           Nsample=buffer_size-buffer_free
           print*,'Number of galaxies above flux limit',Nsample

           if (Nsample==0) then
              !skip resize and output catalogue if no object is found
              deallocate(radioflux,samplex,lums,mstar)
              goto 500 
           endif
           ! resize Radioflux to the final sample size
           if (buffer_free /=0) then
              print*,'resize final flux ' 
              allocate(radioflux_copy(Nfreq,buffer_size),samplex_copy(buffer_size),&
                   &lums_copy(buffer_size),mstar_copy(buffer_size),stat=iostat)
              if (iostat /=0) then
                 print*,'Allocation error'
                 stop
              endif
              radioflux_copy=radioflux
              samplex_copy=samplex
              lums_copy=lums
              mstar_copy=mstar
              deallocate(radioflux,samplex,lums,mstar,stat=iostat)
              if (iostat /=0) then
                 print*,'Deallocation error'
                 stop
              endif

              Nsample=buffer_size-buffer_free
              allocate(radioflux(Nfreq,Nsample),samplex(nsample),lums(nsample),mstar(nsample),stat=iostat)
              if (iostat /=0) then
                 print*,'Allocation error'
                 stop
              endif
              radioflux(:,:)=radioflux_copy(:,1:Nsample)
              samplex(:)=samplex_copy(1:Nsample)
              lums(:)=lums_copy(1:Nsample)
              mstar(:)=mstar_copy(1:Nsample)
              deallocate(radioflux_copy,samplex_copy,lums_copy,mstar_copy,stat=iostat)
              if (iostat /=0) then
                 print*,'Deallocation error'
                 stop
              endif
           endif


           print*,'SFRs and flux generated'


           !vectors for polarization model
           allocate(Polaflux(Nfreq,Nsample),inclinations(Nsample)) !polarized flux and view angle 

           ! generating view angle with sin(i) distribution
           do i=1,Nsample
              cos_i=rand()
              inclinations(i)=acos(cos_i)*180./pi
           enddo

           ! polarized flux
           call pola_SFGS(inclinations,frequencies,radioflux,polaflux)


           ! compute halo mass from sfr 
           allocate(Darkmass(Nsample),himass(Nsample)&
                &,latitudes(Nsample),longitudes(Nsample),z_gals(Nsample),&
                &sizes(Nsample),ellipticity1(Nsample),ellipticity2(Nsample)&
                &,thetas(Nsample),optclass(Nsample),qrat(Nsample),ids(Nsample),stat=iostat)
           if (iostat/=0) then
              print*,'sampler: allocation error'
              stop
           endif


           himass(:)=-100.
           optclass(:)=1
           if (ii ==1) optclass(:)=2 !spheroids -> elliptical, late-type -> spiral


           do i=1,Nsample
              dm_model=interpol(lums(i),L14,masses_L14,Nfunction)

              Darkmass(i)=dm_model ! dark halo mass 


              latitudes(i)=(rand()-0.5)*sim_side
              longitudes(i)=(rand()-0.5)*sim_side
              z_gals(i)=rand()*(zhigh-zlow)+zlow

              !if (z_gals(i) <=0.525) himass(i)=(0.9-z_gals(i)*0.4)*samplex(i)+9.15+0.075*z_gals(i)+random_normal()*0.3 !derived by comparing mass distribution of hi and continuum.
              ! obtained with aboundance-matching of M_HI with evolution as in Paul et al. 2023
              a_sfr=1.1-2.46*z_gals(i)+2.06*z_gals(i)**2.
              b_sfr=9.5-z_gals(i)
              if (z_gals(i) <=0.525) himass(i)=samplex(i)*a_sfr+b_sfr+random_normal()*0.2

           enddo


           !associate size to objects
           do i=1,Nsample
              z_i=dble(z_gals(i))
              dim=size_SF(Darkmass(i),z_i,ii)   ! physical size
              sizes(i)=theta_p(dim,z_i)           ! apparent size
              q=sqrt(squarealphas(ii)+(cos(inclinations(i)))**2.*(1.-squarealphas(ii))) ! axis ratio, linked to inclination and sub-population
              qrat(i)=q
              ellipticity1(i)=sqrt(sizes(i)**2./q) ! apparent bmaj
              ellipticity2(i)=q*ellipticity1(i)     ! apparent bmin
              Thetas(i)=rand()*360. !PA in degs
              current_count=current_count+1.
              ids(i)=current_count
           enddo



           !preparing to output the data in catalogue format
           jstart=Nsample_old+1
           if (.not. allocated(catout)) then

              !start filling the output catalogue array
              allocate(catout(Ncat_common,Nsample),stat=iostat)
              if (iostat/=0) then
                 print*,'sampler: allocation error'
                 stop
              endif

           else

              !copy the old part of the output catalogue and add the new part

              allocate(catout_copy(Ncat_common,Nsample_old))
              catout_copy=catout
              deallocate(catout)

              allocate(catout(Ncat_common,Nsample_old+Nsample))
              catout(:,1:Nsample_old)=catout_copy
              deallocate(catout_copy)
           endif


           catout(1,jstart:jstart+Nsample-1)=ids !Unique ID
           catout(2,jstart:jstart+Nsample-1)=-100. !Lum1400
           catout(3,jstart:jstart+Nsample-1)=samplex(1:Nsample)        !logSFR
           catout(4:nfreq_out+3,jstart:jstart+Nsample-1)=radioflux(4:Nfreq,:)  ! total intensity
           catout(nfreq_out+4:2*nfreq_out+3,jstart:jstart+Nsample-1)=polaflux(4:Nfreq,:)  !polarization
           catout(2*nfreq_out+4,jstart:jstart+Nsample-1)=darkmass      ! dark mass
           catout(2*nfreq_out+5,jstart:jstart+Nsample-1)=mstar      ! stellar mass
           catout(2*nfreq_out+6,jstart:jstart+Nsample-1)=himass      ! HI mass
           catout(2*nfreq_out+7,jstart:jstart+Nsample-1)=latitudes     !cartesian coordinates - to be projected on the sphere by wrapper
           catout(2*nfreq_out+8,jstart:jstart+Nsample-1)=longitudes
           catout(2*nfreq_out+9,jstart:jstart+Nsample-1)=0. ! spherical coordinates - to be filled by wrapper
           catout(2*nfreq_out+10,jstart:jstart+Nsample-1)=0. 
           catout(2*nfreq_out+11,jstart:jstart+Nsample-1)=z_gals
           catout(2*nfreq_out+12,jstart:jstart+Nsample-1)=sizes  
           catout(2*nfreq_out+13,jstart:jstart+Nsample-1)=inclinations
           catout(2*nfreq_out+14,jstart:jstart+Nsample-1)=qrat
           catout(2*nfreq_out+15,jstart:jstart+Nsample-1)=ellipticity1  !bmaj
           catout(2*nfreq_out+16,jstart:jstart+Nsample-1)=ellipticity2  !bmin
           catout(2*nfreq_out+17,jstart:jstart+Nsample-1)=thetas !PA
           catout(2*nfreq_out+18,jstart:jstart+Nsample-1)=-100. !spot dist
           catout(2*nfreq_out+19,jstart:jstart+Nsample-1)=ii !SGSs
           catout(2*nfreq_out+20,jstart:jstart+Nsample-1)=optclass !early/late-type

           ! compute intrinsic luminosities and store them in the catalogue if requested
           if (save_lums ==1) then




              allocate(lums_save(Nfreq,Nsample),freq_rest(Nfreq),lums_i(Nfreq))

              do i=1,Nsample
                 zlum=z_gals(i)

                 conv_flux=1./(4.*pi*(lumr(zlum)*Mpc)**2)*1.e26*(1.+zlum) !L [erg/s/Hz] ->S [mJy]

                 lums_i=log10(radioflux(:,i)/conv_flux) !luminosities at the rest frequencies 



                 freq_rest=log10(frequencies*(1.+z_gals(i))) ! rest frequencies


                 !interpolate the luminosity as a function of frequency from rest frequencies to frequencies


                 do iii=4,Nfreq
                    nu_iii=log10(frequencies(iii))
                    l_iii=interpol_nosort(nu_iii,freq_rest,lums_i,Nfreq)
                    lums_save(iii,i)=l_iii

                 enddo


              enddo


              catout(2*nfreq_out+21:3*nfreq_out+20,jstart:jstart+Nsample-1)=lums_save(4:Nfreq,:)

              deallocate(lums_save,lums_i,freq_rest)
           endif

           Nsample_old=Nsample_old+Nsample ! size of the catalogue filled so far

           deallocate(radioflux,polaflux,inclinations,samplex,darkmass,himass,&
                &latitudes,longitudes,z_gals,sizes,lums,mstar,ellipticity1,&
                &ellipticity2,optclass,qrat,ids,thetas,stat=iostat)
           if (iostat/=0) then
              print*,'sampler: deallocation error'

           endif
500        continue
        enddo ! loop on SFR populations


        ! writing catalogue to disk
        if (Nsample_old /=0) then 
           cat_filename=rawdir(1:l_rawdir)//'/catalogue_continuum_z'//zstr//'.fits'
           call write_catalogue_new(cat_filename,catout,Ncat_common,tagnames,tunit,tform)
           !write in summary file:
           write(42,*)zstr
           ! catalogue written
           print*,'done'
           deallocate(catout)
        endif
        deallocate(x,px,data)
        deallocate(masses_L14,l14)

     endif

  enddo

  ! close summary file (input data for wrapper)
  close(42)

  ! free memory
  deallocate(poladistr,x3,px3,alphascaling_48_20,alphascaling_14_48,alphascaling_015_14)
  deallocate(redshifts,redshift_names)
  deallocate(dustsed,dif)
  deallocate(Lsyn,Lfree,Ld)
  deallocate(frequencies,frequencies_rest,tagnames,tunit,tform)

  ! compute time for the run
  call date_and_time(values = values_time(:,2))

  values_time(:,1) = values_time(:,2) - values_time(:,1)
  clock_time =  (  (values_time(3,1)*24 &
       &           + values_time(5,1))*60. &
       &           + values_time(6,1))*60. &
       &           + values_time(7,1) &
       &           + values_time(8,1)/1000.
  PRINT*,"Total clock time [m]:",clock_time/60.
  PRINT*,"               "//code//" > Normal completion"
end program sampler
