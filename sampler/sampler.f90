!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program is part of the T-RECS code
! Author: A. Bonaldi
! see Bonaldi et al. (2018) MNRAS for more details
! generate samples of radio galaxies from models
! save outputs to files per redshift bin and per galaxy sub-population
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program sampler

  use healpix_types
  USE fitstools
  USE utilities
  USE pix_tools
  USE random_tools
  use random,only:random_poisson
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
  character(LEN=filenamelen)::paramfile,description,dummy,outdir
  real(dp),allocatable::frequencies(:),frequencies_rest(:),Lsyn(:),Lfree(:),Ld(:)
  real(sp),allocatable::spec(:)
  real(dp)::nu,deltanu,sfr,mn,mx,volume,integ,volumetot,fom,fom_old
  character(LEN=filenamelen)::SFR_filename,SFR2Mh_filename,filename,freq_filename
  character(LEN=filenamelen)::filename1,filename2,filename3,filename4,filename5,filename6
  character(LEN=filenamelen)::filename7,filename8,filename9,filename10
  character(LEN=filenamelen)::AGN_filename,LERG_filename,HERG_filename
  character(LEN=filenamelen)::cone_filename,chline,filestat,cat_filename
  character(LEN=10),allocatable::tagnames(:)
  character(LEN=10)::names(3)
  CHARACTER(LEN=5) :: output,output2,tag
  CHARACTER(LEN=10) ::output3
  CHARACTER(LEN=80), DIMENSION(1:120) :: header
  CHARACTER(LEN=7)::line,col_str,zstr_long
  CHARACTER(LEN=4),allocatable::redshift_names(:)
  CHARACTER(LEN=7),allocatable::redshift_names_lum(:)
  CHARACTER(LEN=4) ::zstr

  !single precision variables
  !real(sp)::center_lat,center_lon
  real(sp)::z_min,z_max
  real(sp),allocatable::delta(:),work(:,:),diff(:)
  real(sp),allocatable::Radioflux(:,:),catout(:,:),samplex(:),lums(:),z_gals(:)
  real(sp),allocatable::Polaflux(:,:),inclinations(:),ellipticity1(:),ellipticity2(:),polafracs(:),lums_save(:,:)
  real(sp),allocatable::Radioflux_copy(:,:),samplex_copy(:),sizes(:),radioflux_slice(:,:),lums_slice(:)
  real(sp),allocatable::angles(:),sizes_3d(:),lums_copy(:),spot_dist(:)

  real(sp),allocatable::darkmass(:),darkmass_halo(:),latitudes(:),longitudes(:),cone(:,:),cone_old(:,:)
  real(sp),allocatable::samplex_run(:),redshifts(:),L14(:)
  real(sp),allocatable::masses_L14(:),redshifts_lum(:),samplex_slice(:),satellite_flag(:)

  real(sp),allocatable::masses_lerg(:),masses_herg(:),lerg_p(:),herg_p(:)
  REAL(SP) :: clock_time,coo_max,coo_min,dlat,dlon,dz,delta_perp,delta_parall
  REAL(SP) ::deltatheta,deltaz,mu,flux14,dmtol,mgal,ncone_binned
  real(sp)::dm_model,dm_model2,dm_best,lat,lon,maxflux(1),redshift,minmass,minmass_cone
  real(sp)::dim,lum_threshold,theta_high,theta_low,sin_i,logs,alpha_low
  real(sp)::alpha_high,minfreq,norm_f,dlogl,col,b,a,arg1,arg2,th,freq_norm
  real(sp)::minm(Nmbins),maxm(Nmbins),mc(1),halfmax(1),mw(1),sg_m(1)
  real(sp),parameter::sigma_alpha=2.35660,mean_alpha=0.609847 

  !double precision variables
  real(dp)::sim_area,skyfrac
  real(dp)::sim_side
  real(dp)::fluxlim,fluxlim_freq
  real(dp),allocatable::data(:,:),x(:),px(:),dustsed(:,:),dif(:)
  real(dp),allocatable::data2(:,:),x2(:),px2(:),x3(:),px3(:)
  real(dp),allocatable::alphascaling_lowf(:,:),alphascaling_highf(:,:)
  real(dp),allocatable::poladistr(:,:),alphascaling_lowf_pola(:,:),alphascaling_highf_pola(:,:)
  real(dp)::alphamed(3),z,alpha,conv_flux,zlow,zhigh,nfullsky,it,Ngen_db,z_i
  real(dp)::dx,xmin,xmax,val,test,thr,norm,over,dx2,lumlim
  real(dp)::zlow_lum,zhigh_lum,zlum,par_2,sigma_lin,L_stop_2

  !integer variables
  integer::no_AGN,no_SFG,save_lums,do_clustering,Nsample_binned,istart(1),iend(1),istart_i,iend_i,iii
  integer*8::Nsample,Nsample_lerg,Nsample_herg,join,k,n_hist,nside,ntiles,t
  integer::iseed,Nfreq,Nfreq_out,i,ii,nrows,Ncolumns,Ncolumns_lf,nskip,Ngen,Ngen2,jj,kk
  integer::buffer_size,buffer_free,jbuf,buffer_free_old,buffer_size_old,l_outdir
  integer::count,Ncat_sfr,Ncat_agn,nrows_lerg,nrows_herg,seed(34),nrows_old,nrows_lf
  integer::Nsample_old,Nfunction,Nsample_surv,Nhaloes,Nsample_mass,nreds_out
  integer::l,j,nreds,nreds2,nreds_cone,zi,nloop,ist,iostat,nreds_lum,islice
  integer(4) :: ic4, crate4, cmax4,ni,sum_plus,p14(1),i14,ilim,i48,p(1),p_i,try
  INTEGER, DIMENSION(8,2) :: values_time
  integer::N,N3,Nran,reason,iunit
  integer,allocatable::indx(:,:),indx2(:,:),poisson_numbers(:)
  logical::first=.true.
  !types variables
  TYPE(paramfile_handle) :: handle

  save iseed


  call date_and_time(values = values_time(:,1))
  call system_clock(count=ic4, count_rate=crate4, count_max=cmax4)
  iseed=ic4

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
  !************************************                                                                                        

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
  fluxlim = parse_double(handle, 'fluxlim', default=10.d-9, descr=description)

  description = concatnl( &
       & " Enter the frequency at which the fluxlimit is imposed [MHz]: ")
  fluxlim_freq = parse_double(handle, 'fluxlim_freq', default=1400.d0, descr=description)

  description = concatnl( &
       & " Do you want to output the luminosities (0=no, 1=yes)? ")
  save_lums = parse_int(handle, 'save_lums', default=0, vmax=1, descr=description)

  description = concatnl( &
       & " Do you want to simulate clustering (only for a max size=5, 0=no, 1=yes)?")
  do_clustering = parse_int(handle, 'do_clustering', default=0, vmax=1, descr=description)

  description = concatnl( &
       & " Enter the minimum redhift to consider: ")
  z_min = parse_real(handle, 'z_min', default=0., vmin=0., descr=description)

  description = concatnl( &
       & " Enter the maximum redshift to consider: ")
  z_max = parse_real(handle, 'z_max', default=8., vmin=0., descr=description)

  description = concatnl( &
       & " Enter the name of the the output file directory:")
  outdir=  parse_string(handle, 'outdir', default='.', descr=description)

  ! end reading input parameters


  !***************************************************
  !Global filenames for the input data to be read:
  !
  !    AGNs
  !
  !effective spectral indices for AGNs
  filename1='../../TRECS_Inputs/alphaeff/alphascaling_flat_1.4_4.8.txt' ! 1.4-4.8 GHz flat 
  filename2='../../TRECS_Inputs/alphaeff/alphascaling_steep_1.4_4.8.txt' ! 1.4-4.8 GHz steep 
  filename3='../../TRECS_Inputs/alphaeff/alphascaling_flat_4.8_20.0.txt' ! 4.8-20 GHz flat 
  filename4='../../TRECS_Inputs/alphaeff/alphascaling_steep_4.8_20.0.txt' ! 4.8-20 GHz steep 
  filename5='../../TRECS_Inputs/alphaeff/alphascaling_flat_pol_1.4_4.8.txt' ! 1.4-4.8 GHz flat polarization
  filename6='../../TRECS_Inputs/alphaeff/alphascaling_flat_pol_4.8_20.0.txt' ! 4.8-20 GHz flat polarization 

  !polarization fractions, Galluzzi et l. and Hales et al. 
  filename7='../../TRECS_Inputs/AGN_polafraction/Polfrac_AGNs_1e4_hales.dat' 

  !characteristic luminosity: redshift bins
  filename8='../../TRECS_Inputs/LF/CharLum/Characteristic_Luminosity/z_bins.dat'

  !Intrinsic AGN size distribution from DiPompeo et al.
  filename9='../../TRECS_Inputs/LF/AGN_sizes.dat'  

  !   SFGs
  !
  !dust SED
  filename10='../../TRECS_Inputs/LF/SEDdust.dat'!nu, UVgal, spheroids, lens_spheroids

  ! end filenames
  !**********************************************************

  !'Seeding random number generators'
  !getting seed for random number generation from the clock
  call system_clock(count=ic4, count_rate=crate4, count_max=cmax4)
  do i=1,12
     seed(i)=ic4*i/12.+iseed  ! combining two clock readings to fill 12 seeds
  enddo

  call random_seed(PUT=seed)  ! feeding to Poisson generator 


  !string formatting: eliminate spaces between path anf file name
  outdir=ADJUSTL(outdir)
  l_outdir=LEN_TRIM(outdir)


  !simulation area 
  sim_area=sim_side*sim_side
  skyfrac=sim_area/fullsky_area

  print*,'simulation area [square degrees]=',sim_area

  ! the lightcone for the clustering simulation if 5X5 degs. Stopping if asked to do clustering simulation bigger than this
  if ((sim_side > 5.) .and. (do_clustering==1)) then
     print*,'Error:'
     print*,'Clustering not supported for simulation size bigger than 5X5'
     stop
  endif


  !read input files with effective spectral indices for the AGN frequency behaviour
  !1.4-4.8 GHz range, flat spectrum

  Ncolumns=4
  n=rows_number(filename1,1,nskip)
  allocate(data(n,ncolumns),alphascaling_lowf(n,4)) !flux, alphaflat,alphaflat,alphasteep
  call read_columns(filename1,2,n,Ncolumns,nskip,data)

  do i=1,n
     alphascaling_lowf(i,1)=data(i,1)           !flux at 1.4 GHz
     alphascaling_lowf(i,2)=data(i,3)-data(i,4) !alpha+delta_alpha
     alphascaling_lowf(i,3)=data(i,3)-data(i,4) !alpha+delta_alpha
  enddo
  deallocate(data)

  !1.4-4.8 GHz range, steep spectrum
  Ncolumns=4
  n=rows_number(filename2,1,nskip)
  allocate(data(n,ncolumns)) !flux, alphaflat,alphasteep
  call read_columns(filename2,2,n,Ncolumns,nskip,data)

  do i=1,n
     if (alphascaling_lowf(i,1)/=data(i,1)) then
        print*,'Error: fluxes do not match!'
        stop
     endif           !flux at 1.4 GHz
     alphascaling_lowf(i,4)=data(i,3)-data(i,4) !alpha+delta_alpha
  enddo
  deallocate(data)

  !4.8-20 GHz range, flat spectrum
  Ncolumns=4
  n=rows_number(filename3,1,nskip)
  allocate(data(n,ncolumns),alphascaling_highf(n,4)) !flux, alphaflat,alphaflat,alphasteep
  call read_columns(filename3,2,n,Ncolumns,nskip,data)

  do i=1,n
     alphascaling_highf(i,1)=data(i,1)           !log flux at 4.8 GHz (Jy)
     alphascaling_highf(i,2)=data(i,3)-data(i,4) !alpha+delta_alpha
     alphascaling_highf(i,3)=data(i,3)-data(i,4) !alpha+delta_alpha
  enddo
  deallocate(data)

  !4.8-20 GHz range, steep spectrum
  Ncolumns=4
  n=rows_number(filename4,1,nskip)
  allocate(data(n,ncolumns)) !flux, alphaflat,alphasteep
  call read_columns(filename4,2,n,Ncolumns,nskip,data)

  do i=1,n
     if (alphascaling_highf(i,1)/=data(i,1)) then
        print*,'Error: fluxes do not match!'
        stop
     endif
     alphascaling_highf(i,4)=data(i,3)-data(i,4)!alpha+delta_alpha
  enddo
  deallocate(data)

  ! polarization effective indices
  !1.4-4.8 GHz, flat
  Ncolumns=4
  n=rows_number(filename5,1,nskip)
  allocate(data(n,ncolumns),alphascaling_lowf_pola(n,4)) !flux, alphaflat,alphaflat,alphasteep
  call read_columns(filename5,2,n,Ncolumns,nskip,data)

  do i=1,n
     alphascaling_lowf_pola(i,1)=data(i,1)   
     alphascaling_lowf_pola(i,2)=-0.1-data(i,4)
     alphascaling_lowf_pola(i,3)=-0.1-data(i,4)
  enddo
  deallocate(data)

  alphascaling_lowf_pola(i,4)=alphascaling_lowf(i,4) ! for steep spectrum use the total intensity frequecy spectrum in polarization

  !4.8-20 GHz, flat
  Ncolumns=4
  n=rows_number(filename6,1,nskip)
  allocate(data(n,ncolumns),alphascaling_highf_pola(n,4)) 
  call read_columns(filename6,2,n,Ncolumns,nskip,data)

  do i=1,n
     alphascaling_highf_pola(i,1)=data(i,1) 
     alphascaling_highf_pola(i,2)=-0.1-data(i,4)
     alphascaling_highf_pola(i,3)=-0.1-data(i,4)
  enddo
  deallocate(data)

  alphascaling_highf_pola(i,4)=alphascaling_highf(i,4) ! for steep spectrum use the total intensity frequecy spectrum in polarization



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
  !1400 and 4800 MHz and the frequency at which the flux cut is done, whether they are already 
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
  
  !AGN simulation follows. If no AGN simulation is needed this is skipped
  if (no_AGN /=0) goto 0101


  !*************************
  !AGN simulation starts
  !*************************


  !structure of the catalogue
  Ncat_agn=2*Nfreq_out+12 !number of catalogue columns: flux in total intensity and polarization, 
  !plus all other parameters
  !the additional items are:
  !lum,mass,x,y,lat,lon,redshift,size,angle,proj_size,rs,popflag
  if (save_lums ==1) Ncat_agn=3*Nfreq_out+12 !optionally saving intrinsic luminosities


  !creating tag names for the catalogue
  allocate(tagnames(Ncat_agn))
  j=1
  tagnames(j)='Lum1400'
  j=2
  do i=1,Nfreq_out
     write(output,"(i5)")int(frequencies(i+3))
     output=ADJUSTL(output)
     l=LEN_TRIM(output)
     dummy = 'I'//output(:l)
     tagnames(j)=dummy
     j=j+1
  enddo

  !polarised intensity
  do i=1,Nfreq_out
     write(output,"(i5)")int(frequencies(i+3))
     output=ADJUSTL(output)
     l=LEN_TRIM(output)
     dummy = 'P'//output(:l)
     tagnames(j)=dummy
     j=j+1
  enddo

  tagnames(j)='Mh'
  j=j+1
  tagnames(j)='x_coord'
  j=j+1
  tagnames(j)='y_coord'
  j=j+1
  tagnames(j)='latitude'
  j=j+1
  tagnames(j)='longitude'
  j=j+1
  tagnames(j)='redshift'
  j=j+1
  tagnames(j)='phys size'
  j=j+1
  tagnames(j)='angle'
  j=j+1
  tagnames(j)='size'
  j=j+1
  tagnames(j)='Rs'
  j=j+1
  tagnames(j)='PopFlag'
  j=j+1
  if (save_lums ==1) then
     !optional tag names for the luminosity
     !Luminosity
     do i=1,Nfreq_out
        write(output,"(i5)")int(frequencies(i+3))
        output=ADJUSTL(output)
        l=LEN_TRIM(output)
        dummy = 'L'//output(:l)
        tagnames(j)=dummy
        j=j+1
     enddo
  endif



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

  !reading the first characteristic luminosity file to get dimension of fiels for later
  AGN_filename='../../TRECS_Inputs/LF/CharLum/Characteristic_Luminosity/luminosity_0.01000.dat' 
  Ncolumns=6
  Ncolumns_lf=Ncolumns
  nrows=rows_number(AGN_filename,1,nskip)
  Nrows_lf=nrows
  n_hist=nrows


  !redshift slices for AGNS
  ! this are the central redhifts at which the Bonato et al. 1.4 GHz luminosity functions are provided.  
  ! the lightcone for AGNs has been binned in the same way. 
  !each redshift slice is processed independently. execution can be parallelised in redhift bins

  nreds=57 !SF gals until z=8 (where I have cone)
  nreds_out=nreds

  allocate(redshifts(nreds),redshift_names(nreds))

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
  names(1)='FSRQ'
  names(2)='BLLac'
  names(3)='SS_AGN'
  alphamed=(/-0.1,-0.1,-0.8/) ! median spectral index for the three populations

  do zi=1,nreds_out-1   ! Main redshift loop 

     z=redshifts(zi)

     if ((z_min <= z) .and. (z_max > z)) then    ! redshift slice with center z is processed

        !getting the size of the redshift slice
        delta_perp=1. 
        deltatheta=theta(delta_perp,z)
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
        LERG_filename='../../TRECS_Inputs/AbundanceMatch/results_AGNs/AGNprob_LERG_z'//zstr//'.txt' 
        HERG_filename='../../TRECS_Inputs/AbundanceMatch/results_AGNs/AGNprob_HERG_z'//zstr//'.txt' 

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

        if (do_clustering==1) then 
           ! preprocessing for AGN clustering:
           ! reading the slice of lightcone with dark halo masses and their position

           cone_filename='../../TRECS_Inputs/cones/cone_5X5_z'//zstr//'.txt_sort' 
           Ncolumns=4
           nrows=rows_number(cone_filename,1,nskip)
           Nhaloes=nrows

           if (nrows ==0) then 
              print*,'no file found!'
              stop
           endif

           allocate(cone(nrows,ncolumns))
           call read_columns_s(cone_filename,2,nrows,Ncolumns,nskip,cone)
           !reading done


           !check that the size of the cone is big enough
           if ((minval(cone(:,3)) >  -0.95*sim_side/2.) .or.  (maxval(cone(:,3)) <0.95*sim_side/2.  ) &
                &.or. (minval(cone(:,4)) >  -0.95*sim_side/2.) .or.  (maxval(cone(:,4)) <0.95*sim_side/2.)) then 
              print*,'Error: cone size too small!'
              print*,cone_filename
              stop
           endif

           ! if the cone is too big, extract a smaller cone from its centre
           if ((minval(cone(:,3)) <  -1.*sim_side/2.) .or.  (maxval(cone(:,3)) >sim_side/2.  ) .or. &
                &(minval(cone(:,4)) <  -1.*sim_side/2.) .or.  (maxval(cone(:,4)) >sim_side/2.)) then
              print*,'resizing the cone to size',sim_side

              !copy cone to old cone
              nrows_old=nrows
              allocate(cone_old(nrows_old,ncolumns))
              cone_old=cone
              deallocate(cone)
              ! count number of objects for the new cone
              nrows=0
              do i=1,nrows_old
                 if ((abs(cone_old(i,3)) <= sim_side/2.) .and. (abs(cone_old(i,4)) <= sim_side/2.)) nrows=nrows+1
              enddo
              Nhaloes=nrows
              print*,'new number of haloes',nrows

              allocate(cone(nrows,ncolumns))
              ii=0
              do i=1,nrows_old
                 if ((abs(cone_old(i,3)) <= sim_side/2.) .and. (abs(cone_old(i,4)) <= sim_side/2.)) then
                    ii=ii+1
                    cone(ii,:)=cone_old(i,:)
                 endif
              enddo

              deallocate(cone_old)
           endif
           ! end resize cone if too big

           ! printing info on the cone
           print*,'************'
           print*,'cone masses',minval(cone(:,1)),maxval(cone(:,1))
           print*,'cone redshifts',minval(cone(:,2)),maxval(cone(:,2))
           print*,'cone lats',minval(cone(:,3)),maxval(cone(:,3))
           print*,'cone lons',minval(cone(:,4)),maxval(cone(:,4))
           print*,'************'



           !indexing cone by mass to make mass matching quicker
           !creating Nmbins mass bins for the cone
           ncone_binned=real(Nhaloes)/real(Nmbins)
           nsample_binned=int(ncone_binned)+1  ! 

           ! computing min/max mass for each mass bin based on the mass distribution of the cone 
           !similar number of objects in each mass bin, for efficiency
           do i=1,Nmbins 
              istart=maxval((/(i-1)*nsample_binned,1/))
              iend=minval((/i*nsample_binned,Nhaloes/))
              istart_i=istart(1)
              iend_i=iend(1)
              minm(i)=cone(istart_i,1)
              maxm(i)=cone(iend_i,1)
              if (iend_i > Nhaloes) goto 112
           enddo

112        continue


        endif ! end if do_clustering=1


        do ii=1,3 !loop on AGN populations  
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
              theta_low=sin(5.*pi/180.)
              theta_high=sin(90.*pi/180.)
           case default
              theta_high=sin(5.*pi/180.)
              theta_low=sin(0.)
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

                 print*,'slice',redshifts_lum(islice)

                 zstr_long=redshift_names_lum(islice) !string tag 

                 conv_flux=1./(4.*pi*(lumr(zlum)*Mpc)**2)*1.e26*(1+zlum) !L [erg/s/Hz] ->S [mJy] 
                 frequencies_rest=frequencies*(1.+zlum) !freq in the rest frame (for K correction)

                 !read luminosities 
                 AGN_filename='../../TRECS_Inputs/LF/CharLum/Characteristic_Luminosity/luminosity_'//zstr_long//'.dat'  

                 call read_columns(AGN_filename,2,nrows_lf,Ncolumns_lf,nskip,data)

                 x=data(:,2) ! log of characteristic luminosity
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
                       samplex_slice(iii)=x(i) 
                       iii=iii+1
                    enddo

                 enddo

                 alpha=alphamed(ii) ! starting point for the spectral index 
                 !computing flux and counting number of objects above flux cut

                 Nsample_surv=0
                 do i=1,Nsample
                    !converting luminosity to flux, scale flux in frequency
                    freq_norm=1400.*(1.+randgauss_boxmuller(iseed)/1000.) ! to introduce some scatter in flux14 from characteristic luminosity
                    flux14=(10.**samplex_slice(i)*(frequencies_rest(i14)/freq_norm)**alpha)*conv_flux !no scatter

                    Radioflux_slice(:,i)=flux14*(frequencies/1400.)**alpha
                    spec=radioflux_slice(:,i)
                    !refine frequency scaling with the "effective indics" and adding scatter
                    call effective_index(i14,i48,ii,alpha,alphascaling_lowf&
                         &,alphascaling_highf,spec,frequencies)
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
           print*,'done'

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

           !poisson-sample from the polarization PDF but on the same Nsample objects already extracted
           call poisson_constrained(iseed,x3,Px3,0.,100.,polafracs,Nsample)

           do i=1,Nfreq
              polaflux(i,:)=polafracs/100.*radioflux(i,:)
           enddo

           alpha=alphamed(ii)
           !frequency scaling in polarization
           do i=1,Nsample
              spec=polaflux(:,i)
              call effective_index(i14,i48,ii,alpha,alphascaling_lowf_pola&
                   &,alphascaling_highf_pola,spec,frequencies)
              polaflux(:,i)=spec
           enddo
           deallocate(polafracs)
           !end polarization model

           !allocate arrays to store the other results
           allocate(Darkmass(Nsample),latitudes(Nsample),longitudes(Nsample),&
                &z_gals(Nsample),sizes_3d(Nsample),sizes(Nsample),angles(Nsample),spot_dist(Nsample))

           ! size and view angle from precomputed distributions
           spot_dist(:)=0. !Distance between two bright posts! 
           !initialised to 0 because BLLAC and FSRQ will not have two lobes

           !AGN size model
           !load file with intrinsic size distribution
           !from DiPompeo et al.

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
           dx2=abs(x2(2)-x2(1))
           N=Nrows

!!$           allocate(poisson_numbers(N),stat=iostat)
!!$if (iostat /= 0) then
!!$              print*,'sizes allocation failed 2'
!!$              stop
!!$           endif


           ! sample from size distribution
           call poisson_constrained(iseed,x2,Px2,0.,2010.,sizes_3d,Nsample)
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
                    Darkmass(iii)=ran_mwc(iseed)*(xmax-xmin)+xmin
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
                 if (Darkmass(iii)==0.) Darkmass(iii)=ran_mwc(iseed)*(xmax-xmin)+xmin
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
                    Darkmass(iii)=ran_mwc(iseed)*(xmax-xmin)+xmin
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
                 if (Darkmass(iii)==0.) Darkmass(iii)=ran_mwc(iseed)*(xmax-xmin)+xmin
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


              !              print*,'Lum range',minval(samplex),maxval(samplex)

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
                       ! gnenerate distance between the bright spots
                       spot_dist(iii)=randgauss_boxmuller(iseed)*0.18+0.62 !FRII
                       Darkmass(iii)=ran_mwc(iseed)*(xmax-xmin)+xmin
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
                       spot_dist(iii)=randgauss_boxmuller(iseed)*0.11+0.17 !FRI
                       Darkmass(iii)=ran_mwc(iseed)*(xmax-xmin)+xmin
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
                    Darkmass(iii)=ran_mwc(iseed)*(xmax-xmin)+xmin
                    ! these sources are morphologically described as FRI
                    ! gnenerate distance between the bright spots
                    spot_dist(iii)=randgauss_boxmuller(iseed)*0.11+0.17 !FRI
                 endif

                 if ((spot_dist(iii) <=0.) .or. (spot_dist(iii)>1.)) spot_dist(iii)=ran_mwc(iseed)  ! correct for out of bounds values - unif distribution

              enddo

           end select
           print*,'dark mass ends'

           !start coordinates and redhifts
           !initialize coordinate and redshifts 
           do i=1,Nsample
              latitudes(i)=(ran_mwc(iseed)-0.5)*sim_side
              longitudes(i)=(ran_mwc(iseed)-0.5)*sim_side
              z_gals(i)=ran_mwc(iseed)*(zhigh-zlow)+zlow
           enddo

           ! clustering: associate galaxies to dark haloes of a lightcone from a cosmo simulation
           if (do_clustering==1) then
              print*,'clustering starts'
              minmass=minval(Darkmass)
              print*,'number of haloes',Nhaloes
              print*,'number of galaxies',Nsample

              do i=1,Nsample
                 mgal=darkmass(i)

                 if ((mgal >=minm(1)) .and. (mgal<=maxm(Nmbins))) then ! otherwise a random position has been assigned already
                    !compare mass of the galaxy with DH masses to find a matching one

                    if (Nhaloes < 1000) then 
                       ! if this redhist slice as a few DHs, I don't bin them in mass
                       istart_i=0
                       iend_i=Nhaloes

                    else
                       ! if this redhist slice as a few DHs, I bin them in mass to speed-up the mass matching

                       do iii=1,Nmbins   ! determine in which bin mass the galaxy falls
                          if ((mgal >=minm(iii)) .and. (mgal <maxm(iii))) then
                             istart=maxval((/(iii-1)*nsample_binned,1/))
                             iend=minval((/iii*nsample_binned,Nhaloes/))
                             istart_i=istart(1)
                             iend_i=iend(1)
                          endif
                       enddo

                    endif

                    ! assign the galaxy to the closest halo mass in the mass bin
                    fom_old=abs(mgal-cone(istart_i,1)) !initialise distance betweem model mass and halo mass
                    p_i=istart_i 

                    do iii=istart_i+1,iend_i
                       fom=abs(mgal-cone(iii,1))

                       if (fom<fom_old) then
                          fom_old=fom
                          p_i=iii
                       endif
                    enddo
                    ! once galaxy is associeted to halo give it the redshift and the coordinates of the halo.
                    dm_best=cone(p_i,1)
                    if (dm_best /=flagvalue) then 
                       latitudes(i)=cone(p_i,3)
                       longitudes(i)=cone(p_i,4)
                       z_gals(i)=cone(p_i,2)
                       cone(p_i,:)=flagvalue
                    endif
                 endif
              enddo

              print*,'clustering ends'


           endif ! end clustering 

           ! Convert intrinsic size to projected angular size
           ! taking into account view angle and redshift
           do i=1,Nsample
              sin_i=ran_mwc(iseed)*(theta_high-theta_low)+theta_low
              angles(i)=asin(sin_i)
              z_i=dble(z_gals(i))
              dim=sizes_3d(i)*sin_i/1000. ! size in Mpc corrected for view angle
              sizes(i)=theta(dim,z_i)     ! angular size
           enddo

           print*,'coordinates'
           print*,minval(latitudes),maxval(latitudes)
           print*,minval(longitudes),maxval(longitudes)

           !preparing to output the data in catalogue format
           allocate(catout(Ncat_agn,Nsample),stat=iostat)
           if (iostat/=0) then
              print*,'sampler: allocation error'
              stop
           endif


           !store the results in the catalogue


!!$           do i=1,Nsample
!!$              catout(1,i)=samplex(i)        !lum_1.4 GHz
!!$              catout(2:nfreq_out+1,i)=radioflux(4:Nfreq,i)  ! total intensity
!!$              catout(nfreq_out+2:2*nfreq_out+1,i)=polaflux(4:Nfreq,i)  !polarization
!!$              catout(2*nfreq_out+2,i)=darkmass(i)      ! mass
!!$              catout(2*nfreq_out+3,i)=latitudes(i)      !cartesian coordinates - to be projected on the sphere by wrapper
!!$              catout(2*nfreq_out+4,i)=longitudes(i) 
!!$              catout(2*nfreq_out+5,i)=0. ! spherical coordinates - to be filled by wrapper
!!$              catout(2*nfreq_out+6,i)=0. 
!!$              catout(2*nfreq_out+7,i)=z_gals(i)        ! redshift
!!$              catout(2*nfreq_out+8,i)=sizes_3d(i)      !intrinsic size
!!$              catout(2*nfreq_out+9,i)=angles(i)        !view angle
!!$              catout(2*nfreq_out+10,i)=sizes(i)        ! projected angular size
!!$              catout(2*nfreq_out+11,i)=spot_dist(i)    ! distance between bright spots (Rs)
!!$              catout(2*nfreq_out+12,i)=ii+3        ! flag to identify population
!!$           enddo
           catout(1,:)=samplex(1:Nsample)        !lum_1.4 GHz
           catout(2:nfreq_out+1,:)=radioflux(4:Nfreq,:)  ! total intensity
           catout(nfreq_out+2:2*nfreq_out+1,:)=polaflux(4:Nfreq,:)  !polarization
           catout(2*nfreq_out+2,:)=darkmass      ! mass
           catout(2*nfreq_out+3,:)=latitudes     !cartesian coordinates - to be projected on the sphere by wrapper
           catout(2*nfreq_out+4,:)=longitudes
           catout(2*nfreq_out+5,:)=0. ! spherical coordinates - to be filled by wrapper
           catout(2*nfreq_out+6,:)=0. 
           catout(2*nfreq_out+7,:)=z_gals       ! redshift
           catout(2*nfreq_out+8,:)=sizes_3d     !intrinsic size
           catout(2*nfreq_out+9,:)=angles       !view angle
           catout(2*nfreq_out+10,:)=sizes       ! projected angular size
           catout(2*nfreq_out+11,:)=spot_dist   ! distance between bright spots (Rs)
           catout(2*nfreq_out+12,:)=ii+3        ! flag to identify population

           ! compute intrinsic luminosities and store them in the catalogue if requested
           if (save_lums ==1) then
              allocate(lums_save(Nfreq,Nsample))
              do i=1,Nsample
                 zlum=z_gals(i)
                 conv_flux=1./(4.*pi*(lumr(zlum)*Mpc)**2)*1.e26*(1+zlum) !L [erg/s/Hz] ->S [mJy]
                 lums_save(:,i)=log10(radioflux(:,i)/conv_flux)
              enddo

              catout(2*nfreq_out+13:3*nfreq_out+12,:)=lums_save(4:Nfreq,:)
              deallocate(lums_save)

           endif

           
           
           write(output,"(i5)")t
           output=ADJUSTL(output)
           l=LEN_TRIM(output)

           ! writing catalogue
           cat_filename=outdir(1:l_outdir)//'/catalogue_z'//zstr//'_'//trim(names(ii))//'.fits'
           
           call write_catalogue(cat_filename,catout,Ncat_agn,tagnames)
                      !free memory
         
           deallocate(catout,radioflux,polaflux,samplex,darkmass,latitudes,&
                &longitudes,z_gals,sizes,sizes_3d,angles,spot_dist,stat=iostat)
           if (iostat/=0) then
              print*,'sampler: deallocation error'
              stop
           endif
           

400        continue

        enddo
        deallocate(masses_lerg,masses_herg,lerg_p,herg_p)
        if (do_clustering==1) deallocate(cone)
     endif

     
  enddo

  ! end of AGN modelling
  print*,'finished AGNs'
  ! free memory
  deallocate(poladistr,x3,px3,alphascaling_highf,alphascaling_lowf)
  deallocate(alphascaling_highf_pola,alphascaling_lowf_pola)
  deallocate(redshifts,redshift_names)
  deallocate(tagnames)


0101 continue

  ! SFG model begins
  if (no_SFG /=0) goto 0202  ! skip SFG if needed
  !  nreds=67  !SF gals
  !  if (do_clustering==1)  

  nreds=57 !SF gals until z=8 (where I have cone)
  nreds_out=nreds

  allocate(redshifts(nreds),redshift_names(nreds))

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

  ! setting number of columns for the catalogue
  Ncat_sfr=2*Nfreq_out+11 !total intensity and polarisation
  if (save_lums ==1) Ncat_sfr=3*Nfreq_out+11 !total intensity, polarisation and luminisity
  !The additional fields are:
  !sfr,mass,x,y,lat,lon,redshift,size,e1,e2,popflag


  ! Assign tag names to all columns in the catalogue
  allocate(tagnames(Ncat_sfr))
  j=1
  tagnames(j)='logSFR'
  j=2
  do i=1,Nfreq_out
     write(output,"(i5)")int(frequencies(i+3))
     output=ADJUSTL(output)
     l=LEN_TRIM(output)
     dummy = 'I'//output(:l)
     tagnames(j)=dummy
     j=j+1
  enddo
  do i=1,Nfreq_out
     write(output,"(i5)")int(frequencies(i+3))
     output=ADJUSTL(output)
     l=LEN_TRIM(output)
     dummy = 'P'//output(:l)
     tagnames(j)=dummy
     j=j+1
  enddo
  tagnames(j)='Mh'
  j=j+1
  tagnames(j)='x_coord'
  j=j+1
  tagnames(j)='y_coord'
  j=j+1
  tagnames(j)='latitude'
  j=j+1
  tagnames(j)='longitude'
  j=j+1
  tagnames(j)='redshift'
  j=j+1
  tagnames(j)='size'
  j=j+1
  tagnames(j)='e1'
  j=j+1
  tagnames(j)='e2'
  j=j+1
  tagnames(j)='PopFlag'
  j=j+1
  if (save_lums ==1) then
     !Luminosity
     do i=1,Nfreq_out
        !write(output,"(i3)")i
        !write(output,"(f10.5)")frequencies(i)
        write(output,"(i5)")int(frequencies(i+3))
        output=ADJUSTL(output)
        l=LEN_TRIM(output)
        dummy = 'L'//output(:l)
        tagnames(j)=dummy
        j=j+1
     enddo
  endif

  tagnames(1)='logSFR'

  !read dust SEDs
  Ncolumns=4
  nrows=rows_number(filename10,1,nskip)
  allocate(dustsed(nrows,ncolumns),dif(nrows))
  call read_columns(filename10,2,nrows,Ncolumns,nskip,dustsed)


  dustsed(:,2)=10.**dustsed(:,2)
  dustsed(:,3)=10.**dustsed(:,3)
  dustsed(:,4)=10.**dustsed(:,4)

  ! main redshift loop
  do zi=1,nreds_out-1

     z=redshifts(zi)

     if ((z_min <= z) .and. (z_max > z)) then    ! redshift slice with center z is processed

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
        print*,'SFGs: Processing redshift',z,' Range',zlow,zhigh
        print*,'************************'

        zstr=redshift_names(zi) 

        conv_flux=1./(4.*pi*(lumr(z)*Mpc)**2)*1.e26*(1+z) !L [erg/s/Hz] ->S [mJy] 
        frequencies_rest=frequencies*(1.+z) !freq in the rest frame (for K correction)

        !z-distance conversion at this redshift
        deltaz=zhigh-zlow
        delta_parall=r(zhigh)-r(zlow)

        !z-angular distance conversion at this redshift
        delta_perp=1. 
        deltatheta=theta(delta_perp,z)

        volumetot=4.*pi/3.*(r(zhigh)**3-r(zlow)**3)
        volume=volumetot*skyfrac  !volume corresponding to FoV

        !relation between L14 and mass of dark halo, from abundance matching
        ! reading from a file
        SFR2Mh_filename='../../TRECS_Inputs/AbundanceMatch/results_LSFR/L2mh_z'//zstr//'.txt' 
        Ncolumns=2
        nrows=rows_number(SFR2Mh_filename,1,nskip)
        Nfunction=nrows

        allocate(data(nrows,ncolumns),L14(nrows),masses_L14(nrows))

        call read_columns(SFR2Mh_filename,2,nrows,Ncolumns,nskip,data)
        L14=data(:,1)
        masses_L14=data(:,2)
        deallocate(data)
        minmass_cone=9.2 ! this is the minimum halo mass in the lightcone

        ! input lightcone data for the clustering model 
        if (do_clustering==1) then 

           ! reading data for the lightcone redhsift slice 
           cone_filename='../../TRECS_Inputs/cones/cone_5X5_z'//zstr//'.txt_sort' 
           Ncolumns=4

           nrows=rows_number(cone_filename,1,nskip)
           Nhaloes=nrows

           if (nrows ==0) then 
              print*,'no file found!'
              stop
           endif

           allocate(cone(nrows,ncolumns))
           call read_columns_s(cone_filename,2,nrows,Ncolumns,nskip,cone)
           ! reading done


           !check that cone is big enough for the specified FoV
           if ((minval(cone(:,3)) >  -0.95*sim_side/2.) .or.  (maxval(cone(:,3)) <0.95*sim_side/2.  ) &
                &.or. (minval(cone(:,4)) >  -0.95*sim_side/2.) .or.  (maxval(cone(:,4)) <0.95*sim_side/2.)) then 
              print*,'Error: cone size too small!'
              print*,cone_filename
              stop
           endif

           ! if cone is bigger than the FoV resize it by selectinv a catout from the centre
           if ((minval(cone(:,3)) <  -1.*sim_side/2.) .or.  (maxval(cone(:,3)) >sim_side/2.  ) .or. &
                &(minval(cone(:,4)) <  -1.*sim_side/2.) .or.  (maxval(cone(:,4)) >sim_side/2.)) then

              print*,'resizing the cone to size',sim_side

              !copy cone to old cone
              nrows_old=nrows

              allocate(cone_old(nrows_old,ncolumns))
              cone_old=cone
              deallocate(cone)


              ! count number of objects for the new cone
              nrows=0
              do i=1,nrows_old
                 if ((abs(cone_old(i,3)) <= sim_side/2.) .and. (abs(cone_old(i,4)) <= sim_side/2.)) nrows=nrows+1
              enddo
              Nhaloes=nrows
              print*,'new number of haloes',nrows

              allocate(cone(nrows,ncolumns))
              ii=0
              do i=1,nrows_old
                 if ((abs(cone_old(i,3)) <= sim_side/2.) .and. (abs(cone_old(i,4)) <= sim_side/2.)) then
                    ii=ii+1
                    cone(ii,:)=cone_old(i,:)
                 endif

              enddo

              deallocate(cone_old)
           endif ! end resizing cone if too big

           print*,'************'
           print*,'cone masses',minval(cone(:,1)),maxval(cone(:,1))
           print*,'cone redshifts',minval(cone(:,2)),maxval(cone(:,2))
           print*,'cone lats',minval(cone(:,3)),maxval(cone(:,3))
           print*,'cone lons',minval(cone(:,4)),maxval(cone(:,4))
           print*,'************'

           !binning cone in mass to speed-up mass matching 
           ncone_binned=real(Nhaloes)/real(Nmbins)
           nsample_binned=int(ncone_binned)+1

           !           print*,'rounding',ncone_binned,nsample_binned

           do i=1,Nmbins 
              istart=maxval((/(i-1)*nsample_binned,1/))
              iend=minval((/i*nsample_binned,Nhaloes/))
              istart_i=istart(1)
              iend_i=iend(1)
              minm(i)=cone(istart_i,1)
              maxm(i)=cone(iend_i,1)
              if (iend_i > Nhaloes) goto 111
           enddo

111        continue

        endif ! end input lightcone file for the clustering model


        !Reading SFR rate functins, from Mancuso et al. 
        SFR_filename='../../TRECS_Inputs/SFRF/SFRF_z'//zstr//'.dat'  

        Ncolumns=5
        nrows=rows_number(SFR_filename,1,nskip)
        Nrows_lf=nrows

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




        do ii=1,3 !loop on SFR populations  
           ! information for limiting ram usage
           ! processing long files in chuncks of lenght buffer_size
           buffer_size=1000
           buffer_free=buffer_size
           jbuf=0 ! index to fill buffer 

           !print*,'beginning loop'
           px=(10.**data(:,ii+2)) !phi(logSFR)


           ! pezzo in test comincia qui
           ! quick calculation of max radioflux at the selection frequency to avoid generating galaxies which will be discarded
           call Ldust(frequencies_rest,dustsed,dif,ii,Ld) 
           do iii=1,nrows
              sfr=10.**x(iii)
              !flux= synchrotron + free-free+dust
              call Lsynch(frequencies_rest,sfr,Lsyn) !Lsyn is average value 
              call Lff(frequencies_rest,sfr,Lfree) 
              !syn+ff with a scatter (evolution relation by Magnelli et al. )
              test=3.5 !max value for the gaussian random number, giving maximum possible flux for the source

              delta=10.0000**(log10(Lsyn+Lfree)+test*0.4000+2.3500*(1.0000 -(1.0000 +z)**(-0.1200))) 
              if (minval(delta)<0.) delta(:)=0. 
              test=(delta(ilim)+Ld(ilim)*sfr)*conv_flux ! add dust SED
              if (test < fluxlim*1000.) px(iii)=0.
           enddo

           integ=trapint(x,px) ! integral with trapezium rule for the PDF, to give the number of galaxies 

           print*,'number of '//names(ii)//'galaxies for Mpc**3',integ

           norm=volume*integ/sum(px) ! normalisation for histogram 
           Ngen_db=volume*integ
           N=Nrows_lf 

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

           allocate(samplex(buffer_size),radioflux(Nfreq,buffer_size),lums(buffer_size))
           if (iostat/=0) then
              print*,'sampler: allocation error'
              stop
           endif

           samplex(:)=0.

           dx=abs(x(2)-x(1))

           ! generating SFR values for the galaxies from the PDF

           !           call Ldust(frequencies_rest,dustsed,dif,ii,Ld) 
           ! this if the same for the same redshift slice. need to be multiplied by sfr 
           Nsample_surv=0

           !print*,poisson_numbers,sum(poisson_numbers)

           do i=1,N
              Nran=poisson_numbers(i)
              xmax=x(i)+dx/2
              xmin=x(i)-dx/2

              allocate(radioflux_slice(Nfreq,Nran),samplex_slice(Nran),lums_slice(Nran))

              do j=1,Nran
                 samplex_slice(j)=ran_mwc(iseed)*(xmax-xmin)+xmin

                 sfr=10.**samplex_slice(j)
                 !flux= synchrotron + free-free+dust
                 call Lsynch(frequencies_rest,sfr,Lsyn) !Lsyn is average value 
                 call Lff(frequencies_rest,sfr,Lfree) 
                 !syn+ff with a scatter (evolution relation by Magnelli et al. )
                 delta=10.0000**(log10(Lsyn+Lfree)+randgauss_boxmuller(iseed)*0.4000+2.3500*(1.0000 -(1.0000 +z)**(-0.1200))) 
                 if (minval(delta)<0.) delta(:)=0. 
                 Lums_slice(j)=dlog10(delta(i14)+Ld(i14)*sfr) 
                 Radioflux_slice(:,j)=(delta+Ld*sfr)*conv_flux ! add dust SED
                 if (Radioflux_slice(ilim,j)>= fluxlim*1000.) Nsample_surv= Nsample_surv+1  ! implement flux threshold
              enddo

              !              print*,'Number of galaxies above flux limit',Nsample_surv


              if (buffer_free < Nsample_surv) then 
                 ! expand buffer
                 Print*,'resize buffer'
                 buffer_size_old=buffer_size
                 buffer_free_old=buffer_free
                 buffer_size=buffer_size_old+(buffer_size_old+Nsample_surv)*10
                 buffer_free=buffer_free+(buffer_size_old+Nsample_surv)*10

                 allocate(Radioflux_copy(Nfreq,buffer_size_old),samplex_copy(buffer_size_old),lums_copy(buffer_size_old))
                 radioflux_copy=radioflux
                 samplex_copy=samplex
                 lums_copy=lums 
                 deallocate(radioflux,samplex,lums)
                 allocate(Radioflux(Nfreq,buffer_size),samplex(buffer_size),lums(buffer_size))
                 Radioflux(:,1:buffer_size_old)=Radioflux_copy(:,:)
                 samplex(1:buffer_size_old)=samplex_copy(:)
                 lums(1:buffer_size_old)=lums_copy(:)
                 deallocate(lums_copy,radioflux_copy,samplex_copy)
              endif

              ! fill buffer
              do j=1,Nran
                 if (Radioflux_slice(ilim,j)>= fluxlim*1000.) then
                    jbuf=jbuf+1
                    Radioflux(:,jbuf)=Radioflux_slice(:,j)
                    samplex(jbuf)=samplex_slice(j)
                    lums(jbuf)=lums_slice(j)
                 endif
              enddo
              buffer_free=buffer_size-jbuf

              deallocate(samplex_slice,radioflux_slice,lums_slice)

           enddo


           deallocate(poisson_numbers)



           Nsample=buffer_size-buffer_free
           print*,'Number of galaxies above flux limit',Nsample
           if (Nsample==0) then
              !skip resize and output catalogue if no object is found
              deallocate(radioflux,samplex,lums)
              goto 500 
           endif
           ! resize Radioflux to the final sample size
           if (buffer_free /=0) then
              print*,'resize final flux ' 
              allocate(radioflux_copy(Nfreq,buffer_size),samplex_copy(buffer_size),lums_copy(buffer_size),stat=iostat)
              if (iostat /=0) then
                 print*,'Allocation error'
                 stop
              endif
              radioflux_copy=radioflux
              samplex_copy=samplex
              lums_copy=lums
              deallocate(radioflux,samplex,lums,stat=iostat)
              if (iostat /=0) then
                 print*,'Deallocation error'
                 stop
              endif

              Nsample=buffer_size-buffer_free
              allocate(radioflux(Nfreq,Nsample),samplex(nsample),lums(nsample),stat=iostat)
              if (iostat /=0) then
                 print*,'Allocation error'
                 stop
              endif
              radioflux(:,:)=radioflux_copy(:,1:Nsample)
              samplex(:)=samplex_copy(1:Nsample)
              lums(:)=lums_copy(1:Nsample)
              deallocate(radioflux_copy,samplex_copy,lums_copy,stat=iostat)
              if (iostat /=0) then
                 print*,'Deallocation error'
                 stop
              endif
           endif

           !           print*,radioflux(ilim,:)
           !           stop
           print*,'SFRs and flux generated'
           !stop
           !vectors for polarization model
           allocate(Polaflux(Nfreq,Nsample),inclinations(Nsample)) !polarized flux and view angle 

           ! generating view angle with sin(i) distribution
           do i=1,Nsample
              sin_i=ran_mwc(iseed)
              inclinations(i)=asin(sin_i)*180./pi
           enddo

           ! polarized flux
           call pola_SFGS(inclinations,frequencies,radioflux,polaflux)


           ! compute halo mass from sfr 
           allocate(Darkmass(Nsample),Darkmass_halo(Nsample),latitudes(Nsample),longitudes(Nsample),z_gals(Nsample),sizes(Nsample),&
                &ellipticity1(Nsample),ellipticity2(Nsample),satellite_flag(Nsample),stat=iostat)
           if (iostat/=0) then
              print*,'sampler: allocation error'
              stop
           endif

           satellite_flag(:)=1 ! by default SFGs are as satellites




           ! galaxies as satellites to big haloes:
           ! read mass probability distribution

           filename='../../TRECS_Inputs/AbundanceMatch/results_satellites/Satelliteprob_z'//zstr//'.txt' 

           Ncolumns=2
           nrows=rows_number(filename,1,nskip)

           allocate(data2(nrows,ncolumns),x2(nrows),px2(nrows)) 

           call read_columns(filename,2,nrows,Ncolumns,nskip,data2)
           x2=data2(:,1)
           px2=data2(:,2)


           ! initialise mass for all galaxies as satellite galaxies
           call poisson_constrained(iseed,x2,Px2,5.,15.,Darkmass_halo,Nsample)
           !mass is that of the halo hosting a satellite
           deallocate(x2,px2,data2)



           do i=1,Nsample
              dm_model=interpol(lums(i),L14,masses_L14,Nfunction)
              satellite_flag(i)=dm_model/darkmass_halo(i) ! galaxy/halo mass ratio
              Darkmass(i)=dm_model ! halo mass to associate with BGC instead of satellite         
              satellite_flag(i)=dm_model/darkmass_halo(i) ! galaxy/halo mass ratio  
              if (dm_model >=minmass_cone) satellite_flag(i)=0 ! only galaxies with mass smaller that minimum halo mass in the lightcone are kept as satellites

              latitudes(i)=(ran_mwc(iseed)-0.5)*sim_side
              longitudes(i)=(ran_mwc(iseed)-0.5)*sim_side
              z_gals(i)=ran_mwc(iseed)*(zhigh-zlow)+zlow
           enddo


           Nsample_mass=sum(satellite_flag)
           print*,'number of satellite galaxies',Nsample_mass

           if (do_clustering==1) then ! clustering model
              print*,'clustering starts'

              minmass=minval(Darkmass)
              print*,'number of haloes',Nhaloes
              print*,'number of galaxies',Nsample



              do i=1,Nsample
                 mgal=darkmass(i)

                 ! mass matching if the galaxy's mass is in the range of masses of the lighcone
                 ! galaxies outside this range have been assigned random positions already
                 ! galaxies that are successfully matched to a halo will be associated the same position of the halo
                 if ((mgal >=minm(1)) .and. (mgal<=maxm(Nmbins))) then 

                    ! two possible implementations of the association of dark haloes to SFGs
                    ! case(0): galaxy always associated to the closest halo mass, first come first served
                    ! this means some galaxies end up being associated with halo masses quite different
                    ! case(1): galaxy associated to halo mass no more distant than dmtol
                    ! this is more consistent for all galaxies, but typically gets galaxies associated with smaller halo masses, due to the shape of the mass function
                    ! mass_scatter parameter controlling this is defined in valiable declaration and set to 1. 
                    ! change it to 0 to use the alternative implementation

                    dmtol=2.

                    select case(mass_scatter) 
                    case(0)
                       !first mass association implementation
                       ! if the redshift slice contains a few haloes, no binning in mass is necessary
                       if (Nhaloes < 1000) then 
                          ! no mass binning for a few haloes
                          istart_i=0
                          iend_i=Nhaloes

                       else
                          ! for big redshift slices, use mass bins to speed-up associating the galaxy to a halo of similar mass

                          do iii=1,Nmbins   ! determine in which bin mass the galaxy falls
                             if ((mgal >=minm(iii)) .and. (mgal <maxm(iii))) then
                                istart=maxval((/(iii-1)*nsample_binned,1/))
                                iend=minval((/iii*nsample_binned,Nhaloes/))
                                istart_i=istart(1)
                                iend_i=iend(1)

                             endif
                          enddo
                       endif


!!$                    ! assign the galaxy to the closest halo mass in the mass bin
                       fom_old=abs(mgal-cone(istart_i,1)) !initialise distance betweem model mass and halo mass
                       ! associate galaxy to the closest available mass
                       p_i=istart_i 
                       do iii=istart_i+1,iend_i
                          fom=abs(mgal-cone(iii,1))
                          if (fom<fom_old) then
                             fom_old=fom
                             p_i=iii
                          endif
                       enddo

                       dm_best=cone(p_i,1)
                       if (dm_best /=flagvalue) then 
                          latitudes(i)=cone(p_i,3)
                          longitudes(i)=cone(p_i,4)
                          z_gals(i)=cone(p_i,2)
                          cone(p_i,:)=flagvalue
                       endif


                    case(1)
                       ! second mass association implementation

                       ! allow some scatter between halo mass and galaxy mass

                       istart_i=0   ! no mass binning, considering the whole cone
                       iend_i=Nhaloes
                       fom=1.e20    ! variable initialization
                       count=0.
                       if (maxval(cone(:,1)) > flagvalue) then  ! this means not all haloes are used 

                          do while (fom > dmtol)
                             try=int(ran_mwc(iseed)*(iend_i-istart_i-1)+istart_i+1)

                             if ((try >= istart_i+1) .and. (try <=iend_i )) then
                                fom=abs(mgal-cone(try,1))
                                count=count+1 
                             endif

                             if (count > Nhaloes) goto 888 
                             ! abort mass association if a suitable halo is not found
                             ! in this case random coordinates have been already assigned
                          enddo
                       endif

                       p_i=try
                       dm_best=cone(p_i,1)

                       if ((fom<=dmtol) .and. (dm_best /=flagvalue)) then 

                          latitudes(i)=cone(p_i,3)
                          longitudes(i)=cone(p_i,4)
                          z_gals(i)=cone(p_i,2)
                          cone(p_i,:)=flagvalue

                       endif
888                    continue



                    end select
                 endif
              enddo

              print*,'clustering ends'


           endif ! end clustering 


           deallocate(satellite_flag,darkmass_halo)



           !associate size to objects
           do i=1,Nsample
              z_i=dble(z_gals(i))
              dim=size_SF(Darkmass(i),z_i,ii)   ! physical size
              sizes(i)=theta(dim,z_i)           ! apparent size
           enddo

           ! ellipticities (Tunbridge et al. 2016)
           N=100
           allocate(x2(N),px2(N)) ! distribution for |e|

           do i=1,N
              x2(i)=i*1./real(N)  !|e| ranging 0 to 1 sampled by N points
              b=0.19
              a=0.58
              arg1=pi*x2(i)/2.
              arg2=(2.*x2(i)/b)**a
              px2(i)=x2(i)*((cos(arg1))**2 *exp(-1.*arg2))
           enddo

           call poisson_constrained(iseed,x2,Px2,0.,1.,ellipticity1,Nsample)
           deallocate(x2,px2)


           do i=1,Nsample
              th=ran_mwc(iseed)*2.*pi*2. !2theta
              ellipticity2(i)=ellipticity1(i)*sin(th)
              ellipticity1(i)=ellipticity1(i)*cos(th)
           enddo
           ! end ellipticities

           !writing the catalogue
           allocate(catout(Ncat_sfr,Nsample))


           !saving the fields to catalogue columns
           catout(1,:)=samplex(1:Nsample)
           catout(2:nfreq_out+1,:)=radioflux(4:Nfreq,:)
           catout(nfreq_out+2:2*nfreq_out+1,:)=polaflux(4:Nfreq,:)
           catout(2*nfreq_out+2,:)=darkmass
           catout(2*nfreq_out+3,:)=latitudes
           catout(2*nfreq_out+4,:)=longitudes
           catout(2*nfreq_out+5,:)=latitudes*0.
           catout(2*nfreq_out+6,:)=longitudes*0.
           catout(2*nfreq_out+7,:)=z_gals
           catout(2*nfreq_out+8,:)=sizes  
           catout(2*nfreq_out+9,:)=ellipticity1  
           catout(2*nfreq_out+10,:)=ellipticity2  
           catout(2*nfreq_out+11,:)=ii !SGSs
           if (save_lums ==1) then
              allocate(lums_save(Nfreq,Nsample))

              do i=1,Nsample
                 zlum=z_gals(i)
                 conv_flux=1./(4.*pi*(lumr(zlum)*Mpc)**2)*1.e26*(1+zlum) !L [erg/s/Hz] ->S [mJy] 
                 lums_save(:,i)=log10(radioflux(:,i)/conv_flux)
              enddo
              catout(2*nfreq_out+12:3*nfreq_out+11,:)=lums_save(4:Nfreq,:)
              deallocate(lums_save)

           endif
           ! end saving fields to catalogue columns

           ! wrinting catalogue to disk
           cat_filename=outdir(1:l_outdir)//'/catalogue_z'//zstr//'_'//trim(names(ii))//'.fits'

           call write_catalogue(cat_filename,catout,Ncat_sfr,tagnames)
           ! catalogue written
           print*,'done'
           deallocate(catout,radioflux,polaflux,inclinations,samplex,darkmass,&
                &latitudes,longitudes,z_gals,sizes,lums,ellipticity1,ellipticity2,stat=iostat)
           if (iostat/=0) then
              print*,'sampler: deallocation error'

           endif
500        continue
        enddo ! loop on SFR populations

        deallocate(x,px,data)
        deallocate(masses_L14,l14)
        if (do_clustering==1) deallocate(cone)
        !stop

     endif

  enddo !loop on redshifts

  !free all memory
  deallocate(redshifts,redshift_names)
  deallocate(dustsed,dif)
  deallocate(Lsyn,Lfree,Ld)
  deallocate(frequencies,frequencies_rest,tagnames)

0202 continue

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
