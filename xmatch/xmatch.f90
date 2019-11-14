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



program xmatch

  use healpix_types
  USE fitstools
  USE utilities
  USE pix_tools
  USE random_tools
  use sampler_io
  USE paramfile_io, ONLY : paramfile_handle, parse_init, parse_int, &
       parse_string, parse_double, parse_lgt, concatnl
  USE extension, ONLY : getEnvironment, getArgument, nArguments
  implicit none

  real, parameter::sterad_deg=3283. 
  character(LEN=filenamelen)::paramfile,description,dummy
  character(LEN=filenamelen)::filename1,filename2,filenameX
  !  character(LEN=filenamelen),allocatable::catfiles(:)
  character*16,allocatable::tagnames1(:),tform1(:),tunit1(:)
  character*16,allocatable::tagnames2(:),tform2(:),tunit2(:)
  character*16,allocatable::tagnamesX(:),tformX(:),tunitX(:)
  character*16::tagname,tag1,tag2
  character*16::lat_tag,lon_tag,x_tag,y_tag
  real(sp),allocatable::data(:,:),column(:)!,sample1(:),sample2(:)
  real(sp),allocatable::data1(:,:),data2(:,:),dataX(:,:),binmins(:),binmaxs(:)
  real(sp)::scatter
  real(dp)::side1,side2
  CHARACTER(LEN=5) :: output,output2,tag
  CHARACTER(LEN=10) ::output3
  CHARACTER(LEN=80), DIMENSION(1:120) :: header
  CHARACTER(LEN=80)::line
  integer::selection,nmatched,j,ntotal
  integer(I4B),allocatable::data_indices1(:),flags1(:),data_indices2(:),flags2(:),sample1(:),sample2(:)
  integer(I4B),allocatable::populations1(:),populations2(:),matches(:),matching_sources(:)
  integer(I4B) :: nest,test(1)
  integer(I4B):: Nrows,Nrows1,Nrows2
  integer::Nfiles,i,Ncols1,Ncols2,Ncols,nbins
  integer::iostat,l,ii,n1,n2
  integer::i_tag1,i_tag2
  integer status,unit,rows,nrows_i,readwrite,blocksize,hdutype,tfields
  integer::bitpix,naxis,naxes ,jj,do_clustering
  integer(I4B), dimension(:),  allocatable :: listpix,listpix_copy(:)
  integer(I4B)::nside,nlist,ipix,ipix2
  INTEGER, DIMENSION(8,2) :: values_time
  REAL(SP) :: clock_time
  REAL(SP) :: lat_target_rad,lon_target_rad,lat_rad,lon_rad
  integer(4) :: ic4, crate4, cmax4,ni,iseed
  CHARACTER(LEN=*), PARAMETER :: code ="Xmatch"
  TYPE(paramfile_handle) :: handle
  character extname*16
  integer varidat, colnum,frow,felem

  !double precision ispline
  save iseed

  !getting seed from the clock
  call system_clock(count=ic4, count_rate=crate4, count_max=cmax4)
  iseed=ic4

  call date_and_time(values = values_time(:,1))

  !1)input file:
  !  -bins in redshift
  !  -number of objects


  if (nArguments() == 0) then
     paramfile=''
  else if (nArguments() == 1) then
     call getArgument(1,paramfile)
  else 
     print '("Usage: exe: [parameter file name]")'
     stop 1
  endif



!!$cat1=/home/a.bonaldi/local2/scratch/Bonaldi/Radio_srccnt/TESTS_newtrecs/HIrun/catalogue_HI
!!$cat2=/home/a.bonaldi/local2/scratch/Bonaldi/Radio_srccnt/TESTS_newtrecs/continuumrun/catalogue_continuum
!!$Ncross=1
!!$tag1=MHI
!!$catX=/home/a.bonaldi/local2/scratch/Bonaldi/Radio_srccnt/TESTS_newtrecs/cross_HIcontinuum/catalogue_HI_continuum
!!$side1=1.
!!$side2=1.
!!$selection=0
!!$scatter=0.1


  !!dopo il cross aggiungeve una colonna che dava provenienza del calalogo. era usata per fare xmatch dei soli oggetti non matchati in precedenza. in realta' avevo fatto una nuova opzione che invece di fare il cross faceva semplicemente il merge portando agn e sfg a stesso formato. potrei anche cambiare in sampler_continuum e dare ad AGN ed SFG lo stesso numero di colonne 



  handle = parse_init(paramfile)
  description = concatnl( &
       & " Enter the name of the first catalogue ")
  filename1=parse_string(handle,'cat1', default='', descr=description)

  handle = parse_init(paramfile)
  description = concatnl( &
       & " Enter the name of the second catalogue ")
  filename2=parse_string(handle,'cat2', default='', descr=description)

  handle = parse_init(paramfile)
  description = concatnl( &
       & " Enter the name of the second catalogue ")
  filenameX=parse_string(handle,'catX', default='', descr=description)

  handle = parse_init(paramfile)
  description = concatnl( &
       & " Enter the tag name of the first catalogue ")
  tag1=parse_string(handle,'tag1', default='', descr=description)

  handle = parse_init(paramfile)
  description = concatnl( &
       & " Enter the tag name of the first catalogue ")
  tag2=parse_string(handle,'tag2', default='', descr=description)


  description = concatnl( &
       & " Enter the simulation size: ")
  side1 = parse_double(handle, 'side1', default=5.d0, vmin=0.d0, descr=description)


  description = concatnl( &
       & " Enter the simulation size: ")
  side2 = parse_double(handle, 'side2', default=5.d0, vmin=0.d0, descr=description)


  description = concatnl( &
       & " Enter the scatter for the match: ")
  scatter = parse_double(handle, 'scatter', default=1.d-2, vmin=1.d-5, descr=description)

  description = concatnl( &
       & " which objects in the final catalogue?")
  selection = parse_int(handle, 'selection', default=0, vmax=2, descr=description)


  !Explore filename1 content
  call columns_catalogue(filename1,Ncols1,status)
  if (status ==0) then 
     print*,'number of columns in file '//trim(filename1),Ncols1
  else
     print*,'Problem reading file '//trim(filename1)
     stop
  endif



  !  Ncols1=9
  allocate(tagnames1(Ncols1),tunit1(Ncols1),tform1(Ncols1),stat=iostat)
  !tunit1(:)=''
  !tform1(:)='1E'


  call rows_catalogue(filename1,Ncols1,nrows1,tagnames1,tunit1)
!!$print*,tagnames1
!!$stop

  i_tag1=-1
  do ii=1,Ncols1
     tagname=tagnames1(ii)
     if (tagname ==tag1) i_tag1=ii
  enddo
  if (i_tag1 ==-1) then
     print*,'Error: tag '//tag1//'not found in '//filename1
     stop  
  else
     print*,'Tag '//tag1//' found in '//trim(filename1)//' as Column',i_tag1

  endif


  !Explore filename2 content
  call columns_catalogue(filename2,Ncols2,status)

  if (status ==0) then 
     print*,'number of columns in file '//trim(filename2),Ncols2
  else
     print*,'Problem reading file '//trim(filename2)
     stop
  endif

  allocate(tagnames2(Ncols2),tunit2(Ncols2),tform2(Ncols2),stat=iostat)
  tunit2(:)=''
  tform2(:)='1E'


  call rows_catalogue(filename2,Ncols2,nrows2,tagnames2,tunit2)

  i_tag2=-1
  do ii=1,Ncols2
     tagname=tagnames2(ii)
     if (tagname ==tag2) i_tag2=ii
  enddo
  if (i_tag2 ==-1) then
     print*,'Error: tag '//tag2//'not found in '//filename2
     stop  
  else
     print*,'Tag '//tag2//' found in '//trim(filename2)//' as Column',i_tag2
  endif



  ! set up Xmatched catalogue columns
  Ncols=Ncols1+Ncols2
  allocate(tagnamesX(Ncols),tunitX(Ncols),tformX(Ncols))



  do i=1,Ncols2
     tagname=tagnames2(i)
     do ii=1,Ncols1
        if (tagname ==tagnames1(ii)) tagnames2(i)=trim(tagnames2(i))//'_2'
     enddo
  enddo

  tagnamesX(1:Ncols1)=tagnames1
  tagnamesX(Ncols1+1:Ncols)=tagnames2
  tunitX(1:Ncols1)=tunit1
  !  tformX(1:Ncols1)=tform1
  tunitX(Ncols1+1:Ncols)=tunit2
  !  tformX(Ncols1+1:Ncols)=tform2

  tformX(:)='1E'

  allocate(data1(nrows1,Ncols1),data2(nrows2,Ncols2),stat=iostat)
  if (iostat /=0) then
     print*,'Error allocating filenames'
     stop
  endif

  call read_catalogue(filename1,Ncols1,Nrows1,data1)
  call read_catalogue(filename2,Ncols2,Nrows2,data2)

  !print*,Nrows1,Nrows2
  !stop
  !TODO:here reduce size if they are different

  if (Nrows1<Nrows2) then
     if(nrows1 <=10) then
        nbins=1
        allocate(binmins(1),binmaxs(1))
        binmins=minval(data1(:,i_tag1))-scatter
        binmaxs=maxval(data1(:,i_tag1))+scatter
     else
        !      nbins=int(nrows1/10.)
        nbins=10
        allocate(column(nrows1),binmins(nbins),binmaxs(nbins))
        column=data1(:,i_tag1)
        call get_bins(column,nbins,binmins,binmaxs)
        deallocate(column)   
     endif
  else

     if(nrows2 <=10) then
        nbins=1
        allocate(binmins(1),binmaxs(1))
        binmins=minval(data2(:,i_tag2))-scatter
        binmaxs=maxval(data2(:,i_tag2))+scatter
     else
        !   nbins=int(nrows2/10.)
        nbins=10
        print*,nbins
        !stop
        allocate(column(nrows2),binmins(nbins),binmaxs(nbins))
        column=data2(:,i_tag2)

        print*,minval(column),maxval(column)
        !stop
        call get_bins(column,nbins,binmins,binmaxs)

        deallocate(column)
     endif
  endif

  !bin both catalogues
  allocate(column(nrows1),flags1(nrows1),data_indices1(nrows1),populations1(nbins))
  column=data1(:,i_tag1)
  flags1(:)=0
  call assign_bins(column,flags1,binmins,binmaxs,nrows1,nbins,data_indices1,populations1)
  deallocate(column)

  allocate(column(nrows2),flags2(nrows2),data_indices2(nrows2),populations2(nbins))
  column=data2(:,i_tag2)
  flags2(:)=0
  call assign_bins(column,flags2,binmins,binmaxs,nrows2,nbins,data_indices2,populations2)
  deallocate(column)


  !print*,sum(populations1),sum(populations2)
  !stop


  allocate(matching_sources(nrows2))
  matching_sources(:)=-100.



  do i=1,nbins

     n1=populations1(i)
     n2=populations2(i)

     if ((n1/=0) .and. (n2 /=0)) then
        allocate(sample1(n1),sample2(n2),matches(n2))
        matches(:)=-100.
        call extract_subsample(data_indices1,nrows1,i,sample1)
        call extract_subsample(data_indices2,nrows2,i,sample2)


        call cross(data1(sample1,i_tag1),data2(sample2,i_tag2),matches,n1,n2,scatter)

        do ii=1,n2
           if (matches(ii) /=-100) matching_sources(sample2(ii))=sample1(matches(ii))
        enddo

        deallocate(sample1,sample2,matches)
     endif

  enddo


  nmatched=0
  do ii=1,nrows2
     if (matching_sources(ii) /=-100) nmatched=nmatched+1
  enddo


  print*,'Matched objects',nmatched,'out of',nrows2



  select case(selection)

  case(0)
     ! write content of both catalogues
     ntotal=nrows1+nrows2-nmatched
     print*,'Total objects',ntotal

     allocate(data(Ncols,ntotal))
     data(:,:)=-100.
     data(1:Ncols1,1:nrows1)=transpose(data1)
     !print*,data(1:nrows1,i_tag1)

     j=nrows1+1

     do i=1,nrows2
        if (matching_sources(i) /=-100) then
           data(Ncols1+1:Ntotal,matching_sources(i))=data2(i,1:Ncols2) 
        else
           data(Ncols1+1:Ntotal,j)=data2(i,1:Ncols2)
           j=j+1
        endif
     enddo

  case(1)
     ! write only content of first catalogue, with matches
     ntotal=nrows1
     print*,'Total objects',ntotal

     allocate(data(Ncols,ntotal))
     data(:,:)=-100.
     data(1:Ncols1,1:nrows1)=transpose(data1)
     !print*,data(1:nrows1,i_tag1)

     !  j=nrows1+1

     do i=1,nrows2
        if (matching_sources(i) /=-100) then
           data(Ncols1+1:Ntotal,matching_sources(i))=data2(i,1:Ncols2) 
!!$      else
!!$         data(Ncols1+1:Ntotal,j)=data2(i,1:Ncols2)
!!$         j=j+1
        endif
     enddo

  case(2)
     ! write only content  of second catalogue, with matches

     ntotal=nrows2
     print*,'Total objects',ntotal

     allocate(data(Ncols,ntotal))
     data(:,:)=-100.
     data(Ncols1+1:Ntotal,1:nrows2)=transpose(data2)
     !print*,data(1:nrows1,i_tag1)

     j=nrows1+1

     do i=1,nrows2
        if (matching_sources(i) /=-100) then
           data(1:Ncols1+1,i)=data1(matching_sources(i),1:Ncols1)
!!$      else
!!$         data(Ncols1+1:Ntotal,j)=data2(i,1:Ncols2)
!!$         j=j+1
        endif
     enddo

  end select



  call write_catalogue_new(filenameX,data,Ncols,tagnamesX,tunitX,tformX)


  call date_and_time(values = values_time(:,2))

  values_time(:,1) = values_time(:,2) - values_time(:,1)
  clock_time =  (  (values_time(3,1)*24 &
       &           + values_time(5,1))*60. &
       &           + values_time(6,1))*60. &
       &           + values_time(7,1) &
       &           + values_time(8,1)/1000.
  PRINT*,"Total clock time [m]:",clock_time/60.
  PRINT*,"               "//code//" > Normal completion"
end program xmatch
