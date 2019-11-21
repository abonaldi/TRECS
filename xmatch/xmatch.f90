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


  character(LEN=filenamelen)::paramfile,description,dummy
  character(LEN=filenamelen)::filename1,filename2,filenameX
  character(LEN=filenamelen)::filename1_root,filename2_root,filenameX_root
  character*16,allocatable::tagnames1(:),tform1(:),tunit1(:)
  character*16,allocatable::tagnames2(:),tform2(:),tunit2(:)
  character*16,allocatable::tagnamesX(:),tformX(:),tunitX(:)
  character*16::redshift_names(67)
  character*16::tagname,tag1,tag2
  real(sp),allocatable::data(:,:),column(:),data_old(:,:)
  real(sp),allocatable::data1(:,:),data1_dummy(:,:),data2(:,:),dataX(:,:),binmins(:),binmaxs(:)
  real(sp)::scatter,init
  real(dp)::side1,side2
  integer::selection,nmatched,j,ntotal,i_x,i_y,nrows_old
  integer(I4B),allocatable::data_indices1(:),flags1(:),data_indices2(:),flags2(:),sample1(:),sample2(:)
  integer(I4B),allocatable::populations1(:),populations2(:),matches(:),matching_sources(:)
  integer(I4B) :: nest,test(1)
  integer(I4B):: Nrows,Nrows1,Nrows2
  integer::Nfiles,i,Ncols1,Ncols2,Ncols,nbins
  integer::iostat,l,ii,n1,n2,n1_actual,n2_actual,jj
  integer::i_tag1,i_tag2,status,status1,status2,nreds_out
  INTEGER, DIMENSION(8,2) :: values_time
  REAL(SP) :: clock_time
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
  init=rand(iseed) ! initalise random number generator seed

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


  handle = parse_init(paramfile)
  description = concatnl( &
       & " Enter the name of the first catalogue ")
  filename1_root=parse_string(handle,'cat1', default='', descr=description)

  handle = parse_init(paramfile)
  description = concatnl( &
       & " Enter the name of the second catalogue ")
  filename2_root=parse_string(handle,'cat2', default='', descr=description)

  handle = parse_init(paramfile)
  description = concatnl( &
       & " Enter the name of the second catalogue ")
  filenameX_root=parse_string(handle,'catX', default='', descr=description)

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




  redshift_names=(/'0.01','0.02','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40','0.45','0.50',&
       '0.55','0.60','0.65','0.70','0.75','0.80','0.85','0.90','0.95','1.00','1.20','1.40','1.60','1.80','2.00'&
       ,'2.20','2.40','2.60','2.80','3.00','3.20','3.40','3.60','3.80','4.00','4.20','4.40','4.60','4.80','5.00'&
       ,'5.20','5.40','5.60','5.80','6.00','6.20','6.40','6.60','6.80','7.00','7.20','7.40','7.60','7.80','8.00',&
       '8.20','8.40','8.60','8.80','9.00','9.20','9.40','9.60','9.80','10.0'/)

  nreds_out=57

  do jj=1,nreds_out 
     ! establish format of the catalogues and set the format of the crossmatched catalogue

     filename1=trim(filename1_root)//'_z'//trim(redshift_names(jj))//'.fits'
     filename2=trim(filename2_root)//'_z'//trim(redshift_names(jj))//'.fits'


     !Explore filename1 content
     call columns_catalogue(filename1,Ncols1,status1)
     call columns_catalogue(filename2,Ncols2,status2)

     if (status1+status2 ==0) then ! both files exist: can set the format for all files in the redshift list

        allocate(tagnames1(Ncols1),tunit1(Ncols1),tform1(Ncols1),stat=iostat)

        call rows_catalogue(filename1,Ncols1,nrows1,tagnames1,tunit1)

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
        tunitX(Ncols1+1:Ncols)=tunit2

        tformX(:)='1E'

        goto 100
     endif
  enddo

  if (status1+status2 /=0) then
     print*,'Error: specified files not found'
     stop
  endif


100 continue
  ! format set, now proceed with the crossmatch for all redshift slices


  do jj=1,nreds_out 
     ! establish format of the catalogues and set the format of the crossmatched catalogue

     filename1=trim(filename1_root)//'_z'//trim(redshift_names(jj))//'.fits'
     filename2=trim(filename2_root)//'_z'//trim(redshift_names(jj))//'.fits'
     filenameX=trim(filenameX_root)//'_z'//trim(redshift_names(jj))//'.fits'

     call columns_catalogue(filename1,Ncols1,status1)
     call columns_catalogue(filename2,Ncols2,status2)

     print*,'Files to Cross-match:'
     print*,filename1
     print*,filename2

     if ((status1 /=0) .and. (status2 /=0)) goto 200! do nothing for this redshift slice if none of the files exist

     if (status1==0) then 
        !if first file exisis, read and eventually resize

        call rows_catalogue(filename1,Ncols1,nrows1,tagnames1,tunit1)

        allocate(data1(nrows1,Ncols1),stat=iostat)
        if (iostat /=0) then
           print*,'Error allocating first catalogue data'
           stop
        endif
        call read_catalogue(filename1,Ncols1,Nrows1,data1)
        if (side2 < side1) then
           ! reduce size of data1

           print*,'Reducing size of first catalogue to',side2
           nrows_old=nrows1
           allocate(data_old(nrows_old,Ncols1))
           data_old=data1
           deallocate(data1)

           !find columns contaning the coordinates
           i_x=-1
           i_y=-1
           do ii=1,Ncols1
              tagname=tagnames1(ii)
              if (tagname =='x_coord') i_x=ii
              if (tagname =='y_coord') i_y=ii
           enddo
           if (i_x ==-1) then
              print*,'Error: tag x_coord not found in '//filename1
              stop  
           else
              print*,'Tag x_coord found in '//trim(filename1)//' as Column',i_x
           endif
           if (i_y ==-1) then
              print*,'Error: tag y_coord not found in '//filename1
              stop  
           else
              print*,'Tag y_coord found in '//trim(filename1)//' as Column',i_y
           endif


           ! count number of objects for the reduced size
           nrows1=0
           do i=1,nrows_old
              if ((abs(data_old(i,i_x)) <= side2/2.) .and. (abs(data_old(i,i_y)) <= side2/2.) ) nrows1=nrows1+1
           enddo


           allocate(data1(nrows1,Ncols1))
           ii=0
           do i=1,nrows_old
              if ((abs(data_old(i,i_x)) <= side2/2.) .and. (abs(data_old(i,i_y)) <= side2/2.)) then
                 ii=ii+1
                 data1(ii,:)=data_old(i,:)
              endif
           enddo
           deallocate(data_old)
           print*,'First catalogue resized to',nrows1
        endif
     endif

     if (status2==0) then 
        !if second file exists, read and eventually resize
        call rows_catalogue(filename2,Ncols2,nrows2,tagnames2,tunit2)
        allocate(data2(nrows2,Ncols2),stat=iostat)
        if (iostat /=0) then
           print*,'Error allocating second catalogue data'
           stop
        endif

        call read_catalogue(filename2,Ncols2,Nrows2,data2)

        if (side1 < side2) then
           ! reduce size of data2

           print*,'Reducing size of second catalogue to',side1
           nrows_old=nrows2
           allocate(data_old(nrows_old,Ncols2))
           data_old=data2
           deallocate(data2)

           !find columns contaning the coordinates
           i_x=-1
           i_y=-1
           do ii=1,Ncols2
              tagname=tagnames2(ii)
              if (tagname =='x_coord') i_x=ii
              if (tagname =='y_coord') i_y=ii
           enddo
           if (i_x ==-1) then
              print*,'Error: tag x_coord not found in '//filename2
              stop  
           else
              print*,'Tag x_coord found in '//trim(filename2)//' as Column',i_x
           endif
           if (i_y ==-1) then
              print*,'Error: tag y_coord not found in '//filename2
              stop  
           else
              print*,'Tag y_coord found in '//trim(filename2)//' as Column',i_y
           endif


           ! count number of objects for the reduced size
           nrows2=0
           do i=1,nrows_old
              if ((abs(data_old(i,i_x)) <= side1/2.) .and. (abs(data_old(i,i_y)) <= side1/2.)) nrows2=nrows2+1
           enddo


           allocate(data2(nrows2,Ncols2))
           ii=0
           do i=1,nrows_old
              if ((abs(data_old(i,i_x)) <= side1/2.) .and. (abs(data_old(i,i_y)) <= side1/2.)) then
                 ii=ii+1
                 data2(ii,:)=data_old(i,:)
              endif
           enddo
           deallocate(data_old)
           print*,'Second catalogue resized to',nrows2
        endif

     endif

  if ((status1 /=0) .and. (status2==0)) then
     ! only one file exis, copy its content on the crossmatched catalogue
     if (selection==1) goto 200 
!!$ print*,Ncols1+1,Ncols,Ncols2
!!$stop
     ! write content of second catalogue to the crossmathed one
     ntotal=nrows2
     print*,'Total objects',ntotal

     allocate(data(Ncols+1,ntotal))
     data(:,:)=-100.
     data(Ncols1+1:Ncols,:)=transpose(data2)

  endif

  if ((status1 ==0) .and. (status2/=0)) then
     ! only one file exis, copy its content on the crossmatched catalogue

     if (selection==2) goto 200 
     ! write content of first catalogue to the crossmathed one
     ntotal=nrows1
     print*,'Total objects',ntotal

     allocate(data(Ncols,ntotal))
     data(:,:)=-100.
     data(1:Ncols1,1:ntotal)=transpose(data1)
  endif


  if (status1+status2==0) then ! both files exist: crssomatch
     ! choices for only one file existing or both or none

     if (Nrows1<Nrows2) then
        if(nrows1 <=10) then
           nbins=1
           allocate(binmins(1),binmaxs(1))
           binmins=minval(data1(:,i_tag1))*(1.-3.*scatter)
           binmaxs=maxval(data1(:,i_tag1))*(1.+3.*scatter)
        else
           !      nbins=int(nrows1/10.)
           nbins=10
           allocate(column(nrows1),binmins(nbins),binmaxs(nbins))
           column=data1(:,i_tag1)
           call get_bins(column,nbins,binmins,binmaxs)
           ! this will create bins overlap
           binmins=binmins*(1.-3.*scatter)
           binmaxs=binmaxs*(1.+3.*scatter)
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
           ! this will create bins overlap
           binmins=binmins*(1.-scatter)
           binmaxs=binmaxs*(1.+scatter)
           deallocate(column)
        endif
     endif

     allocate(matching_sources(nrows2),sample1(nrows1),sample2(nrows2),data1_dummy(nrows1,Ncols1))
     matching_sources(:)=-100.
     data1_dummy=data1

     do i=1,nbins
        !print*,'bin',i

        call extract_subsample_new(data1_dummy(:,i_tag1),nrows1,binmins(i),binmaxs(i),sample1,n1)
        call extract_subsample_new(data2(:,i_tag2),nrows2,binmins(i),binmaxs(i),sample2,n2)


        if ((n1/=0) .and. (n2 /=0)) call cross_new(data1_dummy(:,i_tag1)&
             &,sample1(1:n1),data2(:,i_tag2),sample2(1:n2),matching_sources,n1,n2,scatter,iseed)

     enddo

     deallocate(sample1,sample2,data1_dummy)

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
              data(Ncols1+1:Ncols,matching_sources(i))=data2(i,1:Ncols2) 
           else
              data(Ncols1+1:Ncols,j)=data2(i,1:Ncols2)
              ! TODO: here overwrite coordinates with the random ones?
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
              data(Ncols1+1:Ncols,matching_sources(i))=data2(i,1:Ncols2) 
           endif
        enddo

     case(2)
        ! write only content  of second catalogue, with matches

        ntotal=nrows2
        print*,'Total objects',ntotal

        allocate(data(Ncols,ntotal))
        data(:,:)=-100.
        data(Ncols1+1:Ncols,1:nrows2)=transpose(data2)

        ! TODO: here overwrite coordinates with the random ones?

        j=nrows1+1

        do i=1,nrows2
           if (matching_sources(i) /=-100) then
              data(1:Ncols1+1,i)=data1(matching_sources(i),1:Ncols1)
           endif
        enddo

     end select

     deallocate(binmins,binmaxs,matching_sources)
  endif

  print*,'Writing to:'
  print*,filenameX
  call write_catalogue_new(filenameX,data,Ncols,tagnamesX,tunitX,tformX)
  deallocate(data)

200 continue
  if (allocated(data1)) deallocate(data1)
  if (allocated(data2)) deallocate(data2)
  !deallocate(data1,data2,data,matching_sources,binmins,binmaxs)

enddo

call date_and_time(values = values_time(:,2))

values_time(:,1) = values_time(:,2) - values_time(:,1)
clock_time =  (  (values_time(3,1)*24 &
    &           + values_time(5,1))*60. &
    &           + values_time(6,1))*60. &
    &           + values_time(7,1) &
    &           + values_time(8,1)/1000.
PRINT*,"Total clock time [m]:",clock_time/60.
PRINT*,"               "//code//" > Normal completion"
end program
