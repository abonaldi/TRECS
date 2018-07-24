! module containing all I/O subroutines
! author A. Bonaldi

module sampler_io
  use healpix_types
  USE fitstools
  USE utilities
  USE pix_tools
  USE paramfile_io, ONLY : paramfile_handle, parse_init, parse_int, &
       parse_string, parse_double, parse_lgt, concatnl
contains


 subroutine rows_catalogue(filename,Ncol,irow,tags)
    !enquire length of catalogue file 
    !C     read and print data values from an ASCII or binary table

    integer status,unit,readwrite,blocksize,hdutype,ntable
    integer felem,nelems,nullj,diameter,nfound,irow,colnum
    real nulle,density
    character filename*200,nullstr*1,name*8
    character*16 ttype(Ncol),tform(Ncol),tunit(Ncol),tags(:)!,units(:)
    logical anynull


1   status=0

    !C     Get an unused Logical Unit Number to use to open the FITS file
2   call ftgiou(unit,status)

    !C     open the FITS file previously created by WRITEIMAGE
    readwrite=0
3   call ftopen(unit,filename,readwrite,blocksize,status)



    !C         read the TTYPEn keywords, which give the names of the columns
6   call ftgkns(unit,'TTYPE',1,Ncol,ttype,nfound,status)

    if (nfound==0) then 
       call ftmrhd(unit,1,hdutype,status) ! go to the next extension
       call ftgkns(unit,'TTYPE',1,Ncol,ttype,nfound,status)
    endif
    

    !C         read the data, one row at a time, and print them out
    felem=1
    nelems=1
    nullstr=' '
    nullj=0
    nulle=0.
  
    status=0
    irow=1
    do while (status ==0)
       colnum=1
9      call ftgcve(unit,colnum,irow,felem,nelems,nulle,density,anynull,status)

       irow=irow+1
    end do
    irow=irow-2 ! number of elements in table
    !print*,'Number of rows',irow

    status=0
    !C     close the file and free the unit number
10  call ftclos(unit, status)
    call ftfiou(unit, status)

    !C     check for any error, and if so print out error messages
11  if (status .gt. 0)call printerror(status)
    tags=ttype
!stop
  end subroutine rows_catalogue


  subroutine read_catalogue(filename,Ncol,Nrows,data)

    !C     read and print data values from an ASCII or binary table

    integer status,unit,readwrite,blocksize,hdutype,ntable
    integer felem,nelems,nullj,diameter,nfound,irow,colnum,Nrows
    real nulle,density,density2
    character filename*200,nullstr*1,name*8
    character*16 ttype(Ncol),tform(Ncol),tunit(Ncol)!
    logical anynull
    real(sp)::data(:,:)

1   status=0

    !C     Get an unused Logical Unit Number to use to open the FITS file
2   call ftgiou(unit,status)

    !C     open the FITS file previously created by WRITEIMAGE
    !filename='ATESTFILEZ.FITS'
    readwrite=0
3   call ftopen(unit,filename,readwrite,blocksize,status)



    !C         read the TTYPEn keywords, which give the names of the columns
6   call ftgkns(unit,'TTYPE',1,Ncol,ttype,nfound,status)
    if (nfound==0) then 
       call ftmrhd(unit,1,hdutype,status) ! go to the next extension
       call ftgkns(unit,'TTYPE',1,Ncol,ttype,nfound,status)
    endif

    felem=1
    nelems=1
    nullstr=' '
    nullj=0
    nulle=0.
    do irow=1,Nrows
       do icol=1,Ncol
          call ftgcve(unit,icol,irow,felem,nelems,nulle,density,anynull,status)
          data(irow,icol)=density
       enddo
    end do

    !C     close the file and free the unit number
10  call ftclos(unit, status)
    call ftfiou(unit, status)

    !C     check for any error, and if so print out error messages
11  if (status .gt. 0)call printerror(status)
  end subroutine read_catalogue

  ! writing catalogue in fits binary table format 
  subroutine write_catalogue(filename,data,Ncol,tagnames)
    real(sp),intent(in)::data(:,:)
    character(LEN=10),intent(in)::tagnames(:)
    integer status,unit,readwrite,blocksize,hdutype,tfields,nrows,Ncol
    integer varidat, colnum,frow,felem
    integer bitpix,naxis,naxes
    real density,flag
    character(LEN=filenamelen)::filename,extname
    character*16 ttype(Ncol),tform(Ncol),tunit(Ncol)!,name(6)
    real(sp),allocatable::column(:)


    nrows=size(data(1,:))
    allocate(column(nrows))

    status=0
    !C     Name of the FITS file to append the ASCII table to:

    !C     Get an unused Logical Unit Number to use to open the FITS file
    call ftgiou(unit,status)

    !C     open the FITS file, with write access
    readwrite=1

    blocksize=1
    call ftinit(unit,filename,blocksize,status)
    !write primary header
    bitpix=8
    naxis=0
    naxes=0
    call FTPHPS(unit,bitpix,naxis,naxes,status)

    !C     move to the last (2nd) HDU in the file
    call ftmahd(unit,1,hdutype,status)

    !C     append/create a new empty HDU onto the end of the file and move to it
    call ftcrhd(unit,status)

    !C     define parameters for the binary table (see the above data statements)
    tfields=Ncol
    extname='Catalogue'
    varidat=0

    tform(:)='1E'
    tunit(1:ncol)=(/ ''/)
    ttype=tagnames

    !C     write the required header parameters for the binary table
    call ftphbn(unit,nrows,tfields,ttype,tform,tunit,extname,varidat,status)



    !C     write catalogue column by column 

    frow=1
    felem=1

    i=1
    do while (i <= Ncol)
       column=sngl(data(i,:))
       call ftpcle(unit,i,frow,felem,nrows,column,status)  
       if (status .gt. 0) call printerror(status)
       i=i+1

    end do

    !C     close the FITS file and free the unit number
    call ftclos(unit, status)
    call ftfiou(unit, status)

    !C     check for any error, and if so print out error messages
    if (status .gt. 0)call printerror(status)
    deallocate(column)
  end subroutine write_catalogue

 
  function rows_number(filename,iunit,nskip)
! count number of rows in an ASCI file
    integer:: reason,iunit,nrows,rows_number,nskip
    character (LEN=200)::filename
    real(DP)::col

    close(iunit)
    open (UNIT=iunit,file=filename,status='unknown',form="formatted")
    print*,filename
    nrows=1
    nskip=0
    DO
       READ(iunit,*,IOSTAT=Reason)  col
       IF (Reason > 0)  THEN 
!          print*,'skipping line'
          nskip=nskip+1
       ELSE IF (Reason < 0) THEN
          exit
       ELSE
          nrows=nrows+1
       END IF
    END DO
    rows_number=nrows-1
    CLOSE (unit=iunit)
 !   print*,'Nrows=',rows_number
    return
  end function rows_number

  subroutine read_columns_s(filename,iunit,nrows,Ncolumns,nskip,data)
! read ASCII file column by column, load as single precision
    character (LEN=200)::filename
    integer:: reason,iunit,nrows,Ncolumns,i,nskip
    real(sp):: col1,col2,col3,col4,col5,col6,col7,col8,col9,col10
    real(sp)::data(:,:)

    data(:,:)=0.
    if (Ncolumns >10) then 
       print*,'Number of columns not supported!'
       stop
    endif
    close(iunit)
    open (UNIT=iunit,file=filename,status='unknown',form="formatted")

    do i=1,nskip
       READ(iunit,*,IOSTAT=Reason)  col1
    enddo


    SELECT CASE(Ncolumns)

    CASE(1)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1
          data(i,1)=col1
       ENDDO
    CASE(2)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2
          data(i,1)=col1
          data(i,2)=col2
       ENDDO
    CASE(3)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2,col3
          data(i,1)=col1
          data(i,2)=col2
          data(i,3)=col3
       ENDDO

    CASE(4)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2,col3,col4
          data(i,1)=col1
          data(i,2)=col2
          data(i,3)=col3
          data(i,4)=col4
       ENDDO
    CASE(5)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2,col3,col4,col5
          data(i,1)=col1
          data(i,2)=col2
          data(i,3)=col3
          data(i,4)=col4
          data(i,5)=col5
       ENDDO
    CASE(6)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2,col3,col4,col5,col6
          data(i,1)=col1
          data(i,2)=col2
          data(i,3)=col3
          data(i,4)=col4
          data(i,5)=col5
          data(i,6)=col6
       ENDDO

    CASE(7)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2,col3,col4,col5,col7
          data(i,1)=col1
          data(i,2)=col2
          data(i,3)=col3
          data(i,4)=col4
          data(i,5)=col5
          data(i,6)=col6
          data(i,7)=col7
       ENDDO
    CASE(8)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2,col3,col4,col5,col7,col8
          data(i,1)=col1
          data(i,2)=col2
          data(i,3)=col3
          data(i,4)=col4
          data(i,5)=col5
          data(i,6)=col6
          data(i,7)=col7
          data(i,8)=col8
       ENDDO

    CASE(9)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2,col3,col4,col5,col7,col8,col9
          data(i,1)=col1
          data(i,2)=col2
          data(i,3)=col3
          data(i,4)=col4
          data(i,5)=col5
          data(i,6)=col6
          data(i,7)=col7
          data(i,8)=col8
          data(i,9)=col9
       ENDDO

    CASE(10)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2,col3,col4,col5,col7,col8,col9,col10
          data(i,1)=col1
          data(i,2)=col2
          data(i,3)=col3
          data(i,4)=col4
          data(i,5)=col5
          data(i,6)=col6
          data(i,7)=col7
          data(i,8)=col8
          data(i,9)=col9
          data(i,10)=col10
       ENDDO
    END SELECT
    CLOSE (unit=iunit)
  end subroutine read_columns_s


 subroutine read_columns(filename,iunit,nrows,Ncolumns,nskip,data)
! read ASCII file column by column, load as double precision 
   character (LEN=200)::filename
    integer:: reason,iunit,nrows,Ncolumns,i,nskip
    real(dp):: col1,col2,col3,col4,col5,col6,col7,col8,col9,col10
    real(dp)::data(:,:)

    data(:,:)=0.
    if (Ncolumns >10) then 
       print*,'Number of columns not supported!'
       stop
    endif
    close(iunit)
    open (UNIT=iunit,file=filename,status='unknown',form="formatted")

    do i=1,nskip
       READ(iunit,*,IOSTAT=Reason)  col1
    enddo


    SELECT CASE(Ncolumns)

    CASE(1)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1
          data(i,1)=col1
       ENDDO
    CASE(2)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2
          data(i,1)=col1
          data(i,2)=col2
       ENDDO
    CASE(3)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2,col3
          data(i,1)=col1
          data(i,2)=col2
          data(i,3)=col3
       ENDDO

    CASE(4)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2,col3,col4
          data(i,1)=col1
          data(i,2)=col2
          data(i,3)=col3
          data(i,4)=col4
       ENDDO
    CASE(5)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2,col3,col4,col5
          data(i,1)=col1
          data(i,2)=col2
          data(i,3)=col3
          data(i,4)=col4
          data(i,5)=col5
       ENDDO
    CASE(6)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2,col3,col4,col5,col6
          data(i,1)=col1
          data(i,2)=col2
          data(i,3)=col3
          data(i,4)=col4
          data(i,5)=col5
          data(i,6)=col6
       ENDDO

    CASE(7)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2,col3,col4,col5,col7
          data(i,1)=col1
          data(i,2)=col2
          data(i,3)=col3
          data(i,4)=col4
          data(i,5)=col5
          data(i,6)=col6
          data(i,7)=col7
       ENDDO
    CASE(8)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2,col3,col4,col5,col7,col8
          data(i,1)=col1
          data(i,2)=col2
          data(i,3)=col3
          data(i,4)=col4
          data(i,5)=col5
          data(i,6)=col6
          data(i,7)=col7
          data(i,8)=col8
       ENDDO

    CASE(9)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2,col3,col4,col5,col7,col8,col9
          data(i,1)=col1
          data(i,2)=col2
          data(i,3)=col3
          data(i,4)=col4
          data(i,5)=col5
          data(i,6)=col6
          data(i,7)=col7
          data(i,8)=col8
          data(i,9)=col9
       ENDDO

    CASE(10)
       DO i=1,nrows
          READ(iunit,*,IOSTAT=Reason)  col1,col2,col3,col4,col5,col7,col8,col9,col10
          data(i,1)=col1
          data(i,2)=col2
          data(i,3)=col3
          data(i,4)=col4
          data(i,5)=col5
          data(i,6)=col6
          data(i,7)=col7
          data(i,8)=col8
          data(i,9)=col9
          data(i,10)=col10
       ENDDO
    END SELECT
    CLOSE (unit=iunit)
  end subroutine read_columns


end module sampler_io
