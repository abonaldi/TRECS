
module interpolations
use healpix_types
contains




  function interpol(val,x,y,N)
    implicit none
    real(sp),intent(IN)::x(:),y(:)
    real(sp),intent(in)::val
    real(sp)::interpol
    integer::N,N2,i,ii,p(1),ny
    real(sp)::dif(N),dist(1)
    real(dp)::x0,x1,x2,y0,y1,y2,A,B,C,tiny,p0,p1,p2,yold

    tiny=1.e-3
    dif=abs(val-x)
    p=minloc(dif)
    dist=minval(dif)
    if (dist(1)<tiny) then
       interpol=y(p(1))
       !no interpolatin needed
    else 
       !quadratic interpolation. find closest three points
       if (p(1)==1) then 
          ii=p(1)+1
       else if (p(1)==N2) then
          ii=p(1)-1
       else
          ii=p(1)
       endif
       x0=x(ii-1) 
       x1=x(ii) 
       x2=x(ii+1) 
       y0=y(ii-1) 
       y1=y(ii) 
       y2=y(ii+1) 
  
       call quadratic(x0,x1,x2,y0,y1,y2,A,B,C)

       interpol=A*val**2.+B*val+C
  endif

end function interpol


  function interpol_nosort(val,x,y,N)
!does not assome the values to be unique nor sorted
    implicit none
    real(sp),intent(IN)::x(:),y(:)
    real(sp),intent(in)::val
    real(sp)::interpol_nosort
    integer::N,N2,i,ii,p(1),ny,p1,p2,p3
    real(sp)::dif(N),dist(1),flagval
    real(dp)::x0,x1,x2,y0,y1,y2,A,B,C,tiny!,p0,p1,p2,yold

    tiny=1.e-3
    dif=abs(val-x)
    p=minloc(dif)
    dist=minval(dif)
    flagval=maxval(dif)*2.
    !print*,dist
    !stop
    if (dist(1)<tiny) then
       interpol_nosort=y(p(1))
       !no interpolatin needed
    else

       !find the three closest points recursively
       p1=p(1)
       do i=1,N
          if (dif(i) ==dist(1)) dif(i)=flagval
       enddo
       p=minloc(dif)
       dist=minval(dif)
       p2=p(1)

       do i=1,N
          if (dif(i) ==dist(1)) dif(i)=flagval
       enddo
       p=minloc(dif)
       p3=p(1)



       !quadratic interpolation. find closest three points
       !if (p(1)==1) then 
       !   ii=p(1)+1
       !else if (p(1)==N2) then
       !   ii=p(1)-1
       !else
       !   ii=p(1)
       !endif
       x0=x(p1) 
       x1=x(p2) 
       x2=x(p3) 
       y0=y(p1) 
       y1=y(p2) 
       y2=y(p3) 

!!$       print*,10.**x0,y0
!!$       print*,10.**x1,y1
!!$       print*,10.**x2,y2
!!$       stop
       call quadratic(x0,x1,x2,y0,y1,y2,A,B,C)

       interpol_nosort=A*val**2.+B*val+C
  endif

     end function interpol_nosort

subroutine quadratic(x0,x1,x2,y0,y1,y2,A,B,C)
implicit none
!Computes quadratic coefficients 
real(dp),intent(IN)::x0,x1,x2,y0,y1,y2
real(dp),intent(OUT)::A,B,C
real(dp):: a0,a1,a2

a0=y0/(x0-x1)/(x0-x2)
a1=y1/(x1-x0)/(x1-x2)
a2=y2/(x2-x0)/(x2-x1)

A=a0+a1+a2
B=-(a0*(x1+x2)+a1*(x0+x2)+a2*(x0+x1))
C=a0*x1*x2+a1*x0*x2+a2*x0*x1

end subroutine quadratic

subroutine quadratic_interpolation(x,y,xnew,ynew,N,N2)
  implicit none
  !interpolate y quadratically y->ynew 
  real(dp),intent(in)::x(:),y(:),xnew(:)
  real(dp),intent(out)::ynew(:)
  integer::N,N2,i,ii,p(1),ny
  real(dp)::dif(N),dist(1)
  real(dp)::x0,x1,x2,y0,y1,y2,A,B,C,tiny,p0,p1,p2,yold

  tiny=1.e-3

  do i=1,N2 
     dif=abs(xnew(i)-x)
     p=minloc(dif)
     dist=minval(dif)
     if (dist(1)<tiny) then
        ynew(i)=y(p(1))
        !no interpolatin needed
     else 
        !quadratic interpolation. find closest three points
        if (p(1)==1) then 
           ii=p(1)+1
        else if (p(1)==N2) then
           ii=p(1)-1
        else
           ii=p(1)
        endif
        x0=x(ii-1) 
        x1=x(ii) 
        x2=x(ii+1) 
        y0=y(ii-1) 
        y1=y(ii) 
        y2=y(ii+1) 

        call quadratic(x0,x1,x2,y0,y1,y2,A,B,C)

        ynew(i)=A*xnew(i)**2.+B*xnew(i)+C

     endif
     enddo

end subroutine quadratic_interpolation


   subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
implicit none
integer n
double precision x(n), y(n), b(n), c(n), d(n)
integer i, j, gap
double precision h

gap = n-1
! check input
if ( n < 2 ) return
if ( n < 3 ) then
  b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end if
!
! step 1: preparation
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do i = 2, gap
  d(i) = x(i+1) - x(i)
  b(i) = 2.0*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
end do
!
! step 2: end conditions 
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.0
c(n) = 0.0
if(n /= 3) then
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
end if
!
! step 3: forward elimination 
!
do i = 2, n
  h = d(i-1)/b(i-1)
  b(i) = b(i) - h*d(i-1)
  c(i) = c(i) - h*c(i-1)
end do
!
! step 4: back substitution
!
c(n) = c(n)/b(n)
do j = 1, gap
  i = n-j
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do
!
! step 5: compute spline coefficients
!
b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
do i = 1, gap
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
end do
c(n) = 3.0*c(n)
d(n) = d(n-1)
end subroutine spline

  function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
implicit none
double precision ispline
integer n
double precision  u, x(n), y(n), b(n), c(n), d(n)
integer i, j, k
double precision dx

! if u is ouside the x() interval take a boundary value (left or right)
if(u <= x(1)) then
  ispline = y(1)
  return
end if
if(u >= x(n)) then
  ispline = y(n)
  return
end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
i = 1
j = n+1
do while (j > i+1)
  k = (i+j)/2
  if(u < x(k)) then
    j=k
    else
    i=k
   end if
end do
!*
!  evaluate spline interpolation
!*
dx = u - x(i)
ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function ispline



end module interpolations


module random_tools

  use healpix_types
  use interpolations
  use random

contains



  function rms(vector)

    real(SP),intent(in):: vector(:)
    real(SP)::rms,avg
    integer(I4B)::n,i

    n=size(vector)

    avg=0.
    do i=1,N
       avg=avg+vector(i)
    enddo
    avg=avg/real(N)

    rms=0.
    do i=1,N
       rms=rms+(vector(i)-avg)**2
    enddo
    rms=rms/real(N)
    rms=sqrt(rms)
    return
  end function rms



  
  subroutine reordersample(iseed,sample)
    !reordering a vector in a random order
    implicit none
    real(sp)::sample(:)
    real(sp),allocatable::samplecopy(:),ordersample(:)
    integer::N,i,p(1),iseed


    N=size(sample)

    allocate(samplecopy(N),ordersample(N))

    samplecopy=sample

    do i=1,N
       ordersample(i)=rand()
    enddo


    do i=1,N
       P=minloc(ordersample)
       ordersample(P(1))=2. 
       sample(i)=samplecopy(P(1))
    enddo
    deallocate(samplecopy,ordersample)



  end subroutine reordersample

  subroutine poisson_constrained(iseed,x,Px,xmin_prior,xmax_prior,sample,Nsample)
    !extract poisson samples from a distribution with the constraint of the total number of objects 
    implicit none
    real(sp)::sample(:),mu,norm,xmin_prior,xmax_prior,sg
    !    real(sp)::x_m(1),x_w1(1),x_w2(1),halfmax(1),sg1(1),sg2(1)
    real(dp)::px(:),x(:),dx,xmin,xmax,peakval
    integer::i,N,j,iii,p(1),iseed,Nran,test,jjj
    integer*8::Nsample,Ngen,iostat,p_m(1)
    logical::first=.true.

    ! this does not work very well with very small samples (<10). Consider havong a miminum sample size and extracting a subsample if necessary


    sample(:)=0.
    dx=abs(x(2)-x(1))
    sg=dx/3. ! scatter for the additional samples
    N=size(x)

    norm=real(Nsample)/sum(px) 


    !print*,Nsample,norm

    !poisson sampling the original distribution. the number of objects could be larger or smaller than the number required
    iii=1
    j=1 ! position of peack distribution
    peakval=0.d0
    Ngen=0 !initialise
    do i=1,N
       mu=real(Px(i)*norm)
       Nran=random_Poisson(mu, first)! integer from poisson distribution with mean mu
       xmax=x(i)+dx/2
       xmin=x(i)-dx/2
       if (xmin < xmin_prior) xmin=xmin_prior
       if (xmax > xmax_prior) xmax=xmax_prior
       if (px(i) > peakval) then 
          j=i
          peakval=px(i)
       endif
       do j=1,Nran
          sample(iii)=rand()*(xmax-xmin)+xmin
          Ngen=iii
          if (Ngen == Nsample) goto 249 
          iii=iii+1
       enddo
    enddo

    ! if the sample is empty after the poisson run, fill the first element with value of peak of distribution
    if (Ngen==0) then 
       xmax=x(j)+dx/2
       xmin=x(j)-dx/2
       if (xmin < xmin_prior) xmin=xmin_prior
       if (xmax > xmax_prior) xmax=xmax_prior
       sample(1)=rand()*(xmax-xmin)+xmin
       Ngen=1
    endif

    !filling the remaining with replicas of the sample plus small scatter
    do iii=Ngen,Nsample
       test=nint(rand()*(ngen-1)+1)
       sample(iii)=sample(test)+random_normal()*sg
       if (sample(iii) < xmin_prior) sample(iii)=xmin_prior
       if (sample(iii) > xmax_prior) sample(iii)=xmax_prior
    enddo




249 continue 

  end subroutine poisson_constrained


  subroutine histogram(data,nbins,binsize,histo)
    implicit none
    real(dp)::histo(:,:)
    real(dp)::data(:)
    integer::ndata,i,nbins,ii,nbins_ini
    real(dp)::avg,st,binsize,test,tot,mn,mx

    ndata=size(data)
    if (nbins==-1) then
       !determine bin size
       avg=sum(data)/dble(ndata)

       test=0.
       do i=1,ndata
          test=test+(data(i)-avg)**2 ! mean is 0.
       enddo
       st=sqrt(test/dble(ndata))
       mn=minval(data)!-st !extrapolate values
       mx=maxval(data)!+st !extrapolate values

       binsize=st/3. ! three bins inside 1 sigma
       nbins=0
       test=mn
       do while (test <=mx)
          nbins=nbins+1
          test=test+binsize
       enddo
       nbins=nbins+2 ! 0 elements on first and last bin 
    else 
       if (sum(histo(:,1))==0.) then !compute binning if not provided
          mn=minval(data)
          mx=maxval(data)

          ! if number of bins is set proceed with histogram
          test=mn+binsize/2.-binsize !centre of bin, one bin more added at the beginning
          do i=1,nbins
             histo(i,1)=test
             test=test+binsize
          enddo
       endif
       do i=1,nbins
          tot=0.
          do ii=1,ndata
             if ((data(ii) >=histo(i,1)-binsize/2.) .and. (data(ii) <histo(i,1)+binsize/2.)) tot=tot+1.
          enddo
          histo(i,2)=tot
       enddo

    endif
  end subroutine histogram


 subroutine get_bins(data,nbins,binmins,binmaxs)
    ! compute bin limits to have similar number of elements in each bin
    implicit none
    real(sp)::data(:),binmaxs(:),binmins(:)
    real(dp)::binsize
    integer::i,nbins,nrows


    nrows=size(data)
  
    binmins(:)=-1.
    binmaxs(:)=-1.
    binsize=(maxval(data)-minval(data))/nbins

    binmins(1)=minval(data)
    binmaxs(1)=binmins(1)+binsize
    do i=2,nbins  
       binmins(i)=binmins(i-1)+binsize
       binmaxs(i)=binmaxs(i-1)+binsize
    enddo

  end subroutine get_bins


  subroutine assign_bins(data,flags,binmins,binmaxs,nrows,nbins,data_indices,populations)
    ! compute bin limits to have similar number of elements in each bin
    implicit none
    real(sp)::data(:)
    real(sp)::binmaxs(:),binmins(:)
    integer::data_indices(:),i,nbins,nrows,j,populations(:),flags(:)

    data_indices(:)=-100
    populations(:)=0

    do i=1,nrows
       do j=1,nbins
          if ((data(i) >= binmins(j)) .and. (data(i) < binmaxs(j)).and. (flags(i)==0)) then
             data_indices(i)=j
             populations(j)=populations(j)+1
          endif
       enddo

    enddo


  end subroutine assign_bins

!!$!this is the old version without the possibility to exclude sources with flagging
!!$  subroutine assign_bins_old(data,binmins,binmaxs,nrows,nbins,data_indices,populations)
!!$    ! compute bin limits to have similar number of elements in each bin
!!$    implicit none
!!$    real(sp)::data(:)
!!$    real(sp)::binmaxs(:),binmins(:)
!!$    integer::data_indices(:),i,nbins,nrows,j,populations(:)
!!$
!!$    data_indices(:)=100
!!$    populations(:)=0
!!$
!!$    do i=1,nrows
!!$       if (data(i) < binmaxs(1)) then 
!!$          data_indices(i)=1
!!$          populations(1)=populations(1)+1
!!$       endif
!!$       do j=2,nbins-1
!!$          if ((data(i) >= binmins(j)) .and. (data(i) < binmaxs(j))) then
!!$             data_indices(i)=j
!!$             populations(j)=populations(j)+1
!!$          endif
!!$       enddo
!!$       if (data(i) >= binmins(nbins)) then
!!$          data_indices(i)=nbins
!!$          populations(nbins)=populations(nbins)+1
!!$       endif
!!$       if (data_indices(i) ==100) then
!!$          print*,data(i),'problema'
!!$          stop
!!$       endif
!!$    enddo
!!$
!!$  end subroutine assign_bins_old
!!$

  subroutine extract_subsample(indices,nrows,i,subsample,n)
! added the extra parameter required for overlapping bins
    implicit none
    integer::indices(:),subsample(:),i,j,nrows,ii,n

    ii=1
    do j=1,nrows
       if (indices(j)==i) then 
          subsample(ii)=j
          ii=ii+1
       endif
    enddo
    n=ii-1 ! actual number of elements I filled
 end subroutine extract_subsample


 subroutine extract_subsample_new(data,nrows,bmin,bmax,subsample,N)
    implicit none
    integer::subsample(:),i,j,nrows,ii,N
    real(sp)::data(:),bmin,bmax

    subsample(:)=0
    ii=1
    do j=1,nrows
       if ((data(j) >=bmin) .and. (data(j)<bmax) .and. data(j) /=-100.) then 
          subsample(ii)=j
          ii=ii+1
       endif
    enddo
    N=ii-1

  end subroutine extract_subsample_new




  subroutine cross_new(data_big,sample_big,data_small,sample_small,matches,nbig,nsmall,scatter,iseed)
    implicit none
    integer::matches(:),nbig,nsmall
    integer::i,j,jstart,iseed,sample_big(:),sample_small(:)
    real(sp)::scatter,data_big(:),data_small(:),s



    ! scroll between the big catalogue starting from a random position. This will give more independent results if the correlation is repeated

!!$print*,iseed
!!$stop
    jstart=int(rand()*nbig)
    if (jstart <1) jstart=1

    do i=1,nsmall
       if (nbig ==0) goto 200 ! job done if all objects have been matched
       if (matches(sample_small(i)) ==-100) then ! this object has not been matched already
          s=abs(random_normal()) ! random number for the scatter
          do j=jstart,nbig

             if((data_big(sample_big(j))/=-100) .and. &
                  &(abs(data_small(sample_small(i)))-data_big(sample_big(j))<=scatter*data_small(sample_small(i))*s  ) ) then 


                matches(sample_small(i))=sample_big(j)
                data_big(sample_big(j))=-100.
                nbig=nbig-1 
                goto 100
             endif
          enddo

          if (jstart >1) then 
             do j=1,jstart-1
                if (nbig ==0) goto 200 ! job done if all objects have been matched
                if((data_big(sample_big(j)) /=-100.) .and. &
                     &(abs(data_small(sample_small(i))-data_big(sample_big(j))) <= scatter*data_small(sample_small(i))*s)) then
                   
                   matches(sample_small(i))=sample_big(j)
                   data_big(sample_big(j))=-100. 
                   nbig=nbig-1
                   goto 100
                endif
             enddo
          endif
100       continue
       endif
    enddo


200 continue
  end subroutine cross_new







  subroutine cross(data_big,data_small,matches,nbig,nsmall,scatter)
    implicit none
    integer::matches(:),nbig,nsmall
    integer::i,j,jstart,iseed
    real(sp)::scatter,data_big(:),data_small(:),s



    ! scroll between the big catalogue starting from a random position. This will give more independent results if the correlation is repeated


    jstart=int(rand()*nbig)
    if (jstart <1) jstart=1

    do i=1,nsmall
       if (nbig ==0) goto 200 ! job done if all objects have been matched
       if (matches(i) ==-100) then ! this object has not been matched already
          s=abs(random_normal()) ! random number for the scatter
           do j=jstart,nbig
             if((data_big(j) /=-100.) .and. (abs(data_small(i)-data_big(j)) <= scatter*data_small(i)*s)) then
                
                matches(i)=j
                data_big(j)=-100.
                nbig=nbig-1 
                goto 100
             endif
          enddo

          if (jstart >1) then 
             do j=1,jstart-1
                if (nbig ==0) goto 200 ! job done if all objects have been matched
                if((data_big(j) /=-100.) .and. (abs(data_small(i)-data_big(j)) <= scatter*data_small(i)*s)) then
                   
                   matches(i)=j
                   data_big(j)=-100. 
                   nbig=nbig-1
                   goto 100
                endif
             enddo
          endif
100       continue
       endif
    enddo


200 continue
  end subroutine cross

end module random_tools


module trecs_mods
! subroutine specifically writeen for T-RECS
! author A. Bonaldi
use healpix_types
use random_tools
use random
contains



subroutine pola_SFGS(inclinations,frequencies,flux,pola_flux)
! polarization for star-forming galaxies
! see Bonaldi et al. (2018) for more details

implicit none
real(sp)::inclinations(:),flux(:,:),pola_flux(:,:)
real(dp),intent(in)::frequencies(:)
!fits to points form Fig 9 Sun & Reich 2012. 
real(sp),parameter::c0=0.207937,c1=0.0326215,c2=-0.000130592,c3=-1.17845e-6 !1.4 GHz
real(sp),parameter::d0=0.795238,d1=-0.0971140,d2=0.00480952,d3=-3.83838e-5 !4.8 GHz
real(sp),parameter::e0=-0.488889,e1=0.0603199,e2=0.000585858,e3=-1.12795e-05 !2.7 GHz
real(sp),parameter::f0= 0.317460,f1=-0.0470298,f2=0.00298124,f3=-2.20539e-05 !8.4 GHz
real(sp),parameter::g0= 0.115873,g1=-0.0125156,g2=0.00110967,g3=-6.73400e-06 !22 GHz
integer::i,i14,i48,ii,j
real(sp)::pola_flux_14,pola_flux_48,pola_flux_27,pola_flux_84,pola_flux_220,y,x
real(sp)::coeff0,coeff1,coeff2,coeff3,coeff4
integer::Nsample,Nfrequencies


Nsample=size(inclinations)
Nfrequencies=size(frequencies)


do i=1,Nsample

pola_flux_14=c0+inclinations(i)*c1+inclinations(i)**2.*c2+inclinations(i)**3.*c3 !pola_fraction 1.4 GHz
pola_flux_48=d0+inclinations(i)*d1+inclinations(i)**2.*d2+inclinations(i)**3.*d3 !pola_fraction 4.8 GHz
pola_flux_27=e0+inclinations(i)*e1+inclinations(i)**2.*e2+inclinations(i)**3.*e3 !pola_fraction 2.7 GHz
pola_flux_84=f0+inclinations(i)*f1+inclinations(i)**2.*f2+inclinations(i)**3.*f3 !pola_fraction 8.4 GHz
pola_flux_220=g0+inclinations(i)*g1+inclinations(i)**2.*g2+inclinations(i)**3.*g3 !pola_fraction 22 GHz

do j=1,Nfrequencies
   x=frequencies(j)
   y= interpol(x,(/1400.,2700.,4800.,8400.,22000./),(/pola_flux_14,pola_flux_27,pola_flux_48,pola_flux_84,pola_flux_220/),5)
   if (y<0.) y=0.
   if (y>100.) y=100.
   pola_flux(j,i)=y
   enddo

pola_flux(:,i)=pola_flux(:,i)/100.*flux(:,i)


enddo



end subroutine pola_SFGS


  ! functions to compute radio luminosity from SFR 
  ! synchrotron emission
subroutine Lsynch_old(nu,sfr,L)
    implicit none
    ! synchrotron luminosity for one galaxy at several frequencies
    real(dp),intent(out)::L(:)
    real(dp),intent(in)::nu(:),sfr
    real(dp),parameter::normstar=1.68e28,norm=1.9e28,beta=3


    !Mancuso et al. 2015 eq 3
!    L=(nu/1000.)**(-0.85)*(1.+(nu/20000.)**0.5)**(-1.) !baseline Mancuso model

    !L=(nu/1000.)**(-0.85-0.07*log(nu/1000.))/1.25 !curved spectrum consistent with Mancuso at nu>1.4

    L=(nu/1000.)**(-0.85-0.07*log(nu/1000.))/1.5 !curved spectrum consistent with lower normalization wrt Mancuso - good agreement @ 150 MHz
   
    L=L*normstar/((normstar/norm/sfr)**beta+(normstar/norm/sfr))  ! erg/s/Hz

  end subroutine Lsynch_old



  ! New model from Smith et al. Linear relation between SFR L, mstar dependence
  ! New frequency dependence as a curved power law. Parameters adjusted to give agreement with observed counts from 150 MHz to 20 GHz.
  
  ! synchrotron emission
subroutine Lsynch(nu,sfr,mstar,L)
    implicit none
    ! synchrotron luminosity for one galaxy at several frequencies
    real(dp),intent(out)::L(:)
    real(dp),intent(in)::nu(:),sfr
    real(sp),intent(in)::mstar
    !real(dp),parameter::normstar=1.68e28,norm=1.9e28,beta=3
    real(sp),parameter::logL0=22.1,beta=0.85,gamma=0.402,dlogL0=-0.82 !Smith et al. 2020 eq 2 with normalization correction
    real(sp),parameter::nu0=1600.,si=-0.85,dsi=-0.1
    real(sp)::logL,L14

    !frequency dependence - Bonaldi et al. 2022
    L=(nu/nu0)**(si+dsi*log(nu/nu0)) !curved spectrum consistent with lower normalization wrt Mancuso - good agreement @ 150 MHz
        
    !Smith et al. 2020 eq (2) with normalization correction
    logL=(logL0+dlogL0 + random_normal()*4.e-3)+(beta + random_normal()*5.e-3)*log10(sfr)+(gamma+random_normal()*5.e-3)*(mstar-10)

    L=L*1.e7*10.**(logL) !1e7 is the conversion from W/Hz to erg/s/Hz

  end subroutine Lsynch




  
  ! free-free emission
  subroutine Lff(nu,sfr,L)
    ! free-free luminosity for one galaxy at several frequencies
    implicit none
    real(dp),intent(out)::L(:)
    real(dp),intent(in)::nu(:),sfr
    real(dp),parameter::Te=1.e4_dp,h=6.6260755e-34_dp,k=1.380658e-23_dp!,pi=

    L=dlog(nu/1000.*(Te/1.e4)**(-1.5))
    L=dlog(exp(5.960-0.551329*L)+2.71828)
    L=3.75e26*sfr*(Te/1.e4_dp)**(-0.5)*exp(-h*nu*1.e6/k/Te)*L

  end subroutine Lff

! dust emission (loading precomputed and selecting the closest frequency)
  subroutine Ldust(nu,tab_sed,dist,ii,L)
    ! free-free luminosity for one galaxy at several frequencies
    implicit none
    real(dp),intent(out)::L(:)
    real(dp),intent(in)::nu(:),tab_sed(:,:)
    integer,intent(in)::ii
    real(dp)::L_ir_simul,dist(:)
    integer::p(1),i
    real(sp)::norm


    norm=2.6 !early type 
    if (ii == 1) norm=3. !late type

    L_ir_simul=norm!*sfr

    do i=1,size(nu)
       dist=abs(tab_sed(:,1)-nu(i))
       p=minloc(dist)
       L(i)=L_ir_simul*tab_sed(p(1),ii+1)
    enddo

  end subroutine Ldust


  !computing frequency spectra with effective spectral indices depending on the flux and frequency
  ! this method is to enforce consistency with counts over a broad frequency range, 150 MHz to 20 GHz. 
  subroutine effective_index(i14,i48,ii,alpha,alphascaling_midf,alphascaling_highf,alphascaling_lowf,radioflux_i,frequencies)
    implicit none
    integer::p_low(1),p_mid(1),p_high(1)
    integer::i14,i48,i,iseed,jj,ii,Nfreq
    real(sp)::radioflux_i(:),alpha_mid,alpha_high,alpha_low,norm_f,logs,alphascat
    real(dp),intent(in)::alphascaling_lowf(:,:),alphascaling_highf(:,:),alphascaling_midf(:,:),frequencies(:),alpha

   
    Nfreq=size(frequencies)
    

    
    !match the source flux with those in the effective index arrays
    logs=log10(radioflux_i(i14))-3.  !1.4 GHz log flux Jy

    !flux at 1.4 GHz for flat sources
    p_mid=minloc(abs(alphascaling_midf(:,1)-logs))
    p_low=minloc(abs(alphascaling_lowf(:,1)-logs))

    alphascat=random_normal()*0.25 !scatter
    alpha_mid=alphascaling_midf(p_mid(1),ii+1)+alphascat
    alpha_low=alphascaling_lowf(p_low(1),ii+1)+alphascat
    

    do jj=1,Nfreq
       if (frequencies(jj) > 1400.) then
          Radioflux_i(jj)=(Radioflux_i(i14)*(frequencies(jj)/1400.)**alpha_mid)
       else
          !Radioflux_i(jj)=(Radioflux_i(i14)*(frequencies(jj)/1400.)**(alpha+alphascat)) !baseline
          Radioflux_i(jj)=(Radioflux_i(i14)*(frequencies(jj)/1400.)**(alpha_low)) !effective below 1.4GHz
       endif

      ! print*,frequencies(jj),radioflux_i(jj)
       
    enddo


  !now up to 4.8 GHz flux is correct
    logs=log10(radioflux_i(i48))-3.  !4.8 GHz log flux Jy
    p_high=minloc(abs(alphascaling_highf(:,1)-logs))
    alpha_high=alphascaling_highf(p_high(1),ii+1)+random_normal()*0.25



    do jj=1,Nfreq 
       if (frequencies(jj) >= 4800.) then
          Radioflux_i(jj)=(radioflux_i(i48)*(frequencies(jj)/4800.)**alpha_high)
       endif

    enddo



  end subroutine effective_index

!size model

! size of star-forming galaxies
function size_SF(logM,z,type)
! Mh DM halo mass
implicit none
real(sp),intent(in)::logM
real(dp),intent(in)::z
integer,intent(in)::type
integer::seed
real(sp)::logMstar,logR,chi,chi2,chi3,logNz,logMbz,M1,M2,alpha,omega,Mstar
real(sp)::sigma_lnR,scat
real(sp)::size_SF

!Aversa et al. 2015 Table 2 M*-Mh with evolution
real(sp),parameter::alpha0=-2.2,alpha1=-1.9,alpha2=-1.6,alpha3=4.7
real(sp),parameter::omega0=-0.75,omega1=-0.30,omega2=-1.8,omega3=2.6
real(sp),parameter::norm0=10.4,norm1=-0.8,norm2=0.8,norm3=-0.2
real(sp),parameter::mb0=11.5,mb1=-0.0,mb2=0.0,mb3=-0.0
real(sp),parameter::theta=-1.

!re-fitted parameters from Shen et al. late-type model 
!real(sp),parameter::alpha_1=0.115114,beta_1= 0.898240,gamma_1=0.198678,M0_1=3.01609e10 ! updated after calibrating from gauss fwhm to exp scale rad.
real(sp),parameter::alpha_1=0.115114,beta_1= 0.898240,gamma_1=0.1,M0_1=3.01609e10 ! updated after using physical size instead of comoving
real(sp),parameter::sigma1=0.47,sigma2=0.34,M0_sigma=4.e10 ! shen et al. table 1 for scatter


save seed

!Aversa et al. 2015 Table 2 M*-Mh with evolution
chi=log10((1.+z)/1.1)
chi2=chi**2.
chi3=chi**3.

alpha=alpha0+alpha1*chi+alpha2*chi2+alpha3*chi3
omega=omega0+omega1*chi+omega2*chi2+omega3*chi3
logMbz=mb0+mb1*chi+mb2*chi2+mb3*chi3
logNz=norm0+norm1*chi+norm2*chi2+norm3*chi3



M1=10.**(alpha*(logM-logMbz))
M2=10.**(omega*(logM-logMbz))
logMstar=logNz+theta*log10(M1+M2)


Mstar=10.**logMstar

size_SF=(gamma_1*(Mstar)**alpha_1*(1.+Mstar/M0_1)**(beta_1-alpha_1))/1000.
sigma_lnR=sigma2+(sigma1-sigma2)/(1.+(Mstar/M0_sigma)**2.)
scat=random_normal()*sigma_lnR
size_SF=log(size_SF)+scat
size_SF=exp(size_SF)

end function size_SF

! size of star-forming galaxies
function starmass(logM,z)
! Mh DM halo mass
implicit none
real(sp),intent(in)::logM
real(sp),intent(in)::z
integer::seed
real(sp)::logMstar,logR,chi,chi2,chi3,logNz,logMbz,M1,M2,alpha,omega
real(sp)::sigma_lnR,scat
real(sp)::starmass

!Aversa et al. 2015 Table 2 M*-Mh with evolution
real(sp),parameter::alpha0=-2.2,alpha1=-1.9,alpha2=-1.6,alpha3=4.7
real(sp),parameter::omega0=-0.75,omega1=-0.30,omega2=-1.8,omega3=2.6
real(sp),parameter::norm0=10.4,norm1=-0.8,norm2=0.8,norm3=-0.2
real(sp),parameter::mb0=11.5,mb1=-0.0,mb2=0.0,mb3=-0.0
real(sp),parameter::theta=-1.

!re-fitted parameters from Shen et al. late-type model 
!real(sp),parameter::alpha_1=0.115114,beta_1= 0.898240,gamma_1=0.198678,M0_1=3.01609e10 ! updated after calibrating from gauss fwhm to exp scale rad. 
!real(sp),parameter::sigma1=0.47,sigma2=0.34,M0_sigma=4.e10 ! shen et al. table 1 for scatter


save seed

!Aversa et al. 2015 Table 2 M*-Mh with evolution
chi=log10((1.+z)/1.1)
chi2=chi**2.
chi3=chi**3.

alpha=alpha0+alpha1*chi+alpha2*chi2+alpha3*chi3
omega=omega0+omega1*chi+omega2*chi2+omega3*chi3
logMbz=mb0+mb1*chi+mb2*chi2+mb3*chi3
logNz=norm0+norm1*chi+norm2*chi2+norm3*chi3



M1=10.**(alpha*(logM-logMbz))
M2=10.**(omega*(logM-logMbz))
logMstar=logNz+theta*log10(M1+M2)-2.


starmass=logMstar*(1.+random_normal()/10.)


end function starmass


end module trecs_mods

module misc
use healpix_types
contains


  function trapint(x,y)
    !integration with trapezium rule                                                                   
    implicit none
    real(dp),intent(in)::x(:),y(:)
    real(dp)::trapint,sum,h
    integer::nx,i

    nx=size(x)

    sum=0.
    do i=1,nx-1
       sum=sum+0.5*abs(x(i+1)-x(i))*(y(i)+y(i+1))
    enddo
    trapint=sum

  end function trapint
end module misc



module cosmology
  ! angular diameter distance, volume element, luminosity distance, cosmological parameters, etc.

  use healpix_types
  use misc
  implicit none
  real(dp),parameter::H0=67.0d0,omegam=0.32d0,omegav=0.68d0,c = 2.99792d8
  real(dp)::h=H0/100.
contains

  !-----------------------------------------------------------------------

  function Eh(z)
    ! E=H(z)/H0
    real(dp) :: Eh
    real(dp), intent(in) :: z
    real(dp) ::  a
    a = 1.0/(1.0+z)
    Eh=sqrt(omegam*a**(-3)+omegav)
    RETURN
  end function Eh

  !-----------------------------------------------------------------------

  function r(z)

    ! comoving coordinate distance , in units of Mpc

    real(dp), INTENT(IN) :: z
    real(dp) :: r
    real(dp)::rint(1000),zetas(1000)
    real(dp) ::dz,z0
    integer::i
    dz=z/(1000.-1.)
    z0=0.
    do i=1,999
       rint(i)=1.d0/Eh(z0)
       zetas(i)=z0
       z0=z0+dz
    enddo
    zetas(1000)=z
    rint(1000)=1.d0/Eh(z)

    r=c*trapint(zetas,rint)/1.0d5/h
    return
  end function r


  !-----------------------------------------------------------------------

  function dA(z)
    ! angular diameter distance in units of Mpc

    real(dp), INTENT(IN) :: z
    real(dp) :: dA
    dA = r(z)/(1._dp+z)
    RETURN
  end function dA
  !-----------------------------------------------------------------------

  function dVdzdO(z)
    ! volume element in units of Mpc^3
    real(dp), INTENT(IN) :: z
    real(dp) :: dVdzdO

    dVdzdO = c/1.0d5*r(z)**2/Eh(z)
  end function dVdzdO

  function lumr(z)
    ! luminosity distance in Mpc
    real(dp), INTENT(IN) :: z
    real(dp) :: lumr

    lumr=r(z)*(1._dp+z)

  end function lumr


  !-----------------------------------------------------------------------

  function theta_c(size,z)
    ! angular size of a source given comoving size and redshift
    ! size in Mpc (comoving)
    ! theta in arcsec
    real(sp),intent(in)::size
    real(dp),intent(in)::z
    real(sp)::theta_c      
    theta_c=size/r(z)*180./PI*60.*60. !arcsec
  end function theta_c

  function theta_p(size,z)
    ! angular size of a source given physical size and redshift
    ! size in Mpc (comoving)
    ! theta in arcsec
    real(sp),intent(in)::size
    real(dp),intent(in)::z
    real(sp)::theta_p      
    theta_p=size/da(z)*180./PI*60.*60. !arcsec
  end function theta_p

  
  function size_c(fov,z)
    ! comoving size of a source given angular size and redshift
    ! fov in degrees
    ! size in MPc (comoving)
    real(dp),intent(in)::fov,z
    real(dp)::size_c
    size_c=fov*PI/180.*r(z)
  end function size_c
  
  function size_p(fov,z)
    ! physical size of a source given angular and redshift
    ! fov in degrees
    ! size in MPc (physical)
    real(dp),intent(in)::fov,z
    real(dp)::size_p
    size_p=fov*PI/180.*da(z)
  end function size_p


end module cosmology

