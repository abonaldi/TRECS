
MODULE random
! A module for random number generation from the following distributions:
!
!     Distribution                    Function/subroutine name
!
!     Normal (Gaussian)               random_normal
!     Gamma                           random_gamma
!     Chi-squared                     random_chisq
!     Exponential                     random_exponential
!     Weibull                         random_Weibull
!     Beta                            random_beta
!     t                               random_t
!     Multivariate normal             random_mvnorm
!     Generalized inverse Gaussian    random_inv_gauss
!     Poisson                         random_Poisson
!     Binomial                        random_binomial1   *
!                                     random_binomial2   *
!     Negative binomial               random_neg_binomial
!     von Mises                       random_von_Mises
!     Cauchy                          random_Cauchy
!
!  Generate a random ordering of the integers 1 .. N
!                                     random_order
!     Initialize (seed) the uniform random number generator for ANY compiler
!                                     seed_random_number

!  ** Two functions are provided for the binomial distribution.
!  If the parameter values remain constant, it is recommended that the
!  first function is used (random_binomial1).   If one or both of the
!  parameters change, use the second function (random_binomial2).

! The compilers own random number generator, SUBROUTINE RANDOM_NUMBER(r),
! is used to provide a source of uniformly distributed random numbers.

! N.B. At this stage, only one random number is generated at each call to
!      one of the functions above.

! The module uses the following functions which are included here:
! bin_prob to calculate a single binomial probability
! lngamma  to calculate the logarithm to base e of the gamma function

! Some of the code is adapted from Dagpunar's book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
!
! In most of Dagpunar's routines, there is a test to see whether the value
! of one or two floating-point parameters has changed since the last call.
! These tests have been replaced by using a logical variable FIRST.
! This should be set to .TRUE. on the first call using new values of the
! parameters, and .FALSE. if the parameter values are the same as for the
! previous call.

! Version 1.10, 1 August 1996
! Changes from version 1.01
! 1. The random_order, random_Poisson & random_binomial routines have been
!    replaced with more efficient routines.
! 2. A routine, seed_random_number, has been added to seed the uniform random
!    number generator.   This requires input of the required number of seeds
!    for the particular compiler from a specified I/O unit such as a keyboard.

!     Author: Alan Miller
!             CSIRO Division of Mathematics & Statistics
!             Private Bag 10, Clayton South MDC
!             Clayton 3169, Victoria, Australia
!     Phone: (+61) 3 9545-8036      Fax: (+61) 3 9545-8080
!     e-mail: Alan.Miller @ mel.dms.csiro.au
!     WWW-page:  http://www.mel.dms.csiro/~alan

REAL, PRIVATE             :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0
PRIVATE                   :: integral


CONTAINS



REAL FUNCTION random_normal()

! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

IMPLICIT NONE

!     Local variables
REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q

!     Generate P = (u,v) uniform in rectangle enclosing acceptance region

DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156 * (v - half)

!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)

!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
END DO

!     Return ratio of P's coordinates as the normal deviate
random_normal = v/u
RETURN

END FUNCTION random_normal

!<
! random_gamma
!
! USAGE
! real random_gamma(s,first)
!
! INPUT PARAMETERS
! s  -  real shape parameter for distribution (0< s)
! first- logical whether or not this is the first call to random_gamma
!
! RETURN VALUES
! x= random number from a gamma distribution
!
! DESCRIPTION
! Generates a random number from a gamma distribution
! 
! CATEGORY
! Lib:Math:Random
!
! COMPILE LEVEL
! LOCAL
! 
!>


REAL FUNCTION random_gamma(s, first)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

!     FUNCTION GENERATES A RANDOM GAMMA VARIATE.
!     CALLS EITHER random_gamma1 (S > 1.0)
!     OR random_exponential (S = 1.0)
!     OR random_gamma2 (S < 1.0).

!     S = SHAPE PARAMETER OF DISTRIBUTION (0 < REAL).

IMPLICIT NONE
REAL, INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first

IF (s <= zero) THEN
  PRINT *, 'SHAPE PARAMETER VALUE MUST BE POSITIVE'
  STOP
END IF

IF (s > one) THEN
  random_gamma = random_gamma1(s, first)
ELSE IF (s < one) THEN
  random_gamma = random_gamma2(s, first)
ELSE
  random_gamma = random_exponential()
END IF

RETURN
END FUNCTION random_gamma



REAL FUNCTION random_gamma1(s, first)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO GAMMA**(S-1)*EXP(-GAMMA),
! BASED UPON BEST'S T DISTRIBUTION METHOD

!     S = SHAPE PARAMETER OF DISTRIBUTION
!          (1.0 < REAL)

IMPLICIT NONE

REAL, INTENT(IN)     :: s
LOGICAL, INTENT(IN)  :: first

!     Local variables
REAL    :: b, h, d, r, g, f, x, sixty4 = 64.0, three = 3.0, pt75 = 0.75
SAVE b, h

IF (s <= one) THEN
  PRINT *, 'IMPERMISSIBLE SHAPE PARAMETER VALUE'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  b = s - one
  h = SQRT(three*s - pt75)
END IF

DO
  CALL RANDOM_NUMBER(r)
  g = r - r*r
  IF (g <= zero) CYCLE
  f = (r - half)*h/SQRT(g)
  x = b + f
  IF (x <= zero) CYCLE
  CALL RANDOM_NUMBER(r)
  d = sixty4*g*(r*g)**2
  IF (d <= zero) EXIT
  IF (d*x < x - two*f*f) EXIT
  IF (LOG(d) < two*(b*LOG(x/b) - f)) EXIT
END DO
random_gamma1 = x

RETURN
END FUNCTION random_gamma1



REAL FUNCTION random_gamma2(s, first)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO
! GAMMA2**(S-1) * EXP(-GAMMA2),
! USING A SWITCHING METHOD.

!    S = SHAPE PARAMETER OF DISTRIBUTION
!          (REAL < 1.0)

IMPLICIT NONE

REAL, INTENT(IN)    :: s
LOGICAL, INTENT(IN) :: first

!     Local variables
REAL                :: a, p, c, uf, vr, d, r, x, w
real ::         vsmall, vlarge


SAVE a, p, c, uf, vr, d

vsmall = TINY(1.0); vlarge = HUGE(1.0)
IF (s <= zero .OR. s >= one) THEN
  PRINT *, 'SHAPE PARAMETER VALUE OUTSIDE PERMITTED RANGE'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  a = one - s
  p = a/(a + s*EXP(-a))
  IF (s < vsmall) THEN
    PRINT *, 'SHAPE PARAMETER VALUE TOO SMALL'
    STOP
  END IF
  c = one/s
  uf = p*(vsmall/a)**s
  vr = one - vsmall
  d = a*LOG(a)
END IF

DO
  CALL RANDOM_NUMBER(r)
  IF (r >= vr) THEN
    CYCLE
  ELSE IF (r > p) THEN
    x = a - LOG((one - r)/(one - p))
    w = a*LOG(x)-d
  ELSE IF (r > uf) THEN
    x = a*(r/p)**c
    w = x
  ELSE
    random_gamma2 = zero
    RETURN
  END IF

  CALL RANDOM_NUMBER(r)
  IF (one-r <= w .AND. r > zero) THEN
    IF (r*(w + one) >= one) CYCLE
    IF (-LOG(r) <= w) CYCLE
  END IF
  EXIT
END DO

random_gamma2 = x
RETURN

END FUNCTION random_gamma2

!<
! random_chisq
!
! USAGE
! real random_chisq(ndf,first)
!
! INPUT PARAMETERS
! ndf  -  integer number of degrees of freedom
! first- logical whether or not this is the first call to random_chisq
!
! RETURN VALUES
! x= random number, chi-sq distribution
!
! DESCRIPTION
! Generates a random number from a chi-sq distribution
! 
! CATEGORY
! Lib:Math:Random
!
! COMPILE LEVEL
! LOCAL
! 
!>



REAL FUNCTION random_chisq(ndf, first)

!     Generates a random variate from the chi-squared distribution with
!     ndf degrees of freedom

IMPLICIT NONE

INTEGER, INTENT(IN) :: ndf
LOGICAL, INTENT(IN) :: first

random_chisq = two * random_gamma(half*FLOAT(ndf), first)
RETURN

END FUNCTION random_chisq


!<
! random_exponential
!
! USAGE
! real random_exponential()
!
!
! RETURN VALUES
! x= random number, exponential distribution
!
! DESCRIPTION
! Generates a random number from an exponential distribution
! 
! CATEGORY
! Lib:Math:Random
!
! COMPILE LEVEL
! LOCAL
! 
!>


REAL FUNCTION random_exponential( )

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
! TO EXP(-random_exponential), USING INVERSION.

IMPLICIT NONE

!     Local variable
REAL  :: r

DO
  CALL RANDOM_NUMBER(r)
  IF (r > zero) EXIT
END DO

random_exponential = -LOG(r)
RETURN

END FUNCTION random_exponential


!<
! random_Weibull
!
! USAGE
! real random_Weibull(a)
!
! INPUT PARAMETERS
! a  -   real probablility distribution
!
! RETURN VALUES
! x= random number, Weibull distribution
!
! DESCRIPTION
! Generates a random number from a Weibull distribution
! 
! CATEGORY
! Lib:Math:Random
!
! COMPILE LEVEL
! LOCAL
! 
!>


REAL FUNCTION random_Weibull(a)

!     Generates a random variate from the Weibull distribution with
!     probability density:
!                      a
!               a-1  -x
!     f(x) = a.x    e

IMPLICIT NONE

REAL, INTENT(IN) :: a

!     For speed, there is no checking that a is not zero or very small.

random_Weibull = random_exponential() ** (one/a)
RETURN

END FUNCTION random_Weibull

!<
! random_beta
!
! USAGE
! real random_beta(aa,bb,first)
!
! INPUT PARAMETERS
! aa     -  real    shape parameter
! bb     -  real    shape parameter
! first  -  logical whether or not this the first time routine has been called
!
! RETURN VALUES
! x= random number, beta distribution
!
! DESCRIPTION
! Generates a random number from a beta distribution
! 
! CATEGORY
! Lib:Math:Random
!
! COMPILE LEVEL
! LOCAL
! 
!>


REAL FUNCTION random_beta(aa, bb, first)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,1]
! FROM A BETA DISTRIBUTION WITH DENSITY
! PROPORTIONAL TO BETA**(AA-1) * (1-BETA)**(BB-1).
! USING CHENG'S LOG LOGISTIC METHOD.

!     AA = SHAPE PARAMETER FROM DISTRIBUTION (0 < REAL)
!     BB = SHAPE PARAMETER FROM DISTRIBUTION (0 < REAL)

IMPLICIT NONE
REAL, INTENT(IN)    :: aa, bb
LOGICAL, INTENT(IN) :: first

!     Local variables
REAL, PARAMETER     :: aln4 = 1.3862944
REAL                :: a, b, d, f, h, t, c, g, r, s, x, y, z
LOGICAL             :: swap
real ::         vsmall, vlarge


SAVE d, f, h, t, c, swap

vsmall = TINY(1.0); vlarge = HUGE(1.0)
IF (aa <= zero .OR. bb <= zero) THEN
  PRINT *, 'IMPERMISSIBLE SHAPE PARAMETER VALUE(S)'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  a = aa
  b = bb
  swap = b > a
  IF (swap) THEN
    g = b
    b = a
    a = g
  END IF
  d = a/b
  f = a+b
  IF (b > one) THEN
    h = SQRT((two*a*b - f)/(f - two))
    t = one
  ELSE
    h = b
    t = one/(one + (a/(vlarge*b))**b)
  END IF
  c = a+h
END IF

DO
  CALL RANDOM_NUMBER(r)
  CALL RANDOM_NUMBER(x)
  s = r*r*x
  IF (r < vsmall .OR. s <= zero) CYCLE
  IF (r < t) THEN
    x = LOG(r/(one - r))/h
    y = d*EXP(x)
    z = c*x + f*LOG((one + d)/(one + y)) - aln4
    IF (s - one > z) THEN
      IF (s - s*z > one) CYCLE
      IF (LOG(s) > z) CYCLE
    END IF
    random_beta = y/(one + y)
  ELSE
    IF (4.0*s > (one + one/d)**f) CYCLE
    random_beta = one
  END IF
  EXIT
END DO

  IF (swap) random_beta = one - random_beta
  RETURN
END FUNCTION random_beta

!<
! random_t
!
! USAGE
! real random_t(m)
!
! INPUT PARAMETERS
! m  -  integer number of degrees of freedom
!
! RETURN VALUES
! x= random number, t distribution
!
! DESCRIPTION
! Generates a random number from a t distribution
! 
! CATEGORY
! Lib:Math:Random
!
! COMPILE LEVEL
! LOCAL
! 
!>


REAL FUNCTION random_t(m)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE FROM A
! T DISTRIBUTION USING KINDERMAN AND MONAHAN'S RATIO METHOD.

!     M = DEGREES OF FREEDOM OF DISTRIBUTION
!           (1 <= 1NTEGER)

IMPLICIT NONE

INTEGER, INTENT(IN) :: m

!     Local variables
REAL      :: s, c, a, f, g, r, x, v, three = 3.0, four = 4.0, quart = 0.25,   &
             five = 5.0, sixteen = 16.0
INTEGER   :: mm = 0
SAVE s, c, a, f, g

IF (m < 1) THEN
  PRINT *, 'IMPERMISSIBLE DEGREES OF FREEDOM'
  STOP
END IF

IF (m .NE. mm) THEN                    ! Initialization, if necessary
  s = m
  c = -quart*(s + one)
  a = four/(one + one/s)**c
  f = sixteen/a
  IF (m > 1) THEN
    g = s - one
    g = ((s + one)/g)**c*SQRT((s+s)/g)
  ELSE
    g = one
  END IF
  mm = m
END IF

DO
  CALL RANDOM_NUMBER(r)
  IF (r <= zero) CYCLE
  CALL RANDOM_NUMBER(v)
  x = (two*v - one)*g/r
  v = x*x
  IF (v > five - a*r) THEN
    IF (m >= 1 .AND. r*(v + three) > f) CYCLE
    IF (r > (one + v/s)**c) CYCLE
  END IF
  EXIT
END DO

random_t = x
RETURN
END FUNCTION random_t

!<
! random_mvnorm
!
! USAGE
! call random_mvnorm(n,h,d,first,x,ier)
!
! INPUT PARAMETERS
! n      -  integer number of variates in vector 
! h      -  real    J'th elemetn of vector of means
! first  -  logical whether or not this is the first time the routine has 
!                    been called
! d      -  real (I,J)th element of input variance matrix
!
! OUTPUT PARAMETERS
! x      -  real J'th element of delivered vector
! f      -  real (I,J)th element of output variance matrix
!
! DESCRIPTION
! Generates a N variate random normal vector using cholesky decomposition
! 
! CATEGORY
! Lib:Math:Random
!
! COMPILE LEVEL
! LOCAL
! 
!>


SUBROUTINE random_mvnorm(n, h, d, f, first, x, ier)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! N.B. An extra argument, ier, has been added to Dagpunar's routine

!     SUBROUTINE GENERATES AN N VARIATE RANDOM NORMAL
!     VECTOR USING A CHOLESKY DECOMPOSITION.

! ARGUMENTS:
!        N = NUMBER OF VARIATES IN VECTOR
!           (INPUT,INTEGER >= 1)
!     H(J) = J'TH ELEMENT OF VECTOR OF MEANS
!           (INPUT,REAL)
!     X(J) = J'TH ELEMENT OF DELIVERED VECTOR
!           (OUTPUT,REAL)
!
!    D(J*(J-1)/2+I) = (I,J)'TH ELEMENT OF VARIANCE MATRIX (J> = I)
!            (INPUT,REAL)
!    F((J-1)*(2*N-J)/2+I) = (I,J)'TH ELEMENT OF LOWER TRIANGULAR
!           DECOMPOSITION OF VARIANCE MATRIX (J <= I)
!            (OUTPUT,REAL)

!    FIRST = .TRUE. IF THIS IS THE FIRST CALL OF THE ROUTINE
!    OR IF THE DISTRIBUTION HAS CHANGED SINCE THE LAST CALL OF THE ROUTINE.
!    OTHERWISE SET TO .FALSE.
!            (INPUT,LOGICAL)

!    ier = 1 if the input covariance matrix is not +ve definite
!        = 0 otherwise

IMPLICIT NONE

INTEGER, INTENT(IN)  :: n
REAL, INTENT(IN)     :: h(n), d(n*(n+1)/2)
REAL, INTENT(OUT)    :: f(n*(n+1)/2), x(n)
LOGICAL, INTENT(IN)  :: first
INTEGER, INTENT(OUT) :: ier

!     Local variables
INTEGER  :: n2, j, i, m
REAL     :: y, v

SAVE n2

IF (n < 1) THEN
  PRINT *, 'SIZE OF VECTOR IS NON POSITIVE'
  STOP
END IF

ier = 0
IF (first) THEN                        ! Initialization, if necessary
  n2 = 2*n
  IF (d(1) < zero) THEN
    ier = 1
    RETURN
  END IF

  f(1) = SQRT(d(1))
  y = one/f(1)
  DO j = 2,n
    f(j) = d(1+j*(j-1)/2) * y
  END DO

  DO i = 2,n
    v = d(i*(i-1)/2+i)
    DO m = 1,i-1
      v = v - f((m-1)*(n2-m)/2+i)**2
    END DO

    IF (v < zero) THEN
      ier = 1
      RETURN
    END IF

    v = SQRT(v)
    y = one/v
    f((i-1)*(n2-i)/2+i) = v
    DO j = i+1,n
      v = d(j*(j-1)/2+i)
      DO m = 1,i-1
        v = v - f((m-1)*(n2-m)/2+i)*f((m-1)*(n2-m)/2 + j)
      END DO ! m = 1,i-1
      f((i-1)*(n2-i)/2 + j) = v*y
    END DO ! j = i+1,n
  END DO ! i = 2,n
END IF

x(1:n) = h(1:n)
DO j = 1,n
  y = random_normal()
  DO i = j,n
    x(i) = x(i) + f((j-1)*(n2-j)/2 + i) * y
  END DO ! i = j,n
END DO ! j = 1,n

END SUBROUTINE random_mvnorm



REAL FUNCTION random_inv_gauss(h, b, first)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY] FROM
! A REPARAMETERISED GENERALISED INVERSE GAUSSIAN (GIG) DISTRIBUTION
! WITH DENSITY PROPORTIONAL TO  GIG**(H-1) * EXP(-0.5*B*(GIG+1/GIG))
! USING A RATIO METHOD.

!     H = PARAMETER OF DISTRIBUTION (0 <= REAL)
!     B = PARAMETER OF DISTRIBUTION (0 < REAL)

IMPLICIT NONE
REAL, INTENT(IN)    :: h, b
LOGICAL, INTENT(IN) :: first

real ::         vsmall, vlarge

!     Local variables
REAL            :: a, c, d, e, ym, xm, r, w, r1, r2, x, quart = 0.25
SAVE a, c, d, e
vsmall = TINY(1.0); vlarge = HUGE(1.0)

IF (h < zero .OR. b <= zero) THEN
  PRINT *, 'IMPERMISSIBLE DISTRIBUTION PARAMETER VALUES'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  IF (h > quart*b*SQRT(vlarge)) THEN
    PRINT *, 'THE RATIO H:B IS TOO SMALL'
    STOP
  END IF
  e = b*b
  d = h + one
  ym = (-d + SQRT(d*d + e))/b
  IF (ym < vsmall) THEN
    PRINT *, 'THE VALUE OF B IS TOO SMALL'
    STOP
  END IF

  d = h - one
  xm = (d + SQRT(d*d + e))/b
  d = half*d
  e = -quart*b
  r = xm + one/xm
  w = xm*ym
  a = w**(-half*h) * SQRT(xm/ym) * EXP(-e*(r - ym - one/ym))
  IF (a < vsmall) THEN
    PRINT *, 'THE VALUE OF H IS TOO LARGE'
    STOP
  END IF
  c = -d*LOG(xm) - e*r
END IF

DO
  CALL RANDOM_NUMBER(r1)
  IF (r1 <= zero) CYCLE
  CALL RANDOM_NUMBER(r2)
  x = a*r2/r1
  IF (x <= zero) CYCLE
  IF (LOG(r1) < d*LOG(x) + e*(x + one/x) + c) EXIT
END DO

random_inv_gauss = x

RETURN
END FUNCTION random_inv_gauss



INTEGER FUNCTION random_Poisson(mu, first)
!**********************************************************************
!     Translated to Fortran 90 by Alan Miller from:
!                           RANLIB
!
!     Library of Fortran Routines for Random Number Generation
!
!                    Compiled and Written by:
!
!                         Barry W. Brown
!                          James Lovato
!
!             Department of Biomathematics, Box 237
!             The University of Texas, M.D. Anderson Cancer Center
!             1515 Holcombe Boulevard
!             Houston, TX      77030
!
! This work was supported by grant CA-16672 from the National Cancer Institute.

!                    GENerate POIsson random deviate

!                            Function

! Generates a single random deviate from a Poisson distribution with mean mu.

!                            Arguments

!     mu --> The mean of the Poisson distribution from which
!            a random deviate is to be generated.
!                              REAL mu

!                              Method

!     For details see:

!               Ahrens, J.H. and Dieter, U.
!               Computer Generation of Poisson Deviates
!               From Modified Normal Distributions.
!               ACM Trans. Math. Software, 8, 2
!               (June 1982),163-179

IMPLICIT NONE

!     TABLES: COEFFICIENTS A0-A7 FOR STEP F. FACTORIALS FACT
!     COEFFICIENTS A(K) - FOR PX = FK*V*V*SUM(A(K)*V**K)-DEL

!     SEPARATION OF CASES A AND B

!     .. Scalar Arguments ..
REAL, INTENT(IN)    :: mu
LOGICAL, INTENT(IN) :: first
!     ..
!     .. Local Scalars ..
REAL a0, a1, a2, a3, a4, a5, a6, a7, b1, b2, c, c0, c1, c2, c3, d, del,      &
     difmuk, e, fk, fx, fy, g, omega, p, p0, px, py, q, s, t, u, v, x, xx
INTEGER j, k, kflag, l, m
LOGICAL :: full_init
!     ..
!     .. Local Arrays ..
REAL fact(10), pp(35)
!     ..
!     .. Intrinsic Functions ..
INTRINSIC ABS, LOG, EXP, FLOAT, IFIX, MAX, MIN, SIGN, SQRT
!     ..
!     .. Data statements ..
DATA a0, a1, a2, a3, a4, a5, a6, a7/-.5, .3333333, -.2500068, .2000118, &
                                    -.1661269, .1421878, -.1384794, .1250060/
DATA fact/1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880./

SAVE pp, s, d, l, m, p, q, p0, full_init

!     ..
!     .. Executable Statements ..
IF (mu > 10.0) THEN
!     C A S E  A. (RECALCULATION OF S, D, L IF MU HAS CHANGED)

  IF (first) THEN
    s = SQRT(mu)
    d = 6.0*mu*mu

!             THE POISSON PROBABILITIES PK EXCEED THE DISCRETE NORMAL
!             PROBABILITIES FK WHENEVER K >= M(MU). L=IFIX(MU-1.1484)
!             IS AN UPPER BOUND TO M(MU) FOR ALL MU >= 10 .

    l = IFIX(mu-1.1484)
    full_init = .false.
  END IF


!     STEP N. NORMAL SAMPLE - random_normal() FOR STANDARD NORMAL DEVIATE

  g = mu + s*random_normal()
  IF (g > 0.0) THEN
    random_Poisson = IFIX(g)

!     STEP I. IMMEDIATE ACCEPTANCE IF random_Poisson IS LARGE ENOUGH

    IF (random_Poisson>=l) RETURN

!     STEP S. SQUEEZE ACCEPTANCE - SAMPLE U

    fk = FLOAT(random_Poisson)
    difmuk = mu - fk
    CALL RANDOM_NUMBER(u)
    IF (d*u >= difmuk*difmuk*difmuk) RETURN
  END IF

!     STEP P. PREPARATIONS FOR STEPS Q AND H.
!             (RECALCULATIONS OF PARAMETERS IF NECESSARY)
!             .3989423=(2*PI)**(-.5)  .416667E-1=1./24.  .1428571=1./7.
!             THE QUANTITIES B1, B2, C3, C2, C1, C0 ARE FOR THE HERMITE
!             APPROXIMATIONS TO THE DISCRETE NORMAL PROBABILITIES FK.
!             C=.1069/MU GUARANTEES MAJORIZATION BY THE 'HAT'-FUNCTION.

  IF (.NOT. full_init) THEN
    omega = .3989423/s
    b1 = .4166667E-1/mu
    b2 = .3*b1*b1
    c3 = .1428571*b1*b2
    c2 = b2 - 15.*c3
    c1 = b1 - 6.*b2 + 45.*c3
    c0 = 1. - b1 + 3.*b2 - 15.*c3
    c = .1069/mu
    full_init = .true.
  END IF

  IF (g<0.0) GO TO 50

!             'SUBROUTINE' F IS CALLED (KFLAG=0 FOR CORRECT RETURN)

  kflag = 0
  GO TO 70

!     STEP Q. QUOTIENT ACCEPTANCE (RARE CASE)

  40 IF (fy-u*fy <= py*EXP(px-fx)) RETURN

!     STEP E. EXPONENTIAL SAMPLE - random_exponential() FOR STANDARD EXPONENTIAL
!             DEVIATE E AND SAMPLE T FROM THE LAPLACE 'HAT'
!             (IF T <= -.6744 THEN PK < FK FOR ALL MU >= 10.)

  50 e = random_exponential()
  CALL RANDOM_NUMBER(u)
  u = u + u - one
  t = 1.8 + SIGN(e, u)
  IF (t <= (-.6744)) GO TO 50
  random_Poisson = IFIX(mu+s*t)
  fk = FLOAT(random_Poisson)
  difmuk = mu - fk

!             'SUBROUTINE' F IS CALLED (KFLAG=1 FOR CORRECT RETURN)

  kflag = 1
  GO TO 70

!     STEP H. HAT ACCEPTANCE (E IS REPEATED ON REJECTION)

  60 IF (c*ABS(u) > py*EXP(px+e) - fy*EXP(fx+e)) GO TO 50
  RETURN

!     STEP F. 'SUBROUTINE' F. CALCULATION OF PX, PY, FX, FY.
!             CASE random_Poisson < 10 USES FACTORIALS FROM TABLE FACT

  70 IF (random_Poisson>=10) GO TO 80
  px = -mu
  py = mu**random_Poisson/fact(random_Poisson+1)
  GO TO 110

!             CASE random_Poisson >= 10 USES POLYNOMIAL APPROXIMATION
!             A0-A7 FOR ACCURACY WHEN ADVISABLE
!             .8333333E-1=1./12.  .3989423=(2*PI)**(-.5)

  80 del = .8333333E-1/fk
  del = del - 4.8*del*del*del
  v = difmuk/fk
  IF (ABS(v)>0.25) THEN
    px = fk*LOG(one + v) - difmuk - del
  ELSE
    px = fk*v*v* (((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0) - del
  END IF
  py = .3989423/SQRT(fk)
  110 x = (half - difmuk)/s
  xx = x*x
  fx = -half*xx
  fy = omega* (((c3*xx + c2)*xx + c1)*xx + c0)
  IF (kflag) 40, 40, 60

!---------------------------------------------------------------------------
!     C A S E  B.    mu < 10
!     START NEW TABLE AND CALCULATE P0 IF NECESSARY

ELSE
  IF (first) THEN
    m = MAX(1, IFIX(mu))
    l = 0
    p = EXP(-mu)
    q = p
    p0 = p
  END IF

!     STEP U. UNIFORM SAMPLE FOR INVERSION METHOD

  DO
    CALL RANDOM_NUMBER(u)
    random_Poisson = 0
    IF (u <= p0) RETURN

!     STEP T. TABLE COMPARISON UNTIL THE END PP(L) OF THE
!             PP-TABLE OF CUMULATIVE POISSON PROBABILITIES
!             (0.458=PP(9) FOR MU=10)

    IF (l == 0) GO TO 150
    j = 1
    IF (u > 0.458) j = MIN(l, m)
    DO k = j, l
      IF (u <= pp(k)) GO TO 180
    END DO
    IF (l == 35) CYCLE

!     STEP C. CREATION OF NEW POISSON PROBABILITIES P
!             AND THEIR CUMULATIVES Q=PP(K)

    150 l = l + 1
    DO k = l, 35
      p = p*mu/FLOAT(k)
      q = q + p
      pp(k) = q
      IF (u <= q) GO TO 170
    END DO
    l = 35
  END DO

  170 l = k
  180 random_Poisson = k
  RETURN
END IF

END FUNCTION random_Poisson



INTEGER FUNCTION random_binomial1(n, p, first)

! FUNCTION GENERATES A RANDOM BINOMIAL VARIATE USING C.D.Kemp's method.
! This algorithm is suitable when many random variates are required
! with the SAME parameter values for n & p.

!    P = BERNOULLI SUCCESS PROBABILITY
!           (0 <= REAL <= 1)
!    N = NUMBER OF BERNOULLI TRIALS
!           (1 <= INTEGER)
!    FIRST = .TRUE. for the first call using the current parameter values
!          = .FALSE. if the values of (n,p) are unchanged from last call

! Reference: Kemp, C.D. (1986). `A modal method for generating binomial
!            variables', Commun. Statist. - Theor. Meth. 15(3), 805-813.

IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL, INTENT(IN)    :: p
LOGICAL, INTENT(IN) :: first

!     Local variables

INTEGER             :: r0, ru, rd
REAL                :: odds_ratio, u, p_r, pd, pu, zero = 0.0, one = 1.0

SAVE odds_ratio, r0, p_r

IF (first) THEN
  r0 = (n+1)*p
  p_r = bin_prob(n, p, r0)
  odds_ratio = p / (one - p)
END IF

CALL RANDOM_NUMBER(u)
u = u - p_r
IF (u < zero) THEN
  random_binomial1 = r0
  RETURN
END IF

pu = p_r
ru = r0
pd = p_r
rd = r0
DO
  rd = rd - 1
  IF (rd >= 0) THEN
    pd = pd * (rd+1) / (odds_ratio * FLOAT(n-rd))
    u = u - pd
    IF (u < zero) THEN
      random_binomial1 = rd
      RETURN
    END IF
  END IF

  ru = ru + 1
  IF (ru <= n) THEN
    pu = pu * (n-ru+1) * odds_ratio / FLOAT(ru)
    u = u - pu
    IF (u < zero) THEN
      random_binomial1 = ru
      RETURN
    END IF
  END IF
END DO

!     This point should not be reached, but just in case:

random_binomial1 = r0
RETURN

END FUNCTION random_binomial1



REAL FUNCTION bin_prob(n, p, r)
!     Calculate a binomial probability

IMPLICIT NONE

INTEGER, INTENT(IN) :: n, r
REAL, INTENT(IN)    :: p

!     Local variable
REAL                :: one = 1.0

bin_prob = EXP( lngamma(DBLE(n+1)) - lngamma(DBLE(r+1)) - lngamma(DBLE(n-r+1)) &
                + r*LOG(p) + (n-r)*LOG(one - p) )
RETURN

END FUNCTION bin_prob



DOUBLE PRECISION FUNCTION lngamma(x)

! Logarithm to base e of the gamma function.
!
! Accurate to about 1.e-14.
! Programmer: Alan Miller

! Latest revision of Fortran 77 version - 28 February 1988

IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: x

!       Local variables

DOUBLE PRECISION :: a1 = -4.166666666554424D-02, a2 = 2.430554511376954D-03,  &
                    a3 = -7.685928044064347D-04, a4 = 5.660478426014386D-04,  &
                    temp, arg, product, lnrt2pi = 9.189385332046727D-1,       &
                    pi = 3.141592653589793D0
LOGICAL          :: reflect

!       lngamma is not defined if x = 0 or a negative integer.

IF (x > 0.d0) GO TO 10
IF (x /= INT(x)) GO TO 10
lngamma = 0.d0
RETURN

!       If x < 0, use the reflection formula:
!               gamma(x) * gamma(1-x) = pi * cosec(pi.x)

10 reflect = (x < 0.d0)
IF (reflect) THEN
  arg = 1.d0 - x
ELSE
  arg = x
END IF

!       Increase the argument, if necessary, to make it > 10.

product = 1.d0
20 IF (arg <= 10.d0) THEN
  product = product * arg
  arg = arg + 1.d0
  GO TO 20
END IF

!  Use a polynomial approximation to Stirling's formula.
!  N.B. The real Stirling's formula is used here, not the simpler, but less
!       accurate formula given by De Moivre in a letter to Stirling, which
!       is the one usually quoted.

arg = arg - 0.5D0
temp = 1.d0/arg**2
lngamma = lnrt2pi + arg * (LOG(arg) - 1.d0 + &
                  (((a4*temp + a3)*temp + a2)*temp + a1)*temp) - LOG(product)
IF (reflect) THEN
  temp = SIN(pi * x)
  lngamma = LOG(pi/temp) - lngamma
END IF
RETURN
END FUNCTION lngamma



INTEGER FUNCTION random_binomial2(n, pp, first)
!**********************************************************************
!     Translated to Fortran 90 by Alan Miller from:
!                              RANLIB
!
!     Library of Fortran Routines for Random Number Generation
!
!                      Compiled and Written by:
!
!                           Barry W. Brown
!                            James Lovato
!
!               Department of Biomathematics, Box 237
!               The University of Texas, M.D. Anderson Cancer Center
!               1515 Holcombe Boulevard
!               Houston, TX      77030
!
! This work was supported by grant CA-16672 from the National Cancer Institute.

!                    GENerate BINomial random deviate

!                              Function

!     Generates a single random deviate from a binomial
!     distribution whose number of trials is N and whose
!     probability of an event in each trial is P.

!                              Arguments

!     N  --> The number of trials in the binomial distribution
!            from which a random deviate is to be generated.
!                              INTEGER N

!     P  --> The probability of an event in each trial of the
!            binomial distribution from which a random deviate
!            is to be generated.
!                              REAL P

!     FIRST --> Set FIRST = .TRUE. for the first call to perform initialization
!               the set FIRST = .FALSE. for further calls using the same pair
!               of parameter values (N, P).
!                              LOGICAL FIRST

!     random_binomial2 <-- A random deviate yielding the number of events
!                from N independent trials, each of which has
!                a probability of event P.
!                              INTEGER random_binomial

!                              Method

!     This is algorithm BTPE from:

!         Kachitvichyanukul, V. and Schmeiser, B. W.
!         Binomial Random Variate Generation.
!         Communications of the ACM, 31, 2 (February, 1988) 216.

!**********************************************************************

IMPLICIT NONE

!*****DETERMINE APPROPRIATE ALGORITHM AND WHETHER SETUP IS NECESSARY

!     ..
!     .. Scalar Arguments ..
REAL, INTENT(IN)    :: pp
INTEGER, INTENT(IN) :: n
LOGICAL, INTENT(IN) :: first
!     ..
!     .. Local Scalars ..
REAL    :: al, alv, amaxp, c, f, f1, f2, ffm, fm, g, p, p1, p2, p3, p4, q,  &
            qn, r, u, v, w, w2, x, x1, x2, xl, xll, xlr, xm, xnp, xnpq, xr, &
            ynorm, z, z2
REAL    :: zero = 0.0, half = 0.5, one = 1.0
INTEGER :: i, ix, ix1, k, m, mp
!     ..
!     .. Intrinsic Functions ..
INTRINSIC ABS, LOG, MIN, INT, SQRT

SAVE p, q, xnp, ffm, m, fm, xnpq, p1, xm, xl, xr, c, al, xll, xlr,  &
     p2, p3, p4, qn, r, g

!     ..
!     .. Executable Statements ..

!*****SETUP, PERFORM ONLY WHEN PARAMETERS CHANGE

IF (first) THEN
  p = MIN(pp, one-pp)
  q = one - p
  xnp = n * p
END IF

IF (xnp > 30.) THEN
  IF (first) THEN
    ffm = xnp + p
    m = ffm
    fm = m
    xnpq = xnp * q
    p1 = INT(2.195*SQRT(xnpq) - 4.6*q) + half
    xm = fm + half
    xl = xm - p1
    xr = xm + p1
    c = 0.134 + 20.5 / (15.3 + fm)
    al = (ffm-xl) / (ffm - xl*p)
    xll = al * (one + half*al)
    al = (xr - ffm) / (xr*q)
    xlr = al * (one + half*al)
    p2 = p1 * (one + c + c)
    p3 = p2 + c / xll
    p4 = p3 + c / xlr
  END IF

!*****GENERATE VARIATE, Binomial mean at least 30.

  20 CALL RANDOM_NUMBER(u)
  u = u * p4
  CALL RANDOM_NUMBER(v)

!     TRIANGULAR REGION

  IF (u <= p1) THEN
    ix = xm - p1 * v + u
    GO TO 110
  END IF

!     PARALLELOGRAM REGION

  IF (u <= p2) THEN
    x = xl + (u-p1) / c
    v = v * c + one - ABS(xm-x) / p1
    IF (v > one .OR. v <= zero) GO TO 20
    ix = x
  ELSE

!     LEFT TAIL

    IF (u <= p3) THEN
      ix = xl + LOG(v) / xll
      IF (ix < 0) GO TO 20
      v = v * (u-p2) * xll
    ELSE

!     RIGHT TAIL

      ix = xr - LOG(v) / xlr
      IF (ix > n) GO TO 20
      v = v * (u-p3) * xlr
    END IF
  END IF

!*****DETERMINE APPROPRIATE WAY TO PERFORM ACCEPT/REJECT TEST

  k = ABS(ix-m)
  IF (k <= 20 .OR. k >= xnpq/2-1) THEN

!     EXPLICIT EVALUATION

    f = one
    r = p / q
    g = (n+1) * r
    IF (m < ix) THEN
      mp = m + 1
      DO i = mp, ix
        f = f * (g/i-r)
      END DO

    ELSE IF (m > ix) THEN
      ix1 = ix + 1
      DO i = ix1, m
        f = f / (g/i-r)
      END DO
    END IF

    IF (v > f) THEN
      GO TO 20
    ELSE
      GO TO 110
    END IF
  END IF

!     SQUEEZING USING UPPER AND LOWER BOUNDS ON LOG(F(X))

  amaxp = (k/xnpq) * ((k*(k/3. + .625) + .1666666666666)/xnpq + half)
  ynorm = -k * k / (2.*xnpq)
  alv = LOG(v)
  IF (alv<ynorm - amaxp) GO TO 110
  IF (alv>ynorm + amaxp) GO TO 20

!     STIRLING'S (actually de Moivre's) FORMULA TO MACHINE ACCURACY FOR
!     THE FINAL ACCEPTANCE/REJECTION TEST

  x1 = ix + 1
  f1 = fm + one
  z = n + 1 - fm
  w = n - ix + one
  z2 = z * z
  x2 = x1 * x1
  f2 = f1 * f1
  w2 = w * w
  IF (alv - (xm*LOG(f1/x1) + (n-m+half)*LOG(z/w) + (ix-m)*LOG(w*p/(x1*q)) +    &
      (13860.-(462.-(132.-(99.-140./f2)/f2)/f2)/f2)/f1/166320. +               &
      (13860.-(462.-(132.-(99.-140./z2)/z2)/z2)/z2)/z/166320. +                &
      (13860.-(462.-(132.-(99.-140./x2)/x2)/x2)/x2)/x1/166320. +               &
      (13860.-(462.-(132.-(99.-140./w2)/w2)/w2)/w2)/w/166320.) > zero) THEN
    GO TO 20
  ELSE
    GO TO 110
  END IF

ELSE
!     INVERSE CDF LOGIC FOR MEAN LESS THAN 30
  IF (first) THEN
    qn = q ** n
    r = p / q
    g = r * (n+1)
  END IF

  90 ix = 0
  f = qn
  CALL RANDOM_NUMBER(u)
  100 IF (u >= f) THEN
    IF (ix > 110) GO TO 90
    u = u - f
    ix = ix + 1
    f = f * (g/ix - r)
    GO TO 100
  END IF
END IF

110 IF (pp > half) ix = n - ix
random_binomial2 = ix
RETURN

END FUNCTION random_binomial2




INTEGER FUNCTION random_neg_binomial(sk, p)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM NEGATIVE BINOMIAL VARIATE USING UNSTORED
! INVERSION AND/OR THE REPRODUCTIVE PROPERTY.

!    SK = NUMBER OF FAILURES REQUIRED (Dagpunar's words!)
!       = the `power' parameter of the negative binomial
!           (0 < REAL)
!    P = BERNOULLI SUCCESS PROBABILITY
!           (0 < REAL < 1)

! THE PARAMETER H IS SET SO THAT UNSTORED INVERSION ONLY IS USED WHEN P <= H,
! OTHERWISE A COMBINATION OF UNSTORED INVERSION AND
! THE REPRODUCTIVE PROPERTY IS USED.

IMPLICIT NONE
REAL, INTENT(IN)   :: sk, p

!     Local variables
! THE PARAMETER ULN = -LOG(MACHINE'S SMALLEST REAL NUMBER).

REAL, PARAMETER    :: h = 0.7
REAL               :: q, x, st, uln, v, r, s, y, g
INTEGER            :: k, i, n
real ::         vsmall, vlarge

vsmall = TINY(1.0); vlarge = HUGE(1.0)

IF (sk <= zero .OR. p <= zero .OR. p >= one) THEN
  PRINT *, 'IMPERMISSIBLE DISTRIBUTION PARAMETER VALUES'
  STOP
END IF

q = one - p
x = zero
st = sk
IF (p > h) THEN
  v = one/LOG(p)
  k = st
  DO i = 1,k
    DO
      CALL RANDOM_NUMBER(r)
      IF (r > zero) EXIT
    END DO
    n = v*LOG(r)
    x = x + n
  END DO
  st = st - k
END IF

s = zero
uln = -LOG(vsmall)
IF (st > -uln/LOG(q)) THEN
  PRINT *, ' P IS TOO LARGE FOR THIS VALUE OF SK'
  STOP
END IF

y = q**st
g = st
CALL RANDOM_NUMBER(r)
DO
  IF (y > r) EXIT
  r = r - y
  s = s + one
  y = y*p*g/s
  g = g + one
END DO

random_neg_binomial = x + s + half
RETURN
END FUNCTION random_neg_binomial



REAL FUNCTION random_von_Mises(k, first, ier)

!     Algorithm VMD from:
!     Dagpunar, J.S. (1990) `Sampling from the von Mises distribution
!     via a comparison of random numbers', J. of Appl. Statist., 17,
!     165-168.

!     Fortran 90 code by Alan Miller
!     CSIRO Division of Mathematics & Statistics

!     Arguments:
!     k (real)        parameter of the von Mises distribution.
!     first (logical) set to .TRUE. the first time that the function
!                     is called, or the first time with a new value
!                     for k.   When first = .TRUE., the function sets
!                     up starting values and may be very much slower.
!     ier (integer)   error indicator
!                     = 0 if k is acceptable
!                     = -1 if k < 0.
!                     = +1 if k > =  10.

IMPLICIT NONE
REAL, INTENT(IN)     :: k
LOGICAL, INTENT(IN)  :: first
INTEGER, INTENT(OUT) :: ier

!     Local variables

INTEGER              :: nk, j, n
REAL, PARAMETER      :: pi = 3.14159265
REAL                 :: p(20), theta(0:20), sump, r, th, lambda, rlast
DOUBLE PRECISION     :: dk

SAVE nk, p, theta

IF (first) THEN                        ! Initialization, if necessary
  ier = 0
  IF (k < zero) THEN
    ier = -1
    RETURN
  END IF

  nk = k + k + one
  IF (nk > 20) THEN
    ier = +1
    RETURN
  END IF

  dk = k
  theta(0) = zero
  IF (k > half) THEN

!     Set up array p of probabilities.

    sump = zero
    DO j = 1, nk
      IF (j < nk) THEN
        theta(j) = ACOS(one - j/k)
      ELSE
        theta(nk) = pi
      END IF

!     Numerical integration of e^[k.cos(x)] from theta(j-1) to theta(j)

      CALL integral(theta(j-1), theta(j), p(j), dk)
      sump = sump + p(j)
    END DO ! j = 1, nk
    p(1:nk) = p(1:nk) / sump
  ELSE
    p(1) = one
    theta(1) = pi
  END IF                         ! if k > 0.5
END IF                           ! if first

CALL RANDOM_NUMBER(r)
DO j = 1, nk
  r = r - p(j)
  IF (r < zero) EXIT
END DO
r = -r/p(j)

DO
  th = theta(j-1) + r*(theta(j) - theta(j-1))
  lambda = k - j + one - k*COS(th)
  n = 1
  rlast = lambda

  DO
    CALL RANDOM_NUMBER(r)
    IF (r > rlast) EXIT
    n = n + 1
    rlast = r
  END DO

  IF (n .NE. 2*(n/2)) EXIT         ! is n even?
  CALL RANDOM_NUMBER(r)
END DO

random_von_Mises = SIGN(th, (r - rlast)/(one - rlast) - half)
RETURN
END FUNCTION random_von_Mises



SUBROUTINE integral(a, b, result, dk)

!     Gaussian integration of exp(k.cosx) from a to b.

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: dk
REAL, INTENT(IN)             :: a, b
REAL, INTENT(OUT)            :: result

!     Local variables

DOUBLE PRECISION  :: xmid, range, x1, x2,                                    &
   x(3) = (/0.238619186083197D0, 0.661209386466265D0, 0.932469514203152D0/), &
   w(3) = (/0.467913934572691D0, 0.360761573048139D0, 0.171324492379170D0/)
INTEGER           :: i

xmid = (a + b)/2.d0
range = (b - a)/2.d0

result = 0.d0
DO i = 1, 3
  x1 = xmid + x(i)*range
  x2 = xmid - x(i)*range
  result = result + w(i)*(EXP(dk*COS(x1)) + EXP(dk*COS(x2)))
END DO

result = result * range
RETURN
END SUBROUTINE integral



REAL FUNCTION random_Cauchy()
IMPLICIT NONE
real ::         vsmall, vlarge


!     Generate a random deviate from the standard Cauchy distribution


!     Local variables
REAL     :: v(2)

vsmall = TINY(1.0); vlarge = HUGE(1.0)
DO
  CALL RANDOM_NUMBER(v)
  v = two*(v - half)
  IF (ABS(v(2)) < vsmall) CYCLE               ! Test for zero
  IF (v(1)**2 + v(2)**2 < one) EXIT
END DO
random_Cauchy = v(1) / v(2)

END FUNCTION random_Cauchy



SUBROUTINE random_order(order, n)

!     Generate a random ordering of the integers 1 ... n.

IMPLICIT NONE

INTEGER, INTENT(IN)  :: n
INTEGER, INTENT(OUT) :: order(n)

!     Local variables

INTEGER :: i, j, k
REAL    :: wk

DO i = 1, n
  order(i) = i
END DO

!     Starting at the end, swap the current last indicator with one
!     randomly chosen from those preceeding it.

DO i = n, 2, -1
  CALL RANDOM_NUMBER(wk)
  j = 1 + i * wk
  IF (j < i) THEN
    k = order(i)
    order(i) = order(j)
    order(j) = k
  END IF
END DO

RETURN
END SUBROUTINE random_order



SUBROUTINE seed_random_number(iounit)
IMPLICIT NONE

INTEGER, INTENT(IN)  :: iounit

! Local variables

INTEGER              :: k
INTEGER, ALLOCATABLE :: seed(:)

CALL RANDOM_SEED(SIZE=k)
ALLOCATE( seed(k) )

WRITE(*, '(1x, a, i2, a)')'Enter ', k, ' integers for random no. seeds: '
READ(*, *) seed
WRITE(iounit, '(1x, a, (7i10))') 'Random no. seeds: ', seed
CALL RANDOM_SEED(PUT=seed)

DEALLOCATE( seed )

END SUBROUTINE seed_random_number


END MODULE random

