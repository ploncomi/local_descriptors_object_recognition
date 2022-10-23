#ifndef LAPACK_FEM_SELECTED_HPP
#define LAPACK_FEM_SELECTED_HPP

#include <fem.hpp> // Fortran EMulation library of fable module
#include <sstream>
#include <string>

namespace lapack_fem {

using namespace fem::major_types;

struct common :
  fem::common
{
  fem::cmn_sve dlamc1_sve;
  fem::cmn_sve dlamc2_sve;
  fem::cmn_sve dlamch_sve;
};

// Fortran file: /net/sigma/raid1/rwgk/dist/lapack_fem/xerbla.f
inline
void
xerbla(
  str_cref srname,
  int const& info)
{
  //C
  //C  -- LAPACK auxiliary routine (preliminary version) --
  //C     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  XERBLA  is an error handler for the LAPACK routines.
  //C  It is called by an LAPACK routine if an input parameter has an
  //C  invalid value.  A message is printed and execution stops.
  //C
  //C  Installers may consider modifying the STOP statement in order to
  //C  call system-specific exception-handling facilities.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  SRNAME  (input) CHARACTER*(*)
  //C          The name of the routine which called XERBLA.
  //C
  //C  INFO    (input) INTEGER
  //C          The position of the invalid parameter in the parameter list
  //C          of the calling routine.
  //C
  //C =====================================================================
  //C
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //CFEM  WRITE( *, FMT = 9999 )SRNAME( 1:LEN_TRIM( SRNAME ) ), INFO
  //C
  std::ostringstream o;
  o << srname.elems() << ": illegal value for parameter " << info;
  throw std::runtime_error(o.str());
  //C
  //C     End of XERBLA
  //C
}

//C
//C***********************************************************************
//C
// Fortran file: /net/sigma/raid1/rwgk/dist/lapack_fem/dlamch.f
inline
double
dlamc3(
  double const& a,
  double const& b)
{
  double return_value = fem::double0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
  //C  the addition of  A  and  B ,  for use in situations where optimizers
  //C  might hold one of these in a register.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  A       (input) DOUBLE PRECISION
  //C  B       (input) DOUBLE PRECISION
  //C          The values A and B.
  //C
  //C =====================================================================
  //C
  //C     .. Executable Statements ..
  //C
  // force optimizer to not keep a and b in registers
  union double_char {
    double v;
    char s[sizeof(double)];
  };
  double_char u[2];
  u[0].v = a;
  u[1].v = b;
  for(unsigned i=0,j=1,k=0;i<sizeof(double);i++,j=1-j,k=1-k) {
    char buffer = u[k].s[i];
    u[k].s[i] = u[j].s[i];
    u[j].s[i] = buffer;
  }
  return_value = u[0].v + u[1].v;
  //C
  return return_value;
  //C
  //C     End of DLAMC3
  //C
}

struct dlamc1_save
{
  bool first;
  int lbeta;
  bool lieee1;
  bool lrnd;
  int lt;

  dlamc1_save() :
    first(fem::bool0),
    lbeta(fem::int0),
    lieee1(fem::bool0),
    lrnd(fem::bool0),
    lt(fem::int0)
  {}
};

//C
//C***********************************************************************
//C
// Fortran file: /net/sigma/raid1/rwgk/dist/lapack_fem/dlamch.f
inline
void
dlamc1(
  common& cmn,
  int& beta,
  int& t,
  bool& rnd,
  bool& ieee1)
{
  FEM_CMN_SVE(dlamc1);
  bool& first = sve.first;
  int& lbeta = sve.lbeta;
  bool& lieee1 = sve.lieee1;
  bool& lrnd = sve.lrnd;
  int& lt = sve.lt;
  if (is_called_first_time) {
    first = true;
  }
  double one = fem::double0;
  double a = fem::double0;
  double c = fem::double0;
  int rac = fem::int0;
  double b = fem::double0;
  double qtr = fem::double0;
  double savec = fem::double0;
  double f = fem::double0;
  double t1 = fem::double0;
  double t2 = fem::double0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLAMC1 determines the machine parameters given by BETA, T, RND, and
  //C  IEEE1.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  BETA    (output) INTEGER
  //C          The base of the machine.
  //C
  //C  T       (output) INTEGER
  //C          The number of ( BETA ) digits in the mantissa.
  //C
  //C  RND     (output) LOGICAL
  //C          Specifies whether proper rounding  ( RND = .TRUE. )  or
  //C          chopping  ( RND = .FALSE. )  occurs in addition. This may not
  //C          be a reliable guide to the way in which the machine performs
  //C          its arithmetic.
  //C
  //C  IEEE1   (output) LOGICAL
  //C          Specifies whether rounding appears to be done in the IEEE
  //C          'round to nearest' style.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  The routine is based on the routine  ENVRON  by Malcolm and
  //C  incorporates suggestions by Gentleman and Marovich. See
  //C
  //C     Malcolm M. A. (1972) Algorithms to reveal properties of
  //C        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
  //C
  //C     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
  //C        that reveal properties of floating point arithmetic units.
  //C        Comms. of the ACM, 17, 276-277.
  //C
  //C =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Save statement ..
  //C     ..
  //C     .. Data statements ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  if (first) {
    one = 1;
    //C
    //C        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
    //C        IEEE1, T and RND.
    //C
    //C        Throughout this routine  we use the function  DLAMC3  to ensure
    //C        that relevant values are  stored and not held in registers,  or
    //C        are not affected by optimizers.
    //C
    //C        Compute  a = 2.0**m  with the  smallest positive integer m such
    //C        that
    //C
    //C           fl( a + 1.0 ) = a.
    //C
    a = 1;
    c = 1;
    //C
    //C+       WHILE( C.EQ.ONE )LOOP
    rac = 0;
    statement_10:
    rac++;
    if (rac == 256) {
      FEM_STOP("DLAMC1 failure (10)");
    }
    if (c == one) {
      a = 2 * a;
      c = dlamc3(a, one);
      c = dlamc3(c, -a);
      goto statement_10;
    }
    //C+       END WHILE
    //C
    //C        Now compute  b = 2.0**m  with the smallest positive integer m
    //C        such that
    //C
    //C           fl( a + b ) .gt. a.
    //C
    b = 1;
    c = dlamc3(a, b);
    //C
    //C+       WHILE( C.EQ.A )LOOP
    rac = 0;
    statement_20:
    rac++;
    if (rac == 256) {
      FEM_STOP("DLAMC1 failure (20)");
    }
    if (c == a) {
      b = 2 * b;
      c = dlamc3(a, b);
      goto statement_20;
    }
    //C+       END WHILE
    //C
    //C        Now compute the base.  a and c  are neighbouring floating point
    //C        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
    //C        their difference is beta. Adding 0.25 to c is to ensure that it
    //C        is truncated to beta and not ( beta - 1 ).
    //C
    qtr = one / 4;
    savec = c;
    c = dlamc3(c, -a);
    lbeta = static_cast<int>(c + qtr);
    if (lbeta < 2) {
      lbeta = 2;
    }
    //C
    //C        Now determine whether rounding or chopping occurs,  by adding a
    //C        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
    //C
    b = lbeta;
    f = dlamc3(b / 2, -b / 100);
    c = dlamc3(f, a);
    if (c == a) {
      lrnd = true;
    }
    else {
      lrnd = false;
    }
    f = dlamc3(b / 2, b / 100);
    c = dlamc3(f, a);
    if ((lrnd) && (c == a)) {
      lrnd = false;
    }
    //C
    //C        Try and decide whether rounding is done in the  IEEE  'round to
    //C        nearest' style. B/2 is half a unit in the last place of the two
    //C        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
    //C        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
    //C        A, but adding B/2 to SAVEC should change SAVEC.
    //C
    t1 = dlamc3(b / 2, a);
    t2 = dlamc3(b / 2, savec);
    lieee1 = (t1 == a) && (t2 > savec) && lrnd;
    //C
    //C        Now find  the  mantissa, t.  It should  be the  integer part of
    //C        log to the base beta of a,  however it is safer to determine  t
    //C        by powering.  So we find t as the smallest positive integer for
    //C        which
    //C
    //C           fl( beta**t + 1.0 ) = 1.0.
    //C
    lt = 0;
    a = 1;
    c = 1;
    //C
    //C+       WHILE( C.EQ.ONE )LOOP
    rac = 0;
    statement_30:
    rac++;
    if (rac == 256) {
      FEM_STOP("DLAMC1 failure (30)");
    }
    if (c == one) {
      lt++;
      a = a * lbeta;
      c = dlamc3(a, one);
      c = dlamc3(c, -a);
      goto statement_30;
    }
    //C+       END WHILE
    //C
  }
  //C
  beta = lbeta;
  t = lt;
  rnd = lrnd;
  ieee1 = lieee1;
  first = false;
  //C
  //C     End of DLAMC1
  //C
}

//C
//C***********************************************************************
//C
// Fortran file: /net/sigma/raid1/rwgk/dist/lapack_fem/dlamch.f
inline
void
dlamc4(
  int& emin,
  double const& start,
  int const& base)
{
  double a = fem::double0;
  double one = fem::double0;
  double rbase = fem::double0;
  double zero = fem::double0;
  double b1 = fem::double0;
  double c1 = fem::double0;
  double c2 = fem::double0;
  double d1 = fem::double0;
  double d2 = fem::double0;
  int i = fem::int0;
  double b2 = fem::double0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLAMC4 is a service routine for DLAMC2.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  EMIN    (output) INTEGER
  //C          The minimum exponent before (gradual) underflow, computed by
  //C          setting A = START and dividing by BASE until the previous A
  //C          can not be recovered.
  //C
  //C  START   (input) DOUBLE PRECISION
  //C          The starting point for determining EMIN.
  //C
  //C  BASE    (input) INTEGER
  //C          The base of the machine.
  //C
  //C =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  a = start;
  one = 1;
  rbase = one / base;
  zero = 0;
  emin = 1;
  b1 = dlamc3(a * rbase, zero);
  c1 = a;
  c2 = a;
  d1 = a;
  d2 = a;
  //C+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
  //C    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
  statement_10:
  if ((c1 == a) && (c2 == a) && (d1 == a) && (d2 == a)) {
    emin = emin - 1;
    a = b1;
    b1 = dlamc3(a / base, zero);
    c1 = dlamc3(b1 * base, zero);
    d1 = zero;
    FEM_DO(i, 1, base) {
      d1 += b1;
    }
    b2 = dlamc3(a * rbase, zero);
    c2 = dlamc3(b2 / rbase, zero);
    d2 = zero;
    FEM_DO(i, 1, base) {
      d2 += b2;
    }
    goto statement_10;
  }
  //C+    END WHILE
  //C
  //C     End of DLAMC4
  //C
}

//C
//C***********************************************************************
//C
// Fortran file: /net/sigma/raid1/rwgk/dist/lapack_fem/dlamch.f
inline
void
dlamc5(
  int const& beta,
  int const& p,
  int const& emin,
  bool const& ieee,
  int& emax,
  double& rmax)
{
  int lexp = fem::int0;
  int exbits = fem::int0;
  int variable_try = fem::int0;
  int uexp = fem::int0;
  int expsum = fem::int0;
  int nbits = fem::int0;
  const double one = 1.0e0;
  double recbas = fem::double0;
  double z = fem::double0;
  const double zero = 0.0e0;
  double y = fem::double0;
  int i = fem::int0;
  double oldy = fem::double0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLAMC5 attempts to compute RMAX, the largest machine floating-point
  //C  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
  //C  approximately to a power of 2.  It will fail on machines where this
  //C  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
  //C  EMAX = 28718).  It will also fail if the value supplied for EMIN is
  //C  too large (i.e. too close to zero), probably with overflow.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  BETA    (input) INTEGER
  //C          The base of floating-point arithmetic.
  //C
  //C  P       (input) INTEGER
  //C          The number of base BETA digits in the mantissa of a
  //C          floating-point value.
  //C
  //C  EMIN    (input) INTEGER
  //C          The minimum exponent before (gradual) underflow.
  //C
  //C  IEEE    (input) LOGICAL
  //C          A logical flag specifying whether or not the arithmetic
  //C          system is thought to comply with the IEEE standard.
  //C
  //C  EMAX    (output) INTEGER
  //C          The largest exponent before overflow
  //C
  //C  RMAX    (output) DOUBLE PRECISION
  //C          The largest machine floating-point number.
  //C
  //C =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     First compute LEXP and UEXP, two powers of 2 that bound
  //C     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
  //C     approximately to the bound that is closest to abs(EMIN).
  //C     (EMAX is the exponent of the required number RMAX).
  //C
  lexp = 1;
  exbits = 1;
  statement_10:
  variable_try = lexp * 2;
  if (variable_try <= (-emin)) {
    lexp = variable_try;
    exbits++;
    goto statement_10;
  }
  if (lexp ==  - emin) {
    uexp = lexp;
  }
  else {
    uexp = variable_try;
    exbits++;
  }
  //C
  //C     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
  //C     than or equal to EMIN. EXBITS is the number of bits needed to
  //C     store the exponent.
  //C
  if ((uexp + emin) > (-lexp - emin)) {
    expsum = 2 * lexp;
  }
  else {
    expsum = 2 * uexp;
  }
  //C
  //C     EXPSUM is the exponent range, approximately equal to
  //C     EMAX - EMIN + 1 .
  //C
  emax = expsum + emin - 1;
  nbits = 1 + exbits + p;
  //C
  //C     NBITS is the total number of bits needed to store a
  //C     floating-point number.
  //C
  if ((fem::mod(nbits, 2) == 1) && (beta == 2)) {
    //C
    //C        Either there are an odd number of bits used to store a
    //C        floating-point number, which is unlikely, or some bits are
    //C        not used in the representation of numbers, which is possible,
    //C        (e.g. Cray machines) or the mantissa has an implicit bit,
    //C        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
    //C        most likely. We have to assume the last alternative.
    //C        If this is true, then we need to reduce EMAX by one because
    //C        there must be some way of representing zero in an implicit-bit
    //C        system. On machines like Cray, we are reducing EMAX by one
    //C        unnecessarily.
    //C
    emax = emax - 1;
  }
  //C
  if (ieee) {
    //C
    //C        Assume we are on an IEEE machine which reserves one exponent
    //C        for infinity and NaN.
    //C
    emax = emax - 1;
  }
  //C
  //C     Now create RMAX, the largest machine number, which should
  //C     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
  //C
  //C     First compute 1.0 - BETA**(-P), being careful that the
  //C     result is less than 1.0 .
  //C
  recbas = one / beta;
  z = beta - one;
  y = zero;
  FEM_DO(i, 1, p) {
    z = z * recbas;
    if (y < one) {
      oldy = y;
    }
    y = dlamc3(y, z);
  }
  if (y >= one) {
    y = oldy;
  }
  //C
  //C     Now multiply by BETA**EMAX to get RMAX.
  //C
  FEM_DO(i, 1, emax) {
    y = dlamc3(y * beta, zero);
  }
  //C
  rmax = y;
  //C
  //C     End of DLAMC5
  //C
}

struct dlamc2_save
{
  bool first;
  bool iwarn;
  int lbeta;
  int lemax;
  int lemin;
  double leps;
  double lrmax;
  double lrmin;
  int lt;

  dlamc2_save() :
    first(fem::bool0),
    iwarn(fem::bool0),
    lbeta(fem::int0),
    lemax(fem::int0),
    lemin(fem::int0),
    leps(fem::double0),
    lrmax(fem::double0),
    lrmin(fem::double0),
    lt(fem::int0)
  {}
};

//C
//C***********************************************************************
//C
// Fortran file: /net/sigma/raid1/rwgk/dist/lapack_fem/dlamch.f
inline
void
dlamc2(
  common& cmn,
  int& beta,
  int& t,
  bool& rnd,
  double& eps,
  int& emin,
  double& rmin,
  int& emax,
  double& rmax)
{
  FEM_CMN_SVE(dlamc2);
  bool& first = sve.first;
  bool& iwarn = sve.iwarn;
  int& lbeta = sve.lbeta;
  int& lemax = sve.lemax;
  int& lemin = sve.lemin;
  double& leps = sve.leps;
  double& lrmax = sve.lrmax;
  double& lrmin = sve.lrmin;
  int& lt = sve.lt;
  if (is_called_first_time) {
    first = true;
    iwarn = false;
  }
  double zero = fem::double0;
  double one = fem::double0;
  double two = fem::double0;
  bool lrnd = fem::bool0;
  bool lieee1 = fem::bool0;
  double b = fem::double0;
  double a = fem::double0;
  double half = fem::double0;
  double sixth = fem::double0;
  double third = fem::double0;
  double c = fem::double0;
  double rbase = fem::double0;
  double smalll = fem::double0;
  int i = fem::int0;
  int ngpmin = fem::int0;
  int ngnmin = fem::int0;
  int gpmin = fem::int0;
  int gnmin = fem::int0;
  bool ieee = fem::bool0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLAMC2 determines the machine parameters specified in its argument
  //C  list.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  BETA    (output) INTEGER
  //C          The base of the machine.
  //C
  //C  T       (output) INTEGER
  //C          The number of ( BETA ) digits in the mantissa.
  //C
  //C  RND     (output) LOGICAL
  //C          Specifies whether proper rounding  ( RND = .TRUE. )  or
  //C          chopping  ( RND = .FALSE. )  occurs in addition. This may not
  //C          be a reliable guide to the way in which the machine performs
  //C          its arithmetic.
  //C
  //C  EPS     (output) DOUBLE PRECISION
  //C          The smallest positive number such that
  //C
  //C             fl( 1.0 - EPS ) .LT. 1.0,
  //C
  //C          where fl denotes the computed value.
  //C
  //C  EMIN    (output) INTEGER
  //C          The minimum exponent before (gradual) underflow occurs.
  //C
  //C  RMIN    (output) DOUBLE PRECISION
  //C          The smallest normalized number for the machine, given by
  //C          BASE**( EMIN - 1 ), where  BASE  is the floating point value
  //C          of BETA.
  //C
  //C  EMAX    (output) INTEGER
  //C          The maximum exponent before overflow occurs.
  //C
  //C  RMAX    (output) DOUBLE PRECISION
  //C          The largest positive number for the machine, given by
  //C          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
  //C          value of BETA.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  The computation of  EPS  is based on a routine PARANOIA by
  //C  W. Kahan of the University of California at Berkeley.
  //C
  //C =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Save statement ..
  //C     ..
  //C     .. Data statements ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  if (first) {
    zero = 0;
    one = 1;
    two = 2;
    //C
    //C        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
    //C        BETA, T, RND, EPS, EMIN and RMIN.
    //C
    //C        Throughout this routine  we use the function  DLAMC3  to ensure
    //C        that relevant values are stored  and not held in registers,  or
    //C        are not affected by optimizers.
    //C
    //C        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
    //C
    dlamc1(cmn, lbeta, lt, lrnd, lieee1);
    //C
    //C        Start to find EPS.
    //C
    b = lbeta;
    a = fem::pow(b, (-lt));
    leps = a;
    //C
    //C        Try some tricks to see whether or not this is the correct  EPS.
    //C
    b = two / 3;
    half = one / 2;
    sixth = dlamc3(b, -half);
    third = dlamc3(sixth, sixth);
    b = dlamc3(third, -half);
    b = dlamc3(b, sixth);
    b = fem::abs(b);
    if (b < leps) {
      b = leps;
    }
    //C
    leps = 1;
    //C
    //C+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
    statement_10:
    if ((leps > b) && (b > zero)) {
      leps = b;
      c = dlamc3(half * leps, (fem::pow(two, 5)) * (fem::pow2(leps)));
      c = dlamc3(half, -c);
      b = dlamc3(half, c);
      c = dlamc3(half, -b);
      b = dlamc3(half, c);
      goto statement_10;
    }
    //C+       END WHILE
    //C
    if (a < leps) {
      leps = a;
    }
    //C
    //C        Computation of EPS complete.
    //C
    //C        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
    //C        Keep dividing  A by BETA until (gradual) underflow occurs. This
    //C        is detected when we cannot recover the previous A.
    //C
    rbase = one / lbeta;
    smalll = one;
    FEM_DO(i, 1, 3) {
      smalll = dlamc3(smalll * rbase, zero);
    }
    a = dlamc3(one, smalll);
    dlamc4(ngpmin, one, lbeta);
    dlamc4(ngnmin, -one, lbeta);
    dlamc4(gpmin, a, lbeta);
    dlamc4(gnmin, -a, lbeta);
    ieee = false;
    //C
    if ((ngpmin == ngnmin) && (gpmin == gnmin)) {
      if (ngpmin == gpmin) {
        lemin = ngpmin;
        //C            ( Non twos-complement machines, no gradual underflow;
        //C              e.g.,  VAX )
      }
      else if ((gpmin - ngpmin) == 3) {
        lemin = ngpmin - 1 + lt;
        ieee = true;
        //C            ( Non twos-complement machines, with gradual underflow;
        //C              e.g., IEEE standard followers )
      }
      else {
        lemin = fem::min(ngpmin, gpmin);
        //C            ( A guess; no known machine )
        iwarn = true;
      }
      //C
    }
    else if ((ngpmin == gpmin) && (ngnmin == gnmin)) {
      if (fem::abs(ngpmin - ngnmin) == 1) {
        lemin = fem::max(ngpmin, ngnmin);
        //C            ( Twos-complement machines, no gradual underflow;
        //C              e.g., CYBER 205 )
      }
      else {
        lemin = fem::min(ngpmin, ngnmin);
        //C            ( A guess; no known machine )
        iwarn = true;
      }
      //C
    }
    else if ((fem::abs(ngpmin - ngnmin) == 1) && (gpmin == gnmin)) {
      if ((gpmin - fem::min(ngpmin, ngnmin)) == 3) {
        lemin = fem::max(ngpmin, ngnmin) - 1 + lt;
        //C            ( Twos-complement machines with gradual underflow;
        //C              no known machine )
      }
      else {
        lemin = fem::min(ngpmin, ngnmin);
        //C            ( A guess; no known machine )
        iwarn = true;
      }
      //C
    }
    else {
      lemin = fem::min(ngpmin, ngnmin, gpmin, gnmin);
      //C         ( A guess; no known machine )
      iwarn = true;
    }
    first = false;
    //C**
    //C Comment out this if block if EMIN is ok
    //CFEM     IF( IWARN ) THEN
    //CFEM        FIRST = .TRUE.
    //CFEM        WRITE( 6, FMT = 9999 )LEMIN
    //CFEM     END IF
    //C**
    //C
    //C        Assume IEEE arithmetic if we found denormalised  numbers above,
    //C        or if arithmetic seems to round in the  IEEE style,  determined
    //C        in routine DLAMC1. A true IEEE machine should have both  things
    //C        true; however, faulty machines may have one or the other.
    //C
    ieee = ieee || lieee1;
    //C
    //C        Compute  RMIN by successive division by  BETA. We could compute
    //C        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
    //C        this computation.
    //C
    lrmin = 1;
    {
      int fem_do_last = 1 - lemin;
      FEM_DO(i, 1, fem_do_last) {
        lrmin = dlamc3(lrmin * rbase, zero);
      }
    }
    //C
    //C        Finally, call DLAMC5 to compute EMAX and RMAX.
    //C
    dlamc5(lbeta, lt, lemin, ieee, lemax, lrmax);
  }
  //C
  beta = lbeta;
  t = lt;
  rnd = lrnd;
  eps = leps;
  emin = lemin;
  rmin = lrmin;
  emax = lemax;
  rmax = lrmax;
  //C
  //C     End of DLAMC2
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/INSTALL/lsame.f
inline
bool
lsame(
  str_cref ca,
  str_cref cb)
{
  bool return_value = fem::bool0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  LSAME returns .TRUE. if CA is the same letter as CB regardless of
  //C  case.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  CA      (input) CHARACTER*1
  //C  CB      (input) CHARACTER*1
  //C          CA and CB specify the single characters to be compared.
  //C
  //C =====================================================================
  //C
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test if the characters are equal
  //C
  return_value = ca == cb;
  if (return_value) {
    return return_value;
  }
  //C
  //C     Now test for equivalence if both characters are alphabetic.
  //C
  int zcode = fem::ichar("Z");
  //C
  //C     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
  //C     machines, on which ICHAR returns a value with bit 8 set.
  //C     ICHAR('A') on Prime machines returns 193 which is the same as
  //C     ICHAR('A') on an EBCDIC machine.
  //C
  int inta = fem::ichar(ca);
  int intb = fem::ichar(cb);
  //C
  if (zcode == 90 || zcode == 122) {
    //C
    //C        ASCII is assumed - ZCODE is the ASCII code of either lower or
    //C        upper case 'Z'.
    //C
    if (inta >= 97 && inta <= 122) {
      inta = inta - 32;
    }
    if (intb >= 97 && intb <= 122) {
      intb = intb - 32;
    }
    //C
  }
  else if (zcode == 233 || zcode == 169) {
    //C
    //C        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
    //C        upper case 'Z'.
    //C
    if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 ||
        inta >= 162 && inta <= 169) {
      inta += 64;
    }
    if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 ||
        intb >= 162 && intb <= 169) {
      intb += 64;
    }
    //C
  }
  else if (zcode == 218 || zcode == 250) {
    //C
    //C        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
    //C        plus 128 of either lower or upper case 'Z'.
    //C
    if (inta >= 225 && inta <= 250) {
      inta = inta - 32;
    }
    if (intb >= 225 && intb <= 250) {
      intb = intb - 32;
    }
  }
  return_value = inta == intb;
  return return_value;
  //C
  //C     RETURN
  //C
  //C     End of LSAME
  //C
}

struct dlamch_save
{
  double base;
  double emax;
  double emin;
  double eps;
  bool first;
  double prec;
  double rmax;
  double rmin;
  double rnd;
  double sfmin;
  double t;

  dlamch_save() :
    base(fem::double0),
    emax(fem::double0),
    emin(fem::double0),
    eps(fem::double0),
    first(fem::bool0),
    prec(fem::double0),
    rmax(fem::double0),
    rmin(fem::double0),
    rnd(fem::double0),
    sfmin(fem::double0),
    t(fem::double0)
  {}
};

// Fortran file: /net/sigma/raid1/rwgk/dist/lapack_fem/dlamch.f
inline
double
dlamch(
  common& cmn,
  str_cref cmach)
{
  double return_value = fem::double0;
  FEM_CMN_SVE(dlamch);
  // SAVE
  double& base = sve.base;
  double& emax = sve.emax;
  double& emin = sve.emin;
  double& eps = sve.eps;
  bool& first = sve.first;
  double& prec = sve.prec;
  double& rmax = sve.rmax;
  double& rmin = sve.rmin;
  double& rnd = sve.rnd;
  double& sfmin = sve.sfmin;
  double& t = sve.t;
  //
  if (is_called_first_time) {
    first = true;
  }
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLAMCH determines double precision machine parameters.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  CMACH   (input) CHARACTER*1
  //C          Specifies the value to be returned by DLAMCH:
  //C          = 'E' or 'e',   DLAMCH := eps
  //C          = 'S' or 's ,   DLAMCH := sfmin
  //C          = 'B' or 'b',   DLAMCH := base
  //C          = 'P' or 'p',   DLAMCH := eps*base
  //C          = 'N' or 'n',   DLAMCH := t
  //C          = 'R' or 'r',   DLAMCH := rnd
  //C          = 'M' or 'm',   DLAMCH := emin
  //C          = 'U' or 'u',   DLAMCH := rmin
  //C          = 'L' or 'l',   DLAMCH := emax
  //C          = 'O' or 'o',   DLAMCH := rmax
  //C
  //C          where
  //C
  //C          eps   = relative machine precision
  //C          sfmin = safe minimum, such that 1/sfmin does not overflow
  //C          base  = base of the machine
  //C          prec  = eps*base
  //C          t     = number of (base) digits in the mantissa
  //C          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
  //C          emin  = minimum exponent before (gradual) underflow
  //C          rmin  = underflow threshold - base**(emin-1)
  //C          emax  = largest exponent before overflow
  //C          rmax  = overflow threshold  - (base**emax)*(1-eps)
  //C
  //C =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Save statement ..
  //C     ..
  //C     .. Data statements ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  int beta = fem::int0;
  int it = fem::int0;
  bool lrnd = fem::bool0;
  int imin = fem::int0;
  int imax = fem::int0;
  const double one = 1.0e+0;
  const double zero = 0.0e+0;
  double smalll = fem::double0;
  if (first) {
    dlamc2(cmn, beta, it, lrnd, eps, imin, rmin, imax, rmax);
    base = beta;
    t = it;
    if (lrnd) {
      rnd = one;
      eps = (fem::pow(base, (1 - it))) / 2;
    }
    else {
      rnd = zero;
      eps = fem::pow(base, (1 - it));
    }
    prec = eps * base;
    emin = imin;
    emax = imax;
    sfmin = rmin;
    smalll = one / rmax;
    if (smalll >= sfmin) {
      //C
      //C           Use SMALL plus a bit, to avoid the possibility of rounding
      //C           causing overflow when computing  1/sfmin.
      //C
      sfmin = smalll * (one + eps);
    }
  }
  //C
  double rmach = fem::double0;
  if (lsame(cmach, "E")) {
    rmach = eps;
  }
  else if (lsame(cmach, "S")) {
    rmach = sfmin;
  }
  else if (lsame(cmach, "B")) {
    rmach = base;
  }
  else if (lsame(cmach, "P")) {
    rmach = prec;
  }
  else if (lsame(cmach, "N")) {
    rmach = t;
  }
  else if (lsame(cmach, "R")) {
    rmach = rnd;
  }
  else if (lsame(cmach, "M")) {
    rmach = emin;
  }
  else if (lsame(cmach, "U")) {
    rmach = rmin;
  }
  else if (lsame(cmach, "L")) {
    rmach = emax;
  }
  else if (lsame(cmach, "O")) {
    rmach = rmax;
  }
  //C
  return_value = rmach;
  first = false;
  return return_value;
  //C
  //C     End of DLAMCH
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/BLAS/SRC/daxpy.f
inline
void
daxpy(
  int const& n,
  double const& da,
  arr_cref<double> dx,
  int const& incx,
  arr_ref<double> dy,
  int const& incy)
{
  dx(dimension(star));
  dy(dimension(star));
  int ix = fem::int0;
  int iy = fem::int0;
  int i = fem::int0;
  int m = fem::int0;
  int mp1 = fem::int0;
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C     DAXPY constant times a vector plus a vector.
  //C     uses unrolled loops for increments equal to one.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C     jack dongarra, linpack, 3/11/78.
  //C     modified 12/3/93, array(1) declarations changed to array(*)
  //C
  //C  =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  if (n <= 0) {
    return;
  }
  if (da == 0.0e0) {
    return;
  }
  if (incx == 1 && incy == 1) {
    goto statement_20;
  }
  //C
  //C        code for unequal increments or equal increments
  //C          not equal to 1
  //C
  ix = 1;
  iy = 1;
  if (incx < 0) {
    ix = (-n + 1) * incx + 1;
  }
  if (incy < 0) {
    iy = (-n + 1) * incy + 1;
  }
  FEM_DO(i, 1, n) {
    dy(iy) += da * dx(ix);
    ix += incx;
    iy += incy;
  }
  return;
  //C
  //C        code for both increments equal to 1
  //C
  //C        clean-up loop
  //C
  statement_20:
  m = fem::mod(n, 4);
  if (m == 0) {
    goto statement_40;
  }
  FEM_DO(i, 1, m) {
    dy(i) += da * dx(i);
  }
  if (n < 4) {
    return;
  }
  statement_40:
  mp1 = m + 1;
  FEM_DOSTEP(i, mp1, n, 4) {
    dy(i) += da * dx(i);
    dy(i + 1) += da * dx(i + 1);
    dy(i + 2) += da * dx(i + 2);
    dy(i + 3) += da * dx(i + 3);
  }
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/BLAS/SRC/dcopy.f
inline
void
dcopy(
  int const& n,
  arr_cref<double> dx,
  int const& incx,
  arr_ref<double> dy,
  int const& incy)
{
  dx(dimension(star));
  dy(dimension(star));
  int ix = fem::int0;
  int iy = fem::int0;
  int i = fem::int0;
  int m = fem::int0;
  int mp1 = fem::int0;
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C     DCOPY copies a vector, x, to a vector, y.
  //C     uses unrolled loops for increments equal to one.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C     jack dongarra, linpack, 3/11/78.
  //C     modified 12/3/93, array(1) declarations changed to array(*)
  //C
  //C  =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  if (n <= 0) {
    return;
  }
  if (incx == 1 && incy == 1) {
    goto statement_20;
  }
  //C
  //C        code for unequal increments or equal increments
  //C          not equal to 1
  //C
  ix = 1;
  iy = 1;
  if (incx < 0) {
    ix = (-n + 1) * incx + 1;
  }
  if (incy < 0) {
    iy = (-n + 1) * incy + 1;
  }
  FEM_DO(i, 1, n) {
    dy(iy) = dx(ix);
    ix += incx;
    iy += incy;
  }
  return;
  //C
  //C        code for both increments equal to 1
  //C
  //C        clean-up loop
  //C
  statement_20:
  m = fem::mod(n, 7);
  if (m == 0) {
    goto statement_40;
  }
  FEM_DO(i, 1, m) {
    dy(i) = dx(i);
  }
  if (n < 7) {
    return;
  }
  statement_40:
  mp1 = m + 1;
  FEM_DOSTEP(i, mp1, n, 7) {
    dy(i) = dx(i);
    dy(i + 1) = dx(i + 1);
    dy(i + 2) = dx(i + 2);
    dy(i + 3) = dx(i + 3);
    dy(i + 4) = dx(i + 4);
    dy(i + 5) = dx(i + 5);
    dy(i + 6) = dx(i + 6);
  }
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/BLAS/SRC/ddot.f
inline
double
ddot(
  int const& n,
  arr_cref<double> dx,
  int const& incx,
  arr_cref<double> dy,
  int const& incy)
{
  double return_value = fem::double0;
  dx(dimension(star));
  dy(dimension(star));
  double dtemp = fem::double0;
  int ix = fem::int0;
  int iy = fem::int0;
  int i = fem::int0;
  int m = fem::int0;
  int mp1 = fem::int0;
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C     DDOT forms the dot product of two vectors.
  //C     uses unrolled loops for increments equal to one.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C     jack dongarra, linpack, 3/11/78.
  //C     modified 12/3/93, array(1) declarations changed to array(*)
  //C
  //C  =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  return_value = 0.0e0;
  dtemp = 0.0e0;
  if (n <= 0) {
    return return_value;
  }
  if (incx == 1 && incy == 1) {
    goto statement_20;
  }
  //C
  //C        code for unequal increments or equal increments
  //C          not equal to 1
  //C
  ix = 1;
  iy = 1;
  if (incx < 0) {
    ix = (-n + 1) * incx + 1;
  }
  if (incy < 0) {
    iy = (-n + 1) * incy + 1;
  }
  FEM_DO(i, 1, n) {
    dtemp += dx(ix) * dy(iy);
    ix += incx;
    iy += incy;
  }
  return_value = dtemp;
  return return_value;
  //C
  //C        code for both increments equal to 1
  //C
  //C        clean-up loop
  //C
  statement_20:
  m = fem::mod(n, 5);
  if (m == 0) {
    goto statement_40;
  }
  FEM_DO(i, 1, m) {
    dtemp += dx(i) * dy(i);
  }
  if (n < 5) {
    goto statement_60;
  }
  statement_40:
  mp1 = m + 1;
  FEM_DOSTEP(i, mp1, n, 5) {
    dtemp += dx(i) * dy(i) + dx(i + 1) * dy(i + 1) + dx(i + 2) * dy(
      i + 2) + dx(i + 3) * dy(i + 3) + dx(i + 4) * dy(i + 4);
  }
  statement_60:
  return_value = dtemp;
  return return_value;
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/BLAS/SRC/dgemm.f
inline
void
dgemm(
  str_cref transa,
  str_cref transb,
  int const& m,
  int const& n,
  int const& k,
  double const& alpha,
  arr_cref<double, 2> a,
  int const& lda,
  arr_cref<double, 2> b,
  int const& ldb,
  double const& beta,
  arr_ref<double, 2> c,
  int const& ldc)
{
  a(dimension(lda, star));
  b(dimension(ldb, star));
  c(dimension(ldc, star));
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DGEMM  performs one of the matrix-matrix operations
  //C
  //C     C := alpha*op( A )*op( B ) + beta*C,
  //C
  //C  where  op( X ) is one of
  //C
  //C     op( X ) = X   or   op( X ) = X',
  //C
  //C  alpha and beta are scalars, and A, B and C are matrices, with op( A )
  //C  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  //C
  //C  Arguments
  //C  ==========
  //C
  //C  TRANSA - CHARACTER*1.
  //C           On entry, TRANSA specifies the form of op( A ) to be used in
  //C           the matrix multiplication as follows:
  //C
  //C              TRANSA = 'N' or 'n',  op( A ) = A.
  //C
  //C              TRANSA = 'T' or 't',  op( A ) = A'.
  //C
  //C              TRANSA = 'C' or 'c',  op( A ) = A'.
  //C
  //C           Unchanged on exit.
  //C
  //C  TRANSB - CHARACTER*1.
  //C           On entry, TRANSB specifies the form of op( B ) to be used in
  //C           the matrix multiplication as follows:
  //C
  //C              TRANSB = 'N' or 'n',  op( B ) = B.
  //C
  //C              TRANSB = 'T' or 't',  op( B ) = B'.
  //C
  //C              TRANSB = 'C' or 'c',  op( B ) = B'.
  //C
  //C           Unchanged on exit.
  //C
  //C  M      - INTEGER.
  //C           On entry,  M  specifies  the number  of rows  of the  matrix
  //C           op( A )  and of the  matrix  C.  M  must  be at least  zero.
  //C           Unchanged on exit.
  //C
  //C  N      - INTEGER.
  //C           On entry,  N  specifies the number  of columns of the matrix
  //C           op( B ) and the number of columns of the matrix C. N must be
  //C           at least zero.
  //C           Unchanged on exit.
  //C
  //C  K      - INTEGER.
  //C           On entry,  K  specifies  the number of columns of the matrix
  //C           op( A ) and the number of rows of the matrix op( B ). K must
  //C           be at least  zero.
  //C           Unchanged on exit.
  //C
  //C  ALPHA  - DOUBLE PRECISION.
  //C           On entry, ALPHA specifies the scalar alpha.
  //C           Unchanged on exit.
  //C
  //C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
  //C           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
  //C           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
  //C           part of the array  A  must contain the matrix  A,  otherwise
  //C           the leading  k by m  part of the array  A  must contain  the
  //C           matrix A.
  //C           Unchanged on exit.
  //C
  //C  LDA    - INTEGER.
  //C           On entry, LDA specifies the first dimension of A as declared
  //C           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
  //C           LDA must be at least  max( 1, m ), otherwise  LDA must be at
  //C           least  max( 1, k ).
  //C           Unchanged on exit.
  //C
  //C  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
  //C           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
  //C           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
  //C           part of the array  B  must contain the matrix  B,  otherwise
  //C           the leading  n by k  part of the array  B  must contain  the
  //C           matrix B.
  //C           Unchanged on exit.
  //C
  //C  LDB    - INTEGER.
  //C           On entry, LDB specifies the first dimension of B as declared
  //C           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
  //C           LDB must be at least  max( 1, k ), otherwise  LDB must be at
  //C           least  max( 1, n ).
  //C           Unchanged on exit.
  //C
  //C  BETA   - DOUBLE PRECISION.
  //C           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
  //C           supplied as zero then C need not be set on input.
  //C           Unchanged on exit.
  //C
  //C  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
  //C           Before entry, the leading  m by n  part of the array  C must
  //C           contain the matrix  C,  except when  beta  is zero, in which
  //C           case C need not be set on entry.
  //C           On exit, the array  C  is overwritten by the  m by n  matrix
  //C           ( alpha*op( A )*op( B ) + beta*C ).
  //C
  //C  LDC    - INTEGER.
  //C           On entry, LDC specifies the first dimension of C as declared
  //C           in  the  calling  (sub)  program.   LDC  must  be  at  least
  //C           max( 1, m ).
  //C           Unchanged on exit.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Level 3 Blas routine.
  //C
  //C  -- Written on 8-February-1989.
  //C     Jack Dongarra, Argonne National Laboratory.
  //C     Iain Duff, AERE Harwell.
  //C     Jeremy Du Croz, Numerical Algorithms Group Ltd.
  //C     Sven Hammarling, Numerical Algorithms Group Ltd.
  //C
  //C  =====================================================================
  //C
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Parameters ..
  //C     ..
  //C
  //C     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
  //C     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
  //C     and  columns of  A  and the  number of  rows  of  B  respectively.
  //C
  bool nota = lsame(transa, "N");
  bool notb = lsame(transb, "N");
  int nrowa = fem::int0;
  int ncola = fem::int0;
  if (nota) {
    nrowa = m;
    ncola = k;
  }
  else {
    nrowa = k;
    ncola = m;
  }
  int nrowb = fem::int0;
  if (notb) {
    nrowb = k;
  }
  else {
    nrowb = n;
  }
  //C
  //C     Test the input parameters.
  //C
  int info = 0;
  if ((!nota) && (!lsame(transa, "C")) && (!lsame(transa, "T"))) {
    info = 1;
  }
  else if ((!notb) && (!lsame(transb, "C")) && (!lsame(transb, "T"))) {
    info = 2;
  }
  else if (m < 0) {
    info = 3;
  }
  else if (n < 0) {
    info = 4;
  }
  else if (k < 0) {
    info = 5;
  }
  else if (lda < fem::max(1, nrowa)) {
    info = 8;
  }
  else if (ldb < fem::max(1, nrowb)) {
    info = 10;
  }
  else if (ldc < fem::max(1, m)) {
    info = 13;
  }
  if (info != 0) {
    xerbla("DGEMM ", info);
    return;
  }
  //C
  //C     Quick return if possible.
  //C
  const double zero = 0.0e+0;
  const double one = 1.0e+0;
  if ((m == 0) || (n == 0) || (((alpha == zero) || (k == 0)) && (
      beta == one))) {
    return;
  }
  //C
  //C     And if  alpha.eq.zero.
  //C
  int j = fem::int0;
  int i = fem::int0;
  if (alpha == zero) {
    if (beta == zero) {
      FEM_DO(j, 1, n) {
        FEM_DO(i, 1, m) {
          c(i, j) = zero;
        }
      }
    }
    else {
      FEM_DO(j, 1, n) {
        FEM_DO(i, 1, m) {
          c(i, j) = beta * c(i, j);
        }
      }
    }
    return;
  }
  //C
  //C     Start the operations.
  //C
  int l = fem::int0;
  double temp = fem::double0;
  if (notb) {
    if (nota) {
      //C
      //C           Form  C := alpha*A*B + beta*C.
      //C
      FEM_DO(j, 1, n) {
        if (beta == zero) {
          FEM_DO(i, 1, m) {
            c(i, j) = zero;
          }
        }
        else if (beta != one) {
          FEM_DO(i, 1, m) {
            c(i, j) = beta * c(i, j);
          }
        }
        FEM_DO(l, 1, k) {
          if (b(l, j) != zero) {
            temp = alpha * b(l, j);
            FEM_DO(i, 1, m) {
              c(i, j) += temp * a(i, l);
            }
          }
        }
      }
    }
    else {
      //C
      //C           Form  C := alpha*A'*B + beta*C
      //C
      FEM_DO(j, 1, n) {
        FEM_DO(i, 1, m) {
          temp = zero;
          FEM_DO(l, 1, k) {
            temp += a(l, i) * b(l, j);
          }
          if (beta == zero) {
            c(i, j) = alpha * temp;
          }
          else {
            c(i, j) = alpha * temp + beta * c(i, j);
          }
        }
      }
    }
  }
  else {
    if (nota) {
      //C
      //C           Form  C := alpha*A*B' + beta*C
      //C
      FEM_DO(j, 1, n) {
        if (beta == zero) {
          FEM_DO(i, 1, m) {
            c(i, j) = zero;
          }
        }
        else if (beta != one) {
          FEM_DO(i, 1, m) {
            c(i, j) = beta * c(i, j);
          }
        }
        FEM_DO(l, 1, k) {
          if (b(j, l) != zero) {
            temp = alpha * b(j, l);
            FEM_DO(i, 1, m) {
              c(i, j) += temp * a(i, l);
            }
          }
        }
      }
    }
    else {
      //C
      //C           Form  C := alpha*A'*B' + beta*C
      //C
      FEM_DO(j, 1, n) {
        FEM_DO(i, 1, m) {
          temp = zero;
          FEM_DO(l, 1, k) {
            temp += a(l, i) * b(j, l);
          }
          if (beta == zero) {
            c(i, j) = alpha * temp;
          }
          else {
            c(i, j) = alpha * temp + beta * c(i, j);
          }
        }
      }
    }
  }
  //C
  //C     End of DGEMM .
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/BLAS/SRC/dgemv.f
inline
void
dgemv(
  str_cref trans,
  int const& m,
  int const& n,
  double const& alpha,
  arr_cref<double, 2> a,
  int const& lda,
  arr_cref<double> x,
  int const& incx,
  double const& beta,
  arr_ref<double> y,
  int const& incy)
{
  a(dimension(lda, star));
  x(dimension(star));
  y(dimension(star));
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DGEMV  performs one of the matrix-vector operations
  //C
  //C     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
  //C
  //C  where alpha and beta are scalars, x and y are vectors and A is an
  //C  m by n matrix.
  //C
  //C  Arguments
  //C  ==========
  //C
  //C  TRANS  - CHARACTER*1.
  //C           On entry, TRANS specifies the operation to be performed as
  //C           follows:
  //C
  //C              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
  //C
  //C              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
  //C
  //C              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
  //C
  //C           Unchanged on exit.
  //C
  //C  M      - INTEGER.
  //C           On entry, M specifies the number of rows of the matrix A.
  //C           M must be at least zero.
  //C           Unchanged on exit.
  //C
  //C  N      - INTEGER.
  //C           On entry, N specifies the number of columns of the matrix A.
  //C           N must be at least zero.
  //C           Unchanged on exit.
  //C
  //C  ALPHA  - DOUBLE PRECISION.
  //C           On entry, ALPHA specifies the scalar alpha.
  //C           Unchanged on exit.
  //C
  //C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
  //C           Before entry, the leading m by n part of the array A must
  //C           contain the matrix of coefficients.
  //C           Unchanged on exit.
  //C
  //C  LDA    - INTEGER.
  //C           On entry, LDA specifies the first dimension of A as declared
  //C           in the calling (sub) program. LDA must be at least
  //C           max( 1, m ).
  //C           Unchanged on exit.
  //C
  //C  X      - DOUBLE PRECISION array of DIMENSION at least
  //C           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
  //C           and at least
  //C           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
  //C           Before entry, the incremented array X must contain the
  //C           vector x.
  //C           Unchanged on exit.
  //C
  //C  INCX   - INTEGER.
  //C           On entry, INCX specifies the increment for the elements of
  //C           X. INCX must not be zero.
  //C           Unchanged on exit.
  //C
  //C  BETA   - DOUBLE PRECISION.
  //C           On entry, BETA specifies the scalar beta. When BETA is
  //C           supplied as zero then Y need not be set on input.
  //C           Unchanged on exit.
  //C
  //C  Y      - DOUBLE PRECISION array of DIMENSION at least
  //C           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
  //C           and at least
  //C           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
  //C           Before entry with BETA non-zero, the incremented array Y
  //C           must contain the vector y. On exit, Y is overwritten by the
  //C           updated vector y.
  //C
  //C  INCY   - INTEGER.
  //C           On entry, INCY specifies the increment for the elements of
  //C           Y. INCY must not be zero.
  //C           Unchanged on exit.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Level 2 Blas routine.
  //C
  //C  -- Written on 22-October-1986.
  //C     Jack Dongarra, Argonne National Lab.
  //C     Jeremy Du Croz, Nag Central Office.
  //C     Sven Hammarling, Nag Central Office.
  //C     Richard Hanson, Sandia National Labs.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C
  //C     Test the input parameters.
  //C
  int info = 0;
  if (!lsame(trans, "N") && !lsame(trans, "T") && !lsame(trans, "C")) {
    info = 1;
  }
  else if (m < 0) {
    info = 2;
  }
  else if (n < 0) {
    info = 3;
  }
  else if (lda < fem::max(1, m)) {
    info = 6;
  }
  else if (incx == 0) {
    info = 8;
  }
  else if (incy == 0) {
    info = 11;
  }
  if (info != 0) {
    xerbla("DGEMV ", info);
    return;
  }
  //C
  //C     Quick return if possible.
  //C
  const double zero = 0.0e+0;
  const double one = 1.0e+0;
  if ((m == 0) || (n == 0) || ((alpha == zero) && (beta == one))) {
    return;
  }
  //C
  //C     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
  //C     up the start points in  X  and  Y.
  //C
  int lenx = fem::int0;
  int leny = fem::int0;
  if (lsame(trans, "N")) {
    lenx = n;
    leny = m;
  }
  else {
    lenx = m;
    leny = n;
  }
  int kx = fem::int0;
  if (incx > 0) {
    kx = 1;
  }
  else {
    kx = 1 - (lenx - 1) * incx;
  }
  int ky = fem::int0;
  if (incy > 0) {
    ky = 1;
  }
  else {
    ky = 1 - (leny - 1) * incy;
  }
  //C
  //C     Start the operations. In this version the elements of A are
  //C     accessed sequentially with one pass through A.
  //C
  //C     First form  y := beta*y.
  //C
  int i = fem::int0;
  int iy = fem::int0;
  if (beta != one) {
    if (incy == 1) {
      if (beta == zero) {
        FEM_DO(i, 1, leny) {
          y(i) = zero;
        }
      }
      else {
        FEM_DO(i, 1, leny) {
          y(i) = beta * y(i);
        }
      }
    }
    else {
      iy = ky;
      if (beta == zero) {
        FEM_DO(i, 1, leny) {
          y(iy) = zero;
          iy += incy;
        }
      }
      else {
        FEM_DO(i, 1, leny) {
          y(iy) = beta * y(iy);
          iy += incy;
        }
      }
    }
  }
  if (alpha == zero) {
    return;
  }
  int jx = fem::int0;
  int j = fem::int0;
  double temp = fem::double0;
  int jy = fem::int0;
  int ix = fem::int0;
  if (lsame(trans, "N")) {
    //C
    //C        Form  y := alpha*A*x + y.
    //C
    jx = kx;
    if (incy == 1) {
      FEM_DO(j, 1, n) {
        if (x(jx) != zero) {
          temp = alpha * x(jx);
          FEM_DO(i, 1, m) {
            y(i) += temp * a(i, j);
          }
        }
        jx += incx;
      }
    }
    else {
      FEM_DO(j, 1, n) {
        if (x(jx) != zero) {
          temp = alpha * x(jx);
          iy = ky;
          FEM_DO(i, 1, m) {
            y(iy) += temp * a(i, j);
            iy += incy;
          }
        }
        jx += incx;
      }
    }
  }
  else {
    //C
    //C        Form  y := alpha*A'*x + y.
    //C
    jy = ky;
    if (incx == 1) {
      FEM_DO(j, 1, n) {
        temp = zero;
        FEM_DO(i, 1, m) {
          temp += a(i, j) * x(i);
        }
        y(jy) += alpha * temp;
        jy += incy;
      }
    }
    else {
      FEM_DO(j, 1, n) {
        temp = zero;
        ix = kx;
        FEM_DO(i, 1, m) {
          temp += a(i, j) * x(ix);
          ix += incx;
        }
        y(jy) += alpha * temp;
        jy += incy;
      }
    }
  }
  //C
  //C     End of DGEMV .
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/BLAS/SRC/dger.f
inline
void
dger(
  int const& m,
  int const& n,
  double const& alpha,
  arr_cref<double> x,
  int const& incx,
  arr_cref<double> y,
  int const& incy,
  arr_ref<double, 2> a,
  int const& lda)
{
  x(dimension(star));
  y(dimension(star));
  a(dimension(lda, star));
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DGER   performs the rank 1 operation
  //C
  //C     A := alpha*x*y' + A,
  //C
  //C  where alpha is a scalar, x is an m element vector, y is an n element
  //C  vector and A is an m by n matrix.
  //C
  //C  Arguments
  //C  ==========
  //C
  //C  M      - INTEGER.
  //C           On entry, M specifies the number of rows of the matrix A.
  //C           M must be at least zero.
  //C           Unchanged on exit.
  //C
  //C  N      - INTEGER.
  //C           On entry, N specifies the number of columns of the matrix A.
  //C           N must be at least zero.
  //C           Unchanged on exit.
  //C
  //C  ALPHA  - DOUBLE PRECISION.
  //C           On entry, ALPHA specifies the scalar alpha.
  //C           Unchanged on exit.
  //C
  //C  X      - DOUBLE PRECISION array of dimension at least
  //C           ( 1 + ( m - 1 )*abs( INCX ) ).
  //C           Before entry, the incremented array X must contain the m
  //C           element vector x.
  //C           Unchanged on exit.
  //C
  //C  INCX   - INTEGER.
  //C           On entry, INCX specifies the increment for the elements of
  //C           X. INCX must not be zero.
  //C           Unchanged on exit.
  //C
  //C  Y      - DOUBLE PRECISION array of dimension at least
  //C           ( 1 + ( n - 1 )*abs( INCY ) ).
  //C           Before entry, the incremented array Y must contain the n
  //C           element vector y.
  //C           Unchanged on exit.
  //C
  //C  INCY   - INTEGER.
  //C           On entry, INCY specifies the increment for the elements of
  //C           Y. INCY must not be zero.
  //C           Unchanged on exit.
  //C
  //C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
  //C           Before entry, the leading m by n part of the array A must
  //C           contain the matrix of coefficients. On exit, A is
  //C           overwritten by the updated matrix.
  //C
  //C  LDA    - INTEGER.
  //C           On entry, LDA specifies the first dimension of A as declared
  //C           in the calling (sub) program. LDA must be at least
  //C           max( 1, m ).
  //C           Unchanged on exit.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Level 2 Blas routine.
  //C
  //C  -- Written on 22-October-1986.
  //C     Jack Dongarra, Argonne National Lab.
  //C     Jeremy Du Croz, Nag Central Office.
  //C     Sven Hammarling, Nag Central Office.
  //C     Richard Hanson, Sandia National Labs.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C
  //C     Test the input parameters.
  //C
  int info = 0;
  if (m < 0) {
    info = 1;
  }
  else if (n < 0) {
    info = 2;
  }
  else if (incx == 0) {
    info = 5;
  }
  else if (incy == 0) {
    info = 7;
  }
  else if (lda < fem::max(1, m)) {
    info = 9;
  }
  if (info != 0) {
    xerbla("DGER  ", info);
    return;
  }
  //C
  //C     Quick return if possible.
  //C
  const double zero = 0.0e+0;
  if ((m == 0) || (n == 0) || (alpha == zero)) {
    return;
  }
  //C
  //C     Start the operations. In this version the elements of A are
  //C     accessed sequentially with one pass through A.
  //C
  int jy = fem::int0;
  if (incy > 0) {
    jy = 1;
  }
  else {
    jy = 1 - (n - 1) * incy;
  }
  int j = fem::int0;
  double temp = fem::double0;
  int i = fem::int0;
  int kx = fem::int0;
  int ix = fem::int0;
  if (incx == 1) {
    FEM_DO(j, 1, n) {
      if (y(jy) != zero) {
        temp = alpha * y(jy);
        FEM_DO(i, 1, m) {
          a(i, j) += x(i) * temp;
        }
      }
      jy += incy;
    }
  }
  else {
    if (incx > 0) {
      kx = 1;
    }
    else {
      kx = 1 - (m - 1) * incx;
    }
    FEM_DO(j, 1, n) {
      if (y(jy) != zero) {
        temp = alpha * y(jy);
        ix = kx;
        FEM_DO(i, 1, m) {
          a(i, j) += x(ix) * temp;
          ix += incx;
        }
      }
      jy += incy;
    }
  }
  //C
  //C     End of DGER  .
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/BLAS/SRC/dnrm2.f
inline
double
dnrm2(
  int const& n,
  arr_cref<double> x,
  int const& incx)
{
  double return_value = fem::double0;
  x(dimension(star));
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DNRM2 returns the euclidean norm of a vector via the function
  //C  name, so that
  //C
  //C     DNRM2 := sqrt( x'*x )
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  -- This version written on 25-October-1982.
  //C     Modified on 14-October-1993 to inline the call to DLASSQ.
  //C     Sven Hammarling, Nag Ltd.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  const double zero = 0.0e+0;
  double norm = fem::double0;
  double scale = fem::double0;
  const double one = 1.0e+0;
  double ssq = fem::double0;
  int ix = fem::int0;
  double absxi = fem::double0;
  if (n < 1 || incx < 1) {
    norm = zero;
  }
  else if (n == 1) {
    norm = fem::abs(x(1));
  }
  else {
    scale = zero;
    ssq = one;
    //C        The following loop is equivalent to this call to the LAPACK
    //C        auxiliary routine:
    //C        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
    //C
    FEM_DOSTEP(ix, 1, 1 + (n - 1) * incx, incx) {
      if (x(ix) != zero) {
        absxi = fem::abs(x(ix));
        if (scale < absxi) {
          ssq = one + ssq * fem::pow2((scale / absxi));
          scale = absxi;
        }
        else {
          ssq += fem::pow2((absxi / scale));
        }
      }
    }
    norm = scale * fem::sqrt(ssq);
  }
  //C
  return_value = norm;
  return return_value;
  //C
  //C     End of DNRM2.
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/BLAS/SRC/drot.f
inline
void
drot(
  int const& n,
  arr_ref<double> dx,
  int const& incx,
  arr_ref<double> dy,
  int const& incy,
  double const& c,
  double const& s)
{
  dx(dimension(star));
  dy(dimension(star));
  int ix = fem::int0;
  int iy = fem::int0;
  int i = fem::int0;
  double dtemp = fem::double0;
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C     DROT applies a plane rotation.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C     jack dongarra, linpack, 3/11/78.
  //C     modified 12/3/93, array(1) declarations changed to array(*)
  //C
  //C  =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  if (n <= 0) {
    return;
  }
  if (incx == 1 && incy == 1) {
    goto statement_20;
  }
  //C
  //C       code for unequal increments or equal increments not equal
  //C         to 1
  //C
  ix = 1;
  iy = 1;
  if (incx < 0) {
    ix = (-n + 1) * incx + 1;
  }
  if (incy < 0) {
    iy = (-n + 1) * incy + 1;
  }
  FEM_DO(i, 1, n) {
    dtemp = c * dx(ix) + s * dy(iy);
    dy(iy) = c * dy(iy) - s * dx(ix);
    dx(ix) = dtemp;
    ix += incx;
    iy += incy;
  }
  return;
  //C
  //C       code for both increments equal to 1
  //C
  statement_20:
  FEM_DO(i, 1, n) {
    dtemp = c * dx(i) + s * dy(i);
    dy(i) = c * dy(i) - s * dx(i);
    dx(i) = dtemp;
  }
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/BLAS/SRC/dscal.f
inline
void
dscal(
  int const& n,
  double const& da,
  arr_ref<double> dx,
  int const& incx)
{
  dx(dimension(star));
  int nincx = fem::int0;
  int i = fem::int0;
  int m = fem::int0;
  int mp1 = fem::int0;
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C     DSCAL scales a vector by a constant.
  //C     uses unrolled loops for increment equal to one.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C     jack dongarra, linpack, 3/11/78.
  //C     modified 3/93 to return if incx .le. 0.
  //C     modified 12/3/93, array(1) declarations changed to array(*)
  //C
  //C  =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  if (n <= 0 || incx <= 0) {
    return;
  }
  if (incx == 1) {
    goto statement_20;
  }
  //C
  //C        code for increment not equal to 1
  //C
  nincx = n * incx;
  FEM_DOSTEP(i, 1, nincx, incx) {
    dx(i) = da * dx(i);
  }
  return;
  //C
  //C        code for increment equal to 1
  //C
  //C        clean-up loop
  //C
  statement_20:
  m = fem::mod(n, 5);
  if (m == 0) {
    goto statement_40;
  }
  FEM_DO(i, 1, m) {
    dx(i) = da * dx(i);
  }
  if (n < 5) {
    return;
  }
  statement_40:
  mp1 = m + 1;
  FEM_DOSTEP(i, mp1, n, 5) {
    dx(i) = da * dx(i);
    dx(i + 1) = da * dx(i + 1);
    dx(i + 2) = da * dx(i + 2);
    dx(i + 3) = da * dx(i + 3);
    dx(i + 4) = da * dx(i + 4);
  }
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/BLAS/SRC/dswap.f
inline
void
dswap(
  int const& n,
  arr_ref<double> dx,
  int const& incx,
  arr_ref<double> dy,
  int const& incy)
{
  dx(dimension(star));
  dy(dimension(star));
  int ix = fem::int0;
  int iy = fem::int0;
  int i = fem::int0;
  double dtemp = fem::double0;
  int m = fem::int0;
  int mp1 = fem::int0;
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C     interchanges two vectors.
  //C     uses unrolled loops for increments equal one.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C     jack dongarra, linpack, 3/11/78.
  //C     modified 12/3/93, array(1) declarations changed to array(*)
  //C
  //C  =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  if (n <= 0) {
    return;
  }
  if (incx == 1 && incy == 1) {
    goto statement_20;
  }
  //C
  //C       code for unequal increments or equal increments not equal
  //C         to 1
  //C
  ix = 1;
  iy = 1;
  if (incx < 0) {
    ix = (-n + 1) * incx + 1;
  }
  if (incy < 0) {
    iy = (-n + 1) * incy + 1;
  }
  FEM_DO(i, 1, n) {
    dtemp = dx(ix);
    dx(ix) = dy(iy);
    dy(iy) = dtemp;
    ix += incx;
    iy += incy;
  }
  return;
  //C
  //C       code for both increments equal to 1
  //C
  //C       clean-up loop
  //C
  statement_20:
  m = fem::mod(n, 3);
  if (m == 0) {
    goto statement_40;
  }
  FEM_DO(i, 1, m) {
    dtemp = dx(i);
    dx(i) = dy(i);
    dy(i) = dtemp;
  }
  if (n < 3) {
    return;
  }
  statement_40:
  mp1 = m + 1;
  FEM_DOSTEP(i, mp1, n, 3) {
    dtemp = dx(i);
    dx(i) = dy(i);
    dy(i) = dtemp;
    dtemp = dx(i + 1);
    dx(i + 1) = dy(i + 1);
    dy(i + 1) = dtemp;
    dtemp = dx(i + 2);
    dx(i + 2) = dy(i + 2);
    dy(i + 2) = dtemp;
  }
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/BLAS/SRC/dsymv.f
inline
void
dsymv(
  str_cref uplo,
  int const& n,
  double const& alpha,
  arr_cref<double, 2> a,
  int const& lda,
  arr_cref<double> x,
  int const& incx,
  double const& beta,
  arr_ref<double> y,
  int const& incy)
{
  a(dimension(lda, star));
  x(dimension(star));
  y(dimension(star));
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DSYMV  performs the matrix-vector  operation
  //C
  //C     y := alpha*A*x + beta*y,
  //C
  //C  where alpha and beta are scalars, x and y are n element vectors and
  //C  A is an n by n symmetric matrix.
  //C
  //C  Arguments
  //C  ==========
  //C
  //C  UPLO   - CHARACTER*1.
  //C           On entry, UPLO specifies whether the upper or lower
  //C           triangular part of the array A is to be referenced as
  //C           follows:
  //C
  //C              UPLO = 'U' or 'u'   Only the upper triangular part of A
  //C                                  is to be referenced.
  //C
  //C              UPLO = 'L' or 'l'   Only the lower triangular part of A
  //C                                  is to be referenced.
  //C
  //C           Unchanged on exit.
  //C
  //C  N      - INTEGER.
  //C           On entry, N specifies the order of the matrix A.
  //C           N must be at least zero.
  //C           Unchanged on exit.
  //C
  //C  ALPHA  - DOUBLE PRECISION.
  //C           On entry, ALPHA specifies the scalar alpha.
  //C           Unchanged on exit.
  //C
  //C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
  //C           Before entry with  UPLO = 'U' or 'u', the leading n by n
  //C           upper triangular part of the array A must contain the upper
  //C           triangular part of the symmetric matrix and the strictly
  //C           lower triangular part of A is not referenced.
  //C           Before entry with UPLO = 'L' or 'l', the leading n by n
  //C           lower triangular part of the array A must contain the lower
  //C           triangular part of the symmetric matrix and the strictly
  //C           upper triangular part of A is not referenced.
  //C           Unchanged on exit.
  //C
  //C  LDA    - INTEGER.
  //C           On entry, LDA specifies the first dimension of A as declared
  //C           in the calling (sub) program. LDA must be at least
  //C           max( 1, n ).
  //C           Unchanged on exit.
  //C
  //C  X      - DOUBLE PRECISION array of dimension at least
  //C           ( 1 + ( n - 1 )*abs( INCX ) ).
  //C           Before entry, the incremented array X must contain the n
  //C           element vector x.
  //C           Unchanged on exit.
  //C
  //C  INCX   - INTEGER.
  //C           On entry, INCX specifies the increment for the elements of
  //C           X. INCX must not be zero.
  //C           Unchanged on exit.
  //C
  //C  BETA   - DOUBLE PRECISION.
  //C           On entry, BETA specifies the scalar beta. When BETA is
  //C           supplied as zero then Y need not be set on input.
  //C           Unchanged on exit.
  //C
  //C  Y      - DOUBLE PRECISION array of dimension at least
  //C           ( 1 + ( n - 1 )*abs( INCY ) ).
  //C           Before entry, the incremented array Y must contain the n
  //C           element vector y. On exit, Y is overwritten by the updated
  //C           vector y.
  //C
  //C  INCY   - INTEGER.
  //C           On entry, INCY specifies the increment for the elements of
  //C           Y. INCY must not be zero.
  //C           Unchanged on exit.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Level 2 Blas routine.
  //C
  //C  -- Written on 22-October-1986.
  //C     Jack Dongarra, Argonne National Lab.
  //C     Jeremy Du Croz, Nag Central Office.
  //C     Sven Hammarling, Nag Central Office.
  //C     Richard Hanson, Sandia National Labs.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C
  //C     Test the input parameters.
  //C
  int info = 0;
  if (!lsame(uplo, "U") && !lsame(uplo, "L")) {
    info = 1;
  }
  else if (n < 0) {
    info = 2;
  }
  else if (lda < fem::max(1, n)) {
    info = 5;
  }
  else if (incx == 0) {
    info = 7;
  }
  else if (incy == 0) {
    info = 10;
  }
  if (info != 0) {
    xerbla("DSYMV ", info);
    return;
  }
  //C
  //C     Quick return if possible.
  //C
  const double zero = 0.0e+0;
  const double one = 1.0e+0;
  if ((n == 0) || ((alpha == zero) && (beta == one))) {
    return;
  }
  //C
  //C     Set up the start points in  X  and  Y.
  //C
  int kx = fem::int0;
  if (incx > 0) {
    kx = 1;
  }
  else {
    kx = 1 - (n - 1) * incx;
  }
  int ky = fem::int0;
  if (incy > 0) {
    ky = 1;
  }
  else {
    ky = 1 - (n - 1) * incy;
  }
  //C
  //C     Start the operations. In this version the elements of A are
  //C     accessed sequentially with one pass through the triangular part
  //C     of A.
  //C
  //C     First form  y := beta*y.
  //C
  int i = fem::int0;
  int iy = fem::int0;
  if (beta != one) {
    if (incy == 1) {
      if (beta == zero) {
        FEM_DO(i, 1, n) {
          y(i) = zero;
        }
      }
      else {
        FEM_DO(i, 1, n) {
          y(i) = beta * y(i);
        }
      }
    }
    else {
      iy = ky;
      if (beta == zero) {
        FEM_DO(i, 1, n) {
          y(iy) = zero;
          iy += incy;
        }
      }
      else {
        FEM_DO(i, 1, n) {
          y(iy) = beta * y(iy);
          iy += incy;
        }
      }
    }
  }
  if (alpha == zero) {
    return;
  }
  int j = fem::int0;
  double temp1 = fem::double0;
  double temp2 = fem::double0;
  int jx = fem::int0;
  int jy = fem::int0;
  int ix = fem::int0;
  if (lsame(uplo, "U")) {
    //C
    //C        Form  y  when A is stored in upper triangle.
    //C
    if ((incx == 1) && (incy == 1)) {
      FEM_DO(j, 1, n) {
        temp1 = alpha * x(j);
        temp2 = zero;
        {
          int fem_do_last = j - 1;
          FEM_DO(i, 1, fem_do_last) {
            y(i) += temp1 * a(i, j);
            temp2 += a(i, j) * x(i);
          }
        }
        y(j) += temp1 * a(j, j) + alpha * temp2;
      }
    }
    else {
      jx = kx;
      jy = ky;
      FEM_DO(j, 1, n) {
        temp1 = alpha * x(jx);
        temp2 = zero;
        ix = kx;
        iy = ky;
        {
          int fem_do_last = j - 1;
          FEM_DO(i, 1, fem_do_last) {
            y(iy) += temp1 * a(i, j);
            temp2 += a(i, j) * x(ix);
            ix += incx;
            iy += incy;
          }
        }
        y(jy) += temp1 * a(j, j) + alpha * temp2;
        jx += incx;
        jy += incy;
      }
    }
  }
  else {
    //C
    //C        Form  y  when A is stored in lower triangle.
    //C
    if ((incx == 1) && (incy == 1)) {
      FEM_DO(j, 1, n) {
        temp1 = alpha * x(j);
        temp2 = zero;
        y(j) += temp1 * a(j, j);
        FEM_DO(i, j + 1, n) {
          y(i) += temp1 * a(i, j);
          temp2 += a(i, j) * x(i);
        }
        y(j) += alpha * temp2;
      }
    }
    else {
      jx = kx;
      jy = ky;
      FEM_DO(j, 1, n) {
        temp1 = alpha * x(jx);
        temp2 = zero;
        y(jy) += temp1 * a(j, j);
        ix = jx;
        iy = jy;
        FEM_DO(i, j + 1, n) {
          ix += incx;
          iy += incy;
          y(iy) += temp1 * a(i, j);
          temp2 += a(i, j) * x(ix);
        }
        y(jy) += alpha * temp2;
        jx += incx;
        jy += incy;
      }
    }
  }
  //C
  //C     End of DSYMV .
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/BLAS/SRC/dsyr2.f
inline
void
dsyr2(
  str_cref uplo,
  int const& n,
  double const& alpha,
  arr_cref<double> x,
  int const& incx,
  arr_cref<double> y,
  int const& incy,
  arr_ref<double, 2> a,
  int const& lda)
{
  x(dimension(star));
  y(dimension(star));
  a(dimension(lda, star));
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DSYR2  performs the symmetric rank 2 operation
  //C
  //C     A := alpha*x*y' + alpha*y*x' + A,
  //C
  //C  where alpha is a scalar, x and y are n element vectors and A is an n
  //C  by n symmetric matrix.
  //C
  //C  Arguments
  //C  ==========
  //C
  //C  UPLO   - CHARACTER*1.
  //C           On entry, UPLO specifies whether the upper or lower
  //C           triangular part of the array A is to be referenced as
  //C           follows:
  //C
  //C              UPLO = 'U' or 'u'   Only the upper triangular part of A
  //C                                  is to be referenced.
  //C
  //C              UPLO = 'L' or 'l'   Only the lower triangular part of A
  //C                                  is to be referenced.
  //C
  //C           Unchanged on exit.
  //C
  //C  N      - INTEGER.
  //C           On entry, N specifies the order of the matrix A.
  //C           N must be at least zero.
  //C           Unchanged on exit.
  //C
  //C  ALPHA  - DOUBLE PRECISION.
  //C           On entry, ALPHA specifies the scalar alpha.
  //C           Unchanged on exit.
  //C
  //C  X      - DOUBLE PRECISION array of dimension at least
  //C           ( 1 + ( n - 1 )*abs( INCX ) ).
  //C           Before entry, the incremented array X must contain the n
  //C           element vector x.
  //C           Unchanged on exit.
  //C
  //C  INCX   - INTEGER.
  //C           On entry, INCX specifies the increment for the elements of
  //C           X. INCX must not be zero.
  //C           Unchanged on exit.
  //C
  //C  Y      - DOUBLE PRECISION array of dimension at least
  //C           ( 1 + ( n - 1 )*abs( INCY ) ).
  //C           Before entry, the incremented array Y must contain the n
  //C           element vector y.
  //C           Unchanged on exit.
  //C
  //C  INCY   - INTEGER.
  //C           On entry, INCY specifies the increment for the elements of
  //C           Y. INCY must not be zero.
  //C           Unchanged on exit.
  //C
  //C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
  //C           Before entry with  UPLO = 'U' or 'u', the leading n by n
  //C           upper triangular part of the array A must contain the upper
  //C           triangular part of the symmetric matrix and the strictly
  //C           lower triangular part of A is not referenced. On exit, the
  //C           upper triangular part of the array A is overwritten by the
  //C           upper triangular part of the updated matrix.
  //C           Before entry with UPLO = 'L' or 'l', the leading n by n
  //C           lower triangular part of the array A must contain the lower
  //C           triangular part of the symmetric matrix and the strictly
  //C           upper triangular part of A is not referenced. On exit, the
  //C           lower triangular part of the array A is overwritten by the
  //C           lower triangular part of the updated matrix.
  //C
  //C  LDA    - INTEGER.
  //C           On entry, LDA specifies the first dimension of A as declared
  //C           in the calling (sub) program. LDA must be at least
  //C           max( 1, n ).
  //C           Unchanged on exit.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Level 2 Blas routine.
  //C
  //C  -- Written on 22-October-1986.
  //C     Jack Dongarra, Argonne National Lab.
  //C     Jeremy Du Croz, Nag Central Office.
  //C     Sven Hammarling, Nag Central Office.
  //C     Richard Hanson, Sandia National Labs.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C
  //C     Test the input parameters.
  //C
  int info = 0;
  if (!lsame(uplo, "U") && !lsame(uplo, "L")) {
    info = 1;
  }
  else if (n < 0) {
    info = 2;
  }
  else if (incx == 0) {
    info = 5;
  }
  else if (incy == 0) {
    info = 7;
  }
  else if (lda < fem::max(1, n)) {
    info = 9;
  }
  if (info != 0) {
    xerbla("DSYR2 ", info);
    return;
  }
  //C
  //C     Quick return if possible.
  //C
  const double zero = 0.0e+0;
  if ((n == 0) || (alpha == zero)) {
    return;
  }
  //C
  //C     Set up the start points in X and Y if the increments are not both
  //C     unity.
  //C
  int kx = fem::int0;
  int ky = fem::int0;
  int jx = fem::int0;
  int jy = fem::int0;
  if ((incx != 1) || (incy != 1)) {
    if (incx > 0) {
      kx = 1;
    }
    else {
      kx = 1 - (n - 1) * incx;
    }
    if (incy > 0) {
      ky = 1;
    }
    else {
      ky = 1 - (n - 1) * incy;
    }
    jx = kx;
    jy = ky;
  }
  //C
  //C     Start the operations. In this version the elements of A are
  //C     accessed sequentially with one pass through the triangular part
  //C     of A.
  //C
  int j = fem::int0;
  double temp1 = fem::double0;
  double temp2 = fem::double0;
  int i = fem::int0;
  int ix = fem::int0;
  int iy = fem::int0;
  if (lsame(uplo, "U")) {
    //C
    //C        Form  A  when A is stored in the upper triangle.
    //C
    if ((incx == 1) && (incy == 1)) {
      FEM_DO(j, 1, n) {
        if ((x(j) != zero) || (y(j) != zero)) {
          temp1 = alpha * y(j);
          temp2 = alpha * x(j);
          FEM_DO(i, 1, j) {
            a(i, j) += x(i) * temp1 + y(i) * temp2;
          }
        }
      }
    }
    else {
      FEM_DO(j, 1, n) {
        if ((x(jx) != zero) || (y(jy) != zero)) {
          temp1 = alpha * y(jy);
          temp2 = alpha * x(jx);
          ix = kx;
          iy = ky;
          FEM_DO(i, 1, j) {
            a(i, j) += x(ix) * temp1 + y(iy) * temp2;
            ix += incx;
            iy += incy;
          }
        }
        jx += incx;
        jy += incy;
      }
    }
  }
  else {
    //C
    //C        Form  A  when A is stored in the lower triangle.
    //C
    if ((incx == 1) && (incy == 1)) {
      FEM_DO(j, 1, n) {
        if ((x(j) != zero) || (y(j) != zero)) {
          temp1 = alpha * y(j);
          temp2 = alpha * x(j);
          FEM_DO(i, j, n) {
            a(i, j) += x(i) * temp1 + y(i) * temp2;
          }
        }
      }
    }
    else {
      FEM_DO(j, 1, n) {
        if ((x(jx) != zero) || (y(jy) != zero)) {
          temp1 = alpha * y(jy);
          temp2 = alpha * x(jx);
          ix = jx;
          iy = jy;
          FEM_DO(i, j, n) {
            a(i, j) += x(ix) * temp1 + y(iy) * temp2;
            ix += incx;
            iy += incy;
          }
        }
        jx += incx;
        jy += incy;
      }
    }
  }
  //C
  //C     End of DSYR2 .
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/BLAS/SRC/dsyr2k.f
inline
void
dsyr2k(
  str_cref uplo,
  str_cref trans,
  int const& n,
  int const& k,
  double const& alpha,
  arr_cref<double, 2> a,
  int const& lda,
  arr_cref<double, 2> b,
  int const& ldb,
  double const& beta,
  arr_ref<double, 2> c,
  int const& ldc)
{
  a(dimension(lda, star));
  b(dimension(ldb, star));
  c(dimension(ldc, star));
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DSYR2K  performs one of the symmetric rank 2k operations
  //C
  //C     C := alpha*A*B' + alpha*B*A' + beta*C,
  //C
  //C  or
  //C
  //C     C := alpha*A'*B + alpha*B'*A + beta*C,
  //C
  //C  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
  //C  and  A and B  are  n by k  matrices  in the  first  case  and  k by n
  //C  matrices in the second case.
  //C
  //C  Arguments
  //C  ==========
  //C
  //C  UPLO   - CHARACTER*1.
  //C           On  entry,   UPLO  specifies  whether  the  upper  or  lower
  //C           triangular  part  of the  array  C  is to be  referenced  as
  //C           follows:
  //C
  //C              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
  //C                                  is to be referenced.
  //C
  //C              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
  //C                                  is to be referenced.
  //C
  //C           Unchanged on exit.
  //C
  //C  TRANS  - CHARACTER*1.
  //C           On entry,  TRANS  specifies the operation to be performed as
  //C           follows:
  //C
  //C              TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +
  //C                                        beta*C.
  //C
  //C              TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +
  //C                                        beta*C.
  //C
  //C              TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +
  //C                                        beta*C.
  //C
  //C           Unchanged on exit.
  //C
  //C  N      - INTEGER.
  //C           On entry,  N specifies the order of the matrix C.  N must be
  //C           at least zero.
  //C           Unchanged on exit.
  //C
  //C  K      - INTEGER.
  //C           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
  //C           of  columns  of the  matrices  A and B,  and on  entry  with
  //C           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
  //C           of rows of the matrices  A and B.  K must be at least  zero.
  //C           Unchanged on exit.
  //C
  //C  ALPHA  - DOUBLE PRECISION.
  //C           On entry, ALPHA specifies the scalar alpha.
  //C           Unchanged on exit.
  //C
  //C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
  //C           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
  //C           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
  //C           part of the array  A  must contain the matrix  A,  otherwise
  //C           the leading  k by n  part of the array  A  must contain  the
  //C           matrix A.
  //C           Unchanged on exit.
  //C
  //C  LDA    - INTEGER.
  //C           On entry, LDA specifies the first dimension of A as declared
  //C           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
  //C           then  LDA must be at least  max( 1, n ), otherwise  LDA must
  //C           be at least  max( 1, k ).
  //C           Unchanged on exit.
  //C
  //C  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
  //C           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
  //C           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
  //C           part of the array  B  must contain the matrix  B,  otherwise
  //C           the leading  k by n  part of the array  B  must contain  the
  //C           matrix B.
  //C           Unchanged on exit.
  //C
  //C  LDB    - INTEGER.
  //C           On entry, LDB specifies the first dimension of B as declared
  //C           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
  //C           then  LDB must be at least  max( 1, n ), otherwise  LDB must
  //C           be at least  max( 1, k ).
  //C           Unchanged on exit.
  //C
  //C  BETA   - DOUBLE PRECISION.
  //C           On entry, BETA specifies the scalar beta.
  //C           Unchanged on exit.
  //C
  //C  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
  //C           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
  //C           upper triangular part of the array C must contain the upper
  //C           triangular part  of the  symmetric matrix  and the strictly
  //C           lower triangular part of C is not referenced.  On exit, the
  //C           upper triangular part of the array  C is overwritten by the
  //C           upper triangular part of the updated matrix.
  //C           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
  //C           lower triangular part of the array C must contain the lower
  //C           triangular part  of the  symmetric matrix  and the strictly
  //C           upper triangular part of C is not referenced.  On exit, the
  //C           lower triangular part of the array  C is overwritten by the
  //C           lower triangular part of the updated matrix.
  //C
  //C  LDC    - INTEGER.
  //C           On entry, LDC specifies the first dimension of C as declared
  //C           in  the  calling  (sub)  program.   LDC  must  be  at  least
  //C           max( 1, n ).
  //C           Unchanged on exit.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Level 3 Blas routine.
  //C
  //C  -- Written on 8-February-1989.
  //C     Jack Dongarra, Argonne National Laboratory.
  //C     Iain Duff, AERE Harwell.
  //C     Jeremy Du Croz, Numerical Algorithms Group Ltd.
  //C     Sven Hammarling, Numerical Algorithms Group Ltd.
  //C
  //C  =====================================================================
  //C
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Parameters ..
  //C     ..
  //C
  //C     Test the input parameters.
  //C
  int nrowa = fem::int0;
  if (lsame(trans, "N")) {
    nrowa = n;
  }
  else {
    nrowa = k;
  }
  bool upper = lsame(uplo, "U");
  //C
  int info = 0;
  if ((!upper) && (!lsame(uplo, "L"))) {
    info = 1;
  }
  else if ((!lsame(trans, "N")) && (!lsame(trans, "T")) && (!lsame(trans,
    "C"))) {
    info = 2;
  }
  else if (n < 0) {
    info = 3;
  }
  else if (k < 0) {
    info = 4;
  }
  else if (lda < fem::max(1, nrowa)) {
    info = 7;
  }
  else if (ldb < fem::max(1, nrowa)) {
    info = 9;
  }
  else if (ldc < fem::max(1, n)) {
    info = 12;
  }
  if (info != 0) {
    xerbla("DSYR2K", info);
    return;
  }
  //C
  //C     Quick return if possible.
  //C
  const double zero = 0.0e+0;
  const double one = 1.0e+0;
  if ((n == 0) || (((alpha == zero) || (k == 0)) && (beta == one))) {
    return;
  }
  //C
  //C     And when  alpha.eq.zero.
  //C
  int j = fem::int0;
  int i = fem::int0;
  if (alpha == zero) {
    if (upper) {
      if (beta == zero) {
        FEM_DO(j, 1, n) {
          FEM_DO(i, 1, j) {
            c(i, j) = zero;
          }
        }
      }
      else {
        FEM_DO(j, 1, n) {
          FEM_DO(i, 1, j) {
            c(i, j) = beta * c(i, j);
          }
        }
      }
    }
    else {
      if (beta == zero) {
        FEM_DO(j, 1, n) {
          FEM_DO(i, j, n) {
            c(i, j) = zero;
          }
        }
      }
      else {
        FEM_DO(j, 1, n) {
          FEM_DO(i, j, n) {
            c(i, j) = beta * c(i, j);
          }
        }
      }
    }
    return;
  }
  //C
  //C     Start the operations.
  //C
  int l = fem::int0;
  double temp1 = fem::double0;
  double temp2 = fem::double0;
  if (lsame(trans, "N")) {
    //C
    //C        Form  C := alpha*A*B' + alpha*B*A' + C.
    //C
    if (upper) {
      FEM_DO(j, 1, n) {
        if (beta == zero) {
          FEM_DO(i, 1, j) {
            c(i, j) = zero;
          }
        }
        else if (beta != one) {
          FEM_DO(i, 1, j) {
            c(i, j) = beta * c(i, j);
          }
        }
        FEM_DO(l, 1, k) {
          if ((a(j, l) != zero) || (b(j, l) != zero)) {
            temp1 = alpha * b(j, l);
            temp2 = alpha * a(j, l);
            FEM_DO(i, 1, j) {
              c(i, j) += a(i, l) * temp1 + b(i, l) * temp2;
            }
          }
        }
      }
    }
    else {
      FEM_DO(j, 1, n) {
        if (beta == zero) {
          FEM_DO(i, j, n) {
            c(i, j) = zero;
          }
        }
        else if (beta != one) {
          FEM_DO(i, j, n) {
            c(i, j) = beta * c(i, j);
          }
        }
        FEM_DO(l, 1, k) {
          if ((a(j, l) != zero) || (b(j, l) != zero)) {
            temp1 = alpha * b(j, l);
            temp2 = alpha * a(j, l);
            FEM_DO(i, j, n) {
              c(i, j) += a(i, l) * temp1 + b(i, l) * temp2;
            }
          }
        }
      }
    }
  }
  else {
    //C
    //C        Form  C := alpha*A'*B + alpha*B'*A + C.
    //C
    if (upper) {
      FEM_DO(j, 1, n) {
        FEM_DO(i, 1, j) {
          temp1 = zero;
          temp2 = zero;
          FEM_DO(l, 1, k) {
            temp1 += a(l, i) * b(l, j);
            temp2 += b(l, i) * a(l, j);
          }
          if (beta == zero) {
            c(i, j) = alpha * temp1 + alpha * temp2;
          }
          else {
            c(i, j) = beta * c(i, j) + alpha * temp1 + alpha * temp2;
          }
        }
      }
    }
    else {
      FEM_DO(j, 1, n) {
        FEM_DO(i, j, n) {
          temp1 = zero;
          temp2 = zero;
          FEM_DO(l, 1, k) {
            temp1 += a(l, i) * b(l, j);
            temp2 += b(l, i) * a(l, j);
          }
          if (beta == zero) {
            c(i, j) = alpha * temp1 + alpha * temp2;
          }
          else {
            c(i, j) = beta * c(i, j) + alpha * temp1 + alpha * temp2;
          }
        }
      }
    }
  }
  //C
  //C     End of DSYR2K.
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/BLAS/SRC/dtrmm.f
inline
void
dtrmm(
  str_cref side,
  str_cref uplo,
  str_cref transa,
  str_cref diag,
  int const& m,
  int const& n,
  double const& alpha,
  arr_cref<double, 2> a,
  int const& lda,
  arr_ref<double, 2> b,
  int const& ldb)
{
  a(dimension(lda, star));
  b(dimension(ldb, star));
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DTRMM  performs one of the matrix-matrix operations
  //C
  //C     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
  //C
  //C  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
  //C  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
  //C
  //C     op( A ) = A   or   op( A ) = A'.
  //C
  //C  Arguments
  //C  ==========
  //C
  //C  SIDE   - CHARACTER*1.
  //C           On entry,  SIDE specifies whether  op( A ) multiplies B from
  //C           the left or right as follows:
  //C
  //C              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
  //C
  //C              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
  //C
  //C           Unchanged on exit.
  //C
  //C  UPLO   - CHARACTER*1.
  //C           On entry, UPLO specifies whether the matrix A is an upper or
  //C           lower triangular matrix as follows:
  //C
  //C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
  //C
  //C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
  //C
  //C           Unchanged on exit.
  //C
  //C  TRANSA - CHARACTER*1.
  //C           On entry, TRANSA specifies the form of op( A ) to be used in
  //C           the matrix multiplication as follows:
  //C
  //C              TRANSA = 'N' or 'n'   op( A ) = A.
  //C
  //C              TRANSA = 'T' or 't'   op( A ) = A'.
  //C
  //C              TRANSA = 'C' or 'c'   op( A ) = A'.
  //C
  //C           Unchanged on exit.
  //C
  //C  DIAG   - CHARACTER*1.
  //C           On entry, DIAG specifies whether or not A is unit triangular
  //C           as follows:
  //C
  //C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
  //C
  //C              DIAG = 'N' or 'n'   A is not assumed to be unit
  //C                                  triangular.
  //C
  //C           Unchanged on exit.
  //C
  //C  M      - INTEGER.
  //C           On entry, M specifies the number of rows of B. M must be at
  //C           least zero.
  //C           Unchanged on exit.
  //C
  //C  N      - INTEGER.
  //C           On entry, N specifies the number of columns of B.  N must be
  //C           at least zero.
  //C           Unchanged on exit.
  //C
  //C  ALPHA  - DOUBLE PRECISION.
  //C           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
  //C           zero then  A is not referenced and  B need not be set before
  //C           entry.
  //C           Unchanged on exit.
  //C
  //C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
  //C           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
  //C           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
  //C           upper triangular part of the array  A must contain the upper
  //C           triangular matrix  and the strictly lower triangular part of
  //C           A is not referenced.
  //C           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
  //C           lower triangular part of the array  A must contain the lower
  //C           triangular matrix  and the strictly upper triangular part of
  //C           A is not referenced.
  //C           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
  //C           A  are not referenced either,  but are assumed to be  unity.
  //C           Unchanged on exit.
  //C
  //C  LDA    - INTEGER.
  //C           On entry, LDA specifies the first dimension of A as declared
  //C           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
  //C           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
  //C           then LDA must be at least max( 1, n ).
  //C           Unchanged on exit.
  //C
  //C  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
  //C           Before entry,  the leading  m by n part of the array  B must
  //C           contain the matrix  B,  and  on exit  is overwritten  by the
  //C           transformed matrix.
  //C
  //C  LDB    - INTEGER.
  //C           On entry, LDB specifies the first dimension of B as declared
  //C           in  the  calling  (sub)  program.   LDB  must  be  at  least
  //C           max( 1, m ).
  //C           Unchanged on exit.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Level 3 Blas routine.
  //C
  //C  -- Written on 8-February-1989.
  //C     Jack Dongarra, Argonne National Laboratory.
  //C     Iain Duff, AERE Harwell.
  //C     Jeremy Du Croz, Numerical Algorithms Group Ltd.
  //C     Sven Hammarling, Numerical Algorithms Group Ltd.
  //C
  //C  =====================================================================
  //C
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Parameters ..
  //C     ..
  //C
  //C     Test the input parameters.
  //C
  bool lside = lsame(side, "L");
  int nrowa = fem::int0;
  if (lside) {
    nrowa = m;
  }
  else {
    nrowa = n;
  }
  bool nounit = lsame(diag, "N");
  bool upper = lsame(uplo, "U");
  //C
  int info = 0;
  if ((!lside) && (!lsame(side, "R"))) {
    info = 1;
  }
  else if ((!upper) && (!lsame(uplo, "L"))) {
    info = 2;
  }
  else if ((!lsame(transa, "N")) && (!lsame(transa, "T")) && (!lsame(transa,
    "C"))) {
    info = 3;
  }
  else if ((!lsame(diag, "U")) && (!lsame(diag, "N"))) {
    info = 4;
  }
  else if (m < 0) {
    info = 5;
  }
  else if (n < 0) {
    info = 6;
  }
  else if (lda < fem::max(1, nrowa)) {
    info = 9;
  }
  else if (ldb < fem::max(1, m)) {
    info = 11;
  }
  if (info != 0) {
    xerbla("DTRMM ", info);
    return;
  }
  //C
  //C     Quick return if possible.
  //C
  if (m == 0 || n == 0) {
    return;
  }
  //C
  //C     And when  alpha.eq.zero.
  //C
  const double zero = 0.0e+0;
  int j = fem::int0;
  int i = fem::int0;
  if (alpha == zero) {
    FEM_DO(j, 1, n) {
      FEM_DO(i, 1, m) {
        b(i, j) = zero;
      }
    }
    return;
  }
  //C
  //C     Start the operations.
  //C
  int k = fem::int0;
  double temp = fem::double0;
  const double one = 1.0e+0;
  if (lside) {
    if (lsame(transa, "N")) {
      //C
      //C           Form  B := alpha*A*B.
      //C
      if (upper) {
        FEM_DO(j, 1, n) {
          FEM_DO(k, 1, m) {
            if (b(k, j) != zero) {
              temp = alpha * b(k, j);
              {
                int fem_do_last = k - 1;
                FEM_DO(i, 1, fem_do_last) {
                  b(i, j) += temp * a(i, k);
                }
              }
              if (nounit) {
                temp = temp * a(k, k);
              }
              b(k, j) = temp;
            }
          }
        }
      }
      else {
        FEM_DO(j, 1, n) {
          FEM_DOSTEP(k, m, 1, -1) {
            if (b(k, j) != zero) {
              temp = alpha * b(k, j);
              b(k, j) = temp;
              if (nounit) {
                b(k, j) = b(k, j) * a(k, k);
              }
              FEM_DO(i, k + 1, m) {
                b(i, j) += temp * a(i, k);
              }
            }
          }
        }
      }
    }
    else {
      //C
      //C           Form  B := alpha*A'*B.
      //C
      if (upper) {
        FEM_DO(j, 1, n) {
          FEM_DOSTEP(i, m, 1, -1) {
            temp = b(i, j);
            if (nounit) {
              temp = temp * a(i, i);
            }
            {
              int fem_do_last = i - 1;
              FEM_DO(k, 1, fem_do_last) {
                temp += a(k, i) * b(k, j);
              }
            }
            b(i, j) = alpha * temp;
          }
        }
      }
      else {
        FEM_DO(j, 1, n) {
          FEM_DO(i, 1, m) {
            temp = b(i, j);
            if (nounit) {
              temp = temp * a(i, i);
            }
            FEM_DO(k, i + 1, m) {
              temp += a(k, i) * b(k, j);
            }
            b(i, j) = alpha * temp;
          }
        }
      }
    }
  }
  else {
    if (lsame(transa, "N")) {
      //C
      //C           Form  B := alpha*B*A.
      //C
      if (upper) {
        FEM_DOSTEP(j, n, 1, -1) {
          temp = alpha;
          if (nounit) {
            temp = temp * a(j, j);
          }
          FEM_DO(i, 1, m) {
            b(i, j) = temp * b(i, j);
          }
          {
            int fem_do_last = j - 1;
            FEM_DO(k, 1, fem_do_last) {
              if (a(k, j) != zero) {
                temp = alpha * a(k, j);
                FEM_DO(i, 1, m) {
                  b(i, j) += temp * b(i, k);
                }
              }
            }
          }
        }
      }
      else {
        FEM_DO(j, 1, n) {
          temp = alpha;
          if (nounit) {
            temp = temp * a(j, j);
          }
          FEM_DO(i, 1, m) {
            b(i, j) = temp * b(i, j);
          }
          FEM_DO(k, j + 1, n) {
            if (a(k, j) != zero) {
              temp = alpha * a(k, j);
              FEM_DO(i, 1, m) {
                b(i, j) += temp * b(i, k);
              }
            }
          }
        }
      }
    }
    else {
      //C
      //C           Form  B := alpha*B*A'.
      //C
      if (upper) {
        FEM_DO(k, 1, n) {
          {
            int fem_do_last = k - 1;
            FEM_DO(j, 1, fem_do_last) {
              if (a(j, k) != zero) {
                temp = alpha * a(j, k);
                FEM_DO(i, 1, m) {
                  b(i, j) += temp * b(i, k);
                }
              }
            }
          }
          temp = alpha;
          if (nounit) {
            temp = temp * a(k, k);
          }
          if (temp != one) {
            FEM_DO(i, 1, m) {
              b(i, k) = temp * b(i, k);
            }
          }
        }
      }
      else {
        FEM_DOSTEP(k, n, 1, -1) {
          FEM_DO(j, k + 1, n) {
            if (a(j, k) != zero) {
              temp = alpha * a(j, k);
              FEM_DO(i, 1, m) {
                b(i, j) += temp * b(i, k);
              }
            }
          }
          temp = alpha;
          if (nounit) {
            temp = temp * a(k, k);
          }
          if (temp != one) {
            FEM_DO(i, 1, m) {
              b(i, k) = temp * b(i, k);
            }
          }
        }
      }
    }
  }
  //C
  //C     End of DTRMM .
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/BLAS/SRC/dtrmv.f
inline
void
dtrmv(
  str_cref uplo,
  str_cref trans,
  str_cref diag,
  int const& n,
  arr_cref<double, 2> a,
  int const& lda,
  arr_ref<double> x,
  int const& incx)
{
  a(dimension(lda, star));
  x(dimension(star));
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DTRMV  performs one of the matrix-vector operations
  //C
  //C     x := A*x,   or   x := A'*x,
  //C
  //C  where x is an n element vector and  A is an n by n unit, or non-unit,
  //C  upper or lower triangular matrix.
  //C
  //C  Arguments
  //C  ==========
  //C
  //C  UPLO   - CHARACTER*1.
  //C           On entry, UPLO specifies whether the matrix is an upper or
  //C           lower triangular matrix as follows:
  //C
  //C              UPLO = 'U' or 'u'   A is an upper triangular matrix.
  //C
  //C              UPLO = 'L' or 'l'   A is a lower triangular matrix.
  //C
  //C           Unchanged on exit.
  //C
  //C  TRANS  - CHARACTER*1.
  //C           On entry, TRANS specifies the operation to be performed as
  //C           follows:
  //C
  //C              TRANS = 'N' or 'n'   x := A*x.
  //C
  //C              TRANS = 'T' or 't'   x := A'*x.
  //C
  //C              TRANS = 'C' or 'c'   x := A'*x.
  //C
  //C           Unchanged on exit.
  //C
  //C  DIAG   - CHARACTER*1.
  //C           On entry, DIAG specifies whether or not A is unit
  //C           triangular as follows:
  //C
  //C              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
  //C
  //C              DIAG = 'N' or 'n'   A is not assumed to be unit
  //C                                  triangular.
  //C
  //C           Unchanged on exit.
  //C
  //C  N      - INTEGER.
  //C           On entry, N specifies the order of the matrix A.
  //C           N must be at least zero.
  //C           Unchanged on exit.
  //C
  //C  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
  //C           Before entry with  UPLO = 'U' or 'u', the leading n by n
  //C           upper triangular part of the array A must contain the upper
  //C           triangular matrix and the strictly lower triangular part of
  //C           A is not referenced.
  //C           Before entry with UPLO = 'L' or 'l', the leading n by n
  //C           lower triangular part of the array A must contain the lower
  //C           triangular matrix and the strictly upper triangular part of
  //C           A is not referenced.
  //C           Note that when  DIAG = 'U' or 'u', the diagonal elements of
  //C           A are not referenced either, but are assumed to be unity.
  //C           Unchanged on exit.
  //C
  //C  LDA    - INTEGER.
  //C           On entry, LDA specifies the first dimension of A as declared
  //C           in the calling (sub) program. LDA must be at least
  //C           max( 1, n ).
  //C           Unchanged on exit.
  //C
  //C  X      - DOUBLE PRECISION array of dimension at least
  //C           ( 1 + ( n - 1 )*abs( INCX ) ).
  //C           Before entry, the incremented array X must contain the n
  //C           element vector x. On exit, X is overwritten with the
  //C           tranformed vector x.
  //C
  //C  INCX   - INTEGER.
  //C           On entry, INCX specifies the increment for the elements of
  //C           X. INCX must not be zero.
  //C           Unchanged on exit.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Level 2 Blas routine.
  //C
  //C  -- Written on 22-October-1986.
  //C     Jack Dongarra, Argonne National Lab.
  //C     Jeremy Du Croz, Nag Central Office.
  //C     Sven Hammarling, Nag Central Office.
  //C     Richard Hanson, Sandia National Labs.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C
  //C     Test the input parameters.
  //C
  int info = 0;
  if (!lsame(uplo, "U") && !lsame(uplo, "L")) {
    info = 1;
  }
  else if (!lsame(trans, "N") && !lsame(trans, "T") && !lsame(trans, "C")) {
    info = 2;
  }
  else if (!lsame(diag, "U") && !lsame(diag, "N")) {
    info = 3;
  }
  else if (n < 0) {
    info = 4;
  }
  else if (lda < fem::max(1, n)) {
    info = 6;
  }
  else if (incx == 0) {
    info = 8;
  }
  if (info != 0) {
    xerbla("DTRMV ", info);
    return;
  }
  //C
  //C     Quick return if possible.
  //C
  if (n == 0) {
    return;
  }
  //C
  bool nounit = lsame(diag, "N");
  //C
  //C     Set up the start point in X if the increment is not unity. This
  //C     will be  ( N - 1 )*INCX  too small for descending loops.
  //C
  int kx = fem::int0;
  if (incx <= 0) {
    kx = 1 - (n - 1) * incx;
  }
  else if (incx != 1) {
    kx = 1;
  }
  //C
  //C     Start the operations. In this version the elements of A are
  //C     accessed sequentially with one pass through A.
  //C
  int j = fem::int0;
  const double zero = 0.0e+0;
  double temp = fem::double0;
  int i = fem::int0;
  int jx = fem::int0;
  int ix = fem::int0;
  if (lsame(trans, "N")) {
    //C
    //C        Form  x := A*x.
    //C
    if (lsame(uplo, "U")) {
      if (incx == 1) {
        FEM_DO(j, 1, n) {
          if (x(j) != zero) {
            temp = x(j);
            {
              int fem_do_last = j - 1;
              FEM_DO(i, 1, fem_do_last) {
                x(i) += temp * a(i, j);
              }
            }
            if (nounit) {
              x(j) = x(j) * a(j, j);
            }
          }
        }
      }
      else {
        jx = kx;
        FEM_DO(j, 1, n) {
          if (x(jx) != zero) {
            temp = x(jx);
            ix = kx;
            {
              int fem_do_last = j - 1;
              FEM_DO(i, 1, fem_do_last) {
                x(ix) += temp * a(i, j);
                ix += incx;
              }
            }
            if (nounit) {
              x(jx) = x(jx) * a(j, j);
            }
          }
          jx += incx;
        }
      }
    }
    else {
      if (incx == 1) {
        FEM_DOSTEP(j, n, 1, -1) {
          if (x(j) != zero) {
            temp = x(j);
            FEM_DOSTEP(i, n, j + 1, -1) {
              x(i) += temp * a(i, j);
            }
            if (nounit) {
              x(j) = x(j) * a(j, j);
            }
          }
        }
      }
      else {
        kx += (n - 1) * incx;
        jx = kx;
        FEM_DOSTEP(j, n, 1, -1) {
          if (x(jx) != zero) {
            temp = x(jx);
            ix = kx;
            FEM_DOSTEP(i, n, j + 1, -1) {
              x(ix) += temp * a(i, j);
              ix = ix - incx;
            }
            if (nounit) {
              x(jx) = x(jx) * a(j, j);
            }
          }
          jx = jx - incx;
        }
      }
    }
  }
  else {
    //C
    //C        Form  x := A'*x.
    //C
    if (lsame(uplo, "U")) {
      if (incx == 1) {
        FEM_DOSTEP(j, n, 1, -1) {
          temp = x(j);
          if (nounit) {
            temp = temp * a(j, j);
          }
          FEM_DOSTEP(i, j - 1, 1, -1) {
            temp += a(i, j) * x(i);
          }
          x(j) = temp;
        }
      }
      else {
        jx = kx + (n - 1) * incx;
        FEM_DOSTEP(j, n, 1, -1) {
          temp = x(jx);
          ix = jx;
          if (nounit) {
            temp = temp * a(j, j);
          }
          FEM_DOSTEP(i, j - 1, 1, -1) {
            ix = ix - incx;
            temp += a(i, j) * x(ix);
          }
          x(jx) = temp;
          jx = jx - incx;
        }
      }
    }
    else {
      if (incx == 1) {
        FEM_DO(j, 1, n) {
          temp = x(j);
          if (nounit) {
            temp = temp * a(j, j);
          }
          FEM_DO(i, j + 1, n) {
            temp += a(i, j) * x(i);
          }
          x(j) = temp;
        }
      }
      else {
        jx = kx;
        FEM_DO(j, 1, n) {
          temp = x(jx);
          ix = jx;
          if (nounit) {
            temp = temp * a(j, j);
          }
          FEM_DO(i, j + 1, n) {
            ix += incx;
            temp += a(i, j) * x(ix);
          }
          x(jx) = temp;
          jx += incx;
        }
      }
    }
  }
  //C
  //C     End of DTRMV .
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/ieeeck.f
inline
int
ieeeck(
  int const& ispec,
  float const& zero,
  float const& one)
{
  int return_value = fem::int0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  IEEECK is called from the ILAENV to verify that Infinity and
  //C  possibly NaN arithmetic is safe (i.e. will not trap).
  //C
  //C  Arguments
  //C  =========
  //C
  //C  ISPEC   (input) INTEGER
  //C          Specifies whether to test just for inifinity arithmetic
  //C          or whether to test for infinity and NaN arithmetic.
  //C          = 0: Verify infinity arithmetic only.
  //C          = 1: Verify infinity and NaN arithmetic.
  //C
  //C  ZERO    (input) REAL
  //C          Must contain the value 0.0
  //C          This is passed to prevent the compiler from optimizing
  //C          away this code.
  //C
  //C  ONE     (input) REAL
  //C          Must contain the value 1.0
  //C          This is passed to prevent the compiler from optimizing
  //C          away this code.
  //C
  //C  RETURN VALUE:  INTEGER
  //C          = 0:  Arithmetic failed to produce the correct answers
  //C          = 1:  Arithmetic produced the correct answers
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Executable Statements ..
  return_value = 1;
  //C
  float posinf = one / zero;
  if (posinf <= one) {
    return_value = 0;
    return return_value;
  }
  //C
  float neginf = -one / zero;
  if (neginf >= zero) {
    return_value = 0;
    return return_value;
  }
  //C
  float negzro = one / (neginf + one);
  if (negzro != zero) {
    return_value = 0;
    return return_value;
  }
  //C
  neginf = one / negzro;
  if (neginf >= zero) {
    return_value = 0;
    return return_value;
  }
  //C
  float newzro = negzro + zero;
  if (newzro != zero) {
    return_value = 0;
    return return_value;
  }
  //C
  posinf = one / newzro;
  if (posinf <= one) {
    return_value = 0;
    return return_value;
  }
  //C
  neginf = neginf * posinf;
  if (neginf >= zero) {
    return_value = 0;
    return return_value;
  }
  //C
  posinf = posinf * posinf;
  if (posinf <= one) {
    return_value = 0;
    return return_value;
  }
  //C
  //C     Return if we were only asked to check infinity arithmetic
  //C
  if (ispec == 0) {
    return return_value;
  }
  //C
  float nan1 = posinf + neginf;
  //C
  float nan2 = posinf / neginf;
  //C
  float nan3 = posinf / posinf;
  //C
  float nan4 = posinf * zero;
  //C
  float nan5 = neginf * negzro;
  //C
  float nan6 = nan5 * 0.0f;
  //C
  if (nan1 == nan1) {
    return_value = 0;
    return return_value;
  }
  //C
  if (nan2 == nan2) {
    return_value = 0;
    return return_value;
  }
  //C
  if (nan3 == nan3) {
    return_value = 0;
    return return_value;
  }
  //C
  if (nan4 == nan4) {
    return_value = 0;
    return return_value;
  }
  //C
  if (nan5 == nan5) {
    return_value = 0;
    return return_value;
  }
  //C
  if (nan6 == nan6) {
    return_value = 0;
    return return_value;
  }
  //C
  return return_value;
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/iladlc.f
inline
int
iladlc(
  int const& m,
  int const& n,
  arr_cref<double, 2> a,
  int const& lda)
{
  int return_value = fem::int0;
  a(dimension(lda, star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2.1)                        --
  //C
  //C  -- April 2009                                                      --
  //C
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  ILADLC scans A for its last non-zero column.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix A.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix A.
  //C
  //C  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
  //C          The m by n matrix A.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A. LDA >= max(1,M).
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Quick test for the common case where one corner is non-zero.
  const double zero = 0.0e+0;
  int i = fem::int0;
  if (n == 0) {
    return_value = n;
  }
  else if (a(1, n) != zero || a(m, n) != zero) {
    return_value = n;
  }
  else {
    //C     Now scan each column from the end, returning with the first non-zero.
    FEM_DOSTEP(return_value, n, 1, -1) {
      FEM_DO(i, 1, m) {
        if (a(i, return_value) != zero) {
          return return_value;
        }
      }
    }
  }
  return return_value;
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/iladlr.f
inline
int
iladlr(
  int const& m,
  int const& n,
  arr_cref<double, 2> a,
  int const& lda)
{
  int return_value = fem::int0;
  a(dimension(lda, star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2.1)                        --
  //C
  //C  -- April 2009                                                      --
  //C
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  ILADLR scans A for its last non-zero row.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix A.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix A.
  //C
  //C  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
  //C          The m by n matrix A.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A. LDA >= max(1,M).
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Quick test for the common case where one corner is non-zero.
  const double zero = 0.0e+0;
  int j = fem::int0;
  int i = fem::int0;
  if (m == 0) {
    return_value = m;
  }
  else if (a(m, 1) != zero || a(m, n) != zero) {
    return_value = m;
  }
  else {
    //C     Scan up each column tracking the last zero row seen.
    return_value = 0;
    FEM_DO(j, 1, n) {
      FEM_DOSTEP(i, m, 1, -1) {
        if (a(i, j) != zero) {
          break;
        }
      }
      return_value = fem::max(return_value, i);
    }
  }
  return return_value;
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/iparmq.f
inline
int
iparmq(
  int const& ispec,
  str_cref /* name */,
  str_cref /* opts */,
  int const& /* n */,
  int const& ilo,
  int const& ihi,
  int const& /* lwork */)
{
  int return_value = fem::int0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C       This program sets problem and machine dependent parameters
  //C       useful for xHSEQR and its subroutines. It is called whenever
  //C       ILAENV is called with 12 <= ISPEC <= 16
  //C
  //C  Arguments
  //C  =========
  //C
  //C       ISPEC  (input) integer scalar
  //C              ISPEC specifies which tunable parameter IPARMQ should
  //C              return.
  //C
  //C              ISPEC=12: (INMIN)  Matrices of order nmin or less
  //C                        are sent directly to xLAHQR, the implicit
  //C                        double shift QR algorithm.  NMIN must be
  //C                        at least 11.
  //C
  //C              ISPEC=13: (INWIN)  Size of the deflation window.
  //C                        This is best set greater than or equal to
  //C                        the number of simultaneous shifts NS.
  //C                        Larger matrices benefit from larger deflation
  //C                        windows.
  //C
  //C              ISPEC=14: (INIBL) Determines when to stop nibbling and
  //C                        invest in an (expensive) multi-shift QR sweep.
  //C                        If the aggressive early deflation subroutine
  //C                        finds LD converged eigenvalues from an order
  //C                        NW deflation window and LD.GT.(NW*NIBBLE)/100,
  //C                        then the next QR sweep is skipped and early
  //C                        deflation is applied immediately to the
  //C                        remaining active diagonal block.  Setting
  //C                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
  //C                        multi-shift QR sweep whenever early deflation
  //C                        finds a converged eigenvalue.  Setting
  //C                        IPARMQ(ISPEC=14) greater than or equal to 100
  //C                        prevents TTQRE from skipping a multi-shift
  //C                        QR sweep.
  //C
  //C              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
  //C                        a multi-shift QR iteration.
  //C
  //C              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
  //C                        following meanings.
  //C                        0:  During the multi-shift QR sweep,
  //C                            xLAQR5 does not accumulate reflections and
  //C                            does not use matrix-matrix multiply to
  //C                            update the far-from-diagonal matrix
  //C                            entries.
  //C                        1:  During the multi-shift QR sweep,
  //C                            xLAQR5 and/or xLAQRaccumulates reflections and uses
  //C                            matrix-matrix multiply to update the
  //C                            far-from-diagonal matrix entries.
  //C                        2:  During the multi-shift QR sweep.
  //C                            xLAQR5 accumulates reflections and takes
  //C                            advantage of 2-by-2 block structure during
  //C                            matrix-matrix multiplies.
  //C                        (If xTRMM is slower than xGEMM, then
  //C                        IPARMQ(ISPEC=16)=1 may be more efficient than
  //C                        IPARMQ(ISPEC=16)=2 despite the greater level of
  //C                        arithmetic work implied by the latter choice.)
  //C
  //C       NAME    (input) character string
  //C               Name of the calling subroutine
  //C
  //C       OPTS    (input) character string
  //C               This is a concatenation of the string arguments to
  //C               TTQRE.
  //C
  //C       N       (input) integer scalar
  //C               N is the order of the Hessenberg matrix H.
  //C
  //C       ILO     (input) INTEGER
  //C       IHI     (input) INTEGER
  //C               It is assumed that H is already upper triangular
  //C               in rows and columns 1:ILO-1 and IHI+1:N.
  //C
  //C       LWORK   (input) integer scalar
  //C               The amount of workspace available.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C       Little is known about how best to choose these parameters.
  //C       It is possible to use different values of the parameters
  //C       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
  //C
  //C       It is probably best to choose different parameters for
  //C       different matrices and different parameters at different
  //C       times during the iteration, but this has not been
  //C       implemented --- yet.
  //C
  //C       The best choices of most of the parameters depend
  //C       in an ill-understood way on the relative execution
  //C       rate of xLAQR3 and xLAQR5 and on the nature of each
  //C       particular eigenvalue problem.  Experiment may be the
  //C       only practical way to determine which choices are most
  //C       effective.
  //C
  //C       Following is a list of default values supplied by IPARMQ.
  //C       These defaults may be adjusted in order to attain better
  //C       performance in any particular computational environment.
  //C
  //C       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
  //C                        Default: 75. (Must be at least 11.)
  //C
  //C       IPARMQ(ISPEC=13) Recommended deflation window size.
  //C                        This depends on ILO, IHI and NS, the
  //C                        number of simultaneous shifts returned
  //C                        by IPARMQ(ISPEC=15).  The default for
  //C                        (IHI-ILO+1).LE.500 is NS.  The default
  //C                        for (IHI-ILO+1).GT.500 is 3*NS/2.
  //C
  //C       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
  //C
  //C       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
  //C                        a multi-shift QR iteration.
  //C
  //C                        If IHI-ILO+1 is ...
  //C
  //C                        greater than      ...but less    ... the
  //C                        or equal to ...      than        default is
  //C
  //C                                0               30       NS =   2+
  //C                               30               60       NS =   4+
  //C                               60              150       NS =  10
  //C                              150              590       NS =  **
  //C                              590             3000       NS =  64
  //C                             3000             6000       NS = 128
  //C                             6000             infinity   NS = 256
  //C
  //C                    (+)  By default matrices of this order are
  //C                         passed to the implicit double shift routine
  //C                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
  //C                         values of NS are used only in case of a rare
  //C                         xLAHQR failure.
  //C
  //C                    (**) The asterisks (**) indicate an ad-hoc
  //C                         function increasing from 10 to 64.
  //C
  //C       IPARMQ(ISPEC=16) Select structured matrix multiply.
  //C                        (See ISPEC=16 above for details.)
  //C                        Default: 3.
  //C
  //C     ================================================================
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  const int ishfts = 15;
  const int inwin = 13;
  const int iacc22 = 16;
  int nh = fem::int0;
  int ns = fem::int0;
  const float two = 2.0f;
  if ((ispec == ishfts) || (ispec == inwin) || (ispec == iacc22)) {
    //C
    //C        ==== Set the number simultaneous shifts ====
    //C
    nh = ihi - ilo + 1;
    ns = 2;
    if (nh >= 30) {
      ns = 4;
    }
    if (nh >= 60) {
      ns = 10;
    }
    if (nh >= 150) {
      ns = fem::max(10, nh / fem::nint(fem::log(fem::real(nh)) /
        fem::log(two)));
    }
    if (nh >= 590) {
      ns = 64;
    }
    if (nh >= 3000) {
      ns = 128;
    }
    if (nh >= 6000) {
      ns = 256;
    }
    ns = fem::max(2, ns - fem::mod(ns, 2));
  }
  //C
  const int inmin = 12;
  const int nmin = 75;
  const int inibl = 14;
  const int nibble = 14;
  const int knwswp = 500;
  const int kacmin = 14;
  const int k22min = 14;
  if (ispec == inmin) {
    //C
    //C        ===== Matrices of order smaller than NMIN get sent
    //C        .     to xLAHQR, the classic double shift algorithm.
    //C        .     This must be at least 11. ====
    //C
    return_value = nmin;
    //C
  }
  else if (ispec == inibl) {
    //C
    //C        ==== INIBL: skip a multi-shift qr iteration and
    //C        .    whenever aggressive early deflation finds
    //C        .    at least (NIBBLE*(window size)/100) deflations. ====
    //C
    return_value = nibble;
    //C
  }
  else if (ispec == ishfts) {
    //C
    //C        ==== NSHFTS: The number of simultaneous shifts =====
    //C
    return_value = ns;
    //C
  }
  else if (ispec == inwin) {
    //C
    //C        ==== NW: deflation window size.  ====
    //C
    if (nh <= knwswp) {
      return_value = ns;
    }
    else {
      return_value = 3 * ns / 2;
    }
    //C
  }
  else if (ispec == iacc22) {
    //C
    //C        ==== IACC22: Whether to accumulate reflections
    //C        .     before updating the far-from-diagonal elements
    //C        .     and whether to use 2-by-2 block structure while
    //C        .     doing it.  A small amount of work could be saved
    //C        .     by making this choice dependent also upon the
    //C        .     NH=IHI-ILO+1.
    //C
    return_value = 0;
    if (ns >= kacmin) {
      return_value = 1;
    }
    if (ns >= k22min) {
      return_value = 2;
    }
    //C
  }
  else {
    //C        ===== invalid value of ispec =====
    return_value = -1;
    //C
  }
  return return_value;
  //C
  //C     ==== End of IPARMQ ====
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/ilaenv.f
inline
int
ilaenv(
  int const& ispec,
  str_cref name,
  str_cref opts,
  int const& n1,
  int const& n2,
  int const& n3,
  int const& n4)
{
  int return_value = fem::int0;
  fem::str<6> subnam = fem::char0;
  int ic = fem::int0;
  int iz = fem::int0;
  int i = fem::int0;
  fem::str<1> c1 = fem::char0;
  bool sname = fem::bool0;
  bool cname = fem::bool0;
  fem::str<2> c2 = fem::char0;
  fem::str<3> c3 = fem::char0;
  fem::str<2> c4 = fem::char0;
  int nb = fem::int0;
  int nbmin = fem::int0;
  int nx = fem::int0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2.1)                        --
  //C
  //C  -- April 2009                                                      --
  //C
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  ILAENV is called from the LAPACK routines to choose problem-dependent
  //C  parameters for the local environment.  See ISPEC for a description of
  //C  the parameters.
  //C
  //C  ILAENV returns an INTEGER
  //C  if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
  //C  if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
  //C
  //C  This version provides a set of parameters which should give good,
  //C  but not optimal, performance on many of the currently available
  //C  computers.  Users are encouraged to modify this subroutine to set
  //C  the tuning parameters for their particular machine using the option
  //C  and problem size information in the arguments.
  //C
  //C  This routine will not function correctly if it is converted to all
  //C  lower case.  Converting it to all upper case is allowed.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  ISPEC   (input) INTEGER
  //C          Specifies the parameter to be returned as the value of
  //C          ILAENV.
  //C          = 1: the optimal blocksize; if this value is 1, an unblocked
  //C               algorithm will give the best performance.
  //C          = 2: the minimum block size for which the block routine
  //C               should be used; if the usable block size is less than
  //C               this value, an unblocked routine should be used.
  //C          = 3: the crossover point (in a block routine, for N less
  //C               than this value, an unblocked routine should be used)
  //C          = 4: the number of shifts, used in the nonsymmetric
  //C               eigenvalue routines (DEPRECATED)
  //C          = 5: the minimum column dimension for blocking to be used;
  //C               rectangular blocks must have dimension at least k by m,
  //C               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
  //C          = 6: the crossover point for the SVD (when reducing an m by n
  //C               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
  //C               this value, a QR factorization is used first to reduce
  //C               the matrix to a triangular form.)
  //C          = 7: the number of processors
  //C          = 8: the crossover point for the multishift QR method
  //C               for nonsymmetric eigenvalue problems (DEPRECATED)
  //C          = 9: maximum size of the subproblems at the bottom of the
  //C               computation tree in the divide-and-conquer algorithm
  //C               (used by xGELSD and xGESDD)
  //C          =10: ieee NaN arithmetic can be trusted not to trap
  //C          =11: infinity arithmetic can be trusted not to trap
  //C          12 <= ISPEC <= 16:
  //C               xHSEQR or one of its subroutines,
  //C               see IPARMQ for detailed explanation
  //C
  //C  NAME    (input) CHARACTER*(*)
  //C          The name of the calling subroutine, in either upper case or
  //C          lower case.
  //C
  //C  OPTS    (input) CHARACTER*(*)
  //C          The character options to the subroutine NAME, concatenated
  //C          into a single character string.  For example, UPLO = 'U',
  //C          TRANS = 'T', and DIAG = 'N' for a triangular routine would
  //C          be specified as OPTS = 'UTN'.
  //C
  //C  N1      (input) INTEGER
  //C  N2      (input) INTEGER
  //C  N3      (input) INTEGER
  //C  N4      (input) INTEGER
  //C          Problem dimensions for the subroutine NAME; these may not all
  //C          be required.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  The following conventions have been used when calling ILAENV from the
  //C  LAPACK routines:
  //C  1)  OPTS is a concatenation of all of the character options to
  //C      subroutine NAME, in the same order that they appear in the
  //C      argument list for NAME, even if they are not used in determining
  //C      the value of the parameter specified by ISPEC.
  //C  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
  //C      that they appear in the argument list for NAME.  N1 is used
  //C      first, N2 second, and so on, and unused problem dimensions are
  //C      passed a value of -1.
  //C  3)  The parameter value returned by ILAENV is checked for validity in
  //C      the calling subroutine.  For example, ILAENV is used to retrieve
  //C      the optimal blocksize for STRTRI as follows:
  //C
  //C      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
  //C      IF( NB.LE.1 ) NB = MAX( 1, N )
  //C
  //C  =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  switch (ispec) {
    case 1: goto statement_10;
    case 2: goto statement_10;
    case 3: goto statement_10;
    case 4: goto statement_80;
    case 5: goto statement_90;
    case 6: goto statement_100;
    case 7: goto statement_110;
    case 8: goto statement_120;
    case 9: goto statement_130;
    case 10: goto statement_140;
    case 11: goto statement_150;
    case 12: goto statement_160;
    case 13: goto statement_160;
    case 14: goto statement_160;
    case 15: goto statement_160;
    case 16: goto statement_160;
    default: break;
  }
  //C
  //C     Invalid value for ISPEC
  //C
  return_value = -1;
  return return_value;
  //C
  statement_10:
  //C
  //C     Convert NAME to upper case if the first character is lower case.
  //C
  return_value = 1;
  subnam = name;
  ic = fem::ichar(subnam(1, 1));
  iz = fem::ichar("Z");
  if (iz == 90 || iz == 122) {
    //C
    //C        ASCII character set
    //C
    if (ic >= 97 && ic <= 122) {
      subnam(1, 1) = fem::fchar(ic - 32);
      FEM_DO(i, 2, 6) {
        ic = fem::ichar(subnam(i, i));
        if (ic >= 97 && ic <= 122) {
          subnam(i, i) = fem::fchar(ic - 32);
        }
      }
    }
    //C
  }
  else if (iz == 233 || iz == 169) {
    //C
    //C        EBCDIC character set
    //C
    if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (
        ic >= 162 && ic <= 169)) {
      subnam(1, 1) = fem::fchar(ic + 64);
      FEM_DO(i, 2, 6) {
        ic = fem::ichar(subnam(i, i));
        if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (
            ic >= 162 && ic <= 169)) {
          subnam(i, i) = fem::fchar(ic + 64);
        }
      }
    }
    //C
  }
  else if (iz == 218 || iz == 250) {
    //C
    //C        Prime machines:  ASCII+128
    //C
    if (ic >= 225 && ic <= 250) {
      subnam(1, 1) = fem::fchar(ic - 32);
      FEM_DO(i, 2, 6) {
        ic = fem::ichar(subnam(i, i));
        if (ic >= 225 && ic <= 250) {
          subnam(i, i) = fem::fchar(ic - 32);
        }
      }
    }
  }
  //C
  c1 = subnam(1, 1);
  sname = c1 == "S" || c1 == "D";
  cname = c1 == "C" || c1 == "Z";
  if (!(cname || sname)) {
    return return_value;
  }
  c2 = subnam(2, 3);
  c3 = subnam(4, 6);
  c4 = c3(2, 3);
  //C
  switch (ispec) {
    case 1: goto statement_50;
    case 2: goto statement_60;
    case 3: goto statement_70;
    default: break;
  }
  //C
  statement_50:
  //C
  //C     ISPEC = 1:  block size
  //C
  //C     In these examples, separate code is provided for setting NB for
  //C     real and complex.  We assume that NB will take the same value in
  //C     single or double precision.
  //C
  nb = 1;
  //C
  if (c2 == "GE") {
    if (c3 == "TRF") {
      if (sname) {
        nb = 64;
      }
      else {
        nb = 64;
      }
    }
    else if (c3 == "QRF" || c3 == "RQF" || c3 == "LQF" || c3 == "QLF") {
      if (sname) {
        nb = 32;
      }
      else {
        nb = 32;
      }
    }
    else if (c3 == "HRD") {
      if (sname) {
        nb = 32;
      }
      else {
        nb = 32;
      }
    }
    else if (c3 == "BRD") {
      if (sname) {
        nb = 32;
      }
      else {
        nb = 32;
      }
    }
    else if (c3 == "TRI") {
      if (sname) {
        nb = 64;
      }
      else {
        nb = 64;
      }
    }
  }
  else if (c2 == "PO") {
    if (c3 == "TRF") {
      if (sname) {
        nb = 64;
      }
      else {
        nb = 64;
      }
    }
  }
  else if (c2 == "SY") {
    if (c3 == "TRF") {
      if (sname) {
        nb = 64;
      }
      else {
        nb = 64;
      }
    }
    else if (sname && c3 == "TRD") {
      nb = 32;
    }
    else if (sname && c3 == "GST") {
      nb = 64;
    }
  }
  else if (cname && c2 == "HE") {
    if (c3 == "TRF") {
      nb = 64;
    }
    else if (c3 == "TRD") {
      nb = 32;
    }
    else if (c3 == "GST") {
      nb = 64;
    }
  }
  else if (sname && c2 == "OR") {
    if (c3(1, 1) == "G") {
      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" ||
          c4 == "HR" || c4 == "TR" || c4 == "BR") {
        nb = 32;
      }
    }
    else if (c3(1, 1) == "M") {
      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" ||
          c4 == "HR" || c4 == "TR" || c4 == "BR") {
        nb = 32;
      }
    }
  }
  else if (cname && c2 == "UN") {
    if (c3(1, 1) == "G") {
      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" ||
          c4 == "HR" || c4 == "TR" || c4 == "BR") {
        nb = 32;
      }
    }
    else if (c3(1, 1) == "M") {
      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" ||
          c4 == "HR" || c4 == "TR" || c4 == "BR") {
        nb = 32;
      }
    }
  }
  else if (c2 == "GB") {
    if (c3 == "TRF") {
      if (sname) {
        if (n4 <= 64) {
          nb = 1;
        }
        else {
          nb = 32;
        }
      }
      else {
        if (n4 <= 64) {
          nb = 1;
        }
        else {
          nb = 32;
        }
      }
    }
  }
  else if (c2 == "PB") {
    if (c3 == "TRF") {
      if (sname) {
        if (n2 <= 64) {
          nb = 1;
        }
        else {
          nb = 32;
        }
      }
      else {
        if (n2 <= 64) {
          nb = 1;
        }
        else {
          nb = 32;
        }
      }
    }
  }
  else if (c2 == "TR") {
    if (c3 == "TRI") {
      if (sname) {
        nb = 64;
      }
      else {
        nb = 64;
      }
    }
  }
  else if (c2 == "LA") {
    if (c3 == "UUM") {
      if (sname) {
        nb = 64;
      }
      else {
        nb = 64;
      }
    }
  }
  else if (sname && c2 == "ST") {
    if (c3 == "EBZ") {
      nb = 1;
    }
  }
  return_value = nb;
  return return_value;
  //C
  statement_60:
  //C
  //C     ISPEC = 2:  minimum block size
  //C
  nbmin = 2;
  if (c2 == "GE") {
    if (c3 == "QRF" || c3 == "RQF" || c3 == "LQF" || c3 == "QLF") {
      if (sname) {
        nbmin = 2;
      }
      else {
        nbmin = 2;
      }
    }
    else if (c3 == "HRD") {
      if (sname) {
        nbmin = 2;
      }
      else {
        nbmin = 2;
      }
    }
    else if (c3 == "BRD") {
      if (sname) {
        nbmin = 2;
      }
      else {
        nbmin = 2;
      }
    }
    else if (c3 == "TRI") {
      if (sname) {
        nbmin = 2;
      }
      else {
        nbmin = 2;
      }
    }
  }
  else if (c2 == "SY") {
    if (c3 == "TRF") {
      if (sname) {
        nbmin = 8;
      }
      else {
        nbmin = 8;
      }
    }
    else if (sname && c3 == "TRD") {
      nbmin = 2;
    }
  }
  else if (cname && c2 == "HE") {
    if (c3 == "TRD") {
      nbmin = 2;
    }
  }
  else if (sname && c2 == "OR") {
    if (c3(1, 1) == "G") {
      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" ||
          c4 == "HR" || c4 == "TR" || c4 == "BR") {
        nbmin = 2;
      }
    }
    else if (c3(1, 1) == "M") {
      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" ||
          c4 == "HR" || c4 == "TR" || c4 == "BR") {
        nbmin = 2;
      }
    }
  }
  else if (cname && c2 == "UN") {
    if (c3(1, 1) == "G") {
      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" ||
          c4 == "HR" || c4 == "TR" || c4 == "BR") {
        nbmin = 2;
      }
    }
    else if (c3(1, 1) == "M") {
      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" ||
          c4 == "HR" || c4 == "TR" || c4 == "BR") {
        nbmin = 2;
      }
    }
  }
  return_value = nbmin;
  return return_value;
  //C
  statement_70:
  //C
  //C     ISPEC = 3:  crossover point
  //C
  nx = 0;
  if (c2 == "GE") {
    if (c3 == "QRF" || c3 == "RQF" || c3 == "LQF" || c3 == "QLF") {
      if (sname) {
        nx = 128;
      }
      else {
        nx = 128;
      }
    }
    else if (c3 == "HRD") {
      if (sname) {
        nx = 128;
      }
      else {
        nx = 128;
      }
    }
    else if (c3 == "BRD") {
      if (sname) {
        nx = 128;
      }
      else {
        nx = 128;
      }
    }
  }
  else if (c2 == "SY") {
    if (sname && c3 == "TRD") {
      nx = 32;
    }
  }
  else if (cname && c2 == "HE") {
    if (c3 == "TRD") {
      nx = 32;
    }
  }
  else if (sname && c2 == "OR") {
    if (c3(1, 1) == "G") {
      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" ||
          c4 == "HR" || c4 == "TR" || c4 == "BR") {
        nx = 128;
      }
    }
  }
  else if (cname && c2 == "UN") {
    if (c3(1, 1) == "G") {
      if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" ||
          c4 == "HR" || c4 == "TR" || c4 == "BR") {
        nx = 128;
      }
    }
  }
  return_value = nx;
  return return_value;
  //C
  statement_80:
  //C
  //C     ISPEC = 4:  number of shifts (used by xHSEQR)
  //C
  return_value = 6;
  return return_value;
  //C
  statement_90:
  //C
  //C     ISPEC = 5:  minimum column dimension (not used)
  //C
  return_value = 2;
  return return_value;
  //C
  statement_100:
  //C
  //C     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
  //C
  return_value = fem::fint(fem::real(fem::min(n1, n2)) * 1.6e0f);
  return return_value;
  //C
  statement_110:
  //C
  //C     ISPEC = 7:  number of processors (not used)
  //C
  return_value = 1;
  return return_value;
  //C
  statement_120:
  //C
  //C     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
  //C
  return_value = 50;
  return return_value;
  //C
  statement_130:
  //C
  //C     ISPEC = 9:  maximum size of the subproblems at the bottom of the
  //C                 computation tree in the divide-and-conquer algorithm
  //C                 (used by xGELSD and xGESDD)
  //C
  return_value = 25;
  return return_value;
  //C
  statement_140:
  //C
  //C     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
  //C
  //C     ILAENV = 0
  return_value = 1;
  if (return_value == 1) {
    return_value = ieeeck(1, 0.0f, 1.0f);
  }
  return return_value;
  //C
  statement_150:
  //C
  //C     ISPEC = 11: infinity arithmetic can be trusted not to trap
  //C
  //C     ILAENV = 0
  return_value = 1;
  if (return_value == 1) {
    return_value = ieeeck(0, 0.0f, 1.0f);
  }
  return return_value;
  //C
  statement_160:
  //C
  //C     12 <= ISPEC <= 16: xHSEQR or one of its subroutines.
  //C
  return_value = iparmq(ispec, name, opts, n1, n2, n3, n4);
  return return_value;
  //C
  //C     End of ILAENV
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlassq.f
inline
void
dlassq(
  int const& n,
  arr_cref<double> x,
  int const& incx,
  double& scale,
  double& sumsq)
{
  x(dimension(star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASSQ  returns the values  scl  and  smsq  such that
  //C
  //C     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
  //C
  //C  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
  //C  assumed to be non-negative and  scl  returns the value
  //C
  //C     scl = max( scale, abs( x( i ) ) ).
  //C
  //C  scale and sumsq must be supplied in SCALE and SUMSQ and
  //C  scl and smsq are overwritten on SCALE and SUMSQ respectively.
  //C
  //C  The routine makes only one pass through the vector x.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  N       (input) INTEGER
  //C          The number of elements to be used from the vector X.
  //C
  //C  X       (input) DOUBLE PRECISION array, dimension (N)
  //C          The vector for which a scaled sum of squares is computed.
  //C             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
  //C
  //C  INCX    (input) INTEGER
  //C          The increment between successive values of the vector X.
  //C          INCX > 0.
  //C
  //C  SCALE   (input/output) DOUBLE PRECISION
  //C          On entry, the value  scale  in the equation above.
  //C          On exit, SCALE is overwritten with  scl , the scaling factor
  //C          for the sum of squares.
  //C
  //C  SUMSQ   (input/output) DOUBLE PRECISION
  //C          On entry, the value  sumsq  in the equation above.
  //C          On exit, SUMSQ is overwritten with  smsq , the basic sum of
  //C          squares from which  scl  has been factored out.
  //C
  //C =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  int ix = fem::int0;
  const double zero = 0.0e+0;
  double absxi = fem::double0;
  if (n > 0) {
    FEM_DOSTEP(ix, 1, 1 + (n - 1) * incx, incx) {
      if (x(ix) != zero) {
        absxi = fem::abs(x(ix));
        if (scale < absxi) {
          sumsq = 1 + sumsq * fem::pow2((scale / absxi));
          scale = absxi;
        }
        else {
          sumsq += fem::pow2((absxi / scale));
        }
      }
    }
  }
  //C
  //C     End of DLASSQ
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlanst.f
inline
double
dlanst(
  str_cref norm,
  int const& n,
  arr_cref<double> d,
  arr_cref<double> e)
{
  double return_value = fem::double0;
  d(dimension(star));
  e(dimension(star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLANST  returns the value of the one norm,  or the Frobenius norm, or
  //C  the  infinity norm,  or the  element of  largest absolute value  of a
  //C  real symmetric tridiagonal matrix A.
  //C
  //C  Description
  //C  ===========
  //C
  //C  DLANST returns the value
  //C
  //C     DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'
  //C              (
  //C              ( norm1(A),         NORM = '1', 'O' or 'o'
  //C              (
  //C              ( normI(A),         NORM = 'I' or 'i'
  //C              (
  //C              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
  //C
  //C  where  norm1  denotes the  one norm of a matrix (maximum column sum),
  //C  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
  //C  normF  denotes the  Frobenius norm of a matrix (square root of sum of
  //C  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  NORM    (input) CHARACTER*1
  //C          Specifies the value to be returned in DLANST as described
  //C          above.
  //C
  //C  N       (input) INTEGER
  //C          The order of the matrix A.  N >= 0.  When N = 0, DLANST is
  //C          set to zero.
  //C
  //C  D       (input) DOUBLE PRECISION array, dimension (N)
  //C          The diagonal elements of A.
  //C
  //C  E       (input) DOUBLE PRECISION array, dimension (N-1)
  //C          The (n-1) sub-diagonal or super-diagonal elements of A.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  const double zero = 0.0e+0;
  double anorm = fem::double0;
  int i = fem::int0;
  double scale = fem::double0;
  const double one = 1.0e+0;
  double sum = fem::double0;
  if (n <= 0) {
    anorm = zero;
  }
  else if (lsame(norm, "M")) {
    //C
    //C        Find max(abs(A(i,j))).
    //C
    anorm = fem::abs(d(n));
    {
      int fem_do_last = n - 1;
      FEM_DO(i, 1, fem_do_last) {
        anorm = fem::max(anorm, fem::abs(d(i)));
        anorm = fem::max(anorm, fem::abs(e(i)));
      }
    }
  }
  else if (lsame(norm, "O") || norm == "1" || lsame(norm, "I")) {
    //C
    //C        Find norm1(A).
    //C
    if (n == 1) {
      anorm = fem::abs(d(1));
    }
    else {
      anorm = fem::max(fem::abs(d(1)) + fem::abs(e(1)), fem::abs(e(
        n - 1)) + fem::abs(d(n)));
      {
        int fem_do_last = n - 1;
        FEM_DO(i, 2, fem_do_last) {
          anorm = fem::max(anorm, fem::abs(d(i)) + fem::abs(e(i)) +
            fem::abs(e(i - 1)));
        }
      }
    }
  }
  else if ((lsame(norm, "F")) || (lsame(norm, "E"))) {
    //C
    //C        Find normF(A).
    //C
    scale = zero;
    sum = one;
    if (n > 1) {
      dlassq(n - 1, e, 1, scale, sum);
      sum = 2 * sum;
    }
    dlassq(n, d, 1, scale, sum);
    anorm = scale * fem::sqrt(sum);
  }
  //C
  return_value = anorm;
  return return_value;
  //C
  //C     End of DLANST
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlartg.f
inline
void
dlartg(
  common& cmn,
  double const& f,
  double const& g,
  double& cs,
  double& sn,
  double& r)
{
  double safmin = fem::double0;
  double eps = fem::double0;
  const double two = 2.0e0;
  double safmn2 = fem::double0;
  const double one = 1.0e0;
  double safmx2 = fem::double0;
  const double zero = 0.0e0;
  double f1 = fem::double0;
  double g1 = fem::double0;
  double scale = fem::double0;
  int count = fem::int0;
  int i = fem::int0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLARTG generate a plane rotation so that
  //C
  //C     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
  //C     [ -SN  CS  ]     [ G ]     [ 0 ]
  //C
  //C  This is a slower, more accurate version of the BLAS1 routine DROTG,
  //C  with the following other differences:
  //C     F and G are unchanged on return.
  //C     If G=0, then CS=1 and SN=0.
  //C     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
  //C        floating point operations (saves work in DBDSQR when
  //C        there are zeros on the diagonal).
  //C
  //C  If F exceeds G in magnitude, CS will be positive.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  F       (input) DOUBLE PRECISION
  //C          The first component of vector to be rotated.
  //C
  //C  G       (input) DOUBLE PRECISION
  //C          The second component of vector to be rotated.
  //C
  //C  CS      (output) DOUBLE PRECISION
  //C          The cosine of the rotation.
  //C
  //C  SN      (output) DOUBLE PRECISION
  //C          The sine of the rotation.
  //C
  //C  R       (output) DOUBLE PRECISION
  //C          The nonzero component of the rotated vector.
  //C
  //C  This version has a few statements commented out for thread safety
  //C  (machine parameters are computed on each entry). 10 feb 03, SJH.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     LOGICAL            FIRST
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Save statement ..
  //C     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
  //C     ..
  //C     .. Data statements ..
  //C     DATA               FIRST / .TRUE. /
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     IF( FIRST ) THEN
  safmin = dlamch(cmn, "S");
  eps = dlamch(cmn, "E");
  safmn2 = fem::pow(dlamch(cmn, "B"), fem::fint(fem::log(safmin /
    eps) / fem::log(dlamch(cmn, "B")) / two));
  safmx2 = one / safmn2;
  //C        FIRST = .FALSE.
  //C     END IF
  if (g == zero) {
    cs = one;
    sn = zero;
    r = f;
  }
  else if (f == zero) {
    cs = zero;
    sn = one;
    r = g;
  }
  else {
    f1 = f;
    g1 = g;
    scale = fem::max(fem::abs(f1), fem::abs(g1));
    if (scale >= safmx2) {
      count = 0;
      statement_10:
      count++;
      f1 = f1 * safmn2;
      g1 = g1 * safmn2;
      scale = fem::max(fem::abs(f1), fem::abs(g1));
      if (scale >= safmx2) {
        goto statement_10;
      }
      r = fem::sqrt(fem::pow2(f1) + fem::pow2(g1));
      cs = f1 / r;
      sn = g1 / r;
      FEM_DO(i, 1, count) {
        r = r * safmx2;
      }
    }
    else if (scale <= safmn2) {
      count = 0;
      statement_30:
      count++;
      f1 = f1 * safmx2;
      g1 = g1 * safmx2;
      scale = fem::max(fem::abs(f1), fem::abs(g1));
      if (scale <= safmn2) {
        goto statement_30;
      }
      r = fem::sqrt(fem::pow2(f1) + fem::pow2(g1));
      cs = f1 / r;
      sn = g1 / r;
      FEM_DO(i, 1, count) {
        r = r * safmn2;
      }
    }
    else {
      r = fem::sqrt(fem::pow2(f1) + fem::pow2(g1));
      cs = f1 / r;
      sn = g1 / r;
    }
    if (fem::abs(f) > fem::abs(g) && cs < zero) {
      cs = -cs;
      sn = -sn;
      r = -r;
    }
  }
  //C
  //C     End of DLARTG
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlaisnan.f
inline
bool
dlaisnan(
  double const& din1,
  double const& din2)
{
  bool return_value = fem::bool0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  This routine is not for general use.  It exists solely to avoid
  //C  over-optimization in DISNAN.
  //C
  //C  DLAISNAN checks for NaNs by comparing its two arguments for
  //C  inequality.  NaN is the only floating-point value where NaN != NaN
  //C  returns .TRUE.  To check for NaNs, pass the same variable as both
  //C  arguments.
  //C
  //C  A compiler must assume that the two arguments are
  //C  not the same variable, and the test will not be optimized away.
  //C  Interprocedural or whole-program optimization may delete this
  //C  test.  The ISNAN functions will be replaced by the correct
  //C  Fortran 03 intrinsic once the intrinsic is widely available.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  DIN1     (input) DOUBLE PRECISION
  //C  DIN2     (input) DOUBLE PRECISION
  //C          Two numbers to compare for inequality.
  //C
  //C  =====================================================================
  //C
  //C  .. Executable Statements ..
  return_value = (din1 != din2);
  return return_value;
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/disnan.f
inline
bool
disnan(
  double const& din)
{
  bool return_value = fem::bool0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
  //C  otherwise.  To be replaced by the Fortran 2003 intrinsic in the
  //C  future.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  DIN      (input) DOUBLE PRECISION
  //C          Input to test for NaN.
  //C
  //C  =====================================================================
  //C
  //C  .. External Functions ..
  //C  ..
  //C  .. Executable Statements ..
  return_value = dlaisnan(din, din);
  return return_value;
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlascl.f
inline
void
dlascl(
  common& cmn,
  str_cref type,
  int const& kl,
  int const& ku,
  double const& cfrom,
  double const& cto,
  int const& m,
  int const& n,
  arr_ref<double, 2> a,
  int const& lda,
  int& info)
{
  a(dimension(lda, star));
  int itype = fem::int0;
  const double zero = 0.0e0;
  double smlnum = fem::double0;
  const double one = 1.0e0;
  double bignum = fem::double0;
  double cfromc = fem::double0;
  double ctoc = fem::double0;
  double cfrom1 = fem::double0;
  double mul = fem::double0;
  bool done = fem::bool0;
  double cto1 = fem::double0;
  int j = fem::int0;
  int i = fem::int0;
  int k3 = fem::int0;
  int k4 = fem::int0;
  int k1 = fem::int0;
  int k2 = fem::int0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASCL multiplies the M by N real matrix A by the real scalar
  //C  CTO/CFROM.  This is done without over/underflow as long as the final
  //C  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
  //C  A may be full, upper triangular, lower triangular, upper Hessenberg,
  //C  or banded.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  TYPE    (input) CHARACTER*1
  //C          TYPE indices the storage type of the input matrix.
  //C          = 'G':  A is a full matrix.
  //C          = 'L':  A is a lower triangular matrix.
  //C          = 'U':  A is an upper triangular matrix.
  //C          = 'H':  A is an upper Hessenberg matrix.
  //C          = 'B':  A is a symmetric band matrix with lower bandwidth KL
  //C                  and upper bandwidth KU and with the only the lower
  //C                  half stored.
  //C          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
  //C                  and upper bandwidth KU and with the only the upper
  //C                  half stored.
  //C          = 'Z':  A is a band matrix with lower bandwidth KL and upper
  //C                  bandwidth KU.
  //C
  //C  KL      (input) INTEGER
  //C          The lower bandwidth of A.  Referenced only if TYPE = 'B',
  //C          'Q' or 'Z'.
  //C
  //C  KU      (input) INTEGER
  //C          The upper bandwidth of A.  Referenced only if TYPE = 'B',
  //C          'Q' or 'Z'.
  //C
  //C  CFROM   (input) DOUBLE PRECISION
  //C  CTO     (input) DOUBLE PRECISION
  //C          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
  //C          without over/underflow if the final result CTO*A(I,J)/CFROM
  //C          can be represented without over/underflow.  CFROM must be
  //C          nonzero.
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix A.  M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix A.  N >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
  //C          storage type.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(1,M).
  //C
  //C  INFO    (output) INTEGER
  //C          0  - successful exit
  //C          <0 - if INFO = -i, the i-th argument had an illegal value.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  //C
  if (lsame(type, "G")) {
    itype = 0;
  }
  else if (lsame(type, "L")) {
    itype = 1;
  }
  else if (lsame(type, "U")) {
    itype = 2;
  }
  else if (lsame(type, "H")) {
    itype = 3;
  }
  else if (lsame(type, "B")) {
    itype = 4;
  }
  else if (lsame(type, "Q")) {
    itype = 5;
  }
  else if (lsame(type, "Z")) {
    itype = 6;
  }
  else {
    itype = -1;
  }
  //C
  if (itype ==  - 1) {
    info = -1;
  }
  else if (cfrom == zero || disnan(cfrom)) {
    info = -4;
  }
  else if (disnan(cto)) {
    info = -5;
  }
  else if (m < 0) {
    info = -6;
  }
  else if (n < 0 || (itype == 4 && n != m) || (itype == 5 && n != m)) {
    info = -7;
  }
  else if (itype <= 3 && lda < fem::max(1, m)) {
    info = -9;
  }
  else if (itype >= 4) {
    if (kl < 0 || kl > fem::max(m - 1, 0)) {
      info = -2;
    }
    else if (ku < 0 || ku > fem::max(n - 1, 0) || ((itype == 4 ||
      itype == 5) && kl != ku)) {
      info = -3;
    }
    else if ((itype == 4 && lda < kl + 1) || (itype == 5 &&
      lda < ku + 1) || (itype == 6 && lda < 2 * kl + ku + 1)) {
      info = -9;
    }
  }
  //C
  if (info != 0) {
    xerbla("DLASCL", -info);
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (n == 0 || m == 0) {
    return;
  }
  //C
  //C     Get machine parameters
  //C
  smlnum = dlamch(cmn, "S");
  bignum = one / smlnum;
  //C
  cfromc = cfrom;
  ctoc = cto;
  //C
  statement_10:
  cfrom1 = cfromc * smlnum;
  if (cfrom1 == cfromc) {
    //C        CFROMC is an inf.  Multiply by a correctly signed zero for
    //C        finite CTOC, or a NaN if CTOC is infinite.
    mul = ctoc / cfromc;
    done = true;
    cto1 = ctoc;
  }
  else {
    cto1 = ctoc / bignum;
    if (cto1 == ctoc) {
      //C           CTOC is either 0 or an inf.  In both cases, CTOC itself
      //C           serves as the correct multiplication factor.
      mul = ctoc;
      done = true;
      cfromc = one;
    }
    else if (fem::abs(cfrom1) > fem::abs(ctoc) && ctoc != zero) {
      mul = smlnum;
      done = false;
      cfromc = cfrom1;
    }
    else if (fem::abs(cto1) > fem::abs(cfromc)) {
      mul = bignum;
      done = false;
      ctoc = cto1;
    }
    else {
      mul = ctoc / cfromc;
      done = true;
    }
  }
  //C
  if (itype == 0) {
    //C
    //C        Full matrix
    //C
    FEM_DO(j, 1, n) {
      FEM_DO(i, 1, m) {
        a(i, j) = a(i, j) * mul;
      }
    }
    //C
  }
  else if (itype == 1) {
    //C
    //C        Lower triangular matrix
    //C
    FEM_DO(j, 1, n) {
      FEM_DO(i, j, m) {
        a(i, j) = a(i, j) * mul;
      }
    }
    //C
  }
  else if (itype == 2) {
    //C
    //C        Upper triangular matrix
    //C
    FEM_DO(j, 1, n) {
      {
        int fem_do_last = fem::min(j, m);
        FEM_DO(i, 1, fem_do_last) {
          a(i, j) = a(i, j) * mul;
        }
      }
    }
    //C
  }
  else if (itype == 3) {
    //C
    //C        Upper Hessenberg matrix
    //C
    FEM_DO(j, 1, n) {
      {
        int fem_do_last = fem::min(j + 1, m);
        FEM_DO(i, 1, fem_do_last) {
          a(i, j) = a(i, j) * mul;
        }
      }
    }
    //C
  }
  else if (itype == 4) {
    //C
    //C        Lower half of a symmetric band matrix
    //C
    k3 = kl + 1;
    k4 = n + 1;
    FEM_DO(j, 1, n) {
      {
        int fem_do_last = fem::min(k3, k4 - j);
        FEM_DO(i, 1, fem_do_last) {
          a(i, j) = a(i, j) * mul;
        }
      }
    }
    //C
  }
  else if (itype == 5) {
    //C
    //C        Upper half of a symmetric band matrix
    //C
    k1 = ku + 2;
    k3 = ku + 1;
    FEM_DO(j, 1, n) {
      FEM_DO(i, fem::max(k1 - j, 1), k3) {
        a(i, j) = a(i, j) * mul;
      }
    }
    //C
  }
  else if (itype == 6) {
    //C
    //C        Band matrix
    //C
    k1 = kl + ku + 2;
    k2 = kl + 1;
    k3 = 2 * kl + ku + 1;
    k4 = kl + ku + 1 + m;
    FEM_DO(j, 1, n) {
      {
        int fem_do_last = fem::min(k3, k4 - j);
        FEM_DO(i, fem::max(k1 - j, k2), fem_do_last) {
          a(i, j) = a(i, j) * mul;
        }
      }
    }
    //C
  }
  //C
  if (!done) {
    goto statement_10;
  }
  //C
  //C     End of DLASCL
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlamrg.f
inline
void
dlamrg(
  int const& n1,
  int const& n2,
  arr_cref<double> a,
  int const& dtrd1,
  int const& dtrd2,
  arr_ref<int> index)
{
  a(dimension(star));
  index(dimension(star));
  int n1sv = fem::int0;
  int n2sv = fem::int0;
  int ind1 = fem::int0;
  int ind2 = fem::int0;
  int i = fem::int0;
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLAMRG will create a permutation list which will merge the elements
  //C  of A (which is composed of two independently sorted sets) into a
  //C  single set which is sorted in ascending order.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  N1     (input) INTEGER
  //C  N2     (input) INTEGER
  //C         These arguements contain the respective lengths of the two
  //C         sorted lists to be merged.
  //C
  //C  A      (input) DOUBLE PRECISION array, dimension (N1+N2)
  //C         The first N1 elements of A contain a list of numbers which
  //C         are sorted in either ascending or descending order.  Likewise
  //C         for the final N2 elements.
  //C
  //C  DTRD1  (input) INTEGER
  //C  DTRD2  (input) INTEGER
  //C         These are the strides to be taken through the array A.
  //C         Allowable strides are 1 and -1.  They indicate whether a
  //C         subset of A is sorted in ascending (DTRDx = 1) or descending
  //C         (DTRDx = -1) order.
  //C
  //C  INDEX  (output) INTEGER array, dimension (N1+N2)
  //C         On exit this array will contain a permutation such that
  //C         if B( I ) = A( INDEX( I ) ) for I=1,N1+N2, then B will be
  //C         sorted in ascending order.
  //C
  //C  =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  n1sv = n1;
  n2sv = n2;
  if (dtrd1 > 0) {
    ind1 = 1;
  }
  else {
    ind1 = n1;
  }
  if (dtrd2 > 0) {
    ind2 = 1 + n1;
  }
  else {
    ind2 = n1 + n2;
  }
  i = 1;
  //C     while ( (N1SV > 0) & (N2SV > 0) )
  statement_10:
  if (n1sv > 0 && n2sv > 0) {
    if (a(ind1) <= a(ind2)) {
      index(i) = ind1;
      i++;
      ind1 += dtrd1;
      n1sv = n1sv - 1;
    }
    else {
      index(i) = ind2;
      i++;
      ind2 += dtrd2;
      n2sv = n2sv - 1;
    }
    goto statement_10;
  }
  //C     end while
  if (n1sv == 0) {
    FEM_DO(n1sv, 1, n2sv) {
      index(i) = ind2;
      i++;
      ind2 += dtrd2;
    }
  }
  else {
    //C     N2SV .EQ. 0
    FEM_DO(n2sv, 1, n1sv) {
      index(i) = ind1;
      i++;
      ind1 += dtrd1;
    }
  }
  //C
  //C     End of DLAMRG
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlacpy.f
inline
void
dlacpy(
  str_cref uplo,
  int const& m,
  int const& n,
  arr_cref<double, 2> a,
  int const& lda,
  arr_ref<double, 2> b,
  int const& ldb)
{
  a(dimension(lda, star));
  b(dimension(ldb, star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLACPY copies all or part of a two-dimensional matrix A to another
  //C  matrix B.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  UPLO    (input) CHARACTER*1
  //C          Specifies the part of the matrix A to be copied to B.
  //C          = 'U':      Upper triangular part
  //C          = 'L':      Lower triangular part
  //C          Otherwise:  All of the matrix A
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix A.  M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix A.  N >= 0.
  //C
  //C  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
  //C          The m by n matrix A.  If UPLO = 'U', only the upper triangle
  //C          or trapezoid is accessed; if UPLO = 'L', only the lower
  //C          triangle or trapezoid is accessed.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(1,M).
  //C
  //C  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
  //C          On exit, B = A in the locations specified by UPLO.
  //C
  //C  LDB     (input) INTEGER
  //C          The leading dimension of the array B.  LDB >= max(1,M).
  //C
  //C  =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  int j = fem::int0;
  int i = fem::int0;
  if (lsame(uplo, "U")) {
    FEM_DO(j, 1, n) {
      {
        int fem_do_last = fem::min(j, m);
        FEM_DO(i, 1, fem_do_last) {
          b(i, j) = a(i, j);
        }
      }
    }
  }
  else if (lsame(uplo, "L")) {
    FEM_DO(j, 1, n) {
      FEM_DO(i, j, m) {
        b(i, j) = a(i, j);
      }
    }
  }
  else {
    FEM_DO(j, 1, n) {
      FEM_DO(i, 1, m) {
        b(i, j) = a(i, j);
      }
    }
  }
  //C
  //C     End of DLACPY
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlapy2.f
inline
double
dlapy2(
  double const& x,
  double const& y)
{
  double return_value = fem::double0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
  //C  overflow.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  X       (input) DOUBLE PRECISION
  //C  Y       (input) DOUBLE PRECISION
  //C          X and Y specify the values x and y.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  double xabs = fem::abs(x);
  double yabs = fem::abs(y);
  double w = fem::max(xabs, yabs);
  double z = fem::min(xabs, yabs);
  const double zero = 0.0e0;
  const double one = 1.0e0;
  if (z == zero) {
    return_value = w;
  }
  else {
    return_value = w * fem::sqrt(one + fem::pow2((z / w)));
  }
  return return_value;
  //C
  //C     End of DLAPY2
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlaset.f
inline
void
dlaset(
  str_cref uplo,
  int const& m,
  int const& n,
  double const& alpha,
  double const& beta,
  arr_ref<double, 2> a,
  int const& lda)
{
  a(dimension(lda, star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASET initializes an m-by-n matrix A to BETA on the diagonal and
  //C  ALPHA on the offdiagonals.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  UPLO    (input) CHARACTER*1
  //C          Specifies the part of the matrix A to be set.
  //C          = 'U':      Upper triangular part is set; the strictly lower
  //C                      triangular part of A is not changed.
  //C          = 'L':      Lower triangular part is set; the strictly upper
  //C                      triangular part of A is not changed.
  //C          Otherwise:  All of the matrix A is set.
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix A.  M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix A.  N >= 0.
  //C
  //C  ALPHA   (input) DOUBLE PRECISION
  //C          The constant to which the offdiagonal elements are to be set.
  //C
  //C  BETA    (input) DOUBLE PRECISION
  //C          The constant to which the diagonal elements are to be set.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On exit, the leading m-by-n submatrix of A is set as follows:
  //C
  //C          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
  //C          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
  //C          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
  //C
  //C          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(1,M).
  //C
  //C =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  int j = fem::int0;
  int i = fem::int0;
  if (lsame(uplo, "U")) {
    //C
    //C        Set the strictly upper triangular or trapezoidal part of the
    //C        array to ALPHA.
    //C
    FEM_DO(j, 2, n) {
      {
        int fem_do_last = fem::min(j - 1, m);
        FEM_DO(i, 1, fem_do_last) {
          a(i, j) = alpha;
        }
      }
    }
    //C
  }
  else if (lsame(uplo, "L")) {
    //C
    //C        Set the strictly lower triangular or trapezoidal part of the
    //C        array to ALPHA.
    //C
    {
      int fem_do_last = fem::min(m, n);
      FEM_DO(j, 1, fem_do_last) {
        FEM_DO(i, j + 1, m) {
          a(i, j) = alpha;
        }
      }
    }
    //C
  }
  else {
    //C
    //C        Set the leading m-by-n submatrix to ALPHA.
    //C
    FEM_DO(j, 1, n) {
      FEM_DO(i, 1, m) {
        a(i, j) = alpha;
      }
    }
  }
  //C
  //C     Set the first min(M,N) diagonal elements to BETA.
  //C
  {
    int fem_do_last = fem::min(m, n);
    FEM_DO(i, 1, fem_do_last) {
      a(i, i) = beta;
    }
  }
  //C
  //C     End of DLASET
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasd2.f
inline
void
dlasd2(
  common& cmn,
  int const& nl,
  int const& nr,
  int const& sqre,
  int& k,
  arr_ref<double> d,
  arr_ref<double> z,
  double const& alpha,
  double const& beta,
  arr_ref<double, 2> u,
  int const& ldu,
  arr_ref<double, 2> vt,
  int const& ldvt,
  arr_ref<double> dsigma,
  arr_ref<double, 2> u2,
  int const& ldu2,
  arr_ref<double, 2> vt2,
  int const& ldvt2,
  arr_ref<int> idxp,
  arr_ref<int> idx,
  arr_ref<int> idxc,
  arr_ref<int> idxq,
  arr_ref<int> coltyp,
  int& info)
{
  d(dimension(star));
  z(dimension(star));
  u(dimension(ldu, star));
  vt(dimension(ldvt, star));
  dsigma(dimension(star));
  u2(dimension(ldu2, star));
  vt2(dimension(ldvt2, star));
  idxp(dimension(star));
  idx(dimension(star));
  idxc(dimension(star));
  idxq(dimension(star));
  coltyp(dimension(star));
  int n = fem::int0;
  int m = fem::int0;
  int nlp1 = fem::int0;
  int nlp2 = fem::int0;
  double z1 = fem::double0;
  int i = fem::int0;
  int idxi = fem::int0;
  double eps = fem::double0;
  double tol = fem::double0;
  const double eight = 8.0e+0;
  int k2 = fem::int0;
  int j = fem::int0;
  int jprev = fem::int0;
  double s = fem::double0;
  double c = fem::double0;
  double tau = fem::double0;
  const double zero = 0.0e+0;
  int idxjp = fem::int0;
  int idxj = fem::int0;
  arr_1d<4, int> ctot(fem::fill0);
  int ct = fem::int0;
  arr_1d<4, int> psm(fem::fill0);
  int jp = fem::int0;
  const double two = 2.0e+0;
  double hlftol = fem::double0;
  const double one = 1.0e+0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASD2 merges the two sets of singular values together into a single
  //C  sorted set.  Then it tries to deflate the size of the problem.
  //C  There are two ways in which deflation can occur:  when two or more
  //C  singular values are close together or if there is a tiny entry in the
  //C  Z vector.  For each such occurrence the order of the related secular
  //C  equation problem is reduced by one.
  //C
  //C  DLASD2 is called from DLASD1.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  NL     (input) INTEGER
  //C         The row dimension of the upper block.  NL >= 1.
  //C
  //C  NR     (input) INTEGER
  //C         The row dimension of the lower block.  NR >= 1.
  //C
  //C  SQRE   (input) INTEGER
  //C         = 0: the lower block is an NR-by-NR square matrix.
  //C         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
  //C
  //C         The bidiagonal matrix has N = NL + NR + 1 rows and
  //C         M = N + SQRE >= N columns.
  //C
  //C  K      (output) INTEGER
  //C         Contains the dimension of the non-deflated matrix,
  //C         This is the order of the related secular equation. 1 <= K <=N.
  //C
  //C  D      (input/output) DOUBLE PRECISION array, dimension(N)
  //C         On entry D contains the singular values of the two submatrices
  //C         to be combined.  On exit D contains the trailing (N-K) updated
  //C         singular values (those which were deflated) sorted into
  //C         increasing order.
  //C
  //C  Z      (output) DOUBLE PRECISION array, dimension(N)
  //C         On exit Z contains the updating row vector in the secular
  //C         equation.
  //C
  //C  ALPHA  (input) DOUBLE PRECISION
  //C         Contains the diagonal element associated with the added row.
  //C
  //C  BETA   (input) DOUBLE PRECISION
  //C         Contains the off-diagonal element associated with the added
  //C         row.
  //C
  //C  U      (input/output) DOUBLE PRECISION array, dimension(LDU,N)
  //C         On entry U contains the left singular vectors of two
  //C         submatrices in the two square blocks with corners at (1,1),
  //C         (NL, NL), and (NL+2, NL+2), (N,N).
  //C         On exit U contains the trailing (N-K) updated left singular
  //C         vectors (those which were deflated) in its last N-K columns.
  //C
  //C  LDU    (input) INTEGER
  //C         The leading dimension of the array U.  LDU >= N.
  //C
  //C  VT     (input/output) DOUBLE PRECISION array, dimension(LDVT,M)
  //C         On entry VT' contains the right singular vectors of two
  //C         submatrices in the two square blocks with corners at (1,1),
  //C         (NL+1, NL+1), and (NL+2, NL+2), (M,M).
  //C         On exit VT' contains the trailing (N-K) updated right singular
  //C         vectors (those which were deflated) in its last N-K columns.
  //C         In case SQRE =1, the last row of VT spans the right null
  //C         space.
  //C
  //C  LDVT   (input) INTEGER
  //C         The leading dimension of the array VT.  LDVT >= M.
  //C
  //C  DSIGMA (output) DOUBLE PRECISION array, dimension (N)
  //C         Contains a copy of the diagonal elements (K-1 singular values
  //C         and one zero) in the secular equation.
  //C
  //C  U2     (output) DOUBLE PRECISION array, dimension(LDU2,N)
  //C         Contains a copy of the first K-1 left singular vectors which
  //C         will be used by DLASD3 in a matrix multiply (DGEMM) to solve
  //C         for the new left singular vectors. U2 is arranged into four
  //C         blocks. The first block contains a column with 1 at NL+1 and
  //C         zero everywhere else; the second block contains non-zero
  //C         entries only at and above NL; the third contains non-zero
  //C         entries only below NL+1; and the fourth is dense.
  //C
  //C  LDU2   (input) INTEGER
  //C         The leading dimension of the array U2.  LDU2 >= N.
  //C
  //C  VT2    (output) DOUBLE PRECISION array, dimension(LDVT2,N)
  //C         VT2' contains a copy of the first K right singular vectors
  //C         which will be used by DLASD3 in a matrix multiply (DGEMM) to
  //C         solve for the new right singular vectors. VT2 is arranged into
  //C         three blocks. The first block contains a row that corresponds
  //C         to the special 0 diagonal element in SIGMA; the second block
  //C         contains non-zeros only at and before NL +1; the third block
  //C         contains non-zeros only at and after  NL +2.
  //C
  //C  LDVT2  (input) INTEGER
  //C         The leading dimension of the array VT2.  LDVT2 >= M.
  //C
  //C  IDXP   (workspace) INTEGER array dimension(N)
  //C         This will contain the permutation used to place deflated
  //C         values of D at the end of the array. On output IDXP(2:K)
  //C         points to the nondeflated D-values and IDXP(K+1:N)
  //C         points to the deflated singular values.
  //C
  //C  IDX    (workspace) INTEGER array dimension(N)
  //C         This will contain the permutation used to sort the contents of
  //C         D into ascending order.
  //C
  //C  IDXC   (output) INTEGER array dimension(N)
  //C         This will contain the permutation used to arrange the columns
  //C         of the deflated U matrix into three groups:  the first group
  //C         contains non-zero entries only at and above NL, the second
  //C         contains non-zero entries only below NL+2, and the third is
  //C         dense.
  //C
  //C  IDXQ   (input/output) INTEGER array dimension(N)
  //C         This contains the permutation which separately sorts the two
  //C         sub-problems in D into ascending order.  Note that entries in
  //C         the first hlaf of this permutation must first be moved one
  //C         position backward; and entries in the second half
  //C         must first have NL+1 added to their values.
  //C
  //C  COLTYP (workspace/output) INTEGER array dimension(N)
  //C         As workspace, this will contain a label which will indicate
  //C         which of the following types a column in the U2 matrix or a
  //C         row in the VT2 matrix is:
  //C         1 : non-zero in the upper half only
  //C         2 : non-zero in the lower half only
  //C         3 : dense
  //C         4 : deflated
  //C
  //C         On exit, it is an array of dimension 4, with COLTYP(I) being
  //C         the dimension of the I-th type columns.
  //C
  //C  INFO   (output) INTEGER
  //C          = 0:  successful exit.
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Based on contributions by
  //C     Ming Gu and Huan Ren, Computer Science Division, University of
  //C     California at Berkeley, USA
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Arrays ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters.
  //C
  info = 0;
  //C
  if (nl < 1) {
    info = -1;
  }
  else if (nr < 1) {
    info = -2;
  }
  else if ((sqre != 1) && (sqre != 0)) {
    info = -3;
  }
  //C
  n = nl + nr + 1;
  m = n + sqre;
  //C
  if (ldu < n) {
    info = -10;
  }
  else if (ldvt < m) {
    info = -12;
  }
  else if (ldu2 < n) {
    info = -15;
  }
  else if (ldvt2 < m) {
    info = -17;
  }
  if (info != 0) {
    xerbla("DLASD2", -info);
    return;
  }
  //C
  nlp1 = nl + 1;
  nlp2 = nl + 2;
  //C
  //C     Generate the first part of the vector Z; and move the singular
  //C     values in the first part of D one position backward.
  //C
  z1 = alpha * vt(nlp1, nlp1);
  z(1) = z1;
  FEM_DOSTEP(i, nl, 1, -1) {
    z(i + 1) = alpha * vt(i, nlp1);
    d(i + 1) = d(i);
    idxq(i + 1) = idxq(i) + 1;
  }
  //C
  //C     Generate the second part of the vector Z.
  //C
  FEM_DO(i, nlp2, m) {
    z(i) = beta * vt(i, nlp2);
  }
  //C
  //C     Initialize some reference arrays.
  //C
  FEM_DO(i, 2, nlp1) {
    coltyp(i) = 1;
  }
  FEM_DO(i, nlp2, n) {
    coltyp(i) = 2;
  }
  //C
  //C     Sort the singular values into increasing order
  //C
  FEM_DO(i, nlp2, n) {
    idxq(i) += nlp1;
  }
  //C
  //C     DSIGMA, IDXC, IDXC, and the first column of U2
  //C     are used as storage space.
  //C
  FEM_DO(i, 2, n) {
    dsigma(i) = d(idxq(i));
    u2(i, 1) = z(idxq(i));
    idxc(i) = coltyp(idxq(i));
  }
  //C
  dlamrg(nl, nr, dsigma(2), 1, 1, idx(2));
  //C
  FEM_DO(i, 2, n) {
    idxi = 1 + idx(i);
    d(i) = dsigma(idxi);
    z(i) = u2(idxi, 1);
    coltyp(i) = idxc(idxi);
  }
  //C
  //C     Calculate the allowable deflation tolerance
  //C
  eps = dlamch(cmn, "Epsilon");
  tol = fem::max(fem::abs(alpha), fem::abs(beta));
  tol = eight * eps * fem::max(fem::abs(d(n)), tol);
  //C
  //C     There are 2 kinds of deflation -- first a value in the z-vector
  //C     is small, second two (or more) singular values are very close
  //C     together (their difference is small).
  //C
  //C     If the value in the z-vector is small, we simply permute the
  //C     array so that the corresponding singular value is moved to the
  //C     end.
  //C
  //C     If two values in the D-vector are close, we perform a two-sided
  //C     rotation designed to make one of the corresponding z-vector
  //C     entries zero, and then permute the array so that the deflated
  //C     singular value is moved to the end.
  //C
  //C     If there are multiple singular values then the problem deflates.
  //C     Here the number of equal singular values are found.  As each equal
  //C     singular value is found, an elementary reflector is computed to
  //C     rotate the corresponding singular subspace so that the
  //C     corresponding components of Z are zero in this new basis.
  //C
  k = 1;
  k2 = n + 1;
  FEM_DO(j, 2, n) {
    if (fem::abs(z(j)) <= tol) {
      //C
      //C           Deflate due to small z component.
      //C
      k2 = k2 - 1;
      idxp(k2) = j;
      coltyp(j) = 4;
      if (j == n) {
        goto statement_120;
      }
    }
    else {
      jprev = j;
      goto statement_90;
    }
  }
  statement_90:
  j = jprev;
  statement_100:
  j++;
  if (j > n) {
    goto statement_110;
  }
  if (fem::abs(z(j)) <= tol) {
    //C
    //C        Deflate due to small z component.
    //C
    k2 = k2 - 1;
    idxp(k2) = j;
    coltyp(j) = 4;
  }
  else {
    //C
    //C        Check if singular values are close enough to allow deflation.
    //C
    if (fem::abs(d(j) - d(jprev)) <= tol) {
      //C
      //C           Deflation is possible.
      //C
      s = z(jprev);
      c = z(j);
      //C
      //C           Find sqrt(a**2+b**2) without overflow or
      //C           destructive underflow.
      //C
      tau = dlapy2(c, s);
      c = c / tau;
      s = -s / tau;
      z(j) = tau;
      z(jprev) = zero;
      //C
      //C           Apply back the Givens rotation to the left and right
      //C           singular vector matrices.
      //C
      idxjp = idxq(idx(jprev) + 1);
      idxj = idxq(idx(j) + 1);
      if (idxjp <= nlp1) {
        idxjp = idxjp - 1;
      }
      if (idxj <= nlp1) {
        idxj = idxj - 1;
      }
      drot(n, u(1, idxjp), 1, u(1, idxj), 1, c, s);
      drot(m, vt(idxjp, 1), ldvt, vt(idxj, 1), ldvt, c, s);
      if (coltyp(j) != coltyp(jprev)) {
        coltyp(j) = 3;
      }
      coltyp(jprev) = 4;
      k2 = k2 - 1;
      idxp(k2) = jprev;
      jprev = j;
    }
    else {
      k++;
      u2(k, 1) = z(jprev);
      dsigma(k) = d(jprev);
      idxp(k) = jprev;
      jprev = j;
    }
  }
  goto statement_100;
  statement_110:
  //C
  //C     Record the last singular value.
  //C
  k++;
  u2(k, 1) = z(jprev);
  dsigma(k) = d(jprev);
  idxp(k) = jprev;
  //C
  statement_120:
  //C
  //C     Count up the total number of the various types of columns, then
  //C     form a permutation which positions the four column types into
  //C     four groups of uniform structure (although one or more of these
  //C     groups may be empty).
  //C
  FEM_DO(j, 1, 4) {
    ctot(j) = 0;
  }
  FEM_DO(j, 2, n) {
    ct = coltyp(j);
    ctot(ct)++;
  }
  //C
  //C     PSM(*) = Position in SubMatrix (of types 1 through 4)
  //C
  psm(1) = 2;
  psm(2) = 2 + ctot(1);
  psm(3) = psm(2) + ctot(2);
  psm(4) = psm(3) + ctot(3);
  //C
  //C     Fill out the IDXC array so that the permutation which it induces
  //C     will place all type-1 columns first, all type-2 columns next,
  //C     then all type-3's, and finally all type-4's, starting from the
  //C     second column. This applies similarly to the rows of VT.
  //C
  FEM_DO(j, 2, n) {
    jp = idxp(j);
    ct = coltyp(jp);
    idxc(psm(ct)) = j;
    psm(ct)++;
  }
  //C
  //C     Sort the singular values and corresponding singular vectors into
  //C     DSIGMA, U2, and VT2 respectively.  The singular values/vectors
  //C     which were not deflated go into the first K slots of DSIGMA, U2,
  //C     and VT2 respectively, while those which were deflated go into the
  //C     last N - K slots, except that the first column/row will be treated
  //C     separately.
  //C
  FEM_DO(j, 2, n) {
    jp = idxp(j);
    dsigma(j) = d(jp);
    idxj = idxq(idx(idxp(idxc(j))) + 1);
    if (idxj <= nlp1) {
      idxj = idxj - 1;
    }
    dcopy(n, u(1, idxj), 1, u2(1, j), 1);
    dcopy(m, vt(idxj, 1), ldvt, vt2(j, 1), ldvt2);
  }
  //C
  //C     Determine DSIGMA(1), DSIGMA(2) and Z(1)
  //C
  dsigma(1) = zero;
  hlftol = tol / two;
  if (fem::abs(dsigma(2)) <= hlftol) {
    dsigma(2) = hlftol;
  }
  if (m > n) {
    z(1) = dlapy2(z1, z(m));
    if (z(1) <= tol) {
      c = one;
      s = zero;
      z(1) = tol;
    }
    else {
      c = z1 / z(1);
      s = z(m) / z(1);
    }
  }
  else {
    if (fem::abs(z1) <= tol) {
      z(1) = tol;
    }
    else {
      z(1) = z1;
    }
  }
  //C
  //C     Move the rest of the updating row to Z.
  //C
  dcopy(k - 1, u2(2, 1), 1, z(2), 1);
  //C
  //C     Determine the first column of U2, the first row of VT2 and the
  //C     last row of VT.
  //C
  dlaset("A", n, 1, zero, zero, u2, ldu2);
  u2(nlp1, 1) = one;
  if (m > n) {
    FEM_DO(i, 1, nlp1) {
      vt(m, i) = -s * vt(nlp1, i);
      vt2(1, i) = c * vt(nlp1, i);
    }
    FEM_DO(i, nlp2, m) {
      vt2(1, i) = s * vt(m, i);
      vt(m, i) = c * vt(m, i);
    }
  }
  else {
    dcopy(m, vt(nlp1, 1), ldvt, vt2(1, 1), ldvt2);
  }
  if (m > n) {
    dcopy(m, vt(m, 1), ldvt, vt2(m, 1), ldvt2);
  }
  //C
  //C     The deflated singular values and their corresponding vectors go
  //C     into the back of D, U, and V respectively.
  //C
  if (n > k) {
    dcopy(n - k, dsigma(k + 1), 1, d(k + 1), 1);
    dlacpy("A", n, n - k, u2(1, k + 1), ldu2, u(1, k + 1), ldu);
    dlacpy("A", n - k, m, vt2(k + 1, 1), ldvt2, vt(k + 1, 1), ldvt);
  }
  //C
  //C     Copy CTOT into COLTYP for referencing in DLASD3.
  //C
  FEM_DO(j, 1, 4) {
    coltyp(j) = ctot(j);
  }
  //C
  //C     End of DLASD2
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlaed6.f
inline
void
dlaed6(
  common& cmn,
  int const& kniter,
  bool const& orgati,
  double const& rho,
  arr_cref<double> d,
  arr_cref<double> z,
  double const& finit,
  double& tau,
  int& info)
{
  d(dimension(3));
  z(dimension(3));
  double lbd = fem::double0;
  double ubd = fem::double0;
  const double zero = 0.0e0;
  int niter = fem::int0;
  const double two = 2.0e0;
  double temp = fem::double0;
  double c = fem::double0;
  double a = fem::double0;
  double b = fem::double0;
  const double four = 4.0e0;
  double eps = fem::double0;
  double base = fem::double0;
  const double three = 3.0e0;
  double small1 = fem::double0;
  const double one = 1.0e0;
  double sminv1 = fem::double0;
  double small2 = fem::double0;
  double sminv2 = fem::double0;
  bool scale = fem::bool0;
  double sclfac = fem::double0;
  double sclinv = fem::double0;
  int i = fem::int0;
  arr_1d<3, double> dscale(fem::fill0);
  arr_1d<3, double> zscale(fem::fill0);
  double fc = fem::double0;
  double df = fem::double0;
  double ddf = fem::double0;
  double temp1 = fem::double0;
  double temp2 = fem::double0;
  double temp3 = fem::double0;
  double f = fem::double0;
  int iter = fem::int0;
  const int maxit = 40;
  double eta = fem::double0;
  double erretm = fem::double0;
  double temp4 = fem::double0;
  const double eight = 8.0e0;
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     February 2007
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLAED6 computes the positive or negative root (closest to the origin)
  //C  of
  //C                   z(1)        z(2)        z(3)
  //C  f(x) =   rho + --------- + ---------- + ---------
  //C                  d(1)-x      d(2)-x      d(3)-x
  //C
  //C  It is assumed that
  //C
  //C        if ORGATI = .true. the root is between d(2) and d(3);
  //C        otherwise it is between d(1) and d(2)
  //C
  //C  This routine will be called by DLAED4 when necessary. In most cases,
  //C  the root sought is the smallest in magnitude, though it might not be
  //C  in some extremely rare situations.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  KNITER       (input) INTEGER
  //C               Refer to DLAED4 for its significance.
  //C
  //C  ORGATI       (input) LOGICAL
  //C               If ORGATI is true, the needed root is between d(2) and
  //C               d(3); otherwise it is between d(1) and d(2).  See
  //C               DLAED4 for further details.
  //C
  //C  RHO          (input) DOUBLE PRECISION
  //C               Refer to the equation f(x) above.
  //C
  //C  D            (input) DOUBLE PRECISION array, dimension (3)
  //C               D satisfies d(1) < d(2) < d(3).
  //C
  //C  Z            (input) DOUBLE PRECISION array, dimension (3)
  //C               Each of the elements in z must be positive.
  //C
  //C  FINIT        (input) DOUBLE PRECISION
  //C               The value of f at 0. It is more accurate than the one
  //C               evaluated inside this routine (if someone wants to do
  //C               so).
  //C
  //C  TAU          (output) DOUBLE PRECISION
  //C               The root of the equation f(x).
  //C
  //C  INFO         (output) INTEGER
  //C               = 0: successful exit
  //C               > 0: if INFO = 1, failure to converge
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  30/06/99: Based on contributions by
  //C     Ren-Cang Li, Computer Science Division, University of California
  //C     at Berkeley, USA
  //C
  //C  10/02/03: This version has a few statements commented out for thread
  //C  safety (machine parameters are computed on each entry). SJH.
  //C
  //C  05/10/06: Modified from a new version of Ren-Cang Li, use
  //C     Gragg-Thornton-Warner cubic convergent scheme for better stability.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Local Arrays ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  info = 0;
  //C
  if (orgati) {
    lbd = d(2);
    ubd = d(3);
  }
  else {
    lbd = d(1);
    ubd = d(2);
  }
  if (finit < zero) {
    lbd = zero;
  }
  else {
    ubd = zero;
  }
  //C
  niter = 1;
  tau = zero;
  if (kniter == 2) {
    if (orgati) {
      temp = (d(3) - d(2)) / two;
      c = rho + z(1) / ((d(1) - d(2)) - temp);
      a = c * (d(2) + d(3)) + z(2) + z(3);
      b = c * d(2) * d(3) + z(2) * d(3) + z(3) * d(2);
    }
    else {
      temp = (d(1) - d(2)) / two;
      c = rho + z(3) / ((d(3) - d(2)) - temp);
      a = c * (d(1) + d(2)) + z(1) + z(2);
      b = c * d(1) * d(2) + z(1) * d(2) + z(2) * d(1);
    }
    temp = fem::max(fem::abs(a), fem::abs(b), fem::abs(c));
    a = a / temp;
    b = b / temp;
    c = c / temp;
    if (c == zero) {
      tau = b / a;
    }
    else if (a <= zero) {
      tau = (a - fem::sqrt(fem::abs(a * a - four * b * c))) / (two * c);
    }
    else {
      tau = two * b / (a + fem::sqrt(fem::abs(a * a - four * b * c)));
    }
    if (tau < lbd || tau > ubd) {
      tau = (lbd + ubd) / two;
    }
    if (d(1) == tau || d(2) == tau || d(3) == tau) {
      tau = zero;
    }
    else {
      temp = finit + tau * z(1) / (d(1) * (d(1) - tau)) + tau * z(
        2) / (d(2) * (d(2) - tau)) + tau * z(3) / (d(3) * (d(3) -
        tau));
      if (temp <= zero) {
        lbd = tau;
      }
      else {
        ubd = tau;
      }
      if (fem::abs(finit) <= fem::abs(temp)) {
        tau = zero;
      }
    }
  }
  //C
  //C     get machine parameters for possible scaling to avoid overflow
  //C
  //C     modified by Sven: parameters SMALL1, SMINV1, SMALL2,
  //C     SMINV2, EPS are not SAVEd anymore between one call to the
  //C     others but recomputed at each call
  //C
  eps = dlamch(cmn, "Epsilon");
  base = dlamch(cmn, "Base");
  small1 = fem::pow(base, (fem::fint(fem::log(dlamch(cmn,
    "SafMin")) / fem::log(base) / three)));
  sminv1 = one / small1;
  small2 = small1 * small1;
  sminv2 = sminv1 * sminv1;
  //C
  //C     Determine if scaling of inputs necessary to avoid overflow
  //C     when computing 1/TEMP**3
  //C
  if (orgati) {
    temp = fem::min(fem::abs(d(2) - tau), fem::abs(d(3) - tau));
  }
  else {
    temp = fem::min(fem::abs(d(1) - tau), fem::abs(d(2) - tau));
  }
  scale = false;
  if (temp <= small1) {
    scale = true;
    if (temp <= small2) {
      //C
      //C        Scale up by power of radix nearest 1/SAFMIN**(2/3)
      //C
      sclfac = sminv2;
      sclinv = small2;
    }
    else {
      //C
      //C        Scale up by power of radix nearest 1/SAFMIN**(1/3)
      //C
      sclfac = sminv1;
      sclinv = small1;
    }
    //C
    //C        Scaling up safe because D, Z, TAU scaled elsewhere to be O(1)
    //C
    FEM_DO(i, 1, 3) {
      dscale(i) = d(i) * sclfac;
      zscale(i) = z(i) * sclfac;
    }
    tau = tau * sclfac;
    lbd = lbd * sclfac;
    ubd = ubd * sclfac;
  }
  else {
    //C
    //C        Copy D and Z to DSCALE and ZSCALE
    //C
    FEM_DO(i, 1, 3) {
      dscale(i) = d(i);
      zscale(i) = z(i);
    }
  }
  //C
  fc = zero;
  df = zero;
  ddf = zero;
  FEM_DO(i, 1, 3) {
    temp = one / (dscale(i) - tau);
    temp1 = zscale(i) * temp;
    temp2 = temp1 * temp;
    temp3 = temp2 * temp;
    fc += temp1 / dscale(i);
    df += temp2;
    ddf += temp3;
  }
  f = finit + tau * fc;
  //C
  if (fem::abs(f) <= zero) {
    goto statement_60;
  }
  if (f <= zero) {
    lbd = tau;
  }
  else {
    ubd = tau;
  }
  //C
  //C        Iteration begins -- Use Gragg-Thornton-Warner cubic convergent
  //C                            scheme
  //C
  //C     It is not hard to see that
  //C
  //C           1) Iterations will go up monotonically
  //C              if FINIT < 0;
  //C
  //C           2) Iterations will go down monotonically
  //C              if FINIT > 0.
  //C
  iter = niter + 1;
  //C
  FEM_DO(niter, iter, maxit) {
    //C
    if (orgati) {
      temp1 = dscale(2) - tau;
      temp2 = dscale(3) - tau;
    }
    else {
      temp1 = dscale(1) - tau;
      temp2 = dscale(2) - tau;
    }
    a = (temp1 + temp2) * f - temp1 * temp2 * df;
    b = temp1 * temp2 * f;
    c = f - (temp1 + temp2) * df + temp1 * temp2 * ddf;
    temp = fem::max(fem::abs(a), fem::abs(b), fem::abs(c));
    a = a / temp;
    b = b / temp;
    c = c / temp;
    if (c == zero) {
      eta = b / a;
    }
    else if (a <= zero) {
      eta = (a - fem::sqrt(fem::abs(a * a - four * b * c))) / (two * c);
    }
    else {
      eta = two * b / (a + fem::sqrt(fem::abs(a * a - four * b * c)));
    }
    if (f * eta >= zero) {
      eta = -f / df;
    }
    //C
    tau += eta;
    if (tau < lbd || tau > ubd) {
      tau = (lbd + ubd) / two;
    }
    //C
    fc = zero;
    erretm = zero;
    df = zero;
    ddf = zero;
    FEM_DO(i, 1, 3) {
      temp = one / (dscale(i) - tau);
      temp1 = zscale(i) * temp;
      temp2 = temp1 * temp;
      temp3 = temp2 * temp;
      temp4 = temp1 / dscale(i);
      fc += temp4;
      erretm += fem::abs(temp4);
      df += temp2;
      ddf += temp3;
    }
    f = finit + tau * fc;
    erretm = eight * (fem::abs(finit) + fem::abs(tau) * erretm) +
      fem::abs(tau) * df;
    if (fem::abs(f) <= eps * erretm) {
      goto statement_60;
    }
    if (f <= zero) {
      lbd = tau;
    }
    else {
      ubd = tau;
    }
  }
  info = 1;
  statement_60:
  //C
  //C     Undo scaling
  //C
  if (scale) {
    tau = tau * sclinv;
  }
  //C
  //C     End of DLAED6
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasd5.f
inline
void
dlasd5(
  int const& i,
  arr_cref<double> d,
  arr_cref<double> z,
  arr_ref<double> delta,
  double const& rho,
  double& dsigma,
  arr_ref<double> work)
{
  d(dimension(2));
  z(dimension(2));
  delta(dimension(2));
  work(dimension(2));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  This subroutine computes the square root of the I-th eigenvalue
  //C  of a positive symmetric rank-one modification of a 2-by-2 diagonal
  //C  matrix
  //C
  //C             diag( D ) * diag( D ) +  RHO *  Z * transpose(Z) .
  //C
  //C  The diagonal entries in the array D are assumed to satisfy
  //C
  //C             0 <= D(i) < D(j)  for  i < j .
  //C
  //C  We also assume RHO > 0 and that the Euclidean norm of the vector
  //C  Z is one.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  I      (input) INTEGER
  //C         The index of the eigenvalue to be computed.  I = 1 or I = 2.
  //C
  //C  D      (input) DOUBLE PRECISION array, dimension ( 2 )
  //C         The original eigenvalues.  We assume 0 <= D(1) < D(2).
  //C
  //C  Z      (input) DOUBLE PRECISION array, dimension ( 2 )
  //C         The components of the updating vector.
  //C
  //C  DELTA  (output) DOUBLE PRECISION array, dimension ( 2 )
  //C         Contains (D(j) - sigma_I) in its  j-th component.
  //C         The vector DELTA contains the information necessary
  //C         to construct the eigenvectors.
  //C
  //C  RHO    (input) DOUBLE PRECISION
  //C         The scalar in the symmetric updating formula.
  //C
  //C  DSIGMA (output) DOUBLE PRECISION
  //C         The computed sigma_I, the I-th updated eigenvalue.
  //C
  //C  WORK   (workspace) DOUBLE PRECISION array, dimension ( 2 )
  //C         WORK contains (D(j) + sigma_I) in its  j-th component.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Based on contributions by
  //C     Ren-Cang Li, Computer Science Division, University of California
  //C     at Berkeley, USA
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  double del = d(2) - d(1);
  double delsq = del * (d(2) + d(1));
  const double one = 1.0e+0;
  const double four = 4.0e+0;
  const double three = 3.0e+0;
  double w = fem::double0;
  const double zero = 0.0e+0;
  double b = fem::double0;
  double c = fem::double0;
  const double two = 2.0e+0;
  double tau = fem::double0;
  if (i == 1) {
    w = one + four * rho * (z(2) * z(2) / (d(1) + three * d(2)) - z(
      1) * z(1) / (three * d(1) + d(2))) / del;
    if (w > zero) {
      b = delsq + rho * (z(1) * z(1) + z(2) * z(2));
      c = rho * z(1) * z(1) * delsq;
      //C
      //C           B > ZERO, always
      //C
      //C           The following TAU is DSIGMA * DSIGMA - D( 1 ) * D( 1 )
      //C
      tau = two * c / (b + fem::sqrt(fem::abs(b * b - four * c)));
      //C
      //C           The following TAU is DSIGMA - D( 1 )
      //C
      tau = tau / (d(1) + fem::sqrt(d(1) * d(1) + tau));
      dsigma = d(1) + tau;
      delta(1) = -tau;
      delta(2) = del - tau;
      work(1) = two * d(1) + tau;
      work(2) = (d(1) + tau) + d(2);
      //C           DELTA( 1 ) = -Z( 1 ) / TAU
      //C           DELTA( 2 ) = Z( 2 ) / ( DEL-TAU )
    }
    else {
      b = -delsq + rho * (z(1) * z(1) + z(2) * z(2));
      c = rho * z(2) * z(2) * delsq;
      //C
      //C           The following TAU is DSIGMA * DSIGMA - D( 2 ) * D( 2 )
      //C
      if (b > zero) {
        tau = -two * c / (b + fem::sqrt(b * b + four * c));
      }
      else {
        tau = (b - fem::sqrt(b * b + four * c)) / two;
      }
      //C
      //C           The following TAU is DSIGMA - D( 2 )
      //C
      tau = tau / (d(2) + fem::sqrt(fem::abs(d(2) * d(2) + tau)));
      dsigma = d(2) + tau;
      delta(1) = -(del + tau);
      delta(2) = -tau;
      work(1) = d(1) + tau + d(2);
      work(2) = two * d(2) + tau;
      //C           DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
      //C           DELTA( 2 ) = -Z( 2 ) / TAU
    }
    //C        TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
    //C        DELTA( 1 ) = DELTA( 1 ) / TEMP
    //C        DELTA( 2 ) = DELTA( 2 ) / TEMP
  }
  else {
    //C
    //C        Now I=2
    //C
    b = -delsq + rho * (z(1) * z(1) + z(2) * z(2));
    c = rho * z(2) * z(2) * delsq;
    //C
    //C        The following TAU is DSIGMA * DSIGMA - D( 2 ) * D( 2 )
    //C
    if (b > zero) {
      tau = (b + fem::sqrt(b * b + four * c)) / two;
    }
    else {
      tau = two * c / (-b + fem::sqrt(b * b + four * c));
    }
    //C
    //C        The following TAU is DSIGMA - D( 2 )
    //C
    tau = tau / (d(2) + fem::sqrt(d(2) * d(2) + tau));
    dsigma = d(2) + tau;
    delta(1) = -(del + tau);
    delta(2) = -tau;
    work(1) = d(1) + tau + d(2);
    work(2) = two * d(2) + tau;
    //C        DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
    //C        DELTA( 2 ) = -Z( 2 ) / TAU
    //C        TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
    //C        DELTA( 1 ) = DELTA( 1 ) / TEMP
    //C        DELTA( 2 ) = DELTA( 2 ) / TEMP
  }
  //C
  //C     End of DLASD5
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasd4.f
inline
void
dlasd4(
  common& cmn,
  int const& n,
  int const& i,
  arr_cref<double> d,
  arr_cref<double> z,
  arr_ref<double> delta,
  double const& rho,
  double& sigma,
  arr_ref<double> work,
  int& info)
{
  d(dimension(star));
  z(dimension(star));
  delta(dimension(star));
  work(dimension(star));
  const double one = 1.0e+0;
  double eps = fem::double0;
  double rhoinv = fem::double0;
  int ii = fem::int0;
  int niter = fem::int0;
  const double two = 2.0e+0;
  double temp = fem::double0;
  double temp1 = fem::double0;
  int j = fem::int0;
  const double zero = 0.0e+0;
  double psi = fem::double0;
  double c = fem::double0;
  double w = fem::double0;
  double tau = fem::double0;
  double delsq = fem::double0;
  double a = fem::double0;
  double b = fem::double0;
  const double four = 4.0e+0;
  double eta = fem::double0;
  double dpsi = fem::double0;
  double erretm = fem::double0;
  double phi = fem::double0;
  double dphi = fem::double0;
  const double eight = 8.0e+0;
  double dtnsq1 = fem::double0;
  double dtnsq = fem::double0;
  int iter = fem::int0;
  const int maxit = 20;
  int ip1 = fem::int0;
  double delsq2 = fem::double0;
  bool orgati = fem::bool0;
  double sg2lb = fem::double0;
  double sg2ub = fem::double0;
  int iim1 = fem::int0;
  int iip1 = fem::int0;
  bool swtch3 = fem::bool0;
  double dw = fem::double0;
  const double three = 3.0e+0;
  double dtipsq = fem::double0;
  double dtisq = fem::double0;
  double dtiim = fem::double0;
  double dtiip = fem::double0;
  arr_1d<3, double> zz(fem::fill0);
  arr_1d<3, double> dd(fem::fill0);
  double prew = fem::double0;
  bool swtch = fem::bool0;
  const double ten = 10.0e+0;
  double temp2 = fem::double0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  This subroutine computes the square root of the I-th updated
  //C  eigenvalue of a positive symmetric rank-one modification to
  //C  a positive diagonal matrix whose entries are given as the squares
  //C  of the corresponding entries in the array d, and that
  //C
  //C         0 <= D(i) < D(j)  for  i < j
  //C
  //C  and that RHO > 0. This is arranged by the calling routine, and is
  //C  no loss in generality.  The rank-one modified system is thus
  //C
  //C         diag( D ) * diag( D ) +  RHO *  Z * Z_transpose.
  //C
  //C  where we assume the Euclidean norm of Z is 1.
  //C
  //C  The method consists of approximating the rational functions in the
  //C  secular equation by simpler interpolating rational functions.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  N      (input) INTEGER
  //C         The length of all arrays.
  //C
  //C  I      (input) INTEGER
  //C         The index of the eigenvalue to be computed.  1 <= I <= N.
  //C
  //C  D      (input) DOUBLE PRECISION array, dimension ( N )
  //C         The original eigenvalues.  It is assumed that they are in
  //C         order, 0 <= D(I) < D(J)  for I < J.
  //C
  //C  Z      (input) DOUBLE PRECISION array, dimension ( N )
  //C         The components of the updating vector.
  //C
  //C  DELTA  (output) DOUBLE PRECISION array, dimension ( N )
  //C         If N .ne. 1, DELTA contains (D(j) - sigma_I) in its  j-th
  //C         component.  If N = 1, then DELTA(1) = 1.  The vector DELTA
  //C         contains the information necessary to construct the
  //C         (singular) eigenvectors.
  //C
  //C  RHO    (input) DOUBLE PRECISION
  //C         The scalar in the symmetric updating formula.
  //C
  //C  SIGMA  (output) DOUBLE PRECISION
  //C         The computed sigma_I, the I-th updated eigenvalue.
  //C
  //C  WORK   (workspace) DOUBLE PRECISION array, dimension ( N )
  //C         If N .ne. 1, WORK contains (D(j) + sigma_I) in its  j-th
  //C         component.  If N = 1, then WORK( 1 ) = 1.
  //C
  //C  INFO   (output) INTEGER
  //C         = 0:  successful exit
  //C         > 0:  if INFO = 1, the updating process failed.
  //C
  //C  Internal Parameters
  //C  ===================
  //C
  //C  Logical variable ORGATI (origin-at-i?) is used for distinguishing
  //C  whether D(i) or D(i+1) is treated as the origin.
  //C
  //C            ORGATI = .true.    origin at i
  //C            ORGATI = .false.   origin at i+1
  //C
  //C  Logical variable SWTCH3 (switch-for-3-poles?) is for noting
  //C  if we are working with THREE poles!
  //C
  //C  MAXIT is the maximum number of iterations allowed for each
  //C  eigenvalue.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Based on contributions by
  //C     Ren-Cang Li, Computer Science Division, University of California
  //C     at Berkeley, USA
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Local Arrays ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Since this routine is called in an inner loop, we do no argument
  //C     checking.
  //C
  //C     Quick return for N=1 and 2.
  //C
  info = 0;
  if (n == 1) {
    //C
    //C        Presumably, I=1 upon entry
    //C
    sigma = fem::sqrt(d(1) * d(1) + rho * z(1) * z(1));
    delta(1) = one;
    work(1) = one;
    return;
  }
  if (n == 2) {
    dlasd5(i, d, z, delta, rho, sigma, work);
    return;
  }
  //C
  //C     Compute machine epsilon
  //C
  eps = dlamch(cmn, "Epsilon");
  rhoinv = one / rho;
  //C
  //C     The case I = N
  //C
  if (i == n) {
    //C
    //C        Initialize some basic variables
    //C
    ii = n - 1;
    niter = 1;
    //C
    //C        Calculate initial guess
    //C
    temp = rho / two;
    //C
    //C        If ||Z||_2 is not one, then TEMP should be set to
    //C        RHO * ||Z||_2^2 / TWO
    //C
    temp1 = temp / (d(n) + fem::sqrt(d(n) * d(n) + temp));
    FEM_DO(j, 1, n) {
      work(j) = d(j) + d(n) + temp1;
      delta(j) = (d(j) - d(n)) - temp1;
    }
    //C
    psi = zero;
    {
      int fem_do_last = n - 2;
      FEM_DO(j, 1, fem_do_last) {
        psi += z(j) * z(j) / (delta(j) * work(j));
      }
    }
    //C
    c = rhoinv + psi;
    w = c + z(ii) * z(ii) / (delta(ii) * work(ii)) + z(n) * z(n) / (
      delta(n) * work(n));
    //C
    if (w <= zero) {
      temp1 = fem::sqrt(d(n) * d(n) + rho);
      temp = z(n - 1) * z(n - 1) / ((d(n - 1) + temp1) * (d(n) - d(
        n - 1) + rho / (d(n) + temp1))) + z(n) * z(n) / rho;
      //C
      //C           The following TAU is to approximate
      //C           SIGMA_n^2 - D( N )*D( N )
      //C
      if (c <= temp) {
        tau = rho;
      }
      else {
        delsq = (d(n) - d(n - 1)) * (d(n) + d(n - 1));
        a = -c * delsq + z(n - 1) * z(n - 1) + z(n) * z(n);
        b = z(n) * z(n) * delsq;
        if (a < zero) {
          tau = two * b / (fem::sqrt(a * a + four * b * c) - a);
        }
        else {
          tau = (a + fem::sqrt(a * a + four * b * c)) / (two * c);
        }
      }
      //C
      //C           It can be proved that
      //C               D(N)^2+RHO/2 <= SIGMA_n^2 < D(N)^2+TAU <= D(N)^2+RHO
      //C
    }
    else {
      delsq = (d(n) - d(n - 1)) * (d(n) + d(n - 1));
      a = -c * delsq + z(n - 1) * z(n - 1) + z(n) * z(n);
      b = z(n) * z(n) * delsq;
      //C
      //C           The following TAU is to approximate
      //C           SIGMA_n^2 - D( N )*D( N )
      //C
      if (a < zero) {
        tau = two * b / (fem::sqrt(a * a + four * b * c) - a);
      }
      else {
        tau = (a + fem::sqrt(a * a + four * b * c)) / (two * c);
      }
      //C
      //C           It can be proved that
      //C           D(N)^2 < D(N)^2+TAU < SIGMA(N)^2 < D(N)^2+RHO/2
      //C
    }
    //C
    //C        The following ETA is to approximate SIGMA_n - D( N )
    //C
    eta = tau / (d(n) + fem::sqrt(d(n) * d(n) + tau));
    //C
    sigma = d(n) + eta;
    FEM_DO(j, 1, n) {
      delta(j) = (d(j) - d(i)) - eta;
      work(j) = d(j) + d(i) + eta;
    }
    //C
    //C        Evaluate PSI and the derivative DPSI
    //C
    dpsi = zero;
    psi = zero;
    erretm = zero;
    FEM_DO(j, 1, ii) {
      temp = z(j) / (delta(j) * work(j));
      psi += z(j) * temp;
      dpsi += temp * temp;
      erretm += psi;
    }
    erretm = fem::abs(erretm);
    //C
    //C        Evaluate PHI and the derivative DPHI
    //C
    temp = z(n) / (delta(n) * work(n));
    phi = z(n) * temp;
    dphi = temp * temp;
    erretm = eight * (-phi - psi) + erretm - phi + rhoinv + fem::abs(
      tau) * (dpsi + dphi);
    //C
    w = rhoinv + phi + psi;
    //C
    //C        Test for convergence
    //C
    if (fem::abs(w) <= eps * erretm) {
      goto statement_240;
    }
    //C
    //C        Calculate the new step
    //C
    niter++;
    dtnsq1 = work(n - 1) * delta(n - 1);
    dtnsq = work(n) * delta(n);
    c = w - dtnsq1 * dpsi - dtnsq * dphi;
    a = (dtnsq + dtnsq1) * w - dtnsq * dtnsq1 * (dpsi + dphi);
    b = dtnsq * dtnsq1 * w;
    if (c < zero) {
      c = fem::abs(c);
    }
    if (c == zero) {
      eta = rho - sigma * sigma;
    }
    else if (a >= zero) {
      eta = (a + fem::sqrt(fem::abs(a * a - four * b * c))) / (two * c);
    }
    else {
      eta = two * b / (a - fem::sqrt(fem::abs(a * a - four * b * c)));
    }
    //C
    //C        Note, eta should be positive if w is negative, and
    //C        eta should be negative otherwise. However,
    //C        if for some reason caused by roundoff, eta*w > 0,
    //C        we simply use one Newton step instead. This way
    //C        will guarantee eta*w < 0.
    //C
    if (w * eta > zero) {
      eta = -w / (dpsi + dphi);
    }
    temp = eta - dtnsq;
    if (temp > rho) {
      eta = rho + dtnsq;
    }
    //C
    tau += eta;
    eta = eta / (sigma + fem::sqrt(eta + sigma * sigma));
    FEM_DO(j, 1, n) {
      delta(j) = delta(j) - eta;
      work(j) += eta;
    }
    //C
    sigma += eta;
    //C
    //C        Evaluate PSI and the derivative DPSI
    //C
    dpsi = zero;
    psi = zero;
    erretm = zero;
    FEM_DO(j, 1, ii) {
      temp = z(j) / (work(j) * delta(j));
      psi += z(j) * temp;
      dpsi += temp * temp;
      erretm += psi;
    }
    erretm = fem::abs(erretm);
    //C
    //C        Evaluate PHI and the derivative DPHI
    //C
    temp = z(n) / (work(n) * delta(n));
    phi = z(n) * temp;
    dphi = temp * temp;
    erretm = eight * (-phi - psi) + erretm - phi + rhoinv + fem::abs(
      tau) * (dpsi + dphi);
    //C
    w = rhoinv + phi + psi;
    //C
    //C        Main loop to update the values of the array   DELTA
    //C
    iter = niter + 1;
    //C
    FEM_DO(niter, iter, maxit) {
      //C
      //C           Test for convergence
      //C
      if (fem::abs(w) <= eps * erretm) {
        goto statement_240;
      }
      //C
      //C           Calculate the new step
      //C
      dtnsq1 = work(n - 1) * delta(n - 1);
      dtnsq = work(n) * delta(n);
      c = w - dtnsq1 * dpsi - dtnsq * dphi;
      a = (dtnsq + dtnsq1) * w - dtnsq1 * dtnsq * (dpsi + dphi);
      b = dtnsq1 * dtnsq * w;
      if (a >= zero) {
        eta = (a + fem::sqrt(fem::abs(a * a - four * b * c))) / (two * c);
      }
      else {
        eta = two * b / (a - fem::sqrt(fem::abs(a * a - four * b * c)));
      }
      //C
      //C           Note, eta should be positive if w is negative, and
      //C           eta should be negative otherwise. However,
      //C           if for some reason caused by roundoff, eta*w > 0,
      //C           we simply use one Newton step instead. This way
      //C           will guarantee eta*w < 0.
      //C
      if (w * eta > zero) {
        eta = -w / (dpsi + dphi);
      }
      temp = eta - dtnsq;
      if (temp <= zero) {
        eta = eta / two;
      }
      //C
      tau += eta;
      eta = eta / (sigma + fem::sqrt(eta + sigma * sigma));
      FEM_DO(j, 1, n) {
        delta(j) = delta(j) - eta;
        work(j) += eta;
      }
      //C
      sigma += eta;
      //C
      //C           Evaluate PSI and the derivative DPSI
      //C
      dpsi = zero;
      psi = zero;
      erretm = zero;
      FEM_DO(j, 1, ii) {
        temp = z(j) / (work(j) * delta(j));
        psi += z(j) * temp;
        dpsi += temp * temp;
        erretm += psi;
      }
      erretm = fem::abs(erretm);
      //C
      //C           Evaluate PHI and the derivative DPHI
      //C
      temp = z(n) / (work(n) * delta(n));
      phi = z(n) * temp;
      dphi = temp * temp;
      erretm = eight * (-phi - psi) + erretm - phi + rhoinv +
        fem::abs(tau) * (dpsi + dphi);
      //C
      w = rhoinv + phi + psi;
    }
    //C
    //C        Return with INFO = 1, NITER = MAXIT and not converged
    //C
    info = 1;
    goto statement_240;
    //C
    //C        End for the case I = N
    //C
  }
  else {
    //C
    //C        The case for I < N
    //C
    niter = 1;
    ip1 = i + 1;
    //C
    //C        Calculate initial guess
    //C
    delsq = (d(ip1) - d(i)) * (d(ip1) + d(i));
    delsq2 = delsq / two;
    temp = delsq2 / (d(i) + fem::sqrt(d(i) * d(i) + delsq2));
    FEM_DO(j, 1, n) {
      work(j) = d(j) + d(i) + temp;
      delta(j) = (d(j) - d(i)) - temp;
    }
    //C
    psi = zero;
    {
      int fem_do_last = i - 1;
      FEM_DO(j, 1, fem_do_last) {
        psi += z(j) * z(j) / (work(j) * delta(j));
      }
    }
    //C
    phi = zero;
    FEM_DOSTEP(j, n, i + 2, -1) {
      phi += z(j) * z(j) / (work(j) * delta(j));
    }
    c = rhoinv + psi + phi;
    w = c + z(i) * z(i) / (work(i) * delta(i)) + z(ip1) * z(ip1) / (
      work(ip1) * delta(ip1));
    //C
    if (w > zero) {
      //C
      //C           d(i)^2 < the ith sigma^2 < (d(i)^2+d(i+1)^2)/2
      //C
      //C           We choose d(i) as origin.
      //C
      orgati = true;
      sg2lb = zero;
      sg2ub = delsq2;
      a = c * delsq + z(i) * z(i) + z(ip1) * z(ip1);
      b = z(i) * z(i) * delsq;
      if (a > zero) {
        tau = two * b / (a + fem::sqrt(fem::abs(a * a - four * b * c)));
      }
      else {
        tau = (a - fem::sqrt(fem::abs(a * a - four * b * c))) / (two * c);
      }
      //C
      //C           TAU now is an estimation of SIGMA^2 - D( I )^2. The
      //C           following, however, is the corresponding estimation of
      //C           SIGMA - D( I ).
      //C
      eta = tau / (d(i) + fem::sqrt(d(i) * d(i) + tau));
    }
    else {
      //C
      //C           (d(i)^2+d(i+1)^2)/2 <= the ith sigma^2 < d(i+1)^2/2
      //C
      //C           We choose d(i+1) as origin.
      //C
      orgati = false;
      sg2lb = -delsq2;
      sg2ub = zero;
      a = c * delsq - z(i) * z(i) - z(ip1) * z(ip1);
      b = z(ip1) * z(ip1) * delsq;
      if (a < zero) {
        tau = two * b / (a - fem::sqrt(fem::abs(a * a + four * b * c)));
      }
      else {
        tau = -(a + fem::sqrt(fem::abs(a * a + four * b * c))) / (two * c);
      }
      //C
      //C           TAU now is an estimation of SIGMA^2 - D( IP1 )^2. The
      //C           following, however, is the corresponding estimation of
      //C           SIGMA - D( IP1 ).
      //C
      eta = tau / (d(ip1) + fem::sqrt(fem::abs(d(ip1) * d(ip1) + tau)));
    }
    //C
    if (orgati) {
      ii = i;
      sigma = d(i) + eta;
      FEM_DO(j, 1, n) {
        work(j) = d(j) + d(i) + eta;
        delta(j) = (d(j) - d(i)) - eta;
      }
    }
    else {
      ii = i + 1;
      sigma = d(ip1) + eta;
      FEM_DO(j, 1, n) {
        work(j) = d(j) + d(ip1) + eta;
        delta(j) = (d(j) - d(ip1)) - eta;
      }
    }
    iim1 = ii - 1;
    iip1 = ii + 1;
    //C
    //C        Evaluate PSI and the derivative DPSI
    //C
    dpsi = zero;
    psi = zero;
    erretm = zero;
    FEM_DO(j, 1, iim1) {
      temp = z(j) / (work(j) * delta(j));
      psi += z(j) * temp;
      dpsi += temp * temp;
      erretm += psi;
    }
    erretm = fem::abs(erretm);
    //C
    //C        Evaluate PHI and the derivative DPHI
    //C
    dphi = zero;
    phi = zero;
    FEM_DOSTEP(j, n, iip1, -1) {
      temp = z(j) / (work(j) * delta(j));
      phi += z(j) * temp;
      dphi += temp * temp;
      erretm += phi;
    }
    //C
    w = rhoinv + phi + psi;
    //C
    //C        W is the value of the secular function with
    //C        its ii-th element removed.
    //C
    swtch3 = false;
    if (orgati) {
      if (w < zero) {
        swtch3 = true;
      }
    }
    else {
      if (w > zero) {
        swtch3 = true;
      }
    }
    if (ii == 1 || ii == n) {
      swtch3 = false;
    }
    //C
    temp = z(ii) / (work(ii) * delta(ii));
    dw = dpsi + dphi + temp * temp;
    temp = z(ii) * temp;
    w += temp;
    erretm = eight * (phi - psi) + erretm + two * rhoinv + three *
      fem::abs(temp) + fem::abs(tau) * dw;
    //C
    //C        Test for convergence
    //C
    if (fem::abs(w) <= eps * erretm) {
      goto statement_240;
    }
    //C
    if (w <= zero) {
      sg2lb = fem::max(sg2lb, tau);
    }
    else {
      sg2ub = fem::min(sg2ub, tau);
    }
    //C
    //C        Calculate the new step
    //C
    niter++;
    if (!swtch3) {
      dtipsq = work(ip1) * delta(ip1);
      dtisq = work(i) * delta(i);
      if (orgati) {
        c = w - dtipsq * dw + delsq * fem::pow2((z(i) / dtisq));
      }
      else {
        c = w - dtisq * dw - delsq * fem::pow2((z(ip1) / dtipsq));
      }
      a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
      b = dtipsq * dtisq * w;
      if (c == zero) {
        if (a == zero) {
          if (orgati) {
            a = z(i) * z(i) + dtipsq * dtipsq * (dpsi + dphi);
          }
          else {
            a = z(ip1) * z(ip1) + dtisq * dtisq * (dpsi + dphi);
          }
        }
        eta = b / a;
      }
      else if (a <= zero) {
        eta = (a - fem::sqrt(fem::abs(a * a - four * b * c))) / (two * c);
      }
      else {
        eta = two * b / (a + fem::sqrt(fem::abs(a * a - four * b * c)));
      }
    }
    else {
      //C
      //C           Interpolation using THREE most relevant poles
      //C
      dtiim = work(iim1) * delta(iim1);
      dtiip = work(iip1) * delta(iip1);
      temp = rhoinv + psi + phi;
      if (orgati) {
        temp1 = z(iim1) / dtiim;
        temp1 = temp1 * temp1;
        c = (temp - dtiip * (dpsi + dphi)) - (d(iim1) - d(iip1)) * (d(
          iim1) + d(iip1)) * temp1;
        zz(1) = z(iim1) * z(iim1);
        if (dpsi < temp1) {
          zz(3) = dtiip * dtiip * dphi;
        }
        else {
          zz(3) = dtiip * dtiip * ((dpsi - temp1) + dphi);
        }
      }
      else {
        temp1 = z(iip1) / dtiip;
        temp1 = temp1 * temp1;
        c = (temp - dtiim * (dpsi + dphi)) - (d(iip1) - d(iim1)) * (d(
          iim1) + d(iip1)) * temp1;
        if (dphi < temp1) {
          zz(1) = dtiim * dtiim * dpsi;
        }
        else {
          zz(1) = dtiim * dtiim * (dpsi + (dphi - temp1));
        }
        zz(3) = z(iip1) * z(iip1);
      }
      zz(2) = z(ii) * z(ii);
      dd(1) = dtiim;
      dd(2) = delta(ii) * work(ii);
      dd(3) = dtiip;
      dlaed6(cmn, niter, orgati, c, dd, zz, w, eta, info);
      if (info != 0) {
        goto statement_240;
      }
    }
    //C
    //C        Note, eta should be positive if w is negative, and
    //C        eta should be negative otherwise. However,
    //C        if for some reason caused by roundoff, eta*w > 0,
    //C        we simply use one Newton step instead. This way
    //C        will guarantee eta*w < 0.
    //C
    if (w * eta >= zero) {
      eta = -w / dw;
    }
    if (orgati) {
      temp1 = work(i) * delta(i);
      temp = eta - temp1;
    }
    else {
      temp1 = work(ip1) * delta(ip1);
      temp = eta - temp1;
    }
    if (temp > sg2ub || temp < sg2lb) {
      if (w < zero) {
        eta = (sg2ub - tau) / two;
      }
      else {
        eta = (sg2lb - tau) / two;
      }
    }
    //C
    tau += eta;
    eta = eta / (sigma + fem::sqrt(sigma * sigma + eta));
    //C
    prew = w;
    //C
    sigma += eta;
    FEM_DO(j, 1, n) {
      work(j) += eta;
      delta(j) = delta(j) - eta;
    }
    //C
    //C        Evaluate PSI and the derivative DPSI
    //C
    dpsi = zero;
    psi = zero;
    erretm = zero;
    FEM_DO(j, 1, iim1) {
      temp = z(j) / (work(j) * delta(j));
      psi += z(j) * temp;
      dpsi += temp * temp;
      erretm += psi;
    }
    erretm = fem::abs(erretm);
    //C
    //C        Evaluate PHI and the derivative DPHI
    //C
    dphi = zero;
    phi = zero;
    FEM_DOSTEP(j, n, iip1, -1) {
      temp = z(j) / (work(j) * delta(j));
      phi += z(j) * temp;
      dphi += temp * temp;
      erretm += phi;
    }
    //C
    temp = z(ii) / (work(ii) * delta(ii));
    dw = dpsi + dphi + temp * temp;
    temp = z(ii) * temp;
    w = rhoinv + phi + psi + temp;
    erretm = eight * (phi - psi) + erretm + two * rhoinv + three *
      fem::abs(temp) + fem::abs(tau) * dw;
    //C
    if (w <= zero) {
      sg2lb = fem::max(sg2lb, tau);
    }
    else {
      sg2ub = fem::min(sg2ub, tau);
    }
    //C
    swtch = false;
    if (orgati) {
      if (-w > fem::abs(prew) / ten) {
        swtch = true;
      }
    }
    else {
      if (w > fem::abs(prew) / ten) {
        swtch = true;
      }
    }
    //C
    //C        Main loop to update the values of the array   DELTA and WORK
    //C
    iter = niter + 1;
    //C
    FEM_DO(niter, iter, maxit) {
      //C
      //C           Test for convergence
      //C
      if (fem::abs(w) <= eps * erretm) {
        goto statement_240;
      }
      //C
      //C           Calculate the new step
      //C
      if (!swtch3) {
        dtipsq = work(ip1) * delta(ip1);
        dtisq = work(i) * delta(i);
        if (!swtch) {
          if (orgati) {
            c = w - dtipsq * dw + delsq * fem::pow2((z(i) / dtisq));
          }
          else {
            c = w - dtisq * dw - delsq * fem::pow2((z(ip1) / dtipsq));
          }
        }
        else {
          temp = z(ii) / (work(ii) * delta(ii));
          if (orgati) {
            dpsi += temp * temp;
          }
          else {
            dphi += temp * temp;
          }
          c = w - dtisq * dpsi - dtipsq * dphi;
        }
        a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
        b = dtipsq * dtisq * w;
        if (c == zero) {
          if (a == zero) {
            if (!swtch) {
              if (orgati) {
                a = z(i) * z(i) + dtipsq * dtipsq * (dpsi + dphi);
              }
              else {
                a = z(ip1) * z(ip1) + dtisq * dtisq * (dpsi + dphi);
              }
            }
            else {
              a = dtisq * dtisq * dpsi + dtipsq * dtipsq * dphi;
            }
          }
          eta = b / a;
        }
        else if (a <= zero) {
          eta = (a - fem::sqrt(fem::abs(a * a - four * b * c))) / (two * c);
        }
        else {
          eta = two * b / (a + fem::sqrt(fem::abs(a * a - four * b * c)));
        }
      }
      else {
        //C
        //C              Interpolation using THREE most relevant poles
        //C
        dtiim = work(iim1) * delta(iim1);
        dtiip = work(iip1) * delta(iip1);
        temp = rhoinv + psi + phi;
        if (swtch) {
          c = temp - dtiim * dpsi - dtiip * dphi;
          zz(1) = dtiim * dtiim * dpsi;
          zz(3) = dtiip * dtiip * dphi;
        }
        else {
          if (orgati) {
            temp1 = z(iim1) / dtiim;
            temp1 = temp1 * temp1;
            temp2 = (d(iim1) - d(iip1)) * (d(iim1) + d(iip1)) * temp1;
            c = temp - dtiip * (dpsi + dphi) - temp2;
            zz(1) = z(iim1) * z(iim1);
            if (dpsi < temp1) {
              zz(3) = dtiip * dtiip * dphi;
            }
            else {
              zz(3) = dtiip * dtiip * ((dpsi - temp1) + dphi);
            }
          }
          else {
            temp1 = z(iip1) / dtiip;
            temp1 = temp1 * temp1;
            temp2 = (d(iip1) - d(iim1)) * (d(iim1) + d(iip1)) * temp1;
            c = temp - dtiim * (dpsi + dphi) - temp2;
            if (dphi < temp1) {
              zz(1) = dtiim * dtiim * dpsi;
            }
            else {
              zz(1) = dtiim * dtiim * (dpsi + (dphi - temp1));
            }
            zz(3) = z(iip1) * z(iip1);
          }
        }
        dd(1) = dtiim;
        dd(2) = delta(ii) * work(ii);
        dd(3) = dtiip;
        dlaed6(cmn, niter, orgati, c, dd, zz, w, eta, info);
        if (info != 0) {
          goto statement_240;
        }
      }
      //C
      //C           Note, eta should be positive if w is negative, and
      //C           eta should be negative otherwise. However,
      //C           if for some reason caused by roundoff, eta*w > 0,
      //C           we simply use one Newton step instead. This way
      //C           will guarantee eta*w < 0.
      //C
      if (w * eta >= zero) {
        eta = -w / dw;
      }
      if (orgati) {
        temp1 = work(i) * delta(i);
        temp = eta - temp1;
      }
      else {
        temp1 = work(ip1) * delta(ip1);
        temp = eta - temp1;
      }
      if (temp > sg2ub || temp < sg2lb) {
        if (w < zero) {
          eta = (sg2ub - tau) / two;
        }
        else {
          eta = (sg2lb - tau) / two;
        }
      }
      //C
      tau += eta;
      eta = eta / (sigma + fem::sqrt(sigma * sigma + eta));
      //C
      sigma += eta;
      FEM_DO(j, 1, n) {
        work(j) += eta;
        delta(j) = delta(j) - eta;
      }
      //C
      prew = w;
      //C
      //C           Evaluate PSI and the derivative DPSI
      //C
      dpsi = zero;
      psi = zero;
      erretm = zero;
      FEM_DO(j, 1, iim1) {
        temp = z(j) / (work(j) * delta(j));
        psi += z(j) * temp;
        dpsi += temp * temp;
        erretm += psi;
      }
      erretm = fem::abs(erretm);
      //C
      //C           Evaluate PHI and the derivative DPHI
      //C
      dphi = zero;
      phi = zero;
      FEM_DOSTEP(j, n, iip1, -1) {
        temp = z(j) / (work(j) * delta(j));
        phi += z(j) * temp;
        dphi += temp * temp;
        erretm += phi;
      }
      //C
      temp = z(ii) / (work(ii) * delta(ii));
      dw = dpsi + dphi + temp * temp;
      temp = z(ii) * temp;
      w = rhoinv + phi + psi + temp;
      erretm = eight * (phi - psi) + erretm + two * rhoinv + three *
        fem::abs(temp) + fem::abs(tau) * dw;
      if (w * prew > zero && fem::abs(w) > fem::abs(prew) / ten) {
        swtch = !swtch;
      }
      //C
      if (w <= zero) {
        sg2lb = fem::max(sg2lb, tau);
      }
      else {
        sg2ub = fem::min(sg2ub, tau);
      }
      //C
    }
    //C
    //C        Return with INFO = 1, NITER = MAXIT and not converged
    //C
    info = 1;
    //C
  }
  //C
  statement_240:;
  //C
  //C     End of DLASD4
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasd3.f
inline
void
dlasd3(
  common& cmn,
  int const& nl,
  int const& nr,
  int const& sqre,
  int const& k,
  arr_ref<double> d,
  arr_ref<double, 2> q,
  int const& ldq,
  arr_ref<double> dsigma,
  arr_ref<double, 2> u,
  int const& ldu,
  arr_cref<double, 2> u2,
  int const& ldu2,
  arr_ref<double, 2> vt,
  int const& ldvt,
  arr_ref<double, 2> vt2,
  int const& ldvt2,
  arr_cref<int> idxc,
  arr_cref<int> ctot,
  arr_ref<double> z,
  int& info)
{
  d(dimension(star));
  q(dimension(ldq, star));
  dsigma(dimension(star));
  u(dimension(ldu, star));
  u2(dimension(ldu2, star));
  vt(dimension(ldvt, star));
  vt2(dimension(ldvt2, star));
  idxc(dimension(star));
  ctot(dimension(star));
  z(dimension(star));
  int n = fem::int0;
  int m = fem::int0;
  int nlp1 = fem::int0;
  int nlp2 = fem::int0;
  const double zero = 0.0e+0;
  int i = fem::int0;
  double rho = fem::double0;
  const double one = 1.0e+0;
  int j = fem::int0;
  const double negone = -1.0e+0;
  double temp = fem::double0;
  int jc = fem::int0;
  int ktemp = fem::int0;
  int ctemp = fem::int0;
  int nrp1 = fem::int0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASD3 finds all the square roots of the roots of the secular
  //C  equation, as defined by the values in D and Z.  It makes the
  //C  appropriate calls to DLASD4 and then updates the singular
  //C  vectors by matrix multiplication.
  //C
  //C  This code makes very mild assumptions about floating point
  //C  arithmetic. It will work on machines with a guard digit in
  //C  add/subtract, or on those binary machines without guard digits
  //C  which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.
  //C  It could conceivably fail on hexadecimal or decimal machines
  //C  without guard digits, but we know of none.
  //C
  //C  DLASD3 is called from DLASD1.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  NL     (input) INTEGER
  //C         The row dimension of the upper block.  NL >= 1.
  //C
  //C  NR     (input) INTEGER
  //C         The row dimension of the lower block.  NR >= 1.
  //C
  //C  SQRE   (input) INTEGER
  //C         = 0: the lower block is an NR-by-NR square matrix.
  //C         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
  //C
  //C         The bidiagonal matrix has N = NL + NR + 1 rows and
  //C         M = N + SQRE >= N columns.
  //C
  //C  K      (input) INTEGER
  //C         The size of the secular equation, 1 =< K = < N.
  //C
  //C  D      (output) DOUBLE PRECISION array, dimension(K)
  //C         On exit the square roots of the roots of the secular equation,
  //C         in ascending order.
  //C
  //C  Q      (workspace) DOUBLE PRECISION array,
  //C                     dimension at least (LDQ,K).
  //C
  //C  LDQ    (input) INTEGER
  //C         The leading dimension of the array Q.  LDQ >= K.
  //C
  //C  DSIGMA (input) DOUBLE PRECISION array, dimension(K)
  //C         The first K elements of this array contain the old roots
  //C         of the deflated updating problem.  These are the poles
  //C         of the secular equation.
  //C
  //C  U      (output) DOUBLE PRECISION array, dimension (LDU, N)
  //C         The last N - K columns of this matrix contain the deflated
  //C         left singular vectors.
  //C
  //C  LDU    (input) INTEGER
  //C         The leading dimension of the array U.  LDU >= N.
  //C
  //C  U2     (input/output) DOUBLE PRECISION array, dimension (LDU2, N)
  //C         The first K columns of this matrix contain the non-deflated
  //C         left singular vectors for the split problem.
  //C
  //C  LDU2   (input) INTEGER
  //C         The leading dimension of the array U2.  LDU2 >= N.
  //C
  //C  VT     (output) DOUBLE PRECISION array, dimension (LDVT, M)
  //C         The last M - K columns of VT' contain the deflated
  //C         right singular vectors.
  //C
  //C  LDVT   (input) INTEGER
  //C         The leading dimension of the array VT.  LDVT >= N.
  //C
  //C  VT2    (input/output) DOUBLE PRECISION array, dimension (LDVT2, N)
  //C         The first K columns of VT2' contain the non-deflated
  //C         right singular vectors for the split problem.
  //C
  //C  LDVT2  (input) INTEGER
  //C         The leading dimension of the array VT2.  LDVT2 >= N.
  //C
  //C  IDXC   (input) INTEGER array, dimension ( N )
  //C         The permutation used to arrange the columns of U (and rows of
  //C         VT) into three groups:  the first group contains non-zero
  //C         entries only at and above (or before) NL +1; the second
  //C         contains non-zero entries only at and below (or after) NL+2;
  //C         and the third is dense. The first column of U and the row of
  //C         VT are treated separately, however.
  //C
  //C         The rows of the singular vectors found by DLASD4
  //C         must be likewise permuted before the matrix multiplies can
  //C         take place.
  //C
  //C  CTOT   (input) INTEGER array, dimension ( 4 )
  //C         A count of the total number of the various types of columns
  //C         in U (or rows in VT), as described in IDXC. The fourth column
  //C         type is any column which has been deflated.
  //C
  //C  Z      (input) DOUBLE PRECISION array, dimension (K)
  //C         The first K elements of this array contain the components
  //C         of the deflation-adjusted updating row vector.
  //C
  //C  INFO   (output) INTEGER
  //C         = 0:  successful exit.
  //C         < 0:  if INFO = -i, the i-th argument had an illegal value.
  //C         > 0:  if INFO = 1, an singular value did not converge
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Based on contributions by
  //C     Ming Gu and Huan Ren, Computer Science Division, University of
  //C     California at Berkeley, USA
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters.
  //C
  info = 0;
  //C
  if (nl < 1) {
    info = -1;
  }
  else if (nr < 1) {
    info = -2;
  }
  else if ((sqre != 1) && (sqre != 0)) {
    info = -3;
  }
  //C
  n = nl + nr + 1;
  m = n + sqre;
  nlp1 = nl + 1;
  nlp2 = nl + 2;
  //C
  if ((k < 1) || (k > n)) {
    info = -4;
  }
  else if (ldq < k) {
    info = -7;
  }
  else if (ldu < n) {
    info = -10;
  }
  else if (ldu2 < n) {
    info = -12;
  }
  else if (ldvt < m) {
    info = -14;
  }
  else if (ldvt2 < m) {
    info = -16;
  }
  if (info != 0) {
    xerbla("DLASD3", -info);
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (k == 1) {
    d(1) = fem::abs(z(1));
    dcopy(m, vt2(1, 1), ldvt2, vt(1, 1), ldvt);
    if (z(1) > zero) {
      dcopy(n, u2(1, 1), 1, u(1, 1), 1);
    }
    else {
      FEM_DO(i, 1, n) {
        u(i, 1) = -u2(i, 1);
      }
    }
    return;
  }
  //C
  //C     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
  //C     be computed with high relative accuracy (barring over/underflow).
  //C     This is a problem on machines without a guard digit in
  //C     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
  //C     The following code replaces DSIGMA(I) by 2*DSIGMA(I)-DSIGMA(I),
  //C     which on any of these machines zeros out the bottommost
  //C     bit of DSIGMA(I) if it is 1; this makes the subsequent
  //C     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
  //C     occurs. On binary machines with a guard digit (almost all
  //C     machines) it does not change DSIGMA(I) at all. On hexadecimal
  //C     and decimal machines with a guard digit, it slightly
  //C     changes the bottommost bits of DSIGMA(I). It does not account
  //C     for hexadecimal or decimal machines without guard digits
  //C     (we know of none). We use a subroutine call to compute
  //C     2*DSIGMA(I) to prevent optimizing compilers from eliminating
  //C     this code.
  //C
  FEM_DO(i, 1, k) {
    dsigma(i) = dlamc3(dsigma(i), dsigma(i)) - dsigma(i);
  }
  //C
  //C     Keep a copy of Z.
  //C
  dcopy(k, z, 1, q, 1);
  //C
  //C     Normalize Z.
  //C
  rho = dnrm2(k, z, 1);
  dlascl(cmn, "G", 0, 0, rho, one, k, 1, z, k, info);
  rho = rho * rho;
  //C
  //C     Find the new singular values.
  //C
  FEM_DO(j, 1, k) {
    dlasd4(cmn, k, j, dsigma, z, u(1, j), rho, d(j), vt(1, j), info);
    //C
    //C        If the zero finder fails, the computation is terminated.
    //C
    if (info != 0) {
      return;
    }
  }
  //C
  //C     Compute updated Z.
  //C
  FEM_DO(i, 1, k) {
    z(i) = u(i, k) * vt(i, k);
    {
      int fem_do_last = i - 1;
      FEM_DO(j, 1, fem_do_last) {
        z(i) = z(i) * (u(i, j) * vt(i, j) / (dsigma(i) - dsigma(j)) /
          (dsigma(i) + dsigma(j)));
      }
    }
    {
      int fem_do_last = k - 1;
      FEM_DO(j, i, fem_do_last) {
        z(i) = z(i) * (u(i, j) * vt(i, j) / (dsigma(i) - dsigma(j +
          1)) / (dsigma(i) + dsigma(j + 1)));
      }
    }
    z(i) = fem::sign(fem::sqrt(fem::abs(z(i))), q(i, 1));
  }
  //C
  //C     Compute left singular vectors of the modified diagonal matrix,
  //C     and store related information for the right singular vectors.
  //C
  FEM_DO(i, 1, k) {
    vt(1, i) = z(1) / u(1, i) / vt(1, i);
    u(1, i) = negone;
    FEM_DO(j, 2, k) {
      vt(j, i) = z(j) / u(j, i) / vt(j, i);
      u(j, i) = dsigma(j) * vt(j, i);
    }
    temp = dnrm2(k, u(1, i), 1);
    q(1, i) = u(1, i) / temp;
    FEM_DO(j, 2, k) {
      jc = idxc(j);
      q(j, i) = u(jc, i) / temp;
    }
  }
  //C
  //C     Update the left singular vector matrix.
  //C
  if (k == 2) {
    dgemm("N", "N", n, k, k, one, u2, ldu2, q, ldq, zero, u, ldu);
    goto statement_100;
  }
  if (ctot(1) > 0) {
    dgemm("N", "N", nl, k, ctot(1), one, u2(1, 2), ldu2, q(2, 1),
      ldq, zero, u(1, 1), ldu);
    if (ctot(3) > 0) {
      ktemp = 2 + ctot(1) + ctot(2);
      dgemm("N", "N", nl, k, ctot(3), one, u2(1, ktemp), ldu2, q(ktemp,
        1), ldq, one, u(1, 1), ldu);
    }
  }
  else if (ctot(3) > 0) {
    ktemp = 2 + ctot(1) + ctot(2);
    dgemm("N", "N", nl, k, ctot(3), one, u2(1, ktemp), ldu2, q(ktemp,
      1), ldq, zero, u(1, 1), ldu);
  }
  else {
    dlacpy("F", nl, k, u2, ldu2, u, ldu);
  }
  dcopy(k, q(1, 1), ldq, u(nlp1, 1), ldu);
  ktemp = 2 + ctot(1);
  ctemp = ctot(2) + ctot(3);
  dgemm("N", "N", nr, k, ctemp, one, u2(nlp2, ktemp), ldu2, q(ktemp,
    1), ldq, zero, u(nlp2, 1), ldu);
  //C
  //C     Generate the right singular vectors.
  //C
  statement_100:
  FEM_DO(i, 1, k) {
    temp = dnrm2(k, vt(1, i), 1);
    q(i, 1) = vt(1, i) / temp;
    FEM_DO(j, 2, k) {
      jc = idxc(j);
      q(i, j) = vt(jc, i) / temp;
    }
  }
  //C
  //C     Update the right singular vector matrix.
  //C
  if (k == 2) {
    dgemm("N", "N", k, m, k, one, q, ldq, vt2, ldvt2, zero, vt, ldvt);
    return;
  }
  ktemp = 1 + ctot(1);
  dgemm("N", "N", k, nlp1, ktemp, one, q(1, 1), ldq, vt2(1, 1),
    ldvt2, zero, vt(1, 1), ldvt);
  ktemp = 2 + ctot(1) + ctot(2);
  if (ktemp <= ldvt2) {
    dgemm("N", "N", k, nlp1, ctot(3), one, q(1, ktemp), ldq, vt2(ktemp,
      1), ldvt2, one, vt(1, 1), ldvt);
  }
  //C
  ktemp = ctot(1) + 1;
  nrp1 = nr + sqre;
  if (ktemp > 1) {
    FEM_DO(i, 1, k) {
      q(i, ktemp) = q(i, 1);
    }
    FEM_DO(i, nlp2, m) {
      vt2(ktemp, i) = vt2(1, i);
    }
  }
  ctemp = 1 + ctot(2) + ctot(3);
  dgemm("N", "N", k, nrp1, ctemp, one, q(1, ktemp), ldq, vt2(ktemp,
    nlp2), ldvt2, zero, vt(1, nlp2), ldvt);
  //C
  //C     End of DLASD3
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasd1.f
inline
void
dlasd1(
  common& cmn,
  int const& nl,
  int const& nr,
  int const& sqre,
  arr_ref<double> d,
  double& alpha,
  double& beta,
  arr_ref<double, 2> u,
  int const& ldu,
  arr_ref<double, 2> vt,
  int const& ldvt,
  arr_ref<int> idxq,
  arr_ref<int> iwork,
  arr_ref<double> work,
  int& info)
{
  d(dimension(star));
  u(dimension(ldu, star));
  vt(dimension(ldvt, star));
  idxq(dimension(star));
  iwork(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASD1 computes the SVD of an upper bidiagonal N-by-M matrix B,
  //C  where N = NL + NR + 1 and M = N + SQRE. DLASD1 is called from DLASD0.
  //C
  //C  A related subroutine DLASD7 handles the case in which the singular
  //C  values (and the singular vectors in factored form) are desired.
  //C
  //C  DLASD1 computes the SVD as follows:
  //C
  //C                ( D1(in)  0    0     0 )
  //C    B = U(in) * (   Z1'   a   Z2'    b ) * VT(in)
  //C                (   0     0   D2(in) 0 )
  //C
  //C      = U(out) * ( D(out) 0) * VT(out)
  //C
  //C  where Z' = (Z1' a Z2' b) = u' VT', and u is a vector of dimension M
  //C  with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros
  //C  elsewhere; and the entry b is empty if SQRE = 0.
  //C
  //C  The left singular vectors of the original matrix are stored in U, and
  //C  the transpose of the right singular vectors are stored in VT, and the
  //C  singular values are in D.  The algorithm consists of three stages:
  //C
  //C     The first stage consists of deflating the size of the problem
  //C     when there are multiple singular values or when there are zeros in
  //C     the Z vector.  For each such occurence the dimension of the
  //C     secular equation problem is reduced by one.  This stage is
  //C     performed by the routine DLASD2.
  //C
  //C     The second stage consists of calculating the updated
  //C     singular values. This is done by finding the square roots of the
  //C     roots of the secular equation via the routine DLASD4 (as called
  //C     by DLASD3). This routine also calculates the singular vectors of
  //C     the current problem.
  //C
  //C     The final stage consists of computing the updated singular vectors
  //C     directly using the updated singular values.  The singular vectors
  //C     for the current problem are multiplied with the singular vectors
  //C     from the overall problem.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  NL     (input) INTEGER
  //C         The row dimension of the upper block.  NL >= 1.
  //C
  //C  NR     (input) INTEGER
  //C         The row dimension of the lower block.  NR >= 1.
  //C
  //C  SQRE   (input) INTEGER
  //C         = 0: the lower block is an NR-by-NR square matrix.
  //C         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
  //C
  //C         The bidiagonal matrix has row dimension N = NL + NR + 1,
  //C         and column dimension M = N + SQRE.
  //C
  //C  D      (input/output) DOUBLE PRECISION array,
  //C                        dimension (N = NL+NR+1).
  //C         On entry D(1:NL,1:NL) contains the singular values of the
  //C         upper block; and D(NL+2:N) contains the singular values of
  //C         the lower block. On exit D(1:N) contains the singular values
  //C         of the modified matrix.
  //C
  //C  ALPHA  (input/output) DOUBLE PRECISION
  //C         Contains the diagonal element associated with the added row.
  //C
  //C  BETA   (input/output) DOUBLE PRECISION
  //C         Contains the off-diagonal element associated with the added
  //C         row.
  //C
  //C  U      (input/output) DOUBLE PRECISION array, dimension(LDU,N)
  //C         On entry U(1:NL, 1:NL) contains the left singular vectors of
  //C         the upper block; U(NL+2:N, NL+2:N) contains the left singular
  //C         vectors of the lower block. On exit U contains the left
  //C         singular vectors of the bidiagonal matrix.
  //C
  //C  LDU    (input) INTEGER
  //C         The leading dimension of the array U.  LDU >= max( 1, N ).
  //C
  //C  VT     (input/output) DOUBLE PRECISION array, dimension(LDVT,M)
  //C         where M = N + SQRE.
  //C         On entry VT(1:NL+1, 1:NL+1)' contains the right singular
  //C         vectors of the upper block; VT(NL+2:M, NL+2:M)' contains
  //C         the right singular vectors of the lower block. On exit
  //C         VT' contains the right singular vectors of the
  //C         bidiagonal matrix.
  //C
  //C  LDVT   (input) INTEGER
  //C         The leading dimension of the array VT.  LDVT >= max( 1, M ).
  //C
  //C  IDXQ  (output) INTEGER array, dimension(N)
  //C         This contains the permutation which will reintegrate the
  //C         subproblem just solved back into sorted order, i.e.
  //C         D( IDXQ( I = 1, N ) ) will be in ascending order.
  //C
  //C  IWORK  (workspace) INTEGER array, dimension( 4 * N )
  //C
  //C  WORK   (workspace) DOUBLE PRECISION array, dimension( 3*M**2 + 2*M )
  //C
  //C  INFO   (output) INTEGER
  //C          = 0:  successful exit.
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value.
  //C          > 0:  if INFO = 1, an singular value did not converge
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Based on contributions by
  //C     Ming Gu and Huan Ren, Computer Science Division, University of
  //C     California at Berkeley, USA
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters.
  //C
  info = 0;
  //C
  if (nl < 1) {
    info = -1;
  }
  else if (nr < 1) {
    info = -2;
  }
  else if ((sqre < 0) || (sqre > 1)) {
    info = -3;
  }
  if (info != 0) {
    xerbla("DLASD1", -info);
    return;
  }
  //C
  int n = nl + nr + 1;
  int m = n + sqre;
  //C
  //C     The following values are for bookkeeping purposes only.  They are
  //C     integer pointers which indicate the portion of the workspace
  //C     used by a particular array in DLASD2 and DLASD3.
  //C
  int ldu2 = n;
  int ldvt2 = m;
  //C
  int iz = 1;
  int isigma = iz + m;
  int iu2 = isigma + n;
  int ivt2 = iu2 + ldu2 * n;
  int iq = ivt2 + ldvt2 * m;
  //C
  int idx = 1;
  int idxc = idx + n;
  int coltyp = idxc + n;
  int idxp = coltyp + n;
  //C
  //C     Scale.
  //C
  double orgnrm = fem::max(fem::abs(alpha), fem::abs(beta));
  const double zero = 0.0e+0;
  d(nl + 1) = zero;
  int i = fem::int0;
  FEM_DO(i, 1, n) {
    if (fem::abs(d(i)) > orgnrm) {
      orgnrm = fem::abs(d(i));
    }
  }
  const double one = 1.0e+0;
  dlascl(cmn, "G", 0, 0, orgnrm, one, n, 1, d, n, info);
  alpha = alpha / orgnrm;
  beta = beta / orgnrm;
  //C
  //C     Deflate singular values.
  //C
  int k = fem::int0;
  dlasd2(cmn, nl, nr, sqre, k, d, work(iz), alpha, beta, u, ldu, vt,
    ldvt, work(isigma), work(iu2), ldu2, work(ivt2), ldvt2, iwork(idxp),
    iwork(idx), iwork(idxc), idxq, iwork(coltyp), info);
  //C
  //C     Solve Secular Equation and update singular vectors.
  //C
  int ldq = k;
  dlasd3(cmn, nl, nr, sqre, k, d, work(iq), ldq, work(isigma), u,
    ldu, work(iu2), ldu2, vt, ldvt, work(ivt2), ldvt2, iwork(idxc),
    iwork(coltyp), work(iz), info);
  if (info != 0) {
    return;
  }
  //C
  //C     Unscale.
  //C
  dlascl(cmn, "G", 0, 0, one, orgnrm, n, 1, d, n, info);
  //C
  //C     Prepare the IDXQ sorting permutation.
  //C
  int n1 = k;
  int n2 = n - k;
  dlamrg(n1, n2, d, 1, -1, idxq);
  //C
  //C     End of DLASD1
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlas2.f
inline
void
dlas2(
  double const& f,
  double const& g,
  double const& h,
  double& ssmin,
  double& ssmax)
{
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLAS2  computes the singular values of the 2-by-2 matrix
  //C     [  F   G  ]
  //C     [  0   H  ].
  //C  On return, SSMIN is the smaller singular value and SSMAX is the
  //C  larger singular value.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  F       (input) DOUBLE PRECISION
  //C          The (1,1) element of the 2-by-2 matrix.
  //C
  //C  G       (input) DOUBLE PRECISION
  //C          The (1,2) element of the 2-by-2 matrix.
  //C
  //C  H       (input) DOUBLE PRECISION
  //C          The (2,2) element of the 2-by-2 matrix.
  //C
  //C  SSMIN   (output) DOUBLE PRECISION
  //C          The smaller singular value.
  //C
  //C  SSMAX   (output) DOUBLE PRECISION
  //C          The larger singular value.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Barring over/underflow, all output quantities are correct to within
  //C  a few units in the last place (ulps), even in the absence of a guard
  //C  digit in addition/subtraction.
  //C
  //C  In IEEE arithmetic, the code works correctly if one matrix element is
  //C  infinite.
  //C
  //C  Overflow will not occur unless the largest singular value itself
  //C  overflows, or is within a few ulps of overflow. (On machines with
  //C  partial overflow, like the Cray, overflow may occur if the largest
  //C  singular value is within a factor of 2 of overflow.)
  //C
  //C  Underflow is harmless if underflow is gradual. Otherwise, results
  //C  may correspond to a matrix modified by perturbations of size near
  //C  the underflow threshold.
  //C
  //C  ====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  double fa = fem::abs(f);
  double ga = fem::abs(g);
  double ha = fem::abs(h);
  double fhmn = fem::min(fa, ha);
  double fhmx = fem::max(fa, ha);
  const double zero = 0.0e0;
  const double one = 1.0e0;
  double as = fem::double0;
  double at = fem::double0;
  double au = fem::double0;
  const double two = 2.0e0;
  double c = fem::double0;
  if (fhmn == zero) {
    ssmin = zero;
    if (fhmx == zero) {
      ssmax = ga;
    }
    else {
      ssmax = fem::max(fhmx, ga) * fem::sqrt(one + fem::pow2((fem::min(fhmx,
        ga) / fem::max(fhmx, ga))));
    }
  }
  else {
    if (ga < fhmx) {
      as = one + fhmn / fhmx;
      at = (fhmx - fhmn) / fhmx;
      au = fem::pow2((ga / fhmx));
      c = two / (fem::sqrt(as * as + au) + fem::sqrt(at * at + au));
      ssmin = fhmn * c;
      ssmax = fhmx / c;
    }
    else {
      au = fhmx / ga;
      if (au == zero) {
        //C
        //C              Avoid possible harmful underflow if exponent range
        //C              asymmetric (true SSMIN may not underflow even if
        //C              AU underflows)
        //C
        ssmin = (fhmn * fhmx) / ga;
        ssmax = ga;
      }
      else {
        as = one + fhmn / fhmx;
        at = (fhmx - fhmn) / fhmx;
        c = one / (fem::sqrt(one + fem::pow2((as * au))) + fem::sqrt(
          one + fem::pow2((at * au))));
        ssmin = (fhmn * c) * au;
        ssmin += ssmin;
        ssmax = ga / (c + c);
      }
    }
  }
  //C
  //C     End of DLAS2
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasq4.f
inline
void
dlasq4(
  int const& i0,
  int const& n0,
  arr_cref<double> z,
  int const& pp,
  int const& n0in,
  double const& dmin,
  double const& dmin1,
  double const& dmin2,
  double const& dn,
  double const& dn1,
  double const& dn2,
  double& tau,
  int& ttype,
  double& g)
{
  z(dimension(star));
  const double zero = 0.0e0;
  int nn = fem::int0;
  double b1 = fem::double0;
  double b2 = fem::double0;
  double a2 = fem::double0;
  const double qurtr = 0.250e0;
  double gap2 = fem::double0;
  double gap1 = fem::double0;
  const double half = 0.50e0;
  double s = fem::double0;
  const double third = 0.3330e0;
  double gam = fem::double0;
  int np = fem::int0;
  int i4 = fem::int0;
  const double hundrd = 100.0e0;
  const double cnst1 = 0.5630e0;
  const double cnst3 = 1.050e0;
  const double one = 1.0e0;
  const double cnst2 = 1.010e0;
  const double two = 2.0e0;
  //C
  //C  -- LAPACK routine (version 3.2)                                    --
  //C
  //C  -- Contributed by Osni Marques of the Lawrence Berkeley National   --
  //C  -- Laboratory and Beresford Parlett of the Univ. of California at  --
  //C  -- Berkeley                                                        --
  //C  -- November 2008                                                   --
  //C
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASQ4 computes an approximation TAU to the smallest eigenvalue
  //C  using values of d from the previous transform.
  //C
  //C  I0    (input) INTEGER
  //C        First index.
  //C
  //C  N0    (input) INTEGER
  //C        Last index.
  //C
  //C  Z     (input) DOUBLE PRECISION array, dimension ( 4*N )
  //C        Z holds the qd array.
  //C
  //C  PP    (input) INTEGER
  //C        PP=0 for ping, PP=1 for pong.
  //C
  //C  NOIN  (input) INTEGER
  //C        The value of N0 at start of EIGTEST.
  //C
  //C  DMIN  (input) DOUBLE PRECISION
  //C        Minimum value of d.
  //C
  //C  DMIN1 (input) DOUBLE PRECISION
  //C        Minimum value of d, excluding D( N0 ).
  //C
  //C  DMIN2 (input) DOUBLE PRECISION
  //C        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
  //C
  //C  DN    (input) DOUBLE PRECISION
  //C        d(N)
  //C
  //C  DN1   (input) DOUBLE PRECISION
  //C        d(N-1)
  //C
  //C  DN2   (input) DOUBLE PRECISION
  //C        d(N-2)
  //C
  //C  TAU   (output) DOUBLE PRECISION
  //C        This is the shift.
  //C
  //C  TTYPE (output) INTEGER
  //C        Shift type.
  //C
  //C  G     (input/output) REAL
  //C        G is passed as an argument in order to save its value between
  //C        calls to DLASQ4.
  //C
  //C  Further Details
  //C  ===============
  //C  CNST1 = 9/16
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     A negative DMIN forces the shift to take that absolute value
  //C     TTYPE records the type of shift.
  //C
  if (dmin <= zero) {
    tau = -dmin;
    ttype = -1;
    return;
  }
  //C
  nn = 4 * n0 + pp;
  if (n0in == n0) {
    //C
    //C        No eigenvalues deflated.
    //C
    if (dmin == dn || dmin == dn1) {
      //C
      b1 = fem::sqrt(z(nn - 3)) * fem::sqrt(z(nn - 5));
      b2 = fem::sqrt(z(nn - 7)) * fem::sqrt(z(nn - 9));
      a2 = z(nn - 7) + z(nn - 5);
      //C
      //C           Cases 2 and 3.
      //C
      if (dmin == dn && dmin1 == dn1) {
        gap2 = dmin2 - a2 - dmin2 * qurtr;
        if (gap2 > zero && gap2 > b2) {
          gap1 = a2 - dn - (b2 / gap2) * b2;
        }
        else {
          gap1 = a2 - dn - (b1 + b2);
        }
        if (gap1 > zero && gap1 > b1) {
          s = fem::max(dn - (b1 / gap1) * b1, half * dmin);
          ttype = -2;
        }
        else {
          s = zero;
          if (dn > b1) {
            s = dn - b1;
          }
          if (a2 > (b1 + b2)) {
            s = fem::min(s, a2 - (b1 + b2));
          }
          s = fem::max(s, third * dmin);
          ttype = -3;
        }
      }
      else {
        //C
        //C              Case 4.
        //C
        ttype = -4;
        s = qurtr * dmin;
        if (dmin == dn) {
          gam = dn;
          a2 = zero;
          if (z(nn - 5) > z(nn - 7)) {
            return;
          }
          b2 = z(nn - 5) / z(nn - 7);
          np = nn - 9;
        }
        else {
          np = nn - 2 * pp;
          b2 = z(np - 2);
          gam = dn1;
          if (z(np - 4) > z(np - 2)) {
            return;
          }
          a2 = z(np - 4) / z(np - 2);
          if (z(nn - 9) > z(nn - 11)) {
            return;
          }
          b2 = z(nn - 9) / z(nn - 11);
          np = nn - 13;
        }
        //C
        //C              Approximate contribution to norm squared from I < NN-1.
        //C
        a2 += b2;
        FEM_DOSTEP(i4, np, 4 * i0 - 1 + pp, -4) {
          if (b2 == zero) {
            goto statement_20;
          }
          b1 = b2;
          if (z(i4) > z(i4 - 2)) {
            return;
          }
          b2 = b2 * (z(i4) / z(i4 - 2));
          a2 += b2;
          if (hundrd * fem::max(b2, b1) < a2 || cnst1 < a2) {
            goto statement_20;
          }
        }
        statement_20:
        a2 = cnst3 * a2;
        //C
        //C              Rayleigh quotient residual bound.
        //C
        if (a2 < cnst1) {
          s = gam * (one - fem::sqrt(a2)) / (one + a2);
        }
      }
    }
    else if (dmin == dn2) {
      //C
      //C           Case 5.
      //C
      ttype = -5;
      s = qurtr * dmin;
      //C
      //C           Compute contribution to norm squared from I > NN-2.
      //C
      np = nn - 2 * pp;
      b1 = z(np - 2);
      b2 = z(np - 6);
      gam = dn2;
      if (z(np - 8) > b2 || z(np - 4) > b1) {
        return;
      }
      a2 = (z(np - 8) / b2) * (one + z(np - 4) / b1);
      //C
      //C           Approximate contribution to norm squared from I < NN-2.
      //C
      if (n0 - i0 > 2) {
        b2 = z(nn - 13) / z(nn - 15);
        a2 += b2;
        FEM_DOSTEP(i4, nn - 17, 4 * i0 - 1 + pp, -4) {
          if (b2 == zero) {
            goto statement_40;
          }
          b1 = b2;
          if (z(i4) > z(i4 - 2)) {
            return;
          }
          b2 = b2 * (z(i4) / z(i4 - 2));
          a2 += b2;
          if (hundrd * fem::max(b2, b1) < a2 || cnst1 < a2) {
            goto statement_40;
          }
        }
        statement_40:
        a2 = cnst3 * a2;
      }
      //C
      if (a2 < cnst1) {
        s = gam * (one - fem::sqrt(a2)) / (one + a2);
      }
    }
    else {
      //C
      //C           Case 6, no information to guide us.
      //C
      if (ttype ==  - 6) {
        g += third * (one - g);
      }
      else if (ttype ==  - 18) {
        g = qurtr * third;
      }
      else {
        g = qurtr;
      }
      s = g * dmin;
      ttype = -6;
    }
    //C
  }
  else if (n0in == (n0 + 1)) {
    //C
    //C        One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN.
    //C
    if (dmin1 == dn1 && dmin2 == dn2) {
      //C
      //C           Cases 7 and 8.
      //C
      ttype = -7;
      s = third * dmin1;
      if (z(nn - 5) > z(nn - 7)) {
        return;
      }
      b1 = z(nn - 5) / z(nn - 7);
      b2 = b1;
      if (b2 == zero) {
        goto statement_60;
      }
      FEM_DOSTEP(i4, 4 * n0 - 9 + pp, 4 * i0 - 1 + pp, -4) {
        a2 = b1;
        if (z(i4) > z(i4 - 2)) {
          return;
        }
        b1 = b1 * (z(i4) / z(i4 - 2));
        b2 += b1;
        if (hundrd * fem::max(b1, a2) < b2) {
          goto statement_60;
        }
      }
      statement_60:
      b2 = fem::sqrt(cnst3 * b2);
      a2 = dmin1 / (one + fem::pow2(b2));
      gap2 = half * dmin2 - a2;
      if (gap2 > zero && gap2 > b2 * a2) {
        s = fem::max(s, a2 * (one - cnst2 * a2 * (b2 / gap2) * b2));
      }
      else {
        s = fem::max(s, a2 * (one - cnst2 * b2));
        ttype = -8;
      }
    }
    else {
      //C
      //C           Case 9.
      //C
      s = qurtr * dmin1;
      if (dmin1 == dn1) {
        s = half * dmin1;
      }
      ttype = -9;
    }
    //C
  }
  else if (n0in == (n0 + 2)) {
    //C
    //C        Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN.
    //C
    //C        Cases 10 and 11.
    //C
    if (dmin2 == dn2 && two * z(nn - 5) < z(nn - 7)) {
      ttype = -10;
      s = third * dmin2;
      if (z(nn - 5) > z(nn - 7)) {
        return;
      }
      b1 = z(nn - 5) / z(nn - 7);
      b2 = b1;
      if (b2 == zero) {
        goto statement_80;
      }
      FEM_DOSTEP(i4, 4 * n0 - 9 + pp, 4 * i0 - 1 + pp, -4) {
        if (z(i4) > z(i4 - 2)) {
          return;
        }
        b1 = b1 * (z(i4) / z(i4 - 2));
        b2 += b1;
        if (hundrd * b1 < b2) {
          goto statement_80;
        }
      }
      statement_80:
      b2 = fem::sqrt(cnst3 * b2);
      a2 = dmin2 / (one + fem::pow2(b2));
      gap2 = z(nn - 7) + z(nn - 9) - fem::sqrt(z(nn - 11)) *
        fem::sqrt(z(nn - 9)) - a2;
      if (gap2 > zero && gap2 > b2 * a2) {
        s = fem::max(s, a2 * (one - cnst2 * a2 * (b2 / gap2) * b2));
      }
      else {
        s = fem::max(s, a2 * (one - cnst2 * b2));
      }
    }
    else {
      s = qurtr * dmin2;
      ttype = -11;
    }
  }
  else if (n0in > (n0 + 2)) {
    //C
    //C        Case 12, more than two eigenvalues deflated. No information.
    //C
    s = zero;
    ttype = -12;
  }
  //C
  tau = s;
  //C
  //C     End of DLASQ4
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasq5.f
inline
void
dlasq5(
  int const& i0,
  int const& n0,
  arr_ref<double> z,
  int const& pp,
  double const& tau,
  double& dmin,
  double& dmin1,
  double& dmin2,
  double& dn,
  double& dnm1,
  double& dnm2,
  bool const& ieee)
{
  z(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2)                                    --
  //C
  //C  -- Contributed by Osni Marques of the Lawrence Berkeley National   --
  //C  -- Laboratory and Beresford Parlett of the Univ. of California at  --
  //C  -- Berkeley                                                        --
  //C  -- November 2008                                                   --
  //C
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASQ5 computes one dqds transform in ping-pong form, one
  //C  version for IEEE machines another for non IEEE machines.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  I0    (input) INTEGER
  //C        First index.
  //C
  //C  N0    (input) INTEGER
  //C        Last index.
  //C
  //C  Z     (input) DOUBLE PRECISION array, dimension ( 4*N )
  //C        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid
  //C        an extra argument.
  //C
  //C  PP    (input) INTEGER
  //C        PP=0 for ping, PP=1 for pong.
  //C
  //C  TAU   (input) DOUBLE PRECISION
  //C        This is the shift.
  //C
  //C  DMIN  (output) DOUBLE PRECISION
  //C        Minimum value of d.
  //C
  //C  DMIN1 (output) DOUBLE PRECISION
  //C        Minimum value of d, excluding D( N0 ).
  //C
  //C  DMIN2 (output) DOUBLE PRECISION
  //C        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
  //C
  //C  DN    (output) DOUBLE PRECISION
  //C        d(N0), the last value of d.
  //C
  //C  DNM1  (output) DOUBLE PRECISION
  //C        d(N0-1).
  //C
  //C  DNM2  (output) DOUBLE PRECISION
  //C        d(N0-2).
  //C
  //C  IEEE  (input) LOGICAL
  //C        Flag for IEEE or non IEEE arithmetic.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameter ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  if ((n0 - i0 - 1) <= 0) {
    return;
  }
  //C
  int j4 = 4 * i0 + pp - 3;
  double emin = z(j4 + 4);
  double d = z(j4) - tau;
  dmin = d;
  dmin1 = -z(j4);
  //C
  double temp = fem::double0;
  int j4p2 = fem::int0;
  const double zero = 0.0e0;
  if (ieee) {
    //C
    //C        Code for IEEE arithmetic.
    //C
    if (pp == 0) {
      FEM_DOSTEP(j4, 4 * i0, 4 * (n0 - 3), 4) {
        z(j4 - 2) = d + z(j4 - 1);
        temp = z(j4 + 1) / z(j4 - 2);
        d = d * temp - tau;
        dmin = fem::min(dmin, d);
        z(j4) = z(j4 - 1) * temp;
        emin = fem::min(z(j4), emin);
      }
    }
    else {
      FEM_DOSTEP(j4, 4 * i0, 4 * (n0 - 3), 4) {
        z(j4 - 3) = d + z(j4);
        temp = z(j4 + 2) / z(j4 - 3);
        d = d * temp - tau;
        dmin = fem::min(dmin, d);
        z(j4 - 1) = z(j4) * temp;
        emin = fem::min(z(j4 - 1), emin);
      }
    }
    //C
    //C        Unroll last two steps.
    //C
    dnm2 = d;
    dmin2 = dmin;
    j4 = 4 * (n0 - 2) - pp;
    j4p2 = j4 + 2 * pp - 1;
    z(j4 - 2) = dnm2 + z(j4p2);
    z(j4) = z(j4p2 + 2) * (z(j4p2) / z(j4 - 2));
    dnm1 = z(j4p2 + 2) * (dnm2 / z(j4 - 2)) - tau;
    dmin = fem::min(dmin, dnm1);
    //C
    dmin1 = dmin;
    j4 += 4;
    j4p2 = j4 + 2 * pp - 1;
    z(j4 - 2) = dnm1 + z(j4p2);
    z(j4) = z(j4p2 + 2) * (z(j4p2) / z(j4 - 2));
    dn = z(j4p2 + 2) * (dnm1 / z(j4 - 2)) - tau;
    dmin = fem::min(dmin, dn);
    //C
  }
  else {
    //C
    //C        Code for non IEEE arithmetic.
    //C
    if (pp == 0) {
      FEM_DOSTEP(j4, 4 * i0, 4 * (n0 - 3), 4) {
        z(j4 - 2) = d + z(j4 - 1);
        if (d < zero) {
          return;
        }
        else {
          z(j4) = z(j4 + 1) * (z(j4 - 1) / z(j4 - 2));
          d = z(j4 + 1) * (d / z(j4 - 2)) - tau;
        }
        dmin = fem::min(dmin, d);
        emin = fem::min(emin, z(j4));
      }
    }
    else {
      FEM_DOSTEP(j4, 4 * i0, 4 * (n0 - 3), 4) {
        z(j4 - 3) = d + z(j4);
        if (d < zero) {
          return;
        }
        else {
          z(j4 - 1) = z(j4 + 2) * (z(j4) / z(j4 - 3));
          d = z(j4 + 2) * (d / z(j4 - 3)) - tau;
        }
        dmin = fem::min(dmin, d);
        emin = fem::min(emin, z(j4 - 1));
      }
    }
    //C
    //C        Unroll last two steps.
    //C
    dnm2 = d;
    dmin2 = dmin;
    j4 = 4 * (n0 - 2) - pp;
    j4p2 = j4 + 2 * pp - 1;
    z(j4 - 2) = dnm2 + z(j4p2);
    if (dnm2 < zero) {
      return;
    }
    else {
      z(j4) = z(j4p2 + 2) * (z(j4p2) / z(j4 - 2));
      dnm1 = z(j4p2 + 2) * (dnm2 / z(j4 - 2)) - tau;
    }
    dmin = fem::min(dmin, dnm1);
    //C
    dmin1 = dmin;
    j4 += 4;
    j4p2 = j4 + 2 * pp - 1;
    z(j4 - 2) = dnm1 + z(j4p2);
    if (dnm1 < zero) {
      return;
    }
    else {
      z(j4) = z(j4p2 + 2) * (z(j4p2) / z(j4 - 2));
      dn = z(j4p2 + 2) * (dnm1 / z(j4 - 2)) - tau;
    }
    dmin = fem::min(dmin, dn);
    //C
  }
  //C
  z(j4 + 2) = dn;
  z(4 * n0 - pp) = emin;
  //C
  //C     End of DLASQ5
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasq6.f
inline
void
dlasq6(
  common& cmn,
  int const& i0,
  int const& n0,
  arr_ref<double> z,
  int const& pp,
  double& dmin,
  double& dmin1,
  double& dmin2,
  double& dn,
  double& dnm1,
  double& dnm2)
{
  z(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2)                                    --
  //C
  //C  -- Contributed by Osni Marques of the Lawrence Berkeley National   --
  //C  -- Laboratory and Beresford Parlett of the Univ. of California at  --
  //C  -- Berkeley                                                        --
  //C  -- November 2008                                                   --
  //C
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASQ6 computes one dqd (shift equal to zero) transform in
  //C  ping-pong form, with protection against underflow and overflow.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  I0    (input) INTEGER
  //C        First index.
  //C
  //C  N0    (input) INTEGER
  //C        Last index.
  //C
  //C  Z     (input) DOUBLE PRECISION array, dimension ( 4*N )
  //C        Z holds the qd array. EMIN is stored in Z(4*N0) to avoid
  //C        an extra argument.
  //C
  //C  PP    (input) INTEGER
  //C        PP=0 for ping, PP=1 for pong.
  //C
  //C  DMIN  (output) DOUBLE PRECISION
  //C        Minimum value of d.
  //C
  //C  DMIN1 (output) DOUBLE PRECISION
  //C        Minimum value of d, excluding D( N0 ).
  //C
  //C  DMIN2 (output) DOUBLE PRECISION
  //C        Minimum value of d, excluding D( N0 ) and D( N0-1 ).
  //C
  //C  DN    (output) DOUBLE PRECISION
  //C        d(N0), the last value of d.
  //C
  //C  DNM1  (output) DOUBLE PRECISION
  //C        d(N0-1).
  //C
  //C  DNM2  (output) DOUBLE PRECISION
  //C        d(N0-2).
  //C
  //C  =====================================================================
  //C
  //C     .. Parameter ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Function ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  if ((n0 - i0 - 1) <= 0) {
    return;
  }
  //C
  double safmin = dlamch(cmn, "Safe minimum");
  int j4 = 4 * i0 + pp - 3;
  double emin = z(j4 + 4);
  double d = z(j4);
  dmin = d;
  //C
  const double zero = 0.0e0;
  double temp = fem::double0;
  if (pp == 0) {
    FEM_DOSTEP(j4, 4 * i0, 4 * (n0 - 3), 4) {
      z(j4 - 2) = d + z(j4 - 1);
      if (z(j4 - 2) == zero) {
        z(j4) = zero;
        d = z(j4 + 1);
        dmin = d;
        emin = zero;
      }
      else if (safmin * z(j4 + 1) < z(j4 - 2) && safmin * z(j4 -
        2) < z(j4 + 1)) {
        temp = z(j4 + 1) / z(j4 - 2);
        z(j4) = z(j4 - 1) * temp;
        d = d * temp;
      }
      else {
        z(j4) = z(j4 + 1) * (z(j4 - 1) / z(j4 - 2));
        d = z(j4 + 1) * (d / z(j4 - 2));
      }
      dmin = fem::min(dmin, d);
      emin = fem::min(emin, z(j4));
    }
  }
  else {
    FEM_DOSTEP(j4, 4 * i0, 4 * (n0 - 3), 4) {
      z(j4 - 3) = d + z(j4);
      if (z(j4 - 3) == zero) {
        z(j4 - 1) = zero;
        d = z(j4 + 2);
        dmin = d;
        emin = zero;
      }
      else if (safmin * z(j4 + 2) < z(j4 - 3) && safmin * z(j4 -
        3) < z(j4 + 2)) {
        temp = z(j4 + 2) / z(j4 - 3);
        z(j4 - 1) = z(j4) * temp;
        d = d * temp;
      }
      else {
        z(j4 - 1) = z(j4 + 2) * (z(j4) / z(j4 - 3));
        d = z(j4 + 2) * (d / z(j4 - 3));
      }
      dmin = fem::min(dmin, d);
      emin = fem::min(emin, z(j4 - 1));
    }
  }
  //C
  //C     Unroll last two steps.
  //C
  dnm2 = d;
  dmin2 = dmin;
  j4 = 4 * (n0 - 2) - pp;
  int j4p2 = j4 + 2 * pp - 1;
  z(j4 - 2) = dnm2 + z(j4p2);
  if (z(j4 - 2) == zero) {
    z(j4) = zero;
    dnm1 = z(j4p2 + 2);
    dmin = dnm1;
    emin = zero;
  }
  else if (safmin * z(j4p2 + 2) < z(j4 - 2) && safmin * z(j4 - 2) < z(
    j4p2 + 2)) {
    temp = z(j4p2 + 2) / z(j4 - 2);
    z(j4) = z(j4p2) * temp;
    dnm1 = dnm2 * temp;
  }
  else {
    z(j4) = z(j4p2 + 2) * (z(j4p2) / z(j4 - 2));
    dnm1 = z(j4p2 + 2) * (dnm2 / z(j4 - 2));
  }
  dmin = fem::min(dmin, dnm1);
  //C
  dmin1 = dmin;
  j4 += 4;
  j4p2 = j4 + 2 * pp - 1;
  z(j4 - 2) = dnm1 + z(j4p2);
  if (z(j4 - 2) == zero) {
    z(j4) = zero;
    dn = z(j4p2 + 2);
    dmin = dn;
    emin = zero;
  }
  else if (safmin * z(j4p2 + 2) < z(j4 - 2) && safmin * z(j4 - 2) < z(
    j4p2 + 2)) {
    temp = z(j4p2 + 2) / z(j4 - 2);
    z(j4) = z(j4p2) * temp;
    dn = dnm1 * temp;
  }
  else {
    z(j4) = z(j4p2 + 2) * (z(j4p2) / z(j4 - 2));
    dn = z(j4p2 + 2) * (dnm1 / z(j4 - 2));
  }
  dmin = fem::min(dmin, dn);
  //C
  z(j4 + 2) = dn;
  z(4 * n0 - pp) = emin;
  //C
  //C     End of DLASQ6
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasq3.f
inline
void
dlasq3(
  common& cmn,
  int const& i0,
  int& n0,
  arr_ref<double> z,
  int& pp,
  double& dmin,
  double& sigma,
  double& desig,
  double& qmax,
  int& nfail,
  int& iter,
  int& ndiv,
  bool const& ieee,
  int& ttype,
  double& dmin1,
  double& dmin2,
  double& dn,
  double& dn1,
  double& dn2,
  double& g,
  double& tau)
{
  z(dimension(star));
  int n0in = fem::int0;
  double eps = fem::double0;
  const double hundrd = 100.0e0;
  double tol = fem::double0;
  double tol2 = fem::double0;
  int nn = fem::int0;
  double s = fem::double0;
  const double half = 0.5e0;
  double t = fem::double0;
  const double one = 1.0e0;
  const double zero = 0.0e0;
  const double cbias = 1.50e0;
  int ipn4 = fem::int0;
  int j4 = fem::int0;
  double temp = fem::double0;
  const double two = 2.0e0;
  const double qurtr = 0.250e0;
  //C
  //C  -- LAPACK routine (version 3.2)                                    --
  //C
  //C  -- Contributed by Osni Marques of the Lawrence Berkeley National   --
  //C  -- Laboratory and Beresford Parlett of the Univ. of California at  --
  //C  -- Berkeley                                                        --
  //C  -- November 2008                                                   --
  //C
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASQ3 checks for deflation, computes a shift (TAU) and calls dqds.
  //C  In case of failure it changes shifts, and tries again until output
  //C  is positive.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  I0     (input) INTEGER
  //C         First index.
  //C
  //C  N0     (input) INTEGER
  //C         Last index.
  //C
  //C  Z      (input) DOUBLE PRECISION array, dimension ( 4*N )
  //C         Z holds the qd array.
  //C
  //C  PP     (input/output) INTEGER
  //C         PP=0 for ping, PP=1 for pong.
  //C         PP=2 indicates that flipping was applied to the Z array
  //C         and that the initial tests for deflation should not be
  //C         performed.
  //C
  //C  DMIN   (output) DOUBLE PRECISION
  //C         Minimum value of d.
  //C
  //C  SIGMA  (output) DOUBLE PRECISION
  //C         Sum of shifts used in current segment.
  //C
  //C  DESIG  (input/output) DOUBLE PRECISION
  //C         Lower order part of SIGMA
  //C
  //C  QMAX   (input) DOUBLE PRECISION
  //C         Maximum value of q.
  //C
  //C  NFAIL  (output) INTEGER
  //C         Number of times shift was too big.
  //C
  //C  ITER   (output) INTEGER
  //C         Number of iterations.
  //C
  //C  NDIV   (output) INTEGER
  //C         Number of divisions.
  //C
  //C  IEEE   (input) LOGICAL
  //C         Flag for IEEE or non IEEE arithmetic (passed to DLASQ5).
  //C
  //C  TTYPE  (input/output) INTEGER
  //C         Shift type.
  //C
  //C  DMIN1, DMIN2, DN, DN1, DN2, G, TAU (input/output) DOUBLE PRECISION
  //C         These are passed as arguments in order to save their values
  //C         between calls to DLASQ3.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. External Function ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  n0in = n0;
  eps = dlamch(cmn, "Precision");
  tol = eps * hundrd;
  tol2 = fem::pow2(tol);
  //C
  //C     Check for deflation.
  //C
  statement_10:
  //C
  if (n0 < i0) {
    return;
  }
  if (n0 == i0) {
    goto statement_20;
  }
  nn = 4 * n0 + pp;
  if (n0 == (i0 + 1)) {
    goto statement_40;
  }
  //C
  //C     Check whether E(N0-1) is negligible, 1 eigenvalue.
  //C
  if (z(nn - 5) > tol2 * (sigma + z(nn - 3)) && z(nn - 2 * pp -
      4) > tol2 * z(nn - 7)) {
    goto statement_30;
  }
  //C
  statement_20:
  //C
  z(4 * n0 - 3) = z(4 * n0 + pp - 3) + sigma;
  n0 = n0 - 1;
  goto statement_10;
  //C
  //C     Check  whether E(N0-2) is negligible, 2 eigenvalues.
  //C
  statement_30:
  //C
  if (z(nn - 9) > tol2 * sigma && z(nn - 2 * pp - 8) > tol2 * z(nn - 11)) {
    goto statement_50;
  }
  //C
  statement_40:
  //C
  if (z(nn - 3) > z(nn - 7)) {
    s = z(nn - 3);
    z(nn - 3) = z(nn - 7);
    z(nn - 7) = s;
  }
  if (z(nn - 5) > z(nn - 3) * tol2) {
    t = half * ((z(nn - 7) - z(nn - 3)) + z(nn - 5));
    s = z(nn - 3) * (z(nn - 5) / t);
    if (s <= t) {
      s = z(nn - 3) * (z(nn - 5) / (t * (one + fem::sqrt(one + s / t))));
    }
    else {
      s = z(nn - 3) * (z(nn - 5) / (t + fem::sqrt(t) * fem::sqrt(t + s)));
    }
    t = z(nn - 7) + (s + z(nn - 5));
    z(nn - 3) = z(nn - 3) * (z(nn - 7) / t);
    z(nn - 7) = t;
  }
  z(4 * n0 - 7) = z(nn - 7) + sigma;
  z(4 * n0 - 3) = z(nn - 3) + sigma;
  n0 = n0 - 2;
  goto statement_10;
  //C
  statement_50:
  if (pp == 2) {
    pp = 0;
  }
  //C
  //C     Reverse the qd-array, if warranted.
  //C
  if (dmin <= zero || n0 < n0in) {
    if (cbias * z(4 * i0 + pp - 3) < z(4 * n0 + pp - 3)) {
      ipn4 = 4 * (i0 + n0);
      FEM_DOSTEP(j4, 4 * i0, 2 * (i0 + n0 - 1), 4) {
        temp = z(j4 - 3);
        z(j4 - 3) = z(ipn4 - j4 - 3);
        z(ipn4 - j4 - 3) = temp;
        temp = z(j4 - 2);
        z(j4 - 2) = z(ipn4 - j4 - 2);
        z(ipn4 - j4 - 2) = temp;
        temp = z(j4 - 1);
        z(j4 - 1) = z(ipn4 - j4 - 5);
        z(ipn4 - j4 - 5) = temp;
        temp = z(j4);
        z(j4) = z(ipn4 - j4 - 4);
        z(ipn4 - j4 - 4) = temp;
      }
      if (n0 - i0 <= 4) {
        z(4 * n0 + pp - 1) = z(4 * i0 + pp - 1);
        z(4 * n0 - pp) = z(4 * i0 - pp);
      }
      dmin2 = fem::min(dmin2, z(4 * n0 + pp - 1));
      z(4 * n0 + pp - 1) = fem::min(z(4 * n0 + pp - 1), z(4 * i0 + pp - 1),
        z(4 * i0 + pp + 3));
      z(4 * n0 - pp) = fem::min(z(4 * n0 - pp), z(4 * i0 - pp), z(4 *
        i0 - pp + 4));
      qmax = fem::max(qmax, z(4 * i0 + pp - 3), z(4 * i0 + pp + 1));
      dmin = -zero;
    }
  }
  //C
  //C     Choose a shift.
  //C
  dlasq4(i0, n0, z, pp, n0in, dmin, dmin1, dmin2, dn, dn1, dn2, tau, ttype, g);
  //C
  //C     Call dqds until DMIN > 0.
  //C
  statement_70:
  //C
  dlasq5(i0, n0, z, pp, tau, dmin, dmin1, dmin2, dn, dn1, dn2, ieee);
  //C
  ndiv += (n0 - i0 + 2);
  iter++;
  //C
  //C     Check status.
  //C
  if (dmin >= zero && dmin1 > zero) {
    //C
    //C        Success.
    //C
    goto statement_90;
    //C
  }
  else if (dmin < zero && dmin1 > zero && z(4 * (n0 - 1) -
    pp) < tol * (sigma + dn1) && fem::abs(dn) < tol * sigma) {
    //C
    //C        Convergence hidden by negative DN.
    //C
    z(4 * (n0 - 1) - pp + 2) = zero;
    dmin = zero;
    goto statement_90;
  }
  else if (dmin < zero) {
    //C
    //C        TAU too big. Select new TAU and try again.
    //C
    nfail++;
    if (ttype <  - 22) {
      //C
      //C           Failed twice. Play it safe.
      //C
      tau = zero;
    }
    else if (dmin1 > zero) {
      //C
      //C           Late failure. Gives excellent shift.
      //C
      tau = (tau + dmin) * (one - two * eps);
      ttype = ttype - 11;
    }
    else {
      //C
      //C           Early failure. Divide by 4.
      //C
      tau = qurtr * tau;
      ttype = ttype - 12;
    }
    goto statement_70;
  }
  else if (disnan(dmin)) {
    //C
    //C        NaN.
    //C
    if (tau == zero) {
      goto statement_80;
    }
    else {
      tau = zero;
      goto statement_70;
    }
  }
  else {
    //C
    //C        Possible underflow. Play it safe.
    //C
    goto statement_80;
  }
  //C
  //C     Risk of underflow.
  //C
  statement_80:
  dlasq6(cmn, i0, n0, z, pp, dmin, dmin1, dmin2, dn, dn1, dn2);
  ndiv += (n0 - i0 + 2);
  iter++;
  tau = zero;
  //C
  statement_90:
  if (tau < sigma) {
    desig += tau;
    t = sigma + desig;
    desig = desig - (t - sigma);
  }
  else {
    t = sigma + tau;
    desig += sigma - (t - tau);
  }
  sigma = t;
  //C
  //C     End of DLASQ3
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasrt.f
inline
void
dlasrt(
  str_cref id,
  int const& n,
  arr_ref<double> d,
  int& info)
{
  d(dimension(star));
  int dir = fem::int0;
  int stkpnt = fem::int0;
  arr_2d<2, 32, int> stack(fem::fill0);
  int start = fem::int0;
  int endd = fem::int0;
  const int select = 20;
  int i = fem::int0;
  int j = fem::int0;
  double dmnmx = fem::double0;
  double d1 = fem::double0;
  double d2 = fem::double0;
  double d3 = fem::double0;
  double tmp = fem::double0;
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  Sort the numbers in D in increasing order (if ID = 'I') or
  //C  in decreasing order (if ID = 'D' ).
  //C
  //C  Use Quick Sort, reverting to Insertion sort on arrays of
  //C  size <= 20. Dimension of STACK limits N to about 2**32.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  ID      (input) CHARACTER*1
  //C          = 'I': sort D in increasing order;
  //C          = 'D': sort D in decreasing order.
  //C
  //C  N       (input) INTEGER
  //C          The length of the array D.
  //C
  //C  D       (input/output) DOUBLE PRECISION array, dimension (N)
  //C          On entry, the array to be sorted.
  //C          On exit, D has been sorted into increasing order
  //C          (D(1) <= ... <= D(N) ) or into decreasing order
  //C          (D(1) >= ... >= D(N) ), depending on ID.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Local Arrays ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input paramters.
  //C
  info = 0;
  dir = -1;
  if (lsame(id, "D")) {
    dir = 0;
  }
  else if (lsame(id, "I")) {
    dir = 1;
  }
  if (dir ==  - 1) {
    info = -1;
  }
  else if (n < 0) {
    info = -2;
  }
  if (info != 0) {
    xerbla("DLASRT", -info);
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (n <= 1) {
    return;
  }
  //C
  stkpnt = 1;
  stack(1, 1) = 1;
  stack(2, 1) = n;
  statement_10:
  start = stack(1, stkpnt);
  endd = stack(2, stkpnt);
  stkpnt = stkpnt - 1;
  if (endd - start <= select && endd - start > 0) {
    //C
    //C        Do Insertion sort on D( START:ENDD )
    //C
    if (dir == 0) {
      //C
      //C           Sort into decreasing order
      //C
      FEM_DO(i, start + 1, endd) {
        FEM_DOSTEP(j, i, start + 1, -1) {
          if (d(j) > d(j - 1)) {
            dmnmx = d(j);
            d(j) = d(j - 1);
            d(j - 1) = dmnmx;
          }
          else {
            goto statement_30;
          }
        }
        statement_30:;
      }
      //C
    }
    else {
      //C
      //C           Sort into increasing order
      //C
      FEM_DO(i, start + 1, endd) {
        FEM_DOSTEP(j, i, start + 1, -1) {
          if (d(j) < d(j - 1)) {
            dmnmx = d(j);
            d(j) = d(j - 1);
            d(j - 1) = dmnmx;
          }
          else {
            goto statement_50;
          }
        }
        statement_50:;
      }
      //C
    }
    //C
  }
  else if (endd - start > select) {
    //C
    //C        Partition D( START:ENDD ) and stack parts, largest one first
    //C
    //C        Choose partition entry as median of 3
    //C
    d1 = d(start);
    d2 = d(endd);
    i = (start + endd) / 2;
    d3 = d(i);
    if (d1 < d2) {
      if (d3 < d1) {
        dmnmx = d1;
      }
      else if (d3 < d2) {
        dmnmx = d3;
      }
      else {
        dmnmx = d2;
      }
    }
    else {
      if (d3 < d2) {
        dmnmx = d2;
      }
      else if (d3 < d1) {
        dmnmx = d3;
      }
      else {
        dmnmx = d1;
      }
    }
    //C
    if (dir == 0) {
      //C
      //C           Sort into decreasing order
      //C
      i = start - 1;
      j = endd + 1;
      statement_60:
      statement_70:
      j = j - 1;
      if (d(j) < dmnmx) {
        goto statement_70;
      }
      statement_80:
      i++;
      if (d(i) > dmnmx) {
        goto statement_80;
      }
      if (i < j) {
        tmp = d(i);
        d(i) = d(j);
        d(j) = tmp;
        goto statement_60;
      }
      if (j - start > endd - j - 1) {
        stkpnt++;
        stack(1, stkpnt) = start;
        stack(2, stkpnt) = j;
        stkpnt++;
        stack(1, stkpnt) = j + 1;
        stack(2, stkpnt) = endd;
      }
      else {
        stkpnt++;
        stack(1, stkpnt) = j + 1;
        stack(2, stkpnt) = endd;
        stkpnt++;
        stack(1, stkpnt) = start;
        stack(2, stkpnt) = j;
      }
    }
    else {
      //C
      //C           Sort into increasing order
      //C
      i = start - 1;
      j = endd + 1;
      statement_90:
      statement_100:
      j = j - 1;
      if (d(j) > dmnmx) {
        goto statement_100;
      }
      statement_110:
      i++;
      if (d(i) < dmnmx) {
        goto statement_110;
      }
      if (i < j) {
        tmp = d(i);
        d(i) = d(j);
        d(j) = tmp;
        goto statement_90;
      }
      if (j - start > endd - j - 1) {
        stkpnt++;
        stack(1, stkpnt) = start;
        stack(2, stkpnt) = j;
        stkpnt++;
        stack(1, stkpnt) = j + 1;
        stack(2, stkpnt) = endd;
      }
      else {
        stkpnt++;
        stack(1, stkpnt) = j + 1;
        stack(2, stkpnt) = endd;
        stkpnt++;
        stack(1, stkpnt) = start;
        stack(2, stkpnt) = j;
      }
    }
  }
  if (stkpnt > 0) {
    goto statement_10;
  }
  //C
  //C     End of DLASRT
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasq2.f
inline
void
dlasq2(
  common& cmn,
  int const& n,
  arr_ref<double> z,
  int& info)
{
  z(dimension(star));
  double eps = fem::double0;
  double safmin = fem::double0;
  const double hundrd = 100.0e0;
  double tol = fem::double0;
  double tol2 = fem::double0;
  const double zero = 0.0e0;
  double d = fem::double0;
  const double half = 0.5e0;
  double t = fem::double0;
  double s = fem::double0;
  const double one = 1.0e0;
  double emin = fem::double0;
  double qmax = fem::double0;
  double zmax = fem::double0;
  double e = fem::double0;
  int k = fem::int0;
  int iinfo = fem::int0;
  double trace = fem::double0;
  bool ieee = fem::bool0;
  int i0 = fem::int0;
  int n0 = fem::int0;
  const double cbias = 1.50e0;
  int ipn4 = fem::int0;
  int i4 = fem::int0;
  double temp = fem::double0;
  int pp = fem::int0;
  int ttype = fem::int0;
  double dmin1 = fem::double0;
  double dmin2 = fem::double0;
  double dn = fem::double0;
  double dn1 = fem::double0;
  double dn2 = fem::double0;
  double g = fem::double0;
  double tau = fem::double0;
  int iter = fem::int0;
  int nfail = fem::int0;
  int ndiv = fem::int0;
  int iwhila = fem::int0;
  double desig = fem::double0;
  double sigma = fem::double0;
  double emax = fem::double0;
  double qmin = fem::double0;
  const double four = 4.0e0;
  double dee = fem::double0;
  double deemin = fem::double0;
  int kmin = fem::int0;
  const double two = 2.0e0;
  double dmin = fem::double0;
  int nbig = fem::int0;
  int iwhilb = fem::int0;
  int splt = fem::int0;
  double oldemn = fem::double0;
  //C
  //C  -- LAPACK routine (version 3.2)                                    --
  //C
  //C  -- Contributed by Osni Marques of the Lawrence Berkeley National   --
  //C  -- Laboratory and Beresford Parlett of the Univ. of California at  --
  //C  -- Berkeley                                                        --
  //C  -- November 2008                                                   --
  //C
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASQ2 computes all the eigenvalues of the symmetric positive
  //C  definite tridiagonal matrix associated with the qd array Z to high
  //C  relative accuracy are computed to high relative accuracy, in the
  //C  absence of denormalization, underflow and overflow.
  //C
  //C  To see the relation of Z to the tridiagonal matrix, let L be a
  //C  unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and
  //C  let U be an upper bidiagonal matrix with 1's above and diagonal
  //C  Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the
  //C  symmetric tridiagonal to which it is similar.
  //C
  //C  Note : DLASQ2 defines a logical variable, IEEE, which is true
  //C  on machines which follow ieee-754 floating-point standard in their
  //C  handling of infinities and NaNs, and false otherwise. This variable
  //C  is passed to DLASQ3.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  N     (input) INTEGER
  //C        The number of rows and columns in the matrix. N >= 0.
  //C
  //C  Z     (input/output) DOUBLE PRECISION array, dimension ( 4*N )
  //C        On entry Z holds the qd array. On exit, entries 1 to N hold
  //C        the eigenvalues in decreasing order, Z( 2*N+1 ) holds the
  //C        trace, and Z( 2*N+2 ) holds the sum of the eigenvalues. If
  //C        N > 2, then Z( 2*N+3 ) holds the iteration count, Z( 2*N+4 )
  //C        holds NDIVS/NIN^2, and Z( 2*N+5 ) holds the percentage of
  //C        shifts that failed.
  //C
  //C  INFO  (output) INTEGER
  //C        = 0: successful exit
  //C        < 0: if the i-th argument is a scalar and had an illegal
  //C             value, then INFO = -i, if the i-th argument is an
  //C             array and the j-entry had an illegal value, then
  //C             INFO = -(i*100+j)
  //C        > 0: the algorithm failed
  //C              = 1, a split was marked by a positive value in E
  //C              = 2, current block of Z not diagonalized after 30*N
  //C                   iterations (in inner while loop)
  //C              = 3, termination criterion of outer while loop not met
  //C                   (program created more than N unreduced blocks)
  //C
  //C  Further Details
  //C  ===============
  //C  Local Variables: I0:N0 defines a current unreduced segment of Z.
  //C  The shifts are accumulated in SIGMA. Iteration count is in ITER.
  //C  Ping-pong is controlled by PP (alternates between 0 and 1).
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments.
  //C     (in case DLASQ2 is not called by DLASQ1)
  //C
  info = 0;
  eps = dlamch(cmn, "Precision");
  safmin = dlamch(cmn, "Safe minimum");
  tol = eps * hundrd;
  tol2 = fem::pow2(tol);
  //C
  if (n < 0) {
    info = -1;
    xerbla("DLASQ2", 1);
    return;
  }
  else if (n == 0) {
    return;
  }
  else if (n == 1) {
    //C
    //C        1-by-1 case.
    //C
    if (z(1) < zero) {
      info = -201;
      xerbla("DLASQ2", 2);
    }
    return;
  }
  else if (n == 2) {
    //C
    //C        2-by-2 case.
    //C
    if (z(2) < zero || z(3) < zero) {
      info = -2;
      xerbla("DLASQ2", 2);
      return;
    }
    else if (z(3) > z(1)) {
      d = z(3);
      z(3) = z(1);
      z(1) = d;
    }
    z(5) = z(1) + z(2) + z(3);
    if (z(2) > z(3) * tol2) {
      t = half * ((z(1) - z(3)) + z(2));
      s = z(3) * (z(2) / t);
      if (s <= t) {
        s = z(3) * (z(2) / (t * (one + fem::sqrt(one + s / t))));
      }
      else {
        s = z(3) * (z(2) / (t + fem::sqrt(t) * fem::sqrt(t + s)));
      }
      t = z(1) + (s + z(2));
      z(3) = z(3) * (z(1) / t);
      z(1) = t;
    }
    z(2) = z(3);
    z(6) = z(2) + z(1);
    return;
  }
  //C
  //C     Check for negative data and compute sums of q's and e's.
  //C
  z(2 * n) = zero;
  emin = z(2);
  qmax = zero;
  zmax = zero;
  d = zero;
  e = zero;
  //C
  FEM_DOSTEP(k, 1, 2 * (n - 1), 2) {
    if (z(k) < zero) {
      info = -(200 + k);
      xerbla("DLASQ2", 2);
      return;
    }
    else if (z(k + 1) < zero) {
      info = -(200 + k + 1);
      xerbla("DLASQ2", 2);
      return;
    }
    d += z(k);
    e += z(k + 1);
    qmax = fem::max(qmax, z(k));
    emin = fem::min(emin, z(k + 1));
    zmax = fem::max(qmax, zmax, z(k + 1));
  }
  if (z(2 * n - 1) < zero) {
    info = -(200 + 2 * n - 1);
    xerbla("DLASQ2", 2);
    return;
  }
  d += z(2 * n - 1);
  qmax = fem::max(qmax, z(2 * n - 1));
  zmax = fem::max(qmax, zmax);
  //C
  //C     Check for diagonality.
  //C
  if (e == zero) {
    FEM_DO(k, 2, n) {
      z(k) = z(2 * k - 1);
    }
    dlasrt("D", n, z, iinfo);
    z(2 * n - 1) = d;
    return;
  }
  //C
  trace = d + e;
  //C
  //C     Check for zero data.
  //C
  if (trace == zero) {
    z(2 * n - 1) = zero;
    return;
  }
  //C
  //C     Check whether the machine is IEEE conformable.
  //C
  ieee = ilaenv(10, "DLASQ2", "N", 1, 2, 3, 4) == 1 && ilaenv(11,
    "DLASQ2", "N", 1, 2, 3, 4) == 1;
  //C
  //C     Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).
  //C
  FEM_DOSTEP(k, 2 * n, 2, -2) {
    z(2 * k) = zero;
    z(2 * k - 1) = z(k);
    z(2 * k - 2) = zero;
    z(2 * k - 3) = z(k - 1);
  }
  //C
  i0 = 1;
  n0 = n;
  //C
  //C     Reverse the qd-array, if warranted.
  //C
  if (cbias * z(4 * i0 - 3) < z(4 * n0 - 3)) {
    ipn4 = 4 * (i0 + n0);
    FEM_DOSTEP(i4, 4 * i0, 2 * (i0 + n0 - 1), 4) {
      temp = z(i4 - 3);
      z(i4 - 3) = z(ipn4 - i4 - 3);
      z(ipn4 - i4 - 3) = temp;
      temp = z(i4 - 1);
      z(i4 - 1) = z(ipn4 - i4 - 5);
      z(ipn4 - i4 - 5) = temp;
    }
  }
  //C
  //C     Initial split checking via dqd and Li's test.
  //C
  pp = 0;
  //C
  FEM_DO(k, 1, 2) {
    //C
    d = z(4 * n0 + pp - 3);
    FEM_DOSTEP(i4, 4 * (n0 - 1) + pp, 4 * i0 + pp, -4) {
      if (z(i4 - 1) <= tol2 * d) {
        z(i4 - 1) = -zero;
        d = z(i4 - 3);
      }
      else {
        d = z(i4 - 3) * (d / (d + z(i4 - 1)));
      }
    }
    //C
    //C        dqd maps Z to ZZ plus Li's test.
    //C
    emin = z(4 * i0 + pp + 1);
    d = z(4 * i0 + pp - 3);
    FEM_DOSTEP(i4, 4 * i0 + pp, 4 * (n0 - 1) + pp, 4) {
      z(i4 - 2 * pp - 2) = d + z(i4 - 1);
      if (z(i4 - 1) <= tol2 * d) {
        z(i4 - 1) = -zero;
        z(i4 - 2 * pp - 2) = d;
        z(i4 - 2 * pp) = zero;
        d = z(i4 + 1);
      }
      else if (safmin * z(i4 + 1) < z(i4 - 2 * pp - 2) && safmin * z(
        i4 - 2 * pp - 2) < z(i4 + 1)) {
        temp = z(i4 + 1) / z(i4 - 2 * pp - 2);
        z(i4 - 2 * pp) = z(i4 - 1) * temp;
        d = d * temp;
      }
      else {
        z(i4 - 2 * pp) = z(i4 + 1) * (z(i4 - 1) / z(i4 - 2 * pp - 2));
        d = z(i4 + 1) * (d / z(i4 - 2 * pp - 2));
      }
      emin = fem::min(emin, z(i4 - 2 * pp));
    }
    z(4 * n0 - pp - 2) = d;
    //C
    //C        Now find qmax.
    //C
    qmax = z(4 * i0 - pp - 2);
    FEM_DOSTEP(i4, 4 * i0 - pp + 2, 4 * n0 - pp - 2, 4) {
      qmax = fem::max(qmax, z(i4));
    }
    //C
    //C        Prepare for the next iteration on K.
    //C
    pp = 1 - pp;
  }
  //C
  //C     Initialise variables to pass to DLASQ3.
  //C
  ttype = 0;
  dmin1 = zero;
  dmin2 = zero;
  dn = zero;
  dn1 = zero;
  dn2 = zero;
  g = zero;
  tau = zero;
  //C
  iter = 2;
  nfail = 0;
  ndiv = 2 * (n0 - i0);
  //C
  {
    int fem_do_last = n + 1;
    FEM_DO(iwhila, 1, fem_do_last) {
      if (n0 < 1) {
        goto statement_170;
      }
      //C
      //C        While array unfinished do
      //C
      //C        E(N0) holds the value of SIGMA when submatrix in I0:N0
      //C        splits from the rest of the array, but is negated.
      //C
      desig = zero;
      if (n0 == n) {
        sigma = zero;
      }
      else {
        sigma = -z(4 * n0 - 1);
      }
      if (sigma < zero) {
        info = 1;
        return;
      }
      //C
      //C        Find last unreduced submatrix's top index I0, find QMAX and
      //C        EMIN. Find Gershgorin-type bound if Q's much greater than E's.
      //C
      emax = zero;
      if (n0 > i0) {
        emin = fem::abs(z(4 * n0 - 5));
      }
      else {
        emin = zero;
      }
      qmin = z(4 * n0 - 3);
      qmax = qmin;
      FEM_DOSTEP(i4, 4 * n0, 8, -4) {
        if (z(i4 - 5) <= zero) {
          goto statement_100;
        }
        if (qmin >= four * emax) {
          qmin = fem::min(qmin, z(i4 - 3));
          emax = fem::max(emax, z(i4 - 5));
        }
        qmax = fem::max(qmax, z(i4 - 7) + z(i4 - 5));
        emin = fem::min(emin, z(i4 - 5));
      }
      i4 = 4;
      //C
      statement_100:
      i0 = i4 / 4;
      pp = 0;
      //C
      if (n0 - i0 > 1) {
        dee = z(4 * i0 - 3);
        deemin = dee;
        kmin = i0;
        FEM_DOSTEP(i4, 4 * i0 + 1, 4 * n0 - 3, 4) {
          dee = z(i4) * (dee / (dee + z(i4 - 2)));
          if (dee <= deemin) {
            deemin = dee;
            kmin = (i4 + 3) / 4;
          }
        }
        if ((kmin - i0) * 2 < n0 - kmin && deemin <= half * z(4 * n0 - 3)) {
          ipn4 = 4 * (i0 + n0);
          pp = 2;
          FEM_DOSTEP(i4, 4 * i0, 2 * (i0 + n0 - 1), 4) {
            temp = z(i4 - 3);
            z(i4 - 3) = z(ipn4 - i4 - 3);
            z(ipn4 - i4 - 3) = temp;
            temp = z(i4 - 2);
            z(i4 - 2) = z(ipn4 - i4 - 2);
            z(ipn4 - i4 - 2) = temp;
            temp = z(i4 - 1);
            z(i4 - 1) = z(ipn4 - i4 - 5);
            z(ipn4 - i4 - 5) = temp;
            temp = z(i4);
            z(i4) = z(ipn4 - i4 - 4);
            z(ipn4 - i4 - 4) = temp;
          }
        }
      }
      //C
      //C        Put -(initial shift) into DMIN.
      //C
      dmin = -fem::max(zero, qmin - two * fem::sqrt(qmin) * fem::sqrt(emax));
      //C
      //C        Now I0:N0 is unreduced.
      //C        PP = 0 for ping, PP = 1 for pong.
      //C        PP = 2 indicates that flipping was applied to the Z array and
      //C               and that the tests for deflation upon entry in DLASQ3
      //C               should not be performed.
      //C
      nbig = 30 * (n0 - i0 + 1);
      FEM_DO(iwhilb, 1, nbig) {
        if (i0 > n0) {
          goto statement_150;
        }
        //C
        //C           While submatrix unfinished take a good dqds step.
        //C
        dlasq3(cmn, i0, n0, z, pp, dmin, sigma, desig, qmax, nfail,
          iter, ndiv, ieee, ttype, dmin1, dmin2, dn, dn1, dn2, g,
          tau);
        //C
        pp = 1 - pp;
        //C
        //C           When EMIN is very small check for splits.
        //C
        if (pp == 0 && n0 - i0 >= 3) {
          if (z(4 * n0) <= tol2 * qmax || z(4 * n0 - 1) <= tol2 * sigma) {
            splt = i0 - 1;
            qmax = z(4 * i0 - 3);
            emin = z(4 * i0 - 1);
            oldemn = z(4 * i0);
            FEM_DOSTEP(i4, 4 * i0, 4 * (n0 - 3), 4) {
              if (z(i4) <= tol2 * z(i4 - 3) || z(i4 - 1) <= tol2 * sigma) {
                z(i4 - 1) = -sigma;
                splt = i4 / 4;
                qmax = zero;
                emin = z(i4 + 3);
                oldemn = z(i4 + 4);
              }
              else {
                qmax = fem::max(qmax, z(i4 + 1));
                emin = fem::min(emin, z(i4 - 1));
                oldemn = fem::min(oldemn, z(i4));
              }
            }
            z(4 * n0 - 1) = emin;
            z(4 * n0) = oldemn;
            i0 = splt + 1;
          }
        }
        //C
      }
      //C
      info = 2;
      return;
      //C
      //C        end IWHILB
      //C
      statement_150:;
      //C
    }
  }
  //C
  info = 3;
  return;
  //C
  //C     end IWHILA
  //C
  statement_170:
  //C
  //C     Move q's to the front.
  //C
  FEM_DO(k, 2, n) {
    z(k) = z(4 * k - 3);
  }
  //C
  //C     Sort and compute sum of eigenvalues.
  //C
  dlasrt("D", n, z, iinfo);
  //C
  e = zero;
  FEM_DOSTEP(k, n, 1, -1) {
    e += z(k);
  }
  //C
  //C     Store trace, sum(eigenvalues) and information on performance.
  //C
  z(2 * n + 1) = trace;
  z(2 * n + 2) = e;
  z(2 * n + 3) = fem::dble(iter);
  z(2 * n + 4) = fem::dble(ndiv) / fem::dble(fem::pow2(n));
  z(2 * n + 5) = hundrd * nfail / fem::dble(iter);
  //C
  //C     End of DLASQ2
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasq1.f
inline
void
dlasq1(
  common& cmn,
  int const& n,
  arr_ref<double> d,
  arr_cref<double> e,
  arr_ref<double> work,
  int& info)
{
  d(dimension(star));
  e(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2)                                    --
  //C
  //C  -- Contributed by Osni Marques of the Lawrence Berkeley National   --
  //C  -- Laboratory and Beresford Parlett of the Univ. of California at  --
  //C  -- Berkeley                                                        --
  //C  -- November 2008                                                   --
  //C
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASQ1 computes the singular values of a real N-by-N bidiagonal
  //C  matrix with diagonal D and off-diagonal E. The singular values
  //C  are computed to high relative accuracy, in the absence of
  //C  denormalization, underflow and overflow. The algorithm was first
  //C  presented in
  //C
  //C  "Accurate singular values and differential qd algorithms" by K. V.
  //C  Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,
  //C  1994,
  //C
  //C  and the present implementation is described in "An implementation of
  //C  the dqds Algorithm (Positive Case)", LAPACK Working Note.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  N     (input) INTEGER
  //C        The number of rows and columns in the matrix. N >= 0.
  //C
  //C  D     (input/output) DOUBLE PRECISION array, dimension (N)
  //C        On entry, D contains the diagonal elements of the
  //C        bidiagonal matrix whose SVD is desired. On normal exit,
  //C        D contains the singular values in decreasing order.
  //C
  //C  E     (input/output) DOUBLE PRECISION array, dimension (N)
  //C        On entry, elements E(1:N-1) contain the off-diagonal elements
  //C        of the bidiagonal matrix whose SVD is desired.
  //C        On exit, E is overwritten.
  //C
  //C  WORK  (workspace) DOUBLE PRECISION array, dimension (4*N)
  //C
  //C  INFO  (output) INTEGER
  //C        = 0: successful exit
  //C        < 0: if INFO = -i, the i-th argument had an illegal value
  //C        > 0: the algorithm failed
  //C             = 1, a split was marked by a positive value in E
  //C             = 2, current block of Z not diagonalized after 30*N
  //C                  iterations (in inner while loop)
  //C             = 3, termination criterion of outer while loop not met
  //C                  (program created more than N unreduced blocks)
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  info = 0;
  double sigmn = fem::double0;
  double sigmx = fem::double0;
  if (n < 0) {
    info = -2;
    xerbla("DLASQ1", -info);
    return;
  }
  else if (n == 0) {
    return;
  }
  else if (n == 1) {
    d(1) = fem::abs(d(1));
    return;
  }
  else if (n == 2) {
    dlas2(d(1), e(1), d(2), sigmn, sigmx);
    d(1) = sigmx;
    d(2) = sigmn;
    return;
  }
  //C
  //C     Estimate the largest singular value.
  //C
  const double zero = 0.0e0;
  sigmx = zero;
  int i = fem::int0;
  {
    int fem_do_last = n - 1;
    FEM_DO(i, 1, fem_do_last) {
      d(i) = fem::abs(d(i));
      sigmx = fem::max(sigmx, fem::abs(e(i)));
    }
  }
  d(n) = fem::abs(d(n));
  //C
  //C     Early return if SIGMX is zero (matrix is already diagonal).
  //C
  int iinfo = fem::int0;
  if (sigmx == zero) {
    dlasrt("D", n, d, iinfo);
    return;
  }
  //C
  FEM_DO(i, 1, n) {
    sigmx = fem::max(sigmx, d(i));
  }
  //C
  //C     Copy D and E into WORK (in the Z format) and scale (squaring the
  //C     input data makes scaling by a power of the radix pointless).
  //C
  double eps = dlamch(cmn, "Precision");
  double safmin = dlamch(cmn, "Safe minimum");
  double scale = fem::sqrt(eps / safmin);
  dcopy(n, d, 1, work(1), 2);
  dcopy(n - 1, e, 1, work(2), 2);
  dlascl(cmn, "G", 0, 0, sigmx, scale, 2 * n - 1, 1, work, 2 * n - 1, iinfo);
  //C
  //C     Compute the q's and e's.
  //C
  {
    int fem_do_last = 2 * n - 1;
    FEM_DO(i, 1, fem_do_last) {
      work(i) = fem::pow2(work(i));
    }
  }
  work(2 * n) = zero;
  //C
  dlasq2(cmn, n, work, info);
  //C
  if (info == 0) {
    FEM_DO(i, 1, n) {
      d(i) = fem::sqrt(work(i));
    }
    dlascl(cmn, "G", 0, 0, scale, sigmx, n, 1, d, n, iinfo);
  }
  //C
  //C     End of DLASQ1
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasr.f
inline
void
dlasr(
  str_cref side,
  str_cref pivot,
  str_cref direct,
  int const& m,
  int const& n,
  arr_cref<double> c,
  arr_cref<double> s,
  arr_ref<double, 2> a,
  int const& lda)
{
  c(dimension(star));
  s(dimension(star));
  a(dimension(lda, star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASR applies a sequence of plane rotations to a real matrix A,
  //C  from either the left or the right.
  //C
  //C  When SIDE = 'L', the transformation takes the form
  //C
  //C     A := P*A
  //C
  //C  and when SIDE = 'R', the transformation takes the form
  //C
  //C     A := A*P**T
  //C
  //C  where P is an orthogonal matrix consisting of a sequence of z plane
  //C  rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
  //C  and P**T is the transpose of P.
  //C
  //C  When DIRECT = 'F' (Forward sequence), then
  //C
  //C     P = P(z-1) * ... * P(2) * P(1)
  //C
  //C  and when DIRECT = 'B' (Backward sequence), then
  //C
  //C     P = P(1) * P(2) * ... * P(z-1)
  //C
  //C  where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
  //C
  //C     R(k) = (  c(k)  s(k) )
  //C          = ( -s(k)  c(k) ).
  //C
  //C  When PIVOT = 'V' (Variable pivot), the rotation is performed
  //C  for the plane (k,k+1), i.e., P(k) has the form
  //C
  //C     P(k) = (  1                                            )
  //C            (       ...                                     )
  //C            (              1                                )
  //C            (                   c(k)  s(k)                  )
  //C            (                  -s(k)  c(k)                  )
  //C            (                                1              )
  //C            (                                     ...       )
  //C            (                                            1  )
  //C
  //C  where R(k) appears as a rank-2 modification to the identity matrix in
  //C  rows and columns k and k+1.
  //C
  //C  When PIVOT = 'T' (Top pivot), the rotation is performed for the
  //C  plane (1,k+1), so P(k) has the form
  //C
  //C     P(k) = (  c(k)                    s(k)                 )
  //C            (         1                                     )
  //C            (              ...                              )
  //C            (                     1                         )
  //C            ( -s(k)                    c(k)                 )
  //C            (                                 1             )
  //C            (                                      ...      )
  //C            (                                             1 )
  //C
  //C  where R(k) appears in rows and columns 1 and k+1.
  //C
  //C  Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
  //C  performed for the plane (k,z), giving P(k) the form
  //C
  //C     P(k) = ( 1                                             )
  //C            (      ...                                      )
  //C            (             1                                 )
  //C            (                  c(k)                    s(k) )
  //C            (                         1                     )
  //C            (                              ...              )
  //C            (                                     1         )
  //C            (                 -s(k)                    c(k) )
  //C
  //C  where R(k) appears in rows and columns k and z.  The rotations are
  //C  performed without ever forming P(k) explicitly.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  SIDE    (input) CHARACTER*1
  //C          Specifies whether the plane rotation matrix P is applied to
  //C          A on the left or the right.
  //C          = 'L':  Left, compute A := P*A
  //C          = 'R':  Right, compute A:= A*P**T
  //C
  //C  PIVOT   (input) CHARACTER*1
  //C          Specifies the plane for which P(k) is a plane rotation
  //C          matrix.
  //C          = 'V':  Variable pivot, the plane (k,k+1)
  //C          = 'T':  Top pivot, the plane (1,k+1)
  //C          = 'B':  Bottom pivot, the plane (k,z)
  //C
  //C  DIRECT  (input) CHARACTER*1
  //C          Specifies whether P is a forward or backward sequence of
  //C          plane rotations.
  //C          = 'F':  Forward, P = P(z-1)*...*P(2)*P(1)
  //C          = 'B':  Backward, P = P(1)*P(2)*...*P(z-1)
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix A.  If m <= 1, an immediate
  //C          return is effected.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix A.  If n <= 1, an
  //C          immediate return is effected.
  //C
  //C  C       (input) DOUBLE PRECISION array, dimension
  //C                  (M-1) if SIDE = 'L'
  //C                  (N-1) if SIDE = 'R'
  //C          The cosines c(k) of the plane rotations.
  //C
  //C  S       (input) DOUBLE PRECISION array, dimension
  //C                  (M-1) if SIDE = 'L'
  //C                  (N-1) if SIDE = 'R'
  //C          The sines s(k) of the plane rotations.  The 2-by-2 plane
  //C          rotation part of the matrix P(k), R(k), has the form
  //C          R(k) = (  c(k)  s(k) )
  //C                 ( -s(k)  c(k) ).
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          The M-by-N matrix A.  On exit, A is overwritten by P*A if
  //C          SIDE = 'R' or by A*P**T if SIDE = 'L'.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(1,M).
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters
  //C
  int info = 0;
  if (!(lsame(side, "L") || lsame(side, "R"))) {
    info = 1;
  }
  else if (!(lsame(pivot, "V") || lsame(pivot, "T") || lsame(pivot, "B"))) {
    info = 2;
  }
  else if (!(lsame(direct, "F") || lsame(direct, "B"))) {
    info = 3;
  }
  else if (m < 0) {
    info = 4;
  }
  else if (n < 0) {
    info = 5;
  }
  else if (lda < fem::max(1, m)) {
    info = 9;
  }
  if (info != 0) {
    xerbla("DLASR ", info);
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if ((m == 0) || (n == 0)) {
    return;
  }
  int j = fem::int0;
  double ctemp = fem::double0;
  double stemp = fem::double0;
  const double one = 1.0e+0;
  const double zero = 0.0e+0;
  int i = fem::int0;
  double temp = fem::double0;
  if (lsame(side, "L")) {
    //C
    //C        Form  P * A
    //C
    if (lsame(pivot, "V")) {
      if (lsame(direct, "F")) {
        {
          int fem_do_last = m - 1;
          FEM_DO(j, 1, fem_do_last) {
            ctemp = c(j);
            stemp = s(j);
            if ((ctemp != one) || (stemp != zero)) {
              FEM_DO(i, 1, n) {
                temp = a(j + 1, i);
                a(j + 1, i) = ctemp * temp - stemp * a(j, i);
                a(j, i) = stemp * temp + ctemp * a(j, i);
              }
            }
          }
        }
      }
      else if (lsame(direct, "B")) {
        FEM_DOSTEP(j, m - 1, 1, -1) {
          ctemp = c(j);
          stemp = s(j);
          if ((ctemp != one) || (stemp != zero)) {
            FEM_DO(i, 1, n) {
              temp = a(j + 1, i);
              a(j + 1, i) = ctemp * temp - stemp * a(j, i);
              a(j, i) = stemp * temp + ctemp * a(j, i);
            }
          }
        }
      }
    }
    else if (lsame(pivot, "T")) {
      if (lsame(direct, "F")) {
        FEM_DO(j, 2, m) {
          ctemp = c(j - 1);
          stemp = s(j - 1);
          if ((ctemp != one) || (stemp != zero)) {
            FEM_DO(i, 1, n) {
              temp = a(j, i);
              a(j, i) = ctemp * temp - stemp * a(1, i);
              a(1, i) = stemp * temp + ctemp * a(1, i);
            }
          }
        }
      }
      else if (lsame(direct, "B")) {
        FEM_DOSTEP(j, m, 2, -1) {
          ctemp = c(j - 1);
          stemp = s(j - 1);
          if ((ctemp != one) || (stemp != zero)) {
            FEM_DO(i, 1, n) {
              temp = a(j, i);
              a(j, i) = ctemp * temp - stemp * a(1, i);
              a(1, i) = stemp * temp + ctemp * a(1, i);
            }
          }
        }
      }
    }
    else if (lsame(pivot, "B")) {
      if (lsame(direct, "F")) {
        {
          int fem_do_last = m - 1;
          FEM_DO(j, 1, fem_do_last) {
            ctemp = c(j);
            stemp = s(j);
            if ((ctemp != one) || (stemp != zero)) {
              FEM_DO(i, 1, n) {
                temp = a(j, i);
                a(j, i) = stemp * a(m, i) + ctemp * temp;
                a(m, i) = ctemp * a(m, i) - stemp * temp;
              }
            }
          }
        }
      }
      else if (lsame(direct, "B")) {
        FEM_DOSTEP(j, m - 1, 1, -1) {
          ctemp = c(j);
          stemp = s(j);
          if ((ctemp != one) || (stemp != zero)) {
            FEM_DO(i, 1, n) {
              temp = a(j, i);
              a(j, i) = stemp * a(m, i) + ctemp * temp;
              a(m, i) = ctemp * a(m, i) - stemp * temp;
            }
          }
        }
      }
    }
  }
  else if (lsame(side, "R")) {
    //C
    //C        Form A * P'
    //C
    if (lsame(pivot, "V")) {
      if (lsame(direct, "F")) {
        {
          int fem_do_last = n - 1;
          FEM_DO(j, 1, fem_do_last) {
            ctemp = c(j);
            stemp = s(j);
            if ((ctemp != one) || (stemp != zero)) {
              FEM_DO(i, 1, m) {
                temp = a(i, j + 1);
                a(i, j + 1) = ctemp * temp - stemp * a(i, j);
                a(i, j) = stemp * temp + ctemp * a(i, j);
              }
            }
          }
        }
      }
      else if (lsame(direct, "B")) {
        FEM_DOSTEP(j, n - 1, 1, -1) {
          ctemp = c(j);
          stemp = s(j);
          if ((ctemp != one) || (stemp != zero)) {
            FEM_DO(i, 1, m) {
              temp = a(i, j + 1);
              a(i, j + 1) = ctemp * temp - stemp * a(i, j);
              a(i, j) = stemp * temp + ctemp * a(i, j);
            }
          }
        }
      }
    }
    else if (lsame(pivot, "T")) {
      if (lsame(direct, "F")) {
        FEM_DO(j, 2, n) {
          ctemp = c(j - 1);
          stemp = s(j - 1);
          if ((ctemp != one) || (stemp != zero)) {
            FEM_DO(i, 1, m) {
              temp = a(i, j);
              a(i, j) = ctemp * temp - stemp * a(i, 1);
              a(i, 1) = stemp * temp + ctemp * a(i, 1);
            }
          }
        }
      }
      else if (lsame(direct, "B")) {
        FEM_DOSTEP(j, n, 2, -1) {
          ctemp = c(j - 1);
          stemp = s(j - 1);
          if ((ctemp != one) || (stemp != zero)) {
            FEM_DO(i, 1, m) {
              temp = a(i, j);
              a(i, j) = ctemp * temp - stemp * a(i, 1);
              a(i, 1) = stemp * temp + ctemp * a(i, 1);
            }
          }
        }
      }
    }
    else if (lsame(pivot, "B")) {
      if (lsame(direct, "F")) {
        {
          int fem_do_last = n - 1;
          FEM_DO(j, 1, fem_do_last) {
            ctemp = c(j);
            stemp = s(j);
            if ((ctemp != one) || (stemp != zero)) {
              FEM_DO(i, 1, m) {
                temp = a(i, j);
                a(i, j) = stemp * a(i, n) + ctemp * temp;
                a(i, n) = ctemp * a(i, n) - stemp * temp;
              }
            }
          }
        }
      }
      else if (lsame(direct, "B")) {
        FEM_DOSTEP(j, n - 1, 1, -1) {
          ctemp = c(j);
          stemp = s(j);
          if ((ctemp != one) || (stemp != zero)) {
            FEM_DO(i, 1, m) {
              temp = a(i, j);
              a(i, j) = stemp * a(i, n) + ctemp * temp;
              a(i, n) = ctemp * a(i, n) - stemp * temp;
            }
          }
        }
      }
    }
  }
  //C
  //C     End of DLASR
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasv2.f
inline
void
dlasv2(
  common& cmn,
  double const& f,
  double const& g,
  double const& h,
  double& ssmin,
  double& ssmax,
  double& snr,
  double& csr,
  double& snl,
  double& csl)
{
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASV2 computes the singular value decomposition of a 2-by-2
  //C  triangular matrix
  //C     [  F   G  ]
  //C     [  0   H  ].
  //C  On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the
  //C  smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and
  //C  right singular vectors for abs(SSMAX), giving the decomposition
  //C
  //C     [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]
  //C     [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].
  //C
  //C  Arguments
  //C  =========
  //C
  //C  F       (input) DOUBLE PRECISION
  //C          The (1,1) element of the 2-by-2 matrix.
  //C
  //C  G       (input) DOUBLE PRECISION
  //C          The (1,2) element of the 2-by-2 matrix.
  //C
  //C  H       (input) DOUBLE PRECISION
  //C          The (2,2) element of the 2-by-2 matrix.
  //C
  //C  SSMIN   (output) DOUBLE PRECISION
  //C          abs(SSMIN) is the smaller singular value.
  //C
  //C  SSMAX   (output) DOUBLE PRECISION
  //C          abs(SSMAX) is the larger singular value.
  //C
  //C  SNL     (output) DOUBLE PRECISION
  //C  CSL     (output) DOUBLE PRECISION
  //C          The vector (CSL, SNL) is a unit left singular vector for the
  //C          singular value abs(SSMAX).
  //C
  //C  SNR     (output) DOUBLE PRECISION
  //C  CSR     (output) DOUBLE PRECISION
  //C          The vector (CSR, SNR) is a unit right singular vector for the
  //C          singular value abs(SSMAX).
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Any input parameter may be aliased with any output parameter.
  //C
  //C  Barring over/underflow and assuming a guard digit in subtraction, all
  //C  output quantities are correct to within a few units in the last
  //C  place (ulps).
  //C
  //C  In IEEE arithmetic, the code works correctly if one matrix element is
  //C  infinite.
  //C
  //C  Overflow will not occur unless the largest singular value itself
  //C  overflows or is within a few ulps of overflow. (On machines with
  //C  partial overflow, like the Cray, overflow may occur if the largest
  //C  singular value is within a factor of 2 of overflow.)
  //C
  //C  Underflow is harmless if underflow is gradual. Otherwise, results
  //C  may correspond to a matrix modified by perturbations of size near
  //C  the underflow threshold.
  //C
  //C =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  double ft = f;
  double fa = fem::abs(ft);
  double ht = h;
  double ha = fem::abs(h);
  //C
  //C     PMAX points to the maximum absolute element of matrix
  //C       PMAX = 1 if F largest in absolute values
  //C       PMAX = 2 if G largest in absolute values
  //C       PMAX = 3 if H largest in absolute values
  //C
  int pmax = 1;
  bool swap = (ha > fa);
  double temp = fem::double0;
  if (swap) {
    pmax = 3;
    temp = ft;
    ft = ht;
    ht = temp;
    temp = fa;
    fa = ha;
    ha = temp;
    //C
    //C        Now FA .ge. HA
    //C
  }
  double gt = g;
  double ga = fem::abs(gt);
  const double zero = 0.0e0;
  const double one = 1.0e0;
  double clt = fem::double0;
  double crt = fem::double0;
  double slt = fem::double0;
  double srt = fem::double0;
  bool gasmal = fem::bool0;
  double d = fem::double0;
  double l = fem::double0;
  double m = fem::double0;
  const double two = 2.0e0;
  double t = fem::double0;
  double mm = fem::double0;
  double tt = fem::double0;
  double s = fem::double0;
  double r = fem::double0;
  const double half = 0.5e0;
  double a = fem::double0;
  const double four = 4.0e0;
  if (ga == zero) {
    //C
    //C        Diagonal matrix
    //C
    ssmin = ha;
    ssmax = fa;
    clt = one;
    crt = one;
    slt = zero;
    srt = zero;
  }
  else {
    gasmal = true;
    if (ga > fa) {
      pmax = 2;
      if ((fa / ga) < dlamch(cmn, "EPS")) {
        //C
        //C              Case of very large GA
        //C
        gasmal = false;
        ssmax = ga;
        if (ha > one) {
          ssmin = fa / (ga / ha);
        }
        else {
          ssmin = (fa / ga) * ha;
        }
        clt = one;
        slt = ht / gt;
        srt = one;
        crt = ft / gt;
      }
    }
    if (gasmal) {
      //C
      //C           Normal case
      //C
      d = fa - ha;
      if (d == fa) {
        //C
        //C              Copes with infinite F or H
        //C
        l = one;
      }
      else {
        l = d / fa;
      }
      //C
      //C           Note that 0 .le. L .le. 1
      //C
      m = gt / ft;
      //C
      //C           Note that abs(M) .le. 1/macheps
      //C
      t = two - l;
      //C
      //C           Note that T .ge. 1
      //C
      mm = m * m;
      tt = t * t;
      s = fem::sqrt(tt + mm);
      //C
      //C           Note that 1 .le. S .le. 1 + 1/macheps
      //C
      if (l == zero) {
        r = fem::abs(m);
      }
      else {
        r = fem::sqrt(l * l + mm);
      }
      //C
      //C           Note that 0 .le. R .le. 1 + 1/macheps
      //C
      a = half * (s + r);
      //C
      //C           Note that 1 .le. A .le. 1 + abs(M)
      //C
      ssmin = ha / a;
      ssmax = fa * a;
      if (mm == zero) {
        //C
        //C              Note that M is very tiny
        //C
        if (l == zero) {
          t = fem::sign(two, ft) * fem::sign(one, gt);
        }
        else {
          t = gt / fem::sign(d, ft) + m / t;
        }
      }
      else {
        t = (m / (s + t) + m / (r + l)) * (one + a);
      }
      l = fem::sqrt(t * t + four);
      crt = two / l;
      srt = t / l;
      clt = (crt + srt * m) / a;
      slt = (ht / ft) * srt / a;
    }
  }
  if (swap) {
    csl = srt;
    snl = crt;
    csr = slt;
    snr = clt;
  }
  else {
    csl = clt;
    snl = slt;
    csr = crt;
    snr = srt;
  }
  //C
  //C     Correct signs of SSMAX and SSMIN
  //C
  double tsign = fem::double0;
  if (pmax == 1) {
    tsign = fem::sign(one, csr) * fem::sign(one, csl) * fem::sign(one, f);
  }
  if (pmax == 2) {
    tsign = fem::sign(one, snr) * fem::sign(one, csl) * fem::sign(one, g);
  }
  if (pmax == 3) {
    tsign = fem::sign(one, snr) * fem::sign(one, snl) * fem::sign(one, h);
  }
  ssmax = fem::sign(ssmax, tsign);
  ssmin = fem::sign(ssmin, tsign * fem::sign(one, f) * fem::sign(one, h));
  //C
  //C     End of DLASV2
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dbdsqr.f
inline
void
dbdsqr(
  common& cmn,
  str_cref uplo,
  int const& n,
  int const& ncvt,
  int const& nru,
  int const& ncc,
  arr_ref<double> d,
  arr_ref<double> e,
  arr_ref<double, 2> vt,
  int const& ldvt,
  arr_ref<double, 2> u,
  int const& ldu,
  arr_ref<double, 2> c,
  int const& ldc,
  arr_ref<double> work,
  int& info)
{
  d(dimension(star));
  e(dimension(star));
  vt(dimension(ldvt, star));
  u(dimension(ldu, star));
  c(dimension(ldc, star));
  work(dimension(star));
  bool lower = fem::bool0;
  bool rotate = fem::bool0;
  int nm1 = fem::int0;
  int nm12 = fem::int0;
  int nm13 = fem::int0;
  int idir = fem::int0;
  double eps = fem::double0;
  double unfl = fem::double0;
  int i = fem::int0;
  double cs = fem::double0;
  double sn = fem::double0;
  double r = fem::double0;
  const double ten = 10.0e0;
  const double hndrd = 100.0e0;
  const double meigth = -0.125e0;
  double tolmul = fem::double0;
  double tol = fem::double0;
  const double zero = 0.0e0;
  double smax = fem::double0;
  double sminl = fem::double0;
  double sminoa = fem::double0;
  double mu = fem::double0;
  const int maxitr = 6;
  double thresh = fem::double0;
  int maxit = fem::int0;
  int iter = fem::int0;
  int oldll = fem::int0;
  int oldm = fem::int0;
  int m = fem::int0;
  double smin = fem::double0;
  int lll = fem::int0;
  int ll = fem::int0;
  double abss = fem::double0;
  double abse = fem::double0;
  double sigmn = fem::double0;
  double sigmx = fem::double0;
  double sinr = fem::double0;
  double cosr = fem::double0;
  double sinl = fem::double0;
  double cosl = fem::double0;
  const double hndrth = 0.01e0;
  double shift = fem::double0;
  double sll = fem::double0;
  const double one = 1.0e0;
  double oldcs = fem::double0;
  double oldsn = fem::double0;
  double h = fem::double0;
  double f = fem::double0;
  double g = fem::double0;
  const double negone = -1.0e0;
  int isub = fem::int0;
  int j = fem::int0;
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     January 2007
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DBDSQR computes the singular values and, optionally, the right and/or
  //C  left singular vectors from the singular value decomposition (SVD) of
  //C  a real N-by-N (upper or lower) bidiagonal matrix B using the implicit
  //C  zero-shift QR algorithm.  The SVD of B has the form
  //C
  //C     B = Q * S * P**T
  //C
  //C  where S is the diagonal matrix of singular values, Q is an orthogonal
  //C  matrix of left singular vectors, and P is an orthogonal matrix of
  //C  right singular vectors.  If left singular vectors are requested, this
  //C  subroutine actually returns U*Q instead of Q, and, if right singular
  //C  vectors are requested, this subroutine returns P**T*VT instead of
  //C  P**T, for given real input matrices U and VT.  When U and VT are the
  //C  orthogonal matrices that reduce a general matrix A to bidiagonal
  //C  form:  A = U*B*VT, as computed by DGEBRD, then
  //C
  //C     A = (U*Q) * S * (P**T*VT)
  //C
  //C  is the SVD of A.  Optionally, the subroutine may also compute Q**T*C
  //C  for a given real input matrix C.
  //C
  //C  See "Computing  Small Singular Values of Bidiagonal Matrices With
  //C  Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
  //C  LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
  //C  no. 5, pp. 873-912, Sept 1990) and
  //C  "Accurate singular values and differential qd algorithms," by
  //C  B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
  //C  Department, University of California at Berkeley, July 1992
  //C  for a detailed description of the algorithm.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  UPLO    (input) CHARACTER*1
  //C          = 'U':  B is upper bidiagonal;
  //C          = 'L':  B is lower bidiagonal.
  //C
  //C  N       (input) INTEGER
  //C          The order of the matrix B.  N >= 0.
  //C
  //C  NCVT    (input) INTEGER
  //C          The number of columns of the matrix VT. NCVT >= 0.
  //C
  //C  NRU     (input) INTEGER
  //C          The number of rows of the matrix U. NRU >= 0.
  //C
  //C  NCC     (input) INTEGER
  //C          The number of columns of the matrix C. NCC >= 0.
  //C
  //C  D       (input/output) DOUBLE PRECISION array, dimension (N)
  //C          On entry, the n diagonal elements of the bidiagonal matrix B.
  //C          On exit, if INFO=0, the singular values of B in decreasing
  //C          order.
  //C
  //C  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
  //C          On entry, the N-1 offdiagonal elements of the bidiagonal
  //C          matrix B.
  //C          On exit, if INFO = 0, E is destroyed; if INFO > 0, D and E
  //C          will contain the diagonal and superdiagonal elements of a
  //C          bidiagonal matrix orthogonally equivalent to the one given
  //C          as input.
  //C
  //C  VT      (input/output) DOUBLE PRECISION array, dimension (LDVT, NCVT)
  //C          On entry, an N-by-NCVT matrix VT.
  //C          On exit, VT is overwritten by P**T * VT.
  //C          Not referenced if NCVT = 0.
  //C
  //C  LDVT    (input) INTEGER
  //C          The leading dimension of the array VT.
  //C          LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.
  //C
  //C  U       (input/output) DOUBLE PRECISION array, dimension (LDU, N)
  //C          On entry, an NRU-by-N matrix U.
  //C          On exit, U is overwritten by U * Q.
  //C          Not referenced if NRU = 0.
  //C
  //C  LDU     (input) INTEGER
  //C          The leading dimension of the array U.  LDU >= max(1,NRU).
  //C
  //C  C       (input/output) DOUBLE PRECISION array, dimension (LDC, NCC)
  //C          On entry, an N-by-NCC matrix C.
  //C          On exit, C is overwritten by Q**T * C.
  //C          Not referenced if NCC = 0.
  //C
  //C  LDC     (input) INTEGER
  //C          The leading dimension of the array C.
  //C          LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0.
  //C
  //C  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  If INFO = -i, the i-th argument had an illegal value
  //C          > 0:
  //C             if NCVT = NRU = NCC = 0,
  //C                = 1, a split was marked by a positive value in E
  //C                = 2, current block of Z not diagonalized after 30*N
  //C                     iterations (in inner while loop)
  //C                = 3, termination criterion of outer while loop not met
  //C                     (program created more than N unreduced blocks)
  //C             else NCVT = NRU = NCC = 0,
  //C                   the algorithm did not converge; D and E contain the
  //C                   elements of a bidiagonal matrix which is orthogonally
  //C                   similar to the input matrix B;  if INFO = i, i
  //C                   elements of E have not converged to zero.
  //C
  //C  Internal Parameters
  //C  ===================
  //C
  //C  TOLMUL  DOUBLE PRECISION, default = max(10,min(100,EPS**(-1/8)))
  //C          TOLMUL controls the convergence criterion of the QR loop.
  //C          If it is positive, TOLMUL*EPS is the desired relative
  //C             precision in the computed singular values.
  //C          If it is negative, abs(TOLMUL*EPS*sigma_max) is the
  //C             desired absolute accuracy in the computed singular
  //C             values (corresponds to relative accuracy
  //C             abs(TOLMUL*EPS) in the largest singular value.
  //C          abs(TOLMUL) should be between 1 and 1/EPS, and preferably
  //C             between 10 (for fast convergence) and .1/EPS
  //C             (for there to be some accuracy in the results).
  //C          Default is to lose at either one eighth or 2 of the
  //C             available decimal digits in each computed singular value
  //C             (whichever is smaller).
  //C
  //C  MAXITR  INTEGER, default = 6
  //C          MAXITR controls the maximum number of passes of the
  //C          algorithm through its inner loop. The algorithms stops
  //C          (and so fails to converge) if the number of passes
  //C          through the inner loop exceeds MAXITR*N**2.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters.
  //C
  info = 0;
  lower = lsame(uplo, "L");
  if (!lsame(uplo, "U") && !lower) {
    info = -1;
  }
  else if (n < 0) {
    info = -2;
  }
  else if (ncvt < 0) {
    info = -3;
  }
  else if (nru < 0) {
    info = -4;
  }
  else if (ncc < 0) {
    info = -5;
  }
  else if ((ncvt == 0 && ldvt < 1) || (ncvt > 0 && ldvt < fem::max(1, n))) {
    info = -9;
  }
  else if (ldu < fem::max(1, nru)) {
    info = -11;
  }
  else if ((ncc == 0 && ldc < 1) || (ncc > 0 && ldc < fem::max(1, n))) {
    info = -13;
  }
  if (info != 0) {
    xerbla("DBDSQR", -info);
    return;
  }
  if (n == 0) {
    return;
  }
  if (n == 1) {
    goto statement_160;
  }
  //C
  //C     ROTATE is true if any singular vectors desired, false otherwise
  //C
  rotate = (ncvt > 0) || (nru > 0) || (ncc > 0);
  //C
  //C     If no singular vectors desired, use qd algorithm
  //C
  if (!rotate) {
    dlasq1(cmn, n, d, e, work, info);
    return;
  }
  //C
  nm1 = n - 1;
  nm12 = nm1 + nm1;
  nm13 = nm12 + nm1;
  idir = 0;
  //C
  //C     Get machine constants
  //C
  eps = dlamch(cmn, "Epsilon");
  unfl = dlamch(cmn, "Safe minimum");
  //C
  //C     If matrix lower bidiagonal, rotate to be upper bidiagonal
  //C     by applying Givens rotations on the left
  //C
  if (lower) {
    {
      int fem_do_last = n - 1;
      FEM_DO(i, 1, fem_do_last) {
        dlartg(cmn, d(i), e(i), cs, sn, r);
        d(i) = r;
        e(i) = sn * d(i + 1);
        d(i + 1) = cs * d(i + 1);
        work(i) = cs;
        work(nm1 + i) = sn;
      }
    }
    //C
    //C        Update singular vectors if desired
    //C
    if (nru > 0) {
      dlasr("R", "V", "F", nru, n, work(1), work(n), u, ldu);
    }
    if (ncc > 0) {
      dlasr("L", "V", "F", n, ncc, work(1), work(n), c, ldc);
    }
  }
  //C
  //C     Compute singular values to relative accuracy TOL
  //C     (By setting TOL to be negative, algorithm will compute
  //C     singular values to absolute accuracy ABS(TOL)*norm(input matrix))
  //C
  tolmul = fem::max(ten, fem::min(hndrd, fem::pow(eps, meigth)));
  tol = tolmul * eps;
  //C
  //C     Compute approximate maximum, minimum singular values
  //C
  smax = zero;
  FEM_DO(i, 1, n) {
    smax = fem::max(smax, fem::abs(d(i)));
  }
  {
    int fem_do_last = n - 1;
    FEM_DO(i, 1, fem_do_last) {
      smax = fem::max(smax, fem::abs(e(i)));
    }
  }
  sminl = zero;
  if (tol >= zero) {
    //C
    //C        Relative accuracy desired
    //C
    sminoa = fem::abs(d(1));
    if (sminoa == zero) {
      goto statement_50;
    }
    mu = sminoa;
    FEM_DO(i, 2, n) {
      mu = fem::abs(d(i)) * (mu / (mu + fem::abs(e(i - 1))));
      sminoa = fem::min(sminoa, mu);
      if (sminoa == zero) {
        goto statement_50;
      }
    }
    statement_50:
    sminoa = sminoa / fem::sqrt(fem::dble(n));
    thresh = fem::max(tol * sminoa, maxitr * n * n * unfl);
  }
  else {
    //C
    //C        Absolute accuracy desired
    //C
    thresh = fem::max(fem::abs(tol) * smax, maxitr * n * n * unfl);
  }
  //C
  //C     Prepare for main iteration loop for the singular values
  //C     (MAXIT is the maximum number of passes through the inner
  //C     loop permitted before nonconvergence signalled.)
  //C
  maxit = maxitr * n * n;
  iter = 0;
  oldll = -1;
  oldm = -1;
  //C
  //C     M points to last element of unconverged part of matrix
  //C
  m = n;
  //C
  //C     Begin main iteration loop
  //C
  statement_60:
  //C
  //C     Check for convergence or exceeding iteration count
  //C
  if (m <= 1) {
    goto statement_160;
  }
  if (iter > maxit) {
    goto statement_200;
  }
  //C
  //C     Find diagonal block of matrix to work on
  //C
  if (tol < zero && fem::abs(d(m)) <= thresh) {
    d(m) = zero;
  }
  smax = fem::abs(d(m));
  smin = smax;
  {
    int fem_do_last = m - 1;
    FEM_DO(lll, 1, fem_do_last) {
      ll = m - lll;
      abss = fem::abs(d(ll));
      abse = fem::abs(e(ll));
      if (tol < zero && abss <= thresh) {
        d(ll) = zero;
      }
      if (abse <= thresh) {
        goto statement_80;
      }
      smin = fem::min(smin, abss);
      smax = fem::max(smax, abss, abse);
    }
  }
  ll = 0;
  goto statement_90;
  statement_80:
  e(ll) = zero;
  //C
  //C     Matrix splits since E(LL) = 0
  //C
  if (ll == m - 1) {
    //C
    //C        Convergence of bottom singular value, return to top of loop
    //C
    m = m - 1;
    goto statement_60;
  }
  statement_90:
  ll++;
  //C
  //C     E(LL) through E(M-1) are nonzero, E(LL-1) is zero
  //C
  if (ll == m - 1) {
    //C
    //C        2 by 2 block, handle separately
    //C
    dlasv2(cmn, d(m - 1), e(m - 1), d(m), sigmn, sigmx, sinr, cosr, sinl, cosl);
    d(m - 1) = sigmx;
    e(m - 1) = zero;
    d(m) = sigmn;
    //C
    //C        Compute singular vectors, if desired
    //C
    if (ncvt > 0) {
      drot(ncvt, vt(m - 1, 1), ldvt, vt(m, 1), ldvt, cosr, sinr);
    }
    if (nru > 0) {
      drot(nru, u(1, m - 1), 1, u(1, m), 1, cosl, sinl);
    }
    if (ncc > 0) {
      drot(ncc, c(m - 1, 1), ldc, c(m, 1), ldc, cosl, sinl);
    }
    m = m - 2;
    goto statement_60;
  }
  //C
  //C     If working on new submatrix, choose shift direction
  //C     (from larger end diagonal element towards smaller)
  //C
  if (ll > oldm || m < oldll) {
    if (fem::abs(d(ll)) >= fem::abs(d(m))) {
      //C
      //C           Chase bulge from top (big end) to bottom (small end)
      //C
      idir = 1;
    }
    else {
      //C
      //C           Chase bulge from bottom (big end) to top (small end)
      //C
      idir = 2;
    }
  }
  //C
  //C     Apply convergence tests
  //C
  if (idir == 1) {
    //C
    //C        Run convergence test in forward direction
    //C        First apply standard test to bottom of matrix
    //C
    if (fem::abs(e(m - 1)) <= fem::abs(tol) * fem::abs(d(m)) || (
        tol < zero && fem::abs(e(m - 1)) <= thresh)) {
      e(m - 1) = zero;
      goto statement_60;
    }
    //C
    if (tol >= zero) {
      //C
      //C           If relative accuracy desired,
      //C           apply convergence criterion forward
      //C
      mu = fem::abs(d(ll));
      sminl = mu;
      {
        int fem_do_last = m - 1;
        FEM_DO(lll, ll, fem_do_last) {
          if (fem::abs(e(lll)) <= tol * mu) {
            e(lll) = zero;
            goto statement_60;
          }
          mu = fem::abs(d(lll + 1)) * (mu / (mu + fem::abs(e(lll))));
          sminl = fem::min(sminl, mu);
        }
      }
    }
    //C
  }
  else {
    //C
    //C        Run convergence test in backward direction
    //C        First apply standard test to top of matrix
    //C
    if (fem::abs(e(ll)) <= fem::abs(tol) * fem::abs(d(ll)) || (
        tol < zero && fem::abs(e(ll)) <= thresh)) {
      e(ll) = zero;
      goto statement_60;
    }
    //C
    if (tol >= zero) {
      //C
      //C           If relative accuracy desired,
      //C           apply convergence criterion backward
      //C
      mu = fem::abs(d(m));
      sminl = mu;
      FEM_DOSTEP(lll, m - 1, ll, -1) {
        if (fem::abs(e(lll)) <= tol * mu) {
          e(lll) = zero;
          goto statement_60;
        }
        mu = fem::abs(d(lll)) * (mu / (mu + fem::abs(e(lll))));
        sminl = fem::min(sminl, mu);
      }
    }
  }
  oldll = ll;
  oldm = m;
  //C
  //C     Compute shift.  First, test if shifting would ruin relative
  //C     accuracy, and if so set the shift to zero.
  //C
  if (tol >= zero && n * tol * (sminl / smax) <= fem::max(eps, hndrth * tol)) {
    //C
    //C        Use a zero shift to avoid loss of relative accuracy
    //C
    shift = zero;
  }
  else {
    //C
    //C        Compute the shift from 2-by-2 block at end of matrix
    //C
    if (idir == 1) {
      sll = fem::abs(d(ll));
      dlas2(d(m - 1), e(m - 1), d(m), shift, r);
    }
    else {
      sll = fem::abs(d(m));
      dlas2(d(ll), e(ll), d(ll + 1), shift, r);
    }
    //C
    //C        Test if shift negligible, and if so set to zero
    //C
    if (sll > zero) {
      if (fem::pow2((shift / sll)) < eps) {
        shift = zero;
      }
    }
  }
  //C
  //C     Increment iteration count
  //C
  iter += m - ll;
  //C
  //C     If SHIFT = 0, do simplified QR iteration
  //C
  if (shift == zero) {
    if (idir == 1) {
      //C
      //C           Chase bulge from top to bottom
      //C           Save cosines and sines for later singular vector updates
      //C
      cs = one;
      oldcs = one;
      {
        int fem_do_last = m - 1;
        FEM_DO(i, ll, fem_do_last) {
          dlartg(cmn, d(i) * cs, e(i), cs, sn, r);
          if (i > ll) {
            e(i - 1) = oldsn * r;
          }
          dlartg(cmn, oldcs * r, d(i + 1) * sn, oldcs, oldsn, d(i));
          work(i - ll + 1) = cs;
          work(i - ll + 1 + nm1) = sn;
          work(i - ll + 1 + nm12) = oldcs;
          work(i - ll + 1 + nm13) = oldsn;
        }
      }
      h = d(m) * cs;
      d(m) = h * oldcs;
      e(m - 1) = h * oldsn;
      //C
      //C           Update singular vectors
      //C
      if (ncvt > 0) {
        dlasr("L", "V", "F", m - ll + 1, ncvt, work(1), work(n), vt(ll,
          1), ldvt);
      }
      if (nru > 0) {
        dlasr("R", "V", "F", nru, m - ll + 1, work(nm12 + 1), work(nm13 + 1),
          u(1, ll), ldu);
      }
      if (ncc > 0) {
        dlasr("L", "V", "F", m - ll + 1, ncc, work(nm12 + 1), work(nm13 + 1),
          c(ll, 1), ldc);
      }
      //C
      //C           Test convergence
      //C
      if (fem::abs(e(m - 1)) <= thresh) {
        e(m - 1) = zero;
      }
      //C
    }
    else {
      //C
      //C           Chase bulge from bottom to top
      //C           Save cosines and sines for later singular vector updates
      //C
      cs = one;
      oldcs = one;
      FEM_DOSTEP(i, m, ll + 1, -1) {
        dlartg(cmn, d(i) * cs, e(i - 1), cs, sn, r);
        if (i < m) {
          e(i) = oldsn * r;
        }
        dlartg(cmn, oldcs * r, d(i - 1) * sn, oldcs, oldsn, d(i));
        work(i - ll) = cs;
        work(i - ll + nm1) = -sn;
        work(i - ll + nm12) = oldcs;
        work(i - ll + nm13) = -oldsn;
      }
      h = d(ll) * cs;
      d(ll) = h * oldcs;
      e(ll) = h * oldsn;
      //C
      //C           Update singular vectors
      //C
      if (ncvt > 0) {
        dlasr("L", "V", "B", m - ll + 1, ncvt, work(nm12 + 1), work(nm13 + 1),
          vt(ll, 1), ldvt);
      }
      if (nru > 0) {
        dlasr("R", "V", "B", nru, m - ll + 1, work(1), work(n), u(1, ll), ldu);
      }
      if (ncc > 0) {
        dlasr("L", "V", "B", m - ll + 1, ncc, work(1), work(n), c(ll, 1), ldc);
      }
      //C
      //C           Test convergence
      //C
      if (fem::abs(e(ll)) <= thresh) {
        e(ll) = zero;
      }
    }
  }
  else {
    //C
    //C        Use nonzero shift
    //C
    if (idir == 1) {
      //C
      //C           Chase bulge from top to bottom
      //C           Save cosines and sines for later singular vector updates
      //C
      f = (fem::abs(d(ll)) - shift) * (fem::sign(one, d(ll)) + shift / d(ll));
      g = e(ll);
      {
        int fem_do_last = m - 1;
        FEM_DO(i, ll, fem_do_last) {
          dlartg(cmn, f, g, cosr, sinr, r);
          if (i > ll) {
            e(i - 1) = r;
          }
          f = cosr * d(i) + sinr * e(i);
          e(i) = cosr * e(i) - sinr * d(i);
          g = sinr * d(i + 1);
          d(i + 1) = cosr * d(i + 1);
          dlartg(cmn, f, g, cosl, sinl, r);
          d(i) = r;
          f = cosl * e(i) + sinl * d(i + 1);
          d(i + 1) = cosl * d(i + 1) - sinl * e(i);
          if (i < m - 1) {
            g = sinl * e(i + 1);
            e(i + 1) = cosl * e(i + 1);
          }
          work(i - ll + 1) = cosr;
          work(i - ll + 1 + nm1) = sinr;
          work(i - ll + 1 + nm12) = cosl;
          work(i - ll + 1 + nm13) = sinl;
        }
      }
      e(m - 1) = f;
      //C
      //C           Update singular vectors
      //C
      if (ncvt > 0) {
        dlasr("L", "V", "F", m - ll + 1, ncvt, work(1), work(n), vt(ll,
          1), ldvt);
      }
      if (nru > 0) {
        dlasr("R", "V", "F", nru, m - ll + 1, work(nm12 + 1), work(nm13 + 1),
          u(1, ll), ldu);
      }
      if (ncc > 0) {
        dlasr("L", "V", "F", m - ll + 1, ncc, work(nm12 + 1), work(nm13 + 1),
          c(ll, 1), ldc);
      }
      //C
      //C           Test convergence
      //C
      if (fem::abs(e(m - 1)) <= thresh) {
        e(m - 1) = zero;
      }
      //C
    }
    else {
      //C
      //C           Chase bulge from bottom to top
      //C           Save cosines and sines for later singular vector updates
      //C
      f = (fem::abs(d(m)) - shift) * (fem::sign(one, d(m)) + shift / d(m));
      g = e(m - 1);
      FEM_DOSTEP(i, m, ll + 1, -1) {
        dlartg(cmn, f, g, cosr, sinr, r);
        if (i < m) {
          e(i) = r;
        }
        f = cosr * d(i) + sinr * e(i - 1);
        e(i - 1) = cosr * e(i - 1) - sinr * d(i);
        g = sinr * d(i - 1);
        d(i - 1) = cosr * d(i - 1);
        dlartg(cmn, f, g, cosl, sinl, r);
        d(i) = r;
        f = cosl * e(i - 1) + sinl * d(i - 1);
        d(i - 1) = cosl * d(i - 1) - sinl * e(i - 1);
        if (i > ll + 1) {
          g = sinl * e(i - 2);
          e(i - 2) = cosl * e(i - 2);
        }
        work(i - ll) = cosr;
        work(i - ll + nm1) = -sinr;
        work(i - ll + nm12) = cosl;
        work(i - ll + nm13) = -sinl;
      }
      e(ll) = f;
      //C
      //C           Test convergence
      //C
      if (fem::abs(e(ll)) <= thresh) {
        e(ll) = zero;
      }
      //C
      //C           Update singular vectors if desired
      //C
      if (ncvt > 0) {
        dlasr("L", "V", "B", m - ll + 1, ncvt, work(nm12 + 1), work(nm13 + 1),
          vt(ll, 1), ldvt);
      }
      if (nru > 0) {
        dlasr("R", "V", "B", nru, m - ll + 1, work(1), work(n), u(1, ll), ldu);
      }
      if (ncc > 0) {
        dlasr("L", "V", "B", m - ll + 1, ncc, work(1), work(n), c(ll, 1), ldc);
      }
    }
  }
  //C
  //C     QR iteration finished, go back and check convergence
  //C
  goto statement_60;
  //C
  //C     All singular values converged, so make them positive
  //C
  statement_160:
  FEM_DO(i, 1, n) {
    if (d(i) < zero) {
      d(i) = -d(i);
      //C
      //C           Change sign of singular vectors, if desired
      //C
      if (ncvt > 0) {
        dscal(ncvt, negone, vt(i, 1), ldvt);
      }
    }
  }
  //C
  //C     Sort the singular values into decreasing order (insertion sort on
  //C     singular values, but only one transposition per singular vector)
  //C
  {
    int fem_do_last = n - 1;
    FEM_DO(i, 1, fem_do_last) {
      //C
      //C        Scan for smallest D(I)
      //C
      isub = 1;
      smin = d(1);
      {
        int fem_do_last = n + 1 - i;
        FEM_DO(j, 2, fem_do_last) {
          if (d(j) <= smin) {
            isub = j;
            smin = d(j);
          }
        }
      }
      if (isub != n + 1 - i) {
        //C
        //C           Swap singular values and vectors
        //C
        d(isub) = d(n + 1 - i);
        d(n + 1 - i) = smin;
        if (ncvt > 0) {
          dswap(ncvt, vt(isub, 1), ldvt, vt(n + 1 - i, 1), ldvt);
        }
        if (nru > 0) {
          dswap(nru, u(1, isub), 1, u(1, n + 1 - i), 1);
        }
        if (ncc > 0) {
          dswap(ncc, c(isub, 1), ldc, c(n + 1 - i, 1), ldc);
        }
      }
    }
  }
  goto statement_220;
  //C
  //C     Maximum number of iterations exceeded, failure to converge
  //C
  statement_200:
  info = 0;
  {
    int fem_do_last = n - 1;
    FEM_DO(i, 1, fem_do_last) {
      if (e(i) != zero) {
        info++;
      }
    }
  }
  statement_220:;
  //C
  //C     End of DBDSQR
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasdq.f
inline
void
dlasdq(
  common& cmn,
  str_cref uplo,
  int const& sqre,
  int const& n,
  int const& ncvt,
  int const& nru,
  int const& ncc,
  arr_ref<double> d,
  arr_ref<double> e,
  arr_ref<double, 2> vt,
  int const& ldvt,
  arr_ref<double, 2> u,
  int const& ldu,
  arr_ref<double, 2> c,
  int const& ldc,
  arr_ref<double> work,
  int& info)
{
  d(dimension(star));
  e(dimension(star));
  vt(dimension(ldvt, star));
  u(dimension(ldu, star));
  c(dimension(ldc, star));
  work(dimension(star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASDQ computes the singular value decomposition (SVD) of a real
  //C  (upper or lower) bidiagonal matrix with diagonal D and offdiagonal
  //C  E, accumulating the transformations if desired. Letting B denote
  //C  the input bidiagonal matrix, the algorithm computes orthogonal
  //C  matrices Q and P such that B = Q * S * P' (P' denotes the transpose
  //C  of P). The singular values S are overwritten on D.
  //C
  //C  The input matrix U  is changed to U  * Q  if desired.
  //C  The input matrix VT is changed to P' * VT if desired.
  //C  The input matrix C  is changed to Q' * C  if desired.
  //C
  //C  See "Computing  Small Singular Values of Bidiagonal Matrices With
  //C  Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
  //C  LAPACK Working Note #3, for a detailed description of the algorithm.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  UPLO  (input) CHARACTER*1
  //C        On entry, UPLO specifies whether the input bidiagonal matrix
  //C        is upper or lower bidiagonal, and wether it is square are
  //C        not.
  //C           UPLO = 'U' or 'u'   B is upper bidiagonal.
  //C           UPLO = 'L' or 'l'   B is lower bidiagonal.
  //C
  //C  SQRE  (input) INTEGER
  //C        = 0: then the input matrix is N-by-N.
  //C        = 1: then the input matrix is N-by-(N+1) if UPLU = 'U' and
  //C             (N+1)-by-N if UPLU = 'L'.
  //C
  //C        The bidiagonal matrix has
  //C        N = NL + NR + 1 rows and
  //C        M = N + SQRE >= N columns.
  //C
  //C  N     (input) INTEGER
  //C        On entry, N specifies the number of rows and columns
  //C        in the matrix. N must be at least 0.
  //C
  //C  NCVT  (input) INTEGER
  //C        On entry, NCVT specifies the number of columns of
  //C        the matrix VT. NCVT must be at least 0.
  //C
  //C  NRU   (input) INTEGER
  //C        On entry, NRU specifies the number of rows of
  //C        the matrix U. NRU must be at least 0.
  //C
  //C  NCC   (input) INTEGER
  //C        On entry, NCC specifies the number of columns of
  //C        the matrix C. NCC must be at least 0.
  //C
  //C  D     (input/output) DOUBLE PRECISION array, dimension (N)
  //C        On entry, D contains the diagonal entries of the
  //C        bidiagonal matrix whose SVD is desired. On normal exit,
  //C        D contains the singular values in ascending order.
  //C
  //C  E     (input/output) DOUBLE PRECISION array.
  //C        dimension is (N-1) if SQRE = 0 and N if SQRE = 1.
  //C        On entry, the entries of E contain the offdiagonal entries
  //C        of the bidiagonal matrix whose SVD is desired. On normal
  //C        exit, E will contain 0. If the algorithm does not converge,
  //C        D and E will contain the diagonal and superdiagonal entries
  //C        of a bidiagonal matrix orthogonally equivalent to the one
  //C        given as input.
  //C
  //C  VT    (input/output) DOUBLE PRECISION array, dimension (LDVT, NCVT)
  //C        On entry, contains a matrix which on exit has been
  //C        premultiplied by P', dimension N-by-NCVT if SQRE = 0
  //C        and (N+1)-by-NCVT if SQRE = 1 (not referenced if NCVT=0).
  //C
  //C  LDVT  (input) INTEGER
  //C        On entry, LDVT specifies the leading dimension of VT as
  //C        declared in the calling (sub) program. LDVT must be at
  //C        least 1. If NCVT is nonzero LDVT must also be at least N.
  //C
  //C  U     (input/output) DOUBLE PRECISION array, dimension (LDU, N)
  //C        On entry, contains a  matrix which on exit has been
  //C        postmultiplied by Q, dimension NRU-by-N if SQRE = 0
  //C        and NRU-by-(N+1) if SQRE = 1 (not referenced if NRU=0).
  //C
  //C  LDU   (input) INTEGER
  //C        On entry, LDU  specifies the leading dimension of U as
  //C        declared in the calling (sub) program. LDU must be at
  //C        least max( 1, NRU ) .
  //C
  //C  C     (input/output) DOUBLE PRECISION array, dimension (LDC, NCC)
  //C        On entry, contains an N-by-NCC matrix which on exit
  //C        has been premultiplied by Q'  dimension N-by-NCC if SQRE = 0
  //C        and (N+1)-by-NCC if SQRE = 1 (not referenced if NCC=0).
  //C
  //C  LDC   (input) INTEGER
  //C        On entry, LDC  specifies the leading dimension of C as
  //C        declared in the calling (sub) program. LDC must be at
  //C        least 1. If NCC is nonzero, LDC must also be at least N.
  //C
  //C  WORK  (workspace) DOUBLE PRECISION array, dimension (4*N)
  //C        Workspace. Only referenced if one of NCVT, NRU, or NCC is
  //C        nonzero, and if N is at least 2.
  //C
  //C  INFO  (output) INTEGER
  //C        On exit, a value of 0 indicates a successful exit.
  //C        If INFO < 0, argument number -INFO is illegal.
  //C        If INFO > 0, the algorithm did not converge, and INFO
  //C        specifies how many superdiagonals did not converge.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Based on contributions by
  //C     Ming Gu and Huan Ren, Computer Science Division, University of
  //C     California at Berkeley, USA
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters.
  //C
  info = 0;
  int iuplo = 0;
  if (lsame(uplo, "U")) {
    iuplo = 1;
  }
  if (lsame(uplo, "L")) {
    iuplo = 2;
  }
  if (iuplo == 0) {
    info = -1;
  }
  else if ((sqre < 0) || (sqre > 1)) {
    info = -2;
  }
  else if (n < 0) {
    info = -3;
  }
  else if (ncvt < 0) {
    info = -4;
  }
  else if (nru < 0) {
    info = -5;
  }
  else if (ncc < 0) {
    info = -6;
  }
  else if ((ncvt == 0 && ldvt < 1) || (ncvt > 0 && ldvt < fem::max(1, n))) {
    info = -10;
  }
  else if (ldu < fem::max(1, nru)) {
    info = -12;
  }
  else if ((ncc == 0 && ldc < 1) || (ncc > 0 && ldc < fem::max(1, n))) {
    info = -14;
  }
  if (info != 0) {
    xerbla("DLASDQ", -info);
    return;
  }
  if (n == 0) {
    return;
  }
  //C
  //C     ROTATE is true if any singular vectors desired, false otherwise
  //C
  bool rotate = (ncvt > 0) || (nru > 0) || (ncc > 0);
  int np1 = n + 1;
  int sqre1 = sqre;
  //C
  //C     If matrix non-square upper bidiagonal, rotate to be lower
  //C     bidiagonal.  The rotations are on the right.
  //C
  int i = fem::int0;
  double cs = fem::double0;
  double sn = fem::double0;
  double r = fem::double0;
  const double zero = 0.0e+0;
  if ((iuplo == 1) && (sqre1 == 1)) {
    {
      int fem_do_last = n - 1;
      FEM_DO(i, 1, fem_do_last) {
        dlartg(cmn, d(i), e(i), cs, sn, r);
        d(i) = r;
        e(i) = sn * d(i + 1);
        d(i + 1) = cs * d(i + 1);
        if (rotate) {
          work(i) = cs;
          work(n + i) = sn;
        }
      }
    }
    dlartg(cmn, d(n), e(n), cs, sn, r);
    d(n) = r;
    e(n) = zero;
    if (rotate) {
      work(n) = cs;
      work(n + n) = sn;
    }
    iuplo = 2;
    sqre1 = 0;
    //C
    //C        Update singular vectors if desired.
    //C
    if (ncvt > 0) {
      dlasr("L", "V", "F", np1, ncvt, work(1), work(np1), vt, ldvt);
    }
  }
  //C
  //C     If matrix lower bidiagonal, rotate to be upper bidiagonal
  //C     by applying Givens rotations on the left.
  //C
  if (iuplo == 2) {
    {
      int fem_do_last = n - 1;
      FEM_DO(i, 1, fem_do_last) {
        dlartg(cmn, d(i), e(i), cs, sn, r);
        d(i) = r;
        e(i) = sn * d(i + 1);
        d(i + 1) = cs * d(i + 1);
        if (rotate) {
          work(i) = cs;
          work(n + i) = sn;
        }
      }
    }
    //C
    //C        If matrix (N+1)-by-N lower bidiagonal, one additional
    //C        rotation is needed.
    //C
    if (sqre1 == 1) {
      dlartg(cmn, d(n), e(n), cs, sn, r);
      d(n) = r;
      if (rotate) {
        work(n) = cs;
        work(n + n) = sn;
      }
    }
    //C
    //C        Update singular vectors if desired.
    //C
    if (nru > 0) {
      if (sqre1 == 0) {
        dlasr("R", "V", "F", nru, n, work(1), work(np1), u, ldu);
      }
      else {
        dlasr("R", "V", "F", nru, np1, work(1), work(np1), u, ldu);
      }
    }
    if (ncc > 0) {
      if (sqre1 == 0) {
        dlasr("L", "V", "F", n, ncc, work(1), work(np1), c, ldc);
      }
      else {
        dlasr("L", "V", "F", np1, ncc, work(1), work(np1), c, ldc);
      }
    }
  }
  //C
  //C     Call DBDSQR to compute the SVD of the reduced real
  //C     N-by-N upper bidiagonal matrix.
  //C
  dbdsqr(cmn, "U", n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc,
    work, info);
  //C
  //C     Sort the singular values into ascending order (insertion sort on
  //C     singular values, but only one transposition per singular vector)
  //C
  int isub = fem::int0;
  double smin = fem::double0;
  int j = fem::int0;
  FEM_DO(i, 1, n) {
    //C
    //C        Scan for smallest D(I).
    //C
    isub = i;
    smin = d(i);
    FEM_DO(j, i + 1, n) {
      if (d(j) < smin) {
        isub = j;
        smin = d(j);
      }
    }
    if (isub != i) {
      //C
      //C           Swap singular values and vectors.
      //C
      d(isub) = d(i);
      d(i) = smin;
      if (ncvt > 0) {
        dswap(ncvt, vt(isub, 1), ldvt, vt(i, 1), ldvt);
      }
      if (nru > 0) {
        dswap(nru, u(1, isub), 1, u(1, i), 1);
      }
      if (ncc > 0) {
        dswap(ncc, c(isub, 1), ldc, c(i, 1), ldc);
      }
    }
  }
  //C
  //C     End of DLASDQ
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasdt.f
inline
void
dlasdt(
  int const& n,
  int& lvl,
  int& nd,
  arr_ref<int> inode,
  arr_ref<int> ndiml,
  arr_ref<int> ndimr,
  int const& msub)
{
  inode(dimension(star));
  ndiml(dimension(star));
  ndimr(dimension(star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASDT creates a tree of subproblems for bidiagonal divide and
  //C  conquer.
  //C
  //C  Arguments
  //C  =========
  //C
  //C   N      (input) INTEGER
  //C          On entry, the number of diagonal elements of the
  //C          bidiagonal matrix.
  //C
  //C   LVL    (output) INTEGER
  //C          On exit, the number of levels on the computation tree.
  //C
  //C   ND     (output) INTEGER
  //C          On exit, the number of nodes on the tree.
  //C
  //C   INODE  (output) INTEGER array, dimension ( N )
  //C          On exit, centers of subproblems.
  //C
  //C   NDIML  (output) INTEGER array, dimension ( N )
  //C          On exit, row dimensions of left children.
  //C
  //C   NDIMR  (output) INTEGER array, dimension ( N )
  //C          On exit, row dimensions of right children.
  //C
  //C   MSUB   (input) INTEGER.
  //C          On entry, the maximum row dimension each subproblem at the
  //C          bottom of the tree can be of.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Based on contributions by
  //C     Ming Gu and Huan Ren, Computer Science Division, University of
  //C     California at Berkeley, USA
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Find the number of levels on the tree.
  //C
  int maxn = fem::max(1, n);
  const double two = 2.0e+0;
  double temp = fem::log(fem::dble(maxn) / fem::dble(msub + 1)) / fem::log(two);
  lvl = fem::fint(temp) + 1;
  //C
  int i = n / 2;
  inode(1) = i + 1;
  ndiml(1) = i;
  ndimr(1) = n - i - 1;
  int il = 0;
  int ir = 1;
  int llst = 1;
  int nlvl = fem::int0;
  int ncrnt = fem::int0;
  {
    int fem_do_last = lvl - 1;
    FEM_DO(nlvl, 1, fem_do_last) {
      //C
      //C        Constructing the tree at (NLVL+1)-st level. The number of
      //C        nodes created on this level is LLST * 2.
      //C
      {
        int fem_do_last = llst - 1;
        FEM_DO(i, 0, fem_do_last) {
          il += 2;
          ir += 2;
          ncrnt = llst + i;
          ndiml(il) = ndiml(ncrnt) / 2;
          ndimr(il) = ndiml(ncrnt) - ndiml(il) - 1;
          inode(il) = inode(ncrnt) - ndimr(il) - 1;
          ndiml(ir) = ndimr(ncrnt) / 2;
          ndimr(ir) = ndimr(ncrnt) - ndiml(ir) - 1;
          inode(ir) = inode(ncrnt) + ndiml(ir) + 1;
        }
      }
      llst = llst * 2;
    }
  }
  nd = llst * 2 - 1;
  //C
  //C     End of DLASDT
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasd0.f
inline
void
dlasd0(
  common& cmn,
  int const& n,
  int const& sqre,
  arr_ref<double> d,
  arr_ref<double> e,
  arr_ref<double, 2> u,
  int const& ldu,
  arr_ref<double, 2> vt,
  int const& ldvt,
  int const& smlsiz,
  arr_ref<int> iwork,
  arr_ref<double> work,
  int& info)
{
  d(dimension(star));
  e(dimension(star));
  u(dimension(ldu, star));
  vt(dimension(ldvt, star));
  iwork(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  Using a divide and conquer approach, DLASD0 computes the singular
  //C  value decomposition (SVD) of a real upper bidiagonal N-by-M
  //C  matrix B with diagonal D and offdiagonal E, where M = N + SQRE.
  //C  The algorithm computes orthogonal matrices U and VT such that
  //C  B = U * S * VT. The singular values S are overwritten on D.
  //C
  //C  A related subroutine, DLASDA, computes only the singular values,
  //C  and optionally, the singular vectors in compact form.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  N      (input) INTEGER
  //C         On entry, the row dimension of the upper bidiagonal matrix.
  //C         This is also the dimension of the main diagonal array D.
  //C
  //C  SQRE   (input) INTEGER
  //C         Specifies the column dimension of the bidiagonal matrix.
  //C         = 0: The bidiagonal matrix has column dimension M = N;
  //C         = 1: The bidiagonal matrix has column dimension M = N+1;
  //C
  //C  D      (input/output) DOUBLE PRECISION array, dimension (N)
  //C         On entry D contains the main diagonal of the bidiagonal
  //C         matrix.
  //C         On exit D, if INFO = 0, contains its singular values.
  //C
  //C  E      (input) DOUBLE PRECISION array, dimension (M-1)
  //C         Contains the subdiagonal entries of the bidiagonal matrix.
  //C         On exit, E has been destroyed.
  //C
  //C  U      (output) DOUBLE PRECISION array, dimension at least (LDQ, N)
  //C         On exit, U contains the left singular vectors.
  //C
  //C  LDU    (input) INTEGER
  //C         On entry, leading dimension of U.
  //C
  //C  VT     (output) DOUBLE PRECISION array, dimension at least (LDVT, M)
  //C         On exit, VT' contains the right singular vectors.
  //C
  //C  LDVT   (input) INTEGER
  //C         On entry, leading dimension of VT.
  //C
  //C  SMLSIZ (input) INTEGER
  //C         On entry, maximum size of the subproblems at the
  //C         bottom of the computation tree.
  //C
  //C  IWORK  (workspace) INTEGER work array.
  //C         Dimension must be at least (8 * N)
  //C
  //C  WORK   (workspace) DOUBLE PRECISION work array.
  //C         Dimension must be at least (3 * M**2 + 2 * M)
  //C
  //C  INFO   (output) INTEGER
  //C          = 0:  successful exit.
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value.
  //C          > 0:  if INFO = 1, an singular value did not converge
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Based on contributions by
  //C     Ming Gu and Huan Ren, Computer Science Division, University of
  //C     California at Berkeley, USA
  //C
  //C  =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters.
  //C
  info = 0;
  //C
  if (n < 0) {
    info = -1;
  }
  else if ((sqre < 0) || (sqre > 1)) {
    info = -2;
  }
  //C
  int m = n + sqre;
  //C
  if (ldu < n) {
    info = -6;
  }
  else if (ldvt < m) {
    info = -8;
  }
  else if (smlsiz < 3) {
    info = -9;
  }
  if (info != 0) {
    xerbla("DLASD0", -info);
    return;
  }
  //C
  //C     If the input matrix is too small, call DLASDQ to find the SVD.
  //C
  if (n <= smlsiz) {
    dlasdq(cmn, "U", sqre, n, m, n, 0, d, e, vt, ldvt, u, ldu, u,
      ldu, work, info);
    return;
  }
  //C
  //C     Set up the computation tree.
  //C
  int inode = 1;
  int ndiml = inode + n;
  int ndimr = ndiml + n;
  int idxq = ndimr + n;
  int iwk = idxq + n;
  int nlvl = fem::int0;
  int nd = fem::int0;
  dlasdt(n, nlvl, nd, iwork(inode), iwork(ndiml), iwork(ndimr), smlsiz);
  //C
  //C     For the nodes on bottom level of the tree, solve
  //C     their subproblems by DLASDQ.
  //C
  int ndb1 = (nd + 1) / 2;
  int ncc = 0;
  int i = fem::int0;
  int i1 = fem::int0;
  int ic = fem::int0;
  int nl = fem::int0;
  int nlp1 = fem::int0;
  int nr = fem::int0;
  int nrp1 = fem::int0;
  int nlf = fem::int0;
  int nrf = fem::int0;
  int sqrei = fem::int0;
  int itemp = fem::int0;
  int j = fem::int0;
  FEM_DO(i, ndb1, nd) {
    //C
    //C     IC : center row of each node
    //C     NL : number of rows of left  subproblem
    //C     NR : number of rows of right subproblem
    //C     NLF: starting row of the left   subproblem
    //C     NRF: starting row of the right  subproblem
    //C
    i1 = i - 1;
    ic = iwork(inode + i1);
    nl = iwork(ndiml + i1);
    nlp1 = nl + 1;
    nr = iwork(ndimr + i1);
    nrp1 = nr + 1;
    nlf = ic - nl;
    nrf = ic + 1;
    sqrei = 1;
    dlasdq(cmn, "U", sqrei, nl, nlp1, nl, ncc, d(nlf), e(nlf), vt(nlf,
      nlf), ldvt, u(nlf, nlf), ldu, u(nlf, nlf), ldu, work, info);
    if (info != 0) {
      return;
    }
    itemp = idxq + nlf - 2;
    FEM_DO(j, 1, nl) {
      iwork(itemp + j) = j;
    }
    if (i == nd) {
      sqrei = sqre;
    }
    else {
      sqrei = 1;
    }
    nrp1 = nr + sqrei;
    dlasdq(cmn, "U", sqrei, nr, nrp1, nr, ncc, d(nrf), e(nrf), vt(nrf,
      nrf), ldvt, u(nrf, nrf), ldu, u(nrf, nrf), ldu, work, info);
    if (info != 0) {
      return;
    }
    itemp = idxq + ic;
    FEM_DO(j, 1, nr) {
      iwork(itemp + j - 1) = j;
    }
  }
  //C
  //C     Now conquer each subproblem bottom-up.
  //C
  int lvl = fem::int0;
  int lf = fem::int0;
  int ll = fem::int0;
  int im1 = fem::int0;
  int idxqc = fem::int0;
  double alpha = fem::double0;
  double beta = fem::double0;
  FEM_DOSTEP(lvl, nlvl, 1, -1) {
    //C
    //C        Find the first node LF and last node LL on the
    //C        current level LVL.
    //C
    if (lvl == 1) {
      lf = 1;
      ll = 1;
    }
    else {
      lf = fem::pow(2, (lvl - 1));
      ll = 2 * lf - 1;
    }
    FEM_DO(i, lf, ll) {
      im1 = i - 1;
      ic = iwork(inode + im1);
      nl = iwork(ndiml + im1);
      nr = iwork(ndimr + im1);
      nlf = ic - nl;
      if ((sqre == 0) && (i == ll)) {
        sqrei = sqre;
      }
      else {
        sqrei = 1;
      }
      idxqc = idxq + nlf - 1;
      alpha = d(ic);
      beta = e(ic);
      dlasd1(cmn, nl, nr, sqrei, d(nlf), alpha, beta, u(nlf, nlf),
        ldu, vt(nlf, nlf), ldvt, iwork(idxqc), iwork(iwk), work,
        info);
      if (info != 0) {
        return;
      }
    }
  }
  //C
  //C     End of DLASD0
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasd7.f
inline
void
dlasd7(
  common& cmn,
  int const& icompq,
  int const& nl,
  int const& nr,
  int const& sqre,
  int& k,
  arr_ref<double> d,
  arr_ref<double> z,
  arr_ref<double> zw,
  arr_ref<double> vf,
  arr_ref<double> vfw,
  arr_ref<double> vl,
  arr_ref<double> vlw,
  double const& alpha,
  double const& beta,
  arr_ref<double> dsigma,
  arr_ref<int> idx,
  arr_ref<int> idxp,
  arr_ref<int> idxq,
  arr_ref<int> perm,
  int& givptr,
  arr_ref<int, 2> givcol,
  int const& ldgcol,
  arr_ref<double, 2> givnum,
  int const& ldgnum,
  double& c,
  double& s,
  int& info)
{
  d(dimension(star));
  z(dimension(star));
  zw(dimension(star));
  vf(dimension(star));
  vfw(dimension(star));
  vl(dimension(star));
  vlw(dimension(star));
  dsigma(dimension(star));
  idx(dimension(star));
  idxp(dimension(star));
  idxq(dimension(star));
  perm(dimension(star));
  givcol(dimension(ldgcol, star));
  givnum(dimension(ldgnum, star));
  int n = fem::int0;
  int m = fem::int0;
  int nlp1 = fem::int0;
  int nlp2 = fem::int0;
  double z1 = fem::double0;
  const double zero = 0.0e+0;
  double tau = fem::double0;
  int i = fem::int0;
  int idxi = fem::int0;
  double eps = fem::double0;
  double tol = fem::double0;
  const double eight = 8.0e+0;
  int k2 = fem::int0;
  int j = fem::int0;
  int jprev = fem::int0;
  int idxjp = fem::int0;
  int idxj = fem::int0;
  int jp = fem::int0;
  const double two = 2.0e+0;
  double hlftol = fem::double0;
  const double one = 1.0e+0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASD7 merges the two sets of singular values together into a single
  //C  sorted set. Then it tries to deflate the size of the problem. There
  //C  are two ways in which deflation can occur:  when two or more singular
  //C  values are close together or if there is a tiny entry in the Z
  //C  vector. For each such occurrence the order of the related
  //C  secular equation problem is reduced by one.
  //C
  //C  DLASD7 is called from DLASD6.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  ICOMPQ  (input) INTEGER
  //C          Specifies whether singular vectors are to be computed
  //C          in compact form, as follows:
  //C          = 0: Compute singular values only.
  //C          = 1: Compute singular vectors of upper
  //C               bidiagonal matrix in compact form.
  //C
  //C  NL     (input) INTEGER
  //C         The row dimension of the upper block. NL >= 1.
  //C
  //C  NR     (input) INTEGER
  //C         The row dimension of the lower block. NR >= 1.
  //C
  //C  SQRE   (input) INTEGER
  //C         = 0: the lower block is an NR-by-NR square matrix.
  //C         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
  //C
  //C         The bidiagonal matrix has
  //C         N = NL + NR + 1 rows and
  //C         M = N + SQRE >= N columns.
  //C
  //C  K      (output) INTEGER
  //C         Contains the dimension of the non-deflated matrix, this is
  //C         the order of the related secular equation. 1 <= K <=N.
  //C
  //C  D      (input/output) DOUBLE PRECISION array, dimension ( N )
  //C         On entry D contains the singular values of the two submatrices
  //C         to be combined. On exit D contains the trailing (N-K) updated
  //C         singular values (those which were deflated) sorted into
  //C         increasing order.
  //C
  //C  Z      (output) DOUBLE PRECISION array, dimension ( M )
  //C         On exit Z contains the updating row vector in the secular
  //C         equation.
  //C
  //C  ZW     (workspace) DOUBLE PRECISION array, dimension ( M )
  //C         Workspace for Z.
  //C
  //C  VF     (input/output) DOUBLE PRECISION array, dimension ( M )
  //C         On entry, VF(1:NL+1) contains the first components of all
  //C         right singular vectors of the upper block; and VF(NL+2:M)
  //C         contains the first components of all right singular vectors
  //C         of the lower block. On exit, VF contains the first components
  //C         of all right singular vectors of the bidiagonal matrix.
  //C
  //C  VFW    (workspace) DOUBLE PRECISION array, dimension ( M )
  //C         Workspace for VF.
  //C
  //C  VL     (input/output) DOUBLE PRECISION array, dimension ( M )
  //C         On entry, VL(1:NL+1) contains the  last components of all
  //C         right singular vectors of the upper block; and VL(NL+2:M)
  //C         contains the last components of all right singular vectors
  //C         of the lower block. On exit, VL contains the last components
  //C         of all right singular vectors of the bidiagonal matrix.
  //C
  //C  VLW    (workspace) DOUBLE PRECISION array, dimension ( M )
  //C         Workspace for VL.
  //C
  //C  ALPHA  (input) DOUBLE PRECISION
  //C         Contains the diagonal element associated with the added row.
  //C
  //C  BETA   (input) DOUBLE PRECISION
  //C         Contains the off-diagonal element associated with the added
  //C         row.
  //C
  //C  DSIGMA (output) DOUBLE PRECISION array, dimension ( N )
  //C         Contains a copy of the diagonal elements (K-1 singular values
  //C         and one zero) in the secular equation.
  //C
  //C  IDX    (workspace) INTEGER array, dimension ( N )
  //C         This will contain the permutation used to sort the contents of
  //C         D into ascending order.
  //C
  //C  IDXP   (workspace) INTEGER array, dimension ( N )
  //C         This will contain the permutation used to place deflated
  //C         values of D at the end of the array. On output IDXP(2:K)
  //C         points to the nondeflated D-values and IDXP(K+1:N)
  //C         points to the deflated singular values.
  //C
  //C  IDXQ   (input) INTEGER array, dimension ( N )
  //C         This contains the permutation which separately sorts the two
  //C         sub-problems in D into ascending order.  Note that entries in
  //C         the first half of this permutation must first be moved one
  //C         position backward; and entries in the second half
  //C         must first have NL+1 added to their values.
  //C
  //C  PERM   (output) INTEGER array, dimension ( N )
  //C         The permutations (from deflation and sorting) to be applied
  //C         to each singular block. Not referenced if ICOMPQ = 0.
  //C
  //C  GIVPTR (output) INTEGER
  //C         The number of Givens rotations which took place in this
  //C         subproblem. Not referenced if ICOMPQ = 0.
  //C
  //C  GIVCOL (output) INTEGER array, dimension ( LDGCOL, 2 )
  //C         Each pair of numbers indicates a pair of columns to take place
  //C         in a Givens rotation. Not referenced if ICOMPQ = 0.
  //C
  //C  LDGCOL (input) INTEGER
  //C         The leading dimension of GIVCOL, must be at least N.
  //C
  //C  GIVNUM (output) DOUBLE PRECISION array, dimension ( LDGNUM, 2 )
  //C         Each number indicates the C or S value to be used in the
  //C         corresponding Givens rotation. Not referenced if ICOMPQ = 0.
  //C
  //C  LDGNUM (input) INTEGER
  //C         The leading dimension of GIVNUM, must be at least N.
  //C
  //C  C      (output) DOUBLE PRECISION
  //C         C contains garbage if SQRE =0 and the C-value of a Givens
  //C         rotation related to the right null space if SQRE = 1.
  //C
  //C  S      (output) DOUBLE PRECISION
  //C         S contains garbage if SQRE =0 and the S-value of a Givens
  //C         rotation related to the right null space if SQRE = 1.
  //C
  //C  INFO   (output) INTEGER
  //C         = 0:  successful exit.
  //C         < 0:  if INFO = -i, the i-th argument had an illegal value.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Based on contributions by
  //C     Ming Gu and Huan Ren, Computer Science Division, University of
  //C     California at Berkeley, USA
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters.
  //C
  info = 0;
  n = nl + nr + 1;
  m = n + sqre;
  //C
  if ((icompq < 0) || (icompq > 1)) {
    info = -1;
  }
  else if (nl < 1) {
    info = -2;
  }
  else if (nr < 1) {
    info = -3;
  }
  else if ((sqre < 0) || (sqre > 1)) {
    info = -4;
  }
  else if (ldgcol < n) {
    info = -22;
  }
  else if (ldgnum < n) {
    info = -24;
  }
  if (info != 0) {
    xerbla("DLASD7", -info);
    return;
  }
  //C
  nlp1 = nl + 1;
  nlp2 = nl + 2;
  if (icompq == 1) {
    givptr = 0;
  }
  //C
  //C     Generate the first part of the vector Z and move the singular
  //C     values in the first part of D one position backward.
  //C
  z1 = alpha * vl(nlp1);
  vl(nlp1) = zero;
  tau = vf(nlp1);
  FEM_DOSTEP(i, nl, 1, -1) {
    z(i + 1) = alpha * vl(i);
    vl(i) = zero;
    vf(i + 1) = vf(i);
    d(i + 1) = d(i);
    idxq(i + 1) = idxq(i) + 1;
  }
  vf(1) = tau;
  //C
  //C     Generate the second part of the vector Z.
  //C
  FEM_DO(i, nlp2, m) {
    z(i) = beta * vf(i);
    vf(i) = zero;
  }
  //C
  //C     Sort the singular values into increasing order
  //C
  FEM_DO(i, nlp2, n) {
    idxq(i) += nlp1;
  }
  //C
  //C     DSIGMA, IDXC, IDXC, and ZW are used as storage space.
  //C
  FEM_DO(i, 2, n) {
    dsigma(i) = d(idxq(i));
    zw(i) = z(idxq(i));
    vfw(i) = vf(idxq(i));
    vlw(i) = vl(idxq(i));
  }
  //C
  dlamrg(nl, nr, dsigma(2), 1, 1, idx(2));
  //C
  FEM_DO(i, 2, n) {
    idxi = 1 + idx(i);
    d(i) = dsigma(idxi);
    z(i) = zw(idxi);
    vf(i) = vfw(idxi);
    vl(i) = vlw(idxi);
  }
  //C
  //C     Calculate the allowable deflation tolerence
  //C
  eps = dlamch(cmn, "Epsilon");
  tol = fem::max(fem::abs(alpha), fem::abs(beta));
  tol = eight * eight * eps * fem::max(fem::abs(d(n)), tol);
  //C
  //C     There are 2 kinds of deflation -- first a value in the z-vector
  //C     is small, second two (or more) singular values are very close
  //C     together (their difference is small).
  //C
  //C     If the value in the z-vector is small, we simply permute the
  //C     array so that the corresponding singular value is moved to the
  //C     end.
  //C
  //C     If two values in the D-vector are close, we perform a two-sided
  //C     rotation designed to make one of the corresponding z-vector
  //C     entries zero, and then permute the array so that the deflated
  //C     singular value is moved to the end.
  //C
  //C     If there are multiple singular values then the problem deflates.
  //C     Here the number of equal singular values are found.  As each equal
  //C     singular value is found, an elementary reflector is computed to
  //C     rotate the corresponding singular subspace so that the
  //C     corresponding components of Z are zero in this new basis.
  //C
  k = 1;
  k2 = n + 1;
  FEM_DO(j, 2, n) {
    if (fem::abs(z(j)) <= tol) {
      //C
      //C           Deflate due to small z component.
      //C
      k2 = k2 - 1;
      idxp(k2) = j;
      if (j == n) {
        goto statement_100;
      }
    }
    else {
      jprev = j;
      goto statement_70;
    }
  }
  statement_70:
  j = jprev;
  statement_80:
  j++;
  if (j > n) {
    goto statement_90;
  }
  if (fem::abs(z(j)) <= tol) {
    //C
    //C        Deflate due to small z component.
    //C
    k2 = k2 - 1;
    idxp(k2) = j;
  }
  else {
    //C
    //C        Check if singular values are close enough to allow deflation.
    //C
    if (fem::abs(d(j) - d(jprev)) <= tol) {
      //C
      //C           Deflation is possible.
      //C
      s = z(jprev);
      c = z(j);
      //C
      //C           Find sqrt(a**2+b**2) without overflow or
      //C           destructive underflow.
      //C
      tau = dlapy2(c, s);
      z(j) = tau;
      z(jprev) = zero;
      c = c / tau;
      s = -s / tau;
      //C
      //C           Record the appropriate Givens rotation
      //C
      if (icompq == 1) {
        givptr++;
        idxjp = idxq(idx(jprev) + 1);
        idxj = idxq(idx(j) + 1);
        if (idxjp <= nlp1) {
          idxjp = idxjp - 1;
        }
        if (idxj <= nlp1) {
          idxj = idxj - 1;
        }
        givcol(givptr, 2) = idxjp;
        givcol(givptr, 1) = idxj;
        givnum(givptr, 2) = c;
        givnum(givptr, 1) = s;
      }
      drot(1, vf(jprev), 1, vf(j), 1, c, s);
      drot(1, vl(jprev), 1, vl(j), 1, c, s);
      k2 = k2 - 1;
      idxp(k2) = jprev;
      jprev = j;
    }
    else {
      k++;
      zw(k) = z(jprev);
      dsigma(k) = d(jprev);
      idxp(k) = jprev;
      jprev = j;
    }
  }
  goto statement_80;
  statement_90:
  //C
  //C     Record the last singular value.
  //C
  k++;
  zw(k) = z(jprev);
  dsigma(k) = d(jprev);
  idxp(k) = jprev;
  //C
  statement_100:
  //C
  //C     Sort the singular values into DSIGMA. The singular values which
  //C     were not deflated go into the first K slots of DSIGMA, except
  //C     that DSIGMA(1) is treated separately.
  //C
  FEM_DO(j, 2, n) {
    jp = idxp(j);
    dsigma(j) = d(jp);
    vfw(j) = vf(jp);
    vlw(j) = vl(jp);
  }
  if (icompq == 1) {
    FEM_DO(j, 2, n) {
      jp = idxp(j);
      perm(j) = idxq(idx(jp) + 1);
      if (perm(j) <= nlp1) {
        perm(j) = perm(j) - 1;
      }
    }
  }
  //C
  //C     The deflated singular values go back into the last N - K slots of
  //C     D.
  //C
  dcopy(n - k, dsigma(k + 1), 1, d(k + 1), 1);
  //C
  //C     Determine DSIGMA(1), DSIGMA(2), Z(1), VF(1), VL(1), VF(M), and
  //C     VL(M).
  //C
  dsigma(1) = zero;
  hlftol = tol / two;
  if (fem::abs(dsigma(2)) <= hlftol) {
    dsigma(2) = hlftol;
  }
  if (m > n) {
    z(1) = dlapy2(z1, z(m));
    if (z(1) <= tol) {
      c = one;
      s = zero;
      z(1) = tol;
    }
    else {
      c = z1 / z(1);
      s = -z(m) / z(1);
    }
    drot(1, vf(m), 1, vf(1), 1, c, s);
    drot(1, vl(m), 1, vl(1), 1, c, s);
  }
  else {
    if (fem::abs(z1) <= tol) {
      z(1) = tol;
    }
    else {
      z(1) = z1;
    }
  }
  //C
  //C     Restore Z, VF, and VL.
  //C
  dcopy(k - 1, zw(2), 1, z(2), 1);
  dcopy(n - 1, vfw(2), 1, vf(2), 1);
  dcopy(n - 1, vlw(2), 1, vl(2), 1);
  //C
  //C     End of DLASD7
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasd8.f
inline
void
dlasd8(
  common& cmn,
  int const& icompq,
  int const& k,
  arr_ref<double> d,
  arr_ref<double> z,
  arr_ref<double> vf,
  arr_ref<double> vl,
  arr_ref<double> difl,
  arr_ref<double, 2> difr,
  int const& lddifr,
  arr_ref<double> dsigma,
  arr_ref<double> work,
  int& info)
{
  d(dimension(star));
  z(dimension(star));
  vf(dimension(star));
  vl(dimension(star));
  difl(dimension(star));
  difr(dimension(lddifr, star));
  dsigma(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     October 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASD8 finds the square roots of the roots of the secular equation,
  //C  as defined by the values in DSIGMA and Z. It makes the appropriate
  //C  calls to DLASD4, and stores, for each  element in D, the distance
  //C  to its two nearest poles (elements in DSIGMA). It also updates
  //C  the arrays VF and VL, the first and last components of all the
  //C  right singular vectors of the original bidiagonal matrix.
  //C
  //C  DLASD8 is called from DLASD6.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  ICOMPQ  (input) INTEGER
  //C          Specifies whether singular vectors are to be computed in
  //C          factored form in the calling routine:
  //C          = 0: Compute singular values only.
  //C          = 1: Compute singular vectors in factored form as well.
  //C
  //C  K       (input) INTEGER
  //C          The number of terms in the rational function to be solved
  //C          by DLASD4.  K >= 1.
  //C
  //C  D       (output) DOUBLE PRECISION array, dimension ( K )
  //C          On output, D contains the updated singular values.
  //C
  //C  Z       (input/output) DOUBLE PRECISION array, dimension ( K )
  //C          On entry, the first K elements of this array contain the
  //C          components of the deflation-adjusted updating row vector.
  //C          On exit, Z is updated.
  //C
  //C  VF      (input/output) DOUBLE PRECISION array, dimension ( K )
  //C          On entry, VF contains  information passed through DBEDE8.
  //C          On exit, VF contains the first K components of the first
  //C          components of all right singular vectors of the bidiagonal
  //C          matrix.
  //C
  //C  VL      (input/output) DOUBLE PRECISION array, dimension ( K )
  //C          On entry, VL contains  information passed through DBEDE8.
  //C          On exit, VL contains the first K components of the last
  //C          components of all right singular vectors of the bidiagonal
  //C          matrix.
  //C
  //C  DIFL    (output) DOUBLE PRECISION array, dimension ( K )
  //C          On exit, DIFL(I) = D(I) - DSIGMA(I).
  //C
  //C  DIFR    (output) DOUBLE PRECISION array,
  //C                   dimension ( LDDIFR, 2 ) if ICOMPQ = 1 and
  //C                   dimension ( K ) if ICOMPQ = 0.
  //C          On exit, DIFR(I,1) = D(I) - DSIGMA(I+1), DIFR(K,1) is not
  //C          defined and will not be referenced.
  //C
  //C          If ICOMPQ = 1, DIFR(1:K,2) is an array containing the
  //C          normalizing factors for the right singular vector matrix.
  //C
  //C  LDDIFR  (input) INTEGER
  //C          The leading dimension of DIFR, must be at least K.
  //C
  //C  DSIGMA  (input/output) DOUBLE PRECISION array, dimension ( K )
  //C          On entry, the first K elements of this array contain the old
  //C          roots of the deflated updating problem.  These are the poles
  //C          of the secular equation.
  //C          On exit, the elements of DSIGMA may be very slightly altered
  //C          in value.
  //C
  //C  WORK    (workspace) DOUBLE PRECISION array, dimension at least 3 * K
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit.
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value.
  //C          > 0:  if INFO = 1, an singular value did not converge
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Based on contributions by
  //C     Ming Gu and Huan Ren, Computer Science Division, University of
  //C     California at Berkeley, USA
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters.
  //C
  info = 0;
  //C
  if ((icompq < 0) || (icompq > 1)) {
    info = -1;
  }
  else if (k < 1) {
    info = -2;
  }
  else if (lddifr < k) {
    info = -9;
  }
  if (info != 0) {
    xerbla("DLASD8", -info);
    return;
  }
  //C
  //C     Quick return if possible
  //C
  const double one = 1.0e+0;
  if (k == 1) {
    d(1) = fem::abs(z(1));
    difl(1) = d(1);
    if (icompq == 1) {
      difl(2) = one;
      difr(1, 2) = one;
    }
    return;
  }
  //C
  //C     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
  //C     be computed with high relative accuracy (barring over/underflow).
  //C     This is a problem on machines without a guard digit in
  //C     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
  //C     The following code replaces DSIGMA(I) by 2*DSIGMA(I)-DSIGMA(I),
  //C     which on any of these machines zeros out the bottommost
  //C     bit of DSIGMA(I) if it is 1; this makes the subsequent
  //C     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
  //C     occurs. On binary machines with a guard digit (almost all
  //C     machines) it does not change DSIGMA(I) at all. On hexadecimal
  //C     and decimal machines with a guard digit, it slightly
  //C     changes the bottommost bits of DSIGMA(I). It does not account
  //C     for hexadecimal or decimal machines without guard digits
  //C     (we know of none). We use a subroutine call to compute
  //C     2*DLAMBDA(I) to prevent optimizing compilers from eliminating
  //C     this code.
  //C
  int i = fem::int0;
  FEM_DO(i, 1, k) {
    dsigma(i) = dlamc3(dsigma(i), dsigma(i)) - dsigma(i);
  }
  //C
  //C     Book keeping.
  //C
  int iwk1 = 1;
  int iwk2 = iwk1 + k;
  int iwk3 = iwk2 + k;
  int iwk2i = iwk2 - 1;
  int iwk3i = iwk3 - 1;
  //C
  //C     Normalize Z.
  //C
  double rho = dnrm2(k, z, 1);
  dlascl(cmn, "G", 0, 0, rho, one, k, 1, z, k, info);
  rho = rho * rho;
  //C
  //C     Initialize WORK(IWK3).
  //C
  dlaset("A", k, 1, one, one, work(iwk3), k);
  //C
  //C     Compute the updated singular values, the arrays DIFL, DIFR,
  //C     and the updated Z.
  //C
  int j = fem::int0;
  FEM_DO(j, 1, k) {
    dlasd4(cmn, k, j, dsigma, z, work(iwk1), rho, d(j), work(iwk2), info);
    //C
    //C        If the root finder fails, the computation is terminated.
    //C
    if (info != 0) {
      return;
    }
    work(iwk3i + j) = work(iwk3i + j) * work(j) * work(iwk2i + j);
    difl(j) = -work(j);
    difr(j, 1) = -work(j + 1);
    {
      int fem_do_last = j - 1;
      FEM_DO(i, 1, fem_do_last) {
        work(iwk3i + i) = work(iwk3i + i) * work(i) * work(iwk2i +
          i) / (dsigma(i) - dsigma(j)) / (dsigma(i) + dsigma(j));
      }
    }
    FEM_DO(i, j + 1, k) {
      work(iwk3i + i) = work(iwk3i + i) * work(i) * work(iwk2i + i) /
        (dsigma(i) - dsigma(j)) / (dsigma(i) + dsigma(j));
    }
  }
  //C
  //C     Compute updated Z.
  //C
  FEM_DO(i, 1, k) {
    z(i) = fem::sign(fem::sqrt(fem::abs(work(iwk3i + i))), z(i));
  }
  //C
  //C     Update VF and VL.
  //C
  double diflj = fem::double0;
  double dj = fem::double0;
  double dsigj = fem::double0;
  double difrj = fem::double0;
  double dsigjp = fem::double0;
  double temp = fem::double0;
  FEM_DO(j, 1, k) {
    diflj = difl(j);
    dj = d(j);
    dsigj = -dsigma(j);
    if (j < k) {
      difrj = -difr(j, 1);
      dsigjp = -dsigma(j + 1);
    }
    work(j) = -z(j) / diflj / (dsigma(j) + dj);
    {
      int fem_do_last = j - 1;
      FEM_DO(i, 1, fem_do_last) {
        work(i) = z(i) / (dlamc3(dsigma(i), dsigj) - diflj) / (dsigma(i) + dj);
      }
    }
    FEM_DO(i, j + 1, k) {
      work(i) = z(i) / (dlamc3(dsigma(i), dsigjp) + difrj) / (dsigma(i) + dj);
    }
    temp = dnrm2(k, work, 1);
    work(iwk2i + j) = ddot(k, work, 1, vf, 1) / temp;
    work(iwk3i + j) = ddot(k, work, 1, vl, 1) / temp;
    if (icompq == 1) {
      difr(j, 2) = temp;
    }
  }
  //C
  dcopy(k, work(iwk2), 1, vf, 1);
  dcopy(k, work(iwk3), 1, vl, 1);
  //C
  //C     End of DLASD8
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasd6.f
inline
void
dlasd6(
  common& cmn,
  int const& icompq,
  int const& nl,
  int const& nr,
  int const& sqre,
  arr_ref<double> d,
  arr_ref<double> vf,
  arr_ref<double> vl,
  double& alpha,
  double& beta,
  arr_ref<int> idxq,
  arr_ref<int> perm,
  int& givptr,
  arr_ref<int, 2> givcol,
  int const& ldgcol,
  arr_ref<double, 2> givnum,
  int const& ldgnum,
  arr_ref<double, 2> poles,
  arr_ref<double> difl,
  arr_ref<double> difr,
  arr_ref<double> z,
  int& k,
  double& c,
  double& s,
  arr_ref<double> work,
  arr_ref<int> iwork,
  int& info)
{
  d(dimension(star));
  vf(dimension(star));
  vl(dimension(star));
  idxq(dimension(star));
  perm(dimension(star));
  givcol(dimension(ldgcol, star));
  givnum(dimension(ldgnum, star));
  poles(dimension(ldgnum, star));
  difl(dimension(star));
  difr(dimension(star));
  z(dimension(star));
  work(dimension(star));
  iwork(dimension(star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLASD6 computes the SVD of an updated upper bidiagonal matrix B
  //C  obtained by merging two smaller ones by appending a row. This
  //C  routine is used only for the problem which requires all singular
  //C  values and optionally singular vector matrices in factored form.
  //C  B is an N-by-M matrix with N = NL + NR + 1 and M = N + SQRE.
  //C  A related subroutine, DLASD1, handles the case in which all singular
  //C  values and singular vectors of the bidiagonal matrix are desired.
  //C
  //C  DLASD6 computes the SVD as follows:
  //C
  //C                ( D1(in)  0    0     0 )
  //C    B = U(in) * (   Z1'   a   Z2'    b ) * VT(in)
  //C                (   0     0   D2(in) 0 )
  //C
  //C      = U(out) * ( D(out) 0) * VT(out)
  //C
  //C  where Z' = (Z1' a Z2' b) = u' VT', and u is a vector of dimension M
  //C  with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros
  //C  elsewhere; and the entry b is empty if SQRE = 0.
  //C
  //C  The singular values of B can be computed using D1, D2, the first
  //C  components of all the right singular vectors of the lower block, and
  //C  the last components of all the right singular vectors of the upper
  //C  block. These components are stored and updated in VF and VL,
  //C  respectively, in DLASD6. Hence U and VT are not explicitly
  //C  referenced.
  //C
  //C  The singular values are stored in D. The algorithm consists of two
  //C  stages:
  //C
  //C        The first stage consists of deflating the size of the problem
  //C        when there are multiple singular values or if there is a zero
  //C        in the Z vector. For each such occurence the dimension of the
  //C        secular equation problem is reduced by one. This stage is
  //C        performed by the routine DLASD7.
  //C
  //C        The second stage consists of calculating the updated
  //C        singular values. This is done by finding the roots of the
  //C        secular equation via the routine DLASD4 (as called by DLASD8).
  //C        This routine also updates VF and VL and computes the distances
  //C        between the updated singular values and the old singular
  //C        values.
  //C
  //C  DLASD6 is called from DLASDA.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  ICOMPQ (input) INTEGER
  //C         Specifies whether singular vectors are to be computed in
  //C         factored form:
  //C         = 0: Compute singular values only.
  //C         = 1: Compute singular vectors in factored form as well.
  //C
  //C  NL     (input) INTEGER
  //C         The row dimension of the upper block.  NL >= 1.
  //C
  //C  NR     (input) INTEGER
  //C         The row dimension of the lower block.  NR >= 1.
  //C
  //C  SQRE   (input) INTEGER
  //C         = 0: the lower block is an NR-by-NR square matrix.
  //C         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
  //C
  //C         The bidiagonal matrix has row dimension N = NL + NR + 1,
  //C         and column dimension M = N + SQRE.
  //C
  //C  D      (input/output) DOUBLE PRECISION array, dimension ( NL+NR+1 ).
  //C         On entry D(1:NL,1:NL) contains the singular values of the
  //C         upper block, and D(NL+2:N) contains the singular values
  //C         of the lower block. On exit D(1:N) contains the singular
  //C         values of the modified matrix.
  //C
  //C  VF     (input/output) DOUBLE PRECISION array, dimension ( M )
  //C         On entry, VF(1:NL+1) contains the first components of all
  //C         right singular vectors of the upper block; and VF(NL+2:M)
  //C         contains the first components of all right singular vectors
  //C         of the lower block. On exit, VF contains the first components
  //C         of all right singular vectors of the bidiagonal matrix.
  //C
  //C  VL     (input/output) DOUBLE PRECISION array, dimension ( M )
  //C         On entry, VL(1:NL+1) contains the  last components of all
  //C         right singular vectors of the upper block; and VL(NL+2:M)
  //C         contains the last components of all right singular vectors of
  //C         the lower block. On exit, VL contains the last components of
  //C         all right singular vectors of the bidiagonal matrix.
  //C
  //C  ALPHA  (input/output) DOUBLE PRECISION
  //C         Contains the diagonal element associated with the added row.
  //C
  //C  BETA   (input/output) DOUBLE PRECISION
  //C         Contains the off-diagonal element associated with the added
  //C         row.
  //C
  //C  IDXQ   (output) INTEGER array, dimension ( N )
  //C         This contains the permutation which will reintegrate the
  //C         subproblem just solved back into sorted order, i.e.
  //C         D( IDXQ( I = 1, N ) ) will be in ascending order.
  //C
  //C  PERM   (output) INTEGER array, dimension ( N )
  //C         The permutations (from deflation and sorting) to be applied
  //C         to each block. Not referenced if ICOMPQ = 0.
  //C
  //C  GIVPTR (output) INTEGER
  //C         The number of Givens rotations which took place in this
  //C         subproblem. Not referenced if ICOMPQ = 0.
  //C
  //C  GIVCOL (output) INTEGER array, dimension ( LDGCOL, 2 )
  //C         Each pair of numbers indicates a pair of columns to take place
  //C         in a Givens rotation. Not referenced if ICOMPQ = 0.
  //C
  //C  LDGCOL (input) INTEGER
  //C         leading dimension of GIVCOL, must be at least N.
  //C
  //C  GIVNUM (output) DOUBLE PRECISION array, dimension ( LDGNUM, 2 )
  //C         Each number indicates the C or S value to be used in the
  //C         corresponding Givens rotation. Not referenced if ICOMPQ = 0.
  //C
  //C  LDGNUM (input) INTEGER
  //C         The leading dimension of GIVNUM and POLES, must be at least N.
  //C
  //C  POLES  (output) DOUBLE PRECISION array, dimension ( LDGNUM, 2 )
  //C         On exit, POLES(1,*) is an array containing the new singular
  //C         values obtained from solving the secular equation, and
  //C         POLES(2,*) is an array containing the poles in the secular
  //C         equation. Not referenced if ICOMPQ = 0.
  //C
  //C  DIFL   (output) DOUBLE PRECISION array, dimension ( N )
  //C         On exit, DIFL(I) is the distance between I-th updated
  //C         (undeflated) singular value and the I-th (undeflated) old
  //C         singular value.
  //C
  //C  DIFR   (output) DOUBLE PRECISION array,
  //C                  dimension ( LDGNUM, 2 ) if ICOMPQ = 1 and
  //C                  dimension ( N ) if ICOMPQ = 0.
  //C         On exit, DIFR(I, 1) is the distance between I-th updated
  //C         (undeflated) singular value and the I+1-th (undeflated) old
  //C         singular value.
  //C
  //C         If ICOMPQ = 1, DIFR(1:K,2) is an array containing the
  //C         normalizing factors for the right singular vector matrix.
  //C
  //C         See DLASD8 for details on DIFL and DIFR.
  //C
  //C  Z      (output) DOUBLE PRECISION array, dimension ( M )
  //C         The first elements of this array contain the components
  //C         of the deflation-adjusted updating row vector.
  //C
  //C  K      (output) INTEGER
  //C         Contains the dimension of the non-deflated matrix,
  //C         This is the order of the related secular equation. 1 <= K <=N.
  //C
  //C  C      (output) DOUBLE PRECISION
  //C         C contains garbage if SQRE =0 and the C-value of a Givens
  //C         rotation related to the right null space if SQRE = 1.
  //C
  //C  S      (output) DOUBLE PRECISION
  //C         S contains garbage if SQRE =0 and the S-value of a Givens
  //C         rotation related to the right null space if SQRE = 1.
  //C
  //C  WORK   (workspace) DOUBLE PRECISION array, dimension ( 4 * M )
  //C
  //C  IWORK  (workspace) INTEGER array, dimension ( 3 * N )
  //C
  //C  INFO   (output) INTEGER
  //C          = 0:  successful exit.
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value.
  //C          > 0:  if INFO = 1, an singular value did not converge
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Based on contributions by
  //C     Ming Gu and Huan Ren, Computer Science Division, University of
  //C     California at Berkeley, USA
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters.
  //C
  info = 0;
  int n = nl + nr + 1;
  int m = n + sqre;
  //C
  if ((icompq < 0) || (icompq > 1)) {
    info = -1;
  }
  else if (nl < 1) {
    info = -2;
  }
  else if (nr < 1) {
    info = -3;
  }
  else if ((sqre < 0) || (sqre > 1)) {
    info = -4;
  }
  else if (ldgcol < n) {
    info = -14;
  }
  else if (ldgnum < n) {
    info = -16;
  }
  if (info != 0) {
    xerbla("DLASD6", -info);
    return;
  }
  //C
  //C     The following values are for bookkeeping purposes only.  They are
  //C     integer pointers which indicate the portion of the workspace
  //C     used by a particular array in DLASD7 and DLASD8.
  //C
  int isigma = 1;
  int iw = isigma + n;
  int ivfw = iw + m;
  int ivlw = ivfw + m;
  //C
  int idx = 1;
  int idxc = idx + n;
  int idxp = idxc + n;
  //C
  //C     Scale.
  //C
  double orgnrm = fem::max(fem::abs(alpha), fem::abs(beta));
  const double zero = 0.0e+0;
  d(nl + 1) = zero;
  int i = fem::int0;
  FEM_DO(i, 1, n) {
    if (fem::abs(d(i)) > orgnrm) {
      orgnrm = fem::abs(d(i));
    }
  }
  const double one = 1.0e+0;
  dlascl(cmn, "G", 0, 0, orgnrm, one, n, 1, d, n, info);
  alpha = alpha / orgnrm;
  beta = beta / orgnrm;
  //C
  //C     Sort and Deflate singular values.
  //C
  dlasd7(cmn, icompq, nl, nr, sqre, k, d, z, work(iw), vf, work(ivfw),
    vl, work(ivlw), alpha, beta, work(isigma), iwork(idx), iwork(idxp),
    idxq, perm, givptr, givcol, ldgcol, givnum, ldgnum, c, s, info);
  //C
  //C     Solve Secular Equation, compute DIFL, DIFR, and update VF, VL.
  //C
  dlasd8(cmn, icompq, k, d, z, vf, vl, difl, difr, ldgnum, work(isigma),
    work(iw), info);
  //C
  //C     Save the poles if ICOMPQ = 1.
  //C
  if (icompq == 1) {
    dcopy(k, d, 1, poles(1, 1), 1);
    dcopy(k, work(isigma), 1, poles(1, 2), 1);
  }
  //C
  //C     Unscale.
  //C
  dlascl(cmn, "G", 0, 0, one, orgnrm, n, 1, d, n, info);
  //C
  //C     Prepare the IDXQ sorting permutation.
  //C
  int n1 = k;
  int n2 = n - k;
  dlamrg(n1, n2, d, 1, -1, idxq);
  //C
  //C     End of DLASD6
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlasda.f
inline
void
dlasda(
  common& cmn,
  int const& icompq,
  int const& smlsiz,
  int const& n,
  int const& sqre,
  arr_ref<double> d,
  arr_ref<double> e,
  arr_ref<double, 2> u,
  int const& ldu,
  arr_ref<double, 2> vt,
  arr_ref<int> k,
  arr_ref<double, 2> difl,
  arr_ref<double, 2> difr,
  arr_ref<double, 2> z,
  arr_ref<double, 2> poles,
  arr_ref<int> givptr,
  arr_ref<int, 2> givcol,
  int const& ldgcol,
  arr_ref<int, 2> perm,
  arr_ref<double, 2> givnum,
  arr_ref<double> c,
  arr_ref<double> s,
  arr_ref<double> work,
  arr_ref<int> iwork,
  int& info)
{
  d(dimension(star));
  e(dimension(star));
  u(dimension(ldu, star));
  vt(dimension(ldu, star));
  k(dimension(star));
  difl(dimension(ldu, star));
  difr(dimension(ldu, star));
  z(dimension(ldu, star));
  poles(dimension(ldu, star));
  givptr(dimension(star));
  givcol(dimension(ldgcol, star));
  perm(dimension(ldgcol, star));
  givnum(dimension(ldu, star));
  c(dimension(star));
  s(dimension(star));
  work(dimension(star));
  iwork(dimension(star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  Using a divide and conquer approach, DLASDA computes the singular
  //C  value decomposition (SVD) of a real upper bidiagonal N-by-M matrix
  //C  B with diagonal D and offdiagonal E, where M = N + SQRE. The
  //C  algorithm computes the singular values in the SVD B = U * S * VT.
  //C  The orthogonal matrices U and VT are optionally computed in
  //C  compact form.
  //C
  //C  A related subroutine, DLASD0, computes the singular values and
  //C  the singular vectors in explicit form.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  ICOMPQ (input) INTEGER
  //C         Specifies whether singular vectors are to be computed
  //C         in compact form, as follows
  //C         = 0: Compute singular values only.
  //C         = 1: Compute singular vectors of upper bidiagonal
  //C              matrix in compact form.
  //C
  //C  SMLSIZ (input) INTEGER
  //C         The maximum size of the subproblems at the bottom of the
  //C         computation tree.
  //C
  //C  N      (input) INTEGER
  //C         The row dimension of the upper bidiagonal matrix. This is
  //C         also the dimension of the main diagonal array D.
  //C
  //C  SQRE   (input) INTEGER
  //C         Specifies the column dimension of the bidiagonal matrix.
  //C         = 0: The bidiagonal matrix has column dimension M = N;
  //C         = 1: The bidiagonal matrix has column dimension M = N + 1.
  //C
  //C  D      (input/output) DOUBLE PRECISION array, dimension ( N )
  //C         On entry D contains the main diagonal of the bidiagonal
  //C         matrix. On exit D, if INFO = 0, contains its singular values.
  //C
  //C  E      (input) DOUBLE PRECISION array, dimension ( M-1 )
  //C         Contains the subdiagonal entries of the bidiagonal matrix.
  //C         On exit, E has been destroyed.
  //C
  //C  U      (output) DOUBLE PRECISION array,
  //C         dimension ( LDU, SMLSIZ ) if ICOMPQ = 1, and not referenced
  //C         if ICOMPQ = 0. If ICOMPQ = 1, on exit, U contains the left
  //C         singular vector matrices of all subproblems at the bottom
  //C         level.
  //C
  //C  LDU    (input) INTEGER, LDU = > N.
  //C         The leading dimension of arrays U, VT, DIFL, DIFR, POLES,
  //C         GIVNUM, and Z.
  //C
  //C  VT     (output) DOUBLE PRECISION array,
  //C         dimension ( LDU, SMLSIZ+1 ) if ICOMPQ = 1, and not referenced
  //C         if ICOMPQ = 0. If ICOMPQ = 1, on exit, VT' contains the right
  //C         singular vector matrices of all subproblems at the bottom
  //C         level.
  //C
  //C  K      (output) INTEGER array,
  //C         dimension ( N ) if ICOMPQ = 1 and dimension 1 if ICOMPQ = 0.
  //C         If ICOMPQ = 1, on exit, K(I) is the dimension of the I-th
  //C         secular equation on the computation tree.
  //C
  //C  DIFL   (output) DOUBLE PRECISION array, dimension ( LDU, NLVL ),
  //C         where NLVL = floor(log_2 (N/SMLSIZ))).
  //C
  //C  DIFR   (output) DOUBLE PRECISION array,
  //C                  dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1 and
  //C                  dimension ( N ) if ICOMPQ = 0.
  //C         If ICOMPQ = 1, on exit, DIFL(1:N, I) and DIFR(1:N, 2 * I - 1)
  //C         record distances between singular values on the I-th
  //C         level and singular values on the (I -1)-th level, and
  //C         DIFR(1:N, 2 * I ) contains the normalizing factors for
  //C         the right singular vector matrix. See DLASD8 for details.
  //C
  //C  Z      (output) DOUBLE PRECISION array,
  //C                  dimension ( LDU, NLVL ) if ICOMPQ = 1 and
  //C                  dimension ( N ) if ICOMPQ = 0.
  //C         The first K elements of Z(1, I) contain the components of
  //C         the deflation-adjusted updating row vector for subproblems
  //C         on the I-th level.
  //C
  //C  POLES  (output) DOUBLE PRECISION array,
  //C         dimension ( LDU, 2 * NLVL ) if ICOMPQ = 1, and not referenced
  //C         if ICOMPQ = 0. If ICOMPQ = 1, on exit, POLES(1, 2*I - 1) and
  //C         POLES(1, 2*I) contain  the new and old singular values
  //C         involved in the secular equations on the I-th level.
  //C
  //C  GIVPTR (output) INTEGER array,
  //C         dimension ( N ) if ICOMPQ = 1, and not referenced if
  //C         ICOMPQ = 0. If ICOMPQ = 1, on exit, GIVPTR( I ) records
  //C         the number of Givens rotations performed on the I-th
  //C         problem on the computation tree.
  //C
  //C  GIVCOL (output) INTEGER array,
  //C         dimension ( LDGCOL, 2 * NLVL ) if ICOMPQ = 1, and not
  //C         referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I,
  //C         GIVCOL(1, 2 *I - 1) and GIVCOL(1, 2 *I) record the locations
  //C         of Givens rotations performed on the I-th level on the
  //C         computation tree.
  //C
  //C  LDGCOL (input) INTEGER, LDGCOL = > N.
  //C         The leading dimension of arrays GIVCOL and PERM.
  //C
  //C  PERM   (output) INTEGER array,
  //C         dimension ( LDGCOL, NLVL ) if ICOMPQ = 1, and not referenced
  //C         if ICOMPQ = 0. If ICOMPQ = 1, on exit, PERM(1, I) records
  //C         permutations done on the I-th level of the computation tree.
  //C
  //C  GIVNUM (output) DOUBLE PRECISION array,
  //C         dimension ( LDU,  2 * NLVL ) if ICOMPQ = 1, and not
  //C         referenced if ICOMPQ = 0. If ICOMPQ = 1, on exit, for each I,
  //C         GIVNUM(1, 2 *I - 1) and GIVNUM(1, 2 *I) record the C- and S-
  //C         values of Givens rotations performed on the I-th level on
  //C         the computation tree.
  //C
  //C  C      (output) DOUBLE PRECISION array,
  //C         dimension ( N ) if ICOMPQ = 1, and dimension 1 if ICOMPQ = 0.
  //C         If ICOMPQ = 1 and the I-th subproblem is not square, on exit,
  //C         C( I ) contains the C-value of a Givens rotation related to
  //C         the right null space of the I-th subproblem.
  //C
  //C  S      (output) DOUBLE PRECISION array, dimension ( N ) if
  //C         ICOMPQ = 1, and dimension 1 if ICOMPQ = 0. If ICOMPQ = 1
  //C         and the I-th subproblem is not square, on exit, S( I )
  //C         contains the S-value of a Givens rotation related to
  //C         the right null space of the I-th subproblem.
  //C
  //C  WORK   (workspace) DOUBLE PRECISION array, dimension
  //C         (6 * N + (SMLSIZ + 1)*(SMLSIZ + 1)).
  //C
  //C  IWORK  (workspace) INTEGER array.
  //C         Dimension must be at least (7 * N).
  //C
  //C  INFO   (output) INTEGER
  //C          = 0:  successful exit.
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value.
  //C          > 0:  if INFO = 1, an singular value did not converge
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Based on contributions by
  //C     Ming Gu and Huan Ren, Computer Science Division, University of
  //C     California at Berkeley, USA
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters.
  //C
  info = 0;
  //C
  if ((icompq < 0) || (icompq > 1)) {
    info = -1;
  }
  else if (smlsiz < 3) {
    info = -2;
  }
  else if (n < 0) {
    info = -3;
  }
  else if ((sqre < 0) || (sqre > 1)) {
    info = -4;
  }
  else if (ldu < (n + sqre)) {
    info = -8;
  }
  else if (ldgcol < n) {
    info = -17;
  }
  if (info != 0) {
    xerbla("DLASDA", -info);
    return;
  }
  //C
  int m = n + sqre;
  //C
  //C     If the input matrix is too small, call DLASDQ to find the SVD.
  //C
  if (n <= smlsiz) {
    if (icompq == 0) {
      dlasdq(cmn, "U", sqre, n, 0, 0, 0, d, e, vt, ldu, u, ldu, u,
        ldu, work, info);
    }
    else {
      dlasdq(cmn, "U", sqre, n, m, n, 0, d, e, vt, ldu, u, ldu, u,
        ldu, work, info);
    }
    return;
  }
  //C
  //C     Book-keeping and  set up the computation tree.
  //C
  int inode = 1;
  int ndiml = inode + n;
  int ndimr = ndiml + n;
  int idxq = ndimr + n;
  int iwk = idxq + n;
  //C
  int ncc = 0;
  int nru = 0;
  //C
  int smlszp = smlsiz + 1;
  int vf = 1;
  int vl = vf + m;
  int nwork1 = vl + m;
  int nwork2 = nwork1 + smlszp * smlszp;
  //C
  int nlvl = fem::int0;
  int nd = fem::int0;
  dlasdt(n, nlvl, nd, iwork(inode), iwork(ndiml), iwork(ndimr), smlsiz);
  //C
  //C     for the nodes on bottom level of the tree, solve
  //C     their subproblems by DLASDQ.
  //C
  int ndb1 = (nd + 1) / 2;
  int i = fem::int0;
  int i1 = fem::int0;
  int ic = fem::int0;
  int nl = fem::int0;
  int nlp1 = fem::int0;
  int nr = fem::int0;
  int nlf = fem::int0;
  int nrf = fem::int0;
  int idxqi = fem::int0;
  int vfi = fem::int0;
  int vli = fem::int0;
  int sqrei = fem::int0;
  const double zero = 0.0e+0;
  const double one = 1.0e+0;
  int itemp = fem::int0;
  int j = fem::int0;
  int nrp1 = fem::int0;
  FEM_DO(i, ndb1, nd) {
    //C
    //C        IC : center row of each node
    //C        NL : number of rows of left  subproblem
    //C        NR : number of rows of right subproblem
    //C        NLF: starting row of the left   subproblem
    //C        NRF: starting row of the right  subproblem
    //C
    i1 = i - 1;
    ic = iwork(inode + i1);
    nl = iwork(ndiml + i1);
    nlp1 = nl + 1;
    nr = iwork(ndimr + i1);
    nlf = ic - nl;
    nrf = ic + 1;
    idxqi = idxq + nlf - 2;
    vfi = vf + nlf - 1;
    vli = vl + nlf - 1;
    sqrei = 1;
    if (icompq == 0) {
      dlaset("A", nlp1, nlp1, zero, one, work(nwork1), smlszp);
      dlasdq(cmn, "U", sqrei, nl, nlp1, nru, ncc, d(nlf), e(nlf),
        work(nwork1), smlszp, work(nwork2), nl, work(nwork2), nl,
        work(nwork2), info);
      itemp = nwork1 + nl * smlszp;
      dcopy(nlp1, work(nwork1), 1, work(vfi), 1);
      dcopy(nlp1, work(itemp), 1, work(vli), 1);
    }
    else {
      dlaset("A", nl, nl, zero, one, u(nlf, 1), ldu);
      dlaset("A", nlp1, nlp1, zero, one, vt(nlf, 1), ldu);
      dlasdq(cmn, "U", sqrei, nl, nlp1, nl, ncc, d(nlf), e(nlf), vt(nlf,
        1), ldu, u(nlf, 1), ldu, u(nlf, 1), ldu, work(nwork1), info);
      dcopy(nlp1, vt(nlf, 1), 1, work(vfi), 1);
      dcopy(nlp1, vt(nlf, nlp1), 1, work(vli), 1);
    }
    if (info != 0) {
      return;
    }
    FEM_DO(j, 1, nl) {
      iwork(idxqi + j) = j;
    }
    if ((i == nd) && (sqre == 0)) {
      sqrei = 0;
    }
    else {
      sqrei = 1;
    }
    idxqi += nlp1;
    vfi += nlp1;
    vli += nlp1;
    nrp1 = nr + sqrei;
    if (icompq == 0) {
      dlaset("A", nrp1, nrp1, zero, one, work(nwork1), smlszp);
      dlasdq(cmn, "U", sqrei, nr, nrp1, nru, ncc, d(nrf), e(nrf),
        work(nwork1), smlszp, work(nwork2), nr, work(nwork2), nr,
        work(nwork2), info);
      itemp = nwork1 + (nrp1 - 1) * smlszp;
      dcopy(nrp1, work(nwork1), 1, work(vfi), 1);
      dcopy(nrp1, work(itemp), 1, work(vli), 1);
    }
    else {
      dlaset("A", nr, nr, zero, one, u(nrf, 1), ldu);
      dlaset("A", nrp1, nrp1, zero, one, vt(nrf, 1), ldu);
      dlasdq(cmn, "U", sqrei, nr, nrp1, nr, ncc, d(nrf), e(nrf), vt(nrf,
        1), ldu, u(nrf, 1), ldu, u(nrf, 1), ldu, work(nwork1), info);
      dcopy(nrp1, vt(nrf, 1), 1, work(vfi), 1);
      dcopy(nrp1, vt(nrf, nrp1), 1, work(vli), 1);
    }
    if (info != 0) {
      return;
    }
    FEM_DO(j, 1, nr) {
      iwork(idxqi + j) = j;
    }
  }
  //C
  //C     Now conquer each subproblem bottom-up.
  //C
  j = fem::pow(2, nlvl);
  int lvl = fem::int0;
  int lvl2 = fem::int0;
  int lf = fem::int0;
  int ll = fem::int0;
  int im1 = fem::int0;
  double alpha = fem::double0;
  double beta = fem::double0;
  FEM_DOSTEP(lvl, nlvl, 1, -1) {
    lvl2 = lvl * 2 - 1;
    //C
    //C        Find the first node LF and last node LL on
    //C        the current level LVL.
    //C
    if (lvl == 1) {
      lf = 1;
      ll = 1;
    }
    else {
      lf = fem::pow(2, (lvl - 1));
      ll = 2 * lf - 1;
    }
    FEM_DO(i, lf, ll) {
      im1 = i - 1;
      ic = iwork(inode + im1);
      nl = iwork(ndiml + im1);
      nr = iwork(ndimr + im1);
      nlf = ic - nl;
      nrf = ic + 1;
      if (i == ll) {
        sqrei = sqre;
      }
      else {
        sqrei = 1;
      }
      vfi = vf + nlf - 1;
      vli = vl + nlf - 1;
      idxqi = idxq + nlf - 1;
      alpha = d(ic);
      beta = e(ic);
      if (icompq == 0) {
        dlasd6(cmn, icompq, nl, nr, sqrei, d(nlf), work(vfi), work(vli),
          alpha, beta, iwork(idxqi), perm, givptr(1), givcol, ldgcol,
          givnum, ldu, poles, difl, difr, z, k(1), c(1), s(1), work(nwork1),
          iwork(iwk), info);
      }
      else {
        j = j - 1;
        dlasd6(cmn, icompq, nl, nr, sqrei, d(nlf), work(vfi), work(vli),
          alpha, beta, iwork(idxqi), perm(nlf, lvl), givptr(j),
          givcol(nlf, lvl2), ldgcol, givnum(nlf, lvl2), ldu, poles(nlf,
          lvl2), difl(nlf, lvl), difr(nlf, lvl2), z(nlf, lvl), k(j),
          c(j), s(j), work(nwork1), iwork(iwk), info);
      }
      if (info != 0) {
        return;
      }
    }
  }
  //C
  //C     End of DLASDA
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dbdsdc.f
inline
void
dbdsdc(
  common& cmn,
  str_cref uplo,
  str_cref compq,
  int const& n,
  arr_ref<double> d,
  arr_ref<double> e,
  arr_ref<double, 2> u,
  int const& ldu,
  arr_ref<double, 2> vt,
  int const& ldvt,
  arr_ref<double> q,
  arr_ref<int> iq,
  arr_ref<double> work,
  arr_ref<int> iwork,
  int& info)
{
  d(dimension(star));
  e(dimension(star));
  u(dimension(ldu, star));
  vt(dimension(ldvt, star));
  q(dimension(star));
  iq(dimension(star));
  work(dimension(star));
  iwork(dimension(star));
  int iuplo = fem::int0;
  int icompq = fem::int0;
  int smlsiz = fem::int0;
  const double one = 1.0e+0;
  int nm1 = fem::int0;
  int wstart = fem::int0;
  int qstart = fem::int0;
  int i = fem::int0;
  double cs = fem::double0;
  double sn = fem::double0;
  double r = fem::double0;
  const double zero = 0.0e+0;
  int iu = fem::int0;
  int ivt = fem::int0;
  double orgnrm = fem::double0;
  int ierr = fem::int0;
  double eps = fem::double0;
  const double two = 2.0e+0;
  int mlvl = fem::int0;
  int smlszp = fem::int0;
  int difl = fem::int0;
  int difr = fem::int0;
  int z = fem::int0;
  int ic = fem::int0;
  int is = fem::int0;
  int poles = fem::int0;
  int givnum = fem::int0;
  int k = fem::int0;
  int givptr = fem::int0;
  int perm = fem::int0;
  int givcol = fem::int0;
  int start = fem::int0;
  int sqre = fem::int0;
  int nsize = fem::int0;
  int ii = fem::int0;
  int kk = fem::int0;
  double p = fem::double0;
  int j = fem::int0;
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DBDSDC computes the singular value decomposition (SVD) of a real
  //C  N-by-N (upper or lower) bidiagonal matrix B:  B = U * S * VT,
  //C  using a divide and conquer method, where S is a diagonal matrix
  //C  with non-negative diagonal elements (the singular values of B), and
  //C  U and VT are orthogonal matrices of left and right singular vectors,
  //C  respectively. DBDSDC can be used to compute all singular values,
  //C  and optionally, singular vectors or singular vectors in compact form.
  //C
  //C  This code makes very mild assumptions about floating point
  //C  arithmetic. It will work on machines with a guard digit in
  //C  add/subtract, or on those binary machines without guard digits
  //C  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
  //C  It could conceivably fail on hexadecimal or decimal machines
  //C  without guard digits, but we know of none.  See DLASD3 for details.
  //C
  //C  The code currently calls DLASDQ if singular values only are desired.
  //C  However, it can be slightly modified to compute singular values
  //C  using the divide and conquer method.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  UPLO    (input) CHARACTER*1
  //C          = 'U':  B is upper bidiagonal.
  //C          = 'L':  B is lower bidiagonal.
  //C
  //C  COMPQ   (input) CHARACTER*1
  //C          Specifies whether singular vectors are to be computed
  //C          as follows:
  //C          = 'N':  Compute singular values only;
  //C          = 'P':  Compute singular values and compute singular
  //C                  vectors in compact form;
  //C          = 'I':  Compute singular values and singular vectors.
  //C
  //C  N       (input) INTEGER
  //C          The order of the matrix B.  N >= 0.
  //C
  //C  D       (input/output) DOUBLE PRECISION array, dimension (N)
  //C          On entry, the n diagonal elements of the bidiagonal matrix B.
  //C          On exit, if INFO=0, the singular values of B.
  //C
  //C  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
  //C          On entry, the elements of E contain the offdiagonal
  //C          elements of the bidiagonal matrix whose SVD is desired.
  //C          On exit, E has been destroyed.
  //C
  //C  U       (output) DOUBLE PRECISION array, dimension (LDU,N)
  //C          If  COMPQ = 'I', then:
  //C             On exit, if INFO = 0, U contains the left singular vectors
  //C             of the bidiagonal matrix.
  //C          For other values of COMPQ, U is not referenced.
  //C
  //C  LDU     (input) INTEGER
  //C          The leading dimension of the array U.  LDU >= 1.
  //C          If singular vectors are desired, then LDU >= max( 1, N ).
  //C
  //C  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
  //C          If  COMPQ = 'I', then:
  //C             On exit, if INFO = 0, VT' contains the right singular
  //C             vectors of the bidiagonal matrix.
  //C          For other values of COMPQ, VT is not referenced.
  //C
  //C  LDVT    (input) INTEGER
  //C          The leading dimension of the array VT.  LDVT >= 1.
  //C          If singular vectors are desired, then LDVT >= max( 1, N ).
  //C
  //C  Q       (output) DOUBLE PRECISION array, dimension (LDQ)
  //C          If  COMPQ = 'P', then:
  //C             On exit, if INFO = 0, Q and IQ contain the left
  //C             and right singular vectors in a compact form,
  //C             requiring O(N log N) space instead of 2*N**2.
  //C             In particular, Q contains all the DOUBLE PRECISION data in
  //C             LDQ >= N*(11 + 2*SMLSIZ + 8*INT(LOG_2(N/(SMLSIZ+1))))
  //C             words of memory, where SMLSIZ is returned by ILAENV and
  //C             is equal to the maximum size of the subproblems at the
  //C             bottom of the computation tree (usually about 25).
  //C          For other values of COMPQ, Q is not referenced.
  //C
  //C  IQ      (output) INTEGER array, dimension (LDIQ)
  //C          If  COMPQ = 'P', then:
  //C             On exit, if INFO = 0, Q and IQ contain the left
  //C             and right singular vectors in a compact form,
  //C             requiring O(N log N) space instead of 2*N**2.
  //C             In particular, IQ contains all INTEGER data in
  //C             LDIQ >= N*(3 + 3*INT(LOG_2(N/(SMLSIZ+1))))
  //C             words of memory, where SMLSIZ is returned by ILAENV and
  //C             is equal to the maximum size of the subproblems at the
  //C             bottom of the computation tree (usually about 25).
  //C          For other values of COMPQ, IQ is not referenced.
  //C
  //C  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //C          If COMPQ = 'N' then LWORK >= (4 * N).
  //C          If COMPQ = 'P' then LWORK >= (6 * N).
  //C          If COMPQ = 'I' then LWORK >= (3 * N**2 + 4 * N).
  //C
  //C  IWORK   (workspace) INTEGER array, dimension (8*N)
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit.
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value.
  //C          > 0:  The algorithm failed to compute an singular value.
  //C                The update process of divide and conquer failed.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Based on contributions by
  //C     Ming Gu and Huan Ren, Computer Science Division, University of
  //C     California at Berkeley, USA
  //C
  //C  =====================================================================
  //C  Changed dimension statement in comment describing E from (N) to
  //C  (N-1).  Sven, 17 Feb 05.
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters.
  //C
  info = 0;
  //C
  iuplo = 0;
  if (lsame(uplo, "U")) {
    iuplo = 1;
  }
  if (lsame(uplo, "L")) {
    iuplo = 2;
  }
  if (lsame(compq, "N")) {
    icompq = 0;
  }
  else if (lsame(compq, "P")) {
    icompq = 1;
  }
  else if (lsame(compq, "I")) {
    icompq = 2;
  }
  else {
    icompq = -1;
  }
  if (iuplo == 0) {
    info = -1;
  }
  else if (icompq < 0) {
    info = -2;
  }
  else if (n < 0) {
    info = -3;
  }
  else if ((ldu < 1) || ((icompq == 2) && (ldu < n))) {
    info = -7;
  }
  else if ((ldvt < 1) || ((icompq == 2) && (ldvt < n))) {
    info = -9;
  }
  if (info != 0) {
    xerbla("DBDSDC", -info);
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (n == 0) {
    return;
  }
  smlsiz = ilaenv(9, "DBDSDC", " ", 0, 0, 0, 0);
  if (n == 1) {
    if (icompq == 1) {
      q(1) = fem::sign(one, d(1));
      q(1 + smlsiz * n) = one;
    }
    else if (icompq == 2) {
      u(1, 1) = fem::sign(one, d(1));
      vt(1, 1) = one;
    }
    d(1) = fem::abs(d(1));
    return;
  }
  nm1 = n - 1;
  //C
  //C     If matrix lower bidiagonal, rotate to be upper bidiagonal
  //C     by applying Givens rotations on the left
  //C
  wstart = 1;
  qstart = 3;
  if (icompq == 1) {
    dcopy(n, d, 1, q(1), 1);
    dcopy(n - 1, e, 1, q(n + 1), 1);
  }
  if (iuplo == 2) {
    qstart = 5;
    wstart = 2 * n - 1;
    {
      int fem_do_last = n - 1;
      FEM_DO(i, 1, fem_do_last) {
        dlartg(cmn, d(i), e(i), cs, sn, r);
        d(i) = r;
        e(i) = sn * d(i + 1);
        d(i + 1) = cs * d(i + 1);
        if (icompq == 1) {
          q(i + 2 * n) = cs;
          q(i + 3 * n) = sn;
        }
        else if (icompq == 2) {
          work(i) = cs;
          work(nm1 + i) = -sn;
        }
      }
    }
  }
  //C
  //C     If ICOMPQ = 0, use DLASDQ to compute the singular values.
  //C
  if (icompq == 0) {
    dlasdq(cmn, "U", 0, n, 0, 0, 0, d, e, vt, ldvt, u, ldu, u, ldu,
      work(wstart), info);
    goto statement_40;
  }
  //C
  //C     If N is smaller than the minimum divide size SMLSIZ, then solve
  //C     the problem with another solver.
  //C
  if (n <= smlsiz) {
    if (icompq == 2) {
      dlaset("A", n, n, zero, one, u, ldu);
      dlaset("A", n, n, zero, one, vt, ldvt);
      dlasdq(cmn, "U", 0, n, n, n, 0, d, e, vt, ldvt, u, ldu, u, ldu,
        work(wstart), info);
    }
    else if (icompq == 1) {
      iu = 1;
      ivt = iu + n;
      dlaset("A", n, n, zero, one, q(iu + (qstart - 1) * n), n);
      dlaset("A", n, n, zero, one, q(ivt + (qstart - 1) * n), n);
      dlasdq(cmn, "U", 0, n, n, n, 0, d, e, q(ivt + (qstart - 1) * n),
        n, q(iu + (qstart - 1) * n), n, q(iu + (qstart - 1) * n), n,
        work(wstart), info);
    }
    goto statement_40;
  }
  //C
  if (icompq == 2) {
    dlaset("A", n, n, zero, one, u, ldu);
    dlaset("A", n, n, zero, one, vt, ldvt);
  }
  //C
  //C     Scale.
  //C
  orgnrm = dlanst("M", n, d, e);
  if (orgnrm == zero) {
    return;
  }
  dlascl(cmn, "G", 0, 0, orgnrm, one, n, 1, d, n, ierr);
  dlascl(cmn, "G", 0, 0, orgnrm, one, nm1, 1, e, nm1, ierr);
  //C
  eps = dlamch(cmn, "Epsilon");
  //C
  mlvl = fem::fint(fem::log(fem::dble(n) / fem::dble(smlsiz + 1)) /
    fem::log(two)) + 1;
  smlszp = smlsiz + 1;
  //C
  if (icompq == 1) {
    iu = 1;
    ivt = 1 + smlsiz;
    difl = ivt + smlszp;
    difr = difl + mlvl;
    z = difr + mlvl * 2;
    ic = z + mlvl;
    is = ic + 1;
    poles = is + 1;
    givnum = poles + 2 * mlvl;
    //C
    k = 1;
    givptr = 2;
    perm = 3;
    givcol = perm + mlvl;
  }
  //C
  FEM_DO(i, 1, n) {
    if (fem::abs(d(i)) < eps) {
      d(i) = fem::sign(eps, d(i));
    }
  }
  //C
  start = 1;
  sqre = 0;
  //C
  FEM_DO(i, 1, nm1) {
    if ((fem::abs(e(i)) < eps) || (i == nm1)) {
      //C
      //C        Subproblem found. First determine its size and then
      //C        apply divide and conquer on it.
      //C
      if (i < nm1) {
        //C
        //C        A subproblem with E(I) small for I < NM1.
        //C
        nsize = i - start + 1;
      }
      else if (fem::abs(e(i)) >= eps) {
        //C
        //C        A subproblem with E(NM1) not too small but I = NM1.
        //C
        nsize = n - start + 1;
      }
      else {
        //C
        //C        A subproblem with E(NM1) small. This implies an
        //C        1-by-1 subproblem at D(N). Solve this 1-by-1 problem
        //C        first.
        //C
        nsize = i - start + 1;
        if (icompq == 2) {
          u(n, n) = fem::sign(one, d(n));
          vt(n, n) = one;
        }
        else if (icompq == 1) {
          q(n + (qstart - 1) * n) = fem::sign(one, d(n));
          q(n + (smlsiz + qstart - 1) * n) = one;
        }
        d(n) = fem::abs(d(n));
      }
      if (icompq == 2) {
        dlasd0(cmn, nsize, sqre, d(start), e(start), u(start, start),
          ldu, vt(start, start), ldvt, smlsiz, iwork, work(wstart),
          info);
      }
      else {
        dlasda(cmn, icompq, smlsiz, nsize, sqre, d(start), e(start),
          q(start + (iu + qstart - 2) * n), n, q(start + (ivt +
          qstart - 2) * n), iq(start + k * n), q(start + (difl +
          qstart - 2) * n), q(start + (difr + qstart - 2) * n), q(
          start + (z + qstart - 2) * n), q(start + (poles + qstart - 2) * n),
          iq(start + givptr * n), iq(start + givcol * n), n, iq(
          start + perm * n), q(start + (givnum + qstart - 2) * n), q(
          start + (ic + qstart - 2) * n), q(start + (is + qstart - 2) * n),
          work(wstart), iwork, info);
        if (info != 0) {
          return;
        }
      }
      start = i + 1;
    }
  }
  //C
  //C     Unscale
  //C
  dlascl(cmn, "G", 0, 0, one, orgnrm, n, 1, d, n, ierr);
  statement_40:
  //C
  //C     Use Selection Sort to minimize swaps of singular vectors
  //C
  FEM_DO(ii, 2, n) {
    i = ii - 1;
    kk = i;
    p = d(i);
    FEM_DO(j, ii, n) {
      if (d(j) > p) {
        kk = j;
        p = d(j);
      }
    }
    if (kk != i) {
      d(kk) = d(i);
      d(i) = p;
      if (icompq == 1) {
        iq(i) = kk;
      }
      else if (icompq == 2) {
        dswap(n, u(1, i), 1, u(1, kk), 1);
        dswap(n, vt(i, 1), ldvt, vt(kk, 1), ldvt);
      }
    }
    else if (icompq == 1) {
      iq(i) = i;
    }
  }
  //C
  //C     If ICOMPQ = 1, use IQ(N,1) as the indicator for UPLO
  //C
  if (icompq == 1) {
    if (iuplo == 1) {
      iq(n) = 1;
    }
    else {
      iq(n) = 0;
    }
  }
  //C
  //C     If B is lower bidiagonal, update U by those Givens rotations
  //C     which rotated B to be upper bidiagonal
  //C
  if ((iuplo == 2) && (icompq == 2)) {
    dlasr("L", "V", "B", n, n, work(1), work(n), u, ldu);
  }
  //C
  //C     End of DBDSDC
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlarf.f
inline
void
dlarf(
  str_cref side,
  int const& m,
  int const& n,
  arr_cref<double> v,
  int const& incv,
  double const& tau,
  arr_ref<double, 2> c,
  int const& ldc,
  arr_ref<double> work)
{
  v(dimension(star));
  c(dimension(ldc, star));
  work(dimension(star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLARF applies a real elementary reflector H to a real m by n matrix
  //C  C, from either the left or the right. H is represented in the form
  //C
  //C        H = I - tau * v * v'
  //C
  //C  where tau is a real scalar and v is a real vector.
  //C
  //C  If tau = 0, then H is taken to be the unit matrix.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  SIDE    (input) CHARACTER*1
  //C          = 'L': form  H * C
  //C          = 'R': form  C * H
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix C.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix C.
  //C
  //C  V       (input) DOUBLE PRECISION array, dimension
  //C                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
  //C                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
  //C          The vector v in the representation of H. V is not used if
  //C          TAU = 0.
  //C
  //C  INCV    (input) INTEGER
  //C          The increment between elements of v. INCV <> 0.
  //C
  //C  TAU     (input) DOUBLE PRECISION
  //C          The value tau in the representation of H.
  //C
  //C  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
  //C          On entry, the m by n matrix C.
  //C          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
  //C          or C * H if SIDE = 'R'.
  //C
  //C  LDC     (input) INTEGER
  //C          The leading dimension of the array C. LDC >= max(1,M).
  //C
  //C  WORK    (workspace) DOUBLE PRECISION array, dimension
  //C                         (N) if SIDE = 'L'
  //C                      or (M) if SIDE = 'R'
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  bool applyleft = lsame(side, "L");
  int lastv = 0;
  int lastc = 0;
  const double zero = 0.0e+0;
  int i = fem::int0;
  if (tau != zero) {
    //C     Set up variables for scanning V.  LASTV begins pointing to the end
    //C     of V.
    if (applyleft) {
      lastv = m;
    }
    else {
      lastv = n;
    }
    if (incv > 0) {
      i = 1 + (lastv - 1) * incv;
    }
    else {
      i = 1;
    }
    //C     Look for the last non-zero row in V.
    while (lastv > 0 && v(i) == zero) {
      lastv = lastv - 1;
      i = i - incv;
    }
    if (applyleft) {
      //C     Scan for the last non-zero column in C(1:lastv,:).
      lastc = iladlc(lastv, n, c, ldc);
    }
    else {
      //C     Scan for the last non-zero row in C(:,1:lastv).
      lastc = iladlr(m, lastv, c, ldc);
    }
  }
  //C     Note that lastc.eq.0 renders the BLAS operations null; no special
  //C     case is needed at this level.
  const double one = 1.0e+0;
  if (applyleft) {
    //C
    //C        Form  H * C
    //C
    if (lastv > 0) {
      //C
      //C           w(1:lastc,1) := C(1:lastv,1:lastc)' * v(1:lastv,1)
      //C
      dgemv("Transpose", lastv, lastc, one, c, ldc, v, incv, zero, work, 1);
      //C
      //C           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)'
      //C
      dger(lastv, lastc, -tau, v, incv, work, 1, c, ldc);
    }
  }
  else {
    //C
    //C        Form  C * H
    //C
    if (lastv > 0) {
      //C
      //C           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
      //C
      dgemv("No transpose", lastc, lastv, one, c, ldc, v, incv, zero, work, 1);
      //C
      //C           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)'
      //C
      dger(lastc, lastv, -tau, work, 1, v, incv, c, ldc);
    }
  }
  //C
  //C     End of DLARF
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlarfg.f
inline
void
dlarfg(
  common& cmn,
  int const& n,
  double& alpha,
  arr_ref<double> x,
  int const& incx,
  double& tau)
{
  x(dimension(star));
  const double zero = 0.0e+0;
  double xnorm = fem::double0;
  double beta = fem::double0;
  double safmin = fem::double0;
  int knt = fem::int0;
  const double one = 1.0e+0;
  double rsafmn = fem::double0;
  int j = fem::int0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLARFG generates a real elementary reflector H of order n, such
  //C  that
  //C
  //C        H * ( alpha ) = ( beta ),   H' * H = I.
  //C            (   x   )   (   0  )
  //C
  //C  where alpha and beta are scalars, and x is an (n-1)-element real
  //C  vector. H is represented in the form
  //C
  //C        H = I - tau * ( 1 ) * ( 1 v' ) ,
  //C                      ( v )
  //C
  //C  where tau is a real scalar and v is a real (n-1)-element
  //C  vector.
  //C
  //C  If the elements of x are all zero, then tau = 0 and H is taken to be
  //C  the unit matrix.
  //C
  //C  Otherwise  1 <= tau <= 2.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  N       (input) INTEGER
  //C          The order of the elementary reflector.
  //C
  //C  ALPHA   (input/output) DOUBLE PRECISION
  //C          On entry, the value alpha.
  //C          On exit, it is overwritten with the value beta.
  //C
  //C  X       (input/output) DOUBLE PRECISION array, dimension
  //C                         (1+(N-2)*abs(INCX))
  //C          On entry, the vector x.
  //C          On exit, it is overwritten with the vector v.
  //C
  //C  INCX    (input) INTEGER
  //C          The increment between elements of X. INCX > 0.
  //C
  //C  TAU     (output) DOUBLE PRECISION
  //C          The value tau.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  if (n <= 1) {
    tau = zero;
    return;
  }
  //C
  xnorm = dnrm2(n - 1, x, incx);
  //C
  if (xnorm == zero) {
    //C
    //C        H  =  I
    //C
    tau = zero;
  }
  else {
    //C
    //C        general case
    //C
    beta = -fem::sign(dlapy2(alpha, xnorm), alpha);
    safmin = dlamch(cmn, "S") / dlamch(cmn, "E");
    knt = 0;
    if (fem::abs(beta) < safmin) {
      //C
      //C           XNORM, BETA may be inaccurate; scale X and recompute them
      //C
      rsafmn = one / safmin;
      statement_10:
      knt++;
      dscal(n - 1, rsafmn, x, incx);
      beta = beta * rsafmn;
      alpha = alpha * rsafmn;
      if (fem::abs(beta) < safmin) {
        goto statement_10;
      }
      //C
      //C           New BETA is at most 1, at least SAFMIN
      //C
      xnorm = dnrm2(n - 1, x, incx);
      beta = -fem::sign(dlapy2(alpha, xnorm), alpha);
    }
    tau = (beta - alpha) / beta;
    dscal(n - 1, one / (alpha - beta), x, incx);
    //C
    //C        If ALPHA is subnormal, it may lose relative accuracy
    //C
    FEM_DO(j, 1, knt) {
      beta = beta * safmin;
    }
    alpha = beta;
  }
  //C
  //C     End of DLARFG
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dgebd2.f
inline
void
dgebd2(
  common& cmn,
  int const& m,
  int const& n,
  arr_ref<double, 2> a,
  int const& lda,
  arr_ref<double> d,
  arr_ref<double> e,
  arr_ref<double> tauq,
  arr_ref<double> taup,
  arr_ref<double> work,
  int& info)
{
  a(dimension(lda, star));
  d(dimension(star));
  e(dimension(star));
  tauq(dimension(star));
  taup(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DGEBD2 reduces a real general m by n matrix A to upper or lower
  //C  bidiagonal form B by an orthogonal transformation: Q' * A * P = B.
  //C
  //C  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  M       (input) INTEGER
  //C          The number of rows in the matrix A.  M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns in the matrix A.  N >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the m by n general matrix to be reduced.
  //C          On exit,
  //C          if m >= n, the diagonal and the first superdiagonal are
  //C            overwritten with the upper bidiagonal matrix B; the
  //C            elements below the diagonal, with the array TAUQ, represent
  //C            the orthogonal matrix Q as a product of elementary
  //C            reflectors, and the elements above the first superdiagonal,
  //C            with the array TAUP, represent the orthogonal matrix P as
  //C            a product of elementary reflectors;
  //C          if m < n, the diagonal and the first subdiagonal are
  //C            overwritten with the lower bidiagonal matrix B; the
  //C            elements below the first subdiagonal, with the array TAUQ,
  //C            represent the orthogonal matrix Q as a product of
  //C            elementary reflectors, and the elements above the diagonal,
  //C            with the array TAUP, represent the orthogonal matrix P as
  //C            a product of elementary reflectors.
  //C          See Further Details.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(1,M).
  //C
  //C  D       (output) DOUBLE PRECISION array, dimension (min(M,N))
  //C          The diagonal elements of the bidiagonal matrix B:
  //C          D(i) = A(i,i).
  //C
  //C  E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)
  //C          The off-diagonal elements of the bidiagonal matrix B:
  //C          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
  //C          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
  //C
  //C  TAUQ    (output) DOUBLE PRECISION array dimension (min(M,N))
  //C          The scalar factors of the elementary reflectors which
  //C          represent the orthogonal matrix Q. See Further Details.
  //C
  //C  TAUP    (output) DOUBLE PRECISION array, dimension (min(M,N))
  //C          The scalar factors of the elementary reflectors which
  //C          represent the orthogonal matrix P. See Further Details.
  //C
  //C  WORK    (workspace) DOUBLE PRECISION array, dimension (max(M,N))
  //C
  //C  INFO    (output) INTEGER
  //C          = 0: successful exit.
  //C          < 0: if INFO = -i, the i-th argument had an illegal value.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  The matrices Q and P are represented as products of elementary
  //C  reflectors:
  //C
  //C  If m >= n,
  //C
  //C     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
  //C
  //C  Each H(i) and G(i) has the form:
  //C
  //C     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
  //C
  //C  where tauq and taup are real scalars, and v and u are real vectors;
  //C  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
  //C  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
  //C  tauq is stored in TAUQ(i) and taup in TAUP(i).
  //C
  //C  If m < n,
  //C
  //C     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
  //C
  //C  Each H(i) and G(i) has the form:
  //C
  //C     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
  //C
  //C  where tauq and taup are real scalars, and v and u are real vectors;
  //C  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
  //C  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
  //C  tauq is stored in TAUQ(i) and taup in TAUP(i).
  //C
  //C  The contents of A on exit are illustrated by the following examples:
  //C
  //C  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
  //C
  //C    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
  //C    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
  //C    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
  //C    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
  //C    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
  //C    (  v1  v2  v3  v4  v5 )
  //C
  //C  where d and e denote diagonal and off-diagonal elements of B, vi
  //C  denotes an element of the vector defining H(i), and ui an element of
  //C  the vector defining G(i).
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters
  //C
  info = 0;
  if (m < 0) {
    info = -1;
  }
  else if (n < 0) {
    info = -2;
  }
  else if (lda < fem::max(1, m)) {
    info = -4;
  }
  if (info < 0) {
    xerbla("DGEBD2", -info);
    return;
  }
  //C
  int i = fem::int0;
  const double one = 1.0e+0;
  const double zero = 0.0e+0;
  if (m >= n) {
    //C
    //C        Reduce to upper bidiagonal form
    //C
    FEM_DO(i, 1, n) {
      //C
      //C           Generate elementary reflector H(i) to annihilate A(i+1:m,i)
      //C
      dlarfg(cmn, m - i + 1, a(i, i), a(fem::min(i + 1, m), i), 1, tauq(i));
      d(i) = a(i, i);
      a(i, i) = one;
      //C
      //C           Apply H(i) to A(i:m,i+1:n) from the left
      //C
      if (i < n) {
        dlarf("Left", m - i + 1, n - i, a(i, i), 1, tauq(i), a(i, i + 1),
          lda, work);
      }
      a(i, i) = d(i);
      //C
      if (i < n) {
        //C
        //C              Generate elementary reflector G(i) to annihilate
        //C              A(i,i+2:n)
        //C
        dlarfg(cmn, n - i, a(i, i + 1), a(i, fem::min(i + 2, n)), lda, taup(i));
        e(i) = a(i, i + 1);
        a(i, i + 1) = one;
        //C
        //C              Apply G(i) to A(i+1:m,i+1:n) from the right
        //C
        dlarf("Right", m - i, n - i, a(i, i + 1), lda, taup(i), a(i + 1,
          i + 1), lda, work);
        a(i, i + 1) = e(i);
      }
      else {
        taup(i) = zero;
      }
    }
  }
  else {
    //C
    //C        Reduce to lower bidiagonal form
    //C
    FEM_DO(i, 1, m) {
      //C
      //C           Generate elementary reflector G(i) to annihilate A(i,i+1:n)
      //C
      dlarfg(cmn, n - i + 1, a(i, i), a(i, fem::min(i + 1, n)), lda, taup(i));
      d(i) = a(i, i);
      a(i, i) = one;
      //C
      //C           Apply G(i) to A(i+1:m,i:n) from the right
      //C
      if (i < m) {
        dlarf("Right", m - i, n - i + 1, a(i, i), lda, taup(i), a(i + 1,
          i), lda, work);
      }
      a(i, i) = d(i);
      //C
      if (i < m) {
        //C
        //C              Generate elementary reflector H(i) to annihilate
        //C              A(i+2:m,i)
        //C
        dlarfg(cmn, m - i, a(i + 1, i), a(fem::min(i + 2, m), i), 1, tauq(i));
        e(i) = a(i + 1, i);
        a(i + 1, i) = one;
        //C
        //C              Apply H(i) to A(i+1:m,i+1:n) from the left
        //C
        dlarf("Left", m - i, n - i, a(i + 1, i), 1, tauq(i), a(i + 1,
          i + 1), lda, work);
        a(i + 1, i) = e(i);
      }
      else {
        tauq(i) = zero;
      }
    }
  }
  //C
  //C     End of DGEBD2
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlabrd.f
inline
void
dlabrd(
  common& cmn,
  int const& m,
  int const& n,
  int const& nb,
  arr_ref<double, 2> a,
  int const& lda,
  arr_ref<double> d,
  arr_ref<double> e,
  arr_ref<double> tauq,
  arr_ref<double> taup,
  arr_ref<double, 2> x,
  int const& ldx,
  arr_ref<double, 2> y,
  int const& ldy)
{
  a(dimension(lda, star));
  d(dimension(star));
  e(dimension(star));
  tauq(dimension(star));
  taup(dimension(star));
  x(dimension(ldx, star));
  y(dimension(ldy, star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLABRD reduces the first NB rows and columns of a real general
  //C  m by n matrix A to upper or lower bidiagonal form by an orthogonal
  //C  transformation Q' * A * P, and returns the matrices X and Y which
  //C  are needed to apply the transformation to the unreduced part of A.
  //C
  //C  If m >= n, A is reduced to upper bidiagonal form; if m < n, to lower
  //C  bidiagonal form.
  //C
  //C  This is an auxiliary routine called by DGEBRD
  //C
  //C  Arguments
  //C  =========
  //C
  //C  M       (input) INTEGER
  //C          The number of rows in the matrix A.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns in the matrix A.
  //C
  //C  NB      (input) INTEGER
  //C          The number of leading rows and columns of A to be reduced.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the m by n general matrix to be reduced.
  //C          On exit, the first NB rows and columns of the matrix are
  //C          overwritten; the rest of the array is unchanged.
  //C          If m >= n, elements on and below the diagonal in the first NB
  //C            columns, with the array TAUQ, represent the orthogonal
  //C            matrix Q as a product of elementary reflectors; and
  //C            elements above the diagonal in the first NB rows, with the
  //C            array TAUP, represent the orthogonal matrix P as a product
  //C            of elementary reflectors.
  //C          If m < n, elements below the diagonal in the first NB
  //C            columns, with the array TAUQ, represent the orthogonal
  //C            matrix Q as a product of elementary reflectors, and
  //C            elements on and above the diagonal in the first NB rows,
  //C            with the array TAUP, represent the orthogonal matrix P as
  //C            a product of elementary reflectors.
  //C          See Further Details.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(1,M).
  //C
  //C  D       (output) DOUBLE PRECISION array, dimension (NB)
  //C          The diagonal elements of the first NB rows and columns of
  //C          the reduced matrix.  D(i) = A(i,i).
  //C
  //C  E       (output) DOUBLE PRECISION array, dimension (NB)
  //C          The off-diagonal elements of the first NB rows and columns of
  //C          the reduced matrix.
  //C
  //C  TAUQ    (output) DOUBLE PRECISION array dimension (NB)
  //C          The scalar factors of the elementary reflectors which
  //C          represent the orthogonal matrix Q. See Further Details.
  //C
  //C  TAUP    (output) DOUBLE PRECISION array, dimension (NB)
  //C          The scalar factors of the elementary reflectors which
  //C          represent the orthogonal matrix P. See Further Details.
  //C
  //C  X       (output) DOUBLE PRECISION array, dimension (LDX,NB)
  //C          The m-by-nb matrix X required to update the unreduced part
  //C          of A.
  //C
  //C  LDX     (input) INTEGER
  //C          The leading dimension of the array X. LDX >= M.
  //C
  //C  Y       (output) DOUBLE PRECISION array, dimension (LDY,NB)
  //C          The n-by-nb matrix Y required to update the unreduced part
  //C          of A.
  //C
  //C  LDY     (input) INTEGER
  //C          The leading dimension of the array Y. LDY >= N.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  The matrices Q and P are represented as products of elementary
  //C  reflectors:
  //C
  //C     Q = H(1) H(2) . . . H(nb)  and  P = G(1) G(2) . . . G(nb)
  //C
  //C  Each H(i) and G(i) has the form:
  //C
  //C     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
  //C
  //C  where tauq and taup are real scalars, and v and u are real vectors.
  //C
  //C  If m >= n, v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored on exit in
  //C  A(i:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+1:n) is stored on exit in
  //C  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
  //C
  //C  If m < n, v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored on exit in
  //C  A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i:n) is stored on exit in
  //C  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
  //C
  //C  The elements of the vectors v and u together form the m-by-nb matrix
  //C  V and the nb-by-n matrix U' which are needed, with X and Y, to apply
  //C  the transformation to the unreduced part of the matrix, using a block
  //C  update of the form:  A := A - V*Y' - X*U'.
  //C
  //C  The contents of A on exit are illustrated by the following examples
  //C  with nb = 2:
  //C
  //C  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
  //C
  //C    (  1   1   u1  u1  u1 )           (  1   u1  u1  u1  u1  u1 )
  //C    (  v1  1   1   u2  u2 )           (  1   1   u2  u2  u2  u2 )
  //C    (  v1  v2  a   a   a  )           (  v1  1   a   a   a   a  )
  //C    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
  //C    (  v1  v2  a   a   a  )
  //C
  //C  where a denotes an element of the original matrix which is unchanged,
  //C  vi denotes an element of the vector defining H(i), and ui an element
  //C  of the vector defining G(i).
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Quick return if possible
  //C
  if (m <= 0 || n <= 0) {
    return;
  }
  //C
  int i = fem::int0;
  const double one = 1.0e0;
  const double zero = 0.0e0;
  if (m >= n) {
    //C
    //C        Reduce to upper bidiagonal form
    //C
    FEM_DO(i, 1, nb) {
      //C
      //C           Update A(i:m,i)
      //C
      dgemv("No transpose", m - i + 1, i - 1, -one, a(i, 1), lda, y(i,
        1), ldy, one, a(i, i), 1);
      dgemv("No transpose", m - i + 1, i - 1, -one, x(i, 1), ldx, a(1,
        i), 1, one, a(i, i), 1);
      //C
      //C           Generate reflection Q(i) to annihilate A(i+1:m,i)
      //C
      dlarfg(cmn, m - i + 1, a(i, i), a(fem::min(i + 1, m), i), 1, tauq(i));
      d(i) = a(i, i);
      if (i < n) {
        a(i, i) = one;
        //C
        //C              Compute Y(i+1:n,i)
        //C
        dgemv("Transpose", m - i + 1, n - i, one, a(i, i + 1), lda, a(i,
          i), 1, zero, y(i + 1, i), 1);
        dgemv("Transpose", m - i + 1, i - 1, one, a(i, 1), lda, a(i,
          i), 1, zero, y(1, i), 1);
        dgemv("No transpose", n - i, i - 1, -one, y(i + 1, 1), ldy, y(1,
          i), 1, one, y(i + 1, i), 1);
        dgemv("Transpose", m - i + 1, i - 1, one, x(i, 1), ldx, a(i,
          i), 1, zero, y(1, i), 1);
        dgemv("Transpose", i - 1, n - i, -one, a(1, i + 1), lda, y(1,
          i), 1, one, y(i + 1, i), 1);
        dscal(n - i, tauq(i), y(i + 1, i), 1);
        //C
        //C              Update A(i,i+1:n)
        //C
        dgemv("No transpose", n - i, i, -one, y(i + 1, 1), ldy, a(i,
          1), lda, one, a(i, i + 1), lda);
        dgemv("Transpose", i - 1, n - i, -one, a(1, i + 1), lda, x(i,
          1), ldx, one, a(i, i + 1), lda);
        //C
        //C              Generate reflection P(i) to annihilate A(i,i+2:n)
        //C
        dlarfg(cmn, n - i, a(i, i + 1), a(i, fem::min(i + 2, n)), lda, taup(i));
        e(i) = a(i, i + 1);
        a(i, i + 1) = one;
        //C
        //C              Compute X(i+1:m,i)
        //C
        dgemv("No transpose", m - i, n - i, one, a(i + 1, i + 1),
          lda, a(i, i + 1), lda, zero, x(i + 1, i), 1);
        dgemv("Transpose", n - i, i, one, y(i + 1, 1), ldy, a(i, i + 1),
          lda, zero, x(1, i), 1);
        dgemv("No transpose", m - i, i, -one, a(i + 1, 1), lda, x(1,
          i), 1, one, x(i + 1, i), 1);
        dgemv("No transpose", i - 1, n - i, one, a(1, i + 1), lda, a(i,
          i + 1), lda, zero, x(1, i), 1);
        dgemv("No transpose", m - i, i - 1, -one, x(i + 1, 1), ldx, x(1,
          i), 1, one, x(i + 1, i), 1);
        dscal(m - i, taup(i), x(i + 1, i), 1);
      }
    }
  }
  else {
    //C
    //C        Reduce to lower bidiagonal form
    //C
    FEM_DO(i, 1, nb) {
      //C
      //C           Update A(i,i:n)
      //C
      dgemv("No transpose", n - i + 1, i - 1, -one, y(i, 1), ldy, a(i,
        1), lda, one, a(i, i), lda);
      dgemv("Transpose", i - 1, n - i + 1, -one, a(1, i), lda, x(i,
        1), ldx, one, a(i, i), lda);
      //C
      //C           Generate reflection P(i) to annihilate A(i,i+1:n)
      //C
      dlarfg(cmn, n - i + 1, a(i, i), a(i, fem::min(i + 1, n)), lda, taup(i));
      d(i) = a(i, i);
      if (i < m) {
        a(i, i) = one;
        //C
        //C              Compute X(i+1:m,i)
        //C
        dgemv("No transpose", m - i, n - i + 1, one, a(i + 1, i),
          lda, a(i, i), lda, zero, x(i + 1, i), 1);
        dgemv("Transpose", n - i + 1, i - 1, one, y(i, 1), ldy, a(i,
          i), lda, zero, x(1, i), 1);
        dgemv("No transpose", m - i, i - 1, -one, a(i + 1, 1), lda, x(1,
          i), 1, one, x(i + 1, i), 1);
        dgemv("No transpose", i - 1, n - i + 1, one, a(1, i), lda, a(i,
          i), lda, zero, x(1, i), 1);
        dgemv("No transpose", m - i, i - 1, -one, x(i + 1, 1), ldx, x(1,
          i), 1, one, x(i + 1, i), 1);
        dscal(m - i, taup(i), x(i + 1, i), 1);
        //C
        //C              Update A(i+1:m,i)
        //C
        dgemv("No transpose", m - i, i - 1, -one, a(i + 1, 1), lda, y(i,
          1), ldy, one, a(i + 1, i), 1);
        dgemv("No transpose", m - i, i, -one, x(i + 1, 1), ldx, a(1,
          i), 1, one, a(i + 1, i), 1);
        //C
        //C              Generate reflection Q(i) to annihilate A(i+2:m,i)
        //C
        dlarfg(cmn, m - i, a(i + 1, i), a(fem::min(i + 2, m), i), 1, tauq(i));
        e(i) = a(i + 1, i);
        a(i + 1, i) = one;
        //C
        //C              Compute Y(i+1:n,i)
        //C
        dgemv("Transpose", m - i, n - i, one, a(i + 1, i + 1), lda, a(i + 1,
          i), 1, zero, y(i + 1, i), 1);
        dgemv("Transpose", m - i, i - 1, one, a(i + 1, 1), lda, a(i + 1,
          i), 1, zero, y(1, i), 1);
        dgemv("No transpose", n - i, i - 1, -one, y(i + 1, 1), ldy, y(1,
          i), 1, one, y(i + 1, i), 1);
        dgemv("Transpose", m - i, i, one, x(i + 1, 1), ldx, a(i + 1,
          i), 1, zero, y(1, i), 1);
        dgemv("Transpose", i, n - i, -one, a(1, i + 1), lda, y(1, i),
          1, one, y(i + 1, i), 1);
        dscal(n - i, tauq(i), y(i + 1, i), 1);
      }
    }
  }
  //C
  //C     End of DLABRD
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dgebrd.f
inline
void
dgebrd(
  common& cmn,
  int const& m,
  int const& n,
  arr_ref<double, 2> a,
  int const& lda,
  arr_ref<double> d,
  arr_ref<double> e,
  arr_ref<double> tauq,
  arr_ref<double> taup,
  arr_ref<double> work,
  int const& lwork,
  int& info)
{
  a(dimension(lda, star));
  d(dimension(star));
  e(dimension(star));
  tauq(dimension(star));
  taup(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DGEBRD reduces a general real M-by-N matrix A to upper or lower
  //C  bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
  //C
  //C  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  M       (input) INTEGER
  //C          The number of rows in the matrix A.  M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns in the matrix A.  N >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the M-by-N general matrix to be reduced.
  //C          On exit,
  //C          if m >= n, the diagonal and the first superdiagonal are
  //C            overwritten with the upper bidiagonal matrix B; the
  //C            elements below the diagonal, with the array TAUQ, represent
  //C            the orthogonal matrix Q as a product of elementary
  //C            reflectors, and the elements above the first superdiagonal,
  //C            with the array TAUP, represent the orthogonal matrix P as
  //C            a product of elementary reflectors;
  //C          if m < n, the diagonal and the first subdiagonal are
  //C            overwritten with the lower bidiagonal matrix B; the
  //C            elements below the first subdiagonal, with the array TAUQ,
  //C            represent the orthogonal matrix Q as a product of
  //C            elementary reflectors, and the elements above the diagonal,
  //C            with the array TAUP, represent the orthogonal matrix P as
  //C            a product of elementary reflectors.
  //C          See Further Details.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(1,M).
  //C
  //C  D       (output) DOUBLE PRECISION array, dimension (min(M,N))
  //C          The diagonal elements of the bidiagonal matrix B:
  //C          D(i) = A(i,i).
  //C
  //C  E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)
  //C          The off-diagonal elements of the bidiagonal matrix B:
  //C          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
  //C          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
  //C
  //C  TAUQ    (output) DOUBLE PRECISION array dimension (min(M,N))
  //C          The scalar factors of the elementary reflectors which
  //C          represent the orthogonal matrix Q. See Further Details.
  //C
  //C  TAUP    (output) DOUBLE PRECISION array, dimension (min(M,N))
  //C          The scalar factors of the elementary reflectors which
  //C          represent the orthogonal matrix P. See Further Details.
  //C
  //C  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  //C
  //C  LWORK   (input) INTEGER
  //C          The length of the array WORK.  LWORK >= max(1,M,N).
  //C          For optimum performance LWORK >= (M+N)*NB, where NB
  //C          is the optimal blocksize.
  //C
  //C          If LWORK = -1, then a workspace query is assumed; the routine
  //C          only calculates the optimal size of the WORK array, returns
  //C          this value as the first entry of the WORK array, and no error
  //C          message related to LWORK is issued by XERBLA.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  The matrices Q and P are represented as products of elementary
  //C  reflectors:
  //C
  //C  If m >= n,
  //C
  //C     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
  //C
  //C  Each H(i) and G(i) has the form:
  //C
  //C     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
  //C
  //C  where tauq and taup are real scalars, and v and u are real vectors;
  //C  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
  //C  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
  //C  tauq is stored in TAUQ(i) and taup in TAUP(i).
  //C
  //C  If m < n,
  //C
  //C     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
  //C
  //C  Each H(i) and G(i) has the form:
  //C
  //C     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
  //C
  //C  where tauq and taup are real scalars, and v and u are real vectors;
  //C  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
  //C  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
  //C  tauq is stored in TAUQ(i) and taup in TAUP(i).
  //C
  //C  The contents of A on exit are illustrated by the following examples:
  //C
  //C  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
  //C
  //C    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
  //C    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
  //C    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
  //C    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
  //C    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
  //C    (  v1  v2  v3  v4  v5 )
  //C
  //C  where d and e denote diagonal and off-diagonal elements of B, vi
  //C  denotes an element of the vector defining H(i), and ui an element of
  //C  the vector defining G(i).
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters
  //C
  info = 0;
  int nb = fem::max(1, ilaenv(1, "DGEBRD", " ", m, n, -1, -1));
  int lwkopt = (m + n) * nb;
  work(1) = fem::dble(lwkopt);
  bool lquery = (lwork ==  - 1);
  if (m < 0) {
    info = -1;
  }
  else if (n < 0) {
    info = -2;
  }
  else if (lda < fem::max(1, m)) {
    info = -4;
  }
  else if (lwork < fem::max(1, m, n) && !lquery) {
    info = -10;
  }
  if (info < 0) {
    xerbla("DGEBRD", -info);
    return;
  }
  else if (lquery) {
    return;
  }
  //C
  //C     Quick return if possible
  //C
  int minmn = fem::min(m, n);
  if (minmn == 0) {
    work(1) = 1;
    return;
  }
  //C
  double ws = fem::max(m, n);
  int ldwrkx = m;
  int ldwrky = n;
  //C
  int nx = fem::int0;
  int nbmin = fem::int0;
  if (nb > 1 && nb < minmn) {
    //C
    //C        Set the crossover point NX.
    //C
    nx = fem::max(nb, ilaenv(3, "DGEBRD", " ", m, n, -1, -1));
    //C
    //C        Determine when to switch from blocked to unblocked code.
    //C
    if (nx < minmn) {
      ws = (m + n) * nb;
      if (lwork < ws) {
        //C
        //C              Not enough work space for the optimal NB, consider using
        //C              a smaller block size.
        //C
        nbmin = ilaenv(2, "DGEBRD", " ", m, n, -1, -1);
        if (lwork >= (m + n) * nbmin) {
          nb = lwork / (m + n);
        }
        else {
          nb = 1;
          nx = minmn;
        }
      }
    }
  }
  else {
    nx = minmn;
  }
  //C
  int i = fem::int0;
  const double one = 1.0e+0;
  int j = fem::int0;
  FEM_DOSTEP(i, 1, minmn - nx, nb) {
    //C
    //C        Reduce rows and columns i:i+nb-1 to bidiagonal form and return
    //C        the matrices X and Y which are needed to update the unreduced
    //C        part of the matrix
    //C
    dlabrd(cmn, m - i + 1, n - i + 1, nb, a(i, i), lda, d(i), e(i),
      tauq(i), taup(i), work, ldwrkx, work(ldwrkx * nb + 1), ldwrky);
    //C
    //C        Update the trailing submatrix A(i+nb:m,i+nb:n), using an update
    //C        of the form  A := A - V*Y' - X*U'
    //C
    dgemm("No transpose", "Transpose", m - i - nb + 1, n - i - nb + 1,
      nb, -one, a(i + nb, i), lda, work(ldwrkx * nb + nb + 1),
      ldwrky, one, a(i + nb, i + nb), lda);
    dgemm("No transpose", "No transpose", m - i - nb + 1, n - i - nb + 1,
      nb, -one, work(nb + 1), ldwrkx, a(i, i + nb), lda, one, a(i + nb,
      i + nb), lda);
    //C
    //C        Copy diagonal and off-diagonal elements of B back into A
    //C
    if (m >= n) {
      {
        int fem_do_last = i + nb - 1;
        FEM_DO(j, i, fem_do_last) {
          a(j, j) = d(j);
          a(j, j + 1) = e(j);
        }
      }
    }
    else {
      {
        int fem_do_last = i + nb - 1;
        FEM_DO(j, i, fem_do_last) {
          a(j, j) = d(j);
          a(j + 1, j) = e(j);
        }
      }
    }
  }
  //C
  //C     Use unblocked code to reduce the remainder of the matrix
  //C
  int iinfo = fem::int0;
  dgebd2(cmn, m - i + 1, n - i + 1, a(i, i), lda, d(i), e(i), tauq(i),
    taup(i), work, iinfo);
  work(1) = ws;
  //C
  //C     End of DGEBRD
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlarfp.f
inline
void
dlarfp(
  common& cmn,
  int const& n,
  double& alpha,
  arr_ref<double> x,
  int const& incx,
  double& tau)
{
  x(dimension(star));
  const double zero = 0.0e+0;
  double xnorm = fem::double0;
  const double two = 2.0e+0;
  int j = fem::int0;
  double beta = fem::double0;
  double safmin = fem::double0;
  int knt = fem::int0;
  const double one = 1.0e+0;
  double rsafmn = fem::double0;
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLARFP generates a real elementary reflector H of order n, such
  //C  that
  //C
  //C        H * ( alpha ) = ( beta ),   H' * H = I.
  //C            (   x   )   (   0  )
  //C
  //C  where alpha and beta are scalars, beta is non-negative, and x is
  //C  an (n-1)-element real vector.  H is represented in the form
  //C
  //C        H = I - tau * ( 1 ) * ( 1 v' ) ,
  //C                      ( v )
  //C
  //C  where tau is a real scalar and v is a real (n-1)-element
  //C  vector.
  //C
  //C  If the elements of x are all zero, then tau = 0 and H is taken to be
  //C  the unit matrix.
  //C
  //C  Otherwise  1 <= tau <= 2.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  N       (input) INTEGER
  //C          The order of the elementary reflector.
  //C
  //C  ALPHA   (input/output) DOUBLE PRECISION
  //C          On entry, the value alpha.
  //C          On exit, it is overwritten with the value beta.
  //C
  //C  X       (input/output) DOUBLE PRECISION array, dimension
  //C                         (1+(N-2)*abs(INCX))
  //C          On entry, the vector x.
  //C          On exit, it is overwritten with the vector v.
  //C
  //C  INCX    (input) INTEGER
  //C          The increment between elements of X. INCX > 0.
  //C
  //C  TAU     (output) DOUBLE PRECISION
  //C          The value tau.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  if (n <= 0) {
    tau = zero;
    return;
  }
  //C
  xnorm = dnrm2(n - 1, x, incx);
  //C
  if (xnorm == zero) {
    //C
    //C        H  =  [+/-1, 0; I], sign chosen so ALPHA >= 0
    //C
    if (alpha >= zero) {
      //C           When TAU.eq.ZERO, the vector is special-cased to be
      //C           all zeros in the application routines.  We do not need
      //C           to clear it.
      tau = zero;
    }
    else {
      //C           However, the application routines rely on explicit
      //C           zero checks when TAU.ne.ZERO, and we must clear X.
      tau = two;
      {
        int fem_do_last = n - 1;
        FEM_DO(j, 1, fem_do_last) {
          x(1 + (j - 1) * incx) = 0;
        }
      }
      alpha = -alpha;
    }
  }
  else {
    //C
    //C        general case
    //C
    beta = fem::sign(dlapy2(alpha, xnorm), alpha);
    safmin = dlamch(cmn, "S") / dlamch(cmn, "E");
    knt = 0;
    if (fem::abs(beta) < safmin) {
      //C
      //C           XNORM, BETA may be inaccurate; scale X and recompute them
      //C
      rsafmn = one / safmin;
      statement_10:
      knt++;
      dscal(n - 1, rsafmn, x, incx);
      beta = beta * rsafmn;
      alpha = alpha * rsafmn;
      if (fem::abs(beta) < safmin) {
        goto statement_10;
      }
      //C
      //C           New BETA is at most 1, at least SAFMIN
      //C
      xnorm = dnrm2(n - 1, x, incx);
      beta = fem::sign(dlapy2(alpha, xnorm), alpha);
    }
    alpha += beta;
    if (beta < zero) {
      beta = -beta;
      tau = -alpha / beta;
    }
    else {
      alpha = xnorm * (xnorm / alpha);
      tau = alpha / beta;
      alpha = -alpha;
    }
    dscal(n - 1, one / alpha, x, incx);
    //C
    //C        If BETA is subnormal, it may lose relative accuracy
    //C
    FEM_DO(j, 1, knt) {
      beta = beta * safmin;
    }
    alpha = beta;
  }
  //C
  //C     End of DLARFP
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dgelq2.f
inline
void
dgelq2(
  common& cmn,
  int const& m,
  int const& n,
  arr_ref<double, 2> a,
  int const& lda,
  arr_ref<double> tau,
  arr_ref<double> work,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DGELQ2 computes an LQ factorization of a real m by n matrix A:
  //C  A = L * Q.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix A.  M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix A.  N >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the m by n matrix A.
  //C          On exit, the elements on and below the diagonal of the array
  //C          contain the m by min(m,n) lower trapezoidal matrix L (L is
  //C          lower triangular if m <= n); the elements above the diagonal,
  //C          with the array TAU, represent the orthogonal matrix Q as a
  //C          product of elementary reflectors (see Further Details).
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(1,M).
  //C
  //C  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
  //C          The scalar factors of the elementary reflectors (see Further
  //C          Details).
  //C
  //C  WORK    (workspace) DOUBLE PRECISION array, dimension (M)
  //C
  //C  INFO    (output) INTEGER
  //C          = 0: successful exit
  //C          < 0: if INFO = -i, the i-th argument had an illegal value
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  The matrix Q is represented as a product of elementary reflectors
  //C
  //C     Q = H(k) . . . H(2) H(1), where k = min(m,n).
  //C
  //C  Each H(i) has the form
  //C
  //C     H(i) = I - tau * v * v'
  //C
  //C  where tau is a real scalar, and v is a real vector with
  //C  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
  //C  and tau in TAU(i).
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  if (m < 0) {
    info = -1;
  }
  else if (n < 0) {
    info = -2;
  }
  else if (lda < fem::max(1, m)) {
    info = -4;
  }
  if (info != 0) {
    xerbla("DGELQ2", -info);
    return;
  }
  //C
  int k = fem::min(m, n);
  //C
  int i = fem::int0;
  double aii = fem::double0;
  const double one = 1.0e+0;
  FEM_DO(i, 1, k) {
    //C
    //C        Generate elementary reflector H(i) to annihilate A(i,i+1:n)
    //C
    dlarfp(cmn, n - i + 1, a(i, i), a(i, fem::min(i + 1, n)), lda, tau(i));
    if (i < m) {
      //C
      //C           Apply H(i) to A(i+1:m,i:n) from the right
      //C
      aii = a(i, i);
      a(i, i) = one;
      dlarf("Right", m - i, n - i + 1, a(i, i), lda, tau(i), a(i + 1,
        i), lda, work);
      a(i, i) = aii;
    }
  }
  //C
  //C     End of DGELQ2
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlarfb.f
inline
void
dlarfb(
  str_cref side,
  str_cref trans,
  str_cref direct,
  str_cref storev,
  int const& m,
  int const& n,
  int const& k,
  arr_cref<double, 2> v,
  int const& ldv,
  arr_cref<double, 2> t,
  int const& ldt,
  arr_ref<double, 2> c,
  int const& ldc,
  arr_ref<double, 2> work,
  int const& ldwork)
{
  v(dimension(ldv, star));
  t(dimension(ldt, star));
  c(dimension(ldc, star));
  work(dimension(ldwork, star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLARFB applies a real block reflector H or its transpose H' to a
  //C  real m by n matrix C, from either the left or the right.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  SIDE    (input) CHARACTER*1
  //C          = 'L': apply H or H' from the Left
  //C          = 'R': apply H or H' from the Right
  //C
  //C  TRANS   (input) CHARACTER*1
  //C          = 'N': apply H (No transpose)
  //C          = 'T': apply H' (Transpose)
  //C
  //C  DIRECT  (input) CHARACTER*1
  //C          Indicates how H is formed from a product of elementary
  //C          reflectors
  //C          = 'F': H = H(1) H(2) . . . H(k) (Forward)
  //C          = 'B': H = H(k) . . . H(2) H(1) (Backward)
  //C
  //C  STOREV  (input) CHARACTER*1
  //C          Indicates how the vectors which define the elementary
  //C          reflectors are stored:
  //C          = 'C': Columnwise
  //C          = 'R': Rowwise
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix C.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix C.
  //C
  //C  K       (input) INTEGER
  //C          The order of the matrix T (= the number of elementary
  //C          reflectors whose product defines the block reflector).
  //C
  //C  V       (input) DOUBLE PRECISION array, dimension
  //C                                (LDV,K) if STOREV = 'C'
  //C                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
  //C                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
  //C          The matrix V. See further details.
  //C
  //C  LDV     (input) INTEGER
  //C          The leading dimension of the array V.
  //C          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
  //C          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
  //C          if STOREV = 'R', LDV >= K.
  //C
  //C  T       (input) DOUBLE PRECISION array, dimension (LDT,K)
  //C          The triangular k by k matrix T in the representation of the
  //C          block reflector.
  //C
  //C  LDT     (input) INTEGER
  //C          The leading dimension of the array T. LDT >= K.
  //C
  //C  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
  //C          On entry, the m by n matrix C.
  //C          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
  //C
  //C  LDC     (input) INTEGER
  //C          The leading dimension of the array C. LDA >= max(1,M).
  //C
  //C  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)
  //C
  //C  LDWORK  (input) INTEGER
  //C          The leading dimension of the array WORK.
  //C          If SIDE = 'L', LDWORK >= max(1,N);
  //C          if SIDE = 'R', LDWORK >= max(1,M).
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Quick return if possible
  //C
  if (m <= 0 || n <= 0) {
    return;
  }
  //C
  fem::str<1> transt = fem::char0;
  if (lsame(trans, "N")) {
    transt = "T";
  }
  else {
    transt = "N";
  }
  //C
  int lastv = fem::int0;
  int lastc = fem::int0;
  int j = fem::int0;
  const double one = 1.0e+0;
  int i = fem::int0;
  if (lsame(storev, "C")) {
    //C
    if (lsame(direct, "F")) {
      //C
      //C           Let  V =  ( V1 )    (first K rows)
      //C                     ( V2 )
      //C           where  V1  is unit lower triangular.
      //C
      if (lsame(side, "L")) {
        //C
        //C              Form  H * C  or  H' * C  where  C = ( C1 )
        //C                                                  ( C2 )
        //C
        lastv = fem::max(k, iladlr(m, k, v, ldv));
        lastc = iladlc(lastv, n, c, ldc);
        //C
        //C              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
        //C
        //C              W := C1'
        //C
        FEM_DO(j, 1, k) {
          dcopy(lastc, c(j, 1), ldc, work(1, j), 1);
        }
        //C
        //C              W := W * V1
        //C
        dtrmm("Right", "Lower", "No transpose", "Unit", lastc, k,
          one, v, ldv, work, ldwork);
        if (lastv > k) {
          //C
          //C                 W := W + C2'*V2
          //C
          dgemm("Transpose", "No transpose", lastc, k, lastv - k,
            one, c(k + 1, 1), ldc, v(k + 1, 1), ldv, one, work,
            ldwork);
        }
        //C
        //C              W := W * T'  or  W * T
        //C
        dtrmm("Right", "Upper", transt, "Non-unit", lastc, k, one, t,
          ldt, work, ldwork);
        //C
        //C              C := C - V * W'
        //C
        if (lastv > k) {
          //C
          //C                 C2 := C2 - V2 * W'
          //C
          dgemm("No transpose", "Transpose", lastv - k, lastc, k,
            -one, v(k + 1, 1), ldv, work, ldwork, one, c(k + 1, 1),
            ldc);
        }
        //C
        //C              W := W * V1'
        //C
        dtrmm("Right", "Lower", "Transpose", "Unit", lastc, k, one,
          v, ldv, work, ldwork);
        //C
        //C              C1 := C1 - W'
        //C
        FEM_DO(j, 1, k) {
          FEM_DO(i, 1, lastc) {
            c(j, i) = c(j, i) - work(i, j);
          }
        }
        //C
      }
      else if (lsame(side, "R")) {
        //C
        //C              Form  C * H  or  C * H'  where  C = ( C1  C2 )
        //C
        lastv = fem::max(k, iladlr(n, k, v, ldv));
        lastc = iladlr(m, lastv, c, ldc);
        //C
        //C              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
        //C
        //C              W := C1
        //C
        FEM_DO(j, 1, k) {
          dcopy(lastc, c(1, j), 1, work(1, j), 1);
        }
        //C
        //C              W := W * V1
        //C
        dtrmm("Right", "Lower", "No transpose", "Unit", lastc, k,
          one, v, ldv, work, ldwork);
        if (lastv > k) {
          //C
          //C                 W := W + C2 * V2
          //C
          dgemm("No transpose", "No transpose", lastc, k, lastv - k,
            one, c(1, k + 1), ldc, v(k + 1, 1), ldv, one, work,
            ldwork);
        }
        //C
        //C              W := W * T  or  W * T'
        //C
        dtrmm("Right", "Upper", trans, "Non-unit", lastc, k, one, t,
          ldt, work, ldwork);
        //C
        //C              C := C - W * V'
        //C
        if (lastv > k) {
          //C
          //C                 C2 := C2 - W * V2'
          //C
          dgemm("No transpose", "Transpose", lastc, lastv - k, k,
            -one, work, ldwork, v(k + 1, 1), ldv, one, c(1, k + 1),
            ldc);
        }
        //C
        //C              W := W * V1'
        //C
        dtrmm("Right", "Lower", "Transpose", "Unit", lastc, k, one,
          v, ldv, work, ldwork);
        //C
        //C              C1 := C1 - W
        //C
        FEM_DO(j, 1, k) {
          FEM_DO(i, 1, lastc) {
            c(i, j) = c(i, j) - work(i, j);
          }
        }
      }
      //C
    }
    else {
      //C
      //C           Let  V =  ( V1 )
      //C                     ( V2 )    (last K rows)
      //C           where  V2  is unit upper triangular.
      //C
      if (lsame(side, "L")) {
        //C
        //C              Form  H * C  or  H' * C  where  C = ( C1 )
        //C                                                  ( C2 )
        //C
        lastv = fem::max(k, iladlr(m, k, v, ldv));
        lastc = iladlc(lastv, n, c, ldc);
        //C
        //C              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
        //C
        //C              W := C2'
        //C
        FEM_DO(j, 1, k) {
          dcopy(lastc, c(lastv - k + j, 1), ldc, work(1, j), 1);
        }
        //C
        //C              W := W * V2
        //C
        dtrmm("Right", "Upper", "No transpose", "Unit", lastc, k,
          one, v(lastv - k + 1, 1), ldv, work, ldwork);
        if (lastv > k) {
          //C
          //C                 W := W + C1'*V1
          //C
          dgemm("Transpose", "No transpose", lastc, k, lastv - k,
            one, c, ldc, v, ldv, one, work, ldwork);
        }
        //C
        //C              W := W * T'  or  W * T
        //C
        dtrmm("Right", "Lower", transt, "Non-unit", lastc, k, one, t,
          ldt, work, ldwork);
        //C
        //C              C := C - V * W'
        //C
        if (lastv > k) {
          //C
          //C                 C1 := C1 - V1 * W'
          //C
          dgemm("No transpose", "Transpose", lastv - k, lastc, k,
            -one, v, ldv, work, ldwork, one, c, ldc);
        }
        //C
        //C              W := W * V2'
        //C
        dtrmm("Right", "Upper", "Transpose", "Unit", lastc, k, one, v(
          lastv - k + 1, 1), ldv, work, ldwork);
        //C
        //C              C2 := C2 - W'
        //C
        FEM_DO(j, 1, k) {
          FEM_DO(i, 1, lastc) {
            c(lastv - k + j, i) = c(lastv - k + j, i) - work(i, j);
          }
        }
        //C
      }
      else if (lsame(side, "R")) {
        //C
        //C              Form  C * H  or  C * H'  where  C = ( C1  C2 )
        //C
        lastv = fem::max(k, iladlr(n, k, v, ldv));
        lastc = iladlr(m, lastv, c, ldc);
        //C
        //C              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
        //C
        //C              W := C2
        //C
        FEM_DO(j, 1, k) {
          dcopy(lastc, c(1, n - k + j), 1, work(1, j), 1);
        }
        //C
        //C              W := W * V2
        //C
        dtrmm("Right", "Upper", "No transpose", "Unit", lastc, k,
          one, v(lastv - k + 1, 1), ldv, work, ldwork);
        if (lastv > k) {
          //C
          //C                 W := W + C1 * V1
          //C
          dgemm("No transpose", "No transpose", lastc, k, lastv - k,
            one, c, ldc, v, ldv, one, work, ldwork);
        }
        //C
        //C              W := W * T  or  W * T'
        //C
        dtrmm("Right", "Lower", trans, "Non-unit", lastc, k, one, t,
          ldt, work, ldwork);
        //C
        //C              C := C - W * V'
        //C
        if (lastv > k) {
          //C
          //C                 C1 := C1 - W * V1'
          //C
          dgemm("No transpose", "Transpose", lastc, lastv - k, k,
            -one, work, ldwork, v, ldv, one, c, ldc);
        }
        //C
        //C              W := W * V2'
        //C
        dtrmm("Right", "Upper", "Transpose", "Unit", lastc, k, one, v(
          lastv - k + 1, 1), ldv, work, ldwork);
        //C
        //C              C2 := C2 - W
        //C
        FEM_DO(j, 1, k) {
          FEM_DO(i, 1, lastc) {
            c(i, lastv - k + j) = c(i, lastv - k + j) - work(i, j);
          }
        }
      }
    }
    //C
  }
  else if (lsame(storev, "R")) {
    //C
    if (lsame(direct, "F")) {
      //C
      //C           Let  V =  ( V1  V2 )    (V1: first K columns)
      //C           where  V1  is unit upper triangular.
      //C
      if (lsame(side, "L")) {
        //C
        //C              Form  H * C  or  H' * C  where  C = ( C1 )
        //C                                                  ( C2 )
        //C
        lastv = fem::max(k, iladlc(k, m, v, ldv));
        lastc = iladlc(lastv, n, c, ldc);
        //C
        //C              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
        //C
        //C              W := C1'
        //C
        FEM_DO(j, 1, k) {
          dcopy(lastc, c(j, 1), ldc, work(1, j), 1);
        }
        //C
        //C              W := W * V1'
        //C
        dtrmm("Right", "Upper", "Transpose", "Unit", lastc, k, one,
          v, ldv, work, ldwork);
        if (lastv > k) {
          //C
          //C                 W := W + C2'*V2'
          //C
          dgemm("Transpose", "Transpose", lastc, k, lastv - k, one, c(k + 1,
            1), ldc, v(1, k + 1), ldv, one, work, ldwork);
        }
        //C
        //C              W := W * T'  or  W * T
        //C
        dtrmm("Right", "Upper", transt, "Non-unit", lastc, k, one, t,
          ldt, work, ldwork);
        //C
        //C              C := C - V' * W'
        //C
        if (lastv > k) {
          //C
          //C                 C2 := C2 - V2' * W'
          //C
          dgemm("Transpose", "Transpose", lastv - k, lastc, k, -one,
            v(1, k + 1), ldv, work, ldwork, one, c(k + 1, 1), ldc);
        }
        //C
        //C              W := W * V1
        //C
        dtrmm("Right", "Upper", "No transpose", "Unit", lastc, k,
          one, v, ldv, work, ldwork);
        //C
        //C              C1 := C1 - W'
        //C
        FEM_DO(j, 1, k) {
          FEM_DO(i, 1, lastc) {
            c(j, i) = c(j, i) - work(i, j);
          }
        }
        //C
      }
      else if (lsame(side, "R")) {
        //C
        //C              Form  C * H  or  C * H'  where  C = ( C1  C2 )
        //C
        lastv = fem::max(k, iladlc(k, n, v, ldv));
        lastc = iladlr(m, lastv, c, ldc);
        //C
        //C              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
        //C
        //C              W := C1
        //C
        FEM_DO(j, 1, k) {
          dcopy(lastc, c(1, j), 1, work(1, j), 1);
        }
        //C
        //C              W := W * V1'
        //C
        dtrmm("Right", "Upper", "Transpose", "Unit", lastc, k, one,
          v, ldv, work, ldwork);
        if (lastv > k) {
          //C
          //C                 W := W + C2 * V2'
          //C
          dgemm("No transpose", "Transpose", lastc, k, lastv - k,
            one, c(1, k + 1), ldc, v(1, k + 1), ldv, one, work,
            ldwork);
        }
        //C
        //C              W := W * T  or  W * T'
        //C
        dtrmm("Right", "Upper", trans, "Non-unit", lastc, k, one, t,
          ldt, work, ldwork);
        //C
        //C              C := C - W * V
        //C
        if (lastv > k) {
          //C
          //C                 C2 := C2 - W * V2
          //C
          dgemm("No transpose", "No transpose", lastc, lastv - k, k,
            -one, work, ldwork, v(1, k + 1), ldv, one, c(1, k + 1),
            ldc);
        }
        //C
        //C              W := W * V1
        //C
        dtrmm("Right", "Upper", "No transpose", "Unit", lastc, k,
          one, v, ldv, work, ldwork);
        //C
        //C              C1 := C1 - W
        //C
        FEM_DO(j, 1, k) {
          FEM_DO(i, 1, lastc) {
            c(i, j) = c(i, j) - work(i, j);
          }
        }
        //C
      }
      //C
    }
    else {
      //C
      //C           Let  V =  ( V1  V2 )    (V2: last K columns)
      //C           where  V2  is unit lower triangular.
      //C
      if (lsame(side, "L")) {
        //C
        //C              Form  H * C  or  H' * C  where  C = ( C1 )
        //C                                                  ( C2 )
        //C
        lastv = fem::max(k, iladlc(k, m, v, ldv));
        lastc = iladlc(lastv, n, c, ldc);
        //C
        //C              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
        //C
        //C              W := C2'
        //C
        FEM_DO(j, 1, k) {
          dcopy(lastc, c(lastv - k + j, 1), ldc, work(1, j), 1);
        }
        //C
        //C              W := W * V2'
        //C
        dtrmm("Right", "Lower", "Transpose", "Unit", lastc, k, one, v(1,
          lastv - k + 1), ldv, work, ldwork);
        if (lastv > k) {
          //C
          //C                 W := W + C1'*V1'
          //C
          dgemm("Transpose", "Transpose", lastc, k, lastv - k, one,
            c, ldc, v, ldv, one, work, ldwork);
        }
        //C
        //C              W := W * T'  or  W * T
        //C
        dtrmm("Right", "Lower", transt, "Non-unit", lastc, k, one, t,
          ldt, work, ldwork);
        //C
        //C              C := C - V' * W'
        //C
        if (lastv > k) {
          //C
          //C                 C1 := C1 - V1' * W'
          //C
          dgemm("Transpose", "Transpose", lastv - k, lastc, k, -one,
            v, ldv, work, ldwork, one, c, ldc);
        }
        //C
        //C              W := W * V2
        //C
        dtrmm("Right", "Lower", "No transpose", "Unit", lastc, k,
          one, v(1, lastv - k + 1), ldv, work, ldwork);
        //C
        //C              C2 := C2 - W'
        //C
        FEM_DO(j, 1, k) {
          FEM_DO(i, 1, lastc) {
            c(lastv - k + j, i) = c(lastv - k + j, i) - work(i, j);
          }
        }
        //C
      }
      else if (lsame(side, "R")) {
        //C
        //C              Form  C * H  or  C * H'  where  C = ( C1  C2 )
        //C
        lastv = fem::max(k, iladlc(k, n, v, ldv));
        lastc = iladlr(m, lastv, c, ldc);
        //C
        //C              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
        //C
        //C              W := C2
        //C
        FEM_DO(j, 1, k) {
          dcopy(lastc, c(1, lastv - k + j), 1, work(1, j), 1);
        }
        //C
        //C              W := W * V2'
        //C
        dtrmm("Right", "Lower", "Transpose", "Unit", lastc, k, one, v(1,
          lastv - k + 1), ldv, work, ldwork);
        if (lastv > k) {
          //C
          //C                 W := W + C1 * V1'
          //C
          dgemm("No transpose", "Transpose", lastc, k, lastv - k,
            one, c, ldc, v, ldv, one, work, ldwork);
        }
        //C
        //C              W := W * T  or  W * T'
        //C
        dtrmm("Right", "Lower", trans, "Non-unit", lastc, k, one, t,
          ldt, work, ldwork);
        //C
        //C              C := C - W * V
        //C
        if (lastv > k) {
          //C
          //C                 C1 := C1 - W * V1
          //C
          dgemm("No transpose", "No transpose", lastc, lastv - k, k,
            -one, work, ldwork, v, ldv, one, c, ldc);
        }
        //C
        //C              W := W * V2
        //C
        dtrmm("Right", "Lower", "No transpose", "Unit", lastc, k,
          one, v(1, lastv - k + 1), ldv, work, ldwork);
        //C
        //C              C1 := C1 - W
        //C
        FEM_DO(j, 1, k) {
          FEM_DO(i, 1, lastc) {
            c(i, lastv - k + j) = c(i, lastv - k + j) - work(i, j);
          }
        }
        //C
      }
      //C
    }
  }
  //C
  //C     End of DLARFB
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlarft.f
inline
void
dlarft(
  str_cref direct,
  str_cref storev,
  int const& n,
  int const& k,
  arr_ref<double, 2> v,
  int const& ldv,
  arr_cref<double> tau,
  arr_ref<double, 2> t,
  int const& ldt)
{
  v(dimension(ldv, star));
  tau(dimension(star));
  t(dimension(ldt, star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLARFT forms the triangular factor T of a real block reflector H
  //C  of order n, which is defined as a product of k elementary reflectors.
  //C
  //C  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
  //C
  //C  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
  //C
  //C  If STOREV = 'C', the vector which defines the elementary reflector
  //C  H(i) is stored in the i-th column of the array V, and
  //C
  //C     H  =  I - V * T * V'
  //C
  //C  If STOREV = 'R', the vector which defines the elementary reflector
  //C  H(i) is stored in the i-th row of the array V, and
  //C
  //C     H  =  I - V' * T * V
  //C
  //C  Arguments
  //C  =========
  //C
  //C  DIRECT  (input) CHARACTER*1
  //C          Specifies the order in which the elementary reflectors are
  //C          multiplied to form the block reflector:
  //C          = 'F': H = H(1) H(2) . . . H(k) (Forward)
  //C          = 'B': H = H(k) . . . H(2) H(1) (Backward)
  //C
  //C  STOREV  (input) CHARACTER*1
  //C          Specifies how the vectors which define the elementary
  //C          reflectors are stored (see also Further Details):
  //C          = 'C': columnwise
  //C          = 'R': rowwise
  //C
  //C  N       (input) INTEGER
  //C          The order of the block reflector H. N >= 0.
  //C
  //C  K       (input) INTEGER
  //C          The order of the triangular factor T (= the number of
  //C          elementary reflectors). K >= 1.
  //C
  //C  V       (input/output) DOUBLE PRECISION array, dimension
  //C                               (LDV,K) if STOREV = 'C'
  //C                               (LDV,N) if STOREV = 'R'
  //C          The matrix V. See further details.
  //C
  //C  LDV     (input) INTEGER
  //C          The leading dimension of the array V.
  //C          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
  //C
  //C  TAU     (input) DOUBLE PRECISION array, dimension (K)
  //C          TAU(i) must contain the scalar factor of the elementary
  //C          reflector H(i).
  //C
  //C  T       (output) DOUBLE PRECISION array, dimension (LDT,K)
  //C          The k by k triangular factor T of the block reflector.
  //C          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
  //C          lower triangular. The rest of the array is not used.
  //C
  //C  LDT     (input) INTEGER
  //C          The leading dimension of the array T. LDT >= K.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  The shape of the matrix V and the storage of the vectors which define
  //C  the H(i) is best illustrated by the following example with n = 5 and
  //C  k = 3. The elements equal to 1 are not stored; the corresponding
  //C  array elements are modified but restored on exit. The rest of the
  //C  array is not used.
  //C
  //C  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
  //C
  //C               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
  //C                   ( v1  1    )                     (     1 v2 v2 v2 )
  //C                   ( v1 v2  1 )                     (        1 v3 v3 )
  //C                   ( v1 v2 v3 )
  //C
  //C  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
  //C
  //C               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
  //C                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
  //C                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
  //C                   (     1 v3 )
  //C                   (        1 )
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Quick return if possible
  //C
  if (n == 0) {
    return;
  }
  //C
  int prevlastv = fem::int0;
  int i = fem::int0;
  const double zero = 0.0e+0;
  int j = fem::int0;
  double vii = fem::double0;
  const double one = 1.0e+0;
  int lastv = fem::int0;
  if (lsame(direct, "F")) {
    prevlastv = n;
    FEM_DO(i, 1, k) {
      prevlastv = fem::max(i, prevlastv);
      if (tau(i) == zero) {
        //C
        //C              H(i)  =  I
        //C
        FEM_DO(j, 1, i) {
          t(j, i) = zero;
        }
      }
      else {
        //C
        //C              general case
        //C
        vii = v(i, i);
        v(i, i) = one;
        if (lsame(storev, "C")) {
          //C                 Skip any trailing zeros.
          FEM_DOSTEP(lastv, n, i + 1, -1) {
            if (v(lastv, i) != zero) {
              break;
            }
          }
          j = fem::min(lastv, prevlastv);
          //C
          //C                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)' * V(i:j,i)
          //C
          dgemv("Transpose", j - i + 1, i - 1, -tau(i), v(i, 1), ldv,
            v(i, i), 1, zero, t(1, i), 1);
        }
        else {
          //C                 Skip any trailing zeros.
          FEM_DOSTEP(lastv, n, i + 1, -1) {
            if (v(i, lastv) != zero) {
              break;
            }
          }
          j = fem::min(lastv, prevlastv);
          //C
          //C                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)'
          //C
          dgemv("No transpose", i - 1, j - i + 1, -tau(i), v(1, i),
            ldv, v(i, i), ldv, zero, t(1, i), 1);
        }
        v(i, i) = vii;
        //C
        //C              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
        //C
        dtrmv("Upper", "No transpose", "Non-unit", i - 1, t, ldt, t(1, i), 1);
        t(i, i) = tau(i);
        if (i > 1) {
          prevlastv = fem::max(prevlastv, lastv);
        }
        else {
          prevlastv = lastv;
        }
      }
    }
  }
  else {
    prevlastv = 1;
    FEM_DOSTEP(i, k, 1, -1) {
      if (tau(i) == zero) {
        //C
        //C              H(i)  =  I
        //C
        FEM_DO(j, i, k) {
          t(j, i) = zero;
        }
      }
      else {
        //C
        //C              general case
        //C
        if (i < k) {
          if (lsame(storev, "C")) {
            vii = v(n - k + i, i);
            v(n - k + i, i) = one;
            //C                    Skip any leading zeros.
            {
              int fem_do_last = i - 1;
              FEM_DO(lastv, 1, fem_do_last) {
                if (v(lastv, i) != zero) {
                  break;
                }
              }
            }
            j = fem::max(lastv, prevlastv);
            //C
            //C                    T(i+1:k,i) :=
            //C                            - tau(i) * V(j:n-k+i,i+1:k)' * V(j:n-k+i,i)
            //C
            dgemv("Transpose", n - k + i - j + 1, k - i, -tau(i), v(j,
              i + 1), ldv, v(j, i), 1, zero, t(i + 1, i), 1);
            v(n - k + i, i) = vii;
          }
          else {
            vii = v(i, n - k + i);
            v(i, n - k + i) = one;
            //C                    Skip any leading zeros.
            {
              int fem_do_last = i - 1;
              FEM_DO(lastv, 1, fem_do_last) {
                if (v(i, lastv) != zero) {
                  break;
                }
              }
            }
            j = fem::max(lastv, prevlastv);
            //C
            //C                    T(i+1:k,i) :=
            //C                            - tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)'
            //C
            dgemv("No transpose", k - i, n - k + i - j + 1, -tau(i),
              v(i + 1, j), ldv, v(i, j), ldv, zero, t(i + 1, i), 1);
            v(i, n - k + i) = vii;
          }
          //C
          //C                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
          //C
          dtrmv("Lower", "No transpose", "Non-unit", k - i, t(i + 1,
            i + 1), ldt, t(i + 1, i), 1);
          if (i > 1) {
            prevlastv = fem::min(prevlastv, lastv);
          }
          else {
            prevlastv = lastv;
          }
        }
        t(i, i) = tau(i);
      }
    }
  }
  //C
  //C     End of DLARFT
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dgelqf.f
inline
void
dgelqf(
  common& cmn,
  int const& m,
  int const& n,
  arr_ref<double, 2> a,
  int const& lda,
  arr_ref<double> tau,
  arr_ref<double> work,
  int const& lwork,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DGELQF computes an LQ factorization of a real M-by-N matrix A:
  //C  A = L * Q.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix A.  M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix A.  N >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the M-by-N matrix A.
  //C          On exit, the elements on and below the diagonal of the array
  //C          contain the m-by-min(m,n) lower trapezoidal matrix L (L is
  //C          lower triangular if m <= n); the elements above the diagonal,
  //C          with the array TAU, represent the orthogonal matrix Q as a
  //C          product of elementary reflectors (see Further Details).
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(1,M).
  //C
  //C  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
  //C          The scalar factors of the elementary reflectors (see Further
  //C          Details).
  //C
  //C  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  //C
  //C  LWORK   (input) INTEGER
  //C          The dimension of the array WORK.  LWORK >= max(1,M).
  //C          For optimum performance LWORK >= M*NB, where NB is the
  //C          optimal blocksize.
  //C
  //C          If LWORK = -1, then a workspace query is assumed; the routine
  //C          only calculates the optimal size of the WORK array, returns
  //C          this value as the first entry of the WORK array, and no error
  //C          message related to LWORK is issued by XERBLA.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  The matrix Q is represented as a product of elementary reflectors
  //C
  //C     Q = H(k) . . . H(2) H(1), where k = min(m,n).
  //C
  //C  Each H(i) has the form
  //C
  //C     H(i) = I - tau * v * v'
  //C
  //C  where tau is a real scalar, and v is a real vector with
  //C  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
  //C  and tau in TAU(i).
  //C
  //C  =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  int nb = ilaenv(1, "DGELQF", " ", m, n, -1, -1);
  int lwkopt = m * nb;
  work(1) = lwkopt;
  bool lquery = (lwork ==  - 1);
  if (m < 0) {
    info = -1;
  }
  else if (n < 0) {
    info = -2;
  }
  else if (lda < fem::max(1, m)) {
    info = -4;
  }
  else if (lwork < fem::max(1, m) && !lquery) {
    info = -7;
  }
  if (info != 0) {
    xerbla("DGELQF", -info);
    return;
  }
  else if (lquery) {
    return;
  }
  //C
  //C     Quick return if possible
  //C
  int k = fem::min(m, n);
  if (k == 0) {
    work(1) = 1;
    return;
  }
  //C
  int nbmin = 2;
  int nx = 0;
  int iws = m;
  int ldwork = fem::int0;
  if (nb > 1 && nb < k) {
    //C
    //C        Determine when to cross over from blocked to unblocked code.
    //C
    nx = fem::max(0, ilaenv(3, "DGELQF", " ", m, n, -1, -1));
    if (nx < k) {
      //C
      //C           Determine if workspace is large enough for blocked code.
      //C
      ldwork = m;
      iws = ldwork * nb;
      if (lwork < iws) {
        //C
        //C              Not enough workspace to use optimal NB:  reduce NB and
        //C              determine the minimum value of NB.
        //C
        nb = lwork / ldwork;
        nbmin = fem::max(2, ilaenv(2, "DGELQF", " ", m, n, -1, -1));
      }
    }
  }
  //C
  int i = fem::int0;
  int ib = fem::int0;
  int iinfo = fem::int0;
  if (nb >= nbmin && nb < k && nx < k) {
    //C
    //C        Use blocked code initially
    //C
    FEM_DOSTEP(i, 1, k - nx, nb) {
      ib = fem::min(k - i + 1, nb);
      //C
      //C           Compute the LQ factorization of the current block
      //C           A(i:i+ib-1,i:n)
      //C
      dgelq2(cmn, ib, n - i + 1, a(i, i), lda, tau(i), work, iinfo);
      if (i + ib <= m) {
        //C
        //C              Form the triangular factor of the block reflector
        //C              H = H(i) H(i+1) . . . H(i+ib-1)
        //C
        dlarft("Forward", "Rowwise", n - i + 1, ib, a(i, i), lda, tau(i),
          work, ldwork);
        //C
        //C              Apply H to A(i+ib:m,i:n) from the right
        //C
        dlarfb("Right", "No transpose", "Forward", "Rowwise", m - i - ib + 1,
          n - i + 1, ib, a(i, i), lda, work, ldwork, a(i + ib, i),
          lda, work(ib + 1), ldwork);
      }
    }
  }
  else {
    i = 1;
  }
  //C
  //C     Use unblocked code to factor the last or only block.
  //C
  if (i <= k) {
    dgelq2(cmn, m - i + 1, n - i + 1, a(i, i), lda, tau(i), work, iinfo);
  }
  //C
  work(1) = iws;
  //C
  //C     End of DGELQF
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dgeqr2.f
inline
void
dgeqr2(
  common& cmn,
  int const& m,
  int const& n,
  arr_ref<double, 2> a,
  int const& lda,
  arr_ref<double> tau,
  arr_ref<double> work,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DGEQR2 computes a QR factorization of a real m by n matrix A:
  //C  A = Q * R.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix A.  M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix A.  N >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the m by n matrix A.
  //C          On exit, the elements on and above the diagonal of the array
  //C          contain the min(m,n) by n upper trapezoidal matrix R (R is
  //C          upper triangular if m >= n); the elements below the diagonal,
  //C          with the array TAU, represent the orthogonal matrix Q as a
  //C          product of elementary reflectors (see Further Details).
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(1,M).
  //C
  //C  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
  //C          The scalar factors of the elementary reflectors (see Further
  //C          Details).
  //C
  //C  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
  //C
  //C  INFO    (output) INTEGER
  //C          = 0: successful exit
  //C          < 0: if INFO = -i, the i-th argument had an illegal value
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  The matrix Q is represented as a product of elementary reflectors
  //C
  //C     Q = H(1) H(2) . . . H(k), where k = min(m,n).
  //C
  //C  Each H(i) has the form
  //C
  //C     H(i) = I - tau * v * v'
  //C
  //C  where tau is a real scalar, and v is a real vector with
  //C  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
  //C  and tau in TAU(i).
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  if (m < 0) {
    info = -1;
  }
  else if (n < 0) {
    info = -2;
  }
  else if (lda < fem::max(1, m)) {
    info = -4;
  }
  if (info != 0) {
    xerbla("DGEQR2", -info);
    return;
  }
  //C
  int k = fem::min(m, n);
  //C
  int i = fem::int0;
  double aii = fem::double0;
  const double one = 1.0e+0;
  FEM_DO(i, 1, k) {
    //C
    //C        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
    //C
    dlarfp(cmn, m - i + 1, a(i, i), a(fem::min(i + 1, m), i), 1, tau(i));
    if (i < n) {
      //C
      //C           Apply H(i) to A(i:m,i+1:n) from the left
      //C
      aii = a(i, i);
      a(i, i) = one;
      dlarf("Left", m - i + 1, n - i, a(i, i), 1, tau(i), a(i, i + 1),
        lda, work);
      a(i, i) = aii;
    }
  }
  //C
  //C     End of DGEQR2
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dgeqrf.f
inline
void
dgeqrf(
  common& cmn,
  int const& m,
  int const& n,
  arr_ref<double, 2> a,
  int const& lda,
  arr_ref<double> tau,
  arr_ref<double> work,
  int const& lwork,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DGEQRF computes a QR factorization of a real M-by-N matrix A:
  //C  A = Q * R.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix A.  M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix A.  N >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the M-by-N matrix A.
  //C          On exit, the elements on and above the diagonal of the array
  //C          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
  //C          upper triangular if m >= n); the elements below the diagonal,
  //C          with the array TAU, represent the orthogonal matrix Q as a
  //C          product of min(m,n) elementary reflectors (see Further
  //C          Details).
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(1,M).
  //C
  //C  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
  //C          The scalar factors of the elementary reflectors (see Further
  //C          Details).
  //C
  //C  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  //C
  //C  LWORK   (input) INTEGER
  //C          The dimension of the array WORK.  LWORK >= max(1,N).
  //C          For optimum performance LWORK >= N*NB, where NB is
  //C          the optimal blocksize.
  //C
  //C          If LWORK = -1, then a workspace query is assumed; the routine
  //C          only calculates the optimal size of the WORK array, returns
  //C          this value as the first entry of the WORK array, and no error
  //C          message related to LWORK is issued by XERBLA.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  The matrix Q is represented as a product of elementary reflectors
  //C
  //C     Q = H(1) H(2) . . . H(k), where k = min(m,n).
  //C
  //C  Each H(i) has the form
  //C
  //C     H(i) = I - tau * v * v'
  //C
  //C  where tau is a real scalar, and v is a real vector with
  //C  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
  //C  and tau in TAU(i).
  //C
  //C  =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  int nb = ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
  int lwkopt = n * nb;
  work(1) = lwkopt;
  bool lquery = (lwork ==  - 1);
  if (m < 0) {
    info = -1;
  }
  else if (n < 0) {
    info = -2;
  }
  else if (lda < fem::max(1, m)) {
    info = -4;
  }
  else if (lwork < fem::max(1, n) && !lquery) {
    info = -7;
  }
  if (info != 0) {
    xerbla("DGEQRF", -info);
    return;
  }
  else if (lquery) {
    return;
  }
  //C
  //C     Quick return if possible
  //C
  int k = fem::min(m, n);
  if (k == 0) {
    work(1) = 1;
    return;
  }
  //C
  int nbmin = 2;
  int nx = 0;
  int iws = n;
  int ldwork = fem::int0;
  if (nb > 1 && nb < k) {
    //C
    //C        Determine when to cross over from blocked to unblocked code.
    //C
    nx = fem::max(0, ilaenv(3, "DGEQRF", " ", m, n, -1, -1));
    if (nx < k) {
      //C
      //C           Determine if workspace is large enough for blocked code.
      //C
      ldwork = n;
      iws = ldwork * nb;
      if (lwork < iws) {
        //C
        //C              Not enough workspace to use optimal NB:  reduce NB and
        //C              determine the minimum value of NB.
        //C
        nb = lwork / ldwork;
        nbmin = fem::max(2, ilaenv(2, "DGEQRF", " ", m, n, -1, -1));
      }
    }
  }
  //C
  int i = fem::int0;
  int ib = fem::int0;
  int iinfo = fem::int0;
  if (nb >= nbmin && nb < k && nx < k) {
    //C
    //C        Use blocked code initially
    //C
    FEM_DOSTEP(i, 1, k - nx, nb) {
      ib = fem::min(k - i + 1, nb);
      //C
      //C           Compute the QR factorization of the current block
      //C           A(i:m,i:i+ib-1)
      //C
      dgeqr2(cmn, m - i + 1, ib, a(i, i), lda, tau(i), work, iinfo);
      if (i + ib <= n) {
        //C
        //C              Form the triangular factor of the block reflector
        //C              H = H(i) H(i+1) . . . H(i+ib-1)
        //C
        dlarft("Forward", "Columnwise", m - i + 1, ib, a(i, i), lda,
          tau(i), work, ldwork);
        //C
        //C              Apply H' to A(i:m,i+ib:n) from the left
        //C
        dlarfb("Left", "Transpose", "Forward", "Columnwise", m - i + 1,
          n - i - ib + 1, ib, a(i, i), lda, work, ldwork, a(i, i + ib),
          lda, work(ib + 1), ldwork);
      }
    }
  }
  else {
    i = 1;
  }
  //C
  //C     Use unblocked code to factor the last or only block.
  //C
  if (i <= k) {
    dgeqr2(cmn, m - i + 1, n - i + 1, a(i, i), lda, tau(i), work, iinfo);
  }
  //C
  work(1) = iws;
  //C
  //C     End of DGEQRF
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlange.f
inline
double
dlange(
  str_cref norm,
  int const& m,
  int const& n,
  arr_cref<double, 2> a,
  int const& lda,
  arr_ref<double> work)
{
  double return_value = fem::double0;
  a(dimension(lda, star));
  work(dimension(star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLANGE  returns the value of the one norm,  or the Frobenius norm, or
  //C  the  infinity norm,  or the  element of  largest absolute value  of a
  //C  real matrix A.
  //C
  //C  Description
  //C  ===========
  //C
  //C  DLANGE returns the value
  //C
  //C     DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
  //C              (
  //C              ( norm1(A),         NORM = '1', 'O' or 'o'
  //C              (
  //C              ( normI(A),         NORM = 'I' or 'i'
  //C              (
  //C              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
  //C
  //C  where  norm1  denotes the  one norm of a matrix (maximum column sum),
  //C  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
  //C  normF  denotes the  Frobenius norm of a matrix (square root of sum of
  //C  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  NORM    (input) CHARACTER*1
  //C          Specifies the value to be returned in DLANGE as described
  //C          above.
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix A.  M >= 0.  When M = 0,
  //C          DLANGE is set to zero.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix A.  N >= 0.  When N = 0,
  //C          DLANGE is set to zero.
  //C
  //C  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
  //C          The m by n matrix A.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(M,1).
  //C
  //C  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
  //C          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
  //C          referenced.
  //C
  //C =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  const double zero = 0.0e+0;
  double value = fem::double0;
  int j = fem::int0;
  int i = fem::int0;
  double sum = fem::double0;
  double scale = fem::double0;
  const double one = 1.0e+0;
  if (fem::min(m, n) == 0) {
    value = zero;
  }
  else if (lsame(norm, "M")) {
    //C
    //C        Find max(abs(A(i,j))).
    //C
    value = zero;
    FEM_DO(j, 1, n) {
      FEM_DO(i, 1, m) {
        value = fem::max(value, fem::abs(a(i, j)));
      }
    }
  }
  else if ((lsame(norm, "O")) || (norm == "1")) {
    //C
    //C        Find norm1(A).
    //C
    value = zero;
    FEM_DO(j, 1, n) {
      sum = zero;
      FEM_DO(i, 1, m) {
        sum += fem::abs(a(i, j));
      }
      value = fem::max(value, sum);
    }
  }
  else if (lsame(norm, "I")) {
    //C
    //C        Find normI(A).
    //C
    FEM_DO(i, 1, m) {
      work(i) = zero;
    }
    FEM_DO(j, 1, n) {
      FEM_DO(i, 1, m) {
        work(i) += fem::abs(a(i, j));
      }
    }
    value = zero;
    FEM_DO(i, 1, m) {
      value = fem::max(value, work(i));
    }
  }
  else if ((lsame(norm, "F")) || (lsame(norm, "E"))) {
    //C
    //C        Find normF(A).
    //C
    scale = zero;
    sum = one;
    FEM_DO(j, 1, n) {
      dlassq(m, a(1, j), 1, scale, sum);
    }
    value = scale * fem::sqrt(sum);
  }
  //C
  return_value = value;
  return return_value;
  //C
  //C     End of DLANGE
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dorgl2.f
inline
void
dorgl2(
  int const& m,
  int const& n,
  int const& k,
  arr_ref<double, 2> a,
  int const& lda,
  arr_cref<double> tau,
  arr_ref<double> work,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DORGL2 generates an m by n real matrix Q with orthonormal rows,
  //C  which is defined as the first m rows of a product of k elementary
  //C  reflectors of order n
  //C
  //C        Q  =  H(k) . . . H(2) H(1)
  //C
  //C  as returned by DGELQF.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix Q. M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix Q. N >= M.
  //C
  //C  K       (input) INTEGER
  //C          The number of elementary reflectors whose product defines the
  //C          matrix Q. M >= K >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the i-th row must contain the vector which defines
  //C          the elementary reflector H(i), for i = 1,2,...,k, as returned
  //C          by DGELQF in the first k rows of its array argument A.
  //C          On exit, the m-by-n matrix Q.
  //C
  //C  LDA     (input) INTEGER
  //C          The first dimension of the array A. LDA >= max(1,M).
  //C
  //C  TAU     (input) DOUBLE PRECISION array, dimension (K)
  //C          TAU(i) must contain the scalar factor of the elementary
  //C          reflector H(i), as returned by DGELQF.
  //C
  //C  WORK    (workspace) DOUBLE PRECISION array, dimension (M)
  //C
  //C  INFO    (output) INTEGER
  //C          = 0: successful exit
  //C          < 0: if INFO = -i, the i-th argument has an illegal value
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  if (m < 0) {
    info = -1;
  }
  else if (n < m) {
    info = -2;
  }
  else if (k < 0 || k > m) {
    info = -3;
  }
  else if (lda < fem::max(1, m)) {
    info = -5;
  }
  if (info != 0) {
    xerbla("DORGL2", -info);
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (m <= 0) {
    return;
  }
  //C
  int j = fem::int0;
  int l = fem::int0;
  const double zero = 0.0e+0;
  const double one = 1.0e+0;
  if (k < m) {
    //C
    //C        Initialise rows k+1:m to rows of the unit matrix
    //C
    FEM_DO(j, 1, n) {
      FEM_DO(l, k + 1, m) {
        a(l, j) = zero;
      }
      if (j > k && j <= m) {
        a(j, j) = one;
      }
    }
  }
  //C
  int i = fem::int0;
  FEM_DOSTEP(i, k, 1, -1) {
    //C
    //C        Apply H(i) to A(i:m,i:n) from the right
    //C
    if (i < n) {
      if (i < m) {
        a(i, i) = one;
        dlarf("Right", m - i, n - i + 1, a(i, i), lda, tau(i), a(i + 1,
          i), lda, work);
      }
      dscal(n - i, -tau(i), a(i, i + 1), lda);
    }
    a(i, i) = one - tau(i);
    //C
    //C        Set A(i,1:i-1) to zero
    //C
    {
      int fem_do_last = i - 1;
      FEM_DO(l, 1, fem_do_last) {
        a(i, l) = zero;
      }
    }
  }
  //C
  //C     End of DORGL2
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dorglq.f
inline
void
dorglq(
  int const& m,
  int const& n,
  int const& k,
  arr_ref<double, 2> a,
  int const& lda,
  arr_cref<double> tau,
  arr_ref<double> work,
  int const& lwork,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DORGLQ generates an M-by-N real matrix Q with orthonormal rows,
  //C  which is defined as the first M rows of a product of K elementary
  //C  reflectors of order N
  //C
  //C        Q  =  H(k) . . . H(2) H(1)
  //C
  //C  as returned by DGELQF.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix Q. M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix Q. N >= M.
  //C
  //C  K       (input) INTEGER
  //C          The number of elementary reflectors whose product defines the
  //C          matrix Q. M >= K >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the i-th row must contain the vector which defines
  //C          the elementary reflector H(i), for i = 1,2,...,k, as returned
  //C          by DGELQF in the first k rows of its array argument A.
  //C          On exit, the M-by-N matrix Q.
  //C
  //C  LDA     (input) INTEGER
  //C          The first dimension of the array A. LDA >= max(1,M).
  //C
  //C  TAU     (input) DOUBLE PRECISION array, dimension (K)
  //C          TAU(i) must contain the scalar factor of the elementary
  //C          reflector H(i), as returned by DGELQF.
  //C
  //C  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  //C
  //C  LWORK   (input) INTEGER
  //C          The dimension of the array WORK. LWORK >= max(1,M).
  //C          For optimum performance LWORK >= M*NB, where NB is
  //C          the optimal blocksize.
  //C
  //C          If LWORK = -1, then a workspace query is assumed; the routine
  //C          only calculates the optimal size of the WORK array, returns
  //C          this value as the first entry of the WORK array, and no error
  //C          message related to LWORK is issued by XERBLA.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument has an illegal value
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  int nb = ilaenv(1, "DORGLQ", " ", m, n, k, -1);
  int lwkopt = fem::max(1, m) * nb;
  work(1) = lwkopt;
  bool lquery = (lwork ==  - 1);
  if (m < 0) {
    info = -1;
  }
  else if (n < m) {
    info = -2;
  }
  else if (k < 0 || k > m) {
    info = -3;
  }
  else if (lda < fem::max(1, m)) {
    info = -5;
  }
  else if (lwork < fem::max(1, m) && !lquery) {
    info = -8;
  }
  if (info != 0) {
    xerbla("DORGLQ", -info);
    return;
  }
  else if (lquery) {
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (m <= 0) {
    work(1) = 1;
    return;
  }
  //C
  int nbmin = 2;
  int nx = 0;
  int iws = m;
  int ldwork = fem::int0;
  if (nb > 1 && nb < k) {
    //C
    //C        Determine when to cross over from blocked to unblocked code.
    //C
    nx = fem::max(0, ilaenv(3, "DORGLQ", " ", m, n, k, -1));
    if (nx < k) {
      //C
      //C           Determine if workspace is large enough for blocked code.
      //C
      ldwork = m;
      iws = ldwork * nb;
      if (lwork < iws) {
        //C
        //C              Not enough workspace to use optimal NB:  reduce NB and
        //C              determine the minimum value of NB.
        //C
        nb = lwork / ldwork;
        nbmin = fem::max(2, ilaenv(2, "DORGLQ", " ", m, n, k, -1));
      }
    }
  }
  //C
  int ki = fem::int0;
  int kk = fem::int0;
  int j = fem::int0;
  int i = fem::int0;
  const double zero = 0.0e+0;
  if (nb >= nbmin && nb < k && nx < k) {
    //C
    //C        Use blocked code after the last block.
    //C        The first kk rows are handled by the block method.
    //C
    ki = ((k - nx - 1) / nb) * nb;
    kk = fem::min(k, ki + nb);
    //C
    //C        Set A(kk+1:m,1:kk) to zero.
    //C
    FEM_DO(j, 1, kk) {
      FEM_DO(i, kk + 1, m) {
        a(i, j) = zero;
      }
    }
  }
  else {
    kk = 0;
  }
  //C
  //C     Use unblocked code for the last or only block.
  //C
  int iinfo = fem::int0;
  if (kk < m) {
    dorgl2(m - kk, n - kk, k - kk, a(kk + 1, kk + 1), lda, tau(kk + 1),
      work, iinfo);
  }
  //C
  int ib = fem::int0;
  int l = fem::int0;
  if (kk > 0) {
    //C
    //C        Use blocked code
    //C
    FEM_DOSTEP(i, ki + 1, 1, -nb) {
      ib = fem::min(nb, k - i + 1);
      if (i + ib <= m) {
        //C
        //C              Form the triangular factor of the block reflector
        //C              H = H(i) H(i+1) . . . H(i+ib-1)
        //C
        dlarft("Forward", "Rowwise", n - i + 1, ib, a(i, i), lda, tau(i),
          work, ldwork);
        //C
        //C              Apply H' to A(i+ib:m,i:n) from the right
        //C
        dlarfb("Right", "Transpose", "Forward", "Rowwise", m - i - ib + 1,
          n - i + 1, ib, a(i, i), lda, work, ldwork, a(i + ib, i),
          lda, work(ib + 1), ldwork);
      }
      //C
      //C           Apply H' to columns i:n of current block
      //C
      dorgl2(ib, n - i + 1, ib, a(i, i), lda, tau(i), work, iinfo);
      //C
      //C           Set columns 1:i-1 of current block to zero
      //C
      {
        int fem_do_last = i - 1;
        FEM_DO(j, 1, fem_do_last) {
          {
            int fem_do_last = i + ib - 1;
            FEM_DO(l, i, fem_do_last) {
              a(l, j) = zero;
            }
          }
        }
      }
    }
  }
  //C
  work(1) = iws;
  //C
  //C     End of DORGLQ
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dorg2r.f
inline
void
dorg2r(
  int const& m,
  int const& n,
  int const& k,
  arr_ref<double, 2> a,
  int const& lda,
  arr_cref<double> tau,
  arr_ref<double> work,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DORG2R generates an m by n real matrix Q with orthonormal columns,
  //C  which is defined as the first n columns of a product of k elementary
  //C  reflectors of order m
  //C
  //C        Q  =  H(1) H(2) . . . H(k)
  //C
  //C  as returned by DGEQRF.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix Q. M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix Q. M >= N >= 0.
  //C
  //C  K       (input) INTEGER
  //C          The number of elementary reflectors whose product defines the
  //C          matrix Q. N >= K >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the i-th column must contain the vector which
  //C          defines the elementary reflector H(i), for i = 1,2,...,k, as
  //C          returned by DGEQRF in the first k columns of its array
  //C          argument A.
  //C          On exit, the m-by-n matrix Q.
  //C
  //C  LDA     (input) INTEGER
  //C          The first dimension of the array A. LDA >= max(1,M).
  //C
  //C  TAU     (input) DOUBLE PRECISION array, dimension (K)
  //C          TAU(i) must contain the scalar factor of the elementary
  //C          reflector H(i), as returned by DGEQRF.
  //C
  //C  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
  //C
  //C  INFO    (output) INTEGER
  //C          = 0: successful exit
  //C          < 0: if INFO = -i, the i-th argument has an illegal value
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  if (m < 0) {
    info = -1;
  }
  else if (n < 0 || n > m) {
    info = -2;
  }
  else if (k < 0 || k > n) {
    info = -3;
  }
  else if (lda < fem::max(1, m)) {
    info = -5;
  }
  if (info != 0) {
    xerbla("DORG2R", -info);
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (n <= 0) {
    return;
  }
  //C
  //C     Initialise columns k+1:n to columns of the unit matrix
  //C
  int j = fem::int0;
  int l = fem::int0;
  const double zero = 0.0e+0;
  const double one = 1.0e+0;
  FEM_DO(j, k + 1, n) {
    FEM_DO(l, 1, m) {
      a(l, j) = zero;
    }
    a(j, j) = one;
  }
  //C
  int i = fem::int0;
  FEM_DOSTEP(i, k, 1, -1) {
    //C
    //C        Apply H(i) to A(i:m,i:n) from the left
    //C
    if (i < n) {
      a(i, i) = one;
      dlarf("Left", m - i + 1, n - i, a(i, i), 1, tau(i), a(i, i + 1),
        lda, work);
    }
    if (i < m) {
      dscal(m - i, -tau(i), a(i + 1, i), 1);
    }
    a(i, i) = one - tau(i);
    //C
    //C        Set A(1:i-1,i) to zero
    //C
    {
      int fem_do_last = i - 1;
      FEM_DO(l, 1, fem_do_last) {
        a(l, i) = zero;
      }
    }
  }
  //C
  //C     End of DORG2R
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dorgqr.f
inline
void
dorgqr(
  int const& m,
  int const& n,
  int const& k,
  arr_ref<double, 2> a,
  int const& lda,
  arr_cref<double> tau,
  arr_ref<double> work,
  int const& lwork,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DORGQR generates an M-by-N real matrix Q with orthonormal columns,
  //C  which is defined as the first N columns of a product of K elementary
  //C  reflectors of order M
  //C
  //C        Q  =  H(1) H(2) . . . H(k)
  //C
  //C  as returned by DGEQRF.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix Q. M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix Q. M >= N >= 0.
  //C
  //C  K       (input) INTEGER
  //C          The number of elementary reflectors whose product defines the
  //C          matrix Q. N >= K >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the i-th column must contain the vector which
  //C          defines the elementary reflector H(i), for i = 1,2,...,k, as
  //C          returned by DGEQRF in the first k columns of its array
  //C          argument A.
  //C          On exit, the M-by-N matrix Q.
  //C
  //C  LDA     (input) INTEGER
  //C          The first dimension of the array A. LDA >= max(1,M).
  //C
  //C  TAU     (input) DOUBLE PRECISION array, dimension (K)
  //C          TAU(i) must contain the scalar factor of the elementary
  //C          reflector H(i), as returned by DGEQRF.
  //C
  //C  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  //C
  //C  LWORK   (input) INTEGER
  //C          The dimension of the array WORK. LWORK >= max(1,N).
  //C          For optimum performance LWORK >= N*NB, where NB is the
  //C          optimal blocksize.
  //C
  //C          If LWORK = -1, then a workspace query is assumed; the routine
  //C          only calculates the optimal size of the WORK array, returns
  //C          this value as the first entry of the WORK array, and no error
  //C          message related to LWORK is issued by XERBLA.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument has an illegal value
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  int nb = ilaenv(1, "DORGQR", " ", m, n, k, -1);
  int lwkopt = fem::max(1, n) * nb;
  work(1) = lwkopt;
  bool lquery = (lwork ==  - 1);
  if (m < 0) {
    info = -1;
  }
  else if (n < 0 || n > m) {
    info = -2;
  }
  else if (k < 0 || k > n) {
    info = -3;
  }
  else if (lda < fem::max(1, m)) {
    info = -5;
  }
  else if (lwork < fem::max(1, n) && !lquery) {
    info = -8;
  }
  if (info != 0) {
    xerbla("DORGQR", -info);
    return;
  }
  else if (lquery) {
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (n <= 0) {
    work(1) = 1;
    return;
  }
  //C
  int nbmin = 2;
  int nx = 0;
  int iws = n;
  int ldwork = fem::int0;
  if (nb > 1 && nb < k) {
    //C
    //C        Determine when to cross over from blocked to unblocked code.
    //C
    nx = fem::max(0, ilaenv(3, "DORGQR", " ", m, n, k, -1));
    if (nx < k) {
      //C
      //C           Determine if workspace is large enough for blocked code.
      //C
      ldwork = n;
      iws = ldwork * nb;
      if (lwork < iws) {
        //C
        //C              Not enough workspace to use optimal NB:  reduce NB and
        //C              determine the minimum value of NB.
        //C
        nb = lwork / ldwork;
        nbmin = fem::max(2, ilaenv(2, "DORGQR", " ", m, n, k, -1));
      }
    }
  }
  //C
  int ki = fem::int0;
  int kk = fem::int0;
  int j = fem::int0;
  int i = fem::int0;
  const double zero = 0.0e+0;
  if (nb >= nbmin && nb < k && nx < k) {
    //C
    //C        Use blocked code after the last block.
    //C        The first kk columns are handled by the block method.
    //C
    ki = ((k - nx - 1) / nb) * nb;
    kk = fem::min(k, ki + nb);
    //C
    //C        Set A(1:kk,kk+1:n) to zero.
    //C
    FEM_DO(j, kk + 1, n) {
      FEM_DO(i, 1, kk) {
        a(i, j) = zero;
      }
    }
  }
  else {
    kk = 0;
  }
  //C
  //C     Use unblocked code for the last or only block.
  //C
  int iinfo = fem::int0;
  if (kk < n) {
    dorg2r(m - kk, n - kk, k - kk, a(kk + 1, kk + 1), lda, tau(kk + 1),
      work, iinfo);
  }
  //C
  int ib = fem::int0;
  int l = fem::int0;
  if (kk > 0) {
    //C
    //C        Use blocked code
    //C
    FEM_DOSTEP(i, ki + 1, 1, -nb) {
      ib = fem::min(nb, k - i + 1);
      if (i + ib <= n) {
        //C
        //C              Form the triangular factor of the block reflector
        //C              H = H(i) H(i+1) . . . H(i+ib-1)
        //C
        dlarft("Forward", "Columnwise", m - i + 1, ib, a(i, i), lda,
          tau(i), work, ldwork);
        //C
        //C              Apply H to A(i:m,i+ib:n) from the left
        //C
        dlarfb("Left", "No transpose", "Forward", "Columnwise", m - i + 1,
          n - i - ib + 1, ib, a(i, i), lda, work, ldwork, a(i, i + ib),
          lda, work(ib + 1), ldwork);
      }
      //C
      //C           Apply H to rows i:m of current block
      //C
      dorg2r(m - i + 1, ib, ib, a(i, i), lda, tau(i), work, iinfo);
      //C
      //C           Set rows 1:i-1 of current block to zero
      //C
      {
        int fem_do_last = i + ib - 1;
        FEM_DO(j, i, fem_do_last) {
          {
            int fem_do_last = i - 1;
            FEM_DO(l, 1, fem_do_last) {
              a(l, j) = zero;
            }
          }
        }
      }
    }
  }
  //C
  work(1) = iws;
  //C
  //C     End of DORGQR
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dorgbr.f
inline
void
dorgbr(
  str_cref vect,
  int const& m,
  int const& n,
  int const& k,
  arr_ref<double, 2> a,
  int const& lda,
  arr_cref<double> tau,
  arr_ref<double> work,
  int const& lwork,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DORGBR generates one of the real orthogonal matrices Q or P**T
  //C  determined by DGEBRD when reducing a real matrix A to bidiagonal
  //C  form: A = Q * B * P**T.  Q and P**T are defined as products of
  //C  elementary reflectors H(i) or G(i) respectively.
  //C
  //C  If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
  //C  is of order M:
  //C  if m >= k, Q = H(1) H(2) . . . H(k) and DORGBR returns the first n
  //C  columns of Q, where m >= n >= k;
  //C  if m < k, Q = H(1) H(2) . . . H(m-1) and DORGBR returns Q as an
  //C  M-by-M matrix.
  //C
  //C  If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T
  //C  is of order N:
  //C  if k < n, P**T = G(k) . . . G(2) G(1) and DORGBR returns the first m
  //C  rows of P**T, where n >= m >= k;
  //C  if k >= n, P**T = G(n-1) . . . G(2) G(1) and DORGBR returns P**T as
  //C  an N-by-N matrix.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  VECT    (input) CHARACTER*1
  //C          Specifies whether the matrix Q or the matrix P**T is
  //C          required, as defined in the transformation applied by DGEBRD:
  //C          = 'Q':  generate Q;
  //C          = 'P':  generate P**T.
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix Q or P**T to be returned.
  //C          M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix Q or P**T to be returned.
  //C          N >= 0.
  //C          If VECT = 'Q', M >= N >= min(M,K);
  //C          if VECT = 'P', N >= M >= min(N,K).
  //C
  //C  K       (input) INTEGER
  //C          If VECT = 'Q', the number of columns in the original M-by-K
  //C          matrix reduced by DGEBRD.
  //C          If VECT = 'P', the number of rows in the original K-by-N
  //C          matrix reduced by DGEBRD.
  //C          K >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the vectors which define the elementary reflectors,
  //C          as returned by DGEBRD.
  //C          On exit, the M-by-N matrix Q or P**T.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A. LDA >= max(1,M).
  //C
  //C  TAU     (input) DOUBLE PRECISION array, dimension
  //C                                (min(M,K)) if VECT = 'Q'
  //C                                (min(N,K)) if VECT = 'P'
  //C          TAU(i) must contain the scalar factor of the elementary
  //C          reflector H(i) or G(i), which determines Q or P**T, as
  //C          returned by DGEBRD in its array argument TAUQ or TAUP.
  //C
  //C  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  //C
  //C  LWORK   (input) INTEGER
  //C          The dimension of the array WORK. LWORK >= max(1,min(M,N)).
  //C          For optimum performance LWORK >= min(M,N)*NB, where NB
  //C          is the optimal blocksize.
  //C
  //C          If LWORK = -1, then a workspace query is assumed; the routine
  //C          only calculates the optimal size of the WORK array, returns
  //C          this value as the first entry of the WORK array, and no error
  //C          message related to LWORK is issued by XERBLA.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  bool wantq = lsame(vect, "Q");
  int mn = fem::min(m, n);
  bool lquery = (lwork ==  - 1);
  if (!wantq && !lsame(vect, "P")) {
    info = -1;
  }
  else if (m < 0) {
    info = -2;
  }
  else if (n < 0 || (wantq && (n > m || n < fem::min(m, k))) || (
    !wantq && (m > n || m < fem::min(n, k)))) {
    info = -3;
  }
  else if (k < 0) {
    info = -4;
  }
  else if (lda < fem::max(1, m)) {
    info = -6;
  }
  else if (lwork < fem::max(1, mn) && !lquery) {
    info = -9;
  }
  //C
  int nb = fem::int0;
  int lwkopt = fem::int0;
  if (info == 0) {
    if (wantq) {
      nb = ilaenv(1, "DORGQR", " ", m, n, k, -1);
    }
    else {
      nb = ilaenv(1, "DORGLQ", " ", m, n, k, -1);
    }
    lwkopt = fem::max(1, mn) * nb;
    work(1) = lwkopt;
  }
  //C
  if (info != 0) {
    xerbla("DORGBR", -info);
    return;
  }
  else if (lquery) {
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (m == 0 || n == 0) {
    work(1) = 1;
    return;
  }
  //C
  int iinfo = fem::int0;
  int j = fem::int0;
  const double zero = 0.0e+0;
  int i = fem::int0;
  const double one = 1.0e+0;
  if (wantq) {
    //C
    //C        Form Q, determined by a call to DGEBRD to reduce an m-by-k
    //C        matrix
    //C
    if (m >= k) {
      //C
      //C           If m >= k, assume m >= n >= k
      //C
      dorgqr(m, n, k, a, lda, tau, work, lwork, iinfo);
      //C
    }
    else {
      //C
      //C           If m < k, assume m = n
      //C
      //C           Shift the vectors which define the elementary reflectors one
      //C           column to the right, and set the first row and column of Q
      //C           to those of the unit matrix
      //C
      FEM_DOSTEP(j, m, 2, -1) {
        a(1, j) = zero;
        FEM_DO(i, j + 1, m) {
          a(i, j) = a(i, j - 1);
        }
      }
      a(1, 1) = one;
      FEM_DO(i, 2, m) {
        a(i, 1) = zero;
      }
      if (m > 1) {
        //C
        //C              Form Q(2:m,2:m)
        //C
        dorgqr(m - 1, m - 1, m - 1, a(2, 2), lda, tau, work, lwork, iinfo);
      }
    }
  }
  else {
    //C
    //C        Form P', determined by a call to DGEBRD to reduce a k-by-n
    //C        matrix
    //C
    if (k < n) {
      //C
      //C           If k < n, assume k <= m <= n
      //C
      dorglq(m, n, k, a, lda, tau, work, lwork, iinfo);
      //C
    }
    else {
      //C
      //C           If k >= n, assume m = n
      //C
      //C           Shift the vectors which define the elementary reflectors one
      //C           row downward, and set the first row and column of P' to
      //C           those of the unit matrix
      //C
      a(1, 1) = one;
      FEM_DO(i, 2, n) {
        a(i, 1) = zero;
      }
      FEM_DO(j, 2, n) {
        FEM_DOSTEP(i, j - 1, 2, -1) {
          a(i, j) = a(i - 1, j);
        }
        a(1, j) = zero;
      }
      if (n > 1) {
        //C
        //C              Form P'(2:n,2:n)
        //C
        dorglq(n - 1, n - 1, n - 1, a(2, 2), lda, tau, work, lwork, iinfo);
      }
    }
  }
  work(1) = lwkopt;
  //C
  //C     End of DORGBR
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dorml2.f
inline
void
dorml2(
  str_cref side,
  str_cref trans,
  int const& m,
  int const& n,
  int const& k,
  arr_ref<double, 2> a,
  int const& lda,
  arr_cref<double> tau,
  arr_ref<double, 2> c,
  int const& ldc,
  arr_ref<double> work,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  c(dimension(ldc, star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DORML2 overwrites the general real m by n matrix C with
  //C
  //C        Q * C  if SIDE = 'L' and TRANS = 'N', or
  //C
  //C        Q'* C  if SIDE = 'L' and TRANS = 'T', or
  //C
  //C        C * Q  if SIDE = 'R' and TRANS = 'N', or
  //C
  //C        C * Q' if SIDE = 'R' and TRANS = 'T',
  //C
  //C  where Q is a real orthogonal matrix defined as the product of k
  //C  elementary reflectors
  //C
  //C        Q = H(k) . . . H(2) H(1)
  //C
  //C  as returned by DGELQF. Q is of order m if SIDE = 'L' and of order n
  //C  if SIDE = 'R'.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  SIDE    (input) CHARACTER*1
  //C          = 'L': apply Q or Q' from the Left
  //C          = 'R': apply Q or Q' from the Right
  //C
  //C  TRANS   (input) CHARACTER*1
  //C          = 'N': apply Q  (No transpose)
  //C          = 'T': apply Q' (Transpose)
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix C. M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix C. N >= 0.
  //C
  //C  K       (input) INTEGER
  //C          The number of elementary reflectors whose product defines
  //C          the matrix Q.
  //C          If SIDE = 'L', M >= K >= 0;
  //C          if SIDE = 'R', N >= K >= 0.
  //C
  //C  A       (input) DOUBLE PRECISION array, dimension
  //C                               (LDA,M) if SIDE = 'L',
  //C                               (LDA,N) if SIDE = 'R'
  //C          The i-th row must contain the vector which defines the
  //C          elementary reflector H(i), for i = 1,2,...,k, as returned by
  //C          DGELQF in the first k rows of its array argument A.
  //C          A is modified by the routine but restored on exit.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A. LDA >= max(1,K).
  //C
  //C  TAU     (input) DOUBLE PRECISION array, dimension (K)
  //C          TAU(i) must contain the scalar factor of the elementary
  //C          reflector H(i), as returned by DGELQF.
  //C
  //C  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
  //C          On entry, the m by n matrix C.
  //C          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
  //C
  //C  LDC     (input) INTEGER
  //C          The leading dimension of the array C. LDC >= max(1,M).
  //C
  //C  WORK    (workspace) DOUBLE PRECISION array, dimension
  //C                                   (N) if SIDE = 'L',
  //C                                   (M) if SIDE = 'R'
  //C
  //C  INFO    (output) INTEGER
  //C          = 0: successful exit
  //C          < 0: if INFO = -i, the i-th argument had an illegal value
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  bool left = lsame(side, "L");
  bool notran = lsame(trans, "N");
  //C
  //C     NQ is the order of Q
  //C
  int nq = fem::int0;
  if (left) {
    nq = m;
  }
  else {
    nq = n;
  }
  if (!left && !lsame(side, "R")) {
    info = -1;
  }
  else if (!notran && !lsame(trans, "T")) {
    info = -2;
  }
  else if (m < 0) {
    info = -3;
  }
  else if (n < 0) {
    info = -4;
  }
  else if (k < 0 || k > nq) {
    info = -5;
  }
  else if (lda < fem::max(1, k)) {
    info = -7;
  }
  else if (ldc < fem::max(1, m)) {
    info = -10;
  }
  if (info != 0) {
    xerbla("DORML2", -info);
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (m == 0 || n == 0 || k == 0) {
    return;
  }
  //C
  int i1 = fem::int0;
  int i2 = fem::int0;
  int i3 = fem::int0;
  if ((left && notran) || (!left && !notran)) {
    i1 = 1;
    i2 = k;
    i3 = 1;
  }
  else {
    i1 = k;
    i2 = 1;
    i3 = -1;
  }
  //C
  int ni = fem::int0;
  int jc = fem::int0;
  int mi = fem::int0;
  int ic = fem::int0;
  if (left) {
    ni = n;
    jc = 1;
  }
  else {
    mi = m;
    ic = 1;
  }
  //C
  int i = fem::int0;
  double aii = fem::double0;
  const double one = 1.0e+0;
  FEM_DOSTEP(i, i1, i2, i3) {
    if (left) {
      //C
      //C           H(i) is applied to C(i:m,1:n)
      //C
      mi = m - i + 1;
      ic = i;
    }
    else {
      //C
      //C           H(i) is applied to C(1:m,i:n)
      //C
      ni = n - i + 1;
      jc = i;
    }
    //C
    //C        Apply H(i)
    //C
    aii = a(i, i);
    a(i, i) = one;
    dlarf(side, mi, ni, a(i, i), lda, tau(i), c(ic, jc), ldc, work);
    a(i, i) = aii;
  }
  //C
  //C     End of DORML2
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dormlq.f
inline
void
dormlq(
  str_cref side,
  str_cref trans,
  int const& m,
  int const& n,
  int const& k,
  arr_ref<double, 2> a,
  int const& lda,
  arr_cref<double> tau,
  arr_ref<double, 2> c,
  int const& ldc,
  arr_ref<double> work,
  int const& lwork,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  c(dimension(ldc, star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DORMLQ overwrites the general real M-by-N matrix C with
  //C
  //C                  SIDE = 'L'     SIDE = 'R'
  //C  TRANS = 'N':      Q * C          C * Q
  //C  TRANS = 'T':      Q**T * C       C * Q**T
  //C
  //C  where Q is a real orthogonal matrix defined as the product of k
  //C  elementary reflectors
  //C
  //C        Q = H(k) . . . H(2) H(1)
  //C
  //C  as returned by DGELQF. Q is of order M if SIDE = 'L' and of order N
  //C  if SIDE = 'R'.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  SIDE    (input) CHARACTER*1
  //C          = 'L': apply Q or Q**T from the Left;
  //C          = 'R': apply Q or Q**T from the Right.
  //C
  //C  TRANS   (input) CHARACTER*1
  //C          = 'N':  No transpose, apply Q;
  //C          = 'T':  Transpose, apply Q**T.
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix C. M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix C. N >= 0.
  //C
  //C  K       (input) INTEGER
  //C          The number of elementary reflectors whose product defines
  //C          the matrix Q.
  //C          If SIDE = 'L', M >= K >= 0;
  //C          if SIDE = 'R', N >= K >= 0.
  //C
  //C  A       (input) DOUBLE PRECISION array, dimension
  //C                               (LDA,M) if SIDE = 'L',
  //C                               (LDA,N) if SIDE = 'R'
  //C          The i-th row must contain the vector which defines the
  //C          elementary reflector H(i), for i = 1,2,...,k, as returned by
  //C          DGELQF in the first k rows of its array argument A.
  //C          A is modified by the routine but restored on exit.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A. LDA >= max(1,K).
  //C
  //C  TAU     (input) DOUBLE PRECISION array, dimension (K)
  //C          TAU(i) must contain the scalar factor of the elementary
  //C          reflector H(i), as returned by DGELQF.
  //C
  //C  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
  //C          On entry, the M-by-N matrix C.
  //C          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
  //C
  //C  LDC     (input) INTEGER
  //C          The leading dimension of the array C. LDC >= max(1,M).
  //C
  //C  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  //C
  //C  LWORK   (input) INTEGER
  //C          The dimension of the array WORK.
  //C          If SIDE = 'L', LWORK >= max(1,N);
  //C          if SIDE = 'R', LWORK >= max(1,M).
  //C          For optimum performance LWORK >= N*NB if SIDE = 'L', and
  //C          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
  //C          blocksize.
  //C
  //C          If LWORK = -1, then a workspace query is assumed; the routine
  //C          only calculates the optimal size of the WORK array, returns
  //C          this value as the first entry of the WORK array, and no error
  //C          message related to LWORK is issued by XERBLA.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Local Arrays ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  bool left = lsame(side, "L");
  bool notran = lsame(trans, "N");
  bool lquery = (lwork ==  - 1);
  //C
  //C     NQ is the order of Q and NW is the minimum dimension of WORK
  //C
  int nq = fem::int0;
  int nw = fem::int0;
  if (left) {
    nq = m;
    nw = n;
  }
  else {
    nq = n;
    nw = m;
  }
  if (!left && !lsame(side, "R")) {
    info = -1;
  }
  else if (!notran && !lsame(trans, "T")) {
    info = -2;
  }
  else if (m < 0) {
    info = -3;
  }
  else if (n < 0) {
    info = -4;
  }
  else if (k < 0 || k > nq) {
    info = -5;
  }
  else if (lda < fem::max(1, k)) {
    info = -7;
  }
  else if (ldc < fem::max(1, m)) {
    info = -10;
  }
  else if (lwork < fem::max(1, nw) && !lquery) {
    info = -12;
  }
  //C
  const int nbmax = 64;
  int nb = fem::int0;
  int lwkopt = fem::int0;
  if (info == 0) {
    //C
    //C        Determine the block size.  NB may be at most NBMAX, where NBMAX
    //C        is used to define the local array T.
    //C
    nb = fem::min(nbmax, ilaenv(1, "DORMLQ", side + trans, m, n, k, -1));
    lwkopt = fem::max(1, nw) * nb;
    work(1) = lwkopt;
  }
  //C
  if (info != 0) {
    xerbla("DORMLQ", -info);
    return;
  }
  else if (lquery) {
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (m == 0 || n == 0 || k == 0) {
    work(1) = 1;
    return;
  }
  //C
  int nbmin = 2;
  int ldwork = nw;
  int iws = fem::int0;
  if (nb > 1 && nb < k) {
    iws = nw * nb;
    if (lwork < iws) {
      nb = lwork / ldwork;
      nbmin = fem::max(2, ilaenv(2, "DORMLQ", side + trans, m, n, k, -1));
    }
  }
  else {
    iws = nw;
  }
  //C
  int iinfo = fem::int0;
  int i1 = fem::int0;
  int i2 = fem::int0;
  int i3 = fem::int0;
  int ni = fem::int0;
  int jc = fem::int0;
  int mi = fem::int0;
  int ic = fem::int0;
  fem::str<1> transt = fem::char0;
  int i = fem::int0;
  int ib = fem::int0;
  const int ldt = nbmax + 1;
  arr<double, 2> t(dimension(ldt, nbmax), fem::fill0);
  if (nb < nbmin || nb >= k) {
    //C
    //C        Use unblocked code
    //C
    dorml2(side, trans, m, n, k, a, lda, tau, c, ldc, work, iinfo);
  }
  else {
    //C
    //C        Use blocked code
    //C
    if ((left && notran) || (!left && !notran)) {
      i1 = 1;
      i2 = k;
      i3 = nb;
    }
    else {
      i1 = ((k - 1) / nb) * nb + 1;
      i2 = 1;
      i3 = -nb;
    }
    //C
    if (left) {
      ni = n;
      jc = 1;
    }
    else {
      mi = m;
      ic = 1;
    }
    //C
    if (notran) {
      transt = "T";
    }
    else {
      transt = "N";
    }
    //C
    FEM_DOSTEP(i, i1, i2, i3) {
      ib = fem::min(nb, k - i + 1);
      //C
      //C           Form the triangular factor of the block reflector
      //C           H = H(i) H(i+1) . . . H(i+ib-1)
      //C
      dlarft("Forward", "Rowwise", nq - i + 1, ib, a(i, i), lda, tau(i),
        t, ldt);
      if (left) {
        //C
        //C              H or H' is applied to C(i:m,1:n)
        //C
        mi = m - i + 1;
        ic = i;
      }
      else {
        //C
        //C              H or H' is applied to C(1:m,i:n)
        //C
        ni = n - i + 1;
        jc = i;
      }
      //C
      //C           Apply H or H'
      //C
      dlarfb(side, transt, "Forward", "Rowwise", mi, ni, ib, a(i, i),
        lda, t, ldt, c(ic, jc), ldc, work, ldwork);
    }
  }
  work(1) = lwkopt;
  //C
  //C     End of DORMLQ
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dorm2r.f
inline
void
dorm2r(
  str_cref side,
  str_cref trans,
  int const& m,
  int const& n,
  int const& k,
  arr_ref<double, 2> a,
  int const& lda,
  arr_cref<double> tau,
  arr_ref<double, 2> c,
  int const& ldc,
  arr_ref<double> work,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  c(dimension(ldc, star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DORM2R overwrites the general real m by n matrix C with
  //C
  //C        Q * C  if SIDE = 'L' and TRANS = 'N', or
  //C
  //C        Q'* C  if SIDE = 'L' and TRANS = 'T', or
  //C
  //C        C * Q  if SIDE = 'R' and TRANS = 'N', or
  //C
  //C        C * Q' if SIDE = 'R' and TRANS = 'T',
  //C
  //C  where Q is a real orthogonal matrix defined as the product of k
  //C  elementary reflectors
  //C
  //C        Q = H(1) H(2) . . . H(k)
  //C
  //C  as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
  //C  if SIDE = 'R'.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  SIDE    (input) CHARACTER*1
  //C          = 'L': apply Q or Q' from the Left
  //C          = 'R': apply Q or Q' from the Right
  //C
  //C  TRANS   (input) CHARACTER*1
  //C          = 'N': apply Q  (No transpose)
  //C          = 'T': apply Q' (Transpose)
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix C. M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix C. N >= 0.
  //C
  //C  K       (input) INTEGER
  //C          The number of elementary reflectors whose product defines
  //C          the matrix Q.
  //C          If SIDE = 'L', M >= K >= 0;
  //C          if SIDE = 'R', N >= K >= 0.
  //C
  //C  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
  //C          The i-th column must contain the vector which defines the
  //C          elementary reflector H(i), for i = 1,2,...,k, as returned by
  //C          DGEQRF in the first k columns of its array argument A.
  //C          A is modified by the routine but restored on exit.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.
  //C          If SIDE = 'L', LDA >= max(1,M);
  //C          if SIDE = 'R', LDA >= max(1,N).
  //C
  //C  TAU     (input) DOUBLE PRECISION array, dimension (K)
  //C          TAU(i) must contain the scalar factor of the elementary
  //C          reflector H(i), as returned by DGEQRF.
  //C
  //C  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
  //C          On entry, the m by n matrix C.
  //C          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
  //C
  //C  LDC     (input) INTEGER
  //C          The leading dimension of the array C. LDC >= max(1,M).
  //C
  //C  WORK    (workspace) DOUBLE PRECISION array, dimension
  //C                                   (N) if SIDE = 'L',
  //C                                   (M) if SIDE = 'R'
  //C
  //C  INFO    (output) INTEGER
  //C          = 0: successful exit
  //C          < 0: if INFO = -i, the i-th argument had an illegal value
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  bool left = lsame(side, "L");
  bool notran = lsame(trans, "N");
  //C
  //C     NQ is the order of Q
  //C
  int nq = fem::int0;
  if (left) {
    nq = m;
  }
  else {
    nq = n;
  }
  if (!left && !lsame(side, "R")) {
    info = -1;
  }
  else if (!notran && !lsame(trans, "T")) {
    info = -2;
  }
  else if (m < 0) {
    info = -3;
  }
  else if (n < 0) {
    info = -4;
  }
  else if (k < 0 || k > nq) {
    info = -5;
  }
  else if (lda < fem::max(1, nq)) {
    info = -7;
  }
  else if (ldc < fem::max(1, m)) {
    info = -10;
  }
  if (info != 0) {
    xerbla("DORM2R", -info);
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (m == 0 || n == 0 || k == 0) {
    return;
  }
  //C
  int i1 = fem::int0;
  int i2 = fem::int0;
  int i3 = fem::int0;
  if ((left && !notran) || (!left && notran)) {
    i1 = 1;
    i2 = k;
    i3 = 1;
  }
  else {
    i1 = k;
    i2 = 1;
    i3 = -1;
  }
  //C
  int ni = fem::int0;
  int jc = fem::int0;
  int mi = fem::int0;
  int ic = fem::int0;
  if (left) {
    ni = n;
    jc = 1;
  }
  else {
    mi = m;
    ic = 1;
  }
  //C
  int i = fem::int0;
  double aii = fem::double0;
  const double one = 1.0e+0;
  FEM_DOSTEP(i, i1, i2, i3) {
    if (left) {
      //C
      //C           H(i) is applied to C(i:m,1:n)
      //C
      mi = m - i + 1;
      ic = i;
    }
    else {
      //C
      //C           H(i) is applied to C(1:m,i:n)
      //C
      ni = n - i + 1;
      jc = i;
    }
    //C
    //C        Apply H(i)
    //C
    aii = a(i, i);
    a(i, i) = one;
    dlarf(side, mi, ni, a(i, i), 1, tau(i), c(ic, jc), ldc, work);
    a(i, i) = aii;
  }
  //C
  //C     End of DORM2R
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dormqr.f
inline
void
dormqr(
  str_cref side,
  str_cref trans,
  int const& m,
  int const& n,
  int const& k,
  arr_ref<double, 2> a,
  int const& lda,
  arr_cref<double> tau,
  arr_ref<double, 2> c,
  int const& ldc,
  arr_ref<double> work,
  int const& lwork,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  c(dimension(ldc, star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DORMQR overwrites the general real M-by-N matrix C with
  //C
  //C                  SIDE = 'L'     SIDE = 'R'
  //C  TRANS = 'N':      Q * C          C * Q
  //C  TRANS = 'T':      Q**T * C       C * Q**T
  //C
  //C  where Q is a real orthogonal matrix defined as the product of k
  //C  elementary reflectors
  //C
  //C        Q = H(1) H(2) . . . H(k)
  //C
  //C  as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
  //C  if SIDE = 'R'.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  SIDE    (input) CHARACTER*1
  //C          = 'L': apply Q or Q**T from the Left;
  //C          = 'R': apply Q or Q**T from the Right.
  //C
  //C  TRANS   (input) CHARACTER*1
  //C          = 'N':  No transpose, apply Q;
  //C          = 'T':  Transpose, apply Q**T.
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix C. M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix C. N >= 0.
  //C
  //C  K       (input) INTEGER
  //C          The number of elementary reflectors whose product defines
  //C          the matrix Q.
  //C          If SIDE = 'L', M >= K >= 0;
  //C          if SIDE = 'R', N >= K >= 0.
  //C
  //C  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
  //C          The i-th column must contain the vector which defines the
  //C          elementary reflector H(i), for i = 1,2,...,k, as returned by
  //C          DGEQRF in the first k columns of its array argument A.
  //C          A is modified by the routine but restored on exit.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.
  //C          If SIDE = 'L', LDA >= max(1,M);
  //C          if SIDE = 'R', LDA >= max(1,N).
  //C
  //C  TAU     (input) DOUBLE PRECISION array, dimension (K)
  //C          TAU(i) must contain the scalar factor of the elementary
  //C          reflector H(i), as returned by DGEQRF.
  //C
  //C  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
  //C          On entry, the M-by-N matrix C.
  //C          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
  //C
  //C  LDC     (input) INTEGER
  //C          The leading dimension of the array C. LDC >= max(1,M).
  //C
  //C  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  //C
  //C  LWORK   (input) INTEGER
  //C          The dimension of the array WORK.
  //C          If SIDE = 'L', LWORK >= max(1,N);
  //C          if SIDE = 'R', LWORK >= max(1,M).
  //C          For optimum performance LWORK >= N*NB if SIDE = 'L', and
  //C          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
  //C          blocksize.
  //C
  //C          If LWORK = -1, then a workspace query is assumed; the routine
  //C          only calculates the optimal size of the WORK array, returns
  //C          this value as the first entry of the WORK array, and no error
  //C          message related to LWORK is issued by XERBLA.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Local Arrays ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  bool left = lsame(side, "L");
  bool notran = lsame(trans, "N");
  bool lquery = (lwork ==  - 1);
  //C
  //C     NQ is the order of Q and NW is the minimum dimension of WORK
  //C
  int nq = fem::int0;
  int nw = fem::int0;
  if (left) {
    nq = m;
    nw = n;
  }
  else {
    nq = n;
    nw = m;
  }
  if (!left && !lsame(side, "R")) {
    info = -1;
  }
  else if (!notran && !lsame(trans, "T")) {
    info = -2;
  }
  else if (m < 0) {
    info = -3;
  }
  else if (n < 0) {
    info = -4;
  }
  else if (k < 0 || k > nq) {
    info = -5;
  }
  else if (lda < fem::max(1, nq)) {
    info = -7;
  }
  else if (ldc < fem::max(1, m)) {
    info = -10;
  }
  else if (lwork < fem::max(1, nw) && !lquery) {
    info = -12;
  }
  //C
  const int nbmax = 64;
  int nb = fem::int0;
  int lwkopt = fem::int0;
  if (info == 0) {
    //C
    //C        Determine the block size.  NB may be at most NBMAX, where NBMAX
    //C        is used to define the local array T.
    //C
    nb = fem::min(nbmax, ilaenv(1, "DORMQR", side + trans, m, n, k, -1));
    lwkopt = fem::max(1, nw) * nb;
    work(1) = lwkopt;
  }
  //C
  if (info != 0) {
    xerbla("DORMQR", -info);
    return;
  }
  else if (lquery) {
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (m == 0 || n == 0 || k == 0) {
    work(1) = 1;
    return;
  }
  //C
  int nbmin = 2;
  int ldwork = nw;
  int iws = fem::int0;
  if (nb > 1 && nb < k) {
    iws = nw * nb;
    if (lwork < iws) {
      nb = lwork / ldwork;
      nbmin = fem::max(2, ilaenv(2, "DORMQR", side + trans, m, n, k, -1));
    }
  }
  else {
    iws = nw;
  }
  //C
  int iinfo = fem::int0;
  int i1 = fem::int0;
  int i2 = fem::int0;
  int i3 = fem::int0;
  int ni = fem::int0;
  int jc = fem::int0;
  int mi = fem::int0;
  int ic = fem::int0;
  int i = fem::int0;
  int ib = fem::int0;
  const int ldt = nbmax + 1;
  arr<double, 2> t(dimension(ldt, nbmax), fem::fill0);
  if (nb < nbmin || nb >= k) {
    //C
    //C        Use unblocked code
    //C
    dorm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, iinfo);
  }
  else {
    //C
    //C        Use blocked code
    //C
    if ((left && !notran) || (!left && notran)) {
      i1 = 1;
      i2 = k;
      i3 = nb;
    }
    else {
      i1 = ((k - 1) / nb) * nb + 1;
      i2 = 1;
      i3 = -nb;
    }
    //C
    if (left) {
      ni = n;
      jc = 1;
    }
    else {
      mi = m;
      ic = 1;
    }
    //C
    FEM_DOSTEP(i, i1, i2, i3) {
      ib = fem::min(nb, k - i + 1);
      //C
      //C           Form the triangular factor of the block reflector
      //C           H = H(i) H(i+1) . . . H(i+ib-1)
      //C
      dlarft("Forward", "Columnwise", nq - i + 1, ib, a(i, i), lda,
        tau(i), t, ldt);
      if (left) {
        //C
        //C              H or H' is applied to C(i:m,1:n)
        //C
        mi = m - i + 1;
        ic = i;
      }
      else {
        //C
        //C              H or H' is applied to C(1:m,i:n)
        //C
        ni = n - i + 1;
        jc = i;
      }
      //C
      //C           Apply H or H'
      //C
      dlarfb(side, trans, "Forward", "Columnwise", mi, ni, ib, a(i,
        i), lda, t, ldt, c(ic, jc), ldc, work, ldwork);
    }
  }
  work(1) = lwkopt;
  //C
  //C     End of DORMQR
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dormbr.f
inline
void
dormbr(
  str_cref vect,
  str_cref side,
  str_cref trans,
  int const& m,
  int const& n,
  int const& k,
  arr_ref<double, 2> a,
  int const& lda,
  arr_cref<double> tau,
  arr_ref<double, 2> c,
  int const& ldc,
  arr_ref<double> work,
  int const& lwork,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  c(dimension(ldc, star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  If VECT = 'Q', DORMBR overwrites the general real M-by-N matrix C
  //C  with
  //C                  SIDE = 'L'     SIDE = 'R'
  //C  TRANS = 'N':      Q * C          C * Q
  //C  TRANS = 'T':      Q**T * C       C * Q**T
  //C
  //C  If VECT = 'P', DORMBR overwrites the general real M-by-N matrix C
  //C  with
  //C                  SIDE = 'L'     SIDE = 'R'
  //C  TRANS = 'N':      P * C          C * P
  //C  TRANS = 'T':      P**T * C       C * P**T
  //C
  //C  Here Q and P**T are the orthogonal matrices determined by DGEBRD when
  //C  reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and
  //C  P**T are defined as products of elementary reflectors H(i) and G(i)
  //C  respectively.
  //C
  //C  Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
  //C  order of the orthogonal matrix Q or P**T that is applied.
  //C
  //C  If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
  //C  if nq >= k, Q = H(1) H(2) . . . H(k);
  //C  if nq < k, Q = H(1) H(2) . . . H(nq-1).
  //C
  //C  If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
  //C  if k < nq, P = G(1) G(2) . . . G(k);
  //C  if k >= nq, P = G(1) G(2) . . . G(nq-1).
  //C
  //C  Arguments
  //C  =========
  //C
  //C  VECT    (input) CHARACTER*1
  //C          = 'Q': apply Q or Q**T;
  //C          = 'P': apply P or P**T.
  //C
  //C  SIDE    (input) CHARACTER*1
  //C          = 'L': apply Q, Q**T, P or P**T from the Left;
  //C          = 'R': apply Q, Q**T, P or P**T from the Right.
  //C
  //C  TRANS   (input) CHARACTER*1
  //C          = 'N':  No transpose, apply Q  or P;
  //C          = 'T':  Transpose, apply Q**T or P**T.
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix C. M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix C. N >= 0.
  //C
  //C  K       (input) INTEGER
  //C          If VECT = 'Q', the number of columns in the original
  //C          matrix reduced by DGEBRD.
  //C          If VECT = 'P', the number of rows in the original
  //C          matrix reduced by DGEBRD.
  //C          K >= 0.
  //C
  //C  A       (input) DOUBLE PRECISION array, dimension
  //C                                (LDA,min(nq,K)) if VECT = 'Q'
  //C                                (LDA,nq)        if VECT = 'P'
  //C          The vectors which define the elementary reflectors H(i) and
  //C          G(i), whose products determine the matrices Q and P, as
  //C          returned by DGEBRD.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.
  //C          If VECT = 'Q', LDA >= max(1,nq);
  //C          if VECT = 'P', LDA >= max(1,min(nq,K)).
  //C
  //C  TAU     (input) DOUBLE PRECISION array, dimension (min(nq,K))
  //C          TAU(i) must contain the scalar factor of the elementary
  //C          reflector H(i) or G(i) which determines Q or P, as returned
  //C          by DGEBRD in the array argument TAUQ or TAUP.
  //C
  //C  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
  //C          On entry, the M-by-N matrix C.
  //C          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q
  //C          or P*C or P**T*C or C*P or C*P**T.
  //C
  //C  LDC     (input) INTEGER
  //C          The leading dimension of the array C. LDC >= max(1,M).
  //C
  //C  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  //C
  //C  LWORK   (input) INTEGER
  //C          The dimension of the array WORK.
  //C          If SIDE = 'L', LWORK >= max(1,N);
  //C          if SIDE = 'R', LWORK >= max(1,M).
  //C          For optimum performance LWORK >= N*NB if SIDE = 'L', and
  //C          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
  //C          blocksize.
  //C
  //C          If LWORK = -1, then a workspace query is assumed; the routine
  //C          only calculates the optimal size of the WORK array, returns
  //C          this value as the first entry of the WORK array, and no error
  //C          message related to LWORK is issued by XERBLA.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value
  //C
  //C  =====================================================================
  //C
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  bool applyq = lsame(vect, "Q");
  bool left = lsame(side, "L");
  bool notran = lsame(trans, "N");
  bool lquery = (lwork ==  - 1);
  //C
  //C     NQ is the order of Q or P and NW is the minimum dimension of WORK
  //C
  int nq = fem::int0;
  int nw = fem::int0;
  if (left) {
    nq = m;
    nw = n;
  }
  else {
    nq = n;
    nw = m;
  }
  if (!applyq && !lsame(vect, "P")) {
    info = -1;
  }
  else if (!left && !lsame(side, "R")) {
    info = -2;
  }
  else if (!notran && !lsame(trans, "T")) {
    info = -3;
  }
  else if (m < 0) {
    info = -4;
  }
  else if (n < 0) {
    info = -5;
  }
  else if (k < 0) {
    info = -6;
  }
  else if ((applyq && lda < fem::max(1, nq)) || (!applyq && lda < fem::max(1,
    fem::min(nq, k)))) {
    info = -8;
  }
  else if (ldc < fem::max(1, m)) {
    info = -11;
  }
  else if (lwork < fem::max(1, nw) && !lquery) {
    info = -13;
  }
  //C
  int nb = fem::int0;
  int lwkopt = fem::int0;
  if (info == 0) {
    if (applyq) {
      if (left) {
        nb = ilaenv(1, "DORMQR", side + trans, m - 1, n, m - 1, -1);
      }
      else {
        nb = ilaenv(1, "DORMQR", side + trans, m, n - 1, n - 1, -1);
      }
    }
    else {
      if (left) {
        nb = ilaenv(1, "DORMLQ", side + trans, m - 1, n, m - 1, -1);
      }
      else {
        nb = ilaenv(1, "DORMLQ", side + trans, m, n - 1, n - 1, -1);
      }
    }
    lwkopt = fem::max(1, nw) * nb;
    work(1) = lwkopt;
  }
  //C
  if (info != 0) {
    xerbla("DORMBR", -info);
    return;
  }
  else if (lquery) {
    return;
  }
  //C
  //C     Quick return if possible
  //C
  work(1) = 1;
  if (m == 0 || n == 0) {
    return;
  }
  //C
  int iinfo = fem::int0;
  int mi = fem::int0;
  int ni = fem::int0;
  int i1 = fem::int0;
  int i2 = fem::int0;
  fem::str<1> transt = fem::char0;
  if (applyq) {
    //C
    //C        Apply Q
    //C
    if (nq >= k) {
      //C
      //C           Q was determined by a call to DGEBRD with nq >= k
      //C
      dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, iinfo);
    }
    else if (nq > 1) {
      //C
      //C           Q was determined by a call to DGEBRD with nq < k
      //C
      if (left) {
        mi = m - 1;
        ni = n;
        i1 = 2;
        i2 = 1;
      }
      else {
        mi = m;
        ni = n - 1;
        i1 = 1;
        i2 = 2;
      }
      dormqr(side, trans, mi, ni, nq - 1, a(2, 1), lda, tau, c(i1,
        i2), ldc, work, lwork, iinfo);
    }
  }
  else {
    //C
    //C        Apply P
    //C
    if (notran) {
      transt = "T";
    }
    else {
      transt = "N";
    }
    if (nq > k) {
      //C
      //C           P was determined by a call to DGEBRD with nq > k
      //C
      dormlq(side, transt, m, n, k, a, lda, tau, c, ldc, work, lwork, iinfo);
    }
    else if (nq > 1) {
      //C
      //C           P was determined by a call to DGEBRD with nq <= k
      //C
      if (left) {
        mi = m - 1;
        ni = n;
        i1 = 2;
        i2 = 1;
      }
      else {
        mi = m;
        ni = n - 1;
        i1 = 1;
        i2 = 2;
      }
      dormlq(side, transt, mi, ni, nq - 1, a(1, 2), lda, tau, c(i1,
        i2), ldc, work, lwork, iinfo);
    }
  }
  work(1) = lwkopt;
  //C
  //C     End of DORMBR
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dgesdd.f
inline
void
dgesdd(
  common& cmn,
  str_cref jobz,
  int const& m,
  int const& n,
  arr_ref<double, 2> a,
  int const& lda,
  arr_ref<double> s,
  arr_ref<double, 2> u,
  int const& ldu,
  arr_ref<double, 2> vt,
  int const& ldvt,
  arr_ref<double> work,
  int const& lwork,
  arr_ref<int> iwork,
  int& info)
{
  a(dimension(lda, star));
  s(dimension(star));
  u(dimension(ldu, star));
  vt(dimension(ldvt, star));
  work(dimension(star));
  iwork(dimension(star));
  //C
  //C  -- LAPACK driver routine (version 3.2.1)                                  --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     March 2009
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DGESDD computes the singular value decomposition (SVD) of a real
  //C  M-by-N matrix A, optionally computing the left and right singular
  //C  vectors.  If singular vectors are desired, it uses a
  //C  divide-and-conquer algorithm.
  //C
  //C  The SVD is written
  //C
  //C       A = U * SIGMA * transpose(V)
  //C
  //C  where SIGMA is an M-by-N matrix which is zero except for its
  //C  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
  //C  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
  //C  are the singular values of A; they are real and non-negative, and
  //C  are returned in descending order.  The first min(m,n) columns of
  //C  U and V are the left and right singular vectors of A.
  //C
  //C  Note that the routine returns VT = V**T, not V.
  //C
  //C  The divide and conquer algorithm makes very mild assumptions about
  //C  floating point arithmetic. It will work on machines with a guard
  //C  digit in add/subtract, or on those binary machines without guard
  //C  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
  //C  Cray-2. It could conceivably fail on hexadecimal or decimal machines
  //C  without guard digits, but we know of none.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  JOBZ    (input) CHARACTER*1
  //C          Specifies options for computing all or part of the matrix U:
  //C          = 'A':  all M columns of U and all N rows of V**T are
  //C                  returned in the arrays U and VT;
  //C          = 'S':  the first min(M,N) columns of U and the first
  //C                  min(M,N) rows of V**T are returned in the arrays U
  //C                  and VT;
  //C          = 'O':  If M >= N, the first N columns of U are overwritten
  //C                  on the array A and all rows of V**T are returned in
  //C                  the array VT;
  //C                  otherwise, all columns of U are returned in the
  //C                  array U and the first M rows of V**T are overwritten
  //C                  in the array A;
  //C          = 'N':  no columns of U or rows of V**T are computed.
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the input matrix A.  M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the input matrix A.  N >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the M-by-N matrix A.
  //C          On exit,
  //C          if JOBZ = 'O',  A is overwritten with the first N columns
  //C                          of U (the left singular vectors, stored
  //C                          columnwise) if M >= N;
  //C                          A is overwritten with the first M rows
  //C                          of V**T (the right singular vectors, stored
  //C                          rowwise) otherwise.
  //C          if JOBZ .ne. 'O', the contents of A are destroyed.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(1,M).
  //C
  //C  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
  //C          The singular values of A, sorted so that S(i) >= S(i+1).
  //C
  //C  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
  //C          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;
  //C          UCOL = min(M,N) if JOBZ = 'S'.
  //C          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M
  //C          orthogonal matrix U;
  //C          if JOBZ = 'S', U contains the first min(M,N) columns of U
  //C          (the left singular vectors, stored columnwise);
  //C          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.
  //C
  //C  LDU     (input) INTEGER
  //C          The leading dimension of the array U.  LDU >= 1; if
  //C          JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.
  //C
  //C  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
  //C          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
  //C          N-by-N orthogonal matrix V**T;
  //C          if JOBZ = 'S', VT contains the first min(M,N) rows of
  //C          V**T (the right singular vectors, stored rowwise);
  //C          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.
  //C
  //C  LDVT    (input) INTEGER
  //C          The leading dimension of the array VT.  LDVT >= 1; if
  //C          JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;
  //C          if JOBZ = 'S', LDVT >= min(M,N).
  //C
  //C  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
  //C
  //C  LWORK   (input) INTEGER
  //C          The dimension of the array WORK. LWORK >= 1.
  //C          If JOBZ = 'N',
  //C            LWORK >= 3*min(M,N) + max(max(M,N),7*min(M,N)).
  //C          If JOBZ = 'O',
  //C            LWORK >= 3*min(M,N) +
  //C                     max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)).
  //C          If JOBZ = 'S' or 'A'
  //C            LWORK >= 3*min(M,N) +
  //C                     max(max(M,N),4*min(M,N)*min(M,N)+4*min(M,N)).
  //C          For good performance, LWORK should generally be larger.
  //C          If LWORK = -1 but other input arguments are legal, WORK(1)
  //C          returns the optimal LWORK.
  //C
  //C  IWORK   (workspace) INTEGER array, dimension (8*min(M,N))
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit.
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value.
  //C          > 0:  DBDSDC did not converge, updating process failed.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  Based on contributions by
  //C     Ming Gu and Huan Ren, Computer Science Division, University of
  //C     California at Berkeley, USA
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Local Arrays ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  int minmn = fem::min(m, n);
  bool wntqa = lsame(jobz, "A");
  bool wntqs = lsame(jobz, "S");
  bool wntqas = wntqa || wntqs;
  bool wntqo = lsame(jobz, "O");
  bool wntqn = lsame(jobz, "N");
  bool lquery = (lwork ==  - 1);
  //C
  if (!(wntqa || wntqs || wntqo || wntqn)) {
    info = -1;
  }
  else if (m < 0) {
    info = -2;
  }
  else if (n < 0) {
    info = -3;
  }
  else if (lda < fem::max(1, m)) {
    info = -5;
  }
  else if (ldu < 1 || (wntqas && ldu < m) || (wntqo && m < n && ldu < m)) {
    info = -8;
  }
  else if (ldvt < 1 || (wntqa && ldvt < n) || (wntqs &&
    ldvt < minmn) || (wntqo && m >= n && ldvt < n)) {
    info = -10;
  }
  //C
  //C     Compute workspace
  //C      (Note: Comments in the code beginning "Workspace:" describe the
  //C       minimal amount of workspace needed at that point in the code,
  //C       as well as the preferred amount for good performance.
  //C       NB refers to the optimal block size for the immediately
  //C       following subroutine, as returned by ILAENV.)
  //C
  int minwrk = fem::int0;
  int maxwrk = fem::int0;
  int mnthr = fem::int0;
  int bdspac = fem::int0;
  int wrkbl = fem::int0;
  if (info == 0) {
    minwrk = 1;
    maxwrk = 1;
    if (m >= n && minmn > 0) {
      //C
      //C           Compute space needed for DBDSDC
      //C
      mnthr = fem::fint(minmn * 11.0e0 / 6.0e0);
      if (wntqn) {
        bdspac = 7 * n;
      }
      else {
        bdspac = 3 * n * n + 4 * n;
      }
      if (m >= mnthr) {
        if (wntqn) {
          //C
          //C                 Path 1 (M much larger than N, JOBZ='N')
          //C
          wrkbl = n + n * ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, 3 * n + 2 * n * ilaenv(1, "DGEBRD",
            " ", n, n, -1, -1));
          maxwrk = fem::max(wrkbl, bdspac + n);
          minwrk = bdspac + n;
        }
        else if (wntqo) {
          //C
          //C                 Path 2 (M much larger than N, JOBZ='O')
          //C
          wrkbl = n + n * ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, n + n * ilaenv(1, "DORGQR", " ", m,
            n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + 2 * n * ilaenv(1, "DGEBRD",
            " ", n, n, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORMBR",
            "QLN", n, n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORMBR",
            "PRT", n, n, n, -1));
          wrkbl = fem::max(wrkbl, bdspac + 3 * n);
          maxwrk = wrkbl + 2 * n * n;
          minwrk = bdspac + 2 * n * n + 3 * n;
        }
        else if (wntqs) {
          //C
          //C                 Path 3 (M much larger than N, JOBZ='S')
          //C
          wrkbl = n + n * ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, n + n * ilaenv(1, "DORGQR", " ", m,
            n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + 2 * n * ilaenv(1, "DGEBRD",
            " ", n, n, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORMBR",
            "QLN", n, n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORMBR",
            "PRT", n, n, n, -1));
          wrkbl = fem::max(wrkbl, bdspac + 3 * n);
          maxwrk = wrkbl + n * n;
          minwrk = bdspac + n * n + 3 * n;
        }
        else if (wntqa) {
          //C
          //C                 Path 4 (M much larger than N, JOBZ='A')
          //C
          wrkbl = n + n * ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, n + m * ilaenv(1, "DORGQR", " ", m,
            m, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + 2 * n * ilaenv(1, "DGEBRD",
            " ", n, n, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORMBR",
            "QLN", n, n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORMBR",
            "PRT", n, n, n, -1));
          wrkbl = fem::max(wrkbl, bdspac + 3 * n);
          maxwrk = wrkbl + n * n;
          minwrk = bdspac + n * n + 3 * n;
        }
      }
      else {
        //C
        //C              Path 5 (M at least N, but not much larger)
        //C
        wrkbl = 3 * n + (m + n) * ilaenv(1, "DGEBRD", " ", m, n, -1, -1);
        if (wntqn) {
          maxwrk = fem::max(wrkbl, bdspac + 3 * n);
          minwrk = 3 * n + fem::max(m, bdspac);
        }
        else if (wntqo) {
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORMBR",
            "QLN", m, n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORMBR",
            "PRT", n, n, n, -1));
          wrkbl = fem::max(wrkbl, bdspac + 3 * n);
          maxwrk = wrkbl + m * n;
          minwrk = 3 * n + fem::max(m, n * n + bdspac);
        }
        else if (wntqs) {
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORMBR",
            "QLN", m, n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORMBR",
            "PRT", n, n, n, -1));
          maxwrk = fem::max(wrkbl, bdspac + 3 * n);
          minwrk = 3 * n + fem::max(m, bdspac);
        }
        else if (wntqa) {
          wrkbl = fem::max(wrkbl, 3 * n + m * ilaenv(1, "DORMBR",
            "QLN", m, m, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORMBR",
            "PRT", n, n, n, -1));
          maxwrk = fem::max(maxwrk, bdspac + 3 * n);
          minwrk = 3 * n + fem::max(m, bdspac);
        }
      }
    }
    else if (minmn > 0) {
      //C
      //C           Compute space needed for DBDSDC
      //C
      mnthr = fem::fint(minmn * 11.0e0 / 6.0e0);
      if (wntqn) {
        bdspac = 7 * m;
      }
      else {
        bdspac = 3 * m * m + 4 * m;
      }
      if (n >= mnthr) {
        if (wntqn) {
          //C
          //C                 Path 1t (N much larger than M, JOBZ='N')
          //C
          wrkbl = m + m * ilaenv(1, "DGELQF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, 3 * m + 2 * m * ilaenv(1, "DGEBRD",
            " ", m, m, -1, -1));
          maxwrk = fem::max(wrkbl, bdspac + m);
          minwrk = bdspac + m;
        }
        else if (wntqo) {
          //C
          //C                 Path 2t (N much larger than M, JOBZ='O')
          //C
          wrkbl = m + m * ilaenv(1, "DGELQF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, m + m * ilaenv(1, "DORGLQ", " ", m,
            n, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + 2 * m * ilaenv(1, "DGEBRD",
            " ", m, m, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORMBR",
            "QLN", m, m, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORMBR",
            "PRT", m, m, m, -1));
          wrkbl = fem::max(wrkbl, bdspac + 3 * m);
          maxwrk = wrkbl + 2 * m * m;
          minwrk = bdspac + 2 * m * m + 3 * m;
        }
        else if (wntqs) {
          //C
          //C                 Path 3t (N much larger than M, JOBZ='S')
          //C
          wrkbl = m + m * ilaenv(1, "DGELQF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, m + m * ilaenv(1, "DORGLQ", " ", m,
            n, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + 2 * m * ilaenv(1, "DGEBRD",
            " ", m, m, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORMBR",
            "QLN", m, m, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORMBR",
            "PRT", m, m, m, -1));
          wrkbl = fem::max(wrkbl, bdspac + 3 * m);
          maxwrk = wrkbl + m * m;
          minwrk = bdspac + m * m + 3 * m;
        }
        else if (wntqa) {
          //C
          //C                 Path 4t (N much larger than M, JOBZ='A')
          //C
          wrkbl = m + m * ilaenv(1, "DGELQF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, m + n * ilaenv(1, "DORGLQ", " ", n,
            n, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + 2 * m * ilaenv(1, "DGEBRD",
            " ", m, m, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORMBR",
            "QLN", m, m, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORMBR",
            "PRT", m, m, m, -1));
          wrkbl = fem::max(wrkbl, bdspac + 3 * m);
          maxwrk = wrkbl + m * m;
          minwrk = bdspac + m * m + 3 * m;
        }
      }
      else {
        //C
        //C              Path 5t (N greater than M, but not much larger)
        //C
        wrkbl = 3 * m + (m + n) * ilaenv(1, "DGEBRD", " ", m, n, -1, -1);
        if (wntqn) {
          maxwrk = fem::max(wrkbl, bdspac + 3 * m);
          minwrk = 3 * m + fem::max(n, bdspac);
        }
        else if (wntqo) {
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORMBR",
            "QLN", m, m, n, -1));
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORMBR",
            "PRT", m, n, m, -1));
          wrkbl = fem::max(wrkbl, bdspac + 3 * m);
          maxwrk = wrkbl + m * n;
          minwrk = 3 * m + fem::max(n, m * m + bdspac);
        }
        else if (wntqs) {
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORMBR",
            "QLN", m, m, n, -1));
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORMBR",
            "PRT", m, n, m, -1));
          maxwrk = fem::max(wrkbl, bdspac + 3 * m);
          minwrk = 3 * m + fem::max(n, bdspac);
        }
        else if (wntqa) {
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORMBR",
            "QLN", m, m, n, -1));
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORMBR",
            "PRT", n, n, m, -1));
          maxwrk = fem::max(wrkbl, bdspac + 3 * m);
          minwrk = 3 * m + fem::max(n, bdspac);
        }
      }
    }
    maxwrk = fem::max(maxwrk, minwrk);
    work(1) = maxwrk;
    //C
    if (lwork < minwrk && !lquery) {
      info = -12;
    }
  }
  //C
  if (info != 0) {
    xerbla("DGESDD", -info);
    return;
  }
  else if (lquery) {
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (m == 0 || n == 0) {
    return;
  }
  //C
  //C     Get machine constants
  //C
  double eps = dlamch(cmn, "P");
  double smlnum = fem::sqrt(dlamch(cmn, "S")) / eps;
  const double one = 1.0e0;
  double bignum = one / smlnum;
  //C
  //C     Scale A if max element outside range [SMLNUM,BIGNUM]
  //C
  arr_1d<1, double> dum(fem::fill0);
  double anrm = dlange("M", m, n, a, lda, dum);
  int iscl = 0;
  const double zero = 0.0e0;
  int ierr = fem::int0;
  if (anrm > zero && anrm < smlnum) {
    iscl = 1;
    dlascl(cmn, "G", 0, 0, anrm, smlnum, m, n, a, lda, ierr);
  }
  else if (anrm > bignum) {
    iscl = 1;
    dlascl(cmn, "G", 0, 0, anrm, bignum, m, n, a, lda, ierr);
  }
  //C
  int itau = fem::int0;
  int nwork = fem::int0;
  int ie = fem::int0;
  int itauq = fem::int0;
  int itaup = fem::int0;
  arr_1d<1, int> idum(fem::fill0);
  int ir = fem::int0;
  int ldwrkr = fem::int0;
  int iu = fem::int0;
  int i = fem::int0;
  int chunk = fem::int0;
  int ldwrku = fem::int0;
  int ivt = fem::int0;
  int il = fem::int0;
  int ldwrkl = fem::int0;
  int blk = fem::int0;
  int ldwkvt = fem::int0;
  if (m >= n) {
    //C
    //C        A has at least as many rows as columns. If A has sufficiently
    //C        more rows than columns, first reduce using the QR
    //C        decomposition (if sufficient workspace available)
    //C
    if (m >= mnthr) {
      //C
      if (wntqn) {
        //C
        //C              Path 1 (M much larger than N, JOBZ='N')
        //C              No singular vectors to be computed
        //C
        itau = 1;
        nwork = itau + n;
        //C
        //C              Compute A=Q*R
        //C              (Workspace: need 2*N, prefer N+N*NB)
        //C
        dgeqrf(cmn, m, n, a, lda, work(itau), work(nwork), lwork - nwork + 1,
          ierr);
        //C
        //C              Zero out below R
        //C
        dlaset("L", n - 1, n - 1, zero, zero, a(2, 1), lda);
        ie = 1;
        itauq = ie + n;
        itaup = itauq + n;
        nwork = itaup + n;
        //C
        //C              Bidiagonalize R in A
        //C              (Workspace: need 4*N, prefer 3*N+2*N*NB)
        //C
        dgebrd(cmn, n, n, a, lda, s, work(ie), work(itauq), work(itaup),
          work(nwork), lwork - nwork + 1, ierr);
        nwork = ie + n;
        //C
        //C              Perform bidiagonal SVD, computing singular values only
        //C              (Workspace: need N+BDSPAC)
        //C
        dbdsdc(cmn, "U", "N", n, s, work(ie), dum, 1, dum, 1, dum,
          idum, work(nwork), iwork, info);
        //C
      }
      else if (wntqo) {
        //C
        //C              Path 2 (M much larger than N, JOBZ = 'O')
        //C              N left singular vectors to be overwritten on A and
        //C              N right singular vectors to be computed in VT
        //C
        ir = 1;
        //C
        //C              WORK(IR) is LDWRKR by N
        //C
        if (lwork >= lda * n + n * n + 3 * n + bdspac) {
          ldwrkr = lda;
        }
        else {
          ldwrkr = (lwork - n * n - 3 * n - bdspac) / n;
        }
        itau = ir + ldwrkr * n;
        nwork = itau + n;
        //C
        //C              Compute A=Q*R
        //C              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
        //C
        dgeqrf(cmn, m, n, a, lda, work(itau), work(nwork), lwork - nwork + 1,
          ierr);
        //C
        //C              Copy R to WORK(IR), zeroing out below it
        //C
        dlacpy("U", n, n, a, lda, work(ir), ldwrkr);
        dlaset("L", n - 1, n - 1, zero, zero, work(ir + 1), ldwrkr);
        //C
        //C              Generate Q in A
        //C              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
        //C
        dorgqr(m, n, n, a, lda, work(itau), work(nwork), lwork - nwork + 1,
          ierr);
        ie = itau;
        itauq = ie + n;
        itaup = itauq + n;
        nwork = itaup + n;
        //C
        //C              Bidiagonalize R in VT, copying result to WORK(IR)
        //C              (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
        //C
        dgebrd(cmn, n, n, work(ir), ldwrkr, s, work(ie), work(itauq),
          work(itaup), work(nwork), lwork - nwork + 1, ierr);
        //C
        //C              WORK(IU) is N by N
        //C
        iu = nwork;
        nwork = iu + n * n;
        //C
        //C              Perform bidiagonal SVD, computing left singular vectors
        //C              of bidiagonal matrix in WORK(IU) and computing right
        //C              singular vectors of bidiagonal matrix in VT
        //C              (Workspace: need N+N*N+BDSPAC)
        //C
        dbdsdc(cmn, "U", "I", n, s, work(ie), work(iu), n, vt, ldvt,
          dum, idum, work(nwork), iwork, info);
        //C
        //C              Overwrite WORK(IU) by left singular vectors of R
        //C              and VT by right singular vectors of R
        //C              (Workspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB)
        //C
        dormbr("Q", "L", "N", n, n, n, work(ir), ldwrkr, work(itauq),
          work(iu), n, work(nwork), lwork - nwork + 1, ierr);
        dormbr("P", "R", "T", n, n, n, work(ir), ldwrkr, work(itaup),
          vt, ldvt, work(nwork), lwork - nwork + 1, ierr);
        //C
        //C              Multiply Q in A by left singular vectors of R in
        //C              WORK(IU), storing result in WORK(IR) and copying to A
        //C              (Workspace: need 2*N*N, prefer N*N+M*N)
        //C
        FEM_DOSTEP(i, 1, m, ldwrkr) {
          chunk = fem::min(m - i + 1, ldwrkr);
          dgemm("N", "N", chunk, n, n, one, a(i, 1), lda, work(iu),
            n, zero, work(ir), ldwrkr);
          dlacpy("F", chunk, n, work(ir), ldwrkr, a(i, 1), lda);
        }
        //C
      }
      else if (wntqs) {
        //C
        //C              Path 3 (M much larger than N, JOBZ='S')
        //C              N left singular vectors to be computed in U and
        //C              N right singular vectors to be computed in VT
        //C
        ir = 1;
        //C
        //C              WORK(IR) is N by N
        //C
        ldwrkr = n;
        itau = ir + ldwrkr * n;
        nwork = itau + n;
        //C
        //C              Compute A=Q*R
        //C              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
        //C
        dgeqrf(cmn, m, n, a, lda, work(itau), work(nwork), lwork - nwork + 1,
          ierr);
        //C
        //C              Copy R to WORK(IR), zeroing out below it
        //C
        dlacpy("U", n, n, a, lda, work(ir), ldwrkr);
        dlaset("L", n - 1, n - 1, zero, zero, work(ir + 1), ldwrkr);
        //C
        //C              Generate Q in A
        //C              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
        //C
        dorgqr(m, n, n, a, lda, work(itau), work(nwork), lwork - nwork + 1,
          ierr);
        ie = itau;
        itauq = ie + n;
        itaup = itauq + n;
        nwork = itaup + n;
        //C
        //C              Bidiagonalize R in WORK(IR)
        //C              (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
        //C
        dgebrd(cmn, n, n, work(ir), ldwrkr, s, work(ie), work(itauq),
          work(itaup), work(nwork), lwork - nwork + 1, ierr);
        //C
        //C              Perform bidiagonal SVD, computing left singular vectors
        //C              of bidiagoal matrix in U and computing right singular
        //C              vectors of bidiagonal matrix in VT
        //C              (Workspace: need N+BDSPAC)
        //C
        dbdsdc(cmn, "U", "I", n, s, work(ie), u, ldu, vt, ldvt, dum,
          idum, work(nwork), iwork, info);
        //C
        //C              Overwrite U by left singular vectors of R and VT
        //C              by right singular vectors of R
        //C              (Workspace: need N*N+3*N, prefer N*N+2*N+N*NB)
        //C
        dormbr("Q", "L", "N", n, n, n, work(ir), ldwrkr, work(itauq),
          u, ldu, work(nwork), lwork - nwork + 1, ierr);
        //C
        dormbr("P", "R", "T", n, n, n, work(ir), ldwrkr, work(itaup),
          vt, ldvt, work(nwork), lwork - nwork + 1, ierr);
        //C
        //C              Multiply Q in A by left singular vectors of R in
        //C              WORK(IR), storing result in U
        //C              (Workspace: need N*N)
        //C
        dlacpy("F", n, n, u, ldu, work(ir), ldwrkr);
        dgemm("N", "N", m, n, n, one, a, lda, work(ir), ldwrkr, zero, u, ldu);
        //C
      }
      else if (wntqa) {
        //C
        //C              Path 4 (M much larger than N, JOBZ='A')
        //C              M left singular vectors to be computed in U and
        //C              N right singular vectors to be computed in VT
        //C
        iu = 1;
        //C
        //C              WORK(IU) is N by N
        //C
        ldwrku = n;
        itau = iu + ldwrku * n;
        nwork = itau + n;
        //C
        //C              Compute A=Q*R, copying result to U
        //C              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
        //C
        dgeqrf(cmn, m, n, a, lda, work(itau), work(nwork), lwork - nwork + 1,
          ierr);
        dlacpy("L", m, n, a, lda, u, ldu);
        //C
        //C              Generate Q in U
        //C              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
        dorgqr(m, m, n, u, ldu, work(itau), work(nwork), lwork - nwork + 1,
          ierr);
        //C
        //C              Produce R in A, zeroing out other entries
        //C
        dlaset("L", n - 1, n - 1, zero, zero, a(2, 1), lda);
        ie = itau;
        itauq = ie + n;
        itaup = itauq + n;
        nwork = itaup + n;
        //C
        //C              Bidiagonalize R in A
        //C              (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
        //C
        dgebrd(cmn, n, n, a, lda, s, work(ie), work(itauq), work(itaup),
          work(nwork), lwork - nwork + 1, ierr);
        //C
        //C              Perform bidiagonal SVD, computing left singular vectors
        //C              of bidiagonal matrix in WORK(IU) and computing right
        //C              singular vectors of bidiagonal matrix in VT
        //C              (Workspace: need N+N*N+BDSPAC)
        //C
        dbdsdc(cmn, "U", "I", n, s, work(ie), work(iu), n, vt, ldvt,
          dum, idum, work(nwork), iwork, info);
        //C
        //C              Overwrite WORK(IU) by left singular vectors of R and VT
        //C              by right singular vectors of R
        //C              (Workspace: need N*N+3*N, prefer N*N+2*N+N*NB)
        //C
        dormbr("Q", "L", "N", n, n, n, a, lda, work(itauq), work(iu),
          ldwrku, work(nwork), lwork - nwork + 1, ierr);
        dormbr("P", "R", "T", n, n, n, a, lda, work(itaup), vt, ldvt,
          work(nwork), lwork - nwork + 1, ierr);
        //C
        //C              Multiply Q in U by left singular vectors of R in
        //C              WORK(IU), storing result in A
        //C              (Workspace: need N*N)
        //C
        dgemm("N", "N", m, n, n, one, u, ldu, work(iu), ldwrku, zero, a, lda);
        //C
        //C              Copy left singular vectors of A from A to U
        //C
        dlacpy("F", m, n, a, lda, u, ldu);
        //C
      }
      //C
    }
    else {
      //C
      //C           M .LT. MNTHR
      //C
      //C           Path 5 (M at least N, but not much larger)
      //C           Reduce to bidiagonal form without QR decomposition
      //C
      ie = 1;
      itauq = ie + n;
      itaup = itauq + n;
      nwork = itaup + n;
      //C
      //C           Bidiagonalize A
      //C           (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB)
      //C
      dgebrd(cmn, m, n, a, lda, s, work(ie), work(itauq), work(itaup),
        work(nwork), lwork - nwork + 1, ierr);
      if (wntqn) {
        //C
        //C              Perform bidiagonal SVD, only computing singular values
        //C              (Workspace: need N+BDSPAC)
        //C
        dbdsdc(cmn, "U", "N", n, s, work(ie), dum, 1, dum, 1, dum,
          idum, work(nwork), iwork, info);
      }
      else if (wntqo) {
        iu = nwork;
        if (lwork >= m * n + 3 * n + bdspac) {
          //C
          //C                 WORK( IU ) is M by N
          //C
          ldwrku = m;
          nwork = iu + ldwrku * n;
          dlaset("F", m, n, zero, zero, work(iu), ldwrku);
        }
        else {
          //C
          //C                 WORK( IU ) is N by N
          //C
          ldwrku = n;
          nwork = iu + ldwrku * n;
          //C
          //C                 WORK(IR) is LDWRKR by N
          //C
          ir = nwork;
          ldwrkr = (lwork - n * n - 3 * n) / n;
        }
        nwork = iu + ldwrku * n;
        //C
        //C              Perform bidiagonal SVD, computing left singular vectors
        //C              of bidiagonal matrix in WORK(IU) and computing right
        //C              singular vectors of bidiagonal matrix in VT
        //C              (Workspace: need N+N*N+BDSPAC)
        //C
        dbdsdc(cmn, "U", "I", n, s, work(ie), work(iu), ldwrku, vt,
          ldvt, dum, idum, work(nwork), iwork, info);
        //C
        //C              Overwrite VT by right singular vectors of A
        //C              (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
        //C
        dormbr("P", "R", "T", n, n, n, a, lda, work(itaup), vt, ldvt,
          work(nwork), lwork - nwork + 1, ierr);
        //C
        if (lwork >= m * n + 3 * n + bdspac) {
          //C
          //C                 Overwrite WORK(IU) by left singular vectors of A
          //C                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
          //C
          dormbr("Q", "L", "N", m, n, n, a, lda, work(itauq), work(iu),
            ldwrku, work(nwork), lwork - nwork + 1, ierr);
          //C
          //C                 Copy left singular vectors of A from WORK(IU) to A
          //C
          dlacpy("F", m, n, work(iu), ldwrku, a, lda);
        }
        else {
          //C
          //C                 Generate Q in A
          //C                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
          //C
          dorgbr("Q", m, n, n, a, lda, work(itauq), work(nwork),
            lwork - nwork + 1, ierr);
          //C
          //C                 Multiply Q in A by left singular vectors of
          //C                 bidiagonal matrix in WORK(IU), storing result in
          //C                 WORK(IR) and copying to A
          //C                 (Workspace: need 2*N*N, prefer N*N+M*N)
          //C
          FEM_DOSTEP(i, 1, m, ldwrkr) {
            chunk = fem::min(m - i + 1, ldwrkr);
            dgemm("N", "N", chunk, n, n, one, a(i, 1), lda, work(iu),
              ldwrku, zero, work(ir), ldwrkr);
            dlacpy("F", chunk, n, work(ir), ldwrkr, a(i, 1), lda);
          }
        }
        //C
      }
      else if (wntqs) {
        //C
        //C              Perform bidiagonal SVD, computing left singular vectors
        //C              of bidiagonal matrix in U and computing right singular
        //C              vectors of bidiagonal matrix in VT
        //C              (Workspace: need N+BDSPAC)
        //C
        dlaset("F", m, n, zero, zero, u, ldu);
        dbdsdc(cmn, "U", "I", n, s, work(ie), u, ldu, vt, ldvt, dum,
          idum, work(nwork), iwork, info);
        //C
        //C              Overwrite U by left singular vectors of A and VT
        //C              by right singular vectors of A
        //C              (Workspace: need 3*N, prefer 2*N+N*NB)
        //C
        dormbr("Q", "L", "N", m, n, n, a, lda, work(itauq), u, ldu,
          work(nwork), lwork - nwork + 1, ierr);
        dormbr("P", "R", "T", n, n, n, a, lda, work(itaup), vt, ldvt,
          work(nwork), lwork - nwork + 1, ierr);
      }
      else if (wntqa) {
        //C
        //C              Perform bidiagonal SVD, computing left singular vectors
        //C              of bidiagonal matrix in U and computing right singular
        //C              vectors of bidiagonal matrix in VT
        //C              (Workspace: need N+BDSPAC)
        //C
        dlaset("F", m, m, zero, zero, u, ldu);
        dbdsdc(cmn, "U", "I", n, s, work(ie), u, ldu, vt, ldvt, dum,
          idum, work(nwork), iwork, info);
        //C
        //C              Set the right corner of U to identity matrix
        //C
        if (m > n) {
          dlaset("F", m - n, m - n, zero, one, u(n + 1, n + 1), ldu);
        }
        //C
        //C              Overwrite U by left singular vectors of A and VT
        //C              by right singular vectors of A
        //C              (Workspace: need N*N+2*N+M, prefer N*N+2*N+M*NB)
        //C
        dormbr("Q", "L", "N", m, m, n, a, lda, work(itauq), u, ldu,
          work(nwork), lwork - nwork + 1, ierr);
        dormbr("P", "R", "T", n, n, m, a, lda, work(itaup), vt, ldvt,
          work(nwork), lwork - nwork + 1, ierr);
      }
      //C
    }
    //C
  }
  else {
    //C
    //C        A has more columns than rows. If A has sufficiently more
    //C        columns than rows, first reduce using the LQ decomposition (if
    //C        sufficient workspace available)
    //C
    if (n >= mnthr) {
      //C
      if (wntqn) {
        //C
        //C              Path 1t (N much larger than M, JOBZ='N')
        //C              No singular vectors to be computed
        //C
        itau = 1;
        nwork = itau + m;
        //C
        //C              Compute A=L*Q
        //C              (Workspace: need 2*M, prefer M+M*NB)
        //C
        dgelqf(cmn, m, n, a, lda, work(itau), work(nwork), lwork - nwork + 1,
          ierr);
        //C
        //C              Zero out above L
        //C
        dlaset("U", m - 1, m - 1, zero, zero, a(1, 2), lda);
        ie = 1;
        itauq = ie + m;
        itaup = itauq + m;
        nwork = itaup + m;
        //C
        //C              Bidiagonalize L in A
        //C              (Workspace: need 4*M, prefer 3*M+2*M*NB)
        //C
        dgebrd(cmn, m, m, a, lda, s, work(ie), work(itauq), work(itaup),
          work(nwork), lwork - nwork + 1, ierr);
        nwork = ie + m;
        //C
        //C              Perform bidiagonal SVD, computing singular values only
        //C              (Workspace: need M+BDSPAC)
        //C
        dbdsdc(cmn, "U", "N", m, s, work(ie), dum, 1, dum, 1, dum,
          idum, work(nwork), iwork, info);
        //C
      }
      else if (wntqo) {
        //C
        //C              Path 2t (N much larger than M, JOBZ='O')
        //C              M right singular vectors to be overwritten on A and
        //C              M left singular vectors to be computed in U
        //C
        ivt = 1;
        //C
        //C              IVT is M by M
        //C
        il = ivt + m * m;
        if (lwork >= m * n + m * m + 3 * m + bdspac) {
          //C
          //C                 WORK(IL) is M by N
          //C
          ldwrkl = m;
          chunk = n;
        }
        else {
          ldwrkl = m;
          chunk = (lwork - m * m) / m;
        }
        itau = il + ldwrkl * m;
        nwork = itau + m;
        //C
        //C              Compute A=L*Q
        //C              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
        //C
        dgelqf(cmn, m, n, a, lda, work(itau), work(nwork), lwork - nwork + 1,
          ierr);
        //C
        //C              Copy L to WORK(IL), zeroing about above it
        //C
        dlacpy("L", m, m, a, lda, work(il), ldwrkl);
        dlaset("U", m - 1, m - 1, zero, zero, work(il + ldwrkl), ldwrkl);
        //C
        //C              Generate Q in A
        //C              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
        //C
        dorglq(m, n, m, a, lda, work(itau), work(nwork), lwork - nwork + 1,
          ierr);
        ie = itau;
        itauq = ie + m;
        itaup = itauq + m;
        nwork = itaup + m;
        //C
        //C              Bidiagonalize L in WORK(IL)
        //C              (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
        //C
        dgebrd(cmn, m, m, work(il), ldwrkl, s, work(ie), work(itauq),
          work(itaup), work(nwork), lwork - nwork + 1, ierr);
        //C
        //C              Perform bidiagonal SVD, computing left singular vectors
        //C              of bidiagonal matrix in U, and computing right singular
        //C              vectors of bidiagonal matrix in WORK(IVT)
        //C              (Workspace: need M+M*M+BDSPAC)
        //C
        dbdsdc(cmn, "U", "I", m, s, work(ie), u, ldu, work(ivt), m,
          dum, idum, work(nwork), iwork, info);
        //C
        //C              Overwrite U by left singular vectors of L and WORK(IVT)
        //C              by right singular vectors of L
        //C              (Workspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB)
        //C
        dormbr("Q", "L", "N", m, m, m, work(il), ldwrkl, work(itauq),
          u, ldu, work(nwork), lwork - nwork + 1, ierr);
        dormbr("P", "R", "T", m, m, m, work(il), ldwrkl, work(itaup),
          work(ivt), m, work(nwork), lwork - nwork + 1, ierr);
        //C
        //C              Multiply right singular vectors of L in WORK(IVT) by Q
        //C              in A, storing result in WORK(IL) and copying to A
        //C              (Workspace: need 2*M*M, prefer M*M+M*N)
        //C
        FEM_DOSTEP(i, 1, n, chunk) {
          blk = fem::min(n - i + 1, chunk);
          dgemm("N", "N", m, blk, m, one, work(ivt), m, a(1, i), lda,
            zero, work(il), ldwrkl);
          dlacpy("F", m, blk, work(il), ldwrkl, a(1, i), lda);
        }
        //C
      }
      else if (wntqs) {
        //C
        //C              Path 3t (N much larger than M, JOBZ='S')
        //C              M right singular vectors to be computed in VT and
        //C              M left singular vectors to be computed in U
        //C
        il = 1;
        //C
        //C              WORK(IL) is M by M
        //C
        ldwrkl = m;
        itau = il + ldwrkl * m;
        nwork = itau + m;
        //C
        //C              Compute A=L*Q
        //C              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
        //C
        dgelqf(cmn, m, n, a, lda, work(itau), work(nwork), lwork - nwork + 1,
          ierr);
        //C
        //C              Copy L to WORK(IL), zeroing out above it
        //C
        dlacpy("L", m, m, a, lda, work(il), ldwrkl);
        dlaset("U", m - 1, m - 1, zero, zero, work(il + ldwrkl), ldwrkl);
        //C
        //C              Generate Q in A
        //C              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
        //C
        dorglq(m, n, m, a, lda, work(itau), work(nwork), lwork - nwork + 1,
          ierr);
        ie = itau;
        itauq = ie + m;
        itaup = itauq + m;
        nwork = itaup + m;
        //C
        //C              Bidiagonalize L in WORK(IU), copying result to U
        //C              (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
        //C
        dgebrd(cmn, m, m, work(il), ldwrkl, s, work(ie), work(itauq),
          work(itaup), work(nwork), lwork - nwork + 1, ierr);
        //C
        //C              Perform bidiagonal SVD, computing left singular vectors
        //C              of bidiagonal matrix in U and computing right singular
        //C              vectors of bidiagonal matrix in VT
        //C              (Workspace: need M+BDSPAC)
        //C
        dbdsdc(cmn, "U", "I", m, s, work(ie), u, ldu, vt, ldvt, dum,
          idum, work(nwork), iwork, info);
        //C
        //C              Overwrite U by left singular vectors of L and VT
        //C              by right singular vectors of L
        //C              (Workspace: need M*M+3*M, prefer M*M+2*M+M*NB)
        //C
        dormbr("Q", "L", "N", m, m, m, work(il), ldwrkl, work(itauq),
          u, ldu, work(nwork), lwork - nwork + 1, ierr);
        dormbr("P", "R", "T", m, m, m, work(il), ldwrkl, work(itaup),
          vt, ldvt, work(nwork), lwork - nwork + 1, ierr);
        //C
        //C              Multiply right singular vectors of L in WORK(IL) by
        //C              Q in A, storing result in VT
        //C              (Workspace: need M*M)
        //C
        dlacpy("F", m, m, vt, ldvt, work(il), ldwrkl);
        dgemm("N", "N", m, n, m, one, work(il), ldwrkl, a, lda, zero, vt, ldvt);
        //C
      }
      else if (wntqa) {
        //C
        //C              Path 4t (N much larger than M, JOBZ='A')
        //C              N right singular vectors to be computed in VT and
        //C              M left singular vectors to be computed in U
        //C
        ivt = 1;
        //C
        //C              WORK(IVT) is M by M
        //C
        ldwkvt = m;
        itau = ivt + ldwkvt * m;
        nwork = itau + m;
        //C
        //C              Compute A=L*Q, copying result to VT
        //C              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
        //C
        dgelqf(cmn, m, n, a, lda, work(itau), work(nwork), lwork - nwork + 1,
          ierr);
        dlacpy("U", m, n, a, lda, vt, ldvt);
        //C
        //C              Generate Q in VT
        //C              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
        //C
        dorglq(n, n, m, vt, ldvt, work(itau), work(nwork), lwork - nwork + 1,
          ierr);
        //C
        //C              Produce L in A, zeroing out other entries
        //C
        dlaset("U", m - 1, m - 1, zero, zero, a(1, 2), lda);
        ie = itau;
        itauq = ie + m;
        itaup = itauq + m;
        nwork = itaup + m;
        //C
        //C              Bidiagonalize L in A
        //C              (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
        //C
        dgebrd(cmn, m, m, a, lda, s, work(ie), work(itauq), work(itaup),
          work(nwork), lwork - nwork + 1, ierr);
        //C
        //C              Perform bidiagonal SVD, computing left singular vectors
        //C              of bidiagonal matrix in U and computing right singular
        //C              vectors of bidiagonal matrix in WORK(IVT)
        //C              (Workspace: need M+M*M+BDSPAC)
        //C
        dbdsdc(cmn, "U", "I", m, s, work(ie), u, ldu, work(ivt),
          ldwkvt, dum, idum, work(nwork), iwork, info);
        //C
        //C              Overwrite U by left singular vectors of L and WORK(IVT)
        //C              by right singular vectors of L
        //C              (Workspace: need M*M+3*M, prefer M*M+2*M+M*NB)
        //C
        dormbr("Q", "L", "N", m, m, m, a, lda, work(itauq), u, ldu,
          work(nwork), lwork - nwork + 1, ierr);
        dormbr("P", "R", "T", m, m, m, a, lda, work(itaup), work(ivt),
          ldwkvt, work(nwork), lwork - nwork + 1, ierr);
        //C
        //C              Multiply right singular vectors of L in WORK(IVT) by
        //C              Q in VT, storing result in A
        //C              (Workspace: need M*M)
        //C
        dgemm("N", "N", m, n, m, one, work(ivt), ldwkvt, vt, ldvt,
          zero, a, lda);
        //C
        //C              Copy right singular vectors of A from A to VT
        //C
        dlacpy("F", m, n, a, lda, vt, ldvt);
        //C
      }
      //C
    }
    else {
      //C
      //C           N .LT. MNTHR
      //C
      //C           Path 5t (N greater than M, but not much larger)
      //C           Reduce to bidiagonal form without LQ decomposition
      //C
      ie = 1;
      itauq = ie + m;
      itaup = itauq + m;
      nwork = itaup + m;
      //C
      //C           Bidiagonalize A
      //C           (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
      //C
      dgebrd(cmn, m, n, a, lda, s, work(ie), work(itauq), work(itaup),
        work(nwork), lwork - nwork + 1, ierr);
      if (wntqn) {
        //C
        //C              Perform bidiagonal SVD, only computing singular values
        //C              (Workspace: need M+BDSPAC)
        //C
        dbdsdc(cmn, "L", "N", m, s, work(ie), dum, 1, dum, 1, dum,
          idum, work(nwork), iwork, info);
      }
      else if (wntqo) {
        ldwkvt = m;
        ivt = nwork;
        if (lwork >= m * n + 3 * m + bdspac) {
          //C
          //C                 WORK( IVT ) is M by N
          //C
          dlaset("F", m, n, zero, zero, work(ivt), ldwkvt);
          nwork = ivt + ldwkvt * n;
        }
        else {
          //C
          //C                 WORK( IVT ) is M by M
          //C
          nwork = ivt + ldwkvt * m;
          il = nwork;
          //C
          //C                 WORK(IL) is M by CHUNK
          //C
          chunk = (lwork - m * m - 3 * m) / m;
        }
        //C
        //C              Perform bidiagonal SVD, computing left singular vectors
        //C              of bidiagonal matrix in U and computing right singular
        //C              vectors of bidiagonal matrix in WORK(IVT)
        //C              (Workspace: need M*M+BDSPAC)
        //C
        dbdsdc(cmn, "L", "I", m, s, work(ie), u, ldu, work(ivt),
          ldwkvt, dum, idum, work(nwork), iwork, info);
        //C
        //C              Overwrite U by left singular vectors of A
        //C              (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
        //C
        dormbr("Q", "L", "N", m, m, n, a, lda, work(itauq), u, ldu,
          work(nwork), lwork - nwork + 1, ierr);
        //C
        if (lwork >= m * n + 3 * m + bdspac) {
          //C
          //C                 Overwrite WORK(IVT) by left singular vectors of A
          //C                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
          //C
          dormbr("P", "R", "T", m, n, m, a, lda, work(itaup), work(ivt),
            ldwkvt, work(nwork), lwork - nwork + 1, ierr);
          //C
          //C                 Copy right singular vectors of A from WORK(IVT) to A
          //C
          dlacpy("F", m, n, work(ivt), ldwkvt, a, lda);
        }
        else {
          //C
          //C                 Generate P**T in A
          //C                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
          //C
          dorgbr("P", m, n, m, a, lda, work(itaup), work(nwork),
            lwork - nwork + 1, ierr);
          //C
          //C                 Multiply Q in A by right singular vectors of
          //C                 bidiagonal matrix in WORK(IVT), storing result in
          //C                 WORK(IL) and copying to A
          //C                 (Workspace: need 2*M*M, prefer M*M+M*N)
          //C
          FEM_DOSTEP(i, 1, n, chunk) {
            blk = fem::min(n - i + 1, chunk);
            dgemm("N", "N", m, blk, m, one, work(ivt), ldwkvt, a(1,
              i), lda, zero, work(il), m);
            dlacpy("F", m, blk, work(il), m, a(1, i), lda);
          }
        }
      }
      else if (wntqs) {
        //C
        //C              Perform bidiagonal SVD, computing left singular vectors
        //C              of bidiagonal matrix in U and computing right singular
        //C              vectors of bidiagonal matrix in VT
        //C              (Workspace: need M+BDSPAC)
        //C
        dlaset("F", m, n, zero, zero, vt, ldvt);
        dbdsdc(cmn, "L", "I", m, s, work(ie), u, ldu, vt, ldvt, dum,
          idum, work(nwork), iwork, info);
        //C
        //C              Overwrite U by left singular vectors of A and VT
        //C              by right singular vectors of A
        //C              (Workspace: need 3*M, prefer 2*M+M*NB)
        //C
        dormbr("Q", "L", "N", m, m, n, a, lda, work(itauq), u, ldu,
          work(nwork), lwork - nwork + 1, ierr);
        dormbr("P", "R", "T", m, n, m, a, lda, work(itaup), vt, ldvt,
          work(nwork), lwork - nwork + 1, ierr);
      }
      else if (wntqa) {
        //C
        //C              Perform bidiagonal SVD, computing left singular vectors
        //C              of bidiagonal matrix in U and computing right singular
        //C              vectors of bidiagonal matrix in VT
        //C              (Workspace: need M+BDSPAC)
        //C
        dlaset("F", n, n, zero, zero, vt, ldvt);
        dbdsdc(cmn, "L", "I", m, s, work(ie), u, ldu, vt, ldvt, dum,
          idum, work(nwork), iwork, info);
        //C
        //C              Set the right corner of VT to identity matrix
        //C
        if (n > m) {
          dlaset("F", n - m, n - m, zero, one, vt(m + 1, m + 1), ldvt);
        }
        //C
        //C              Overwrite U by left singular vectors of A and VT
        //C              by right singular vectors of A
        //C              (Workspace: need 2*M+N, prefer 2*M+N*NB)
        //C
        dormbr("Q", "L", "N", m, m, n, a, lda, work(itauq), u, ldu,
          work(nwork), lwork - nwork + 1, ierr);
        dormbr("P", "R", "T", n, n, m, a, lda, work(itaup), vt, ldvt,
          work(nwork), lwork - nwork + 1, ierr);
      }
      //C
    }
    //C
  }
  //C
  //C     Undo scaling if necessary
  //C
  if (iscl == 1) {
    if (anrm > bignum) {
      dlascl(cmn, "G", 0, 0, bignum, anrm, minmn, 1, s, minmn, ierr);
    }
    if (anrm < smlnum) {
      dlascl(cmn, "G", 0, 0, smlnum, anrm, minmn, 1, s, minmn, ierr);
    }
  }
  //C
  //C     Return optimal workspace in WORK(1)
  //C
  work(1) = maxwrk;
  //C
  //C     End of DGESDD
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dgesvd.f
inline
void
dgesvd(
  common& cmn,
  str_cref jobu,
  str_cref jobvt,
  int const& m,
  int const& n,
  arr_ref<double, 2> a,
  int const& lda,
  arr_ref<double> s,
  arr_ref<double, 2> u,
  int const& ldu,
  arr_ref<double, 2> vt,
  int const& ldvt,
  arr_ref<double> work,
  int const& lwork,
  int& info)
{
  a(dimension(lda, star));
  s(dimension(star));
  u(dimension(ldu, star));
  vt(dimension(ldvt, star));
  work(dimension(star));
  //C
  //C  -- LAPACK driver routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DGESVD computes the singular value decomposition (SVD) of a real
  //C  M-by-N matrix A, optionally computing the left and/or right singular
  //C  vectors. The SVD is written
  //C
  //C       A = U * SIGMA * transpose(V)
  //C
  //C  where SIGMA is an M-by-N matrix which is zero except for its
  //C  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
  //C  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
  //C  are the singular values of A; they are real and non-negative, and
  //C  are returned in descending order.  The first min(m,n) columns of
  //C  U and V are the left and right singular vectors of A.
  //C
  //C  Note that the routine returns V**T, not V.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  JOBU    (input) CHARACTER*1
  //C          Specifies options for computing all or part of the matrix U:
  //C          = 'A':  all M columns of U are returned in array U:
  //C          = 'S':  the first min(m,n) columns of U (the left singular
  //C                  vectors) are returned in the array U;
  //C          = 'O':  the first min(m,n) columns of U (the left singular
  //C                  vectors) are overwritten on the array A;
  //C          = 'N':  no columns of U (no left singular vectors) are
  //C                  computed.
  //C
  //C  JOBVT   (input) CHARACTER*1
  //C          Specifies options for computing all or part of the matrix
  //C          V**T:
  //C          = 'A':  all N rows of V**T are returned in the array VT;
  //C          = 'S':  the first min(m,n) rows of V**T (the right singular
  //C                  vectors) are returned in the array VT;
  //C          = 'O':  the first min(m,n) rows of V**T (the right singular
  //C                  vectors) are overwritten on the array A;
  //C          = 'N':  no rows of V**T (no right singular vectors) are
  //C                  computed.
  //C
  //C          JOBVT and JOBU cannot both be 'O'.
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the input matrix A.  M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the input matrix A.  N >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the M-by-N matrix A.
  //C          On exit,
  //C          if JOBU = 'O',  A is overwritten with the first min(m,n)
  //C                          columns of U (the left singular vectors,
  //C                          stored columnwise);
  //C          if JOBVT = 'O', A is overwritten with the first min(m,n)
  //C                          rows of V**T (the right singular vectors,
  //C                          stored rowwise);
  //C          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
  //C                          are destroyed.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(1,M).
  //C
  //C  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
  //C          The singular values of A, sorted so that S(i) >= S(i+1).
  //C
  //C  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
  //C          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
  //C          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
  //C          if JOBU = 'S', U contains the first min(m,n) columns of U
  //C          (the left singular vectors, stored columnwise);
  //C          if JOBU = 'N' or 'O', U is not referenced.
  //C
  //C  LDU     (input) INTEGER
  //C          The leading dimension of the array U.  LDU >= 1; if
  //C          JOBU = 'S' or 'A', LDU >= M.
  //C
  //C  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
  //C          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
  //C          V**T;
  //C          if JOBVT = 'S', VT contains the first min(m,n) rows of
  //C          V**T (the right singular vectors, stored rowwise);
  //C          if JOBVT = 'N' or 'O', VT is not referenced.
  //C
  //C  LDVT    (input) INTEGER
  //C          The leading dimension of the array VT.  LDVT >= 1; if
  //C          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
  //C
  //C  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
  //C          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
  //C          superdiagonal elements of an upper bidiagonal matrix B
  //C          whose diagonal is in S (not necessarily sorted). B
  //C          satisfies A = U * B * VT, so it has the same singular values
  //C          as A, and singular vectors related by U and VT.
  //C
  //C  LWORK   (input) INTEGER
  //C          The dimension of the array WORK.
  //C          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).
  //C          For good performance, LWORK should generally be larger.
  //C
  //C          If LWORK = -1, then a workspace query is assumed; the routine
  //C          only calculates the optimal size of the WORK array, returns
  //C          this value as the first entry of the WORK array, and no error
  //C          message related to LWORK is issued by XERBLA.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit.
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value.
  //C          > 0:  if DBDSQR did not converge, INFO specifies how many
  //C                superdiagonals of an intermediate bidiagonal form B
  //C                did not converge to zero. See the description of WORK
  //C                above for details.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Local Arrays ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  int minmn = fem::min(m, n);
  bool wntua = lsame(jobu, "A");
  bool wntus = lsame(jobu, "S");
  bool wntuas = wntua || wntus;
  bool wntuo = lsame(jobu, "O");
  bool wntun = lsame(jobu, "N");
  bool wntva = lsame(jobvt, "A");
  bool wntvs = lsame(jobvt, "S");
  bool wntvas = wntva || wntvs;
  bool wntvo = lsame(jobvt, "O");
  bool wntvn = lsame(jobvt, "N");
  bool lquery = (lwork ==  - 1);
  //C
  if (!(wntua || wntus || wntuo || wntun)) {
    info = -1;
  }
  else if (!(wntva || wntvs || wntvo || wntvn) || (wntvo && wntuo)) {
    info = -2;
  }
  else if (m < 0) {
    info = -3;
  }
  else if (n < 0) {
    info = -4;
  }
  else if (lda < fem::max(1, m)) {
    info = -6;
  }
  else if (ldu < 1 || (wntuas && ldu < m)) {
    info = -9;
  }
  else if (ldvt < 1 || (wntva && ldvt < n) || (wntvs && ldvt < minmn)) {
    info = -11;
  }
  //C
  //C     Compute workspace
  //C      (Note: Comments in the code beginning "Workspace:" describe the
  //C       minimal amount of workspace needed at that point in the code,
  //C       as well as the preferred amount for good performance.
  //C       NB refers to the optimal block size for the immediately
  //C       following subroutine, as returned by ILAENV.)
  //C
  int minwrk = fem::int0;
  int maxwrk = fem::int0;
  int mnthr = fem::int0;
  int bdspac = fem::int0;
  int wrkbl = fem::int0;
  if (info == 0) {
    minwrk = 1;
    maxwrk = 1;
    if (m >= n && minmn > 0) {
      //C
      //C           Compute space needed for DBDSQR
      //C
      mnthr = ilaenv(6, "DGESVD", jobu + jobvt, m, n, 0, 0);
      bdspac = 5 * n;
      if (m >= mnthr) {
        if (wntun) {
          //C
          //C                 Path 1 (M much larger than N, JOBU='N')
          //C
          maxwrk = n + n * ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
          maxwrk = fem::max(maxwrk, 3 * n + 2 * n * ilaenv(1,
            "DGEBRD", " ", n, n, -1, -1));
          if (wntvo || wntvas) {
            maxwrk = fem::max(maxwrk, 3 * n + (n - 1) * ilaenv(1,
              "DORGBR", "P", n, n, n, -1));
          }
          maxwrk = fem::max(maxwrk, bdspac);
          minwrk = fem::max(4 * n, bdspac);
        }
        else if (wntuo && wntvn) {
          //C
          //C                 Path 2 (M much larger than N, JOBU='O', JOBVT='N')
          //C
          wrkbl = n + n * ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, n + n * ilaenv(1, "DORGQR", " ", m,
            n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + 2 * n * ilaenv(1, "DGEBRD",
            " ", n, n, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORGBR",
            "Q", n, n, n, -1));
          wrkbl = fem::max(wrkbl, bdspac);
          maxwrk = fem::max(n * n + wrkbl, n * n + m * n + n);
          minwrk = fem::max(3 * n + m, bdspac);
        }
        else if (wntuo && wntvas) {
          //C
          //C                 Path 3 (M much larger than N, JOBU='O', JOBVT='S' or
          //C                 'A')
          //C
          wrkbl = n + n * ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, n + n * ilaenv(1, "DORGQR", " ", m,
            n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + 2 * n * ilaenv(1, "DGEBRD",
            " ", n, n, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORGBR",
            "Q", n, n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + (n - 1) * ilaenv(1,
            "DORGBR", "P", n, n, n, -1));
          wrkbl = fem::max(wrkbl, bdspac);
          maxwrk = fem::max(n * n + wrkbl, n * n + m * n + n);
          minwrk = fem::max(3 * n + m, bdspac);
        }
        else if (wntus && wntvn) {
          //C
          //C                 Path 4 (M much larger than N, JOBU='S', JOBVT='N')
          //C
          wrkbl = n + n * ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, n + n * ilaenv(1, "DORGQR", " ", m,
            n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + 2 * n * ilaenv(1, "DGEBRD",
            " ", n, n, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORGBR",
            "Q", n, n, n, -1));
          wrkbl = fem::max(wrkbl, bdspac);
          maxwrk = n * n + wrkbl;
          minwrk = fem::max(3 * n + m, bdspac);
        }
        else if (wntus && wntvo) {
          //C
          //C                 Path 5 (M much larger than N, JOBU='S', JOBVT='O')
          //C
          wrkbl = n + n * ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, n + n * ilaenv(1, "DORGQR", " ", m,
            n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + 2 * n * ilaenv(1, "DGEBRD",
            " ", n, n, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORGBR",
            "Q", n, n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + (n - 1) * ilaenv(1,
            "DORGBR", "P", n, n, n, -1));
          wrkbl = fem::max(wrkbl, bdspac);
          maxwrk = 2 * n * n + wrkbl;
          minwrk = fem::max(3 * n + m, bdspac);
        }
        else if (wntus && wntvas) {
          //C
          //C                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' or
          //C                 'A')
          //C
          wrkbl = n + n * ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, n + n * ilaenv(1, "DORGQR", " ", m,
            n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + 2 * n * ilaenv(1, "DGEBRD",
            " ", n, n, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORGBR",
            "Q", n, n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + (n - 1) * ilaenv(1,
            "DORGBR", "P", n, n, n, -1));
          wrkbl = fem::max(wrkbl, bdspac);
          maxwrk = n * n + wrkbl;
          minwrk = fem::max(3 * n + m, bdspac);
        }
        else if (wntua && wntvn) {
          //C
          //C                 Path 7 (M much larger than N, JOBU='A', JOBVT='N')
          //C
          wrkbl = n + n * ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, n + m * ilaenv(1, "DORGQR", " ", m,
            m, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + 2 * n * ilaenv(1, "DGEBRD",
            " ", n, n, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORGBR",
            "Q", n, n, n, -1));
          wrkbl = fem::max(wrkbl, bdspac);
          maxwrk = n * n + wrkbl;
          minwrk = fem::max(3 * n + m, bdspac);
        }
        else if (wntua && wntvo) {
          //C
          //C                 Path 8 (M much larger than N, JOBU='A', JOBVT='O')
          //C
          wrkbl = n + n * ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, n + m * ilaenv(1, "DORGQR", " ", m,
            m, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + 2 * n * ilaenv(1, "DGEBRD",
            " ", n, n, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORGBR",
            "Q", n, n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + (n - 1) * ilaenv(1,
            "DORGBR", "P", n, n, n, -1));
          wrkbl = fem::max(wrkbl, bdspac);
          maxwrk = 2 * n * n + wrkbl;
          minwrk = fem::max(3 * n + m, bdspac);
        }
        else if (wntua && wntvas) {
          //C
          //C                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' or
          //C                 'A')
          //C
          wrkbl = n + n * ilaenv(1, "DGEQRF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, n + m * ilaenv(1, "DORGQR", " ", m,
            m, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + 2 * n * ilaenv(1, "DGEBRD",
            " ", n, n, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * n + n * ilaenv(1, "DORGBR",
            "Q", n, n, n, -1));
          wrkbl = fem::max(wrkbl, 3 * n + (n - 1) * ilaenv(1,
            "DORGBR", "P", n, n, n, -1));
          wrkbl = fem::max(wrkbl, bdspac);
          maxwrk = n * n + wrkbl;
          minwrk = fem::max(3 * n + m, bdspac);
        }
      }
      else {
        //C
        //C              Path 10 (M at least N, but not much larger)
        //C
        maxwrk = 3 * n + (m + n) * ilaenv(1, "DGEBRD", " ", m, n, -1, -1);
        if (wntus || wntuo) {
          maxwrk = fem::max(maxwrk, 3 * n + n * ilaenv(1, "DORGBR",
            "Q", m, n, n, -1));
        }
        if (wntua) {
          maxwrk = fem::max(maxwrk, 3 * n + m * ilaenv(1, "DORGBR",
            "Q", m, m, n, -1));
        }
        if (!wntvn) {
          maxwrk = fem::max(maxwrk, 3 * n + (n - 1) * ilaenv(1,
            "DORGBR", "P", n, n, n, -1));
        }
        maxwrk = fem::max(maxwrk, bdspac);
        minwrk = fem::max(3 * n + m, bdspac);
      }
    }
    else if (minmn > 0) {
      //C
      //C           Compute space needed for DBDSQR
      //C
      mnthr = ilaenv(6, "DGESVD", jobu + jobvt, m, n, 0, 0);
      bdspac = 5 * m;
      if (n >= mnthr) {
        if (wntvn) {
          //C
          //C                 Path 1t(N much larger than M, JOBVT='N')
          //C
          maxwrk = m + m * ilaenv(1, "DGELQF", " ", m, n, -1, -1);
          maxwrk = fem::max(maxwrk, 3 * m + 2 * m * ilaenv(1,
            "DGEBRD", " ", m, m, -1, -1));
          if (wntuo || wntuas) {
            maxwrk = fem::max(maxwrk, 3 * m + m * ilaenv(1, "DORGBR",
              "Q", m, m, m, -1));
          }
          maxwrk = fem::max(maxwrk, bdspac);
          minwrk = fem::max(4 * m, bdspac);
        }
        else if (wntvo && wntun) {
          //C
          //C                 Path 2t(N much larger than M, JOBU='N', JOBVT='O')
          //C
          wrkbl = m + m * ilaenv(1, "DGELQF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, m + m * ilaenv(1, "DORGLQ", " ", m,
            n, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + 2 * m * ilaenv(1, "DGEBRD",
            " ", m, m, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * m + (m - 1) * ilaenv(1,
            "DORGBR", "P", m, m, m, -1));
          wrkbl = fem::max(wrkbl, bdspac);
          maxwrk = fem::max(m * m + wrkbl, m * m + m * n + m);
          minwrk = fem::max(3 * m + n, bdspac);
        }
        else if (wntvo && wntuas) {
          //C
          //C                 Path 3t(N much larger than M, JOBU='S' or 'A',
          //C                 JOBVT='O')
          //C
          wrkbl = m + m * ilaenv(1, "DGELQF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, m + m * ilaenv(1, "DORGLQ", " ", m,
            n, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + 2 * m * ilaenv(1, "DGEBRD",
            " ", m, m, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * m + (m - 1) * ilaenv(1,
            "DORGBR", "P", m, m, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORGBR",
            "Q", m, m, m, -1));
          wrkbl = fem::max(wrkbl, bdspac);
          maxwrk = fem::max(m * m + wrkbl, m * m + m * n + m);
          minwrk = fem::max(3 * m + n, bdspac);
        }
        else if (wntvs && wntun) {
          //C
          //C                 Path 4t(N much larger than M, JOBU='N', JOBVT='S')
          //C
          wrkbl = m + m * ilaenv(1, "DGELQF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, m + m * ilaenv(1, "DORGLQ", " ", m,
            n, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + 2 * m * ilaenv(1, "DGEBRD",
            " ", m, m, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * m + (m - 1) * ilaenv(1,
            "DORGBR", "P", m, m, m, -1));
          wrkbl = fem::max(wrkbl, bdspac);
          maxwrk = m * m + wrkbl;
          minwrk = fem::max(3 * m + n, bdspac);
        }
        else if (wntvs && wntuo) {
          //C
          //C                 Path 5t(N much larger than M, JOBU='O', JOBVT='S')
          //C
          wrkbl = m + m * ilaenv(1, "DGELQF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, m + m * ilaenv(1, "DORGLQ", " ", m,
            n, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + 2 * m * ilaenv(1, "DGEBRD",
            " ", m, m, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * m + (m - 1) * ilaenv(1,
            "DORGBR", "P", m, m, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORGBR",
            "Q", m, m, m, -1));
          wrkbl = fem::max(wrkbl, bdspac);
          maxwrk = 2 * m * m + wrkbl;
          minwrk = fem::max(3 * m + n, bdspac);
        }
        else if (wntvs && wntuas) {
          //C
          //C                 Path 6t(N much larger than M, JOBU='S' or 'A',
          //C                 JOBVT='S')
          //C
          wrkbl = m + m * ilaenv(1, "DGELQF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, m + m * ilaenv(1, "DORGLQ", " ", m,
            n, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + 2 * m * ilaenv(1, "DGEBRD",
            " ", m, m, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * m + (m - 1) * ilaenv(1,
            "DORGBR", "P", m, m, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORGBR",
            "Q", m, m, m, -1));
          wrkbl = fem::max(wrkbl, bdspac);
          maxwrk = m * m + wrkbl;
          minwrk = fem::max(3 * m + n, bdspac);
        }
        else if (wntva && wntun) {
          //C
          //C                 Path 7t(N much larger than M, JOBU='N', JOBVT='A')
          //C
          wrkbl = m + m * ilaenv(1, "DGELQF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, m + n * ilaenv(1, "DORGLQ", " ", n,
            n, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + 2 * m * ilaenv(1, "DGEBRD",
            " ", m, m, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * m + (m - 1) * ilaenv(1,
            "DORGBR", "P", m, m, m, -1));
          wrkbl = fem::max(wrkbl, bdspac);
          maxwrk = m * m + wrkbl;
          minwrk = fem::max(3 * m + n, bdspac);
        }
        else if (wntva && wntuo) {
          //C
          //C                 Path 8t(N much larger than M, JOBU='O', JOBVT='A')
          //C
          wrkbl = m + m * ilaenv(1, "DGELQF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, m + n * ilaenv(1, "DORGLQ", " ", n,
            n, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + 2 * m * ilaenv(1, "DGEBRD",
            " ", m, m, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * m + (m - 1) * ilaenv(1,
            "DORGBR", "P", m, m, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORGBR",
            "Q", m, m, m, -1));
          wrkbl = fem::max(wrkbl, bdspac);
          maxwrk = 2 * m * m + wrkbl;
          minwrk = fem::max(3 * m + n, bdspac);
        }
        else if (wntva && wntuas) {
          //C
          //C                 Path 9t(N much larger than M, JOBU='S' or 'A',
          //C                 JOBVT='A')
          //C
          wrkbl = m + m * ilaenv(1, "DGELQF", " ", m, n, -1, -1);
          wrkbl = fem::max(wrkbl, m + n * ilaenv(1, "DORGLQ", " ", n,
            n, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + 2 * m * ilaenv(1, "DGEBRD",
            " ", m, m, -1, -1));
          wrkbl = fem::max(wrkbl, 3 * m + (m - 1) * ilaenv(1,
            "DORGBR", "P", m, m, m, -1));
          wrkbl = fem::max(wrkbl, 3 * m + m * ilaenv(1, "DORGBR",
            "Q", m, m, m, -1));
          wrkbl = fem::max(wrkbl, bdspac);
          maxwrk = m * m + wrkbl;
          minwrk = fem::max(3 * m + n, bdspac);
        }
      }
      else {
        //C
        //C              Path 10t(N greater than M, but not much larger)
        //C
        maxwrk = 3 * m + (m + n) * ilaenv(1, "DGEBRD", " ", m, n, -1, -1);
        if (wntvs || wntvo) {
          maxwrk = fem::max(maxwrk, 3 * m + m * ilaenv(1, "DORGBR",
            "P", m, n, m, -1));
        }
        if (wntva) {
          maxwrk = fem::max(maxwrk, 3 * m + n * ilaenv(1, "DORGBR",
            "P", n, n, m, -1));
        }
        if (!wntun) {
          maxwrk = fem::max(maxwrk, 3 * m + (m - 1) * ilaenv(1,
            "DORGBR", "Q", m, m, m, -1));
        }
        maxwrk = fem::max(maxwrk, bdspac);
        minwrk = fem::max(3 * m + n, bdspac);
      }
    }
    maxwrk = fem::max(maxwrk, minwrk);
    work(1) = maxwrk;
    //C
    if (lwork < minwrk && !lquery) {
      info = -13;
    }
  }
  //C
  if (info != 0) {
    xerbla("DGESVD", -info);
    return;
  }
  else if (lquery) {
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (m == 0 || n == 0) {
    return;
  }
  //C
  //C     Get machine constants
  //C
  double eps = dlamch(cmn, "P");
  double smlnum = fem::sqrt(dlamch(cmn, "S")) / eps;
  const double one = 1.0e0;
  double bignum = one / smlnum;
  //C
  //C     Scale A if max element outside range [SMLNUM,BIGNUM]
  //C
  arr_1d<1, double> dum(fem::fill0);
  double anrm = dlange("M", m, n, a, lda, dum);
  int iscl = 0;
  const double zero = 0.0e0;
  int ierr = fem::int0;
  if (anrm > zero && anrm < smlnum) {
    iscl = 1;
    dlascl(cmn, "G", 0, 0, anrm, smlnum, m, n, a, lda, ierr);
  }
  else if (anrm > bignum) {
    iscl = 1;
    dlascl(cmn, "G", 0, 0, anrm, bignum, m, n, a, lda, ierr);
  }
  //C
  int itau = fem::int0;
  int iwork = fem::int0;
  int ie = fem::int0;
  int itauq = fem::int0;
  int itaup = fem::int0;
  int ncvt = fem::int0;
  int ir = fem::int0;
  int ldwrku = fem::int0;
  int ldwrkr = fem::int0;
  int iu = fem::int0;
  int i = fem::int0;
  int chunk = fem::int0;
  int ncu = fem::int0;
  int nru = fem::int0;
  int blk = fem::int0;
  int nrvt = fem::int0;
  if (m >= n) {
    //C
    //C        A has at least as many rows as columns. If A has sufficiently
    //C        more rows than columns, first reduce using the QR
    //C        decomposition (if sufficient workspace available)
    //C
    if (m >= mnthr) {
      //C
      if (wntun) {
        //C
        //C              Path 1 (M much larger than N, JOBU='N')
        //C              No left singular vectors to be computed
        //C
        itau = 1;
        iwork = itau + n;
        //C
        //C              Compute A=Q*R
        //C              (Workspace: need 2*N, prefer N+N*NB)
        //C
        dgeqrf(cmn, m, n, a, lda, work(itau), work(iwork), lwork - iwork + 1,
          ierr);
        //C
        //C              Zero out below R
        //C
        dlaset("L", n - 1, n - 1, zero, zero, a(2, 1), lda);
        ie = 1;
        itauq = ie + n;
        itaup = itauq + n;
        iwork = itaup + n;
        //C
        //C              Bidiagonalize R in A
        //C              (Workspace: need 4*N, prefer 3*N+2*N*NB)
        //C
        dgebrd(cmn, n, n, a, lda, s, work(ie), work(itauq), work(itaup),
          work(iwork), lwork - iwork + 1, ierr);
        ncvt = 0;
        if (wntvo || wntvas) {
          //C
          //C                 If right singular vectors desired, generate P'.
          //C                 (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
          //C
          dorgbr("P", n, n, n, a, lda, work(itaup), work(iwork),
            lwork - iwork + 1, ierr);
          ncvt = n;
        }
        iwork = ie + n;
        //C
        //C              Perform bidiagonal QR iteration, computing right
        //C              singular vectors of A in A if desired
        //C              (Workspace: need BDSPAC)
        //C
        dbdsqr(cmn, "U", n, ncvt, 0, 0, s, work(ie), a, lda, dum, 1,
          dum, 1, work(iwork), info);
        //C
        //C              If right singular vectors desired in VT, copy them there
        //C
        if (wntvas) {
          dlacpy("F", n, n, a, lda, vt, ldvt);
        }
        //C
      }
      else if (wntuo && wntvn) {
        //C
        //C              Path 2 (M much larger than N, JOBU='O', JOBVT='N')
        //C              N left singular vectors to be overwritten on A and
        //C              no right singular vectors to be computed
        //C
        if (lwork >= n * n + fem::max(4 * n, bdspac)) {
          //C
          //C                 Sufficient workspace for a fast algorithm
          //C
          ir = 1;
          if (lwork >= fem::max(wrkbl, lda * n + n) + lda * n) {
            //C
            //C                    WORK(IU) is LDA by N, WORK(IR) is LDA by N
            //C
            ldwrku = lda;
            ldwrkr = lda;
          }
          else if (lwork >= fem::max(wrkbl, lda * n + n) + n * n) {
            //C
            //C                    WORK(IU) is LDA by N, WORK(IR) is N by N
            //C
            ldwrku = lda;
            ldwrkr = n;
          }
          else {
            //C
            //C                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N
            //C
            ldwrku = (lwork - n * n - n) / n;
            ldwrkr = n;
          }
          itau = ir + ldwrkr * n;
          iwork = itau + n;
          //C
          //C                 Compute A=Q*R
          //C                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
          //C
          dgeqrf(cmn, m, n, a, lda, work(itau), work(iwork), lwork - iwork + 1,
            ierr);
          //C
          //C                 Copy R to WORK(IR) and zero out below it
          //C
          dlacpy("U", n, n, a, lda, work(ir), ldwrkr);
          dlaset("L", n - 1, n - 1, zero, zero, work(ir + 1), ldwrkr);
          //C
          //C                 Generate Q in A
          //C                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
          //C
          dorgqr(m, n, n, a, lda, work(itau), work(iwork), lwork - iwork + 1,
            ierr);
          ie = itau;
          itauq = ie + n;
          itaup = itauq + n;
          iwork = itaup + n;
          //C
          //C                 Bidiagonalize R in WORK(IR)
          //C                 (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
          //C
          dgebrd(cmn, n, n, work(ir), ldwrkr, s, work(ie), work(itauq),
            work(itaup), work(iwork), lwork - iwork + 1, ierr);
          //C
          //C                 Generate left vectors bidiagonalizing R
          //C                 (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)
          //C
          dorgbr("Q", n, n, n, work(ir), ldwrkr, work(itauq), work(iwork),
            lwork - iwork + 1, ierr);
          iwork = ie + n;
          //C
          //C                 Perform bidiagonal QR iteration, computing left
          //C                 singular vectors of R in WORK(IR)
          //C                 (Workspace: need N*N+BDSPAC)
          //C
          dbdsqr(cmn, "U", n, 0, n, 0, s, work(ie), dum, 1, work(ir),
            ldwrkr, dum, 1, work(iwork), info);
          iu = ie + n;
          //C
          //C                 Multiply Q in A by left singular vectors of R in
          //C                 WORK(IR), storing result in WORK(IU) and copying to A
          //C                 (Workspace: need N*N+2*N, prefer N*N+M*N+N)
          //C
          FEM_DOSTEP(i, 1, m, ldwrku) {
            chunk = fem::min(m - i + 1, ldwrku);
            dgemm("N", "N", chunk, n, n, one, a(i, 1), lda, work(ir),
              ldwrkr, zero, work(iu), ldwrku);
            dlacpy("F", chunk, n, work(iu), ldwrku, a(i, 1), lda);
          }
          //C
        }
        else {
          //C
          //C                 Insufficient workspace for a fast algorithm
          //C
          ie = 1;
          itauq = ie + n;
          itaup = itauq + n;
          iwork = itaup + n;
          //C
          //C                 Bidiagonalize A
          //C                 (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB)
          //C
          dgebrd(cmn, m, n, a, lda, s, work(ie), work(itauq), work(itaup),
            work(iwork), lwork - iwork + 1, ierr);
          //C
          //C                 Generate left vectors bidiagonalizing A
          //C                 (Workspace: need 4*N, prefer 3*N+N*NB)
          //C
          dorgbr("Q", m, n, n, a, lda, work(itauq), work(iwork),
            lwork - iwork + 1, ierr);
          iwork = ie + n;
          //C
          //C                 Perform bidiagonal QR iteration, computing left
          //C                 singular vectors of A in A
          //C                 (Workspace: need BDSPAC)
          //C
          dbdsqr(cmn, "U", n, 0, m, 0, s, work(ie), dum, 1, a, lda,
            dum, 1, work(iwork), info);
          //C
        }
        //C
      }
      else if (wntuo && wntvas) {
        //C
        //C              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A')
        //C              N left singular vectors to be overwritten on A and
        //C              N right singular vectors to be computed in VT
        //C
        if (lwork >= n * n + fem::max(4 * n, bdspac)) {
          //C
          //C                 Sufficient workspace for a fast algorithm
          //C
          ir = 1;
          if (lwork >= fem::max(wrkbl, lda * n + n) + lda * n) {
            //C
            //C                    WORK(IU) is LDA by N and WORK(IR) is LDA by N
            //C
            ldwrku = lda;
            ldwrkr = lda;
          }
          else if (lwork >= fem::max(wrkbl, lda * n + n) + n * n) {
            //C
            //C                    WORK(IU) is LDA by N and WORK(IR) is N by N
            //C
            ldwrku = lda;
            ldwrkr = n;
          }
          else {
            //C
            //C                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N
            //C
            ldwrku = (lwork - n * n - n) / n;
            ldwrkr = n;
          }
          itau = ir + ldwrkr * n;
          iwork = itau + n;
          //C
          //C                 Compute A=Q*R
          //C                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
          //C
          dgeqrf(cmn, m, n, a, lda, work(itau), work(iwork), lwork - iwork + 1,
            ierr);
          //C
          //C                 Copy R to VT, zeroing out below it
          //C
          dlacpy("U", n, n, a, lda, vt, ldvt);
          if (n > 1) {
            dlaset("L", n - 1, n - 1, zero, zero, vt(2, 1), ldvt);
          }
          //C
          //C                 Generate Q in A
          //C                 (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
          //C
          dorgqr(m, n, n, a, lda, work(itau), work(iwork), lwork - iwork + 1,
            ierr);
          ie = itau;
          itauq = ie + n;
          itaup = itauq + n;
          iwork = itaup + n;
          //C
          //C                 Bidiagonalize R in VT, copying result to WORK(IR)
          //C                 (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
          //C
          dgebrd(cmn, n, n, vt, ldvt, s, work(ie), work(itauq), work(itaup),
            work(iwork), lwork - iwork + 1, ierr);
          dlacpy("L", n, n, vt, ldvt, work(ir), ldwrkr);
          //C
          //C                 Generate left vectors bidiagonalizing R in WORK(IR)
          //C                 (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)
          //C
          dorgbr("Q", n, n, n, work(ir), ldwrkr, work(itauq), work(iwork),
            lwork - iwork + 1, ierr);
          //C
          //C                 Generate right vectors bidiagonalizing R in VT
          //C                 (Workspace: need N*N+4*N-1, prefer N*N+3*N+(N-1)*NB)
          //C
          dorgbr("P", n, n, n, vt, ldvt, work(itaup), work(iwork),
            lwork - iwork + 1, ierr);
          iwork = ie + n;
          //C
          //C                 Perform bidiagonal QR iteration, computing left
          //C                 singular vectors of R in WORK(IR) and computing right
          //C                 singular vectors of R in VT
          //C                 (Workspace: need N*N+BDSPAC)
          //C
          dbdsqr(cmn, "U", n, n, n, 0, s, work(ie), vt, ldvt, work(ir),
            ldwrkr, dum, 1, work(iwork), info);
          iu = ie + n;
          //C
          //C                 Multiply Q in A by left singular vectors of R in
          //C                 WORK(IR), storing result in WORK(IU) and copying to A
          //C                 (Workspace: need N*N+2*N, prefer N*N+M*N+N)
          //C
          FEM_DOSTEP(i, 1, m, ldwrku) {
            chunk = fem::min(m - i + 1, ldwrku);
            dgemm("N", "N", chunk, n, n, one, a(i, 1), lda, work(ir),
              ldwrkr, zero, work(iu), ldwrku);
            dlacpy("F", chunk, n, work(iu), ldwrku, a(i, 1), lda);
          }
          //C
        }
        else {
          //C
          //C                 Insufficient workspace for a fast algorithm
          //C
          itau = 1;
          iwork = itau + n;
          //C
          //C                 Compute A=Q*R
          //C                 (Workspace: need 2*N, prefer N+N*NB)
          //C
          dgeqrf(cmn, m, n, a, lda, work(itau), work(iwork), lwork - iwork + 1,
            ierr);
          //C
          //C                 Copy R to VT, zeroing out below it
          //C
          dlacpy("U", n, n, a, lda, vt, ldvt);
          if (n > 1) {
            dlaset("L", n - 1, n - 1, zero, zero, vt(2, 1), ldvt);
          }
          //C
          //C                 Generate Q in A
          //C                 (Workspace: need 2*N, prefer N+N*NB)
          //C
          dorgqr(m, n, n, a, lda, work(itau), work(iwork), lwork - iwork + 1,
            ierr);
          ie = itau;
          itauq = ie + n;
          itaup = itauq + n;
          iwork = itaup + n;
          //C
          //C                 Bidiagonalize R in VT
          //C                 (Workspace: need 4*N, prefer 3*N+2*N*NB)
          //C
          dgebrd(cmn, n, n, vt, ldvt, s, work(ie), work(itauq), work(itaup),
            work(iwork), lwork - iwork + 1, ierr);
          //C
          //C                 Multiply Q in A by left vectors bidiagonalizing R
          //C                 (Workspace: need 3*N+M, prefer 3*N+M*NB)
          //C
          dormbr("Q", "R", "N", m, n, n, vt, ldvt, work(itauq), a,
            lda, work(iwork), lwork - iwork + 1, ierr);
          //C
          //C                 Generate right vectors bidiagonalizing R in VT
          //C                 (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
          //C
          dorgbr("P", n, n, n, vt, ldvt, work(itaup), work(iwork),
            lwork - iwork + 1, ierr);
          iwork = ie + n;
          //C
          //C                 Perform bidiagonal QR iteration, computing left
          //C                 singular vectors of A in A and computing right
          //C                 singular vectors of A in VT
          //C                 (Workspace: need BDSPAC)
          //C
          dbdsqr(cmn, "U", n, n, m, 0, s, work(ie), vt, ldvt, a, lda,
            dum, 1, work(iwork), info);
          //C
        }
        //C
      }
      else if (wntus) {
        //C
        if (wntvn) {
          //C
          //C                 Path 4 (M much larger than N, JOBU='S', JOBVT='N')
          //C                 N left singular vectors to be computed in U and
          //C                 no right singular vectors to be computed
          //C
          if (lwork >= n * n + fem::max(4 * n, bdspac)) {
            //C
            //C                    Sufficient workspace for a fast algorithm
            //C
            ir = 1;
            if (lwork >= wrkbl + lda * n) {
              //C
              //C                       WORK(IR) is LDA by N
              //C
              ldwrkr = lda;
            }
            else {
              //C
              //C                       WORK(IR) is N by N
              //C
              ldwrkr = n;
            }
            itau = ir + ldwrkr * n;
            iwork = itau + n;
            //C
            //C                    Compute A=Q*R
            //C                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
            //C
            dgeqrf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Copy R to WORK(IR), zeroing out below it
            //C
            dlacpy("U", n, n, a, lda, work(ir), ldwrkr);
            dlaset("L", n - 1, n - 1, zero, zero, work(ir + 1), ldwrkr);
            //C
            //C                    Generate Q in A
            //C                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
            //C
            dorgqr(m, n, n, a, lda, work(itau), work(iwork), lwork - iwork + 1,
              ierr);
            ie = itau;
            itauq = ie + n;
            itaup = itauq + n;
            iwork = itaup + n;
            //C
            //C                    Bidiagonalize R in WORK(IR)
            //C                    (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
            //C
            dgebrd(cmn, n, n, work(ir), ldwrkr, s, work(ie), work(itauq),
              work(itaup), work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Generate left vectors bidiagonalizing R in WORK(IR)
            //C                    (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)
            //C
            dorgbr("Q", n, n, n, work(ir), ldwrkr, work(itauq), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + n;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of R in WORK(IR)
            //C                    (Workspace: need N*N+BDSPAC)
            //C
            dbdsqr(cmn, "U", n, 0, n, 0, s, work(ie), dum, 1, work(ir),
              ldwrkr, dum, 1, work(iwork), info);
            //C
            //C                    Multiply Q in A by left singular vectors of R in
            //C                    WORK(IR), storing result in U
            //C                    (Workspace: need N*N)
            //C
            dgemm("N", "N", m, n, n, one, a, lda, work(ir), ldwrkr,
              zero, u, ldu);
            //C
          }
          else {
            //C
            //C                    Insufficient workspace for a fast algorithm
            //C
            itau = 1;
            iwork = itau + n;
            //C
            //C                    Compute A=Q*R, copying result to U
            //C                    (Workspace: need 2*N, prefer N+N*NB)
            //C
            dgeqrf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("L", m, n, a, lda, u, ldu);
            //C
            //C                    Generate Q in U
            //C                    (Workspace: need 2*N, prefer N+N*NB)
            //C
            dorgqr(m, n, n, u, ldu, work(itau), work(iwork), lwork - iwork + 1,
              ierr);
            ie = itau;
            itauq = ie + n;
            itaup = itauq + n;
            iwork = itaup + n;
            //C
            //C                    Zero out below R in A
            //C
            dlaset("L", n - 1, n - 1, zero, zero, a(2, 1), lda);
            //C
            //C                    Bidiagonalize R in A
            //C                    (Workspace: need 4*N, prefer 3*N+2*N*NB)
            //C
            dgebrd(cmn, n, n, a, lda, s, work(ie), work(itauq), work(itaup),
              work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Multiply Q in U by left vectors bidiagonalizing R
            //C                    (Workspace: need 3*N+M, prefer 3*N+M*NB)
            //C
            dormbr("Q", "R", "N", m, n, n, a, lda, work(itauq), u,
              ldu, work(iwork), lwork - iwork + 1, ierr);
            iwork = ie + n;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of A in U
            //C                    (Workspace: need BDSPAC)
            //C
            dbdsqr(cmn, "U", n, 0, m, 0, s, work(ie), dum, 1, u, ldu,
              dum, 1, work(iwork), info);
            //C
          }
          //C
        }
        else if (wntvo) {
          //C
          //C                 Path 5 (M much larger than N, JOBU='S', JOBVT='O')
          //C                 N left singular vectors to be computed in U and
          //C                 N right singular vectors to be overwritten on A
          //C
          if (lwork >= 2 * n * n + fem::max(4 * n, bdspac)) {
            //C
            //C                    Sufficient workspace for a fast algorithm
            //C
            iu = 1;
            if (lwork >= wrkbl + 2 * lda * n) {
              //C
              //C                       WORK(IU) is LDA by N and WORK(IR) is LDA by N
              //C
              ldwrku = lda;
              ir = iu + ldwrku * n;
              ldwrkr = lda;
            }
            else if (lwork >= wrkbl + (lda + n) * n) {
              //C
              //C                       WORK(IU) is LDA by N and WORK(IR) is N by N
              //C
              ldwrku = lda;
              ir = iu + ldwrku * n;
              ldwrkr = n;
            }
            else {
              //C
              //C                       WORK(IU) is N by N and WORK(IR) is N by N
              //C
              ldwrku = n;
              ir = iu + ldwrku * n;
              ldwrkr = n;
            }
            itau = ir + ldwrkr * n;
            iwork = itau + n;
            //C
            //C                    Compute A=Q*R
            //C                    (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
            //C
            dgeqrf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Copy R to WORK(IU), zeroing out below it
            //C
            dlacpy("U", n, n, a, lda, work(iu), ldwrku);
            dlaset("L", n - 1, n - 1, zero, zero, work(iu + 1), ldwrku);
            //C
            //C                    Generate Q in A
            //C                    (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
            //C
            dorgqr(m, n, n, a, lda, work(itau), work(iwork), lwork - iwork + 1,
              ierr);
            ie = itau;
            itauq = ie + n;
            itaup = itauq + n;
            iwork = itaup + n;
            //C
            //C                    Bidiagonalize R in WORK(IU), copying result to
            //C                    WORK(IR)
            //C                    (Workspace: need 2*N*N+4*N,
            //C                                prefer 2*N*N+3*N+2*N*NB)
            //C
            dgebrd(cmn, n, n, work(iu), ldwrku, s, work(ie), work(itauq),
              work(itaup), work(iwork), lwork - iwork + 1, ierr);
            dlacpy("U", n, n, work(iu), ldwrku, work(ir), ldwrkr);
            //C
            //C                    Generate left bidiagonalizing vectors in WORK(IU)
            //C                    (Workspace: need 2*N*N+4*N, prefer 2*N*N+3*N+N*NB)
            //C
            dorgbr("Q", n, n, n, work(iu), ldwrku, work(itauq), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Generate right bidiagonalizing vectors in WORK(IR)
            //C                    (Workspace: need 2*N*N+4*N-1,
            //C                                prefer 2*N*N+3*N+(N-1)*NB)
            //C
            dorgbr("P", n, n, n, work(ir), ldwrkr, work(itaup), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + n;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of R in WORK(IU) and computing
            //C                    right singular vectors of R in WORK(IR)
            //C                    (Workspace: need 2*N*N+BDSPAC)
            //C
            dbdsqr(cmn, "U", n, n, n, 0, s, work(ie), work(ir),
              ldwrkr, work(iu), ldwrku, dum, 1, work(iwork), info);
            //C
            //C                    Multiply Q in A by left singular vectors of R in
            //C                    WORK(IU), storing result in U
            //C                    (Workspace: need N*N)
            //C
            dgemm("N", "N", m, n, n, one, a, lda, work(iu), ldwrku,
              zero, u, ldu);
            //C
            //C                    Copy right singular vectors of R to A
            //C                    (Workspace: need N*N)
            //C
            dlacpy("F", n, n, work(ir), ldwrkr, a, lda);
            //C
          }
          else {
            //C
            //C                    Insufficient workspace for a fast algorithm
            //C
            itau = 1;
            iwork = itau + n;
            //C
            //C                    Compute A=Q*R, copying result to U
            //C                    (Workspace: need 2*N, prefer N+N*NB)
            //C
            dgeqrf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("L", m, n, a, lda, u, ldu);
            //C
            //C                    Generate Q in U
            //C                    (Workspace: need 2*N, prefer N+N*NB)
            //C
            dorgqr(m, n, n, u, ldu, work(itau), work(iwork), lwork - iwork + 1,
              ierr);
            ie = itau;
            itauq = ie + n;
            itaup = itauq + n;
            iwork = itaup + n;
            //C
            //C                    Zero out below R in A
            //C
            dlaset("L", n - 1, n - 1, zero, zero, a(2, 1), lda);
            //C
            //C                    Bidiagonalize R in A
            //C                    (Workspace: need 4*N, prefer 3*N+2*N*NB)
            //C
            dgebrd(cmn, n, n, a, lda, s, work(ie), work(itauq), work(itaup),
              work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Multiply Q in U by left vectors bidiagonalizing R
            //C                    (Workspace: need 3*N+M, prefer 3*N+M*NB)
            //C
            dormbr("Q", "R", "N", m, n, n, a, lda, work(itauq), u,
              ldu, work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Generate right vectors bidiagonalizing R in A
            //C                    (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
            //C
            dorgbr("P", n, n, n, a, lda, work(itaup), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + n;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of A in U and computing right
            //C                    singular vectors of A in A
            //C                    (Workspace: need BDSPAC)
            //C
            dbdsqr(cmn, "U", n, n, m, 0, s, work(ie), a, lda, u, ldu,
              dum, 1, work(iwork), info);
            //C
          }
          //C
        }
        else if (wntvas) {
          //C
          //C                 Path 6 (M much larger than N, JOBU='S', JOBVT='S'
          //C                         or 'A')
          //C                 N left singular vectors to be computed in U and
          //C                 N right singular vectors to be computed in VT
          //C
          if (lwork >= n * n + fem::max(4 * n, bdspac)) {
            //C
            //C                    Sufficient workspace for a fast algorithm
            //C
            iu = 1;
            if (lwork >= wrkbl + lda * n) {
              //C
              //C                       WORK(IU) is LDA by N
              //C
              ldwrku = lda;
            }
            else {
              //C
              //C                       WORK(IU) is N by N
              //C
              ldwrku = n;
            }
            itau = iu + ldwrku * n;
            iwork = itau + n;
            //C
            //C                    Compute A=Q*R
            //C                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
            //C
            dgeqrf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Copy R to WORK(IU), zeroing out below it
            //C
            dlacpy("U", n, n, a, lda, work(iu), ldwrku);
            dlaset("L", n - 1, n - 1, zero, zero, work(iu + 1), ldwrku);
            //C
            //C                    Generate Q in A
            //C                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
            //C
            dorgqr(m, n, n, a, lda, work(itau), work(iwork), lwork - iwork + 1,
              ierr);
            ie = itau;
            itauq = ie + n;
            itaup = itauq + n;
            iwork = itaup + n;
            //C
            //C                    Bidiagonalize R in WORK(IU), copying result to VT
            //C                    (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
            //C
            dgebrd(cmn, n, n, work(iu), ldwrku, s, work(ie), work(itauq),
              work(itaup), work(iwork), lwork - iwork + 1, ierr);
            dlacpy("U", n, n, work(iu), ldwrku, vt, ldvt);
            //C
            //C                    Generate left bidiagonalizing vectors in WORK(IU)
            //C                    (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)
            //C
            dorgbr("Q", n, n, n, work(iu), ldwrku, work(itauq), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Generate right bidiagonalizing vectors in VT
            //C                    (Workspace: need N*N+4*N-1,
            //C                                prefer N*N+3*N+(N-1)*NB)
            //C
            dorgbr("P", n, n, n, vt, ldvt, work(itaup), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + n;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of R in WORK(IU) and computing
            //C                    right singular vectors of R in VT
            //C                    (Workspace: need N*N+BDSPAC)
            //C
            dbdsqr(cmn, "U", n, n, n, 0, s, work(ie), vt, ldvt, work(iu),
              ldwrku, dum, 1, work(iwork), info);
            //C
            //C                    Multiply Q in A by left singular vectors of R in
            //C                    WORK(IU), storing result in U
            //C                    (Workspace: need N*N)
            //C
            dgemm("N", "N", m, n, n, one, a, lda, work(iu), ldwrku,
              zero, u, ldu);
            //C
          }
          else {
            //C
            //C                    Insufficient workspace for a fast algorithm
            //C
            itau = 1;
            iwork = itau + n;
            //C
            //C                    Compute A=Q*R, copying result to U
            //C                    (Workspace: need 2*N, prefer N+N*NB)
            //C
            dgeqrf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("L", m, n, a, lda, u, ldu);
            //C
            //C                    Generate Q in U
            //C                    (Workspace: need 2*N, prefer N+N*NB)
            //C
            dorgqr(m, n, n, u, ldu, work(itau), work(iwork), lwork - iwork + 1,
              ierr);
            //C
            //C                    Copy R to VT, zeroing out below it
            //C
            dlacpy("U", n, n, a, lda, vt, ldvt);
            if (n > 1) {
              dlaset("L", n - 1, n - 1, zero, zero, vt(2, 1), ldvt);
            }
            ie = itau;
            itauq = ie + n;
            itaup = itauq + n;
            iwork = itaup + n;
            //C
            //C                    Bidiagonalize R in VT
            //C                    (Workspace: need 4*N, prefer 3*N+2*N*NB)
            //C
            dgebrd(cmn, n, n, vt, ldvt, s, work(ie), work(itauq),
              work(itaup), work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Multiply Q in U by left bidiagonalizing vectors
            //C                    in VT
            //C                    (Workspace: need 3*N+M, prefer 3*N+M*NB)
            //C
            dormbr("Q", "R", "N", m, n, n, vt, ldvt, work(itauq), u,
              ldu, work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Generate right bidiagonalizing vectors in VT
            //C                    (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
            //C
            dorgbr("P", n, n, n, vt, ldvt, work(itaup), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + n;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of A in U and computing right
            //C                    singular vectors of A in VT
            //C                    (Workspace: need BDSPAC)
            //C
            dbdsqr(cmn, "U", n, n, m, 0, s, work(ie), vt, ldvt, u,
              ldu, dum, 1, work(iwork), info);
            //C
          }
          //C
        }
        //C
      }
      else if (wntua) {
        //C
        if (wntvn) {
          //C
          //C                 Path 7 (M much larger than N, JOBU='A', JOBVT='N')
          //C                 M left singular vectors to be computed in U and
          //C                 no right singular vectors to be computed
          //C
          if (lwork >= n * n + fem::max(n + m, 4 * n, bdspac)) {
            //C
            //C                    Sufficient workspace for a fast algorithm
            //C
            ir = 1;
            if (lwork >= wrkbl + lda * n) {
              //C
              //C                       WORK(IR) is LDA by N
              //C
              ldwrkr = lda;
            }
            else {
              //C
              //C                       WORK(IR) is N by N
              //C
              ldwrkr = n;
            }
            itau = ir + ldwrkr * n;
            iwork = itau + n;
            //C
            //C                    Compute A=Q*R, copying result to U
            //C                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
            //C
            dgeqrf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("L", m, n, a, lda, u, ldu);
            //C
            //C                    Copy R to WORK(IR), zeroing out below it
            //C
            dlacpy("U", n, n, a, lda, work(ir), ldwrkr);
            dlaset("L", n - 1, n - 1, zero, zero, work(ir + 1), ldwrkr);
            //C
            //C                    Generate Q in U
            //C                    (Workspace: need N*N+N+M, prefer N*N+N+M*NB)
            //C
            dorgqr(m, m, n, u, ldu, work(itau), work(iwork), lwork - iwork + 1,
              ierr);
            ie = itau;
            itauq = ie + n;
            itaup = itauq + n;
            iwork = itaup + n;
            //C
            //C                    Bidiagonalize R in WORK(IR)
            //C                    (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
            //C
            dgebrd(cmn, n, n, work(ir), ldwrkr, s, work(ie), work(itauq),
              work(itaup), work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Generate left bidiagonalizing vectors in WORK(IR)
            //C                    (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)
            //C
            dorgbr("Q", n, n, n, work(ir), ldwrkr, work(itauq), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + n;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of R in WORK(IR)
            //C                    (Workspace: need N*N+BDSPAC)
            //C
            dbdsqr(cmn, "U", n, 0, n, 0, s, work(ie), dum, 1, work(ir),
              ldwrkr, dum, 1, work(iwork), info);
            //C
            //C                    Multiply Q in U by left singular vectors of R in
            //C                    WORK(IR), storing result in A
            //C                    (Workspace: need N*N)
            //C
            dgemm("N", "N", m, n, n, one, u, ldu, work(ir), ldwrkr,
              zero, a, lda);
            //C
            //C                    Copy left singular vectors of A from A to U
            //C
            dlacpy("F", m, n, a, lda, u, ldu);
            //C
          }
          else {
            //C
            //C                    Insufficient workspace for a fast algorithm
            //C
            itau = 1;
            iwork = itau + n;
            //C
            //C                    Compute A=Q*R, copying result to U
            //C                    (Workspace: need 2*N, prefer N+N*NB)
            //C
            dgeqrf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("L", m, n, a, lda, u, ldu);
            //C
            //C                    Generate Q in U
            //C                    (Workspace: need N+M, prefer N+M*NB)
            //C
            dorgqr(m, m, n, u, ldu, work(itau), work(iwork), lwork - iwork + 1,
              ierr);
            ie = itau;
            itauq = ie + n;
            itaup = itauq + n;
            iwork = itaup + n;
            //C
            //C                    Zero out below R in A
            //C
            dlaset("L", n - 1, n - 1, zero, zero, a(2, 1), lda);
            //C
            //C                    Bidiagonalize R in A
            //C                    (Workspace: need 4*N, prefer 3*N+2*N*NB)
            //C
            dgebrd(cmn, n, n, a, lda, s, work(ie), work(itauq), work(itaup),
              work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Multiply Q in U by left bidiagonalizing vectors
            //C                    in A
            //C                    (Workspace: need 3*N+M, prefer 3*N+M*NB)
            //C
            dormbr("Q", "R", "N", m, n, n, a, lda, work(itauq), u,
              ldu, work(iwork), lwork - iwork + 1, ierr);
            iwork = ie + n;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of A in U
            //C                    (Workspace: need BDSPAC)
            //C
            dbdsqr(cmn, "U", n, 0, m, 0, s, work(ie), dum, 1, u, ldu,
              dum, 1, work(iwork), info);
            //C
          }
          //C
        }
        else if (wntvo) {
          //C
          //C                 Path 8 (M much larger than N, JOBU='A', JOBVT='O')
          //C                 M left singular vectors to be computed in U and
          //C                 N right singular vectors to be overwritten on A
          //C
          if (lwork >= 2 * n * n + fem::max(n + m, 4 * n, bdspac)) {
            //C
            //C                    Sufficient workspace for a fast algorithm
            //C
            iu = 1;
            if (lwork >= wrkbl + 2 * lda * n) {
              //C
              //C                       WORK(IU) is LDA by N and WORK(IR) is LDA by N
              //C
              ldwrku = lda;
              ir = iu + ldwrku * n;
              ldwrkr = lda;
            }
            else if (lwork >= wrkbl + (lda + n) * n) {
              //C
              //C                       WORK(IU) is LDA by N and WORK(IR) is N by N
              //C
              ldwrku = lda;
              ir = iu + ldwrku * n;
              ldwrkr = n;
            }
            else {
              //C
              //C                       WORK(IU) is N by N and WORK(IR) is N by N
              //C
              ldwrku = n;
              ir = iu + ldwrku * n;
              ldwrkr = n;
            }
            itau = ir + ldwrkr * n;
            iwork = itau + n;
            //C
            //C                    Compute A=Q*R, copying result to U
            //C                    (Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
            //C
            dgeqrf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("L", m, n, a, lda, u, ldu);
            //C
            //C                    Generate Q in U
            //C                    (Workspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB)
            //C
            dorgqr(m, m, n, u, ldu, work(itau), work(iwork), lwork - iwork + 1,
              ierr);
            //C
            //C                    Copy R to WORK(IU), zeroing out below it
            //C
            dlacpy("U", n, n, a, lda, work(iu), ldwrku);
            dlaset("L", n - 1, n - 1, zero, zero, work(iu + 1), ldwrku);
            ie = itau;
            itauq = ie + n;
            itaup = itauq + n;
            iwork = itaup + n;
            //C
            //C                    Bidiagonalize R in WORK(IU), copying result to
            //C                    WORK(IR)
            //C                    (Workspace: need 2*N*N+4*N,
            //C                                prefer 2*N*N+3*N+2*N*NB)
            //C
            dgebrd(cmn, n, n, work(iu), ldwrku, s, work(ie), work(itauq),
              work(itaup), work(iwork), lwork - iwork + 1, ierr);
            dlacpy("U", n, n, work(iu), ldwrku, work(ir), ldwrkr);
            //C
            //C                    Generate left bidiagonalizing vectors in WORK(IU)
            //C                    (Workspace: need 2*N*N+4*N, prefer 2*N*N+3*N+N*NB)
            //C
            dorgbr("Q", n, n, n, work(iu), ldwrku, work(itauq), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Generate right bidiagonalizing vectors in WORK(IR)
            //C                    (Workspace: need 2*N*N+4*N-1,
            //C                                prefer 2*N*N+3*N+(N-1)*NB)
            //C
            dorgbr("P", n, n, n, work(ir), ldwrkr, work(itaup), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + n;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of R in WORK(IU) and computing
            //C                    right singular vectors of R in WORK(IR)
            //C                    (Workspace: need 2*N*N+BDSPAC)
            //C
            dbdsqr(cmn, "U", n, n, n, 0, s, work(ie), work(ir),
              ldwrkr, work(iu), ldwrku, dum, 1, work(iwork), info);
            //C
            //C                    Multiply Q in U by left singular vectors of R in
            //C                    WORK(IU), storing result in A
            //C                    (Workspace: need N*N)
            //C
            dgemm("N", "N", m, n, n, one, u, ldu, work(iu), ldwrku,
              zero, a, lda);
            //C
            //C                    Copy left singular vectors of A from A to U
            //C
            dlacpy("F", m, n, a, lda, u, ldu);
            //C
            //C                    Copy right singular vectors of R from WORK(IR) to A
            //C
            dlacpy("F", n, n, work(ir), ldwrkr, a, lda);
            //C
          }
          else {
            //C
            //C                    Insufficient workspace for a fast algorithm
            //C
            itau = 1;
            iwork = itau + n;
            //C
            //C                    Compute A=Q*R, copying result to U
            //C                    (Workspace: need 2*N, prefer N+N*NB)
            //C
            dgeqrf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("L", m, n, a, lda, u, ldu);
            //C
            //C                    Generate Q in U
            //C                    (Workspace: need N+M, prefer N+M*NB)
            //C
            dorgqr(m, m, n, u, ldu, work(itau), work(iwork), lwork - iwork + 1,
              ierr);
            ie = itau;
            itauq = ie + n;
            itaup = itauq + n;
            iwork = itaup + n;
            //C
            //C                    Zero out below R in A
            //C
            dlaset("L", n - 1, n - 1, zero, zero, a(2, 1), lda);
            //C
            //C                    Bidiagonalize R in A
            //C                    (Workspace: need 4*N, prefer 3*N+2*N*NB)
            //C
            dgebrd(cmn, n, n, a, lda, s, work(ie), work(itauq), work(itaup),
              work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Multiply Q in U by left bidiagonalizing vectors
            //C                    in A
            //C                    (Workspace: need 3*N+M, prefer 3*N+M*NB)
            //C
            dormbr("Q", "R", "N", m, n, n, a, lda, work(itauq), u,
              ldu, work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Generate right bidiagonalizing vectors in A
            //C                    (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
            //C
            dorgbr("P", n, n, n, a, lda, work(itaup), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + n;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of A in U and computing right
            //C                    singular vectors of A in A
            //C                    (Workspace: need BDSPAC)
            //C
            dbdsqr(cmn, "U", n, n, m, 0, s, work(ie), a, lda, u, ldu,
              dum, 1, work(iwork), info);
            //C
          }
          //C
        }
        else if (wntvas) {
          //C
          //C                 Path 9 (M much larger than N, JOBU='A', JOBVT='S'
          //C                         or 'A')
          //C                 M left singular vectors to be computed in U and
          //C                 N right singular vectors to be computed in VT
          //C
          if (lwork >= n * n + fem::max(n + m, 4 * n, bdspac)) {
            //C
            //C                    Sufficient workspace for a fast algorithm
            //C
            iu = 1;
            if (lwork >= wrkbl + lda * n) {
              //C
              //C                       WORK(IU) is LDA by N
              //C
              ldwrku = lda;
            }
            else {
              //C
              //C                       WORK(IU) is N by N
              //C
              ldwrku = n;
            }
            itau = iu + ldwrku * n;
            iwork = itau + n;
            //C
            //C                    Compute A=Q*R, copying result to U
            //C                    (Workspace: need N*N+2*N, prefer N*N+N+N*NB)
            //C
            dgeqrf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("L", m, n, a, lda, u, ldu);
            //C
            //C                    Generate Q in U
            //C                    (Workspace: need N*N+N+M, prefer N*N+N+M*NB)
            //C
            dorgqr(m, m, n, u, ldu, work(itau), work(iwork), lwork - iwork + 1,
              ierr);
            //C
            //C                    Copy R to WORK(IU), zeroing out below it
            //C
            dlacpy("U", n, n, a, lda, work(iu), ldwrku);
            dlaset("L", n - 1, n - 1, zero, zero, work(iu + 1), ldwrku);
            ie = itau;
            itauq = ie + n;
            itaup = itauq + n;
            iwork = itaup + n;
            //C
            //C                    Bidiagonalize R in WORK(IU), copying result to VT
            //C                    (Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
            //C
            dgebrd(cmn, n, n, work(iu), ldwrku, s, work(ie), work(itauq),
              work(itaup), work(iwork), lwork - iwork + 1, ierr);
            dlacpy("U", n, n, work(iu), ldwrku, vt, ldvt);
            //C
            //C                    Generate left bidiagonalizing vectors in WORK(IU)
            //C                    (Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)
            //C
            dorgbr("Q", n, n, n, work(iu), ldwrku, work(itauq), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Generate right bidiagonalizing vectors in VT
            //C                    (Workspace: need N*N+4*N-1,
            //C                                prefer N*N+3*N+(N-1)*NB)
            //C
            dorgbr("P", n, n, n, vt, ldvt, work(itaup), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + n;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of R in WORK(IU) and computing
            //C                    right singular vectors of R in VT
            //C                    (Workspace: need N*N+BDSPAC)
            //C
            dbdsqr(cmn, "U", n, n, n, 0, s, work(ie), vt, ldvt, work(iu),
              ldwrku, dum, 1, work(iwork), info);
            //C
            //C                    Multiply Q in U by left singular vectors of R in
            //C                    WORK(IU), storing result in A
            //C                    (Workspace: need N*N)
            //C
            dgemm("N", "N", m, n, n, one, u, ldu, work(iu), ldwrku,
              zero, a, lda);
            //C
            //C                    Copy left singular vectors of A from A to U
            //C
            dlacpy("F", m, n, a, lda, u, ldu);
            //C
          }
          else {
            //C
            //C                    Insufficient workspace for a fast algorithm
            //C
            itau = 1;
            iwork = itau + n;
            //C
            //C                    Compute A=Q*R, copying result to U
            //C                    (Workspace: need 2*N, prefer N+N*NB)
            //C
            dgeqrf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("L", m, n, a, lda, u, ldu);
            //C
            //C                    Generate Q in U
            //C                    (Workspace: need N+M, prefer N+M*NB)
            //C
            dorgqr(m, m, n, u, ldu, work(itau), work(iwork), lwork - iwork + 1,
              ierr);
            //C
            //C                    Copy R from A to VT, zeroing out below it
            //C
            dlacpy("U", n, n, a, lda, vt, ldvt);
            if (n > 1) {
              dlaset("L", n - 1, n - 1, zero, zero, vt(2, 1), ldvt);
            }
            ie = itau;
            itauq = ie + n;
            itaup = itauq + n;
            iwork = itaup + n;
            //C
            //C                    Bidiagonalize R in VT
            //C                    (Workspace: need 4*N, prefer 3*N+2*N*NB)
            //C
            dgebrd(cmn, n, n, vt, ldvt, s, work(ie), work(itauq),
              work(itaup), work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Multiply Q in U by left bidiagonalizing vectors
            //C                    in VT
            //C                    (Workspace: need 3*N+M, prefer 3*N+M*NB)
            //C
            dormbr("Q", "R", "N", m, n, n, vt, ldvt, work(itauq), u,
              ldu, work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Generate right bidiagonalizing vectors in VT
            //C                    (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
            //C
            dorgbr("P", n, n, n, vt, ldvt, work(itaup), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + n;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of A in U and computing right
            //C                    singular vectors of A in VT
            //C                    (Workspace: need BDSPAC)
            //C
            dbdsqr(cmn, "U", n, n, m, 0, s, work(ie), vt, ldvt, u,
              ldu, dum, 1, work(iwork), info);
            //C
          }
          //C
        }
        //C
      }
      //C
    }
    else {
      //C
      //C           M .LT. MNTHR
      //C
      //C           Path 10 (M at least N, but not much larger)
      //C           Reduce to bidiagonal form without QR decomposition
      //C
      ie = 1;
      itauq = ie + n;
      itaup = itauq + n;
      iwork = itaup + n;
      //C
      //C           Bidiagonalize A
      //C           (Workspace: need 3*N+M, prefer 3*N+(M+N)*NB)
      //C
      dgebrd(cmn, m, n, a, lda, s, work(ie), work(itauq), work(itaup),
        work(iwork), lwork - iwork + 1, ierr);
      if (wntuas) {
        //C
        //C              If left singular vectors desired in U, copy result to U
        //C              and generate left bidiagonalizing vectors in U
        //C              (Workspace: need 3*N+NCU, prefer 3*N+NCU*NB)
        //C
        dlacpy("L", m, n, a, lda, u, ldu);
        if (wntus) {
          ncu = n;
        }
        if (wntua) {
          ncu = m;
        }
        dorgbr("Q", m, ncu, n, u, ldu, work(itauq), work(iwork),
          lwork - iwork + 1, ierr);
      }
      if (wntvas) {
        //C
        //C              If right singular vectors desired in VT, copy result to
        //C              VT and generate right bidiagonalizing vectors in VT
        //C              (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
        //C
        dlacpy("U", n, n, a, lda, vt, ldvt);
        dorgbr("P", n, n, n, vt, ldvt, work(itaup), work(iwork),
          lwork - iwork + 1, ierr);
      }
      if (wntuo) {
        //C
        //C              If left singular vectors desired in A, generate left
        //C              bidiagonalizing vectors in A
        //C              (Workspace: need 4*N, prefer 3*N+N*NB)
        //C
        dorgbr("Q", m, n, n, a, lda, work(itauq), work(iwork),
          lwork - iwork + 1, ierr);
      }
      if (wntvo) {
        //C
        //C              If right singular vectors desired in A, generate right
        //C              bidiagonalizing vectors in A
        //C              (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
        //C
        dorgbr("P", n, n, n, a, lda, work(itaup), work(iwork),
          lwork - iwork + 1, ierr);
      }
      iwork = ie + n;
      if (wntuas || wntuo) {
        nru = m;
      }
      if (wntun) {
        nru = 0;
      }
      if (wntvas || wntvo) {
        ncvt = n;
      }
      if (wntvn) {
        ncvt = 0;
      }
      if ((!wntuo) && (!wntvo)) {
        //C
        //C              Perform bidiagonal QR iteration, if desired, computing
        //C              left singular vectors in U and computing right singular
        //C              vectors in VT
        //C              (Workspace: need BDSPAC)
        //C
        dbdsqr(cmn, "U", n, ncvt, nru, 0, s, work(ie), vt, ldvt, u,
          ldu, dum, 1, work(iwork), info);
      }
      else if ((!wntuo) && wntvo) {
        //C
        //C              Perform bidiagonal QR iteration, if desired, computing
        //C              left singular vectors in U and computing right singular
        //C              vectors in A
        //C              (Workspace: need BDSPAC)
        //C
        dbdsqr(cmn, "U", n, ncvt, nru, 0, s, work(ie), a, lda, u,
          ldu, dum, 1, work(iwork), info);
      }
      else {
        //C
        //C              Perform bidiagonal QR iteration, if desired, computing
        //C              left singular vectors in A and computing right singular
        //C              vectors in VT
        //C              (Workspace: need BDSPAC)
        //C
        dbdsqr(cmn, "U", n, ncvt, nru, 0, s, work(ie), vt, ldvt, a,
          lda, dum, 1, work(iwork), info);
      }
      //C
    }
    //C
  }
  else {
    //C
    //C        A has more columns than rows. If A has sufficiently more
    //C        columns than rows, first reduce using the LQ decomposition (if
    //C        sufficient workspace available)
    //C
    if (n >= mnthr) {
      //C
      if (wntvn) {
        //C
        //C              Path 1t(N much larger than M, JOBVT='N')
        //C              No right singular vectors to be computed
        //C
        itau = 1;
        iwork = itau + m;
        //C
        //C              Compute A=L*Q
        //C              (Workspace: need 2*M, prefer M+M*NB)
        //C
        dgelqf(cmn, m, n, a, lda, work(itau), work(iwork), lwork - iwork + 1,
          ierr);
        //C
        //C              Zero out above L
        //C
        dlaset("U", m - 1, m - 1, zero, zero, a(1, 2), lda);
        ie = 1;
        itauq = ie + m;
        itaup = itauq + m;
        iwork = itaup + m;
        //C
        //C              Bidiagonalize L in A
        //C              (Workspace: need 4*M, prefer 3*M+2*M*NB)
        //C
        dgebrd(cmn, m, m, a, lda, s, work(ie), work(itauq), work(itaup),
          work(iwork), lwork - iwork + 1, ierr);
        if (wntuo || wntuas) {
          //C
          //C                 If left singular vectors desired, generate Q
          //C                 (Workspace: need 4*M, prefer 3*M+M*NB)
          //C
          dorgbr("Q", m, m, m, a, lda, work(itauq), work(iwork),
            lwork - iwork + 1, ierr);
        }
        iwork = ie + m;
        nru = 0;
        if (wntuo || wntuas) {
          nru = m;
        }
        //C
        //C              Perform bidiagonal QR iteration, computing left singular
        //C              vectors of A in A if desired
        //C              (Workspace: need BDSPAC)
        //C
        dbdsqr(cmn, "U", m, 0, nru, 0, s, work(ie), dum, 1, a, lda,
          dum, 1, work(iwork), info);
        //C
        //C              If left singular vectors desired in U, copy them there
        //C
        if (wntuas) {
          dlacpy("F", m, m, a, lda, u, ldu);
        }
        //C
      }
      else if (wntvo && wntun) {
        //C
        //C              Path 2t(N much larger than M, JOBU='N', JOBVT='O')
        //C              M right singular vectors to be overwritten on A and
        //C              no left singular vectors to be computed
        //C
        if (lwork >= m * m + fem::max(4 * m, bdspac)) {
          //C
          //C                 Sufficient workspace for a fast algorithm
          //C
          ir = 1;
          if (lwork >= fem::max(wrkbl, lda * n + m) + lda * m) {
            //C
            //C                    WORK(IU) is LDA by N and WORK(IR) is LDA by M
            //C
            ldwrku = lda;
            chunk = n;
            ldwrkr = lda;
          }
          else if (lwork >= fem::max(wrkbl, lda * n + m) + m * m) {
            //C
            //C                    WORK(IU) is LDA by N and WORK(IR) is M by M
            //C
            ldwrku = lda;
            chunk = n;
            ldwrkr = m;
          }
          else {
            //C
            //C                    WORK(IU) is M by CHUNK and WORK(IR) is M by M
            //C
            ldwrku = m;
            chunk = (lwork - m * m - m) / m;
            ldwrkr = m;
          }
          itau = ir + ldwrkr * m;
          iwork = itau + m;
          //C
          //C                 Compute A=L*Q
          //C                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
          //C
          dgelqf(cmn, m, n, a, lda, work(itau), work(iwork), lwork - iwork + 1,
            ierr);
          //C
          //C                 Copy L to WORK(IR) and zero out above it
          //C
          dlacpy("L", m, m, a, lda, work(ir), ldwrkr);
          dlaset("U", m - 1, m - 1, zero, zero, work(ir + ldwrkr), ldwrkr);
          //C
          //C                 Generate Q in A
          //C                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
          //C
          dorglq(m, n, m, a, lda, work(itau), work(iwork), lwork - iwork + 1,
            ierr);
          ie = itau;
          itauq = ie + m;
          itaup = itauq + m;
          iwork = itaup + m;
          //C
          //C                 Bidiagonalize L in WORK(IR)
          //C                 (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
          //C
          dgebrd(cmn, m, m, work(ir), ldwrkr, s, work(ie), work(itauq),
            work(itaup), work(iwork), lwork - iwork + 1, ierr);
          //C
          //C                 Generate right vectors bidiagonalizing L
          //C                 (Workspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB)
          //C
          dorgbr("P", m, m, m, work(ir), ldwrkr, work(itaup), work(iwork),
            lwork - iwork + 1, ierr);
          iwork = ie + m;
          //C
          //C                 Perform bidiagonal QR iteration, computing right
          //C                 singular vectors of L in WORK(IR)
          //C                 (Workspace: need M*M+BDSPAC)
          //C
          dbdsqr(cmn, "U", m, m, 0, 0, s, work(ie), work(ir), ldwrkr,
            dum, 1, dum, 1, work(iwork), info);
          iu = ie + m;
          //C
          //C                 Multiply right singular vectors of L in WORK(IR) by Q
          //C                 in A, storing result in WORK(IU) and copying to A
          //C                 (Workspace: need M*M+2*M, prefer M*M+M*N+M)
          //C
          FEM_DOSTEP(i, 1, n, chunk) {
            blk = fem::min(n - i + 1, chunk);
            dgemm("N", "N", m, blk, m, one, work(ir), ldwrkr, a(1,
              i), lda, zero, work(iu), ldwrku);
            dlacpy("F", m, blk, work(iu), ldwrku, a(1, i), lda);
          }
          //C
        }
        else {
          //C
          //C                 Insufficient workspace for a fast algorithm
          //C
          ie = 1;
          itauq = ie + m;
          itaup = itauq + m;
          iwork = itaup + m;
          //C
          //C                 Bidiagonalize A
          //C                 (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
          //C
          dgebrd(cmn, m, n, a, lda, s, work(ie), work(itauq), work(itaup),
            work(iwork), lwork - iwork + 1, ierr);
          //C
          //C                 Generate right vectors bidiagonalizing A
          //C                 (Workspace: need 4*M, prefer 3*M+M*NB)
          //C
          dorgbr("P", m, n, m, a, lda, work(itaup), work(iwork),
            lwork - iwork + 1, ierr);
          iwork = ie + m;
          //C
          //C                 Perform bidiagonal QR iteration, computing right
          //C                 singular vectors of A in A
          //C                 (Workspace: need BDSPAC)
          //C
          dbdsqr(cmn, "L", m, n, 0, 0, s, work(ie), a, lda, dum, 1,
            dum, 1, work(iwork), info);
          //C
        }
        //C
      }
      else if (wntvo && wntuas) {
        //C
        //C              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O')
        //C              M right singular vectors to be overwritten on A and
        //C              M left singular vectors to be computed in U
        //C
        if (lwork >= m * m + fem::max(4 * m, bdspac)) {
          //C
          //C                 Sufficient workspace for a fast algorithm
          //C
          ir = 1;
          if (lwork >= fem::max(wrkbl, lda * n + m) + lda * m) {
            //C
            //C                    WORK(IU) is LDA by N and WORK(IR) is LDA by M
            //C
            ldwrku = lda;
            chunk = n;
            ldwrkr = lda;
          }
          else if (lwork >= fem::max(wrkbl, lda * n + m) + m * m) {
            //C
            //C                    WORK(IU) is LDA by N and WORK(IR) is M by M
            //C
            ldwrku = lda;
            chunk = n;
            ldwrkr = m;
          }
          else {
            //C
            //C                    WORK(IU) is M by CHUNK and WORK(IR) is M by M
            //C
            ldwrku = m;
            chunk = (lwork - m * m - m) / m;
            ldwrkr = m;
          }
          itau = ir + ldwrkr * m;
          iwork = itau + m;
          //C
          //C                 Compute A=L*Q
          //C                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
          //C
          dgelqf(cmn, m, n, a, lda, work(itau), work(iwork), lwork - iwork + 1,
            ierr);
          //C
          //C                 Copy L to U, zeroing about above it
          //C
          dlacpy("L", m, m, a, lda, u, ldu);
          dlaset("U", m - 1, m - 1, zero, zero, u(1, 2), ldu);
          //C
          //C                 Generate Q in A
          //C                 (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
          //C
          dorglq(m, n, m, a, lda, work(itau), work(iwork), lwork - iwork + 1,
            ierr);
          ie = itau;
          itauq = ie + m;
          itaup = itauq + m;
          iwork = itaup + m;
          //C
          //C                 Bidiagonalize L in U, copying result to WORK(IR)
          //C                 (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
          //C
          dgebrd(cmn, m, m, u, ldu, s, work(ie), work(itauq), work(itaup),
            work(iwork), lwork - iwork + 1, ierr);
          dlacpy("U", m, m, u, ldu, work(ir), ldwrkr);
          //C
          //C                 Generate right vectors bidiagonalizing L in WORK(IR)
          //C                 (Workspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB)
          //C
          dorgbr("P", m, m, m, work(ir), ldwrkr, work(itaup), work(iwork),
            lwork - iwork + 1, ierr);
          //C
          //C                 Generate left vectors bidiagonalizing L in U
          //C                 (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB)
          //C
          dorgbr("Q", m, m, m, u, ldu, work(itauq), work(iwork),
            lwork - iwork + 1, ierr);
          iwork = ie + m;
          //C
          //C                 Perform bidiagonal QR iteration, computing left
          //C                 singular vectors of L in U, and computing right
          //C                 singular vectors of L in WORK(IR)
          //C                 (Workspace: need M*M+BDSPAC)
          //C
          dbdsqr(cmn, "U", m, m, m, 0, s, work(ie), work(ir), ldwrkr,
            u, ldu, dum, 1, work(iwork), info);
          iu = ie + m;
          //C
          //C                 Multiply right singular vectors of L in WORK(IR) by Q
          //C                 in A, storing result in WORK(IU) and copying to A
          //C                 (Workspace: need M*M+2*M, prefer M*M+M*N+M))
          //C
          FEM_DOSTEP(i, 1, n, chunk) {
            blk = fem::min(n - i + 1, chunk);
            dgemm("N", "N", m, blk, m, one, work(ir), ldwrkr, a(1,
              i), lda, zero, work(iu), ldwrku);
            dlacpy("F", m, blk, work(iu), ldwrku, a(1, i), lda);
          }
          //C
        }
        else {
          //C
          //C                 Insufficient workspace for a fast algorithm
          //C
          itau = 1;
          iwork = itau + m;
          //C
          //C                 Compute A=L*Q
          //C                 (Workspace: need 2*M, prefer M+M*NB)
          //C
          dgelqf(cmn, m, n, a, lda, work(itau), work(iwork), lwork - iwork + 1,
            ierr);
          //C
          //C                 Copy L to U, zeroing out above it
          //C
          dlacpy("L", m, m, a, lda, u, ldu);
          dlaset("U", m - 1, m - 1, zero, zero, u(1, 2), ldu);
          //C
          //C                 Generate Q in A
          //C                 (Workspace: need 2*M, prefer M+M*NB)
          //C
          dorglq(m, n, m, a, lda, work(itau), work(iwork), lwork - iwork + 1,
            ierr);
          ie = itau;
          itauq = ie + m;
          itaup = itauq + m;
          iwork = itaup + m;
          //C
          //C                 Bidiagonalize L in U
          //C                 (Workspace: need 4*M, prefer 3*M+2*M*NB)
          //C
          dgebrd(cmn, m, m, u, ldu, s, work(ie), work(itauq), work(itaup),
            work(iwork), lwork - iwork + 1, ierr);
          //C
          //C                 Multiply right vectors bidiagonalizing L by Q in A
          //C                 (Workspace: need 3*M+N, prefer 3*M+N*NB)
          //C
          dormbr("P", "L", "T", m, n, m, u, ldu, work(itaup), a, lda,
            work(iwork), lwork - iwork + 1, ierr);
          //C
          //C                 Generate left vectors bidiagonalizing L in U
          //C                 (Workspace: need 4*M, prefer 3*M+M*NB)
          //C
          dorgbr("Q", m, m, m, u, ldu, work(itauq), work(iwork),
            lwork - iwork + 1, ierr);
          iwork = ie + m;
          //C
          //C                 Perform bidiagonal QR iteration, computing left
          //C                 singular vectors of A in U and computing right
          //C                 singular vectors of A in A
          //C                 (Workspace: need BDSPAC)
          //C
          dbdsqr(cmn, "U", m, n, m, 0, s, work(ie), a, lda, u, ldu,
            dum, 1, work(iwork), info);
          //C
        }
        //C
      }
      else if (wntvs) {
        //C
        if (wntun) {
          //C
          //C                 Path 4t(N much larger than M, JOBU='N', JOBVT='S')
          //C                 M right singular vectors to be computed in VT and
          //C                 no left singular vectors to be computed
          //C
          if (lwork >= m * m + fem::max(4 * m, bdspac)) {
            //C
            //C                    Sufficient workspace for a fast algorithm
            //C
            ir = 1;
            if (lwork >= wrkbl + lda * m) {
              //C
              //C                       WORK(IR) is LDA by M
              //C
              ldwrkr = lda;
            }
            else {
              //C
              //C                       WORK(IR) is M by M
              //C
              ldwrkr = m;
            }
            itau = ir + ldwrkr * m;
            iwork = itau + m;
            //C
            //C                    Compute A=L*Q
            //C                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
            //C
            dgelqf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Copy L to WORK(IR), zeroing out above it
            //C
            dlacpy("L", m, m, a, lda, work(ir), ldwrkr);
            dlaset("U", m - 1, m - 1, zero, zero, work(ir + ldwrkr), ldwrkr);
            //C
            //C                    Generate Q in A
            //C                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
            //C
            dorglq(m, n, m, a, lda, work(itau), work(iwork), lwork - iwork + 1,
              ierr);
            ie = itau;
            itauq = ie + m;
            itaup = itauq + m;
            iwork = itaup + m;
            //C
            //C                    Bidiagonalize L in WORK(IR)
            //C                    (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
            //C
            dgebrd(cmn, m, m, work(ir), ldwrkr, s, work(ie), work(itauq),
              work(itaup), work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Generate right vectors bidiagonalizing L in
            //C                    WORK(IR)
            //C                    (Workspace: need M*M+4*M, prefer M*M+3*M+(M-1)*NB)
            //C
            dorgbr("P", m, m, m, work(ir), ldwrkr, work(itaup), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + m;
            //C
            //C                    Perform bidiagonal QR iteration, computing right
            //C                    singular vectors of L in WORK(IR)
            //C                    (Workspace: need M*M+BDSPAC)
            //C
            dbdsqr(cmn, "U", m, m, 0, 0, s, work(ie), work(ir),
              ldwrkr, dum, 1, dum, 1, work(iwork), info);
            //C
            //C                    Multiply right singular vectors of L in WORK(IR) by
            //C                    Q in A, storing result in VT
            //C                    (Workspace: need M*M)
            //C
            dgemm("N", "N", m, n, m, one, work(ir), ldwrkr, a, lda,
              zero, vt, ldvt);
            //C
          }
          else {
            //C
            //C                    Insufficient workspace for a fast algorithm
            //C
            itau = 1;
            iwork = itau + m;
            //C
            //C                    Compute A=L*Q
            //C                    (Workspace: need 2*M, prefer M+M*NB)
            //C
            dgelqf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Copy result to VT
            //C
            dlacpy("U", m, n, a, lda, vt, ldvt);
            //C
            //C                    Generate Q in VT
            //C                    (Workspace: need 2*M, prefer M+M*NB)
            //C
            dorglq(m, n, m, vt, ldvt, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            ie = itau;
            itauq = ie + m;
            itaup = itauq + m;
            iwork = itaup + m;
            //C
            //C                    Zero out above L in A
            //C
            dlaset("U", m - 1, m - 1, zero, zero, a(1, 2), lda);
            //C
            //C                    Bidiagonalize L in A
            //C                    (Workspace: need 4*M, prefer 3*M+2*M*NB)
            //C
            dgebrd(cmn, m, m, a, lda, s, work(ie), work(itauq), work(itaup),
              work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Multiply right vectors bidiagonalizing L by Q in VT
            //C                    (Workspace: need 3*M+N, prefer 3*M+N*NB)
            //C
            dormbr("P", "L", "T", m, n, m, a, lda, work(itaup), vt,
              ldvt, work(iwork), lwork - iwork + 1, ierr);
            iwork = ie + m;
            //C
            //C                    Perform bidiagonal QR iteration, computing right
            //C                    singular vectors of A in VT
            //C                    (Workspace: need BDSPAC)
            //C
            dbdsqr(cmn, "U", m, n, 0, 0, s, work(ie), vt, ldvt, dum,
              1, dum, 1, work(iwork), info);
            //C
          }
          //C
        }
        else if (wntuo) {
          //C
          //C                 Path 5t(N much larger than M, JOBU='O', JOBVT='S')
          //C                 M right singular vectors to be computed in VT and
          //C                 M left singular vectors to be overwritten on A
          //C
          if (lwork >= 2 * m * m + fem::max(4 * m, bdspac)) {
            //C
            //C                    Sufficient workspace for a fast algorithm
            //C
            iu = 1;
            if (lwork >= wrkbl + 2 * lda * m) {
              //C
              //C                       WORK(IU) is LDA by M and WORK(IR) is LDA by M
              //C
              ldwrku = lda;
              ir = iu + ldwrku * m;
              ldwrkr = lda;
            }
            else if (lwork >= wrkbl + (lda + m) * m) {
              //C
              //C                       WORK(IU) is LDA by M and WORK(IR) is M by M
              //C
              ldwrku = lda;
              ir = iu + ldwrku * m;
              ldwrkr = m;
            }
            else {
              //C
              //C                       WORK(IU) is M by M and WORK(IR) is M by M
              //C
              ldwrku = m;
              ir = iu + ldwrku * m;
              ldwrkr = m;
            }
            itau = ir + ldwrkr * m;
            iwork = itau + m;
            //C
            //C                    Compute A=L*Q
            //C                    (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
            //C
            dgelqf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Copy L to WORK(IU), zeroing out below it
            //C
            dlacpy("L", m, m, a, lda, work(iu), ldwrku);
            dlaset("U", m - 1, m - 1, zero, zero, work(iu + ldwrku), ldwrku);
            //C
            //C                    Generate Q in A
            //C                    (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
            //C
            dorglq(m, n, m, a, lda, work(itau), work(iwork), lwork - iwork + 1,
              ierr);
            ie = itau;
            itauq = ie + m;
            itaup = itauq + m;
            iwork = itaup + m;
            //C
            //C                    Bidiagonalize L in WORK(IU), copying result to
            //C                    WORK(IR)
            //C                    (Workspace: need 2*M*M+4*M,
            //C                                prefer 2*M*M+3*M+2*M*NB)
            //C
            dgebrd(cmn, m, m, work(iu), ldwrku, s, work(ie), work(itauq),
              work(itaup), work(iwork), lwork - iwork + 1, ierr);
            dlacpy("L", m, m, work(iu), ldwrku, work(ir), ldwrkr);
            //C
            //C                    Generate right bidiagonalizing vectors in WORK(IU)
            //C                    (Workspace: need 2*M*M+4*M-1,
            //C                                prefer 2*M*M+3*M+(M-1)*NB)
            //C
            dorgbr("P", m, m, m, work(iu), ldwrku, work(itaup), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Generate left bidiagonalizing vectors in WORK(IR)
            //C                    (Workspace: need 2*M*M+4*M, prefer 2*M*M+3*M+M*NB)
            //C
            dorgbr("Q", m, m, m, work(ir), ldwrkr, work(itauq), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + m;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of L in WORK(IR) and computing
            //C                    right singular vectors of L in WORK(IU)
            //C                    (Workspace: need 2*M*M+BDSPAC)
            //C
            dbdsqr(cmn, "U", m, m, m, 0, s, work(ie), work(iu),
              ldwrku, work(ir), ldwrkr, dum, 1, work(iwork), info);
            //C
            //C                    Multiply right singular vectors of L in WORK(IU) by
            //C                    Q in A, storing result in VT
            //C                    (Workspace: need M*M)
            //C
            dgemm("N", "N", m, n, m, one, work(iu), ldwrku, a, lda,
              zero, vt, ldvt);
            //C
            //C                    Copy left singular vectors of L to A
            //C                    (Workspace: need M*M)
            //C
            dlacpy("F", m, m, work(ir), ldwrkr, a, lda);
            //C
          }
          else {
            //C
            //C                    Insufficient workspace for a fast algorithm
            //C
            itau = 1;
            iwork = itau + m;
            //C
            //C                    Compute A=L*Q, copying result to VT
            //C                    (Workspace: need 2*M, prefer M+M*NB)
            //C
            dgelqf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("U", m, n, a, lda, vt, ldvt);
            //C
            //C                    Generate Q in VT
            //C                    (Workspace: need 2*M, prefer M+M*NB)
            //C
            dorglq(m, n, m, vt, ldvt, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            ie = itau;
            itauq = ie + m;
            itaup = itauq + m;
            iwork = itaup + m;
            //C
            //C                    Zero out above L in A
            //C
            dlaset("U", m - 1, m - 1, zero, zero, a(1, 2), lda);
            //C
            //C                    Bidiagonalize L in A
            //C                    (Workspace: need 4*M, prefer 3*M+2*M*NB)
            //C
            dgebrd(cmn, m, m, a, lda, s, work(ie), work(itauq), work(itaup),
              work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Multiply right vectors bidiagonalizing L by Q in VT
            //C                    (Workspace: need 3*M+N, prefer 3*M+N*NB)
            //C
            dormbr("P", "L", "T", m, n, m, a, lda, work(itaup), vt,
              ldvt, work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Generate left bidiagonalizing vectors of L in A
            //C                    (Workspace: need 4*M, prefer 3*M+M*NB)
            //C
            dorgbr("Q", m, m, m, a, lda, work(itauq), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + m;
            //C
            //C                    Perform bidiagonal QR iteration, compute left
            //C                    singular vectors of A in A and compute right
            //C                    singular vectors of A in VT
            //C                    (Workspace: need BDSPAC)
            //C
            dbdsqr(cmn, "U", m, n, m, 0, s, work(ie), vt, ldvt, a,
              lda, dum, 1, work(iwork), info);
            //C
          }
          //C
        }
        else if (wntuas) {
          //C
          //C                 Path 6t(N much larger than M, JOBU='S' or 'A',
          //C                         JOBVT='S')
          //C                 M right singular vectors to be computed in VT and
          //C                 M left singular vectors to be computed in U
          //C
          if (lwork >= m * m + fem::max(4 * m, bdspac)) {
            //C
            //C                    Sufficient workspace for a fast algorithm
            //C
            iu = 1;
            if (lwork >= wrkbl + lda * m) {
              //C
              //C                       WORK(IU) is LDA by N
              //C
              ldwrku = lda;
            }
            else {
              //C
              //C                       WORK(IU) is LDA by M
              //C
              ldwrku = m;
            }
            itau = iu + ldwrku * m;
            iwork = itau + m;
            //C
            //C                    Compute A=L*Q
            //C                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
            //C
            dgelqf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Copy L to WORK(IU), zeroing out above it
            //C
            dlacpy("L", m, m, a, lda, work(iu), ldwrku);
            dlaset("U", m - 1, m - 1, zero, zero, work(iu + ldwrku), ldwrku);
            //C
            //C                    Generate Q in A
            //C                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
            //C
            dorglq(m, n, m, a, lda, work(itau), work(iwork), lwork - iwork + 1,
              ierr);
            ie = itau;
            itauq = ie + m;
            itaup = itauq + m;
            iwork = itaup + m;
            //C
            //C                    Bidiagonalize L in WORK(IU), copying result to U
            //C                    (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
            //C
            dgebrd(cmn, m, m, work(iu), ldwrku, s, work(ie), work(itauq),
              work(itaup), work(iwork), lwork - iwork + 1, ierr);
            dlacpy("L", m, m, work(iu), ldwrku, u, ldu);
            //C
            //C                    Generate right bidiagonalizing vectors in WORK(IU)
            //C                    (Workspace: need M*M+4*M-1,
            //C                                prefer M*M+3*M+(M-1)*NB)
            //C
            dorgbr("P", m, m, m, work(iu), ldwrku, work(itaup), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Generate left bidiagonalizing vectors in U
            //C                    (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB)
            //C
            dorgbr("Q", m, m, m, u, ldu, work(itauq), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + m;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of L in U and computing right
            //C                    singular vectors of L in WORK(IU)
            //C                    (Workspace: need M*M+BDSPAC)
            //C
            dbdsqr(cmn, "U", m, m, m, 0, s, work(ie), work(iu),
              ldwrku, u, ldu, dum, 1, work(iwork), info);
            //C
            //C                    Multiply right singular vectors of L in WORK(IU) by
            //C                    Q in A, storing result in VT
            //C                    (Workspace: need M*M)
            //C
            dgemm("N", "N", m, n, m, one, work(iu), ldwrku, a, lda,
              zero, vt, ldvt);
            //C
          }
          else {
            //C
            //C                    Insufficient workspace for a fast algorithm
            //C
            itau = 1;
            iwork = itau + m;
            //C
            //C                    Compute A=L*Q, copying result to VT
            //C                    (Workspace: need 2*M, prefer M+M*NB)
            //C
            dgelqf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("U", m, n, a, lda, vt, ldvt);
            //C
            //C                    Generate Q in VT
            //C                    (Workspace: need 2*M, prefer M+M*NB)
            //C
            dorglq(m, n, m, vt, ldvt, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Copy L to U, zeroing out above it
            //C
            dlacpy("L", m, m, a, lda, u, ldu);
            dlaset("U", m - 1, m - 1, zero, zero, u(1, 2), ldu);
            ie = itau;
            itauq = ie + m;
            itaup = itauq + m;
            iwork = itaup + m;
            //C
            //C                    Bidiagonalize L in U
            //C                    (Workspace: need 4*M, prefer 3*M+2*M*NB)
            //C
            dgebrd(cmn, m, m, u, ldu, s, work(ie), work(itauq), work(itaup),
              work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Multiply right bidiagonalizing vectors in U by Q
            //C                    in VT
            //C                    (Workspace: need 3*M+N, prefer 3*M+N*NB)
            //C
            dormbr("P", "L", "T", m, n, m, u, ldu, work(itaup), vt,
              ldvt, work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Generate left bidiagonalizing vectors in U
            //C                    (Workspace: need 4*M, prefer 3*M+M*NB)
            //C
            dorgbr("Q", m, m, m, u, ldu, work(itauq), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + m;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of A in U and computing right
            //C                    singular vectors of A in VT
            //C                    (Workspace: need BDSPAC)
            //C
            dbdsqr(cmn, "U", m, n, m, 0, s, work(ie), vt, ldvt, u,
              ldu, dum, 1, work(iwork), info);
            //C
          }
          //C
        }
        //C
      }
      else if (wntva) {
        //C
        if (wntun) {
          //C
          //C                 Path 7t(N much larger than M, JOBU='N', JOBVT='A')
          //C                 N right singular vectors to be computed in VT and
          //C                 no left singular vectors to be computed
          //C
          if (lwork >= m * m + fem::max(n + m, 4 * m, bdspac)) {
            //C
            //C                    Sufficient workspace for a fast algorithm
            //C
            ir = 1;
            if (lwork >= wrkbl + lda * m) {
              //C
              //C                       WORK(IR) is LDA by M
              //C
              ldwrkr = lda;
            }
            else {
              //C
              //C                       WORK(IR) is M by M
              //C
              ldwrkr = m;
            }
            itau = ir + ldwrkr * m;
            iwork = itau + m;
            //C
            //C                    Compute A=L*Q, copying result to VT
            //C                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
            //C
            dgelqf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("U", m, n, a, lda, vt, ldvt);
            //C
            //C                    Copy L to WORK(IR), zeroing out above it
            //C
            dlacpy("L", m, m, a, lda, work(ir), ldwrkr);
            dlaset("U", m - 1, m - 1, zero, zero, work(ir + ldwrkr), ldwrkr);
            //C
            //C                    Generate Q in VT
            //C                    (Workspace: need M*M+M+N, prefer M*M+M+N*NB)
            //C
            dorglq(n, n, m, vt, ldvt, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            ie = itau;
            itauq = ie + m;
            itaup = itauq + m;
            iwork = itaup + m;
            //C
            //C                    Bidiagonalize L in WORK(IR)
            //C                    (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
            //C
            dgebrd(cmn, m, m, work(ir), ldwrkr, s, work(ie), work(itauq),
              work(itaup), work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Generate right bidiagonalizing vectors in WORK(IR)
            //C                    (Workspace: need M*M+4*M-1,
            //C                                prefer M*M+3*M+(M-1)*NB)
            //C
            dorgbr("P", m, m, m, work(ir), ldwrkr, work(itaup), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + m;
            //C
            //C                    Perform bidiagonal QR iteration, computing right
            //C                    singular vectors of L in WORK(IR)
            //C                    (Workspace: need M*M+BDSPAC)
            //C
            dbdsqr(cmn, "U", m, m, 0, 0, s, work(ie), work(ir),
              ldwrkr, dum, 1, dum, 1, work(iwork), info);
            //C
            //C                    Multiply right singular vectors of L in WORK(IR) by
            //C                    Q in VT, storing result in A
            //C                    (Workspace: need M*M)
            //C
            dgemm("N", "N", m, n, m, one, work(ir), ldwrkr, vt, ldvt,
              zero, a, lda);
            //C
            //C                    Copy right singular vectors of A from A to VT
            //C
            dlacpy("F", m, n, a, lda, vt, ldvt);
            //C
          }
          else {
            //C
            //C                    Insufficient workspace for a fast algorithm
            //C
            itau = 1;
            iwork = itau + m;
            //C
            //C                    Compute A=L*Q, copying result to VT
            //C                    (Workspace: need 2*M, prefer M+M*NB)
            //C
            dgelqf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("U", m, n, a, lda, vt, ldvt);
            //C
            //C                    Generate Q in VT
            //C                    (Workspace: need M+N, prefer M+N*NB)
            //C
            dorglq(n, n, m, vt, ldvt, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            ie = itau;
            itauq = ie + m;
            itaup = itauq + m;
            iwork = itaup + m;
            //C
            //C                    Zero out above L in A
            //C
            dlaset("U", m - 1, m - 1, zero, zero, a(1, 2), lda);
            //C
            //C                    Bidiagonalize L in A
            //C                    (Workspace: need 4*M, prefer 3*M+2*M*NB)
            //C
            dgebrd(cmn, m, m, a, lda, s, work(ie), work(itauq), work(itaup),
              work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Multiply right bidiagonalizing vectors in A by Q
            //C                    in VT
            //C                    (Workspace: need 3*M+N, prefer 3*M+N*NB)
            //C
            dormbr("P", "L", "T", m, n, m, a, lda, work(itaup), vt,
              ldvt, work(iwork), lwork - iwork + 1, ierr);
            iwork = ie + m;
            //C
            //C                    Perform bidiagonal QR iteration, computing right
            //C                    singular vectors of A in VT
            //C                    (Workspace: need BDSPAC)
            //C
            dbdsqr(cmn, "U", m, n, 0, 0, s, work(ie), vt, ldvt, dum,
              1, dum, 1, work(iwork), info);
            //C
          }
          //C
        }
        else if (wntuo) {
          //C
          //C                 Path 8t(N much larger than M, JOBU='O', JOBVT='A')
          //C                 N right singular vectors to be computed in VT and
          //C                 M left singular vectors to be overwritten on A
          //C
          if (lwork >= 2 * m * m + fem::max(n + m, 4 * m, bdspac)) {
            //C
            //C                    Sufficient workspace for a fast algorithm
            //C
            iu = 1;
            if (lwork >= wrkbl + 2 * lda * m) {
              //C
              //C                       WORK(IU) is LDA by M and WORK(IR) is LDA by M
              //C
              ldwrku = lda;
              ir = iu + ldwrku * m;
              ldwrkr = lda;
            }
            else if (lwork >= wrkbl + (lda + m) * m) {
              //C
              //C                       WORK(IU) is LDA by M and WORK(IR) is M by M
              //C
              ldwrku = lda;
              ir = iu + ldwrku * m;
              ldwrkr = m;
            }
            else {
              //C
              //C                       WORK(IU) is M by M and WORK(IR) is M by M
              //C
              ldwrku = m;
              ir = iu + ldwrku * m;
              ldwrkr = m;
            }
            itau = ir + ldwrkr * m;
            iwork = itau + m;
            //C
            //C                    Compute A=L*Q, copying result to VT
            //C                    (Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
            //C
            dgelqf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("U", m, n, a, lda, vt, ldvt);
            //C
            //C                    Generate Q in VT
            //C                    (Workspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB)
            //C
            dorglq(n, n, m, vt, ldvt, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Copy L to WORK(IU), zeroing out above it
            //C
            dlacpy("L", m, m, a, lda, work(iu), ldwrku);
            dlaset("U", m - 1, m - 1, zero, zero, work(iu + ldwrku), ldwrku);
            ie = itau;
            itauq = ie + m;
            itaup = itauq + m;
            iwork = itaup + m;
            //C
            //C                    Bidiagonalize L in WORK(IU), copying result to
            //C                    WORK(IR)
            //C                    (Workspace: need 2*M*M+4*M,
            //C                                prefer 2*M*M+3*M+2*M*NB)
            //C
            dgebrd(cmn, m, m, work(iu), ldwrku, s, work(ie), work(itauq),
              work(itaup), work(iwork), lwork - iwork + 1, ierr);
            dlacpy("L", m, m, work(iu), ldwrku, work(ir), ldwrkr);
            //C
            //C                    Generate right bidiagonalizing vectors in WORK(IU)
            //C                    (Workspace: need 2*M*M+4*M-1,
            //C                                prefer 2*M*M+3*M+(M-1)*NB)
            //C
            dorgbr("P", m, m, m, work(iu), ldwrku, work(itaup), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Generate left bidiagonalizing vectors in WORK(IR)
            //C                    (Workspace: need 2*M*M+4*M, prefer 2*M*M+3*M+M*NB)
            //C
            dorgbr("Q", m, m, m, work(ir), ldwrkr, work(itauq), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + m;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of L in WORK(IR) and computing
            //C                    right singular vectors of L in WORK(IU)
            //C                    (Workspace: need 2*M*M+BDSPAC)
            //C
            dbdsqr(cmn, "U", m, m, m, 0, s, work(ie), work(iu),
              ldwrku, work(ir), ldwrkr, dum, 1, work(iwork), info);
            //C
            //C                    Multiply right singular vectors of L in WORK(IU) by
            //C                    Q in VT, storing result in A
            //C                    (Workspace: need M*M)
            //C
            dgemm("N", "N", m, n, m, one, work(iu), ldwrku, vt, ldvt,
              zero, a, lda);
            //C
            //C                    Copy right singular vectors of A from A to VT
            //C
            dlacpy("F", m, n, a, lda, vt, ldvt);
            //C
            //C                    Copy left singular vectors of A from WORK(IR) to A
            //C
            dlacpy("F", m, m, work(ir), ldwrkr, a, lda);
            //C
          }
          else {
            //C
            //C                    Insufficient workspace for a fast algorithm
            //C
            itau = 1;
            iwork = itau + m;
            //C
            //C                    Compute A=L*Q, copying result to VT
            //C                    (Workspace: need 2*M, prefer M+M*NB)
            //C
            dgelqf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("U", m, n, a, lda, vt, ldvt);
            //C
            //C                    Generate Q in VT
            //C                    (Workspace: need M+N, prefer M+N*NB)
            //C
            dorglq(n, n, m, vt, ldvt, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            ie = itau;
            itauq = ie + m;
            itaup = itauq + m;
            iwork = itaup + m;
            //C
            //C                    Zero out above L in A
            //C
            dlaset("U", m - 1, m - 1, zero, zero, a(1, 2), lda);
            //C
            //C                    Bidiagonalize L in A
            //C                    (Workspace: need 4*M, prefer 3*M+2*M*NB)
            //C
            dgebrd(cmn, m, m, a, lda, s, work(ie), work(itauq), work(itaup),
              work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Multiply right bidiagonalizing vectors in A by Q
            //C                    in VT
            //C                    (Workspace: need 3*M+N, prefer 3*M+N*NB)
            //C
            dormbr("P", "L", "T", m, n, m, a, lda, work(itaup), vt,
              ldvt, work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Generate left bidiagonalizing vectors in A
            //C                    (Workspace: need 4*M, prefer 3*M+M*NB)
            //C
            dorgbr("Q", m, m, m, a, lda, work(itauq), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + m;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of A in A and computing right
            //C                    singular vectors of A in VT
            //C                    (Workspace: need BDSPAC)
            //C
            dbdsqr(cmn, "U", m, n, m, 0, s, work(ie), vt, ldvt, a,
              lda, dum, 1, work(iwork), info);
            //C
          }
          //C
        }
        else if (wntuas) {
          //C
          //C                 Path 9t(N much larger than M, JOBU='S' or 'A',
          //C                         JOBVT='A')
          //C                 N right singular vectors to be computed in VT and
          //C                 M left singular vectors to be computed in U
          //C
          if (lwork >= m * m + fem::max(n + m, 4 * m, bdspac)) {
            //C
            //C                    Sufficient workspace for a fast algorithm
            //C
            iu = 1;
            if (lwork >= wrkbl + lda * m) {
              //C
              //C                       WORK(IU) is LDA by M
              //C
              ldwrku = lda;
            }
            else {
              //C
              //C                       WORK(IU) is M by M
              //C
              ldwrku = m;
            }
            itau = iu + ldwrku * m;
            iwork = itau + m;
            //C
            //C                    Compute A=L*Q, copying result to VT
            //C                    (Workspace: need M*M+2*M, prefer M*M+M+M*NB)
            //C
            dgelqf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("U", m, n, a, lda, vt, ldvt);
            //C
            //C                    Generate Q in VT
            //C                    (Workspace: need M*M+M+N, prefer M*M+M+N*NB)
            //C
            dorglq(n, n, m, vt, ldvt, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Copy L to WORK(IU), zeroing out above it
            //C
            dlacpy("L", m, m, a, lda, work(iu), ldwrku);
            dlaset("U", m - 1, m - 1, zero, zero, work(iu + ldwrku), ldwrku);
            ie = itau;
            itauq = ie + m;
            itaup = itauq + m;
            iwork = itaup + m;
            //C
            //C                    Bidiagonalize L in WORK(IU), copying result to U
            //C                    (Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
            //C
            dgebrd(cmn, m, m, work(iu), ldwrku, s, work(ie), work(itauq),
              work(itaup), work(iwork), lwork - iwork + 1, ierr);
            dlacpy("L", m, m, work(iu), ldwrku, u, ldu);
            //C
            //C                    Generate right bidiagonalizing vectors in WORK(IU)
            //C                    (Workspace: need M*M+4*M, prefer M*M+3*M+(M-1)*NB)
            //C
            dorgbr("P", m, m, m, work(iu), ldwrku, work(itaup), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Generate left bidiagonalizing vectors in U
            //C                    (Workspace: need M*M+4*M, prefer M*M+3*M+M*NB)
            //C
            dorgbr("Q", m, m, m, u, ldu, work(itauq), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + m;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of L in U and computing right
            //C                    singular vectors of L in WORK(IU)
            //C                    (Workspace: need M*M+BDSPAC)
            //C
            dbdsqr(cmn, "U", m, m, m, 0, s, work(ie), work(iu),
              ldwrku, u, ldu, dum, 1, work(iwork), info);
            //C
            //C                    Multiply right singular vectors of L in WORK(IU) by
            //C                    Q in VT, storing result in A
            //C                    (Workspace: need M*M)
            //C
            dgemm("N", "N", m, n, m, one, work(iu), ldwrku, vt, ldvt,
              zero, a, lda);
            //C
            //C                    Copy right singular vectors of A from A to VT
            //C
            dlacpy("F", m, n, a, lda, vt, ldvt);
            //C
          }
          else {
            //C
            //C                    Insufficient workspace for a fast algorithm
            //C
            itau = 1;
            iwork = itau + m;
            //C
            //C                    Compute A=L*Q, copying result to VT
            //C                    (Workspace: need 2*M, prefer M+M*NB)
            //C
            dgelqf(cmn, m, n, a, lda, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            dlacpy("U", m, n, a, lda, vt, ldvt);
            //C
            //C                    Generate Q in VT
            //C                    (Workspace: need M+N, prefer M+N*NB)
            //C
            dorglq(n, n, m, vt, ldvt, work(itau), work(iwork),
              lwork - iwork + 1, ierr);
            //C
            //C                    Copy L to U, zeroing out above it
            //C
            dlacpy("L", m, m, a, lda, u, ldu);
            dlaset("U", m - 1, m - 1, zero, zero, u(1, 2), ldu);
            ie = itau;
            itauq = ie + m;
            itaup = itauq + m;
            iwork = itaup + m;
            //C
            //C                    Bidiagonalize L in U
            //C                    (Workspace: need 4*M, prefer 3*M+2*M*NB)
            //C
            dgebrd(cmn, m, m, u, ldu, s, work(ie), work(itauq), work(itaup),
              work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Multiply right bidiagonalizing vectors in U by Q
            //C                    in VT
            //C                    (Workspace: need 3*M+N, prefer 3*M+N*NB)
            //C
            dormbr("P", "L", "T", m, n, m, u, ldu, work(itaup), vt,
              ldvt, work(iwork), lwork - iwork + 1, ierr);
            //C
            //C                    Generate left bidiagonalizing vectors in U
            //C                    (Workspace: need 4*M, prefer 3*M+M*NB)
            //C
            dorgbr("Q", m, m, m, u, ldu, work(itauq), work(iwork),
              lwork - iwork + 1, ierr);
            iwork = ie + m;
            //C
            //C                    Perform bidiagonal QR iteration, computing left
            //C                    singular vectors of A in U and computing right
            //C                    singular vectors of A in VT
            //C                    (Workspace: need BDSPAC)
            //C
            dbdsqr(cmn, "U", m, n, m, 0, s, work(ie), vt, ldvt, u,
              ldu, dum, 1, work(iwork), info);
            //C
          }
          //C
        }
        //C
      }
      //C
    }
    else {
      //C
      //C           N .LT. MNTHR
      //C
      //C           Path 10t(N greater than M, but not much larger)
      //C           Reduce to bidiagonal form without LQ decomposition
      //C
      ie = 1;
      itauq = ie + m;
      itaup = itauq + m;
      iwork = itaup + m;
      //C
      //C           Bidiagonalize A
      //C           (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
      //C
      dgebrd(cmn, m, n, a, lda, s, work(ie), work(itauq), work(itaup),
        work(iwork), lwork - iwork + 1, ierr);
      if (wntuas) {
        //C
        //C              If left singular vectors desired in U, copy result to U
        //C              and generate left bidiagonalizing vectors in U
        //C              (Workspace: need 4*M-1, prefer 3*M+(M-1)*NB)
        //C
        dlacpy("L", m, m, a, lda, u, ldu);
        dorgbr("Q", m, m, n, u, ldu, work(itauq), work(iwork),
          lwork - iwork + 1, ierr);
      }
      if (wntvas) {
        //C
        //C              If right singular vectors desired in VT, copy result to
        //C              VT and generate right bidiagonalizing vectors in VT
        //C              (Workspace: need 3*M+NRVT, prefer 3*M+NRVT*NB)
        //C
        dlacpy("U", m, n, a, lda, vt, ldvt);
        if (wntva) {
          nrvt = n;
        }
        if (wntvs) {
          nrvt = m;
        }
        dorgbr("P", nrvt, n, m, vt, ldvt, work(itaup), work(iwork),
          lwork - iwork + 1, ierr);
      }
      if (wntuo) {
        //C
        //C              If left singular vectors desired in A, generate left
        //C              bidiagonalizing vectors in A
        //C              (Workspace: need 4*M-1, prefer 3*M+(M-1)*NB)
        //C
        dorgbr("Q", m, m, n, a, lda, work(itauq), work(iwork),
          lwork - iwork + 1, ierr);
      }
      if (wntvo) {
        //C
        //C              If right singular vectors desired in A, generate right
        //C              bidiagonalizing vectors in A
        //C              (Workspace: need 4*M, prefer 3*M+M*NB)
        //C
        dorgbr("P", m, n, m, a, lda, work(itaup), work(iwork),
          lwork - iwork + 1, ierr);
      }
      iwork = ie + m;
      if (wntuas || wntuo) {
        nru = m;
      }
      if (wntun) {
        nru = 0;
      }
      if (wntvas || wntvo) {
        ncvt = n;
      }
      if (wntvn) {
        ncvt = 0;
      }
      if ((!wntuo) && (!wntvo)) {
        //C
        //C              Perform bidiagonal QR iteration, if desired, computing
        //C              left singular vectors in U and computing right singular
        //C              vectors in VT
        //C              (Workspace: need BDSPAC)
        //C
        dbdsqr(cmn, "L", m, ncvt, nru, 0, s, work(ie), vt, ldvt, u,
          ldu, dum, 1, work(iwork), info);
      }
      else if ((!wntuo) && wntvo) {
        //C
        //C              Perform bidiagonal QR iteration, if desired, computing
        //C              left singular vectors in U and computing right singular
        //C              vectors in A
        //C              (Workspace: need BDSPAC)
        //C
        dbdsqr(cmn, "L", m, ncvt, nru, 0, s, work(ie), a, lda, u,
          ldu, dum, 1, work(iwork), info);
      }
      else {
        //C
        //C              Perform bidiagonal QR iteration, if desired, computing
        //C              left singular vectors in A and computing right singular
        //C              vectors in VT
        //C              (Workspace: need BDSPAC)
        //C
        dbdsqr(cmn, "L", m, ncvt, nru, 0, s, work(ie), vt, ldvt, a,
          lda, dum, 1, work(iwork), info);
      }
      //C
    }
    //C
  }
  //C
  //C     If DBDSQR failed to converge, copy unconverged superdiagonals
  //C     to WORK( 2:MINMN )
  //C
  if (info != 0) {
    if (ie > 2) {
      {
        int fem_do_last = minmn - 1;
        FEM_DO(i, 1, fem_do_last) {
          work(i + 1) = work(i + ie - 1);
        }
      }
    }
    if (ie < 2) {
      FEM_DOSTEP(i, minmn - 1, 1, -1) {
        work(i + 1) = work(i + ie - 1);
      }
    }
  }
  //C
  //C     Undo scaling if necessary
  //C
  if (iscl == 1) {
    if (anrm > bignum) {
      dlascl(cmn, "G", 0, 0, bignum, anrm, minmn, 1, s, minmn, ierr);
    }
    if (info != 0 && anrm > bignum) {
      dlascl(cmn, "G", 0, 0, bignum, anrm, minmn - 1, 1, work(2), minmn, ierr);
    }
    if (anrm < smlnum) {
      dlascl(cmn, "G", 0, 0, smlnum, anrm, minmn, 1, s, minmn, ierr);
    }
    if (info != 0 && anrm < smlnum) {
      dlascl(cmn, "G", 0, 0, smlnum, anrm, minmn - 1, 1, work(2), minmn, ierr);
    }
  }
  //C
  //C     Return optimal workspace in WORK(1)
  //C
  work(1) = maxwrk;
  //C
  //C     End of DGESVD
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlae2.f
inline
void
dlae2(
  double const& a,
  double const& b,
  double const& c,
  double& rt1,
  double& rt2)
{
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
  //C     [  A   B  ]
  //C     [  B   C  ].
  //C  On return, RT1 is the eigenvalue of larger absolute value, and RT2
  //C  is the eigenvalue of smaller absolute value.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  A       (input) DOUBLE PRECISION
  //C          The (1,1) element of the 2-by-2 matrix.
  //C
  //C  B       (input) DOUBLE PRECISION
  //C          The (1,2) and (2,1) elements of the 2-by-2 matrix.
  //C
  //C  C       (input) DOUBLE PRECISION
  //C          The (2,2) element of the 2-by-2 matrix.
  //C
  //C  RT1     (output) DOUBLE PRECISION
  //C          The eigenvalue of larger absolute value.
  //C
  //C  RT2     (output) DOUBLE PRECISION
  //C          The eigenvalue of smaller absolute value.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  RT1 is accurate to a few ulps barring over/underflow.
  //C
  //C  RT2 may be inaccurate if there is massive cancellation in the
  //C  determinant A*C-B*B; higher precision or correctly rounded or
  //C  correctly truncated arithmetic would be needed to compute RT2
  //C  accurately in all cases.
  //C
  //C  Overflow is possible only if RT1 is within a factor of 5 of overflow.
  //C  Underflow is harmless if the input data is 0 or exceeds
  //C     underflow_threshold / macheps.
  //C
  //C =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Compute the eigenvalues
  //C
  double sm = a + c;
  double df = a - c;
  double adf = fem::abs(df);
  double tb = b + b;
  double ab = fem::abs(tb);
  double acmx = fem::double0;
  double acmn = fem::double0;
  if (fem::abs(a) > fem::abs(c)) {
    acmx = a;
    acmn = c;
  }
  else {
    acmx = c;
    acmn = a;
  }
  const double one = 1.0e0;
  double rt = fem::double0;
  const double two = 2.0e0;
  if (adf > ab) {
    rt = adf * fem::sqrt(one + fem::pow2((ab / adf)));
  }
  else if (adf < ab) {
    rt = ab * fem::sqrt(one + fem::pow2((adf / ab)));
  }
  else {
    //C
    //C        Includes case AB=ADF=0
    //C
    rt = ab * fem::sqrt(two);
  }
  const double zero = 0.0e0;
  const double half = 0.5e0;
  if (sm < zero) {
    rt1 = half * (sm - rt);
    //C
    //C        Order of execution important.
    //C        To get fully accurate smaller eigenvalue,
    //C        next line needs to be executed in higher precision.
    //C
    rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
  }
  else if (sm > zero) {
    rt1 = half * (sm + rt);
    //C
    //C        Order of execution important.
    //C        To get fully accurate smaller eigenvalue,
    //C        next line needs to be executed in higher precision.
    //C
    rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
  }
  else {
    //C
    //C        Includes case RT1 = RT2 = 0
    //C
    rt1 = half * rt;
    rt2 = -half * rt;
  }
  //C
  //C     End of DLAE2
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlaev2.f
inline
void
dlaev2(
  double const& a,
  double const& b,
  double const& c,
  double& rt1,
  double& rt2,
  double& cs1,
  double& sn1)
{
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
  //C     [  A   B  ]
  //C     [  B   C  ].
  //C  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
  //C  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
  //C  eigenvector for RT1, giving the decomposition
  //C
  //C     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
  //C     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
  //C
  //C  Arguments
  //C  =========
  //C
  //C  A       (input) DOUBLE PRECISION
  //C          The (1,1) element of the 2-by-2 matrix.
  //C
  //C  B       (input) DOUBLE PRECISION
  //C          The (1,2) element and the conjugate of the (2,1) element of
  //C          the 2-by-2 matrix.
  //C
  //C  C       (input) DOUBLE PRECISION
  //C          The (2,2) element of the 2-by-2 matrix.
  //C
  //C  RT1     (output) DOUBLE PRECISION
  //C          The eigenvalue of larger absolute value.
  //C
  //C  RT2     (output) DOUBLE PRECISION
  //C          The eigenvalue of smaller absolute value.
  //C
  //C  CS1     (output) DOUBLE PRECISION
  //C  SN1     (output) DOUBLE PRECISION
  //C          The vector (CS1, SN1) is a unit right eigenvector for RT1.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  RT1 is accurate to a few ulps barring over/underflow.
  //C
  //C  RT2 may be inaccurate if there is massive cancellation in the
  //C  determinant A*C-B*B; higher precision or correctly rounded or
  //C  correctly truncated arithmetic would be needed to compute RT2
  //C  accurately in all cases.
  //C
  //C  CS1 and SN1 are accurate to a few ulps barring over/underflow.
  //C
  //C  Overflow is possible only if RT1 is within a factor of 5 of overflow.
  //C  Underflow is harmless if the input data is 0 or exceeds
  //C     underflow_threshold / macheps.
  //C
  //C =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Compute the eigenvalues
  //C
  double sm = a + c;
  double df = a - c;
  double adf = fem::abs(df);
  double tb = b + b;
  double ab = fem::abs(tb);
  double acmx = fem::double0;
  double acmn = fem::double0;
  if (fem::abs(a) > fem::abs(c)) {
    acmx = a;
    acmn = c;
  }
  else {
    acmx = c;
    acmn = a;
  }
  const double one = 1.0e0;
  double rt = fem::double0;
  const double two = 2.0e0;
  if (adf > ab) {
    rt = adf * fem::sqrt(one + fem::pow2((ab / adf)));
  }
  else if (adf < ab) {
    rt = ab * fem::sqrt(one + fem::pow2((adf / ab)));
  }
  else {
    //C
    //C        Includes case AB=ADF=0
    //C
    rt = ab * fem::sqrt(two);
  }
  const double zero = 0.0e0;
  const double half = 0.5e0;
  int sgn1 = fem::int0;
  if (sm < zero) {
    rt1 = half * (sm - rt);
    sgn1 = -1;
    //C
    //C        Order of execution important.
    //C        To get fully accurate smaller eigenvalue,
    //C        next line needs to be executed in higher precision.
    //C
    rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
  }
  else if (sm > zero) {
    rt1 = half * (sm + rt);
    sgn1 = 1;
    //C
    //C        Order of execution important.
    //C        To get fully accurate smaller eigenvalue,
    //C        next line needs to be executed in higher precision.
    //C
    rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
  }
  else {
    //C
    //C        Includes case RT1 = RT2 = 0
    //C
    rt1 = half * rt;
    rt2 = -half * rt;
    sgn1 = 1;
  }
  //C
  //C     Compute the eigenvector
  //C
  double cs = fem::double0;
  int sgn2 = fem::int0;
  if (df >= zero) {
    cs = df + rt;
    sgn2 = 1;
  }
  else {
    cs = df - rt;
    sgn2 = -1;
  }
  double acs = fem::abs(cs);
  double ct = fem::double0;
  double tn = fem::double0;
  if (acs > ab) {
    ct = -tb / cs;
    sn1 = one / fem::sqrt(one + ct * ct);
    cs1 = ct * sn1;
  }
  else {
    if (ab == zero) {
      cs1 = one;
      sn1 = zero;
    }
    else {
      tn = -cs / tb;
      cs1 = one / fem::sqrt(one + tn * tn);
      sn1 = tn * cs1;
    }
  }
  if (sgn1 == sgn2) {
    tn = cs1;
    cs1 = -sn1;
    sn1 = tn;
  }
  //C
  //C     End of DLAEV2
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlansy.f
inline
double
dlansy(
  str_cref norm,
  str_cref uplo,
  int const& n,
  arr_cref<double, 2> a,
  int const& lda,
  arr_ref<double> work)
{
  double return_value = fem::double0;
  a(dimension(lda, star));
  work(dimension(star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLANSY  returns the value of the one norm,  or the Frobenius norm, or
  //C  the  infinity norm,  or the  element of  largest absolute value  of a
  //C  real symmetric matrix A.
  //C
  //C  Description
  //C  ===========
  //C
  //C  DLANSY returns the value
  //C
  //C     DLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'
  //C              (
  //C              ( norm1(A),         NORM = '1', 'O' or 'o'
  //C              (
  //C              ( normI(A),         NORM = 'I' or 'i'
  //C              (
  //C              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
  //C
  //C  where  norm1  denotes the  one norm of a matrix (maximum column sum),
  //C  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
  //C  normF  denotes the  Frobenius norm of a matrix (square root of sum of
  //C  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  NORM    (input) CHARACTER*1
  //C          Specifies the value to be returned in DLANSY as described
  //C          above.
  //C
  //C  UPLO    (input) CHARACTER*1
  //C          Specifies whether the upper or lower triangular part of the
  //C          symmetric matrix A is to be referenced.
  //C          = 'U':  Upper triangular part of A is referenced
  //C          = 'L':  Lower triangular part of A is referenced
  //C
  //C  N       (input) INTEGER
  //C          The order of the matrix A.  N >= 0.  When N = 0, DLANSY is
  //C          set to zero.
  //C
  //C  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
  //C          The symmetric matrix A.  If UPLO = 'U', the leading n by n
  //C          upper triangular part of A contains the upper triangular part
  //C          of the matrix A, and the strictly lower triangular part of A
  //C          is not referenced.  If UPLO = 'L', the leading n by n lower
  //C          triangular part of A contains the lower triangular part of
  //C          the matrix A, and the strictly upper triangular part of A is
  //C          not referenced.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(N,1).
  //C
  //C  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
  //C          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
  //C          WORK is not referenced.
  //C
  //C =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  const double zero = 0.0e+0;
  double value = fem::double0;
  int j = fem::int0;
  int i = fem::int0;
  double sum = fem::double0;
  double absa = fem::double0;
  double scale = fem::double0;
  const double one = 1.0e+0;
  if (n == 0) {
    value = zero;
  }
  else if (lsame(norm, "M")) {
    //C
    //C        Find max(abs(A(i,j))).
    //C
    value = zero;
    if (lsame(uplo, "U")) {
      FEM_DO(j, 1, n) {
        FEM_DO(i, 1, j) {
          value = fem::max(value, fem::abs(a(i, j)));
        }
      }
    }
    else {
      FEM_DO(j, 1, n) {
        FEM_DO(i, j, n) {
          value = fem::max(value, fem::abs(a(i, j)));
        }
      }
    }
  }
  else if ((lsame(norm, "I")) || (lsame(norm, "O")) || (norm == "1")) {
    //C
    //C        Find normI(A) ( = norm1(A), since A is symmetric).
    //C
    value = zero;
    if (lsame(uplo, "U")) {
      FEM_DO(j, 1, n) {
        sum = zero;
        {
          int fem_do_last = j - 1;
          FEM_DO(i, 1, fem_do_last) {
            absa = fem::abs(a(i, j));
            sum += absa;
            work(i) += absa;
          }
        }
        work(j) = sum + fem::abs(a(j, j));
      }
      FEM_DO(i, 1, n) {
        value = fem::max(value, work(i));
      }
    }
    else {
      FEM_DO(i, 1, n) {
        work(i) = zero;
      }
      FEM_DO(j, 1, n) {
        sum = work(j) + fem::abs(a(j, j));
        FEM_DO(i, j + 1, n) {
          absa = fem::abs(a(i, j));
          sum += absa;
          work(i) += absa;
        }
        value = fem::max(value, sum);
      }
    }
  }
  else if ((lsame(norm, "F")) || (lsame(norm, "E"))) {
    //C
    //C        Find normF(A).
    //C
    scale = zero;
    sum = one;
    if (lsame(uplo, "U")) {
      FEM_DO(j, 2, n) {
        dlassq(j - 1, a(1, j), 1, scale, sum);
      }
    }
    else {
      {
        int fem_do_last = n - 1;
        FEM_DO(j, 1, fem_do_last) {
          dlassq(n - j, a(j + 1, j), 1, scale, sum);
        }
      }
    }
    sum = 2 * sum;
    dlassq(n, a, lda + 1, scale, sum);
    value = scale * fem::sqrt(sum);
  }
  //C
  return_value = value;
  return return_value;
  //C
  //C     End of DLANSY
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dlatrd.f
inline
void
dlatrd(
  common& cmn,
  str_cref uplo,
  int const& n,
  int const& nb,
  arr_ref<double, 2> a,
  int const& lda,
  arr_ref<double> e,
  arr_ref<double> tau,
  arr_ref<double, 2> w,
  int const& ldw)
{
  a(dimension(lda, star));
  e(dimension(star));
  tau(dimension(star));
  w(dimension(ldw, star));
  //C
  //C  -- LAPACK auxiliary routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DLATRD reduces NB rows and columns of a real symmetric matrix A to
  //C  symmetric tridiagonal form by an orthogonal similarity
  //C  transformation Q' * A * Q, and returns the matrices V and W which are
  //C  needed to apply the transformation to the unreduced part of A.
  //C
  //C  If UPLO = 'U', DLATRD reduces the last NB rows and columns of a
  //C  matrix, of which the upper triangle is supplied;
  //C  if UPLO = 'L', DLATRD reduces the first NB rows and columns of a
  //C  matrix, of which the lower triangle is supplied.
  //C
  //C  This is an auxiliary routine called by DSYTRD.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  UPLO    (input) CHARACTER*1
  //C          Specifies whether the upper or lower triangular part of the
  //C          symmetric matrix A is stored:
  //C          = 'U': Upper triangular
  //C          = 'L': Lower triangular
  //C
  //C  N       (input) INTEGER
  //C          The order of the matrix A.
  //C
  //C  NB      (input) INTEGER
  //C          The number of rows and columns to be reduced.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  //C          n-by-n upper triangular part of A contains the upper
  //C          triangular part of the matrix A, and the strictly lower
  //C          triangular part of A is not referenced.  If UPLO = 'L', the
  //C          leading n-by-n lower triangular part of A contains the lower
  //C          triangular part of the matrix A, and the strictly upper
  //C          triangular part of A is not referenced.
  //C          On exit:
  //C          if UPLO = 'U', the last NB columns have been reduced to
  //C            tridiagonal form, with the diagonal elements overwriting
  //C            the diagonal elements of A; the elements above the diagonal
  //C            with the array TAU, represent the orthogonal matrix Q as a
  //C            product of elementary reflectors;
  //C          if UPLO = 'L', the first NB columns have been reduced to
  //C            tridiagonal form, with the diagonal elements overwriting
  //C            the diagonal elements of A; the elements below the diagonal
  //C            with the array TAU, represent the  orthogonal matrix Q as a
  //C            product of elementary reflectors.
  //C          See Further Details.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= (1,N).
  //C
  //C  E       (output) DOUBLE PRECISION array, dimension (N-1)
  //C          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
  //C          elements of the last NB columns of the reduced matrix;
  //C          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
  //C          the first NB columns of the reduced matrix.
  //C
  //C  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
  //C          The scalar factors of the elementary reflectors, stored in
  //C          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
  //C          See Further Details.
  //C
  //C  W       (output) DOUBLE PRECISION array, dimension (LDW,NB)
  //C          The n-by-nb matrix W required to update the unreduced part
  //C          of A.
  //C
  //C  LDW     (input) INTEGER
  //C          The leading dimension of the array W. LDW >= max(1,N).
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  If UPLO = 'U', the matrix Q is represented as a product of elementary
  //C  reflectors
  //C
  //C     Q = H(n) H(n-1) . . . H(n-nb+1).
  //C
  //C  Each H(i) has the form
  //C
  //C     H(i) = I - tau * v * v'
  //C
  //C  where tau is a real scalar, and v is a real vector with
  //C  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
  //C  and tau in TAU(i-1).
  //C
  //C  If UPLO = 'L', the matrix Q is represented as a product of elementary
  //C  reflectors
  //C
  //C     Q = H(1) H(2) . . . H(nb).
  //C
  //C  Each H(i) has the form
  //C
  //C     H(i) = I - tau * v * v'
  //C
  //C  where tau is a real scalar, and v is a real vector with
  //C  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
  //C  and tau in TAU(i).
  //C
  //C  The elements of the vectors v together form the n-by-nb matrix V
  //C  which is needed, with W, to apply the transformation to the unreduced
  //C  part of the matrix, using a symmetric rank-2k update of the form:
  //C  A := A - V*W' - W*V'.
  //C
  //C  The contents of A on exit are illustrated by the following examples
  //C  with n = 5 and nb = 2:
  //C
  //C  if UPLO = 'U':                       if UPLO = 'L':
  //C
  //C    (  a   a   a   v4  v5 )              (  d                  )
  //C    (      a   a   v4  v5 )              (  1   d              )
  //C    (          a   1   v5 )              (  v1  1   a          )
  //C    (              d   1  )              (  v1  v2  a   a      )
  //C    (                  d  )              (  v1  v2  a   a   a  )
  //C
  //C  where d denotes a diagonal element of the reduced matrix, a denotes
  //C  an element of the original matrix that is unchanged, and vi denotes
  //C  an element of the vector defining H(i).
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Quick return if possible
  //C
  if (n <= 0) {
    return;
  }
  //C
  int i = fem::int0;
  int iw = fem::int0;
  const double one = 1.0e+0;
  const double zero = 0.0e+0;
  const double half = 0.5e+0;
  double alpha = fem::double0;
  if (lsame(uplo, "U")) {
    //C
    //C        Reduce last NB columns of upper triangle
    //C
    FEM_DOSTEP(i, n, n - nb + 1, -1) {
      iw = i - n + nb;
      if (i < n) {
        //C
        //C              Update A(1:i,i)
        //C
        dgemv("No transpose", i, n - i, -one, a(1, i + 1), lda, w(i,
          iw + 1), ldw, one, a(1, i), 1);
        dgemv("No transpose", i, n - i, -one, w(1, iw + 1), ldw, a(i,
          i + 1), lda, one, a(1, i), 1);
      }
      if (i > 1) {
        //C
        //C              Generate elementary reflector H(i) to annihilate
        //C              A(1:i-2,i)
        //C
        dlarfg(cmn, i - 1, a(i - 1, i), a(1, i), 1, tau(i - 1));
        e(i - 1) = a(i - 1, i);
        a(i - 1, i) = one;
        //C
        //C              Compute W(1:i-1,i)
        //C
        dsymv("Upper", i - 1, one, a, lda, a(1, i), 1, zero, w(1, iw), 1);
        if (i < n) {
          dgemv("Transpose", i - 1, n - i, one, w(1, iw + 1), ldw, a(1,
            i), 1, zero, w(i + 1, iw), 1);
          dgemv("No transpose", i - 1, n - i, -one, a(1, i + 1), lda,
            w(i + 1, iw), 1, one, w(1, iw), 1);
          dgemv("Transpose", i - 1, n - i, one, a(1, i + 1), lda, a(1,
            i), 1, zero, w(i + 1, iw), 1);
          dgemv("No transpose", i - 1, n - i, -one, w(1, iw + 1),
            ldw, w(i + 1, iw), 1, one, w(1, iw), 1);
        }
        dscal(i - 1, tau(i - 1), w(1, iw), 1);
        alpha = -half * tau(i - 1) * ddot(i - 1, w(1, iw), 1, a(1, i), 1);
        daxpy(i - 1, alpha, a(1, i), 1, w(1, iw), 1);
      }
      //C
    }
  }
  else {
    //C
    //C        Reduce first NB columns of lower triangle
    //C
    FEM_DO(i, 1, nb) {
      //C
      //C           Update A(i:n,i)
      //C
      dgemv("No transpose", n - i + 1, i - 1, -one, a(i, 1), lda, w(i,
        1), ldw, one, a(i, i), 1);
      dgemv("No transpose", n - i + 1, i - 1, -one, w(i, 1), ldw, a(i,
        1), lda, one, a(i, i), 1);
      if (i < n) {
        //C
        //C              Generate elementary reflector H(i) to annihilate
        //C              A(i+2:n,i)
        //C
        dlarfg(cmn, n - i, a(i + 1, i), a(fem::min(i + 2, n), i), 1, tau(i));
        e(i) = a(i + 1, i);
        a(i + 1, i) = one;
        //C
        //C              Compute W(i+1:n,i)
        //C
        dsymv("Lower", n - i, one, a(i + 1, i + 1), lda, a(i + 1, i),
          1, zero, w(i + 1, i), 1);
        dgemv("Transpose", n - i, i - 1, one, w(i + 1, 1), ldw, a(i + 1,
          i), 1, zero, w(1, i), 1);
        dgemv("No transpose", n - i, i - 1, -one, a(i + 1, 1), lda, w(1,
          i), 1, one, w(i + 1, i), 1);
        dgemv("Transpose", n - i, i - 1, one, a(i + 1, 1), lda, a(i + 1,
          i), 1, zero, w(1, i), 1);
        dgemv("No transpose", n - i, i - 1, -one, w(i + 1, 1), ldw, w(1,
          i), 1, one, w(i + 1, i), 1);
        dscal(n - i, tau(i), w(i + 1, i), 1);
        alpha = -half * tau(i) * ddot(n - i, w(i + 1, i), 1, a(i + 1, i), 1);
        daxpy(n - i, alpha, a(i + 1, i), 1, w(i + 1, i), 1);
      }
      //C
    }
  }
  //C
  //C     End of DLATRD
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dorg2l.f
inline
void
dorg2l(
  int const& m,
  int const& n,
  int const& k,
  arr_ref<double, 2> a,
  int const& lda,
  arr_cref<double> tau,
  arr_ref<double> work,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DORG2L generates an m by n real matrix Q with orthonormal columns,
  //C  which is defined as the last n columns of a product of k elementary
  //C  reflectors of order m
  //C
  //C        Q  =  H(k) . . . H(2) H(1)
  //C
  //C  as returned by DGEQLF.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix Q. M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix Q. M >= N >= 0.
  //C
  //C  K       (input) INTEGER
  //C          The number of elementary reflectors whose product defines the
  //C          matrix Q. N >= K >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the (n-k+i)-th column must contain the vector which
  //C          defines the elementary reflector H(i), for i = 1,2,...,k, as
  //C          returned by DGEQLF in the last k columns of its array
  //C          argument A.
  //C          On exit, the m by n matrix Q.
  //C
  //C  LDA     (input) INTEGER
  //C          The first dimension of the array A. LDA >= max(1,M).
  //C
  //C  TAU     (input) DOUBLE PRECISION array, dimension (K)
  //C          TAU(i) must contain the scalar factor of the elementary
  //C          reflector H(i), as returned by DGEQLF.
  //C
  //C  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
  //C
  //C  INFO    (output) INTEGER
  //C          = 0: successful exit
  //C          < 0: if INFO = -i, the i-th argument has an illegal value
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  if (m < 0) {
    info = -1;
  }
  else if (n < 0 || n > m) {
    info = -2;
  }
  else if (k < 0 || k > n) {
    info = -3;
  }
  else if (lda < fem::max(1, m)) {
    info = -5;
  }
  if (info != 0) {
    xerbla("DORG2L", -info);
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (n <= 0) {
    return;
  }
  //C
  //C     Initialise columns 1:n-k to columns of the unit matrix
  //C
  int j = fem::int0;
  int l = fem::int0;
  const double zero = 0.0e+0;
  const double one = 1.0e+0;
  {
    int fem_do_last = n - k;
    FEM_DO(j, 1, fem_do_last) {
      FEM_DO(l, 1, m) {
        a(l, j) = zero;
      }
      a(m - n + j, j) = one;
    }
  }
  //C
  int i = fem::int0;
  int ii = fem::int0;
  FEM_DO(i, 1, k) {
    ii = n - k + i;
    //C
    //C        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
    //C
    a(m - n + ii, ii) = one;
    dlarf("Left", m - n + ii, ii - 1, a(1, ii), 1, tau(i), a, lda, work);
    dscal(m - n + ii - 1, -tau(i), a(1, ii), 1);
    a(m - n + ii, ii) = one - tau(i);
    //C
    //C        Set A(m-k+i+1:m,n-k+i) to zero
    //C
    FEM_DO(l, m - n + ii + 1, m) {
      a(l, ii) = zero;
    }
  }
  //C
  //C     End of DORG2L
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dorgql.f
inline
void
dorgql(
  int const& m,
  int const& n,
  int const& k,
  arr_ref<double, 2> a,
  int const& lda,
  arr_cref<double> tau,
  arr_ref<double> work,
  int const& lwork,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DORGQL generates an M-by-N real matrix Q with orthonormal columns,
  //C  which is defined as the last N columns of a product of K elementary
  //C  reflectors of order M
  //C
  //C        Q  =  H(k) . . . H(2) H(1)
  //C
  //C  as returned by DGEQLF.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  M       (input) INTEGER
  //C          The number of rows of the matrix Q. M >= 0.
  //C
  //C  N       (input) INTEGER
  //C          The number of columns of the matrix Q. M >= N >= 0.
  //C
  //C  K       (input) INTEGER
  //C          The number of elementary reflectors whose product defines the
  //C          matrix Q. N >= K >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the (n-k+i)-th column must contain the vector which
  //C          defines the elementary reflector H(i), for i = 1,2,...,k, as
  //C          returned by DGEQLF in the last k columns of its array
  //C          argument A.
  //C          On exit, the M-by-N matrix Q.
  //C
  //C  LDA     (input) INTEGER
  //C          The first dimension of the array A. LDA >= max(1,M).
  //C
  //C  TAU     (input) DOUBLE PRECISION array, dimension (K)
  //C          TAU(i) must contain the scalar factor of the elementary
  //C          reflector H(i), as returned by DGEQLF.
  //C
  //C  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  //C
  //C  LWORK   (input) INTEGER
  //C          The dimension of the array WORK. LWORK >= max(1,N).
  //C          For optimum performance LWORK >= N*NB, where NB is the
  //C          optimal blocksize.
  //C
  //C          If LWORK = -1, then a workspace query is assumed; the routine
  //C          only calculates the optimal size of the WORK array, returns
  //C          this value as the first entry of the WORK array, and no error
  //C          message related to LWORK is issued by XERBLA.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument has an illegal value
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  bool lquery = (lwork ==  - 1);
  if (m < 0) {
    info = -1;
  }
  else if (n < 0 || n > m) {
    info = -2;
  }
  else if (k < 0 || k > n) {
    info = -3;
  }
  else if (lda < fem::max(1, m)) {
    info = -5;
  }
  //C
  int lwkopt = fem::int0;
  int nb = fem::int0;
  if (info == 0) {
    if (n == 0) {
      lwkopt = 1;
    }
    else {
      nb = ilaenv(1, "DORGQL", " ", m, n, k, -1);
      lwkopt = n * nb;
    }
    work(1) = lwkopt;
    //C
    if (lwork < fem::max(1, n) && !lquery) {
      info = -8;
    }
  }
  //C
  if (info != 0) {
    xerbla("DORGQL", -info);
    return;
  }
  else if (lquery) {
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (n <= 0) {
    return;
  }
  //C
  int nbmin = 2;
  int nx = 0;
  int iws = n;
  int ldwork = fem::int0;
  if (nb > 1 && nb < k) {
    //C
    //C        Determine when to cross over from blocked to unblocked code.
    //C
    nx = fem::max(0, ilaenv(3, "DORGQL", " ", m, n, k, -1));
    if (nx < k) {
      //C
      //C           Determine if workspace is large enough for blocked code.
      //C
      ldwork = n;
      iws = ldwork * nb;
      if (lwork < iws) {
        //C
        //C              Not enough workspace to use optimal NB:  reduce NB and
        //C              determine the minimum value of NB.
        //C
        nb = lwork / ldwork;
        nbmin = fem::max(2, ilaenv(2, "DORGQL", " ", m, n, k, -1));
      }
    }
  }
  //C
  int kk = fem::int0;
  int j = fem::int0;
  int i = fem::int0;
  const double zero = 0.0e+0;
  if (nb >= nbmin && nb < k && nx < k) {
    //C
    //C        Use blocked code after the first block.
    //C        The last kk columns are handled by the block method.
    //C
    kk = fem::min(k, ((k - nx + nb - 1) / nb) * nb);
    //C
    //C        Set A(m-kk+1:m,1:n-kk) to zero.
    //C
    {
      int fem_do_last = n - kk;
      FEM_DO(j, 1, fem_do_last) {
        FEM_DO(i, m - kk + 1, m) {
          a(i, j) = zero;
        }
      }
    }
  }
  else {
    kk = 0;
  }
  //C
  //C     Use unblocked code for the first or only block.
  //C
  int iinfo = fem::int0;
  dorg2l(m - kk, n - kk, k - kk, a, lda, tau, work, iinfo);
  //C
  int ib = fem::int0;
  int l = fem::int0;
  if (kk > 0) {
    //C
    //C        Use blocked code
    //C
    FEM_DOSTEP(i, k - kk + 1, k, nb) {
      ib = fem::min(nb, k - i + 1);
      if (n - k + i > 1) {
        //C
        //C              Form the triangular factor of the block reflector
        //C              H = H(i+ib-1) . . . H(i+1) H(i)
        //C
        dlarft("Backward", "Columnwise", m - k + i + ib - 1, ib, a(1,
          n - k + i), lda, tau(i), work, ldwork);
        //C
        //C              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
        //C
        dlarfb("Left", "No transpose", "Backward", "Columnwise", m -
          k + i + ib - 1, n - k + i - 1, ib, a(1, n - k + i), lda,
          work, ldwork, a, lda, work(ib + 1), ldwork);
      }
      //C
      //C           Apply H to rows 1:m-k+i+ib-1 of current block
      //C
      dorg2l(m - k + i + ib - 1, ib, ib, a(1, n - k + i), lda, tau(i),
        work, iinfo);
      //C
      //C           Set rows m-k+i+ib:m of current block to zero
      //C
      {
        int fem_do_last = n - k + i + ib - 1;
        FEM_DO(j, n - k + i, fem_do_last) {
          FEM_DO(l, m - k + i + ib, m) {
            a(l, j) = zero;
          }
        }
      }
    }
  }
  //C
  work(1) = iws;
  //C
  //C     End of DORGQL
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dorgtr.f
inline
void
dorgtr(
  str_cref uplo,
  int const& n,
  arr_ref<double, 2> a,
  int const& lda,
  arr_cref<double> tau,
  arr_ref<double> work,
  int const& lwork,
  int& info)
{
  a(dimension(lda, star));
  tau(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DORGTR generates a real orthogonal matrix Q which is defined as the
  //C  product of n-1 elementary reflectors of order N, as returned by
  //C  DSYTRD:
  //C
  //C  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
  //C
  //C  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
  //C
  //C  Arguments
  //C  =========
  //C
  //C  UPLO    (input) CHARACTER*1
  //C          = 'U': Upper triangle of A contains elementary reflectors
  //C                 from DSYTRD;
  //C          = 'L': Lower triangle of A contains elementary reflectors
  //C                 from DSYTRD.
  //C
  //C  N       (input) INTEGER
  //C          The order of the matrix Q. N >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the vectors which define the elementary reflectors,
  //C          as returned by DSYTRD.
  //C          On exit, the N-by-N orthogonal matrix Q.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A. LDA >= max(1,N).
  //C
  //C  TAU     (input) DOUBLE PRECISION array, dimension (N-1)
  //C          TAU(i) must contain the scalar factor of the elementary
  //C          reflector H(i), as returned by DSYTRD.
  //C
  //C  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  //C
  //C  LWORK   (input) INTEGER
  //C          The dimension of the array WORK. LWORK >= max(1,N-1).
  //C          For optimum performance LWORK >= (N-1)*NB, where NB is
  //C          the optimal blocksize.
  //C
  //C          If LWORK = -1, then a workspace query is assumed; the routine
  //C          only calculates the optimal size of the WORK array, returns
  //C          this value as the first entry of the WORK array, and no error
  //C          message related to LWORK is issued by XERBLA.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input arguments
  //C
  info = 0;
  bool lquery = (lwork ==  - 1);
  bool upper = lsame(uplo, "U");
  if (!upper && !lsame(uplo, "L")) {
    info = -1;
  }
  else if (n < 0) {
    info = -2;
  }
  else if (lda < fem::max(1, n)) {
    info = -4;
  }
  else if (lwork < fem::max(1, n - 1) && !lquery) {
    info = -7;
  }
  //C
  int nb = fem::int0;
  int lwkopt = fem::int0;
  if (info == 0) {
    if (upper) {
      nb = ilaenv(1, "DORGQL", " ", n - 1, n - 1, n - 1, -1);
    }
    else {
      nb = ilaenv(1, "DORGQR", " ", n - 1, n - 1, n - 1, -1);
    }
    lwkopt = fem::max(1, n - 1) * nb;
    work(1) = lwkopt;
  }
  //C
  if (info != 0) {
    xerbla("DORGTR", -info);
    return;
  }
  else if (lquery) {
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (n == 0) {
    work(1) = 1;
    return;
  }
  //C
  int j = fem::int0;
  int i = fem::int0;
  const double zero = 0.0e+0;
  const double one = 1.0e+0;
  int iinfo = fem::int0;
  if (upper) {
    //C
    //C        Q was determined by a call to DSYTRD with UPLO = 'U'
    //C
    //C        Shift the vectors which define the elementary reflectors one
    //C        column to the left, and set the last row and column of Q to
    //C        those of the unit matrix
    //C
    {
      int fem_do_last = n - 1;
      FEM_DO(j, 1, fem_do_last) {
        {
          int fem_do_last = j - 1;
          FEM_DO(i, 1, fem_do_last) {
            a(i, j) = a(i, j + 1);
          }
        }
        a(n, j) = zero;
      }
    }
    {
      int fem_do_last = n - 1;
      FEM_DO(i, 1, fem_do_last) {
        a(i, n) = zero;
      }
    }
    a(n, n) = one;
    //C
    //C        Generate Q(1:n-1,1:n-1)
    //C
    dorgql(n - 1, n - 1, n - 1, a, lda, tau, work, lwork, iinfo);
    //C
  }
  else {
    //C
    //C        Q was determined by a call to DSYTRD with UPLO = 'L'.
    //C
    //C        Shift the vectors which define the elementary reflectors one
    //C        column to the right, and set the first row and column of Q to
    //C        those of the unit matrix
    //C
    FEM_DOSTEP(j, n, 2, -1) {
      a(1, j) = zero;
      FEM_DO(i, j + 1, n) {
        a(i, j) = a(i, j - 1);
      }
    }
    a(1, 1) = one;
    FEM_DO(i, 2, n) {
      a(i, 1) = zero;
    }
    if (n > 1) {
      //C
      //C           Generate Q(2:n,2:n)
      //C
      dorgqr(n - 1, n - 1, n - 1, a(2, 2), lda, tau, work, lwork, iinfo);
    }
  }
  work(1) = lwkopt;
  //C
  //C     End of DORGTR
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dsteqr.f
inline
void
dsteqr(
  common& cmn,
  str_cref compz,
  int const& n,
  arr_ref<double> d,
  arr_ref<double> e,
  arr_ref<double, 2> z,
  int const& ldz,
  arr_ref<double> work,
  int& info)
{
  d(dimension(star));
  e(dimension(star));
  z(dimension(ldz, star));
  work(dimension(star));
  int icompz = fem::int0;
  const double one = 1.0e0;
  double eps = fem::double0;
  double eps2 = fem::double0;
  double safmin = fem::double0;
  double safmax = fem::double0;
  const double three = 3.0e0;
  double ssfmax = fem::double0;
  double ssfmin = fem::double0;
  const double zero = 0.0e0;
  const int maxit = 30;
  int nmaxit = fem::int0;
  int jtot = fem::int0;
  int l1 = fem::int0;
  int nm1 = fem::int0;
  int m = fem::int0;
  double tst = fem::double0;
  int l = fem::int0;
  int lsv = fem::int0;
  int lend = fem::int0;
  int lendsv = fem::int0;
  double anorm = fem::double0;
  int iscale = fem::int0;
  int lendm1 = fem::int0;
  double p = fem::double0;
  double rt1 = fem::double0;
  double rt2 = fem::double0;
  double c = fem::double0;
  double s = fem::double0;
  const double two = 2.0e0;
  double g = fem::double0;
  double r = fem::double0;
  int mm1 = fem::int0;
  int i = fem::int0;
  double f = fem::double0;
  double b = fem::double0;
  int mm = fem::int0;
  int lendp1 = fem::int0;
  int lm1 = fem::int0;
  int ii = fem::int0;
  int k = fem::int0;
  int j = fem::int0;
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DSTEQR computes all eigenvalues and, optionally, eigenvectors of a
  //C  symmetric tridiagonal matrix using the implicit QL or QR method.
  //C  The eigenvectors of a full or band symmetric matrix can also be found
  //C  if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
  //C  tridiagonal form.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  COMPZ   (input) CHARACTER*1
  //C          = 'N':  Compute eigenvalues only.
  //C          = 'V':  Compute eigenvalues and eigenvectors of the original
  //C                  symmetric matrix.  On entry, Z must contain the
  //C                  orthogonal matrix used to reduce the original matrix
  //C                  to tridiagonal form.
  //C          = 'I':  Compute eigenvalues and eigenvectors of the
  //C                  tridiagonal matrix.  Z is initialized to the identity
  //C                  matrix.
  //C
  //C  N       (input) INTEGER
  //C          The order of the matrix.  N >= 0.
  //C
  //C  D       (input/output) DOUBLE PRECISION array, dimension (N)
  //C          On entry, the diagonal elements of the tridiagonal matrix.
  //C          On exit, if INFO = 0, the eigenvalues in ascending order.
  //C
  //C  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
  //C          On entry, the (n-1) subdiagonal elements of the tridiagonal
  //C          matrix.
  //C          On exit, E has been destroyed.
  //C
  //C  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
  //C          On entry, if  COMPZ = 'V', then Z contains the orthogonal
  //C          matrix used in the reduction to tridiagonal form.
  //C          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
  //C          orthonormal eigenvectors of the original symmetric matrix,
  //C          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
  //C          of the symmetric tridiagonal matrix.
  //C          If COMPZ = 'N', then Z is not referenced.
  //C
  //C  LDZ     (input) INTEGER
  //C          The leading dimension of the array Z.  LDZ >= 1, and if
  //C          eigenvectors are desired, then  LDZ >= max(1,N).
  //C
  //C  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
  //C          If COMPZ = 'N', then WORK is not referenced.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value
  //C          > 0:  the algorithm has failed to find all the eigenvalues in
  //C                a total of 30*N iterations; if INFO = i, then i
  //C                elements of E have not converged to zero; on exit, D
  //C                and E contain the elements of a symmetric tridiagonal
  //C                matrix which is orthogonally similar to the original
  //C                matrix.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters.
  //C
  info = 0;
  //C
  if (lsame(compz, "N")) {
    icompz = 0;
  }
  else if (lsame(compz, "V")) {
    icompz = 1;
  }
  else if (lsame(compz, "I")) {
    icompz = 2;
  }
  else {
    icompz = -1;
  }
  if (icompz < 0) {
    info = -1;
  }
  else if (n < 0) {
    info = -2;
  }
  else if ((ldz < 1) || (icompz > 0 && ldz < fem::max(1, n))) {
    info = -6;
  }
  if (info != 0) {
    xerbla("DSTEQR", -info);
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (n == 0) {
    return;
  }
  //C
  if (n == 1) {
    if (icompz == 2) {
      z(1, 1) = one;
    }
    return;
  }
  //C
  //C     Determine the unit roundoff and over/underflow thresholds.
  //C
  eps = dlamch(cmn, "E");
  eps2 = fem::pow2(eps);
  safmin = dlamch(cmn, "S");
  safmax = one / safmin;
  ssfmax = fem::sqrt(safmax) / three;
  ssfmin = fem::sqrt(safmin) / eps2;
  //C
  //C     Compute the eigenvalues and eigenvectors of the tridiagonal
  //C     matrix.
  //C
  if (icompz == 2) {
    dlaset("Full", n, n, zero, one, z, ldz);
  }
  //C
  nmaxit = n * maxit;
  jtot = 0;
  //C
  //C     Determine where the matrix splits and choose QL or QR iteration
  //C     for each block, according to whether top or bottom diagonal
  //C     element is smaller.
  //C
  l1 = 1;
  nm1 = n - 1;
  //C
  statement_10:
  if (l1 > n) {
    goto statement_160;
  }
  if (l1 > 1) {
    e(l1 - 1) = zero;
  }
  if (l1 <= nm1) {
    FEM_DO(m, l1, nm1) {
      tst = fem::abs(e(m));
      if (tst == zero) {
        goto statement_30;
      }
      if (tst <= (fem::sqrt(fem::abs(d(m))) * fem::sqrt(fem::abs(d(
          m + 1)))) * eps) {
        e(m) = zero;
        goto statement_30;
      }
    }
  }
  m = n;
  //C
  statement_30:
  l = l1;
  lsv = l;
  lend = m;
  lendsv = lend;
  l1 = m + 1;
  if (lend == l) {
    goto statement_10;
  }
  //C
  //C     Scale submatrix in rows and columns L to LEND
  //C
  anorm = dlanst("I", lend - l + 1, d(l), e(l));
  iscale = 0;
  if (anorm == zero) {
    goto statement_10;
  }
  if (anorm > ssfmax) {
    iscale = 1;
    dlascl(cmn, "G", 0, 0, anorm, ssfmax, lend - l + 1, 1, d(l), n, info);
    dlascl(cmn, "G", 0, 0, anorm, ssfmax, lend - l, 1, e(l), n, info);
  }
  else if (anorm < ssfmin) {
    iscale = 2;
    dlascl(cmn, "G", 0, 0, anorm, ssfmin, lend - l + 1, 1, d(l), n, info);
    dlascl(cmn, "G", 0, 0, anorm, ssfmin, lend - l, 1, e(l), n, info);
  }
  //C
  //C     Choose between QL and QR iteration
  //C
  if (fem::abs(d(lend)) < fem::abs(d(l))) {
    lend = lsv;
    l = lendsv;
  }
  //C
  if (lend > l) {
    //C
    //C        QL Iteration
    //C
    //C        Look for small subdiagonal element.
    //C
    statement_40:
    if (l != lend) {
      lendm1 = lend - 1;
      FEM_DO(m, l, lendm1) {
        tst = fem::pow2(fem::abs(e(m)));
        if (tst <= (eps2 * fem::abs(d(m))) * fem::abs(d(m + 1)) + safmin) {
          goto statement_60;
        }
      }
    }
    //C
    m = lend;
    //C
    statement_60:
    if (m < lend) {
      e(m) = zero;
    }
    p = d(l);
    if (m == l) {
      goto statement_80;
    }
    //C
    //C        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
    //C        to compute its eigensystem.
    //C
    if (m == l + 1) {
      if (icompz > 0) {
        dlaev2(d(l), e(l), d(l + 1), rt1, rt2, c, s);
        work(l) = c;
        work(n - 1 + l) = s;
        dlasr("R", "V", "B", n, 2, work(l), work(n - 1 + l), z(1, l), ldz);
      }
      else {
        dlae2(d(l), e(l), d(l + 1), rt1, rt2);
      }
      d(l) = rt1;
      d(l + 1) = rt2;
      e(l) = zero;
      l += 2;
      if (l <= lend) {
        goto statement_40;
      }
      goto statement_140;
    }
    //C
    if (jtot == nmaxit) {
      goto statement_140;
    }
    jtot++;
    //C
    //C        Form shift.
    //C
    g = (d(l + 1) - p) / (two * e(l));
    r = dlapy2(g, one);
    g = d(m) - p + (e(l) / (g + fem::sign(r, g)));
    //C
    s = one;
    c = one;
    p = zero;
    //C
    //C        Inner loop
    //C
    mm1 = m - 1;
    FEM_DOSTEP(i, mm1, l, -1) {
      f = s * e(i);
      b = c * e(i);
      dlartg(cmn, g, f, c, s, r);
      if (i != m - 1) {
        e(i + 1) = r;
      }
      g = d(i + 1) - p;
      r = (d(i) - g) * s + two * c * b;
      p = s * r;
      d(i + 1) = g + p;
      g = c * r - b;
      //C
      //C           If eigenvectors are desired, then save rotations.
      //C
      if (icompz > 0) {
        work(i) = c;
        work(n - 1 + i) = -s;
      }
      //C
    }
    //C
    //C        If eigenvectors are desired, then apply saved rotations.
    //C
    if (icompz > 0) {
      mm = m - l + 1;
      dlasr("R", "V", "B", n, mm, work(l), work(n - 1 + l), z(1, l), ldz);
    }
    //C
    d(l) = d(l) - p;
    e(l) = g;
    goto statement_40;
    //C
    //C        Eigenvalue found.
    //C
    statement_80:
    d(l) = p;
    //C
    l++;
    if (l <= lend) {
      goto statement_40;
    }
    goto statement_140;
    //C
  }
  else {
    //C
    //C        QR Iteration
    //C
    //C        Look for small superdiagonal element.
    //C
    statement_90:
    if (l != lend) {
      lendp1 = lend + 1;
      FEM_DOSTEP(m, l, lendp1, -1) {
        tst = fem::pow2(fem::abs(e(m - 1)));
        if (tst <= (eps2 * fem::abs(d(m))) * fem::abs(d(m - 1)) + safmin) {
          goto statement_110;
        }
      }
    }
    //C
    m = lend;
    //C
    statement_110:
    if (m > lend) {
      e(m - 1) = zero;
    }
    p = d(l);
    if (m == l) {
      goto statement_130;
    }
    //C
    //C        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
    //C        to compute its eigensystem.
    //C
    if (m == l - 1) {
      if (icompz > 0) {
        dlaev2(d(l - 1), e(l - 1), d(l), rt1, rt2, c, s);
        work(m) = c;
        work(n - 1 + m) = s;
        dlasr("R", "V", "F", n, 2, work(m), work(n - 1 + m), z(1, l - 1), ldz);
      }
      else {
        dlae2(d(l - 1), e(l - 1), d(l), rt1, rt2);
      }
      d(l - 1) = rt1;
      d(l) = rt2;
      e(l - 1) = zero;
      l = l - 2;
      if (l >= lend) {
        goto statement_90;
      }
      goto statement_140;
    }
    //C
    if (jtot == nmaxit) {
      goto statement_140;
    }
    jtot++;
    //C
    //C        Form shift.
    //C
    g = (d(l - 1) - p) / (two * e(l - 1));
    r = dlapy2(g, one);
    g = d(m) - p + (e(l - 1) / (g + fem::sign(r, g)));
    //C
    s = one;
    c = one;
    p = zero;
    //C
    //C        Inner loop
    //C
    lm1 = l - 1;
    FEM_DO(i, m, lm1) {
      f = s * e(i);
      b = c * e(i);
      dlartg(cmn, g, f, c, s, r);
      if (i != m) {
        e(i - 1) = r;
      }
      g = d(i) - p;
      r = (d(i + 1) - g) * s + two * c * b;
      p = s * r;
      d(i) = g + p;
      g = c * r - b;
      //C
      //C           If eigenvectors are desired, then save rotations.
      //C
      if (icompz > 0) {
        work(i) = c;
        work(n - 1 + i) = s;
      }
      //C
    }
    //C
    //C        If eigenvectors are desired, then apply saved rotations.
    //C
    if (icompz > 0) {
      mm = l - m + 1;
      dlasr("R", "V", "F", n, mm, work(m), work(n - 1 + m), z(1, m), ldz);
    }
    //C
    d(l) = d(l) - p;
    e(lm1) = g;
    goto statement_90;
    //C
    //C        Eigenvalue found.
    //C
    statement_130:
    d(l) = p;
    //C
    l = l - 1;
    if (l >= lend) {
      goto statement_90;
    }
    goto statement_140;
    //C
  }
  //C
  //C     Undo scaling if necessary
  //C
  statement_140:
  if (iscale == 1) {
    dlascl(cmn, "G", 0, 0, ssfmax, anorm, lendsv - lsv + 1, 1, d(lsv), n, info);
    dlascl(cmn, "G", 0, 0, ssfmax, anorm, lendsv - lsv, 1, e(lsv), n, info);
  }
  else if (iscale == 2) {
    dlascl(cmn, "G", 0, 0, ssfmin, anorm, lendsv - lsv + 1, 1, d(lsv), n, info);
    dlascl(cmn, "G", 0, 0, ssfmin, anorm, lendsv - lsv, 1, e(lsv), n, info);
  }
  //C
  //C     Check for no convergence to an eigenvalue after a total
  //C     of N*MAXIT iterations.
  //C
  if (jtot < nmaxit) {
    goto statement_10;
  }
  {
    int fem_do_last = n - 1;
    FEM_DO(i, 1, fem_do_last) {
      if (e(i) != zero) {
        info++;
      }
    }
  }
  goto statement_190;
  //C
  //C     Order eigenvalues and eigenvectors.
  //C
  statement_160:
  if (icompz == 0) {
    //C
    //C        Use Quick Sort
    //C
    dlasrt("I", n, d, info);
    //C
  }
  else {
    //C
    //C        Use Selection Sort to minimize swaps of eigenvectors
    //C
    FEM_DO(ii, 2, n) {
      i = ii - 1;
      k = i;
      p = d(i);
      FEM_DO(j, ii, n) {
        if (d(j) < p) {
          k = j;
          p = d(j);
        }
      }
      if (k != i) {
        d(k) = d(i);
        d(i) = p;
        dswap(n, z(1, i), 1, z(1, k), 1);
      }
    }
  }
  //C
  statement_190:;
  //C
  //C     End of DSTEQR
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dsterf.f
inline
void
dsterf(
  common& cmn,
  int const& n,
  arr_ref<double> d,
  arr_ref<double> e,
  int& info)
{
  d(dimension(star));
  e(dimension(star));
  double eps = fem::double0;
  double eps2 = fem::double0;
  double safmin = fem::double0;
  const double one = 1.0e0;
  double safmax = fem::double0;
  const double three = 3.0e0;
  double ssfmax = fem::double0;
  double ssfmin = fem::double0;
  const int maxit = 30;
  int nmaxit = fem::int0;
  const double zero = 0.0e0;
  double sigma = fem::double0;
  int jtot = fem::int0;
  int l1 = fem::int0;
  int m = fem::int0;
  int l = fem::int0;
  int lsv = fem::int0;
  int lend = fem::int0;
  int lendsv = fem::int0;
  double anorm = fem::double0;
  int iscale = fem::int0;
  int i = fem::int0;
  double p = fem::double0;
  double rte = fem::double0;
  double rt1 = fem::double0;
  double rt2 = fem::double0;
  const double two = 2.0e0;
  double r = fem::double0;
  double c = fem::double0;
  double s = fem::double0;
  double gamma = fem::double0;
  double bb = fem::double0;
  double oldc = fem::double0;
  double oldgam = fem::double0;
  double alpha = fem::double0;
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
  //C  using the Pal-Walker-Kahan variant of the QL or QR algorithm.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  N       (input) INTEGER
  //C          The order of the matrix.  N >= 0.
  //C
  //C  D       (input/output) DOUBLE PRECISION array, dimension (N)
  //C          On entry, the n diagonal elements of the tridiagonal matrix.
  //C          On exit, if INFO = 0, the eigenvalues in ascending order.
  //C
  //C  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
  //C          On entry, the (n-1) subdiagonal elements of the tridiagonal
  //C          matrix.
  //C          On exit, E has been destroyed.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value
  //C          > 0:  the algorithm failed to find all of the eigenvalues in
  //C                a total of 30*N iterations; if INFO = i, then i
  //C                elements of E have not converged to zero.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters.
  //C
  info = 0;
  //C
  //C     Quick return if possible
  //C
  if (n < 0) {
    info = -1;
    xerbla("DSTERF", -info);
    return;
  }
  if (n <= 1) {
    return;
  }
  //C
  //C     Determine the unit roundoff for this environment.
  //C
  eps = dlamch(cmn, "E");
  eps2 = fem::pow2(eps);
  safmin = dlamch(cmn, "S");
  safmax = one / safmin;
  ssfmax = fem::sqrt(safmax) / three;
  ssfmin = fem::sqrt(safmin) / eps2;
  //C
  //C     Compute the eigenvalues of the tridiagonal matrix.
  //C
  nmaxit = n * maxit;
  sigma = zero;
  jtot = 0;
  //C
  //C     Determine where the matrix splits and choose QL or QR iteration
  //C     for each block, according to whether top or bottom diagonal
  //C     element is smaller.
  //C
  l1 = 1;
  //C
  statement_10:
  if (l1 > n) {
    goto statement_170;
  }
  if (l1 > 1) {
    e(l1 - 1) = zero;
  }
  {
    int fem_do_last = n - 1;
    FEM_DO(m, l1, fem_do_last) {
      if (fem::abs(e(m)) <= (fem::sqrt(fem::abs(d(m))) * fem::sqrt(
          fem::abs(d(m + 1)))) * eps) {
        e(m) = zero;
        goto statement_30;
      }
    }
  }
  m = n;
  //C
  statement_30:
  l = l1;
  lsv = l;
  lend = m;
  lendsv = lend;
  l1 = m + 1;
  if (lend == l) {
    goto statement_10;
  }
  //C
  //C     Scale submatrix in rows and columns L to LEND
  //C
  anorm = dlanst("I", lend - l + 1, d(l), e(l));
  iscale = 0;
  if (anorm > ssfmax) {
    iscale = 1;
    dlascl(cmn, "G", 0, 0, anorm, ssfmax, lend - l + 1, 1, d(l), n, info);
    dlascl(cmn, "G", 0, 0, anorm, ssfmax, lend - l, 1, e(l), n, info);
  }
  else if (anorm < ssfmin) {
    iscale = 2;
    dlascl(cmn, "G", 0, 0, anorm, ssfmin, lend - l + 1, 1, d(l), n, info);
    dlascl(cmn, "G", 0, 0, anorm, ssfmin, lend - l, 1, e(l), n, info);
  }
  //C
  {
    int fem_do_last = lend - 1;
    FEM_DO(i, l, fem_do_last) {
      e(i) = fem::pow2(e(i));
    }
  }
  //C
  //C     Choose between QL and QR iteration
  //C
  if (fem::abs(d(lend)) < fem::abs(d(l))) {
    lend = lsv;
    l = lendsv;
  }
  //C
  if (lend >= l) {
    //C
    //C        QL Iteration
    //C
    //C        Look for small subdiagonal element.
    //C
    statement_50:
    if (l != lend) {
      {
        int fem_do_last = lend - 1;
        FEM_DO(m, l, fem_do_last) {
          if (fem::abs(e(m)) <= eps2 * fem::abs(d(m) * d(m + 1))) {
            goto statement_70;
          }
        }
      }
    }
    m = lend;
    //C
    statement_70:
    if (m < lend) {
      e(m) = zero;
    }
    p = d(l);
    if (m == l) {
      goto statement_90;
    }
    //C
    //C        If remaining matrix is 2 by 2, use DLAE2 to compute its
    //C        eigenvalues.
    //C
    if (m == l + 1) {
      rte = fem::sqrt(e(l));
      dlae2(d(l), rte, d(l + 1), rt1, rt2);
      d(l) = rt1;
      d(l + 1) = rt2;
      e(l) = zero;
      l += 2;
      if (l <= lend) {
        goto statement_50;
      }
      goto statement_150;
    }
    //C
    if (jtot == nmaxit) {
      goto statement_150;
    }
    jtot++;
    //C
    //C        Form shift.
    //C
    rte = fem::sqrt(e(l));
    sigma = (d(l + 1) - p) / (two * rte);
    r = dlapy2(sigma, one);
    sigma = p - (rte / (sigma + fem::sign(r, sigma)));
    //C
    c = one;
    s = zero;
    gamma = d(m) - sigma;
    p = gamma * gamma;
    //C
    //C        Inner loop
    //C
    FEM_DOSTEP(i, m - 1, l, -1) {
      bb = e(i);
      r = p + bb;
      if (i != m - 1) {
        e(i + 1) = s * r;
      }
      oldc = c;
      c = p / r;
      s = bb / r;
      oldgam = gamma;
      alpha = d(i);
      gamma = c * (alpha - sigma) - s * oldgam;
      d(i + 1) = oldgam + (alpha - gamma);
      if (c != zero) {
        p = (gamma * gamma) / c;
      }
      else {
        p = oldc * bb;
      }
    }
    //C
    e(l) = s * p;
    d(l) = sigma + gamma;
    goto statement_50;
    //C
    //C        Eigenvalue found.
    //C
    statement_90:
    d(l) = p;
    //C
    l++;
    if (l <= lend) {
      goto statement_50;
    }
    goto statement_150;
    //C
  }
  else {
    //C
    //C        QR Iteration
    //C
    //C        Look for small superdiagonal element.
    //C
    statement_100:
    FEM_DOSTEP(m, l, lend + 1, -1) {
      if (fem::abs(e(m - 1)) <= eps2 * fem::abs(d(m) * d(m - 1))) {
        goto statement_120;
      }
    }
    m = lend;
    //C
    statement_120:
    if (m > lend) {
      e(m - 1) = zero;
    }
    p = d(l);
    if (m == l) {
      goto statement_140;
    }
    //C
    //C        If remaining matrix is 2 by 2, use DLAE2 to compute its
    //C        eigenvalues.
    //C
    if (m == l - 1) {
      rte = fem::sqrt(e(l - 1));
      dlae2(d(l), rte, d(l - 1), rt1, rt2);
      d(l) = rt1;
      d(l - 1) = rt2;
      e(l - 1) = zero;
      l = l - 2;
      if (l >= lend) {
        goto statement_100;
      }
      goto statement_150;
    }
    //C
    if (jtot == nmaxit) {
      goto statement_150;
    }
    jtot++;
    //C
    //C        Form shift.
    //C
    rte = fem::sqrt(e(l - 1));
    sigma = (d(l - 1) - p) / (two * rte);
    r = dlapy2(sigma, one);
    sigma = p - (rte / (sigma + fem::sign(r, sigma)));
    //C
    c = one;
    s = zero;
    gamma = d(m) - sigma;
    p = gamma * gamma;
    //C
    //C        Inner loop
    //C
    {
      int fem_do_last = l - 1;
      FEM_DO(i, m, fem_do_last) {
        bb = e(i);
        r = p + bb;
        if (i != m) {
          e(i - 1) = s * r;
        }
        oldc = c;
        c = p / r;
        s = bb / r;
        oldgam = gamma;
        alpha = d(i + 1);
        gamma = c * (alpha - sigma) - s * oldgam;
        d(i) = oldgam + (alpha - gamma);
        if (c != zero) {
          p = (gamma * gamma) / c;
        }
        else {
          p = oldc * bb;
        }
      }
    }
    //C
    e(l - 1) = s * p;
    d(l) = sigma + gamma;
    goto statement_100;
    //C
    //C        Eigenvalue found.
    //C
    statement_140:
    d(l) = p;
    //C
    l = l - 1;
    if (l >= lend) {
      goto statement_100;
    }
    goto statement_150;
    //C
  }
  //C
  //C     Undo scaling if necessary
  //C
  statement_150:
  if (iscale == 1) {
    dlascl(cmn, "G", 0, 0, ssfmax, anorm, lendsv - lsv + 1, 1, d(lsv), n, info);
  }
  if (iscale == 2) {
    dlascl(cmn, "G", 0, 0, ssfmin, anorm, lendsv - lsv + 1, 1, d(lsv), n, info);
  }
  //C
  //C     Check for no convergence to an eigenvalue after a total
  //C     of N*MAXIT iterations.
  //C
  if (jtot < nmaxit) {
    goto statement_10;
  }
  {
    int fem_do_last = n - 1;
    FEM_DO(i, 1, fem_do_last) {
      if (e(i) != zero) {
        info++;
      }
    }
  }
  goto statement_180;
  //C
  //C     Sort eigenvalues in increasing order.
  //C
  statement_170:
  dlasrt("I", n, d, info);
  //C
  statement_180:;
  //C
  //C     End of DSTERF
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dsytd2.f
inline
void
dsytd2(
  common& cmn,
  str_cref uplo,
  int const& n,
  arr_ref<double, 2> a,
  int const& lda,
  arr_ref<double> d,
  arr_ref<double> e,
  arr_ref<double> tau,
  int& info)
{
  a(dimension(lda, star));
  d(dimension(star));
  e(dimension(star));
  tau(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal
  //C  form T by an orthogonal similarity transformation: Q' * A * Q = T.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  UPLO    (input) CHARACTER*1
  //C          Specifies whether the upper or lower triangular part of the
  //C          symmetric matrix A is stored:
  //C          = 'U':  Upper triangular
  //C          = 'L':  Lower triangular
  //C
  //C  N       (input) INTEGER
  //C          The order of the matrix A.  N >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  //C          n-by-n upper triangular part of A contains the upper
  //C          triangular part of the matrix A, and the strictly lower
  //C          triangular part of A is not referenced.  If UPLO = 'L', the
  //C          leading n-by-n lower triangular part of A contains the lower
  //C          triangular part of the matrix A, and the strictly upper
  //C          triangular part of A is not referenced.
  //C          On exit, if UPLO = 'U', the diagonal and first superdiagonal
  //C          of A are overwritten by the corresponding elements of the
  //C          tridiagonal matrix T, and the elements above the first
  //C          superdiagonal, with the array TAU, represent the orthogonal
  //C          matrix Q as a product of elementary reflectors; if UPLO
  //C          = 'L', the diagonal and first subdiagonal of A are over-
  //C          written by the corresponding elements of the tridiagonal
  //C          matrix T, and the elements below the first subdiagonal, with
  //C          the array TAU, represent the orthogonal matrix Q as a product
  //C          of elementary reflectors. See Further Details.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(1,N).
  //C
  //C  D       (output) DOUBLE PRECISION array, dimension (N)
  //C          The diagonal elements of the tridiagonal matrix T:
  //C          D(i) = A(i,i).
  //C
  //C  E       (output) DOUBLE PRECISION array, dimension (N-1)
  //C          The off-diagonal elements of the tridiagonal matrix T:
  //C          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
  //C
  //C  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
  //C          The scalar factors of the elementary reflectors (see Further
  //C          Details).
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value.
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  If UPLO = 'U', the matrix Q is represented as a product of elementary
  //C  reflectors
  //C
  //C     Q = H(n-1) . . . H(2) H(1).
  //C
  //C  Each H(i) has the form
  //C
  //C     H(i) = I - tau * v * v'
  //C
  //C  where tau is a real scalar, and v is a real vector with
  //C  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
  //C  A(1:i-1,i+1), and tau in TAU(i).
  //C
  //C  If UPLO = 'L', the matrix Q is represented as a product of elementary
  //C  reflectors
  //C
  //C     Q = H(1) H(2) . . . H(n-1).
  //C
  //C  Each H(i) has the form
  //C
  //C     H(i) = I - tau * v * v'
  //C
  //C  where tau is a real scalar, and v is a real vector with
  //C  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
  //C  and tau in TAU(i).
  //C
  //C  The contents of A on exit are illustrated by the following examples
  //C  with n = 5:
  //C
  //C  if UPLO = 'U':                       if UPLO = 'L':
  //C
  //C    (  d   e   v2  v3  v4 )              (  d                  )
  //C    (      d   e   v3  v4 )              (  e   d              )
  //C    (          d   e   v4 )              (  v1  e   d          )
  //C    (              d   e  )              (  v1  v2  e   d      )
  //C    (                  d  )              (  v1  v2  v3  e   d  )
  //C
  //C  where d and e denote diagonal and off-diagonal elements of T, and vi
  //C  denotes an element of the vector defining H(i).
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters
  //C
  info = 0;
  bool upper = lsame(uplo, "U");
  if (!upper && !lsame(uplo, "L")) {
    info = -1;
  }
  else if (n < 0) {
    info = -2;
  }
  else if (lda < fem::max(1, n)) {
    info = -4;
  }
  if (info != 0) {
    xerbla("DSYTD2", -info);
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (n <= 0) {
    return;
  }
  //C
  int i = fem::int0;
  double taui = fem::double0;
  const double zero = 0.0e0;
  const double one = 1.0e0;
  const double half = 1.0e0 / 2.0e0;
  double alpha = fem::double0;
  if (upper) {
    //C
    //C        Reduce the upper triangle of A
    //C
    FEM_DOSTEP(i, n - 1, 1, -1) {
      //C
      //C           Generate elementary reflector H(i) = I - tau * v * v'
      //C           to annihilate A(1:i-1,i+1)
      //C
      dlarfg(cmn, i, a(i, i + 1), a(1, i + 1), 1, taui);
      e(i) = a(i, i + 1);
      //C
      if (taui != zero) {
        //C
        //C              Apply H(i) from both sides to A(1:i,1:i)
        //C
        a(i, i + 1) = one;
        //C
        //C              Compute  x := tau * A * v  storing x in TAU(1:i)
        //C
        dsymv(uplo, i, taui, a, lda, a(1, i + 1), 1, zero, tau, 1);
        //C
        //C              Compute  w := x - 1/2 * tau * (x'*v) * v
        //C
        alpha = -half * taui * ddot(i, tau, 1, a(1, i + 1), 1);
        daxpy(i, alpha, a(1, i + 1), 1, tau, 1);
        //C
        //C              Apply the transformation as a rank-2 update:
        //C                 A := A - v * w' - w * v'
        //C
        dsyr2(uplo, i, -one, a(1, i + 1), 1, tau, 1, a, lda);
        //C
        a(i, i + 1) = e(i);
      }
      d(i + 1) = a(i + 1, i + 1);
      tau(i) = taui;
    }
    d(1) = a(1, 1);
  }
  else {
    //C
    //C        Reduce the lower triangle of A
    //C
    {
      int fem_do_last = n - 1;
      FEM_DO(i, 1, fem_do_last) {
        //C
        //C           Generate elementary reflector H(i) = I - tau * v * v'
        //C           to annihilate A(i+2:n,i)
        //C
        dlarfg(cmn, n - i, a(i + 1, i), a(fem::min(i + 2, n), i), 1, taui);
        e(i) = a(i + 1, i);
        //C
        if (taui != zero) {
          //C
          //C              Apply H(i) from both sides to A(i+1:n,i+1:n)
          //C
          a(i + 1, i) = one;
          //C
          //C              Compute  x := tau * A * v  storing y in TAU(i:n-1)
          //C
          dsymv(uplo, n - i, taui, a(i + 1, i + 1), lda, a(i + 1, i),
            1, zero, tau(i), 1);
          //C
          //C              Compute  w := x - 1/2 * tau * (x'*v) * v
          //C
          alpha = -half * taui * ddot(n - i, tau(i), 1, a(i + 1, i), 1);
          daxpy(n - i, alpha, a(i + 1, i), 1, tau(i), 1);
          //C
          //C              Apply the transformation as a rank-2 update:
          //C                 A := A - v * w' - w * v'
          //C
          dsyr2(uplo, n - i, -one, a(i + 1, i), 1, tau(i), 1, a(i + 1,
            i + 1), lda);
          //C
          a(i + 1, i) = e(i);
        }
        d(i) = a(i, i);
        tau(i) = taui;
      }
    }
    d(n) = a(n, n);
  }
  //C
  //C     End of DSYTD2
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dsytrd.f
inline
void
dsytrd(
  common& cmn,
  str_cref uplo,
  int const& n,
  arr_ref<double, 2> a,
  int const& lda,
  arr_ref<double> d,
  arr_ref<double> e,
  arr_ref<double> tau,
  arr_ref<double> work,
  int const& lwork,
  int& info)
{
  a(dimension(lda, star));
  d(dimension(star));
  e(dimension(star));
  tau(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DSYTRD reduces a real symmetric matrix A to real symmetric
  //C  tridiagonal form T by an orthogonal similarity transformation:
  //C  Q**T * A * Q = T.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  UPLO    (input) CHARACTER*1
  //C          = 'U':  Upper triangle of A is stored;
  //C          = 'L':  Lower triangle of A is stored.
  //C
  //C  N       (input) INTEGER
  //C          The order of the matrix A.  N >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
  //C          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  //C          N-by-N upper triangular part of A contains the upper
  //C          triangular part of the matrix A, and the strictly lower
  //C          triangular part of A is not referenced.  If UPLO = 'L', the
  //C          leading N-by-N lower triangular part of A contains the lower
  //C          triangular part of the matrix A, and the strictly upper
  //C          triangular part of A is not referenced.
  //C          On exit, if UPLO = 'U', the diagonal and first superdiagonal
  //C          of A are overwritten by the corresponding elements of the
  //C          tridiagonal matrix T, and the elements above the first
  //C          superdiagonal, with the array TAU, represent the orthogonal
  //C          matrix Q as a product of elementary reflectors; if UPLO
  //C          = 'L', the diagonal and first subdiagonal of A are over-
  //C          written by the corresponding elements of the tridiagonal
  //C          matrix T, and the elements below the first subdiagonal, with
  //C          the array TAU, represent the orthogonal matrix Q as a product
  //C          of elementary reflectors. See Further Details.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(1,N).
  //C
  //C  D       (output) DOUBLE PRECISION array, dimension (N)
  //C          The diagonal elements of the tridiagonal matrix T:
  //C          D(i) = A(i,i).
  //C
  //C  E       (output) DOUBLE PRECISION array, dimension (N-1)
  //C          The off-diagonal elements of the tridiagonal matrix T:
  //C          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
  //C
  //C  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
  //C          The scalar factors of the elementary reflectors (see Further
  //C          Details).
  //C
  //C  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  //C
  //C  LWORK   (input) INTEGER
  //C          The dimension of the array WORK.  LWORK >= 1.
  //C          For optimum performance LWORK >= N*NB, where NB is the
  //C          optimal blocksize.
  //C
  //C          If LWORK = -1, then a workspace query is assumed; the routine
  //C          only calculates the optimal size of the WORK array, returns
  //C          this value as the first entry of the WORK array, and no error
  //C          message related to LWORK is issued by XERBLA.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value
  //C
  //C  Further Details
  //C  ===============
  //C
  //C  If UPLO = 'U', the matrix Q is represented as a product of elementary
  //C  reflectors
  //C
  //C     Q = H(n-1) . . . H(2) H(1).
  //C
  //C  Each H(i) has the form
  //C
  //C     H(i) = I - tau * v * v'
  //C
  //C  where tau is a real scalar, and v is a real vector with
  //C  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
  //C  A(1:i-1,i+1), and tau in TAU(i).
  //C
  //C  If UPLO = 'L', the matrix Q is represented as a product of elementary
  //C  reflectors
  //C
  //C     Q = H(1) H(2) . . . H(n-1).
  //C
  //C  Each H(i) has the form
  //C
  //C     H(i) = I - tau * v * v'
  //C
  //C  where tau is a real scalar, and v is a real vector with
  //C  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
  //C  and tau in TAU(i).
  //C
  //C  The contents of A on exit are illustrated by the following examples
  //C  with n = 5:
  //C
  //C  if UPLO = 'U':                       if UPLO = 'L':
  //C
  //C    (  d   e   v2  v3  v4 )              (  d                  )
  //C    (      d   e   v3  v4 )              (  e   d              )
  //C    (          d   e   v4 )              (  v1  e   d          )
  //C    (              d   e  )              (  v1  v2  e   d      )
  //C    (                  d  )              (  v1  v2  v3  e   d  )
  //C
  //C  where d and e denote diagonal and off-diagonal elements of T, and vi
  //C  denotes an element of the vector defining H(i).
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters
  //C
  info = 0;
  bool upper = lsame(uplo, "U");
  bool lquery = (lwork ==  - 1);
  if (!upper && !lsame(uplo, "L")) {
    info = -1;
  }
  else if (n < 0) {
    info = -2;
  }
  else if (lda < fem::max(1, n)) {
    info = -4;
  }
  else if (lwork < 1 && !lquery) {
    info = -9;
  }
  //C
  int nb = fem::int0;
  int lwkopt = fem::int0;
  if (info == 0) {
    //C
    //C        Determine the block size.
    //C
    nb = ilaenv(1, "DSYTRD", uplo, n, -1, -1, -1);
    lwkopt = n * nb;
    work(1) = lwkopt;
  }
  //C
  if (info != 0) {
    xerbla("DSYTRD", -info);
    return;
  }
  else if (lquery) {
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (n == 0) {
    work(1) = 1;
    return;
  }
  //C
  int nx = n;
  int iws = 1;
  int ldwork = fem::int0;
  int nbmin = fem::int0;
  if (nb > 1 && nb < n) {
    //C
    //C        Determine when to cross over from blocked to unblocked code
    //C        (last block is always handled by unblocked code).
    //C
    nx = fem::max(nb, ilaenv(3, "DSYTRD", uplo, n, -1, -1, -1));
    if (nx < n) {
      //C
      //C           Determine if workspace is large enough for blocked code.
      //C
      ldwork = n;
      iws = ldwork * nb;
      if (lwork < iws) {
        //C
        //C              Not enough workspace to use optimal NB:  determine the
        //C              minimum value of NB, and reduce NB or force use of
        //C              unblocked code by setting NX = N.
        //C
        nb = fem::max(lwork / ldwork, 1);
        nbmin = ilaenv(2, "DSYTRD", uplo, n, -1, -1, -1);
        if (nb < nbmin) {
          nx = n;
        }
      }
    }
    else {
      nx = n;
    }
  }
  else {
    nb = 1;
  }
  //C
  int kk = fem::int0;
  int i = fem::int0;
  const double one = 1.0e+0;
  int j = fem::int0;
  int iinfo = fem::int0;
  if (upper) {
    //C
    //C        Reduce the upper triangle of A.
    //C        Columns 1:kk are handled by the unblocked method.
    //C
    kk = n - ((n - nx + nb - 1) / nb) * nb;
    FEM_DOSTEP(i, n - nb + 1, kk + 1, -nb) {
      //C
      //C           Reduce columns i:i+nb-1 to tridiagonal form and form the
      //C           matrix W which is needed to update the unreduced part of
      //C           the matrix
      //C
      dlatrd(cmn, uplo, i + nb - 1, nb, a, lda, e, tau, work, ldwork);
      //C
      //C           Update the unreduced submatrix A(1:i-1,1:i-1), using an
      //C           update of the form:  A := A - V*W' - W*V'
      //C
      dsyr2k(uplo, "No transpose", i - 1, nb, -one, a(1, i), lda,
        work, ldwork, one, a, lda);
      //C
      //C           Copy superdiagonal elements back into A, and diagonal
      //C           elements into D
      //C
      {
        int fem_do_last = i + nb - 1;
        FEM_DO(j, i, fem_do_last) {
          a(j - 1, j) = e(j - 1);
          d(j) = a(j, j);
        }
      }
    }
    //C
    //C        Use unblocked code to reduce the last or only block
    //C
    dsytd2(cmn, uplo, kk, a, lda, d, e, tau, iinfo);
  }
  else {
    //C
    //C        Reduce the lower triangle of A
    //C
    FEM_DOSTEP(i, 1, n - nx, nb) {
      //C
      //C           Reduce columns i:i+nb-1 to tridiagonal form and form the
      //C           matrix W which is needed to update the unreduced part of
      //C           the matrix
      //C
      dlatrd(cmn, uplo, n - i + 1, nb, a(i, i), lda, e(i), tau(i),
        work, ldwork);
      //C
      //C           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
      //C           an update of the form:  A := A - V*W' - W*V'
      //C
      dsyr2k(uplo, "No transpose", n - i - nb + 1, nb, -one, a(i + nb,
        i), lda, work(nb + 1), ldwork, one, a(i + nb, i + nb), lda);
      //C
      //C           Copy subdiagonal elements back into A, and diagonal
      //C           elements into D
      //C
      {
        int fem_do_last = i + nb - 1;
        FEM_DO(j, i, fem_do_last) {
          a(j + 1, j) = e(j);
          d(j) = a(j, j);
        }
      }
    }
    //C
    //C        Use unblocked code to reduce the last or only block
    //C
    dsytd2(cmn, uplo, n - i + 1, a(i, i), lda, d(i), e(i), tau(i), iinfo);
  }
  //C
  work(1) = lwkopt;
  //C
  //C     End of DSYTRD
  //C
}

// Fortran file: /net/sigma/raid1/rwgk/lapack-3.2.1/SRC/dsyev.f
inline
void
dsyev(
  common& cmn,
  str_cref jobz,
  str_cref uplo,
  int const& n,
  arr_ref<double, 2> a,
  int const& lda,
  arr_ref<double> w,
  arr_ref<double> work,
  int const& lwork,
  int& info)
{
  a(dimension(lda, star));
  w(dimension(star));
  work(dimension(star));
  //C
  //C  -- LAPACK driver routine (version 3.2) --
  //C  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  //C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  //C     November 2006
  //C
  //C     .. Scalar Arguments ..
  //C     ..
  //C     .. Array Arguments ..
  //C     ..
  //C
  //C  Purpose
  //C  =======
  //C
  //C  DSYEV computes all eigenvalues and, optionally, eigenvectors of a
  //C  real symmetric matrix A.
  //C
  //C  Arguments
  //C  =========
  //C
  //C  JOBZ    (input) CHARACTER*1
  //C          = 'N':  Compute eigenvalues only;
  //C          = 'V':  Compute eigenvalues and eigenvectors.
  //C
  //C  UPLO    (input) CHARACTER*1
  //C          = 'U':  Upper triangle of A is stored;
  //C          = 'L':  Lower triangle of A is stored.
  //C
  //C  N       (input) INTEGER
  //C          The order of the matrix A.  N >= 0.
  //C
  //C  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
  //C          On entry, the symmetric matrix A.  If UPLO = 'U', the
  //C          leading N-by-N upper triangular part of A contains the
  //C          upper triangular part of the matrix A.  If UPLO = 'L',
  //C          the leading N-by-N lower triangular part of A contains
  //C          the lower triangular part of the matrix A.
  //C          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
  //C          orthonormal eigenvectors of the matrix A.
  //C          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
  //C          or the upper triangle (if UPLO='U') of A, including the
  //C          diagonal, is destroyed.
  //C
  //C  LDA     (input) INTEGER
  //C          The leading dimension of the array A.  LDA >= max(1,N).
  //C
  //C  W       (output) DOUBLE PRECISION array, dimension (N)
  //C          If INFO = 0, the eigenvalues in ascending order.
  //C
  //C  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
  //C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
  //C
  //C  LWORK   (input) INTEGER
  //C          The length of the array WORK.  LWORK >= max(1,3*N-1).
  //C          For optimal efficiency, LWORK >= (NB+2)*N,
  //C          where NB is the blocksize for DSYTRD returned by ILAENV.
  //C
  //C          If LWORK = -1, then a workspace query is assumed; the routine
  //C          only calculates the optimal size of the WORK array, returns
  //C          this value as the first entry of the WORK array, and no error
  //C          message related to LWORK is issued by XERBLA.
  //C
  //C  INFO    (output) INTEGER
  //C          = 0:  successful exit
  //C          < 0:  if INFO = -i, the i-th argument had an illegal value
  //C          > 0:  if INFO = i, the algorithm failed to converge; i
  //C                off-diagonal elements of an intermediate tridiagonal
  //C                form did not converge to zero.
  //C
  //C  =====================================================================
  //C
  //C     .. Parameters ..
  //C     ..
  //C     .. Local Scalars ..
  //C     ..
  //C     .. External Functions ..
  //C     ..
  //C     .. External Subroutines ..
  //C     ..
  //C     .. Intrinsic Functions ..
  //C     ..
  //C     .. Executable Statements ..
  //C
  //C     Test the input parameters.
  //C
  bool wantz = lsame(jobz, "V");
  bool lower = lsame(uplo, "L");
  bool lquery = (lwork ==  - 1);
  //C
  info = 0;
  if (!(wantz || lsame(jobz, "N"))) {
    info = -1;
  }
  else if (!(lower || lsame(uplo, "U"))) {
    info = -2;
  }
  else if (n < 0) {
    info = -3;
  }
  else if (lda < fem::max(1, n)) {
    info = -5;
  }
  //C
  int nb = fem::int0;
  int lwkopt = fem::int0;
  if (info == 0) {
    nb = ilaenv(1, "DSYTRD", uplo, n, -1, -1, -1);
    lwkopt = fem::max(1, (nb + 2) * n);
    work(1) = lwkopt;
    //C
    if (lwork < fem::max(1, 3 * n - 1) && !lquery) {
      info = -8;
    }
  }
  //C
  if (info != 0) {
    xerbla("DSYEV ", -info);
    return;
  }
  else if (lquery) {
    return;
  }
  //C
  //C     Quick return if possible
  //C
  if (n == 0) {
    return;
  }
  //C
  const double one = 1.0e0;
  if (n == 1) {
    w(1) = a(1, 1);
    work(1) = 2;
    if (wantz) {
      a(1, 1) = one;
    }
    return;
  }
  //C
  //C     Get machine constants.
  //C
  double safmin = dlamch(cmn, "Safe minimum");
  double eps = dlamch(cmn, "Precision");
  double smlnum = safmin / eps;
  double bignum = one / smlnum;
  double rmin = fem::sqrt(smlnum);
  double rmax = fem::sqrt(bignum);
  //C
  //C     Scale matrix to allowable range, if necessary.
  //C
  double anrm = dlansy("M", uplo, n, a, lda, work);
  int iscale = 0;
  const double zero = 0.0e0;
  double sigma = fem::double0;
  if (anrm > zero && anrm < rmin) {
    iscale = 1;
    sigma = rmin / anrm;
  }
  else if (anrm > rmax) {
    iscale = 1;
    sigma = rmax / anrm;
  }
  if (iscale == 1) {
    dlascl(cmn, uplo, 0, 0, one, sigma, n, n, a, lda, info);
  }
  //C
  //C     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
  //C
  int inde = 1;
  int indtau = inde + n;
  int indwrk = indtau + n;
  int llwork = lwork - indwrk + 1;
  int iinfo = fem::int0;
  dsytrd(cmn, uplo, n, a, lda, w, work(inde), work(indtau), work(indwrk),
    llwork, iinfo);
  //C
  //C     For eigenvalues only, call DSTERF.  For eigenvectors, first call
  //C     DORGTR to generate the orthogonal matrix, then call DSTEQR.
  //C
  if (!wantz) {
    dsterf(cmn, n, w, work(inde), info);
  }
  else {
    dorgtr(uplo, n, a, lda, work(indtau), work(indwrk), llwork, iinfo);
    dsteqr(cmn, jobz, n, w, work(inde), a, lda, work(indtau), info);
  }
  //C
  //C     If matrix was scaled, then rescale eigenvalues appropriately.
  //C
  int imax = fem::int0;
  if (iscale == 1) {
    if (info == 0) {
      imax = n;
    }
    else {
      imax = info - 1;
    }
    dscal(imax, one / sigma, w, 1);
  }
  //C
  //C     Set WORK(1) to optimal workspace size.
  //C
  work(1) = lwkopt;
  //C
  //C     End of DSYEV
  //C
}

} // namespace lapack_fem

#endif // GUARD
