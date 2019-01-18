!====================================================================
! Ordinary Bessel functions.
!====================================================================
RECURSIVE real*8 function funBessel(alpha,x) result(res)
implicit none
real*8,intent(in)::alpha,x
integer::m,m_max
real*8::pi=3.141592653589793d0
real*8::a_0,a_m,djmuz

  if(alpha.gt.0.d0)then
    res=djmuz(alpha,x)  !==function from DIN for arbitrary alpha>0 and x>0==!
  else
    res=2.d0*(alpha+1.d0)/x*funBessel(alpha+1.d0,x)-funBessel(alpha+2.d0,x)
  end if
  goto 122
  !===============================================!

  !if(x.gt.1.d2*abs(alpha**2-0.25d0))then
  !if(x.gt.1.d2*abs(alpha**2-0.25d0))then
  !  funBessel=sqrt(2.d0/pi/x)*cos(x-alpha*pi/2.d0-pi/4.d0)
  !else
  !  if(abs(alpha-nint(alpha)).le.1.d-4)then
  !    funBessel=BESSEL_JN(nint(alpha),x)    !==Fortran internal function==!
  !  else
  !    m_max=100   !==has to be checked==!
  !    m=0; funBessel=0.d0; a_0=1.d0/gamma(alpha+1.d0)*(x/2.d0)**alpha; a_m=a_0
  !    do while(m.le.m_max)
  !      !funBessel=funBessel+((-1.d0)**m)*((x/2.d0)**(2.d0*m+alpha))/gamma(1.d0*m+1.d0)/gamma(1.d0*m+alpha+1.d0)
  !      funBessel=funBessel+a_m
  !      a_m=a_m*(-1.d0)/(m+1.d0)/(m+1.d0+alpha)*(x/2.d0)**2
  !      m=m+1
  !    end do
  !  end if
  !end if
122 return
end function funBessel


!=====================================================================
! Has to be revritten for integer alpha!
!=====================================================================
RECURSIVE real*8 function funNeumann(alpha,x) result(res)
implicit none
real*8,intent(in)::alpha,x
real*8::funBessel  !==functions==!
real*8::pi=3.141592653589793d0
real*8::delta,a_0,a_k,b_0,b_k
integer::k

  if(alpha.gt.0.d0)then
    if(abs(alpha-nint(alpha)).le.1.d-4)then
      res=BESSEL_YN(nint(alpha),x)    !==Fortran internal function==!
    else
      res=funNeumann_alpha_nonint(alpha-nint(alpha),nint(alpha),x)
    end if
  else
    res=-2.d0*(alpha+1.d0)/x*funNeumann(alpha+1.d0,x)-funNeumann(alpha+2.d0,x)
  end if
  goto 122
  !===========================================!

  !if(x.gt.1.d2*abs(alpha**2-0.25d0))then
  !  funNeumann=sqrt(2.d0/pi/x)*sin(x-alpha*pi/2.d0-pi/4.d0)
  !else
  !  delta=alpha/1-int(alpha/1)
  !  if((delta.le.1.d-6).or.(delta.ge.0.999999d0))then
  !    funNeumann=funNeumann_int(idnint(alpha),x)
  !  else
  !    k=0
  !    a_0=(x/2.d0)**alpha/gamma(1.d0+alpha)*cos(alpha*pi)/sin(alpha*pi)
  !    b_0=-(x/2.d0)**(-alpha)/gamma(1.d0-alpha)/sin(alpha*pi)
  !    funNeumann=a_0+b_0
  !    a_k=a_0; b_k=b_0
  !    do while(k.le.30)
  !      k=k+1; a_k=a_k*(-1.d0)*(x/2.d0)**2/k/(alpha+k); b_k=b_k*(-1.d0)*(x/2.d0)**2/k/(1.d0*k-alpha)
  !      funNeumann=funNeumann+a_k+b_k
  !    end do
  !  end if
  !end if
122 return
contains

  !==positive orders only==!
  real*8 function funNeumann_alpha_nonint(alpha_frac,alpha_int,x)
  implicit none
  real*8,intent(in)::alpha_frac,x
  integer,intent(in)::alpha_int
  real*8:: res
  dimension res(alpha_int+1)
  integer::ncalc
    call rybesl ( x, alpha_frac, alpha_int+1, res, ncalc )   !==we use the subroutine from numerical library "specfun"==!
    funNeumann_alpha_nonint=res(alpha_int+1)
  return
  end function funNeumann_alpha_nonint

  real*8 function funNeumann_int(n,x)
  implicit none
  integer,intent(in)::n
  real*8,intent(in)::x
  real*8::funBessel  !==function==!
  real*8::sum1,sum2,sum3,sum31,sum32
  real*8::C,pi  !==constants
  integer::k,k_max,m
    k_max=40
    C=0.5772156649015325d0  !==Euler constant
    pi=3.141592653589793d0
    funNeumann_int=2.d0/pi*funBessel(1.d0*n,x)*(log(x/2.d0)+C)
    !==============
    k=1; sum1=0.d0
    do while(k.le.n)
      sum1=sum1+1.d0/(k*1.d0); k=k+1
    end do
    sum1=sum1/pi*((x/2.d0)**n)/gamma(1.d0*n+1.d0)
    !==============
    k=1; sum2=0.d0
    do while(k.le.n)
      sum2=sum2+gamma(1.d0*n-1.d0*k+1.d0)/gamma(1.d0*k)*((x/2.d0)**(-1.d0*n+2.d0*k-2.d0))
      k=k+1
    end do
    sum2=sum2/pi
    !==============
    k=1; sum3=0.d0
    do while(k.le.k_max)
      sum31=0.d0; m=1
      do while(m.le.(n+k)); sum31=sum31+1.d0/(1.d0*m); m=m+1; end do
      sum32=0.d0; m=1
      do while(m.le.k); sum32=sum32+1.d0/(1.d0*m); m=m+1; end do
      sum3=sum3+((-1.d0)**k)*((x/2.d0)**(n+2.d0*k))/gamma(k+1.d0)/gamma(k+n+1.d0)*(sum31+sum32)
      k=k+1
    end do
    sum3=sum3/pi
    funNeumann_int=funNeumann_int-sum1-sum2-sum3
  return
  end function funNeumann_int

end function funNeumann
!=============================================================================




!====================================================================
! Modified Bessel functions of the first kind.
!====================================================================
real*8 function modifiedBessel1(alpha,x)
implicit none
real*8,intent(in)::alpha,x
integer::i,k

 if(x.le.(0.013d0*sqrt(alpha+1.d0)))then
   !==aproximation for small parameter==!
   modifiedBessel1=1.d0/gamma(alpha+1.d0)*(x/2.d0)**alpha
   !write(*,*)"use"
 else
    k=100
    i=0; modifiedBessel1=0.d0
    do while(i.le.k)
      modifiedBessel1=modifiedBessel1+((x**2/4.d0)**i)/dgamma(1.d0*i+1.d0)/dgamma(alpha+1.d0*i+1.d0)
      i=i+1
    end do
    modifiedBessel1=modifiedBessel1*((x/2.d0)**alpha)
  end if
return
end function modifiedBessel1

