!==========================================================================================================================
! The function calculates the Green function in the frequency domain. The calculations are based on semi-analytical manner.
! The radial dirivative in the definition of dotM is taken analytically.
! fv_in - frequency in units of viscous frequency at R=1.
!==========================================================================================================================
complex*16 function GreenLip_F2(R,R_,fv_in,R_in,R_out,n,k_i,num)
implicit none
real*8,intent(in)::R,R_,R_in,R_out,fv_in,n
integer::num
real*8::res_Re,res_Im,k,l,x_in,x_out,f
complex*16::c_make  !==function==!
complex*16::fGLip
real*8::pi=3.141592653589793d0
real*8:: k_i
dimension k_i(num)
integer::i

  !f=fv_in*(R_in/R_out)**(n-2.d0)   !==the frequency in units of viscous frequency at R_out, this one is used in further calculations==!
  f=fv_in*(1.d0/R_out)**(n-2.d0)

  l=0.5d0/(2.d0-n)
  x_in=(R_in/R_out)**(1.d0-n/2.d0)
  x_out=1.d0

  GreenLip_F2=(0.d0,0.d0)
  i=1
  do while(i.le.num)
    fGLip=fGreenLip_C(k_i(i),R,R_,R_in,R_out,f,n)

    GreenLip_F2=GreenLip_F2+fGLip
    i=i+1
  end do

  !==construct complex Green function out of the Re and Im parts==!
  GreenLip_F2=GreenLip_F2*6*pi*R**0.5d0*R_**0.25d0*(2.d0-n)
  GreenLip_F2=GreenLip_F2/4.d0
return
contains

  complex*16 function fGreenLip_C(k_i,R,R_,R_in,R_out,f,n)
  implicit none
  real*8::pi=3.141592653589793d0
  real*8,intent(in)::k_i,R,R_,R_in,R_out,f,n
  real*8::l,x,x_,x_in,x_out
  complex*16::res,ii
  real*8::funBessel  !==functions==!
    ii=(0.d0,1.d0)
    l=0.5d0/(2.d0-n)
    x =(R/R_out)**(1.d0-0.5d0*n)
    x_=(R_/R_out)**(1.d0-0.5d0*n)
    x_in=(R_in/R_out)**(1.d0-0.5d0*n)
    x_out=1.d0 !(R_out/R_out)**(1.d0-0.5d0*n)
    res=1.d0/(4.d0*ii*pi*f+k_i**2*(n-2.d0)**2)

    fGreenLip_C=res*Vi(k_i,x_,x_in,l)/(Vi(k_i,x_out,x_in,l))**2 &
                *R**(-0.75d0)*0.25d0*(k_i*(2.-n)*x* &
                (funBessel(-l,k_i*x_in)*(funBessel(l-1.,k_i*x)-funBessel(l+1.,k_i*x))&
                -funBessel(l,k_i*x_in)*(funBessel(-l-1.,k_i*x)-funBessel(-l+1,k_i*x))) &
                +Vi(k_i,x,x_in,l))
  return
  end function fGreenLip_C


  real*8 function fGreenLip_Re(k_i,R,R_,R_in,R_out,f,n)
  implicit none
  real*8::pi=3.141592653589793d0
  real*8,intent(in)::k_i,R,R_,R_in,R_out,f,n
  real*8::l,x,x_,x_in,x_out
  real*8::res_Re
  complex*16::res,ii
  real*8::funBessel  !==functions==!
    ii=(0.d0,1.d0)
    l=0.5d0/(2.d0-n)
    x =(R/R_out)**(1.d0-0.5d0*n)
    x_=(R_/R_out)**(1.d0-0.5d0*n)
    x_in=(R_in/R_out)**(1.d0-0.5d0*n)
    x_out=(R_out/R_out)**(1.d0-0.5d0*n)
    res=1.d0/(4.d0*ii*pi*f+k_i**2*(n-2.d0)**2)
    res_Re=realpart(res)

    fGreenLip_Re=res_Re*Vi(k_i,x_,x_in,l)/(Vi(k_i,x_out,x_in,l))**2 &
                 *R**(-0.75d0)*0.25d0*(k_i*(2.-n)*x* &
                (funBessel(-l,k_i*x_in)*(funBessel(l-1.,k_i*x)-funBessel(l+1.,k_i*x))&
                 -funBessel(l,k_i*x_in)*(funBessel(-l-1.,k_i*x)-funBessel(-l+1,k_i*x))) &
                +Vi(k_i,x,x_in,l))

  return
  end function fGreenLip_Re
  !===========================================================================================


  real*8 function fGreenLip_Im(k_i,R,R_,R_in,R_out,f,n)
  implicit none
  real*8::pi=3.141592653589793d0
  real*8,intent(in)::k_i,R,R_,R_in,R_out,f,n
  real*8::l,x,x_,x_in,x_out
  real*8::res_Im
  complex*16::res,ii
  real*8::funBessel  !==functions==!
    ii=(0.d0,1.d0)
    l=0.5d0/(2.d0-n)
    x =(R/R_out)**(1.d0-0.5d0*n)
    x_=(R_/R_out)**(1.d0-0.5d0*n)
    x_in=(R_in/R_out)**(1.d0-0.5d0*n)
    x_out=(R_out/R_out)**(1.d0-0.5d0*n)
    res=1.d0/(4.d0*ii*pi*f+k_i**2*(n-2.d0)**2)
    res_Im=imagpart(res)

    fGreenLip_Im=res_Im*Vi(k_i,x_,x_in,l)/(Vi(k_i,x_out,x_in,l))**2 &
                  *R**(-0.75d0)*0.25d0*(k_i*(2.-n)*x* &
                  (funBessel(-l,k_i*x_in)*(funBessel(l-1.,k_i*x)-funBessel(l+1.,k_i*x))&
                  -funBessel(l,k_i*x_in)*(funBessel(-l-1.,k_i*x)-funBessel(-l+1,k_i*x))) &
                  +Vi(k_i,x,x_in,l))
  return
  end function fGreenLip_Im
  !========================================================================================

  !==function for GreenLip==!
  real*8 function Vi(k_i,x,x_in,l)
  implicit none
  real*8,intent(in)::k_i,x,x_in,l
  real*8::funBessel  !==functions==!
    Vi=funBessel(l,k_i*x)*funBessel(-l,k_i*x_in)-funBessel(-l,k_i*x)*funBessel(l,k_i*x_in)
  return
  end function Vi
end function GreenLip_F2
!================================================================================================



!==========================================================================================================================
! The function calculates the Green function in the frequency domain. The calculations are based on semi-analytical manner.
! fv_in - frequency in units of viscous frequency at the inner disc radius.
!==========================================================================================================================
complex*16 function GreenLip_F(R,R_,fv_in,R_in,R_out,n,k_i,num)
implicit none
real*8,intent(in)::R,R_,R_in,R_out,fv_in,n
integer::num
real*8::res_Re,res_Im,res_Re_,res_Im_,k,dR,l,x_in,x_out,f
complex*16::c_make
real*8::pi=3.141592653589793d0
real*8:: k_i
dimension k_i(num)
integer::i

  f=fv_in*(R_in/R_out)**(n-2.d0)   !==the frequency in units of viscous frequency at R_out, this one is used in further calculations==!

  dR=1.d-5   !0.05d0  !==this value provides the most accurate resolution of the numerical scheme

  l=0.5d0/(2.d0-n)
  x_in=(R_in/R_out)**(1.d0-n/2.d0)
  x_out=(R_out/R_out)**(1.d0-n/2.d0)
  !call rootsLipunova(k_i,num,l,x_in,x_out)

  res_Re=0.d0
  res_Im=0.d0
  res_Re_=0.d0
  res_Im_=0.d0
  i=1
  do while(i.le.num)
    res_Re=res_Re+fGreenLip_Re(k_i(i),R,R_,R_in,R_out,f,n)
    res_Im=res_Im+fGreenLip_Im(k_i(i),R,R_,R_in,R_out,f,n)

    res_Re_=res_Re_+fGreenLip_Re(k_i(i),R+dR,R_,R_in,R_out,f,n)
    res_Im_=res_Im_+fGreenLip_Im(k_i(i),R+dR,R_,R_in,R_out,f,n)

    i=i+1
  end do
  res_Re=res_Re *R**(-n-0.25d0)*R_**(5.d0/4.d0)*R_out**(n-2.d0) *(2.d0-n)
  res_Im=res_Im *R**(-n-0.25d0)*R_**(5.d0/4.d0)*R_out**(n-2.d0) *(2.d0-n)

  res_Re_=res_Re_ *(R+dR)**(-n-0.25d0)*R_**(5./4.)*R_out**(n-2.d0) *(2.d0-n)
  res_Im_=res_Im_ *(R+dR)**(-n-0.25d0)*R_**(5./4.)*R_out**(n-2.d0) *(2.d0-n)

  !==differentiation==!
  res_Re=(res_Re_*(R+dR)**(n+0.5d0)-res_Re*(R)**(n+0.5d0))/dR
  res_Im=(res_Im_*(R+dR)**(n+0.5d0)-res_Im*(R)**(n+0.5d0))/dR
  !write(*,*)"# ",res_Re,res_Im

  !==construct complex Green function out of the Re and Im parts==!
  GreenLip_F=(1.d0,0.d0)*res_Re+(0.d0,1.d0)*res_Im
  GreenLip_F=GreenLip_F*sqrt(R)*6.d0*pi *(R_out/R_in)**(2.d0-n)
  GreenLip_F=GreenLip_F/4.d0/R_    !=sheck the reason for the 4==!
return
contains

  real*8 function fGreenLip_Re(k_i,R,R_,R_in,R_out,f,n)
  implicit none
  real*8::pi=3.141592653589793d0
  real*8,intent(in)::k_i,R,R_,R_in,R_out,f,n
  real*8::l,x,x_,x_in,x_out,y,y_,y_in
  real*8::res_mod,res_arg,res_Re,res_Im
  complex*16::res,ii
    ii=(0.d0,1.d0)
    l=0.5d0/(2.d0-n)
    x =(R/R_out)**(1.d0-0.5d0*n)
    x_=(R_/R_out)**(1.d0-0.5d0*n)
    x_in=(R_in/R_out)**(1.d0-0.5d0*n)
    x_out=(R_out/R_out)**(1.d0-0.5d0*n)
    res=1.d0/(4*ii*pi*f+k_i**2*(n-2.d0)**2)
    res_Re=realpart(res)
    fGreenLip_Re=res_Re*Vi(k_i,x_,x_in,l)*Vi(k_i,x,x_in,l)/(Vi(k_i,x_out,x_in,l))**2
    !fGreenLip_Re=res_Re !Vi(k_i,x_,x_in,l)*Vi(k_i,x,x_in,l)/(Vi(k_i,x_out,x_in,l))**2

    !fGreenLip_Re=res_Re*Vi(k_i,x_,x_in,l)/(Vi(k_i,x_out,x_in,l))**2 * ((n+1.d0)*Vi(k_i,x,x_in,l)+ k_i*x/4.d0/l)
  return
  end function fGreenLip_Re
  !===========================================================================================


  real*8 function fGreenLip_Im(k_i,R,R_,R_in,R_out,f,n)
  implicit none
  real*8::pi=3.141592653589793d0
  real*8,intent(in)::k_i,R,R_,R_in,R_out,f,n
  real*8::l,x,x_,x_in,x_out,y,y_,y_in
  real*8::res_mod,res_arg,res_Re,res_Im,res_2
  complex*16::res,ii

    ii=(0.d0,1.d0)
    l=0.5d0/(2.d0-n)
    x =(R/R_out)**(1.d0-0.5d0*n)
    x_=(R_/R_out)**(1.d0-0.5d0*n)
    x_in=(R_in/R_out)**(1.d0-0.5d0*n)
    x_out=(R_out/R_out)**(1.d0-0.5d0*n)
    res=1.d0/(4.d0*ii*pi*f+k_i**2*(n-2.d0)**2)
    res_Im=imagpart(res)
    fGreenLip_Im=res_Im*Vi(k_i,x_,x_in,l)*Vi(k_i,x,x_in,l)/(Vi(k_i,x_out,x_in,l))**2
    !fGreenLip_Im=res_Im!Vi(k_i,x_,x_in,l)*Vi(k_i,x,x_in,l)/(Vi(k_i,x_out,x_in,l))**2

    !fGreenLip_Im=res_Im*Vi(k_i,x_,x_in,l)/(Vi(k_i,x_out,x_in,l))**2 * ((n+1.d0)*Vi(k_i,x,x_in,l)+ k_i*x/4.d0/l)
  return
  end function fGreenLip_Im
  !========================================================================================

  !==function for GreenLip==!
  real*8 function Vi(k_i,x,x_in,l)
  implicit none
  real*8,intent(in)::k_i,x,x_in,l
  real*8::funBessel  !==functions==!
    Vi=funBessel(l,k_i*x)*funBessel(-l,k_i*x_in)-funBessel(-l,k_i*x)*funBessel(l,k_i*x_in)
  return
  end function Vi
end function GreenLip_F
!================================================================================================




!================================================================================================
! The subroutine controlls the mass in accretion disc described by Lipunovas Green functions.
!================================================================================================
real*8 function GreenLipDicsMass(R_,t,R_in,R_out,n,k_i,num)
implicit none
real*8,intent(in)::R_,t,R_in,R_out,n,k_i
integer,intent(in)::num
dimension k_i(num)
real*8::R,dR,res,Sigma,Mdot
real*8::pi=3.141592653589793d0
integer::nn
  nn=2000
  dR=(R_out-R_in)/nn
  R=R_in
  res=0.d0
  do while(R.le.R_out)
    call GreenLip(R,R_,t,R_in,R_out,n,k_i,num,Sigma,Mdot)
    !Sigma=Sigma/R**(n+0.5d0)
    !write(*,*)R,Sigma,res
    !read(*,*)
    res=res+2.d0*pi*R*dR*Sigma
    R=R+dR
  end do
  GreenLipDicsMass=res
return
end function GreenLipDicsMass





!================================================================================================
! The subroutine calculates the surface density and mass accretion rate according to solutions given
! by Lipunova.
! R - radial coordiante of interest; R_ - the radial coordinate of the initial perturbation
! R_in - inner radius, R_out - outher radius
! n - power describing kinematic viscosity
! k_i(num) - array with the roots of Lipunova's equation
! num - the number of roots which we use. 
! tv_in - time-scale in units of viscous time at R=1
!!!!!! still numerical differentiation!!!!!
!================================================================================================
subroutine GreenLip(R,R_,tv_in,R_in,R_out,n,k_i,num,Sigma,Mdot)
implicit none
real*8,intent(in)::R,R_,R_in,R_out,tv_in,n
integer::num
real*8::t,res,res_,k,dR,l,x_in,x_out,resMdot,Sigma,Mdot
real*8::pi=3.141592653589793d0
real*8:: k_i
dimension k_i(num)
integer::i

  t=tv_in*R_out**(n-2.d0)

  dR=1.d-4    !==this value provides the most accurate resolution of the numerical scheme

  l=0.5d0/(2.d0-n)
  x_in=(R_in/R_out)**(1.d0-n/2.d0)
  x_out=(R_out/R_out)**(1.d0-n/2.d0)

  res=0.d0
  res_=0.d0
  i=1
  do while(i.le.num)
    res =res +fGreenLip(k_i(i),R,R_,R_in,R_out,t,n)
    res_=res_+fGreenLip(k_i(i),R+dR,R_,R_in,R_out,t,n)
    i=i+1
  end do
  !==getting the surface density==!
  res =res *R**(-n-0.25d0)*R_**(5.d0/4.d0)*R_out**(n-2.d0) *(2.d0-n)
  res_=res_*(R+dR)**(-n-0.25d0)*R_**(5.d0/4.d0)*R_out**(n-2.d0) *(2.d0-n)

  !==differentiation and getting dotM==!
  resMdot=(res_*(R+dR)**(n+0.5d0)-res*(R)**(n+0.5d0))/dR
  resMdot=resMdot*sqrt(R)*6.d0*pi

  Sigma=res
  Mdot=resMdot
return
contains

  real*8 function fGreenLip(k_i,R,R_,R_in,R_out,t,n)
  implicit none
  real*8::pi=3.141592653589793d0
  real*8,intent(in)::k_i,R,R_,R_in,R_out,t,n
  real*8::l,x,x_,x_in,x_out,y,y_,y_in
  real*8::res
    l=0.5d0/(2.d0-n)
    x =(R/R_out)**(1.d0-0.5d0*n)
    x_=(R_/R_out)**(1.d0-0.5d0*n)
    x_in=(R_in/R_out)**(1.d0-0.5d0*n)
    x_out=(R_out/R_out)**(1.d0-0.5d0*n)
    res=exp(-2.d0*(1.d0-n/2.d0)**2*k_i**2*t)
    fGreenLip=res*Vi(k_i,x_,x_in,l)*Vi(k_i,x,x_in,l)/(Vi(k_i,x_out,x_in,l))**2
  return
  end function fGreenLip

  !==function for GreenLip==!
  real*8 function Vi(k_i,x,x_in,l)
  implicit none
  real*8,intent(in)::k_i,x,x_in,l
  real*8::funBessel  !==functions==!
    Vi=funBessel(l,k_i*x)*funBessel(-l,k_i*x_in)-funBessel(-l,k_i*x)*funBessel(l,k_i*x_in)
  return
  end function Vi
end subroutine GreenLip



!=========================================================================================================
! The subroutine calculates the first n roots of the equation...
!=========================================================================================================
subroutine rootsLipunova(root,n,l,x_in,x_out)
implicit none
integer,intent(in)::n
real*8,intent(in)::l,x_in,x_out
real*8::root,root_down,root_up
dimension root(n),root_down(n,2),root_up(n,2)
real*8::k,dk,res1,res2,eps,k1,k2,res
integer::i

  dk=1.d-1      !==step: should be shousen somehow==!
  k=0.d0
  !root_down(1,1)=0.d0;  root_down(1,2)=0.d0
  res2=eq_Lip(k,l,x_in,x_out)
  i=1
  do while(i.le.n)
    res1=res2
    res2=eq_Lip(k+dk,l,x_in,x_out)
    !write(*,*)k,res1
    if((res1*res2).le.0.d0)then
      root_down(i,1)=k;  root_down(i,2)=res1
      root_up(i,1)=k+dk; root_up(i,2)=res2
      !write(*,*)i,root_down(i,2),root_up(i,2)
      i=i+1
    end if
    !write(*,*)k,eq_Lip(k,l,x_in,x_out)
    k=k+dk
  end do
  eps=1.e-5 !min(root_down(2,1)-root_down(1,1),root_down(n,1)-root_down(n-1,1))/1.d4
  !write(*,*)eps

  i=1
  do while(i.le.n)
    k1=root_down(i,1); k2=root_up(i,1)
    res1=root_down(i,2); res2=root_up(i,2)
    do while(abs(k2-k1).gt.eps)
      k=(k2+k1)/2.d0
      res=eq_Lip(k,l,x_in,x_out)
      if((res*res1).gt.0.d0)then
        res1=res; k1=k
      else
        res2=res; k2=k
      end if
    end do
    root(i)=(k1+k2)/2.d0
    !write(*,*)i,root(i),root(i)-root(max(i-1,1))
    !read(*,*)
    i=i+1
  end do

return
contains
  real*8 function eq_Lip(k,l,x_in,x_out)
  implicit none
  real*8,intent(in)::k,l,x_in,x_out
  real*8::funBessel,mo,mi
    mo=k*x_out
    mi=k*x_in
    !eq_Lip=l/x_out*((1.d0-k)*funBessel(l,k*x_out)*funBessel(-l,k*x_in)-(1.d0+k)*funBessel(-l,k*x_out)*funBessel(l,k*x_in))&
    !      +k**2*(funBessel(l-1.d0,k*x_out)*funBessel(-l,k*x_in)-funBessel(-l-1.d0,k*x_out)*funBessel(l,k*x_in))
    !eq_Lip=-2.d0*l/x_out*funBessel(-l,k*x_out)*funBessel(l,k*x_in)&
    !       +k*(funBessel(l-1.d0,k*x_out)*funBessel(-l,k*x_in)-funBessel(-l-1.d0,k*x_out)*funBessel(l,k*x_in))
    eq_Lip=funBessel(-l,mi)*(l/x_out*funBessel(l,mo)+k/2.d0*(funBessel(l-1.d0,mo)-funBessel(l+1.d0,mo)))&
          +funBessel(+l,mi)*(k/2.d0*(funBessel(-l+1.d0,mo)-funBessel(-l-1.d0,mo))-l/x_out*funBessel(-l,mo))
  return
  end function eq_Lip
end subroutine rootsLipunova


