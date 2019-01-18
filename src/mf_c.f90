!=======================================
! making complex conjugation
!=======================================
complex*16 function c_link(c)
implicit none
complex*16,intent(in)::c
complex*16::ii
  ii=(0.d00,1.d00)
  c_link=c-2*ii*aimag(c)
  !if(c.ne.0.d00)then
  !  write(*,*)"==>",c,c_link
  !end if
return
end function c_link


!============================================================
! Making a complex numer out of absolute value and argument.
!============================================================
complex*16 function c_make(r,fi)
implicit none
real*8,intent(in)::r,fi
complex*16::ii
  ii=(0.d0,1.d0)
  c_make=r*cos(fi)+ii*r*sin(fi)
return
end function c_make


!========================================================
! An absolute value of a complex number q.
!========================================================
real*8 function c_mod(q)
implicit none
complex*16,intent(in)::q
real*8::x,y
  x=q
  y=aimag(q)
  c_mod=sqrt(x*x+y*y)
return
end function c_mod


!================================================
! An argument of a complex number q.
!================================================
real*8 function arg(q)
implicit none
complex*16,intent(in)::q
real*8::q1,q2,l,dop
real*8::pi
  pi=3.141592653589793d0
  if(q.eq.(0.d0,0.d0))then
    arg=0.d0; goto 11
  end if
  q1=q
  q2=aimag(q)
  l=sqrt(q1**2+q2**2)
  arg=acos(q1/l)
  dop=acos(q1/l)
  if(q2.le.0) then
    arg=2.d0*pi-dop
  end if
11 return
end function arg



complex*16 function f_rc(x)  !создание комплексного числа из вещественного
implicit none
real*8::x
complex*16::ii
  ii=(0.d00,1.d00)
  f_rc=x+0.d00*ii
return
end function f_rc
