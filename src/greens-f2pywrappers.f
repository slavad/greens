C     -*- fortran -*-
C     This file is autogenerated with f2py (version:2)
C     It contains Fortran 77 wrappers to fortran functions.

      subroutine f2pywrapgreenlip_f2 (greenlip_f2f2pywrap, r, r_, 
     &fv_in, r_in, r_out, n, k_i, num)
      external greenlip_f2
      real*8 r
      real*8 r_
      real*8 fv_in
      real*8 r_in
      real*8 r_out
      real*8 n
      integer num
      real*8 k_i(num)
      complex*16 greenlip_f2f2pywrap, greenlip_f2
      greenlip_f2f2pywrap = greenlip_f2(r, r_, fv_in, r_in, r_out,
     & n, k_i, num)
      end


      subroutine f2pywrapgreenlip_f (greenlip_ff2pywrap, r, r_, fv
     &_in, r_in, r_out, n, k_i, num)
      external greenlip_f
      real*8 r
      real*8 r_
      real*8 fv_in
      real*8 r_in
      real*8 r_out
      real*8 n
      integer num
      real*8 k_i(num)
      complex*16 greenlip_ff2pywrap, greenlip_f
      greenlip_ff2pywrap = greenlip_f(r, r_, fv_in, r_in, r_out, n
     &, k_i, num)
      end


      subroutine f2pywrapgreenlipdicsmass (greenlipdicsmassf2pywra
     &p, r_, t, r_in, r_out, n, k_i, num)
      external greenlipdicsmass
      real*8 r_
      real*8 t
      real*8 r_in
      real*8 r_out
      real*8 n
      integer num
      real*8 k_i(num)
      real*8 greenlipdicsmassf2pywrap, greenlipdicsmass
      greenlipdicsmassf2pywrap = greenlipdicsmass(r_, t, r_in, r_o
     &ut, n, k_i, num)
      end

