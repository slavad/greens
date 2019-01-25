program green
  implicit none
  real*8::R_in,R_out,R_,n,t
  integer::n_R, n_root
  R_in = 1.d0
  R_out = 10.d0
  R_ = 2.d0
  n = 0.75d0
  t = 1.1d0
  n_root = 2000
  n_R = 200
  call TestGreenFun01(R_in,R_out,R_,n,t,n_root,n_R)
end program green