c the following 20 subprograms calculates Bessel function J_mu(z)
c for arbitrary values of argument z>=0 and index mu>=0
c program calculates the Bessel function J_mu(z)
c by series, asymptotic and Cherry method
c djmuz chooses the suitable subprogram for given z and mu (z and dmu)
c djser calculates the function with series, djasp --- with asymptotics,
c djchr --- by Cherry method
c subroutines ph01ln, ph01tg, dartgm, darthm are auxilary for djchr
c dcheb0 calculates Bessel functions with expansions on Chebyshev polinomials
c di13ml, di13pl, dj13ml, dj13pl, dk13gs and yj12gs calculate Bessel
c functions with indices 1/3 and -1/3, dgam calculates gamma-function
c dj0, dj0l, yjog and yj1g calculate the Bessel function J_0(z)
c value of j_mu(z) can be obtained with the command a=djmuz(mu,z)
 
      function djmuz(dmu,z)
      implicit real*8(a-h,o-z)
      data zlw,zup,dmul,dmup /1d1,1d2,3.5d0,4.5d0/
      if((z.eq.0d0).or.(dmu.eq.0d0)) goto 100
      if((z.le.zlw).and.(dmu.le.zlw)) then
          djmuz=djser(dmu,z)
        else
          if((z.ge.zup).and.(dmu.le.dmul)) then
             djmuz=djasp(dmu,z)
           else
             if(z.le.zlw) dml=1d1
             if((z.gt.zlw).and.(z.le.12.5d0)) dml=dmup
             if(z.gt.12.5d0) dml=z-8d0
             if(dmu.ge.dml) then
                djmuz=djchr(dmu,z)
               else
                imu=idint(dml-dmu)+1
                dmn=dmu+1d0*imu
                djp1=djchr(dmn+1d0,z)
                djo=djchr(dmn,z)
                do nd=1,imu
                   djm1=-djp1+2d0*(dmn-1d0*(nd-1))/z*djo
                   djp1=djo
                   djo=djm1
                enddo
                djmuz=djo
             endif
             return
          endif
      endif
      return
  100 continue
      if(dmu.gt.0d0) then
         djmuz=0d0
        else
         djmuz=dj0(z)
      endif
      return
      end
 
c subprogram calculates the Bessel function J_mu(z) through series,
c J_mu(z)=(1/gamma(mu+1))(z/2)^mu
c sum_{n=0}^infty prod_{j=1}^n(-1)(z/2)^2/j(mu+j), 0<mu<=10, 0<=z<=5
      function djser(dmu,z)
      implicit real*8(a-h,o-z)
      data eps /1d-9/
      if(z.eq.0d0) djser=0d0
      if(z.eq.0d0) return
      zs=0.5d0*z
      nmu=idint(dmu)
      ddmu=dmu-nmu
      if(ddmu.eq.0d0) then
        zsmu=1d0
       else
        zsmu=zs**ddmu
      endif
        dgmu=dgam(1d0+dmu)
      if(nmu.eq.0) goto 30
      do 20 nn=1,nmu
         zsmu=zs*zsmu
   20 continue
   30 continue
      dan=1d3
      dmz=1d0
      dn=1d0
      zz=-zs*zs
      an=1d0
      do while (dan.gt.eps)
         an=an*zz/dn/(dmu+dn)
         dmz=dmz+an
         dan=dabs(an)
         dn=dn+1d0
      enddo
      djser=zsmu/dgmu*dmz
      return
      end
 
c subprogram calculates the Bessel function J_mu(z) through asymptotics,
c J_mu(z)=sqrt(2/pi z)[cos(z-pi mu/2-pi/4)sum_{n=0}^infty a_{2n}+
c sin(z-pi mu/2-pi/4)sum_{n=0}^infty a_{2n+1}]
c a_n=(-1)^{[n/2]}prod_{k=1}^n[mu^2-(2k-1)^2/4]/2zk
c  calculation for z>=10, mu=<10
      function djasp(dmu,z)
      implicit real*8(a-h,o-z)
      data eps /1d-9/
      data pis,pif /1.570796326794d0,0.7853981633974483d0/
      zt=0.5d0/z
      dsqz=dsqrt(zt/pif)
      dmu2=dmu*dmu
      zmu=z-pis*dmu-pif
      czm=dcos(zmu)
      szm=dsin(zmu)
      sume=0d0
      sumo=0d0
      dn=1d0
      an=1d0
      dann=1d3
      dan=1d0
  300 continue
      do while ((dan.gt.eps).and.(dan.lt.dann))
         sume=sume+an
         dno=dn-0.5d0
         an=zt/dn*(dmu2-dno*dno)*an
         sumo=sumo+an
         dn=dn+1d0
         dne=dn-0.5d0
         an=-an*zt/dn*(dmu2-dne*dne)
         dann=dan
         dan=dabs(an)
      enddo
      djasp=dsqz*(czm*sume-szm*sumo)
 2001 format(1x,6d12.5)
      end
 
c subprogram calculates the Bessel function J_mu(z) by Cherry method:
c          /K_{1/3}(nu*xi), v_{+}=arth(zeta)/zeta, 0=<zeta<1,
c J_mu(z)=|     xi=nu*phi(v), phi(v)=v+phi_1(v)/mu^2+phi_2/mu^4,
c          \J_{1/3}(nu*xi)+J_{-1/3}(nu*xi), v_{-}=arctg(zeta)/zeta, z>0.
c z=lambda*sqrt(1-zeta^2) for upper row, z=lambda*sqrt(1+zeta^2) for low
c nu^2/mu^2=1+c_1/mu^2+c_2/mu^4, lambda^2/mu^2=1+gamma_1/mu^2+gamma_2/mu
      function djchr(dmu,z)
      implicit real*8(a-h,o-z)
      data c1,c2 /0.43809523809523809d-01,-0.481826817173756d-02/
      data g1,g2 /0.285714285714d-01,-0.217017268445840d-02/
      data pif,pit /0.785398163397448d0,0.628318530717959d1/
      data pis /1.570796326794897d0/
      data c13,sq13 /0.333333333333333d0,0.577350269189626d0/
      dm2=dmu*dmu  ! mu^2
      dmr2=1d0/dm2 ! 1/mu^2
      dl2=dm2*(1d0+dmr2*(g1+g2*dmr2)) ! [1+(g_1+g_2/mu^2)/mu^2]/mu^2
      dlb=dsqrt(dl2) ! sqroot(dl2)
      dnu=dmu*dsqrt(1d0+dmr2*(c1+c2*dmr2)) ! mu sqroot(1+(c_1+c_2/mu^2)/mu^2)
      z2l2=z*z/dl2 ! z^2/dl2
      if(z2l2.gt.1d0) goto 200
      zet=dsqrt(1d0-z2l2)
      zet1=z2l2/(1d0+zet)
      call ph01ln(zet,zet1,v,phi10n,phi11n,phi20n,phi21n)
      phi=v+dmr2*(phi10n+phi20n*dmr2)
      phi1=1d0+dmr2*(phi11n+phi21n*dmr2)
      xi=dnu*phi
      if(xi.lt.6d0) then
          xi13=xi**c13
          wstar=(di13ml(xi)/xi13-di13pl(xi)*xi13)*sq13
        else
          wstar=dk13gs(xi)/dsqrt(pit*xi)*dexp(-xi)
      endif
      djchr=dsqrt(phi/zet/phi1)*wstar
      return
  200 continue
      zet=dsqrt(z2l2-1d0)
      call ph01tg(zet,v,phi10g,phi11g,phi20g,phi21g)
      phi=v+dmr2*(phi10g+phi20g*dmr2)
      phi1=1d0+dmr2*(phi11g+phi21g*dmr2)
      xi=dnu*phi
      if(xi.lt.6d0) then
          xi13=xi**c13
          dj13m=dj13ml(xi)
          dj13p=dj13pl(xi)
          wstar=(dj13ml(xi)/xi13+dj13pl(xi)*xi13)*sq13
        else
          call yj13gs(xi,yjr,yji)
          x=xi-pif
          wstar=(dcos(x)*yjr-dsin(x)*yji)/dsqrt(pis*xi)
      endif
      djchr=dsqrt(phi/zet/phi1)*wstar
      return
      end
 
c programm calculates functions phi_1(v) and phi_2(v) and their derivatives
c for v_{+}, 0=<z<1, phi_j=y^{5/2}*p_j(y),phi_j'=y*p_{j1}(y),y=zeta^2
      subroutine ph01ln(z,z1,v,phi10n,phi11n,phi20n,phi21n)
      implicit real*8(a-h,o-z)
      common /warh/ arthm(30)
      data c13 /0.3333333333333333d0/
      data zlim /0.9d0/
 
      data aa1 /-0.438095238095238d-02/ ! -23/5250
      data aa2 /-0.612244897959184d-02/ ! -3/490
      data aa3 /-0.416666666666667d-01/ ! -1/24
      data aa4 / 0.568181818181818d-01/ ! 5/88
 
      data ab10 / 0.131428571428571d-01/  ! 23/1750
      data ab11 /-0.457142857142857d-02/  ! 4/875
      data ab20 / 0.183673469387755d-01/  ! 9/490
      data ab21 /-0.596938775510204d-01/  ! -117/1960
      data ab30 / 0.125000000000000d+00/  ! 1/8
      data ab31 / 0.250000000000000d+00/  ! 1/4
      data ab40 /-0.170454545454545d+00/  ! -15/88
 
      data ad10 /-0.789051184492001d-02/ ! -5224793159/662161500000
      data ad11 /-0.451342199750363d-03/ ! -24905119/55180125000
      data ad12 / 0.625772395405049d-03/ ! 3139109/5016375000
      data ad20 /-0.740895804875397d-02/ ! -21804119/2942940000
      data ad21 /-0.692858162245917d-04/ ! -2124/30655625
      data ad22 / 0.642879569410182d-03/ ! 14333/22295000
      data ad30 /-0.810192056620628d-02/ ! -471629/58212000
      data ad31 / 0.627705627705628d-03/ ! 29/46200
      data ad32 / 0.154195011337869d-02/ ! 17/11025
      data ad40 /-0.104130591630592d-01/ ! -5773/554400
      data ad41 / 0.194444444444444d-02/ ! 7/3600
      data ad42 /-0.219223484848485d-01/ ! -463/21120
      data ad50 /-0.172390109890110d-01/ ! -251/14560
      data ad51 / 0.493303571428571d-01/ ! 221/4480
      data ad52 / 0.193526785714286d+00/ ! 867/4480
      data ad60 /-0.130729166666667d+00/ ! -251/1920
      data ad61 /-0.230208333333333d+00/ ! -221/960
      data ad62 /-0.345312500000000d+00/ ! -221/640
      data ad70 / 0.192248774509804d+00/ ! 1255/6528
      data ad71 / 0.169270833333333d+00/ ! 65/384
 
      data ae10 / 0.710146066042801d-01/ ! 5224793159/73573500000
      data ae11 / 0.430123547201098d-02/ ! 79114237/18393375000
      data ae12 /-0.419736981385961d-02/ ! -2339509/557375000
      data ae13 / 0.977285789040891d-03/ ! 204268/209015625
      data ae20 / 0.666806224387857d-01/ ! 65412357/980980000
      data ae21 / 0.636783624538727d-03/ ! 19521/30655625
      data ae22 /-0.381988215865767d-02/ ! -936807/245245000
      data ae23 / 0.262538685803992d-02/ ! 58533/22295000
      data ae30 / 0.729172850958565d-01/ ! 471629/6468000
      data ae31 /-0.609461966604824d-02/ ! -657/107800
      data ae32 /-0.337043908472480d-02/ ! -109/32340
      data ae33 /-0.819238945578231d-01/ ! -38537/470400
      data ae40 / 0.937175324675325d-01/ ! 5773/61600
      data ae41 /-0.225324675324675d-01/ ! -347/15400
      data ae42 / 0.182366071428571d+00/ ! 817/4480
      data ae43 / 0.120933441558442d+01/ ! 14899/12320
      data ae50 / 0.155151098901099d+00/ ! 2259/14560
      data ae51 /-0.464062500000000d+00/ ! -297/640
      data ae52 /-0.174174107142857d+01/ ! -7803/4480
      data ae53 /-0.375669642857143d+01/ ! -1683/448
      data ae60 / 0.117656250000000d+01/ ! 753/640
      data ae61 / 0.216562500000000d+01/ ! 693/320
      data ae62 / 0.310781250000000d+01/ ! 1989/640
      data ae63 / 0.414375000000000d+01/ ! 663/160
      data ae70 /-0.173023897058824d+01/ ! -3765/2176
      data ae71 /-0.159237132352941d+01/ ! -3465/2176
      data ae72 /-0.152343750000000d+01/ ! -195/128
 
      data aav1 /-0.219047619047619d-01/ ! -23d0/1050d0
      data aav2 /-0.694444444444445d-01/ ! -5d0/72d0
      data aaz1 / 0.208333333333333d+00/ ! 5d0/24d0
      data aaz2 /-0.125000000000000d+00/ ! -1d0/8d0
      data aaz3 /-0.142857142857143d-01/ ! -1d0/70d0
 
      data abv1 /-0.761904761904762d-02/ ! -4d0/525d0
      data abv2 / 0.694444444444445d-01/ ! 5d0/72d0
      data abz1 /-0.625000000000000d+00/ ! -5d0/8d0
      data abz2 / 0.750000000000000d+00/ ! 3d0/4d0
      data abz3 /-0.139285714285714d+00/ ! -39d0/280d0
 
      data adv1 / 0.312886197702524d-02/ ! 3139109d0/1003275000d0
      data adv2 / 0.152116402116402d-02/ ! 23d0/15120d0
      data adv3 /-0.992063492063492d-03/ ! -1d0/1008d0
      data adv4 /-0.868055555555556d-02/ ! -5d0/576d0
      data adv5 / 0.144675925925926d-01/ ! 25d0/1728d0
      data adv6 /-0.403485082304527d-01/ ! -1255d0/31104d0
      data adz1 / 0.150005232862376d-02/ ! 14333d0/9555000d0
      data adz2 / 0.462585034013606d-02/ ! 17d0/3675d0
      data adz3 /-0.803819444444445d-01/ ! -463d0/5760d0
      data adz4 / 0.838616071428571d+00/ ! 3757d0/4480d0
      data adz5 /-0.172656250000000d+01/ ! -221d0/128d0
      data adz6 / 0.959201388888889d+00/ ! 1105d0/1152d0
 
      data aev1 / 0.162880964840149d-02/ ! 204268d0/125409375d0
      data aev2 /-0.529100529100529d-03/ ! -1d0/1890d0
      data aev3 /-0.967261904761905d-02/ ! -13d0/1344d0
      data aev4 / 0.520833333333333d-01/ ! 5d0/96d0
      data aev5 /-0.434027777777778d-01/ ! -25d0/576d0
      data aev6 / 0.198412698412698d-02/ ! 1d0/504d0
      data aev7 / 0.173611111111111d-01/ ! 5d0/288d0
      data aev8 /-0.289351851851852d-01/ ! -25d0/864d0
      data aev9 / 0.121045524691358d+00/ ! 1255d0/10368d0
      data aez1 / 0.612590266875981d-02/ ! 19511d0/3185000d0
      data aez2 /-0.245771683673469d+00/ ! -38537d0/156800d0
      data aez3 / 0.443422619047619d+01/ ! 14899d0/3360d0
      data aez4 /-0.162790178571429d+02/ ! -7293d0/448d0
      data aez5 / 0.207187500000000d+02/ ! 663d0/32d0
      data aez6 /-0.863281250000000d+01/ ! -1105d0/128d0
 
      y=z*z
      if(z.gt.zlim) goto 100
      call darthm(z,z1,9)
      v0=arthm(2)
      v1=arthm(3)
      v2=arthm(4)
      v3=arthm(5)
      v4=arthm(6)
      v5=arthm(7)
      v6=arthm(8)
      v7=arthm(9)
      v=c13*y*z*v0
      v01=v0+1d0
      v02=v0*v0
      v03=v0*v02
      v04=v0*v03
      p10n=aa1*v01*v1+aa2*v2+aa3*v3+aa4*v4
      p11n=(ab10+ab11*v0)*v1+(ab20+ab21*v0)*v2+(ab30+ab31*v0)*v3+
     ,ab40*v01*v4
      p20n=(ad10+ad11*v0+ad12*v02*v01)*v1+
     ,(ad20+ad21*v0+ad22*v02)*v2+
     ,(ad30+ad31*v0+ad32*v02)*v3+
     ,(ad40+ad41*v0+ad42*v02)*v4+
     ,(ad50+ad51*v0+ad52*v02)*v5+
     ,(ad60+ad61*v0+ad62*v02)*v6+
     ,(ad70+ad71*v0*v01)*v7
      p21n=(ae10+ae11*v0+ae12*v02+ae13*v03)*v1+
     ,(ae20+ae21*v0+ae22*v02+ae23*v03)*v2+
     ,(ae30+ae31*v0+ae32*v02+ae33*v03)*v3+
     ,(ae40+ae41*v0+ae42*v02+ae43*v03)*v4+
     ,(ae50+ae51*v0+ae52*v02+ae53*v03)*v5+
     ,(ae60+ae61*v0+ae62*v02+ae63*v03)*v6+
     ,(ae70+ae71*v0+ae72*v02*v01)*v7
      y52=y*y*z
      phi10n=y52*p10n/v0
      phi11n=y*p11n/v02
      phi20n=y52*p20n/v03
      phi21n=y*p21n/v04
      return
  100 continue
      zr=1d0/z
      yr=zr*zr
      yr2=yr*yr
      yr3=yr*yr2
      zyr=zr*yr
      call darthm(z,z1,1)
      v=c13*z*y*arthm(2)
      vr=1d0/v
      vr2=vr*vr
      vr3=vr*vr2
      vr4=vr*vr3
      phi10n=aav1*v+aav2*vr+(aaz1*yr+aaz2+aaz3*y)*zr
      phi11n=abv1+abv2*vr2+((abz1*yr+abz2)*yr+abz3)*yr
      phi20n=adv1*v+adv2*vr+adv6*vr3+z*((adv3+yr*(adv4+yr*adv5))*vr2+
     ,(adz1+yr*(adz2+yr*(adz3+yr*(adz4+yr*(adz5+yr*adz6))))))
      phi21n=aev1+(aev2+yr*(aev3+yr*(aev4+yr*aev5)))*vr2+
     ,z*(aev6+yr*(aev7+yr*aev8))*vr3+aev9*vr4+
     ,yr*(aez1+yr*(aez2+yr*(aez3+yr*(aez4+yr*(aez5+yr*aez6)))))
      return
      end
 
c programm calculates functions phi_1(v) and phi_2(v) and their derivati
c for v_{-}, zeta>=0, phi_j=y^{5/2}*p_j(y), phi_j'=y*p_{j1}(y), y=zeta^2
      subroutine ph01tg(z,v,phi10g,phi11g,phi20g,phi21g)
      implicit real*8(a-h,o-z)
      common /warg/ artgm(30)
      data c13 /0.3333333333333333d0/
      data zlim /0.9d0/
 
      data aa1 / 0.438095238095238d-02/ ! 23/5250
      data aa2 / 0.612244897959184d-02/ ! 3/490
      data aa3 / 0.416666666666667d-01/ ! 1/24
      data aa4 /-0.568181818181818d-01/ ! -5/88
 
      data ab10 /-0.131428571428571d-01/  ! -23/1750
      data ab11 /0.457142857142857d-02/  ! 4/875
      data ab20 /-0.183673469387755d-01/  ! -9/490
      data ab21 / 0.596938775510204d-01/  ! 117/1960
      data ab30 /-0.125000000000000d+00/  ! -1/8
      data ab31 /-0.250000000000000d+00/  ! -1/4
      data ab40 / 0.170454545454545d+00/  ! 15/88
 
      data ad10 / 0.789051184492001d-02/ ! 5224793159/662161500000
      data ad11 / 0.451342199750363d-03/ ! 24905119/55180125000
      data ad12 /-0.625772395405049d-03/ ! -3139109/5016375000
      data ad20 /0.740895804875397d-02/ ! 21804119/2942940000
      data ad21 / 0.692858162245917d-04/ ! 2124/30655625
      data ad22 /-0.642879569410182d-03/ ! -14333/22295000
      data ad30 / 0.810192056620628d-02/ ! 471629/58212000
      data ad31 /-0.627705627705628d-03/ ! -29/46200
      data ad32 /-0.154195011337869d-02/ ! -17/11025
      data ad40 / 0.104130591630592d-01/ ! 5773/554400
      data ad41 /-0.194444444444444d-02/ ! -7/3600
      data ad42 / 0.219223484848485d-01/ ! 463/21120
      data ad50 / 0.172390109890110d-01/ ! 251/14560
      data ad51 /-0.493303571428571d-01/ ! -221/4480
      data ad52 /-0.193526785714286d+00/ ! -867/4480
      data ad60 / 0.130729166666667d+00/ ! 251/1920
      data ad61 / 0.230208333333333d+00/ ! 221/960
      data ad62 / 0.345312500000000d+00/ ! 221/640
      data ad70 /-0.192248774509804d+00/ ! -1255/6528
      data ad71 /-0.169270833333333d+00/ ! -65/384
 
      data ae10 /-0.710146066042801d-01/ ! -5224793159/73573500000
      data ae11 /-0.430123547201098d-02/ ! -79114237/18393375000
      data ae12 / 0.419736981385961d-02/ ! 2339509/557375000
      data ae13 /-0.977285789040891d-03/ ! -204268/209015625
      data ae20 /-0.666806224387857d-01/ ! -65412357/980980000
      data ae21 /-0.636783624538727d-03/ ! -19521/30655625
      data ae22 / 0.381988215865767d-02/ ! 936807/245245000
      data ae23 /-0.262538685803992d-02/ ! -58533/22295000
      data ae30 /-0.729172850958565d-01/ ! -471629/6468000
      data ae31 / 0.609461966604824d-02/ ! 657/107800
      data ae32 / 0.337043908472480d-02/ ! 109/32340
      data ae33 / 0.819238945578231d-01/ ! 38537/470400
      data ae40 /-0.937175324675325d-01/ ! -5773/61600
      data ae41 / 0.225324675324675d-01/ ! 347/15400
      data ae42 /-0.182366071428571d+00/ ! -817/4480
      data ae43 /-0.120933441558442d+01/ ! -14899/12320
      data ae50 /-0.155151098901099d+00/ ! -2259/14560
      data ae51 / 0.464062500000000d+00/ ! 297/640
      data ae52 / 0.174174107142857d+01/ ! 7803/4480
      data ae53 / 0.375669642857143d+01/ ! 1683/448
      data ae60 /-0.117656250000000d+01/ ! -753/640
      data ae61 /-0.216562500000000d+01/ ! -693/320
      data ae62 /-0.310781250000000d+01/ ! -1989/640
      data ae63 /-0.414375000000000d+01/ ! -663/160
      data ae70 / 0.173023897058824d+01/ ! 3765/2176
      data ae71 / 0.159237132352941d+01/ ! 3465/2176
      data ae72 / 0.152343750000000d+01/ ! 195/128
 
      data aav1 /-0.219047619047619d-01/ ! -23d0/1050d0
      data aav2 / 0.694444444444445d-01/ ! 5d0/72d0
      data aaz1 /-0.208333333333333d+00/ ! -5d0/24d0
      data aaz2 /-0.125000000000000d+00/ ! -1d0/8d0
      data aaz3 / 0.142857142857143d-01/ ! 1d0/70d0
 
      data abv1 /-0.761904761904762d-02/ ! -4d0/525d0
      data abv2 /-0.694444444444445d-01/ ! -5d0/72d0
      data abz1 / 0.625000000000000d+00/ ! 5d0/8d0
      data abz2 / 0.750000000000000d+00/ ! 3d0/4d0
      data abz3 / 0.139285714285714d+00/ ! 39d0/280d0
 
      data adv1 / 0.312886197702524d-02/ ! 3139109d0/1003275000d0
      data adv2 /-0.152116402116402d-02/ ! -23d0/15120d0
      data adv3 /-0.992063492063492d-03/ ! -1d0/1008d0
      data adv4 / 0.868055555555556d-02/ ! 5d0/576d0
      data adv5 / 0.144675925925926d-01/ ! 25d0/1728d0
      data adv6 /-0.403485082304527d-01/ ! -1255d0/31104d0
      data adz1 /-0.150005232862376d-02/ ! -14333d0/9555000d0
      data adz2 / 0.462585034013606d-02/ ! 17d0/3675d0
      data adz3 / 0.803819444444445d-01/ ! 463d0/5760d0
      data adz4 / 0.838616071428571d+00/ ! 3757d0/4480d0
      data adz5 / 0.172656250000000d+01/ ! 221d0/128d0
      data adz6 / 0.959201388888889d+00/ ! 1105d0/1152d0
 
      data aev1 / 0.162880964840149d-02/ ! 204268d0/125409375d0
      data aev2 / 0.529100529100529d-03/ ! 1d0/1890d0
      data aev3 /-0.967261904761905d-02/ ! -13d0/1344d0
      data aev4 /-0.520833333333333d-01/ ! -5d0/96d0
      data aev5 /-0.434027777777778d-01/ ! -25d0/576d0
      data aev6 / 0.198412698412698d-02/ ! 1d0/504d0
      data aev7 /-0.173611111111111d-01/ ! -5d0/288d0
      data aev8 /-0.289351851851852d-01/ ! -25d0/864d0
      data aev9 / 0.121045524691358d+00/ ! 1255d0/10368d0
      data aez1 /-0.612590266875981d-02/ ! -19511d0/3185000d0
      data aez2 /-0.245771683673469d+00/ ! -38537d0/156800d0
      data aez3 /-0.443422619047619d+01/ ! -14899d0/3360d0
      data aez4 /-0.162790178571429d+02/ ! -7293d0/448d0
      data aez5 /-0.207187500000000d+02/ ! -663d0/32d0
      data aez6 /-0.863281250000000d+01/ ! -1105d0/128d0
 
      y=z*z
      if(z.gt.zlim) goto 100
      call dartgm(z,9)
      v0=artgm(2)
      v1=artgm(3)
      v2=artgm(4)
      v3=artgm(5)
      v4=artgm(6)
      v5=artgm(7)
      v6=artgm(8)
      v7=artgm(9)
      v=c13*y*z*v0
      v01=v0+1d0
      v02=v0*v0
      v03=v0*v02
      v04=v0*v03
      p10g=aa1*v01*v1+aa2*v2+aa3*v3+aa4*v4
      p11g=(ab10+ab11*v0)*v1+(ab20+ab21*v0)*v2+(ab30+ab31*v0)*v3+
     ,ab40*v01*v4
      p20g=(ad10+ad11*v0+ad12*v02*v01)*v1+
     ,(ad20+ad21*v0+ad22*v02)*v2+
     ,(ad30+ad31*v0+ad32*v02)*v3+
     ,(ad40+ad41*v0+ad42*v02)*v4+
     ,(ad50+ad51*v0+ad52*v02)*v5+
     ,(ad60+ad61*v0+ad62*v02)*v6+
     ,(ad70+ad71*v0*v01)*v7
      p21g=(ae10+ae11*v0+ae12*v02+ae13*v03)*v1+
     ,(ae20+ae21*v0+ae22*v02+ae23*v03)*v2+
     ,(ae30+ae31*v0+ae32*v02+ae33*v03)*v3+
     ,(ae40+ae41*v0+ae42*v02+ae43*v03)*v4+
     ,(ae50+ae51*v0+ae52*v02+ae53*v03)*v5+
     ,(ae60+ae61*v0+ae62*v02+ae63*v03)*v6+
     ,(ae70+ae71*v0+ae72*v02*v01)*v7
      y52=y*y*z
      phi10g=y52*p10g/v0
      phi11g=y*p11g/v02
      phi20g=y52*p20g/v03
      phi21g=y*p21g/v04
      return
  100 continue
      zr=1d0/z
      yr=zr*zr
      yr2=yr*yr
      yr3=yr*yr2
      zyr=zr*yr
      call dartgm(z,1)
      v=c13*z*y*artgm(2)
      vr=1d0/v
      vr2=vr*vr
      vr3=vr*vr2
      vr4=vr*vr3
      vr5=vr*vr4
      phi10g=aav1*v+aav2*vr+(aaz1*yr+aaz2+aaz3*y)*zr
      phi11g=abv1+abv2*vr2+((abz1*yr+abz2)*yr+abz3)*yr
      phi20g=adv1*v+adv2*vr+adv6*vr3+z*((adv3+yr*(adv4+yr*adv5))*vr2+
     ,(adz1+yr*(adz2+yr*(adz3+yr*(adz4+yr*(adz5+yr*adz6))))))
      phi21g=aev1+(aev2+yr*(aev3+yr*(aev4+yr*aev5)))*vr2+
     ,z*(aev6+yr*(aev7+yr*aev8))*vr3+aev9*vr4+
     ,yr*(aez1+yr*(aez2+yr*(aez3+yr*(aez4+yr*(aez5+yr*aez6)))))
      return
      end
 
c programm calculates residuals of artg(x)/x series, v_0=artg(x)/x
c v_m(x)=sum_{n=0}^infty (-1)^n(2m+3)/(2m+2n+3)x^m
c v_m=1+(-1)^m(2m+3)/(2m+5)*x*v_{m+1}
c x>=0, x=x, artgm(m)=v_{m-1}
      subroutine dartgm(x,m)
      implicit real*8(a-h,o-z)
      common /warg/ artgm(30)
      data pis /1.570796326794897d0/
      data eps,zlim /1d-9,0.9d0/
      if(x.lt.0d0) goto 70
      mp1=m+1
      dmt=2d0*m+1d0
      if(x.gt.zlim) goto 40
      xs=-x*x
      vm=1d0
      dn=2d0
      xn=xs
      an=1d0
      do while (dabs(an).ge.eps)
         an=xn*dmt/(dmt+dn)
         vm=vm+an
         xn=xs*xn
         dn=dn+2d0
      enddo
      artgm(mp1)=vm
      if(m.eq.0) return
      do 20 mn=1,m
         mn1=mn-1
         dd=dmt-2d0*mn
         artgm(m-mn1)=1d0+xs*dd/(dd+2d0)*artgm(mp1-mn1)
   20 continue
      return
   40 continue
      xr=1d0/x
      if(x.le.1d0) then
          artgm(1)=datan(x)*xr
        else
          artgm(1)=(pis-datan(xr))*xr
      endif
      if(m.eq.0) return
      do 50 mn=1,m
         dmn=2d0*mn-1d0
         artgm(mn+1)=(1d0-artgm(mn))*(dmn+2d0)*xr*xr/dmn
   50 continue
      return
   70 continue
      write(6,2003) x
      return
 2003 format(1x,'negative x=',3x,d12.5)
      end
 
c programm calculates residuals of arth(x)/x series, v_0=arth(x)/x
c v_m(x)=sum_{n=0}^infty (2m+3)/(2m+2n+3)x^m
c v_m=1+(2m+3)/(2m+5)*x*v_{m+1}
c 0=<x<1, x=x,x1=1-x, arthm(m)=v_{m-1}
      subroutine darthm(x,x1,m)
      implicit real*8(a-h,o-z)
c      common /ww/ eps,zlim,ylim
      common /warh/ arthm(30)
      data eps,zlim /1d-9,0.9d0/
      if((x.lt.0d0).or.(x.ge.1d0)) goto 70
      mp1=m+1
      dmt=2d0*m+1d0
      if(x.gt.zlim) goto 40
      xs=x*x
      vm=1d0
      dn=2d0
      xn=xs
      an=1d0
      do while (an.ge.eps)
         an=xn*dmt/(dmt+dn)
         vm=vm+an
         xn=xs*xn
         dn=dn+2d0
      enddo
      arthm(mp1)=vm
      if(m.eq.0) return
      do 20 mn=1,m
         mn1=mn-1
         dnmt=2d0*mn
         dd=dmt-dnmt
         arthm(m-mn1)=1d0+xs*dd/(dd+2d0)*arthm(mp1-mn1)
   20 continue
      return
   40 continue
      xr=1d0/x
      arthm(1)=0.5d0*dlog((1d0+x)/x1)*xr
      if(m.eq.0) return
      do 50 mn=1,m
         dmn=2d0*mn-1d0
         arthm(mn+1)=(arthm(mn)-1d0)*(dmn+2d0)*xr*xr/dmn
   50 continue
      return
   70 continue
      write(6,2005) x
      return
 2005 format(1x,'must be 0=<x<1',3x,'x=',3x,d12.5)
      end

c subroutine di13ml(x) calculates for given argument x from [-8,8]
c the value of the Bessel function x^{1/3}I_{-1/3}(x) through its
c expansion on the Chebyshev polinomials of even order
c Luke, Special mathematical functions and their approximations, 
c ``MIR'', Moscow, 1980, 608 pp., table 9.18, p. 373
      function di13ml(x)
      implicit real*8(a-h,o-z)
      dimension a(19)
      data nout /6/
      data a / 246.76446 77610 21950 14476d0,
     1         374.90419 63977 61772 38356d0,
     2         168.94690 54665 51934 32342d0,
     3          47.81292 47970 09021 28986d0,
     4           9.02594 18911 35278 81180d0,
     5           1.19852 10167 30412 23176d0,
     6           0.11697 76116 93126 92166d0,
     7           0.00869 83095 97421 34176d0,
     8           0.00050 74091 61720 07577d0,
     9           0.00002 37873 65337 43276d0,
     ,           0.00000 09143 62383 07318d0,
     1           0.00000 00293 09947 86363d0,
     2           0.00000 00007 94846 32582d0,
     3           0.00000 00000 18462 55797d0,
     4           0.00000 00000 00371 27729d0,
     5           0.00000 00000 00006 52493d0,
     6           0.10105d-15,0.139d-17,0.2d-19/
c
      w=x*0.125d0
      di13ml=dcheb0(a,w,18,2,nout)
      return
      end
 
c subroutine di13pl(x) calculates for given argument x from [-8,8]
c the value of the Bessel function x^{-1/3}I_{1/3}(x) through its
c expansion on the Chebyshev polinomials of even order
c Luke, Special mathematical functions and their approximations, 
c ``MIR'', Moscow, 1980, 608 pp., table 9.18, p. 373
      function di13pl(x)
      implicit real*8(a-h,o-z)
      dimension a(18)
      data nout /6/
      data a / 65.17366 54884 16526 95482d0,
     1         95.23190 27040 92163 97860d0,
     2         39.59331 03435 21063 45067d0,
     3         10.20622 18597 97178 58861d0,
     4          1.75552 58978 34803 58155d0,
     5          0.21351 04236 00373 85749d0,
     6          0.01921 26868 62831 23642d0,
     7          0.00132 57975 71063 41284d0,
     8          0.00007 22075 54088 40446d0,
     9          0.00000 31776 13642 90919d0,
     ,          0.00000 01152 07184 65074d0,
     1          0.00000 00034 97935 67902d0,
     2          0.00000 00000 90183 77386d0,
     3          0.00000 00000 01998 07607d0,
     4          0.00000 00000 00038 43793d0,
     5          0.00000 00000 00000 64790d0,
     6          0.965d-17,0.13d-18/
c
      w=x*0.125d0
      di13pl=dcheb0(a,w,17,2,nout)
      return
      end

c subroutine dj13ml(x) calculates for given argument x from [-8,8]
c the value of the Bessel function x^{1/3}J_{-1/3}(x) through its
c expansion on the Chebyshev polinomials of even order
c Luke, Special mathematical functions and their approximations, 
c ``MIR'', Moscow, 1980, 608 pp., table 9.12, p. 365
      function dj13ml(x)
      implicit real*8(a-h,o-z)
      dimension a(19)
      data nout /6/
      data a / 0.15921 18677 83189 60054d0,
     1         0.13935 64007 14896 26069d0,
     2         0.05406 07319 77580 98869d0,
     3        -0.49889 25940 08176 73779d0,
     4         0.27634 28099 41437 09958d0,
     5        -0.06962 17491 14124 72070d0,
     6         0.01050 51479 21580 21916d0,
     7        -0.00107 36049 94092 01433d0,
     8         0.00007 98632 63462 05385d0,
     9        -0.00000 45371 68415 38583d0,
     ,         0.00000 02038 09350 93453d0,
     1        -0.00000 00074 32215 74519d0,
     2         0.00000 00002 24648 92145d0,
     3        -0.00000 00000 05724 18891d0,
     4         0.00000 00000 00124 68853d0,
     5        -0.00000 00000 00002 34948d0,
     6         0.00000 00000 00000 03868d0,
     7        -0.56d-18,0.1d-19/
c
      w=x*0.125d0
      dj13ml=dcheb0(a,w,18,2,nout)
      return
      end


c subroutine dj13pl(x) calculates for given argument x from [-8,8]
c the value of the Bessel function x^{-1/3}J_{1/3}(x) through its
c expansion on the Chebyshev polinomials of even order
c Luke, Special mathematical functions and their approximations, 
c ``MIR'', Moscow, 1980, 608 pp., table 9.12, p. 365
      function dj13pl(x)
      implicit real*8(a-h,o-z)
      dimension a(18)
      data nout /6/
      data a / 0.13437 97942 27756 44764d0,
     1        -0.11984 86974 10732 73261d0,
     2         0.28710 68687 28687 47721d0,
     3        -0.24251 07291 97801 31136d0,
     4         0.08569 82591 68484 14668d0,
     5        -0.01691 33949 03176 58537d0,
     6         0.00215 74934 76006 00765d0,
     7        -0.00019 38490 54788 73423d0,
     8         0.00001 29832 45927 14860d0,
     9        -0.00000 06748 56479 41628d0,
     ,         0.00000 00280 59239 71246d0,
     1        -0.00000 00009 55461 01654d0,
     2         0.00000 00000 27154 60828d0,
     3        -0.00000 00000 00654 22144d0,
     4         0.00000 00000 00013 53678d0,
     5        -0.00000 00000 00000 24324d0,
     6         0.00000 00000 00000 00383d0,
     7        -0.5d-19/
c
      w=x*0.125d0
      dj13pl=dcheb0(a,w,17,2,nout)
      return
      end

c subroutine dk13gs(x) calculates for given argument x>=5
c the value of the Bessel function exp(x)*sqrt(2*x/pi)K_{1/3}(x)
c through its expansion on the Chebyshev displased polinomials
c Luke, Special mathematical functions and their approximations, 
c ``MIR'', Moscow, 1980, 608 pp., table 9.18, p. 374
      function dk13gs(x)
      implicit real*8(a-h,o-z)
      dimension d(21)
      data nout /6/
      data d / 0.99353 64122 76093 38920d0,
     1        -0.00631 44392 60798 63137d0,
     2         0.00014 30095 80961 13131d0,
     3        -0.00000 57870 60592 02472d0,
     4         0.00000 03265 50333 19976d0,
     5        -0.00000 00231 23231 95077d0,
     6         0.00000 00019 39555 14434d0,
     7        -0.00000 00001 85897 88507d0,
     8         0.00000 00000 19868 42439d0,
     9        -0.00000 00000 02326 78966d0,
     ,         0.00000 00000 00294 68313d0,
     1        -0.00000 00000 00039 95293d0,
     2         0.00000 00000 00005 75225d0,
     3        -0.00000 00000 00000 87375d0,
     4  0.13927d-15,-0.2319d-16,0.402d-17,-0.72d-18,
     8  0.13d-18,-0.3d-19,0.1d-19/
c
      w=5d0/x
      dk13gs=dcheb0(d,w,20,3,nout)
      return
      end

c subroutine yj13gs(x) calculates for given argument |x|>5
c the combinations of values of the Bessel functions
c namely the real yjr and imaginary yji parts of the complex sum
c (\pi x/2)^{1/2}[J_{1/3}(x)+i Y_{1/3}(x)]e^{-i(x-\pi/4-\pi\nu/2)}
c through their expansions on the Chebyshev even polinomials
c Luke, Special mathematical functions and their approximations, 
c ``MIR'', Moscow, 1980, 608 pp., table 9.12, p. 365
      subroutine yj13gs(x,yjr,yji)
      implicit real*8(a-h,o-z)
      dimension br(22),bi(22)
      data nout /6/
      data br / 0.99946 50123 35293 52860d0,
     1         -0.00070 77638 36250 14260d0,
     2         -0.00016 87952 26919 56340d0,
     3          0.00000 43855 75455 07034d0,
     4          0.00000 03580 42553 40564d0,
     5         -0.00000 00464 10843 62258d0,
     6          0.00000 00007 21026 33596d0,
     7          0.00000 00004 73021 33710d0,
     8         -0.00000 00000 70556 54580d0,
     9          0.00000 00000 01600 09976d0,
     ,          0.00000 00000 01288 42696d0,
     1         -0.00000 00000 00282 67591d0,
     2          0.00000 00000 00021 03354d0,
     3          0.00000 00000 00004 44013d0,
     4         -0.00000 00000 00001 82827d0,
     5          0.00000 00000 00000 29922d0,
     6         -0.00000 00000 00000 00462d0,
     7         -0.1232d-16,0.398d-17,-0.63d-18,
     ,         -0.1d-19,0.4d-19/
      data bi /-0.00685 69474 27212 69534d0,
     1         -0.00681 43334 21755 53208d0,
     2          0.00005 01456 42019 67809d0,
     3          0.00000 71142 29548 86394d0,
     4         -0.00000 04367 92366 55846d0,
     5         -0.00000 00141 58283 43520d0,
     6          0.00000 00049 24417 94772d0,
     7         -0.00000 00003 454775 277630,
     8         -0.00000 00000 31070 98997d0,
     9          0.00000 00000 11008 31850d0,
     ,         -0.00000 00000 01162 71069d0,
     1         -0.00000 00000 00063 01191d0,
     2          0.00000 00000 00046 74357d0,
     3         -0.00000 00000 00008 46207d0,
     4          0.00000 00000 00000 37472d0,
     5          0.00000 00000 00000 23430d0,
     6         -0.00000 00000 00000 08025d0,
     7          0.1257d-16,0.9d-19,-0.68d-18,
     ,          0.22d-18,-0.4d-19 /
c
      w=5d0/x
      yjr=dcheb0(br,w,21,0,nout)
      yji=dcheb0(bi,w,21,0,nout)
      return
      end
 
c subprogram for calculation of gamma-function gamma(1+x) for real argum
c and gamma(m+x)=(m-1+x)(m-2+x)...(1+x)gamma(1+x)
      function dgam(xx)
      implicit real*8(a-h,o-z)
      dimension a(30)
      data xl,eps /3.2d0,1d-17/
      data a /1.00000 00000 00000 00000d0,
     1         0.57721 56649 01532 86061d0,
     2       -0.65587 80715 20253 88108d0,
     3       -0.04200 26350 34095 23553d0,
     4        0.16653 86113 82291 48950d0,
     5       -0.04219 77345 55544 33675d0,
     6       -0.00962 19715 27876 97356d0,
     7        0.00721 89432 46663 09954d0,
     8       -0.00116 51675 91859 06511d0,
     9       -0.00021 52416 74114 95097d0,
     ,        0.00012 80502 82388 11619d0,
     1       -0.00002 01348 54780 78824d0,
     2       -0.00000 12504 93482 14267d0,
     3        0.00000 11330 27231 98170d0,
     4       -0.00000 02056 33841 69776d0,
     5        0.00000 00061 16095 10448d0,
     6        0.00000 00050 02007 64447d0,
     7       -0.00000 00011 81274 57049d0,
     8        0.00000 00001 04342 67117d0,
     9        0.00000 00000 07782 26344d0,
     ,       -0.00000 00000 03696 80562d0,
     1        0.00000 00000 00510 03703d0,
     2       -0.20 58326d-13,-0.5 34812d-14,
     4        0.1 22678d-14,-0.11813d-15,
     6        0.119d-17,0.141d-17,
     8       -0.23d-18,0.2d-19/
c
      ddg=1d0
      dgam=0d0
      x=xx-1d0
      if(x.lt.xl) goto 300
      nxn=idint(x)
      x=x-nxn*1d0
      do nx=1,nxn
         ddg=ddg*(xx-nx*1d0)
      enddo
  300 continue
c      write(6,2000) nxn
c      write(6,2001) xx,x,ddg
      s=0d0
      xn=1d0
      do n=1,30
         sn=a(n)*xn
         s=s+sn
         if(dabs(sn).lt.eps) goto 330
         xn=x*xn
      enddo
  330 continue
      dgam=ddg/s
      return
 2000 format(1x,4i4)
 2001 format(1x,5d12.5)
      end

c subprogram calculates for given first n+1 elements of the vector a
c coefficients of the expansion on the Chebyshev polinomials and the
c argument x the following value:
c f(n,x)=sum_{k=0}^n a(k) T(l(key),x).
cHere the key
c  0  - subsequent integers,
c  1  -  odd numbers,
c  2  - even numbers,
c  3  - displased
c by using the Klenshoe summing algorithm
c nout --- the number of the canal of messages of the program
c Luke, Special mathematical functions and their approximations, 
      function dcheb0(a,x,n,key,nout)
      implicit real*8 (a-h,o-z)
      dimension a(n)
      data c0,c2,c4 /0d0,2d0,4d0/
      data cmax /1d76/
c control of the maximal number of terms
      ier=0
      if(n) 3,5,7
c the reaction to negative number
    3 write(nout,1000)
      dcheb0=cmax
      ier=1
c control the key 
    5 if((key.ge.0).and.(key.le.3)) goto 7
c reaction to the erroneous value of the key
      write(nout,1001)
      dcheb0=1d76
      ier=1
    7 if(ier.ne.0) return
      n1=n+1
      if(key.eq.0) z=c2*x
      if(key.eq.1) z=c4*x*x-c2
      if(key.eq.2) z=c4*x*x-c2
      if(key.eq.3) z=c4*x-c2
      b0=c0
      b1=c0
      do 10 i=1,n1
        k=n1-i
        b2=b1
        b1=b0
        b0=z*b1-b2+a(k+1)
   10 continue
      if(key.eq.0) dcheb0=b0-x*b1
      if(key.eq.1) dcheb0=x*(b0-b1)
      if(key.eq.2) dcheb0=b0-b1*z/c2
      if(key.eq.3) dcheb0=b0-b1*z/c2
      return
 1000 format(1x,'dcheb0: the third parameter is negative.')
 1001 format(1x,'dcheb0: erroneous value of the key')
      end

c subroutine calculates the function J_0(x) for given x
      function dj0(x)
      implicit real*8(a-h,o-z)
      data pif,pis /0.785398163397448d0,1.570796326794897d0/
      if(x.lt.6d0) then
          dj0=dj0l(x)
        else
          call yj0g(x,yjr0,yji0)
          xpi=x-pif
          dj0=(dcos(xpi)*yjr0-dsin(xpi)*yji0)/dsqrt(pis*x)
      endif
      return
      end

c subroutine dj0l(x) calculates for given argument x <-8=<x=<8
c the values of the Bessel function J_0(x)
c the real yjr and imaginary yji parts of the complex sum                       
c through its expansion on the Chebyshev even polinomials
c Luke, Special mathematical functions and their approximations, 
c ``MIR'', Moscow, 1980, 608 pp., table 9.1, p. 344
      function dj0l(x)
      implicit real*8(a-h,o-z)
      dimension a(18)
      data nout /6/
      data a /   0.15772 79714 74890 11956d0,
     1          -0.00872 34423 52852 22129d0,
     2           0.26517 86132 03336 80987d0,
     3          -0.37009 49938 72649 77903d0,
     4           0.15806 71023 32097 26128d0,
     5          -0.03489 37694 11408 88516d0,
     6           0.00481 91800 69467 60450d0,
     7          -0.00046 06261 66206 27505d0,
     8           0.00003 24603 28821 00508d0,
     9          -0.00000 17619 46907 76215d0,
     ,           0.00000 00760 81635 92419d0,
     1          -0.00000 00026 79253 53056d0,
     2           0.00000 00000 78486 96314d0,
     3          -0.00000 00000 01943 83469d0,
     4           0.00000 00000 00041 25321d0,
     5          -0.00000 00000 00000 75885d0,
     6           0.1222d-16,-0.17d-18 /
c
      w=x*0.125d0
      dj0l=dcheb0(a,w,17,2,nout)
      return
      end

c the subprogram yj0g(x,dre0,dim0) calculates for the given argument x>=5
c the values of the Bessel function:
c dre0=Re([J_0(x) +i*Y_0(x)]*sqrt(pi*x/2)*exp(i*(pi/4-x))) and
c dim0=Im([J_0(x)+i*Y(0,x)]*sqrt(pi*x/2)exp(i*(pi/4-x)))
c through their expansions of the even Chebyshev polinomials
c Luke, Special mathematical functions and their approximations, 
c ``MIR'', Moscow, 1980, 608 pp., table 9.1, p. 344
c if re and im are found then the functions J_0(x) and 
c with the formulae
cJ_0(x)=sqrt(2/(pi*x))*[dre0*cos(pi/4-x)+dim0*sin( pi/4-x)],
cY_0(x)=sqrt\(2/(pi*x))*[dim0*cos(pi/4-x )-dre0 * sin( pi/4-x )].
      subroutine yj0g(x,dre0,dim0)
      implicit real*8(a-h,o-z)
      dimension dr(23),di(22)
      data nout /6/
      data dr/ 0.99898 80898 58965 15390d0,
     1        -0.00133 84285 49971 85578d0,
     2        -0.00031 87898 78061 89289d0,
     3         0.00000 85112 32210 65665d0,
     4         0.00000 06915 42349 13894d0,
     5        -0.00000 00907 70101 53734d0,
     6         0.00000 00014 54928 07929d0,
     7         0.00000 00009 26762 48672d0,
     8        -0.00000 00001 39166 19797d0,
     9         0.00000 00000 03237 97518d0,
     ,         0.00000 00000 02535 35729d0,
     1        -0.00000 00000 00559 09032d0,
     2         0.00000 00000 00041 91896d0,
     3         0.00000 00000 00008 73316d0,
     4        -0.00000 00000 00003 61861d0,
     5 0.59438d-15,-0.964d-17,-0.2436d-16,0.789d-17,
     9-0.125d-17,-0.2d-19,0.8d-19,-0.3d-19 /
      data di/-0.01233 15205 78544 14382d0,
     1        -0.01224 94962 81259 47486d0,
     2         0.00009 64941 84993 42287d0,
     3         0.00001 36555 70490 35682d0,
     4        -0.00000 08518 06644 42635d0,
     5        -0.00000 00272 44053 41355d0,
     6         0.00000 00096 46421 33771d0,
     7        -0.00000 00006 83347 51799d0,
     8        -0.00000 00000 60627 38000d0,
     9         0.00000 00000 21695 71634d0,
     ,        -0.00000 00000 02304 89890d0,
     1        -0.00000 00000 00122 55390d0,
     2         0.00000 00000 00092 31372d0,
     3        -0.00000 00000 00016 77838d0,
     4         0.00000 00000 00000 75375d0,
     5 0.46244d-15,-0.15906d-15, 0.2500d-16,0.15d-18,
     9-0.135d-17, 0.44d-18,-0.7d-19 /
      w=5d0/x
      dre0=dcheb0(dr,w,22,3,nout)
      dim0=dcheb0(di,w,21,3,nout)
      return
      end

c the subprogram yj1g(x,dre1,dim1) calculates for the given argument x>=5
c the values of the Bessel function:
c dre1=Re([J_1(x)+i Y_1(x)]sqrt(pi*x/2)exp(i(pi/4-x))) and
c dim1=Im([J_0(x)+i Y(0,x)]sqrt(pi*x/2)exp(i(pi/4-x)))
c through their expansions of the even Chebyshev polinomials
c Luke, Special mathematical functions and their approximations, 
c ``MIR'', Moscow, 1980, 608 pp., table 9.2, p. 346
c if dre1 and dim1 are found then the functions J_1(x) and Y_1(x)
c can be calculated with the formulae
c J_1(x)=sqrt(2/(pi*x)) [dre1 cos(3 pi/4-x)+dim1*sin(3 pi/4-x)],
c Y_1(x)=sqrt\(2/(pi*x) [dim1*cos(3 pi/4-x)-dre1*sin(3 pi/4-x )].
      subroutine yj1g(x,dre1,dim1)
      implicit real*8(a-h,o-z)
      dimension dr(24),di(24)
      data nout /6/
      data dr/ 1.00170 22348 53820 99565d0,
     1         0.00225 55728 46561 17976d0,
     2         0.00054 32164 87508 01325d0,
     3        -0.00001 11794 61895 40836d0,
     4        -0.00000 09469 01382 39192d0,
     5         0.00000 01110 32677 12082d0,
     6        -0.00000 00012 94398 92684d0,
     7        -0.00000 00011 14905 94420d0,
     8         0.00000 00001 57637 23196d0,
     9        -0.00000 00000 02830 45747d0,
     ,        -0.00000 00000 02932 16857d0,
     1         0.00000 00000 00617 80854d0,
     2        -0.00000 00000 00043 16153d0,
     3        -0.00000 00000 00010 13289d0,
     4         0.00000 00000 00003 97343d0,
     5 -0.63161d-15, 0.575d-17, 0.2696d-16,-0.847d-17,
     9  0.130d-17, 0.3d-19,-0.9d-19, 0.3d-19,-0.1d-19 /
      data di/ 0.03726 17150 00537 65365d0,
     1         0.03714 53224 79807 68994d0,
     2        -0.00013 72632 38201 90679d0,
     3        -0.00001 98512 94687 59687d0,
     4         0.00000 10700 14057 38568d0,
     5         0.00000 00383 05261 71449d0,
     6        -0.00000 00116 28723 27663d0,
     7         0.00000 00007 59733 09244d0,
     8         0.00000 00000 75476 07460d0,
     9        -0.00000 00000 24752 78088d0,
     ,         0.00000 00000 02493 89256d0,
     1         0.00000 00000 00156 19784d0,
     2        -0.00000 00000 00103 38521d0,
     3         0.00000 00000 00018 12876d0,
     4        -0.00000 00000 00000 70876d0,
     5 -0.52042d-15, 0.17235d-15,-0.2625d-16,-0.38d-18,
     9  0.148d-17,-0.47d-18, 0.8d-19,0d0,-0.1d-19/
      w=5d0/x
      dre1=dcheb0(dr,w,23,3,nout)
      dim1=dcheb0(di,w,23,3,nout)
      return
      end