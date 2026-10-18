c      Subroutine capture_s solves classical equations for
c       relative motion coupled to surface deformation modes
c       specified through restoring force and mass parameters
c       Cin(i,lam,n) and Din(i,lam,n) in common/CD/ for projectile
c       (i=1) and target (i=2) for different multipolarities (lam=2,3).
c       The program also includes damping and non-linear restoring 
c       force, which is read or taken by default according to ierr-value.
c       Coulomb excitation is included or excluded according to ipt-value.
c      Reference: Broglia,Dasso and Winther: Proceedings of the 
c        Enrico Fermi International School of Physics 1979; "Nuclear
c        Structure and Heavy-ion Reactions" p.327, North Holland
c        Publ.Co., Amsterdam 1981.
c      To save time ipt = 1,2 for projectile (target) excitation only, 
c       ipt = 0 if both, Coulomb excitation may be neglected by setting
c       ipt = 3 for both(or 4,5 for projectile (target) excitation only).
c       The number of states included is limited by the parameter inelast.
c      If at entry iin = 0 it finds capture angular momentum (alcapt), 
c       i.e. the lowest impact parameter where the surface-surface 
c       distance (s) becomes negative, and capture distance (rcapt).
c       The orbital angular momentum at this point is alatr.
c       At exit ierr = 0 in this case.
c       If this situation does not occur (no capture), the bombarding
c        energy is increased until capture takes place, and
c        Ecm = ecmcap>ecm, r_capture= rcapt,  and l_capt=alcapt= 0.
c       In this case ierr = 3 at exit.
c       The program also evaluates the square of eta (eta2)
c        to be used for barrier distributions.
c       At high energy, where no real capture takes place
c        the program assumes capture when the half-densities
c        overlap during the collision.
c
c      If at entry iin = 1 the program calculates instead for r=rscat
c       (r is unperturbed distance of closest approach) the scattering 
c       angle (thscat), the minimum distance (rcapt) and the minimum
c       surface-surface distance (smin) as well as the total energy 
c       (hvt) and total l (alt) absorbed by each mode (i,lam,n).
c       If successful, then ierr = 0 at exit and alatr = final l.
c       If at exit ierr = 3, it means that the trajectory merged.
c       For iin = 1 no change of parameters is possible.
c
c      If at entry ierr = 0 s-damping is used with relative width 
c        of lowest state = .1, gnon = 0 and rotations are included.
c      If at entry ierr = 1  non-standard parameters may be read.
c-----------------------------------------------------------
      subroutine capture_s(ipt,iin,algraz,rgraz,rscat,rcapt,alatr,smi
     &,alcapt,ecmcapt,eta2,phido,rdud,thscat,alt,hvt
     &,istep,nestep,iout,inelast,ierr)
      implicit real*8(a-h,o-z)
      external funcap,f_out
      common/za/za,aa,zaa,aaa
      common/ecm/ecm
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      common/capt/eb,rbar,i5,intp
      common/giantwid/giantwid(3)
      common/CD/Cin(2,3,2),Din(2,3,2)
      common/char/elemp,elemt
      common/ebe23/eex(2,3,2)
      common/head/head
      common/captureloc/ Y20,Y22,Y31,Y33,am0,vel,alc,rac,raac,sdamp
     &,gamma,wid,gnon,ainert,aainert,ddelp,ee(2),rmini,smini,alfz(2,3,2)
     &,reta(2,3,2,4),aeta(2,3,2,4),rdur(2,3),phid
     &,nmin,nmax,nnn,iptt,ieta,iiin
      dimension hvib(2,3,2),hpot(2,3,2),hkin(2,3,2),alvib(2,3,2)
     &,alt(2,3,2),hvt(2,3,2),altt(8),hvtt(8),eta2(2,3,2,4),rdud(2,3)
      dimension array(200,7),z(200),zdot(200)
      character*5 elemp,elemt
c
      call constant_s
      if(ipt.gt.2)then
         r0c = 0.0d0
         iptt = ipt - 3
      else
         r0c = 1.2d0
         iptt = ipt
      endif
      istand = ierr
      iiin = iin
      ieta = 0
c
c     so that orbital integrals are only calculated for last trajectory
      am0 = aa*aaa/(aa+aaa)
      rrgraz = rgraz+1.d0
      ierr1 = 1
      delalc = alofr(rrgraz,ierr1)-algraz
      if(ierr1.eq.1)goto 2
      almax = 0.d0
      almin = -1.d0
      ecmmax = eb*1.05d0
      ecmmin = ecm
      ee(1) = -1.d0
      ee(2) = 1.d0
      if(inelast.eq.-1)then
         write(*,*)' Capture_s skipped: no inelastic modes'
         ierr = 0
         return
        elseif(inelast.eq.0)then 
         nmin = 2
         nmax = 2
        elseif(inelast.eq.1)then 
         nmin = 1
         nmax = 1
        elseif(inelast.eq.2)then 
         nmin = 1
         nmax = 2
      endif
c
      if(iin.eq.1)then
         ierr = 1
         alc = alofr(rscat,ierr)
         if(ierr.eq.1)goto 2
         goto 57
      endif
c
      if(nestep.eq.1.and.istep.eq.1)then    
         sdamp = -1.d0
         gamma = 0.d0
         wid = .1d0
         gnon = 0.d0
         ddelp = 1.d0
      endif   
c
      write(*,*)' Calculating capture angular momentum'
      if(ecm.lt.eb)then
         ibar = 1
         alc = 0.d0
        else
         ibar = 0
         alc = algraz
      endif
c
   57 continue
      rac = aa3*r0c
      raac = aaa3*r0c
      ecm1 = ecm
      ecmstore = ecm
      timeinit = X05BAF()
      notraj = 0
      ierr1 = 0
      call ylm_s(Y20,2,0,ierr1)
      if(ierr1.eq.1)goto 2
      ierr1 = 0
      call ylm_s(Y22,2,2,ierr1)
      if(ierr1.eq.1)goto 2
      ierr1 = 0
      call ylm_s(Y31,3,1,ierr1)
      if(ierr1.eq.1)goto 2
      ierr1 = 0
      call ylm_s(Y33,3,3,ierr1)
      if(ierr1.eq.1)goto 2
      uup = 1.d0
c  
c  Zeropoint fluctuations:
      do i = 1,2
       do lam=2,3
        do n=nmin,nmax
         if(Cin(i,lam,n).lt.5.d4)then
          alfz(i,lam,n)=sqrt(eex(i,lam,n)/(2.d0*Cin(i,lam,n)+1.d-10))  
         else
          alfz(i,lam,n)=0.d0
         endif
        enddo
       enddo
      enddo
      delw = .2d0
c     (step in w between checks)
      wstop = 3.d0
c     (to stop if distance remains small)
      scut = 0.d0
c     (minimum s to decide whether capture)
      rminn = 100.d0
      sminn = 100.d0
c
c  Here enter for new trajectory
c
   58 continue
c
      notraj = notraj+1
      time1 = X05BAF()
      vel = sqrt(2.d0*ecm1*hc2om/am0)
      eta = za*zaa*esq/vel
      epsilon = sqrt(1.d0+(alc/eta)**2)
      a0 = za*zaa*esq/2.d0/ecm1
c
      rinit = sqrt(eta/accur)*uup/10.d0
c     Calculate initial values including Coulomb deflection at rinit: 
      yy = (rinit-a0)/a0/epsilon
      xxx = yy**2-1.d0
      if(xxx.lt.0.d0)xxx=0.d0
      awinit = log(yy+sqrt(xxx))
      winitin = -dble(int(awinit/delw)+1)*delw
c
      winit = winitin
      tinit = a0/vel*(epsilon*sinh(winit)+winit)
      argg = epsilon**2-1.d0
      if(argg.lt.0.d0)argg = 0.d0
      phiinit =   asin(sqrt(argg)/epsilon)
     &          + asin(sqrt(argg)*sinh(winit)
     &          /(epsilon*cosh(winit)+1.d0))
      ierr11 = 1
      aarg = 2.d0*am0/hc2om*(ecm1-ueff(alc,rinit,ierr11))
      if(aarg.lt.0.d0)aarg = 0.d0
      pr0 = sqrt(aarg)
      if(ierr11.eq.1)goto 2
c
c                                NO!
c
c List of the variables:    
c for trajectory:
c             z(1) =time, z(2) =r, z(3) =pr, z(4) =phi, z(5) =pphi=l,
c no deform.  z(6) =time, z(7) =r, z(8) =pr, z(9) =phi, z(10) =pphi=l,
c for the modes of i=1(projectile):  n=1(low modes)  (alf = alpha)
c   z(11)=alf20,z(12)=RE(alf22),z(13)=IM(alf22),z(14)=0, z(15)=0
c   z(16)=pi20, z(17)=RE(pi22), z(18)=IM(pi22), z(19)=0, z(20)=0
c z(21)=RE(alf31),z(22)=IM(alf31),z(23)=RE(alf33),z(24)=IM(alf3),z(25)=0
c   z(26)=RE(pi31),z(27)=IM(pi31),z(28)=RE(pi33),z(29)=IM(pi33),z(30)=0
c                  i=1        n=2(high modes)
c alf2   z(31)=     , z(32)=     , z(33)=     , z(34)=     , z(35)=
c pi2    z(36)=     , z(37)=     , z(38)=     , z(39)=     , z(40)=
c alf3   z(41)=     , z(42)=     , z(43)=     , z(44)=     , z(45)=
c pi3    z(46)=     , z(47)=     , z(48)=     , z(49)=     , z(50)=
c modes for target i=2       n=1(low modes)
c   z(51)=alf20,z(52)=RE(alf22),z(53)=IM(alf22),z(54)=0, z(55)=0
c   z(56)=pi20, z(57)=RE(pi22), z(58)=IM(pi22), z(59)=0, z(60)=0
c z(61)=RE(alf31),z(62)=IM(alf31),z(63)=RE(alf33),z(64)=IM(alf3),z(65)=0
c   z(66)=RE(pi31),z(67)=IM(pi31),z(68)=RE(pi33),z(69)=IM(pi33),z(70)=0
c                  i=2       n=2(high modes)
c alf2   z(71)=     , z(72)=     , z(73)=     , z(74)=     , z(75)=
c pi2    z(76)=     , z(77)=     , z(78)=     , z(79)=     , z(80)=
c alf3   z(81)=     , z(82)=     , z(83)=     , z(84)=     , z(85)=
c pi3    z(86)=     , z(87)=     , z(88)=     , z(89)=     , z(90)=
c
c       z(91)=la_dissip, z(92)=laa_dissip,  (dissipated angular momentum)
c       z(93)=hviba_dissip, z(94)=hvibaa_dissip, (dissipated energy)
c       z(95)= 0, z(97)=delphiaa, z(96)=delphia   (rotation angle)
c       z(98)= 0, z(99)= 0, z(100) = 0,
c       z(101)=aldis(1,2,1),z(102)=aldis(1,2,2),z(103)=aldis(1,3,1),
c       z(104)=aldis(1,3,2),z(105)=aldis(2,2,1),z(106)=aldis(2,2,2),
c       z(107)=aldis(2,3,1),z(108)=aldis(2,3,2),z(109)= 0,z(110)= 0
c       z(111)=hdis(1,2,1),z(112)=hdis(1,2,2),z(113)=hdis(1,3,1),
c       z(114)=hdis(1,3,2),z(115)=hdis(2,2,1),z(116)=hdis(2,2,2),
c       z(117)=hdis(2,3,1),z(118)=hdis(2,3,2),z(119)= 0,z(120)= 0
c     orbital integrals
c       z(121)=RE(orb(1,2,1,-2),z(122)=IM(orb(1,2,1,-2)),
c       z(123)=RE(orb(1,2,1,0),z(124)=IM(orb(1,2,1,0)),
c       z(125)=RE(orb(1,2,1,2),z(126)=IM(orb(1,2,1,2)),
c       z(127)=0,z(128)=0,z(129)=0,z(130)=0,
c       z(131)=RE(orb(1,3,1,-3),z(132)=IM(orb(1,3,1,-3)),
c       z(133)=RE(orb(1,3,1,-1),z(134)=IM(orb(1,3,1,-1)),
c       z(135)=RE(orb(1,3,1,1),z(136)=IM(orb(1,3,1,1)),
c       z(136)=RE(orb(1,3,1,3),z(138)=IM(orb(1,3,1,3)),
c       z(139)=0,z(140)=0,
c       z(141)=RE(orb(1,2,2,-2),z(142)=IM(orb(1,2,2,-2)),
c       z(143)=RE(orb(1,2,2,0),z(144)=IM(orb(1,2,2,0)),
c       z(145)=RE(orb(1,2,2,2),z(146)=IM(orb(1,2,2,2)),
c       z(147)=0,z(148)=0,z(149)=0,z(150)=0,
c       z(151)=RE(orb(1,3,2,-3),z(152)=IM(orb(1,3,2,-3)),
c       z(153)=RE(orb(1,3,2,-1),z(154)=IM(orb(1,3,2,-1)),
c       z(155)=RE(orb(1,3,2,1),z(156)=IM(orb(1,3,2,1)),
c       z(156)=RE(orb(1,3,2,3),z(158)=IM(orb(1,3,2,3)),
c       z(159)=0,z(160)=0,
c       z(161)=RE(orb(2,2,1,-2),z(162)=IM(orb(2,2,1,-2)),
c       z(163)=RE(orb(2,2,1,0),z(164)=IM(orb(2,2,1,0)),
c       z(165)=RE(orb(2,2,1,2),z(166)=IM(orb(2,2,1,2)),
c       z(167)=0,z(168)=0,z(169)=0,z(170)=0,
c       z(171)=RE(orb(2,3,1,-3),z(172)=IM(orb(2,3,1,-3)),
c       z(173)=RE(orb(2,3,1,-1),z(174)=IM(orb(2,3,1,-1)),
c       z(175)=RE(orb(2,3,1,1),z(176)=IM(orb(2,3,1,1)),
c       z(176)=RE(orb(2,3,1,3),z(178)=IM(orb(2,3,1,3)),
c       z(179)=0,z(180)=0,
c       z(181)=RE(orb(2,2,2,-2),z(182)=IM(orb(2,2,2,-2)),
c       z(183)=RE(orb(2,2,2,0),z(184)=IM(orb(2,2,2,0)),
c       z(185)=RE(orb(2,2,2,2),z(186)=IM(orb(2,2,2,2)),
c       z(187)=0,z(188)=0,z(189)=0,z(190)=0,
c       z(191)=RE(orb(2,3,2,-3),z(192)=IM(orb(2,3,2,-3)),
c       z(193)=RE(orb(2,3,2,-1),z(194)=IM(orb(2,3,2,-1)),
c       z(195)=RE(orb(2,3,2,1),z(196)=IM(orb(2,3,2,1)),
c       z(196)=RE(orb(2,3,2,3),z(198)=IM(orb(2,3,2,3)),
c       z(199)=0,z(200)=0
c
      aainert = .4d0/hc2om*aaa*raa**2
      ainert = .4d0/hc2om*aa*ra**2
cf
      nmo=nmax-nmin+1
      neq=16+2*nmo*(14+2+2+14)
      neq1=16+2*nmo*14+1
cf
      z(1) = tinit
      z(6) = tinit
      z(2)=rinit
      z(3)=-pr0
      z(4)=phiinit
      z(5)=alc
      z(7)=rinit
      z(8)=-pr0
      z(9)=phiinit 
      z(10)=alc
      do n=11,neq
         z(n)=0.d0
         zdot(n) = 0.d0
      enddo
      phidg = phiinit*180.d0/pi
      s0 = z(2)-ra-raa
      smin = s0
      smini = s0
      rmin = rinit
      rmini = rinit
c
c  return here for new step in integration:
c
   29 continue
c
      tol = accur          !/100.d0
      ifail = 0
      wfin = winit+delw
      w = winit
c
c      call D02BAF(winit,wfin,neq,z,tol,funcap,array,ifail)
      call D02BBF(winit,wfin,neq,z,tol,0,funcap,f_out,array,ifail)
c
      do i=1,2
         do lam=2,3
            do n=nmin,nmax
               hpot(i,lam,n) = 0.d0
               hkin(i,lam,n) = 0.d0
               hvib(i,lam,n) = 0.d0
               alvib(i,lam,n) = 0.d0
            enddo
         enddo
      enddo
c
      madd=0
      do i=1,2
         do n=nmin,nmax
            hpot(i,2,n)=Cin(i,2,n)*(z(17+madd)**2/2.d0+z(19+madd)**2
     &                 +z(21+madd)**2)
            hkin(i,2,n)=(z(18+madd)**2/2.d0+z(20+madd)**2
     &                  +z(22+madd)**2)/Din(i,2,n)
            hvib(i,2,n)= hpot(i,2,n)+ hkin(i,2,n)
            alvib(i,2,n)=-4.d0*(z(19+madd)*z(22+madd)
     &                  +z(21+madd)*z(20+madd))
            hpot(i,3,n)=Cin(i,3,n)*(z(23+madd)**2+z(25+madd)**2
     &                 +z(27+madd)**2+z(29+madd)**2)
            hkin(i,3,n)=(z(24+madd)**2+z(26+madd)**2
     &                 +z(28+madd)**2+z(30+madd)**2)/Din(i,3,n)
            hvib(i,3,n)= hpot(i,3,n)+hkin(i,3,n)
            alvib(i,3,n) =-2.d0*(z(23+madd)*z(26+madd)
     &                   +z(25+madd)*z(24+madd))
     &                   -6.d0*(z(27+madd)*z(30+madd)
     &                   +z(29+madd)*z(28+madd))
            madd=madd+14
         enddo
      enddo
c
c     Energy and angular momentum absorbed by each mode:
      mm=neq1
      nn=1
      do i=1,2
         do n=nmin,nmax
            do lam=2,3
               alt(i,lam,n) = alvib(i,lam,n)+z(mm)
               hvt(i,lam,n) = hvib(i,lam,n)+z(mm+2)
c
               altt(nn) = alt(i,lam,n)
               hvtt(nn) = hvt(i,lam,n)
               mm=mm+1
               nn=nn+1
            enddo
            mm=mm+2
         enddo
      enddo
c
c  quantities for printing and check:
      alviba = alvib(1,2,1)+alvib(1,2,2)+alvib(1,3,1)+alvib(1,3,2)
      alvibaa = alvib(2,2,1)+alvib(2,2,2)+alvib(2,3,1)+alvib(2,3,2)
      alvibt = alviba+alvibaa
      aldisaa = z(12)
      aldisa = z(11)
      altot = z(12)+z(11)+z(5)+alvibt
      hpota = hpot(1,2,1)+hpot(1,3,1)+hpot(1,2,2)+hpot(1,3,2)
      hpotaa = hpot(2,2,1)+hpot(2,3,1)+hpot(2,2,2)+hpot(2,3,2)
      if(inelast.ge.1)then 
         hpota = hpota+(hpot(1,2,1)**2/Cin(1,2,1)
     &           +hpot(1,3,1)**2/Cin(1,3,1))*gnon
         hpotaa= hpotaa+(hpot(2,2,1)**2/Cin(2,2,1)
     &           +hpot(2,3,1)**2/Cin(2,3,1))*gnon
      endif
      hkina = hkin(1,2,1)+hkin(1,3,1)+hkin(1,2,2)+hkin(1,3,2)
      hkinaa = hkin(2,2,1)+hkin(2,3,1)+hkin(2,2,2)+hkin(2,3,2)
      hviba = hpota + hkina
      hvibaa = hpotaa + hkinaa
      hvibt = hvibaa + hviba
      hdisaa = z(14)
      hdisa = z(13)
      hdis = hdisaa+hdisa
c
c     compute surface-surface distance
      delphia = z(15)
      delphadg = z(15)*180.d0/pi
      delphiaa = z(16)
      delphaadg = z(16)*180.d0/pi
      phidg = z(4)*180.d0/pi
      s0 = z(2) - ra - raa
      phi = z(4)
      phicora = phi - delphia
      phicoraa = phi - delphiaa
      cospa =cos(phicora)
      sinpa =sin(phicora)
      cos2pa =cos(2.d0*phicora)
      sin2pa =sin(2.d0*phicora)
      cos3pa =cos(3.d0*phicora)
      sin3pa =sin(3.d0*phicora)
      cospaa =cos(phicoraa)
      sinpaa =sin(phicoraa)
      cos2paa =cos(2.d0*phicoraa)
      sin2paa =sin(2.d0*phicoraa)
      cos3paa =cos(3.d0*phicoraa)
      sin3paa =sin(3.d0*phicoraa)
      dr2a = 0.d0
      dr3a = 0.d0
      dr2aa = 0.d0
      dr3aa = 0.d0
      nadd=0
      do n=nmin,nmax
         dr2a = dr2a+z(17+nadd)*Y20+2.d0*Y22*(z(19+nadd)*cos2pa 
     &          +z(21+nadd)*sin2pa)
         dr3a = dr3a+2.d0*Y31*(z(23+nadd)*cospa+z(25+nadd)*sinpa)
     &          +2.d0*Y33*(z(27+nadd)*cos3pa+z(29+nadd)*sin3pa)
         dr2aa = dr2aa+z(17+14*nmo+nadd)*Y20+2.d0*Y22*(z(19+14*nmo+nadd)
     &          *cos2paa+z(21+14*nmo+nadd)*sin2paa)
         dr3aa = dr3aa+2.d0*Y31*(z(23+14*nmo+nadd)*cospaa
     &          +z(25+14*nmo+nadd)*sinpaa)+2.d0*Y33*(z(27+14*nmo+nadd)
     &          *cos3paa+z(29+14*nmo+nadd)*sin3paa)
         nadd=nadd+14
      enddo
      dra = (dr2a - dr3a)*ra
      draa = (dr2aa + dr3aa)*raa
      s = s0-dra-draa 
c
c  calculate H_tot and T_kin
      aex = s*ainv/2.d0
      aexx = (z(7)-ra-raa)*ainv/2.d0
      if(aex.gt.40)then
         dundr = 0.d0
         un = 0.d0
        else
         eexp = exp(aex)
         un = -v0/(1.d0+eexp**2)
         dundr = v0*ainv/(eexp+1.d0/eexp)**2 
      endif
      if(aexx.gt.40)then
         dudrr = 0.d0
         unu = 0.d0
        else
         eeexp = exp(aexx)
         unu = -v0/(1.d0+eeexp**2)
         dundrr = v0*ainv/(eeexp+1.d0/eeexp)**2 
      endif      
      uc = za*zaa*esq/z(2)
      ucc2aa = uc*3.d0/5.d0*(raac/z(2))**2
      ucc3aa = uc*3.d0/7.d0*(raac/z(2))**3
      ucc2a = uc*3.d0/5.d0*(rac/z(2))**2
      ucc3a = uc*3.d0/7.d0*(rac/z(2))**3
      uc2aa = ucc2aa*dr2aa
      uc3aa = ucc3aa*dr3aa
      uc2a = ucc2a*dr2a
      uc3a = ucc3a*dr3a
      ur = un+uc+uc2aa+uc3aa+uc2a-uc3a
      dudr=dundr-(uc+3.d0*(uc2aa+uc2a)+4.d0*(uc3aa-uc3a))/z(2)
      tkin = (z(3)**2 + (z(5)/z(2))**2)/2.d0/am0*hc2om
      ecmf = ur + tkin
      etot = hvibt + ecmf + hdis
c
c     Energy conservation (cannot be checked if rotation is included) 
      if(ddelp.eq.0.d0.and.abs(etot-ecm1)/ecm1.gt.10.d0*accur)then
       write(*,'(a,f8.4,a,f7.4)')
     &' Warning: Bad energy conservation at l = ', alc,' and w = ',wfin
       write(*,*)'      energy = ',etot,' should be ',ecm1     
      endif
      if(ddelp.eq.0.d0.and.abs(altot-alc)/alc.gt.10.d0*accur)then
       write(*,'(a,f8.4,a,f7.4)')
     &' Warning: Bad ang.mom. conservation at l = ',alc,' and w = ',wfin
       write(*,*)'     l = ',altot,' should be ',alc
      endif
      if(s.le.smin)then
         smin = s
         ala = z(5)
      endif
      if(z(2).lt.rmin)rmin = z(2)
      if(iin.eq.0.and.s.gt.smin+1.d0)then
         alatr = ala
         goto 30
        elseif(iin.eq.1.and.z(2).gt.rinit)then 
         alatr = z(5)
         goto 30
c        (there is no capture !)
      endif
      if(s.lt.scut.or.wfin.gt.wstop.and.s.lt.2.d0)then
c        (there is capture!)
         if(iin.eq.1)then
            if(wfin.lt.wstop)goto 59
c           no capture should have occured
            ierr = 3
            smi = smini
            rmi = rmini
            ecm = ecmstore
            return
         endif
         if(ibar.eq.1.and.abs(ecm1-ecmmin).gt.accur)then
            goto 60
           elseif(ibar.eq.1.and.abs(ecm1-ecmmin).lt.accur)then
            ibar = 0
            almin = alc
            if(almax.eq.0.)then
               alc = alc + delalc
               goto 58
              elseif(almax.gt.accur)then
               alc = (almax+alc)/2.d0
               goto 3
            endif
           else
            almin = alc
            if(almax.eq.0.)then
               alc = alc + delalc
               goto 58
              elseif(almax.gt.accur)then
               if(ieta.eq.0) alc = (almax+alc)/2.d0
               goto 3
            endif
         endif
      endif
   59 continue
      time2 = X05BAF()
      timestep = time2 - time1
      winit = wfin
c
      goto 29
c
c
c     here if trajectory finished (no capture):
   30 continue
      if(rmini.lt.rminn)rminn = rmini
c     rmini is smallest r for this scattering trajectory
      if(smini.lt.sminn) sminn = smini
c     smini is smallest s for this trajectory (not at the smallest r!)
c
      if(iin.eq.0)goto 46 
c
c     Calculate scattering angle
      velf = sqrt(2.d0*ecmf*hc2om/am0)
      etaf = za*zaa*esq/velf
      epsilonf = sqrt(1.d0+(z(5)/etaf)**2)
      a0f = za*zaa*esq/2.d0/ecmf
c
      rfin = z(2)
      yy = (rfin-a0f)/a0f/epsilonf
      wfin = log(yy+sqrt(yy**2-1.d0))
      tfin = a0f/velf*(epsilonf*sinh(wfin)+wfin)
      dphif = asin(sqrt(epsilonf**2-1.d0)/epsilonf)
     &        - asin(sqrt(epsilonf**2-1.d0)*sinh(wfin)/
     &          (epsilonf*cosh(wfin)+1.d0))
      rfinu = z(7)
      yyu = (rfinu-a0)/a0/epsilon
      ssq = yyu**2-1.d0
      if(ssq.lt.0.)then
         wfinu = log(yyu)
        else       
         wfinu = log(yyu+sqrt(ssq))
      endif
      tfinu = a0/vel*(epsilon*sinh(wfinu)+wfinu)
      dphifu =  asin(sqrt(epsilon**2-1.d0)/epsilon)
     &       - asin(sqrt(epsilon**2-1.d0)*sinh(wfinu)/
     &       (epsilon*cosh(wfinu)+1.d0))
c
      thscat = pi - (z(4) + dphif)
      thetadg = 180.d0/pi*thscat
      thetaundef = pi - (z(9) + dphifu)
      thetaundefdg = 180.d0/pi*thetaundef
      thetacoul = 2.d0*asin(1.d0/epsilon)
      thetacouldg = 180.d0/pi*thetacoul
c
      rcapt = rmini
      smi = smini
      ierr = 0
      ecm = ecmstore
c
      return
c
   46 continue
      if(ibar.eq.1)then
         ecmmin = ecm1
         ecmmintemp = ecm1
         goto 61
        else
         almax = alc
         almaxtemp = alc
      endif
cnan      write(232,'(a,1p3e12.4)') ' --- ',ecm1,alc,almin
      if(alc.eq.0.d0)then
         ibar = 1
         ecmmin = ecm1
         ecmmintemp = ecm1
         goto 61
      endif
      if(almin.ge.0.d0)then
         alc = (alc+almin)/2.d0
         goto 3
        else
         alc = alc - delalc
         if(alc.lt.0.d0)alc = 0.d0
         goto 58
      endif
   60 continue
      ecmmax = ecm1
   61 continue
      ecm1 = (ecmmin+ecmmax)/2.d0
cc      if((ecmmax-ecmmin)).gt.ecm1*5.d0*accur)then
      if(((ecmmax-ecmmin)/ecm1).gt.accur) then
         goto 58  
        else
         if(ieta.eq.0)then
            ecm1 = ecmmintemp*(1.d0-accur)
            ecmtemp = ecm1
            ieta = 1
            goto 58
           elseif(ieta.eq.1)then
            write(*,*)' calculated eta for emc = ',ecmtemp
            phido = phid
            do i=1,2
               do lam=2,3
                  do n=nmin,nmax
                     do m=1,4
                        eta2(i,lam,n,m) = reta(i,lam,n,m)**2
     &                                   +aeta(i,lam,n,m)**2
                     enddo 
                  enddo
                  rdud(i,lam)=rdur(i,lam)
               enddo
            enddo
         endif
   62    write(*,*)' Energy too low for capture'
         ecmcapt = (ecmmin+ecmmax)/2.d0
         alcapt = .0d0
         rcapt = rmini
         smi = smini
         ecm = ecmstore
         ierr = 3
         goto 10
      endif
c
    3 continue
c
      if(ieta.eq.0)then
         if(almax-almin.gt.max(almax,10.d0)*5.d0*accur) goto 58
cnan         alc = almaxtemp*(1.d0+accur)
         alc = almaxtemp*(1.d0+100*accur)
cnan         write(232, *) alc
         altemp = alc
         ieta = 1
         goto 58
        elseif(ieta.eq.1)then
         write(*,*)' calculated eta for l = ',altemp
cnan         write(232,*)' calculated eta for l = ',altemp
         phido = phid
         do i=1,2
            do lam=2,3
               do n=nmin,nmax
                  do m=1,4
                     eta2(i,lam,n,m) = reta(i,lam,n,m)**2
     &                                +aeta(i,lam,n,m)**2
                  enddo 
               enddo
               rdud(i,lam)=rdur(i,lam)
            enddo
         enddo
      endif
      rcapt = rminn
      smi = sminn
      alcapt = alc
      ecm = ecmstore
      ierr = 0
   10 continue
      timefin = X05BAF()
      timemerg = timefin - timeinit
      write(*,'(a,f7.2,a)')'  CPU time for capture= ',timemerg,' sec.'
      return
c
    2 continue
c
      ierr = 1
      ecm = ecmstore
      return
      end        
c
c     -------------
c
      subroutine funcap(w,z,zdot)
      implicit real*8(a-h,o-z)
      dimension z(200),zdot(200)
      common/za/za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      common/giantwid/giantwid(3)
      common/CD/Cin(2,3,2),Din(2,3,2)
      common/ebe23/eex(2,3,2)
      common/captureloc/ Y20,Y22,Y31,Y33,am0,vel,alc,rac,raac,sdamp
     &,gamma,wid,gnon,ainert,aainert,ddelp,ee(2),rmini,smini,alfz(2,3,2)
     &,reta(2,3,2,4),aeta(2,3,2,4),rdur(2,3),phid
     &,nmin,nmax,nnn,iptt,ieta,iiin
      dimension hvib(2,3,2),hpot(2,3,2),hkin(2,3,2),alvib(2,3,2)
     &,gam(2,3,2),ucc(2,3),rdudr(2),zdotcal(2),zdotch(2),phi(2),rdu(2,3)
     &,sinp(3,2),cosp(3,2)
cf
      nmo=nmax-nmin+1
      neq=16+2*nmo*(14+2+2+14)
      neq1=16+2*nmo*14+1
      neq2=16+2*nmo*(14+2+2)
cf
      iin = iiin
      w = w
      zdot(1) = z(2)/vel
      do mm = neq1,neq
       zdot(mm) = 0.d0
      enddo
      madd=0
      mm=neq1
      do i=1,2
       zdotcal(i)=0.d0                                 
       zdotch(i)=0.d0
       do n=nmin,nmax
        hpot(i,2,n)=Cin(i,2,n)*(z(17+madd)**2/2.d0+z(19+madd)**2
     &+z(21+madd)**2)
        hkin(i,2,n) = (z(18+madd)**2/2.d0+z(20+madd)**2
     &+z(22+madd)**2)/Din(i,2,n)
        hvib(i,2,n)= hpot(i,2,n)+ hkin(i,2,n)
        alvib(i,2,n)=-4.d0*(z(19+madd)*z(22+madd)+z(21+madd)*z(20+madd))
        hpot(i,3,n)=Cin(i,3,n)*(z(23+madd)**2+z(25+madd)**2
     &+z(27+madd)**2+z(29+madd)**2)
        hkin(i,3,n)=(z(24+madd)**2+z(26+madd)**2
     &+z(28+madd)**2+z(30+madd)**2)/Din(i,3,n)
        hvib(i,3,n)= hpot(i,3,n)+hkin(i,3,n)
        alvib(i,3,n) =  -2.d0*(z(23+madd)*z(26+madd)
     &+z(25+madd)*z(24+madd))-6.d0*(z(27+madd)*z(30+madd)
     &+ z(29+madd)*z(28+madd))
        do lam = 2,3
         if(n.eq.2)then
          gam(i,lam,2) = giantwid(lam)
         elseif(sdamp.ne.0..and.gamma.eq.0..and.wid.eq.0.)then
          gam(i,lam,1) = sdamp*hvib(i,lam,1)
         elseif(gamma.ne.0..and.sdamp.eq.0..and.wid.eq.0.)then
          gam(i,lam,1) = gamma
         elseif(sdamp.ne.0..and.gamma.eq.0..and.wid.ne.0.)then
          gam(i,lam,1) = wid*hvib(i,lam,1)*2.d0*Din(i,lam,1)
         elseif(gamma.ne.0..and.sdamp.eq.0..and.wid.ne.0.)then
          gam(i,lam,1) = wid*2.d0*sqrt(Cin(i,lam,1)*Din(i,lam,1))
         elseif(gamma.eq.0..and.sdamp.eq.0.)then
          gam(i,lam,1) = 0.d0
         elseif(gamma.ne.0..and.sdamp.ne.0.)then
          write(*,*)' Error gamma.ne.0..and.sdamp.ne.0.'
         endif
         if(Din(i,lam,n).lt.5.d4)then
          zdot(mm) = gam(i,lam,n)/Din(i,lam,n)*alvib(i,lam,n)
          zdotcal(i) = zdotcal(i)+ zdot(mm)
          zdot(mm+2) = 2.d0*gam(i,lam,n)/Din(i,lam,n)*hkin(i,lam,n)
          zdotch(i) = zdotch(i)+ zdot(mm+2)
         endif
         mm=mm+1
        enddo
        mm=mm+2
        madd=madd+14
       enddo
      enddo
      phi(1) = z(4)-z(15)
      phi(2) = z(4)-z(16)
      do i = 1,2
       do mu = 1,3
        cosp(mu,i) = cos(dble(mu)*phi(i))
        sinp(mu,i) = sin(dble(mu)*phi(i))
       enddo
      enddo
      dr2a = 0.d0
      dr2aa = 0.d0
      dr3a = 0.d0
      dr3aa = 0.d0
      nadd=0
      do n=nmin,nmax
       dr2a = dr2a
     & + Y20*z(17+nadd)+2.d0*Y22*(z(19+nadd)*cosp(2,1)
     & + z(21+nadd)*sinp(2,1))
       dr3a = dr3a
     & + 2.d0*Y31*(z(23+nadd)*cosp(1,1)+z(25+nadd)*sinp(1,1))
     & + 2.d0*Y33*(z(27+nadd)*cosp(3,1)+z(29+nadd)*sinp(3,1))
       dr2aa = dr2aa
     & + z(17+14*nmo+nadd)*Y20+2.d0*Y22*(z(19+14*nmo+nadd)*cosp(2,2) 
     & + z(21+14*nmo+nadd)*sinp(2,2))
       dr3aa = dr3aa
     & + 2.d0*Y31*(z(23+14*nmo+nadd)*cosp(1,2)+z(25+14*nmo+nadd)
     & *sinp(1,2))+2.d0*Y33*(z(27+14*nmo+nadd)*cosp(3,2)
     & +z(29+14*nmo+nadd)*sinp(3,2))
       nadd=nadd+14
      enddo
      s0 = z(2)-ra-raa
      s = z(2)-ra*(1.d0+dr2a-dr3a)-raa*(1.d0+dr2aa+dr3aa)
      aex = s*ainv/2.d0
      aexx = (z(7)-ra-raa)*ainv/2.d0
      if(aex.gt.40)then
       dundr = 0.d0
       un = 0.d0
      else
       eexp = exp(aex)
       un = -v0/(1.d0+eexp**2)
       dundr = v0*ainv/(eexp+1.d0/eexp)**2
      endif
      if(aexx.gt.40)then
       dundrr = 0.d0
       unu = 0.d0
      else
       eeexp = exp(aexx)
       unu = -v0/(1.d0+eeexp**2)
       dundrr = v0*ainv/(eeexp+1.d0/eeexp)**2
      endif
      rdudr(2) = raa*dundr
      rdudr(1) = ra*dundr
      uc = za*zaa*esq/z(2)
      urr = unu+za*zaa*esq/z(7)
      dudrr = dundrr - za*zaa*esq/z(7)**2
      ucc(1,2) = uc*3.d0/5.d0*(rac/z(2))**2
      ucc(1,3) = uc*3.d0/7.d0*(rac/z(2))**3
      ucc(2,2) = uc*3.d0/5.d0*(raac/z(2))**2
      ucc(2,3) = uc*3.d0/7.d0*(raac/z(2))**3
      uc2aa = ucc(2,2)*dr2aa
      uc3aa = ucc(2,3)*dr3aa
      uc2a = ucc(1,2)*dr2a
      uc3a = ucc(1,3)*dr3a
      ur = un+uc+uc2aa+uc3aa+uc2a-uc3a
      dudr=dundr-(uc+3.d0*(uc2aa+uc2a)+4.d0*(uc3aa-uc3a))/z(2)
c
c List of variables:
c trajectory:
c             z(1) =time, z(2) =r, z(3) =pr, z(4) =phi, z(5) =pphi=l,
c no deform.  z(6) =time, z(7) =r, z(8) =pr, z(9) =phi, z(10) =pphi=l,
c
c             z(11)=la_dissip, z(12)=laa_dissip,
c             z(13)=hviba_dissip, z(14)=hvibaa_dissip,
c             z(15)=delphia, z(16)=delphiaa,
c
c modes for i=1(projectile):   lowest n
c   z(17)=alpha20,z(18)=pi20,
c   z(19)=RE(alpha22),z(20)=RE(pi22),z(21)=IM(alpha22),z(22)=IM(pi22),
c   z(23)=RE(alpha31),z(24)=RE(pi31),z(25)=IM(alpha31),z(26)=IM(pi31),
c   z(27)=RE(alpha33),z(28)=RE(pi33),z(29)=IM(alpha33),z(30)=IM(pi33),
c
c modes for i=2(target):       lowest n
c   z(17+nmo*14)=alpha20,z(18+nmo*14)=pi20,
c   z(19+nmo*14)=RE(alpha22),z(20+nmo*14)=RE(pi22),
c   z(21+nmo*14)=IM(alpha22),z(22+nmo*14)=IM(pi22),
c   z(23+nmo*14)=RE(alpha31),z(24+nmo*14)=RE(pi31),
c   z(25+nmo*14)=IM(alpha31),z(26+nmo*14)=IM(pi31),
c   z(27+nmo*14)=RE(alpha33),z(28+nmo*14)=RE(pi33),
c   z(29+nmo*14)=IM(alpha33),z(30+nmo*14)=IM(pi33),
c
c
c      neq1=16+2*nmo*14+1
c
c       z(neq1)=aldis(1,2,lowest n),z(neq1+1)=aldis(1,3,lowest n),
c       z(neq1+2)=hdis(1,2,lowest n),z(neq1+3)=hdis(1,3,lowest n),
c
c       z(neq1+nmo*4)=aldis(2,2,lowest n),z(neq1+nmo*4+1)=aldis(2,3,lowest n),
c       z(neq1+nmo*4+2)=hdis(2,2,lowest n),z(neq1+nmo*4+3)=hdis(2,3,lowest n),
c
c
c      neq2=16+2*nmo*(14+2+2)
c
c     orbital integrals
c       z(neq2+1)=RE(orb(1,2,lowest n,-2),
c       z(neq2+2)=IM(orb(1,2,lowest n,-2),
c       z(neq2+3)=RE(orb(1,2,lowest n,0),
c       z(neq2+4)=IM(orb(1,2,lowest n,0),
c       z(neq2+5)=RE(orb(1,2,lowest n,2),
c       z(neq2+6)=IM(orb(1,2,lowest n,2),
c
c       z(neq2+7)=RE(orb(1,3,lowest n,-3),
c       z(neq2+8)=IM(orb(1,3,lowest n,-3),
c       z(neq2+9)=RE(orb(1,3,lowest n,-1),
c       z(neq2+10)=IM(orb(1,3,lowest n,-1),
c       z(neq2+11)=RE(orb(1,3,lowest n,1),
c       z(neq2+12)=IM(orb(1,3,lowest n,1),
c       z(neq2+13)=RE(orb(1,3,lowest n,3),
c       z(neq2+14)=IM(orb(1,3,lowest n,3),
c
c       z(neq2+14*nmo+1)=RE(orb(2,2,lowest n,-2),
c       z(neq2+14*nmo+2)=IM(orb(2,2,lowest n,-2),
c       z(neq2+14*nmo+3)=RE(orb(2,2,lowest n,0),
c       z(neq2+14*nmo+4)=IM(orb(2,2,lowest n,0),
c       z(neq2+14*nmo+5)=RE(orb(2,2,lowest n,2),
c       z(neq2+14*nmo+6)=IM(orb(2,2,lowest n,2),
c
c       z(neq2+14*nmo+7)=RE(orb(2,3,lowest n,-3),
c       z(neq2+14*nmo+8)=IM(orb(2,3,lowest n,-3),
c       z(neq2+14*nmo+9)=RE(orb(2,3,lowest n,-1),
c       z(neq2+14*nmo+10)=IM(orb(2,3,lowest n,-1),
c       z(neq2+14*nmo+11)=RE(orb(2,3,lowest n,1),
c       z(neq2+14*nmo+12)=IM(orb(2,3,lowest n,1),
c       z(neq2+14*nmo+13)=RE(orb(2,3,lowest n,3),
c       z(neq2+14*nmo+14)=IM(orb(2,3,lowest n,3),
c
      zdot(2) = hc2om/am0*z(3)*zdot(1)
      zdot(3) = (hc2om/am0*z(5)**2/z(2)**3 - dudr)*zdot(1)
      zdot(4) = hc2om/am0*z(5)/z(2)**2*zdot(1)
      zdot(5) = 0.d0
      madd=0
      do i = 1,2
       do n=nmin,nmax
        zdot(5) = zdot(5) +
     & (-(rdudr(i)-ucc(i,2))*(Y22*4.d0*(z(19+madd)*sinp(2,i)
     & -z(21+madd)*cosp(2,i)))
     & -ee(i)*(rdudr(i)-ucc(i,3))*
     & (Y31*2.d0*(z(23+madd)*sinp(1,i)-z(25+madd)*cosp(1,i))
     & +Y33*6.d0*(z(27+madd)*sinp(3,i)-z(29+madd)*cosp(3,i))))*zdot(1)
        madd=madd+14
       enddo
      enddo
      zdot(6) = z(7)/vel
      zdot(7) = hc2om/am0*z(8)*zdot(6)
      zdot(8) = (hc2om/am0*z(10)**2/z(7)**3 - dudrr)*zdot(6)
      zdot(9) = hc2om/am0*z(10)/z(7)**2*zdot(6)
      zdot(10) = 0.d0
      madd=0
      do i=1,2
       if(iptt.ne.0.and.i.ne.iptt)goto 11
       if(i.eq.2) madd=nmo*14
       do n=nmin,nmax
        if(Din(i,2,n).lt.5.d4)then
         zdot(17+madd) = z(18+madd)/Din(i,2,n)*zdot(1)
         zdot(19+madd) = z(20+madd)/Din(i,2,n)*zdot(1)
         zdot(21+madd) = -z(22+madd)/Din(i,2,n)*zdot(1)
        zdot(18+madd) = (-(Cin(i,2,n)+2.d0*gnon*hpot(i,2,n))*z(17+madd)
     & + (rdudr(i)-ucc(i,2))*Y20
     & - gam(i,2,n)*z(18+madd)/Din(i,2,n))*zdot(1)
        zdot(20+madd) = (-(Cin(i,2,n)+2.d0*gnon*hpot(i,2,n))*z(19+madd)
     & + (rdudr(i)-ucc(i,2))*Y22*cosp(2,i)
     & - gam(i,2,n)*z(20+madd)/Din(i,2,n))*zdot(1)
        zdot(22+madd) = ((Cin(i,2,n)+2.d0*gnon*hpot(i,2,n))*z(21+madd)
     & - (rdudr(i)-ucc(i,2))*Y22*sinp(2,i)
     & - gam(i,2,n)*z(22+madd)/Din(i,2,n))*zdot(1)
        endif
        if(Din(i,3,n).lt.5.d4)then
         zdot(23+madd) = z(24+madd)/Din(i,3,n)*zdot(1)
         zdot(25+madd) = -z(26+madd)/Din(i,3,n)*zdot(1)
         zdot(27+madd) = z(28+madd)/Din(i,3,n)*zdot(1)
         zdot(29+madd) = -z(30+madd)/Din(i,3,n)*zdot(1)
        zdot(24+madd) = (-(Cin(i,3,n)+2.d0*gnon*hpot(i,3,n))*z(23+madd)
     & + ee(i)*(rdudr(i)-ucc(i,3))*Y31*cosp(1,i)
     & - gam(i,3,n)*z(24+madd)/Din(i,3,n))*zdot(1)
        zdot(26+madd) = ((Cin(i,3,n)+2.d0*gnon*hpot(i,3,n))*z(25+madd)
     & - ee(i)*(rdudr(i)-ucc(i,3))*Y31*sinp(1,i)
     & - gam(i,3,n)*z(26+madd)/Din(i,3,n))*zdot(1)
        zdot(28+madd) = (-(Cin(i,3,n)+2.d0*gnon*hpot(i,3,n))*z(27+madd)
     & + ee(i)*(rdudr(i)-ucc(i,3))*Y33*cosp(3,i)
     & - gam(i,3,n)*z(28+madd)/Din(i,3,n))*zdot(1)
        zdot(30+madd) = ((Cin(i,3,n)+2.d0*gnon*hpot(i,3,n))*z(29+madd)
     & - ee(i)*(rdudr(i)-ucc(i,3))*Y33*sinp(3,i)
     & - gam(i,3,n)*z(30+madd)/Din(i,3,n))*zdot(1)
        endif
        madd=madd+14
       enddo
11     continue
      enddo
      zdot(11) = zdotcal(1)*zdot(1)
      zdot(12) = zdotcal(2)*zdot(1)
      zdot(13) = zdotch(1)*zdot(1)
      zdot(14) = zdotch(2)*zdot(1)
      zdot(15) = ddelp*z(11)*zdot(1)/ainert
      zdot(16) = ddelp*z(12)*zdot(1)/aainert
      do mm = neq1,neq2
       zdot(mm) = zdot(mm)*zdot(1)
      enddo
c      write(13,'(2(/a,2f13.8))')
c     &' ucc2a, ucc3a =', ucc(1,2),ucc(1,3),
c     &' r, radudr = ',z(2),rdudr(1)    
c      write(13,'(/a,4f9.5))')
c     &' gam2a,gam3a,gam2aa,gam3aa = ',gam(1,2,1),gam(1,3,1)
c     &,gam(2,2,1),gam(2,3,1)
c      write(13,'(/a,4f12.8,/a,2(/a,3f12.8/a,4f12.8)
c     &/a,2(/a,3f12.8/a,4f12.8))')
c     &' r,pr,phi,pphi = ',zdot(2),zdot(3),zdot(4),zdot(5),
c     &' target:',
c     &' Ralf20,Ralf22,Ialf22 = ',zdot(17+nmo*14),zdot(19+nmo*14),
c     &                           zdot(21+nmo*14),
c     &' Ralf31,Ialf31,Ralf33,Ialf33=',zdot(23+nmo*14),zdot(25+nmo*14),
c     &                                zdot(27+nmo*14),zdot(29+nmo*14),
c     &' Rpi20,Rpi22,Ipi22 = ',zdot(18+nmo*14),zdot(20+nmo*14),zdot(22+nmo*14),
c     &' Rpi31,Ipi31,Rpi33,Ipi33=',zdot(24+nmo*14),zdot(26+nmo*14),
c     &                            zdot(28+nmo*14),zdot(30+nmo*14),
c     &' projectile:',
c     &' Ralf20,Ralf22,Ialf22 = ',zdot(17),zdot(19),zdot(21),
c     &' Ralf31,Ialf31,Ralf33,Ialf33=',zdot(23),zdot(25),zdot(27),zdot(29),
c     &' Rpi20,Rpi22,Ipi22 = ',zdot(18),zdot(20),zdot(22),
c     &' Rpi31,Ipi31,Rpi33,Ipi33=',zdot(24),zdot(26),zdot(28),zdot(30)
c      nnn=nnn+1
c      if(nnn.gt.20)stop
c
c  Orbital integrals evaluated if ieta = 1 or iin = 1
      if(ieta.eq.0.and.iin.eq.0)goto 114
      madd=neq2
      do i=1,2
       if(iptt.ne.0.and.i.ne.iptt)goto 111
       if(i.eq.2) madd=neq2+nmo*14 
       rdu(i,2) = ucc(i,2)-rdudr(i)
       rdu(i,3) = ucc(i,3)-rdudr(i)
       do n=nmin,nmax
        zdot(1+madd) = alfz(i,2,n)*Y22*rdu(i,2)
     &*cos(eex(i,2,n)*z(1)+2.d0*z(9))*zdot(1)
        zdot(2+madd) = alfz(i,2,n)*Y22*rdu(i,2)
     &*sin(eex(i,2,n)*z(1)+2.d0*z(9))*zdot(1)
        zdot(3+madd) = alfz(i,2,n)*Y20*rdu(i,2)
     &*cos(eex(i,2,n)*z(1))*zdot(1)
        zdot(4+madd) = alfz(i,2,n)*Y20*rdu(i,2)
     &*sin(eex(i,2,n)*z(1))*zdot(1)
        zdot(5+madd) = alfz(i,2,n)*Y22*rdu(i,2)
     &*cos(eex(i,2,n)*z(1)-2.d0*z(9))*zdot(1)
        zdot(6+madd) = alfz(i,2,n)*Y22*rdu(i,2)
     &*sin(eex(i,2,n)*z(1)-2.d0*z(9))*zdot(1)
        zdot(7+madd) = -(-1)**i*alfz(i,3,n)*Y33*rdu(i,3)
     &*cos(eex(i,3,n)*z(1)+3.d0*z(9))*zdot(1)
        zdot(8+madd) = -(-1)**i*alfz(i,3,n)*Y33*rdu(i,3)
     &*sin(eex(i,3,n)*z(1)+3.d0*z(9))*zdot(1)
        zdot(9+madd) = -(-1)**i*alfz(i,3,n)*Y31*rdu(i,3)
     &*cos(eex(i,3,n)*z(1)+1.d0*z(9))*zdot(1)
        zdot(10+madd) = -(-1)**i*alfz(i,3,n)*Y31*rdu(i,3)
     &*sin(eex(i,3,n)*z(1)+1.d0*z(9))*zdot(1)
        zdot(11+madd) = (-1)**i*alfz(i,3,n)*Y31*rdu(i,3)
     &*cos(eex(i,3,n)*z(1)-1.d0*z(9))*zdot(1)
        zdot(12+madd) = (-1)**i*alfz(i,3,n)*Y31*rdu(i,3)
     &*sin(eex(i,3,n)*z(1)-1.d0*z(9))*zdot(1)
        zdot(13+madd) = (-1)**i*alfz(i,3,n)*Y33*rdu(i,3)
     &*cos(eex(i,3,n)*z(1)-3.d0*z(9))*zdot(1)
        zdot(14+madd) = (-1)**i*alfz(i,3,n)*Y33*rdu(i,3)
     &*sin(eex(i,3,n)*z(1)-3.d0*z(9))*zdot(1)
        madd=madd+14
       enddo
111    continue 
      enddo
114   continue
      if(z(2).lt.rmini)rmini = z(2)
c  The minimal distance during integration is stored in common
cc       if(iin.eq.0) then
cc       write(232,'(2i10,2f10.4)') iin,ieta,s,smini
cc       endif
       if(s.lt.smini) smini = s
c  The minimal surface distance during integration is stored in common
c       hdis = z(13)+z(14)
c       hkina = hkin(1,2,1)+hkin(1,3,1)+hkin(1,2,2)+hkin(1,3,2)
c       hkinaa = hkin(2,2,1)+hkin(2,3,1)+hkin(2,2,2)+hkin(2,3,2)
c       hviba = hpota + hkina
c       hpota = hpot(1,2,1)+hpot(1,3,1)+hpot(1,2,2)+hpot(1,3,2)
c       hpotaa = hpot(2,2,1)+hpot(2,3,1)+hpot(2,2,2)+hpot(2,3,2)
c       hvibaa = hpotaa + hkinaa
c       hvibt = hvibaa + hviba
c       form = dudr +uc/z(2)
c       if(abs(form).gt.accur)eta2 = (hvibt +hdis +ur-un-uc)/form**2
c
c      if(z(2).eq.rmini)then
c
      if(s.eq.smini)then  
c
c  (two possibilities!)
c  also the real and imaginary parts of eta are stored if ieta= 1
       if(ieta.eq.0.or.iin.eq.1)goto 113
cc       write(232,'(a)') ' scrivo orbital integrals'
       phid = zdot(9)/zdot(1)
       madd=neq2
       do i=1,2
        if(iptt.ne.0.and.i.ne.iptt)goto 112
        if(i.eq.2) madd=neq2+nmo*14
        rdur(i,2) = ucc(i,2)-rdudr(i)
        rdur(i,3) = ucc(i,3)-rdudr(i)
        do n=nmin,nmax
         reta(i,2,n,1)= -z(2+madd)- alfz(i,2,n)*Y22*rdur(i,2)
     &*cos(eex(i,2,n)*z(1)+2.d0*z(9))/(eex(i,2,n)+2.d0*phid)
         aeta(i,2,n,1)= +z(1+madd)- alfz(i,2,n)*Y22*rdur(i,2)
     &*sin(eex(i,2,n)*z(1)+2.d0*z(9))/(eex(i,2,n)+2.d0*phid)
         reta(i,2,n,2)= -z(4+madd)- alfz(i,2,n)*Y20*rdur(i,2)
     &*cos(eex(i,2,n)*z(1))/eex(i,2,n)
         aeta(i,2,n,2)= +z(3+madd)- alfz(i,2,n)*Y20*rdur(i,2)
     &*sin(eex(i,2,n)*z(1))/eex(i,2,n)
         reta(i,2,n,3)= -z(6+madd)- alfz(i,2,n)*Y22*rdur(i,2)
     &*cos(eex(i,2,n)*z(1)-2.d0*z(9))/(eex(i,2,n)-2.d0*phid)
         aeta(i,2,n,3)= +z(5+madd)- alfz(i,2,n)*Y22*rdur(i,2)
     &*sin(eex(i,2,n)*z(1)-2.d0*z(9))/(eex(i,2,n)-2.d0*phid)
       reta(i,3,n,1)=-z(8+madd)+(-1)**i*alfz(i,3,n)*Y33*rdur(i,3)
     &*cos(eex(i,3,n)*z(1)+3.d0*z(9))/(eex(i,3,n)+3.d0*phid)
       aeta(i,3,n,1)=+z(7+madd)+(-1)**i*alfz(i,3,n)*Y33*rdur(i,3)
     &*sin(eex(i,3,n)*z(1)+3.d0*z(9))/(eex(i,3,n)+3.d0*phid)
       reta(i,3,n,2)=-z(10+madd)+(-1)**i*alfz(i,3,n)*Y31*rdur(i,3)
     &*cos(eex(i,3,n)*z(1)+1.d0*z(9))/(eex(i,3,n)+phid)
       aeta(i,3,n,2)=+z(9+madd)+(-1)**i*alfz(i,3,n)*Y31*rdur(i,3)
     &*sin(eex(i,3,n)*z(1)+1.d0*z(9))/(eex(i,3,n)+phid)
       reta(i,3,n,3)=-z(12+madd)-(-1)**i*alfz(i,3,n)*Y31*rdur(i,3)
     &*cos(eex(i,3,n)*z(1)-1.d0*z(9))/(eex(i,3,n)-phid)
       aeta(i,3,n,3)=+z(11+madd)-(-1)**i*alfz(i,3,n)*Y31*rdur(i,3)
     &*sin(eex(i,3,n)*z(1)-1.d0*z(9))/(eex(i,3,n)-phid)
       reta(i,3,n,4)=-z(14+madd)-(-1)**i*alfz(i,3,n)*Y33*rdur(i,3)
     &*cos(eex(i,3,n)*z(1)-3.d0*z(9))/(eex(i,3,n)-3.d0*phid)
       aeta(i,3,n,4)=+z(13+madd)-(-1)**i*alfz(i,3,n)*Y33*rdur(i,3)
     &*sin(eex(i,3,n)*z(1)-3.d0*z(9))/(eex(i,3,n)-3.d0*phid)
         madd=madd+14
        enddo
112     continue 
       enddo
      endif
113   continue
      return
      end


