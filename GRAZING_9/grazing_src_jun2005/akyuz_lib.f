c      AKYUZ_LIB contains subroutines and functions to calculate
c      Akyuz potential, effective potential and derivatives
c      Second part contains subroutines for barrier height etc.
c       This part uses nag_shr/opt. 
c      Version where r0 can be changed and accuracy can be set
c      All subroutines work with accuracy = accur/100
c      Most subroutines depend on common blocks /za/,/vdr/ and /ecm/  
c       while /constant/ and /parm/ are calculated in subroutine 
c       constant_s, which is executed in all programs,
c       except if ierr=1 at entry.
c      All subroutines are in double precision.
c      The file contains: 
c       Functions:  unuc, ueff, dunuc, dueff, d2unuc, d2ueff
c       Subroutines: constant_s, unuc_s, ueff_s, ebar_s, 
c       elbar_s, ecrit_s, lgraz_s, scatangle_s, scatangleaver_s,
c       scatangleder_s, scatangleav_s, rainbow_s, rhoftheta_s, turnp_s
c      At entry ierr=1 indicates that /parm/ need not be calculated,
c       ierr=0 is normal entry. 
c      At exit ierr = 0 if normal, ierr=1 if error
c       ierr = -1,2,3,etc means different abnormal conditions.
c      _________________________________________________________________
c
      subroutine constant_s
c
c      Subroutine calculates constants for use in Akyuz
c      Is called before any use of akyuz_lib programs if ierr=0 at entry
c      Results are stored in common/constants/
c

      implicit real*8(a-h,o-z)
      common/za/za,aa,zaa,aaa
      common/ecm/ecm
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      if(za.eq.0.)then
         print*,'Needs data for constant_s'
         stop
      endif
      pi = 3.1415926535897932384d0
      hbarc = 197.33d0
      c2m = 931.5d0
      esq = 1.4397d0
      hc2om = hbarc**2/c2m
      if(aa.lt.1.d0)then
         aa=1.d0
         write(*,*)' Warning aa.lt.1 in constant_s'
      endif
      aa3 = aa**(1./3.d0)
      if(aaa.lt.1.d0)then
         aaa=1.d0
         write(*,*)' Warning aaa.lt.1 in constant_s'
      endif
      aaa3 = aaa**(1./3.d0)
      gam = .95d0*(1.d0-1.8d0*(aa-2.d0*za)*(aaa-2.d0*zaa)/aa/aaa)
      ainv = 1.17d0*(1.d0+.53d0*(1.d0/aa3 + 1.d0/aaa3))
      ra0 = 1.2d0*aa3 - .09d0
      raa0 = 1.2d0*aaa3 - .09d0
      r00 = ra0 + raa0
      ra = ra0*(1.d0+vdr/r00)
      raa = raa0*(1.d0+vdr/r00)
      r0 = ra + raa
      epslab = ecm*(1.d0/aa+1.d0/aaa)
cnan      zeta = exp(-.4d0/50.d0*(epslab-5.d0))
      zeta=1.
      v0 = 16.d0*pi*gam*ra0*raa0/r00/ainv*zeta
c
      return
      end      
c
c     -----------------------------------------------------------
c
      real*8 function unuc(r,ierr)
c
c      Function calculates akyuz potential
c
      implicit real*8(a-h,o-z)
c
      common/za/za, aa, zaa, aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
c
      if(ierr.ne.1) call constant_s
      if (r.lt.0.d0) then 
         print *,' Error  r.lt.0  in unuc '
         ierr = 1
         unuc = 0.d0
        else
         zz = (r-r0)*ainv
         ierr = 0
         if(zz.gt.70.d0)then
            unuc = 0.d0
           else
            unuc = -v0*exp(-zz)/(1.d0+ exp(-zz))
         endif
      endif
c
      return
      end
c
c     ___________________________________________________
c
      real*8 function ueff(al,r,ierr)
c
c     Function calculates the effective potential with Akyuz 
c     for given l and r.  Uses function unuc
c
      implicit real*8(a-h,o-z)
c
      common/za/za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
c
      if(ierr.ne.1) call constant_s
      if(al.lt.0.d0) then
         print*,'error l.lt.0 in ueff'
         ierr=1
        else
         ierr1 = 1
         ueff = za*zaa*esq/r + unuc(r,ierr1)
     1         +al**2/aa/aaa*(aa+aaa)*hc2om/2.d0/r**2
         ierr = 0
         if (ierr1.eq.1)ierr = 1
      endif
      return
      end                                                             
c
c     _____________________________________________________
c
      real*8 function dunuc(r,ierr)
c
c     Function calculates derivative of the Akyuz potential
c
      implicit real*8(a-h,o-z)
      common/za/za, aa, zaa, aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
c
      if(ierr.ne.1) call constant_s
      if(r.lt.0.d0)then
         print*,' Error r.lt.0 in dunuc '
         ierr = 1
         dunuc = 0.d0
        else
         zz = (r-r0)*ainv
         ierr = 0
         if(zz.gt.70.d0)then
            dunuc = 0.d0
           else
            dunuc = v0/4.d0*ainv/(cosh(zz/2.d0))**2
         endif
      endif
c
      return
      end
c
c     ______________________________________________________
c 
      real*8 function dueff(al,r,ierr)
c
c      Function calculates the derivative of the effective 
c      potential with Akyuz for given l and r
c      Uses function dunuc
c
      implicit real*8(a-h,o-z)
c
      common/za/za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
c
      if(ierr.ne.1) call constant_s
      if(al.lt.0.d0)then
         print*,'error l.lt.0 in ueff'
         ierr = 1
         dueff = 0.d0
        else
         ierr1 = 1
         dueff = -za*zaa*esq/r**2+ dunuc(r,ierr1)
     1        -al**2/aa/aaa*(aa+aaa)*hc2om/r**3
         ierr = 0
         if(ierr1.eq.1) ierr = 1
      endif
      return
      end
c
c     -------------------------------------------------------
c
      real*8 function d2unuc(r,ierr)
c
c      Function calculates 2'derivative of the Akyuz potential
c
      implicit real*8(a-h,o-z)
c
      common/za/za, aa, zaa, aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
c
      if(ierr.ne.1) call constant_s
      if (r.lt.0.d0) then
         print*,'Error r.lt.0 in d2unuc '
         ierr = 1
         d2unuc = 0.d0
        else
         arg = (r-r0)*ainv/2.d0
         ierr = 0
         if(arg.gt.20.d0) then
            d2unuc = 0.d0
           else 
            d2unuc = -v0/4.d0*ainv**2/(cosh(arg))**3*sinh(arg)
         endif
      endif
c
      return
      end
c
c     ------------------------------------------------------
c
      real*8 function d3unuc(r,ierr)
c
c     Function calculates 3'derivative of the Akyuz potential
c
      implicit real*8(a-h,o-z)
c
      common/za/za, aa, zaa, aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
c
      if(ierr.ne.1) call constant_s
      if (r.lt.0.d0) then
         print*,'Error r.lt.0 in d3unuc '
         ierr = 1
         d3unuc = 0.d0
        else
         arg = (r-r0)*ainv/2.d0
         ierr = 0
         if(arg.gt.20.d0)then
            d3unuc = 0.d0
           else
            d3unuc = v0/8.d0*ainv**3*(3.d0*(sinh(arg))**2/
     1               (cosh(arg))**4-1.d0/(cosh(arg))**2)
         endif  
      endif
c
      return
      end
c
c     --------------------------------------
c
      real*8 function d2ueff(al,r,ierr)
c
c      Function calculates the 2'derivative of the effective 
c      potential with Akyuz for given l and r
c      Uses function d2unuc
c
      implicit real*8(a-h,o-z)
c
      common/za/za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
c
      if(ierr.ne.1) call constant_s
      if(al.lt.0d0)then
         print*,'error l.lt.0 in d2ueff'
         ierr = 1
        else
         ierr1 = 1
         d2ueff = 2.d0*za*zaa*esq/r**3+ d2unuc(r,ierr1)
     1           +3.d0*al**2/aa/aaa*(aa+aaa)*hc2om/r**4
         ierr = 0
         if (ierr1.eq.1)ierr = 1
      endif
      return
      end
c
c     ------------------------------------------------
c
      subroutine unuc_s(r,pot,dpot,ierr)
c
c     subroutine  Unuc calculates the akyuz potential (pot)
c     and force (dpot)
c
      implicit real*8(a-h,o-z)
c
      common/za/za, aa, zaa, aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
c
      if(ierr.ne.1) call constant_s
c
      if (r.lt.0.d0) then
         print*,'Error r.lt.0 in unuc'
         ierr = 1
        else
         zz = (r-r0)*ainv
         ierr = 0
         if(zz.gt.70.d0)then
            pot = 0.d0
            dpot = 0.d0
           else
            pot = - v0/(1.d0 + exp(zz))
            dpot = v0/4.d0*ainv/(cosh(zz/2.d0))**2
         endif
      endif
c
      return
      end
c
c     ----------------------------------------------------------
c
      subroutine ueff_s(al,r,effpot,ierr)
c
c      subroutine calculates the effective potential with Akyuz 
c      for given l and r
c      Uses subroutine unuc
c
      implicit real*8(a-h,o-z)
c
      common/za/za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
c
      if(ierr.ne.1) call constant_s
      if(al.lt.0.d0)then
         print*,'error l.lt.0 in ueff'
         ierr = 1
        else
         ierr1 = 1
         call unuc_s(r,pot,dpot,ierr1)
         if (ierr1.eq.1)return
         effpot = za*zaa*esq/r+ pot
     1            +al**2/aa/aaa*(aa+aaa)*hc2om/2.d0/r**2
         ierr = 0
      endif
c
      return
      end
c
c     ----------------------------------------------------------
c
      subroutine unuc_c(r,pot,dpot,ierr)
c
c      Subroutine  Unuc_c calculates the akyuz potential (pot)
c      and force (dpot) for complex r
c      complex r,pot,dpot
c
      implicit real*8(a-h,o-z)
c
      complex*16 r,pot,dpot,arg,comcosh
      common/za/  za, aa, zaa, aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
c
      if(ierr.ne.1) call constant_s
c
      rx = dreal(r)
      if(rx.lt.0.d0) then
         print*,'Error rx.lt.0 in unuc_c'
         ierr = 1
        else
         arg = (r-r0)*ainv/2.d0
         ierr = 0
         if(abs(arg).gt.35.d0)then
            pot = 0.d0
            dpot = 0.d0
           else
            pot = - v0/(1.d0 + exp(arg*2.d0))
            comcosh = (exp(arg)+exp(-arg))/2.d0
            dpot = v0/4.d0*ainv/(comcosh)**2
         endif
      endif
c
      return
      end
c
c     -------------------------------------------------------
c
      subroutine ueff_c(al,r,effpot,ierr)
c
c      Subroutine calculates the effective potential with Akyuz 
c      for given complex l and r
c      Uses subroutine unuc_c
c      complex al, r, effpot
c
      implicit real*8(a-h,o-z)
c
      complex*16 al, r, effpot, pot, dpot
      common/za/ za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
c
      if(ierr.ne.1)call constant_s
c
      alx = dreal(al)
      if(alx.lt.0.d0)then 
         print*,'error lx.lt.0 in ueff'
         ierr = 1
        else
         ierr1 = 1
         call unuc_c(r,pot,dpot,ierr1)
         if(ierr1.eq.1) then
            ierr = 1  
           else
            effpot = za*zaa*esq/r+ pot
     1              +al**2/aa/aaa*(aa+aaa)*hc2om/2.d0/r**2
            ierr = 0
         endif
      endif
c
      return
      end
c
c
c--------------------------------------------------------------
c      The second part of library contains subroutines for 
c      trajectory geometry 
c      Uses akyuz_lib and naglib
c--------------------------------------------------------------
c
c
      subroutine ebar_s(eb,rbar,ierr)
c
c      Subroutine Ebar calculates the barrier energy 
c      for l = 0 and its radius 
c      by akyuz potential. Double precision
c      Uses nag_shr/opt subroutine and akyuz_lib
c      At exit ierr = 2 if no barrier(then rbar = r0 and eb = ueff(0,r0))
c
      implicit real*8(a-h,o-z)
      external force,dueff,dunuc,ueff
      common/za/ za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      if(ierr.eq.1)goto 81
      call constant_s
81    continue
      rmin = r0
83    ff = force(rmin)
      if(ff.gt.0.d0)goto 82
      rmin = rmin+.1d0
      if(rmin.gt.r0+5.d0)goto 300
      goto 83
82    rmax = r0 + 20.d0
      ifail = 1
      if(accur.eq.0.d0)goto 200
      acc = accur/100.d0
      acc1 = acc
      call C05ADF(rmin,rmax,acc1,0.d0,force,rbar,ifail)
      if(ifail.ne.0)then
       goto 83
      endif
      ierr = 1
      al = 0.d0
      eb = ueff(al,rbar,ierr)
      if(ierr.eq.1)goto 400
      goto 400
300   rbar = r0
      ierr = 1
      al = 0.d0
      eb = ueff(al,rbar,ierr)
      if(ierr.eq.1)goto 400
      ierr = 2
      goto 100
200   print*,' specify accuracy in akyuz'
      ierr = 1
100   return
400   ierr = 0
      return
      end
      
      real*8 function force(r)
      implicit real*8(a-h,o-z)
      external dueff,dunuc
      ierr = 1
      al = 0.d0
      force = dueff(al,r,ierr)
      return
      end            
c_______________________________________________________________________
c      Subroutine elbar_s calculates the barrier energy 
c      for given l and its radius 
c      by akyuz potential. Double precision
c      Uses nag_shr/opt subroutine and akyuz_lib
c      At exit ierr = 2 if no barrier(then rbar = r0 and eb = ueff(r0))
c--------------------------------------------------------------------------------

      subroutine elbar_s(al,eb,rbar,ierr)
      implicit real*8(a-h,o-z)
      external forcel,dueff,dunuc,ueff
      common/za/ za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      common/elbarloc/all
      if(ierr.eq.1)goto 81
      call constant_s
81    continue
      all = al
      rmax = r0 + 20.d0
      ifail = 1
      if(accur.le.0.d0)goto 200
      acc = accur/100.d0
      acc1 = acc
      call C05ADF(r0,rmax,acc1,0.d0,forcel,rbar,ifail)
      if(ifail.ne.0) goto 300
      ierr = 1
      eb = ueff(al,rbar,ierr)
      if(ierr.eq.1)goto 100
      goto 400
300   rbar = r0
      ierr = 1
      eb = ueff(al,rbar,ierr)
      if(ierr.eq.1)goto 100
      ierr = 2
      goto 100
200   print*,' specify accuracy in akyuz'
      ierr = 1
100   return
400   ierr = 0
      return
      end
      
      real*8 function forcel(r)
      implicit real*8(a-h,o-z)
      external dueff,dunuc
      common/elbarloc/all
      ierr = 1
      forcel = dueff(all,r,ierr)
      return
      end            
c___________________________________________________________________
c      Ecrit_s subroutine calculates the critical energy, 
c      critical l-value and its radius by Akyuz potential
c      Uses nag_shr/opt 
c      External fun
c      At exit ierr = 2 if no critical energy(then rcrit=r0) 
c------------------------------------------------------

      subroutine ecrit_s(ecr,alc,rcrit,ierr)
      implicit real*8(a-h,o-z)
      external fun
      common/za/ za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      if(ierr.eq.1)goto 81
      call constant_s
81    continue
      rmax = r0 + 10.d0
      ifail = 0
      if(accur.eq.0.)goto 200
      acc = accur/100.d0
      call C05ADF(r0,rmax,acc,0.d0,fun,rcrit,ifail)
      if(ifail.eq.1) goto 300
      ierr1 = 1
      em0 = aa*aaa/(aa+aaa)/hc2om
      alg2 = em0*rcrit**3*(dunuc(rcrit,ierr1)-za*zaa*esq/rcrit**2)
      if(ierr1.eq.1)then
       print*,' error1 in ecrit'
       goto 100
      endif
      if(alg2.lt.0.d0)goto 300
      alc = sqrt(alg2)
      ierr1 = 1
      ecr = unuc(rcrit,ierr1)+za*zaa*esq/rcrit+alg2/em0/2.d0/rcrit**2 
      if(ierr1.eq.1)then
       print*,' error2 in ecrit'
       goto 100
      endif
      ierr = 0
      goto 400
300   ierr = 2
      rcrit = r0
      goto 400
200   print*,' specify accuracy in akyuz'
100   ierr = 1
400   return
      end

      real*8 function fun(r)
      implicit real*8(a-h,o-z)
      external d2unuc,dunuc
      common/za/ za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      ierr1 = 1
      ierr2 = 1
      fun = 3.d0/r*dunuc(r,ierr1)+d2unuc(r,ierr2)-2.d0*za*zaa*esq/r**3
      if(ierr1.eq.1.or.ierr2.eq.1)then
       print*,'Error in fun'
       stop
      endif
      return
      end         
c_________________________________________________________________________
c      subroutine lgraz_s calculates grazing l with Akyuz
c      Uses nag_shr/opt,
c      External ful
c      If at exit
c       ierr=0: ecr=critical energy, alc=critical l, eb=Coulomb barrier,
c                 alg=grazing(orbiting) l, rgz=grazing distance
c       ierr=1 if error
c       ierr=2 if below barrier, then eb=Coulomb barrier, rgz = rbar
c       ierr=3 if ecm above Ecritical, then ecr=critical energy, 
c              alc=critical l, eb=Coulomb barrier, rgz=critical distance
c              alg=angular momentum for turning point at r=r0  
c       ierr=4 if no barrier, then ecr = Potential at r0, rgz = r0
c              alg=angular momentum for turning point at r=r0  
c       Unspecified variables are unchanged
c_____________________________________________________________

      subroutine lgraz_s(ecr,alc,eb,alg,rgz,ierr)
      implicit real*8(a-h,o-z)
      external ful
      common/za/ za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/ecm/ecm
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      if(ierr.eq.1)goto 81
      call constant_s
81    continue
      ierr1 = 1
      call ebar_s(eb,rbar,ierr1)
      if(ierr1.eq.1)goto 400
      if(ecm.lt.eb)then
       rgz = rbar
       goto 200
      endif
      if(ierr1.eq.2)goto 600
      ierr1 = 1
      call ecrit_s(ecr,alc,rcrit,ierr1)
      if(ierr1.eq.1)goto 400
      if(ecm.gt.ecr)goto 500
      rmax = rbar
      ifail = 1
      if(accur.eq.0.)goto 222
      acc = accur/100.d0
      acc1 = acc
cnan      eta = 0.0d0
      eta = 1.24242d-5
      call C05ADF(rcrit,rmax,acc1,eta,ful,rgz,ifail)
      if(ifail.ne.0) goto 300
      ierr1 = 1
      alg2 = aa*aaa/(aa+aaa)*rgz**3/hc2om*(dunuc(rgz,ierr1)
     1      -za*zaa*esq/rgz**2)
      if(ierr1.eq.1)goto 400
      if(alg2.lt.0.d0)goto 301
      alg = sqrt(alg2)
      ierr = 0
      goto 100
600   ierr = 4
      rgz = r0
      ierr2 = 1
      tkin = ecm - unuc(r0,ierr2) - za*zaa*esq/r0
      if(tkin.lt.0.d0)then
       write(*,*)' Error in lgraz_s'
       ierr = 1
       return
      endif
      alg = r0*sqrt(aa*aaa/(aa+aaa)*2.d0/hc2om*tkin)
      ierr1 = 1
      ecr = za*zaa*esq/r0+unuc(r0,ierr1)
      if(ierr1.eq.1)then
       write(*,*)' Error in lgraz_s'
       ierr = 1
       return
      endif
      goto 100
300   print*,' lgrazing could not be found'
      print*,' ifail in C05ADF =',ifail
      goto 400
301   print*,' lgraz^2 = ',alg2,' lgraz put = 0'
      algraz = 0.d0
      ierr=0
      return
200   ierr = 2
      goto 100
500   rgz = rcrit
      ierr2 = 1
      tkin = ecm - unuc(r0,ierr2) - za*zaa*esq/r0
      if(tkin.lt.0.d0)then
       write(*,*)' Error in lgraz_s'
       ierr = 1
       return
      endif
      alg = r0*sqrt(aa*aaa/(aa+aaa)*2.d0/hc2om*tkin)
      ierr = 3
      goto 100
222   print*,' specify accuracy in akyuz'
      goto 400
400   ierr = 1
100   return
      end
      
      real*8 function ful(r)
      implicit real*8(a-h,o-z)
      common/za/ za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/ecm/ecm
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      ierr1 = 1
      ierr2 = 1
cnan      ful = r/2.d0*dunuc(r,ierr1)+unuc(r,ierr2)+za*zaa*esq/r/2.d0-ecm
      ful = -ecm+r*dunuc(r,ierr1)/2.d0+unuc(r,ierr2)+za*zaa*esq/r/2.d0
      if(ierr1.eq.1.or.ierr2.eq.1)goto 1
      return
1     print*,'error in ful'
      stop
      end        
c____________________________________________________________________
c      subroutine calculates outermost turning point (rturn) with Akyuz
c      and maximum (emax) and minimum (emin) in effective potential
c      If below barrier in effective potential (ueff) also
c      inner (rinn) and inside (rins) turning points are calculated
c      at exit: ierr = 4 if no inflection point in ueff
c               ierr = 3 if above barrier in ueff
c               ierr = 2 if no barrier in ueff
c               ierr = 1 if error
c               ierr = 0 if one or "three" turning points
c               ierr = -1 if rturn > r0+15 fm 
c      rinn=rins=0  except if ierr = 0 and three turning points
c      Uses nag C05ADF subroutine 
c_____________________________________________________________

      subroutine turnp_s(al,rins,rinn,rturn,emax,emin,ierr)
      implicit real*8(a-h,o-z)
      external fuk,eforce,e2force
      common/za/ za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/ecm/ ecm
      common/turnploc/ aal
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      if(ierr.eq.1)goto 81
      call constant_s
81    continue
      rinn = 0.d0
      rins = 0.d0
      aal=al
      if(aal.lt.0.d0)then
       write(*,*)' l was lt. 0 in turnp_s. Have put l = 0'
       aal = 0.d0
      endif
      if(ecm.le.0.d0)then
       write(*,*)' Error: ecm was le. 0 in turnp_s'
       ierr = 1
       return
      endif
      eta = za*zaa*esq/sqrt(2.d0/aa/aaa*(aa+aaa)*ecm*hc2om)
      ac = za*zaa*esq/2.d0/ecm
c  Coulomb scattering 
      rturnc = ac*(1.d0+sqrt(1.d0+(aal/eta)**2))
      if(rturnc.gt.r0+15.d0)then
       rturn = rturnc
       ierr = -1
       return
      endif
c  calculate the barrier energy and its radius
c  find first,inflection point
      rma = r0 + 1.d0
      rmi = r0
      ifail = 1
      if(accur.eq.0.d0)goto 222
      acc = accur/100.d0
      acc1 = acc
      call C05ADF(rmi,rma,acc1,0.d0,e2force,rinfl1,ifail)
      if(ifail.ne.0) goto 31
      rma = r0 + 50.d0
      rmi = rinfl1+5.d0*acc
      ifail = 1
      acc1 = acc
      call C05ADF(rmi,rma,acc1,0.d0,e2force,rinfl2,ifail)
      if(ifail.ne.0)then
       print*,'error in turnp_s, only one inflection point'
       goto 32
      end if

c  then max and min
      rm = rinfl2
      ri = rinfl1
      ifail = 1
      acc1 = acc
      call C05ADF(ri,rm,acc1,0.d0,eforce,rg,ifail)
      if(ifail.ne.0) goto 30
      rgraz = rg
      ierr1 = 1
      emax = ueff(aal,rgraz,ierr1)
      if(ierr1.eq.1)goto 251
      ifail = 1
      ri = rinfl1
      acc1 = acc
      call C05ADF(.001d0,ri,acc1,0.d0,eforce,rmin,ifail)
      if(ifail.ne.0) goto 300
      ierr1 = 1
      emin = ueff(aal,rmin,ierr1)
      if(ierr1.eq.1)goto 251
      if(ecm.gt.emax) goto 22
      if(ecm.gt.emin)goto 3
      goto 1
30    ierr = 2
      goto 2
31    continue
32    ierr = 4
      goto 2

c  Then find turning point
3     rm = rmin
      ifail = 1
      acc1 = acc
      call C05ADF(.001d0,rm,acc1,0.d0,fuk,rins,ifail)
      if(ifail.ne.0)goto 250
      ri = rmin
      rm = rgraz
      ifail = 1       
      acc1 = acc
      call C05ADF(ri,rm,acc1,0.d0,fuk,rinn,ifail)
      if(ifail.ne.0)goto 250
      goto 1
22    ierr = 3
2     rmax = 50.d0
      ifail = 1
      acc1 = acc
      call C05ADF(.001d0,rmax,acc1,0.d0,fuk,rturn,ifail)
      if(ifail.ne.0)goto 200
      return      
1     ri = rgraz
      rm = 50.d0
      ifail = 1
      acc1 = acc
      call C05ADF(ri,rm,acc1,0.d0,fuk,rturn,ifail)
      if(ifail.ne.0)goto 200
      ierr = 0
      return
200   print*,'turning point outside 50 fm'
      ierr = 1
      return
250   print*,'Error: could not find inner or inside turning point'
      goto 251
222   print*,' specify accuracy in akyuz'
251   ierr = 1
      return
300   print*,'Error: found no minimum'
      return
      end

      real*8 function fuk(r)
      implicit real*8(a-h,o-z)
      external ueff,unuc
      common/za/ za,aa,zaa,aaa
      common/ecm/ ecm
      common/turnploc/ aal
      al=aal
      ierr1 = 1
      fuk = ecm-ueff(al,r,ierr1)
      if(ierr1.eq.1)then
       print*,'Error in fuk'
       stop
      endif
      return      
      end

      real*8 function eforce(r)
      implicit real*8(a-h,o-z)
      external dueff
      common/za/ za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/turnploc/ aal
      al=aal
      ierr1 = 1
      eforce = dueff(al,r,ierr1)
      if(ierr1.eq.1)then
       print*,'Error in eforce'
       stop
      endif
      return
      end            

      real*8 function e2force(r)
      implicit real*8(a-h,o-z)
      external d2ueff
      common/za/ za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/turnploc/ aal
      al=aal
      ierr1 = 1
      e2force = d2ueff(al,r,ierr1)
      if(ierr1.eq.1)then
       print*,'Error in e2force'
       stop
      endif
      return
      end            
c_____________________________________________________________________
c      Subroutine calculates scattering angle (theta) and derivative
c       (dthet) with respect to impact parameter rho by integrating 
c       complex path in double precision (dthet=-dtheta/drho)
c      If at entry ierr = 1 it uses constants in common
c         otherwise it calculates constants
c      If at entry ierr = 2 it stores the turning point in aux
c         otherwise it stores the Coulomb scattering angle in aux
c      At exit: ierr=2 if turning point inside r0
c      Uses subroutines from Akyuz_lib and nag_shr/opt
c       and external d2phidw
c----------------------------------------------------------------

      subroutine scatangle_s(al,aux,theta,dthet,ierr)
      implicit real*8(a-h,o-z)
      complex*16 rturncmpl,rinit,rfin,arg,cce
      external d2phidw,f_out
      common/za/ za,aa,zaa,aaa
      common/ecm/ ecm
      common/vdr/vdr,accur
      common/scatangleloc/rturnx,alph,rho,aal
      dimension phi(4),dphi(4),array(4,7)
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      if(ierr.eq.2)then
       iaux = 1
      elseif(ierr.eq.1)then
       iaux = 0
       goto 81
      else
       iaux = 0
      endif
      call constant_s
81    continue
      aal=al
      if(aal.lt.0.d0)then
       write(*,*)' l was lt. 0 in scatangle_s. Have put l = 0'
       aal = 0.d0
      endif
      alph = .5d0
      acc = accur/100.d0
      wup = .5d0+sqrt(log(alph/acc))
      winit = -wup
      wfin = wup
      rho = aal/sqrt(aa*aaa/(aa+aaa)*ecm*2.d0/hc2om)
c   Calculate turning point with Akyuz
      ierr1 = 1
      call turnp_s(aal,rins,rinn,rturn,emax,emin,ierr1)
      if(ierr1.eq.1)goto 400
      rturnx = rturn
      if(ierr1.ne.0)goto 71
      qu = sqrt((rturn-rinn)/.65d0)
      if(alph.gt.qu)alph = qu
71    elab=ecm*(1.+aa/aaa)
      eta=esq/sqrt(2.d0*hc2om)*za*zaa*sqrt(aa/elab)
      eps=sqrt(1.d0+(aal/eta)**2)
      thetcoul = 2.d0* asin(1.d0/eps)  
      acoul=esq/2.d0*za*zaa/ecm
      ierr2 = 1
      delthet = unuc(rturn,ierr2)/ecm
      if(ierr2.eq.1)goto 400
      if(abs(delthet).lt.acc)then
c  pure Coulomb scattering
       theta = thetcoul
       dthet = 2.d0*acoul/(acoul**2+rho**2)
       if(iaux.eq.0)then
        aux = thetcoul
       elseif(iaux.eq.1)then
        aux = rturnx
       endif
       ierr = 0
       return
      endif
      rturncmpl = dcmplx(rturnx)
      ffa = alph*exp(-winit**2)
      arg = dcmplx(winit,ffa)
      cce = exp(arg)
      rinit = rturncmpl + ((cce+1.d0/cce)/2.d0-1.d0)*.65d0
      rire = dreal(rinit)
      sqtoe = sqrt(1.d0-2.d0*acoul/rire-(rho/rire)**2)
      phi(1) = -asin(1.d0/eps)+asin((1.d0+rho**2/acoul/rire)/eps)
      phi(2) = 0.d0
      phi(3) = -((rire-acoul)/acoul/eps**2-1.d0)/rire/sqtoe
     &+1.d0/acoul/eps**2
      phi(4) = 0.d0
      do n=1,4
       dphi(n) = dphi(n)
      enddo
      ifail = 0
      win = winit
      wfi = wfin
      if(accur.eq.0.)goto 222
      acc1 = acc
c      call D02BAF(win,wfi,4,phi,acc1,d2phidw,array,ifail)
      call D02BBF(win,wfi,4,phi,acc1,0,d2phidw,f_out,array,ifail)
      if(ifail.ne.0)print*,'Warning from integration routine'
      ffa = alph*exp(-wfin**2)
      arg = cmplx(wfin,ffa)
      cce = exp(arg)
      rfin = rturncmpl + ((cce+1.d0/cce)/2.d0-1.d0)*.65d0
      rfre = dreal(rfin)
      sqtoe = sqrt(1.d0-2.d0*acoul/rfre-(rho/rfre)**2)
      phi(1) = phi(1)-asin(1.d0/eps)+asin((1.d0+rho**2/acoul/rfre)/eps)
      theta = pi - phi(1)
      phi(3) = phi(3)-((rfre-acoul)/acoul/eps**2-1.d0)/rfre/sqtoe
     &+1.d0/acoul/eps**2
      dthet = phi(3)
      if(rturnx.lt.r0)then
       ierr = 2
      else
       ierr = 0
      endif 
      if(iaux.eq.0)then
       aux = thetcoul
      elseif(iaux.eq.1)then
       aux = rturnx
      endif
      return
222   print*,' specify accuracy in akyuz'
400   ierr = 1
      return
      end
c
      subroutine d2phidw(w,phi,dphi)
      implicit real*8(a-h,o-z)
      common/za/ za,aa,zaa,aaa
      common/ecm/ ecm
      common/vdr/vdr,accur
      common/scatangleloc/ rturnx,alph,rho,aal
      dimension phi(4),dphi(4)
      complex*16 arg,cce,drdw,fac,dphicom,d2phicom
      complex*16 rturn,rcpl,cal,velov,toe,ueff
      do n=1,4
       phi(n) = phi(n)
      enddo
      cal = dcmplx(aal)
      rturn = dcmplx(rturnx)
      ffa = alph*exp(-w**2)
      arg = dcmplx(w,ffa)
      cce = exp(arg)
      rcpl = rturn + .65d0*((cce+1.d0/cce)/2.d0-1.d0)
      fac = dcmplx(1.d+0,-2.d+0*w*ffa)
      drdw = .65d0*(cce-1.d0/cce)/2.d0*fac
      ierr1 = 1
      call ueff_c(cal,rcpl,ueff,ierr1)
      if(ierr1.eq.1)then
       print*,' Error in d2phidw'
       stop
      endif
      toe=(ecm-ueff)/ecm
      velov = sqrt(toe)
      vely = dimag(velov)
      if (vely.lt.0.d0) velov=-velov
      dphicom = 1.d0/velov/rcpl**2*drdw*rho
      dphi(1) = dreal(dphicom)
      dphi(2) = dimag(dphicom)
      d2phicom = 1.d0/velov**3/rcpl**2*drdw*(toe+(rho/rcpl)**2)
      dphi(3) = dreal(d2phicom)
      dphi(4) = dimag(d2phicom)
      return
      end
c_________________________________________________________________
c      Subroutine calculates up to second derivative of scattering
c       angle (phi=theta-pi) with respect to angular momentum al 
c       by integrating along a complex path in double precision 
c      Uses subroutines from Akyuz_lib and nag_shr/opt
c      external dnphidw
c      At exit: ierr=2 if turning point inside r0
c----------------------------------------------------------------

      subroutine scatangleder_s(al,ph,dph,d2ph,d3ph,ierr)
      implicit real*8(a-h,o-z)
      external dnphidw2,dnphidang,f_out
      common/za/ za,aa,zaa,aaa
      common/ecm/ ecm
      common/vdr/vdr,accur
      common/scatangleloc/rturnx,aal,rho,rad
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      dimension phi(8),array(8,7)
c
      if(ierr.eq.1)goto 81
      call constant_s
81    continue
      alorho = sqrt(aa*aaa/(aa+aaa)*ecm*2.d0/hc2om)
      aal = al
      if(aal.lt.0.d0)then
       write(*,*)' l was lt. 0 in scatangleder_s. Have put l = 0'
       aal = 0.d0
      endif
      rad1 = 1.d0/ainv/2.d0
      rad = rad1
      acc = accur/100.d0
      wup = 1.d0-log(acc)/3.d0
      winit = -wup
      wfin = wup
      rho = aal/alorho
c   Calculate turning point with Akyuz
      ierr1 = 1
      call turnp_s(aal,rins,rinn,rturn,emax,emin,ierr1)
      if(ierr1.eq.1)goto 2
      if(ierr1.eq.3)then
       rrr = r0-rturn
       if(emax-ecm.lt.v0/10.d0.and.rrr.lt.rad)rad = rrr
      endif
      if(ierr1.eq.0)then
       rrr = (rturn - rinn)/2.d0
       if(rrr.lt.rad)rad = rrr
      endif
      acc1 = acc*(rad1/rad)**2
      rturnx = rturn
      elab=ecm*(1.+aa/aaa)
      eta=esq/sqrt(2.d0*hc2om)*za*zaa*sqrt(aa/elab)
      eps=sqrt(1.d0+(aal/eta)**2)
      thetcoul = 2.d0* asin(1.d0/eps)
      acoul=esq/2.d0*za*zaa/ecm
      ierr2 = 1
      delthet = unuc(rturn,ierr2)/ecm
      if(ierr2.eq.1)goto 2
      if(abs(delthet).lt.acc)then
c  pure Coulomb scattering
       ph = pi-thetcoul
       dph = 2.d0*acoul/(acoul**2+rho**2)/alorho
       d2ph = -4.d0*acoul*rho/(acoul**2+rho**2)**2/alorho**2
       d3ph = 4.d0*acoul*(3.d0*rho**2-acoul**2)/(acoul**2+rho**2)**3
     &/alorho**3
       ierr = 0
       return
      endif
      cce = exp(winit)
      rinit = rturnx + rad +((cce+1.d0/cce)/2.d0-1.d0)/ainv
      stkinoe = sqrt(1.d0-2.d0*acoul/rinit-(rho/rinit)**2)
      phi(1) = -asin(1.d0/eps)+asin((1.d0+rho**2/acoul/rinit)/eps)
      phi(2) = 0.d0
      phi(3) = -((rinit-acoul)/acoul/eps**2-1.d0)/rinit/stkinoe
     &+1.d0/acoul/eps**2
      phi(4) = 0.d0
      phi(5) = 0.d0
      phi(6) = 0.d0
      phi(7) = 0.d0
      phi(8) = 0.d0
      win = winit
      wfi = 0.d0
      ifail = 0
c      call D02BAF(win,wfi,8,phi,acc1,dnphidw2,array,ifail)
      call D02BBF(win,wfi,8,phi,acc1,0,dnphidw2,f_out,array,ifail)
      if(ifail.ne.0)write(*,*)
     &' Warning from 1.integration routine in scatangleder'
      pha = phi(1)
      dpha = phi(3)
      d2pha = phi(5)
      d3pha = phi(7)
      angin = 0.d0
      angfi = 2.d0*pi
      ifail = 0
c      call D02BAF(angin,angfi,8,phi,acc1,dnphidang,array,ifail)
      call D02BBF(angin,angfi,8,phi,acc1,0,dnphidang,f_out,array,ifail)
      if(ifail.ne.0)write(*,'(a/a,i4,4e15.5)')
     &' Warning from 2.integration routine in scatangleder',
     &' ifail, phi(n) = ',ifail,phi(2),phi(4),phi(6),phi(8)
      ph = phi(1)+pha
      dph = (phi(3)+dpha)/alorho
      d2ph = (phi(5)+d2pha)/alorho**2
      d3ph = (phi(7)+d3pha)/alorho**3
      if(rturnx.lt.r0)then
       ierr = 2
      else
       ierr = 0
      endif
      return
2     ierr = 1
      return
      end
c
      subroutine dnphidw2(w,phi,dphi)
      implicit real*8(a-h,o-z)
      common/za/ za,aa,zaa,aaa
      common/ecm/ ecm
      common/vdr/vdr,accur
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      common/scatangleloc/rturnx,aal,rho,rad
      dimension phi(8),dphi(8)
      do n=1,8
       phi(n) = phi(n)
      enddo
      cce = exp(w)
      rcpl = rturnx + rad + ((cce+1.d0/cce)/2.d0-1.d0)/ainv
      drdw = -(cce-1.d0/cce)/2.d0/ainv
      ierr1 = 1
      tkoe=(ecm-ueff(aal,rcpl,ierr1))/ecm
      if(ierr1.eq.1)then
       print*,' Error in dnphidw2'
       stop
      endif
      rhoor2 = (rho/rcpl)**2
      tkoead = tkoe+rhoor2
      velov = sqrt(tkoe)
      if(velov.eq.0.d0)velov = 1.d-10 
      dphi(1) = 1.d0/velov/rcpl**2*drdw*rho
      dphi(2) = 0.d0
      dphi(3) = 1.d0/velov**3/rcpl**2*drdw*tkoead
      dphi(4) = 0.d0
      dphi(5) = 3.d0/velov**5/rcpl**4*drdw*tkoead*rho
      dphi(6) = 0.d0
      dphi(7) = 3.d0/velov**7/rcpl**4*drdw*tkoead*(tkoead+4.d0*rhoor2)
      dphi(8) = 0.d0
      return
      end
c
      subroutine dnphidang(ang,phi,dphi)
      implicit real*8(a-h,o-z)
      complex*16 cce,drdang,dphicom,d2phicom,d3phicom,d4phicom
      complex*16 aiunit,rcpl,alc,velov,tkoe,tkoead,rhoor2,ueff
      common/za/ za,aa,zaa,aaa
      common/ecm/ ecm
      common/vdr/vdr,accur
      common/scatangleloc/rturnx,aal,rho,rad
      dimension phi(8),dphi(8)
      do n=1,8
       phi(n) = phi(n)
      enddo
      aiunit = dcmplx(0.d0,1.d0)
      alc = dcmplx(aal)
      cce = rad*exp(aiunit*ang) 
      rcpl = rturnx + cce
      drdang = aiunit*cce
      ierr1 = 1
      call ueff_c(alc,rcpl,ueff,ierr1)
      if(ierr1.eq.1)then
       print*,' Error in dnphidang'
       stop
      endif
      tkoe=(ecm-ueff)/ecm
      rhoor2 = (rho/rcpl)**2
      tkoead = tkoe+rhoor2
      velov = sqrt(tkoe)
      vely = dimag(velov)
      if (vely.ge.0.d0) velov=-velov
      dphicom = 1.d0/velov/rcpl**2*drdang*rho
      dphi(1) = dreal(dphicom)
      dphi(2) = dimag(dphicom)
      d2phicom = 1.d0/velov**3/rcpl**2*drdang*tkoead
      dphi(3) = dreal(d2phicom)
      dphi(4) = dimag(d2phicom)
      d3phicom = 3.d0/velov**5/rcpl**4*drdang*tkoead*rho
      dphi(5) = dreal(d3phicom)
      dphi(6) = dimag(d3phicom)
      d4phicom=3.d0/velov**7/rcpl**4*drdang*tkoead*(tkoead+4.d0*rhoor2)
      dphi(7) = dreal(d4phicom)
      dphi(8) = dimag(d4phicom)
      return
      end

c_______________________________________________________________
c      Subroutine calculates average scattering angle (thetav)
c       for entrance and exit channels and dtheta/drho (dthetav)
c       (negative inside rainbow). rturnex = turning point in exit channel
c      If at entry theta and dthet.ne.0 it uses these values for 
c       the entrance channel. If theta.eq.0 they are calculated.
c      At exit: ierr = 2 if trajectory inside r0 in entrance or exit channel
c      uses akyuz_lib and nag_shr/opt
c----------------------------------------------------------------
      subroutine scatangleav_s(al,theta,dthet,dela,delz,delecm,
     &delal,thetav,dthetav,rturnex,ierr)
      implicit real*8(a-h,o-z)
      common/za/ za,aa,zaa,aaa
      common/ecm/ ecm
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
c
      ecmin=ecm
      zain=za
      zaain=zaa
      aain=aa
      aaain=aaa
      iwarn = 0
      if(abs(theta).lt.1.d-20)then
       ierr = 2
       al1 = al
       if(al1.lt.0.d0)then
        write(*,*)' l1 was lt. 0 in scatangleav_s. Have put l1 = 0'
        al1 = 0.d0
       endif
       call scatangle_s(al1,rturn1,theta1,dthet1,ierr)
       if(ierr.eq.1)then
        write(*,'(a)')' Error in scatangleav_s for al1'
        goto 3
       else
        if(ierr.eq.2)iwarn = 1
        theta = theta1
        dthet = dthet1
       endif
      endif
      ecm=ecm+delecm
      za=za+delz
      zaa=zaa-delz
      aa=aa+dela
      aaa=aaa-dela
      al2=al+delal
      if(al2.lt.0.d0)then
       write(*,*)' l2 was lt. 0 in turnp_s. Have put l2 = 0'
       al2 = 0.d0
      endif
      ierr = 2
      call scatangle_s(al2,rturnex,theta2,dthet2,ierr)
      if(ierr.eq.1)then
       write(*,'(a)')' Error in scatangleav_s for al2'
       goto 2
      elseif(ierr.eq.2)then
       iwarn = 2
      endif
      dthetav=(dthet + dthet2)/2.d0
      thetav=(theta + theta2)/2.d0
      ierr = 0
      if(iwarn.ne.0)ierr = 2
2     ecm=ecmin
      za=zain
      zaa=zaain
      aa=aain
      aaa=aaain
      call constant_s
3     return
      end
c________________________________________________________________
c      Subroutine calculates average scattering angle for entrance 
c       and exit channels when same distance of closest approach
c       and dtheta/drho (negative inside rainbow)
c      At exit: ierr = 2 if rdist below barrier in exit channel
c               ierr = 2 if rdist inside rgraz or rcrit 
c      uses akyuz_lib and nag_shr/opt
c----------------------------------------------------------------
      subroutine scatangleaver_s(rdist,ideln,idelz,delecm,theta
     &,dthet,dell,ierr)
      implicit real*8(a-h,o-z)
      common/za/ za,aa,zaa,aaa
      common/ecm/ ecm
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
c
      ecmin=ecm
      zain=za
      zaain=zaa
      aain=aa
      aaain=aaa
      iwarn = 0
      deln=dble(ideln)
      delz=dble(idelz)
      if(ierr.eq.0)call constant_s
      ierr = 1
      al1 = alofr(rdist,ierr)
      if(ierr.eq.1)goto 3
      call scatangle_s(al1,thetcoul,theta1,dthet1,ierr)
      if(ierr.eq.1)goto 3
      ierr1 = 0
      ierr2 = 0
      drho1 = (1.d0-(za*zaa*esq/rdist-2.d0*unuc(rdist,ierr1)
     &-rdist*dunuc(rdist,ierr2))/2.d0/ecm)*rdist
      if(ierr1.eq.1.or.ierr2.eq.1)then
       ierr = 1
       goto 3
      endif
      ecm=ecm+delecm
      za=za+delz
      zaa=zaa-delz
      aa=aa+deln+delz
      aaa=aaa-deln-delz
      call constant_s
      ierr = 0
      al2=alofr(rdist,ierr)
      if(ierr.eq.2)then
       print*,' Warning: below barrier in exit channel'
       theta = theta1/2.d0
       dell = al1
       dthet = dthet1/2.d0
       ierr = 2
       goto 2
      endif
      ierr = 1
      call scatangle_s(al2,thetcoul,theta2,dthet2,ierr)
      if(ierr.eq.1)goto 2
      ierr1 = 0
      ierr2 = 0
      drho2 = (1.d0-(za*zaa*esq/rdist-2.d0*unuc(rdist,ierr1)
     &-rdist*dunuc(rdist,ierr2))/2.d0/ecm)*rdist
      if(ierr1.eq.1.or.ierr2.eq.1)then
       ierr = 1
       goto 2
      endif 
      dthet=(dthet1 + dthet2*drho2/drho1)/2.d0
      theta=(theta1 + theta2)/2.d0
      dell = al1-al2
      ierr = 0
2     ecm=ecmin
      za=zain
      zaa=zaain
      aa=aain
      aaa=aaain
      call constant_s
      return
3     return
      end      
c______________________________________________________________________
c      Function calculates l+1/2 for given turning point
c      using akyuz potential
c      At exit ierr=1 if error,ierr=2 if r below barrier,
c      and ierr=3 if r inside rgrazing (or rcritical)
c      At entry ierr = 2 avoids error printout for r below barrier
c-----------------------------------------------------------------

      real*8 function alofr(rturn,ierr)
      implicit real*8(a-h,o-z)
      common/za/ za,aa,zaa,aaa
      common/ecm/ecm
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
c
      if(ierr.eq.1)goto 81
      call constant_s
81    continue
      alofr = 0.d0
      ierr1 = 1
      call lgraz_s(ecr,alc,eb,alg,rgraz,ierr1)
      if(ierr1.eq.1)then
       ierr = 1
       return
      elseif(ierr1.ge.2)then
       goto 5
      endif
      if(rturn.le.rgraz)goto 3
5     ierr2 = 1
      tkin = ecm - unuc(rturn,ierr2) - za*zaa*esq/rturn
      if(ierr2.eq.1)return
      if(tkin.lt.0.d0)goto 2
      alofr = rturn*sqrt(aa*aaa/(aa+aaa)*2.d0/hc2om*tkin)
      ierr = 0
      if(ierr1.eq.3.and.rturn.lt.rgraz)ierr = 3
      return
2     if(ierr.ne.2)write(*,*)' rturn below barrier in alofr'
      ierr = 2
      return
3     ierr = 3
      return
      end
c___________________________________________________________________
c      Subroutine calculates rainbow angles
c      uses akyuz_lib and nag_shr/opt
c      At exit: ierr = 1 if error
c        ierr = 0 if one rainbow (then rainbn=0 and alrn=l_grazing+1/2)
c        ierr = 2  no rainbows (rainbn = rainbc = alrn = alrc = 0)
c        ierr = 3  two rainbows
c        ierr = 4  at least two rainbows (but only one outside r0)
c------------------------------------------------------------

      subroutine rainbow_s(rainbc,rainbn,alrc,alrn,ierr)
      implicit real*8(a-h,o-z)
      common/za/ za,aa,zaa,aaa
      common/ecm/ ecm
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      dimension C(11)
      
c
      if(ierr.eq.1)goto 81
      call constant_s
81    ierr1 = 1
      call lgraz_s(ecr,alc,eb,alg,rgr,ierr1)
      if(ierr1.eq.0)goto 50
      if(ierr1.eq.1)goto 52
      if(ierr1.eq.2)then
       ierr = 2
       return
      elseif(ierr1.eq.3)then
       ierr2 = 1
       call turnp_s(0.d0,rins,rinn,rturn,emax,emin,ierr2)
       if(ierr2.eq.1)goto 2
       if(rturn.lt.r0-.6)then
        rmin0 = r0-.6
       else
        rmin0 = rturn+accur
       endif
       goto 55
      elseif(ierr1.eq.4)then
       ierr2 = 1
       call turnp_s(0.d0,rins,rinn,rturn,emax,emin,ierr2)
       if(ierr2.eq.1)goto 2
       if(rturn.gt.r0)then
        ierr = 2
        return
       elseif(rturn.gt.r0-.6)then
        rmin0 = rturn + accur
       else
        rmin0 = r0-.6
       endif
       goto 55
      endif
50    ierr1 = 1
      rmax = rgr+4.d0
      alpos = alofr(rmax,ierr1)
      alneg = alg+accur
      n = 0
      nmax = 100
      rmin = rgr+.6d0
      ierr1 = 1
      alr = alofr(rmin,ierr1)
      ac = accur*alr
501   ierr1 = 1 
      call scatangleder_s(alr,ph,dph,d2ph,d3ph,ierr1)
      if(ierr1.eq.1)goto 52
      ierr2 = 0
      call zero_s(dph,d2ph,alneg,alpos,alr,n,nmax,ac,ierr2)
      if(ierr2.eq.2)goto 501
      if(ierr2.eq.1)goto 52
      rainbn = 0.d0
      alrc = alr
      ierr2 = 1
      call scatangle_s(alrc,thetcoul,rainb,dthet,ierr2)
      if(ierr2.eq.1)goto 52
      if(abs(dthet).gt.10.d0*accur)write(*,*)
     &' dthet = ',dthet,' rainbow not accurate'
      rainbc =rainb
      ierr = 0
      return
55    alposc = 0.d0
      alpos2 = 0.d0
      rmax = r0+3.d0
559   ierr1= 1
      alm = alofr(rmax,ierr1)
      ierr1 = 1 
      call scatangleder_s(alm,ph,dph,d2ph,d3ph,ierr1)
      if(ierr1.eq.1)goto 52
c  d2ph is negative at large l
      if(d2ph.gt.0)then
       rmax = rmax+1.d0
       alpos2 = alm
       if(dph.gt.0.d0)alposc = alpos2
       goto 559
      elseif(d2ph.le.0.d0)then
       almax = alm
       if(alposc.eq.0.d0)alposc = alm
      endif
      rmin = rmin0
541   ierr1 = 1
      almin = alofr(rmin,ierr1)
      ierr1 = 1 
      call scatangleder_s(almin,ph,dph,d2ph,d3ph,ierr1)
      if(ierr1.eq.1)goto 52
      if(dph.lt.0.d0)then
       alneg = almin
       if(alposn.lt.accur)then
        ierr = 4
        goto 58
       else
        goto 53
       endif
      elseif(dph.ge.0.d0.and.d2ph.gt.0.d0)then
       alposn = almin
       if(d3ph.lt.0.d0)then
        rmin = rmin + .3
        goto  541
       elseif(d3ph.ge.0.d0)then
        rmin = rmin + .3
        if(rmin.gt.r0+.3d0)then
         ierr = 2
         return
        else
         goto 541
        endif
       endif
      elseif(dph.ge.0.d0.and.d2ph.le.0.d0)then
       alneg2 = almin
       if(alpos2.gt.accur)then
        goto 57
       endif
c  else searches zero for d2ph between almin and almax
       alr1 = (almin+almax)/2.d0
       ind = -1
       H = (almax-almin)/20.d0
       nn = 0
551    nn = nn + 1
       ierr1 = 1 
       call scatangleder_s(alr1,ph,dph,d2ph,d3ph,ierr1)
       if(ierr1.eq.1)goto 52
       if(d2ph.gt.0.d0)then
        alneg2 = almin
        alpos2 = alr1
        goto 57
       endif
       ifail = -1
       call C05AVF(alr1,d2ph,H,almin,almax,aly,C,ind,ifail)
       if(ifail.eq.1)then
        write(*,'(a,i2,a)')' ifail = ',ifail,' in C05AVF in rainbow_s'
        ierr = 1
        return
       elseif(ifail.eq.4)then
        ierr = 2
        return
       elseif(ifail.gt.0)then
        write(*,'(a,i2,a)')' ifail = ',ifail,' in C05AVF in rainbow_s'
       endif
       if(ind.ne.0)then
        if(nn.gt.20)then
         ierr = 2
         return
        else
         goto 551
        endif
       elseif(ind.eq.0)then
        print*,' Mystery in rainbow'
        alneg2 = almin
        if(d3ph.gt.0.d0)then
         alpos2 = alr1+H
         alr2 = alr1
         goto 571      
        else
         alpos2 = min(alr1,aly)
         goto 57
        endif
       endif
      endif
57    alr2 = (alneg2 +alpos2)/2.d0
571   n = 0
      nmax = 20
c  find inflection between rainbows (d2ph = 0)
552   ierr1 = 1 
      call scatangleder_s(alr2,ph,dph,d2ph,d3ph,ierr1)
      if(ierr1.eq.1)goto 52
      ierr2 = 0
      ac = accur*alr2
      call zero_s(d2ph,d3ph,alneg2,alpos2,alr2,n,nmax,ac,ierr2)
      if(ierr2.eq.2)goto 552
      if(ierr2.eq.1)goto 52
      if(dph.gt.0.d0)then
        ierr = 2
        return
      endif
c  Now rainbows
      alposn = almin
      alneg = alr2
53    alrn = (alneg +alposn)/2.d0
      n = 0
      nmax = 20
553   ierr1 = 1
      call scatangleder_s(alrn,ph,dph,d2ph,d3ph,ierr1)
      if(ierr1.eq.1)goto 52
      ierr2 = 0
      ac = accur*alrn
      call zero_s(dph,d2ph,alneg,alposn,alrn,n,nmax,ac,ierr2)
      if(ierr2.eq.2)goto 553
      if(ierr2.eq.1)goto 52
      rainbn = pi-ph
      n = 0
      nmax = 20
      ierr = 3
58    alrc = (alneg +alposc)/2.d0
554   ierr1 = 1
      call scatangleder_s(alrc,ph,dph,d2ph,d3ph,ierr1)
      if(ierr1.eq.1)goto 52
      ierr2 = 0
      ac = accur*alrc
      call zero_s(dph,d2ph,alneg,alposc,alrc,n,nmax,ac,ierr2)
      if(ierr2.eq.2)goto 554
      if(ierr2.eq.1)goto 52
      rainbc = pi-ph
      return
52    print*,'error in rainbow_s'
2     ierr = 1
      return
      end
c________________________________________________________________
c      Subroutine calculates turning points and l-values leading to 
c      a given scattering angle >0 with the Akyuz potential
c      Maximum number of steps in zero_s is nmax=20
c      At exit: ierr = 1 if error
c      ierr = 0 if one rainbow (two solutions rin=alin=0)
c      ierr = 2 if below barrier(one solution rint=rin=alint=alin=0)
c      ierr = 3 if above E_crit (three solutions) 
c      ierr = 4 if one rainbow (no solutions with theta>0)
c      ierr = 5 if above E_crit (one solution rout=rint=alout=alint=0)
c      ierr = 6 if above E_crit (one solution rin=rint=alin=alint=0)
c      Uses akyuz_lib, zero_s and nag_shr/opt
c-----------------------------------------------------------------
      subroutine rhoftheta_s(theta,rout,rint,rin,alout,alint,alin,ierr)
      implicit real*8(a-h,o-z)
      common/za/ za,aa,zaa,aaa
      common/ecm/ ecm
      common/vdr/vdr,accur
      common/rhofthetaloc/theeta,alorho
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      nmax = 20
      if(ierr.eq.1)goto 88
      call constant_s
88    continue
      alorho = sqrt(2.d0*aa*aaa/(aa+aaa)/hc2om*ecm)
      rhoc = za*zaa*esq/2.d0/ecm/tan(theta/2.d0)
      theeta = theta
      if(accur.le.0.d0)goto 2
      ierr1 = 1
      call rainbow_s(rainbc,rainbn,alrc,alrn,ierr1)
      rainbcdg = rainbc*180.d0/pi
      rainbndg = rainbn*180.d0/pi
      if(ierr1.eq.1)goto 52
c
      if(ierr1.eq.2)then
c  Below barrier
       ierr = 2
       rhoin = rhoc
       tol = rhoin*accur/100.d0
       rhopos = .5d0*rhoin
       rhoneg = 2.d0*rhoin
       n = 0
10     ierr2 = 0
       call tofrho(rhoin,func,dfunc,ierr2)
       if(ierr2.eq.1)goto 52
       call zero(func,dfunc,rhoneg,rhopos,rhoin,n,nmax,accur,ierr3)
       if(ierr3.eq.0)goto 11
       if(ierr3.eq.1)goto 52
       if(ierr3.eq.2)goto 10
11     rint=0.d0
       rin=0.d0
       alint=0.d0
       alin=0.d0
       alout=rhoin*alorho
       ierr2 = 1
       call turnp_s(alout,rins,rinn,rturn,emax,emin,ierr2)
       if(ierr2.eq.1)goto 52
       rout=rturn
       return
c
      elseif(ierr1.eq.3)then 
c  Above E_crit
       if(theta.lt.rainbn)then
        ierr = 5
c  One solution
        rhoin=rhoc
        tol = rhoin*accur/100.d0
        rhopos = alrc/alorho
        rhoneg = rhoc +10.d0*tol
        n = 0
20      ierr2 = 0
        call tofrho(rhoin,func,dfunc,ierr2)
        if(ierr2.eq.1)goto 52
        call zero(func,dfunc,rhoneg,rhopos,rhoin,n,nmax,accur,ierr3)
        if(ierr3.eq.0)goto 21
        if(ierr3.eq.1)goto 52
        if(ierr3.eq.2)goto 20
21      alout=rhoin*alorho
        ierr2 = 1
        call turnp_s(alout,rins,rinn,rturn,emax,emin,ierr2)
        if(ierr2.eq.1)goto 52
        rout=rturn                
       elseif(theta.gt.rainbc)then
        ierr = 6
c  Also one solution
        rhoin = alrn/2.d0/alorho
        tol = rhoin*accur/100.d0
        rhopos = .5d0/alorho
        rhoneg = alrn/alorho-tol
        n = 0
30      ierr2 = 0
        call tofrho(rhoin,func,dfunc,ierr2)
        if(ierr2.eq.1)goto 52
        call zero(func,dfunc,rhoneg,rhopos,rhoin,n,nmax,accur,ierr3)
        if(ierr3.eq.0)goto 31
        if(ierr3.eq.1)goto 52
        if(ierr3.eq.2)goto 30
31      alin=rhoin*alorho
        alint=0.d0
        alout=0.d0
        ierr2 = 1
        call turnp_s(alin,rins,rinn,rturn,emax,emin,ierr2)
        if(ierr2.eq.1)goto 52
        rin=rturn
       else
c  Three solutions
        ierr = 3 
        rhoin = alrn/2.d0/alorho
        tol = rhoin*accur/100.d0
        rhopos = .5d0/alorho
        rhoneg = alrn/alorho-tol
        n = 0
40      ierr2 = 0
        call tofrho(rhoin,func,dfunc,ierr2)
        if(ierr2.eq.1)goto 52
        call zero(func,dfunc,rhoneg,rhopos,rhoin,n,nmax,accur,ierr3)
        if(ierr3.eq.0)goto 41
        if(ierr3.eq.1)goto 52
        if(ierr3.eq.2)goto 40
41      alin=rhoin*alorho
        ierr2 = 1
        call turnp_s(alin,rins,rinn,rturn,emax,emin,ierr2)
        if(ierr2.eq.1)goto 52
        rin=rturn
c  Second
        rhoin = (alrn+alrc)/2.d0/alorho
        tol = rhoin*accur/100.d0
        rhoneg = alrn/alorho+tol
        rhopos = alrc/alorho-tol
        n = 0
50      ierr2 = 0
        call tofrho(rhoin,func,dfunc,ierr2)
        if(ierr2.eq.1)goto 52
        call zero(func,dfunc,rhoneg,rhopos,rhoin,n,nmax,accur,ierr3)
        if(ierr3.eq.0)goto 51
        if(ierr3.eq.1)goto 52
        if(ierr3.eq.2)goto 50
51      alint=rhoin*alorho
        ierr2 = 1
        call turnp_s(alint,rins,rinn,rturn,emax,emin,ierr2)
        if(ierr2.eq.1)goto 52
        rint=rturn
c  Third
        rhoin = rhoc
        tol = rhoin*accur/100.d0
        rhopos = alrc/alorho+tol
        rhoneg = rhoc + 10.d0*tol
        n = 0
60      ierr2 = 0
        call tofrho(rhoin,func,dfunc,ierr2)
        if(ierr2.eq.1)goto 52
        call zero(func,dfunc,rhoneg,rhopos,rhoin,n,nmax,accur,ierr3)
        if(ierr3.eq.0)goto 61
        if(ierr3.eq.1)goto 52
        if(ierr3.eq.2)goto 60
61      alout=rhoin*alorho
        ierr2 = 1
        call turnp_s(alout,rins,rinn,rturn,emax,emin,ierr2)
        if(ierr2.eq.1)goto 52
        rout=rturn
       endif
c
      elseif(ierr1.eq.0)then
c  One rainbow angle. In this case alrn is lgraz
       ierr = 0
       if(theta.lt.rainbc)then
c  Two solutions. 
        rhoin = (alrn+alrc)/2.d0/alorho
        tol = rhoin*accur/100.d0
        rhoneg = alrn/alorho+tol
        rhopos = alrc/alorho-tol
        n = 0
70      ierr2 = 0
        call tofrho(rhoin,func,dfunc,ierr2)
        if(ierr2.eq.1)goto 52
        call zero(func,dfunc,rhoneg,rhopos,rhoin,n,nmax,accur,ierr3)
        if(ierr3.eq.0)goto 71
        if(ierr3.eq.1)goto 52
        if(ierr3.eq.2)goto 70
71      alint=rhoin*alorho
        ierr2 = 1
        call turnp_s(alint,rins,rinn,rturn,emax,emin,ierr2)
        if(ierr2.eq.1)goto 52
        rint=rturn
c  Second solution
        rhoin=rhoc
        tol = rhoin*accur/100.d0
        rhopos = alrc/alorho+tol
        rhoneg = 40.d0
        n = 0
80      ierr2 = 0
        call tofrho(rhoin,func,dfunc,ierr2)
        if(ierr2.eq.1)goto 52
        call zero(func,dfunc,rhoneg,rhopos,rhoin,n,nmax,accur,ierr3)
        if(ierr3.eq.0)goto 81
        if(ierr3.eq.1)goto 52
        if(ierr3.eq.2)goto 80
81      alout=rhoin*alorho
        ierr2 = 1
        call turnp_s(alout,rins,rinn,rturn,emax,emin,ierr2)
        if(ierr2.eq.1)goto 52
        rout=rturn
       elseif(theta.gt.rainbc)then
c  No solution
        ierr=4
        alin=0.d0
        alint=0.d0
        alout=0.d0
        rin=0.d0
        rint=0.d0
        rout=0.d0
       endif 
      endif
      return
52    print*,'Error in tofrho_s'
      ierr = 1
      return
2     print*,'accur.lt.0 in roftheta_s'
      ierr=1
      return
      end
            
      subroutine tofrho(x,func,dfunc,ierr)
      implicit real*8(a-h,o-z)
      common/rhofthetaloc/theeta,alorho
      ierr1 = 1
      al = alorho*x
      call scatangle_s(al,thetcoul,thetl,dthet,ierr1)
      if(ierr1.eq.1)then
       print*,' Error in tofrho'
       ierr = 1
       return
      endif
      func=thetl-theeta
      dfunc = -dthet
      return
      end

c   Subroutine is described under the name zero_s in math_lib
      subroutine zero(func,dfunc,xneg,xpos,x,n,nmax,accur,ierr)
      implicit real*8(a-h,o-z)
      xin = x
      if(abs(xin-xneg).lt.accur.or.abs(xin-xpos).lt.accur)then
       ierr = 0
       return
      endif
      xnext = x - func/dfunc
      if(xpos.lt.xneg)then
       if(x.gt.xneg.or.x.lt.xpos)then
        ierr = 1
        print*,' Error in zero_s: x outside interval'
        return
       endif 
       if(func.gt.0.)xpos=x
       if(func.lt.0.)xneg=x       
       if(xnext.gt.xpos.and.xnext.lt.xneg)then
        x=xnext
       else
        x = (xneg + xpos)/2.d0
       endif
      elseif(xpos.gt.xneg)then
       if(x.lt.xneg.or.x.gt.xpos)then
        ierr = 1
        print*,' Error in zero_s: x outside interval'
        return 
       endif
       if(func.gt.0.d0)xpos=x
       if(func.lt.0.d0)xneg=x       
       if(xnext.lt.xpos.and.xnext.gt.xneg)then
        x=xnext
       else
        x = (xneg + xpos)/2.d0
       endif
      else
       ierr = 1
       print*,' Error in zero_s: xneg=xpos'
      endif
      ierr = 2
      n=n+1
      if(n.gt.nmax)then
       write(*,'(a,i5,a)')
     &' Error in zero_s: zero not found after ',nmax,' iterations' 
       print*,' Check whether func(xneg)<0 and func(xpos)>0'
       ierr = 1
       return
      endif
      if(abs(xin-xnext).lt.accur)ierr = 0
      return
      end
c
c     ----------------------------
c
      subroutine f_out(xsol,y)
      implicit real*8(a-h,o-z)
      dimension y(200)
      xsol=xsol+100.
c      print *,xsol,(y(i),i=1,4)
      return
      end
c
c    -------------------------------------
c
      subroutine f_out4(xsol,y)
c
c     routine introduced for the change from d02baf to d02bbf
c
      implicit real*8(a-h,o-z)
      dimension y(4)
      xsol=xsol+100.d0
      return
      end
