c      Subroutine calculates scattering angle (phi=theta-pi) 
c       for a reaction with energy change dele, angular momentum
c       change delal, ejectile mass and charge change dela and delz
c       The change takes place over a region of width delw 
c       centered at w0 from the turning point
c      Uses subroutines from Akyuz_lib, math_lib and nag_shr/opt
c      external xphidw
c      At exit: ierr=1 if scattering angle could not be found
c----------------------------------------------------------------

      subroutine scatangleavx_s(al,ph,dphidrho,rturnex,dele,delal,dela,
     &delz,w0,delw,ierr)
      implicit real*8(a-h,o-z)
      external dxphidw1,dxphidw2,dxphidang,f_out4
      common/za/ za,aa,zaa,aaa
      common/ecm/ ecm
      common/vdr/vdr,accur
      common/scatangleloc/rturnin,aal,rho,rad,dle,dlal,dla,dlz,w00,dlw
     &,iwf
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      dimension phi(4),array(4,7)
c
      dle = dele
      dlal = delal
      dla = dela
      dlz = delz
      ecmin = ecm
      zain = za
      aain = aa
      zaain = zaa
      aaain = aaa
c   Calculate turning point in exit channel
      za = za + delz
      aa = aa + dela
      zaa = zaa -delz
      aaa = aaa - dela
      ecm = ecm + dele
      aaal = al + delal
      ierr1 = 0
      call turnp_s(aaal,rins,rinn,rturnex,emax,emin,ierr1)
c  reinstalling parameters in common blocks (even if error)
      ecm = ecmin
      za = zain
      aa = aain
      zaa = zaain
      aaa = aaain
      call constant_s
      if(ierr1.eq.1)then
       write(*,*)' Error in turnp in scatangleav'
       goto 2
      endif
      alorho = sqrt(aa*aaa/(aa+aaa)*ecm*2.d0/hc2om)
c  Preparation of first integration
      aal = al
      acc = accur
      alga= -log(acc)
      wup = 1.d0+alga/3.d0
c  ie. upper limits for joining with Coulomb trajectory
      winit = -wup
      wfin = wup
      rho = aal/alorho
c    radius in integration around turning point in second integration:
      rad1 = 1.d0/ainv/2.d0
      rad = rad1
c   Calculate turning point in entrance channel
      ierr1 = 1
      call turnp_s(aal,rins,rinn,rturnin,emax,emin,ierr1)
      if(ierr1.eq.1)then
       write(*,*)' error in second turnp in scatangle'
       goto 2
      elseif(ierr1.eq.3)then
       rrr = r0-rturnin
       if(emax-ecm.lt.v0/10.d0.and.rrr.lt.rad)rad = rrr
c   ie. changing this radius if entrance channel above barrier
      elseif(ierr1.eq.0)then
       rrr = (rturnin - rinn)/2.d0
       if(rrr.lt.rad)rad = rrr
c   ie. changing this radius if turning point close to rinn
      endif
      acc1 = acc*(rad1/rad)**2
c calculate Coulomb tail
      elab=ecm*(1.+aa/aaa)
      eta=esq/sqrt(2.d0*hc2om)*za*zaa*sqrt(aa/elab)
      eps=sqrt(1.d0+(al/eta)**2)
      thetcoul = 2.d0* asin(1.d0/eps)
      acoul=esq/2.d0*za*zaa/ecm
c  now integrating entrance trajectory up to rturnin+rad
      cce = exp(winit)
      rinit = rturnin + rad +((cce+1.d0/cce)/2.d0-1.d0)/ainv
      stkinoe = sqrt(1.d0-2.d0*acoul/rinit-(rho/rinit)**2)
      phi(1) = -asin(1.d0/eps)+asin((1.d0+rho**2/acoul/rinit)/eps)
c  ie. the Coulomb tail of phi
      phi(2) = 0.d0
c  is reserved for imaginary part of phi in 2. integration
      phi(3) = -((rinit-acoul)/acoul/eps**2-1.d0)/rinit/stkinoe
     &+1.d0/acoul/eps**2
c  ie the derivative dphi/drho
      phi(4) = 0.d0
c  and the imaginary part
      win = winit
      wfi = 0.d0
      ifail = 0
      accc = acc1
c      call D02BAF(win,wfi,4,phi,accc,dxphidw1,array,ifail)
      call D02BBF(win,wfi,4,phi,accc,0,dxphidw1,f_out4,array,ifail)
      if(ifail.ne.0)write(*,*)
     &' Warning from 1.integration routine in scatangleavx'
c  now we have integrated along the real axis up to rturnin+rad
c   and perform complex integrating around rturnin with radius rad
      angin = 0.d0
      angfi = 2.d0*pi
      ifail = 0
      accc = acc1
c      call D02BAF(angin,angfi,4,phi,accc,dxphidang,array,ifail)
      call D02BBF(angin,angfi,4,phi,accc,0,dxphidang,f_out4,array,ifail)
      if(ifail.ne.0)write(*,'(a/a,i4,4e15.5)')
     &' Warning from 2.integration routine in scatangleavx',
     &' ifail, phi(n) = ',ifail,phi(1),phi(2),phi(3),phi(4)
c  prepare integration in exit channel
      if(rturnex.ge.rturnin+2.d0*rad)then
       arg = 1.d0+((rturnex-rturnin)/2.d0-rad)*ainv
       wturn = log(arg+sqrt(arg**2-1.d0))
c  ie. the value of w where r(w) = rturnex in the exit channel integration
      else
       wturn = 0
      endif
      iwf = 0
      nn=0
33    call ierfc(1,(w0-wturn)/delw,sw)
      w0 = w0 +delw/2.d0*sw
      nn=nn+1
34    w00 = w0
      dlw = delw
      win = 0.d0
      wfi = wfin
c  upper limit changed if wfin too small
      wwf = w0+delw*sqrt(alga)
      if(wfi.lt.wwf)wfi = wwf
      if(wfi.lt.w00+3.d0*dlw)wfi= w00+3.d0*dlw 
      ifail = 1
      accc = acc1
c      call D02BAF(win,wfi,4,phi,accc,dxphidw2,array,ifail)
      call D02BBF(win,wfi,4,phi,accc,0,dxphidw2,f_out4,array,ifail)
      if(iwf.eq.1)then
c  ie. increase wturn if warning during integration
       iwf = 0
       wturn = wturn*1.1d0
       tloss = -dele
       if(nn.le.10)goto 33
       if(ifail.ne.0)write(*,'(a,i2,a)')
     &' ifail = ',ifail,' in 3.integration in scatangleavx'
       write(*,*)'   nn>10        for tloss = ',tloss
      endif
c  Coulomb tail in exit channel
      ecms = ecm + dele
      als = aal+delal
      zas = za+delz
      aas = aa+dela
      zaas = zaa-delz
      aaas = aaa-dela
      rhos = als/sqrt(2.d0*aas*aaas/(aas+aaas)/hc2om*ecms)
      elabs=ecms*(1.d0+aas/aaas)
      etas=esq/sqrt(2.d0*hc2om)*zas*zaas*sqrt(aas/elabs)
      epss=sqrt(1.d0+(als/etas)**2)
      acouls=esq/2.d0*zas*zaas/ecms
      cce = exp(wfin)
      rfin = rturnin + rad +((cce+1.d0/cce)/2.d0-1.d0)/ainv
      argu = (1.d0+rhos**2/acouls/rfin)/epss
      if(argu.gt.1.d0)then
       wfin = wfin*1.2d0
       goto 34
      endif
      stkino = sqrt(1.d0-2.d0*acouls/rfin-(rhos/rfin)**2)
      ph = phi(1)-asin(1.d0/epss)+asin(argu)
      dphidrho = phi(3)-((rfin-acouls)/acouls/epss**2-1.d0)/rfin/stkino
     &+1.d0/acouls/epss**2
      ierr = 0
      return
2     ierr = 1
      return
      end
c
      subroutine dxphidw1(w,phi,dphi)
      implicit real*8(a-h,o-z)
      common/za/ za,aa,zaa,aaa
      common/ecm/ ecm
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      common/scatangleloc/rturnin,aal,rho,rad,dle,dlal,dla,dlz,w00,dlw
     &,iwf
      dimension phi(4),dphi(4)
      do n=1,4
       phi(n) = phi(n)
      enddo
      cce = exp(w)
      rcpl = rturnin + rad + ((cce+1.d0/cce)/2.d0-1.d0)/ainv
      drdw = -(cce-1.d0/cce)/2.d0/ainv
      ierr1 = 1
      tkoe=(ecm-ueff(aal,rcpl,ierr1))/ecm
      if(ierr1.eq.1)then
       print*,' Error in dxphidw1 in scatangleavx_s'
       stop
      endif
      rhoor2 = (rho/rcpl)**2
      tkoead = tkoe+rhoor2
      velov = sqrt(tkoe)
      dphi(1) = 1.d0/velov/rcpl**2*drdw*rho
      dphi(2) = 0.d0
      dphi(3) = 1.d0/velov**3/rcpl**2*drdw*tkoead
      dphi(4) = 0.d0
      return
      end
c

      subroutine dxphidw2(w,phi,dphi)
      implicit real*8(a-h,o-z)
      common/za/ za,aa,zaa,aaa
      common/ecm/ ecm
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      common/scatangleloc/rturnin,aal,rho,rad,dle,dlal,dla,dlz,w00,dlw
     &,iwf
      dimension phi(4),dphi(4)
      do n=1,4
c  to avoid error message
       phi(n) = phi(n)
      enddo
      if(dlw.lt.accur)then
       if(w.gt.w00)then
        sw=1.d0
       elseif(w.le.w00)then
        sw=0.d0
       endif
      else
       call ierfc(0,(w-w00)/dlw,sw)
       sw1 = (1.d0-sw/2.d0)
       call ierfc(0,(w-dlw-w00/10.d0)*2.d0/dlw,sw)
       sw2 = (1.d0-sw/2.d0)
      endif
      ecms = ecm+dle*sw1
      als = aal+dlal*sw2
      zas = za+dlz*sw2
      aas = aa+dla*sw2
      zaas = zaa-dlz*sw2
      aaas = aaa-dla*sw2
      rhos = als/sqrt(2.d0*aas*aaas/(aas+aaas)/hc2om*ecms)
c
      cce = exp(w)
      rcpl = rturnin + rad + ((cce+1.d0/cce)/2.d0-1.d0)/ainv
      drdw = -(cce-1.d0/cce)/2.d0/ainv
      ierr1 = 1
      ueffs = zas*zaas*esq/rcpl+ unuc(rcpl,ierr1)
     1        +als**2/aas/aaas*(aas+aaas)*hc2om/2.d0/rcpl**2
      tkoe=(ecms-ueffs)/ecms
cnan
cnan      print *, ' tkoe = ', tkoe
cnan      print *, ' rcpl = ', rcpl
cnan
      if(tkoe.le.0.d0)then
       if(iwf.eq.0) write(*,*)
     &              ' Warning: diminished accuracy in scatangleavx_s'
       tkoe = .1d0
       iwf = 1
cnan
       dphi(1) = 0.0   !  1.d0*drdw*rhos/velov/rcpl**2
       dphi(2) = 0.d0
       dphi(3) = 0.0   !  1.d0*drdw*tkoead/velov**3/rcpl**2
       dphi(4) = 0.d0
       return
cnan
      endif
      rhoor2 = (rho/rcpl)**2
      tkoead = tkoe+rhoor2
      velov = -sqrt(tkoe)
      dphi(1) = 1.d0*drdw*rhos/velov/rcpl**2
      dphi(2) = 0.d0
      dphi(3) = 1.d0*drdw*tkoead/velov**3/rcpl**2
      dphi(4) = 0.d0
      return
      end
c
      subroutine dxphidang(ang,phi,dphi)
      implicit real*8(a-h,o-z)
      complex*16 cce,drdang,dphicom,d2phicom,d3phicom,d4phicom
      complex*16 aiunit,rcpl,alc,velov,tkoe,tkoead,rhoor2,ueff
      common/za/ za,aa,zaa,aaa
      common/ecm/ ecm
      common/vdr/vdr,accur
      common/scatangleloc/rturnin,aal,rho,rad,dle,dlal,dla,dlz,w00,dlw
     &,iwf
      dimension phi(4),dphi(4)
      do n=1,4
       phi(n) = phi(n)
      enddo
      aiunit = dcmplx(0.d0,1.d0)
      alc = dcmplx(aal)
      cce = rad*exp(aiunit*ang)
      rcpl = rturnin + cce
      drdang = aiunit*cce
      ierr1 = 1
      call ueff_c(alc,rcpl,ueff,ierr1)
      if(ierr1.eq.1)then
       print*,' Error in dnphidang in scatangleavx_s'
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
      return
      end




