c      Subroutine calculates change vdrcap in radius of nuclear
c       potential to make capture happen at given angular momentum
c       alcap and a given energy ecmcap, where capture means either 
c       alcapt = algraz or distance of closest approach for l=alcap 
c       is equal to r0. Correction factor for mass parameter = corr
c      At exit ierr = 0 if alcapt=algraz for this potential
c       and inverse barrier frequency homegg is given. 
c       Also rcap = capture (barrier) radius is given.
c      If at exit ierr = 2 there is no barrier. homegg = 0, rcap = r0
c--------------------------------------------------------------------
      subroutine barrierpos_s(alcap,ecmcap,corr,vdrcap,rcap,homegg,ierr)
      implicit real*8(a-h,o-z)
      common/za/ za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/ecm/ecm
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
c
      ecmstore = ecm
      vdrstore = vdr
c  now ecm and vdr are adjusted to find vdrcap
      ecm = ecmcap
      acc = accur
      acc100 = accur/100.d0
      alcapt = alcap
      vdrg = vdrstore
      r00 = r0
      n = 0
      vdrmi=-100.d0
      vdrma=100.d0
c  here reenter for new iteration:
5     n = n+1
      if(n.gt.30)then
       write(*,*)' n.gt.30 in barrierpos_s'
       goto 2
      endif
      if(vdrma.lt.100.d0.and.vdrmi.gt.-100.d0)vdrg=(vdrma+vdrmi)/2.d0
      vdr = vdrg
c      write(*,*)' vdr = ',vdr
      call constant_s
      if(r0.lt.1.5d0)goto 2
      ierr = 1
      call lgraz_s(ecr,alc,eba,alg,rgr,ierr)
c
c       print *, ' ----- ',n,alc,eba,alg
c
       if(ierr.eq.1)then
       write(*,*)' Error in lgraz in barrierpos'
       goto 2
      elseif(ierr.eq.0)then
c       write(*,'(a,f16.12,a,f16.12,a)')' New lgrazing+1/2= ',alg,
c     & ' New grazing distance=',rgr,' fm'
       cent=(alg+1.d0)*alg*hc2om/2.d0*(aa+aaa)/aa/aaa/rgr**2
       if(abs(alg-alcapt)/(alg+1.d0).gt.acc.and.cent.gt.acc100)then
        if(alg.gt.alcapt.and.vdrg.lt.vdrma)vdrma=vdrg
        if(alg.lt.alcapt.and.vdrg.gt.vdrmi)vdrmi=vdrg
        ierr1 = 0
        vdrg = vdrg+.5d0*(alg+alcapt+1.d0)*(alcapt-alg)*hc2om/2.d0
     &*(aa+aaa)/aa/aaa/rgr**2/dunuc(rgr,ierr1)
        if(vdrg.lt.-r00)vdrg = 1.d0-r00         
        goto 5
       else
        ierr2 = 0
        homeg2 = -(aa+aaa)/aa/aaa*d2ueff(alcapt,rgr,ierr2)*hc2om*corr
        if(ierr2.eq.1)then
         write(*,*)' Error in barrierpos,homegg2'
         goto 2
        endif
        if(homeg2.gt..1d0)then
         homegg = sqrt(homeg2)
        else
         homegg = 0.d0
c         write(*,*)' No penetrability of barrier'
        endif
        ierr = 0
        vdrcap = vdrg
        alcap = alcapt
        rcap = rgr
        goto 8
       endif
      elseif(ierr.eq.2)then
       ierr2 = 0
       al = alcapt
       call turnp_s(al,rins,rinn,rturn,emax,emin,ierr2)
       if(ierr2.eq.2)then
        emax = eba+alcapt**2*hc2om/2.d0*(aa+aaa)/aa/aaa/rturn**2
       endif
c       write(*,'(a,f16.12,a)')
c     &' Energy below effective barrier =',emax,' MeV'
       if(abs(emax-ecm).lt.acc)then
        ierr2 = 0
        homeg2 = -(aa+aaa)/aa/aaa*d2ueff(alcapt,rgr,ierr2)*hc2om*corr
        if(ierr2.eq.1)then
         write(*,*)' Error in barrierpos,homeg2'
         goto 2
        endif
        if(homeg2.gt..1d0)then
         homegg = sqrt(homeg2)
        else
         homegg = 0.d0
c         write(*,*)' No penetrability of barrier'
        endif
        ierr = 0
        vdrcap = vdrg
        alcap = alcapt
        rcap = rgr
        goto 8
       else
        if(ecm.gt.emax.and.vdrg.lt.vdrma)vdrma=vdrg
        if(ecm.lt.emax.and.vdrg.gt.vdrmi)vdrmi=vdrg
        ierr1 = 1
        vdrg = vdrg + .5d0*(emax-ecm)/dunuc(rgr,ierr1)
        if(vdrg.lt.-r00)vdrg = .1d0-r00 
        goto 5
       endif
      elseif(ierr.eq.3.or.ierr.eq.4)then
c       if(ierr.eq.3)write(*,'(a,f8.4,a,/a,f8.4,a)')
c     &' Energy above critical energy = ',ecr,' MeV'
c     & ,' New critical radius =',rgr,' fm.'
c       if(ierr.eq.4)write(*,'(a,f10.5)')
c     &' No barrier even for l=0, U(r0)= ',ecr
       ierr2 = 1
       al = alcapt
       call turnp_s(al,rins,rinn,rturn,emax,emin,ierr2)
c       write(*,*)' rturn = ',rturn,' r0 = ',r0
       if(ierr.eq.4)then
c  (only if no barrier even for l = 0)
        er0 = ecr + alcapt**2*hc2om/2.d0*(aa+aaa)/aa/aaa/r0**2
c        write(*,'(a,f8.4,a)')
c     &' Energy below potential at r0 =',er0,' MeV'
        if(abs(rturn-r0)/r0.gt.acc)then
         if(ecm.gt.er0.and.vdrg.lt.vdrma)vdrma=vdrg
         if(ecm.lt.er0.and.vdrg.gt.vdrmi)vdrmi=vdrg
         ierr2 = 1
         vdrg = vdrg+(er0-ecm)/dunuc(r0,ierr2)
         if(vdrg.lt.-r00)vdrg = .1d0-r00 
         goto 5
        else
         ierr = 2
         vdrcap = vdrg
         homegg = 0.d0
         rcap = r0
         goto 8
        endif
       endif
       ierr1 = 0
       alr0 = alofr(r0,ierr1)
       if(ierr1.eq.1)then
        write(*,*)' Error in alofr in barrierpos'
        goto 2
       endif
       if(abs(alr0-alcapt)/(alcapt+1.d0).gt.acc)then
        if(r0.gt.rturn.and.vdrg.lt.vdrma)vdrma=vdrg
        if(r0.lt.rturn.and.vdrg.gt.vdrmi)vdrmi=vdrg
        vdrg = vdrg + (rturn-r0)/4.d0
        if(vdrg.lt.-r00)vdrg = .1d0-r00 
        goto 5
       else
        ierr = 2
        vdrcap = vdrg
        homegg = 0.d0
        rcap = r0
        goto 8
       endif
      endif
2     ierr = 1
8     vdrcap = vdrcap-vdrstore
      ecm = ecmstore
      vdr = vdrstore
      call constant_s
      return
      end
     
 
