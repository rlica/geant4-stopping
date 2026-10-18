c      Subroutine calculates the probability for capture for given l-value 
c       using the barrier distribution barprob(481) 
c       and the WKB expression for the transmission through the barrier 
c       It also calculates the mean shift in the nominal barrier height
c      If barrier penetrability cannot be calculated it is assumed small
c----------------------------------------------------------------
      subroutine transm_s(al,ecm,rcapt,alcapt,ecmcapt,ppcapt,descat
     &,ierr)
      implicit real*8(a-h,o-z)
      common/za/za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/trans/rturn,rgraz,algraz,homegg,alatr
      common/probbar/delebar(481),barprob(481),ebarstep
c
c  vdr,rturn,rgraz,algraz and alatr are not used
c 
      pptot  = 0.d0
      ppcapt = 0.d0
      descat = 0.d0
c  Mean barrier in radial energy:
      centr = hc2om/2.d0*(aa+aaa)/aa/aaa/rcapt**2
      capbar= ecmcapt+(al**2-alcapt**2)*centr
      cbar = capbar - ecm
      elimtr = - homegg/2.d0*pi*log(accur)
      do n=1,481
       eer = delebar(n)
       if(cbar+eer.gt.elimtr)then
        pwkb = 0.d0
       elseif(cbar+eer.le.-elimtr)then
        pwkb = 1.d0
       else
        pwkb = 1.d0/(1.d0 + exp(2.d0*pi*(cbar+eer)/homegg))
       endif
       pptot = pptot + barprob(n)*ebarstep
       aux_gzp=barprob(n)*pwkb*ebarstep
       ppcapt = ppcapt + barprob(n)*pwkb*ebarstep
       descat = descat + barprob(n)*(1.d0-pwkb)*ebarstep*eer
c       if(abs(al-68.3428722).lt.accur.and.mod(n,5).eq.0)then
c        dp = barprob(n)*pwkb
c        dd = barprob(n)*(1.d0-pwkb)*eer
c        print*,' eer,dp,dd = ',eer,dp,dd
c       endif
      enddo      
      if(ppcapt.lt.(1.d0-30.d0*accur))then
       descat = descat/(1.d0-ppcapt)
      else
       descat = ecm-capbar+.1d0
      endif
c      print*,' pptot = ',pptot
      if(abs(pptot-1.d0).gt.30.d0*accur)then
       write(*,*)' pptot in transm_s differs from 1'
       ierr = 1
       goto 2
      endif
      ierr = 0  
c      print*,' descat = ',descat
2     return
      end
