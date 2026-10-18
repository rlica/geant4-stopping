c      Subroutine calculates the approximate energy distribution 
c       of the two nuclei including the low-lying collective states. 
c      eavap and wid2ap are the transfer mean energy and squared widths
c       excluding elastic. The latter probability is p0ap
c      imm is calculated if 0 at input
c      iq indicates the projectile(1), target(2) or both(0) 
c       (for the inclusion of inelastic scattering)
c      At exit ierr = 2 indicates warning about accuracy
c------------------------------------------------------------
      subroutine edistapp_s(iq,delta,eavap,wid2ap,p0ap,deleloss,accur,
     &imm,mmmax,mmmin,eloss,ppapp,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
      external fapre,fapim
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/giantwid/giantwid(3)
      common/inel/p0ina,p0inaa,eavina,eavinaa,wid2ina,wid2inaa,
     &avmina,avminaa,eava00,eavaa00,wid2a00,wid2aa00,inelast
      common/edisloc/dela,p0app,en,pn,acc1,f0inq,iqmin,iqmax
      dimension ppapp(ienergy),eloss(ienergy) 
      character*5 version
      common/vers/version

c
      iprep = ierr
      if(version.eq.'vms  '.or.version.eq.'unix* ')then
       open (unit=98,file='grazing_out.dat',status='old')
      endif
      if(iq.eq.0)then
       iqmin = 1
       iqmax = 2
      else
       iqmin = iq
       iqmax = iq        
      endif
      p0app = p0ap
      acc1 = accur/10.d0
      dela = delta
      nn = 8
c   doubling stepwidth each nn's step (must be even)
      do mm= 1,ienergy
       ppapp(mm)= 0.d0
      enddo
      pn = eavap**2/wid2ap
      if(pn.lt.1.d0)then
       write(*,*)' Warning: pn in edistapp_s less than 1'
       iwarn = 1
      endif
      en = eavap/pn
      f0inq = 0.d0
      do ii = iqmin,iqmax
       do lam = 2,3
        do n = 1,2
         f0inq = f0inq + ff0in(ii,lam,n)
        enddo
       enddo
      enddo
      sum = 0.d0
      sume = 0.d0
      sumw = 0.d0
      factor = 1.d0/sqrt(6.28318530717958647d0)/dela
      if(imm.eq.0)imm = int(sqrt(-8.d0*log(acc1))*dela/deleloss)+1
      if(imm.lt.3)imm = 3
c      write(*,*)' imm = ',imm
c  address for eloss = 0
      eloss(imm) = 0.d0
c  initial stepwidth
      del = deleloss/2.d0
c  negative values for eloss
      do mm = imm-1,1,-1
       eloss(mm) = eloss(mm+1) - del
       el = eloss(mm)
       if(1.d0-p0ap.lt.accur)then
        psn=factor*exp(-el**2/2.d0/dela**2)
        goto 11
       endif
c       write(*,*)' mm = ',mm,' eloss = ',el,' MeV'
       accc = accur
26     ierr1 = 0
       call fourier_s(psn,fapre,fapim,el,accc,ierr1)
c       functions fapre and fapim defined below
       if(ierr1.eq.1)then
        iwarn = 1
        if(accc.lt.accur*1.01)then
         accc=10.d0*accc
         goto 26
        elseif(1.d0-p0ap.lt.10.d0*accur)then
         psn=factor*exp(-el**2/2.d0/dela**2)
        endif
       endif
11     ppapp(mm) = psn
c       write(*,*)' ppapp = ',ppapp(mm),'/MeV'
       ifact = 2*(1+mod(imm-mm,2))
       fact = dble(ifact)/3.d0*del
       if(mod(imm-mm,nn).eq.0.or.abs(mm-imm)*2.eq.nn)then
c   double stepwidth every nn steps
        if(iprep.eq.1)goto 3
c        write(*,*)' mm_double = ',mm
        fact = del
        del = del*2.d0
        else
       endif
       addt = fact*ppapp(mm)
       adde = fact*ppapp(mm)*eloss(mm)
       addw = fact*ppapp(mm)*eloss(mm)**2
       sum = sum + addt
       sume = sume + adde
       sumw = sumw + addw
       if(addt.lt.accur/10.d0)goto 435
      enddo
      mm=1
      write(*,'(a)')
     &' Warning: Attempt to calculate below 1 in edistapp_s'
      iwarn = 1
435   mmmin = mm
      do mm=mmmin-1,1,-1
       eloss(mm)=eloss(mm+1)-del
       ppapp(mm) = 0.d0
       if(mod(imm-mm,nn).eq.0.or.abs(mm-imm)*2.eq.nn)del = del*2.d0
      enddo
c
c   positive values for eloss
      del = deleloss/2.d0 
      do mm = imm,ienergy
       if(mm.eq.imm)then
        eloss(mm) = 0.d0
       else        
        eloss(mm) = eloss(mm-1) + del
       endif
       el = eloss(mm)
c       write(*,*)' mm = ',mm,' eloss = ',el,' MeV'
       accc = accur
27     ierr1 = 0
       call fourier_s(psn,fapre,fapim,el,accc,ierr1)
c       functions fapre and fapim defined below
       if(ierr1.eq.1)then
        iwarn = 1
        if(accc.lt.accur*1.01)then
         accc=10.d0*accc
         goto 27
        elseif(1.d0-p0ap.lt.10.d0*accur)then
         psn=factor*exp(-el**2/2.d0/dela**2)
        endif
       endif
       ppapp(mm) = psn 
c       write(*,*)' ppapp = ',ppapp(mm),'/MeV'
       ifact = 2*(1+mod(mm-imm,2))
       fact = dble(ifact)/3.d0*del
       if(imm.ne.mm.and.mod(imm-mm,nn).eq.0.or.abs(mm-imm)*2.eq.nn)then
c    double stepwith every nn steps
c        write(*,*)' mm_double = ',mm
        fact = del
        del = del*2.d0
       endif
       addt = fact*ppapp(mm)
       adde = fact*ppapp(mm)*eloss(mm)
       addw = fact*ppapp(mm)*eloss(mm)**2
       sum = sum + addt
       sume = sume + adde
       sumw = sumw + addw
       if(addt.lt.accur/10.d0.and.eloss(mm).gt.eavap)goto 436
      enddo
      mm=ienergy
      write(*,'(a)')
     &' Warning: Attempt to calculate above ienergy in edistapp_s'
      iwarn = 1
436   mmmax = mm
      do mm=mmmax+1,ienergy
       eloss(mm)=eloss(mm-1)+del
       ppapp(mm) = 0.d0
       if(imm.ne.mm.and.mod(imm-mm,nn).eq.0.or.abs(mm-imm)*2.eq.nn)then
c    double stepwith every nn steps
        del = del*2.d0
       endif
      enddo
      if(abs(1.d0-sum).gt.10.*accur)then
       write(*,*)' Warning: Total probability in edistapp inaccurate'
       write(*,'(a,f9.6)')'  integral probability = ',sum
       iwarn = 1
      endif
      if(iq.eq.1)then
       eavich=eavina
       wid2ich=wid2ina
      elseif(iq.eq.2)then
       eavich=eavinaa
       wid2ich=wid2inaa
      elseif(iq.eq.0)then
       eavich=eavinaa+eavina
       wid2ich=wid2inaa+wid2ina
      endif
      eavch = eavap*(1.d0-p0ap) + eavich
      if(abs(1.d0-eavch/sume).gt.20.*accur)then
       write(*,*)' Warning: Mean energy - eav in edistapp inaccurate'
       write(*,'(2(a,f9.6))')'   Mean energy = ',sume,' , eav = ',eavch
       iwarn = 1
      endif
      wid2ch = (wid2ap+p0ap*eavap**2)*(1.d0-p0ap) + wid2ich + dela**2
      widch = sqrt(wid2ch)
      sig2w = sumw - sume**2
      if(sig2w.gt.0.d0)then
       sigw = sqrt(sig2w)
      else
       sigw = accur
      endif
      if(abs(1.d0-widch/sigw).gt.40.*accur)then
       write(*,*)' Warning: Spread - wide in edistapp inaccurate'
       write(*,'(2(a,f9.6))')'  Spread = ',sigw,' , wide = ',widch
       iwarn = 1
      endif
      ierr = 0
      if(iwarn.eq.1)ierr = 2
      return
3     ierr = 0
      if(iwarn.eq.1)ierr = 2
      return
      end      
c
      real*8 function fapre(x)
      implicit real*8(a-h,o-z)
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/giantwid/giantwid(3)
      common/edisloc/dela,p0app,en,pn,acc1,f0inq,iqmin,iqmax
      complex*16 fff,zzz,aiunit
      aiunit = dcmplx(0.d0,1.d0)
      zzz = exp(-(x*dela)**2/2.d0)
     &*(p0app+(1-p0app)*(1.d0-aiunit*x*en)**(-pn))
      fff = dcmplx(0.d0,0.d0)
      do ii = iqmin,iqmax
       do lam = 2,3
        if(ff0in(ii,lam,1).gt.acc1)then
         fff = fff+ff0in(ii,lam,1)*exp(aiunit*eex(ii,lam,1)*x)
        endif
        if(ff0in(ii,lam,2).gt.acc1)then
         fff = fff+ff0in(ii,lam,2)*exp(aiunit*eex(ii,lam,2)*x
     &-(giantwid(lam)*x)**2/2.d0)
        endif
       enddo
      enddo
      zzz = zzz*exp(fff-f0inq)
      fapre = dreal(zzz)
      return
      end
c
      real*8 function fapim(x)
      implicit real*8(a-h,o-z)
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/giantwid/giantwid(3)
      common/edisloc/dela,p0app,en,pn,acc1,f0inq,iqmin,iqmax
      complex*16 fff,zzz,aiunit
      aiunit = dcmplx(0.d0,1.d0)
      zzz = exp(-(x*dela)**2/2.d0)
     &*(p0app+(1-p0app)*(1.d0-aiunit*x*en)**(-pn))
      fff = dcmplx(0.d0,0.d0)
      do ii = iqmin,iqmax
       do lam = 2,3
        if(ff0in(ii,lam,1).gt.acc1)then
         fff = fff+ff0in(ii,lam,1)*exp(aiunit*eex(ii,lam,1)*x)
        endif
        if(ff0in(ii,lam,2).gt.acc1)then
         fff = fff+ff0in(ii,lam,2)*exp(aiunit*eex(ii,lam,2)*x
     &-(giantwid(lam)*x)**2/2.d0)
        endif
       enddo
      enddo
      zzz = zzz*exp(fff-f0inq)
      fapim = dimag(zzz)
      return
      end
