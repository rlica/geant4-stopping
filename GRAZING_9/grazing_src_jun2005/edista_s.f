c      Subroutine calculates the energy distribution of the two nuclei
c       with finite energy resolution delta. The calculation
c       includes low-lying collective states and giant resonances
c      imm (address for eloss = 0) is calculated if 0 at entry
c      mmmax and mmmin (actual range for eloss and ppe) are calculated
c      At exit ierr = 2 indicates warning about accuracy
c------------------------------------------------------------
      subroutine edista_s(iqq,delta,deleloss,accur,imm,mmmax,mmmin,
     &eloss,ppe,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
      external fare,faim
      dimension ppe(ienergy),eloss(ienergy) 
c      common:main
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/q/qns,qnp,qps,qpp
      common/bhom/bhomn,bhomp
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/giantwid/giantwid(3)
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/eav/eav,eavns,eavnp,eavps,eavpp,eavtr,
     &eav0,eava0,eavaa0
      common/wid/wid2e,wid2ns,wid2np,wid2ps,wid2pp,wid2tr,wid2ea,
     &wid2eaa,wid2e0,wid2ea0,wid2eaa0
      common/edisloc/dela,acc1,f0inq,iq
      character*5 version
      common/vers/version

c
      iprep = ierr
      if(version.eq.'vms  '.or.version.eq.'unix* ')then
       open (unit=98,file='grazing_out.dat',status='old')
      endif
      iq = iqq
      nn = 8
c   doubling stepwidth each nn's step (must be even)
      if(iq.eq.1)then
       eavch = eava0*(1.d0-p0a)
       widech = sqrt(wid2ea+delta**2)
      elseif(iq.eq.2)then
       eavch = eavaa0*(1.d0-p0aa)
       widech = sqrt(wid2eaa+delta**2)
      endif
      f0inq = 0.d0
      do lam = 2,3
       do n = 1,2
        f0inq = f0inq + ff0in(iq,lam,n)
       enddo
      enddo
      dela = delta
      if(abs(1.d0-p0).gt..1d0)then
       ff = abs(1.d0-p0)
      else
       ff = .1d0
      endif
      acc = accur*ff
      acc1 = acc/10.d0
      iwarn = 0
      do mm= 1,ienergy
       ppe(mm)= 0.d0
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
       if(1.d0-p0.lt.accur)then
        psn=factor*exp(-el**2/2.d0/dela**2)
        goto 11
       endif
c       write(*,*)' mm = ',mm,' eloss = ',el,' MeV'
       ierr1 = 0
       accc=acc
25     ierr1 = 0
       call fourier_s(psn,fare,faim,el,accc,ierr1)
c       functions fare and faim defined below
       if(ierr1.eq.1)then
        iwarn = 1
        if(accc.lt.acc*1.01)then
         accc=10.d0*accc
         goto 25
        elseif(1.d0-p0.lt.10.d0*accc)then
         psn=factor*exp(-el**2/2.d0/dela**2)
        endif
       endif
11     ppe(mm) = psn
c       write(*,*)' ppe = ',ppe(mm),'/MeV'
       ifact = 2*(1+mod(imm-mm,2))
       fact = dble(ifact)/3.d0*del
       if(mod(imm-mm,nn).eq.0.or.abs(mm-imm)*2.eq.nn)then
c   double stepwidth every nn steps
        if(iprep.eq.1)goto 3
c        write(*,*)' mm_double = ',mm
        fact = del
        del = del*2.d0
       endif
       addt = fact*ppe(mm)
       adde = fact*ppe(mm)*eloss(mm)
       addw = fact*ppe(mm)*eloss(mm)**2
       sum = sum + addt
       sume = sume + adde
       sumw = sumw + addw
       if(addt.lt.acc)goto 435
      enddo
      mm=1
      write(*,'(a)')
     &' Warning: Attempt to calculate below 1 in edista_s'
      write(98,'(a)')
     &' Warning: Attempt to calculate below 1 in edista_s'
      iwarn = 1
435   mmmin = mm
      do mm=mmmin-1,1,-1
       eloss(mm)=eloss(mm+1)-del
       ppe(mm) = 0.d0
       if(mod(imm-mm,nn).eq.0.or.abs(mm-imm)*2.eq.nn)del = del*2.d0
      enddo
c
c   positive values for eloss
      del = deleloss/2.d0 
      do mm = imm,ienergy
       if(mm.eq.imm)then
        eloss(mm) = 0.d+0
       else        
        eloss(mm) = eloss(mm-1) + del
       endif
       el = eloss(mm)
c       write(*,*)' mm = ',mm,' eloss = ',el,' MeV' 
       accc=acc
26     ierr1 = 0
       call fourier_s(psn,fare,faim,el,accc,ierr1)
c       functions fare and faim defined below
       if(ierr1.eq.1)then
        iwarn = 1
        if(accc.lt.acc*1.01)then
         accc=10.d0*accc
         goto 26
        elseif(1.d0-p0.lt.10.d0*accc)then
         psn=factor*exp(-el**2/2.d0/dela**2)
        endif
       endif
       ppe(mm) = psn 
c       write(*,*)' ppe = ',ppe(mm),'/MeV'
       ifact = 2*(1+mod(mm-imm,2))
       fact = dble(ifact)/3.d0*del
       if(mm.ne.imm.and.mod(imm-mm,nn).eq.0.or.abs(mm-imm)*2.eq.nn)then
c    double stepwith every nn steps
c        write(*,*)' mm_double = ',mm
        fact = del
        del = del*2.d0
       endif
       addt = fact*ppe(mm)
       adde = fact*ppe(mm)*eloss(mm)
       addw = fact*ppe(mm)*eloss(mm)**2
       sum = sum + addt
       sume = sume + adde
       sumw = sumw + addw
       if(abs(1.d0-p0).lt.accur.and.addt.gt.acc/100.d0)goto 434
       if(el.gt.eavch.and.addt.lt.acc)then
        goto 436
       elseif(sum.lt.accur)then
        mmmin = mm
       endif
434    continue
      enddo
      mm=ienergy
      write(*,'(a)')
     &' Warning: Attempt to calculate above ienergy in edista_s'
      iwarn = 1
436   mmmax = mm
      do mm=mmmax+1,ienergy
       eloss(mm)=eloss(mm-1)+del
       ppe(mm) = 0.d0
       if(mm.ne.imm.and.mod(imm-mm,nn).eq.0.or.abs(mm-imm)*2.eq.nn)then
        del = del*2.d0
       endif
      enddo
      if(abs(1.-sum).gt.10.*accur)then
       write(*,*)' Total probability in edista_s inaccurate'
       write(*,'(a,f9.6)')'   integral probability = ',sum
       iwarn = 1
      endif
      if(abs((eavch-sume)/eavch).gt.20.*accur)then
       write(*,*)' Mean energy - eav in edista_s inaccurate'
       write(*,'(2(a,f9.6))')'   mean energy = ',sume,'  eav = ',eavch
       iwarn = 1
      endif
      sig2w = sumw - sume**2
      if(sig2w.gt.0.d0)then
       sigw = sqrt(sig2w)
      else
       sigw = 0.d0
      endif
      if(abs((widech-sigw)/widech).gt.40.*accur)then
       write(*,*)' Spread - wide in edista_s inaccurate'
       write(*,'(2(a,f9.6))')'   spread = ',sigw,' ,   wide = ',widech
       iwarn = 1
      endif
3     ierr = 0
      if(iwarn.eq.1)ierr = 2
      return
      end      
c
      real*8 function fare(x)
      implicit real*8(a-h,o-z)
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/bhom/bhomn,bhomp
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/giantwid/giantwid(3)
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/edisloc/dela,acc1,f0inq,iq
      complex*16 gnns,gnnp,gnps,gnpp,fff,zzz,aiunit,gncmplx,fin
      aiunit = dcmplx(0.d0,1.d0)
      xim = -x*bhomn/2.d0
      xip = -x*bhomp/2.d0
      if(abs(xim).lt.(acc1)**.3333)then
       xxin = xim/2.d0
       xxip = xip/2.d0
       gnns = gncmplx(ejns,xxin,1,acc1,ierr)*aiunit
       gnnp = gncmplx(ejnp,xxin,1,acc1,ierr)*aiunit
       gnps = gncmplx(ejps,xxip,1,acc1,ierr)*aiunit
       gnpp = gncmplx(ejpp,xxip,1,acc1,ierr)*aiunit
       goto 22
      endif
      ierr = 0
      xi=xim
      gnns = gncmplx(ejns,xi,0,acc1,ierr)
      if(ierr.eq.1)then
       print*,'Error in fare: ejns,xim,gnns = ',ejns,xi,gnns
       stop
      endif
      ierr=0
      xi=0.d0
      g0nns = gncmplx(ejns,xi,0,acc1,ierr)
      if(ierr.eq.1)then
       print*,'Error in fare: ejns,xi,gnns = ',ejns,xi,g0nns
       stop
      endif
      gnns=-(gnns-g0nns)/2.d0/xim
      xi=xim
      ierr = 0
      gnnp = gncmplx(ejnp,xi,0,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in fare: ejnp,xim,gnnp = ',ejnp,xi,gnnp
       stop                                   
      endif
      ierr = 0
      xi=0.d0
      g0nnp = gncmplx(ejnp,xi,0,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in fare: ejnp,xi,gnnp = ',ejnp,xi,g0nnp
       stop
      endif
      gnnp=-(gnnp-g0nnp)/2.d0/xim
      xi=xip
      ierr = 0
      gnps = gncmplx(ejps,xi,0,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in fare: ejps,xip,gnps = ',ejps,xi,gnps
       stop
      endif
      ierr = 0
      xi=0.d0
      g0nps = gncmplx(ejps,xi,0,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in fare: ejps,xi,g0nps = ',ejps,xi,g0nps
       stop
      endif
      gnps=-(gnps-g0nps)/2.d0/xip
      xi=xip
      ierr = 0
      gnpp = gncmplx(ejpp,xi,0,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in fare: ejpp,xip,gnpp = ',ejpp,xi,gnpp
       stop
      endif
      ierr = 0
      xi=0.d0
      g0npp = gncmplx(ejpp,xi,0,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in fare: ejpp,xi,g0npp = ',ejpp,xi,g0npp
       stop
      endif
      gnpp=-(gnpp-g0npp)/2.d0/xip
22    fin = dcmplx(0.d0,0.d0)
      do lam = 2,3
       fin = fin+ff0in(iq,lam,1)*exp(aiunit*eex(iq,lam,1)*x)
     &+ff0in(iq,lam,2)*exp(aiunit*eex(iq,lam,2)*x
     &-(giantwid(lam)*x)**2/2.d0)
      enddo
      fff = (c1ns*gnns+c1np*gnnp+c1ps*gnps+c1pp*gnpp)/aiunit
     &+fin-ep0n-ep0p-f0inq-(x*dela)**2/2.d0
      zzz = exp(fff)
      fare = dreal(zzz)
      return
      end

      real*8 function faim(x)
      implicit real*8(a-h,o-z)
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/bhom/bhomn,bhomp
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/giantwid/giantwid(3)
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/edisloc/dela,acc1,f0inq,iq
      complex*16 gnns,gnnp,gnps,gnpp,fff,zzz,aiunit,gncmplx,fin
      aiunit = dcmplx(0.d0,1.d0)
      xim = -x*bhomn/2.d0
      xip = -x*bhomp/2.d0
      if(abs(xim).lt.(acc1)**.3333)then
       xxin = xim/2.d0
       xxip = xip/2.d0
       gnns = gncmplx(ejns,xxin,1,acc1,ierr)*aiunit
       gnnp = gncmplx(ejnp,xxin,1,acc1,ierr)*aiunit
       gnps = gncmplx(ejps,xxip,1,acc1,ierr)*aiunit
       gnpp = gncmplx(ejpp,xxip,1,acc1,ierr)*aiunit
       goto 22
      endif
      ierr = 0
      xi=xim
      gnns = gncmplx(ejns,xi,0,acc1,ierr)
      if(ierr.eq.1)then
       print*,'Error in faim: ejns,xim,gnns = ',ejns,xi,gnns
       stop
      endif
      ierr = 0
      xi=0.d0
      g0nns = gncmplx(ejns,xi,0,acc1,ierr)
      if(ierr.eq.1)then
       print*,'Error in faim: ejns,xi,gnns = ',ejns,xi,g0nns
       stop
      endif
      gnns=-(gnns-g0nns)/2.d0/xim
      xi=xim
      ierr = 0
      gnnp = gncmplx(ejnp,xi,0,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in faim: ejnp,xim,gnnp = ',ejnp,xi,gnnp
       stop                                   
      endif
      ierr = 0
       xi=0.d0
      g0nnp = gncmplx(ejnp,xi,0,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in faim: ejnp,xi,gnnp = ',ejnp,xi,g0nnp
       stop
      endif
      gnnp=-(gnnp-g0nnp)/2.d0/xim
      xi=xip
      ierr = 0
      gnps = gncmplx(ejps,xi,0,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in faim: ejps,xip,gnps = ',ejps,xi,gnps
       stop
      endif
      ierr = 0
      xi=0.d0
      g0nps = gncmplx(ejps,xi,0,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in faim: ejps,xi,g0nps = ',ejps,xi,g0nps
       stop
      endif
      gnps=-(gnps-g0nps)/2.d0/xip
      xi=xip
      ierr = 0
      gnpp = gncmplx(ejpp,xi,0,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in faim: ejpp,xim,gnpp = ',ejpp,xi,gnpp
       stop
      endif
      ierr = 0
      xi=0.d0
      g0npp = gncmplx(ejpp,xi,0,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in faim: ejpp,xi,g0npp = ',ejpp,xi,g0npp
       stop
      endif
      gnpp=-(gnpp-g0npp)/2.d0/xip
22    fin = dcmplx(0.d0,0.d0)
      do lam = 2,3
       fin = fin+ff0in(iq,lam,1)*exp(aiunit*eex(iq,lam,1)*x)
     &+ff0in(iq,lam,2)*exp(aiunit*eex(iq,lam,2)*x
     &-(giantwid(lam)*x)**2/2.d0)
      enddo
      fff = (c1ns*gnns+c1np*gnnp+c1ps*gnps+c1pp*gnpp)/aiunit
     &+fin-ep0n-ep0p-f0inq-(x*dela)**2/2.d0
      zzz = exp(fff)
      faim = dimag(zzz)
      return
      end
