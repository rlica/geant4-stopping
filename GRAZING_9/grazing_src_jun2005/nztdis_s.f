c      Subroutine calculates the total kinetic energy distribution
c       for a given partition with finite energy resolution delta. 
c       The calculation includes low-lying collective states and 
c       giant resonances. Initial step in energy = deleloss
c      imm (address for eloss = 0) is calculated if 0 at entry
c      mmmax and mmmin (actual range for eloss and ppe) are calculated
c      At exit ierr = 2 indicates warning about accuracy
c------------------------------------------------------------
      subroutine nztdis_s(delta,deleloss,accur,imm,mmmax,mmmin
     &,isn,isz,eloss,ppe,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
      external fnzre,fnzim
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/q/qns,qnp,qps,qpp
      common/bhom/bhomn,bhomp
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/giantwid/giantwid(3)
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/inel/p0ina,p0inaa,eavina,eavinaa,wid2ina,wid2inaa,
     &avmina,avminaa,eava00,eavaa00,wid2a00,wid2aa00,inelast
      common/tot/in00,ip00,i0nmin,i0nmax,i0zmin,i0zmax
      common/nzdis/ p0nn(irangen),p0zz(irangep),
     &eavdist(irangen,irangep),eavdist00,avmadist(irangen,irangep),
     &avmadist00,avmdist(irangen,irangep),avmdist00,
     &wid2dist(irangen,irangep),wid2dist00,wid2adist(irangen,irangep)
      common/nzedis/dela,acc1,f0inq,ppnz,in,ip
      dimension ppe(ienergy),eloss(ienergy) 
      character*5 version
      common/vers/version

c
      iprep = ierr
      if(version.eq.'vms  '.or.version.eq.'unix* ')then
       open (unit=98,file='grazing_out.dat',status='old')
      endif
      in = isn
      ip = isz
      ppnz = p0nn(in+in00)*p0zz(ip+ip00)
      dela = delta
      nn = 8
c   doubling stepwidth each nn's step (must be even)
      eavch = eavdist(in+in00,ip+ip00)+eavina+eavinaa
      widech = sqrt(wid2dist(in+in00,ip+ip00)+wid2ina+wid2inaa+dela**2)
      f0inq = 0.d0
      do iq=1,2
       do lam = 2,3
        do n = 1,2
         f0inq = f0inq + ff0in(iq,lam,n)
        enddo
       enddo
      enddo
      acc1 = accur/100.d0
      iwarn = 0
      do mm= 1,ienergy
       ppe(mm)= 0.d0
      enddo
      sum = 0.d0
      sume = 0.d0
      sumw = 0.d0
      factor = 1.d0/sqrt(6.28318530717958647d0)/(dela+accur**2)
      if(imm.eq.0)imm = int(sqrt(-8.d0*log(acc1))*dela/deleloss)+1
      if(imm.lt.3)imm = 3
c  address for eloss = 0
      eloss(imm) = 0.d0
c  initial stepwidth
      del = deleloss/2.d0
c  negative values for eloss
      do mm = imm-1,1,-1
       eloss(mm) = eloss(mm+1) - del
       el = eloss(mm)
c       write(*,*)' mm = ',mm,' eloss = ',el,' MeV'
       if(isn.eq.0.and.isz.eq.0.and.1.d0-p0.lt.accur/10.d0)then
        psn=factor*exp(-el**2/2.d0/dela**2)
        goto 11
       endif
       accc=accur/5.d0
25     ierr1 = 0
       call fourier_s(psn,fnzre,fnzim,el,accc,ierr1)
c       functions fnzre and fnzim defined below
       if(ierr1.ne.0)then
        iwarn = 1
        if(accc.lt.accur)then
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
       if(addt.lt.acc1)goto 435
      enddo
      mm=1
      write(*,'(a)')
     &' Warning: Attempt to calculate below 1 in nztdis_s'
      write(98,'(a)')
     &' Warning: Attempt to calculate below 1 in nztdis_s'
      iwarn = 1
435   mmmin = mm
      do mm=mmmin-1,1,-1
       eloss(mm)=eloss(mm+1)-del
       ppe(mm) = 0.d0
       if(mod(imm-mm,nn).eq.0.or.abs(mm-imm)*2.eq.nn)del = del*2.d0
      enddo
c
c   positive values for eloss
      if(iprep.eq.1)goto 3
      del = deleloss/2.d0 
      do mm = imm,ienergy
       if(mm.eq.imm)then
        eloss(mm) = 0.d+0
       else        
        eloss(mm) = eloss(mm-1) + del
       endif
       el = eloss(mm)
c       write(*,*)' mm = ',mm,' eloss = ',el,' MeV' 
       if(isn.eq.0.and.isz.eq.0.and.1.d0-p0.lt.accur/10.d0)then
        psn=factor*exp(-el**2/2.d0/dela**2)
        goto 12
       endif
       accc=accur/5.d0
26     ierr1 = 0
       call fourier_s(psn,fnzre,fnzim,el,accc,ierr1)
c       functions fnzre and fnzim defined below
       if(ierr1.ne.0)then
        iwarn = 1
        if(accc.lt.accur)then
         accc=10.d0*accc
         goto 26
        elseif(1.d0-p0.lt.10.d0*accc)then
         psn=factor*exp(-el**2/2.d0/dela**2)
        endif
       endif
12     ppe(mm) = psn 
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
       if(abs(1.d0-p0).lt.accur.and.addt.gt.acc1/10.d0)goto 434
       if(el.gt.eavch.and.addt.lt.acc1)then
        goto 436
       elseif(sum.lt.accur)then
        mmmin = mm
       endif
434    continue
      enddo
      mm=ienergy
      write(*,'(a)')
     &' Warning: Attempt to calculate above ienergy in nztdis_s'
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
       write(*,*)' Total probability in nztdis_s inaccurate'
       write(*,'(a,f9.6)')'   integral probability = ',sum
       iwarn = 1
      endif
      if(abs((eavch-sume)/eavch).gt.20.*accur)then
       write(*,*)' Mean energy - eav in nztdis_s inaccurate'
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
       write(*,*)' Spread - wide in nztdis_s inaccurate'
       write(*,'(2(a,f9.6))')'   spread = ',sigw,' ,   wide = ',widech
       iwarn = 1
      endif
3     ierr = 0
      if(iwarn.eq.1)ierr = 2
      return
      end      
c
      real*8 function fnzre(x)
      implicit real*8(a-h,o-z)
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/bhom/bhomn,bhomp
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/giantwid/giantwid(3)
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/nzedis/dela,acc1,f0inq,ppnz,in,ip
      complex*16 gnns,gnnp,gnps,gnpp,fff,aiunit,gncmplx,fin
     &,besn,besp
      aiunit = dcmplx(0.d0,1.d0)
      xim = -x*bhomn/2.d0
      xip = -x*bhomp/2.d0
      ierr = 0
      xi=xim
      gnns = gncmplx(ejns,xi,1,acc1,ierr)
      if(ierr.eq.1)then
       print*,'Error in fnzre: ejns,xim,gnns = ',ejns,xi,gnns
       stop
      endif
      gnns=c1ns*gnns
      xi=xim
      ierr = 0
      gnnp = gncmplx(ejnp,xi,1,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in fnzre: ejnp,xim,gnnp = ',ejnp,xi,gnnp
       stop                                   
      endif
      gnnp=c1np*gnnp
      xi=xip
      ierr = 0
      gnps = gncmplx(ejps,xi,1,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in fnzre: ejps,xip,gnps = ',ejps,xi,gnps
       stop
      endif
      gnps=c1ps*gnps
      xi=xip
      ierr = 0
      gnpp = gncmplx(ejpp,xi,1,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in fnzre: ejpp,xip,gnpp = ',ejpp,xi,gnpp
       stop
      endif
      gnpp=c1pp*gnpp
      ierr=0
      call besseldis_cs(gnnp,gnns,in,besn,ierr)
      if(ierr.eq.1)then
       print*,' Error in fnzre: gnns,gnnp,in,besn=',gnns,gnnp,in,besn
       stop
      endif
      ierr=0
      call besseldis_cs(gnpp,gnps,ip,besp,ierr)
      if(ierr.eq.1)then
       print*,' Error in fnzre: gnps,gnpp,ip,besp=',gnps,gnpp,ip,besp
       stop
      endif
      fin = dcmplx(0.d0,0.d0)
      do iq=1,2
       do lam = 2,3
        fin = fin+ff0in(iq,lam,1)*exp(aiunit*eex(iq,lam,1)*x)
     &+ff0in(iq,lam,2)*exp(aiunit*eex(iq,lam,2)*x
     &-(giantwid(lam)*x)**2/2.d0)
       enddo
      enddo
      fff = exp(gnns+gnnp+gnps+gnpp-ep0n-ep0p)*besn*besp
     &*exp(fin-f0inq-(x*dela)**2/2.d0)
      fnzre = dreal(fff)/ppnz
      return
      end
c
      real*8 function fnzim(x)
      implicit real*8(a-h,o-z)
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/bhom/bhomn,bhomp
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/giantwid/giantwid(3)
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/nzedis/dela,acc1,f0inq,ppnz,in,ip
      complex*16 gnns,gnnp,gnps,gnpp,fff,aiunit,gncmplx,fin
     &,besn,besp
      aiunit = dcmplx(0.d0,1.d0)
      xim = -x*bhomn/2.d0
      xip = -x*bhomp/2.d0
      ierr = 0
      xi=xim
      gnns = gncmplx(ejns,xi,1,acc1,ierr)
      if(ierr.eq.1)then
       print*,'Error in fnzre: ejns,xim,gnns = ',ejns,xi,gnns
       stop
      endif
      gnns=c1ns*gnns
      xi=xim
      ierr = 0
      gnnp = gncmplx(ejnp,xi,1,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in fnzre: ejnp,xim,gnnp = ',ejnp,xi,gnnp
       stop                                   
      endif
      gnnp=c1np*gnnp
      xi=xip
      ierr = 0
      gnps = gncmplx(ejps,xi,1,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in fnzre: ejps,xip,gnps = ',ejps,xi,gnps
       stop
      endif
      gnps=c1ps*gnps
      xi=xip
      ierr = 0
      gnpp = gncmplx(ejpp,xi,1,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in fnzre: ejpp,xip,gnpp = ',ejpp,xi,gnpp
       stop
      endif
      gnpp=c1pp*gnpp
      ierr=0
      call besseldis_cs(gnnp,gnns,in,besn,ierr)
      if(ierr.eq.1)then
       print*,' Error in fnzre: gnns,gnnp,in,besn=',gnns,gnnp,in,besn
       stop
      endif
      ierr=0
      call besseldis_cs(gnpp,gnps,ip,besp,ierr)
      if(ierr.eq.1)then
       print*,' Error in fnzre: gnps,gnpp,ip,besp=',gnps,gnpp,ip,besp
       stop
      endif
      fin = dcmplx(0.d0,0.d0)
      do iq=1,2
       do lam = 2,3
        fin = fin+ff0in(iq,lam,1)*exp(aiunit*eex(iq,lam,1)*x)
     &+ff0in(iq,lam,2)*exp(aiunit*eex(iq,lam,2)*x
     &-(giantwid(lam)*x)**2/2.d0)
       enddo
      enddo
      fff = exp(gnns+gnnp+gnps+gnpp-ep0n-ep0p)*besn*besp
     &*exp(fin-f0inq-(x*dela)**2/2.d0)
      fnzim = imag(fff)/ppnz
      return
      end
