c      Subroutine calculates total kinetic energy distribution 
c       with finite energy resolution delta. The calculation
c       includes low-lying collective states and giant resonances
c      imm (address for eloss = 0) is calculated if 0 at entry
c      mmmax and mmmin (actual range for eloss and ppe) are calculated
c      At entry ierr = 1 indicates preparatory run
c      At exit ierr = 2 indicates warning about accuracy
c------------------------------------------------------------
      subroutine tdistri_s(delta,deleloss,accur,imm,mmmax,mmmin,
     &eloss,ppe,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
      external ftre,ftim
      dimension ppe(ienergy),eloss(ienergy) 
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/q/qns,qnp,qps,qpp
      common/bhom/bhomn,bhomp
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/giantwid/giantwid(3)
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/f0/f0ina,f0inaa,f0ns,f0np,f0ps,f0pp
      common/tdisloc/acc1,dela,f0in
      character*5 version
      common/vers/version
c
      iprep = ierr
      if(version.eq.'vms  '.or.version.eq.'unix* ')then
       open (unit=98,file='grazing_out.dat',status='old')
      endif
      nn = 8
c   doubling stepwidth each nn's step (must be even)
      f0in = 0.d0
      do i = 1,2
       do lam = 2,3
        do n = 1,2
         f0in = f0in + ff0in(i,lam,n)
        enddo
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
      sum = 0.d+0
      factor = 1.d+0/sqrt(6.28318530717958647d0)/dela
      if(imm.eq.0)imm = int(dble(ienergy)/4.)
c  address for eloss = 0
      eloss(imm) = 0.d0
c  initial stepwidth
      del = deleloss/2.d0
cnan      del = deleloss
c  negative values for eloss
      do mm = imm-1,1,-1
       eloss(mm) = eloss(mm+1) - del
       el = eloss(mm)
       if(abs(1.d0-p0).lt.accur)then
        psn=factor*exp(-el**2/2.d0/dela**2)
        goto 11
       endif
c       write(*,*)' teloss = ',eloss(mm) 
       accc=acc
25     ierr1 = 0
       call fourier_s(psn,ftre,ftim,el,accc,ierr1)
c       functions ftre and ftim defined below
       if(ierr1.eq.1)then
        iwarn = 1
        if(accc.lt.acc*1.01)then
         accc=10.d0*accc
         goto 25
        elseif(1.d0-p0.lt.10.d0*accur)then
         psn=factor*exp(-el**2/2.d0/dela**2)
        endif
       endif
11     ppe(mm) = psn 
c       write(*,*)'   ppt = ',ppe(mm),'/MeV'
       ifact = 2*(1+mod(imm-mm,2))
       fact = dble(ifact)/3.d0*del
       if(mod(imm-mm,nn).eq.0.or.abs(mm-imm)*2.eq.nn)then
c   double stepwidth every nn steps
        if(iprep.eq.1)goto 3
        fact = del
        del = del*2.d0
        else
       endif
       addt = fact*ppe(mm)      
       sum = sum + addt
       if(addt.lt.acc)goto 435
      enddo
      mm=1
      write(*,'(a)')
     &' Warning: Attempt to calculate below 1 in tdistri_s'
      write(98,'(a)')
     &' Warning: Attempt to calculate below 1 in tdistri_s'
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
        eloss(imm) = 0.d+0
       else        
        eloss(mm) = eloss(mm-1) + del
       endif
c       write(*,*)' teloss = ',eloss(mm) 
       el = eloss(mm)
       if(abs(1.d0-p0).lt.accur)then
        psn=factor*exp(-el**2/2.d+0/dela**2)
        goto 12
       endif
       accc=acc
26     ierr1 = 0
       call fourier_s(psn,ftre,ftim,el,accc,ierr1)
c       functions ftre and ftim defined below
       if(ierr1.eq.1)then
        iwarn = 1
        if(accc.lt.acc*1.01)then
         accc=10.d0*accc
         goto 26
        elseif(1.d0-p0.lt.10.d0*accur)then
         psn=factor*exp(-el**2/2.d0/dela**2)
        endif
       endif
12     ppe(mm) = psn 
c       write(*,*)'   ppt = ',ppe(mm),'/MeV'
       ifact = 2*(1+mod(mm-imm,2))
       fact = dble(ifact)/3.d0*del
       if(imm.ne.mm.and.mod(imm-mm,nn).eq.0.or.abs(mm-imm)*2.eq.nn)then
c    double stepwith every nn steps
        fact = del
        del = del*2.d0
        else
       endif
       addt = fact*ppe(mm)      
       sum = sum + addt
       if(abs(1.d0-p0).lt.accur.and.addt.gt.acc/100.d0)goto 434
       if(sum.gt..2d0.and.addt.lt.acc)then
        goto 436
       elseif(sum.lt.accur)then
        mmmin = mm
       endif
434    continue
      enddo
      mm=ienergy
      write(*,'(a)')
     &' Warning: Attempt to calculate above ienergy in tdistri_s'
      iwarn = 1
436   mmmax = mm-1
      do mm=mmmax+1,ienergy
       eloss(mm)=eloss(mm-1)+del
       ppe(mm) = 0.d0
       if(imm.ne.mm.and.mod(imm-mm,nn).eq.0.or.abs(mm-imm)*2.eq.nn)then
       del = del*2.d0
       endif
      enddo
      if(abs(1.-sum).gt.10.*accur)then
       write(*,*)' Total probability in tdistri - 1 .gt. 10*accur'
       write(*,'(a,f9.6)')'  integral probability plus elastic =',sum
       iwarn = 1
       sum=1.                 ! messo ad aprile 2002
      endif
3     ierr = 0
      if(iwarn.eq.1)ierr = 2
      return
      end      
c      
      real*8 function ftre(x)
      implicit real*8(a-h,o-z)
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/q/qns,qnp,qps,qpp
      common/bhom/bhomn,bhomp
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/giantwid/giantwid(3)
      common/tdisloc/acc1,dela,f0in
      complex*16 gnns,gnnp,gnps,gnpp,fff,zzz,aiunit,gncmplx,fin
c
      aiunit = dcmplx(0.d0,1.d0)
      xim = -x*bhomn/2.d0
      ierr = 0
      gnns = gncmplx(ejns,xim,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)'Error in ftre: ejns,xim,gnns=',ejns,xim,gnns
       stop
      endif
      ierr = 0
      gnnp = gncmplx(ejnp,xim,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)'Error in ftre: ejnp,xim,gnnp=',ejnp,xim,gnnp
       stop
      endif
      xip = -x*bhomp/2.d0
      ierr = 0
      gnps = gncmplx(ejps,xip,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)'Error in ftre: ejps,xip,gnps=',ejps,xip,gnps
       stop
      endif
      ierr = 0
      gnpp = gncmplx(ejpp,xip,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)'Error in ftre: ejpp,xip,gnpp=',ejpp,xip,gnpp
       stop
      endif
      fin = dcmplx(0.d0,0.d0)
      do i = 1,2
       do lam = 2,3
        fin = fin+ff0in(i,lam,1)*exp(aiunit*eex(i,lam,1)*x)
     &+ff0in(i,lam,2)*exp(aiunit*eex(i,lam,2)*x
     &-(giantwid(lam)*x)**2/2.d0)
       enddo
      enddo
      fff = c1ns*gnns*exp(-aiunit*qns*x)+c1np*gnnp*exp(-aiunit*qnp*x)
     &+c1ps*gnps*exp(-aiunit*qps*x)+c1pp*gnpp*exp(-aiunit*qpp*x)
     &+fin-ep0n-ep0p-f0in-dela**2*x**2/2.d0
      zzz = exp(fff)
      ftre = dreal(zzz)
      return
      end
c
      real*8 function ftim(x)
      implicit real*8(a-h,o-z)
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/q/qns,qnp,qps,qpp
      common/bhom/bhomn,bhomp
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/giantwid/giantwid(3)
      common/tdisloc/acc1,dela,f0in
      complex*16 gnns,gnnp,gnps,gnpp,fff,zzz,aiunit,gncmplx,fin
c
      aiunit = dcmplx(0.d0,1.d0)
      xim = -x*bhomn/2.d0
      ierr = 0
      gnns = gncmplx(ejns,xim,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)'Error in ftim: ejns,xim,gnns=',ejns,xim,gnns
       stop
      endif
      ierr = 0
      gnnp = gncmplx(ejnp,xim,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)'Error in ftim: ejnp,xim,gnnp=',ejnp,xim,gnnp
       stop
      endif
      xip = -x*bhomp/2.d0
      ierr = 0
      gnps = gncmplx(ejps,xip,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)'Error in ftim: ejps,xip,gnps=',ejps,xip,gnps
       stop
      endif
      ierr = 0
      gnpp = gncmplx(ejpp,xip,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)'Error in ftim: ejpp,xip,gnpp=',ejpp,xip,gnpp
       stop
      endif
      fin = dcmplx(0.d0,0.d0)
      do i = 1,2
       do lam = 2,3
        fin = fin+ff0in(i,lam,1)*exp(aiunit*eex(i,lam,1)*x)
     &+ff0in(i,lam,2)*exp(aiunit*eex(i,lam,2)*x
     &-(giantwid(lam)*x)**2/2.d0)
       enddo
      enddo
      fff = c1ns*gnns*exp(-aiunit*qns*x)+c1np*gnnp*exp(-aiunit*qnp*x)
     &+c1ps*gnps*exp(-aiunit*qps*x)+c1pp*gnpp*exp(-aiunit*qpp*x)
     &+fin-ep0n-ep0p-f0in-dela**2*x**2/2.d0
      zzz = exp(fff)
      ftim = dimag(zzz)
      return
      end
