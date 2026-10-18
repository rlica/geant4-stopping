c      Subroutine calculates mean total angular momentum loss
c       for given total kinetic energy with finite energy resolution delta. 
c       The calculation includes giant resonances and low-lying collective
c       states. Initial step in energy = deleloss
c      imm (address for eloss = 0) is calculated if 0 at entry
c      mmmax and mmmin (actual range for eloss and ppe) are calculated
c      At exit ierr = 2 indicates warning about accuracy
c------------------------------------------------------------
      subroutine mavtdist_s(delta,deleloss,accur,imm,mmmax,mmmin,
     &eloss,avmtdis,abserr,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
      external favmtre,favmtim
      common/mav/avma,avmaa,avmans,avmanp,avmaps,avmapp,
     &avmaans,avmaanp,avmaaps,avmaapp,avmatr,avmaatr,
     &avma0,avmaa0,avm0
      common/leg/ecminit,rainbcdg,decm,cs,r,delr,al,theta,
     &cscoul,cslab,thet,thetlabdeg,thetcou,ppcapt,ww,pr0,
     &drifta,driftz,wida,widz,eava,eavaa,tavs,widea,wideaa,avmin,avm,
     &changer,changers,ppscat,totmerg,dthet,
     &deltae,sig0drut,sigreac0,sigreac0lab,cross,alorho,rho,
     &istep,iout,ipopn,norbit,irstep,lab
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/q/qns,qnp,qps,qpp
      common/avm1/avmt1n,avmt1p,avmr1n,avmr1p
      common/bhom/bhomn,bhomp
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/f0muin/ff0muin(2,3,2,4),f0inm(2,3,2)
      common/giantwid/giantwid(3)
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/mavtdisloc/acc1,dela,f0in
      common/mavtloc1/fac1ns,fac1np,fac1ps,fac1pp,
     &fac2ns,fac2np,fac2ps,fac2pp
      dimension avmtdis(ienergy),eloss(ienergy) 
      character*5 version
      common/vers/version
c
      if(version.eq.'vms  '.or.version.eq.'unix* ')then
       open (unit=98,file='grazing_out.dat',status='old')
      endif
      iprep = ierr
      nn = 8
c  i.e. doubling stepwidth each nn's step (must be even)
      dela = delta
      iwarn = 0
c  preparation of factors 
      fac1ns = c1ns*(-avmr1n+2.d0*ejns*avmt1n)
      fac1np = c1np*(avmr1n+2.d0*ejnp*avmt1n)
      fac1ps = c1ps*(-avmr1p+2.d0*ejps*avmt1p)
      fac1pp = c1pp*(avmr1p+2.d0*ejpp*avmt1p)
      fac2ns = c1ns*4.d0*avmt1n
      fac2np = c1np*4.d0*avmt1n
      fac2ps = c1ps*4.d0*avmt1p
      fac2pp = c1pp*4.d0*avmt1p
      f0in = 0.d0
      do n=1,2
       do lam=2,3
        do i=1,2
         f0in = f0in+ff0in(i,lam,n)
        enddo
       enddo
      enddo
      do mm = 1,ienergy
       avmtdis(mm) = 0.d0
      enddo
      sum = 0.d+0
      del = deleloss/2.d0
      acc = accur
      acc1 = acc/100.d0
c  negative values for eloss
      do 438 mm = imm-1,mmmin,-1
       eloss(mm) = eloss(mm+1) - del
       el = eloss(mm)
       if(1.d0-p0.lt.accur)then
        psn=0.d0
        goto 11
       endif
c       write(*,*)' teloss = ',eloss(mm) 
       accc = acc/10.d0
25     ierr1 = 0
       call fourier_s(psn,favmtre,favmtim,el,accc,ierr1)
c       functions favmtre and favmtim defined below
       if(ierr1.eq.1)then
        iwarn = 1
        if(accc.lt.acc)then
         accc=10.d0*accc
         write(*,*)
     &' Accuracy in mavtdist diminished for eloss = ',eloss(mm)
         goto 25
        elseif(1.d0-p0.lt.10.d0*accur)then
         write(*,*)
     &' avm neglected in mavtdist for eloss = ',eloss(mm)
         psn=0.d0
        endif
       endif
11     avmtdis(mm) = psn
c       write(*,*)'   <M> = ',avmtdis(mm),'/MeV'
       ifact = 2*(1+mod(imm-mm,2))
       fact = dble(ifact)/3.d0*del
       if(mod(imm-mm,nn).eq.0.or.abs(mm-imm)*2.eq.nn)then
c    double stepwith every nn steps
        if(iprep.eq.1)goto 3
        fact = del
        del = del*2.d0
        else
       endif
       addt = fact*avmtdis(mm)
       sum = sum + addt
438   continue
c
c   positive values for eloss
      del = deleloss/2.d0 
      do 434, mm = imm,mmmax
       if(mm.eq.imm)then
        eloss(imm) = 0.d0
       else        
        eloss(mm) = eloss(mm-1) + del
       endif
c       write(*,*)' teloss = ',eloss(mm) 
       el = eloss(mm)
       if(1.d0-p0.lt.accur)then
        psn=0.d0
        goto 12
       endif
       accc = acc/10.d0
26     ierr1 = 0
       call fourier_s(psn,favmtre,favmtim,el,accc,ierr1)
c       functions favmtre and favmtim defined below
       if(ierr1.eq.1)then
        iwarn = 1
        if(accc.lt.acc)then
         accc=10.d0*accc
         write(*,*)
     &' Accuracy in mavtdist diminished for eloss = ',eloss(mm)
         goto 26
        elseif(1.d0-p0.lt.10.d0*accur)then
         psn=0.d0
        endif
       endif
12     avmtdis(mm) = psn 
c       write(*,*)'   <M> = ',avmtdis(mm),'/MeV'
       ifact = 2*(1+mod(mm-imm,2))
       fact = dble(ifact)/3.d+0*del
       if(imm.ne.mm.and.mod(imm-mm,nn).eq.0.or.abs(mm-imm)*2.eq.nn)then
c    double stepwith every nn steps
        fact = del
        del = del*2.d0
        else
       endif
       addt = fact*avmtdis(mm)
       sum = sum + addt
434   continue
      abserr = abs(avm-sum)
      if(abserr.gt.10.*accur)then
       write(*,*)' Mean angular mom. not reproduced in mavtdist_s'
       write(*,'(2(a,f9.6))')
     &'  integral probability = ',sum,' avm = ',avm
       iwarn = 1
      endif
3     ierr = 0
      if(iwarn.eq.1)ierr = 2
      return
      end      
c
      real*8 function favmtre(x)
      implicit real*8(a-h,o-z)
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/q/qns,qnp,qps,qpp
      common/bhom/bhomn,bhomp
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/f0muin/ff0muin(2,3,2,4),f0inm(2,3,2)
      common/giantwid/giantwid(3)
      common/mavtdisloc/acc1,dela,f0in
      common/mavtloc1/fac1ns,fac1np,fac1ps,fac1pp,
     &fac2ns,fac2np,fac2ps,fac2pp
      complex*16 gnns,gnnp,gnps,gnpp,fff,zzz,aiunit,func(11),fin,fim
      complex*16 g2ns,g2np,g2ps,g2pp
c
      if(abs(x).gt.100.d0/dela)then
       favmtre = 0.d0      
       return
      endif
      aiunit = dcmplx(0.d0,1.d0)
      xim = -x*bhomn/2.d0
      no = 2
      ierr = 0
      call gn_cs(ejns,xim,no,func,acc1,ierr)
      gnns = func(2)
      g2ns = func(3)
      if(ierr.eq.1)then
       write(*,*)'Error in favmtre: ejns,xim,gnns= ',ejns,xim,func(2)
       stop
      endif
      ierr = 0
      call gn_cs(ejnp,xim,no,func,acc1,ierr)
      gnnp = func(2)
      g2np = func(3)
      if(ierr.eq.1)then
       write(*,*)'Error in favmtre: ejnp,xim,gnnp= ',ejnp,xim,func(2)
       stop
      endif
      xip = -x*bhomp/2.d0
      ierr = 0
      call gn_cs(ejps,xip,no,func,acc1,ierr)
      gnps = func(2)
      g2ps = func(3)
      if(ierr.eq.1)then
       write(*,*)'Error in favmtre: ejps,xip,gnps= ',ejps,xip,func(2)
       stop
      endif
      ierr = 0
      call gn_cs(ejpp,xip,no,func,acc1,ierr)
      gnpp = func(2)
      g2pp = func(3)
      if(ierr.eq.1)then
       write(*,*)'Error in favmtre: ejpp,xip,gnpp= ',ejpp,xip,func(2)
       stop
      endif
      fin = dcmplx(0)
      fim = dcmplx(0)
      do i = 1,2
       do lam = 2,3
        fin = fin+ff0in(i,lam,1)*exp(aiunit*eex(i,lam,1)*x)
     &+ff0in(i,lam,2)*exp(aiunit*eex(i,lam,2)*x
     &-(giantwid(lam)*x)**2/2.d0)
        fim = fim+f0inm(i,lam,1)*exp(aiunit*eex(i,lam,1)*x)
     &+f0inm(i,lam,2)*exp(aiunit*eex(i,lam,2)*x
     &-(giantwid(lam)*x)**2/2.d0)
       enddo
      enddo
      fff = c1ns*gnns*exp(-aiunit*qns*x)+c1np*gnnp*exp(-aiunit*qnp*x)
     &+c1ps*gnps*exp(-aiunit*qps*x)+c1pp*gnpp*exp(-aiunit*qpp*x)
     &+fin-ep0n-ep0p-f0in-dela**2*x**2/2.d0
      zzz = exp(fff)
      zzz = zzz*(exp(-aiunit*qns*x)*(fac1ns*gnns+fac2ns*g2ns)
     &+exp(-aiunit*qnp*x)*(fac1np*gnnp+fac2np*g2np)
     &+exp(-aiunit*qps*x)*(fac1ps*gnps+fac2ps*g2ps)
     &+exp(-aiunit*qpp*x)*(fac1pp*gnpp+fac2pp*g2pp)
     &+fim)
      favmtre = dreal(zzz)
      return
      end
c
      real*8 function favmtim(x)
      implicit real*8(a-h,o-z)
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/q/qns,qnp,qps,qpp
      common/bhom/bhomn,bhomp
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/giantwid/giantwid(3)
      common/f0muin/ff0muin(2,3,2,4),f0inm(2,3,2)
      common/mavtdisloc/acc1,dela,f0in
      common/mavtloc1/fac1ns,fac1np,fac1ps,fac1pp,
     &fac2ns,fac2np,fac2ps,fac2pp
      complex*16 gnns,gnnp,gnps,gnpp,fff,zzz,aiunit,func(11),fin,fim
      complex*16 g2ns,g2np,g2ps,g2pp
c
      if(abs(x).gt.100.d0/dela)then
       favmtim = 0.d0      
       return
      endif
      aiunit = dcmplx(0.d0,1.d0)
      xim = -x*bhomn/2.d+0
      no = 2
      ierr = 0
      call gn_cs(ejns,xim,no,func,acc1,ierr)
      gnns = func(2)
      g2ns = func(3)
      if(ierr.eq.1)then
       write(*,*)'Error in favmtim: ejns,xim,gnns= ',ejns,xim,func(2)
       stop
      endif
      ierr = 0
      call gn_cs(ejnp,xim,no,func,acc1,ierr)
      gnnp = func(2)
      g2np = func(3)
      if(ierr.eq.1)then
       write(*,*)'Error in favmtim: ejnp,xim,gnnp= ',ejnp,xim,func(2)
       stop
      endif
      xip = -x*bhomp/2.d0
      ierr = 0
      call gn_cs(ejps,xip,no,func,acc1,ierr)
      gnps = func(2)
      g2ps = func(3)
      if(ierr.eq.1)then
       write(*,*)'Error in favmtim: ejps,xip,gnps= ',ejps,xip,func(2)
       stop
      endif
      ierr = 0
      call gn_cs(ejpp,xip,no,func,acc1,ierr)
      gnpp = func(2)
      g2pp = func(3)
      if(ierr.eq.1)then
       write(*,*)'Error in favmtim: ejpp,xip,gnpp= ',ejpp,xip,func(2)
       stop
      endif
      fin = dcmplx(0)
      fim = dcmplx(0)
      do i = 1,2
       do lam = 2,3
        fin = fin+ff0in(i,lam,1)*exp(aiunit*eex(i,lam,1)*x)
     &+ff0in(i,lam,2)*exp(aiunit*eex(i,lam,2)*x
     &-(giantwid(lam)*x)**2/2.d0)
        fim = fim+f0inm(i,lam,1)*exp(aiunit*eex(i,lam,1)*x)
     &+f0inm(i,lam,2)*exp(aiunit*eex(i,lam,2)*x
     &-(giantwid(lam)*x)**2/2.d0)
       enddo
      enddo
      fff = c1ns*gnns*exp(-aiunit*qns*x)+c1np*gnnp*exp(-aiunit*qnp*x)
     &+c1ps*gnps*exp(-aiunit*qps*x)+c1pp*gnpp*exp(-aiunit*qpp*x)
     &+fin-ep0n-ep0p-f0in-dela**2*x**2/2.d0
      zzz = exp(fff)
      zzz = zzz*(exp(-aiunit*qns*x)*(fac1ns*gnns+fac2ns*g2ns)
     &+exp(-aiunit*qnp*x)*(fac1np*gnnp+fac2np*g2np)
     &+exp(-aiunit*qps*x)*(fac1ps*gnps+fac2ps*g2ps)
     &+exp(-aiunit*qpp*x)*(fac1pp*gnpp+fac2pp*g2pp)
     &+fim)
      favmtim = dimag(zzz)
      return
      end
