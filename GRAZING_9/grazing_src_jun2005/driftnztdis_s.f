c      Subroutine calculates mean neutron and proton drift for given 
c       total kinetic energy with finite energy resolution delta. 
c       The calculation includes giant resonances and includes 
c       low-lying collective states which are not selected for 
c       special treatment. The selected states are specified 
c       in common/selec/  to be included later
c      delta, deleloss and eloss are input as well as mgmax 
c       (number of giant res.), imm (address for eloss = 0), 
c       and mmmax and mmmin (actual range for eloss and ppe).
c      At exit ierr = 2 indicates warning about accuracy
c------------------------------------------------------------
      subroutine driftnztdis_s(delta,deleloss,accur,imm,mmmax,mmmin,
     &eloss,driftndis,driftzdis,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
      external fdntre,fdntim,fdztre,fdztim 
      common/leg/ecminit,rainbcdg,decm,cs,r,delr,al,theta,
     &cscoul,cslab,thet,thetlabdeg,thetcou,ppcapt,ww,pr0,
     &drifta,driftz,wida,widz,eava,eavaa,tavs,widea,wideaa,avmin,avm,
     &changer,changers,ppscat,totmerg,dthet,
     &deltae,sig0drut,sigreac0,sigreac0lab,cross,alorho,rho,
     &istep,iout,ipopn,norbit,irstep,lab
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/q/qns,qnp,qps,qpp
      common/bhom/bhomn,bhomp
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/giantwid/giantwid(3)
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/tdisloc/acc1,dela,f0in
      dimension driftndis(ienergy),driftzdis(ienergy),eloss(ienergy) 
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
      f0in = 0.d0
      do n=1,2
       do lam=2,3
        do i=1,2
         f0in = f0in+ff0in(i,lam,n)
        enddo
       enddo
      enddo
      iwarn = 0
      do mm= 1,ienergy
       driftndis(mm)= 0.d0
       driftzdis(mm)= 0.d0
      enddo
      sumn = 0.d0
      sumz = 0.d0
      del = deleloss/2.d0 
      acc = accur/10.d0
      acc1 = acc/100.d0
c  negative values of eloss
      do 438, mm = imm-1,mmmin,-1
       el = eloss(mm)
       if(1.d0-p0.lt.accur)then
        psn=0.d0
        psm=0.d0
c        write(*,*)' teloss = ',eloss(mm) 
        goto 11
       endif
c       write(*,*)' teloss = ',eloss(mm) 
       accc = acc
25     ierr1 = 0
       call fourier_s(psn,fdntre,fdntim,el,accc,ierr1)
c       functions fdntre and fdntim defined below
       if(ierr1.eq.1)then
        iwarn = 1
        if(accc.lt.acc*1.01)then
         accc=10.d0*accc
         goto 25
        elseif(1.d0-p0.lt.10.d0*accur)then
         psn=0.d0
        endif
       endif
       accc = acc
26     ierr1 = 0
       call fourier_s(psm,fdztre,fdztim,el,accc,ierr1)
c       functions fdztre and fdztim defined below
       if(ierr1.eq.1)then
        iwarn = 1
        if(accc.lt.acc*1.01)then
         accc=10.d0*accc
         goto 26
        elseif(1.d0-p0.lt.10.d0*accur)then
         psm=0.d0
        endif
       endif
11     driftndis(mm) = psn 
       driftzdis(mm) = psm 
c       write(*,*)'   driftn = ',driftndis(mm),'/MeV'
c       write(*,*)'   driftz = ',driftzdis(mm),'/MeV'
       ifact = 2*(1+mod(imm-mm,2))
       fact = dble(ifact)/3.d0*del
       if(mod(imm-mm,nn).eq.0.or.abs(mm-imm)*2.eq.nn)then
c    double stepwith every nn steps
        if(iprep.eq.1)goto 3
        fact = del
        del = del*2.d0
        else
       endif
       addn = fact*driftndis(mm)      
       addz = fact*driftzdis(mm)      
       sumn = sumn + addn
       sumz = sumz + addz
438   continue
c
c   positive values for eloss
      del = deleloss/2.d0 
      do 434, mm = imm,mmmax
c       write(*,*)' teloss = ',eloss(mm) 
       el = eloss(mm)
       if(1.d0-p0.lt.accur)then
        psn=0.d0
        psm=0.d0
        goto 12
       endif
       accc = acc
27     ierr1 = 0
       call fourier_s(psn,fdntre,fdntim,el,accc,ierr1)
c       functions fdntre and fdntim defined below
       if(ierr1.eq.1)then
        iwarn = 1
        if(accc.lt.acc*1.01)then
         accc=10.d0*accc
         goto 27
        elseif(1.d0-p0.lt.10.d0*accur)then
         psn=0.d0
        endif
       endif
       accc = acc
28     ierr1 = 0
       call fourier_s(psm,fdztre,fdztim,el,accc,ierr1)
c       functions fdztre and fdztim defined below
       if(ierr1.eq.1)then
        iwarn = 1
        if(accc.lt.acc*1.01)then
         accc=10.d0*accc
         goto 28
        elseif(1.d0-p0.lt.10.d0*accur)then
         psm=0.d0
        endif
       endif
12     driftndis(mm) = psn 
       driftzdis(mm) = psm 
c       write(*,*)'   driftn = ',driftndis(mm),'/MeV'
c       write(*,*)'   driftz = ',driftzdis(mm),'/MeV'
       ifact = 2*(1+mod(mm-imm,2))
       fact = dble(ifact)/3.d+0*del
       if(imm.ne.mm.and.mod(imm-mm,nn).eq.0.or.abs(mm-imm)*2.eq.nn)then
c    double stepwith every nn steps
        fact = del
        del = del*2.d0
        else
       endif
       addn = fact*driftndis(mm)      
       addz = fact*driftzdis(mm)      
       sumn = sumn + addn
       sumz = sumz + addz
434   continue
      suma = sumz+sumn
      if(abs(drifta-suma).gt.10.*accur)then
       write(*,*)' Mean drift in A not reproduced in driftnztdis_s'
       write(*,'(2(a,f9.6))')
     &'  integral = ',suma,' drifta = ',drifta
       iwarn = 1
      endif
      if(abs(driftz-sumz).gt.10.*accur)then
       write(*,*)' Mean drift in Z not reproduced in driftnztdis_s'
       write(*,'(2(a,f9.6))')
     &'  integral = ',sumz,' driftz = ',driftz
       iwarn = 1
      endif
3     ierr = 0
      if(iwarn.eq.1)ierr = 2
      return
      end      
c
      real*8 function fdntre(x)
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
       write(*,*)' Error in fdntre: ejns,xim,gnns = ',ejns,xim,gnns
       stop
      endif
      ierr = 0
      gnnp = gncmplx(ejnp,xim,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)' Error in fdntre: ejnp,xim,gnnp = ',ejnp,xim,gnnp
       stop
      endif
      xip = -x*bhomp/2.d0
      ierr = 0
      gnps = gncmplx(ejps,xip,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)' Error in fdntre: ejps,xip,gnps=',ejps,xip,gnps
       stop
      endif
      ierr = 0
      gnpp = gncmplx(ejpp,xip,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)' Error in fdntre: ejpp,xip,gnpp=',ejpp,xip,gnpp
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
      zzz = zzz*(c1np*gnnp-c1ns*gnns)
      fdntre = dreal(zzz)
      return
      end
c
      real*8 function fdntim(x)
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
       write(*,*)' Error in fdntim: ejns,xim,gnns = ',ejns,xim,gnns
       stop
      endif
      ierr = 0
      gnnp = gncmplx(ejnp,xim,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)' Error in fdntim: ejnp,xim,gnnp = ',ejnp,xim,gnnp
       stop
      endif
      xip = -x*bhomp/2.d0
      ierr = 0
      gnps = gncmplx(ejps,xip,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)' Error in fdntim: ejps,xip,gnps = ',ejps,xip,gnps
       stop
      endif
      ierr = 0
      gnpp = gncmplx(ejpp,xip,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)'Error in fdntim: ejpp,xip,gnpp = ',ejpp,xip,gnpp
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
      zzz = zzz*(c1np*gnnp-c1ns*gnns)
      fdntim = dimag(zzz)
      return
      end
c
      real*8 function fdztre(x)
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
       write(*,*)' Error in fdztre: ejns,xim,gnns = ',ejns,xim,gnns
       stop
      endif
      ierr = 0
      gnnp = gncmplx(ejnp,xim,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)' Error in fdztre: ejnp,xim,gnnp = ',ejnp,xim,gnnp
       stop
      endif
      xip = -x*bhomp/2.d0
      ierr = 0
      gnps = gncmplx(ejps,xip,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)' Error in fdztre: ejps,xip,gnps = ',ejps,xip,gnps
       stop
      endif
      ierr = 0
      gnpp = gncmplx(ejpp,xip,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)' Error in fdztre: ejpp,xip,gnpp = ',ejpp,xip,gnpp
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
      zzz = zzz*(c1pp*gnpp-c1ps*gnps)
      fdztre = dreal(zzz)
      return
      end
c
      real*8 function fdztim(x)
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
       write(*,*)' Error in fdztim: ejns,xim,gnns = ',ejns,xim,gnns
       stop
      endif
      ierr = 0
      gnnp = gncmplx(ejnp,xim,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)' Error in fdztim: ejnp,xim,gnnp = ',ejnp,xim,gnnp
       stop
      endif
      xip = -x*bhomp/2.d0
      ierr = 0
      gnps = gncmplx(ejps,xip,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)' Error in fdztim: ejps,xip,gnps = ',ejps,xip,gnps
       stop
      endif
      ierr = 0
      gnpp = gncmplx(ejpp,xip,1,acc1,ierr)
      if(ierr.eq.1)then
       write(*,*)' Error in fdztim: ejpp,xip,gnpp = ',ejpp,xip,gnpp
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
      zzz = zzz*(c1pp*gnpp-c1ps*gnps)
      fdztim = dimag(zzz)
      return
      end
