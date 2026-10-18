c      Subroutine calculates probability for kinetic energy loss 
c      greater than given amount = tres for finite energy 
c      resolution = delta. A continuous energy distribution from the
c      low-lying collective states is included. Giants are included
c      as reactions above the resolution.
c------------------------------------------------------------
      subroutine treslim_s(delta,tlim,accur,preres,ierr)
      implicit real*8(a-h,o-z)
      external ftlimrere,ftlimreim 
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/q/qns,qnp,qps,qpp
      common/bhom/bhomn,bhomp
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/treslimloc/acc1,dela,deltap,pn(8),en(8),nnmax,nn
      nn = 0
      dela = delta
      acc1 = accur/1000.d0
      deltap = accur
      n=0 
      do i=1,2
       do lam=2,3
        if(ff0in(i,lam,1).gt.accur)then
         n=n+1
         pn(n) = ff0in(i,lam,1)
         en(n) = eex(i,lam,1)
        endif
       enddo
      enddo
      nnmax = n
      el = tlim
      ierr1 = 0
      acc = accur
      if(abs(1.-p0).lt.accur)goto 10
      call fourier_s(preres,ftlimrere,ftlimreim,el,acc,ierr1)
c       functions ftlimrere and ftlimreim defined below
      if(ierr1.eq.1)goto 2
      goto 4
10    no = 0
      elim = el/sqrt(2.d+0)/dela
      call ierfc(no,elim,pre)
      preres = pre/2.d+0
4     ierr = 0
      return  
2     ierr = 1
      return
      end      
c
      real*8 function ftlimrere(x)
      implicit real*8(a-h,o-z)
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/q/qns,qnp,qps,qpp
      common/bhom/bhomn,bhomp
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/treslimloc/acc1,dela,deltap,pn(8),en(8),nnmax,nn
      complex*16 gnns,gnnp,gnps,gnpp,fff,zzz,aiunit,gncmplx,eexx
      aiunit = dcmplx(0.d+0,1.d+0)
      xim = -x*bhomn/2.d+0
      nn = nn+1
      ierr = 0
      gnns = gncmplx(ejns,xim,1,acc1,ierr)
      if(ierr.eq.1)then
       print*,'Error in ftlimrere: ejns,xim,gnns=',ejns,xim,gnns
       stop
      else      
      end if
      ierr = 0
      gnnp = gncmplx(ejnp,xim,1,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in ftlimrere: ejnp,xim,gnnp=',ejnp,xim,gnnp
       stop
      else      
      end if
      ierr = 0
      gnps = gncmplx(ejps,xim,1,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in ftlimrere: ejps,xim,gnps=',ejps,xim,gnps
       stop
      else      
      end if
      ierr = 0
      gnpp = gncmplx(ejpp,xim,1,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in ftlimrere: ejpp,xim,gnpp=',ejpp,xim,gnpp
       stop
      else      
      end if
      fff = c1ns*gnns*exp(-aiunit*qns*x)+c1np*gnnp*exp(-aiunit*qnp*x)
     &+c1ps*gnps*exp(-aiunit*qps*x)+c1pp*gnpp*exp(-aiunit*qpp*x)
     &-ep0n-ep0p-dela**2*x**2/2.d+0
      if(dreal(fff).lt.-60.d0)then
       eexx = 0.d0
      else
       eexx = exp(fff)
      endif
c      if(mod(nn,20).eq.0)write(98,*)'x, nn, fff =',x,nn,fff
      do n=1,nnmax
       eexx = eexx*(1.d0-aiunit*x*en(n))**(-pn(n))
      enddo
      zzz = -aiunit*eexx/(x-aiunit*deltap)
      ftlimrere = dreal(zzz)
c      if(mod(nn,20).eq.0)write(98,*)' ftlimrere =',ftlimrere
      return
      end

      real*8 function ftlimreim(x)
      implicit real*8(a-h,o-z)
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/q/qns,qnp,qps,qpp 
      common/bhom/bhomn,bhomp
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/treslimloc/acc1,dela,deltap,pn(8),en(8),nnmax,nn
      complex*16 gnns,gnnp,gnps,gnpp,fff,zzz,aiunit,gncmplx,eexx
      aiunit = dcmplx(0.d+0,1.d+0)
      xim = -x*bhomn/2.d+0
      ierr = 0
       gnns = gncmplx(ejns,xim,1,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in ftlimreim: ejns,xim,gnns=',ejns,xim,gnns
       stop
      else      
      end if
      ierr = 0
      gnnp = gncmplx(ejnp,xim,1,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in ftlimreim: ejnp,xim,gnnp=',ejnp,xim,gnnp
       stop
      else      
      end if
      ierr = 0
      gnps = gncmplx(ejps,xim,1,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error inn ftlimreim: ejps,xim,gnps=',ejps,xim,gnps
       stop
      else      
      end if
      ierr = 0
      gnpp = gncmplx(ejpp,xim,1,acc1,ierr)
      if(ierr.eq.1)then
       print*,' Error in ftlimreim: ejpp,xim,gnpp=',ejpp,xim,gnpp
       stop
      else      
      end if
      fff = c1ns*gnns*exp(-aiunit*qns*x)+c1np*gnnp*exp(-aiunit*qnp*x)
     &+c1ps*gnps*exp(-aiunit*qps*x)+c1pp*gnpp*exp(-aiunit*qpp*x)
     &-ep0n-ep0p-dela**2*x**2/2.d+0
      if(dreal(fff).lt.-60.d0)then
       eexx = 0.d0
      else
       eexx = exp(fff)
      endif
      do n=1,nnmax
       eexx = eexx*(1.d0-aiunit*x*en(n))**(-pn(n))
      enddo
      zzz = -aiunit*eexx/(x - aiunit*deltap)
      ftlimreim = dimag(zzz)
c      if(mod(nn,20).eq.0)write(98,*)' ftlimreim =',ftlimreim
      return
      end
