c      Library containing subroutines to manipulate tabulated 
c      functions of one variable
c      Uses subroutines from nag library
c------------------------------------------------------------------
c------------------------------------------------------------
c      Subroutine interpolates a tabulated functions of one variable
c       with nso value of increasing or decreasing x values in
c       sox(1000) and corresponding y values in soy(1000)
c       to produce an equidistant table with steplength delx
c       in rex(1000), rey(1000) with increasing x.
c------------------------------------------------------------------
      subroutine interpolatex_s(sox,soy,nso,rex,rey,delx,ierr)
      implicit real*8(a-h,o-z)
      dimension x(1000),y(1000),rex(1000),rey(1000),dydx(1000)
     &,sox(1000),soy(1000)
      nre = int(abs(sox(nso)-sox(1))/delx + 1.01)
      if(nre.gt.1000.or.nso.gt.1000)then
       write(*,*)' nre.or.nso.gt.1000 in interpolatex_s'
       goto 2
      endif
      do n=1,nre
       rex(n)=0.d0
       rey(n)=0.d0
      enddo
      if(sox(nso).lt.sox(1))then
       do n=1,nso
        x(n)=sox(nso-n+1)
        y(n)=soy(nso-n+1)
       enddo
      else
       do n=1,nso
        x(n)=sox(n)
         y(n)=soy(n)
       enddo
      endif
      ifail = 0
      call E01BEF(nso,x,y,dydx,ifail)
      if(ifail.ne.0)then
       write(*,*)' problem in E01BEF in interpolatex, ifail = ',ifail
       goto 2
      endif
      do n=1,nre-1
       rex(n)=x(1)+dble(n-1)*delx
      enddo
      rex(nre)=x(nso)
      ifail=0
      call E01BFF(nso,x,y,dydx,nre,rex,rey,ifail)
      if(ifail.ne.0)then
       write(*,*)' problem in E01BFF in interpolatex, ifail = ',ifail
       goto 2
      endif
      ierr = 0
      return
2     ierr = 1
      return
      end
c------------------------------------------------------------
c      Subroutine interpolates a tabulated functions of one variable
c       with nso increasing or decreasing x values in sox(n) and 
c       corresponding y values in soy(n) to produce an equidistant 
c       table with steplength delx > 0, starting at xmin and ending at 
c       xmax. The final nred values of x  are in rex(i). The final (nred)
c       values of y in rey(i) and its derivative dydx(i) are 0 
c       outside the interval where y is originally defined.
c------------------------------------------------------------------
      subroutine interpolateder_s(sox,soy,nso,xmin,xmax,delx,
     &rex,rey,dydx,nred,ierr)
      implicit real*8(a-h,o-z)
      dimension x(1000),y(1000),dd(1000),rex(nred),rey(nred),dydx(nred)
     &,sox(nso),soy(nso),rx(1000),ry(1000),dyd(1000)
      nre = int((xmax-xmin)/delx +1.01)
      if(nre.gt.nred)then
       write(*,*)' nred less than dimension of rex in interpolateder'
       goto 2
      endif       
      if(nre.gt.1000.or.nso.gt.1000)then
       write(*,*)' nre or nso.gt.1000 in interpolateder_s'
      elseif(nre.lt.0.or.nso.lt.0)then
       write(*,*)' nre or nso.lt.0 in interpolateder_s'
      else      
       goto 3
      endif
      goto 2
3     do n=1,nred
       rex(n)=0.d0
       rey(n)=0.d0
      enddo
      if(sox(nso).lt.sox(1))then
       do n=1,nso
        x(n)=sox(nso-n+1)
        y(n)=soy(nso-n+1)
       enddo
      else
       do n=1,nso
        x(n)=sox(n)
        y(n)=soy(n)
       enddo
      endif
      ifail = 1
      call E01BEF(nso,x,y,dd,ifail)
      if(ifail.ne.0)then
       write(*,*)' problem in E01BEF in interpolateder, ifail = ',ifail
       goto 2
      endif
      m=1
      nmin=1
      nmax=nre
      do n=1,nre
       rex(n)=xmin+dble(n-1)*delx
       if(n.eq.nre)rex(nre)=xmax
       if(rex(n).gt.x(1).and.rex(n).lt.x(nso))then
        if(nmin.eq.1)nmin=n
        rx(m) = rex(n)
        m=m+1
       else
        if(rex(n).gt.x(1).and.nmax.eq.nre)nmax=n-1
        rey(n) = 0.d0
        dydx(n) = 0.d0
       endif
      enddo
      nnr = nmax-nmin+1
      ifail=0
      call E01BGF(nso,x,y,dd,nnr,rx,ry,dyd,ifail)
      if(ifail.ne.0)then
       write(*,*)' problem in E01BGF in interpolateder, ifail = ',ifail
       goto 2
      endif
      do n=nmin,nmax
       rey(n) = ry(n-nmin+1)
       dydx(n) = dyd(n-nmin+1)
      enddo
      ierr = 0
      return
2     ierr = 1
      return
      end
