c      Subroutine for distribution in N and Z  and excitation energy 
c       mean values and spreads due to transfer reactions
c       irangen and irangep must be le.ibesrange
c      If at entry in00=ip00=0, in00 and ip00 are calculated
c       otherwise they are input
c      If at entry ierr = 1 the distribution is printed
c      At exit ierr = 2 is warning for accuracy
c-------------------------------------------------------------------
c
      subroutine nzdistri_s(in00,ip00,inmax,inmin,izmax,izmin,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
      parameter (nrange=100,ibesrange=2*nrange+1,nn0=nrange+1)
c      irangen and irangep are defined in grazing.icl
      common/za/za,aa,zaa,aaa
      common/ecm/ecm
      common/vdr/vdr,accur 
      common/q/qns,qnp,qps,qpp
      common/ba/bna,bnaa,bpa,bpaa
      common/f0/f0ina,f0inaa,f0ns,f0np,f0ps,f0pp
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/eav/eav,eavns,eavnp,eavps,eavpp,eavtr,
     &eav0,eava0,eavaa0
      common/mav/avma,avmaa,avmans,avmanp,avmaps,avmapp,
     &avmaans,avmaanp,avmaaps,avmaapp,avmatr,avmaatr,
     &avma0,avmaa0,avm0
      common/wid/wid2e,wid2ns,wid2np,wid2ps,wid2pp,wid2tr,wid2ea,
     &wid2eaa,wid2e0,wid2ea0,wid2eaa0
      common/nzdis/ p0nn(irangen),p0zz(irangep),
     &eavdist(irangen,irangep),eavdist00,avmadist(irangen,irangep),
     &avmadist00,avmdist(irangen,irangep),avmdist00,
     &wid2dist(irangen,irangep),wid2dist00,wid2adist(irangen,irangep)
      dimension besdis(ibesrange),avmaadist(irangen,irangep)

c
      iout = ierr
c      distribution in N and Z
      p0tn = 0.d0
      p0tz = 0.d0
      do n=1,irangen
       p0nn(n) = 0.d0
      enddo
      do m=1,irangep
       p0zz(m) = 0.d0
       do n=1,irangen
        eavdist(n,m) = 0.d0
        avmadist(n,m) = 0.d0
        avmaadist(n,m) = 0.d0
        avmdist(n,m) = 0.d0
        wid2dist(n,m) = 0.d0
        wid2adist(n,m) = 0.d0
       enddo
      enddo       
7     ierr1 = 0
      acc = accur/1000.d0
      call besseldist_s(f0np,f0ns,nmin,nmax,besdis,acc,ierr1)
      if(ierr1.eq.1)goto 2
      acu = accur/100.d0
74    do i = nn0+nmin,nn0
c  (nn0 = (ibesrange+1)/2 position for N=0 is defined as parameter above)
       if(besdis(i).gt.acu)then 
        goto 76
       endif
      enddo
      i = nn0
76    inmin = i-nn0
      if(in00.eq.0.and.inmin.ge.-1)inmin=-1
      do i = nn0+nmax,nn0,-1
       if(besdis(i).gt.acu)then
        goto 75
       endif
      enddo
      i = nn0
75    inmax = i-nn0
      if(in00.eq.0.and.inmax.le.1)inmax=1
      if(inmax-inmin.gt.irangen-5)then
       acu = 2.d0*acu
       ierr = 2
       if(acu.ge.1.d-2)then
        write(*,*)' You must change irangen'
        ierr = 1
        goto 2
       endif
       goto 74
      endif
      if(ierr.eq.2)write(*,*)
     &' Reduced accuracy: inmax-inmin out of range in nzdistri_s'
      if(in00.eq.0)then
       in00 = int(dble(-inmin*irangen)/dble(inmax-inmin))+1
      else
       if(inmin.lt.1-in00)then
        inmin = 1-in00
        write(*,*)' Warning: inmin.lt.previous inmin in nzdistri'
       endif
       if(inmax+in00+2.gt.irangen)then
        inmax = irangen-in00-2
        write(*,*)' Warning: inmax.gt.previous inmax in nzdistri'
       endif        
      endif
      if(iout.eq.1)write(*,'(/)')
      do i=inmin-2,inmax+2
       n=in00+i
       p0nn(n)=besdis(nn0+i)
       if(iout.eq.1)write(*,'(a,i3,a,d12.5)')' P0n(',i,') = ',p0nn(n)
       p0tn=p0tn+p0nn(in00+i)
      enddo
      if(iout.eq.1)
     & write(*,'(/a,f10.7,/)')' Total probability p0tn = ',p0tn
      acc = accur/100.d0
8     ierr1 = 0
      call besseldist_s(f0pp,f0ps,nmin,nmax,besdis,acc,ierr1)
      if(ierr1.eq.1)goto 2
      acu = accur/20.d0
84    do 87 i = nn0+nmin,nn0
       if(besdis(i).gt.acu)then
        goto 86
       endif
87    continue
      i = nn0
86    izmin = i-nn0
      if(ip00.eq.0.and.izmin.gt.-1)izmin=-1
      do 88 i = nn0+nmax,nn0,-1
       if(besdis(i).gt.acu)then
        goto 85
       endif
88    continue
85    izmax = i-nn0
      if(ip00.eq.0.and.izmax.lt.1)izmax=1
      if(izmax-izmin.gt.irangep-5)then
       acu = 2.d0*acu
       if(acc.ge..01d0)then
        write(*,*)' You must change irangep'
        ierr = 1
        goto 2
       endif
       ierr = 3
       goto 84
      endif
      if(ierr.eq.3)then
       write(*,*)
     &' Reduced accuracy: izmax-izmin out of range in nzdistri_s'
       ierr = 2
      endif
      if(ip00.eq.0)then
       ip00 = int(dble(-izmin*irangep)/dble(izmax-izmin))+1
      else
       if(izmin.lt.1-ip00)then
        izmin = 1-ip00
        write(*,*)' Warning: izmin.lt.previous izmin in nzdistri'
       endif
       if(izmax+ip00+2.gt.irangep)then
        izmax = irangep-ip00-2
        write(*,*)' Warning: izmax.gt.previous izmax in nzdistri'
       endif
      endif
      do i=izmin-2,izmax+2
       m=ip00+i
       p0zz(m)=besdis(nn0+i)
       if(iout.eq.1)write(*,'(a,i3,a,d12.5)')' P0z(',i,') = ',p0zz(m)
       p0tz=p0tz+p0zz(ip00+i) 
      enddo
      if(iout.eq.1)
     & write(*,'(/a,f7.4)')'         Total probability p0tz = ',p0tz
c      write(*,'(a,6i3)')
c     &' in00,ip00,inmax,inmin,izmax,izmin = ',
c     & in00,ip00,inmax,inmin,izmax,izmin
      cheav = 0.d0
      che2 = 0.d0
      che2a = 0.d0
      chavma = 0.d0
      chavmaa = 0.d0   
      do 441,i=inmin-1,inmax+1 
       do 441,j=izmin-1,izmax+1 
        n = in00+i
        m = ip00+j
        if(p0nn(n).eq.0.d0)goto 441
        if(p0zz(m).eq.0.d0)goto 441
        eavdist(n,m)=(eavnp*p0nn(n-1)+eavns*p0nn(n+1))/p0nn(n)
     &             + (eavpp*p0zz(m-1)+eavps*p0zz(m+1))/p0zz(m)
        avmadist(n,m)=(avmanp*p0nn(n-1)+avmans*p0nn(n+1))/p0nn(n)
     &             + (avmapp*p0zz(m-1)+avmaps*p0zz(m+1))/p0zz(m)
        avmaadist(n,m)=(avmaanp*p0nn(n-1)+avmaans*p0nn(n+1))/p0nn(n)
     &             + (avmaapp*p0zz(m-1)+avmaaps*p0zz(m+1))/p0zz(m)
        avmdist(n,m)=avmadist(n,m)+avmaadist(n,m)
c      write(*,'(a,2i3,2f8.2)')' i,j, avmadist,avmaadist = ',i,j
c     &,avmadist(n,m),avmaadist(n,m)
c      sumcheck of eav, and avm 
       cheav = cheav + p0nn(n)*p0zz(m)*eavdist(n,m)
       chavma = chavma + p0nn(n)*p0zz(m)*avmadist(n,m)
       chavmaa = chavmaa + p0nn(n)*p0zz(m)*avmaadist(n,m)
441   continue
      do 442,i=inmin,inmax 
       do 442,j=izmin,izmax 
        n = in00+i
        m = ip00+j
        if(p0nn(n).eq.0.d0)goto 442
        if(p0zz(m).eq.0.d0)goto 442
         eavn = (eavnp*p0nn(n-1)+eavns*p0nn(n+1))/p0nn(n)
         eavp = (eavpp*p0zz(m-1)+eavps*p0zz(m+1))/p0zz(m)
         wid2dist(n,m)= (p0nn(n-1)*wid2np+p0nn(n+1)*wid2ns)/p0nn(n)
     &    +(p0nn(n-2)*eavnp**2+p0nn(n+2)*eavns**2)/p0nn(n)
     &    - eavn**2 + 2.d0*eavns*eavnp
     &                + (p0zz(m-1)*wid2pp+p0zz(m+1)*wid2ps)/p0zz(m)
     &    +(p0zz(m-2)*eavpp**2+p0zz(m+2)*eavps**2)/p0zz(m)
     &    - eavp**2 + 2.d0*eavps*eavpp
         wid2adist(n,m)= 
     &    ((p0nn(n-1)*wid2np+p0nn(n+1)*wid2ns)/p0nn(n))/3.d0
     &    +((p0nn(n-2)*eavnp**2+p0nn(n+2)*eavns**2)/p0nn(n)
     &    - eavn**2 + 2.d0*eavns*eavnp)/4.d0
     &  + ((p0zz(m-1)*wid2pp+p0zz(m+1)*wid2ps)/p0zz(m))/3.d0
     &    +((p0zz(m-2)*eavpp**2+p0zz(m+2)*eavps**2)/p0zz(m)
     &    - eavp**2 + 2.d0*eavps*eavpp)/4.d0
c      print*,' n,m,wid2dist(n,m) = ',n,m,wid2dist(n,m)
c  sumcheck of  wid2e
       che2=che2+p0nn(n)*p0zz(m)*(wid2dist(n,m)+eavdist(n,m)**2)
       che2a=che2a+p0nn(n)*p0zz(m)*(wid2adist(n,m)+eavdist(n,m)**2/4.d0)
442   continue
      chwid2 = che2 - cheav**2
      chwid2a = che2a - cheav**2/4.d0
c
      k = in00
      kk = ip00
      denom = (p0nn(k)*p0zz(kk)-p0tr)
      if(denom.eq.0.d0)then
       write(*,*)' eavdist00 can not be calculated'
       ierr = 2
       goto 3
      endif
      fact = p0nn(k)*p0zz(kk)/denom
      eavdist00 = eavdist(k,kk)*fact
      avmadist00 = avmadist(k,kk)*fact
      avmdist00 = avmdist(k,kk)*fact
      wid2dist00 = (wid2dist(k,kk)+eavdist(k,kk)**2)*fact-eavdist00**2
      wid2dist00alt = fact*(wid2dist(k,kk)-p0tr/denom*eavdist(k,kk)**2)
c      print*,' wid2dist00,wid2dist00alt=',wid2dist00,wid2dist00alt 
      if(wid2dist00.lt.0.d0.and.p0.gt.1.d0-10.d0*accur)then
       write(*,*)' Warning: wid2dist00.lt.0 in nzdistri_s'
       ierr = 2
       goto 3
      endif
c
      eavtr = eavns+eavnp+eavps+eavpp
      if(abs(1.d0-cheav/eavtr).gt.1.d1*accur)then
       write(*,'(a,f9.4,a,f9.4)')
     &' Warning: Sumcheck eav =',cheav,' eavtr =',eavtr
       ierr = 2
       goto 3
      endif
      avmatr = avmans + avmanp + avmaps + avmapp 
      if(abs(1.d0-chavma/avmatr).gt.1.d1*accur)then
       write(*,'(a,f9.4,a,f9.4)')
     &' Warning: Sumcheck avma =',chavma,' avmatr =',avmatr
       ierr = 2
       goto 3
      endif 
      avmaatr = avmaans + avmaanp + avmaaps + avmaapp 
      if(abs(1.d0-chavmaa/avmaatr).gt.1.d1*accur)then
       write(*,'(a,f9.4,a,f9.4)')
     &' Warning: Sumcheck avmaa =',chavmaa,' avmaatr =',avmaatr
       ierr = 2
       goto 3
      endif
      wid2etr = wid2ns + wid2np + wid2ps + wid2pp 
      if(abs(1.d0-chwid2/wid2etr).gt.1.d1*accur)then
       write(*,'(a,f9.4,a,f9.4)')
     &' Warning: Sumcheck wid2e =',chwid2,' wid2etr =',wid2etr
       ierr = 2
       goto 3
      endif
      wid2eatr = wid2etr/3.d0
      if(abs(1.d0-chwid2a/wid2eatr).gt.1.d1*accur)then
       write(*,'(a,f9.4,a,f9.4)')
     &' Warning: Sumcheck wid2ea =',chwid2a,' wid2eatr =',wid2eatr
       ierr = 2
       goto 3
      endif
c
      ierr = 0
      return
2     ierr = 1
3     return
      end
