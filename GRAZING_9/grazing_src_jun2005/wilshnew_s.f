c      Subroutine WILSH produces table for Wilshinsky plot 
c      If at entry ierr = lab = 1 it produces the plot in the lab system
c      If at entry ierr = lab = 0 it produces it in the center of mass 
c       la = lab if plot in dsig/dO (mb/sterad)
c       la = lab+2 if plot in dsig/d theta (mb/rad)
c      At exit ierr = 1 for error
c              ierr = 2 for warning
c---------------------------------------------------------------------
      subroutine wilsh_s(thav,thavlab,dsr,rstep,tloss,widang,ecut
     &,angt,angd,la,nfile,minm,maxm,irstmax,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/ecm/ecm
      common/capt/eb,rbar,i5,intp
      dimension thav(ienergy,50),thavlab(ienergy,50),rstep(50)
     &,dsr(ienergy,50),tloss(ienergy),angt(181),angd(181)
      dimension sig(ienergy,181),sigtab(200,91),sigt(ienergy),D(50)
     &,aint(800),the(181),sofe(200),ei(200),dsde(200),tke(ienergy)
     &,rsti(50),sigri(50),theti(50),sigrint(800),thetint(800),rint(800)
     &,irstmi(ienergy),irstma(ienergy),thetmi(ienergy),thetma(ienergy)
      character*5 intp
      character*30 addtext
c
      degorad = 180.d0/pi
      wid = widang/degorad
      lab = ierr
c
      thetmax = -pi
      thetmin = pi
      print*,' maxm,minm,irstmax = ',maxm,minm,irstmax
      do mm=minm,maxm
       irstmi(mm) = 1
       irstma(mm) = 1
       thetma(mm) = -pi
       thetmi(mm) = pi
      enddo
c  For each energy (mm) we determine the minimum and maximum angle
      do mm=minm,maxm
       thetlast = 0.d0
       thetnextlast = 0.d0
       do i = 1,irstmax
        if(lab.eq.0)theto = thav(mm,i)
        if(lab.eq.1)theto = thavlab(mm,i)
c        if(i.eq.2.and.theto.eq.0.d0)then 
c         thetmi(mm) = pi
c         thetlast = 0.d0
c        endif
        if(theto.eq.0.d0)then
         if(thetlast.ne.0.d0.and.thetnextlast.ne.0.d0)then
          irstma(mm) = i-1
          goto 42
         endif
        elseif(theto.ne.0.d0)then
         if(theto.gt.thetmax)thetma(mm)=theto
         if(theto.lt.thetmin)thetmi(mm)=theto
         if(thetlast.eq.0.d0)irstmi(mm) = i
         if(i.eq.irstmax)irstma(mm) = irstmax
        endif
        thetnextlast = thetlast
        thetlast = theto
       enddo
42     continue
       if(irstma(mm).gt.irstmi(mm)+1)then
        if(thetma(mm).gt.thetmax)thetmax=thetma(mm)
        if(thetmi(mm).lt.thetmin)thetmin=thetmi(mm)
       endif
       thettmdg = thetmi(mm)*degorad
       thettxdg = thetma(mm)*degorad
       write(*,'(a,3i5,2f10.3)')
     &' mm,irstmi,irstma,thetmi(mm),thetma(mm), = '
     &,mm,irstmi(mm),irstma(mm),thettmdg,thettxdg
      enddo
25    if(thetmin.lt.-pi)thetmin = -pi
c  Now prepare interpolation in angles and cross sections
c   Precaution for negative angles
      if(thetmin.lt.0.d0)then
       lil = 1
       tdgmin = dble(int(thetmin*degorad/10.d0-1.d0)*10)
       if(abs(thetmin+pi).lt.accur)tdgmin = -180.d0
      else
       lil = 0       
       tdgmin = dble(int(thetmin*degorad/10.d0+1.d0)*10)
      endif
      tdgmax = dble(int(thetmax*degorad/10.d0+2.d0)*10)
      if(tdgmax.gt.179.999d0)tdgmax = 180.d0
      deltdg = 1.d0      
27    nre = int((tdgmax-tdgmin)/deltdg+1.001)
      if(nre.gt.361)then
       deltdg = 5.d0*deltdg
       goto 27
      elseif(nre.gt.181)then
       deltdg = 2.d0*deltdg
       goto 27
      endif
      tdgmax = tdgmin + dble(nre-1)*deltdg
      if(lil.eq.1)then
       n0 = int(-tdgmin/deltdg+1.0001d0)
       nnmax = nre-n0+1
      elseif(lil.eq.0)then
       nnmax = nre
      endif
      delt = deltdg/degorad
      tmin = tdgmin/degorad
      tmax = tdgmax/degorad
c      write(*,*)' tmin,tmax = ',tdgmin,tdgmax
      dela=wid/2.35482d0
      factor = 1.d0/sqrt(6.28318530717958647d0)/dela
c  input for interpolation for each energy (mm)
      do mm=minm,maxm
       kmax = irstma(mm)-irstmi(mm)+1
       do kk=1,800
        thetint(kk) = 0.d0
        sigrint(kk) = 0.d0
        rint(kk) = 0.d0
       enddo
c  drop data if only two points
       if(kmax.le.2)goto 38
       do k=1,kmax
        theti(k) = 0.d0
        sigri(k) = 0.d0
        rsti(k) = 0.d0
       enddo
       kk = 0
        print *,' kmax : ',kmax
       do k=1,kmax
        i = irstmi(mm)+k-1
        rsti(k) = rstep(i)
c  temporary angle and cross section to be interpolated
        if(lab.eq.0)theti(k) = thav(mm,i)
        if(lab.eq.1)theti(k) = thavlab(mm,i)
        sigri(k) = dsr(mm,i)
c      if(tloss(mm).gt.45.d0.and.tloss(mm).lt.70.d0.and.k.le.5)then
c       if(k.eq.1)write(*,'(a,f10.5)')' eloss = ',tloss(mm)
c       thd = theti(k)*degorad
c       write(*,'(a,i3,2f10.3,1pe15.3)')
c     &' k,rstep,thet,sig = ',k,rsti(k),thd,sigri(k)
c      endif
        mult = min(int(abs(theti(k+1)-theti(k))*degorad+1.d0),100)
        dlr = (rsti(k+1)-rsti(k))/float(mult)
cnan
        print *,' mult,dlr : ',mult,dlr
cnan
        do kkk = 1,mult
         kk = kk+1
         rint(kk) = rsti(k)+float(kkk)*dlr
cnan      print*,' kk,rint(kk) = ',kk,rint(kk)
        enddo
       enddo
       kkmax = kk
       if(kkmax.gt.800)then
        write(*,'(a,i5,a,f10.3)')' Error: kk = ',kkmax,
     &'.gt.800 for eloss = ',tloss(mm),' in wilsh_s'
        ierr = 1
        return
       endif
       ifail = 1
       call E01BEF(kmax,rsti,theti,D,ifail)
       ifail = 1
       call E01BFF(kmax,rsti,theti,D,kkmax,rint,thetint,ifail)
       ifail = 1
       call E01BEF(kmax,rsti,sigri,D,ifail)
       ifail = 1
       call E01BFF(kmax,rsti,sigri,D,kkmax,rint,sigrint,ifail)
       do ii = 1,nre
        the(ii) = tmin + dble(ii-1)*delt
        do kk = 1,kkmax
         aint(kk)=abs(sigrint(kk))
     &*exp(-((the(ii)-thetint(kk))/dela)**2/2.d0)
        enddo
        ifail = 1
        call D01GAF(rint,aint,kkmax,ref,err,ifail)
        if(ifail.ne.0)then
         write(*,*)' trouble in wilsh_s, ifail = ',ifail
         write(*,*)' kkmax = ',kkmax
         ref = 0.d0
         stop
        endif    
        sig(mm,ii)=ref*factor
c      thee = 180.d0/pi*the(ii)
c      write(*,'(a,2f10.5,1pe12.4)')
c     &' eloss,the,sig = ',tloss(mm),thee,sig(mm,ii)
       enddo
38     continue
      enddo
      tkemin = dble(int((ecm-tloss(maxm))/10.d0)*10)
      tkemax = dble(int(ecm-tloss(minm))+3)
      dele = 1.d0
29    me = int((tkemax-tkemin)/dele+1.001d0)
      if(me.gt.400)then
       dele = 5.d0*dele
       goto 29
      elseif(me.gt.200)then
       dele = 2.d0*dele
       goto 29
      endif
      mmm = maxm-minm+1
      do m=1,mmm
       tke(m)= ecm-tloss(maxm-m+1)
c      print*,'tke(m) = ', tke(m)
       sigt(m) = 0.d0
      enddo
      do m = 1,me
       do nn = 1,nnmax
        sigtab(m,nn) = 0.d0
       enddo
      enddo
      do n=1,nre
       do m=1,mmm 
        mm=maxm-m+1
        sigt(m)=sig(mm,n)
       enddo
       ierr = 0
       call interpolateder_s(tke,sigt,mmm,tkemin,tkemax,dele,ei
     &,sofe,dsde,me,ierr)
       if(ierr.eq.1)goto 2
c
c  cross section in dsig/dO in mb/sr or dsig/dthet in mb/rad
       thee = tdgmin + dble(n-1)*deltdg
       theee = thee/degorad
       if(la.le.1)then
        sia = abs(sin(theee))
        if(sia.lt.accur)sia = accur**2
        ffsig = 1.d0/sia
       elseif(la.ge.2)then
        ffsig = 2.d0*pi
       endif
       if(lil.eq.1)then
        if(n.le.n0)then
         nn = n0-n+1
        elseif(n.gt.n0)then
         nn = n-n0+1
         angt(nn) = thee
        endif
       else
        nn = n
        angt(nn) = thee
       endif
       angd(nn) = 0.d0
       do m=1,me
        ss = sofe(m)*ffsig
        sigtab(m,nn)=sigtab(m,nn)+ss
        if(ss.lt.1.d-35)sigtab(m,nn)=0.d0
        if(ei(m).lt.ecm-ecut)then
         angd(nn) = angd(nn) + ss*dele
        endif
       enddo
      enddo
      if(lil.eq.1)ttmin = 0.d0
      if(lil.eq.0)ttmin = tdgmin
      iunit=nfile
      write(iunit,'(6e12.4)')ttmin,tdgmax,deltdg,tkemin,tkemax,dele
      do m=1,me
       write(iunit,'(1p7e11.3)')(sigtab(m,nn),nn=1,nnmax)
      enddo
      ierr = 0
      close(nfile)
      return
2     end
