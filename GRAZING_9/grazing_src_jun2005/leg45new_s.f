c      Subroutine LEG 45 for program grazing
c      At exit ierr = 1 for error
c              ierr = 2 for warning
c              ierr = 3 for differential cross section
c              ierr = 4 for total cross section and new value of r
c              ierr = 0 for end or new value of energy
c---------------------------------------------------------------------
      subroutine leg45_s(nestep,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
c      list,ienergy,irangen and irangep are defined in grazing.icl
      common/za/za,aa,zaa,aaa
      common/ecm/ecm
      common/capt/eb,rbar,i5,intp
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      common/iza/iza,iaa,izaa,iaaa
      common/q/qns,qnp,qps,qpp
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/eav/eav,eavns,eavnp,eavps,eavpp,eavtr,
     &eav0,eava0,eavaa0
      common/mav/avma,avmaa,avmans,avmanp,avmaps,avmapp,
     &avmaans,avmaanp,avmaaps,avmaapp,avmatr,avmaatr,
     &avma0,avmaa0,avm0
      common/wid/wid2e,wid2ns,wid2np,wid2ps,wid2pp,wid2tr,wid2ea,
     &wid2eaa,wid2e0,wid2ea0,wid2eaa0
      common/inel/p0ina,p0inaa,eavina,eavinaa,wid2ina,wid2inaa,
     &avmina,avminaa,eava00,eavaa00,wid2a00,wid2aa00,inelast
      common/nzdis/ p0nn(irangen),p0zz(irangep),
     &eavdist(irangen,irangep),eavdist00,avmadist(irangen,irangep),
     &avmadist00,avmdist(irangen,irangep),avmdist00,
     &wid2dist(irangen,irangep),wid2dist00
      common/char/elemp,elemt
      common/head/head
      common/leg/ecminit,rainbcdg,decm,cs,r,delr,al,theta,
     &cscoul,cslab,thet,thetlabdeg,thetcou,ppcapt,ww,pr0,
     &drifta,driftz,wida,widz,eava,eavaa,tavs,widea,wideaa,avmin,avm,
     &changer,changers,ppscat,totmerg,dthet,
     &deltae,sig0drut,sigreac0,sigreac0lab,cross,alorho,rho,
     &istep,iout,ipopn,norbit,irstep,lab
      common/leg2/sigcaptn
      common/leg3/delintloss,rcapt,alcapt
      common/tot/in00,ip00,i0nmin,i0nmax,i0zmin,i0zmax
      common/leg4/dsigz(irangep),dsigzlab(irangep),
     &ppnz(irangen,irangep),iselectn(list),iselectz(list),
     &dsig(irangen,irangep),dsiglab(irangen,irangep),
     &inmin,inmax,izmin,izmax,leg4
      common/nztdis/eav00,wid200
      dimension sig(irangen,irangep)
      dimension gamd(100),tke(ienergy)
      dimension ppa(irangea),tdista(irangea,ienergy)
     &,sigtdista(irangea,ienergy),siga(irangea),sumsiga(irangea)
     &,sigtdistz(irangep,ienergy),sigz(irangep),sumsigz(irangep)
     &,ppz(irangep),tdistz(irangep,ienergy),dista(irangea)
     &,distz(irangep),alsigtdista(irangea,ienergy)
     &,alsigtdistz(irangep,ienergy)
      common/vers/version
      character*5 intp
      character*5 version
      character*60 addtext
      character*40 efilename,efilename1
      character*16 head
c
      iprep = ierr
      sumdista = 0.d0
      sumdistz = 0.d0
      do n=1,irangea
       do i=1,ienergy
        tdista(n,i)=0.d0
       enddo
      enddo
      do n=1,irangep
       ppz(n) = 0.d0
       do i=1,ienergy
        tdistz(n,i)=0.d0
       enddo
      enddo
      if(istep.ne.1)goto 79
      if(nestep.ne.1)goto 77
320   if(intp.eq.'inter')write(*,'(4(/a))')
     &' You must calculate the energy distributions',
     &'  with a finite energy resolution (FWHM).',
     &' Enter width in MeV',
     &'  (if you enter a number <1 MeV, width is set equal to 1 MeV)'
      read(i5,*,err=1101)tres
      icheck = 0
      if(tres.lt.1.)tres = 1.d0
      goto 1111
1101  write(*,*)' you must provide a number'
      icheck = icheck + 1
      if(icheck.gt.5)goto 2
      goto 320
1111  addtext = '    !energy resolution in MeV'
      if(intp.eq.'inter')write(10,'(f7.4,a)')tres,addtext
      deltae = tres/2.35482d0
      dela = deltae
c   Preparation for total cross sections
77    i0amin = -2
      i0amax = 2
      sigdela = 0.d0
      sigdelz = 0.d0
      sigcha = 0.d0
      sigchz = 0.d0
      cross = 0.d0
      crossres = 0.d0
      do n=1,irangen
       do m=1,irangep
        sig(n,m)=0.d+0
       enddo
      enddo
      do n=1,irangea
       siga(n) = 0.d0
       sumsiga(n) = 0.d0
       do i=1,ienergy
        sigtdista(n,i)=0.d0
       enddo
      enddo
      do n=1,irangep
       sigz(n) = 0.d0
       sumsigz(n) = 0.d0
       do i=1,ienergy
        sigtdistz(n,i)=0.d0
       enddo
      enddo
      tkemax = dble(10*int(.1d0*(1.1d0*ecm+1.d0)))
      tkemin = dble(10*int(0.065d0*eb))
      delt = (tkemax-tkemin)/float(ienergy-1)
      tminl = ecm-tkemax
      tmaxl = ecm-tkemin
      do i = 1,ienergy
       tke(i) = tkemin+float(i)*delt
      enddo
      nnum = (i0nmax-i0nmin)*(i0zmax-i0zmin)
      acc = accur/dble(nnum)
c  TKE versus A  and Z 
79    ia00 = in00+ip00
      iamax = 1
      factor = 1.d0/sqrt(6.28318530717958647d0)/dela
      do isa = inmin+izmin,inmax+izmax
       na = ia00 + isa
       ppa(na)=0.d0
       dista(na) = 0.d0
       do isn = inmin,inmax
        in = in00+isn
        isz = isa-isn
        iz = ip00+isz
        if(isz.lt.izmin.or.isz.gt.izmax)goto 453
        pp = ppnz(in,iz)
        ppp = pp*ppscat
        if(ppp.lt.acc)then
         goto 453
        else
         ppa(na) = ppa(na)+ ppp
         ppz(iz) = ppz(iz)+ ppp
         if(isn.eq.0.and.isz.eq.0)then
          eea = eav00
          wwid = wid200
         else
          eea = eavdist(in,iz)+eavina+eavinaa
          wwid = wid2dist(in,iz)+wid2ina+wid2inaa+dela**2
         endif
         if(eea.gt.0.d0.and.wwid.gt.0.d0)then
          eqns = qns
          if(isn.gt.0)eqns=-qnp
          eqps = qps
          if(isz.gt.0)eqps=-qpp
          eminloss = tminl-isn*eqns-isz*eqps
          e0 = eminloss-delt
          emaxloss = tmaxl-isn*eqns-isz*eqps
          ierr = 0
          call probdtab_s(eea,wwid,e0,delt,gamd,accur,ierr)
          if(ierr.eq.1)then
           write(*,*)' Error in leg45_s for isn,isz = ',isn,isz
           return
          endif
          if(isn.eq.0.and.isz.eq.0)then
           do i=1,ienergy
            if(tke(i).gt.ecm)then
             emm = ecm-tke(i-1)
             ierr= 0
             call prob_s(eea,wwid,emm,probdel,ierr)
             if(ierr.eq.1)goto 2
             probdel = (1.d0 - probdel - .5d0*gamd(ienergy-i+1)*delt)*pp
             isel = i
             goto 11
            endif
           enddo
11         continue
          endif
          do i=1,ienergy
           tdista(na,i)=tdista(na,i)+gamd(ienergy-i+1)*pp
           tdistz(iz,i)=tdistz(iz,i)+gamd(ienergy-i+1)*pp
          enddo
         else
          goto 453
         endif
        endif
453     continue
       enddo
       if(ppa(na).lt.acc.and.isa.le.-2)iamin = isa
       if(ppa(na).gt.acc.and.isa.ge.2)iamax = isa
452   continue
      enddo
      if(iamin.lt.i0amin)i0amin = iamin
      if(iamax.gt.i0amax)i0amax = iamax
c
      do isa = iamin,iamax
       na = isa+ia00
       do i=1,ienergy
        dista(na) = dista(na) + tdista(na,i)*delt      
       enddo
       if(isa.eq.0)dista(na)=dista(na)+probdel
       iach = iaa + isa
       if(dista(na).lt.1.d30.and.dista(na).gt.0.d0)then
        distan = dista(na)
       else
        write(*,'(a,i3,a/a,f10.4)')
     &' Error: dista strange for isa = ',isa,
     &' in leg45_s',' for r = ',r
       endif 
       ppan = ppa(na)/ppscat
c       if(isa.eq.0)print*,' p0 = ',p0
c       write(*,'(a,i4,2(a,f8.6))')
c     &' A = ',iach,' dista(A) = ', distan,' ppa(A) = ',ppan
       sumdista = sumdista + dista(na)
      enddo
      do isz = izmin,izmax
       nz = isz + ip00
       if(dista(na).lt.1.d30.and.dista(na).gt.0.d0)then
        distz(nz) = 0.d0
       else
        write(*,'(a,i3,a/a,f10.4)')
     &' Error: distz strange for isz =',isz,
     &' in leg45_s',' for r = ',r
       endif 
       do i=1,ienergy
         distz(nz) = distz(nz) + tdistz(nz,i)*delt      
       enddo
       if(isz.eq.0)distz(nz)=distz(nz)+probdel
       sumdistz = sumdistz + distz(nz)
      enddo
      sumdista = sumdista+p0
      sumdistz = sumdistz + p0
c      print*,' sumdista = ',sumdista,' sumdistz = ',sumdistz
c

c  integration of cross sections (not accurate for irstep.eq.0)
27    if(istep.eq.1.and.totmerg.lt.0.d0)then
       alpr = -1.d-20
      elseif(istep.eq.1.and.totmerg.ge.0.d0)then
       alpr = alorho*sqrt(totmerg/10.d0/pi)
      endif
      if(irstep.eq.0)then
       factor=2.d0*pi*10.d0*al*(al-alpr)/alorho**2
       if(alpr.lt.0.d0)then
        factor = 0.d0
       elseif(alpr.lt.accur)then
        factor = factor/2.d0
       endif
       alpr = al
      elseif(irstep.eq.1)then
       simp1 = delr/3.d0
       factor=2.d0*pi*10.d0*al*(al-alpr)/alorho**2
     &       +2.d0*pi*10.d0*r*(ecm-za*zaa*esq/2.d0/r
     &-unuc(r,ierr1)-r/2.d0*dunuc(r,ierr2))/ecm*simp1
      else
       ifact = 2*(1+mod(irstep-1,2))
       simp = dble(ifact)/3.d0*delr
       if(mod(irstep-1,10).eq.0)then
        simp = delr
        delr = delr*2.d0
       endif
       ierr1 = 1
       ierr2 = 1
       factor=2.d0*pi*10.d0*r*(ecm-za*zaa*esq/2.d0/r
     &-unuc(r,ierr1)-r/2.d0*dunuc(r,ierr2))/ecm*simp
       if(ierr1.eq.1.or.ierr2.eq.1)then
        ierr = 1
        goto 2
       endif
      endif
      cross = cross + factor*(1.d0-p0)*ppscat
      sigcaptn = sigcaptn + factor*(1.d0-ppscat)
      do in=inmin,inmax
       do im=izmin,izmax
        n=in00+in
        m=ip00+im
        sig(n,m)= sig(n,m)+factor*ppnz(n,m)*ppscat
       enddo
      enddo
      do ia = iamin,iamax
       na = ia+ia00
       do i = 1,ienergy
        sigtdista(na,i) = sigtdista(na,i) + 
     &tdista(na,i)*factor*ppscat
       enddo
       siga(na) = siga(na)+factor*ppa(na)
      enddo
      sigdela = sigdela + probdel*factor*ppscat
      do iz = izmin,izmax
       nz = ip00 + iz
       do i = 1,ienergy
        sigtdistz(nz,i) = sigtdistz(nz,i) + 
     &tdistz(nz,i)*factor*ppscat
       enddo
       sigz(nz) = sigz(nz)+factor*ppz(nz)
      enddo
      sigdelz = sigdelz + probdel*factor*ppscat
      if(irstep.ne.0)irstep = irstep+1
      istep = istep+1
      if(iprep.eq.1.and.istep.ge.3)goto 433
c      write(*,'(4(a,i4))')' iamin = ',iamin,' iamax = ',iamax,
c     &' izmin = ',izmin,' izmax = ',izmax
c      write(*,'(4(a,i4))')' i0amin = ',i0amin,' i0amax = ',i0amax,
c     &' i0zmin = ',i0zmin,' i0zmax = ',i0zmax
      if(abs(1.-p0).lt.accur)goto 433
      alpr = al
      ierr = 4
      goto 2
c
c4  make table for contour plotting of TKE versus A
433   do isa = i0amin,i0amax
       na = isa+ia00
       do i = 1,ienergy
        if(sigtdista(na,i).le.0.d0)then
         alsigtdista(na,i)= -10.d0
        else
         alsigtdista(na,i) = log10(sigtdista(na,i))
        endif
        sumsiga(na)=sumsiga(na)+sigtdista(na,i)*delt
       enddo
       if(isa.eq.0)sumsiga(na)=sumsiga(na)+sigdela
      enddo
      do isz = i0zmin,i0zmax
       nz = isz+ip00
       do i = 1,ienergy
        if(sigtdistz(nz,i).le.0.d0)then
         alsigtdistz(nz,i) = -10.d0
        else 
         alsigtdistz(nz,i) = log10(sigtdistz(nz,i))
        endif
        sumsigz(nz)=sumsigz(nz)+sigtdistz(nz,i)*delt
       enddo
       if(isz.eq.0)sumsigz(nz)=sumsigz(nz)+sigdelz
      enddo
      if(version.eq.'vms  '.or.version.eq.'unix* ')then
       efilename1 = 'graz_'//head//'_tdista'
       call condens_s(efilename1,efilename,le)
       write(*,'(a/a,a31,a)')
     &' Table of log_10(sig) as function of TKE versus A'
     &,' is stored in ',efilename,' Coulomb barrier in FOR018'
       write(98,'(a/a,a31,a)')
     &' Table of log_10(sig) as function of TKE versus A'
     &,' is stored in ',efilename,' Coulomb barrier in FOR018'
       open (unit=93,file=efilename,status='unknown')
      elseif(version.eq.'unix ')then
       write(*,'(a/a)')
     &' Table of log_10(sig) as function of TKE versus A'
     &,' is stored in FOR093. Coulomb barrier in FOR018' 
       write(98,'(a/a)')
     &' Table of log_10(sig) as function of TKE versus A'
     &,' is stored in FOR093. Coulomb barrier in FOR018'
      endif
      tkemax = ecm - tminl
      tkemin = ecm - tmaxl
      eiamin = dble(i0amin+iaa)
      eiamax = dble(i0amax+iaa)
      write(93,'(1x,6e12.4)')eiamin,eiamax,1.d0,tkemin,tkemax,delt
      do i = 1,ienergy
       if(alsigtdista(isa+ia00,i).gt.10.d0.or
     &  .alsigtdista(isa+ia00,i).lt.-10.d0)
     &    print*,' sigtdista =  ',sigtdista(isa+ia00,i)
       write(93,'(1x,1p8e10.2)')
     &(alsigtdista(isa+ia00,i),isa=i0amin,i0amax)
      enddo
      write(*,'(a)')'     A     sig(mb)   sumsig(mb)'
      write(98,'(a)')'     A     sig(mb)   sumsig(mb)'
      do isa = i0amin,i0amax
       na = isa +ia00
       ipa = iaa+isa
       write(*,'(i5,3x,2e12.5)')ipa,siga(na),sumsiga(na)
       write(98,'(i5,3x,2e12.5)')ipa,siga(na),sumsiga(na)
       avar = dble(ipa)
       zvar = .487d0*avar/(1.d0+avar**.66666d0/166.d0)      
       ebarrier = zvar*(zaa+za-zvar)*esq/rcapt*(1.d0-1.d0/(ainv*rcapt))
       sigcha = sigcha + sumsiga(na)     
      enddo
c  make table for contour plotting of TKE versus Z
      if(version.eq.'vms  '.or.version.eq.'unix* ')then
       efilename1 = 'graz_'//head//'_tdistz'
       call condens_s(efilename1,efilename,le)
       write(*,'(a/a,a31,a)')
     &' Table of log_10(sig) as function of TKE versus Z'
     &,' is stored in ',efilename,' Coulomb barrier in FOR017'
       write(98,'(a/a,a31,a)')
     &' Table of log_10(sig) as function of TKE versus Z'
     &,' is stored in ',efilename,' Coulomb barrier in FOR017'
       open (unit=94,file=efilename,status='unknown')
      elseif(version.eq.'unix ')then
       write(*,'(a/a)')
     &' Table of log_10(sig) as function of TKE versus Z'
     &,' is stored in FOR094. Coulomb barrier in FOR017'
       write(98,'(a/a)')
     &' Table of log_10(sig) as function of TKE versus Z'
     &,' is stored in FOR094. Coulomb barrier in FOR017'
      endif
      tkemax = ecm - tminl
      tkemin = ecm - tmaxl
      eizmin = dble(i0zmin+iza)
      eizmax = dble(i0zmax+iza)
      write(94,'(1x,6e12.4)')eizmin,eizmax,1.d0,tkemin,tkemax,delt
      do i = 1,ienergy
       write(94,'(1x,1p8e10.2)')
     &(alsigtdistz(isz+ip00,i),isz=i0zmin,i0zmax)
      enddo
      write(*,'(a)')'     Z     sig(mb)   sumsig(mb)'
      write(98,'(a)')'     Z     sig(mb)   sumsig(mb)'
      do isz = i0zmin,i0zmax
       nz = isz +ip00
       ipz = iza+isz
       if(sigz(nz).ne.0.d0)then
        write(*,'(i5,3x,2e12.5)')ipz,sigz(nz),sumsigz(nz)
        write(98,'(i5,3x,2e12.5)')ipz,sigz(nz),sumsigz(nz)
        ebarrier = dble(ipz*(izaa+iza-ipz))*esq/rcapt
     &*(1.d0-1.d0/(ainv*rcapt))
       endif
       sigchz = sigchz + sumsigz(nz)     
      enddo
      write(*,'(a,f13.6,a,f13.6)')
     &' Total cross section (sum Z) = ',sigchz,' (sum A)= ',sigcha
      close(unit=93)
      write(*,'(/a,e12.4,a)')
     &' Total grazing cross section = ',cross,'  mbarn'
      write(98,'(/a,e12.4,a)')
     &' Total grazing cross section = ',cross,'  mbarn'
      if(tres.gt.accur)then
       write(*,'(/a,e12.4,a)')
     &' "Total grazing cross section" = ',crossres,'  mbarn'
       write(*,'(a,f6.2,a/a,f6.1,a)')
     &'       with energy resolution(FWHM) = ',tres,' MeV',
     &'       for scattering angles larger than ',scatmin,' degrees'
       write(98,'(/a,e12.4,a)')
     &' "Total grazing cross section" = ',crossres,'  mbarn'
       write(98,'(a,f6.2,a/a,f6.1,a)')
     &'       with energy resolution(FWHM) = ',tres,' MeV',
     &'       for scattering angles larger than ',scatmin,' degrees'
      endif
      sigtot = sigcaptn + cross
      sigtotlim = sigcaptn + crossres
      write(*,'(a,e12.4,a)')
     &' Total capture cross section = ',sigcaptn,' mbarn'
      write(98,'(a,e12.4,a)')
     &' Total capture cross section = ',sigcaptn,' mbarn'
      write(*,'(a,e12.4,a)')
     &' Total reaction cross section = ', sigtot,' mbarn'
      write(98,'(a,e12.4,a)')
     &' Total reaction cross section = ', sigtot,' mbarn'
c
      ierr = 0
2     return
      end





