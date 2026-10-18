c      Subroutine LEG 41 for program grazing
c      At exit ierr = 1 for error
c              ierr = 2 for warning
c              ierr = 3 for differential cross section
c              ierr = 4 for total cross section and new value of r
c              ierr = 0 for end or new value of energy
c---------------------------------------------------------------------
      subroutine leg41_s(nestep,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
c      include 'common.fil'
c      list,ienergy,irangen and irangep are defined in grazing.icl
      common/za/za,aa,zaa,aaa
      common/ecm/ecm
      common/capt/eb,rbar,i5,intp
c      common:akyuz
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
c      common:qvaletc
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
     &wid2dist(irangen,irangep),wid2dist00,wid2adist(irangen,irangep)
      common/char/elemp,elemt
      common/head/head
      common/leg/ecminit,rainbcdg,decm,cs,r,delr,al,theta,
     &cscoul,cslab,thet,thetlabdeg,thetcou,ppcapt,ww,pr0,
     &drifta,driftz,wida,widz,eava,eavaa,tavs,widea,wideaa,avmin,avm,
     &changer,changers,ppscat,totmerg,dthet,
     &deltae,sig0drut,sigreac0,sigreac0lab,cross,alorho,rho,
     &istep,iout,ipopn,norbit,irstep,lab
      common/tot/in00,ip00,i0nmin,i0nmax,i0zmin,i0zmax
      common/leg2/sigcaptn
      common/leg4/dsigz(irangep),dsigzlab(irangep),
     &ppnz(irangen,irangep),iselectn(list),iselectz(list),
     &dsig(irangen,irangep),dsiglab(irangen,irangep),
     &inmin,inmax,izmin,izmax,leg4
      dimension sig(irangen,irangep)
      dimension pnzcor(irangen,irangep),sigcor(irangen,irangep)
      dimension amass(irangen,irangep),amassapp(irangen,irangep)
     &,aamass(irangen,irangep),aamassapp(irangen,irangep)
     &,evapdis(irangen,irangep,ienergy)
      character*5 intp
      character*16 head
      character*5 elemp
      character*5 elemt
      character*40 addtext
c
      iprep = ierr
      if(istep.ne.1)goto 51
      if(nestep.ne.1)goto 15
       write(*,'(/a/a/a,f7.3,a/a/a)')
     &' You have decided to calculate total cross sections for each N',
     &' and Z integrated over L, but corrected for actual mass ',
     &' Stepwidth in r is delr =',delr,' fm',
     &' Standard output is  sig_reac, sig_merg',
     &' as well as total cross sections sig_NZ.'
48    if(intp.eq.'inter')write(*,'(a/a)')
     &' Do you want energy distribution of ejectiles or recoils?',
     &' Enter: 1 for ejectiles, 2 for recoils' 
      read(i5,*,err=1100)ievap
      icheck = 0
      goto 1110
1100  write(*,*)' you must provide a number'
      icheck = icheck + 1
      if(icheck.gt.5)goto 2
      goto 48
1110  addtext = '  ! 1 = ejectile, 2 = recoil'
      if(ievap.ne.1.and.ievap.ne.2)goto 48
      if(intp.eq.'inter')write(10,'(i5,a40)')ievap,addtext
      write(*,'(/a)')
     &' You are calculating approximate energy distributions'
      if(ievap.eq.1)write(*,'(a)')
     &'    of the ejectiles after correction for actual masses'
      if(ievap.eq.2)write(*,'(a)')
     &'    of the recoils after correction for actual masses'
c
15    inmi = inmin-2
      inma = inmax+2
      izmi = izmin-2
      izma = izmax+2
      ierr = 1
      call emass_s(amass,amassapp,in00,ip00,inmi,inma,izmi,izma,ierr)
      if(ierr.eq.1)then
       goto 2
      elseif(ierr.eq.2)then
       izmin = izmi+2
       i0zmin = izmin
      endif
      ierr = 2
      call emass_s(aamass,aamassapp,in00,ip00,inmi,inma,izmi,izma
     &,ierr)
      if(ierr.eq.1)then
       goto 2
      elseif(ierr.eq.2)then
       izmax = izma-2
       i0zmax = izmax
      endif
c  preparation for integrating cross section
       cross = 0.d0
       sigtar = 0.d0
       sigproj = 0.d0
       do n=1,irangen
        do m=1,irangep
         sig(n,m)=0.d+0
        enddo
       enddo
       do n=1,irangen
        do m=1,irangep
         sigcor(n,m) = 0.d0
        enddo
       enddo
c
51    continue
       iout = 0
      ierr = ievap
       if(inmin.lt.i0nmin)inmin=i0nmin
       if(inmax.gt.i0nmax)inmax=i0nmax
       if(izmin.lt.i0zmin)izmin=i0zmin
       if(izmax.gt.i0zmax)izmax=i0zmax
      if(istep.eq.1)dele = 0.d0
c      print*,' in00,ip00, = ',in00,ip00
c      print*,' inmin,inmax,izmin,izmax =',inmin,inmax,izmin,izmax
      call enzcor_s(amass,aamass,amassapp,aamassapp,in00,ip00
     &,inmin,inmax,izmin,izmax,pnzcor,dele,evapdis,iout,ierr)
      if(ierr.eq.1)goto 2
      if(istep.eq.1.and.totmerg.lt.0.d0)then
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
      sigtar = sigtar + factor*(p0a-p0)*ppscat
      sigproj = sigproj + factor*(p0aa-p0)*ppscat
      do i=inmin,inmax
       do j=izmin,izmax
        n=in00+i
        m=ip00+j
        sigcor(n,m)= sigcor(n,m)+factor*pnzcor(n,m)*ppscat
       enddo
      enddo
      do i=inmin,inmax
       do j=izmin,izmax
        n=in00+i
        m=ip00+j
        sig(n,m)= sig(n,m)+factor*ppnz(n,m)*ppscat
       enddo
      enddo
      if(irstep.ne.0)irstep = irstep+1
      istep = istep+1
c  Check for finished integration
      if(iprep.eq.1.and.istep.ge.3)goto 33
      if(abs(1.-p0).lt.accur)goto 33
      alpr = al
      ierr = 4
      goto 2
c
c  printout of total cross sections 
33    write(*,*)'-----------------------------------------------'
      write(98,*)'-----------------------------------------------'
      write(*,'(a,f10.3,a)')
     &' Total reaction cross section = ',cross,' mbarn'      
      write(98,'(a,f10.3,a)')
     &' Total reaction cross section = ',cross,' mbarn'      
      sigctot = 0.d0
      sigutot = 0.d0 
      if(ievap.eq.1)then
       write(*,'(a)')
     &' Total cross sections for ejectiles after masscorrection in mb'      
       write(98,'(a)')
     &' Total cross sections for ejectiles after masscorrection in mb'      
       do j=i0zmin,i0zmax
        m=ip00+j
        sigcort = 0.d0
        sigt = 0.d0
        nnz = iza + j
        write(*,'(a,i3)')' For Z = ',nnz
        write(98,'(a,i3)')' For Z = ',nnz
        do i=i0nmin-2,i0nmax
         n=in00+i
         naa = iaa + i + j
         write(*,'(a,i3,2(a,f10.3))')' A = ',naa,
     &' Uncorrected: sig = ',sig(n,m),'  sigcor = ',sigcor(n,m)
         write(98,'(a,i3,2(a,f10.3))')' A = ',naa,
     &' Uncorrected: sig = ',sig(n,m),'  sigcor = ',sigcor(n,m)
         if(i.eq.0.and.j.eq.0.and.inelast.ge.0)then
          write(*,'(a,f7.3,a)')
     &' projectile elastic (but target excit.) sig_el = ',sigtar,' mb'
          write(98,'(a,f7.3,a)')
     &' projectile elastic (but target excit.) sig_el = ',sigtar,' mb'
         endif
         sigt = sigt + sig(n,m)
         sigcort = sigcort+sigcor(n,m)
        enddo
        write(*,'(a,i3,2(a,f10.3))')
     &' Total for Z = ',nnz,' sig = ',sigt,' sigcor = ',sigcort
        write(98,'(a,i3,2(a,f10.3))')
     &' Total for Z = ',nnz,' sig = ',sigt,' sigcor = ',sigcort
        sigctot = sigctot + sigcort
        sigutot = sigutot + sigt
       enddo
       write(*,'(2(a,f10.3,a/))')
     &' Grand total: sigcor = ',sigctot,' mbarn',
     &'  projectile elastic (but target excit.) = ',sigtar,' mbarn'
       write(98,'(2(a,f10.3,a/))')
     &' Grand total: sigcor = ',sigctot,' mbarn',
     &'  projectile elastic (but target excit.) = ',sigtar,' mbarn'
      elseif(ievap.eq.2)then
       write(*,'(a)')
     &' Total cross sections for recoils after masscorrection in mbarn'      
       write(98,'(a)')
     &' Total cross sections for recoils after masscorrection in mbarn'      
       do j=i0zmin,i0zmax
        m=ip00+j
        sigcort = 0.d0
        sigt = 0.d0
        nnz = izaa - j
        write(*,'(a,i3)')' For Z = ',nnz
        write(98,'(a,i3)')' For Z = ',nnz
        do i=i0nmin,i0nmax+2
         n=in00+i
         naa = iaaa - i - j
         write(*,'(a,i3,2(a,f10.3))')' A = ',naa,
     &' Uncorrected: sig = ',sig(n,m),'  sigcor = ',sigcor(n,m)
         write(98,'(a,i3,2(a,f10.3))')' A = ',naa,
     &' Uncorrected: sig = ',sig(n,m),'  sigcor = ',sigcor(n,m)
         if(i.eq.0.and.j.eq.0.and.inelast.ge.0)then
          write(*,'(a,f7.3,a)')
     &' target elastic (but projectile excit.) sig_el = ',sigproj,' mb'
          write(98,'(a,f7.3,a)')
     &' target elastic (but projectile excit.) sig_el = ',sigproj,' mb'
         endif
         sigt = sigt + sig(n,m)
         sigcort = sigcort+sigcor(n,m)
        enddo
        write(*,'(a,i3,2(a,f10.3))')
     &' Total for Z = ',nnz,' sig = ',sigt,' sigcor = ',sigcort
        write(98,'(a,i3,2(a,f10.3))')
     &' Total for Z = ',nnz,' sig = ',sigt,' sigcor = ',sigcort
        sigctot = sigctot + sigcort
        sigutot = sigutot + sigt
       enddo
       write(*,'(2(a,f10.3,a/))')
     &' Grand total: sigcor = ',sigctot,' mbarn',
     &'  target elastic (but projectile excit.) = ',sigproj,' mbarn'
       write(98,'(2(a,f10.3,a/))')
     &' Grand total: sigcor = ',sigctot,' mbarn',
     &'  target elastic (but projectile excit.) = ',sigproj,' mbarn'
      endif
      write(*,'(/a,e12.4,a)')
     &' Total grazing cross section = ',cross,'  mbarn'
      write(98,'(/a,e12.4,a)')
     &' Total grazing cross section = ',cross,'  mbarn'
      sigtot = sigcaptn + cross
      write(*,'(a,e12.4,a)')
     &' Total capture cross section = ',sigcaptn,' mbarn'
      write(98,'(a,e12.4,a)')
     &' Total capture cross section = ',sigcaptn,' mbarn'
      write(*,'(a,e12.4,a)')
     &' Total reaction cross section = ', sigtot,' mbarn'
      write(98,'(a,e12.4,a)')
     &' Total reaction cross section = ', sigtot,' mbarn'
      close(unit=43)
      ierr = 0
2     return
      end
