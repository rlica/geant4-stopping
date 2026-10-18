c      Subroutine LEG 42 for program grazing
c      At exit ierr = 1 for error
c              ierr = 2 for warning
c              ierr = 3 for differential cross section
c              ierr = 4 for total cross section and new value of r
c              ierr = 0 for end or new value of energy
c---------------------------------------------------------------------
      subroutine leg42_s(nestep,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
c      include 'common.fil'
c      list,ienergy,irangen and irangep are defined in grazing.icl
      common/za/za,aa,zaa,aaa
      common/ecm/ecm
      common/elab/elab
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
      dimension pnzevap(irangen,irangep),sigevap(irangen,irangep)
     &,sigevapnsum(irangep),sigevappsum(irangen)
      character*5 intp
      character*16 head
      character*5 elemp
      character*5 elemt
      character*34 efilename
cnan
      character*13 headshort,headshort1
      character*40 file
      character*40 filex
cnan
      character*5 version
      common/vers/version
c
      print*,' p0tr,p0ina,p0inaa = ',p0tr,p0ina,p0inaa 
      iprep = ierr
      if(istep.ne.1)goto 51
      if(nestep.ne.1)goto 15
       write(*,'(/a/a/a,f7.3,a/a/a)')
     &' You have decided to calculate total cross sections for each N',
     &' and Z integrated over L, after evaporation',
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
1110  if(ievap.ne.1.and.ievap.ne.2)goto 48
      if(intp.eq.'inter')write(10,*)ievap
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
         sigevap(n,m) = 0.d0
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
      call enzcor_s(amass,aamass,amassapp,aamassapp,in00,ip00
     &,inmin,inmax,izmin,izmax,pnzcor,dele,evapdis,iout,ierr)
      if(ierr.eq.1)goto 2
      ierr = ievap
      call evap_s(amass,amassapp,aamass,aamassapp,in00,ip00
     &,inmin,inmax,izmin,izmax,pnzevap,dele,evapdis,iout,ierr)
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
      do i=inmin-2,inmax+2
       do j=izmin-2,izmax+2
        n=in00+i
        m=ip00+j
        sigevap(n,m)= sigevap(n,m)+factor*pnzevap(n,m)*ppscat
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
c   printout of total cross sections 
33    write(*,*)'------------------------------------------------'
      write(98,*)'------------------------------------------------'
      write(*,'(a,f10.3,a)')
     &' Total reaction cross section = ',cross,' mbarn'
      write(98,'(a,f10.3,a)')
     &' Total reaction cross section = ',cross,' mbarn'
      sigctot = 0.d0
      sigevtot = 0.d0
      if(ievap.eq.1)then
       write(*,'(a)')
     &' Total cross sections for ejectiles after evaporation in mbarn'
       write(98,'(a)')
     &' Total cross sections for ejectiles after evaporation in mbarn'
cnan
      headshort1 = elemp//elemt//'_Z'
      call condens_s(headshort1,headshort,len)
      file(1:len)=headshort(1:len)
cnan
       do j=i0zmin,i0zmax
        m=ip00+j
        sigcort = 0.d0
        sigevapt = 0.d0
        nnz = iza + j
        write(*,'(a,i3)')' For Z = ',nnz
        write(98,'(a,i3)')' For Z = ',nnz
cnan
        if(nnz.le.9)then
           len1=len+1
           len2=len+1
           write(file(len1:len2),'(i1)') nnz
          elseif(nnz.le.99)then
           len1=len+1
           len2=len+2
           write(file(len1:len2),'(i2)')nnz
          elseif(nnz.le.999)then
           len1=len+1
           len2=len+3
           write(file(len1:len2),'(i3)')nnz
        endif
        file(len2+1:len2+1)='_'
        len2=len2+1
        iielab=elab+0.0001
        if(iielab.le.9)then
           len1=len2+1
           len2=len2+1
           write(file(len1:len2),'(i1)') iielab
          elseif(iielab.le.99)then
           len1=len2+1
           len2=len2+2
           write(file(len1:len2),'(i2)')iielab
          elseif(iielab.le.999)then
           len1=len2+1
           len2=len2+3
           write(file(len1:len2),'(i3)')iielab
          elseif(iielab.le.9999)then
           len1=len2+1
           len2=len2+4
           write(file(len1:len2),'(i4)')iielab
        endif
        filex=file(1:len2)//'.dat'
        open(unit=33,file=filex,status='new')
cnan
        do i=i0nmin,i0nmax
         n=in00+i
         naa = iaa + i + j
         write(*,'(a,i3,2(a,1pe12.4))')' A = ',naa,
     &' Before: sigcor = ',sigcor(n,m),'  sigevap = ',sigevap(n,m)
         write(98,'(a,i3,2(a,1pe12.4))')' A = ',naa,
     &' Before: sigcor = ',sigcor(n,m),'  sigevap = ',sigevap(n,m)
cnan
         write(33,'(i6,1p2e12.4)') naa,sigcor(n,m),sigevap(n,m)
cnan
         if(i.eq.0.and.j.eq.0.and.inelast.eq.1)then
          write(*,'(a,f7.3,a)')
     &' projectile elastic (target excitation) sig_el = ',sigtar,' mb'
          write(98,'(a,f7.3,a)')
     &' projectile elastic (target excitation) sig_el = ',sigtar,' mb'
         endif
         sigcort = sigcort+sigcor(n,m)
         sigevapt = sigevapt+sigevap(n,m)
        enddo
cnan
        close(33)
cnan
        write(*,'(a,i3,2(a,f10.3))')' Total for Z = ',nnz,' sigcor = '
     &,sigcort,' sigevap = ',sigevapt
        write(98,'(a,i3,2(a,f10.3))')' Total for Z = ',nnz,' sigcor = '
     &,sigcort,' sigevap = ',sigevapt
        sigctot = sigctot + sigcort
        sigevtot = sigevtot + sigevapt
       enddo
       write(*,'(4(a,f10.3,a/))')
     &' Grand total: sigcor = ',sigctot,' mb',
     &' projectile elastic (target excitation) = ',sigtar,' mb',
     &' sigevap = ',sigevtot,' mb',' sig_merg = ',totmerg,' mb'
       write(98,'(4(a,f10.3,a/))')
     &' Grand total: sigcor = ',sigctot,' mb',
     &' projectile elastic (target excitation) = ',sigtar,' mb',
     &' sigevap = ',sigevtot,' mb',' sig_merg = ',totmerg,' mb'
      elseif(ievap.eq.2)then
       write(*,'(a)')
     &' Total cross sections for recoils after evaporation in mbarn'
       write(98,'(a)')
     &' Total cross sections for recoils after evaporation in mbarn'
cnan
      headshort1 = elemp//elemt//'_Z'
      call condens_s(headshort1,headshort,len)
      file(1:len)=headshort(1:len)
      print '(a,i6)', file,len
cnan
       do j=i0zmin,i0zmax
        m=ip00+j
        sigcort = 0.d0
        sigevapt = 0.d0
        nnz = izaa - j
        write(*,'(a,i3)')' For Z = ',nnz
        write(98,'(a,i3)')' For Z = ',nnz
cnan
        if(nnz.le.9)then
           len1=len+1
           len2=len+1
           write(file(len1:len2),'(i1)')nnz
          elseif(nnz.le.99)then
           len1=len+1
           len2=len+2
           write(file(len1:len2),'(i2)')nnz
          elseif(nnz.le.999)then
           len1=len+1
           len2=len+3
           write(file(len1:len2),'(i3)')nnz
        endif
        file(len2+1:len2+1)='_'
        len2=len2+1
        iielab=elab+0.0001
        if(iielab.le.9)then
           len1=len2+1
           len2=len2+1
           write(file(len1:len2),'(i1)') iielab
          elseif(iielab.le.99)then
           len1=len2+1
           len2=len2+2
           write(file(len1:len2),'(i2)')iielab
          elseif(iielab.le.999)then
           len1=len2+1
           len2=len2+3
           write(file(len1:len2),'(i3)')iielab
          elseif(iielab.le.9999)then
           len1=len2+1
           len2=len2+4
           write(file(len1:len2),'(i4)')iielab
        endif
        filex=file(1:len2)//'.dat'
        open(unit=33,file=filex,status='new')
cnan
        do i=i0nmin,i0nmax
         n=in00+i
         naa = iaaa - i - j
         write(*,'(a,i3,2(a,1pe12.4))') ' A = ',naa,
     &' Before: sigcor = ',sigcor(n,m),'  sigevap = ',sigevap(n,m)
         write(98,'(a,i3,2(a,1pe12.4))')' A = ',naa,
     &' Before: sigcor = ',sigcor(n,m),'  sigevap = ',sigevap(n,m)
cnan
         write(33,'(i6,1p2e12.4)') naa,sigcor(n,m),sigevap(n,m)
cnan  
        if(i.eq.0.and.j.eq.0.and.inelast.eq.1)then
          write(*,'(a,f7.3,a)')
     &' target elastic (projectile excitation) sig_el = ',sigproj,' mb'
          write(98,'(a,f7.3,a)')
     &' target elastic (projectile excitation) sig_el = ',sigproj,' mb'
         endif
         sigcort = sigcort+sigcor(n,m)
         sigevapt = sigevapt+sigevap(n,m)
        enddo
cnan
        close(33)
cnan
        write(*,'(a,i3,2(a,f10.3))')
     &' Total for Z = ',nnz,' sigcor = ',sigcort,' sigevap = ',sigevapt
        write(98,'(a,i3,2(a,f10.3))')
     &' Total for Z = ',nnz,' sigcor = ',sigcort,' sigevap = ',sigevapt
        sigctot = sigctot + sigcort
        sigevtot = sigevtot + sigevapt
       enddo
       write(*,'(4(a,f10.3,a/))')
     &' Grand total: sigcor = ',sigctot,' mb',
     &'  target elastic (projectile excitation) = ',sigproj,' mb',
     &'  sigevap = ',sigevtot,' mb',' sig_merg = ',totmerg,' mb'
       write(98,'(4(a,f10.3,a/))')
     &' Grand total: sigcor = ',sigctot,' mb',
     &'  target elastic (projectile excitation) = ',sigproj,' mb',
     &'  sigevap = ',sigevtot,' mb',' sig_merg = ',totmerg,' mb'
      endif
c
c printout of total cross sections after evaporationfor given N and Z
c   first determine range to be tabulated
      do m=1,irangep
       sigevapnsum(m)= 0.d0
      enddo
      do n=1,irangen
       sigevappsum(n)= 0.d0
      enddo
      do n=i0nmin+in00,i0nmax+in00
       do m = i0zmin+ip00,i0zmax+ip00
        sigevappsum(n)= sigevappsum(n)+abs(sigevap(n,m))
       enddo
      enddo
      do m=i0zmin+ip00,i0zmax+ip00
       do n = i0nmin+in00,i0nmax+in00
        sigevapnsum(m)= sigevapnsum(m)+abs(sigevap(n,m))
       enddo
      enddo
      do n=i0nmin+in00,in00
       if(sigevappsum(n).gt.accur)then
        goto 90
       endif
      enddo
90    i1nmin = n - in00 -2
      do n=in00+1,i0nmax+in00
       if(sigevappsum(n).lt.accur)then
        goto 491
       endif
      enddo
491   i1nmax = n-in00+1
      do m=i0zmin+ip00,ip00
       if(sigevapnsum(m).gt.accur)then
        goto 492
       endif
      enddo
492   i1pmin = m - ip00 - 1
      do m=ip00+1,i0zmax+ip00
       if(sigevapnsum(m).lt.accur)then
        goto 494
       endif
      enddo
494   i1pmax = m-ip00+1
      if(ievap.eq.1)then
       write(98,'(/a/a)')
     &' Total cross sections in mbarn for given N and Z',
     &'    for ejectiles after evaporation'
       ina0 = iaa - iza
       ntable = int((i1nmax-i1nmin)/10)+1
       do nn=1,ntable
        ina = min(i1nmax-i1nmin-(nn-1)*10,9)
        write(98,'(a,10i7)')
     &' Z\N = ',(ina0+i1nmin+(nn-1)*10+i, i=0,ina)
        write(*,'(a,10i7)')
     &' Z\N = ',(ina0+i1nmin+(nn-1)*10+i, i=0,ina)
        do j=i1pmax,i1pmin,-1
         izb = iza+j
         m=ip00+j
         n0 = in00+i1nmin +(nn-1)*10
         write(*,'(i7,10f7.2)')izb,(sigevap(n,m),n=n0,n0+ina)
         write(98,'(i7,10f7.2)')izb,(sigevap(n,m),n=n0,n0+ina)
        enddo
       enddo
       if(inelast.eq.1)then
        write(*,'(a,f7.3,a)')
     &' projectile elastic (target excitation) sig_el = ',sigtar,' mb'
        write(98,'(a,f7.3,a)')
     &' projectile elastic (target excitation) sig_el = ',sigtar,' mb'
       endif
      elseif(ievap.eq.2)then
       write(*,'(/a/a)')
     &' Total cross sections in mbarn for given N and Z',
     &'    for recoils after evaporation'
       ina0 = iaaa - izaa
       ntable = int((i1nmax-i1nmin)/10)+1
        do nn=1,ntable
        ina = min(i1nmax-i1nmin-(nn-1)*10,9)
        write(98,'(a,10i7)')
     &' Z\N = ',(ina0-i1nmin-(nn-1)*10-i, i=0,ina)
        write(*,'(a,10i7)')' Z\N = ',(ina0-i1nmin-(nn-1)*10-i, i=0,ina)
        do j=i1pmax,i1pmin,-1
         izbb = izaa-j
         m=ip00+j
         n0 = in00+i1nmin +(nn-1)*10
         write(*,'(i7,10f7.2)')izbb,(sigevap(n,m),n=n0,n0+ina)
         write(98,'(i7,10f7.2)')izbb,(sigevap(n,m),n=n0,n0+ina)
        enddo
       enddo
       if(inelast.eq.1)then
        write(*,'(a,f7.3,a)')
     &' target elastic (projectile excitation) sig_el = ',sigproj,' mb'
        write(98,'(a,f7.3,a)')
     &' target elastic (projectile excitation) sig_el = ',sigproj,' mb'
       endif
      endif
      chsigevaptot = 0.d0
      do n=i0nmin+in00,i0nmax+in00
       do m = i0zmin+ip00,i0zmax+ip00
         chsigevaptot = chsigevaptot + sigevap(n,m)
       enddo
      enddo
      write(*,'(a,f10.3,/a,f10.3)')' total grazing cross section =',
     & cross,' Sum of cross sections =',chsigevaptot
c
c   make table for contour plotting
      if(version.eq.'vms  '.or.version.eq.'unix*')then
       efilename = 'graz_'//head//'_nzevap.dat'
       open (unit=93,file=efilename,status='unknown')
      else
       efilename = 'FOR093.DAT'
      endif
      write(*,'(a/a,a31)')' Table of NZ distribution after evaporation'
     &,' is stored in ',efilename
      toler = accur
      do  j = ip00,ip00+i0zmax
       if(sigevapnsum(j).lt.toler)goto 477
      enddo
477   jmax=j
      do j = ip00-1,ip00+i0zmin,-1
       if(sigevapnsum(j).lt.toler)goto 495
      enddo
495   jmin=j
      ymin = dble(jmin-ip00)
      ymax = dble(jmax-ip00)
      dy = 1.
      do i = in00,in00+i0nmax
        if(sigevappsum(i).lt.toler)goto 487
      enddo
487   imax=i
      do i = in00-1,in00+i0nmin,-1
       if(sigevappsum(i).lt.toler)goto 493
      enddo
493   imin=i
      xmin = dble(imin-in00)
      xmax = dble(imax-in00)
      dx = 1.
      write(93,'(1x,6e12.4)')xmin,xmax,dx,ymin,ymax,dy
      do jy=jmin,jmax
        write(93,'(1x,1p8e10.2)') (sigevap(ix,jy),ix=imin,imax)
      enddo
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
      close(unit=93)
      close(unit=44)
      ierr = 0
2     return
      end
