c      Subroutine LEG 2 for program grazing
c      At exit ierr = 1 for error
c              ierr = 2 for warning
c              ierr = 3 for differential cross section
c              ierr = 4 for total cross section and new value of r
c              ierr = 0 for end or new value of energy
c      subroutines: convert_s,nzdistri_s,scatangleav_s
c---------------------------------------------------------------------
      subroutine leg2_s(nestep,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
c      list,ienergy,irangen and irangep are defined in grazing.icl
      common/za/za,aa,zaa,aaa
      common/ecm/ecm
      common/capt/eb,rbar,i5,intp
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      common/pop/rwarn,dtann,dtanp,dtaann,dtaanp
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
      common/iza/iza,iaa,izaa,iaaa
      common/q/qns,qnp,qps,qpp
      common/head/head
      common/char/elemp,elemt
      common/leg/ecminit,rainbcdg,decm,cs,r,delr,al,theta,
     &cscoul,cslab,thet,thetlabdeg,thetcou,ppcapt,ww,pr0,
     &drifta,driftz,wida,widz,eava,eavaa,tavs,widea,wideaa,avmin,avm,
     &changer,changers,ppscat,totmerg,dthet,
     &deltae,sig0drut,sigreac0,sigreac0lab,cross,alorho,rho,
     &istep,iout,ipopn,norbit,irstep,lab
      common/leg2/sigcaptn
      common/tot/in00,ip00,i0nmin,i0nmax,i0zmin,i0zmax
      common/scatang/delw,w0
      common/nzdis/ p0nn(irangen),p0zz(irangep),
     &eavdist(irangen,irangep),eavdist00,avmadist(irangen,irangep),
     &avmadist00,avmdist(irangen,irangep),avmdist00,
     &wid2dist(irangen,irangep),wid2dist00,wid2adist(irangen,irangep)
      common/nztdis/eav00,wid200
      dimension iselectn(list),iselectz(list),sig(irangen,irangep),
     &rav(irangen,irangep),ppnz(irangen,irangep),dsig(irangen,irangep),
     &dsiglab(irangen,irangep),widdist(irangen,irangep),
     &signsum(irangep),sigpsum(irangen),sigtota(irangea),
     &sigtotz(irangep),thetaver(irangen,irangep),
     &thetavlab(irangen,irangep),tabpr(20),alosig(irangen,irangep)
      character*5 elem1(list),elem2(list)
      character*5 elemp,elemt
      character*5 elempp,elemtt
      character*5 intp,strn,strp
      character*3 out
      character*7 outleg2(20)
      character*7 outlegg
      character*16 head
      character*50 xfilename1,efilename1
      character*34 xfilename
      character*34 efilenamea,efilenameb,efilenamec,efilenamed
      character*50 addtext
      character*5 version
      common/vers/version
c
      intprep = ierr
c      print*,' istep,irstep = ',istep,irstep
      if(istep.ne.1)goto 11
      if(nestep.ne.1)goto 8
21     write(*,'(3(/a)/a,f7.3,a,4(/a))')
     &' You have decided to calculate total cross sections',
     &'  for each N and Z,  integrated over excitation energy E*', 
     &'  and transferred angular momentum L.',
     &' Stepwidth in r is delr =',delr,' fm',
     &' Standard output is total reaction and merging cross sections',
     &'  and cross sections sig_NZ for each N and Z, for which',
     &'  the total cross section is larger than aminmb millibarns',
     &' Specify aminmb in mb'
       read(i5,*,err=101)aminmbb
       icheck = 0
       goto 111 
101    write(*,*)' you must provide a number'
       icheck = icheck + 1
       if(icheck.gt.5)goto 2
       goto 21 
111    addtext = '     !Minimum total cross section in mb'
       if(intp.eq.'inter')write(10,'(e20.12,a)')aminmbb,addtext
       aminmb = aminmbb/10.d0
c
5      if(version.eq.'vms  '.or.version.eq.'unix*')then
        if(intp.eq.'inter')write(*,'(/a/a/a)')
     &' You may obtain results for each r i.e. angular distributions',
     &'  stored in the grazing_out.dat file.',
     &' If so enter: yes. Enter: no if not.'
       elseif(version.eq.'unix ')then
        if(intp.eq.'inter')write(*,'(/a/a/a)')
     &' You may obtain results for each r i.e. angular distributions',
     &'  stored in the FOR098.DAT file.',
     &' If so enter: yes. Enter: no if not.'
       endif
       read(i5,'(a)')out
       addtext = '       !angular distributions?'
       if(intp.eq.'inter')write(10,'(a,a)')out,addtext
       if(out.eq.'no ')then
        leg2out = 0
        goto 8
       elseif(out.eq.'yes')then
        leg2out = 1
        outleg2(1) = 'dsig'
        if(intp.eq.'inter')write(*,'(4(/a))')
     &' You now obtain for each value of r:',
     &'  The mean scattering angle and differential cross section',
     &'  for each partition as well as the elastic scattering angle',
     &'  and the probability for capture and for elastic scattering.'
6       if(intp.eq.'inter')write(*,'(14(/a))')
     &' You may obtain further variables from the',
     &'  following list (in paranthesis their significance):',
     &' sigdrut (the elastic cross section divided by Rutherford),',
     &' ww (the absorptive potential W(r)),',
     &' tav (mean total kinetic energy loss for each partition),',
     &' eava (mean excitation energy of each ejectile),',
     &' eavaa (mean excitation energy of each recoil),',
     &' widtha (fluctuation in energy for each ejectile),',
     &' widthaa (fluctuation in energy for each recoil),',
     &' width (fluctuation in energy for each partition),',
     &' avma (mean angular momentum of each ejectile),',
     &' avmaa (mean angular momentum of each recoil),',
     &' avm (mean angular momentum loss for each partition).',
     &' Enter one of the names. Enter: no if no more'
        read(i5,'(a)')outlegg
        if(outlegg.ne.'sigdrut'.and.outlegg.ne.'ww'.and.outlegg.ne.'tav'
     &.and.outlegg.ne.'eava'.and.outlegg.ne.'eavaa'
     &.and.outlegg.ne.'widtha'.and.outlegg.ne.'widthaa'
     &.and.outlegg.ne.'width'.and.outlegg.ne.'avma'
     &.and.outlegg.ne.'avmaa'.and.outlegg.ne.'avm'
     &.and.outlegg.ne.'no')then
         write(*,*)' Wrong spelling, try again!'
         goto 6  
        endif
        addtext = '      !quantities for printout'
        if(intp.eq.'inter')write(10,'(a,a)')outlegg,addtext
        if(outlegg.eq.'no')then
         goto 9  
        else
         leg2out = leg2out+1
         outleg2(leg2out) = outlegg
        endif
        if(leg2out.eq.20)then
         write(*,*)' There is only room for 20 quantities!'
         goto 9  
        else
         goto 6  
        endif
       else
        write(*,*)' try again!'
        goto 5  
       endif      
9      if(intp.eq.'inter')write(*,'(a/a)')
     &' If you want a table of r, and the quantities you specified', 
     &'  enter tab.     If not enter: no'
       read(i5,'(a)')out
       addtext = '      !tab = tables of output,no if not'
       if(intp.eq.'inter')write(10,'(a,a)')out,addtext
       if(out.eq.'no ')then
        ntab = 0
       elseif(out.eq.'tab')then
        ntab = 1
       else
        write(*,*)' Try again!'
        goto 9  
       endif
       l = 1
58     if(intp.eq.'inter')write(*,'(a/a/a/a,i3,a)')
     &' Now you must specify the values of delta N and delta Z,',
     &'  i.e. the change in projectile neutron and proton numbers,',
     &'  for which the intermediate results are to be tabulated.',
     &' Enter',l,'. selected delta N and delta Z: iselectn,iselectz',
     &'  no no if no more'
       read(i5,*,err=104)iseln,iselz
       goto 114
104    if(intp.eq.'inter')write(10,'(a)')'no   no   ! more'
       lnzmax=l-1
       write(*,'(a,i2,a)')
     &' You have asked for ',lnzmax,' angular distributions'
       goto 8
114    addtext = '      !deltaN,deltaZ to be printed'
       if(intp.eq.'inter')write(10,'(2i4,a)')iseln,iselz,addtext 
       iselectn(l)=iseln 
       iselectz(l)=iselz
       l=l+1
       goto 58
c
c  preparation for total cross sections
8     continue
       if(version.eq.'vms  '.or.version.eq.'unix* ')then
        efilename1 = 'graz_'//head//'_nztot.dat'
        call condens_s(efilename1,efilenamea,le)
        open (unit=90,file=efilenamea,status='unknown')
        efilename1 = 'graz_'//head//'_adist.dat'
        call condens_s(efilename1,efilenameb,le)
        open (unit=91,file=efilenameb,status='unknown')
        efilename1 = 'graz_'//head//'_zdist.dat'
        call condens_s(efilename1,efilenamec,le)
        open (unit=92,file=efilenamec,status='unknown')
        efilename1 = 'graz_'//head//'_nzlog.dat'
        call condens_s(efilename1,efilenamed,le)
        open (unit=93,file=efilenamed,status='unknown')
       elseif(version.eq.'unix ')then
        efilenamea = 'FOR090.DAT'
        efilenameb = 'FOR091.DAT'
        efilenamec = 'FOR092.DAT'
        efilenamed = 'FOR093.DAT'
       endif
       do l=1,lnzmax
        ll = 19 + l
        iseln = iselectn(l)
        iselz = iselectz(l)
        if(version.eq.'vms  '.or.version.eq.'unix*')then
         call convert_s(iseln,strn,ierr)
         call convert_s(iselz,strp,ierr)
         xfilename1 = 'graz_'//head//'_'//strn//'_'//strp//'.dat'
         call condens_s(xfilename1,xfilename,le)
         open(unit=ll,file=xfilename,status='unknown')
        else
         write(xfilename,'(a4,i2,a4)')'FOR0',ll,'.DAT'
        endif
        ierr1 = 0
        call element_s(iza+iselz,iaa+iseln+iselz,elem1(l),ierr1)
        call element_s(izaa-iselz,iaaa-iseln-iselz,elem2(l),ierr1)
c        write(*,'(a,i3,a,i3,/,a,a5,a,a5,/,a,a40)')
       enddo
       cross = 0.d0
       do n=1,irangen
        do m=1,irangep
         sig(n,m)=0.d0
         rav(n,m)=0.d0
        enddo
       enddo       
c
c  calculate n,z distribution
11    continue
      ierr = 0
      call nzdistri_s(in00,ip00,inmax,inmin,izmax,izmin,ierr)
      if(ierr.eq.1)goto 2
      do i=1,irangen
       do j=1,irangep
        ppnz(i,j) = 0.d0 
       enddo
      enddo
      do i=inmin,inmax
       do j=izmin,izmax
        n=in00+i
        m=ip00+j
        ppnz(n,m)=p0nn(n)*p0zz(m)
        if(i.eq.0.and.j.eq.0)ppnz(n,m)=(p0nn(n)*p0zz(m)-p0)
       enddo
      enddo
      if(leg2out.eq.0)goto 46 
c
c  mean energies etc. including inelastic 
      pfac = p0nn(in00)*p0zz(ip00)
      fact = pfac/(pfac-p0)
      eav00 = (eavdist(in00,ip00)+eavina+eavinaa)*fact
c      print*,' eavina,eavinaa,eavdist = ',eavina,eavinaa
c     &,eavdist(in00,ip00)
      wid200 = (wid2dist(in00,ip00)+wid2ina+wid2inaa
     &+(eavdist(in00,ip00)+eavina+eavinaa)**2)*fact-eav00**2
c      print*,' wid2ina,wid2inaa = ',wid2ina,wid2inaa
      if(wid200.lt.0.d0)then
       write(*,*)' Warning: wid200 = ',wid200
       width00 = 0.d0
      else
       wid00 = sqrt(wid200)
       width00=2.35482d0*wid00
      endif
      facta = pfac/(pfac-p0tr*p0ina)
      eava00 = (eavdist(in00,ip00)/2+eavina)*facta
      wid2a00 = (wid2adist(in00,ip00)+wid2ina
     &+(eavdist(in00,ip00)/2.d0+eavina)**2)*facta-eava00**2
      if(wid2a00.lt.0.d0)then
       write(*,*)' Warning: wid2a00 = ',wid2a00
       widtha00 = 0.d-9
      else
       wida00 = sqrt(wid2a00)
       widtha00=2.35482d0*wida00
      endif
      factaa = pfac/(pfac-p0tr*p0inaa)
      eavaa00 = (eavdist(in00,ip00)/2+eavinaa)*factaa
      wid2aa00 = (wid2adist(in00,ip00)+wid2inaa
     &+(eavdist(in00,ip00)/2.d0+eavinaa)**2)*factaa-eavaa00**2
      if(wid2aa00.lt.0.d0)then
       write(*,*)' Warning: wid2aa00 = ',wid2aa00
       widthaa00 = 0.d-9
      else
       widaa00 = sqrt(wid2aa00)
       widthaa00=2.35482d0*widaa00
      endif
      avma00 = (avmadist(in00,ip00)+avmina)*facta
      avmaa00 = (avmdist(in00,ip00)-avmadist(in00,ip00)+avminaa)*factaa
      avm00 = (avmdist(in00,ip00)+avmin)*fact
c  
c  calculate differential cross sections
      do n=1,irangen
       do m=1,irangep
        dsig(n,m) = 0.d0
        dsiglab(n,m) = 0.d0
       enddo
      enddo
      do i=inmin,inmax 
       do j=izmin,izmax 
        n = in00+i
        m = ip00+j
        dsiglab(n,m)=cslab*ppnz(n,m)
        dsig(n,m)=cs*ppnz(n,m)
       enddo
      enddo
      do i=inmin,inmax
       do j=izmin,izmax
        n = in00+i 
        m = ip00+j
         if(i.gt.i0nmax.or.i.lt.i0nmin)goto 55
         if(j.gt.i0zmax.or.j.lt.i0zmin)goto 54
         do l=1,lnzmax
          if(iselectn(l).eq.i.and.iselectz(l).eq.j)goto 57
         enddo
         goto 54
57      deln = dble(i)
        delz = dble(j)
        dela = deln + delz
        qn1 = -qns
        if(deln.gt.0.d0)qn1=qnp
        qp1 = -qps
        if(delz.gt.0.d0)qp1=qpp
        delecm = -eavdist(n,m)-eavina-eavinaa+deln*qn1+delz*qp1
        dell = -avmdist(n,m)-avmina-avminaa
        w00 = w0
        dellw = delw       
        if(i.eq.0.and.j.eq.0)then
         delecm = -eavdist00-eava00-eavaa00
         dell = -avm00
        endif
        ierr = 0
c        print*,' al,ph,dph = ',al,ph,dph
c        print*,' rtex,delecm,dell ',rtex,delecm,dell
c        print*,' dela,delz,w00 ',dela,delz,w00
        call scatangleavx_s(al,ph,dph,rtex,delecm,dell,dela,delz
     &,w00,dellw,ierr)
        thetav = pi-ph
        dthetav = dph
        if(ierr.eq.1)then
         write(98,'(a,i2,a,i2,a)')
     &' Mean scattering angle for i = ',i,'j = ',j,' could not be found'
        else
         thetaver(n,m) = thetav
         sinth = sin(theta)
         if(sinth.lt.accur)then
          dsig(n,m) = dsig(n,m)*dthet/dthetav
         else
          dsig(n,m) = dsig(n,m)*sinth*dthet/sin(thetav)/dthetav
         endif
c  In laboratory
         taua = sqrt(ecm/(ecm+delecm)*aa/aaa*(aa+dela)/(aaa-dela))
         thetrlab = atan(sin(thetav)/(cos(thetav)+taua))
         if(thetrlab.lt.0.d0)thetrlab=thetrlab + pi
         thetavlab(n,m) = thetrlab
         if(abs(sin(thetav)).lt.accur)then
          dsiglab(n,m) = dsig(n,m)*(1.d0-taua)**2
         else
          dsiglab(n,m) = dsig(n,m)*(sin(thetav)/sin(thetrlab))**2
     &/abs(cos(thetav-thetrlab))
         endif
        endif
54      continue
       enddo
55     continue
      enddo
c
c  printing of results for each r
22    write(98,'(/a,f7.2,a,f6.2,a)')
     &'    Results for l= ',al,' hbar (r = ',r,' fm)' 
      write(*,'(/a,f7.2,a,f6.2,a)')
     &'    Results for l= ',al,' hbar (r = ',r,' fm)' 
      if(dtann.gt..5.or.dtaann.gt..5)
     &  write(98,'(a/a,4f6.3)')
     &' Warning: Maximum change in neutron single particle ',
     &'      populations for projectile and target = ',dtann,dtaann
      if(dtanp.gt..5.or.dtaanp.gt..5)
     &  write(98,'(a/a,4f6.3)')
     &' Warning: Maximum change in proton single particle ',
     &'      populations for projectile and target = ',dtanp,dtaanp
      tabpr(1) = r
      tabpr(2) = al
      if(lab.ne.1)then
       write(98,'(a,f7.2,a)')
     &' Elastic scattering angle = ',thet,' degrees'
       write(*,'(a,f7.2,a)')
     &' Elastic scattering angle = ',thet,' degrees'
      endif
      if(lab.ne.0)then
       write(98,'(a,f7.2,a)')
     &' Laboratory elastic scattering angle = ',thetlabdeg,' degrees'   
       write(*,'(a,f7.2,a)')
     &' Laboratory elastic scattering angle = ', thetlabdeg,' degrees'
      endif
      pp00 = p0*(1.d0-ppcapt)
      write(98,'(a,d9.2)')' Elastic probability = ',pp00
      write(*,'(a,d9.2)')' Elastic probability = ',pp00
      if(ppcapt.gt.accur)then 
       write(98,'(a,d9.2)')' Capture probability = ',ppcapt
       write(*,'(a,d9.2)')' Capture probability = ',ppcapt 
      endif
c
      do  k = 1,leg2out
       if(outleg2(k).eq.'sigdrut')then
        write(98,'(a,d9.2)')
     &' Elastic cross section/Rutherford =',sig0drut
        write(*,'(a,d9.2)')
     &' Elastic cross section/Rutherford =',sig0drut
       endif
       if(outleg2(k).eq.'ww')then
        write(98,'(a,f6.2,a,d9.2,a,d10.3)')' W(',r,') = ',ww,'  MeV'
        write(*,'(a,f6.2,a,d9.2,a,d10.3)')' W(',r,') = ',ww,'  MeV'
       endif
      enddo
c
      wid2in = wid2ina+wid2inaa
      do 45 ,i=inmin,inmax
      do 49 ,j=izmin,izmax
       n = in00+i 
       m = ip00+j
       www = wid2dist(n,m)+wid2in
       if(www.lt.0.d0)www=0.d0
       widdist(n,m)=sqrt(www)
        lpr = 0
        do l=1,lnzmax
         if(iselectn(l).eq.i.and.iselectz(l).eq.j)lpr = l
        enddo
        if(lpr.eq.0)goto 49
       if(lab.ne.1)then
        thedg = thetaver(n,m)*180.d0/pi
        tabpr(3) = thedg
        tabpr(4) = abs(dsig(n,m))
        write(*,'(a,f6.2,a,2i3,a,f8.3,a)')
     &' thetav = ',thedg,',  dsig(',i,j,')= ',dsig(n,m),' mb/sterad'
        write(98,'(a,f6.2,a,2i3,a,f8.3,a)')
     &' thetav = ',thedg,',  dsig(',i,j,')= ',dsig(n,m),' mb/sterad'
        if(i.eq.0.and.j.eq.0)then
         dsigel=cs*p0
         write(*,'(a,f6.2,a,f8.3,a)')
     &' theta = ',thet,' Elastic cross section =',dsigel,' mb/sterad'
         write(98,'(a,f6.2,a,f8.3,a)')
     &' theta = ',thet,' Elastic cross section =',dsigel,' mb/sterad'
        endif
       endif
       if(lab.ne.0)then
        theldg = thetavlab(n,m)*180.d0/pi
        tabpr(3) = theldg
        tabpr(4) = abs(dsiglab(n,m))
        write(*,'(a,f6.2,a,2i3,a,f8.3,a)')
     &' thetlab = ',theldg,' dsiglab(',i,j,')= ',dsiglab(n,m),' mb/st'
        write(98,'(a,f6.2,a,2i3,a,f8.3,a)')
     &' thetlab = ',theldg,' dsiglab(',i,j,')= ',dsiglab(n,m),' mb/st'
       endif
       kmax = leg2out
       itabpr = 4
       do k = 2,kmax      
        itabpr = itabpr+1
        if(outleg2(k).eq.'eava')then
         if(i.eq.0.and.j.eq.0)then
          tabpr(itabpr) = eava00
          write(*,'(a,f6.2,a)')
     &' Mean excitation of projectile excluding elastic =',eava00,' MeV'
          write(98,'(a,f6.2,a)')
     &' Mean excitation of projectile excluding elastic =',eava00,' MeV'
         else
          eavoa = eavdist(n,m)/2.d0+eavina
          write(98,'(a,f6.2,a)')
     &' Mean excitation energy of ejectile = ',eavoa,' MeV'
          write(*,'(a,f6.2,a)')
     &' Mean excitation energy of ejectile = ',eavoa,' MeV'
          tabpr(itabpr) = eavoa
         endif
        endif
        if(outleg2(k).eq.'eavaa')then
         if(i.eq.0.and.j.eq.0)then
          write(*,'(a,f6.2,a)')
     &' Mean excitation of target excluding elastic = ',eavaa00,' MeV'
          write(98,'(a,f6.2,a)')
     &' Mean excitation of target excluding elastic = ',eavaa00,' MeV'
          tabpr(itabpr) = eavaa00
         else
          eavoaa = eavdist(n,m)/2.d0+eavinaa
          write(98,'(a,f6.2,a)')
     &' Mean excitation energy of recoil =',eavoaa,' MeV'
          write(*,'(a,f6.2,a)')
     &' Mean excitation energy of recoil =',eavoaa,' MeV'
          tabpr(itabpr) = eavoaa
         endif
        endif
        if(outleg2(k).eq.'tav')then
         if(i.eq.0.and.j.eq.0)then
          write(*,'(a,f6.2,a)')
     &' Mean total kinetic energy loss excluding elastic =',eav00,' MeV'
          write(98,'(a,f6.2,a)')
     &' Mean total kinetic energy loss excluding elastic =',eav00,' MeV'
          tabpr(itabpr) = eava00
         else
          qn1 = -qns
          if(i.gt.0)qn1=qnp
          qp1 = -qps
          if(j.gt.0)qp1=qpp
          tavo = eavdist(n,m)+eavina+eavinaa-i*qn1-j*qp1
          write(98,'(a,f6.2,a)')
     &' Mean total kinetic energy loss = ',tavo,' MeV'
          write(*,'(a,f6.2,a)')
     &' Mean total kinetic energy loss = ',tavo,' MeV'
          tabpr(itabpr) = tavo
         endif
        endif
        if(outleg2(k).eq.'width')then
         if(i.eq.0.and.j.eq.0)then
          write(*,'(a,f6.2,a)')
     &' Energy fluctuation (FWHM) excluding elastic = ',width00,' MeV'
          write(98,'(a,f6.2,a)')
     &' Energy fluctuation (FWHM) excluding elastic = ',width00,' MeV'
          tabpr(itabpr) = width00
         else
          width=2.35482d0*widdist(n,m)
          write(*,'(a,f6.2,a)')
     &' Fluctuation in excitation energy(FWHM) = ',width,' MeV'
          write(98,'(a,f6.2,a)')
     &' Fluctuation in excitation energy (FWHM) = ',width,' MeV'
          tabpr(itabpr) = width
         endif
        endif
        if(outleg2(k).eq.'widtha')then
         if(i.eq.0.and.j.eq.0)then
          write(*,'(a/a,f6.2,a)')
     &' Excitation energy fluctuation (FWHM) of projectile',
     &'            excluding elastic = ',widtha00,' MeV' 
          write(98,'(a/a,f6.2,a)')
     &' Excitation energy fluctuation (FWHM) of projectile',
     &'            excluding elastic = ',widtha00,' MeV'
          tabpr(itabpr) = widtha00
         else
          width=2.35482d0*sqrt(wid2adist(n,m)+wid2ina)
          write(*,'(a,f6.2,a)')
     &' Fluctuation in excitation energy(FWHM) of ejectile = '
     &,width,' MeV'
          write(98,'(a,f6.2,a)')
     &' Fluctuation in excitation energy(FWHM) of ejectile = '
     &,width,' MeV'
          tabpr(itabpr) = width
         endif
        endif
        if(outleg2(k).eq.'widthaa')then
         if(i.eq.0.and.j.eq.0)then
          write(*,'(a/a,f6.2,a)')
     &' Excitation energy fluctuation (FWHM) of target',
     &'              excluding elastic = ',widthaa00,' MeV'
          write(98,'(a/a,f6.2,a)')
     &' Excitation energy fluctuation (FWHM) of target',
     &'             excluding elastic = ',widthaa00,' MeV'
          tabpr(itabpr) = widthaa00
         else
          width=2.35482d0*sqrt(wid2adist(n,m)+wid2inaa)
          write(*,'(a,f6.2,a)')
     &' Fluctuation in excitation energy(FWHM) of recoil = '
     &,width,' MeV'
          write(98,'(a,f6.2,a)')
     &' Fluctuation in excitation energy (FWHM) of recoil = '
     &,width,' MeV'
          tabpr(itabpr) = width
         endif
        endif
        if(outleg2(k).eq.'avma')then
         if(i.eq.0.and.j.eq.0)then
          write(*,'(a,f6.2,a)')
     &' Mean angular momentum of projectile excluding elastic = '
     &,avma00,' hbar'
          write(98,'(a,f6.2,a)')
     &' Mean angular momentum of projectile excluding elastic = '
     &,avma00,' hbar'
          tabpr(itabpr) = avma00
         else
          avmoa = avmadist(n,m)+avmina
          write(*,'(a,f6.2,a)')
     &' Mean angular momentum of ejectile = ',avmoa,' hbar'
          write(98,'(a,f6.2,a)')
     &' Mean angular momentum of ejectile = ',avmoa,' hbar'
          tabpr(itabpr) = avmoa
         endif
        endif
        if(outleg2(k).eq.'avmaa')then
         if(i.eq.0.and.j.eq.0)then
          write(*,'(a,f6.2,a)')
     &' Mean angular momentum of target excluding elastic = '
     &,avmaa00,' hbar'
          write(98,'(a,f6.2,a)')
     &' Mean angular momentum of target excluding elastic = '
     &,avmaa00,' hbar'
          tabpr(itabpr) = avmaa00
         else
          avmoaa = avmdist(n,m)-avmadist(n,m)+avminaa
          write(*,'(a,f6.2,a)')
     &' Mean angular momentum of recoil = ',avmoaa,' hbar'
          write(98,'(a,f6.2,a)')
     &' Mean angular momentum of recoil = ',avmoaa,' hbar'
          tabpr(itabpr) = avmoaa
         endif
        endif
        if(outleg2(k).eq.'avm')then
         if(i.eq.0.and.j.eq.0)then
          write(*,'(a,f6.2,a)')
     &' Mean total angular momentum  excluding elastic = ',avm00,' hbar'
          write(98,'(a,f6.2,a)')
     &' Mean total angular momentum  excluding elastic = ',avm00,' hbar'
          tabpr(itabpr) = avm00
         else
          avmo = avmdist(n,m)+avmin
          write(*,'(a,f6.2,a)')
     &' Mean total angular momentum = ',avmo,' hbar'
          write(98,'(a,f6.2,a)')
     &' Mean total angular momentum = ',avmo,' hbar'
          tabpr(itabpr) = avmo
         endif
        endif
        if(outleg2(k).eq.'ww')then
         write(*,'(a,f6.2,a)')
     &' Absorptive potential = ',ww,' MeV'
         write(98,'(a,f6.2,a)')
     &' Absorptive potential = ',ww,' MeV'
         tabpr(itabpr) = ww
        endif
       enddo
       itabmax = 3+leg2out
       if(itabpr.ne.itabmax)write(*,'(a,i3,a,i3)')
     &' itabpr = ',itabpr,' different from leg2out+3 = ',itabmax
c
c  produce tables as functions of r
       if(ntab.eq.0)goto 49 
       nn=19+lpr
       if(istep.eq.1)write(nn,'(a,10(2x,a7))')
     &'   r     l    thetav ',(outleg2(l),l=1,leg2out)
       write(nn,'(f6.2,f6.1,f7.1,f10.4,9(f9.2))')
     &(tabpr(itabpr),itabpr = 1,itabmax)
49    continue
45    continue
c
c  integration of cross sections (not accurate for irstep.eq.0)
46    if(istep.eq.1.and.totmerg.lt.0.d0)then
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
      do i=inmin,inmax
       do j=izmin,izmax
        n=in00+i
        m=ip00+j
        sig(n,m)= sig(n,m)+factor*ppnz(n,m)*ppscat
        rav(n,m)= rav(n,m)+r*factor*ppnz(n,m)*ppscat
       enddo
      enddo
      if(irstep.ne.0)irstep = irstep+1
      istep = istep+1
      if(intprep.eq.1.and.istep.ge.3)goto 33 
c  Check for finished integration
      if(abs(1.-p0).lt.accur)goto 33 
      if(abs(1.-p0tr).lt.accur)then
       write(*,*)' Only pure Coulomb excitation'
      endif
      alpr = al
      ierr = 4
      goto 2
c
c  printout of total cross sections
33    write(*,*)' ------------------------------------------'
      write(98,*)' ------------------------------------------'
      write(*,'(//a,f9.4,a)')
     &' Total grazing cross section = ',cross,'  mbarn'
      write(98,'(a,f8.3,a)')
     &' Total grazing cross section = ',cross,'  mbarn'
      sigtr = cross - sig(in00,ip00)
      write(*,'(a,f8.3,a)')
     &' Total transfer cross section = ',sigtr,'  mbarn'
      write(98,'(a,f8.3,a)')
     &' Total transfer cross section = ',sigtr,'  mbarn'
      sigtot = sigcaptn + cross
      write(*,'(a,e12.4,a)')
     &' Total capture cross section = ',sigcaptn,' mbarn'
      write(98,'(a,e12.4,a)')
     &' Total capture cross section = ',sigcaptn,' mbarn'
      write(*,'(a,e12.4,a)')
     &' Total reaction cross section = ', sigtot,' mbarn'
      write(98,'(a,e12.4,a)')
     &' Total reaction cross section = ', sigtot,' mbarn'
c
      write(12,'(2f10.4,1p4e12.4)') ecm,ecm-eb,cross,sigcaptn,sigtr
     1                              ,sigtot
c
      if(leg2out.ne.0)then
       do l=1,lnzmax
        i = iselectn(l)
        if(i.gt.i0nmax.or.i.lt.i0nmin)goto 81 
        j = iselectz(l)
        if(j.gt.i0zmax.or.j.lt.i0zmin)goto 81 
        n = i+in00
        m = j+ip00
        write(98,'(a,i3,a,i3,a,e10.3,a)')
     &' sig(',i,',',j,') = ',sig(n,m),' mb'
        if(sig(n,m).ne.0.d0)then
         rav(n,m) = rav(n,m)/sig(n,m)
         write(98,'(a,f6.3,a)')
     &' Mean distance of closest approach for this partition = '
     &,rav(n,m),' fm'
         ravv = rav(n,m)
         if(ravv.lt.rwarn)write(98,'(a)')
     &' Warning: Depopulation larger than .5 at this distance'
        endif
81      continue
       enddo
      endif
c printout of total cross sections for given N and Z
c   first determine range to be tabulated
      do m=1,irangep
       signsum(m)= 0.d0
      enddo
      do n=1,irangen
       sigpsum(n)= 0.d0 
      enddo
      do n=i0nmin+in00,i0nmax+in00
       do m = i0zmin+ip00,i0zmax+ip00
        sigpsum(n)= sigpsum(n)+abs(sig(n,m))
       enddo
      enddo
      do m=i0zmin+ip00,i0zmax+ip00
       do n = i0nmin+in00,i0nmax+in00
        signsum(m)= signsum(m)+abs(sig(n,m))
       enddo
      enddo
      do n=i0nmin+in00,in00
       if(sigpsum(n).gt.aminmbb)then
        goto 90
       endif
      enddo 
90    i1nmin = n - in00 -2
      do n=in00+1,i0nmax+in00
       if(sigpsum(n).lt.aminmbb)then
        goto 91 
       endif
      enddo
91    i1nmax = n-in00+1
      do m=i0zmin+ip00,ip00
       if(signsum(m).gt.aminmbb)then
        goto 92 
       endif
      enddo
92    i1pmin = m - ip00 - 1
      do m=ip00+1,i0zmax+ip00
       if(signsum(m).lt.aminmbb)then
        goto 94 
       endif
      enddo
94    i1pmax = m-ip00+1
      write(98,'(/a)')' Total cross sections in mbarn for given N and Z'
      write(*,'(/a)')' Total cross sections in mbarn for given N and Z'
      ina0 = iaa - iza
      ntable = int((i1nmax-i1nmin)/10)+1
      do nn=1,ntable
       ina = min(i1nmax-i1nmin-(nn-1)*10,9)
       write(98,'(a,10i7)')' Z/N = ',(ina0+i1nmin+(nn-1)*10+i, i=0,ina)
       write(*,'(a,10i7)')' Z/N = ',(ina0+i1nmin+(nn-1)*10+i, i=0,ina)
       do j=i1pmax,i1pmin,-1
        izb = iza+j
        m=ip00+j
        n0 = in00+i1nmin +(nn-1)*10
        write(*,'(i7,10f7.2)')izb,(sig(n,m),n=n0,n0+ina)
        write(98,'(i7,10f7.2)')izb,(sig(n,m),n=n0,n0+ina)
       enddo
      enddo
      chsigtot = 0.d0
      do n=i0nmin+in00,i0nmax+in00
       do m = i0zmin+ip00,i0zmax+ip00
        chsigtot = chsigtot + sig(n,m)
       enddo
      enddo 
      write(*,'(a,f10.3,/a,f10.3)')' total grazing cross section =',
     & cross,' Sum of cross sections =',chsigtot
c
c   make table for contour plotting
      toler = aminmbb/100.d0
      do  j = ip00,ip00+i0zmax
       if(signsum(j).lt.toler)goto 77 
      enddo
77    jmax=j-1 
      do j = ip00-1,ip00+i0zmin,-1
       if(signsum(j).lt.toler)goto 95 
      enddo
95    jmin=j+1 
      ymin = dble(jmin-ip00)
      ymax = dble(jmax-ip00)
      dy = 1.
      do i = in00,in00+i0nmax
       if(sigpsum(i).lt.toler)goto 87 
      enddo
87    imax=i-1 
      do i = in00-1,in00+i0nmin,-1 
       if(sigpsum(i).lt.toler)goto 93
      enddo
93    imin=i+1 
      xmin = dble(imin-in00)
      xmax = dble(imax-in00)
      dx = 1.
      write(90,'(1x,6e12.4)')xmin,xmax,dx,ymin,ymax,dy
      write(93,'(1x,6e12.4)')xmin,xmax,dx,ymin,ymax,dy
      do jy=jmin,jmax
        write(90,'(1x,1p8e10.2)') (sig(ix,jy),ix=imin,imax)
        do ix=imin,imax
         if(sig(ix,jy).le.0.d0)then
          alosig(ix,jy)= -25.d0
         else
          alosig(ix,jy)=log10(sig(ix,jy))
         endif
        enddo
        write(93,'(1x,8f10.5)') (alosig(ix,jy),ix=imin,imax)
      enddo
c
c  Table of A distribution
      write(98,*)'     A distribution of ejectile in mbarn'
      do i=1,irangea
       sigtota(i) = 0.d0
      enddo
      ia00 = in00+ip00
      do idelaa = i0nmin+i0zmin,i0nmax+i0zmax
       do i=i0nmin,i0nmax
        j = idelaa - i
        if(j.ge.i0zmin.and.j.le.i0zmax)then
         sigtota(ia00+idelaa)= sigtota(ia00+idelaa)+sig(in00+i,ip00+j)
        endif
       enddo
       if(sigtota(ia00+idelaa).gt.toler)then
        write(91,*)idelaa+iaa,sigtota(ia00+idelaa)
        write(98,'(i10,f12.5)')idelaa+iaa,sigtota(ia00+idelaa)
       endif
      enddo
c
c  Table of Z distribution
      write(98,*)'     Z distribution of ejectile in mbarn'
      do i=1,irangep
       sigtotz(i) = 0.d0
      enddo
      do idelzz = i0zmin,i0zmax
       do i=i0nmin,i0nmax
        sigtotz(ip00+idelzz)=sigtotz(ip00+idelzz)+
     &sig(in00+i,ip00+idelzz)
       enddo
       if(sigtotz(ip00+idelzz).gt.toler)then
        write(92,*)idelzz+iza,sigtotz(ip00+idelzz)
        write(98,'(i10,f12.5)')idelzz+iza,sigtotz(ip00+idelzz)
       endif
      enddo
c
      close(unit=90)
      close(unit=91)
      close(unit=92)
      close(unit=93)
      write(*,*)' -------------------------------------------------'
      write(*,'(a,4(/a,a))')
     &' Total cross section distributions for plotting:',
     &'             sig_NZ    are stored in ',efilenamea,
     &'                   A-distribution in ',efilenameb,
     &'              and  Z-distribution in ',efilenamec,
     &'        log_10(sig_NZ) are stored in ',efilenamed
      write(98,*)' -------------------------------------------------'
      write(98,'(a,4(/a,a))')
     &' Total cross section distributions for plotting:',
     &'             sig_NZ    are stored in ',efilenamea,
     &'                   A-distribution in ',efilenameb,
     &'              and  Z-distribution in ',efilenamec,
     &'        log_10(sig_NZ) are stored in ',efilenamed
      if(ntab.eq.1)then
       do l=1,lnzmax
        ll = 19 + l
        iseln = iselectn(l)
        iselz = iselectz(l)
        if(version.eq.'vms  '.or.version.eq.'unix*')then
         call convert_s(iseln,strn,ierr)
         call convert_s(iselz,strp,ierr)
         xfilename1 = 'graz_'//head//'_'//strn//'_'//strp//'.dat'
         call condens_s(xfilename1,xfilename,le)
         close(unit=ll)
        else
         write(xfilename,'(a4,i2,a4)')'FOR0',ll,'.DAT'
         close(unit=ll)
        endif
        ierr1 = 0
        call element_s(iza+iselz,iaa+iseln+iselz,elempp,ierr1)
        elem1(l) = elempp 
        call element_s(izaa-iselz,iaaa-iseln-iselz,elemtt,ierr1)
        elem2(l) = elemtt
        write(*,'(a,i3,a,i3,/,a,a5,a,a5,/,a,a40)')
     &' Angular distributions for delta N = ',iselectn(l),
     &' and delta Z = ',iselectz(l),
     &'  i.e. for the partition ',elem1(l),' + ',elem2(l),
     &'  is stored in ',xfilename
       enddo
      endif


      ierr = 0
2     return
      end
