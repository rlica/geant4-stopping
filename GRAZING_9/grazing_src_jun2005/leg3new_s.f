c      Subroutine LEG 3 for program grazing including inelastic
c      At exit ierr = 1 for error
c              ierr = 2 for warning
c              ierr = 3 for differential cross section
c              ierr = 4 for total cross section and new value of r
c              ierr = 0 for finished integration
c---------------------------------------------------------------------
      subroutine leg3_s(nestep,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
c      list,ienergy,irangen and irangep are defined in grazing.icl
c
      dimension thetav(ienergy,101),thetavlab(ienergy,101),rstep(101)
     &,dsigr(ienergy,101),eloss(ienergy),angtet(181),angdis(181)
      dimension table(10,11,ienergy),tableh(10,11,3)
      dimension sigetot(ienergy),sigeapptot(ienergy),ppe(ienergy),
     &ppapp(ienergy),driftndis(ienergy),
     &driftzdis(ienergy),avmtdis(ienergy),changerdis(ienergy),
     &dsig(ienergy),dsiglab(ienergy)
c      common/char/elemp,elemt
      logical*1 edistout,driftout,mavtout,scatout
      character*5 elemp
      character*5 elemt
      character*5 intp
      character*3 out
      character*7 outleg3(10)
      character*7 outlegg
      character*7 name
      character*36 filename,filename1
      character*36 efilename,efilename1
      character*60 addtext
      character*16 head
      character*5 version
c      common/vers/version
c
      common/za/za,aa,zaa,aaa
      common/ecm/ecm
      common/capt/eb,rbar,i5,intp
c      common:akyuz
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
     &wid2dist(irangen,irangep),wid2dist00,wid2adist(irangen,irangep)
      common/head/head
      common/leg/ecminit,rainbcdg,decm,cs,r,delr,al,theta,
     &cscoul,cslab,thet,thetlabdeg,thetcou,ppcapt,ww,pr0,
     &drifta,driftz,wida,widz,eava,eavaa,tavs,widea,wideaa,avmin,avm,
     &changer,changers,ppscat,totmerg,dthet,
     &deltae,sig0drut,sigreac0,sigreac0lab,cross,alorho,rho,
     &istep,iout,ipopn,norbit,irstep,lab
      common/scatang/delw,w0
      common/angdis/deltet
      common/leg2/sigcaptn
      common/leg3/delintloss,rcapt,alcapt
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/f0muin/ff0muin(2,3,2,4),f0inm(2,3,2)
      common/giantwid/giantwid(3)
      common/vers/version
      common/char/elemp,elemt
c
      iprep = ierr
      rstep(istep) = r
      if(istep.ne.1)goto 12
      if(nestep.ne.1)goto 11
      edistout = .false.
      driftout = .false.
      mavtout = .false.
      scatout = .false.
c
c  Total cross sections
10    continue
      write(*,'(/a/a/a/a/a,f7.3,a)')
     &' You have decided to calculate total cross sections',
     &'  as function of kinetic energy loss or excitation energy,',
     &'  i.e. the cross section integrated over angles',
     &'  summed over N and Z and integrated over transferred L',
     &' Stepwidth in r is delr =',delr,' fm'
303   if(intp.eq.'inter')write(*,'(5(a/))')
     &'  You have the following choices:',
     &' 1. Distribution in total kinetic energy', 
     &' 2. Distribution in excitation energy of projectile',
     &' 3. Distribution in excitation energy of target',
     &' 4. Wilczynski plot',
     &' Enter your choice:'
      read(i5,*,err=1102)leg3
      icheck = 0
      goto 1112
1102  write(*,*)' you must provide a number'
      icheck = icheck + 1
      if(icheck.gt.5)goto 2
      goto 303
1112  addtext = '   !options for leg3, cannot be changed to or from 1!'
      if(intp.eq.'inter')write(10,'(i5,a)')leg3,addtext
      if(leg3.lt.1.or.leg3.gt.4)then
       write(*,*)' your choice must be 1,2,3 or 4!'
       goto 303
      endif
313   if(intp.eq.'inter')write(*,'(a/a/a/a)')
     &' You must calculate the energy distributions',
     &'  with a finite energy resolution (FWHM).',
     &' Enter width in MeV',
     &'  (if you enter a number <1 MeV, width is set equal to 1 MeV)'
      read(i5,*,err=1103)tres
      icheck = 0
      if(tres.lt.1.)tres = 1.d0
      goto 1113
1103  write(*,*)' you must provide a number'
      icheck = icheck + 1
      if(icheck.gt.5)goto 2
      goto 313
1113  addtext = '    !energy resolution in MeV'
      if(intp.eq.'inter')write(10,'(f7.4,a)')tres,addtext
314   if(intp.eq.'inter')write(*,'(a/a/a/a)')
     &' Since the total cross section diverges you should enter',
     &'  the minimum elastic scattering angle to be included', 
     &'  (in degrees) for large impact parameters,',
     &'  if negative, this angle is measured from the rainbow angle'
      read(i5,*,err=1104)scatminn
      icheck = 0
      goto 1114
1104  write(*,*)' you must provide a number'
      icheck = icheck + 1
      if(icheck.gt.5)goto 2
      goto 314
1114  if(rainbcdg.gt.130.d0)then
       scatmin = 100.d0
      else
       if(scatminn.lt.accur)then
        scatmin = rainbcdg + scatminn
       elseif(scatminn.gt.rainbcdg)then
        scatmin = rainbcdg
       else
        scatmin = scatminn
       endif
      endif
      addtext = '   !minimum scattering angle i degrees'
      if(intp.eq.'inter')write(10,'(f8.3,a)')scatmin,addtext
c
c  Choice of output
      if(leg3.eq.2.or.leg3.eq.3)then
       if(version.eq.'vms  '.or.version.eq.'unix* ')then
        if(intp.eq.'inter')write(*,'(a/a)')
     &' You obtain the energy distribution for each trajectory',
     &'   in the grazing_out.dat file.' 
       else
        if(intp.eq.'inter')write(*,'(a/a)')
     &' You obtain the energy distribution for each trajectory',
     &'   in the FOR098.DAT file.' 
       endif
       ntab = 1
       goto 16
      elseif(leg3.eq.4)then
       edistout = .true.
       driftout = .true.
       mavtout = .true.
       scatout = .true.
       leg3out=0
       goto 30
      elseif(leg3.eq.1)then
       if(version.eq.'vms  '.or.version.eq.'unix* ')then
        if(intp.eq.'inter')write(*,'(a/a)')
     &' You may obtain the results for each trajectory',
     &' stored in grazing_out.dat file. Enter yes, or no if not.'
       else
        if(intp.eq.'inter')write(*,'(a/a)')
     &' You may obtain the results for each trajectory',
     &' stored in FOR098.DAT file. Enter yes, or no if not.'
       endif
       read(i5,'(a)')out
       addtext = '    !additional output?'
       if(intp.eq.'inter')write(10,'(a,a)')out,addtext
       if(out.eq.'no ')then
        leg3out = 0
       elseif(out.eq.'yes')then
        leg3out = 1
        edistout = .true.
        outleg3(1) = 'edist'
307     if(intp.eq.'inter')write(*,'(a/a/a/a)')
     &' You now obtain for each value of r:',
     &'  The elastic scattering angle,',
     &'  the elastic cross section divided by Rutherford,',
     &'  and the kinetic energy distribution.'
       endif      
      endif
      if(intp.eq.'inter')write(*,'(a/a)')
     &' Besides the energy distribution (edist)',
     &'  You may calculate variables from the following list'
30    if(leg3.eq.4.and.intp.eq.'inter')write(*,'(a/a)')
     &' Besides the TKE versus angle (Wilczynski plot)',
     &'  You may tabulate variables from the following list'
      if(intp.eq.'inter')write(*,'(6(a/),a)')
     &'  (in paranthesis their significance):',
     &' driftn (mean drift in N for given energy loss)',
     &' driftz (mean drift in Z for given energy loss)',
     &' mavdis (mean ang. momentum loss for given energy loss)',
     &' scatang (mean scattering angle for given energy loss)',
     &' dsigde (double differential cross section)',
     &' Enter name of variable, no if no more'
      read(i5,'(a)')outlegg
      if(outlegg.eq.'driftn')then
       driftout = .true.
      elseif(outlegg.eq.'driftz')then
       driftout = .true.
      elseif(outlegg.eq.'mavdis')then
       mavtout = .true.
      elseif(outlegg.eq.'scatang')then
       scatout = .true.
       if(lab.eq.1)then
        outlegg = 'scatlab'
       elseif(lab.eq.2)then
        leg3out = leg3out+1
        outleg3(leg3out) = 'scatlab'
       endif
      elseif(outlegg.eq.'dsigde')then
       scatout = .true.
       if(lab.eq.1)then
        outlegg = 'dsiglab'
       elseif(lab.eq.2)then
        leg3out = leg3out+1
        outleg3(leg3out) = 'dsiglab'
       endif
      elseif(outlegg.eq.'edist')then
       write(*,*)' edist is automatically chosen'
       goto 30
      elseif(outlegg.ne.'no')then
       write(*,*)' Wrong spelling, try again!'
      endif
      addtext = '      !additional quantities'
      if(intp.eq.'inter')write(10,'(a,a)')outlegg,addtext
      if(outlegg.eq.'no')then
       goto 309
      else
       leg3out= leg3out+1
       outleg3(leg3out) = outlegg
       goto 30
      endif
309   if(leg3out.eq.0)goto 11
      if(leg3.eq.4)then
       edistout = .true.
       leg3out = leg3out+1
       outleg3(leg3out) = 'edist'
       ntab = 1
       goto 16
      endif
      if(intp.eq.'inter')write(*,'(a/a)')
     &' If you want tables of each of the quantities you specified', 
     &'  enter yes. Otherwise enter no'
      read(i5,'(a)')out
      addtext = '    !tables of quantities specified?'
      if(intp.eq.'inter')write(10,'(a,a)')out,addtext
      if(out.eq.'yes')then
       ntab = 1
      elseif(out.eq.'no ')then
       ntab = 0
      else
       write(*,*)' Try again!'
       goto 309
      endif
c
c  Printout for each trajectory found in:
16    if(ntab.eq.1)then
       if(leg3.eq.1.or.leg3.eq.4)then
        do n=1,leg3out
         nn=20+n    
         name = outleg3(n)    
         if(version.eq.'vms  '.or.version.eq.'unix* ')then
          filename1 = 'graz_'//head//'_'//name//'.dat'
          call condens_s(filename1,filename,le)
          open(unit=nn,file=filename,status='unknown')
         else
          write(filename,'(a4,i2,a4)') 'FOR0',nn,'.DAT'
         endif
        enddo
       elseif(leg3.eq.2)then
        leg3out = 1
        if(version.eq.'vms  '.or.version.eq.'unix* ')then
         efilename1 = 'graz_'//head//'_edista.dat'
         call condens_s(efilename1,efilename,le)
         open (unit=21,file=efilename,status='unknown')
        else
         efilename = 'FOR021.DAT'
        endif
       elseif(leg3.eq.3)then
        leg3out = 1
        if(version.eq.'vms  '.or.version.eq.'unix*')then
         efilename1 = 'graz_'//head//'_edistaa.dat'
         call condens_s(efilename1,efilename,le)
         open (unit=21,file=efilename,status='unknown')
        else
         efilename = 'FOR021.DAT'
        endif
       endif
      endif       
c  and for total:
      if(leg3.eq.1.or.leg3.eq.4)then
       if(version.eq.'vms  '.or.version.eq.'unix*')then
        efilename1 = 'graz_'//head//'_tdisttot.dat'
        call condens_s(efilename1,efilename,le)
        open (unit=95,file=efilename,status='unknown')
       else
        efilename = 'FOR095.DAT'
       endif
      endif
      if(leg3.eq.2)then
       if(version.eq.'vms  '.or.version.eq.'unix* ')then
        efilename1 = 'graz_'//head//'_edistatot.dat'
        call condens_s(efilename1,efilename,le)
        open (unit=95,file=efilename,status='unknown')
       else
        efilename = 'FOR095.DAT'
       endif
c       write(*,'(a/,a40)')
c     &' Total energy distribution of ejectile stored in ',efilename
      endif
      if(leg3.eq.3)then
       if(version.eq.'vms  '.or.version.eq.'unix*')then
        efilename1 = 'graz_'//head//'_edistaatot.dat'
        call condens_s(efilename1,efilename,le)
        open (unit=95,file=efilename,status='unknown')
       else
        efilename = 'FOR095.DAT'
       endif
      endif
c
c  preparation for total cross sections
11    continue
      cross = 0.d0
c     crossres = 0.d0
      do  mm=1,ienergy
       sigetot(mm) = 0.d0
       sigeapptot(mm)=0.d0
      enddo
      mmaxmax = 0
      mminmin = ienergy
      imm = 0
c
c   calculate energy distributions
12    ierr1 = 1
      ierr2 = 1
      dsigdr=10.d0*r*(ecm-za*zaa*esq/2.d0/r
     &-unuc(r,ierr1)-r/2.d0*dunuc(r,ierr2))/ecm
      if(ierr1.eq.1.or.ierr2.eq.1)then
       write(*,*)' error in leg3 '
       ierr = 1
       goto 2
      endif
      write(*,'(/a,f7.4)')
     & '  Calculating distribution in energy for r = ',r
      timeedist1 = X05BAF()
      if(leg3.eq.1.or.leg3.eq.4)then
       deleloss = .5d0*min(deleloss,tres)
       deleloss = dble(1.d0+int(10.*deleloss))/10.d0
       if(ppscat.lt.accur/10.d0)then
        do mm=1,ienergy
         ppe(mm) = 0.d0
        enddo
        goto 29
       else
        deltae = tres/2.35482d0
        ecut = 5.d0*tres
        ierr = iprep
cnan        acc = accur/10.d0
        acc = accur
        call tdistri_s(deltae,deleloss,acc,imm,mmmax,mmmin,
     &eloss,ppe,ierr)
        if(ierr.eq.1)goto 2
        if(ierr.gt.1)write(*,'(a,f10.7)')'    for r = ',r
       endif
       do mm = 1,ienergy
        dsigr(mm,istep) = 0.d0
       enddo
       do mm=mmmin,mmmax
        dsigr(mm,istep) = dsigdr*ppe(mm)*ppscat
       enddo
       if(driftout.or.scatout)then
        write(*,*)
     &'   now finding average drift in N and Z for each energy'
        ierr = iprep
        call driftnztdis_s(deltae,deleloss,accur,imm,mmmax,mmmin,
     &eloss,driftndis,driftzdis,ierr)
        if(ierr.eq.1)goto 2
        if(ierr.gt.1)write(*,'(a,f10.7)')'     for r = ',r
        do mm=mmmin,mmmax
         if(ppe(mm).ne.0.d0)then
          dn = driftndis(mm)/ppe(mm)
          driftndis(mm) = dn
          dz = driftzdis(mm)/ppe(mm)
          driftzdis(mm) = dz
          if(aa+dn.lt.1.d0.or.aaa-dn.lt.1.d0.or.aa+dn-za-dz.lt.1.d0.or.
     &za+dz.lt.1.d0.or.zaa-dz.lt.1.d0.or.aaa-dn-zaa+dz.lt.1.d0)then
           write(*,*)' Warning: too large drifts for eloss = ',eloss(mm)
           if(eloss(mm).lt.tavs)then
            mmmin = mm+1
            goto 24
           elseif(eloss(mm).gt.tavs)then
            mmmax = mm-1
            goto 23
           endif
          endif
         endif
24       continue         
        enddo
23      continue
       endif 
       if(mavtout.or.scatout)then
        write(*,*)
     &'   now finding average angular momentum for each energy'
        ierr = iprep
cnan        acc = accur/10.d0
        acc = accur
        call mavtdist_s(deltae,deleloss,acc,imm,mmmax,mmmin,
     &eloss,avmtdis,abserr,ierr)
       if(ierr.eq.1)goto 2
        if(ierr.gt.1)write(*,'(a,f10.7)')'     for r = ',r
        do mm=mmmin,mmmax
         if(ppe(mm).ne.0.d0)then
          avd = avmtdis(mm)/ppe(mm)
          avmtdis(mm) = avd
          if(al-avd.lt.0.1d0.or.abserr/al.gt.ppe(mm))then
           write(*,*)
     &' Warning: too large error in avm for eloss = ',eloss(mm)
           if(eloss(mm).lt.tavs)then
            mmmin = mm+1
            goto 25
           elseif(eloss(mm).gt.tavs)then
            mmmax = mm-1
            goto 26
           endif
          endif
         endif
25       continue           
        enddo
26      continue
       endif 
       if(scatout)then
        if(iprep.eq.1)then
         mmmin = imm -3
         mmmax = imm
        endif
        do ii = 1,ienergy
         thetav(ii,istep) = 0.d0
         thetavlab(ii,istep) = 0.d0
        enddo
        do mm = mmmin,mmmax
         ierr = 0
         dela = driftndis(mm)+driftzdis(mm)
         delz = driftzdis(mm)
         delecm = - eloss(mm)
         delal = - avmtdis(mm)
         w00 = w0
         dellw = delw
         call scatangleavx_s(al,ph,dph,rturnex,delecm,delal,dela,delz
     &,w00,dellw,ierr)
         if(ierr.eq.1)goto 2
         changerdis(mm) = rturnex-r
         thetaver = pi-ph
cnan
cnan         write(99,*) al,thetaver*180/pi
cnan
         dthetav = dph
         if(ph.eq.0.d0)then
          thetav(mm,istep) = 0.d0
          dsig(mm) = 0.d0
          dsiglab(mm) = 0.d0
          goto 27
         endif
         thetav(mm,istep) = thetaver
         sinth = sin(theta)
         if(sinth.gt.accur)then
          dsig(mm) = cs*ppe(mm)*sinth*dthet/sin(thetaver)/dthetav
         else
          dsig(mm) = cs*ppe(mm)*dthet/dthetav
         endif
c  In laboratory
         taua = sqrt(ecm/(ecm+delecm)*aa/aaa*(aa+dela)/(aaa-dela))
         thetrlab = atan(sin(thetaver)/(cos(thetaver)+taua))
         if(thetrlab.lt.0.d0)thetrlab=thetrlab + pi
         thetavlab(mm,istep) = thetrlab
         if(abs(sin(thetaver)).lt.accur)then
          dsiglab(mm) = dsig(mm)*(1.d0-taua)**2
         else
          dsiglab(mm) = dsig(mm)*(sin(thetaver)
     &/sin(thetrlab))**2/abs(cos(thetaver-thetrlab))
         endif
c         endif
27       continue         
        enddo
       endif
      endif
29    continue 
      if(leg3.ge.2.and.leg3.le.3)then
       if(leg3.eq.2)iq = 1
       if(leg3.eq.3)iq = 2
       deleloss = dble(1.d0+int(10.d0*deleloss))/10.d0
       if(ppcapt.gt.1.d0-5.d0*accur)then
        do mm=1,ienergy
         ppe(mm)= 0.d0
         ppapp(mm) = 0.d0
        enddo
       else
        deltae = tres/2.35482d0
        ierr = iprep
        call edista_s(iq,deltae,deleloss,accur,imm,mmmax,mmmin,
     &eloss,ppe,ierr)
        if(ierr.gt.1)write(*,*)'     for r = ',r
        ierr = iprep
        eavap = eavtr/2.d0/(1.d0-p0tr)
        wid2ap = wid2tr/3.d0/(1.d0-p0tr)-p0tr*eavap**2
        p0ap = p0tr
        call edistapp_s(iq,deltae,eavap,wid2ap,p0ap,deleloss,accur,
     &imm,mmmax,mmmin,eloss,ppapp,ierr)
        if(ierr.gt.1)write(*,*)'     for r = ',r
        do mm = mmmin,mmmax
         write(*,'(a,f8.3,2(a,f7.5))')
     &' eloss = ',eloss(mm),' ppe = ',ppe(mm),' ppapp = ',ppapp(mm)
        enddo
       endif
      endif
      if(mmmax.gt.mmaxmax)mmaxmax = mmmax
      if(mmmin.lt.mminmin)mminmin = mmmin
c  mmmax = maximum of mm calculated
      timeedist2 = X05BAF()
      timeedist = timeedist2 - timeedist1
      write(*,'(a,f7.2,a)')' CPU time for energy distribution = '
     &,timeedist,' sec'
c
c  Printout for each trajectory
      if(iprep.eq.1)then
cnan       goto 2
         goto 22
      endif
      if(leg3out.eq.0)then
       goto 15 
      else
       write(*,'(/a,f7.2,a,f5.2,a)')
     &'   Results for l = ',al,' (r = ', r,' fm)'
       write(98,'(/a,f7.2,a,f5.2,a)')
     &'   Results for l = ',al,' (r = ', r,' fm)'
       if(lab.ne.1)then
        write(98,'(a,f7.2,a)')
     &' Elastic scattering angle = ',thet,' degrees'
        write(*,'(a,f7.2,a)')
     &' Elastic scattering angle = ',thet,' degrees'
       elseif(lab.ne.0)then
        write(98,'(a,f7.2,a)')
     &' Elastic laboratory scattering angle = ',thetlabdeg,' degrees'   
        write(*,'(a,f7.2,a)')
     &' Elastic laboratory scattering angle = ',thetlabdeg,' degrees'
       endif
       write(98,'(a,d11.4)')
     &' Elastic cross section/Rutherford = ',sig0drut
       write(*,'(a,d11.4)')
     &' Elastic cross section/Rutherford = ',sig0drut
       write(98,'(a,d11.4)')
     &' Capture probability = ',ppcapt
       write(*,'(a,d11.4)')
     &' Capture probability = ',ppcapt
       if(istep.eq.1.and.leg3.eq.2)then
        write(98,'(a/a/a,f7.3,a/a)')
     &' Distribution in projectile excitation energy',
     &'  integrated over mass and charge drift',
     &'  with energy  resolution (FWHM) =',tres,' MeV',
     &'  E_loss              Prob             Prob_app'
       elseif(istep.eq.1.and.leg3.eq.3)then
        write(98,'(a/a/a,f7.3,a/a)')
     &' Distribution in target excitation energy',
     &'  integrated over mass and charge drift',
     &'  with energy  resolution (FWHM) =',tres,' MeV',
     &'  E_loss             Prob               Prob_app'
       elseif(istep.eq.1.and.leg3.eq.1.or.istep.eq.1.and.leg3.eq.4)then
        if(lab.eq.0)then
         write(98,'(a/a/a,f7.3,a)')
     &' Distribution in total kinetic energy loss in center of mass',
     &'  integrated over mass and charge drift',
     &'  with energy  resolution (FWHM) =',tres,' MeV'
         write(98,'(a,a)')
     &'  T_loss  Probability  <del_N>   <del_Z>     <L>    theta  ',   
     &'   dsig/dO/dT'
        elseif(lab.eq.1)then
         write(98,'(a/a/a,f7.3,a)')
     &' Distribution in total kinetic energy loss in lab',
     &'  integrated over mass and charge drift',
     &'  with energy  resolution (FWHM) =',tres,' MeV'
         write(98,'(a,a)')
     &'  T_loss  Probability  <del_N>   <del_Z>     <L>    thet_lab',   
     &' dsig_lab/dO/dT'
        endif
       endif
       do 317,mm = mmmin,mmmax
        if(leg3.eq.2.or.leg3.eq.3)then
         write(98,'(f16.2,2(f16.5))') eloss(mm),ppe(mm),ppapp(mm)
        elseif(leg3.eq.1.or.leg3.eq.4)then
         thedg = thetav(mm,istep)*180.d0/pi
         thelbdg = thetavlab(mm,istep)*180.d0/pi
         tke=ecm-eloss(mm)
         if(lab.eq.0)then
          write(98,'(a,f5.1,a,f6.4,4(a,f6.2),a,e10.3)')
     & ' P(',eloss(mm),') = ',ppe(mm),'/Mev  <del_N> = ',driftndis(mm),
     &' <del_Z> = ',driftzdis(mm),' <M>= ',avmtdis(mm),
     &' thet= ',thedg,' dsig= ',dsig(mm)
         elseif(lab.ne.0)then
          write(98,'(a,f5.1,a,f6.4,4(a,f6.2),a,e10.3)')
     & ' P(',eloss(mm),') = ',ppe(mm),'/Mev  <del_N> = ',driftndis(mm),
     &' <del_Z> = ',driftzdis(mm),' <M>= ',avmtdis(mm),
     &' thet= ',thelbdg,' dsig= ',dsiglab(mm)
         endif
        endif
317    continue   
      endif
c
c  Printout of tables
15    if(ntab.eq.1)then
       if(mod(istep,10).eq.1) then
        do m=1,11
         do leo = 1,leg3out
          tableh(leo,m,1)=0.d0
          tableh(leo,m,2)=0.d0
          tableh(leo,m,3)=0.d0
          do n=1,ienergy
           table(leo,m,n)=0.d0
          enddo
         enddo
        enddo
       endif
       do leo = 1,leg3out
        nn = 20+leo
        if(istep.eq.1) then
         if(leg3.eq.2)then
          write(nn,'(a/a/a)')
     &' Double differential cross sections for projectile excitation',
     &'  as functions of excitation energy and angle',
     &'  (negative inside rainbow)'
          write(nn,'(a,f7.3,a/a,f8.2,a)')
     &'     with ',tres,' MeV energy resolution (FWHM)'
     &,'     and scattering angles larger than ',scatmin,' degrees' 
         elseif(leg3.eq.3)then
          write(nn,'(a/a/a)')
     &' Double differential cross sections for target excitation',
     &'  as functions of excitation energy and angle',
     &'  (negative inside rainbow)'
          write(nn,'(a,f7.3,a/a,f8.2,a)')
     &'     with ',tres,' MeV energy resolution (FWHM)'
     &,'     and scattering angles larger than ',scatmin,' degrees' 
         elseif(leg3.eq.1.or.leg3.eq.4)then
          if(outleg3(leo).eq.'edist')then
      print*,' nn = ',nn
           write(*,'(a)')
     &' Distributions in total kinetic energy loss for each trajectory'
           write(nn,'(a)')
     &' Distributions in total kinetic energy loss for each trajectory'
          elseif(outleg3(leo).eq.'driftn')then
           write(nn,'(a)')
     &' Distributions in drift of neutrons for each trajectory'
          elseif(outleg3(leo).eq.'driftz')then
           write(nn,'(a)')
     &' Distributions in drift of protons for each trajectory'
          elseif(outleg3(leo).eq.'mavdis')then
           write(nn,'(a)')
     &' Distributions in transferred ang. momentum for each trajectory'
          elseif(outleg3(leo).eq.'scatang')then
           write(nn,'(a)')
     &' Distributions in cm scattering angle for each trajectory'
          elseif(outleg3(leo).eq.'scatlab')then
           write(nn,'(a)')
     &' Distributions in lab scattering angle for each trajectory'
          elseif(outleg3(leo).eq.'dsigde')then
           write(nn,'(a/a)')
     &' Double differential cross section for each trajectory',
     &'  in mb/MeV/sr in center of mass' 
          elseif(outleg3(leo).eq.'dsiglab')then
           write(nn,'(a/a)')
     &' Double differential cross section for each trajectory',
     &'  in mb/MeV/sr in laboratory system' 
          endif
          write(nn,'(a,f7.3,a/a,f7.2,a)')
     &'     with ',tres,' MeV energy resolution (FWHM)',
     &'     and elastic scattering angles > ',scatmin,' degrees' 
          write(nn,'(a)')
     &'----------------------------------------------------------------'
         endif
        endif
        jcol = mod(istep-1,10)+2 
        tableh(leo,jcol,1)=r
        tableh(leo,jcol,2)=thet
        if(outleg3(leo).eq.'scatlab'.or.outleg3(leo).eq.'dsiglab')
     &tableh(leo,jcol,2)=thetlabdeg
        tableh(leo,jcol,3)=sigreac0
        do mm=1,mmaxmax
         table(leo,1,mm)= eloss(mm)
        enddo
        do mm=mmmin,mmmax 
         if(leg3.eq.2.or.leg3.eq.3)then
          table(leo,jcol,mm)= ppe(mm)*cs
         endif
         if(leg3.eq.1.or.leg3.eq.4)then
          if(outleg3(leo).eq.'edist')then
           table(leo,jcol,mm)= ppe(mm)
          elseif(outleg3(leo).eq.'driftn')then
           table(leo,jcol,mm)= driftndis(mm)
          elseif(outleg3(leo).eq.'driftz')then
           table(leo,jcol,mm)= driftzdis(mm)
          elseif(outleg3(leo).eq.'mavdis')then
           table(leo,jcol,mm)= avmtdis(mm)
          elseif(outleg3(leo).eq.'scatang')then
           thedg = thetav(mm,istep)*180.d0/pi
           table(leo,jcol,mm)= thedg
          elseif(outleg3(leo).eq.'scatlab')then
           thelbdg = thetavlab(mm,istep)*180.d0/pi
           table(leo,jcol,mm)= thelbdg
          elseif(outleg3(leo).eq.'dsigde')then
           table(leo,jcol,mm)= dsig(mm)
          elseif(outleg3(leo).eq.'dsiglab')then
           table(leo,jcol,mm)= dsiglab(mm)
          endif
         endif
        enddo
        if(jcol.eq.11)then
         write(nn,'(a,10f11.2)')'      r   ',(tableh(leo,n,1),n=2,11)
         write(nn,'(a,10f11.2)')'    theta ',(tableh(leo,n,2),n=2,11)
         write(nn,'(a,10e11.3)')'  dsig_r/dl',(tableh(leo,n,3),n=2,11)
         write(nn,'(/)')
         if(outleg3(leo).eq.'dsigde'.or.outleg3(leo).eq.'dsiglab')then
          do mm=1,mmaxmax
           write(nn,'(f10.3,10e11.3)')
     &table(leo,1,mm),(table(leo,n,mm),n=2,11)
          enddo
         else
          do mm=1,mmaxmax
          write(nn,'(f10.3,10f11.5)')
     &table(leo,1,mm),(table(leo,n,mm),n=2,11)
          enddo
         endif
         write(nn,*)
     &'----------------------------------------------------------------'
         do m=1,11
          do n=1,ienergy
           table(leo,m,n)=0.d0
          enddo
          if(leo.eq.leg3out)then
           tableh(leo,m,1)=0.d0
           tableh(leo,m,2)=0.d0
           tableh(leo,m,3)=0.d0
          endif
         enddo
        endif
       enddo
      endif 
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
      do n=mmmin,mmmax
       sigetot(n)=sigetot(n)+factor*ppe(n)*ppscat
       if(leg3.eq.2.or.leg3.eq.3)then
        sigeapptot(n)=sigeapptot(n)+factor*ppapp(n)*ppscat
       endif
      enddo
cnan
cnan  before it was 50 but with delal=2 was to smal
      if(istep+1.gt.100)then
cnan
       write(*,*)' istep+1.gt.50 in leg3'
       goto 14
      endif
      if(irstep.ne.0)irstep = irstep+1
      istep = istep+1
      if(iprep.eq.1)then
       ierr = 0
       goto 2
      endif
c  Check for finished integration
      if(tres.gt.accur.and.thet.lt.scatmin.and.dthet.gt.0.)then
       goto 14
      elseif(leg3.eq.4.and.thet.gt.scatmin)then
       ierr = 4
       goto 2
      elseif(abs(1.-p0).lt.accur)then
       goto 14
      else
       alpr = al
       ierr = 4
       goto 2
      endif
14    continue
      istmax=istep
c
      if(ntab.eq.1.and.jcol.ne.11)then
c  printout of remaining differential
       do leo = 1,leg3out
        nn = 20+leo
        write(nn,'(a,10f11.2)')'      r   ',(tableh(leo,n,1),n=2,jcol)
        write(nn,'(a,10f11.2)')'    theta ',(tableh(leo,n,2),n=2,jcol)
        write(nn,'(a,10e11.3)')'  dsig_r/dl',(tableh(leo,n,3),n=2,jcol)
        write(nn,'(/)')
        if(outleg3(leo).eq.'dsigde'.or.outleg3(leo).eq.'dsiglab')then
         do mm=1,mmaxmax
          write(nn,'(f10.3,10e11.3)')
     &table(leo,1,mm),(table(leo,n,mm),n=2,jcol)
         enddo
        else
         do mm=1,mmaxmax
         write(nn,'(f10.3,10f11.5)')
     &table(leo,1,mm),(table(leo,n,mm),n=2,jcol)
         enddo
        endif
       enddo
      endif
c
c  Printout of total cross section distribution
      write(*,*)' --------------------------------------------'
      write(98,*)' --------------------------------------------'
      if(leg3.eq.2)then
       write(95,'(a/a,f7.3,a/a,f8.2,a/a)')
     &' Excitation energy total cross section of projectile', 
     &'     with ',tres,' MeV energy resolution (FWHM)',
     &'     and scattering angles larger than ',scatmin,' degrees',
     &' eloss MeV  sig mb/MeV'
       write(98,'(a/a,f7.3,a/a,f8.2,a/a)')
     &' Excitation energy total cross section of projectile', 
     &'     with ',tres,' MeV energy resolution (FWHM)',
     &'     and scattering angles larger than ',scatmin,' degrees',
     &' eloss MeV  sig mb/MeV'
       write(*,'(a/a,f7.3,a/a,f8.2,a/a)')
     &' Excitation energy total cross section of projectile', 
     &'     with ',tres,' MeV energy resolution (FWHM)',
     &'     and scattering angles larger than ',scatmin,' degrees',
     &' eloss MeV  sig mb/MeV  sig_app mb/MeV'
       do n=mminmin,mmaxmax
        write(95,'(f8.1,2d15.4)') eloss(n),sigetot(n)
        write(98,'(f8.1,2d15.4)') eloss(n),sigetot(n)
        write(*,'(f8.1,3d15.4)') eloss(n),sigetot(n),sigeapptot(n)
       enddo
      elseif(leg3.eq.3)then
       write(95,'(a/a,f7.3,a/a,f8.2,a/a)')
     &' Excitation energy total cross section of target',
     &'     with ',tres,' MeV energy resolution (FWHM)',
     &'     and scattering angles larger than ',scatmin,' degrees',
     &' eloss MeV   sig mb/MeV'
       write(98,'(a/a,f7.3,a/a,f8.2,a/a)')
     &' Excitation energy total cross section of target',
     &'     with ',tres,' MeV energy resolution (FWHM)',
     &'     and scattering angles larger than ',scatmin,' degrees',
     &' eloss MeV   sig mb/MeV'
       write(*,'(a/a,f7.3,a/a,f8.2,a/a)')
     &' Excitation energy total cross section of target',
     &'     with ',tres,' MeV energy resolution (FWHM)',
     &'     and scattering angles larger than ',scatmin,' degrees',
     &' eloss MeV   sig mb/MeV   sig_app mb/MeV'
       do n=mminmin,mmaxmax
        write(95,'(f8.1,2d15.4)') eloss(n),sigetot(n)
        write(98,'(f8.1,2d15.4)') eloss(n),sigetot(n)
        write(*,'(f8.1,3d15.4)') eloss(n),sigetot(n),sigeapptot(n)
       enddo
      elseif(leg3.eq.1.or.leg3.eq.4)then
       write(95,'(a/a,f7.3,a/a,f8.2,a/a)')
     &' Total kinetic energy loss cross section',
     &'     with ',tres,' MeV energy resolution (FWHM)',
     &'     and scattering angles larger than ',scatmin,' degrees',
     &' teloss MeV      sig mb/MeV'
       write(98,'(a/a,f7.3,a/a,f8.2,a/a)')
     &' Total kinetic energy loss cross section',
     &'     with ',tres,' MeV energy resolution (FWHM)',
     &'     and scattering angles larger than ',scatmin,' degrees',
     &' teloss MeV      sig mb/MeV'
       write(*,'(a/a,f7.3,a/a,f8.2,a/a)')
     &' Total kinetic energy loss cross section',
     &'     with ',tres,' MeV energy resolution (FWHM)',
     &'     and scattering angles larger than ',scatmin,' degrees',
     &' teloss MeV      sig mb/MeV'
       do n=mminmin,mmaxmax
        if(sigetot(n).gt.0.d0)then
         write(95,'(f7.1,1pd15.4)') eloss(n),sigetot(n)
         write(98,'(f7.1,1pd15.4)') eloss(n),sigetot(n)
         write(*,'(f7.1,1pd15.4)') eloss(n),sigetot(n)
        endif
       enddo
      endif
22    continue
      if(leg3.eq.4)then
       widang = deltet
       if(lab.ne.1)then
        ierr = 0
        nnfile = 57
        call wilsh_s(thetav,thetavlab,dsigr,rstep,eloss,widang,ecut
     &,angtet,angdis,la,nnfile,mminmin,mmaxmax,istmax,ierr)
        if(ierr.eq.1)goto 2
        if(version.eq.'unix '.or.version.eq.'unix*')then 
         write(95,'(a/a/a,f7.3,a/a,f8.2,a/a)')
     &' The differential cross section in center of mass as function',
     &'  of total kinetic energy loss and angle (Wilczynski plot)',
     &'     with ',tres,' MeV energy resolution (FWHM)',
     &'     and scattering angles larger than ',scatmin,' degrees',
     &'  is stored in FOR057.dat ready for plotting'
        elseif(version.eq.'vms  ')then
         ierr = la
         call wilsh_plot_s(ierr)
        endif
        write(*,'(a,f6.1/a)')
     &' Angular distribution integrated over energy loss >',ecut,
     &'     theta_cm      dsig_cm/dO (mb/sr)'
        write(98,'(a,f6.1/a)')
     &' Angular distribution integrated over energy loss >',ecut,
     &'     theta_cm      dsig_cm/dO (mb/sr)'
        angtot = 0.d0
        do nn=1,181
         ang = angtet(nn)*pi/180.d0
         if(nn.eq.1)then
          angp = 0.d0
         else
          angp = angtet(nn-1)*pi/180.d0
         endif
         if(angdis(nn).gt.0.d0)then
          write(*,'(f14.2,1pe23.3)')angtet(nn),angdis(nn)
          write(98,'(f14.2,1pe23.3)')angtet(nn),angdis(nn)
         endif
         angtot = angtot + angdis(nn)*2.d0*pi*sin(ang)*(ang-angp) 
        enddo
        write(*,'(a,1pd12.4,a)')' Total cross section = ',angtot,' mb'
        write(98,'(a,1pd12.4,a)')' Total cross section = ',angtot,' mb'
       endif        
       if(lab.ne.0)then
        ierr = 1
        nnfile = 58
        call wilsh_s(thetav,thetavlab,dsigr,rstep,eloss,widang,ecut
     &,angtet,angdis,la,nnfile,mminmin,mmaxmax,istmax,ierr)
        if(ierr.eq.1)goto 2
        if(version.eq.'unix '.or.version.eq.'unix*')then 
         write(95,'(a/a/a,f7.3,a/a,f8.2,a/a)')
     &' The differential cross section in the laboratory as function',
     &'  of total kinetic energy loss and angle (Wilczynski plot)',
     &'     with ',tres,' MeV energy resolution (FWHM)',
     &'     and scattering angles larger than ',scatmin,' degrees',
     &'  is stored in FOR058.dat ready for plotting'
        elseif(version.eq.'vms  ')then
         ierr = la
         call wilsh_plot_s(ierr)
        endif
        write(*,'(a,f6.1/a)')
     &' Angular distribution integrated over energy loss >',ecut,
     &'     theta_lab    dsig_lab/dO (mb/sr)'
        write(98,'(a,f6.1/a)')
     &' Angular distribution integrated over energy loss >',ecut,
     &'     theta_lab    dsig_lab/dO (mb/sr)'
        angtot = 0.d0
        do nn=1,181
         ang = angtet(nn)*pi/180.d0
         if(nn.eq.1)then
          angp = 0.d0
         else
          angp = angtet(nn-1)*pi/180.d0
         endif
         if(angdis(nn).gt.0.d0)then
          write(*,'(f13.2,1pe23.3)')angtet(nn),angdis(nn)
          write(98,'(f13.2,1pe23.3)')angtet(nn),angdis(nn)
         endif
         angtot = angtot + angdis(nn)*2.d0*pi*sin(ang)*(ang-angp) 
        enddo
        write(*,'(a,1pd12.4,a)')' Total cross section = ',angtot,' mb'
        write(98,'(a,1pd12.4,a)')' Total cross section = ',angtot,' mb'
       endif
      endif

c  printout of total cross sections
      write(*,'(/a,e12.4,a)')
     &' Total grazing cross section = ',cross,'  mbarn'
      write(98,'(/a,e12.4,a)')
     &' Total grazing cross section = ',cross,'  mbarn'
c      if(tres.gt.accur)then
c       write(*,'(/a,e12.4,a)')
c     &' "Total grazing cross section" = ',crossres,'  mbarn'
c       write(*,'(a,f6.2,a/a,f6.1,a)')
c     &'       and with energy resolution(FWHM) = ',tres,' MeV',
c     &'       for scattering angles larger than ',scatmin,' degrees'
c       write(98,'(/a,e12.4,a)')
c     &' "Total grazing cross section" = ',crossres,'  mbarn'
c       write(98,'(a,f6.2,a/a,f6.1,a)')
c     &'       and with energy resolution(FWHM) = ',tres,' MeV',
c     &'       for scattering angles larger than ',scatmin,' degrees'
c      endif
      sigtot = sigcaptn + cross
c      sigtotlim = sigcaptn + crossres
      write(*,'(a,e12.4,a)')
     &' Total capture cross section = ',sigcaptn,' mbarn'
      write(98,'(a,e12.4,a)')
     &' Total capture cross section = ',sigcaptn,' mbarn'
      write(*,'(a,e12.4,a)')
     &' Total reaction cross section = ', sigtot,' mbarn'
      write(98,'(a,e12.4,a)')
     &' Total reaction cross section = ', sigtot,' mbarn'
c      if(tres.gt.accur)then
c       write(*,'(/a,e12.4,a)')
c     &' "Total reaction cross section" = ',sigtotlim,'  mbarn'
c       write(*,'(a,f6.2,a/a,f6.1,a)')
c     &' and with energy resolution = ',tres,' MeV',
c     &'       for scattering angles larger than ',scatmin,' degrees'
c       write(98,'(/a,e12.4,a)')
c     &' "Total reaction cross section" = ',sigtotlim,'  mbarn'
c       write(98,'(a,f6.2,a/a,f6.1,a)')
c     &'       and with energy resolution = ',tres,' MeV',
c     &'       for scattering angles larger than ',scatmin,' degrees'
c      endif

c  Printout for each trajectory found in:
      if(ntab.eq.1)then
       write(*,*)
     &' -------------------------------------------------------------'
       write(98,*)
     &' -------------------------------------------------------------'
       if(leg3.eq.1.or.leg3.eq.4)then
        do n=1,leg3out
         nn=20+n    
         name = outleg3(n)
         if(version.eq.'vms  '.or.version.eq.'unix* ')then
          filename1 = 'graz_'//head//'_'//name//'.dat'
          call condens_s(filename1,filename,le)
          close(unit=nn)
         else
          write(filename,'(a4,i2,a4)') 'FOR0',nn,'.DAT'
         endif
         if(outleg3(n).eq.'edist')then
          write(*,'(a/a,a)')
     &' Kinetic energy distributions for each trajectory are stored',
     &'                          in ',filename
          write(98,'(a/a,a)')
     &' Kinetic energy distributions for each trajectory are stored',
     &'                          in ',filename
         elseif(outleg3(n).eq.'mavdis')then
          write(*,'(a/a,a)')
     &' Distributions in transferred angular momentum for each',
     &'    trajectory are stored in ',filename
          write(98,'(a/a,a)')
     &' Distributions in transferred angular momentum for each',
     &'    trajectory are stored in ',filename
         elseif(outleg3(n).eq.'driftn')then
          write(*,'(a/a,a)')
     &' Distributions in N for each trajectory are stored in',
     &'                             ',filename
          write(98,'(a/a,a)')
     &' Distributions in N for each trajectory are stored in',
     &'                             ',filename
         elseif(outleg3(n).eq.'driftz')then
          write(*,'(a/a,a)')
     &' Distributions in Z for each trajectory are stored in',
     &'                             ',filename
          write(98,'(a/a,a)')
     &' Distributions in Z for each trajectory are stored in',
     &'                             ',filename
         elseif(outleg3(n).eq.'scatang')then
          write(*,'(a/a,a)')
     &' Angular distribution for each trajectory in center of mass',
     &'   are stored in             ',filename
          write(98,'(a/a,a)')
     &' Angular distribution for each trajectory in center of mass',
     &'   are stored in             ',filename
         elseif(outleg3(n).eq.'scatlab')then
          write(*,'(a/a,a)')
     &' Angular distribution for each trajectory in laboratory system',
     &'   are stored in             ',filename
          write(98,'(a/a,a)')
     &' Angular distribution for each trajectory in laboratory system',
     &'   are stored in             ',filename
         elseif(outleg3(n).eq.'dsigde')then
          write(*,'(a/a,a)')
     &' Double differential cross sections for each trajectory',
     &'   in center of mass are stored in ',filename
          write(98,'(a/a,a)')
     &' Double differential cross sections for each trajectory',
     &'   in center of mass are stored in ',filename
         elseif(outleg3(n).eq.'dsiglab')then
          write(*,'(a/a,a)')
     &' Double differential cross sections for each trajectory',
     &'   in laboratory system are stored in ',filename
          write(98,'(a/a,a)')
     &' Double differential cross sections for each trajectory',
     &'   in laboratory system are stored in ',filename
         endif
        enddo
       elseif(leg3.eq.2)then
        leg3out = 1
        if(version.eq.'vms  '.or.version.eq.'unix* ')then
         efilename1 = 'graz_'//head//'_edista.dat'
         call condens_s(efilename1,efilename,le)
         close(unit=21)
        else
         efilename = 'FOR021.DAT'
        endif
        write(*,'(a/a,a40)')
     &' Excitation energy distribution of ejectile stored in ',
     &'            ',efilename
        write(98,'(a/a,a40)')
     &' Excitation energy distribution of ejectile stored in ',
     &'            ',efilename
       elseif(leg3.eq.3)then
        leg3out = 1
        if(version.eq.'vms  '.or.version.eq.'unix*')then
         efilename1 = 'graz_'//head//'_edistaa.dat'
         call condens_s(efilename1,efilename,le)
         close (unit=21)
        else
         efilename = 'FOR021.DAT'
        endif
        write(*,'(a/a,a40)')
     &' Excitation energy distribution of recoil stored in ',
     &'             ',efilename
        write(98,'(a/a,a40)')
     &' Excitation energy distribution of recoil stored in ',
     &'             ',efilename
       endif
      endif       
c  and for total:
      if(leg3.eq.1.or.leg3.eq.4)then
       if(version.eq.'vms  '.or.version.eq.'unix*')then
        efilename1 = 'graz_'//head//'_tdisttot.dat'
        call condens_s(efilename1,efilename,le)
        close (unit=95)
       else
        efilename = 'FOR095.DAT'
       endif
       write(*,'(a/,a40)')
     &' Total kinetic energy distribution stored in ',efilename
       write(98,'(a/,a40)')
     &' Total kinetic energy distribution stored in ',efilename
      endif
      if(leg3.eq.2)then
       if(version.eq.'vms  '.or.version.eq.'unix* ')then
        efilename1 = 'graz_'//head//'_edistatot.dat'
        call condens_s(efilename1,efilename,le)
        close (unit=95)
       else
        efilename = 'FOR095.DAT'
       endif
       write(*,'(a/,a40)')
     &' Total energy distribution of ejectile stored in ',efilename
       write(98,'(a/,a40)')
     &' Total energy distribution of ejectile stored in ',efilename
      endif
      if(leg3.eq.3)then
       if(version.eq.'vms  '.or.version.eq.'unix*')then
        efilename1 = 'graz_'//head//'_edistaatot.dat'
        call condens_s(efilename1,efilename,le)
        close (unit=95)
       else
        efilename = 'FOR095.DAT'
       endif
       write(*,'(a/,a40)')
     &' Total energy distribution of recoil stored in ',efilename
       write(98,'(a/,a40)')
     &' Total energy distribution of recoil stored in ',efilename
      endif
c 
      ierr = 0
2     continue
      return
      end
