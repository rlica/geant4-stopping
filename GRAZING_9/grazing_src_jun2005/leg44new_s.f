c      Subroutine LEG 44 for program grazing including inelastic
c      At exit ierr = 1 for error
c              ierr = 2 for warning
c              ierr = 3 for differential cross section
c              ierr = 4 for total cross section and new value of r
c              ierr = 0 for finished integration
c---------------------------------------------------------------------
      subroutine leg44_s(nestep,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
c
c      list,ienergy,irangen and irangep are defined in grazing.icl
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
      common/tot/in00,ip00,i0nmin,i0nmax,i0zmin,i0zmax
      common/scatang/delw,w0
      common/angdis/deltet
      common/leg2/sigcaptn
      common/leg3/delintloss,rcapt,alcapt
      common/ebe23/eex(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/f0muin/ff0muin(2,3,2,4),f0inm(2,3,2)
      common/giantwid/giantwid(3)
      common/leg4/dsigz(irangep),dsigzlab(irangep),
     &ppnz(irangen,irangep),iselectn(list),iselectz(list),
     &dsig(irangen,irangep),dsiglab(irangen,irangep),
     &inmin,inmax,izmin,izmax,leg4
      common/char/elemp,elemt
c
       dimension thetnz(list,ienergy,101),thetnzlab(list,ienergy,101)
     &,dsigr(list,ienergy,101),rstep(101),eloss(ienergy),tloss(ienergy)
     &,tnzloss(list,ienergy),delel(list),mnzmin(list),mnzmax(list)
     &,mmaxmax(list),mminmin(list),iml(list),istore(5)
     &,irstmax(list)
      dimension thav(ienergy,101),thavlab(ienergy,101),dsr(ienergy,101)
     &,angt(181),angd(181),angtet(list,181),angdis(list,181)
      dimension tablenz(10,list,11,ienergy),tableh(10,list,11,3)
     &,jma(list)
      dimension sig(list)
     &,signzetot(list,ienergy),signzeapptot(list,ienergy)
     &,ppe(ienergy),pnze(list,ienergy),pnzeapp(list,ienergy)
     &,ppapp(ienergy),avmtdis(ienergy),avmnztdis(list,ienergy)
     &,changerdis(ienergy),dsignz(list,ienergy),dsignzlab(list,ienergy)
c      common/char/elemp,elemt
      logical*1 edistout,mavtout,scatout,wilczynski
      character*5 elemp,elemt
      character*5 elempp,elemtt
      character*5 intp
      character*3 out,what
      character*7 outleg4(10)
      character*7 outlegg
      character*7 name
      character*7 xname(5)
      character*40 filename,filename1
      character*40 efilename,efilename1
      character*60 addtext
      character*16 head
      character*5 version
      common/vers/version
c
      iprep = ierr
      if(istep.ne.1)goto 12
      if(nestep.ne.1)goto 11
      ecut = -20.d0
      wilczynski = .false.
      edistout = .false.
      mavtout = .false.
      scatout = .false.
c
c  Total cross sections
10    continue 
       write(*,'(/a/a/a/a/a,f7.3,a)')
     &' You have decided to calculate total cross sections',
     &'  for selected partitions as functions of total kinetic energy',
     &'  i.e. cross sections integrated over angles',
     &'  (and integrated over transferred L)',
     &' Stepwidth in r is delr =',delr,' fm'
       if(intp.eq.'inter')write(*,'(a/a/a)')
     &' You may obtain a table to make a Wilczynski plot',
     &'  for each partition (to be specified later)',
     &' Answer yes or no'
       read(i5,'(a)')out
       if(out.eq.'yes')then
        wilczynski = .true.
        scatout = .true.
       elseif(out.ne.'no ')then
        write(*,*)' Try again!'
        goto 10
       endif
       addtext = '       !wilczynski plot?'
       if(intp.eq.'inter')write(10,'(a,a)')out,addtext
       if(version.eq.'vms  '.or.version.eq.'unix ')then
        if(intp.eq.'inter')write(*,'(a/a)')
     &' You may also obtain the results for each trajectory',
     &' stored in grazing_out.dat file. Enter yes, or no if not.'
       else
        if(intp.eq.'inter')write(*,'(a/a)')
     &' You may also obtain the results for each trajectory',
     &' stored in FOR098.DAT file. Enter yes, or no if not.'
       endif
       read(i5,'(a)')out
       addtext = '    !additional output?'
       if(intp.eq.'inter')write(10,'(a,a)')out,addtext
       if(out.eq.'no ')then
        leg4out = 0
       elseif(out.eq.'yes')then
        leg4out = 1
        edistout = .true.
        outleg4(1) = 'edist'
307     if(intp.eq.'inter')write(*,'(a/a/a/a)')
     &' You now obtain for each value of r:',
     &'  The elastic scattering angle,',
     &'  the elastic cross section divided by Rutherford,',
     &'  and the kinetic energy distribution.'
       endif      
30     if(intp.eq.'inter')write(*,'(7(a/))')
     &' Besides the energy distribution (edist)',
     &'  You may calculate, and display for each r, variables from',
     &'  the following list (in paranthesis their significance):',
     &' mavdis (mean ang. momentum loss for given energy loss)',
     &' scatang (mean scattering angle for given energy loss)',
     &' dsigde (double differential cross section)',
     &' Enter name of variable, no if no more'
       read(i5,'(a)')outlegg
       if(outlegg.eq.'mavdis')then
        mavtout = .true.
       elseif(outlegg.eq.'scatang')then
        scatout = .true.
        if(lab.eq.1)then
         outlegg = 'scatlab'
        elseif(lab.eq.2)then
         leg4out = leg4out+1
         outleg4(leg4out) = 'scatlab'
        endif
       elseif(outlegg.eq.'dsigde')then
        scatout = .true.
        if(lab.eq.1)then
         outlegg = 'dsiglab'
        elseif(lab.eq.2)then
         leg4out = leg4out+1
         outleg4(leg4out) = 'dsiglab'
        endif
       elseif(outlegg.ne.'no')then
        write(*,*)' Wrong spelling, try again!'
       endif
       addtext = '      !additional quantities'
       if(intp.eq.'inter')write(10,'(a,a)')outlegg,addtext
       if(outlegg.eq.'no')then
        goto 309
       else
        leg4out= leg4out+1
        outleg4(leg4out) = outlegg
        goto 30
       endif
309    if(leg4out.eq.0)goto 72
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
       if(ntab.eq.1)then
        lout = leg4out
        do lo = 1,leg4out
         istore(lo) = 20 + (lo-1)*10
         xname(lo) = outleg4(lo)
        enddo
       else
        lout = 0
       endif
72     continue
       ierr = 0
       call selectnz_s(lnzmax,istore,xname,lout,accur,ierr)
       if(ierr.eq.1)goto 2
305    if(intp.eq.'inter')write(*,'(4(/a))')
     &' You may calculate the energy distributions with a finite',
     &'  energy resolution tres(FWHM) and cut-off.',
     &'  The cut-off may be negative (default = -20MeV)',
     &' Enter tres and ecut in MeV (0,-20 if not)'
       read(i5,*,err=1301)tres,ecut
       icheck = 0
       goto 1311
1301   write(*,*)' you must provide two numbers'
       icheck = icheck + 1
       if(icheck.gt.5)goto 2
       goto 305
1311   addtext = '    !energy resolution(FWHM)and cutoff in MeV'
       if(intp.eq.'inter')write(10,'(2f7.4,a)')tres,ecut,addtext
       if(tres.lt.0.d0)then
        write(*,*)' Try again'
        goto 305
       endif
       if(tres.gt.accur)write(*,'(a,f7.3,/a,f7.3)')
     &' Energy distributions are calculated with resolution = ',tres,
     &'  and cut-off in kinetic energy loss ',ecut
       do l=1,lnzmax
        if(iselectn(l).eq.0.and.iselectz(l).eq.0)then
         if(tres.lt..2d0)then
         if(intp.eq.'inter')write(*,'(2(/a))')
     &' Since you have selected delta_N=delta_Z=0,',
     &'   tres will be set equal to .2 for this partition'
         endif
314      if(intp.eq.'inter')write(*,'(a/a/a/a)')
     &' Since the total cross section diverges you should enter',
     &'  the minimum elastic scattering angle to be included', 
     &'  (in degrees) for large impact parameters,',
     &'  if negative, this angle is measured from the rainbow angle'
         read(i5,*,err=1104)scatminn
         icheck = 0
         goto 1114
1104     write(*,*)' you must provide a number'
         icheck = icheck + 1
         if(icheck.gt.5)goto 2
         goto 314
1114     if(rainbcdg.gt.130.d0)then
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
        endif
       enddo
c
c  printout for total in:
       if(version.eq.'vms  ')then
        efilename1 = 'graz_'//head//'_tdisttot.dat'
        call condens_s(efilename1,efilename,le)
        open (unit=95,file=efilename,status='unknown')
       else
        efilename = 'FOR095.DAT'
       endif
c
c  preparation for total cross sections
11    continue
      cross = 0.d0
      crossres = 0.d0
      do l=1,list
         irstmax(l)=0.
      enddo
      do l=1,lnzmax
        sig(l) = 0.d0
        do  mm=1,ienergy
         signzetot(l,mm) = 0.d0
         signzeapptot(l,mm)=0.d0
        enddo
       enddo
       do l=1,lnzmax
        irstmax(l)= 50
        mmaxmax(l) = 1
        mminmin(l) = ienergy
       enddo
c      endif
c
c   calculate energy distributions
12    continue
      rstep(istep) = r
      ierr1 = 1
      ierr2 = 1
      dsigdr=10.d0*r*(ecm-za*zaa*esq/2.d0/r
     &-unuc(r,ierr1)-r/2.d0*dunuc(r,ierr2))/ecm
      if(ierr1.eq.1.or.ierr2.eq.1)then
       write(*,*)' error in leg44_s '
       ierr = 1
       goto 2
      endif
      write(*,'(a,f7.2,a)')'  Calculating for r = ',r,' fm'
      do l = 1,lnzmax
       timeedist1 = X05BAF()
       if(istep.gt.irstmax(l))goto 441
       isn = iselectn(l)
       isz = iselectz(l)
       eqns = qns
       if(isn.gt.0)eqns=-qnp
       eqps = qps
       if(isz.gt.0)eqps=-qpp
       write(*,'(a,i3,a,i3)')
     &'    distribution in energy for deln = ',isn,' and delz = ',isz
       delta = tres/2.35482d0
       if(isn.eq.0.and.isz.eq.0.and.tres.lt..2d0)then
        tres = .2d0
        delta = .2/2.35482d0
       endif
       do mi = 1,ienergy
        ppe(mi) = 0.d0
        ppapp(mi) = 0.d0
        pnze(l,mi) = 0.d0
        pnzeapp(l,mi) = 0.d0
        avmnztdis(l,mi) = 0.d0
       enddo
       eavpr = eavdist(in00+isn,ip00+isz)/2.d0
       wid2epr = wid2adist(in00+isn,ip00+isz)
       pnz = ppnz(in00+isn,ip00+isz)
       if(isn.eq.0.and.isz.eq.0)pnz=pnz+p0
       if(pnz.lt.accur**2)then
        irstmax(l) = istep-1
        goto 441
       endif
       widepr = sqrt(wid2epr)
       if(istep.eq.1)then
        deleloss = 2.d0*(widepr+delta)*accur**(.25d0)
        delel(l) = dble(1.d0+int(10.*deleloss))/10.d0
        iml(l) = 0
       endif
       dele=delel(l)
       imm = iml(l)
       ierr = iprep
       call nztdis_s(delta,dele,accur,imm,mmmax,mmmin
     &,isn,isz,eloss,ppe,ierr)
       iml(l)=imm
       if(ierr.eq.1)goto 2
       if(ierr.gt.1)write(*,'(a,f10.7)')'     for r = ',r
       ierr = iprep
       eavap = eavdist(isn+in00,isz+ip00)
       wid2ap = wid2dist(isn+in00,isz+ip00)
       p0ap = 0.d0
       if(isn.eq.0.and.isz.eq.0)then
        eavap = eavdist00
        wid2ap = wid2dist00
        p0ap = p0tr
       endif
       deltae = delta
       iqq = 0
       call edistapp_s(iqq,deltae,eavap,wid2ap,p0ap,dele,accur,
     &imm,mmmax,mmmin,eloss,ppapp,ierr)
       if(ierr.gt.1)write(*,*)'     for r = ',r
       do mm = mmmin,mmmax
        write(*,'(a,f8.3,2(a,f7.5))')
     &' eloss = ',eloss(mm),' ppe = ',ppe(mm),' ppapp = ',ppapp(mm)
       enddo
       do ii = 1,ienergy
        tnzloss(l,ii) = eloss(ii)+isn*eqns+isz*eqps
        pnze(l,ii) = ppe(ii)
        pnzeapp(l,ii) = ppapp(ii)
       enddo
c
       do mm = 1,ienergy
        if(istep.eq.1)dsigr(l,mm,1) = 0.d0
        dsigr(l,mm,istep) = 0.d0
       enddo
       do mm=mmmin,mmmax
        dsigr(l,mm,istep) = dsigdr*pnze(l,mm)*pnz*ppscat
       enddo
       if(mavtout.or.scatout)then
       write(*,*)' calculating mean angular momentum for each energy'
        ierr = iprep
        acc = accur
        call mavnzdis_s(delta,dele,acc,imm,mmmax,mmmin
     &,isn,isz,eloss,avmtdis,abserr,ierr)
        if(ierr.eq.1)goto 2
        if(ierr.gt.1)write(*,'(a,f10.7)')'     for r = ',r
        do mm=mmmin,mmmax
         if(pnze(l,mm).ne.0.d0)then
          avd = avmtdis(mm)/pnze(l,mm)
          avmnztdis(l,mm) = avd
          if(al-avd.lt.0.1d0.or.abserr/(1.d0+al).gt.pnze(l,mm))then
           eavp = eavdist(isn+in00,isz+ip00)+eavina+eavinaa
           if(eloss(mm).lt.eavp)then
            mmmin = mm+1
            do mi=1,mmmin-1
             avmnztdis(l,mi) = 0.d0
            enddo
            goto 25
           elseif(eloss(mm).gt.eavp)then
            mmmax = mm-1
            do mi=mmmax+1,ienergy
             avmnztdis(l,mi) = 0.d0
            enddo
            goto 26
           endif
          endif
         endif
25       continue           
        enddo
26      continue
       endif 
       if(scatout)then
        write(*,*)' calculating scattering angle for each energy'
        if(iprep.eq.1)then
         mmmin = imm -3
         mmmax = imm
        endif
        do ii = 1,ienergy
         thetnz(l,ii,istep) = 0.d0
         thetnzlab(l,ii,istep) = 0.d0
        enddo
        do mm = mmmin,mmmax
         ierr = 0
         dela = float(isn+isz)
         delz = float(isz)
         delecm = -tnzloss(l,mm)
         delal = - avmnztdis(l,mm)
         w00 = w0
         dellw = delw
         call scatangleavx_s(al,ph,dph,rturnex,delecm,delal,dela,delz
     &,w00,dellw,ierr)
         if(ierr.eq.1)write(*,*)'  Error for tloss = ',tnzloss(l,mm)
         changerdis(mm) = rturnex-r
         thetaver = pi-ph
         dthetav = dph
         thetnz(l,mm,istep) = thetaver
         sinth = sin(theta)
         if(sinth.gt.accur)then
          dsignz(l,mm) = cs*pnze(l,mm)*sinth*dthet/sin(thetaver)/dthetav
         else
          dsignz(l,mm) = cs*pnze(l,mm)*dthet/dthetav
         endif
c  In laboratory
         taua = sqrt(ecm/(ecm+delecm)*aa/aaa*(aa+dela)/(aaa-dela))
         thetrlab = atan(sin(thetaver)/(cos(thetaver)+taua))
         if(thetrlab.lt.0.d0)thetrlab=thetrlab + pi
         thetnzlab(l,mm,istep) = thetrlab
         if(abs(sin(thetaver)).lt.accur)then
          dsignzlab(l,mm) = dsignz(l,mm)*(1.d0-taua)**2
         else
          dsignzlab(l,mm) = dsignz(l,mm)*(sin(thetaver)
     &/sin(thetrlab))**2/abs(cos(thetaver-thetrlab))
         endif
27       continue         
        enddo
       endif
       mnzmin(l)=mmmin
       mnzmax(l)=mmmax
c  mnzmax(l) = mmmax for given l for printout
       if(mnzmax(l).gt.mmaxmax(l))mmaxmax(l) = mnzmax(l)
       if(mnzmin(l).lt.mminmin(l))mminmin(l) = mnzmin(l)
c  mmaxmax(l) = maximum mnzmax for all r-values and given l for total
441    continue
       timeedist2 = X05BAF()
       timeedist = timeedist2 - timeedist1
       write(*,'(a,i3,a,f7.2,a/a,i3,a)')
     &' CPU time for',l,'. energy distribution =',timeedist,' sec.',
     &'   You ordered',lnzmax,' distributions'
      enddo
      if(leg4out.eq.0)then
       goto 15 
      else
       write(*,'(/a,f7.2,a,f5.2,a)')
     &'   Results for l = ',al,' ( r = ', r,' fm)'
       write(98,'(/a,f7.2,a,f5.2,a)')
     &'   Results for l = ',al,' ( r = ', r,' fm)'
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
       do l=1,lnzmax
c      print*,' irstmax(l) = ',l,irstmax(l)
        if(istep.gt.irstmax(l))goto 112
        isn = iselectn(l)
        isz = iselectz(l)
        if(lab.eq.0)then
         write(98,'(a/a,i3,a,i3/a,f7.3,a)')
     &' Distribution in total kinetic energy loss, etc.',
     &'  in center of mass  for deln = ',isn,' and delz = ',isz,
     &'  with energy  resolution (FWHM) =',tres,' MeV'
         write(98,'(a)')
     &'  T_loss  Probability    <L>    theta   dsig/dO/dT '
         write(*,'(a/a,i3,a,i3/a,f7.3,a)')
     &' Distribution in total kinetic energy loss, etc.',
     &'  in center of mass  for deln = ',isn,' and delz = ',isz,
     &'  with energy  resolution (FWHM) =',tres,' MeV'
         write(*,'(a)')
     &'  T_loss  Probability    <L>    theta   dsig/dO/dT '
        elseif(lab.eq.1)then
         write(98,'(a/a,i3,a,i3/a,f7.3,a)')
     &' Distribution in total kinetic energy loss, etc.',
     &'  in lab system for deln = ',isn,' and delz = ',isz,
     &'  with energy  resolution (FWHM) =',tres,' MeV'
         write(98,'(a)')
     &'  T_loss  Probability    <L>  thet_lab dsig_lab/dO/dT '
         write(*,'(a/a,i3,a,i3/a,f7.3,a)')
     &' Distribution in total kinetic energy loss, etc.',
     &'  in lab system for deln = ',isn,' and delz = ',isz,
     &'  with energy  resolution (FWHM) =',tres,' MeV'
         write(*,'(a)')
     &'  T_loss  Probability    <L>  thet_lab dsig_lab/dO/dT '
        endif
      print*,' mnzmin(l),mnzmax(l) = ',mnzmin(l),mnzmax(l)
      print*,' mminmin(l),mmaxmax(l) = ',mminmin(l),mmaxmax(l)
        do mm = mnzmin(l),mnzmax(l)
         if(thetnz(l,mm,istep).eq.0.d0.or.
     &pnze(l,mm).lt.accur/100.d0)goto 111
         thedg = thetnz(l,mm,istep)*180.d0/pi
         thelbdg = thetnzlab(l,mm,istep)*180.d0/pi
         tnze=tnzloss(l,mm)
         if(lab.eq.0)then
          write(98,'(f9.1,f10.4,f10.2,f10.2,e12.3)')
     & tnze,pnze(l,mm),avmnztdis(l,mm),thedg,dsignz(l,mm)
          write(*,'(f9.1,f10.4,f10.2,f10.2,e12.3)')
     & tnze,pnze(l,mm),avmnztdis(l,mm),thedg,dsignz(l,mm)
         elseif(lab.ne.0)then
          write(98,'(f9.1,f10.4,f10.2,f10.2,e12.3)')
     & tnze,pnze(l,mm),avmnztdis(l,mm),thelbdg,dsignzlab(l,mm)
          write(*,'(f9.1,f10.4,f10.2,f10.2,e12.3)')
     & tnze,pnze(l,mm),avmnztdis(l,mm),thelbdg,dsignzlab(l,mm)
         endif
111      continue
        enddo
112     continue
       enddo
      endif
c
c  Printout of tables
15    if(ntab.eq.1)then
       if(mod(istep,10).eq.1)then
        do m=1,11
         do leo = 1,leg4out
          do l = 1,lnzmax
           tableh(leo,l,m,1)=0.d0
           tableh(leo,l,m,2)=0.d0
           tableh(leo,l,m,3)=0.d0
           do n=1,ienergy
            tablenz(leo,l,m,n)=0.d0
           enddo
          enddo
         enddo
        enddo
       endif
       do l = 1,lnzmax
        isn = iselectn(l)
        isz = iselectz(l)
        do leo = 1,leg4out
         nn = 20+10*(leo-1)+l
         if(istep.eq.1) then
          if(outleg4(leo).eq.'edist')then
           write(nn,'(a/a,i3,a,i3)')
     &' Distributions in total kinetic energy loss for each trajectory',
     &'  for deln = ',isn,' and delz = ',isz
          elseif(outleg4(leo).eq.'mavdis')then
           write(nn,'(a/a,i3,a,i3)')
     &' Distributions in transferred ang. momentum for each trajectory',
     &'  for deln = ',isn,' and delz = ',isz
          elseif(outleg4(leo).eq.'scatang')then
           write(nn,'(a/a,i3,a,i3)')
     &' Distributions in cm scattering angle for each trajectory',
     &'  for deln = ',isn,' and delz = ',isz
          elseif(outleg4(leo).eq.'scatlab')then
           write(nn,'(a/a,i3,a,i3)')
     &' Distributions in lab scattering angle for each trajectory',
     &'  for deln = ',isn,' and delz = ',isz
          elseif(outleg4(leo).eq.'dsigde')then
           write(nn,'(a/a/a,i3,a,i3)')
     &' Double differential cross section for each trajectory',
     &'  in mb/MeV/sr in center of mass',
     &'  for deln = ',isn,' and delz = ',isz
          elseif(outleg4(leo).eq.'dsiglab')then
           write(nn,'(a/a/a,i3,a,i3)')
     &' Double differential cross section for each trajectory',
     &'  in mb/MeV/sr in laboratory system',
     &'  for deln = ',isn,' and delz = ',isz
          endif
          write(nn,'(a,f7.3,a)')
     &'     with ',tres,' MeV energy resolution (FWHM)'
          write(nn,'(a)')
     &'----------------------------------------------------------------'
         endif
         jcol = mod(istep-1,10)+2
         if(istep.le.irstmax(l)-10)then
          jma(l) = 11
         elseif(istep.gt.irstmax(l)-10.and.istep.le.irstmax(l))then
          jma(l) = mod(irstmax(l),10)+1 
         else
          goto 113
         endif
         tableh(leo,l,jcol,1)=r
         tableh(leo,l,jcol,2)=thet
         if(outleg4(leo).eq.'scatlab'.or.outleg4(leo).eq.'dsiglab')
     &tableh(leo,l,jcol,2)=thetlabdeg
         tableh(leo,l,jcol,3)=sigreac0*ppnz(in00+isn,ip00+isz)
         do mm=mminmin(l),mmaxmax(l)
          tablenz(leo,l,1,mm)= tnzloss(l,mm)
          if(outleg4(leo).eq.'edist')then
           tablenz(leo,l,jcol,mm)= pnze(l,mm)
          elseif(outleg4(leo).eq.'mavdis')then
           tablenz(leo,l,jcol,mm)= avmnztdis(l,mm)
          elseif(outleg4(leo).eq.'scatang')then
           thedg = thetnz(l,mm,istep)*180.d0/pi
           tablenz(leo,l,jcol,mm)= thedg
          elseif(outleg4(leo).eq.'scatlab')then
           thelbdg = thetnzlab(l,mm,istep)*180.d0/pi
           tablenz(leo,l,jcol,mm)= thelbdg
          elseif(outleg4(leo).eq.'dsigde')then
           tablenz(leo,l,jcol,mm)= dsignz(l,mm)
          elseif(outleg4(leo).eq.'dsiglab')then
           tablenz(leo,l,jcol,mm)= dsignzlab(l,mm)
          endif
         enddo
113      if(jcol.eq.11)then
          jmm = jma(l)
          write(nn,'(a,10f11.2)')
     &'      r   ',(tableh(leo,l,n,1),n=2,jmm)
          write(nn,'(a,10f11.2)')
     &'    theta ',(tableh(leo,l,n,2),n=2,jmm)
          write(nn,'(a,10e11.3)')
     &'  dsig_r/dl',(tableh(leo,l,n,3),n=2,jmm)
          write(nn,'(/)')
          if(outleg4(leo).eq.'dsigde'.or.outleg4(leo).eq.'dsiglab')then
           do mm=mminmin(l),mmaxmax(l)
            write(nn,'(f10.3,10e11.3)')
     &tablenz(leo,l,1,mm),(tablenz(leo,l,n,mm),n=2,jmm)
           enddo
          else
           do mm=mminmin(l),mmaxmax(l)
           write(nn,'(f10.3,10f11.5)')
     &tablenz(leo,l,1,mm),(tablenz(leo,l,n,mm),n=2,jmm)
           enddo
          endif
          write(nn,*)
     &'----------------------------------------------------------------'
         endif
        enddo
       enddo
      endif 
c
c  integration of cross sections (not accurate for irstep.eq.0)
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
      do l=1,lnzmax
       n=in00+iselectn(l)
       m=ip00+iselectz(l)
       pppnz = ppnz(n,m)
       sig(l)= sig(l)+factor*pppnz*ppscat
       do mm=mnzmin(l),mnzmax(l)
        signzetot(l,mm) = signzetot(l,mm) + factor*pppnz*pnze(l,mm)*
     &ppscat
        signzeapptot(l,mm) = signzeapptot(l,mm) + factor*pppnz*
     & pnzeapp(l,mm)*ppscat
       enddo
      enddo
      if(istep+1.gt.100)then
       write(*,*)' istep+1.gt.100 in leg4'
       goto 14
      endif
      if(irstep.ne.0)irstep = irstep+1
      istep = istep+1
c  Check for finished integration
      if(iprep.eq.1.and.istep.ge.3)goto 14
      if(abs(1.-p0).lt.accur)goto 14
      alpr = al
      ierr = 4
      goto 2
14    do l=1,lnzmax
       if(irstmax(l).eq.50)irstmax(l)=istep
      enddo
c
      if(ntab.eq.1.and.jcol.ne.11)then
c  printout of remaining differential
       do l=1,lnzmax
        if(istep.gt.irstmax(l))goto 114
        do leo = 1,leg4out
         nn = 20+10*(leo-1)+l
         write(nn,'(a,10f11.2)')
     &'      r   ',(tableh(leo,l,n,1),n=2,jcol)
         write(nn,'(a,10f11.2)')
     &'    theta ',(tableh(leo,l,n,2),n=2,jcol)
         write(nn,'(a,10e11.3)')
     &'  dsig_r/dl',(tableh(leo,l,n,3),n=2,jcol)
         write(nn,'(/)')
         if(outleg4(leo).eq.'dsigde'.or.outleg4(leo).eq.'dsiglab')then
          do mm=mminmin(l),mmaxmax(l)
           write(nn,'(f10.3,10e11.3)')
     &tablenz(leo,l,1,mm),(tablenz(leo,l,n,mm),n=2,jcol)
          enddo
         else
          do mm=mminmin(l),mmaxmax(l)
           write(nn,'(f10.3,10f11.5)')
     &tablenz(leo,l,1,mm),(tablenz(leo,l,n,mm),n=2,jcol)
          enddo
         endif
        enddo
114     continue
       enddo
      endif
c
c  Printout of total cross section distribution
      write(*,*)'------------------------------------------------'
      write(98,*)'------------------------------------------------'
      do l = 1,lnzmax
       isn = iselectn(l)
       isz = iselectz(l)
       write(98,'(a/a,i3,a,i3/a,f7.3,a/a)')
     &' Total kinetic energy loss cross section',
     &'  for deln = ',isn,' and delz = ',isz,
     &'     with ',tres,' MeV energy resolution (FWHM)',
     &' teloss MeV      sig mb/MeV   sig_appr'
       write(*,'(a/a,i3,a,i3/a,f7.3,a/a)')
     &' Total kinetic energy loss cross section',
     &'  for deln = ',isn,' and delz = ',isz,
     &'     with ',tres,' MeV energy resolution (FWHM)',
     &' teloss MeV      sig mb/MeV   sig_appr'
       if(isn.eq.0.and.isz.eq.0)then
        write(98,'(a,f8.2,a)')
     &'     and scattering angles larger than ',scatmin,' degrees'
        write(*,'(a,f8.2,a)')
     &'     and scattering angles larger than ',scatmin,' degrees'
       endif
       do mm=mminmin(l),mmaxmax(l)
        if(signzetot(l,mm).gt.0.d0)then
         write(95,'(f7.1,2(1pd15.4))')
     & tnzloss(l,mm),signzetot(l,mm),signzeapptot(l,mm)
         write(98,'(f7.1,2(1pd15.4))')
     & tnzloss(l,mm),signzetot(l,mm),signzeapptot(l,mm)
         write(*,'(f7.1,2(1pd15.4))')
     & tnzloss(l,mm),signzetot(l,mm),signzeapptot(l,mm)
        endif
       enddo
       write(98,'(a,1pd12.4,a)')' Total cross section = ',sig(l),' mb'  
       write(*,'(a,1pd12.4,a)')' Total cross section = ',sig(l),' mb'  
      enddo
c
22    if(wilczynski)then
222    if(intp.eq.'inter')write(*,'(a/a/a)')
     &' Do you want Wilczynski plot as dsig/dE/dO in mb/MeV/sterad',
     &'   or as dsig/dE/dtheta in mb/Mev/rad?',
     &'     Enter sr or rad.'
       read(i5,'(a)')what
       addtext = ' Wilczynski plot? sterad or rad'
       if(what.ne.'rad'.and.what.ne.'sr ')then
        write(*,*)' try again!'
        goto 222
       endif
       if(intp.eq.'inter')write(10,'(a,a)') what,addtext
       widang = deltet
       do l = 1,lnzmax
        isn = iselectn(l)
        isz = iselectz(l)
        minm = mminmin(l)
        maxm = mmaxmax(l)
        do mm=1,ienergy
         tloss(mm)=tnzloss(l,mm)
         do ir=1,irstmax(l)
          thav(mm,ir)= thetnz(l,mm,ir)
          thavlab(mm,ir)=thetnzlab(l,mm,ir)
          dsr(mm,ir)=dsigr(l,mm,ir)
         enddo
        enddo
        irstma = irstmax(l)
        if(lab.ne.1)then
         if(what.eq.'sr ')then
          la = 0
         elseif(what.eq.'rad')then
          la = 2
         endif
         ierr = 0
         nnfile = 78+2*l
         call wilsh_s(thav,thavlab,dsr,rstep,tloss,widang,ecut
     &,angt,angd,la,nnfile,minm,maxm,irstma,ierr)
         if(ierr.eq.1)goto 2
         do n=1,181
          angtet(l,n) = angt(n)
          angdis(l,n) = angd(n)
         enddo
         if(version.eq.'unix ')then
          dummy=42424.1616
         elseif(version.eq.'vms  '.or.version.eq.'unix ')then
          ierr = la
          call wilsh_nzplot_s(isn,isz,nnfile,ierr)
         endif
         write(*,'(a,f6.1/a,i3,a,i3/a)')
     &' Angular distribution integrated over energy loss >',ecut,
     &'  for the partition ideln = ',isn,' and idelz = ',isz,
     &'     theta_cm    dsig_cm/dO (mb/sr)'
         write(98,'(a,f6.1/a,i3,a,i3/a)')
     &' Angular distribution integrated over energy loss >',ecut,
     &'  for the partition ideln = ',isn,' and idelz = ',isz,
     &'     theta_cm    dsig_cm/dO (mb/sr)'
        angtot = 0.d0
         do nn=2,181
          ang = angtet(l,nn)*pi/180.d0
          angp = angtet(l,nn-1)*pi/180.d0
          if(ang.gt.0.d0)then
           write(*,'(f13.2,1pe22.3)')angtet(l,nn),angdis(l,nn)
           write(96,'(f13.2,1pe22.3)')angtet(l,nn),angdis(l,nn)
           write(98,'(f13.2,1pe22.3)')angtet(l,nn),angdis(l,nn)
          endif
          angtot = angtot + 2.d0*pi*sin(ang)*angdis(l,nn)*(ang-angp)
         enddo
         write(*,'(a,1pd12.4,a)')' Total cross section = ',angtot,' mb'
         write(98,'(a,1pd12.4,a)')' Total cross section = ',angtot,' mb'
        endif
        if(lab.ne.0)then
         if(what.eq.'sr ')then
          la = 1
         elseif(what.eq.'rad')then
          la = 3
         endif
         ierr = 1
         nnfile = 79+2*l
         call wilsh_s(thav,thavlab,dsr,rstep,tloss,widang,ecut
     &,angt,angd,la,nnfile,minm,maxm,irstma,ierr)
         if(ierr.eq.1)goto 2
         do n=1,181
          angtet(l,n) = angt(n)
          angdis(l,n) = angd(n)
         enddo
         if(version.eq.'unix ')then
          dummy=4242.1616
         elseif(version.eq.'vms  '.or.version.eq.'unix ')then
          ierr = la
          call wilsh_nzplot_s(isn,isz,nnfile,ierr)
         endif
         write(*,'(a,f6.1/a,i3,a,i3/a)')
     &' Angular distribution integrated over energy loss >',ecut,
     &'  for the partition ideln = ',isn,' and idelz = ',isz,
     &'     theta_lab    dsig_lab/dO (mb/sr)'
         write(98,'(a,f6.1/a,i3,a,i3/a)')
     &' Angular distribution integrated over energy loss >',ecut,
     &'  for the partition ideln = ',isn,' and idelz = ',isz,
     &'     theta_lab    dsig_lab/dO (mb/sr)'
         angtot = 0.d0
         do nn=2,181
          ang = angtet(l,nn)*pi/180.d0
          angp = angtet(l,nn-1)*pi/180.d0
          if(ang.gt.0.d0)then
           write(*,'(f13.2,1pe22.3)')angtet(l,nn),angdis(l,nn)
           write(96,'(f13.2,1pe22.3)')angtet(l,nn),angdis(l,nn)
           write(98,'(f13.2,1pe22.3)')angtet(l,nn),angdis(l,nn)
          endif
          angtot = angtot + 2.d0*pi*sin(ang)*angdis(l,nn)*(ang-angp)
         enddo
         write(*,'(a,1pd12.4,a)')' Total cross section = ',angtot,' mb'
         write(98,'(a,1pd12.4,a)')' Total cross section = ',angtot,' mb'
        endif
       enddo
      endif
      write(*,'(/a,e12.4,a)')
     &' Total grazing cross section = ',cross,'  mbarn'
      write(98,'(/a,e12.4,a)')
     &' Total grazing cross section = ',cross,'  mbarn'
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
      write(*,*)' -------------------------------------------------'
      write(98,*)' -------------------------------------------------'
      write(*,'(a,/,a40)')
     &' Total kinetic energy distribution is stored in ',efilename
      write(98,'(a,/,a40)')
     &' Total kinetic energy distribution stored in ',efilename
cn      close(unit=15)
      close(unit=95)
      do l=1,lnzmax
       do lo = 1,lout
        ll = istore(lo) + l
        close(unit=ll)
       enddo
      enddo
      ierr = 0
2     return
      end
