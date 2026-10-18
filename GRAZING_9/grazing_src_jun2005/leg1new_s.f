c      Subroutine LEG 1 for program grazing:
c       Average quantities integrated over excitation energy,
c       transferred angular momentum and mass and charge transfer.
c      At exit ierr = 1 for error
c              ierr = 2 for warning
c              ierr = 3 for differential cross section
c              ierr = 4 for total cross section and new value of r
c              ierr = 0 for total cross section finished
c      At entry ierr = 1 indicates preparation of input file
c              (unreliable output)
c      subroutines: condens_s,treslim_s,population_s
c---------------------------------------------------------------------
      subroutine leg1_s(nestep,ierr)
      implicit real*8(a-h,o-z)
      common/za/za,aa,zaa,aaa
      common/ecm/ecm
      common/capt/eb,rbar,i5,intp
c      common:akyuz
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/eav/eav,eavns,eavnp,eavps,eavpp,eavtr,
     &eav0,eava0,eavaa0
      common/mav/avma,avmaa,avmans,avmanp,avmaps,avmapp,
     &avmaans,avmaanp,avmaaps,avmaapp,avmatr,avmaatr,
     &avma0,avmaa0,avm0
      common/wid/wid2e,wid2ns,wid2np,wid2ps,wid2pp,wid2tr,wid2ea,
     &wid2eaa,wid2e0,wid2ea0,wid2eaa0
      common/f0in/ff0in(2,3,2)
      common/inel/p0ina,p0inaa,eavina,eavinaa,wid2ina,wid2inaa,
     &avmina,avminaa,eava00,eavaa00,wid2a00,wid2aa00,inelast
      common/head/head
      common/leg/ecminit,rainbcdg,decm,cs,r,delr,al,theta,
     &cscoul,cslab,thet,thetlabdeg,thetcou,ppcapt,ww,pr0,
     &drifta,driftz,wida,widz,eava,eavaa,tavs,widea,wideaa,avmin,avm,
     &changer,changers,ppscat,totmerg,dthet,
     &deltae,sig0drut,sigreac0,sigreac0lab,cross,alorho,rho,
     &istep,iout,ipopn,norbit,irstep,lab
      common/leg1/excom,dsigcaptn(11),aln(11),homegan,homegap
      common/leg2/sigcaptn
      common/probc/adia,capfac,ecmcapt,widcbar,sigcest,icapt
      common/leg3/delintloss,rcapt,alcapt
      common/scatang/delw,w0
c      character*4 ecmstr
      character*5 intp
      character*16 head
      character*7 name
      character*8 outl
      character*8 outleg1(20)
      character*34 filename,filename1
      character*6 text1
      character*9 text2
      character*29 text
      character*40 addtext
      character*5 version
      character*3 out
      common/vers/version
c
      intprep = ierr
c      print*,' leg1 entered with istep,irstep = ',istep,irstep
      if(istep.ne.1)goto 41
      if(nestep.ne.1)goto 8
c
c  Menu for total cross section
       write(*,'(3(/a)/a,f7.3,a,2(/a))')
     &' You are calculating total cross sections',
     &'  summed over N and Z and integrated over excitation energy',
     &'  and transferred angular momentum.',
     &' Stepwidth in distance of closest approach r is ',delr,' fm',
     &' Standard output is reaction and capture cross sections',
     &' sig_reac and sig_merg'
18     if(intp.eq.'inter')write(*,'(5(/a))')
     &' You have the option to obtain elastic scattering',
     &'  with finite (Gaussian) energy resolution (tres)',
     &'  and finite kinetic energy cutoff (telim)',
     &' Enter width (FWHM) and cutoff in MeV:',
     &'   both may be 0, but if cutoff.ne.0 then choose width.ne.0'
       read(i5,*,err=101)tres,telim
       icheck = 0
       goto 111
101    write(*,*)' you must provide two numbers'
       icheck = icheck + 1
       if(icheck.gt.5)goto 2
       goto 18
111    if(tres.lt.0.)then
        write(*,*)' Try again!'
        goto 18
       elseif(tres.lt.accur.and.abs(telim).gt.accur)then
        write(*,*)' tres should not be zero when tlim.ne.0'
        goto 18
       endif
       addtext = '     ! energy resolution and cutt-off'
       if(intp.eq.'inter')write(10,'(2f8.3,a)')tres,telim,addtext 
       if(tres.gt.accur)then
19      if(intp.eq.'inter')write(*,'(a/a/a,f6.2)')
     &' Since the total cross section diverges you should enter',
     &'  the minimum angle to be included (in degrees),',
     &'  if negative or 0, it is measured from the rainbow angle = '
     &,rainbcdg
        read(i5,*,err=102)scatminn
        icheck = 0
        goto 112
102     write(*,*)' you must provide a number'
        icheck = icheck + 1
        if(icheck.gt.5)goto 2
        goto 19
112     addtext = '   !minimum scattering angle (- if from rainbow)' 
        if(intp.eq.'inter')write(10,'(f9.3,a)')scatminn,addtext
        if(scatminn.lt.accur)then
         scatmin = rainbcdg + scatminn
        else
         scatmin = scatminn
        endif
       endif
c      print*,' scatminn,scatmin,rainbcdg = ',scatminn,scatmin,rainbcdg
45     if(version.eq.'vms  '.or.version.eq.'unix* ')then
        if(intp.eq.'inter')write(*,'(a/a/a)')
     &' You may obtain intermediate results for the contributing',
     &'  trajectories stored in tables and in the grazing_out.dat file',
     &' Enter yes if yes, no if not'
       else
        if(intp.eq.'inter')write(*,'(a/a/a)')
     &' You may obtain intermediate results for the contributing',
     &'  trajectories stored in tables and in the FOR098.DAT file',
     &' Enter yes if yes, no if not'
       endif
       read(i5,'(a)')out
       addtext = '    ! variable for printout'
       if(intp.eq.'inter')write(10,'(a,a)')out,addtext
       if(out.eq.'no ')then
        leg1out = 0
       elseif(out.eq.'yes')then
        leg1out = 1
107     if(intp.eq.'inter')write(*,'(a,23(/a))')     
     &' You may chose among the variables from the following list',
     &' (in paranthesis its significance):',
     &' thet (elastic scattering angle)',
     &' pelas (elastic probability)',
     &' ww (absorptive potential)',
     &' sigdrut (elastic divided by Rutherford)',
     &' ppcapt (capture probability)',
     &' dsigcap (l-distribution of capture)',
     &' yrast (approximate yrast-line for compound nucleus)',
     &' preac (reaction probability)',
     &' thetav (average scattering angle for reactions)',
     &' sigreac (reaction cross section)',
     &' drift (mass and charge drift)',
     &' wida (mass and charge spread)',
     &' eav (mean excitation energy of recoil and ejectile)',
     &' eavin (inelastic excitation energy of recoil and ejectile)',
     &' tav (mean total kinetic energy)',
     &' wide (fluctuation in excitation energy)',
     &' avma (mean angular momentum of ejectile and recoil)',
     &' avm (mean total angular momentum loss).',
     &' pinelas (inelastic "probabilities").',
     &' changer (mean change in distance of closest approach)',
     &' Enter the name of one of the variables. Enter: no if no more'
        read(i5,'(a)')outl
        if(outl.ne.'thet'
     &.and.outl.ne.'pelas'.and.outl.ne.'ww'
     &.and.outl.ne.'ppcapt'.and.outl.ne.'dsigcap'
     &.and.outl.ne.'sigdrut'.and.outl.ne.'sigreac'.and.outl.ne.'drift'
     &.and.outl.ne.'wida'.and.outl.ne.'eav'.and.outl.ne.'tav'
     &.and.outl.ne.'wide'.and.outl.ne.'thetav'.and.outl.ne.'eavin'
     &.and.outl.ne.'avma'.and.outl.ne.'avm'.and.outl.ne.'changer'
     &.and.outl.ne.'pinelas'.and.outl.ne.'yrast'.and.outl.ne.'no')then
         write(*,*)' Wrong spelling, try again!'
         goto 107
        endif
        outleg1(leg1out)=outl
        if(intp.eq.'inter')write(10,'(a,a)')outleg1(leg1out),addtext
        if(outleg1(leg1out).eq.'no'.or.leg1out.eq.10)then
         leg1out = leg1out-1
         goto 49
        endif
        leg1out = leg1out+1
        goto 107
       else
        write(*,*)' Try again!'
        goto 45
       endif      
49     if(leg1out.eq.0)goto 8
46     if(intp.eq.'inter')write(*,'(a/a)')
     &' If you want tables of the variables which you specified', 
     &' Enter: tab.  If not enter: no'
       read(i5,'(a)')out
       addtext = '     !tab = table, no = no table'
       if(intp.eq.'inter')write(10,'(a,a)')out,addtext
       if(out.eq.'tab')then
        ntab = 1
       elseif(out.eq.'no ')then
        ntab = 0
       else
        write(*,*)' Try again!'
        goto 46
       endif
c
c  Preparation for total cross section
8     continue
       cross = 0.d0
       crossres = 0.d0
       sigcaptn = totmerg                  !   corr. nanni      
       if(ntab.eq.1)then
        do n=1,leg1out
         nn=20+n    
         name = outleg1(n)
         if(version.eq.'vms  '.or.version.eq.'unix* ')then
          filename1 = 'graz_'//head//'_'//name//'.dat'
          call condens_s(filename1,filename,le)
          open(unit=nn,file=filename,status='unknown')
         else
          write(filename,'(a4,i2,a4)') 'FOR0',nn,'.DAT'
         endif
         write(*,'(4a)')
     &' A table of ',name,' will be stored in ',filename
        enddo
       endif
c
c   Calculates probability for kinetic energy > telim
c   p0 (in common) is elastic probability neglecting capture
41    if(tres.lt.accur)then
       pelas = p0
      elseif(tres.ge.accur)then
       write(*,'(a,f7.2,/,a,f6.2)')
     &' Calculates probability for kinetic energy loss >',telim,
     &'     with energy resolution (FWHM) =', tres 
       timeedist1 = X05BAF()
       deltae = tres/2.35482d0
       ierr = 0
       call treslim_s(deltae,telim,accur,prelim,ierr)
       if(ierr.eq.1)goto 2
       timeedist2 = X05BAF()
       timeedist = timeedist2 - timeedist1
       pelas = 1.d0 - prelim
       wweff = log(pelas*ppscat)*(homegan+homegap)/4.d0/sqrt(pi)
c  Could perhaps be done better!
      endif
      widea0 = sqrt(wid2ea0)
      wideaa0 = sqrt(wid2eaa0)
c
c    I added the if statement as in grazing_9.for
c
      if(abs(1.d0-p0).gt.1.d-5) then
         drifta0 = drifta/(1.d0-p0)
         driftz0 = driftz/(1.d0-p0)
         tavs0 = tavs/(1.d0-p0)
        else
         drifta0 = 0.
         driftz0 = 0.
         tavs0 = 0.
      endif
c
c  Differential cross sections. The following quantities are in common:
c   cs (and cslab) are geometric cross sections i.e. differential
c    cross sections (elastic + reaction) with elastic scattering angle
c    but corrected for barrier penetration.
c   cscoul is simple Rutherford cross section
c   sig0drut is elastic cross section divided by Rutherford 
c   theta is elastic scattering angle in radian 
c   thet is elastic scattering angle in degrees
c   thetlabdeg is elastic scattering angle in laboratory in degrees
c   thetcou is elastic scattering angle in degrees
c   dthet = dtheta/drho
c   totmerg is cross section inside lmerg (no scattering)
      preac = 1.d0 - pelas
c   pelas and preac are corrected for finite resolution and cutoff
c    but neglecting capture
      sigelas = cs*pelas
      sigdrut = sigelas/cscoul
      sigreac = cs*preac
      if(lab.ne.0)then
       sigreaclab = cslab*preac 
      endif
c  Now correction for energy and angular momentum loss
      dell = -avm0
      delecm = -tavs0
      w00 = w0
      dellw = delw
      call scatangleavx_s(al,ph,dph,rtex,delecm,dell,drifta0,driftz0
     &,w00,dellw,ierr)
      if(ierr.eq.1)then
       print*,' Could not find scattering angle for l = ',al
       goto 77
      endif
      thetav = pi-ph
      dthetav = dph
      sinth = sin(theta)
      if(sinth.gt.accur)then
       dsigreac = sigreac*sinth*dthet/sin(thetav)/dthetav
      else
       dsigreac = sigreac*dthet/dthetav
      endif
      thetavdg = 180.d0*thetav/pi
c  In laboratory
      taua = sqrt(ecm/(ecm+delecm)*aa/aaa*(aa+drifta0)/(aaa-drifta0))
      thetrlab = atan(sin(thetav)/(cos(thetav)+taua))
      if(thetrlab.lt.0.)thetrlab=thetrlab + pi
      thetrlabdg = 180.d0/pi*thetrlab
      if(abs(sin(thetav)).lt.accur)then
       dsigreaclab = dsigreac*(1.d0-taua)**2
      else
       dsigreaclab = dsigreac*(sin(thetav)/sin(thetrlab))**2
     &/abs(cos(thetav-thetrlab))
      endif
c
77    dsigcapt = ppcapt*2.d0*pi/alorho**2*al*10.d0
      yrastfac = 34.7d0/(aa+aaa)**(5.d0/3.d0)
      yrast = yrastfac*al**2
c
c  Printout of standard output for differential cross sections
      if(leg1out.eq.0)then
       write(*,'(a,f7.2,a,f6.2,a)')
     &'    Calculates for l = ',al,' hbar (r = ',r,' fm)' 
      else
       write(98,'(/a,f7.2,a,f6.2,a)')
     &'    Results for l = ',al,' hbar (r = ',r,' fm)' 
       write(*,'(/a,f7.2,a,f6.2,a)')
     &'    Results for l = ',al,' hbar (r = ',r,' fm)' 
      endif
      nmaxx = max(1,leg1out)
      do n = 1,nmaxx
       if(outleg1(n).eq.'ppcapt')then
        write(98,'(a,d12.5)')' Capture probability = ',ppcapt
        write(*,'(a,d12.5)')' Capture probability = ',ppcapt 
       endif
       if(outleg1(n).eq.'dsigcap')then
        write(98,'(a,d12.5)')
     &' Capture cross section dsig/dl(mb) = ',dsigcapt
        write(*,'(a,d12.5)')
     &' Capture cross section dsig/dl(mb) = ',dsigcapt
       endif
       if(iout.eq.1..or.outleg1(n).eq.'yrast')then
        write(98,'(a,f8.3)')
     &' "Yrast-line" for compound nucleus (MeV) = ',yrast
        write(*,'(a,f8.3)')
     &' "Yrast-line" for compound nucleus (MeV) = ',yrast
       endif
       if(outleg1(n).eq.'thet')then
        if(lab.ne.1)then
         write(98,'(a,f7.2,a)')
     &' Elastic scattering angle = ',thet,' degrees'
         write(*,'(a,f7.2,a)')
     &' Elastic scattering angle = ',thet,' degrees'
        elseif(lab.ne.0)then
         write(98,'(a,f7.2,a)')
     &' Laboratory scattering angle = ',thetlabdeg,' degrees'   
         write(*,'(a,f7.2,a)')
     &' Laboratory scattering angle = ',thetlabdeg,' degrees'
        endif
       endif
       if(outleg1(n).eq.'pelas')then
        ppel = pelas*(1.d0-ppcapt)
        if(tres.gt.accur)then
         write(98,'(a,d12.5)')
     &' "Elastic probability" for emerging trajectory = ',ppel
         write(98,'(a,f8.2,a)')
     &'       with energy resolution =',tres,' MeV'
         write(98,'(a,f8.2,a)')
     &'       including kinetic energy loss less than ',telim,' MeV'
         write(*,'(a,d12.5)')
     &' "Elastic probability" for emerging trajectory = ',ppel
         write(*,'(a,f8.2,a)')
     &'       with energy resolution(FWHM) =',tres,' MeV'
         write(*,'(a,f8.2,a)')
     &'       including kinetic energy loss less than ',telim,' MeV'
        endif
        write(98,'(a,d12.5)')
     &' Elastic probability =',p0
        write(*,'(a,d12.5)')
     &' Elastic probability =',p0
       endif
       if(outleg1(n).eq.'sigdrut')then
        if(tres.gt.accur)then
         write(98,'(a,d12.5)')
     &' "Elastic" cross section/Rutherford = ',sigdrut
         write(98,'(a,f8.3,a)')
     &'       including kinetic energy loss less than = ', telim,' MeV'
         write(98,'(a,f8.3,a)')
     &'       with energy resolution(FWFM) =',tres,' MeV'
         write(*,'(a,d12.5)')
     &' "Elastic" cross section/Rutherford =',sigdrut
         write(*,'(a,f8.3,a)')
     &'       including kinetic energy loss less than = ', telim,' MeV'
         write(*,'(a,f8.3,a)')
     &'       with energy resolution(FWFM) = ',tres,' MeV'
        else
         write(98,'(a,d12.5)')
     &' Elastic cross section/Rutherford =',sigdrut
         write(*,'(a,d12.5)')
     &' Elastic cross section/Rutherford =',sigdrut
        endif
       endif
       if(outleg1(n).eq.'ww')then
        if(tres.gt.accur)then
         write(98,'(a,f6.2,a,d12.5,a)')' W(',r,') = ',wweff,'  MeV'
         write(*,'(a,f6.2,a,d12.5,a)')' W(',r,') = ',wweff,'  MeV'
        else 
         write(98,'(a,f6.2,a,d12.5,a)')' W(',r,') = ',ww,'  MeV'
         write(*,'(a,f6.2,a,d12.5,a)')' W(',r,') = ',ww,'  MeV'
        endif
       endif
c  Reaction cross sections,
       if(outleg1(n).eq.'preac')then
        ppra = preac*(1.d0-ppcapt)
        if(tres.gt.accur)then
         write(98,'(a,d12.5)')
     &' "Reaction probability" for emerging trajectory = ',ppra
         write(98,'(a,f8.2,a)')
     &'       with energy resolution(FWHM) = ',tres,' MeV'
         write(98,'(a,f8.2,a)')
     &'       excluding kinetic energy loss < ',telim,' MeV'
         write(*,'(a,d12.5)')
     &' "Reaction probability" for emerging trajectory =',ppra
          write(*,'(a,f8.2,a)')
     &'       with energy resolution(FWHM) = ',tres,' MeV'
          write(*,'(a,f8.2,a)')
     &'       excluding kinetic energy loss < ',telim,' MeV'
        else
         write(98,'(a,d12.5)')
     &' Reaction probability (without capture) = ',pr0
         write(*,'(a,d12.5)')
     &' Reaction probability (without capture) = ',pr0
        endif
       endif
       if(outleg1(n).eq.'thetav')then
        if(lab.ne.1)then
         write(98,'(a,f7.2,a)')
     &' Mean scattering angle for reactions = ',thetavdg,' degrees'
         write(*,'(a,f7.2,a)')
     &' Mean scattering angle for reactions = ',thetavdg,' degrees'
        elseif(lab.ne.0)then
         write(98,'(a,f7.2,a)')
     &' Lab scattering angle for reactions = ',thetrlabdg,' degrees'   
         write(*,'(a,f7.2,a)')
     &' Lab scattering angle for reactions = ',thetrlabdg,' degrees'
        endif
       endif
       if(outleg1(n).eq.'sigreac')then
        if(tres.gt.accur)then
         write(98,'(a,d13.5,a)')
     &' "Reaction" cross section = ',dsigreac,' mb/sterad'
         write(98,'(a,f8.2,a)')
     &'       with energy resolution(FWHM) = ',tres,' MeV'
         write(98,'(a,f8.2,a)')
     &'       excluding kinetic energy loss < ',telim,' MeV'
         write(*,'(a,d13.5,a)')
     &' "Reaction" cross section = ',dsigreac,' mb/sterad'
         write(*,'(a,f8.2,a)')
     &'       with energy resolution(FWHM) = ',tres,' MeV'
         write(*,'(a,f8.2,a)')
     &'       excluding kinetic energy loss < ',telim,' MeV'
        else
         write(98,'(a,d13.5,a)')
     &' Reaction cross section = ',dsigreac,' mb/sterad'
         write(*,'(a,d13.5,a)')
     &' Reaction cross section = ',dsigreac,' mb/sterad'
        endif
        if(lab.ne.0)then
         write(*,'(a,d12.5,a)')'     In lab = ',dsigreaclab,' mb/sterad'
         write(98,'(a,d12.5,a)')'    In lab = ',dsigreaclab,' mb/sterad'
        endif
       endif
c
c  Printout of drifts and spreads
       if(outleg1(n).eq.'drift')then
        write(98,'(a,f7.3,a,f7.3)')
     &' Mass drift=',drifta,'  Charge drift= ',driftz
        write(*,'(a,f7.3,a,f7.3)')
     &' Mass drift=',drifta,'  Charge drift= ',driftz
       endif
       if(outleg1(n).eq.'wida')then      
        write(98,'(a,f7.3,a,f7.3)')
     &' Mass width=',wida,' Charge width= ',widz
        write(*,'(a,f7.3,a,f7.3)')
     &' Mass width=',wida,'  Charge width= ',widz
       endif
       if(outleg1(n).eq.'eav')then
        write(98,'(a,f8.3,a)')
     &' Mean excitation energy of ejectile = ',eava,' MeV'
        write(98,'(a,f8.3,a)')
     &'         excluding elastic = ',eava0,' MeV'
        write(*,'(a,f8.3,a)')
     &' Mean excitation energy of ejectile = ',eava,' MeV'
        write(*,'(a,f8.3,a)')
     &'         excluding elastic = ',eava0,' MeV'
        write(98,'(a,f8.3,a)')
     &' Mean excitation energy of target = ',eavaa,' MeV'
        write(98,'(a,f8.3,a)')
     &'         excluding elastic = ',eavaa0,' MeV'
        write(*,'(a,f8.3,a)')
     &' Mean excitation energy of target = ',eavaa,' MeV'
        write(*,'(a,f8.3,a)')
     &'         excluding elastic = ',eavaa0,' MeV'
       endif
       if(outleg1(n).eq.'eavin')then
        write(98,'(a,f8.3,a)')
     &' Inelastic excitation energy of ejectile = ',eavina,' MeV'
        write(*,'(a,f8.3,a)')
     &' Inelastic excitation energy of ejectile = ',eavina,' MeV'
        write(98,'(a,f8.3,a)')
     &' Inelastic excitation energy of target = ',eavinaa,' MeV'
        write(*,'(a,f8.3,a)')
     &' Inelastic excitation energy of target = ',eavinaa,' MeV'
       endif
       if(outleg1(n).eq.'tav')then
        write(98,'(a,f8.3,a)')
     &' Mean kinetic energy loss = ',tavs,' MeV'
        write(*,'(a,f8.3,a)')
     &'         excluding elastic = ',tavs0,' MeV'
        write(*,'(a,f8.3,a)')
     &' Mean kinetic energy loss = ',tavs,' MeV'
        write(*,'(a,f8.3,a)')
     &'         excluding elastic = ',tavs0,' MeV'
       endif
       if(outleg1(n).eq.'wide')then
        write(98,'(a,f8.3,a)')
     &' Fluctuation in ejectile energy = ',widea,' MeV'
        write(98,'(a,f8.3,a)')
     &'         excluding elastic = ',widea0,' MeV'
        write(*,'(a,f8.3,a)')
     &' Fluctuation in ejectile energy = ',widea,' MeV'
        write(*,'(a,f8.3,a)')
     &'         excluding elastic = ',widea0,' MeV'
        write(98,'(a,f8.3,a)')
     &' Fluctuation in recoil energy = ',wideaa,' MeV'
        write(98,'(a,f8.3,a)')
     &'         excluding elastic = ',wideaa0,' MeV'
        write(*,'(a,f8.3,a)')
     &' Fluctuation in recoil energy = ',wideaa,' MeV'
        write(*,'(a,f8.3,a)')
     &'         excluding elastic = ',wideaa0,' MeV'
       endif
       if(outleg1(n).eq.'avma')then      
        write(98,'(a,f8.3,a)')
     &' Average angular momentum of ejectile = ',avma,' hbar'
        write(98,'(a,f8.3,a)')
     &'         excluding elastic = ',avma0,' hbar'
        write(*,'(a,f8.3,a)')
     &' Average angular momentum of ejectile = ',avma,' hbar'
        write(*,'(a,f8.3,a)')
     &'         excluding elastic = ',avma0,' hbar'
       endif
       if(outleg1(n).eq.'avma')then      
        write(98,'(a,f8.3,a)')
     &' Average angular momentum of recoil = ',avmaa,' hbar'
        write(98,'(a,f8.3,a)')
     &'         excluding elastic = ',avmaa0,' hbar'
        write(*,'(a,f8.3,a)')
     &' Average angular momentum of recoil = ',avmaa,' hbar'
        write(*,'(a,f8.3,a)')
     &'         excluding elastic = ',avmaa0,' hbar'
       endif
       if(outleg1(n).eq.'avm')then 
        write(98,'(a,f8.3,a)')
     &' Average total angular momentum loss = ',avm,' hbar'
        write(98,'(a,f8.3,a)')
     &'         excluding elastic = ',avm0,' hbar'
        write(*,'(a,f8.3,a)')
     &' Average total angular momentum loss = ',avm,' hbar'
        write(*,'(a,f8.3,a)')
     &'         excluding elastic = ',avm0,' hbar'
       endif
       if(outleg1(n).eq.'pinelas')then
        do nnn=1,2 
         if(nnn.eq.1)then
          text1 = 'lowest'
         elseif(nnn.eq.2)then
          text1 = 'giant'
         endif         
         do lam=2,3
          if(lam.eq.2)then
           text2 = text1//' 2+'
          elseif(lam.eq.3)then
           text2 = text1//' 3-'
          endif         
          do ii = 1,2
           if(ii.eq.1)then
            text = text2//' state in projectile'
           elseif(ii.eq.2)then
            text = text2//' state in target'
           endif         
           if(ff0in(ii,lam,nnn).gt.accur)then
            write(*,'(a,a,a,f9.6)')
     &' Excitation "probability" of ',text,' = ',ff0in(ii,lam,nnn)
            write(98,'(a,a,a,f9.6)')
     &' Excitation "probability" of ',text,' = ',ff0in(ii,lam,nnn)
           endif
          enddo
         enddo
        enddo 
       endif
       if(iout.eq.1.or.outleg1(n).eq.'changer')then      
        changer = rtex-r
        write(98,'(a,f7.4,a)')
     &' Average change in distance of closest approach = ',changer,' fm'
        write(*,'(a,f7.4,a)')
     &' Average change in distance of closest approach = ',changer,' fm'
       endif
      enddo
27    continue
c
c  Produce tables
      if(ntab.eq.0)goto 28
      do n = 1,leg1out
       nn=20+n
       if(outleg1(n).eq.'thet')then
        if(lab.eq.0)then
         if(istep.eq.1)write(nn,'(a)')
     &'    r        l      Elastic scattering angle in  degrees'
         write(nn,'(3f10.3)')r,al,thet
        elseif(lab.eq.1)then
         if(istep.eq.1)write(nn,'(a)')
     &'    r        l      Laboratory elastic scattering angle'
         write(nn,'(3f10.3)')r,al,thetlabdeg   
        elseif(lab.eq.2)then
         if(istep.eq.1)write(nn,'(a)')
     &'     r        l      Elastic scattering angle, Laboratory angle'
         write(nn,'(4f10.3)')r,al,thet,thetlabdeg   
        endif
       endif
       if(outleg1(n).eq.'pelas')then
        ppel = pelas*(1.d0-ppcapt)
        if(tres.gt.accur)then
         if(istep.eq.1)then
          write(nn,'(a)')
     &'     r  , theta  ,    l   ,    p0   , "Elastic" probability' 
          if(tres.gt.accur.and.istep.eq.1)then
           write(nn,'(30x,a,f8.2,a)')
     &' with energy resolution(FWHM) = ',tres,' MeV'
           write(nn,'(30x,a,f8.2,a)')
     &' including kinetic energy loss < ',telim,' MeV'
          endif
         endif
         write(nn,'(3f10.3,2f10.6)')r,thet,al,p0,ppel
        else
         if(istep.eq.1)write(nn,'(a)')
     &'     r    ,   theta   ,   l     ,    p0   ,Elastic probability'
         write(nn,'(3f10.3,2f10.6)')r,thet,al,p0,ppel 
        endif
       endif
       if(outleg1(n).eq.'sigdrut')then
        if(istep.eq.1)write(nn,'(a,/,a,f8.3,a/a,f8.3,a)')
     &'    r     ,  theta   ,"Elastic" cross section/Rutherford ',
     &'                 including kinetic energy loss < ', telim,' MeV',
     &'                      with energy resolution = ',tres,' MeV'
        write(nn,'(2f10.3,e12.3)')r,thet,sigdrut 
       endif
       if(outleg1(n).eq.'ww')then
        if(istep.eq.1)write(nn,'(a)')'     r       l       W(r)  MeV'
        if(tres.gt.accur)then
         write(nn,'(2f10.3,e12.3)')r,al,wweff
        else
         write(nn,'(2f10.3,e12.3)')r,al,ww
        endif
       endif
       if(outleg1(n).eq.'ppcapt')then
        if(istep.eq.1)write(nn,'(a)')
     &'     r        l      capture probability'
        write(nn,'(2f10.3,f10.6)')r,al,ppcapt
       endif
       if(outleg1(n).eq.'dsigcap')then
        if(istep.eq.1)then
         if(icapt.eq.1)then
          write(nn,'(a,f7.2/a)')
     &' Classical capture at l =',alcapt,
     &'     l        dsig_{cap}/dl(mb) dsig_{reac}/dl(mb)    r     '
         elseif(icapt.eq.2)then
          write(nn,'(a,f7.2,a/a)')
     &' Classical capture for l = 0 at ecm = ',ecmcapt,' MeV',
     &'     l        dsig_{cap}/dl(mb) dsig_{reac}/dl(mb)    r     '
         endif
         if(istep.eq.1.and.irstep.eq.0)then
          do j=1,11
           if(j.gt.1.and.dsigcaptn(j).le.1.d-20)goto 71
           if(j.eq.1.and.dsigcaptn(2).le.1.d-20)goto 71
           write(nn,'(f9.3,f18.6)')aln(j),dsigcaptn(j)
          enddo
         endif
        endif 
71      write(nn,'(f9.3,2f18.6,f13.6)')al,dsigcapt,sigreac0,r
       endif
       if(outleg1(n).eq.'yrast')then
        if(istep.eq.1)then
         write(nn,'(a,f7.2,a/a)')
     &'  Excitation energy of compound nucleus = ',excom,' MeV',
     &'    l    ,  "yrast" (MeV)   ,  "yrast" + 8(MeV)     r    '
         if(irstep.eq.1)then
          dlel = al/10.d0
          do j=0,9
           alnn = dble(j)*dlel
           yrastn = yrastfac*alnn**2
           yrastplus = yrastn + 8.d0          
           write(nn,'(f9.3,2f18.6)')alnn,yrastn,yrastplus
          enddo
         elseif(irstep.eq.0)then
          alnn = 0.d0
          yrastn = 0.d0
          yrastplus = yrastn + 8.d0          
          write(nn,'(f9.3,2f18.6)')alnn,yrastn,yrastplus
         endif
        endif 
        yrastplus = yrast + 8.d0          
        write(nn,'(f9.3,2f18.6,f11.6)')al,yrast,yrastplus,r
       endif
       if(outleg1(n).eq.'thetav')then
        if(lab.eq.0)then
         if(istep.eq.1)write(nn,'(/a)')
     &'       r         l       theta     thetav '
         write(nn,'(4f10.3)')r,al,thet,thetavdg
        elseif(lab.ne.0)then 
         if(istep.eq.1)write(nn,'(/a)')
     &'       r         l       theta     thetav   thetavlab'
         write(nn,'(5f10.3)')r,al,thet,thetavdg,thetrlabdg
        endif
       endif
c  Reaction cross sections,
       if(outleg1(n).eq.'sigreac')then
        if(lab.eq.0)then
         if(istep.eq.1)write(nn,'(/a,/,30x,a,f8.3,a,/,30x,a,f8.3,a)')
     &'    r    ,   l  ,  thetav , Reaction cross section in CM (mb/s)',
     &'    with energy resolution(FWHM) = ',tres,' MeV',
     &'    excluding kinetic energy loss < ',telim,' MeV.'
         write(nn,'(3f10.3,e12.3)')r,al,thetavdg,dsigreac
        elseif(lab.ne.0)then 
         if(istep.eq.1)write(nn,'(/a,/,30x,a,f8.3,a/,30x,a,f8.3,a)')
     &'    r   , theta ,thetalab ,React. cross sect. lab.(mb/s), in CM',
     &'    with energy resolution = ',tres,' MeV',
     &'    excluding kinetic energy loss < ',telim,' MeV.'
         write(nn,'(3f10.3,2e12.3)')r,thetavdg,thetrlabdg,dsigreaclab,
     &dsigreac
        endif
       endif
c
c   drifts and spreads
       if(outleg1(n).eq.'drift')then
        if(istep.eq.1)write(nn,'(a/a)')
     &'    Mass and charge drift excluding elastic',
     &'     r         l     Mass drift   Charge drift'
        write(nn,'(4f10.3)')r,al,drifta0,driftz0
       endif
       if(outleg1(n).eq.'wida')then      
        if(istep.eq.1)write(nn,'(a)')
     &'    r          l     Mass width, Charge width'
        write(nn,'(4f10.3)')r,al,wida,widz
       endif
       if(outleg1(n).eq.'eav')then
        if(istep.eq.1)write(nn,'(a/a)')
     &'          Mean excitation energy in MeV (excluding elastic)',
     &'       r        l     of ejectile   of recoil     total'
        write(nn,'(2f10.3,3f12.3)')r,al,eava0,eavaa0,eav0
       endif
       if(outleg1(n).eq.'eavin')then
        if(istep.eq.1)write(nn,'(a/a)')
     &'          Inelastic excitation energy in MeV ',
     &'       r        l     of ejectile   of recoil '
        write(nn,'(2f10.3,2f12.3)')r,al,eavina,eavinaa
       endif
       if(outleg1(n).eq.'tav')then
        tke = ecm-tavs0
        if(istep.eq.1)write(nn,'(a,/a)')
     &'     Mean kinetic energy loss in MeV and excluding elastic',
     &'       r        l      T_loss     (T_loss)ex.el.    TKE   '  
        write(nn,'(5f10.3)')r,al,tavs,tavs0,tke
       endif
       if(outleg1(n).eq.'wide')then
        if(istep.eq.1)write(nn,'(a/a)')
     &'        Fluctuation in excitation energy in MeV (elastic excl.)',
     &'    r           l     of ejectile  of recoil'
        write(nn,'(2f10.3,2f12.3)')r,al,widea0,wideaa0
       endif
       if(outleg1(n).eq.'avma')then      
        if(istep.eq.1)write(nn,'(a/a)')
     &' Average angular momentum of ejectile and recoil (0 excl.elast)',
     &'       r        l        avma      avma0     avmaa     avmaa0'
        write(nn,'(6f10.3)')r,al,avma,avma0,avmaa,avmaa0
       endif
       if(outleg1(n).eq.'avm')then      
        if(istep.eq.1)write(nn,'(a/a)')
     &'        Average total angular momentum loss and excl. elastic ',
     &'       r         l        avm       avm0'
        write(nn,'(4f10.3)')r,al,avm,avm0
       endif
       if(outleg1(n).eq.'changer')then      
        if(istep.eq.1)write(nn,'(a/a)')
     &'        Average change in distance of closest approach in fm',
     &'     r       l      changer'
        write(nn,'(3f10.3)')r,al,changer 
       endif
      enddo
c
c  integration of cross sections (not accurate for irstep.eq.0)
28    if(istep.eq.1.and.totmerg.lt.0.d0)then
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
      elseif(irstep.eq.1)then
       simp1 = delr/3.d0
       factor=2.d0*pi*10.d0*al*(al-alpr)/alorho**2
     &       +2.d0*pi*10.d0*r*(ecm-za*zaa*esq/2.d0/r
     &-unuc(r,ierr1)-r/2.d0*dunuc(r,ierr2))/ecm*simp1
      elseif(irstep.gt.1)then
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
c      print*,' factor = ', factor
      cross = cross + factor*(1.d0-p0)*ppscat
      crossres = crossres + factor*(1.d0-pelas)*ppscat
      sigcaptn = sigcaptn + factor*(1.d0-ppscat)
      if(irstep.ne.0)irstep = irstep+1
      istep = istep + 1
      if(intprep.eq.1.and.istep.ge.2)goto 33
c  Check for finished integration
      if(tres.gt.accur.and.thet.lt.scatmin.and.dthet.gt.0.)goto 33
      if(abs(1.-p0).lt.accur)goto 33
      alpr = al
      ierr = 4
      goto 2
c
c  printout of total cross sections 
33    write(*,'(/a,e12.4,a)')
     &' Total grazing cross section = ',cross,'  mbarn'
      write(98,'(/a,e12.4,a)')
     &' Total grazing cross section = ',cross,'  mbarn'
      if(tres.gt.accur)then
       write(*,'(/a,e12.4,a/a,f6.1,a)')
     &' "Total grazing cross section" = ',crossres,'  mbarn',
     &'       exluding energy loss < ', telim,' MeV'
       write(*,'(a,f6.2,a/a,f6.1,a)')
     &'       and with energy resolution(FWHM) = ',tres,' MeV',
     &'       for scattering angles larger than ',scatmin,' degrees'
       write(98,'(/a,e12.4,a/a,f6.1,a)')
     &' "Total grazing cross section" = ',crossres,'  mbarn',
     &'       exluding energy loss < ', telim,' MeV'
       write(98,'(a,f6.2,a/a,f6.1,a)')
     &'       and with energy resolution(FWHM) = ',tres,' MeV',
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
      write(12,'(2f10.4,1p3e12.4)') ecm,ecm-eb,cross,sigcaptn,sigtot
c
      if(tres.gt.accur)then
       write(*,'(/a,e12.4,a/a,f6.2,a)')
     &' "Total reaction cross section" = ',sigtotlim,'  mbarn',
     &'       exluding energy loss < ', telim,' MeV'
       write(*,'(a,f6.2,a/a,f6.1,a)')
     &' and with energy resolution = ',tres,' MeV',
     &'       for scattering angles larger than ',scatmin,' degrees'
       write(98,'(/a,e12.4,a/a,f6.2,a)')
     &' "Total reaction cross section" = ',sigtotlim,'  mbarn',
     &'       exluding energy loss < ', telim,' MeV'
       write(98,'(a,f6.2,a/a,f6.1,a)')
     &'       and with energy resolution = ',tres,' MeV',
     &'       for scattering angles larger than ',scatmin,' degrees'
      endif
      if(ntab.eq.1)then
       do n=1,leg1out
        nn=20+n    
        close(unit=nn)
       enddo
      endif
      ierr = 0
2     return
      end
