c      GRAZING calculates distributions in mass, energy and angular 
c       momentum according to a model which uses the Akyuz potential,
c       average formfactors and level densities and q-values as well 
c       as low-lying collective states. The calculations follow the
c       description (and notation) given in Nucl.Phys. A572(1994),191
c       and A594 (1995),203. A documentation with flow diagrams is
c       available together with the program in tex file.
c      This 9.th version has the date: 1997.05.23
c
c      In the VMS or UNIX* versions of the program
c       (version = vms or unix*) a record of the run is stored in 
c       the grazing_out.dat file.
c       When the program runs interactively the input data are 
c       stored in grazing_in.dat or grazing_in_*.dat for automatic rerun
c       In the UNIX version the input is allways stored in fort.10
c       The record of the run is then in fort.98
c      In the VMS or UNIX* versions the input and output files are (mostly)
c       given appropriate names. In the UNIX version all output is in 
c       numbered files.
c       One changes version setting the variable version = unix,unix* 
c       or vms  (first statement in the program).
c_______________________________________________________________
c
      program grazing_9
      implicit real*8(a-h,o-z)
      character*5 version
      common/vers/version
      character*40 defaultname
      character*40 programname
c  enter here name of default directory:
      data defaultname /'USR$DISK:[WINTHER.PROGRAM.GRAZING.GRAZING_9]'/
      data programname /'grazing_9'/
c  list,ienergy,irangen and irangep are defined in grazing.icl:
      include 'grazing.icl'
      common/za/za,aa,zaa,aaa
      common/iza/iza,iaa,izaa,iaaa
      common/ecm/ecm
      common/elab/elab
      common/iee/iealf1
      common/head/head
      common/leg/ecminit,rainbcdg,decm,cs,r,delr,al,theta,
     &cscoul,cslab,thet,thetlabdeg,thetcou,ppcapt,ww,pr0,
     &drifta,driftz,wida,widz,eava,eavaa,tavs,widea,wideaa,avmin,avm,
     &changer,changers,ppscat,totmerg,dthet,
     &deltae,sig0drut,sigreac0,sigreac0lab,cross,alorho,rho,
     &istep,iout,ipopn,norbit,irstep,lab
      common/tot/in00,ip00,i0nmin,i0nmax,i0zmin,i0zmax
      common/leg1/excom,dsigcaptn(11),aln(11),homegan,homegap
      common/leg2/sigcaptn
      common/leg3/delintloss,rcapt,alcapt
      common/probc/adia,capfac,ecmcapt,widcbar,sigcest,icapt
      common/capt/eb,rbar,i5,intp
      common/trans/rturn,rgraz,algraz,homegg,alatr
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      common/mml/ssl,qvalcom,threshn,threshp
      common/q/qns,qnp,qps,qpp
      common/cap/capn,capp,rr
      common/ba/bna,bnaa,bpa,bpaa
      common/ks/akbarn,akbarp
      common/capr/capred,pow
      common/ebe23/eex(2,3,2)
      common/be23/be23(2,3,2)
      common/giantwid/giantwid(3)
      common/zp23/zp(2,3,2)
      common/CD/Cin(2,3,2),Din(2,3,2)
      common/g_0/g0na,g0naa,g0pa,g0paa
      common/salt/amml,cckn,cckp,cden,cdep,r00
      common/pop/rwarn,dtann,dtanp,dtaann,dtaanp
      common/f0/f0ina,f0inaa,f0ns,f0np,f0ps,f0pp
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/eav/eav,eavns,eavnp,eavps,eavpp,eavtr,
     &eav0,eava0,eavaa0
      common/mav/avma,avmaa,avmans,avmanp,avmaps,avmapp,
     &avmaans,avmaanp,avmaaps,avmaapp,avmatr,avmaatr,
     &avma0,avmaa0,qavm0
      common/wid/wid2e,wid2ns,wid2np,wid2ps,wid2pp,wid2tr,wid2ea,
     &wid2eaa,wid2e0,wid2ea0,wid2eaa0
      common/inel/p0ina,p0inaa,eavina,eavinaa,wid2ina,wid2inaa,
     &avmina,avminaa,eava00,eavaa00,wid2a00,wid2aa00,inelast
      common/selec/ecutod,pcut,isel(4),lamsel(4),msmax
      common/f0in/ff0in(2,3,2)
      common/f0muin/ff0muin(2,3,2,4),f0inm(2,3,2)
      common/widin/wi2in(2,3,2)
      common/widcin/wi2cin(2,3,2)
      common/widcns/wid2cn,wid2cp
      common/probbar/delebar(481),barprob(481),ebarstep
      common/nzdis/ p0nn(irangen),p0zz(irangep),
     &eavdist(irangen,irangep),eavdist00,avmadist(irangen,irangep),
     &avmadist00,avmdist(irangen,irangep),avmdist00,
     &wid2dist(irangen,irangep),wid2dist00,wid2adist(irangen,irangep)
      common/c1/c1ns,c1np,c1ps,c1pp
      common/ej/ejns,ejnp,ejps,ejpp
      common/avm1/avmt1n,avmt1p,avmr1n,avmr1p
      common/scatang/delw,w0
      common/angdis/deltet
      common/bhom/bhomn,bhomp
      common/char/elemp,elemt
      dimension gn(11),form(2,3,2),ppin(2,3,2),dun(2),rdur(2,3)
      dimension eain(2,3,2),avin(2,3,2)
      dimension alt(2,3,2),hvt(2,3,2),et2(2,3,2,4)
      character*2 stand
      character*4 alab,oout
      character*5 elemp,elemt,elmc
      character*6 dum
      character*13 intpex
      character*5 intp
      character*3 out
      character*7 outgen(20)
      character*7 outgenn
      character*60 addtext
      character*11 autoadd1,autoadd
      character*16 head,head1
      character*11 headshort,headshort1
      character*7 addin
      character*15 outfilename
      character*25 infilename
      character*35 infilename1
      character*30 yfilename
      character*35 filename
      character*45 filename1
      character*26 runfile1,runfile
      character*4 eealf
      character*11 ddd
      character*5 ttt
c
c  Here starts program
      version='unix '
      if(version.eq.'vms  ')then
       open (unit=98,file='grazing_out.dat',status='new')
      elseif(version.eq.'unix*')then
       open (unit=98,file='grazing_out.dat',status='unknown')
      endif
      call date_s(ddd)
      call time_s(ttt)
      write(98,'(a,a11,a,a5)')
     &'  This output generated on : ', ddd,' at ',ttt
      write(*,'(a,a11,a,a5)')
     &'  This output generated on : ', ddd,' at ',ttt
      tinit = X05BAF()
1     iealf1=0
      addin = '       '
      intprep = 0
c  (for table headings)
      write(98,'(/5(/a))')
     &' This program calculates distributions in masses, charges,',
     &'  energies and angular momenta of fragments in collisions ',
     &'  between heavy nuclei according to a semiclassical model with',
     &'  average parameters for formfactors, level densities, etc.',
     &' The trajectory is calculated using the Akyuz potential'
      write(*,'(/5(/a))')
     &' This program calculates distributions in masses, charges,',
     &'  energies and angular momenta of fragments in collisions ',
     &'  between heavy nuclei according to a semiclassical model with',
     &'  average parameters for formfactors, level densities, etc.',
     &' The trajectory is calculated using the Akyuz potential'
      if(version.eq.'vms  ')then
       write(*,'(/10(/a))')
     $' You are running the VMS version of the program.',
     &' 1) If you want to run the program interactively, enter inter',
     &' 2) If you only want to prepare an input file, enter prep',
     &' 3) If you have made this allready, or have run the program',
     &'  interactively before, you have created an input file',
     &'  "grazing_in" and you may run automatically by entering: auto',
     &' 4) You may construct several different input files and run',
     &'  them automatically by adding _* to inter, prep or auto.',
     &'  * stands for any 6 character identification, i.e. inter_PB',
     &' 5) You may exit by entering: exit'
       read(5,'(a)',end=2)intpex
       if(intpex(1:5).eq.'inter')then
        intp = 'inter'
        i5=5
        if(intpex(1:6).eq.'inter_')addin = intpex(6:12)
        infilename1 = 'grazing_in'//addin//'.dat'
        call condens_s(infilename1,infilename,lefi)
        open(unit=10,file=infilename,status='unknown')
        autoadd1 = 'auto'//addin
        call condens_s(autoadd1,autoadd,leau)
        write(*,'(/a/a,a21,a/a/a,a15,/a/a/a)')
     &' You are now running interactively.',
     &' Your input is stored in ',infilename,' file',
     &'  so that you may run automatically next time',
     &'  by entering:',autoadd,
     &' The output is partly on the screen, but most of it',
     &'  is also stored in the grazing_out.dat file',
     &'  and in files specified later' 
       elseif(intpex(1:4).eq.'prep')then
        intprep = 1
        intp = 'inter'
        i5=5
        if(intpex(1:5).eq.'prep_')addin = intpex(5:11)
        infilename1 = 'grazing_in'//addin//'.dat'
        call condens_s(infilename1,infilename,lefi)
        open(unit=10,file=infilename,status='unknown')
        autoadd1 = 'auto'//addin
        call condens_s(autoadd1,autoadd,leau)
        write(*,'(/a/a,a21,a/a/a,a11/a))')
     &' You are now preparing in input file.',
     &' Your input is stored in ',infilename,' file',
     &'  so that you may run automatically next time',
     &'  by entering: ',autoadd,
     &' Most of the output is incomplete or wrong!'
       elseif(intpex(1:4).eq.'auto')then
        intp = 'auto'
        i5=10
        if(intpex(1:5).eq.'auto_')addin = intpex(5:11)
        infilename1 = 'grazing_in'//addin//'.dat'
        call condens_s(infilename1,infilename,lefi)
        open(unit=10,file=infilename,status='old')
       elseif(intpex(1:4).eq.'exit')then
        goto 2
       else
        write(*,*)' Try again'
        goto 1
       endif
      elseif(version.eq.'unix*')then
       write(*,'(/10(/a))')
     &' You are running the UNIX* version of the program.',
     &' 1) If you want to run the program interactively, enter inter',
     &' 2) If you only want to prepare an input file, enter prep',
     &' 3) If you have made this allready, or have run the program',
     &'  interactively before, you have created an input file',
     &'  "grazing_in.dat" and you may run with this as input',
     &' 4) You may construct several different input files ',
     &'  by adding _* to inter or prep, where * stands for  ',
     &'  any 6 character identification, i.e. inter_Ni_Pb',
     &' 5) You may exit by entering: exit'
       read(5,'(a)',end=2)intpex
       if(intpex(1:5).eq.'inter')then
        intp = 'inter'
        i5=5
        if(intpex(1:6).eq.'inter_')addin = intpex(6:12)
        infilename1 = 'grazing_in'//addin//'.dat'
        call condens_s(infilename1,infilename,lefi)
        open(unit=10,file=infilename,status='unknown')
        autoadd1 = 'auto'//addin
        call condens_s(autoadd1,autoadd,leau)
        write(10,'(a)')autoadd
        write(*,'(/a/a,a21,a/a/a/a/a)')
     &' You are now running interactively.',
     &' Your input is stored in ',infilename,' file',
     &'  so that you may run automatically next time',
     &' The output is partly on the screen, but most of it',
     &'  is also stored in the grazing_out.dat file',
     &'  and in files specified later'
       elseif(intpex(1:4).eq.'prep')then
        intprep = 1
        intp = 'inter'
        i5=5
        if(intpex(1:5).eq.'prep_')addin = intpex(5:11)
        infilename1 = 'grazing_in'//addin//'.dat'
        call condens_s(infilename1,infilename,lefi)
        open(unit=10,file=infilename,status='unknown')
        autoadd1 = 'auto'//addin
        call condens_s(autoadd1,autoadd,leau)
        write(10,'(a)')autoadd
        write(*,'(/a/a,a21,a/a/a))')
     &' You are now preparing in input file.',
     &' Your input is stored in ',infilename,' file',
     &'  so that you may run automatically next time',
     &' Most of the output is incomplete or wrong!'
       elseif(intpex(1:4).eq.'auto')then
        intp = 'auto'
        i5=5
        if(intpex(1:5).eq.'auto_')addin = intpex(5:11)
        infilename1 = 'grazing_in'//addin//'.dat'
        call condens_s(infilename1,infilename,lefi)
       elseif(intpex(1:4).eq.'exit')then
        goto 2
       else
        write(*,*)' Not an option. Try again!'
        goto 1
       endif
      elseif(version.eq.'unix ')then
       write(*,'(//a/a/a)')
     $' You are running the UNIX version of the program.'
     &, '  Output are stored in numbered files.'
     &, '  Write OK if this is OK'
       read(5,'(a)')intpex
       if(intpex.eq.'OK ')then
        intp = 'inter'
        i5 = 5
        intprep = 0
        write(10,'(a)')'auto'
        write(*,*)' Your input is stored in FOR10 file for rerun'
       elseif(intpex.eq.'auto')then
        intp = 'auto'
        i5 = 5
        intprep = 0
       else
        write(*,*)' wrong input, try again!'
        goto 1
       endif
      endif
c
c  Here re-enter for new run
17    tstart = X05BAF()
      ecmm = 0.d0
      decm = 0.d0
16    if(intp.eq.'inter')write(*,'(/a/a)')
     &' Enter projectile charge and mass numbers Z_a, A_a (integers)',
     &'  and target charge and mass numbers Z_A and A_A (integers):'
      read(i5,*,err=150) iza, iaa, izaa, iaaa
      icheck = 0
      goto 250
150   write(*,*)' your entry was not four numbers'
      icheck = icheck + 1
      if(icheck.gt.5)goto 2
      goto 16
250   addtext = '      ! Z_a,A_a,Z_A,A_A '     
      if(intp.eq.'inter')write(10,'(4i5,a24)') 
     &iza, iaa, izaa, iaaa,addtext
      za = dble(iza)
      aa = dble(iaa)
      zaa = dble(izaa)
      aaa = dble(iaaa)
      ierr1 = 0
      call element_s(iza,iaa,elemp,ierr1)
      ierr2 = 0
      call element_s(izaa,iaaa,elemt,ierr2)
      if(ierr1.eq.1.or.ierr2.eq.1)then
       if(ierr1.eq.1)write(*,*)' The projectile is unknown!'
       if(ierr2.eq.1)write(*,*)' The target is unknown!'
       goto 16
      endif
      write(*,*)' Projectile is ',elemp,', Target is ',elemt
      write(98,'(a,i3,a,i4/a,i3,a,i4)')
     &' Projectile charge and mass no. Za =',iza,' , Aa =',iaa,
     &' Target charge and mass no. ZA =',izaa,' , AA =',iaaa
      write(*,'(a,i3,a,i4/a,i3,a,i4)')
     &' Projectile charge and mass no. Za =',iza,' , Aa =',iaa,
     &' Target charge and mass no. ZA =',izaa,' , AA =',iaaa
c
c  Change parameters
19    if(intp.eq.'inter')write(*,'(5(/a))')
     &' If you want to run calculation with standard parameters',
     &' enter: s for standard.',
     &' If you want to change parameters enter: c for change.',
     &' If you only want to change the accuracy enter: a for accuracy.',
     &' If you have regrets enter: e for exit.'
      read(i5,'(a)')stand
      addtext = '      !  cannot be changed'
      if(intp.eq.'inter')write(10,'(a,a)') stand,addtext
      if(stand.eq.'s')then
       write(98,*)' Standard parameters overall'
       addtext='     !except if you change to a and enter accuracy here' 
       if(intp.eq.'inter')write(10,'(a,a)')addtext
       istand = 0
      else if(stand.eq.'c')then
       write(98,*)' Nonstandard parameters '
       istand = 1
      elseif(stand.eq.'a')then
       write(98,*)' Standard parameters, except accuracy'
       istand = 2
      elseif(stand.eq.'e')then
       goto 2
      else
       write(*,*)' Try again!'
       goto 19
      endif
c  These are the standard parameters:
      accur = .0001d0
      accur0 = accur
      iopt = 2
      ipopn = 0
      vdr = 0.d0
      capfac = 1.15d0
      cden = 8.d0
      cdep = 8.d0      
      adia = 1.d0
      cckn = 1.d0
      cckp = 1.d0
      amml = 1.6d0
      corr = 1.d0
      scale_delta=0.
      iqvalst = 0
      r00 = 1.27d0
      ecutod = 1.d0
      pcut = 5.d0
      delw = .5d0
      w0 = 1.d0
      deltet = 6.d0
      icoul = 0
      if(istand.eq.1)then
172   if(intp.eq.'inter')write(*,'(13(/a))')
     &' You have the possibility to change the following parameters:',
     &'  Prescription for the choice of effective Q-values: Enter Q',
     &'  The ion-ion potential: Enter I',
     &'  The nuclear radius parameter: Enter R',
     &'  Level density parameter a: Enter D',
     &'  Parameters for formfactors: Enter F',
     &'  Parameters for giant resonances: Enter G',
     &'  Parameters for barrier penetration and shift: Enter B',
     &'  Parameters for angular distributions: Enter L',
     &'  Overall accuracy of calculation: Enter A',
     &'  Parameter for special treatment of collective states: Enter P',
     &'  Options for Coulomb excitation: Enter C',
     &'  If exit: Enter E',
     &'  More possibilities appear later.',
     &' If no more changes here: Enter N'
       read(i5,'(a)')stand
       addtext = '      !  parameters to be changed, cannot be changed'
       if(intp.eq.'inter')write(10,'(a,a)') stand,addtext
       if(stand.ne.'Q'.and.stand.ne.'I'.and.stand.ne.'D'.and.stand.
     &ne.'F'.and.stand.ne.'B'.and.stand.ne.'N'.and.stand.ne.'A'
     &.and.stand.ne.'E'.and.stand.ne.'R'.and.stand.ne.'C'.
     &and.stand.ne.'G'.and.stand.ne.'P'.and.stand.ne.'L')then
        write(*,*)' Wrong spelling, try again!'
        goto 172
       endif
373    if(stand.eq.'N')goto 174
       if(stand.eq.'E')goto 17
       if(stand.eq.'C')then
        if(intp.eq.'inter')write(*,'(7(/a))')
     &' In the program Coulomb excitation is included for smaller',
     &'  impact parameters where transfer is important and neglected',
     &'  for larger distances to save time',
     &' You may choose other options: Enter' ,
     &'   1 to include Coulomb excitation for all impact parameters,',
     &'   2 to neglect Coulomb excitation altogether,',
     &'   3 to make no change (default)'
        read(i5,*,err=473) iopt
        icheck = 0
        goto 573
473     write(*,*)' your entry was not a number'   
        icheck = icheck + 1
        if(icheck.gt.5)goto 2
        goto 373
573     addtext = '      ! Coulomb excitation option'
        if(intp.eq.'inter')write(10,'(i5,a)')iopt,addtext
        if(iopt.gt.3.or.iopt.lt.1)goto 373
        if(iopt.eq.1)then
         icoul = 2
         write(*,*)
     &' Coulomb excitation included for all impact parameters'
         write(98,*)
     &' Coulomb excitation included for all impact parameters'
        elseif(iopt.eq.2)then
         icoul = 3
         write(*,*)' Coulomb excitation is neglected'
         write(98,*)' Coulomb excitation is neglected'
        elseif(iopt.eq.3)then
         icoul = 0
         write(*,*)
     &' Coulomb excitation neglected for large impact parameters'
         write(98,*)
     &' Coulomb excitation neglected for large impact parameters'
        endif
       endif
       if(stand.eq.'Q')then
        if(intp.eq.'inter')write(*,'(6(/a))')
     &' In the program the Q-values are taken as difference between',
     &'  Fermi levels in target and projectile (option 2 below)',
     &' You may choose other options for the Q-values: Enter' ,
     &'   1 for average over two steps,' ,
     &'   2 for difference between Fermi levels (default),',
     &'   3 for average between 1 and 2'
        read(i5,*,err=151) iopt
        icheck = 0
        goto 251
151     write(*,*)' your entry was not a number'   
        icheck = icheck + 1
        if(icheck.gt.5)goto 2
        goto 373
251     addtext = '      ! Q-value option'
        if(intp.eq.'inter')write(10,'(i5,a)')iopt,addtext
        if(iopt.gt.4.or.iopt.lt.1)goto 373
       elseif(stand.eq.'R')then
        if(intp.eq.'inter')write(*,'(4(/a))')
     &' The program uses a nuclear radius parameter of r_0 = 1.27 fm',
     &'  to calculate formfactors for transfer.',
     &'  You may change r_0 without changing the ion-ion potential.',
     &' Enter: r_0 (1.27 if no change)'
        read(i5,*,err=152) r00
        icheck = 0
        goto 252
152     write(*,*)' your entry was not a number'   
        icheck = icheck + 1
        if(icheck.gt.5)goto 2
        goto 373
252     addtext = '      ! nuclear radius parameter'
        if(intp.eq.'inter')write(10,'(f10.7,a)')r00,addtext 
       elseif(stand.eq.'I')then
        if(intp.eq.'inter')write(*,'(4(/a))')
     &' The program uses the ion-ion potential of Akyuz-Winther',
     &'  to calculate the trajectory and inelastic formfactors.',
     &'  You may change the radius by an amount vdr',
     &' Enter: vdr in fm. (0 if no change)'
        read(i5,*,err=160) vdr
        icheck = 0
        goto 260
160     write(*,*)' your entry was not a number'   
        icheck = icheck + 1
        if(icheck.gt.5)goto 2
        goto 373
260     addtext = '      ! change in ion-ion potential radius'
        if(intp.eq.'inter')write(10,'(f10.7,a)')vdr,addtext 
       elseif(stand.eq.'D')then
        if(intp.eq.'inter')write(*,'(/a/a/a)')
     &' Level densities are calculated with level density parameter a=8' 
     &,' You may change this by entering a for neutrons and protons',
     &'  i.e. Enter: a_n,a_p (8,8 if no change)'
        read(i5,*,err=153) cden,cdep
        icheck = 0
        goto 253
153     write(*,*)' your entry was not two numbers'   
        icheck = icheck + 1
        if(icheck.gt.5)goto 2
        goto 373
253     addtext = '      !  a_n, a_p level-density parameters'
        if(intp.eq.'inter')write(10,'(2f8.5,a)')cden,cdep,addtext
       elseif(stand.eq.'F')then      
        if(intp.eq.'inter')write(*,'(/a/a/a/a)')
     &' The formfactors for transfer are calculated as in the paper',
     &' You may change the overall normalization by entering',
     &'   a correction factor for the neutron formfactors k_n',
     &'   and the proton formfactors k_p (i.e. enter 1 1 for no change)'
        read(i5,*,err=154) cckn,cckp
        icheck = 0
        goto 254
154     write(*,*)' your entry was not two numbers'   
        icheck = icheck + 1
        if(icheck.gt.5)goto 2
        goto 373
254     addtext = '      ! transfer formfactor corrections to  k_n,k_p' 
        if(intp.eq.'inter')write(10,'(2f6.4,a)')cckn,cckp,addtext
        if(intp.eq.'inter')write(*,'(/a/a/a)')
     &' You may also enter a correction factor for the angular',
     &'   momentum dependence (M_l) in the formfactors',
     &'  (i.e. enter 1 for no change)'
        read(i5,*,err=155) camml
        icheck=0
        goto 255
155     write(*,*)' your entry was not a number'   
        icheck = icheck + 1
        if(icheck.gt.5)goto 2
        goto 373
255     amml = camml*1.6d0
        addtext = '      ! corrections to M_l in formfactor' 
        if(intp.eq.'inter')write(10,'(f6.4,a)')camml,addtext
       elseif(stand.eq.'P')then
        if(intp.eq.'inter')write(*,'(6(/a))')
     &' In the calculation of energy distributions collective',
     &'  states with energy less than ecutod*deltae(energy resolution)',
     &'  or with excitation "probability" larger than pcut are treated',
     &'  as continuos distributions. Standard values are ecutod = 1',
     &'  and pcut = 5. You may change this',
     &'  Enter: ecutod, pcut (i.e. 1, 5 if no change)'
        read(i5,*,err=190) ecutod,pcut
        icheck = 0
        goto 290
190     write(*,*)' your entry was not two numbers'   
        icheck = icheck + 1
        if(icheck.gt.5)goto 2
        goto 373
290     addtext = '      !ecutod and pcut for selection of coll. states'
        if(intp.eq.'inter')write(10,'(2f6.2,a)')ecutod,pcut,addtext
       elseif(stand.eq.'G')then
        iqvalst = 1      
       elseif(stand.eq.'L')then
        if(intp.eq.'inter')write(*,'(6(/a))')
     &' The scattering angle for a given reaction is calculated',
     &'   with a delay in energy loss over an interval delw',
     &'   around the value w0 of the parameter w, which is',
     &'   the time measured in units of the collision time.',     
     &'   Standard values are delw = .5 and w0 = 1.0',
     &'   Enter delw,w0 (i.e. .5, 1.0 if no change)'
        read(i5,*,err=191) delw,w0
        icheck = 0
        goto 291
191     write(*,*)' your entry was not two numbers'   
        icheck = icheck + 1
        if(icheck.gt.5)goto 2
        goto 373
291     addtext = '      !delw and w0 for scattering angle'
        if(intp.eq.'inter')write(10,'(2f6.2,a)')delw,w0,addtext
        if(intp.eq.'inter')write(*,'(4(/a))')
     &' In a Wilczinsky plot the angular distributions are smoothed',
     &'   over an angular range deltet (FWHM)',
     &'   Standard value is deltet = 6 degrees.',
     &'   Enter deltet (i.e. 6 if no change)'
        read(i5,*,err=192) deltet
        icheck = 0
        goto 292
192     write(*,*)' your entry was not a number'   
        icheck = icheck + 1
        if(icheck.gt.5)goto 2
        goto 373
292     addtext = '      !deltet for smoothing scattering angle'
        if(intp.eq.'inter')write(10,'(f6.2,a)')deltet,addtext
       elseif(stand.eq.'B')then
        if(intp.eq.'inter')write(*,'(6(/a))')
     &' Barrier penetration is calculated by WKB-method using the',
     &'   inverse of the reduced mass of the colliding nuclei.',
     &' You may introduce a correction factor (1+sc*delta) to account',
     &'   for less mass involved. Enter correction factor sc',
     &'   (i.e. 0 for no change )'
        read(i5,*,err=157) scale_delta
        icheck = 0
        goto 257
157     write(*,*)' your entry was not a number'   
        icheck = icheck + 1
        if(icheck.gt.5)goto 2
        goto 373
257     addtext = '      ! correction factor for barrier penetration'
        if(intp.eq.'inter')write(10,'(f7.3,a)')scale_delta,addtext
        capff=1.
        capfac = capff*1.15d0
       elseif(stand.eq.'A')then
173     if(intp.eq.'inter')write(*,'(/a/a)')
     &' Standard accuracy is .0001. You may change this by',
     &'   entering new overall accuracy < .01'
        read(i5,*,err=161) accur
        accur0 = accur
        icheck = 0
        goto 261
161     write(*,*)' your entry was not a number'   
        icheck = icheck + 1
        if(icheck.gt.5)goto 2
        goto 173
261     addtext = '      ! overall accuracy'
        if(intp.eq.'inter')write(10,'(f15.13,a)')accur,addtext
        if(accur.gt..01d0.or.accur.lt.0.)then
         write(*,*)' Accuracy must be less than .01'
         goto 173
        endif
       endif
       goto 172
      elseif(istand.eq.2)then
374    if(intp.eq.'inter')write(*,'(/a/a)')
     &' Standard accuracy is .0001. You may change this by',
     &'   entering new overall accuracy < .01'
        read(i5,*,err=162) accur
        accur0 = accur
        icheck = 0
        goto 262
162     write(*,*)' your entry was not a number'   
        icheck = icheck + 1
        if(icheck.gt.5)goto 2
        goto 374
262     addtext = '      ! overall accuracy'
       if(intp.eq.'inter')write(10,'(f15.13,a)')accur,addtext
       if(accur.gt..01d0.or.accur.lt.0.)then
        write(*,*)' Accuracy must be less than .01'
        goto 374
       endif
      endif
c      
174   write(98,'(a,e7.1)')' Overall accuracy = ',accur
      write(*,'(a,e7.1)')' Overall accuracy = ',accur
      if(intprep.eq.1)then
       accur = .001
       accur0 = accur
      endif
c  
77    if(version.eq.'vms  '.or.version.eq.'unix* ')then
       outfilename = 'grazing_out.dat'
      else
       outfilename = 'FOR098.DAT     '
      endif
      if(intp.eq.'inter')write(*,'(2(/a)/a,a15,a,4(/a)))')
     &' Besides the standard output, you may want additional output',
     &'  or a record of various parameters and variables',
     &'  in the ',outfilename,' file. They are specified later',
     &'   Enter yes if you want this.',
     &' You may also choose various special options',
     &'   Enter opt if you want this.',
     &' If you want both: Enter both. If you want neither: Enter no'     
      read(i5,'(a)')oout
      addtext = '     !additional output, special options, don`t change'
      if(intp.eq.'inter')write(10,'(a3,a60)')oout,addtext
      if(oout.eq.'no')then
       iout = 0
       goto 76
      elseif(oout.eq.'both'.or.oout.eq.'opt')then
       if(oout.eq.'both')iout = 1
       if(oout.eq.'opt')iout = 0
       if(intp.eq.'inter')write(*,'(/a/a)')
     &' If you want to run transfer reactions as a diffusion process',
     &' enter an artificial adiabaticity parameter << 1.If no enter 1'
        read(i5,*,err=163) adia
        icheck = 0
        goto 263
163     write(*,*)' your entry was not a number'   
        icheck = icheck + 1
        if(icheck.gt.5)goto 2
        goto 77
263     addtext = '      ! adiabaticity parameter'
       if(intp.eq.'inter')write(10,'(f6.4,a)')adia,addtext 
       if(intp.eq.'inter')write(*,'(/a/a)')
     &' If you want to calculate population of single particle states',
     &'  after the collision, enter yes, otherwise no'
       read(i5,'(a)')out
       if(out.eq.'yes')then
        ipopn = 1
       elseif(out.eq.'no ')then
        ipopn = 0
       else
        write(*,*)' Try again!'
        goto 77
       endif
       addtext = '      ! option to calculate populations'
       if(intp.eq.'inter')write(10,'(a,a)')out,addtext
      elseif(oout.eq.'yes')then
       iout = 1
      else
       write(*,*)' Try again'
       goto 77
      endif
76    ecm = 5.d0/(1.d0/aa+1.d0/aaa)
c  to make zeta = 1
      call constant_s
c  i.e. calculation of standard parameters for potential at low
c   bombarding energies (also hbar, e^2, etc.)
18    if(intp.eq.'inter')write(*,'(9(/a))')
     &' You must include inelastic scattering due to collective modes.',
     &'  in order to obtain realistic estimates for close collisions',
     &'  and for total cross sections',
     &' The low-lying modes are crucial, giant normally unimportant',
     &' You have four options:',
     &'  Enter 1 if only low-lying modes',
     &'  Enter 2 if both low-lying and high-lying (giant) modes',
     &'  Enter 0 if only high-lying (giant) modes',
     &'  Enter -1 if no inelastic modes.'
      read(i5,*,err=164) inelast
      icheck = 0
      goto 264
164   write(*,*)' your entry was not a number'   
      icheck = icheck + 1
      if(icheck.gt.5)goto 2
      goto 18
264   addtext = '   !inelastic options'
      if(intp.eq.'inter')write(10,'(i4,a60)')inelast,addtext
      if(inelast.lt.-1.or.inelast.gt.2)then
       write(*,*)' Try again inelast'
       goto 18
      elseif(inelast.ge.1.and.intp.eq.'auto ')then
       inelas = inelast+2
c  (this will avoid a stop in a batch run if no data are available)
      else
       inelas = inelast
      endif 
30    ierr = iqvalst
      call qvaletc_s(iopt,iout,inelas,ierr)
c  (Here Q-values and collective modes are found)
      if(ierr.eq.1)goto 2
      if(ierr.eq.3)then
       inelast = inelast-2
       write(*,*)
     &'   so low-lying modes not included, and capture is unrealistic'
      endif
c
c      General information
31    write(98,'(/a)')' General information for low energy (zeta=1)'
      write(*,'(/a)')' General information for low energy (zeta=1)'
      poten = 0.d0
      eb = 0.d0
      ierr = 1
      call ebar_s(eb,rbar,ierr)
c  (here information on barrier)
      if(ierr.eq.1)goto 2
      if(ierr.eq.2)then
       write(98,'(a,/a,f6.2,a)')' No Coulomb barrier for this system',
     &'    Sum of nuclear radii =',rbar,' fm'
       write(*,'(a,/a,f6.2,a)')' No Coulomb barrier for this system',
     &'    Sum of nuclear radii =',rbar,' fm'
       ierr1 = 1
       poten = unuc(rbar,ierr1)+za*zaa*esq/rbar
       if(ierr1.eq.1)goto 2
       eb = poten
       write(98,'(a,f8.2,a)')'    Potential energy at this distance ='
     &,poten,' MeV' 
       write(*,'(a,f8.2,a)')'    Potential energy at this distance ='
     &,poten,' MeV' 
      else
       write(98,'(a,f8.2,a)')
     &' Height of nominal Coulomb barrier =',eb,' MeV' 
       write(98,'(a,f6.2,a)')
     &' Radius of nominal Coulomb barrier =',rbar,' fm'
       write(*,'(a,f8.2,a)')
     &' Height of nominal Coulomb barrier =',eb,' MeV' 
       write(*,'(a,f6.2,a)')
     &' Radius of nominal Coulomb barrier =',rbar,' fm'
       ierr = 1
       call ecrit_s(ecr,alc,rcrit,ierr)
       if(ierr.eq.1)goto 2
       write(98,'(a,f8.2,a)')
     &' Critical energy =',ecr,' MeV, where there is no nominal barrier' 
       write(98,'(a,f8.2,a)')
     &' Corresponding critical angular momentum =',alc,' hbar' 
       write(98,'(a,f6.2,a)')
     &'   and critical radius =',rcrit,' fm'
       write(*,'(a,f8.2,a)')
     &' Critical energy =',ecr,' MeV, where there is no nominal barrier' 
       write(*,'(a,f8.2,a)')
     &' Corresponding critical angular momentum =',alc,' hbar' 
       write(*,'(a,f6.2,a)')
     &'   and critical radius =',rcrit,' fm'
      endif
      rdens = r0 - .52
c  since R_dens = R_pot - .26
      write(*,'(a,f5.2,a)')
     &' Sum of nuclear half-density radii = ',rdens,' fm'
      write(98,'(a,f5.2,a)')
     &' Sum of nuclear half-density radii = ',rdens,' fm'
c
c  Enter energy
      nestep = 1
22    if(intp.eq.'inter')write(*,'(/a/a/a)')
     &' Enter center of mass energy in MeV.',
     &'  If you know the laboratory energy or lab energy per nucleon',
     &'  you may first enter 0. Exit by negative number.'
      read(i5,*,err=165) ecm
      icheck = 0
      goto 265
165   write(*,*)' your entry was not a number'   
      icheck = icheck + 1
      if(icheck.gt.5)goto 2
      goto 22
265   addtext = '      ! Center of mass energy'
      if(intp.eq.'inter')write(10,'(f10.5,a)')ecm,addtext
      write(98,'(/)')
      if(ecm.lt.0.d0)goto 2
      if(ecm.gt.accur)then
       elab = ecm*(1.d0+aa/aaa)
       epslab = elab/aa
c  For heading of tables
       iealf1 = 1 
       goto 222
      endif
246   if(intp.eq.'inter')write(*,'(/a/a/a)')
     &' Enter laboratory bombarding energy in MeV.',
     &'  If you know laboratory energy per nucleon enter 0.',
     &'  If you want to repeat the energy input enter -1'
      read(i5,*,err=166) elab
      icheck = 0
      goto 266
166   write(*,*)' your entry was not a number'   
      icheck = icheck + 1
      if(icheck.gt.5)goto 2
      goto 246
266   addtext = '      ! Laboratory energy'
      if(intp.eq.'inter')write(10,'(f10.5,a)')elab,addtext
      if(elab.gt.accur)then
       ecm = elab/(1.+ aa/aaa)           
       epslab = elab/aa
c  For heading of tables
       iealf1 = 2 
       goto 221
      elseif(elab.lt.-accur)then
       goto 22
      endif
247   if(intp.eq.'inter')write(*,'(/a/a)')
     &' Enter laboratory energy per nucleon in MeV.',
     &'  If you want to repeat the energy input enter -1'
      read(i5,*,err=167) epslab
      icheck = 0
      goto 267
167   write(*,*)' your entry was not a number'   
      icheck = icheck + 1
      if(icheck.gt.5)goto 2
      goto 247
267   addtext = '      ! laboratory energy per nucleon'
      if(intp.eq.'inter')write(10,'(f10.5,a)')epslab,addtext
      if(epslab.gt.0)then
       elab = epslab*aa 
       ecm = elab/(1.+ aa/aaa)
c  For heading of tables
       iealf1 = 3 
      else
       goto 22
      endif
      write(98,'(a,f7.1,a)')' E_lab/A =' ,epslab,'  MeV/A'
      write(*,'(a,f7.1,a)')' E_lab/A =' ,epslab,'  MeV/A'      
221   write(98,'(a,f7.1,a)')' E_lab =' ,elab,'  MeV'
      write(*,'(a,f7.1,a)')' E_lab =' ,elab,'  MeV'
222   write(98,'(a,f7.1,a)') ' E_cm =',  ecm, ' Mev'
      write(*,'(a,f7.1,a)') ' E_cm =',  ecm, ' Mev'
      ecminit = ecm
      if(epslab.gt.10.)then
       write(*,*)
     &' Results may be unreliable: Bombarding energy too high!'
       goto 296
      elseif(ecm.lt..7d0*eb.or.ecm.lt..7d0*poten)then
       write(*,*)
     &' Bombarding energy very low (pure Coulomb scattering)'
       goto 296
      else
       goto 297
      endif
296   if(intp.eq.'inter')
     &write(*,*)' Do you really want this? (yes or no)'
      read(i5,'(a)')out
      addtext = '      ! Do you really want this energy?'
      if(intp.eq.'inter')write(10,'(a,a)')out,addtext
      if(out.eq.'no ')then
       goto 22
      elseif(out.eq.'yes')then
       goto 297
      else
       write(*,*)' Try again!'
       goto 296
      endif
297   continue 
1000  format
     &(/,3x,a5,' + ',a5,' at ',f6.1,' Mev lab. Date: ',a11,' at ',a5)
      write(*,1000)elemp,elemt,elab,ddd,ttt
      write(98,1000)elemp,elemt,elab,ddd,ttt
      if(iopt.eq.1)then
       write(98,*)' Q-value by average over two steps'
       write(*,*)' Q-value by average over two steps'
      elseif(iopt.eq.2)then
       write(98,*)' Q-values as difference between Fermi levels'
       write(*,*)' Q-values as difference between Fermi levels'
      elseif(iopt.eq.3)then
       write(98,*)' Q-value average between options 1 and 2'
       write(*,*)' Q-value average between options 1 and 2'
      elseif(iopt.eq.4)then
       write(98,*)' Experimental Q-values'
       write(*,*)' Experimental Q-values'
      endif
c
c  Printout of changed parameters
      vdr0 = vdr
      if(vdr.ne.0.)then
       write(98,'(a,f5.2,a)')
     &' Sum of nuclear radii in potential changed by ',vdr,' fm'
       write(*,'(a,f5.2,a)')
     &' Sum of nuclear radii in potential changed by ',vdr,' fm'
      endif
      if(abs(r00-1.27d0).gt.accur)then
       write(98,'(a,f5.2,a)')
     &' Nuclear radius parameter has been changed to ',r00,' fm'
       write(*,'(a,f5.2,a)')
     &' Nuclear radius parameter has been changed to ',r00,' fm'
      endif
      if(adia.ne.1.)then
       write(98,'(a,f8.6)')
     &' Program running as diffusion with artificial adia = ',adia
       write(*,'(a,f8.6)')
     &' Program running as diffusion with artificial adia = ',adia
      endif
      if(cden.ne.8.or.cdep.ne.8.)then
       write(98,'(a,/,2(f5.2,a))')
     &' The level density parameter a has been changed to',cden,
     &' for neutron and ',cdep,' for protons'
       write(*,'(a,/,2(f5.2,a))')
     &' The level density parameter a has been changed to',cden,
     &' for neutron and ',cdep,' for protons'
      endif
      if(abs(cckn-1.d0).gt.accur.or.abs(cckp-1.d0).gt.accur)then
       write(98,'(a,/,a,2(f5.3,a))')
     &' Entered correction factor for overall strength',
     &'   cckn = ',cckn,' for neutrons  and cckp = ',cckp,' for protons'
       write(*,'(a,/,a,2(f5.3,a))')
     &' Entered correction factor for overall strength',
     &'   cckn = ',cckn,' for neutrons  and cckp = ',cckp,' for protons'
      endif
      if(amml.ne.1.6d0)then
       write(98,'(a/a,f4.2)')
     &' The factor M_l in formfactor has been changed to',
     &'    M_l = ',amml 
       write(*,'(a/a,f4.2)')
     &' The factor M_l in formfactor has been changed to',
     &'    M_l = ',amml 
      endif
c
      write(98,*)' General r-independent information'
      write(*,*)' General r-independent information'
c  Here re-enter with new energy
80    continue
      zeta = exp(-.4d0/50.d0*(epslab-5.d0))
      zeta=1.
      write(98,'(a,f5.3)')
     &' Mean field reduction factor zeta = ',zeta 
      write(*,'(a,f5.3)')
     &' Mean field reduction factor zeta = ',zeta
      call constant_s
      ierr = 1
      call ebar_s(eb,rbar,ierr)
      if(ierr.eq.1) then
         goto 2
      endif
      ierr = 1
      call ecrit_s(ecr,alc,rcrit,ierr)
      if(ierr.eq.1)goto 2              
      if(iealf1.eq.1)iealf = int(ecm+.000001)
      if(iealf1.eq.2)iealf = int(elab+.000001)
      if(iealf1.eq.3)iealf = int(10.d0*epslab+.0000001)
      call convert_s(iealf,eealf,ierr)
      if(ierr.eq.1)goto 1
      head1 = elemp//'_'//elemt//'_'//eealf
      call condens_s(head1,head,len)
      headshort1 = elemp//'_'//elemt
      call condens_s(headshort1,headshort,len)
      alorho = sqrt(2.d0*aa*aaa/(aa+aaa)/hc2om*ecm)
c  Calculate l-grazing etc.
      ierr = 1
      call lgraz_s(ecr,alc,eb,algraz,rgraz,ierr)
      if(ierr.eq.0)then
       norbit = 1
       write(98,'(a,f7.2,a)')' Grazing l-value =',algraz,' hbar'  
       write(98,'(a,f7.2,a)')' Grazing radius =',rgraz,' fm'
       write(*,'(a,f7.2,a)')' Grazing l-value =',algraz,' hbar'  
       write(*,'(a,f6.2,a)')' Grazing radius =',rgraz,' fm'
       rinit = rgraz
      elseif(ierr.eq.1)then
       goto 2
      elseif(ierr.eq.2)then
       norbit = 2
       write(98,*)' Energy below barrier'
       write(*,*)' Energy below barrier'
       algraz = 0.d0
       ierr = 1
       al = 0.d0
       call turnp_s(al,rins,rinn,rturn,emax,emin,ierr)
       if(ierr.eq.1)goto 2
       write(98,'(a,f5.2,a)')
     &'   Turning point for central collisions = ',rturn,' fm'
       write(*,'(a,f5.2,a)')
     &'   Turning point for central collisions =',rturn,' fm'
       rinit = rturn
       rainbcdg = 180.d0
       goto 82
      elseif(ierr.eq.3)then
       norbit = 3
       write(98,*)' Energy above E_critical: No orbiting'
       write(*,*)' Energy above E_critical: No orbiting'
       rgraz = r0
       ierr1 = 1
       algraz = alofr(r0,ierr1)
       if(ierr1.eq.1)goto 2
       write(*,'(a,f5.2,a,f7.2)')
     &' Minimum safe distance = ',r0,' fm, minimum l = ',algraz
       write(98,'(a,f5.2,a,f7.2)')
     &' Minimum safe distance = ',r0,' fm, minimum l = ',algraz
       rinit = r0
      elseif(ierr.eq.4)then
       write(*,'(a,a)') '  ** Warning: Capture cross section not well'
     &                 ,' defined since no real barrier exist'
       potr0 = ecr
       if(ecm.gt.potr0)then
        rgraz = r0
        ierr1 = 1
        algraz = alofr(r0,ierr1)
        if(ierr1.eq.1)goto 2
        write(*,'(a,/a,f5.2,a,f7.2)')
     &' No barrier even for central collisions',
     &' Minimum safe distance = ',r0,' fm, minimum l = ',algraz
        write(98,'(a,/a,f5.2,a,f7.2)')
     &' No barrier even for central collisions',
     &' Minimum safe distance = ',r0,' fm, minimum l = ',algraz
        norbit = 4
        rinit = r0
       elseif(ecm.le.potr0)then
        ierr = 1
        al = 0.d0
        call turnp_s(al,rins,rinn,rturn,emax,emin,ierr)
        if(ierr.eq.1)goto 2
        write(*,'(a,/a,f5.2,a)')' No barrier',
     &' Turning point = ',rturn,' fm outside r0 even for l = 0'
        write(98,'(a,/a,f5.2,a)')' No barrier',
     &' Turning point = ',rturn,' fm outside r0 even for l = 0'
        algraz = 0.d0
        rinit = rturn
        norbit = 5
       endif
      endif
      timerain1 = X05BAF()
      write(*,*)' Calculating rainbow angles'
      ierr = 1
      accur = accur*10.d0
      call rainbow_s(rainbc,rainbn,alrc,alrn,ierr)
      accur = accur0
      timerain2 = X05BAF()
      timerainbow = timerain2-timerain1
      write(*,'(a,f7.2,a)')'  CPU time for rainbow =',timerainbow,
     &'sec.'
      if(ierr.eq.1)then
       write(*,*)' Could not find rainbow angle'
       write(98,*)' Could not find rainbow angle'
       rainbcdg = 180.d0
c (for scatmin later)
      elseif(ierr.eq.2)then
       write(*,*)' No rainbow - below barrier (should not happen!)'
       write(98,*)' No rainbow - below barrier(should not happen!)'
       rainbcdg = 180.d0
c (for scatmin later)
      elseif(ierr.eq.0.or.ierr.eq.3)then
       rainbcdg = rainbc*180.d0/pi
       write(98,'(a,f6.1,a)')' Coulomb rainbow angle =',rainbcdg,
     &' degrees'
        write(*,'(a,f6.1,a)')' Coulomb rainbow angle =',rainbcdg,
     &' degrees'
        write(98,'(a,f7.1,a)')' Rainbow angular momentum =',alrc,
     &' hbar'
        write(*,'(a,f7.1,a)')' Rainbow angular momentum =',alrc,
     &' hbar'
       ierr1 = 1
       call turnp_s(alrc,rins_r,rinn_r,r_rain,emax_r,emin_r,ierr1)
       if(ierr1.eq.1)goto 2
       write(98,'(a,f7.3,a)')' Rainbow distance =',r_rain,' fm'
       write(*,'(a,f7.3,a)')' Rainbow distance =',r_rain,' fm'
       if(ierr.eq.3)then
        rainbndg = rainbn*180.d0/pi
        write(98,'(a,f6.1,a)')' Nuclear rainbow angle =',rainbndg,
     &' degrees'
        write(*,'(a,f6.1,a)')' Nuclear rainbow angle =',rainbndg,
     &' degrees'
        ierr1 = 1
        call turnp_s(alrn,rins_r,rinn_r,r_rain,emax_r,emin_r,ierr1)
        if(ierr1.eq.1)goto 2
        write(98,'(a,f7.3,a)')' Nuclear rainbow distance =',r_rain,' fm'
        write(*,'(a,f7.3,a)')' Nuclear rainbow distance =',r_rain,' fm'
       endif
      endif
82    If(abs(qvalcom).gt.accur)then
       izc = iza+izaa
       iac = iaa + iaaa
       call element_s(izc,iac,elmc,ierr2)
       zcom = dble(izc)
       acom = dble(iac)
       excom = qvalcom + ecm
       thbarp = threshp +zcom*esq/(1.07d0*(acom**.3333d0+1.d0)+2.72d0)
       write(*,'(a,a5,a,f6.2,a/a,f6.2,a/a,f6.2,a)')
     &'  Excitation energy of fused system ',elmc,' is ',excom,' MeV',
     &'   Threshhold for neutron evaporation = ',threshn,' MeV',
     &'   Barrier for proton evaporation = ',thbarp,' MeV'
       write(98,'(a,a5,a,f6.2,a/a,f6.2,a/a,f6.2,a)')
     &'  Excitation energy of fused system ',elmc,' is ',excom,' MeV',
     &'   Threshhold for neutron evaporation = ',threshn,' MeV',
     &'   Barrier for proton evaporation = ',thbarp,' MeV'
      endif
c
c  Correct inelast, set ipt and calculate delta for barrier penetration
      if(inelast.ge.1)then
       if(be23(1,2,1).eq.0.d0.and.be23(1,3,1).eq.0.d0)then
        ipt = 2
c  (only target inelastic)
        if(be23(2,2,1).eq.0.d0.and.be23(2,3,1).eq.0.d0)then
         inelast = inelast-2
c  inelast > 0 was chosen, but no data!
         goto 669
        endif
        delta = 1.d0/4.d0/pi*aa*aaa/(aa+aaa)/hc2om*
     &raa**2*(5.d0/Din(2,2,1)+7.d0/Din(2,3,1))
c  (approximate capture radius (neglecting Coulexc.)rrcapt=rgraz+delta)
       elseif(be23(2,2,1).eq.0.d0.and.be23(2,3,1).eq.0.d0)then
c  (only projectile excitation)
        ipt = 1
        delta = 1.d0/4.d0/pi*aa*aaa/(aa+aaa)/hc2om*ra**2*
     &(5.d0/Din(1,2,1)+7.d0/Din(1,3,1))
       else
c  both target and projectile excitation
        ipt = 0
        delta = 1.d0/4.d0/pi*aa*aaa/(aa+aaa)/hc2om*(ra**2*
     &(5.d0/Din(1,2,1)+7.d0/Din(1,3,1))+raa**2*(5.d0/Din(2,2,1)
     &+7.d0/Din(2,3,1)))
       endif
      endif
669   if(inelast.eq.0)then
       write(*,*)
     &' capture unrealistic: low-lying collective states not included'
       ipt  = 0
       delta = 0.d0
      elseif(inelast.eq.-1)then
       delta = 0.d0
       if(norbit.eq.1)then
        rcapt = rgraz
        alcapt = algraz
        icapt = 1
       elseif(norbit.eq.2)then
        rcapt = rbar
        alcapt = 0.d0
        ecmcapt = eb
        icapt = 2
       elseif(norbit.eq.3)then
        rcapt = rgraz
        alcapt = algraz
        icapt = 1
       elseif(norbit.eq.4)then
        rcapt = r0
        alcapt = algraz
        icapt = 1
       elseif(norbit.eq.5)then
        rcapt = r0
        alcapt = 0.d0
        ecmcapt = potr0
        icapt = 2
       endif
       goto 69
      elseif(inelast.eq.2)then
       ipt = 0
      endif
      corr = 1.d0 + scale_delta*delta
c  i.e. calculates capture
      rgraz0 = rgraz
      rturn0 = rturn
      algraz0 = algraz
      iin = 0
      ierr = istand
      if(norbit.eq.2.or.norbit.eq.5)then
       algrazc = 0.d0
       rgrazc = rturn0
      else      
       algrazc = algraz0
       rgrazc = rgraz0
      endif
      is = 1
      if(accur.gt..0001d0)accur = .0001d0
      if(icoul.eq.3)ipt = ipt + 3
      call capture_s(ipt,iin,algrazc,rgrazc,rscat,rrcapt,alatr,sminn
     &,aalcapt,eecmcapt,et2,phid,rdur,thscat,alt,hvt
     &,is,nestep,iout,inelast,ierr)
      if(sminn.lt.0.) then
         sminn=0.
      endif
      accur = accur0
68    if(ierr.eq.1)then
       goto 2
      elseif(ierr.eq.3)then
       write(*,*)' No classical capture at this low energy'
       write(98,*)' No classical capture at this low energy'
       write(*,'(a,f7.2,a/a,f5.2,a)')
     &'  Lowest capture energy =',eecmcapt,' MeV at l=0',
     &'   with capture distance = ',rrcapt,' fm'
       write(98,'(a,f7.2,a/a,f5.2,a)')
     &'  Lowest capture energy =',eecmcapt,' MeV at l=0',
     &'   with capture distance = ',rrcapt,' fm'
       ecmcapt = eecmcapt 
       alcapt = 0.d0
       icapt = 2
      else
       write(*,'(/a,f6.1,/,a,f5.2,a)')
     &' Capture angular momentum = ',aalcapt,
     &'  capture distance = ',rrcapt,' fm'
       write(*,'(a,f6.1)')
     &'  Orbital angular momentum at capture distance = ',alatr
       write(98,'(/a,f6.1,/,a,f5.2,a)')
     &' Capture angular momentum = ',aalcapt,
     &'  capture distance = ',rrcapt,' fm'
       write(98,'(a,f6.1)')
     &' Orbital angular momentum at capture distance = ',alatr
       ecmcapt = ecm
       alcapt = aalcapt
       icapt = 1
       ierr1 = 1
       call turnp_s(alcapt,rimns,rimnn,runcapt,emmx,emmn,ierr1)
       if(ierr1.eq.1)goto 2
       if(alcapt.gt.algraz)then
        write(*,'(a,f7.4,a,/)')
     &'  Unperturbed turning point for l=l_capt ',runcapt,' fm'
        write(98,'(a,f7.4,a)')
     &'  Unperturbed turning point for l=l_capt ',runcapt,' fm'
       else
        write(*,*)' l_capt is smaller than l_graz'
        write(98,*)' l_capt is smaller than l_graz'
       endif
      endif
c   We now calculate the vdr that would make the mean barrier 
c   for l = alcapt be at an energy ecmcapt and use this to calculate 
c   the inverse barrier frequency, homegg and rcapt used in calculating
c   barrier probability distribution (probbar) and transmission (transm)
c   
69    ierr = 0
      vdr = vdr0
      call barrierpos_s(alcapt,ecmcapt,corr,vdrcap,rcapt,homegg,ierr)
      if(ierr.eq.1)then
       write(*,*)' Error in barrierpos_s - uno-'
       goto 2
      elseif(ierr.eq.0)then
       if(icapt.eq.1)write(*,'(a,f7.4,a)')
     &' Shift of radius to make capture at this l = ',vdrcap,' fm'
       if(icapt.eq.2)write(*,'(a,f7.4,a)')
     &' Shift of radius to make capture at this energy =',vdrcap,' fm'
       write(*,'(2(a,f7.4,a/))')
     &'  Radius for this barrier = ',rcapt,' fm',
     &'  Inverse barrier frequency = ',homegg,' MeV'
       if(icapt.eq.1)write(98,'(a,f7.4,a)')
     &' Shift of radius to make capture at this l = ',vdrcap,' fm'
       if(icapt.eq.2)write(98,'(a,f7.4,a)')
     &' Shift of radius to make capture at this energy =',vdrcap,' fm'
       write(98,'(2(a,f7.4,a/))')
     &'  Radius for this barrier = ',rcapt,' fm',
     &'  Inverse barrier frequency = ',homegg,' MeV'
      elseif(ierr.eq.2)then
       write(*,'(3(a,f7.4,a/))')
     &' No barrier i.e. homegg = ',homegg,'MeV', 
     &'  For l = l_capt capture is at r = ',rcapt,' fm',
     &'  with shift of radius vdrcapt = ',vdrcap,' fm'
       write(98,'(3(a,f7.4,a/))')
     &' No barrier i.e. homegg = ',homegg,'MeV', 
     &'  For l = l_capt capture is at r = ',rcapt,' fm',
     &'  with shift of radius vdrcapt = ',vdrcap,' fm'
      endif
c          
c  Calculate barrier distribution
      vdr = vdr0+vdrcap-accur
      ierr = 1
      smin = sminn
      accur=0.0001
      call probbar_s(inelast,iout,rcapt,alcapt,smin,rdur,phid,et2,ierr)
      accur = accur0
      vdr = vdr0
      if(ierr.eq.1)goto 2
c  Resulting barrierdistribution  is in common (delebar(n),barprob(n))
c   and may be printed out.
c  This finishes preparations that should be made in all runs
c  
c  Now for printout and jobtype:
      if(iout.ne.0.and.nestep.eq.1)then
771    if(intp.eq.'inter')write(*,'(4(/a))')     
     &' You may chose to print',
     &'  inelastic and transfer widths at capture distance',
     &'  and a table of the probability distribution of barriers',
     &' Answer: Yes or no'
       read(i5,'(a)')out
       if(out.eq.'yes')then
        ibartab = 1
       elseif(out.eq.'no ')then
        ibartab = 0
       else
        write(*,*)' Try again'
        goto 771
       endif
       addtext = '      !  printout of barrier distribution'
       if(intp.eq.'inter')write(10,'(a,a)')out,addtext
      endif    
      if(ibartab.eq.1)then
       if(version.eq.'vms  '.or.version.eq.'unix* ')then
        filename1 = 'graz_'//head//'_probbar.dat'
        call condens_s(filename1,filename,lefi)
        open(unit=11,file=filename,status='unknown')
       else
        filename = 'FOR011.DAT'
       endif
       acc10 = accur/10.d0
       do nn = 1,481
        barpro = barprob(nn)
        if(barpro.gt.acc10)
     &write(11,'(2f12.7)')delebar(nn),barprob(nn)
       enddo
       write(*,'(a,a35)')' Barrier distribution stored in ',filename
       widcin2a=sqrt(wi2cin(1,2,1))
       widcin2aa=sqrt(wi2cin(2,2,1))
       widcin3a=sqrt(wi2cin(1,3,1))
       widcin3aa=sqrt(wi2cin(2,3,1))
       write(98,'(a,2(/a,f6.3,a,f6.3))')
     &' Inelastic widths (MeV) at capture for low-lying states:',
     &' 2+ proj. widc_2a = ',widcin2a,' 2+ targ. widc_2A = ',widcin2aa,
     &' 3+ proj. widc_3a = ',widcin3a,' 3+ targ. widc_3A = ',widcin3aa
       widcn=sqrt(wid2cn)
       widcp=sqrt(wid2cp)
       write(98,'(a,/,2(a,f6.3))')
     &' Widths (MeV) at capture for transfer of',
     &'  neutrons widc_n = ',widcn,' and protons widc_p= ',widcp
      endif
12    if(nestep.ne.1)goto 49
71    if(iout.ne.0)then
       if(intp.eq.'inter')write(*,'(13(/a))')     
     &' You may chose to print intermediate variables',
     &'  from the following list (in paranthesis their significance):',
     &'  homega (hbar*omega_n, hbar*omega_p)',
     &'  phidot (hbar*phidot)',
     &'  alpha3 (phidot/omega_n, phidot/omega_p)',
     &'  transk (k_trans)',
     &'  c (k_trans/kappa_n, k_trans/kappa_p)',
     &'  bb (B_n, B_p), dd (D_n, D_p), ak (k_n, k_p)',
     &'  qopt (Q_opt^ns, Q_opt^ps), qs (Qs_ns, Qs_np, Qs_ps, Qs_pp)',
     &'  ej (J_ns,J_np,J_ps,J_pp), c0 (exp(L_n), exp(L_p))',
     &'  f0 (F0_ns,F0_np,F0_ps,F0_pp)',
     &'  dta (dtann,dtanp,dtaann,dtaanp)',
     &' Enter one of the names from the list. Enter: no if no more'
       read(i5,'(a)')outgenn
       if(outgenn.ne.'homega'.and.outgenn.ne.'phidot'.and.outgenn.ne.
     &'alpha3'.and.outgenn.ne.'transk'.and.outgenn.ne.'c'.and.outgenn.
     &ne.'bb'.and.outgenn.ne.'dd'.and.outgenn.ne.'ak'.and.outgenn.ne.
     &'qopt'.and.outgenn.ne.'qs'.and.outgenn.ne.'ej'.and.outgenn.ne.
     &'c0'.and.outgenn.ne.'f0'.and.outgenn.ne.'dta'.and.outgenn.ne.
     &'no')then
        write(*,*)' wrong spelling, try again!'
        goto 71
       endif
       outgen(iout)=outgenn
       addtext = '      !  variable for printout, no = no more'
       if(intp.eq.'inter')write(10,'(a,a)')outgen(iout),addtext
       if(outgen(iout).eq.'no')then
        iout = iout-1
        if(iout.eq.0)iout=1
        goto 21
       endif
       if(iout.ge.20)goto 21
       iout = iout+1
       goto 71
      endif
c      
c      Specify jobtype
21    if(intp.eq.'inter')write(*,'(14(/a))')
     &' You now specify which kind of job you want performed.',
     &' There are 8 types of jobs (legs):',
     &' 0. Only general information (independent of r)',
     &' 1. Sum over N and Z and integrated over E* and transferred L',
     &' 2. N and Z distribution  integrated over E* and transferred L ',
     &'       (cannot include evaporation: If so use 4)',
     &' 3. E* distributions, summed over N and Z and integrated over L',
     &' 4. E* distributions for each N and Z, integrated over L',
     &' 5. L distributions summed over N and Z and integrated over E*', 
     &' 6. L distributions for each N and Z and integrated over E*',
     &' 7. Distributions in E* and transferred L summed over N and Z',
     &' 8. Distributions in E* and transferred L for each N and Z',
     &' Of these 8 types the first 4 are implemented in various ways',
     &'  specified later. Enter which type you want:(exit on leg<0)'
      read(i5,*,err=168) leg
      icheck = 0
      goto 268
168   write(*,*)' your entry was not a number'   
      icheck = icheck + 1
      if(icheck.gt.5)goto 2
      goto 21
268   addtext = '      ! leg,  cannot be changed'
      if(intp.eq.'inter')write(10,'(i4,a)')leg,addtext
49    if(leg.gt.4)then
       write(*,*)' This option is not implemented'
       goto 21
      elseif(leg.eq.0)then
       if(nestep.gt.1)then
        goto 9
       elseif(nestep.eq.1)then
        if(version.eq.'vms  '.or.version.eq.'unix* ')then
         filename1 = 'graz_'//headshort//'_capture.dat'
         call condens_s(filename1,filename,len)
         open (unit=89,file=filename,status='unknown')      
         filename1 = 'graz_'//headshort//'_yrast_etc.dat'
         call condens_s(filename1,yfilename,len)
         open (unit=88,file=yfilename,status='unknown')
         write(*,*)' -------------------------------------------------'
         write(*,'(a/a,a28,a,a30,a)')
     &' You get a printout of quantities relevant for capture',
     &'   in the ',filename,' and ',yfilename,' files'
         write(98,*)' -------------------------------------------------'
         write(98,'(a/a,a28,a,a30,a)')
     &' You get a printout of quantities relevant for capture',
     &'   in the ',filename,' and ',yfilename,' files'
        else
         write(*,'(a/a)')
     &' You get a printout of quantities relevant for capture',
     &'   in the FOR088.DAT and FOR089.DAT files'
         write(98,'(a/a)')
     &' You get a printout of quantities relevant for capture',
     &'   in the FOR088.DAT and FOR089.DAT files'
        endif
        write(88,'(a,a5/a)')
     &' Excitation energy etc. in fused system ',elmc,
     &'      l     E_yrast   E-neut    E_prot '
        write(89,'(a,a5/a)')
     &' Excitation energy etc. in fused system ',elmc,
     &'  E_ex   l_capt  r_capt l_at_r  l_graz  r_graz widcbar sigcest'
        goto 89
       endif
      elseif(leg.lt.0)then
       goto 2
      endif
      descatfac = 1.d0
      if(nestep.gt.1)goto 48
c
c  Total cross section?
  139 continue
       delr = 2.d0*accur**.25d0
       dllr=delr
       if(istand.eq.1)then
        if(intp.eq.'inter')write(*,'(/a,f5.3,a/a)')
     &' To keep stated accuracy you should use delr = ',delr,' fm.',
     &'   Enter 0 if so. Otherwise enter delr in fm'
c
c       achtung read dllr=0.2 by reading 0. problems arize
c 
        read(i5,*,err=169) dllr
        icheck = 0
        goto 269
169     write(*,*)' your entry was not a number'   
        icheck = icheck + 1
        if(icheck.gt.5)goto 2
        goto 139
269     addtext = 
     &' ! stepwidth in r for total cross section (0 = default)'
        if(intp.eq.'inter')write(10,'(f6.4,a)')dllr,addtext
        if(dllr.gt.0.001d0)delr=dllr
       endif
       delrinit = delr
48    continue
      write(98,'(a,f5.3,a)')
     &' Total cross sections are calculated with stepwidth =', delrinit,
     &' fm'
      write(*,'(a,f5.3,a)')
     &' Total cross sections are calculated with stepwidth =', delrinit,
     &' fm'
c
c  Here starts preparations for total cross sections:
c
c  Determine steps in barrierregion (since potential and formfactors
c      are changed) first barrierheight for scattering at capture l 
c      or (below barrier) at l = 0
      ierr = 0
      call transm_s(alcapt,ecm,rcapt,alcapt,ecmcapt,ppcapbar,
     &descat,ierr)
      if(ierr.eq.1)then
       write(*,*)' Error in transm for l = alcapt'
       goto 2
      endif
      write(*,*)' Prepares total cross sections by calculating'
      if(icapt.eq.1)write(*,'(a,f12.7)')
     &,'   Capture probability at capture barrier = ',ppcapbar
      if(icapt.eq.2)write(*,'(a,f12.7)')
     &,'   Capture probability at l = 0 is ',ppcapbar
      vdr = vdr0
      ierr = 0                                               
      call elbar_s(alcapt,elba,rba,ierr)
      if(ierr.eq.1)goto 2
      write(*,'(a,f12.7)')'   Height of barrier  = ',elba
37    ebarscat = elba + descat
      write(*,'(a,f12.7)')
     &,'   Mean height of barrier for scattering = ',ebarscat
c
c  now vdr for mean barrier for scattering at capture l or l = 0
      vdr = vdr0 
      ierr = 0
      call barrierpos_s(alcapt,ebarscat,corr,vdral,rcal,homal,ierr)
      if(ierr.eq.1)then
       write(*,*)' Error in barrierpos_s -due-'
       goto 2
      elseif(ierr.eq.0)then
       write(*,'(a,f7.4,a/a,f7.4,a)')
     &'   Radius for this barrier = ',rcal,' fm',
     &'   Corresponding radiusshift vdral = ',vdral,' fm'
      elseif(ierr.eq.2)then
       write(*,'(a,f7.4,a/a,f7.4,a)')
     &'   No barrier. Capture at r = ',rcal,' fm',
     &'   Corresponding radiusshift vdral = ',vdral,' fm'
      endif
c  to ensure that there is no capture for scattering: 
      if(icapt.eq.1.and.vdral.gt.0.d0)then
       descatfac = descatfac*1.1d0
       descat = descat*descatfac
       write(*,*)' descatfac changed to ',descatfac
       goto 37
      endif
      algraz = alcapt
      vdr = vdr0+vdral
      ierr = 0
      call turnp_s(alcapt,ris,rin,rtur,ema,emi,ierr)
      if(ierr.eq.1)goto 2
      if(iout.ne.0.and.icapt.eq.1)write(*,*)
     &'  Turning point for scattering at barrier = ',rtur
      if(icapt.eq.2)write(*,*)
     &'  Turning point for scattering at l = 0 is ',rtur
c  then step in r
      ierr1 = 0
      rdis = rtur+delrinit
c  translated to step in l at this potential
      vdr = vdr0 + vdral
      ierr1 = 0
      al1 = alofr(rdis,ierr1)
      if(ierr1.eq.1)then
       write(*,*)' Error in delal'
       goto 2
      endif
      delal = (al1-alcapt)*1.d0 
      delal = int((al1-alcapt)*0.25d0) 
      delal = 2.d0
cnan      delal = 4.d0

      write(*,'(a,f8.3)')
     &' Suggested step in l: delal = ',delal
c  We now have the step delal in l
c
c  Prepare energy distribution steps for total cross sections 
c   by estimating distributions at capture distance or for l=0
      al = alcapt
      r = rtur
      vdr = vdr0+vdral
      ieprep = 1
      write(*,'(/a,f8.3,a/)')
     &' calculating for l = ',alcapt,' for preparation:'
c  leaves here knowing al, r and vdr for scattering at barrier or l=0:
      goto 57
c  and returns here to calculate ranges and energy-steps for total cs.
87    ierr = iout
      in00 = 0
      ip00 = 0
      call nzdistri_s(in00,ip00,inmax,inmin,izmax,izmin,ierr)
      if(ierr.eq.1)goto 2
c  Range in N and Z
      i0nmin=inmin
      i0nmax=inmax
      i0zmin=izmin
      i0zmax=izmax
      eavtr0 = eavtr/(1.d0-p0tr)
      wid2tr0 = wid2tr/(1.d0-p0tr)- eavtr0**2
      if(wid2tr0.lt.0.)then
       widtr0 = 0.d0
      else
       widtr0 = sqrt(wid2tr0)
      endif
c  Stepwidth in energy distributions
      delintloss = (eavtr0+widtr0)*accur**(.25d0)
      vdr = vdr0
      call constant_s
c  finished preparation of steps for tables
c
c   Excitation functions?
88    if(nestep.gt.1)goto 20
89    if(iealf1.eq.1)then
       if(intp.eq.'inter')write(*,'(/a/a/a)')
     &' If excitation functions wanted',
     &'  Enter step in center of mass energy dele>0 and Ecm_max,',
     &'  or dele<0 and Ecm_min.        If not enter 0 0'
       read(i5,*,err=170) decm,ecmm
       icheck = 0
       goto 270
170    write(*,*)' your entry was not two numbers'   
       icheck = icheck + 1
       if(icheck.gt.5)goto 2
       goto 89
270    addtext = '      ! Step in E_cm and maximum(or minimum)'
       if(intp.eq.'inter')write(10,'(2f9.4,a)')decm,ecmm,addtext
      elseif(iealf1.eq.2)then
       if(intp.eq.'inter')write(*,'(/a/a/a)')
     &' If excitation functions wanted',
     &'  Enter step in laboratory energy decm>0 and Elab_max,',
     &'  or dele<0 and Elab_min.        If not enter 0 0'
       read(i5,*,err=171) ddecm,eecmm
       icheck = 0
       goto 271
171    write(*,*)' your entry was not two numbers'   
       icheck = icheck + 1
       if(icheck.gt.5)goto 2
       goto 89
271    addtext = '      ! Step in E_lab and maximum(or minimum)'
       if(intp.eq.'inter')write(10,'(2f9.4,a)')ddecm,eecmm,addtext
       decm = ddecm/(1.d0+aa/aaa)
       ecmm = eecmm/(1.d0+aa/aaa)
      elseif(iealf1.eq.3)then
       if(intp.eq.'inter')write(*,'(/a/a/a)')
     &' If excitation functions wanted',
     &'  Enter step in laboratory energy/nucleon dele/A>0 and E/A_max,',
     &'  or dele/A<0 and E/A_min.        If not enter 0 0'
       read(i5,*,err=1172) ddecm,eecmm
       icheck = 0
       goto 272
1172   write(*,*)' your entry was not two numbers'   
       icheck = icheck + 1
       if(icheck.gt.5)goto 2
       goto 89
272    addtext = '      ! Step in E_lab/A and maximum(or minimum)'
       if(intp.eq.'inter')write(10,'(2f9.4,a)')ddecm,eecmm,addtext
       decm = aa*ddecm/(1.d0+aa/aaa)
       ecmm = aa*eecmm/(1.d0+aa/aaa)
      endif
      if(leg.eq.0)goto 9
      if(abs(decm).lt.accur)then
       addtext='    !If you change to Step.ne.0 enter maxtime(min) here' 
       if(intp.eq.'inter')write(10,'(a4,a)')'    ',addtext
       goto 20
      endif
      ecmmm = 10.d0*aa*aaa/(aa+aaa)
      if(decm*(ecmm-ecminit).lt.0.d0)then
       write(*,*)' Try again'
       goto 88
      endif
      neestim = int(1.d0+abs((ecmm-ecminit)/decm))
      if(neestim.gt.50)then
       write(*,*)' Too many steps!'
       goto 88
      endif
90    if(intp.eq.'inter')write(*,'(/a/a)')
     &' Warning: For some options this may take excessive CPU time',
     &' Indicate the maximum CPU time you allow (timemax) in minutes'
      read(i5,*,err=1173) timemax
      icheck = 0
      goto 273
1173  write(*,*)' your entry was not a number'   
      icheck = icheck + 1
      if(icheck.gt.5)goto 2
      goto 90
273   addtext = '      ! Allowed time in min. Delete line if Step is 0.'
      if(intp.eq.'inter')write(10,'(f6.2,a)')timemax,addtext
c  finished possibility for energy steps
c      
c   Determine region of full capture
20    istep = 0
      irstep = 0
      alpr = -1.d0
      al = 0.d0
201   ierr = 0
      if(accur.gt..0001d0)accur = .0001d0
      print *, ' prima di transm_s'
      call transm_s(al,ecm,rcapt,alcapt,ecmcapt,ppcapt,descat,ierr)
      print *, ' dopo di transm_s'
      descat = descat*descatfac
      accur = accur0
      if(ierr.eq.1)then
       write(*,*)' Error in transm for l = al'
       goto 2
      endif
      ppscat = 1.d0-ppcapt
      if(ppscat.lt.0.d0)ppscat = 0.d0
      pplim = 10.d0*accur
      if(pplim.gt..001d0)pplim = .001d0
      print *, ' pplim  ==',pplim
      print *, ' ppscat ==',ppscat
      if(ppscat.lt.pplim)then
c  ie. full capture (for l= alpr or less)
       alpr = al
       al = al+delal
       goto 201
      endif
c
      al = al+delal
      print *, ' al ==',al
c
      if(alpr.ge.0.d0)then
       rhomerg = almerg/alorho
       totmerg = 10.d0*pi*rhomerg**2
      elseif(alpr.lt.0.d0)then
       almerg = 0.d0
       totmerg = -1.d-20
      endif
      if(almerg.gt.0.d0)then
       write(*,'(a,f7.2,a)')'  Full capture for l < ',almerg
       write(98,'(a,f7.2,a)')'  Full capture for l < ',almerg
      endif
c   merging part in ten steps or less for table of dsigcaptn
      nndl = min(int(al/delal+.1),11)-1
cnan
      if(nndl.lt.1) nndl=1
cnan      
      print *,' al,delal,nndl ===',al,delal,nndl
      dlel = al/float(nndl)
      facl = 2.d0*pi/alorho**2*10.d0
      do j = 1,11
       dsigcaptn(j) = 0.d0
      enddo
      do j=0,nndl-1
       aln(j+1) = dble(j)*dlel
       alnn = aln(j+1)
       dsigcaptn(j+1) = facl*alnn
      enddo
      sigcaptn = totmerg
      irstep = 0
      istep =  1
      rwarn = 0.d0
c
       if(istep.eq.1.and.nestep.eq.1)then
c  Laboratory system?
231    if(intp.eq.'inter')write(*,'(/a/a/a)')
     &' Do you want results in the laboratory system?',
     &' Enter: lab (for laboratory), cm (for center of mass),',
     &'        both (for both systems)(not always implemented)'
       read(i5,'(a4)')alab
       if(alab.eq.'lab ')then
        lab=1
       elseif(alab.eq.'cm  ')then
        lab=0
       elseif(alab.eq.'both')then
        lab=2
       else
        write(*,*)' Wrong spelling'
        goto 231
       endif      
       addtext = '  ! lab=lab, cm=center of mass, both=both system'
       if(intp.eq.'inter')write(10,'(a,a)')alab,addtext
      endif
c   Here starts calculation for trajectory knowing only al (and 
c     barrier distribution).
51    ierr = 0
      if(accur.gt..0001d0)accur = .0001d0
      call transm_s(al,ecm,rcapt,alcapt,ecmcapt,ppcapt,descat,ierr)
      descat = descat*descatfac
      accur = accur0
      if(ierr.eq.1)then
       write(*,*)' Error in transm for l = al'
       goto 2
      endif 
      ppscat = 1.d0-ppcapt
      if(ppscat.lt.0.d0)ppscat = 0.d0
cnan      pplim = accur
      pplim = accur/1000.d0                  ! attenzione 
      if(ppscat.lt.pplim)then
       write(*,*)' Full capture for l = ',al
       write(98,*)' Full capture for l = ',al
      endif
      if(ppcapt.lt.pplim.and.irstep.eq.0)irstep = 1
      vdr = vdr0 
      ierr = 0
c  calculate unshifted barrier for vdr=vdr0
      call elbar_s(al,ebaral,rbaral,ierr)
      ebarscat = ebaral + descat
      write(*,'(/,3(a,f12.7/),a,f12.7)')'     For l = ',al
     &,'       Capture probability = ',ppcapt
     &,'       Mean barrier for scattering = ',ebarscat
     &,'       Mean shift in barrier for scattering = ',descat
      ierr = 0
      vdr = vdr0
      call barrierpos_s(al,ebarscat,corr,vdral,rcal,homal,ierr)
      if(ierr.eq.1)then
       write(*,*)' Warning: Error in barrierpos for l = ',al
       goto 2
      elseif(iout.eq.1.and.ierr.eq.0)then
       write(*,'(a,f7.4,a/,a,f7.4,a)')
     &'       Radius for scattering barrier = ',rcal,' fm',
     &'       Barrierposition shift vdral = ',vdral,' fm'
      elseif(iout.eq.1.and.ierr.eq.2)then
       write(*,'(a,f7.4,a/a,f7.4,a)')
     &'       No barrier. Capture at r = ',rcal,' fm',
     &'       Shift of radius vdral = ',vdral,' fm'
      endif
      vdr = vdr0+vdral
      ierr = 0
      call lgraz_s(ecr,alc,eb,algraz,rgraz,ierr)
      if(ierr.eq.1)goto 2
      call turnp_s(al,rins,rinn,rtur,emax,emin,ierr)
      if(ierr.eq.1)goto 2
      write(*,*)'      Actual turning point = ',rtur
      r = rtur
      rwarn = 0   
c
c   Here starts calculations for each trajectory knowing al,r,vdr
c    first recalculating algraz etc if vdr.ne.vdr0
57    continue
      if(irstep.le.2)then
       ierr = 0
       call lgraz_s(ecr,alc,eb,algraz,rgraz,ierr)
       if(ierr.eq.0)then
        norbit = 1
        write(*,'(a,f7.2,a)')' Grazing l-value =',algraz,' hbar'  
        write(*,'(a,f6.2,a)')' Grazing radius =',rgraz,' fm'
       elseif(ierr.eq.1)then
        write(*,*)' Error in lgraz in main'
        goto 2
       elseif(ierr.eq.2)then
        norbit = 2
       write(*,*)' Energy below barrier(for l=0)'
        algraz = 0.d0
        ierr = 1
        all = 0.d0
        call turnp_s(all,rins,rinn,rturn,emax,emin,ierr)
        if(ierr.eq.1)goto 2
c       write(98,'(a,f5.2,a)')
c     &'   Turning point for central collisions = ',rturn,' fm'
       write(*,'(a,f5.2,a)')
     &'   Turning point for central collisions =',rturn,' fm'
c       rinit = rturn
c       goto 82
       elseif(ierr.eq.3)then
        norbit = 3
       write(98,*)' Energy above E_critical: No orbiting'
c       write(*,*)' Energy above E_critical: No orbiting'
cnan
        rgraz = r0 +0.2
        ierr1 = 1
        algraz = alofr(rgraz,ierr1)
cnan
        if(ierr1.eq.1)goto 2
       write(*,'(a,f7.2,a,f5.2)')
     &' No barrier: Minimum l = ',algraz,' for r = ',r0
       elseif(ierr.eq.4)then
        potr0 = ecr
        if(ecm.gt.potr0)then
         rgraz = r0
         ierr1 = 1
         algraz = alofr(r0,ierr1)
         if(ierr1.eq.1)goto 2
        write(*,'(a,/a,f7.2,a,f5.2)')
     &'  No barrier even for central collisions',
     &'  Minimum l = ',algraz,' for r = ',r0
         norbit = 4
c         rinit = r0
        elseif(ecm.le.potr0)then
         ierr = 1
         all = 0.d0
         call turnp_s(all,rins,rinn,rturn,emax,emin,ierr)
         if(ierr.eq.1)goto 2
        write(*,'(a,/a,f5.2,a)')' No barrier',
     &' Turning point = ',rturn,' fm outside r0 even for l = 0'
         algraz = 0.d0
c         rinit = rturn
         norbit = 5
        endif
       endif
      endif
      ierr = 0
      rho = al/sqrt(aa*aaa/(aa+aaa)*ecm*2.d0/hc2om)
c
70    ppscat = 1.d0-ppcapt
      if(ppscat.lt.0.d0)ppscat = 0.d0
c  Here inelastic probabilities including nonlinear and Coulomb excitation:
      nma = 2
      nmi = 1
      if(inelast.eq.0)nmi = 2                                                   
      if(inelast.eq.1)nma = 1
      if(inelast.eq.-1)nma = 0
      iin = 1
      rscat = r
      ierr = istand
      if(icoul.eq.3)ipt = ipt + 3
      ist = 2
c ie. so that no questions asked about changing parameters
      call capture_s(ipt,iin,al,r,rscat,rcaptt,alatrr,smin,
     &alca,ecmca,et2,phid,rdur,thscat,alt,hvt
     &,ist,nestep,iout,inelast,ierr)
      if(ierr.eq.1)goto 2
      if(ierr.eq.3)then
c  increasing descat if capture
       descatfac = descatfac*1.1d0
       write(*,*)' Increasing descat because of capture for l =',al
       goto 51
      endif
      do n=nmi,nma
       do i=1,2
        do lam = 2,3
         f0inm(i,lam,n) = alt(i,lam,n)
         avin(i,lam,n) = alt(i,lam,n)
         eain(i,lam,n) = hvt(i,lam,n)
         wi2in(i,lam,n) = hvt(i,lam,n)*eex(i,lam,n)
         ff0in(i,lam,n) = hvt(i,lam,n)/eex(i,lam,n)
        enddo
       enddo
      enddo
c
      f0ina = 0.d0
      f0inaa = 0.d0
      eavina = 0.d0
      eavinaa = 0.d0
      wid2ina = 0.d0
      wid2inaa = 0.d0
      avmina = 0.d0
      avminaa = 0.d0
      do n = 1,2
       do lam = 2,3
        f0ina = f0ina + ff0in(1,lam,n)
        f0inaa = f0inaa + ff0in(2,lam,n)
        eavina = eavina + eain(1,lam,n)
        eavinaa = eavinaa + eain(2,lam,n)
        wid2ina = wid2ina + wi2in(1,lam,n)
        wid2inaa = wid2inaa + wi2in(2,lam,n)
        avmina = avmina + avin(1,lam,n)
        avminaa = avminaa + avin(2,lam,n)
       enddo
      enddo
      f0in = f0ina + f0inaa
      eavin = eavina + eavinaa
      wid2in = wid2ina + wid2inaa
      avmin = avmina + avminaa
c
c  Here starts calculation of transfer probabilities with smin 
c   from running capture_s in scattering mode above
c  kappa_n = cappan  kappa_p = cappap
      cappan = capn/2.d0
      cappap = capp/2.d0
c  homegan = hbar*omega_n   homegap = hbar*omega_p  for transfer          
      ierr1 = 1
      dur = dueff(al,r,ierr1)
      hom2 = - cappan*hc2om*(aa+aaa)/aa/aaa*dur   
      if(ierr1.eq.1)goto 2
      if(hom2.le.0.)then
       write(*,*)' Error: hom2 = ',hom2
       goto 2
      endif
      homegan = sqrt(hom2)*adia
      homegap = homegan*sqrt(cappap/cappan)
c  homegain = hbar*omega_in
      homegain = sqrt(hom2/cappan*ainv)
c  calculate phidot 
      ierr1 = 1
      phidot2 = 2.d0*hc2om*(aa+aaa)/aa/aaa*
     &         (ecm-unuc(r,ierr1)-za*zaa*esq/r)/r**2
      if(ierr1.eq.1)goto 2
      if(phidot2.lt.-accur)then
       write(*,*)' Error: phidot^2 = ',phidot2
       goto 2
      elseif(phidot2.lt.0.d0)then
       phidot2 = 1.d-12
      endif
      phidot = sqrt(phidot2)
c  alpha3n = phidot/hbar/omega_n  alpha3p = phidot/hbar/omega_p
      alpha3n = phidot/homegan
      alpha3p = phidot/homegap
c  calculate k_trans= transk 
      transk = r*phidot/hc2om
c  cn = k_trans/kappa_n, cp = k_trans/kappa_p
      cn = transk/cappan
      cp = transk/cappap
c  alf32n = (phidot/hbar/omegan)**2  alf32p = (phidot/hbar/omegap)**2
      alf32n = alpha3n**2
      alf32p = alpha3p**2
c  calculate b_a, b_A
      bn = bna+bnaa
      bp = bpa+bpaa
c  B_n = bbn   B_p = bbp 
      bbn = sqrt(1.d0+ 2.d0*bn*alf32n) 
      bbp = sqrt(1.d0+ 2.d0*bp*alf32p) 
c  D_n = ddn   D_p = ddp 
      ddn = 1.d0 + .5d0*bn*ssl
      ddp = 1.d0 + .5d0*bp*ssl
c  (c_n3)**2 = c3n2   (c_p3)**2 = c3p2 
      c3n2 = (1.d0 + 2.d0*bna*bnaa/bn*ssl)/ddn
      c3p2 = (1.d0 + 2.d0*bpa*bpaa/bp*ssl)/ddp
c  alphs3n = phidot/hbar/omega_n/B_n/D_n
      alphs3n = alpha3n/bbn/ddn
      alphs3p = alpha3p/bbp/ddp
      bhomn = bbn*homegan
      bhomp = bbp*homegap
c  calculate k_n and k_p with smin from capture_s
      smiun = r - rr - vdr
      smintr = smin - (r00-1.20)*(aa3+aaa3)-.18
      dsmi = smintr-smiun
      write(*,'(a,f7.4,a)')
     &' Surface distance increased by ',dsmi,' fm'
      eeen = (capn+cn**2/4.d0*cappan)*smintr
      eeep = (capp+cp**2/4.d0*cappap)*smintr
      if(eeen.gt.70..or.eeep.gt.70.)then
       en = 0.d0
       ep = 0.d0
      else
       en = 1.d0/(1.d0+exp(eeen))
       ep = 1.d0/(1.d0+exp(eeep))
      endif
c  k_n(r) = akn  k_p(r) = akp
      akn = adia*akbarn*bbn/sqrt(ddn)*zeta**2*en
      akp = adia*akbarp*bbp/sqrt(ddp)*zeta**2*ep
c  calculate Q_opt
      qoptns = 0.
      qoptps = qoptns + (za-zaa)*esq/r
c  calculate Qs
      qsns = -(qns-qoptns)/bhomn
      qsnp = -(qnp+qoptns)/bhomn
      qsps = -(qps-qoptps)/bhomp
      qspp = -(qpp+qoptps)/bhomp
      acc = accur/1000.d0
c  J^nS = ejns
      ejns = qsns + (bna-bnaa)*cn*alphs3n
c  exp(-(J^nS)**2) = efacns
      if(abs(ejns).gt.8.)then
       efacns = 0.d0
       iins = 1
      else
       efacns = exp(-ejns**2)
       iins = 0
      endif
c  calculate error functions
      if(ejns.lt.-5.)then
       iins = 0
       call ierfc(0,ejns,g0ns)
       call ierfc(1,ejns,g1ns)
       call ierfc(2,ejns,g2ns)
       call ierfc(3,ejns,g3ns)
      else
       ierr = 0
       call gn_s(ejns,3,gn,acc,ierr)
       if(ierr.eq.1)goto 2
       if(ierr.ne.0)write(*,*)' Warning: no convergence in gn_s for ns'
       gg0ns = gn(1)     
       gg1ns = gn(2)     
       gg2ns = gn(3)     
       gg3ns = gn(4)     
       g0ns = gg0ns*efacns
       g1ns = gg1ns*efacns
       g2ns = gg2ns*efacns
       g3ns = gg3ns*efacns
      endif
c  J^nP = ejnP
      ejnp = qsnp + (bnaa-bna)*cn*alphs3n
c  exp(-(J^nP)**2) = efacnp
      if(abs(ejnp).gt.8.)then
       efacnp = 0.d0
       iinp = 1
      else
       efacnp = exp(-ejnp**2)
       iinp = 0
      endif
      if(ejnp.lt.-5.)then
       iinp = 0
       call ierfc(0,ejnp,g0np)
       call ierfc(1,ejnp,g1np)
       call ierfc(2,ejnp,g2np)
       call ierfc(3,ejnp,g3np)
      else
       ierr = 0
       call gn_s(ejnp,3,gn,acc,ierr)
       if(ierr.eq.1)goto 2
       if(ierr.ne.0)write(*,*)' Warning: no convergence in gn_s for np'
       gg0np = gn(1)     
       gg1np = gn(2)     
       gg2np = gn(3)     
       gg3np = gn(4)     
       g0np = gg0np*efacnp
       g1np = gg1np*efacnp
       g2np = gg2np*efacnp
       g3np = gg3np*efacnp
      endif
c  J^pS = ejps
      ejps = qsps + (bpa-bpaa)*cp*alphs3p
c  exp(-(J^pS)**2) = efacps
      if(abs(ejps).gt.8.)then
       efacps = 0.d0
       iips = 1
      else
       efacps = exp(-ejps**2)
       iips = 0
      endif
      if(ejps.lt.-5.)then
       iips = 0
       call ierfc(0,ejps,g0ps)
       call ierfc(1,ejps,g1ps)
       call ierfc(2,ejps,g2ps)
       call ierfc(3,ejps,g3ps)
      else
       ierr = 0
       call gn_s(ejps,3,gn,acc,ierr)
       if(ierr.eq.1)goto 2
       if(ierr.ne.0)write(*,*)' Warning: no convergence in gn_s for ps'
       gg0ps = gn(1)     
       gg1ps = gn(2)     
       gg2ps = gn(3)     
       gg3ps = gn(4)     
       g0ps = gg0ps*efacps
       g1ps = gg1ps*efacps
       g2ps = gg2ps*efacps
       g3ps = gg3ps*efacps
      endif
c  J^pP = ejpp
      ejpp = qspp + (bpaa-bpa)*cp*alphs3p
c  exp(-(J^pP)**2) = efacpp
      if(abs(ejpp).gt.8.)then
       efacpp = 0.d0
       iipp = 1
      else
       efacpp = exp(-ejpp**2)
       iipp = 0
      endif
      if(ejpp.lt.-5.)then
       iipp = 0
       call ierfc(0,ejpp,g0pp)
       call ierfc(1,ejpp,g1pp)
       call ierfc(2,ejpp,g2pp)
       call ierfc(3,ejpp,g3pp)
      else
       ierr = 0
       call gn_s(ejpp,3,gn,acc,ierr)
       if(ierr.eq.1)goto 2
       if(ierr.ne.0)write(*,*)' Warning: no convergence in gn_s for pp'
       gg0pp = gn(1)     
       gg1pp = gn(2)     
       gg2pp = gn(3)     
       gg3pp = gn(4)     
       g0pp = gg0pp*efacpp
       g1pp = gg1pp*efacpp
       g2pp = gg2pp*efacpp
       g3pp = gg3pp*efacpp
      endif
c  exp(L_n^0) = c0n  exp(L_p^0) = c0p
      c0n = exp((bn/2.d0/ddn-cappan/4.d0*(rr+vdr))*cn**2)
      c0p = exp((bp/2.d0/ddp-cappap/4.d0*(rr+vdr))*cp**2)
c   rr is sum of nuclear radii for the densities (not equal to r0)
c calculate F_0's
      f0ns = akn*c0n*g1ns
      f0np = akn*c0n*g1np 
      f0ps = akp*c0p*g1ps
      f0pp = akp*c0p*g1pp
c
c elastic probabilities
      ep0n = f0ns+f0np 
      ep0p = f0ps+f0pp
      if(ep0n.gt.80.) ep0n = 80.d0
      if(ep0p.gt.80.) ep0p = 80.d0
      p0n = exp(-ep0n)
      p0p = exp(-ep0p)
      p0tr = p0n*p0p
c  i.e. probability for no transfer
      p0ina = exp(-f0ina)
      p0inaa = exp(-f0inaa)
      p0in = exp(-f0in)
c  i.e. probability for no inelastic
      p0a = p0tr*p0ina
c  i.e. probability for projectile in ground state
      p0aa = p0tr*p0inaa
c  i.e. probability for target in ground state
      if(p0tr.gt.1.d0-accur/10.d0.and.icoul.eq.0)p0in = 1.d0
c  i.e. neglect Coulomb excitation for large impact parameters
      p0 = p0tr*p0in
c  i.e. probability for elastic scattering
c      
c calculate mean values and spreads  of mass and charge
      driftz = (f0pp - f0ps)
      driftn = (f0np - f0ns)
      drifta = driftn+driftz
      wid2z = (f0pp + f0ps)
      wid2n = (f0np + f0ns)
      wida = sqrt(wid2n+wid2z)
      widz = sqrt(wid2z)
c   of energies
      if(iins.eq.0)then
       eavns = 2.d0*homegan*bbn*f0ns*g2ns/g1ns
       wid2ns = (homegan*bbn)**2*6.d0*f0ns*g3ns/g1ns
      else
       eavns = 2.d0*homegan*bbn*f0ns*gg2ns/gg1ns
       wid2ns = (homegan*bbn)**2*6.d0*f0ns*gg3ns/gg1ns
      endif
      if(iinp.eq.0)then
       eavnp = 2.d0*homegan*bbn*f0np*g2np/g1np
       wid2np = (homegan*bbn)**2*6.d0*f0np*g3np/g1np
      else
       eavnp = 2.d0*homegan*bbn*f0np*gg2np/gg1np
       wid2np = (homegan*bbn)**2*6.d0*f0np*gg3np/gg1np
      endif
      if(iips.eq.0)then
       eavps = 2.d0*homegap*bbp*f0ps*g2ps/g1ps
       wid2ps = (homegap*bbp)**2*6.d0*f0ps*g3ps/g1ps
      else
       eavps = 2.d0*homegap*bbp*f0ps*gg2ps/gg1ps
       wid2ps = (homegap*bbp)**2*6.d0*f0ps*gg3ps/gg1ps
      endif
      if(iipp.eq.0)then
       eavpp = 2.d0*homegap*bbp*f0pp*g2pp/g1pp
       wid2pp = (homegap*bbp)**2*6.d0*f0pp*g3pp/g1pp
      else
       eavpp = 2.d0*homegap*bbp*f0pp*gg2pp/gg1pp
       wid2pp = (homegap*bbp)**2*6.d0*f0pp*gg3pp/gg1pp
      endif
      eavtr = eavns+eavnp+eavps+eavpp
      eava = eavtr/2.d0 + eavina
      eavaa = eavtr/2.d0 + eavinaa
      eav = eavtr + eavin
      qp1 = -qps
      if(driftz.gt.0.d0)qp1=qpp
      qn1 = -qns
      if(driftn.gt.0.d0)qn1=qnp
      tavs = eav - driftz*qp1 - driftn*qn1
c  mean energies excluding elastic
      if(abs(1.d0-p0).gt.1.d-5)then
       eav0 = eav/(1.d0-p0)
       tavs0 = tavs/(1.d0-p0)
      else
       eav0 = 0.d0
       tavs0 = 0.d0
      endif
      if(abs(1.d0-p0a).gt.1.d-5)then
       eava0 = eava/(1.d0-p0a)
      else
       eava0 = 0.d0
      endif
      if(abs(1.d0-p0aa).gt.1.d-5)then
       eavaa0 = eavaa/(1.d0- p0aa)
      else
       eavaa0 = 0.d0
      endif
c
      wid2tr = wid2ns+wid2np+wid2ps+wid2pp
      wid2ea = wid2tr/3.d0+wid2ina
      wid2eaa = wid2tr/3.d0+wid2inaa
      wid2e = wid2tr+wid2in
      wide = sqrt(wid2e)
      widea = sqrt(wid2ea)
      wideaa = sqrt(wid2eaa)
c  energy spreads except elastic
      if(abs(1.d0-p0).gt.1.d-5)then
       wid2e0 = wid2e/(1.d0-p0)-p0*eav0**2
      else
       wid2e0 = 0.d0
      endif
      if(wid2e0.lt.0.d0)wid2e0 = 0
      if(abs(1.d0-p0a).gt.1.d-5)then
       wid2ea0 = wid2ea/(1.d0-p0a)-p0a*eava0**2
      else
       wid2ea0 = 0.d0
      endif
      if(wid2ea0.lt.0.d0)wid2ea0 = 0
      if(abs(1.d0-p0aa).gt.1.d-5)then
       wid2eaa0 = wid2eaa/(1.d0-p0aa)-p0aa*eavaa0**2
      else
       wid2eaa0 = 0.d0
      endif
      if(wid2eaa0.lt.0.d0)wid2eaa0 = 0
c            
c of angular momenta
      dnaa = 1.d0 + bnaa*ssl
      dna = 1.d0 + bna*ssl
      dpaa = 1.d0 + bpaa*ssl
      dpa = 1.d0 + bpa*ssl
      avmans = bna*akn*c0n*(alphs3n*dnaa*g0ns-cn/ddn*g1ns)
      avmanp = bnaa*akn*c0n*(alphs3n*dna*g0np-cn/ddn*g1np)
      avmaps = bpa*akp*c0p*(alphs3p*dpaa*g0ps-cp/ddp*g1ps)
      avmapp = bpaa*akp*c0p*(alphs3p*dpa*g0pp-cp/ddp*g1pp)
      avmaans = bnaa*akn*c0n*(alphs3n*dna*g0ns+cn/ddn*g1ns)
      avmaanp = bna*akn*c0n*(alphs3n*dnaa*g0np+cn/ddn*g1np)
      avmaaps = bpaa*akp*c0p*(alphs3p*dpa*g0ps+cp/ddp*g1ps)
      avmaapp = bpa*akp*c0p*(alphs3p*dpaa*g0pp+cp/ddp*g1pp)
c
      avmatr = avmans + avmanp + avmaps + avmapp 
      avma = avmatr + avmina
      avmaatr = avmaans + avmaanp + avmaaps + avmaapp 
      avmaa = avmaatr + avminaa
      avm = avma + avmaa
      avmtr = avmatr+avmaatr
c  avmt1n = (b_na+b_nA)*phidot/B_n/omega_n*c_3**2 
      avmt1n = bn*alpha3n/bbn*c3n2
      avmt1p = bp*alpha3p/bbp*c3p2
c  avmr1n = (b_an-b_nA)/D_n*k_trans/kappa_n
      avmr1n = cn*(bna-bnaa)/ddn
      avmr1p = cp*(bpa-bpaa)/ddp
      avmch=2.d0*akn*c0n*(avmt1n+avmr1n*ejns)
     &+2.d0*akp*c0p*(avmt1p+avmr1p*ejps)
      if(iopt.le.3.and.abs(avmtr-avmch).gt.accur)then
       write(*,'(2(a,f12.6/a,f10.5))')
     &' Warning: avmtr = ',avmtr,' avmch = ',avmch
     &,'     for r = ',r
      endif
      avmns = avmans+avmaans
      avmps = avmaps+avmaaps
      avmnp = avmanp+avmaanp
      avmpp = avmapp+avmaapp
c  mean angular momenta except elastic
      if(abs(1.d0-p0a).gt.1.d-5)then
       avma0 = avma/(1.d0-p0a)
      else
       avma0 = 0.d0
      endif
      if(abs(1.d0-p0aa).gt.1.d-5)then
       avmaa0 = avmaa/(1.d0-p0aa)
      else
       avmaa0 = 0.d0
      endif
      if(abs(1.d0-p0).gt.1.d-5)then
       avm0 = avm/(1.d0-p0)      
      else
       avm0 = 0.d0
      endif
      if(ieprep.eq.1)then
       ieprep = 0
       goto 87
      endif
c
c  calculation of change in r
      ierr = 1
      forc = - dueff(al,r,ierr)
      if(ierr.eq.1)goto 2
      changer = (tavs + ((za+driftz)*(zaa-driftz)-za*zaa)*esq/r
     &    -phidot*avm)/forc
     &-al*phidot/2.d0*((aa+drifta)*(aaa-drifta)/aa/aaa-1.d0)/forc
      changers = (homegan/bbn*akn*c0n+ homegap/bbp*akp*c0p
     &     +eavin-phidot*avmin)/forc
c  c1ns = k_n(r)*exp(L_n^0-(J_n^0)**2)  etc.
      c1ns = c0n*akn*efacns
      c1np = c0n*akn*efacnp
      c1ps = c0p*akp*efacps
      c1pp = c0p*akp*efacpp
c
c  elastic scattering angle
      call scatangle_s(al,thetcoul,theta,dthet,ierr)
      if(ierr.eq.1)goto 2
      thet = 180.d0/pi*theta
      thetcou = 180.d0/pi*thetcoul
c  Geometric elastic cross section in millibarn
      snit = sin(theta)
      if(snit.gt.accur)then
       gcs = rho/snit/dthet*10.d0
      else
       gcs = 10.d0/dthet**2
      endif
      cscoul = (za*zaa*esq/4./ecm)**2/(sin(thetcoul/2.d0))**4*10.d0
c  geometric elastic cross section corrected for capture 
c      (but not for reaction probability!)
      cs = gcs*ppscat
      pr0 = 1.d0-p0
      if(abs(1.d0-p0).lt.1.d-5)pr0 = 0.d0
c  reaction cross section dsig_reac/dl
      sigreac0 = 2.d0*pi*rho/alorho*10.d0*ppscat*pr0
c  geometric elastic cross section:
      sigelas0 = cs*p0
      sig0drut = sigelas0/cscoul
      if(abs(aa-aaa).lt.1.d-5) then
       thetlab=3.1415
      else
       thetlab = atan(sin(theta)/(cos(theta)+aa/aaa))
      endif
      if(thetlab.lt.0.)thetlab=thetlab + pi
      thetlabdeg = 180.d0/pi*thetlab
      if(abs(sin(theta)).lt.accur)then
       cslab = cs*(1.d0-aa/aaa)**2
      else
       cslab = cs*(sin(theta)/sin(thetlab))**2/abs(cos(theta-thetlab))
      endif
c
c  The absorptive potential:
      ww=-.5d0/sqrt(pi)*((homegan*(f0ns+f0np)+homegap*(f0ps+f0pp))/adia
     &   +homegain*f0in)
c  This is rather meaningless:
      sigreac0lab = cslab*pr0 
c
c  Population of single particle levels     
      dtann=2.d0/homegan/bbn/g0na*akn*c0n
      dtaann=dtann*g0na/g0naa
      dtanp=2.d0/homegap/bbp/g0pa*akp*c0p
      dtaanp=dtanp*g0pa/g0paa
c  Warning if change in population tooo large
      if(dtann.gt..5.or.dtaann.gt..5.or.dtanp.gt..5
     &.or.dtaanp.gt..5)then
       if(r.gt.rwarn)rwarn = r
      endif
      if(ipopn.ne.0)then
       ierr = 1
       call population_s(ierr)
       if(ierr.eq.1)goto 2
       write(*,'(a)')
     &' Table of single particle population stored in FOR019.DAT'
      endif
c
c  Printout of variables
      if(iout.ne.0)then
       write(*,'(/a,2(f7.2,a))')
     &' Calculating for l = ',al,'( r = ',r,' fm)'
       if(dtann.gt..5.or.dtaann.gt..5)
     &  write(98,'(a/a,4f6.3)')
     &' Warning: Maximum change in neutron single particle',
     &'    populations in projectile  and target = ',dtann,dtaann
       if(dtanp.gt..5.or.dtaanp.gt..5)
     &  write(98,'(a/a,4f6.3)')
     &' Warning: Maximum change in proton single particle',
     &'    populations in projectile and target = ',dtanp,dtaanp
       do i=1,iout
        if(outgen(i).eq.'dta')
     &  write(98,'(a/,4f6.3)')
     &' Maximum change in single particle populations ='
     &,dtann,dtaann,dtanp,dtaanp
        if(outgen(i).eq.'homega')
     &  write(98,'(a,f5.2,a,f5.2,a)')
     &' homega_n =',homegan,' MeV,  homega_p =',homegap,' MeV'
        if(outgen(i).eq.'phidot')
     &  write(98,'(a,f7.2,a)')' hbar*phidot =',phidot,' MeV'
        if(outgen(i).eq.'alpha3')
     &  write(98,'(a,f5.2,a,f5.2)')'  phidot/omega_n =',alpha3n,
     &'   phidot/omega_p =',alpha3p  
        if(outgen(i).eq.'transk')
     &  write(98,'(a,f6.3,a)')' k_trans =',transk,' 1/fm'
        if(outgen(i).eq.'c')
     &  write(98,'(a,f5.3,a,f5.3)')' k_trans/kappa_n =',cn,
     & '   k_trans/kappa_p =',cp
        if(outgen(i).eq.'bb')
     &  write(98,'(a,f6.3,a,f6.3)')' B_n = ',bbn,',  B_p = ',bbp
        if(outgen(i).eq.'dd')
     &  write(98,'(a,f6.3,a,f6.3)')' D_n = ',ddn,',  D_p = ',ddp
        if(outgen(i).eq.'ak')
     &  write(98,'(a,f11.8,a,f11.8)')' k_n =',akn,',  k_p =',akp      
        if(outgen(i).eq.'qopt')
     &  write(98,'(a,f6.3,a,f6.3,a)')
     &' Q_opt^ns =',qoptns,' MeV, Q_opt^ps =',qoptps,' MeV'
        if(outgen(i).eq.'qs')
     &  write(98,'(a,4f6.2)')' Qs_ns, Qs_np, Qs_ps, Qs_pp='
     &,qsns,qsnp,qsps,qspp
        if(outgen(i).eq.'ej')
     &  write(98,'(a,4f6.2)')' J_ns,J_np,J_ps,J_pp=',ejns,ejnp,ejps,ejpp
        if(outgen(i).eq.'c0')
     &  write(98,'(a,f7.2,a,f7.2)')' exp(L_n) =',c0n,' exp(L_p) =',c0p
        if(outgen(i).eq.'f0')then
         write(98,'(a,4f8.4)')' F0_ns,F0_np,F0_ps,F0_pp ='
     &,f0ns,f0np,f0ps,f0pp
         if(inelast.eq.1)write(98,'(a,2f7.2)')' F0_in^a, F0_in^A ='
     &,f0ina,f0inaa
        endif
       enddo
      endif
c
c  Branching according to leg
99    ierr = intprep
      if(leg.eq.1)then
       call leg1_s(nestep,ierr)
      elseif(leg.eq.2)then
       call leg2_s(nestep,ierr)
      elseif(leg.eq.3)then
       call leg3_s(nestep,ierr)
      elseif(leg.eq.4)then
       call leg4_s(nestep,ierr)
      endif
      if(ierr.eq.1)goto 2
      if(ierr.eq.0)goto 11
      if(ierr.eq.4)goto 24
      write(*,*)' leg outside allowed range'
      goto 2
c
c  next step in total (istep and irstep determined in the legs):
24    if(irstep.gt.200.or.istep.gt.200)then
       write(*,*)' program stop because istep or irstep >200'
       goto 2
      endif
      if(irstep.eq.0)then
c  i.e. for close collisions where trajectory is defined by l   
       alpr = al
       rpr = r
       al = al + delal
       goto 51
      elseif(irstep.ne.0)then
c  i.e. distant collisions where trajectory defined by r
       vdr = vdr0
       ierr = 0
       if(accur.gt..0001d0)accur = .0001d0
cnan
cnan   attenzione queste sono state spostate da quelle di sotto
ccna       r = r + delr
ccna       al = alofr(r,ierr)
ccna       if(ierr.eq.1)goto 2
ccna       ierr = 0
cnan
       call transm_s(al,ecm,rcapt,alcapt,ecmcapt,ppcapt,descat,ierr)
       accur = accur0
       if(ierr.eq.1)then
        write(*,*)' Error in transm for l = al'
        goto 2
       endif 
       ppscat = 1.d0-ppcapt

       r = r + delr
       ierr = 0
       al = alofr(r,ierr)
       if(ierr.eq.1)goto 2

       goto 57      
      endif
c
c  Printout for leg = 0
9     dum = '------'
      if(norbit.eq.2)then
       if(icapt.eq.2)then
        if(decm.gt.accur.and.ecm+decm.gt.ecmcapt)then
         write(89,'(f6.2,5a8,f8.2,f10.5)')
     &excom,dum,dum,dum,dum,dum,widcbar,sigcest
         ecx=ecmcapt+qvalcom
         alx=0.d0
         alatrx=0.d0
         write(89,'(f6.2,3f8.2,4a8)')
     &ecx,alx,rcapt,alatrx,dum,dum,dum,dum
         if(ecm+decm.gt.eb)then
          alggg = 0.d0
          ecx=eb+qvalcom
          write(89,'(f6.2,3a8,f8.2,3a8)')
     &ecx,dum,dum,dum,alggg,dum,dum,dum
         endif         
        elseif(decm.lt.accur.and.ecm-decm.gt.ecmcapt)then
         if(ecm-decm.gt.eb)then
          alggg = 0.d0
          ecx=eb+qvalcom
          write(89,'(f6.2,3a8,f8.2,3a8)')
     &ecx,dum,dum,dum,alggg,dum,dum,dum
         endif         
         ecx=ecmcapt+qvalcom
         alx=0.d0
         alatrx=0.d0
         write(89,'(f6.2,3f8.2,4a8)')
     &ecx,alx,rcapt,alatrx,dum,dum,dum,dum
         write(89,'(f6.2,5a8,f8.2,f10.5)')
     &excom,dum,dum,dum,dum,dum,widcbar,sigcest
        endif
       elseif(icapt.eq.1)then
        if(decm.gt.accur.and.ecm+decm.gt.eb)then
         write(89,'(f6.2,3f8.2,2a8,f8.2,f10.5)')
     &excom,alcapt,rcapt,alatr,dum,dum,widcbar,sigcest
         alggg = 0.d0
         ecx=eb+qvalcom
         write(89,'(f6.2,3a8,f8.2,3a8)')
     &ecx,dum,dum,dum,alggg,dum,dum,dum
        elseif(decm.lt.accur.and.ecm-decm.gt.eb)then
         alggg = 0.d0
         ecx=eb+qvalcom
         write(89,'(f6.2,3a8,f8.2,3a8)')
     &ecx,dum,dum,dum,alggg,dum,dum,dum
         write(89,'(f6.2,3f8.2,2a8,f8.2,f10.5)')
     &excom,alcapt,rcapt,alatr,dum,dum,widcbar,sigcest
        endif         
       endif 
      else
       write(89,'(f6.2,6f8.2,f10.5)')
     &excom,alcapt,rcapt,alatr,algraz,rgraz,widcbar,sigcest
      endif
      if(decm.lt.-accur)then
       if(abs(ecm-ecminit).lt.accur)then
        goto 13
       else
        goto 11
       endif
      elseif(decm.gt.accur)then
       if(abs(ecm-ecmm).lt.accur)then
        goto 13
       else
        goto 11
       endif
      elseif(abs(decm).le.accur)then
       goto 13
      endif
13    delyrl=sqrt((aa+aaa)**1.666d0*excom/30.d0)/dble(30)
      yefac= 34.7d0/(aa+aaa)**1.6666d0
      do n=0,30
       yrl = dble(n)*delyrl
       yex=yefac*yrl**2
       yexn = yex+threshn
       yexp = yex+thbarp
       write(88,'(4f10.2)')yrl,yex,yexn,yexp
      enddo
11    if(decm.eq.0.d0)then
         goto 81
        else
         ecm = ecm + decm
         delr=dllr
         delrinit=dllr
         elab = ecm*(1.d0+aa/aaa)
         epslab = elab/aa
         if(decm.lt.0.d0.and.ecm.lt.ecmm-accur.or.decm.gt.0.d0.and
     &       .ecm.gt.ecmm+accur)then
           tend = X05BAF()
           cpu = tend-tstart
        write(*,*)' -------------------------------------------------'
        write(98,*)' -------------------------------------------------'
        write(98,'(a,f8.2,a)')' CPU-time for this run was ',cpu,' sec.'
        write(*,'(a,f8.2,a)')' CPU-time for this run was ',cpu,' sec.'
           goto 91
         endif
         write(98,'(/a,f7.1,a)') ' E_cm =',  ecm, ' Mev'
         write(*,'(/a,f7.1,a)') ' E_cm =',  ecm, ' Mev'
         write(98,'(a,f7.1,a)')' E_lab =' ,elab,'  MeV'
         write(*,'(a,f7.1,a)')' E_lab =' ,elab,'  MeV'
         tend = X05BAF()
         cpu = tend-tstart
         if(abs(-ecm+ecminit+decm).lt.accur.and.leg.ne.0)then
           timeest = cpu*dble(neestim)/60.d0
       write(*,'(a,f6.1,a)')' Estimated time for job = ',timeest,' min'
           if(timeest.gt.timemax*.7d0)then
            write(*,*)' The job will presumably exceed the time alloted'
             goto 2
           endif
         endif
         nestep = nestep + 1
         if(intprep.eq.1.and.nestep.eq.3)goto 91
         goto 80
      endif
81    if(intp.eq.'inter')write(*,'(/a)')
     &' Do you want to change the energy? yes or no'
      read(i5,'(a)')out
      addtext = '      ! New energy'
      if(intp.eq.'inter')write(10,'(a,a)')out,addtext
      if(out.eq.'no ')then
       goto 91
      elseif(out.eq.'yes')then
       write(*,*)
     &' --------------------------------------------------------------' 
       write(98,*)
     &' --------------------------------------------------------------' 
       nestep = nestep + 1
       goto 22
      else
       goto 81
      endif
91    if(intp.eq.'inter')write(*,'(/a)')
     &' Do you want make a new run? yes or no'
      read(i5,'(a)')out
      addtext = '      ! new run enters at Z_a,A_aetc'
      if(intp.eq.'inter')write(10,'(a,a)')out,addtext
      if(out.eq.'no ')then
       goto 2
      elseif(out.eq.'yes')then
       write(*,*)
     &' --------------------------------------------------------------' 
       write(*,*)
     &' --------------------------------------------------------------' 
       write(98,*)
     &' --------------------------------------------------------------' 
       write(98,*)
     &' --------------------------------------------------------------' 
       goto 17
      else
       goto 91
      endif
2     tend = X05BAF()
      cpu = tend-tinit
      icpumin = int(cpu/60.d0)
      icpusec = int(mod(cpu,60.d0))
      write(*,*)' ---------------------------------------------------' 
      write(98,*)' ---------------------------------------------------'
      if(cpu.gt..01)then
       write(98,'(a,i3,a,i2,a)')' CPU time for this job was ',icpumin,
     &' min ',icpusec,' sec.'
       write(*,'(a,i3,a,i2,a)')' CPU time for this job was ',icpumin,
     &' min ',icpusec,' sec.'
      endif
      if(intp.eq.'inter'.and.version.eq.'vms  ')then 
       runfile1 = 'grazing_run'//addin//'.com'
       call condens_s(runfile1,runfile,leru)
       open (unit=97,file=runfile,status='unknown')
       write(97,'(a,a)')'$ set def ',defaultname
       write(97,'(a,a)')'$ run ',programname
       write(97,'(a11)')autoadd
       write(97,'(a)')'$ stop'
       if(intprep.eq.0)write(*,'(/a,3(/a,a21))')
     &' You may rerun the program in batch:',
     &'     a) with the same input by running the file:',runfile,
     &'     b) with changed data by modifying the file ',infilename,
     &'        and then running the command file; ',runfile
       if(intprep.eq.1)write(*,'(/a/a,a21/a,2(/a,a21))')
     &' The output that you created was incomplete,',
     &'  but, you have prepared an input file:',infilename,
     &'  and may now run the program, either',
     &'     a) "interactively" by running grazing entering ',autoadd,
     &'     b) in batch: running the command file; ',runfile
      elseif(intp.eq.'inter'.and.version.eq.'unix* ')then
       if(intprep.eq.0)write(*,'(/a,/a,a21)')
     &' You may rerun the program with the same or modified',
     &'  input data by using the input file ',infilename 
       if(intprep.eq.1)write(*,'(/a/a,a21/a/a,a21)')
     &' The output that you created was incomplete,',
     &'  but, you have prepared an input file:',infilename,
     &'  and you may run the program with the same or modified',
     &'  input data by using the input file ',infilename
      endif
      do nn = 5,98
       close(unit=nn)
      enddo 
      end
