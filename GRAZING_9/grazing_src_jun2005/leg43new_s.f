c      Subroutine LEG 43 for program grazing
c      At exit ierr = 1 for error
c              ierr = 2 for warning
c              ierr = 3 for differential cross section
c              ierr = 4 for total cross section and new value of r
c              ierr = 0 for end or new value of energy
c---------------------------------------------------------------------
      subroutine leg43_s(nestep,ierr)
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
      common/selec/ecutod,pcut,isel(4),lamsel(4),msmax
      dimension sigetot(ienergy),sig(irangen,irangep)
     &,sigeapptot(ienergy),signzetot(list,ienergy)
     &,signzeapptot(list,ienergy),mnzmax(list),iml(list)
      dimension thetcorr(list,ienergy),thetlabcorr(list,ienergy)
     &,cscorr(list,ienergy),cslabcorr(list,ienergy)
      dimension enzloss(list,ienergy),
     &pnze(list,ienergy),pnzeapp(list,ienergy),delel(list)
      dimension eloss(ienergy),ppe(ienergy),ppapp(ienergy),istore(5)
      character*5 intp
      character*3 out
      character*7 outleg4(20)
      character*44 nzfilename
      character*16 head
      character*10 xname(5)
      character*50 addtext
      character*5 elem1(list),elem2(list)
      common/element/elem1,elem2
      character*5 version
      common/vers/version
c
      iprep = ierr
      if(istep.ne.1)goto 10 
      if(nestep.ne.1)goto 9  
44    if(intp.eq.'inter')write(*,'(4(/a))')
     &' Do you want the energy distribution ',
     &'  of the ejectile or of the recoil?',
     &' Enter 1 for ejectile or 2 for recoil'
      read(i5,*,err=1100)iq
      icheck = 0
      goto 1110
1100  write(*,*)' you must provide a number'
      icheck = icheck + 1
      if(icheck.gt.5)goto 2
      goto 44
1110  addtext = '     !  1 = ejectile, 2 = recoil'
      if(intp.eq.'inter')write(10,'(i5,a)')iq,addtext
      if(iq.ne.1.and.iq.ne.2)then
       write(*,*)' Try again'
       goto 44
      endif
45    istore(1) = 60
      xname(1) = 'edis'
      lout = 1
      ierr = iq
      call selectnz_s(lnzmax,istore,xname,lout,accur,ierr)
      if(ierr.eq.1)goto 2
405   if(intp.eq.'inter')write(*,'(3(/a))')
     &' You may calculate the energy distributions',
     &'  with a finite energy resolution tres(FWHM).',
     &' Enter tres in MeV'
      read(i5,*,err=1101)tres
      icheck = 0
      goto 1111
1101  write(*,*)' you must provide a number'
      icheck = icheck + 1
      if(icheck.gt.5)goto 2
      goto 405
1111  addtext = '    !energy resolution(FWHM) in MeV'
      if(intp.eq.'inter')write(10,'(f7.4,a)')tres,addtext
      if(tres.lt.0.d0)then
       write(*,*)' Try again'
       goto 405
      elseif(tres.lt..2d0)then
       if(intp.eq.'inter')write(*,'(2(/a))')
     &' If you have selected delta_N=delta_Z=0,',
     &'   tres will be set equal to .2 for this partition'
      endif
c  Here total cross sections
      write(*,'(/a/a/a,f7.3,a/a/a)')
     &' You calculate total cross section distributions in energy',
     &' for selected N and Z integrated over L ',
     &' Stepwidth in r is delr =',delr,' fm',
     &' Standard output is  sig_reac, sig_merg',
     &' as well as total cross sections dsig_NZ/dE* or dsig_NZ/dT.'     
406   if(version.eq.'vms  '.or.version.eq.'unix*')then
       if(intp.eq.'inter')write(*,'(a/a)')
     &' You may obtain intermediate results for each value of r',
     &' stored in grazing_out.dat file. Enter yes, or no if not.'
      elseif(version.eq.'unix ')then
       if(intp.eq.'inter')write(*,'(a/a)')
     &' You may obtain intermediate results for each value of r',
     &' stored in FOR098.DAT file. Enter yes, or no if not.'
      endif
      read(i5,'(a)')out
      addtext = '    !  intermediate results'
      if(intp.eq.'inter')write(10,'(a,a)')out,addtext
      if(out.eq.'no ')then
       leg4out = 0
      elseif(out.ne.'no '.and.out.ne.'yes')then
       goto 406
      elseif(out.eq.'yes')then
       leg4out = 1
407    if(intp.eq.'inter')write(*,'(a/a/a/a/a/a/a/a)')
     &' You now obtain for each value of r:',
     &' The scattering angle, dsig_el/dsig_Ruth, dsig_reac/dO',
     &' as well as cross sections d^2sig_NZ/dO/dE* or d^2sig_NZ/dO/dT',
     &' and the differential reaction cross section.',
     &' In addition you may obtain the variables from',
     &' the following list (in paranthesis its significance):',
     &' Enter name of variable, no if no more',
     &' (not implemented! enter no)'
       read(i5,'(a)')outleg4(leg4out)
       addtext = '    !  not implemented, enter no'
       if(intp.eq.'inter')write(10,'(a,a)')outleg4(leg4out),addtext
       if(outleg4(leg4out).eq.'no')then
        outleg4(leg4out) = 'edist'
        goto 9  
       endif
       leg4out = leg4out+1
       if(leg4out.gt.9)then
        write(*,*)' leg4out outside range'
        ierr = 1
        goto 2
       endif
      endif
c preparation for integrating cross section
9     do l=1,list
       do mm=1,ienergy
        signzetot(l,mm)=0.d0
        signzeapptot(l,mm)=0.d0  
       enddo
      enddo
      cross = 0.d0
      crossres = 0.d0
      sigtar = 0.d0
      sigproj = 0.d0
      do n=1,irangen
       do m=1,irangep
        sig(n,m)=0.d+0
       enddo
      enddo
c calculate energy distribution of projectile for given N and Z
10    do l=1,lnzmax
       isn = iselectn(l)
       isz = iselectz(l)
       if(iq.eq.1)write(*,'(a,a5/,a,i3,a,i3,a)')
     &'  Calculates distribution in energy of ejectile ',elem1(l),
     &' for delta N =',isn,' and delta Z =',isz,' (of projectile)'
       if(iq.eq.2)write(*,'(a,a5/,a,i3,a,i3,a)')
     &'  Calculates distribution in energy of recoil ',elem2(l),
     &' for delta N =',isn,' and delta Z =',isz,' (of projectile)'
       timeedist1 = X05BAF()
       delta = tres/2.35482d0
       if(isn.eq.0.and.isz.eq.0.and.tres.lt..2d0)delta = .2/2.35482d0
       eavpr = eavdist(in00+isn,ip00+isz)/2.d0
       wid2epr = wid2adist(in00+isn,ip00+isz)
       pnz = ppnz(in00+isn,ip00+isz)
       if(isn.eq.0.and.isz.eq.0)pnz=pnz+p0
       if(pnz.lt.accur**2)then
        write(*,*)' Cross section too small for this accuracy'
        goto 441
       endif
       widepr = sqrt(wid2epr)
       if(istep.eq.1)then
        deleloss = 2.d0*(widepr+delta)*accur**(.25d0)
c        if(isn.eq.0.and.isz.eq.0)deleloss = delta
        delel(l) = dble(1.d0+int(10.*deleloss))/10.d0
        iml(l) = 0
       endif
       dele=delel(l)
       imm = iml(l)
       ierr = iprep
       call nzedisa_s(iq,delta,dele,accur,imm,mmmax,mmmin
     &,isn,isz,eloss,ppe,ierr)
c  mmmax/mmmin = maximum/minimum of mm tabulated 
       if(ierr.eq.1)goto 2
       if(ierr.gt.1)write(*,*)'     for r = ',r
       ierr = iprep
       eavap = eavdist(isn+in00,isz+ip00)/2.d0
       wid2ap = wid2adist(isn+in00,isz+ip00)
       if(isn.eq.0.and.isz.eq.0)then
        eavap = eavdist00/2.d0
        wid2ap = wid2dist00/3.d0
       endif
       p0ap = p0tr
       deltae = delta
       call edistapp_s(iq,deltae,eavap,wid2ap,p0ap,dele,accur,
     &imm,mmmax,mmmin,eloss,ppapp,ierr)
       if(ierr.gt.1)write(*,*)'     for r = ',r
       do mm = mmmin,mmmax
        write(*,'(a,f8.3,2(a,f7.5))')
     &' eloss = ',eloss(mm),' ppe = ',ppe(mm),' ppapp = ',ppapp(mm)
       enddo
       do ii = 1,ienergy
        enzloss(l,ii) = eloss(ii)
        pnze(l,ii) = ppe(ii)
        pnzeapp(l,ii) = ppapp(ii)
        mnzmax(l)=mmmax
       enddo
441    continue      
       timeedist2 = X05BAF()
       timeedist = timeedist2 - timeedist1
       write(*,'(a,i3,a,f7.2,a/a,i3,a)')
     &' CPU time for',l,'. energy distribution =',timeedist,' sec.',
     &'   You ordered',lnzmax,' distributions'
      enddo
c
c calculate corrected scattering angle 
      dsigtar = abs(cs)*p0tr*p0ina*(1.d0-p0inaa)
      dsigproj = abs(cs)*p0tr*p0inaa*(1.d0-p0ina)
      if(leg4out.eq.0)goto 27
      write(*,'(/a,f5.2,a)')'   Calculation for r =', r,' fm'
      write(98,'(a,f7.2,a,f5.2,a)')
     &'  Results for l = ',al,' (r = ', r,' fm)'
      write(98,'(/a,f7.2,a,f5.2,a)')
     &'  Results for l = ',al,' (r = ', r,' fm)'
      if(lab.ne.1)then
       write(98,'(a,f7.2,a)')
     &' Elastic scattering angle =',thet,' degrees'
       write(*,'(a,f7.2,a)')
     &' Elastic scattering angle =',thet,' degrees'
      elseif(lab.ne.0)then
       write(98,'(a,f7.2,a)')
     &' Elastic laboratory scattering angle=',thetlabdeg,' degrees'   
       write(*,'(a,f7.2,a)')
     &' Elastic laboratory scattering angle=',thetlabdeg,' degrees'
      endif
      write(98,'(a,d11.4)')
     &' Elastic cross section/Rutherford =',sig0drut
      write(*,'(a,d11.4)')
     &' Elastic cross section/Rutherford =',sig0drut
      if(lab.ne.1)then
       write(98,'(a,d11.4,a)')          
     &' Reaction cross section =',sigreac0,' mb/sterad'
       write(*,'(a,d11.4,a)')
     &' Reaction cross section =',sigreac0,' mb/sterad'
      elseif(lab.ne.0)then
       write(*,'(a,d11.4,a)')
     &' Reaction cross section  i lab = ',sigreac0lab,' mb/sterad'
       write(98,'(a,d11.4,a)')
     &' Reaction cross section in lab = ',sigreac0lab,' mb/sterad'
      endif
c
81    do l=1,lnzmax 
       ll = 60 + l
       isn = iselectn(l)
       n = in00 + isn
       isz = iselectz(l)
       m = ip00 + isz
       if(istep.eq.1)then
        if(iq.eq.1)then
         write(98,'(a,a5)')
     &' Distribution in excitation energy of ejectile ',elem1(l)
         write(ll,'(a,a5)')
     &' Distribution in excitation energy of ejectile ',elem1(l)
        elseif(iq.eq.2)then
         write(98,'(a,a5)')
     &' Distribution in excitation energy of recoil ',elem2(l)
         write(ll,'(a,a5)')
     &' Distribution in excitation energy of recoil ',elem2(l)
        endif
        write(98,'(a,i3,a,i3)')
     &'  with delta N =',isn,' and delta Z =',isz
        if(tres.gt.accur)write(98,'(a,f5.2,a)')
     &'  and energy resolution ',tres,' Mev'
        write(ll,'(a,i3,a,i3)')
     &'  with delta N =',isn,' and delta Z =',isz
        if(tres.gt.accur)write(ll,'(a,f5.2,a)')
     &'  and energy resolution ',tres,' Mev'
        if(lab.eq.1)then
         write(98,*)'   in laboratory system'
         write(ll,*)'   in laboratory system'
        endif
        if(leg4.eq.3)then
         write(98,'(a)') '      E*_b   P(E*_b)   P(E*_b)  dsig/dO/dE'
         write(98,'(a)') '      MeV     /MeV     (appr.)   mb/sr/MeV'
         write(ll,'(a)') '      E*_b   P(E*_b)   P(E*_b)  dsig/dO/dE'
         write(ll,'(a)') '      MeV     /MeV     (appr.)   mb/sr/MeV'
        endif
       endif
       do mm = 1,mnzmax(l)
        if(lab.eq.0)then
         dsigcm = ppnz(n,m)*pnze(l,mm)*cs
         write(98,'(f10.1,2f10.4,e12.3)')enzloss(l,mm),pnze(l,mm)
     &,pnzeapp(l,mm),dsigcm
         write(ll,'(f10.1,2f10.4,e12.3)')enzloss(l,mm),pnze(l,mm)
     &,pnzeapp(l,mm),dsigcm
        elseif(lab.eq.1)then
         dsigl = ppnz(n,m)*pnze(l,mm)*cslab
         write(98,'(f10.1,2f10.4,e12.3)')enzloss(l,mm),pnze(l,mm)
     &,pnzeapp(l,mm),dsigl
         write(ll,'(f10.1,2f10.4,e12.3)')enzloss(l,mm),pnze(l,mm)
     &,pnzeapp(l,mm),dsigl
        endif
       enddo
      enddo
27    continue
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
      sigtar = sigtar + factor*(p0a-p0)*ppscat
      sigproj = sigproj + factor*(p0aa-p0)*ppscat
      do n=1,mmmax
       sigeapptot(n)=sigeapptot(n)+factor*ppapp(n)*ppscat
       sigetot(n)=sigetot(n)+factor*ppe(n)*ppscat
      enddo
      do l=1,lnzmax
       n=in00+iselectn(l)
       m=ip00+iselectz(l)
       pppnz = ppnz(n,m)
       do mm=1,mnzmax(l)
        signzetot(l,mm) = signzetot(l,mm) + factor*pppnz*pnze(l,mm)*
     &ppscat
        signzeapptot(l,mm) = signzeapptot(l,mm) + factor*pppnz*
     & pnzeapp(l,mm)*ppscat
       enddo
      enddo
432   do i=inmin,inmax
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
c printout of total cross sections 
33    write(*,*)'------------------------------------------------'
      write(98,*)'------------------------------------------------'
      do l=1,lnzmax
       isn = iselectn(l)
       isz = iselectz(l)
       nndat = 40 + l
       if(iq.eq.1)then
        write(98,'(a,a5)')
     &' Total excitation energy cross section of ejectile ',elem1(l)
        write(*,'(a,a5)')
     &' Total excitation energy cross section of ejectile ',elem1(l)
        write(nndat,'(a,a5)')
     &' Total excitation energy cross section of ejectile ',elem1(l)
       elseif(iq.eq.2)then
        write(98,'(a,a5)')
     &' Total excitation energy cross section of recoil ',elem2(l)
        write(*,'(a,a5)')
     &' Total excitation energy cross section of recoil ',elem2(l)
        write(nndat,'(a,a5)')
     &' Total excitation energy cross section of recoil ',elem2(l)
       endif
       write(98,'(a,i4,a,i4)')
     &'  with delta_N =',isn,' and delta_Z =',isz
       write(*,'(a,i4,a,i4)')
     &'  with delta_N =',isn,' and delta_Z =',isz
       write(nndat,'(a,i4,a,i4)')
     &'  with delta_N =',isn,' and delta_Z =',isz
       if(tres.gt.accur)then
        write(98,'(a,f5.2,a)')
     &'  and energy resolution ',tres,' Mev'
        write(*,'(a,f5.2,a)')
     &'  and energy resolution ',tres,' Mev'
        write(nndat,'(a,f5.2,a)')
     &'  and energy resolution ',tres,' Mev'
       endif
       write(98,'(a)')
     &'  E*_b MeV  signze mb/MeV  signzeapptot mb/MeV'
       write(*,'(a)')
     &'  E*_b MeV  signze mb/MeV  signzeapptot mb/MeV'
       write(nndat,'(a)')
     &'  E*_b MeV  signze mb/MeV  signzeapptot mb/MeV'
       do n=1,ienergy
        if(n.eq.1.or.signzetot(l,n).gt.accur)then
         write(98,'(f7.1,2d15.4)') enzloss(l,n),signzetot(l,n)
     &,signzeapptot(l,n)
         write(*,'(f7.1,2d15.4)') enzloss(l,n),signzetot(l,n)
     &,signzeapptot(l,n) 
         write(nndat,'(f7.1,2d15.4)') enzloss(l,n),signzetot(l,n)
     &,signzeapptot(l,n)
        endif
       enddo
      enddo
      write(*,'(a)')'----------------------------------------------'
      write(98,'(a)')'----------------------------------------------'
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
c      if(tres.gt.accur)then
c       write(*,'(/a,e12.4,a)')
c     &' "Total reaction cross section" = ',sigtotlim,'  mbarn'
c       write(*,'(a,f6.2,a/a,f6.1,a)')
c     &'       with energy resolution = ',tres,' MeV',
c     &'       for scattering angles larger than ',scatmin,' degrees'
c       write(98,'(/a,e12.4,a)')
c     &' "Total reaction cross section" = ',sigtotlim,'  mbarn'
c       write(98,'(a,f6.2,a/a,f6.1,a)')
c     &'       with energy resolution = ',tres,' MeV',
c     &'       for scattering angles larger than ',scatmin,' degrees'
c      endif
      ierr = 0
2     return
      end










