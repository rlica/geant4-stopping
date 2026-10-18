c      Subroutine LEG 4 for program grazing
c      At exit ierr = 1 for error
c              ierr = 2 for warning
c              ierr = 3 for differential cross section
c              ierr = 4 for total cross section and new value of r
c              ierr = 0 for end or new value of energy
c      subroutines: leg41_s,leg42_s,leg43_s,leg44_s,leg45_s
c---------------------------------------------------------------------
      subroutine leg4_s(nestep,ierr)
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
      common/nztdis/eav00,wid200
      common/char/elemp,elemt
      common/head/head
      common/leg/ecminit,rainbcdg,decm,cs,r,delr,al,theta,
     &cscoul,cslab,thet,thetlabdeg,thetcou,ppcapt,ww,pr0,
     &drifta,driftz,wida,widz,eava,eavaa,tavs,widea,wideaa,avmin,avm,
     &changer,changers,ppscat,totmerg,dthet,
     &deltae,sig0drut,sigreac0,sigreac0lab,cross,alorho,rho,
     &istep,iout,ipopn,norbit,irstep,lab
      common/tot/in00,ip00,i0nmin,i0nmax,i0zmin,i0zmax
      common/leg4/dsigz(irangep),dsigzlab(irangep),
     &ppnz(irangen,irangep),iselectn(list),iselectz(list),
     &dsig(irangen,irangep),dsiglab(irangen,irangep),
     &inmin,inmax,izmin,izmax,leg4
      character*5 elemp
      character*5 elemt
      character*5 intp
      character*16 head
      character*40 addtext
c
c  calculate distribution over N and Z and mean energies and spreads
      intprep = ierr
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
        if(i.eq.0.and.j.eq.0)then
         ppnz(n,m)=(p0nn(n)*p0zz(m)-p0)
c  ppnz(in00,ip00) is probability for i=j=0 minus elastic-elastic
         if(ppnz(n,m).lt.accur**2)ppnz(n,m)=0.d0
        endif
       enddo
      enddo
cnan
      inwrite=min(inmax,8)
      write(65,'(f6.1,1p9e11.3)') al,(p0nn(in00+i),i=0,inwrite)         
      inwrite=max(izmin,-8)
      write(66,'(f6.1,1p9e11.3)') al,(p0zz(ip00+i),i=0,inwrite,-1)         
cnan
c
      pfac = p0nn(in00)*p0zz(ip00)
      if(ppnz(in00,ip00).ne.0.d0)then
       fact = pfac/ppnz(in00,ip00)
       eav00 = (eavdist(in00,ip00)+eavina+eavinaa)*fact
       wid200 = (wid2dist(in00,ip00)+wid2ina+wid2inaa
     &+(eavdist(in00,ip00)+eavina+eavinaa)**2)*fact-eav00**2
       avm00 = (avmdist(in00,ip00)+avmin)*fact
c  i.e. mean energy and width of nonelastic entrance channel
       if(wid200.lt.0.d0)then
        write(*,*)' Warning: wid200.lt.0'
        width00 = 0.d0
       endif
      else
       eav00 = 0.d0
       wid200 = 0.d0
       avm00 = 0.d0
      endif
      wid00 = sqrt(wid200)
      width00=2.35482d0*wid00
      p00a = pfac-p0tr*p0ina
c  i.e. probability for i=j=0 minus ejectile elastic
      if(p00a.ge.accur**2)then
       facta = pfac/p00a
       eava00 = (eavdist(in00,ip00)/2+eavina)*facta
       wid2a00 = (wid2adist(in00,ip00)+wid2ina
     &+(eavdist(in00,ip00)/2.d0+eavina)**2)*facta-eava00**2
       avma00 = (avmadist(in00,ip00)+avmina)*facta
c  i.e. mean energy and width of nonelastic ejectile
       if(wid2a00.lt.0.d0)then
        write(*,*)' Warning: wid2a00.lt.0'
        wid2a00 = 0.d-9
       endif
      else
       eava00 = 0.d0
       wid2a00 = 0.d0
       avma00 = 0.d0
      endif
      wida00 = sqrt(wid2a00)
      widtha00=2.35482d0*wida00
      p00aa = pfac-p0tr*p0inaa
c  i.e. probability for i=j=0 minus recoil elastic
      if(p00aa.ge.accur**2/10.d0)then
       factaa = pfac/p00aa
       eavaa00 = (eavdist(in00,ip00)/2+eavinaa)*factaa
       wid2aa00 = (wid2adist(in00,ip00)+wid2inaa
     &+(eavdist(in00,ip00)/2.d0+eavinaa)**2)*factaa-eavaa00**2
       avmaa00 = (avmdist(in00,ip00)-avmadist(in00,ip00)+avminaa)*factaa
c  i.e. mean energy and width of nonelastic recoil
       if(wid2aa00.lt.0.d0)then
        write(*,*)' Warning: wid2aa00.lt.0'
        wid2aa00 = 0.d-9
       endif
      else
       eava00 = 0.d0
       wid2a00 = 0.d0
       avmaa00 = 0.d0
      endif
      widaa00 = sqrt(wid2aa00)
      widthaa00=2.35482d0*widaa00
c4 calculate cross sections
      do n=1,irangen
       do m=1,irangep
        dsig(n,m) = 0.d0
        dsiglab(n,m) = 0.d0
        dsigz(m) = 0.d0
        dsigzlab(m) = 0.d0
       enddo
      enddo
      do i=inmin,inmax 
       do j=izmin,izmax 
        n = in00+i
        m = ip00+j
        dsiglab(n,m)=cslab*ppnz(n,m)
        dsig(n,m)=cs*ppnz(n,m)
c  cs and cslab include ppscat
       enddo
      enddo
c
      if(nestep.ne.1.or.istep.ne.1)goto 50
   69 continue 
      if(intp.eq.'inter')write(*,'(11(a/))')
     &' You have the following choices:',
     &' 1. Distribution in mass and charge of ejectiles or recoils',
     &'       corrected for actual masses',
     &' 2. Distribution in mass and charge of ejectiles or recoils',
     &'       after evaporation',
     &' 3. Distribution in excitation energy of ejectiles or recoil',
     &'       for selected N and Z, numerically (uncorrected)',
     &' 4. Distribution in total kinetic energy for each partition',
     &' 5. Table for plotting the total kinetic energy',
     &'       versus A of the ejectile',
     &' Enter your choice: (exit on 0)'
       read(i5,*,err=1101)leg4
       icheck = 0
       goto 1111
1101   write(*,*)' you must provide a number'
       icheck = icheck + 1
       if(icheck.gt.5)goto 2
       goto 69
1111   addtext = '  !leg4 cannot be changed!'
       if(intp.eq.'inter')write(10,'(i5,a)')leg4,addtext
       if(leg4.gt.5.or.leg4.lt.0)goto 69
       if(leg4.eq.0)then
        ierr = 1
        goto 2
       endif
cnan      endif
50    ierr = intprep
      if(leg4.eq.1)call leg41_s(nestep,ierr)
      if(leg4.eq.2)call leg42_s(nestep,ierr)
      if(leg4.eq.3)call leg43_s(nestep,ierr)
      if(leg4.eq.4)call leg44_s(nestep,ierr)
      if(leg4.eq.5)call leg45_s(nestep,ierr)
2     return
      end
