c      Subroutine calculates the correction to an energy
c       distribution due to selected inelastic processes
c       as specified in common/selec/
c      imm is input while mmmax and mmmin are changed
c      Correction is done by interpolation in stepsize delx 
c       and shifting the distribution ppe
c      Input distribution from eloss(mm),ppe(mm)
c      Output will be in the same arrays
c      ierr at input operates printout 
c----------------------------------------------------------------
      subroutine distriin_s(eloss,ppe,mmmax,mmmin,delx,ierr)
      include 'grazing.icl'
      implicit real*8(a-h,o-z)
      common/ecm/ecm
      common/za/za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/vers/version
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/ebe23/eex(2,3,2)
      common/widin/wi2in(2,3,2)
      common/f0in/ff0in(2,3,2)
      common/eav/eav,eavns,eavnp,eavps,eavpp,eavtr,
     &eav0,eava0,eavaa0
      common/selec/ecutod,pcut,isel(4),lamsel(4),msmax
      dimension eben(4),eta2(4)
      dimension ppinr(1000),excit(1000),ppein(1000)
     &,eloss(ienergy),ppe(ienergy)
      character*5 version
c
      if(version.eq.'vms  ')then
       open (unit=98,file='grazing_out.dat',status='old')
      else
      endif
      iout = ierr
      if(iout.eq.1)open(unit=16,status='new')
      iwarn = 0
      do ms = 1,4
       if(ms.le.msmax)then
        eben(ms)= eex(isel(ms),lamsel(ms),1)
        eta2(ms)= ff0in(isel(ms),lamsel(ms),1)
        write(*,'(2(a,i2),2(a,f7.4))')
     &' i = ',isel(ms),' lambda = ',lamsel(ms),' E_ex = ',eben(ms),
     &' P = ',eta2(ms)
       else
        eben(ms)= accur/10.d0
        eta2(ms)= accur/10.d0
       endif
      enddo
      ierr = 0
      call interpolate_s(eloss,ppe,mmmax,excit,ppinr,delx,ierr)
      if(ierr.eq.1)goto 2
      nma = int((eloss(mmmax)-eloss(1))/delx+1.01)
c      iim = int((eloss(imm)-eloss(1))/delx+1.01)
      do mm=1,1000
       ppein(mm)=0.d0
      enddo
      n1max = 20
      n2max = 20
      n3max = 20
      n4max = 20
      ex1 = exp(-eta2(1))
      ex2 = exp(-eta2(2))
      ex3 = exp(-eta2(3))
      ex4 = exp(-eta2(4))
      acc=accur/10.d0
      probtot = 0.d0
      fac1 = 1.d0
      do n1 = 0,n1max
       if(n1.ne.0)fac1 = fac1*dble(n1)      
       prob1 = eta2(1)**n1/fac1*ex1
       if(n1.gt.eta2(1).and.prob1.le.acc)goto 11
       shift1 =  eben(1)*dble(n1)
       fac2 = 1.d0
       do n2=0,n2max
        if(n2.ne.0)fac2 = fac2*dble(n2)
        prob2 = eta2(2)**n2/fac2*ex2
        prob12 = prob1*prob2
        if(n2.gt.eta2(2).and.prob12.le.acc)goto 12
        shift2 = eben(2)*dble(n2)
        fac3 = 1.d0
        do n3=0,n3max
         if(n3.ne.0)fac3 = fac3*dble(n3)
         prob3 = eta2(3)**n3/fac3*ex3
         prob123 = prob12*prob3
         if(n3.gt.eta2(3).and.prob123.le.acc)goto 13
         shift3 = eben(3)*dble(n3)
         fac4 = 1.d0
         do n4=0,n4max
          if(n4.ne.0)fac4 = fac4*dble(n4)
          probn = prob123*eta2(4)**n4/fac4*ex4
          probtot = probtot+probn
          if(n4.gt.eta2(4).and.probn.le.acc)goto 14
          shift4 = eben(4)*dble(n4)
          shift = shift1+shift2+shift3+shift4
          anshift = shift/delx
          nshift = int(anshift)
          del = anshift-dble(nshift)
          cdel = 1.d0-del       
c      write(*,'(a,4i3,i6,f9.6)')
c     &' n1,n2,n3,n4,nshift,probn = ',n1,n2,n3,n4,nshift,probn
          ppein(1+nshift) = ppein(1+nshift)+probn*cdel*ppinr(1)
          iimax = min(nma+nshift,1000)
          do ii = 2+nshift,iimax
           ppein(ii) = ppein(ii)+
     &probn*(cdel*ppinr(ii-nshift)+del*ppinr(ii-nshift-1))
          enddo
          if(nma+nshift.gt.1000.and.ppein(1000).gt.acc)then
           write(*,*)' range for ppein in distriin_s too small'
           iwarn = 1
          endif
19        continue
         enddo
         if(n4.gt.n4max)write(*,*)' n4.gt.n4max in distriin'
14       continue  
        enddo
        if(n3.gt.n3max)write(*,*)' n3.gt.n3max in distriin'
13      continue  
       enddo
       if(n2.gt.n2max)write(*,*)' n2.gt.n2max in distriin'
12     continue  
      enddo
      if(n1.gt.n1max)write(*,*)' n1.gt.n1max in distriin'
11    continue
      do ii = nma+1,1000
       if(ppein(ii).lt.acc)goto 18
       excit(ii) = excit(nma)+dble(ii-nma)*delx
      enddo
18    nre = ii-1
      ppeintot = 0.d0
      ppinrtot = 0.d0
      if(iout.eq.1)then
       write(16,*)' eloss (MeV),Probability(/MeV)'
      endif
      do ii = 1,nre
       ppeintot = ppeintot + ppein(ii)*delx
       ppinrtot = ppinrtot + ppinr(ii)*delx
       if(iout.eq.1)then
        write(16,'(f13.3,f13.8)')excit(ii),ppein(ii)
       endif
      enddo
c  Now insert result in eloss(mm) table
      iimi = 1
      do mm=mmmin,ienergy
       do ii = iimi,nre
        if(abs(eloss(mm)-excit(ii)).lt.001d0)then
         ppe(mm)=ppein(ii)
         iimi = ii
         goto 111
        endif
       enddo
111    if(eloss(mm).ge.excit(nre)+.001d0)then
        mmmax = mm - 1
        goto 112
       endif
      enddo       
c      demm = excit(nre)
c      do mm=1,ienergy
c       if(eloss(mm).lt.demm)mmmax = mm+1
c       ppe(mm) = 0.d0
c      enddo
112   if(mm.gt.ienergy)then
       mmmax=ienergy
       write(*,'(a)')
     &' Warning: Attempt to extend table above ienergy in distriin_s'
       write(98,'(a)')
     &' Warning: Attempt to extend table above ienergy in distriin_s'
       iwarn = 1
      endif
      if(abs(ppeintot-1.d0).gt.10.d0*accur)then
       write(*,*)
     &' Warning: |ppeintot - ppinrtot|.gt.accur in distriin_s' 
       write(*,*)'      ppeintot, ppinrtot  = ',ppeintot,ppinrtot 
       write(98,*)
     &' Warning: |ppeintot - ppinrtot|.gt.accur in distriin_s' 
       write(98,*)'      ppeintot, ppinrtot  = ',ppeintot,ppinrtot 
       iwarn = 1
      endif
      ierr = 0
      if(iwarn.eq.1) ierr = 2
2     if(iout.eq.1)close(unit=16)
      return
      end
