c      Subroutine calculates evaporation from mass tables (emasss 
c      and emasssapp) and gammadistribution.
c      In ejectile if ierr = 1 at entry
c      or in recoil if ierr = 2 at entry - In both cases the 
c      data are stored in the position of the corresponding ejectile. 
c      result is in evapdis(n,m,i) with E* = dele*(i-1)
c      Printout according to iout
c---------------------------------------------------------
      subroutine evap_s(amasss,amasssapp,aamasss,aamasssapp,in00,ip00
     &,inmin,inmax,izmin,izmax,pnzevap,dele,evapdis,iout,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
c  irangen and irangep are defined in grazing.icl
      common/za/za,aa,zaa,aaa
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      common/ex/ep0n,ep0p,p0tr,p0a,p0aa,p0in,p0
      common/inel/p0ina,p0inaa,eavina,eavinaa,wid2ina,wid2inaa,
     &avmina,avminaa,eava00,eavaa00,wid2a00,wid2aa00,inelast
      common/nzdis/ p0nn(irangen),p0zz(irangep),
     &eavdist(irangen,irangep),eavdist00,avmadist(irangen,irangep),
     &avmadist00,avmdist(irangen,irangep),avmdist00,
     &wid2dist(irangen,irangep),wid2dist00,wid2adist(irangen,irangep)
      dimension amasss(irangen,irangep),amasssapp(irangen,irangep)
     &,aamasss(irangen,irangep),aamasssapp(irangen,irangep)
     &,pnzevap(irangen,irangep),evapdis(irangen,irangep,ienergy)
     &,thrn(irangen),naa(irangen)
c
      iza = int(za+.0001)
      iaa = int(aa+.0001)
      izaa = int(zaa+.0001)
      iaaa = int(aaa+.0001)
      ievap = ierr
      pnzevaptot = 0.d0
      do m=izmin-2+ip00,izmax+2+ip00
       do n=inmin-2+in00,inmax+2+in00
        pnzevap(n,m)=0.d0
       enddo
      enddo
      if(ievap.eq.1)then
       delp = p0ina
      else
       delp = p0inaa
      endif
      if(ievap.eq.1)then
       do m=izmax+ip00,izmin+ip00,-1
        im = m-ip00
        pnzevapt = 0.d0
        nnmax = int(dble(inmax-inmin+1)/10.d0+.001d0)
        do nn=0,nnmax
         nmax = min(10,inmax-inmin+1-10*nn)
         nmi = inmin+in00+nn*10
         iemax = 0
         lout = 1
      print*,' nnmax,nmax,nmi = ',nnmax,nmax,nmi 
         if(nmax.eq.0)goto 75
         do n=nmi+nmax-1,nmi,-1
          in = n-in00
          emn = amasss(n,m)
          emn1 = amasss(n-1,m)
          thrneu = emn1-emn
          if(emn.eq.0.d0)then
           nnn = iaa-iza+in
           nz = iza+im
           write(*,'(2(a,i3))')' Exact mass not available for N = '
     &,nnn,' Z = ',nz
           thrneu = amasssapp(n-1,m)-amasssapp(n,m)
          endif
          if(emn1.eq.0.d0)then
           nnn = iaa-iza+in-1
           nz = iza+im
           write(*,'(2(a,i3))')' Exact mass not available for N = '
     &,nnn,' Z = ',nz
           thrneu = amasssapp(n-1,m)-amasssapp(n,m)
          endif
          thrn(n) = thrneu
          j = int(thrneu/dele+1.)
          if(j.ge.0)then
           if(j.gt.ienergy)j=ienergy
           do i= 1,ienergy-j
            evapdis(n-1,m,i)=evapdis(n-1,m,i)+evapdis(n,m,i+j)
            evapdis(n,m,i+j)=0.d0
           enddo
           do i = 1,ienergy
            pnzevap(n,m) = pnzevap(n,m)+dele*evapdis(n,m,i)
           enddo
          else
           do i= 1,ienergy+j
            evapdis(n-1,m,i-j)=evapdis(n-1,m,i-j)+evapdis(n,m,i)
            evapdis(n,m,i)=0.d0
           enddo
           pnzevap(n,m) = 0.d0
          endif
          do i=1,ienergy
           if(evapdis(n,m,i).gt.0.d0)then
            if(i.gt.iemax)iemax = i
           endif
          enddo
          if(iout.eq.2)then
           nzz = iza+im
           naa(n) = iaa+im+in
          endif
          pnzevapt = pnzevapt + pnzevap(n,m)
         enddo
         if(iout.eq.2.and.lout.eq.1.and.iemax.gt.0)then
          nrec = nmax
1000      format (/,a5,i3,a5,7f11.2)
1001      format(/,a8,a5,7i11)
          write(34,1000)
     &' Z = ',nzz,' t = ',(thrn(n),n=nmi+nmax-1,nmi,-1)
          write(34,1001)
     &'  E* MeV',' A = ',(naa(n),n=nmi+nmax-1,nmi,-1)
          do i=1,iemax
           ee = dele*dble(i-1)
           write(34,'(f13.2,10(1pe11.3))')
     &ee,(evapdis(n,m,i),n=nmi+nmax-1,nmi,-1)
          enddo
         endif
c  here ends N
        enddo
75      continue
        do i = 1,ienergy
         pnzevap(inmin-2+in00,m) = pnzevap(n,m)
     &+dele*evapdis(inmin-2+in00,m,i)
        enddo
        pnzevapt = pnzevapt + pnzevap(inmin-2+in00,m)
        nzz = iza+im
        print*,' pnzevapt = ',pnzevapt,' for Z = ',nzz
        pnzevaptot = pnzevaptot + pnzevapt
c  here ends Z
       enddo
       pnzevaptot = pnzevaptot + p0tr*p0ina
c
      elseif(ievap.eq.2)then
       do m=izmin+ip00,izmax+ip00
        im = m-ip00
        pnzevapt = 0.d0
        nnmax = int(dble(inmax-inmin+1)/10.d0+.001d0)
        do nn=0,nnmax
         nmax = min(10,inmax-inmin+1-10*nn)
         nmi = inmin+in00+nn*10
         iemax = 0
         lout = 1
         do n=nmi,nmi+nmax-1
          in = n-in00
          emn = aamasss(n,m)
          emn1 = aamasss(n+1,m)
          thrneu = emn1-emn
          if(emn.eq.0.d0)then
           nnn = iaaa-izaa-in
           nz = izaa-im
           write(*,'(2(a,i3))')' Exact mass not available for N = '
     &,nnn,' Z = ',nz
           thrneu = aamasssapp(n+1,m)-aamasssapp(n,m)
          endif
          if(emn1.eq.0.d0)then
           nnn = iaaa-izaa-in-1
           nz = izaa-im
           write(*,'(2(a,i3))')' Exact mass not available for N = '
     &,nnn,' Z = ',nz
           thrneu = aamasssapp(n+1,m)-aamasssapp(n,m)
          endif
          j = int(thrneu/dele+1.)
          thrn(n) = thrneu
          if(j.ge.0)then
           if(j.gt.ienergy)j=ienergy
           do i= 1,ienergy-j
            evapdis(n+1,m,i)=evapdis(n+1,m,i)+evapdis(n,m,i+j)
            evapdis(n,m,i+j)=0.d0
           enddo
           do i = 1,ienergy
            pnzevap(n,m) = pnzevap(n,m)+dele*evapdis(n,m,i)
           enddo
          else
           do i= 1,ienergy+j
            evapdis(n+1,m,i-j)=evapdis(n+1,m,i-j)+evapdis(n,m,i)
            evapdis(n,m,i)=0.d0
           enddo
           pnzevap(n,m) = 0.d0
          endif
          do i=1,ienergy
           if(evapdis(n,m,i).gt.0.d0)then
            if(i.gt.iemax)iemax = i
           endif
          enddo
          if(iout.eq.2)then
           nzz = izaa-im
           naa(n) = iaaa-im-in
          endif
          pnzevapt = pnzevapt + pnzevap(n,m)
         enddo
         if(iout.eq.2.and.lout.eq.1.and.iemax.gt.0)then
          nrec = nmax
1002      format (/,a5,i3,a5,7f11.2)
1003      format(/,a8,a5,7i11)
          write(34,1002)
     &' Z = ',nzz,' t = ',(thrn(n),n=nmi,nmi+nmax-1)
          write(34,1003)
     &'  E* MeV',' A = ',(naa(n),n=nmi,nmi+nmax-1)
          do i=1,iemax
           ee = dele*dble(i-1)
           write(34,'(f13.2,10(1pe11.3))')
     &ee,(evapdis(n,m,i),n=nmi,nmi+nmax-1)
          enddo
         endif
c  here ends N
        enddo
        do i = 1,ienergy
         pnzevap(inmax+2+in00,m) = pnzevap(n,m)
     &+dele*evapdis(inmax+2+in00,m,i)
        enddo
        pnzevapt = pnzevapt + pnzevap(inmax+2+in00,m)
        pnzevaptot = pnzevaptot + pnzevapt
       enddo
       pnzevaptot = pnzevaptot + p0tr*p0inaa
      endif
c
      ierr = 0
      if(abs(pnzevaptot-1.d0).gt.20.d0*accur)then
       write(*,*)' Warning: pnzevaptot = ',pnzevaptot 
       ierr = 2
      endif 
      return
      end
