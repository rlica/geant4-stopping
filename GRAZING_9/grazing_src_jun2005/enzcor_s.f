c      Subroutine calculates corrections to energy distributions
c      from mass tables (amass and aamassapp) and gammadistribution.
c      In ejectile if ierr = 1 at entry or in recoil if ierr = 2 
c      In both cases the data are stored in the position 
c        of the corresponding ejectile. 
c      Result is in evapdis(n,m,i) with E* = dele*(i-1)
c      Printout according to iout
c      If at entry dele = 0, dele is calculated
c---------------------------------------------------------
      subroutine enzcor_s(amass,aamass,amassapp,aamassapp,in00,ip00
     &,inmin,inmax,izmin,izmax,pnzcor,dele,evapdis,iout,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
c  irangen, irangep and ienergy are defined in grazing.icl
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
      dimension amass(irangen,irangep),amassapp(irangen,irangep)
     &,aamass(irangen,irangep),aamassapp(irangen,irangep)
     &,gammd(100),evapdis(irangen,irangep,ienergy)
     &,dee(irangen,irangep),pnzcor(irangen,irangep)
     &,ppc(irangen,irangep),naa(irangen)
c
      iza = int(za+.0001)
      iaa = int(aa+.0001)
      izaa = int(zaa+.0001)
      iaaa = int(aaa+.0001)
      ievap = ierr
      if(ievap.eq.1)then
       delp = p0a
      elseif(ievap.eq.2)then
       delp = p0aa
      else
       write(*,*)' wrong input in enzcor_s, ievap.ne.1 or 2'
       goto 33
      endif
      pnzcortot = 0.d0
      pnzcorttot = 0.d0
      eavcmax = 0.d0
      do n=1,irangen
       do m=1,irangep
        pnzcor(n,m) = 0.d0
        do i=1,ienergy
         evapdis(n,m,i) = 0.d0
        enddo
       enddo
      enddo
      num = (inmax-inmin)*(izmax-izmin)
      nnum = max(num,1)
      pnzcomp = 10.d0*accur/dble(nnum)
c      print*,' pnzcomp = ',pnzcomp 
      do n=inmin+in00,in00+inmax
       in=n-in00
       do m=izmin+ip00,ip00+izmax
        im=m-ip00
        pnzcor(n,m)=p0nn(n)*p0zz(m)
        eavda = .5d0*eavdist(n,m)+eavina
        sig2da = wid2adist(n,m)+wid2ina
        eavdaa = .5d0*eavdist(n,m)+eavinaa
        sig2daa = wid2adist(n,m)+wid2inaa
        if(pnzcor(n,m).gt.pnzcomp)then
c      write(*,'(a,2i3,a,f10.7)')' pnzcor(',in,im,') = ',pnzcor(n,m)
         if(ievap.eq.1)then
          eavd = eavda
          sig2d = sig2da
         else 
          eavd = eavdaa
          sig2d = sig2daa
         endif
         if(in.eq.0.and.im.eq.0)then
          eavda = eava00
          sig2da = wid2a00
          eavdaa = eavaa00
          sig2daa = wid2aa00
          if(ievap.eq.1)then
           eavd = eava00
           sig2d = wid2a00
           pnzcor(n,m) = pnzcor(n,m)-p0a
          elseif(ievap.eq.2)then 
           eavd = eavaa00
           sig2d = wid2aa00
           pnzcor(n,m) = pnzcor(n,m)-p0aa
          endif
          if(eavd.lt.0.d0)eavd = abs(eavd)
          if(sig2d.lt.0.d0)sig2d=eavd
          if(pnzcor(n,m).lt.0.d0)pnzcor(n,m)=abs(pnzcor(n,m))
         endif
c
         dea = amass(n,m)-amassapp(n,m)
         ierr1 = 0
c      print*,' i,j,eavd,sig2d = ',in,im,eavd,sig2d 
         call prob_s(eavda,sig2da,dea,prob,ierr1)
         if(ierr1.eq.1)then
           write(*,*)' For n,m = ',n,m,' in enzcor_s'
          goto 33
         endif
         decra = prob
c =probability that ejectile (in,im) is excited above its ground state
c
         deaa = aamass(n,m)-aamassapp(n,m)
         ierr1 = 0
         call prob_s(eavdaa,sig2daa,deaa,prob,ierr1)
         if(ierr1.eq.1)then
          write(*,*)' For n,m = ',n,m,' in enzcor_s'
          goto 33
         endif
         decraa = prob
c =probability that recoil(in,im) is excited above its ground state
         decr = (decra+decraa)
         if(decr.lt.1.d-10)then
          ppc(n,m) = 0.d0
         else
          ppc(n,m) = 2.d0*decra*decraa/decr
         endif
         if(ievap.eq.1)dee(n,m) = dea
         if(ievap.eq.2)dee(n,m) = deaa
c      write(*,'(a,2i3,a,f10.7)')' dee(',in,im,') = ',dee(n,m)
         eavc=eavd-dee(n,m)
         if(pnzcor(n,m).gt.pnzcomp.and.eavc.gt.eavcmax)then
          eavcmax = eavc
          wwid = sig2d
         endif
        endif
       enddo
      enddo
      do n=inmin+in00,in00-1
       in=n-in00
       do m=izmin+ip00,ip00-1
        im=m-ip00
        if(pnzcor(n,m).gt.pnzcomp)then
         pcorn = ppc(n+1,m)
         pcorp = ppc(n,m+1)
         pcor = pcorn+pcorp
         if(pcor.ne.0.d0)then
          pnzcor(n+1,m)=pcorn/pcor*(1.d0-ppc(n,m))*pnzcor(n,m)
     &+pnzcor(n+1,m)
          pnzcor(n,m+1)=pcorp/pcor*(1.d0-ppc(n,m))*pnzcor(n,m)
     &+pnzcor(n,m+1)
         endif
         pnzcor(n,m)=pnzcor(n,m)*ppc(n,m)
        endif
       enddo
       do m=ip00+izmax,ip00+1,-1
        im=m-ip00
        if(pnzcor(n,m).gt.pnzcomp)then
         pcorn = ppc(n+1,m)
         pcorp = ppc(n,m-1)
         pcor = pcorn+pcorp
         if(pcor.ne.0.d0)then
          pnzcor(n+1,m)=pcorn/pcor*(1.d0-ppc(n,m))*pnzcor(n,m)
     &+pnzcor(n+1,m)
          pnzcor(n,m-1)=pcorp/pcor*(1.d0-ppc(n,m))*pnzcor(n,m)
     &+pnzcor(n,m-1)
         endif
         pnzcor(n,m)=pnzcor(n,m)*ppc(n,m)
        endif
       enddo
      enddo
      do n=in00+inmax,in00+1,-1
       in=n-in00
       do m=izmin+ip00,ip00-1
        im=m-ip00
        if(pnzcor(n,m).gt.pnzcomp)then
         pcorn = ppc(n-1,m)
         pcorp = ppc(n,m+1)
         pcor = pcorn+pcorp
         if(pcor.ne.0.d0)then
          pnzcor(n-1,m)=pcorn/pcor*(1.d0-ppc(n,m))*pnzcor(n,m)
     &+pnzcor(n-1,m)
          pnzcor(n,m+1)=pcorp/pcor*(1.d0-ppc(n,m))*pnzcor(n,m)
     &+pnzcor(n,m+1)
          endif
         pnzcor(n,m)=pnzcor(n,m)*ppc(n,m)
        endif
       enddo
       do m=ip00+izmax,ip00+1,-1
        im=m-ip00
        if(pnzcor(n,m).gt.pnzcomp)then
         pcorn = ppc(n-1,m)
         pcorp = ppc(n,m-1)
         pcor = pcorn+pcorp
         if(pcor.ne.0.d0)then
          pnzcor(n-1,m)=pcorn/pcor*(1.d0-ppc(n,m))*pnzcor(n,m)
     &+pnzcor(n-1,m)
          pnzcor(n,m-1)=pcorp/pcor*(1.d0-ppc(n,m))*pnzcor(n,m)
     &+pnzcor(n,m-1)
         endif
         pnzcor(n,m)=pnzcor(n,m)*ppc(n,m)
        endif
       enddo
      enddo
      m = ip00
      im = 0
      do n=inmin+in00,in00-1
       in=n-in00
       pnzcor(n+1,m)=(1.d0-ppc(n,m))*pnzcor(n,m)+pnzcor(n+1,m)
       pnzcor(n,m)=pnzcor(n,m)*ppc(n,m)
      enddo
      do n=in00+inmax,in00+1,-1
       in=n-in00
       pnzcor(n-1,m)=(1.d0-ppc(n,m))*pnzcor(n,m)+pnzcor(n-1,m)
       pnzcor(n,m)=pnzcor(n,m)*ppc(n,m)
      enddo
      n = in00
      in = 0
      do m=izmin+ip00,ip00-1
       im=m-ip00
       pnzcor(n,m+1)=(1.d0-ppc(n,m))*pnzcor(n,m)+pnzcor(n,m+1)
       pnzcor(n,m)=pnzcor(n,m)*ppc(n,m)
      enddo
      do m=ip00+izmax,ip00+1,-1
       im=m-ip00
       pnzcor(n,m-1)=(1.d0-ppc(n,m))*pnzcor(n,m)+pnzcor(n,m-1)
       pnzcor(n,m)=pnzcor(n,m)*ppc(n,m)
      enddo
c      print*,' eavcmax,wwid = ',eavcmax,wwid 
      if(dele.lt.1.d-10)then
       dele = (eavcmax+3.d0*sqrt(wwid))/dble(ienergy)*2.d0
      endif
      do m=izmin+ip00,ip00+izmax
       im=m-ip00
       nnmax = int(dble(inmax-inmin+1)/10.d0+.001d0)
       do nn=0,nnmax
        nmax = min(10,inmax-inmin+1-10*nn)
        nmi = inmin+in00+nn*10
        lout = 0
        iimax = 0
        do n=nmi,nmi+nmax-1
         in=n-in00
         if(pnzcor(n,m).gt.pnzcomp)then
          if(ievap.eq.1)then
           eavd = .5d0*eavdist(n,m)+eavina
           sig2d = wid2adist(n,m)+wid2ina
          else 
           eavd = .5d0*eavdist(n,m)+eavinaa
           sig2d = wid2adist(n,m)+wid2inaa
          endif
          if(in.eq.0.and.im.eq.0)then
           if(ievap.eq.1)then
            eavd = eava00
            sig2d = wid2a00
           else 
            eavd = eavaa00
            sig2d = wid2aa00
           endif
           if(eavd.lt.0.d0)eavd = abs(eavd)
           if(sig2d.lt.0.d0)sig2d=eavd
           if(pnzcor(n,m).lt.0.d0)pnzcor(n,m)=abs(pnzcor(n,m))
          endif
          ierr1 = 0
          e0 = 0.d0
          call probdtab_s(eavd,sig2d,e0,dele,gammd,accur,ierr1)
          if(ierr1.eq.1)then
           write(*,*)' For in,im = ',in,im,' in enzcor_s'
           do iii = 1,100
            gammd(iii) = 0.d0
           enddo
          elseif(ierr1.eq.3)then
           write(*,*)' Warning reduced accuracy for in,im = ',in,im,
     &' in enzcor_s'
          endif 
          deee = dee(n,m)
          if(deee.ge.dele)then
           imin = int(deee/dele)
           imax = min(imin+ienergy-1,100)
           do i = imin,imax
            evapdis(n,m,i-imin+1) = pnzcor(n,m)*gammd(i)
           enddo
          else
           if(deee.gt.0.d0)then
            imin = 0
           else
            imin = -int(deee/dele)+1
           endif           
           if(imin.ge.ienergy)imin = ienergy-1
           do i = 1,imin+1
            evapdis(n,m,i)= 0.d0
           enddo
           do i = imin+1,ienergy-1
            evapdis(n,m,i+1) = pnzcor(n,m)*gammd(i-imin)
           enddo
           if(imin.lt.ienergy)then
            if(ierr1.eq.2)then
             ierr2 = 0
             call prob_s(eavd,sig2d,dele,probe,ierr2)
             if(ierr2.eq.1)then
              write(*,*)' For n,m = ',n,m,' in enzcor_s'
              goto 33
             endif
             gam0 =.5d0*(1.d0-probe)/dele
             evapdis(n,m,imin+1) = gam0*pnzcor(n,m)
            else
             evapdis(n,m,imin+1) = 0.d0
            endif
           endif
          endif
          if(iout.eq.1)then
           if(ievap.eq.1)then
            nzz = iza+im
            naa(n) = iaa+im+in
           else
            nzz = izaa-im
            naa(n) = iaaa-im-in
           endif
          endif
          pnzcort = 0.d0
c       pnzcort = .5d0*dele*evapdis(n,m,1)
          do i=1,ienergy
           ev =evapdis(n,m,i)
           if(ev.gt.0.d0.and.i.gt.iimax)iimax = i
           pnzcort = pnzcort + dele*ev
          enddo
          pnzcortot= pnzcortot + pnzcor(n,m)
          pnzcorttot = pnzcorttot + pnzcort
          lout = 1
         endif
        enddo
        if(iout.eq.1.and.lout.eq.1)then
         nrec = nmax
1000     format (/,a5,i3,a5,7f11.2)
1001     format (a8,a5,7i11)
         nnn  = nmi+nmax-1
         write(33,1000)
     &' Z = ',nzz,' d = ',(dee(n,m),n=nmi,nnn)
         write(33,1001)
     &'  E* MeV','   A =  ',(naa(n),n=nmi,nnn)
         do i=1,iimax
          ee = dele*dble(i-1)
          write(33,'(f13.2,10(1pe11.3))')
     &ee,(evapdis(n,m,i),n=nmi,nmi+nmax-1)
         enddo
        endif
       enddo
c       pnzcorttot = pnzcorttot + pnzcortt
      enddo
      pnzcortot= pnzcortot + delp
      pnzcorttot= pnzcorttot + delp
      ierr = 0
      if(abs(pnzcortot-1.d0).gt.5.d0*accur)then
       write(*,*)' Warning: pnzcortot = ',pnzcortot
       ierr = 2
      endif
      if(abs(pnzcorttot-1.d0).gt.20.d0*accur)then
       write(*,*)' Warning: pnzcorttot = ',pnzcorttot
       ierr = 2
      endif
      return
33    ierr = 1
      return
      end
