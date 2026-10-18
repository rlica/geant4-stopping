c      Subroutine calculates mass table for evaporation 
c      and corrected mean excitation energies
c      In ejectile if ierr = 1 at entry or in recoil if ierr = 2 
c      In the latter case the data are stored in the position 
c       of the corresponding ejectile. 
c---------------------------------------------------------
      subroutine emass_s(emass,emassapp,in00,ip00
     &,inmin,inmax,izmin,izmax,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
      include 'grazing_file.icl'
c  irangen and irangep are defined in grazing.icl
      common/za/za,aa,zaa,aaa
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/emas/delans,delanp,delaps,delapp,eneut,eprot
     &,delaans,delaanp,delaaps,delaapp,corma,cormaa
      dimension emass(irangen,irangep),emassapp(irangen,irangep)
c      character*80 fnix
c
cccccc      data fnix /'/home/nanni/grazing/mtablex.dat'/
c
      open(unit=17,file=mass_nix,status='old')
c
      irec = ierr
      iza = int(za+.0001)
      iaa = int(aa+.0001)
      izaa = int(zaa+.0001)
      iaaa = int(aaa+.0001)
      if(iza.lt.7.or.izaa.lt.7)then
       write(*,*)' The charge of projectile and target must be >6'
       goto 33
      endif
      if(irec.eq.1.and.iza+izmin.lt.6)then
       izmin = 6-iza
       write(*,*)' Warning: izmi has been set to ',izmin,' in emass_s'
       ierrr = 2
      elseif(irec.eq.2.and.izaa-izmax.lt.6)then
       izmax = izaa - 6
       write(*,*)' Warning: izma has been set to ',izmax,' in emass_s'
       ierrr = 2
      endif
      if(irec.eq.1)then
       do n=1,irangen
        do m=1,irangep
         emass(n,m)=0.d0
        enddo
       enddo
       open(unit=17,file=mass_nix,status='old')
       do n=1,10000
        read(17,'(3i5,f10.3)',err=13)iz,in,ia,exth
        if(iz.ge.iza+izmin.and.iz.le.iza+izmax.and.
     &in.ge.iaa-iza+inmin.and.in.le.iaa-iza+inmax) then
         exd=dble(exth)
         emass(in-iaa+iza+in00,iz-iza+ip00)=exd-dble(in)*eneut
     &-dble(iz)*eprot
        endif
        if(iz.gt.iza+izmax)goto 11
       enddo
13     write(*,*)' Error in reading mass table in emass_s'
        goto 33
11     close(17)
       do n=inmin+in00,inmax+in00
        in=n-in00
        if(in.gt.0)then
         deln = -delanp
        else
         deln = delans
        endif
        do m=izmin+ip00,izmax+ip00
         im=m-ip00
         if(im.gt.0)then
          delp = -delapp
         else
          delp = delaps
         endif
         emassapp(n,m)=dble(in)*deln+dble(im)*delp
     &-dble(in+iaa-iza)*eneut-dble(im+iza)*eprot+corma
c      nz = iza + im
c      nn = iaa - iza + in
c      write(*,'(2(a,i3),2(a,f10.3))')' Z = ',nz,' N = ',nn,
c     &' emass = ',emass(n,m),' emassapp = ',emassapp(n,m)
        enddo
       enddo
      elseif(irec.eq.2)then
       do n=1,irangen
        do m=1,irangep
         emass(n,m)=0.d0
        enddo
       enddo
       open(unit=17,file=mass_nix,status='old')
       do n=1,10000
        read(17,'(3i5,9f10.3,6f10.2)',err=14)iz,in,ia,exth
        if(iz.ge.izaa-izmax.and.iz.le.izaa-izmin.and.
     &in.ge.iaaa-izaa-inmax.and.in.le.iaaa-izaa-inmin) then
         exd=dble(exth)
         emass(-(in-iaaa+izaa)+in00,-(iz-izaa)+ip00)=exd-dble(in)*eneut
     &-dble(iz)*eprot
        endif
        if(iz.gt.izaa-izmin)goto 12
       enddo
14     write(*,*)' Error in reading mass table in emass_s'
        goto 33
12     close(17)
       do n=inmin+in00,inmax+in00
        in=n-in00
        if(in.gt.0)then
         deln = -delaanp
        else
         deln = delaans
        endif
        do m=izmin+ip00,izmax+ip00
         im=m-ip00
         if(im.gt.0)then
          delp = -delaapp
         else
          delp = delaaps
         endif
         emassapp(n,m)=dble(in)*deln+dble(im)*delp
     &-dble(-in+iaaa-izaa)*eneut-dble(-im+izaa)*eprot+cormaa
c      nz = izaa - im
c      nn = iaaa - izaa - in
c      write(*,'(2(a,i3),2(a,f10.3))')' Z = ',nz,' N = ',nn,
c     &' emass = ',emass(n,m),' emassapp = ',emassapp(n,m)
        enddo
       enddo
      endif
      ierr = 0
      if(ierrr.eq.2)ierr = 2
      return
33    ierr = 1
      return
      end




