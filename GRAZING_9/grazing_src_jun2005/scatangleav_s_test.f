c      Tests scatangleavx_s subroutine
c------------------------------------------------------
      implicit real*8(a-h,o-z)
      common/za/ za,aa,zaa,aaa
      common/ecm/ ecm
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      accur = .0001d0
      vdr = 0.d0
      print*,' Enter za,aa,zaa,aaa'
      read*,za,aa,zaa,aaa
6     print*,' Enter ecm (exit if ecm.lt.0)'
      read*,ecm
      if(ecm.lt.0.d0)goto 2
      ierr = 0
      call lgraz_s(ecr,alc,eb,alg,rgz,ierr)
      if(ierr.eq.1)goto 2
      print*,' ebar,algraz = ',eb,alg
      ierr = 1
      call rainbow_s(rainbc,rainbn,alrc,alrn,ierr)
      if(ierr.eq.1)goto 2
      rainbdg = 180.d0/pi*rainbc
      print*,' alrain,rainbdg = ',alrc,rainbdg
4     print*,' Enter al (new ecm if al.lt.0)'
      read*,al
      if(al.lt.0.d0)goto 6
      ierr = 2
      call scatangle_s(al,aux,theta,dthet,ierr)
      if(ierr.eq.1)goto 2
      thetdg = 180.d0/pi*theta
      print*,' turning point in entrance channel = ',aux
      print*,' elastic scattering angle = ',thetdg
5     print*,' enter delw,w0 (new al if delw.lt.0)'
      read*,delw0,w00
      if(delw0.lt.0.d0)goto 4
3     print*,' Now dele,delal,dela,delz (new delw if dele.gt.ecm)'
      read*,dele,delal,dela,delz
      if(dele.gt.ecm)goto 5
      ierr = 0
      theta = 0.d0
      call scatangleav_s(al,theta,dthet,dela,delz,dele,
     &delal,thetav,dthetav,rturnex,ierr)
      if(ierr.eq.1)then
       goto 2
      elseif(ierr.eq.0)then
       scatdg = 180.d0/pi*thetav
       write(*,'(2(a,f10.6,/))')
     &' average scatangle between entrance and exit = ',scatdg,
     &' dthet = ',dthet
      elseif(ierr.eq.2)then
       print*,' rturn inside rgraz in entrance or exit channel'
      endif
      ierr = 1
      w0 = w00
      delw = delw0
      call scatangleavx_s(al,ph,dph,rtex,dele,delal,dela,delz
     &,w0,delw,ierr)
      if(ierr.eq.1)goto 2
      print*,'  w0 increased to ',w0
      scatx =  180.d0/pi*(pi-ph)
      print*,' turning point in exit channel = ',rtex
      write(*,'(2(a,f10.6/))')
     &' scattering angle from scatangleavx_s = ',scatx,
     &' dthet = ',dph
      goto 3
2     end
