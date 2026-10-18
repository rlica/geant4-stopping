c      Subroutine for selecting N and Z for calculation of energy
c       distributions of probability, mean l, scattering angle etc.. 
c      At entry ierr decides output (=1 ejectile, 2= target, 0=partition)
c      Files FOR0//istore(5) to FOR0//istore(5)+ilist will be 
c      used for storing the results 'xname(5)'. lout = number of xnames<6
c---------------------------------------------------------------------
      subroutine selectnz_s(lnzmax,istore,xname,lout,accur,ierr)
      implicit real*8(a-h,o-z)
      include 'grazing.icl'
      common/iza/iza,iaa,izaa,iaaa
      common/capt/eb,rbar,i5,intp
      common/tot/in00,ip00,i0nmin,i0nmax,i0zmin,i0zmax
      common/leg4/dsigz(irangep),dsigzlab(irangep),
     &ppnz(irangen,irangep),iselectn(list),iselectz(list),
     &dsig(irangen,irangep),dsiglab(irangen,irangep),
     &inmin,inmax,izmin,izmax,leg4
      common/head/head
      dimension dsigmax(6),istore(5)
      character*64 xfilename1
      character*44 xfilename
      character*50 addtext
      character*16 head
c      character*16 nanni
      character*7 xna,xname(5)
      character*5 intp
      character*3 out
      character*3 strn,strp
      character*5 elem1(list),elem2(list)
      common/element/elem1,elem2
      character*5 version
      common/vers/version
c
      iq = ierr
      if(intp.eq.'auto')goto 420
419   continue
      if(intp.eq.'inter')write(*,'(a/a/a/a)')
     &' You must now specify the values of delta N and delta Z',
     &'  (of the ejectile) for which you want the energy distribution',
     &' In order to decide you might want to see the five largest',
     &'  cross sections, yes or no!'
      read (i5,'(a)') out
      if(out.eq.'no ')goto 420
      if(out.ne.'yes')goto 419
      dsigmax(6)=1.d+8
      do 445 k=5,1,-1
       dsigmax(k) = 0.d0
       do 445,i=inmin,inmax
        do 445,j=izmin,izmax
         n = in00+i
         m = ip00+j
         if(abs(dsig(n,m)).gt.dsigmax(k).
     &and.abs(dsig(n,m)).lt.dsigmax(k+1))then
          dsigmax(k)=abs(dsig(n,m))
          iselectn(k)=i
          iselectz(k)=j
         endif
445   continue
      do k = 5,1,-1
       i=iselectn(k)
       j=iselectz(k)
       n = in00+i
       m = ip00+j
       if(dsigmax(k).gt.accur)then
        write(*,'(a,2i3,a,f8.3,a)')
     &' dsig(',i,j,') = ',dsig(n,m),' mbarn/sterad'
       else
        kmin = k-1
        goto 420
       endif
      enddo
447   write(*,'(a/a/a/a)')
     & ' You might want to see the next five largest cross sections,',
     &' yes or no'
      read(i5,'(a)')out
      if(out.eq.'no ')goto 420
      if(out.ne.'yes')goto 447
      dsigmax(6)=dsigmax(1)
      do 446 k=5,1,-1
       dsigmax(k) = 0.d0
       do 446,i=inmin,inmax
        do 446,j=izmin,izmax
         n = in00+i
         m = ip00+j
         if(abs(dsig(n,m)).gt.dsigmax(k).
     &and.abs(dsig(n,m)).lt.dsigmax(k+1))then
          dsigmax(k)=abs(dsig(n,m))
          iselectn(k)=i
          iselectz(k)=j
         endif
446   continue
      do k = 5,1,-1
       i=iselectn(k)
       j=iselectz(k)
       n = in00+i
       m = ip00+j
       if(dsigmax(k).gt.accur)then
        write(*,'(a,2i3,a,f8.3,a)')
     &' dsig(',i,j,') = ',dsig(n,m),' mbarn/sterad'
       else
        kmin = k-1
        goto 420
       endif
      enddo
      goto 447
420   continue
      if(intp.eq.'inter')write(*,'(a)')
     &' Now, first enter number of cases to be calculated'
      read(i5,*,err=1100)lnzmax
      icheck = 0
      goto 1110
1100  write(*,*)' you must provide a number'
      icheck = icheck + 1
      if(icheck.gt.5)goto 2
      goto 420
1110  addtext = '    ! number of selected partitions'
      if(intp.eq.'inter')write(10,'(i5,a)')lnzmax,addtext
      if(lnzmax.gt.list)then
       write(*,'(a,i3)')' Maximum number is ',list
       goto 420
      endif
      do l=1,lnzmax
421    if(intp.eq.'inter')write(*,'(a,i3,a)')
     &' Enter',l,'. selected delta N and delta Z (iselectn,iselectz)'
       read(i5,*,err=1101)iseln,iselz
       icheck = 0
       goto 1111
1101   write(*,*)' you must provide two numbers'
       icheck = icheck + 1
       if(icheck.gt.5)goto 2
       goto 421
1111   addtext = '    !  selected partition: delta_N  delta_Z'
       if(intp.eq.'inter')write(10,'(2i4,a)')iseln,iselz,addtext
       iselectn(l)=iseln
       iselectz(l)=iselz
       do lo = 1,lout
        xna = xname(lo)
        ll = istore(lo) + l
        if(version.eq.'vms  '.or.version.eq.'unix* ')then
         call convert_s(iseln,strn,ierr)
         call convert_s(iselz,strp,ierr)
         xfilename1 = 
     &'graz_'//head//'_'//xna//'_'//strn//'_'//strp//'.dat'
         call condens_s(xfilename1,xfilename,le)
         open(unit=ll,file=xfilename,status='unknown')
        else
         write(xfilename,'(a4,i2,a4)')'FOR0',ll,'.DAT'
        endif
        ierr1 = 0
        call element_s(iza+iselz,iaa+iseln+iselz,elem1(l),ierr1)
        call element_s(izaa-iselz,iaaa-iseln-iselz,elem2(l),ierr1)
        if(iq.eq.1)then
         write(*,'(a,a7,a,i3,a,i3,/a,a5,a,a40)')
     &'  ',xna,' distribution for delta N = ',iselectn(l),
     &' and delta Z = ',iselectz(l),
     &'  i.e. for ',elem1(l),' is stored in ',xfilename
         write(15,'(a7,i3,i3,a5,a40)')
     &,xna,iselectn(l),iselectz(l),elem1(l),xfilename
        elseif(iq.eq.2)then
         write(*,'(a,a7,a,i3,a,i3,/a,a5,a,a40)')
     &'  ',xna,' distribution for delta N = ',iselectn(l),
     &' and delta Z = ',iselectz(l),
     &'  i.e. for ',elem2(l),' is stored in ',xfilename
         write(15,'(a7,i3,i3,a5,a40)')
     &,xna,iselectn(l),iselectz(l),elem2(l),xfilename
        else
         write(*,'(a,a7,a,i3,a,i3,/a,a5,a,a5,/a,a40)')
     &'  ',xna,' distribution for delta N = ',iselectn(l),
     &' and delta Z = ',iselectz(l),
     &'  i.e. for the partition ',elem1(l),' + ',elem2(l),
     &'  is stored in ',xfilename
         write(15,'(a7,i3,i3,a5,a5,a40)')
     &,xna,iselectn(l),iselectz(l),elem1(l),elem2(l),xfilename
c  15 is a scratch file to be read later 
        endif
       enddo
      enddo
      ierr = 0
      return
2     ierr = 1
      return
      end
