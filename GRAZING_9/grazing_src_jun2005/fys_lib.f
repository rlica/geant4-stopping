c     FYS_LIB Library contains subroutines
c       for element names, masses, BE2 and BE3 values
c       contains element_s, masse_s and be23c_s subroutines
c       and uses tables of masses (MASSEXP.DAT (experimental) and 
c       MTABLEX.DAT constructed by Nix) and BE values (BE23.DAT)
c--------------------------------------------------------------
c       Element_s subroutine reads for given values of iz and ia the name
c       of the element in the form (ia element) or (element)if ia.le.0
c       At exit ierr = 1 if element is not found
c---------------------------------------------------------------------
c
      subroutine element_s(iz,ia,elem,ierr)
      character elem*5
      character*2 sym(109)
      character*2 symm
      data sym /
     *     'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne'
     1   , 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca'
     2   , 'Sr', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'
     3   , 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr'
     4   , 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn'
     5   , 'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd'
     6   , 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb'
     7   , 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'
     8   , 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th'
     9   , 'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm'
     *   , 'Md', 'No', 'Lr', 'Rf', 'Ha', 'Sg', 'Ns', 'Hs', 'Mt'/
c
      if(ia.lt.iz)then
       write(*,*)' ia less than iz in element_s'
       ierr = 1
       return
      endif      
      if(iz.gt.109)then
       symm = 'Xx'
      elseif(iz.gt.0)then
       symm = sym(iz)
      elseif(iz.lt.0)then
       write(*,*)' iz less than 0 in element_s'
       ierr = 1
       return
      endif      
      if(ia.eq.1.and.iz.eq.0)then
       elem = 'n    '
      elseif(ia.lt.0)then
       write(elem,'(a)') symm
      elseif(ia.le.9)then
       write(elem,'(i1,a)') ia,symm
      else if(ia.le.99) then
       write(elem,'(i2,a)') ia,symm
      else
       write(elem,'(i3,a)') ia,symm
      endif
      ierr=0
      return
      end
c---------------------------------------------------------------
c      Masse_s subroutine reads for given values of iz,ia the
c      mass excess ex measured with respect to 12C and its
c      error dex, both in MeV
c      If at entry ierr = 2 it reads from MASSEXP.DAT (experimental)
c        then at exit ierr = 1 if not found
c      If at entry ierr = 3 it reads from MTABLEX.DAT constructed by Nix
c        then at exit ierr = 1 if not found
c      If at entry ierr = 0 it reads first MASSEXP.DAT and if not available
c        then MTABLEX.DAT. At exit ierr = 2 or 3 respectively
c--------------------------------------------------------------------
      subroutine masse_s(izd,iad,exd,dexd,ierr)
      real*8 exd,dexd
      character*5 ele
c
      include 'grazing_file.icl'
c
      data nmass /257/
c
      if(izd.lt.0 .or. iad.lt.0 .or. iad.lt.izd)then
       write(*,*)' Error: wrong input in masse_s'
       ierr = 1
       return
      endif
      if(izd.eq.0 .and. iad.eq.0) then
       exd=0.d0
       dexd=0.d0
       ierr = 2
       return  
      endif
      if(ierr.eq.3)goto 10
      if(iad.gt.nmass)then
       if(ierr.eq.0)then
        goto 10
       elseif(ierr.eq.2)then
        goto 13
       endif
      endif
20    continue
      open(unit=18,file=mass_exp,status='old')
      do k=1,10000
       read(18,'(2i5,1f12.4)',end=11)iz,ia,exth
       if(ia.eq.iad) then
        if(iz.eq.izd)then
         exd=dble(exth)/1000.d0
         dexd=0.d0
         ierr = 2
         close(18)
         return
        endif
       endif
      enddo
11    if(ierr.eq.0)then
       write(*,*)' Experimental mass could not be found'
       close(unit=18)
       goto 10
      elseif(ierr.eq.2)then
       goto 13
      endif
10    open(unit=17,file=mass_nix,status='old')
      do n=1,20000
c  here is read a condensed version of the nix table:
       read(17,'(3i5,1f10.2)',end=12)iz,in,ia,exth
       if(iz.eq.izd) then
        if(ia.eq.iad)then
         exd=dble(exth)
         dexd=0.d0
         ierr = 3
         close(17)
         return
        endif
       endif
       if(iz.gt.izd)goto 12
      enddo
12    close(17)
      ierr = 0
      call element_s(izd,iad,ele,ierr)
      if(ierr.eq.0)then
       write(*,'(a,a5,a)')
     &' Mass of ',ele,' could not be found in MTABLEX.DAT'
      else
       write(*,'(a,i3,a,i3,a)')
     &' Mass of A,Z = ',iad,',',izd,' could not be found in MTABLEX.DAT'
      endif
      ierr = 1
      return
13    close(18)
      ierr = 0
      call element_s(izd,iad,ele,ierr)
      if(ierr.eq.0)then
       write(*,'(a,a5,a)')
     &' Mass of ',ele,' could not be found in MASSEXP.DAT'
      else
       write(*,'(a,i3,a,i3,a)')
     &' Mass of A,Z = ',iad,',',izd,' could not be found in MASSEXP.DAT'
      endif
      ierr = 1
      return
      end
c------------------------------------------------------------------
c      Subroutine be23c_s finds energy (in MeV) of lowest 2+,
c       belam=B(E2,0-2)(in e^2*b^2),energy (in MeV) of lowest 3-
c       or belam=B(E3,0-3) (in e^2*b^3) for given A, Z, lam=multipole
c       (integers). It furthermore delivers the mass parameter
c       the restoring force and the zeropoint amplitude=alpha.
c       In calculating these quantities it uses the Coulomb radius 
c       parameter r_0c = 1.20 fm.
c      Uses datafile for099, renamed be23.dat.
c      If at entry ierr=0 and data not available: it asks for information
c       and stores for future use. If belam is entered as a negative number
c       it interpolates the value between neighbors and ierr = 3 at exit. 
c       If belam is entered as zero it indicates that the collective mode 
c       is neglected then energy=1,clam=1.d+5,clam=1.d+5,alpha=0
c       will indicate that data are missing and ierr = 2 at exit
c      If be23c_s is run in batch set ierr = 1 at entry:
c       then there will be automatic interpolation (ierr = 3 at exit)
c       and automatic neglect of states with missing information 
c       indicated by belam=0,energy=1,clam=1.d+5,clam=1.d+5,alpha=0
c       and ierr = 2 at exit.
c----------------------------------------------------

      subroutine be23c_s(ize,iae,lam,energy,belam,clam,dlam,alpha,ierr)
      implicit real*8(a-h,o-z)
      parameter (nodat=99)
      character*2 par
      character*3 ok
      character*5 eleme
c
      include 'grazing_file.icl'
c
      open (unit=nodat,file=file_be23,status='old')
      istat = ierr
      iz = ize
      ia = iae
cnan
      if(iz.le.5) goto 6000
cnan
      izodd = mod(ize,2)
      iaodd = mod(iae,2)
      if(izodd.eq.0.and.iaodd.eq.0)then
       interp = 0
       goto 4000
      elseif(izodd.eq.1.and.iaodd.eq.1)then
       iz = ize + 1
       ia = iae + 1
       interp = 1
       goto 4000
      elseif(izodd.eq.0.and.iaodd.eq.1)then
       iz = ize 
       ia = iae + 1
       interp = 1
       goto 4000
      else
       goto 6000      
      endif
4000  rewind nodat
3000  read (nodat,*,end=1000,err=8000) iia,iiz,ilam,Energ,BEla
      if (iia.eq.ia.and.iiz.eq.iz.and.ilam.eq.lam)  goto 2000
      goto 3000
2000  if(interp.eq.0)then
       goto 5000
      elseif(interp.eq.1)then
       enplus = Energ
       belplus = BEla
       if(izodd.eq.1.and.iaodd.eq.1)then
        iz = ize - 1
        ia = iae - 1
       elseif(izodd.eq.0.and.iaodd.eq.1)then
        iz = ize 
        ia = iae - 1
       else
        goto 6000      
       endif
       interp = 2
       goto 4000
      elseif(interp.eq.2)then
       enminus = Energ
       belminus = BEla
       Energ = (enplus+enminus)/2.d0
       BEla = (belplus+belminus)/2.d0
      endif
5000  energy = Energ
      belam = BEla
      if(belam.lt.1.d-7)goto 6000
      belamf = belam*10.**(2*lam)
      fac = 3.d0/4.d0/3.14159d0
      ain = dble(iae)
      zin = dble(ize)
      alpha = sqrt(belamf/dble(2*lam +1))/(fac*zin
     1    *(1.2d0)**lam*ain**(dble(lam)/3.d0))
      beta = sqrt(dble(2*lam + 1))*alpha
      sigma = beta/sqrt(4.d0*3.14159d0)*ain**.333333d0*1.2d0
      clam = energy/2.d0/alpha**2
      dlam = clam/energy**2
      f_lam = energy*belamf/(fac*dble(2*lam + 1)*dble(lam)
     1    *20.8949d0*zin**2/ain*(1.2d0*ain**.33333d0)**(2*lam-2))
      ierr = 0
      if(iaodd.eq.1)ierr = 3
      close(unit=nodat)
      return
1000  if(istat.eq.1)goto 6000
1001  write(*,'(a,i3,a,i4,a,i1,a,i1,/a/a)')' For Z =',iz,' A =',ia,
     &' Enter Energy in MeV , B(E',lam,') in e^2 b^',lam,
     &' If you enter 0 0, data are not entered in the table!',
     &'    and the mode is neglected'
      read(*,*,err=150) Ener,BEl
      goto 250
150   write(*,*)' your entry was not two numbers'
      goto 1001
250   if(BEl.gt.1.d-10)then
       if(lam.eq.2)write(par,'(i1,a)') lam,'+'
       if(lam.eq.3)write(par,'(i1,a)') lam,'-'
       call element_s(iz,ia,eleme,ierr)
       write(*,'(/a/a,a5/a,a2,a,f6.3,a/a,i1,a,e10.3,a,i1/a)')
     &' Before the number is entered into the table',
     &' please check that for the nucleus ',eleme,
     &'  excitation energy of the lowest ',par,' state is ',Ener,' MeV',
     &'  B(E',lam,') value is ',BEl,' e^2 b^',lam,
     &' Is it OK? (yes or no)'
       read(*,'(a)')ok
       if(ok.eq.'yes')then
        write (nodat,'(3i12,f11.6,e19.7)',err=9000)
     & ia,iz,lam,Ener,BEl
       elseif(ok.eq.'no ')then
        goto 1001
       else
        write(*,*)' Try again'
        goto 250
       endif
      else
       goto 6000
      endif
      goto 4000
6000  belam = 0.d0
      energy = 1.d0
      clam = 1.d5
      dlam = 1.d5
      alpha = 0.d0
      ierr = 2
      close(unit=nodat)
      return
8000  print*,' Error by reading be23.dat'
      goto 7000
9000  print*,' Error by writing on be23.dat'
7000  ierr = 1
      close(unit=nodat)
      return
      end








