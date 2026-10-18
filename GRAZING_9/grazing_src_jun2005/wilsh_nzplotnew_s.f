c      This program constructs a command file for producing a
c        Wilczynski plot for given partition isn,isz
c        using the interactive input program [nanni]graph for DISPLA
c      At entry nfile is the number (two digits) of the file where the
c        table with data is found
c      At entry ierr = 0 or 2 indicates cm system ierr = 1 or 3 laboratory.
c        ierr = 0 or 1 in mb/sterad  ierr = 2 or 3 in mb/rad
c-----------------------------------------------------------------
      subroutine wilsh_nzplot_s(isn,isz,nfile,ierr)
      implicit real*8(a-h,o-z)
      common/za/za,aa,zaa,aaa
      common/ecm/ecm
      common/char/elemp,elemt
      common/iee/iealf1
      common/capt/eb,rbar,i5,intp
      character*1 letter
      character*3 strn,strp
      character*5 elemp,elemt
      character*5 intp
      character*8 proj,targ
      character*11 ddd
      character*18 heading
      character*30 addtext
      character*30 xfilename,filename
      character*50 heading1
c
      call date_s(ddd)
      la = ierr
      if(mod(la,2).eq.0)lab = 0
      if(mod(la,2).eq.1)lab = 1
      if(la-2.ge.0)modd = 1
      if(la-2.lt.0)modd = 0
      call convert_s(isn,strn,ierr)
      call convert_s(isz,strp,ierr)
      if(lab.eq.0)then
       xfilename = 'plot_wilsh_cm_'//strn//'_'//strp//'.com'
       call condens_s(xfilename,filename,le)
      elseif(lab.eq.1)then
       xfilename = 'plot_wilsh_lab_'//strn//'_'//strp//'.com'
       call condens_s(xfilename,filename,le)
      endif
      open(unit=11,file=filename,status='unknown')
      write(11,'(a)')'$ set def [WINTHER]'
      write(11,'(a)')'$ run SCR$DISK:[SCRNBINUC.NANNI]graphn'
      write(11,'(a)')'x'
      write(11,'(a)')'18.   23.   14.   14.   3.   3.'
      write(11,'(a)')'0,0'
c The heading is projectile + target 
      proj(1:2) = '^{'
      iflag = 0
      n=2
      do k=1,5
       n=n+1
       letter = elemp(k:k)
       if(letter.eq.'0'.or.letter.eq.'1'.or.letter.eq.'2'
     &.or.letter.eq.'3'.or.letter.eq.'4'.or.letter.eq.'5'
     &.or.letter.eq.'6'.or.letter.eq.'7'.or.letter.eq.'8'
     &.or.letter.eq.'9')then
        proj(n:n) = letter
       elseif(letter.eq.' ')then
        goto 11
       else
        if(iflag.eq.0)then
         proj(n:n) = '}'
         n=n+1
         proj(n:n) = letter
         iflag = 1
        elseif(iflag.eq.1)then
         proj(n:n) = letter
        endif
       endif
      enddo
11    targ(1:2) = '^{'
      iflag = 0
      n=2
      do k=1,5
       n=n+1
       letter = elemt(k:k)
       if(letter.eq.'0'.or.letter.eq.'1'.or.letter.eq.'2'
     &.or.letter.eq.'3'.or.letter.eq.'4'.or.letter.eq.'5'
     &.or.letter.eq.'6'.or.letter.eq.'7'.or.letter.eq.'8'
     &.or.letter.eq.'9')then
        targ(n:n) = letter
       elseif(letter.eq.' ')then
        goto 3
       else
        if(iflag.eq.0)then
         targ(n:n) = '}'
         n=n+1
         targ(n:n) = letter
         iflag = 1
        elseif(iflag.eq.1)then
         targ(n:n) = letter
        endif
       endif
      enddo
3     heading1 = proj//' + '//targ
      call condens_s(heading1,heading,ll)
      if(iealf1.eq.1)then
       ie = int(ecm)
       write(11,'(a,a8,i4,a9)')heading,'at E =',ie,'MeV(cm)'
      elseif(iealf1.eq.2)then
       elab=ecm*(1.d0+aa/aaa)
       ie = int(elab)
       write(11,'(a,a8,i4,a10)')heading,'at E =',ie,'MeV(lab)'
      elseif(iealf1.eq.3)then
       elab=ecm*(1.d0+aa/aaa)
       epslab = elab/aa
       write(11,'(a,a8,f5.2,a6)')heading,'at E =',epslab,'MeV/u'
      endif
      write(11,'(a,i3,a,i3)')
     &' for  \DeltaN =',isn,' and  \DeltaZ =',isz
      read(nfile,'(6e12.4)')tdgmin,tdgmax,deltdg,tkemin,tkemax,dele
      rewind(unit=nfile)
      write(11,'(f5.0,a1,f5.0,a7)')tdgmin,',',tdgmax,',10.,5.'
      if(lab.eq.0)then
       write(11,'(a23)')'\theta _{cm} (degrees)'
      elseif(lab.eq.1)then
       write(11,'(a24)')'\theta _{lab} (degrees)'
      endif
      rewind(unit=nfile)
      write(11,'(f6.0,a1,f6.0,a8)')tkemin,',',tkemax,',20.,10.'
      write(11,'(a9)')'TKE (MeV)'
      write(11,'(a3)')'w,k'
      write(11,'(a4,i2,a4)')'for0',nfile,'.dat'
      write(11,'(21(a/))')
     &'--,5000.,5000.,1.','--,2000.,2000.,1.','.-,1000.,1000.,1.',
     &'--,500.,500.,1.','--,200.,200.,1.','--,100.,100.,1.',
     &'--,50.,50.,1.','--,20.,20.,1.','--,10.,10.,1.','--,5.,5.,1.'
     &,'--,2.,2.,1.',' -,1.,1.,1.','--,0.5,0.5,1.','--,0.2,0.2,1.',
     &'--,0.1,0.1,1.','--,0.05,0.05,1.','--,0.02,0.02,1.',
     &'--,0.01,0.01,1.','--,0.005,0.005,1.','--,0.002,0.002,1.',
     &'..,0.001,0.001,1.'
      write(11,'(a)')'end'
      write(11,'(a)')'w,T .5'
      if(modd.eq.0)write(11,'(a)')
     &'d^{2}\sigma /dE /d\Omega in mb/MeV/sr'
      if(modd.eq.1)write(11,'(a)')
     &'d^{2}\sigma /dE /d\theta in mb/MeV/rad'
      write(11,'(a)')'4.0,17.5'
      write(11,'(a)')'w,T .2'
      write(11,'(a11)')ddd
      write(11,'(a)')'11.0,-1.3'
      write(11,'(2(a/))')'end','$ stop'
2     close(11)
      write(*,'(a/a,a)')
     &' You may make a Wilczynski plot by running',
     &'     the command file ',filename
      return
      end
