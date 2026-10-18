c     ---------------------------------------------------------------
c     GZP2 - 
c     program to produce graph from data_file produced by
c     the user. The program is based on the PGPLOT library. 
c
c     Coded by Nanni (Giovanni Pollarolo)               October  1997
c     -  modified  (histogram and greek letter)         January  1998
c     -  modified  (logx,logy scale)                    February 1998
c     -  modified  (legenda,text_cm)                    February 1998
c     -  modified  (linux version)                      September2000
c
c     NB - If you install the program in your directory/machine
c          please do not change the name (if you do not like the
c          name probably you wan't like the program)
c          Se installate questo programma sulle vostre directory e/o
c          macchine, per favore, non cambiatene il nome (se non vi 
c          piace il nome forse non vi piacera' il programma).
c
c     User Manual (only in Italian, sorry):
c          For an HTML version refer to the page 'software' on 
c          the  WEBpage www.ph.unito.it/~nanni and from there surf  
c          to the gzp2 manual
c
c     Contact:
c          For mistakes, comments, explantions, suggestions, ...  send 
c          messages to  nanni@to.infn.it
c     ---------------------------------------------------------------
c
c     dimension and character arrays                      
      character*72  ttlcr    ,sttlcr   ,txcr      ,tycr
      character*200 ttlc     ,sttlc    ,txc       ,tyc
      character*80  command
      character*80  fileinput
      character*8   dev
      character*12  optionx,optiony
      character*1   curve,col,mycl(0:8)
c
      data mycl /'n', 'w', 'r', 'g', 'b', 'c', 'v',
     1           'y', 'm'/
c
      common /device/ dev
      common /runit/ i5
      common /logsc/ ilogx,ilogy
      common /limit/ xmin,xmax,ymin,ymax
      common /dims/  pgx ,pgy ,arx ,ary ,pox ,poy,   scx,  scy
      common /frame/ dx,dy,itickx,iticky
      common /opt/ optionx,optiony
c
c     some defaults (input unit...)
      i5=5
      inext=0
      scx=1.
      scy=1.
      pox=0.
      poy=0.
      isublab=0
c
c     enter filename with your plotting command
c
      print *,' enter inputfile name'
      read(i5,'(a)') fileinput
      if(fileinput.eq.' ') then
         i5=5
        else
         i5=10
         ll=index(fileinput,'   ')
         open(unit=10,file=fileinput(1:ll),status='old')
      endif
c
c     entering device
c
      print *,' Enter device:' 
      print *,' /PS    (PostScript file, landscape orientation)'
      print *,' /VPS   (PostScript file, portrait orientation)'
      print *,' /CPS   (Colour PostScript file, landscape orientation)'
      print *,' /VCPS  (Colour PostScript file, portrait orientation)'
      print *,' /CC    (LJ250 Sixel Colour Printer file .ccplt)'
      print *,' /XWINDOW (X window window@node:display.screen/xw)'
      read(*,'(a)') dev
      call pgbeg(0,dev,1,1)
c
c     entering plot dimensions (page,area-plot and margins)
c
   77 continue
      ips=index(dev,'ps')
      call pgsci(1)
c      if(ips.ne.0)then
c         call pgsci(0)
c      endif      
c
      print *,' Enter:'
      print *,' pagex,pagey,arx,  ary,  pox,  poy [f6.2]'
      print *,' 18.   26.   10.   14.   4.    5.  [default]'
      read(i5,'(a)') command
      read(command,'(6f6.2)') pgxr,pgyr,arxr,aryr,poxr,poyr
      if(pgxr.lt.0.001) pgxr=18.
      if(pgyr.lt.0.001) pgyr=26.
      if(arxr.lt.0.001) arxr=10.
      if(aryr.lt.0.001) aryr=14.
      if(poxr.lt.0.001.and.inext.eq.0) poxr=4.
      if(poyr.lt.0.001.and.inext.eq.0) poyr=5.  
      if(inext.eq.0) then
      if(dev.eq.'/xwindow'.or.dev.eq.'/XWINDOW'
     1      .or.dev.eq.' ') then
             scy=1.2672
             scx=1.44
             sizexw=35.     ! 35.        ! 25.   ! 30.
             sizeyw=23.     ! 24.75      ! 18.33 ! 22.
             aux1=pgxr/sizexw
             aux2=pgyr/sizeyw
             if(aux1.gt.1.or.aux2.gt.1) then
                aux=max(aux1,aux2)
                scx=scx*aux
                scy=scy*aux
             endif
cc          if(scx.gt.1..or.scy.gt.1.) sc=max(scx,scy)
         endif
         pgx=pgxr/scx
         pgy=pgyr/scy
         arx=arxr/scx
         ary=aryr/scy
         pox=poxr/scx
         poy=poyr/scy
         call pgpap(pgx/2.504,pgy/pgx)
        else
         pox=pox+poxr/scx
         poy=poy+poyr/scy
         arx=arxr/scx
         ary=aryr/scy
      endif
c
c     boxit the full plotting area (only for test)
      xr1=0.0
      xr2=0.99
      yr1=0.
      yr2=0.99
      call pgsvp(xr1,xr2,yr1,yr2)
      call pgswin(0.,1.,0.,1.)
c      call pgbox('bc',1.,0,'bc',1.,0)
c
      xr1=pox/pgx
      xr2=(pox+arx)/pgx
      yr1=poy/pgy
      yr2=(poy+ary)/pgy
      call pgsvp(xr1,xr2,yr1,yr2)
c
c     log scale
      print *, ' Enter type of scale fox x/y: (0,1) (linx,logy)'
      read(i5,*) ilogx,ilogy
      if(isublab.eq.0)then
         optionx='bcntsv'
         optiony='bcntsv'
         if(ilogy.eq.1) optiony='bcnltsv2'
         if(ilogx.eq.1) optionx='bcnltsv2'
        elseif(isublab.eq.2)then
         optionx='bcts'
         optiony='bcntsv'
         if(ilogy.eq.1) optiony='bcnltsv2'
         if(ilogx.eq.1) optionx='bclts2'
        elseif(isublab.eq.1)then
         optionx='bcntsv'
         optiony='bcts'
         if(ilogy.eq.1) optiony='bclts2'
         if(ilogx.eq.1) optionx='bcnltsv2'
        elseif(isublab.eq.3)then
         optionx='bcts'
         optiony='bcts'
         if(ilogy.eq.1) optiony='bclts2'
         if(ilogx.eq.1) optionx='bclts2'
      endif
c      
c     entering the plot title 
      print *,' Enter Headline'
      read(i5,'(a)') ttlcr
      call rosetta(ttlcr,ttlc)
      lttl=index(ttlc,'   ')-1
      print *,' Enter Sub-headline'
      read(i5,'(a)') sttlcr
      call rosetta(sttlcr,sttlc)
      lsttl=index(sttlc,'   ')-1
c
c     reading the frame limits (and texts)
      print *,' Enter xmin,xmax,dx'
      read(i5,'(a)') command
      iout=0
      itickx=1
      read(command,*,err=201,end=201) xmin,xmax,dx,itickx
  201 continue
      print *,' Any text for x-axis?'
      read(i5,'(a)') txcr
      call rosetta(txcr,txc)
      ltx=index(txc,'   ')-1
      if(ltx.eq.0) then
         ltx=1
         txc=' '
      endif
      if(ilogx.eq.1) then
         xmin=log10(xmin)
         xmax=log10(xmax)
      endif         
c
      print *,' Enter ymin,ymax,dy'
      read(i5,'(a)')command
      iticky=1
      read(command,*,err=202,end=202) ymin,ymax,dy,iticky
  202 continue
      print *,' Any text for y-axis?'
      read(i5,'(a)') tycr
      call rosetta(tycr,tyc)
      lty=index(tyc,'   ')-1
      if(lty.eq.0) then
         lty=1
         tyc=' '
      endif
      if(ilogy.eq.1) then
         ymin=log10(ymin)
         ymax=log10(ymax)
      endif         
c
c     ploting the frame
      ifont=2
      call pgscf(ifont)
      call pgslw(2)
      call pgsch(1.)
      call pgswin(xmin,xmax,ymin,ymax)
      call pgbox(optionx,dx,itickx,optiony,dy,iticky)
      call pgslw(3)
      call pgsch(1.2)  !(1.2)
c
      disp=3.
      coord=0.5
      fjust=0.5
      call pgmtxt('B',disp,coord,fjust,txc(1:ltx))
c
      disp=3.5
      coord=0.5
      fjust=0.5
      call pgmtxt('L',disp,coord,fjust,tyc(1:lty))
c
      disp=1.9
      coord=0.5
      fjust=0.5
      call pgsch(1.8)
      call pgmtxt('T',disp,coord,fjust,ttlc(1:lttl))
c
      disp=0.9
      coord=0.5
      fjust=0.5
      call pgsch(1.3)
      call pgmtxt('T',disp,coord,fjust,sttlc(1:lsttl))
c      call pgsch(1.)
c
      call pgsch(1.)
      call pgslw(1)
c
c     starting the plot segments
c
 1000 continue
      print *,' Enter :'
      print *,' col,curve,kind,[icol1,icol2]'
  100 continue
      read(i5,'(a)')command
      if(command(1:1).eq.'!'.or.command(1:1).eq.'%') goto 100
      if(command(1:4).eq.'next'.or.command(1:4).eq.'NEXT')then
         inext=1
         isublab=0
         read(command(5:50),*,end=190,err=190) isublab
  190    continue
         goto 77
      endif
      if(command(1:3).eq.'end'.or.command(1:3).eq.'END') go to 1001
c
      read(command,'(a1,1x,a1)') col,curve
c
c     get the ncol from the letter:
      do n=0,8
         if(col.eq.mycl(n)) goto 102
      enddo
  102 continue
      call pgsci(n)
c
c     drawing the curve
      if(curve.eq.'c') then
         call gr_line(command)
        elseif(curve.eq.'C') then
         call gr_line_clip(command)
        elseif(curve.eq.'p') then
         call gr_point(command)
        elseif(curve.eq.'P') then
         call gr_point_sc(command)
        elseif(curve.eq.'e') then
         call gr_errorv(command)
        elseif(curve.eq.'E') then
         call gr_errorvo(command)
        elseif(curve.eq.'h') then
         call gr_hist(command)
        elseif(curve.eq.'H') then
         call gr_hist_center(command)
        elseif(curve.eq.'t') then
         call gr_text(command)
        elseif(curve.eq.'T') then
         call gr_text_cm(command)
        elseif(curve.eq.'k') then
         call gr_cont(command)
        elseif(curve.eq.'m') then
         call gr_mosaic(command)
        elseif(curve.eq.'l') then
         call gr_legenda(command)
        elseif(curve.eq.'s') then
         call gr_shadow(command)
      endif
      goto 1000
c
 1001 continue
      i5=5
      call pgend
c
      stop
      end
c
c     -------------------------------------------------
c
      subroutine gr_line(command)
c
      character*2   fline
      character*80  filename ,command
      dimension xv(500),yv(500),aaa(20)
c
      common/runit/i5
      common /limit/ xmin,xmax,ymin,ymax
      common /logsc/ ilogx,ilogy
c
      read(command(4:50),'(a2)') fline
      ncol1=1
      ncol2=2
      scale_fac=1.
      shiftx_fac=0.
      shifty_fac=0.
      line_w=1
      read(command(7:50),*,end=200,err=200) ncol1,ncol2,scale_fac
     1,shiftx_fac,shifty_fac,line_w
  200 continue
c
  300 continue
      print *,'Enter file-name to be plotted '
      read(i5,'(a)') filename
      if(filename.ne.'sys$input')then
         iread=69
         open(unit=iread,file=filename,status='old',err=300)
        else
         iread=i5
      endif
      ncolmax=max(ncol1,ncol2)
      npoint=0
      do  n=1,10000
         read(iread,*,end=150,err=150) (aaa(m),m=1,ncolmax)
         xxx=shiftx_fac+aaa(ncol1)
         if(ilogx.eq.1) xxx=log10(abs(xxx)+1.e-23)
         yyy=shifty_fac+scale_fac*aaa(ncol2)
         if(ilogy.eq.1) yyy=log10(abs(yyy)+1.e-23)
         npoint=npoint+1
         xv(npoint)=xxx
         yv(npoint)=yyy
      enddo
  150 continue
      close(unit=69)
c
      call pgslw(line_w)
      if(fline.eq.'--')then
         iline=1
         call pgsls(iline)
         call pgline(npoint,xv,yv)
         iline=1
         call pgsls(iline)
        else if(fline.eq.'- '.or.fline.eq.' -') then
         iline=2
         call pgsls(iline)
         call pgline(npoint,xv,yv)
         iline=1
         call pgsls(iline)
        else if(fline.eq.'-.'.or.fline.eq.'.-') then
         iline=3
         call pgsls(iline)
         call pgline(npoint,xv,yv)
         iline=1
         call pgsls(iline)
        else if(fline.eq.'..') then
         iline=4
         call pgsls(iline)
         call pgline(npoint,xv,yv)
         iline=1
         call pgsls(iline)
      endif
      line_w=1
      call pgslw(line_w)
c
      return
      end
c
c     -------------------------------------------------
c
      subroutine gr_line_clip(command)
c
      character*2   fline
      character*80  filename ,command
      dimension xv(500),yv(500),aaa(20)
c
      common/runit/i5
      common /limit/ xmin,xmax,ymin,ymax
      common /logsc/ ilogx,ilogy
c
      read(command(4:50),'(a2)') fline
      ncol1=1
      ncol2=2
      scale_fac=1.
      shiftx_fac=0.
      shifty_fac=0.
      line_w=1
      read(command(7:50),*,end=200,err=200) ncol1,ncol2,scale_fac
     1,shiftx_fac,shifty_fac,line_w
  200 continue
c
  300 continue
      print *,'Enter file-name to be plotted '
      read(i5,'(a)') filename
      if(filename.ne.'sys$input')then
         iread=69
         open(unit=iread,file=filename,status='old',err=300)
        else
         iread=i5
      endif
      ncolmax=max(ncol1,ncol2)
      npoint=0
      do  n=1,10000
         read(iread,*,end=150,err=150) (aaa(m),m=1,ncolmax)
         xxx=shiftx_fac+aaa(ncol1)
         if(ilogx.eq.1) xxx=log10(abs(xxx)+1.e-23)
         yyy=shifty_fac+scale_fac*aaa(ncol2)
         if(ilogy.eq.1) yyy=log10(abs(yyy)+1.e-23)
         xtmp=(xxx-xmin)*(xxx-xmax)
         ytmp=(yyy-ymin)*(yyy-ymax)
         if(xtmp.le.0.and.ytmp.le.0.)then
            npoint=npoint+1
            xv(npoint)=xxx
            yv(npoint)=yyy
         endif
      enddo
  150 continue
      close(unit=69)
c
      call pgslw(line_w)
      if(fline.eq.'--')then
         iline=1
         call pgsls(iline)
         call pgline(npoint,xv,yv)
         iline=1
         call pgsls(iline)
        else if(fline.eq.'- '.or.fline.eq.' -') then
         iline=2
         call pgsls(iline)
         call pgline(npoint,xv,yv)
         iline=1
         call pgsls(iline)
        else if(fline.eq.'-.'.or.fline.eq.'.-') then
         iline=3
         call pgsls(iline)
         call pgline(npoint,xv,yv)
         iline=1
         call pgsls(iline)
        else if(fline.eq.'..') then
         iline=4
         call pgsls(iline)
         call pgline(npoint,xv,yv)
         iline=1
         call pgsls(iline)
      endif
      line_w=1
      call pgslw(line_w)
c
      return
      end
c
c     ---------------------------------------
c
      subroutine gr_point(command)
c
      common /runit/ i5
      common /logsc/ ilogx,ilogy
c
      character*80  filename  ,command
      dimension aaa(20),xv(2),yv(2)
c
      ncol1=1
      ncol2=2
      scale_fac=1.
      shiftx_fac=0.
      shifty_fac=0.
      kmark=10
      read(command(4:50),*,err=200,end=200) kmark,ncol1,ncol2,scale_fac,
     1    shiftx_fac,shifty_fac
  200 continue
c
  300 continue
      print *,'Enter file-name to be plotted '
      read(i5,'(a80)') filename
      if(filename.ne.'sys$input')then
         iread=69
         open(unit=iread,file=filename,status='old',err=300)
        else
         iread=i5
      endif
      ncolmax=max(ncol1,ncol2)
c
      do  n=1,50000
         read(69,*,end=150,err=150) (aaa(m),m=1,ncolmax)
         xv(1)=shiftx_fac+aaa(ncol1)
         if(ilogx.eq.1)xv(1)=log10(abs(xv(1))+1.e-23)
         yv(1)=shifty_fac+scale_fac*aaa(ncol2)
         if(ilogy.eq.1)yv(1)=log10(abs(yv(1))+1.e-23)
         call pgpt(1,xv,yv,kmark)
      enddo
  150 continue
      close(unit=69)
c
      return
      end
c
c     ---------------------------------------
c
      subroutine gr_point_sc(command)
c
      common /runit/ i5
      common /logsc/ ilogx,ilogy
c
      character*80  filename  ,command
      dimension aaa(20),xv(2),yv(2)
c
      ncol1=1
      ncol2=2
      scale_fac=1.
      shiftx_fac=0.
      shifty_fac=0.
      kmark=10
      scale_ch=1.
      read(command(4:50),*,err=200,end=200) kmark,ncol1,ncol2,scale_fac,
     1    shiftx_fac,shifty_fac,scale_ch
  200 continue
c
  300 continue
      print *,'Enter file-name to be plotted '
      read(i5,'(a80)') filename
      if(filename.ne.'sys$input')then
         iread=69
         open(unit=iread,file=filename,status='old',err=300)
        else
         iread=i5
      endif
      ncolmax=max(ncol1,ncol2)
c
      if(scale_ch.gt.1.000001) then
         line_w=4
         call pgslw(line_w)
         call pgsch(scale_ch)
      endif
c
      do  n=1,50000
         read(69,*,end=150) (aaa(m),m=1,ncolmax)
         xv(1)=shiftx_fac+aaa(ncol1)
         if(ilogx.eq.1)xv(1)=log10(abs(xv(1))+1.e-23)
         yv(1)=shifty_fac+scale_fac*aaa(ncol2)
         if(ilogy.eq.1)yv(1)=log10(abs(yv(1))+1.e-23)
         call pgpt(1,xv,yv,kmark)
      enddo
  150 continue
      close(unit=69)
      line_w=1
      scale_ch=1.
      call pgslw(line_w)
      call pgsch(scale_ch)
c
      return
      end
c
c     ---------------------------------------
c
      subroutine gr_errorv(command)
c
      common /runit/ i5
      common /logsc/ ilogx,ilogy
      common /limit/ xmin,xmax,ymin,ymax
c
      character*80  filename  ,command
      dimension aaa(3),xv(2),yv(2),ev(2)
c
      scale_fac=1.
      shiftx_fac=0.
      shifty_fac=0.
      kmark=10
      scale_ch=1.
      read(command(4:50),*,err=200,end=200) kmark,ncol1,ncol2,scale_fac,
     1    shiftx_fac,shifty_fac,scale_ch
  200 continue
c
  300 continue
      print *,'Enter file-name to be plotted '
      read(i5,'(a80)') filename
c      open(unit=69,file=filename,status='old',err=300)
      if(filename.ne.'sys$input')then
         iread=69
         open(unit=iread,file=filename,status='old',err=300)
        else
         iread=i5
      endif
      ncolmax=3
c
      if(scale_ch.gt.1.000001) then
         line_w=4
         call pgslw(line_w)
         call pgsch(scale_ch)
      endif
c
      do  n=1,50000
         read(iread,*,end=150,err=150) (aaa(m),m=1,ncolmax)
         xv(1)=shiftx_fac+aaa(1)
         if(ilogy.eq.1) then
            aux=shifty_fac+scale_fac*aaa(2)
            yv(1)=log10(abs(aux)+1.e-23)
            ev(1)=log10(abs(aux+scale_fac*aaa(3)))-yv(1)
            call pgpt(1,xv,yv,kmark)
            call pgerrb(2,1,xv,yv,ev,0.0)
            ev(1)=yv(1)-log10(abs(aux-scale_fac*aaa(3)))
            call pgerrb(4,1,xv,yv,ev,0.0)
           else
            aux=shifty_fac+scale_fac*aaa(2)
cc            yv(1)=log10(aux)            ! avevo lasciato sto pasticcio
            yv(1)=aux
            ev(1)=scale_fac*aaa(3)
            call pgpt(1,xv,yv,kmark)
            call pgerrb(6,1,xv,yv,ev,0.0)
         endif
      enddo
  150 continue
      close(unit=69)
c
      line_w=1
      scale_ch=1.
      call pgslw(line_w)
      call pgsch(scale_ch)
c
      return
      end
c
c     ---------------------------------------
c
      subroutine gr_errorvo(command)
c
      common /runit/ i5
c
      character*80  filename  ,command
      dimension aaa(4),xv(2),yv(2),ev(2),eo(2)
c
      scale_fac=1.
      shiftx_fac=0.
      shifty_fac=0.
      kmark=10
      scale_ch=1.
      read(command(4:50),*,err=200,end=200) kmark,ncol1,ncol2,scale_fac,
     1    shiftx_fac,shifty_fac,cale_ch
  200 continue
c
  300 continue
      print *,'Enter file-name to be plotted '
      read(i5,'(a80)') filename
c      open(unit=69,file=filename,status='old',err=300)
      if(filename.ne.'sys$input')then
         iread=69
         open(unit=iread,file=filename,status='old',err=300)
        else
         iread=i5
      endif
      ncolmax=4
c
      if(scale_ch.gt.1.000001) then
         line_w=4
         call pgslw(line_w)
         call pgsch(scale_ch)
      endif
c
      do  n=1,50000
         read(iread,*,end=150,err=150) (aaa(m),m=1,ncolmax)
         xv(1)=shiftx_fac+aaa(1)
         yv(1)=shifty_fac+scale_fac*aaa(2)
         ev(1)=scale_fac*aaa(3)
         eo(1)=scale_fac*aaa(4)

         call pgpt(1,xv,yv,kmark)
         call pgerrb(6,1,xv,yv,ev,0.)
         call pgerrb(5,1,xv,yv,eo,0.)
      enddo
  150 continue
      close(unit=69)
c
      line_w=1
      scale_ch=1.
      call pgslw(line_w)
      call pgsch(scale_ch)
c
      return
      end
c
c     -------------------------------------------------
c
      subroutine gr_hist(command)
c
      character*1   shd
      character*2   fline
      character*80  filename ,command
      dimension xv(800),yv(800),aaa(20)
      dimension xx(4),yy(4)
c
      common/runit/i5
      common /logsc/ ilogx,ilogy
      common /limit/ xmin,xmax,ymin,ymax
      common /frame/ dx,dy,itickx,iticky
c
      read(command(4:50),'(a2)') fline
      kpat=0.
      ncol1=1
      ncol2=2
      scale_fac=1.
      shiftx_fac=0.
      shifty_fac=0.
      line_w=1
c
      read(command(7:50),'(a1)',end=200,err=200) shd
      if(shd.eq.'s'.or.shd.eq.'S') then
         read(command(8:50),*,end=200,err=200) kpat,ncol1,ncol2,
     1        scale_fac,shiftx_fac,shifty_fac,line_w
         else
         read(command(7:50),*,end=200,err=200) ncol1,ncol2,
     1        scale_fac,shiftx_fac,shifty_fac,line_w
      endif
  200 continue
c
  300 continue
      print *,'Enter file-name to be plotted '
      read(i5,'(a)') filename
      open(unit=69,file=filename,status='old',err=300)
      ncolmax=max(ncol1,ncol2)
      do  n=1,1000
         read(69,*,end=150) (aaa(m),m=1,ncolmax)
         xv(n)=shiftx_fac+aaa(ncol1)
         if(ilogx.eq.1) xv(n)=log10(abs(xv(n))+1.e-23)
         yv(n)=shifty_fac+scale_fac*aaa(ncol2)
         if(ilogy.eq.1) then
            if(abs(yv(n)).lt.abs(10.**ymin)) yv(n)=10.**ymin
            yv(n)=log10(abs(yv(n))+1.e-23)
         endif
      enddo
  150 continue
      close(unit=69)
      npoint=n-1
c
c     attenzione queste correzioni devono ancora essere verificate
c     completamente
      do n=npoint,1,-1
         xv(n+1)=xv(n)
         yv(n+1)=yv(n)
      enddo
      xv(1)=xv(2)+(xv(2)-xv(3))
      yv(1)=ymin
      npoint=npoint+1
      xv(npoint+1)=xv(npoint)+(xv(npoint)-xv(npoint-1))
      yv(npoint+1)=ymin
      npoint=npoint+1
c
      call pgslw(line_w)
c           
      call gr_defpattern(kpat,ifs,angle,sepn,phase)
c     
      if(kpat.eq.0) then
         dx=xv(2)-xv(1)
         do n=1,npoint
            xv(n)=xv(n)-dx
         enddo
c
         if(fline.eq.'--')then
            iline=1
            call pgsls(iline)
            call pgbin (npoint, xv, yv, .false.)
            iline=1
            call pgsls(iline)
           else if(fline.eq.'- '.or.fline.eq.' -') then
            iline=2
            call pgsls(iline)
            call pgbin (npoint, xv, yv, .false.)
            iline=1
            call pgsls(iline)
           else if(fline.eq.'-.'.or.fline.eq.'.-') then
            iline=3
            call pgsls(iline)
            call pgbin (npoint, xv, yv, .false.)
            iline=1
            call pgsls(iline)
          else if(fline.eq.'..') then
            iline=4
            call pgsls(iline)
            call pgbin (npoint, xv, yv, .false.)
            iline=1
            call pgsls(iline)
         endif
        else
         call pgshs(angle,sepn,phase)
         nvert=4
         iline=1
         call pgsls(iline)
         call pgsfs(ifs)
         do n=1,npoint-1
cc            if((xv(n)-xmin)*(xv(n)-xmax).lt.0.)then
               dx=xv(n+1)-xv(n)
               xx(1)=xv(n)-dx
               yy(1)=ymin
               xx(2)=xx(1)
               yy(2)=yv(n)
               xx(3)=xv(n)
               yy(3)=yy(2)
               xx(4)=xx(3)
               yy(4)=ymin
               call pgpoly(nvert,xx,yy)
cc            endif
         enddo
c
         n=npoint
         xx(1)=xv(n)-dx
         yy(1)=ymin
         xx(2)=xx(1)
         yy(2)=yv(n)
         xx(3)=xv(n)
         yy(3)=yy(2)
         xx(4)=xx(3)
         yy(4)=ymin
         call pgpoly(nvert,xx,yy)
c
c        drow the borderingline of the hist. box
         do n=1,npoint
            xv(n)=xv(n)-dx
         enddo
c
         if(fline.eq.'--')then
            iline=1
            call pgsls(iline)
            call pgbin (npoint, xv, yv, .false.)
            iline=1
            call pgsls(iline)
           else if(fline.eq.'- '.or.fline.eq.' -') then
            iline=2
            call pgsls(iline)
            call pgbin (npoint, xv, yv, .false.)
            iline=1
            call pgsls(iline)
           else if(fline.eq.'-.'.or.fline.eq.'.-') then
            iline=3
            call pgsls(iline)
            call pgbin (npoint, xv, yv, .false.)
            iline=1
            call pgsls(iline)
           else if(fline.eq.'..') then
            iline=4
            call pgsls(iline)
            call pgbin (npoint, xv, yv, .false.)
            iline=1
            call pgsls(iline)
         endif
      endif
c
      itickx=1
      call pgsci(1)
      call pgbox('bc',dx,itickx,'bcts',dy,iticky)
c
      line_w=1
      call pgslw(line_w)
c
      return
      end
c
c     --------------------------------------------------
c
      subroutine gr_hist_center(command)
c
      character*1   shd
      character*2   fline
      character*80  filename ,command
      dimension xv(800),yv(800),aaa(20)
      dimension xx(4),yy(4)
c
      common/runit/i5
      common /logsc/ ilogx,ilogy
      common /limit/ xmin,xmax,ymin,ymax
      common /frame/ dx,dy,itickx,iticky
c
      read(command(4:50),'(a2)') fline
      kpat=0.
      ncol1=1
      ncol2=2
      scale_fac=1.
      shiftx_fac=0.
      shifty_fac=0.
      line_w=1
c
      read(command(7:50),'(a1)',end=200,err=200) shd
      if(shd.eq.'s'.or.shd.eq.'S') then
         read(command(8:50),*,end=200,err=200) kpat,ncol1,ncol2,
     1        scale_fac,shiftx_fac,shifty_fac,line_w
         else
         read(command(7:50),*,end=200,err=200) ncol1,ncol2,
     1        scale_fac,shiftx_fac,shifty_fac,line_w
      endif
  200 continue
c
  300 continue
      print *,'Enter file-name to be plotted '
      read(i5,'(a)') filename
      open(unit=69,file=filename,status='old',err=300)
      ncolmax=max(ncol1,ncol2)
      do  n=1,1000
         read(69,*,end=150) (aaa(m),m=1,ncolmax)
         xv(n)=shiftx_fac+aaa(ncol1)
         if(ilogx.eq.1) xv(n)=log10(abs(xv(n))+1.e-23)
         yv(n)=shifty_fac+scale_fac*aaa(ncol2)
         if(ilogy.eq.1) then
            if(abs(yv(n)).lt.abs(10.**ymin)) yv(n)=10.**ymin
            yv(n)=log10(abs(yv(n))+1.e-23)
         endif
      enddo
  150 continue
      close(unit=69)
      npoint=n-1
c
c     attenzione queste correzioni devono ancora essere verificate
c     completamente
      do n=npoint,1,-1
         xv(n+1)=xv(n)
         yv(n+1)=yv(n)
      enddo
      xv(1)=xv(2)+(xv(2)-xv(3))
      yv(1)=ymin
      npoint=npoint+1
      xv(npoint+1)=xv(npoint)+(xv(npoint)-xv(npoint-1))
      yv(npoint+1)=ymin
      npoint=npoint+1
c
      call pgslw(line_w)
c
      call gr_defpattern(kpat,ifs,angle,sepn,phase)
c
      if(kpat.eq.0) then
c
         if(fline.eq.'--')then
            iline=1
            call pgsls(iline)
            call pgbin (npoint, xv, yv, .true.)
            iline=1
            call pgsls(iline)
           else if(fline.eq.'- '.or.fline.eq.' -') then
            iline=2
            call pgsls(iline)
            call pgbin (npoint, xv, yv, .true.)
            iline=1
            call pgsls(iline)
           else if(fline.eq.'-.'.or.fline.eq.'.-') then
            iline=3
            call pgsls(iline)
            call pgbin (npoint, xv, yv, .true.)
            iline=1
            call pgsls(iline)
          else if(fline.eq.'..') then
            iline=4
            call pgsls(iline)
            call pgbin (npoint, xv, yv, .true.)
            iline=1
            call pgsls(iline)
         endif
c
        else
c
         call pgshs(angle,sepn,phase)
         nvert=4
         iline=1
         call pgsls(iline)
         call pgsfs(ifs)
         do n=1,npoint-1
cc            if((xv(n)-xmin)*(xv(n)-xmax).lt.0.)then
               dx=xv(n+1)-xv(n)
               xx(1)=xv(n)-dx/2.
               yy(1)=ymin
               xx(2)=xx(1)
               yy(2)=yv(n)
               xx(3)=xv(n)+dx/2.
               yy(3)=yy(2)
               xx(4)=xx(3)
               yy(4)=ymin
               call pgpoly(nvert,xx,yy)
cc            endif
         enddo
c
         n=npoint
         xx(1)=xv(n)-dx/2.
         yy(1)=ymin
         xx(2)=xx(1)
         yy(2)=yv(n)
         xx(3)=xv(n)+dx/2.
         yy(3)=yy(2)
         xx(4)=xx(3)
         yy(4)=ymin
         call pgpoly(nvert,xx,yy)
c
c        drow the borderingline of the hist. box
         do n=1,npoint
            xv(n)=xv(n)-dx/2.
         enddo
c
         if(fline.eq.'--')then
            iline=1
            call pgsls(iline)
            call pgbin (npoint, xv, yv, .false.)
            iline=1
            call pgsls(iline)
           else if(fline.eq.'- '.or.fline.eq.' -') then
            iline=2
            call pgsls(iline)
            call pgbin (npoint, xv, yv, .false.)
            iline=1
            call pgsls(iline)
           else if(fline.eq.'-.'.or.fline.eq.'.-') then
            iline=3
            call pgsls(iline)
            call pgbin (npoint, xv, yv, .false.)
            iline=1
            call pgsls(iline)
           else if(fline.eq.'..') then
            iline=4
            call pgsls(iline)
            call pgbin (npoint, xv, yv, .false.)
            iline=1
            call pgsls(iline)
         endif
      endif
c
      itickx=1
      call pgsci(1)
      call pgbox('bc',dx,itickx,'bcts',dy,iticky)
c
      line_w=1
      call pgslw(line_w)
c
      return
      end
c
c     --------------------------------------------------
c
      subroutine gr_text(command)
c
      character*80   command
      character*78   textr
      character*200  text
c
      common /runit/ i5
      common /logsc/ ilogx,ilogy
      common /limit/ xmin,xmax,ymin,ymax
c
      fchscal=1.
      read(command(5:50),*,err=100,end=100) fchscal
  100 continue
c
      print *,' Enter your text'
      read(i5,'(A)') textr
      call rosetta(textr,text)
      ltext=index(text,'   ')-1
      print *,' Enter [x,y] of text in current axis system'
      read(i5,*) xval,yval
      if(ilogy.eq.1) yval=log10(abs(yval))
      if(ilogx.eq.1) xval=log10(abs(xval))
      angle=0.
      fjust=0.
      call pgsch(fchscal)
      call pgptxt(xval,yval,angle,fjust,text(1:ltext))
      call pgsch(1.)
c
      return
      end
c
c     --------------------------------------------------
c
      subroutine gr_cont(command)
c
      character*2   fline
      character*80  filename
      character*80  command
      parameter (nxdim=301,nydim=301)
      dimension zh(20),zread(nxdim,nydim),tr(6)
c
      common /runit/ i5
      common/limit/xmin,xmax,ymin,ymax
c
      do ny=1,nydim
         do nx=1,nxdim
            zread(nx,ny)=-1.e32
         enddo
      enddo
c
  100 continue
      print *,'Enter file-name to be plotted '
      read(i5,'(a)') filename
      open(unit=69,file=filename,status='old',err=100)
      read (69,*) sxmin,sxmax,sdx,symin,symax,sdy
      nxtot=(sxmax-sxmin)/sdx+1.01               
      nytot=(symax-symin)/sdy+1.01               
      do  ny=1,nytot
         read(69,*) (zread(nx,ny),nx=1,nxtot)
      enddo
      close(unit=69)
c
c     mapping the (x-y) surface (zread) region in the display area
c     defined in the common block /area/
c
      print *,' xmin,xsmin ',xmin,sxmin
      nx_start=(xmin-sxmin)/sdx+1.01
      ny_start=(ymin-symin)/sdy+1.01
      idimx=(xmax-xmin)/sdx+1.01
      idimy=(ymax-ymin)/sdy+1.01
      nx_end=nx_start+idimx-1
      ny_end=ny_start+idimy-1
      if(nx_start.le.0.or.ny_start.le.0) then
         print *,' ------------------------------------------------'
         print *,' The  display area does not agree with the region'
         print *,' of definition of the surface. SORRY, START AGAIN'
         print *,' ------------------------------------------------'
         return
      endif
c
c     definition for the trasformation (see pgplot manual
      tr(1)=sxmin-sdx/2.
      tr(2)=sdx
      tr(3)=0.
      tr(4)=symin-sdy/2.
      tr(5)=0.
      tr(6)=sdy
c
c     define the max/min of the surface  to be plotted
      i=1
      z0=1.e+20
      z1=-1.e+20
      do ny=ny_start,ny_end
         do nx=nx_start,nx_end
            z1=max(zread(nx,ny),z1)
            z0=min(zread(nx,ny),z0)
         enddo
      enddo
      print *,' zmin,zmax: ',z0,z1
c
  200 continue
      print *,' Enter line_stile and contours height [zmin,zmaz,dz]'
      read(i5,'(a)')command
      if(command.eq.'end'.or.command.eq.'END') go to 300
      read(command,'(a2)')fline
      read(command(4:50),*,err=200,end=200)zmin,zmax,dz
c
      nzh=(zmax-zmin)/dz+1
      zh(1)=zmin
      do nz=2,nzh
         zh(nz)=zh(nz-1)+dz
      enddo
c
      if(fline.eq.'--')then
         iline=1
         call pgsls(iline)
        else if(fline.eq.'- '.or.fline.eq.' -') then
         iline=2
         call pgsls(iline)
        else if(fline.eq.'-.'.or.fline.eq.'.-') then
         iline=3
         call pgsls(iline)
        else if(fline.eq.'..') then
         iline=4
         call pgsls(iline)
      endif
c
      call pgcont(zread,nxdim,nydim
     1                 ,nx_start,nx_end
     2                 ,ny_start,ny_end
     3                 ,zh,nzh,tr)
c
      iline=1
      call pgsls(iline)
c
      go to 200
  300 continue
      return
      end
c
c     ---------------------------------------------------
c
      subroutine gr_mosaic(command)
c
      character*2   fline
      character*2   side
      character*80  filename
      character*80  command
      character*8 dev
      parameter (nxdim=1201,nydim=1201)
      dimension zh(20),zread(nxdim,nydim),tr(6)
c
      common /device/ dev
      common /runit/ i5
      common/limit/xmin,xmax,ymin,ymax
      common /frame/ dx,dy,itickx,iticky
c
c
      itype=1
      izlog=0
      read(command(4:50),*,end=99,err=99) itype,izlog
   99 continue
c
  100 continue
      print *,'Enter file-name to be plotted '
      read(i5,'(a)') filename
      open(unit=69,file=filename,status='old',err=100)
      read (69,*) sxmin,sxmax,sdx,symin,symax,sdy
      nxtot=(sxmax-sxmin)/sdx+1.01               
      nytot=(symax-symin)/sdy+1.01               
      do  ny=1,nytot
         read(69,*) (zread(nx,ny),nx=1,nxtot)
      enddo
      close(unit=69)
      if(izlog.eq.1) then
         do ny=1,nytot
            do nx=1,nxtot
               zread(nx,ny)=log10(zread(nx,ny))
            enddo
         enddo
      endif            
c
c     mapping the (x-y) surface (zread) region in the display area
c     defined in the common block /area/
c
      print *,' xmin,xsmin ',xmin,sxmin
      nx_start=(xmin-sxmin)/sdx+1.01
      ny_start=(ymin-symin)/sdy+1.01
      idimx=(xmax-xmin)/sdx+1.01
      idimy=(ymax-ymin)/sdy+1.01
      nx_end=nx_start+idimx-1
      ny_end=ny_start+idimy-1
      if(nx_start.le.0.or.ny_start.le.0) then
         print *,' ------------------------------------------------'
         print *,' The  display area does not agree with the region'
         print *,' of definition of the surface. SORRY, START AGAIN'
         print *,' ------------------------------------------------'
         return
      endif
c
c     definition for the trasformation (see pgplot manual
      tr(1)=sxmin-sdx/2.
      tr(2)=sdx
      tr(3)=0.
      tr(4)=symin-sdy/2.
      tr(5)=0.
      tr(6)=sdy
c
c     define the max/min of the surface  to be plotted
      i=1
      z0=1.e+20
      z1=-1.e+20
      do ny=ny_start,ny_end
         do nx=nx_start,nx_end
            z1=max(zread(nx,ny),z1)
            z0=min(zread(nx,ny),z0)
         enddo
      enddo
      print *,' zmin,zmax: ',z0,z1
      print *,' Enter zmin,zmax'
      read(i5,*)z0,z1
      print *,' Enter annotation'
      read(i5,'(a)') command
      side='no'
      disp=4.
      width=5.
      read(command(1:1),'(a)',err=125,end=125) side(1:1)
  125 continue
      read(command(3:50),*,err=126,end=126) disp,width
  126 continue
c
      if(itype.eq.0) then
         zmin=z0
         zmax=z1
         ips=index(dev,'ps')
         if(ips.ne.0)then
            zmin=z1
            zmax=z0
         endif      
c
         call pggray(zread,nxdim,nydim
     1                    ,nx_start,nx_end
     2                    ,ny_start,ny_end
     3                    ,zmin,zmax,tr)
         if(ips.eq.0) call pgsci(0)
         call pgbox('bcts',dx,itickx,'bcts',dy,iticky)
        else
         zmin=z0
         zmax=z1
         icps=index(dev,'cps')
         if(icps.eq.0)then
            ips=index(dev,'ps')
            if(ips.ne.0)then
               itype=1
               zmin=z1
               zmax=z0
            endif      
         endif
c
         call pgsci(1)
         bright = 0.5
         contra  = 1.0
         call gr_palett(itype, contra, bright)
         call pgimag(zread,nxdim,nydim
     1                    ,nx_start,nx_end
     2                    ,ny_start,ny_end
     3                    ,zmin,zmax,tr)
         call pgbox('bcts',dx,itickx,'bcts',dy,iticky)
         if(side.ne.'no') then
            side=side(1:1)//'I'
            call pgsch(0.6)
            call pgwedg(side,disp,width,zmin,zmax,' ')
            call pgsch(1.)
         endif
      endif
c
      return
      end
c
c     ------------------------------------------------------------
c
      subroutine gr_legenda(command)
c
      character*80  command 
      character*78  textr 
      character*200 text
      character*2   where   ,fline
      character*1   curve   ,col     ,curv  
      character*1   mycl(0:8)
      character*8   dev
c
      dimension xl(2),yl(2),xx(5),yy(5)
c
      data mycl /'n', 'w', 'r', 'g', 'b', 'c', 'v',
     1           'y', 'm'/
c
      common /device/ dev
      common /runit/ i5
      common /logsc/ ilogx,ilogy
      common /limit/ xmin,xmax,ymin,ymax
      common /dims/  pgx ,pgy ,arx ,ary ,pox ,poy,   scx,  scy
      common /frame/ dx,dy,itickx,iticky
c
      where='ru'
      read(command(5:50),'(a2)',end=200,err=200) where
      shift_fac=0.
      read(command(8:50),*,end=200,err=200) shift_fac
      print * , shift_fac
  200 continue
c
      fchscal=1.
c
      scax=(xmax-xmin)/arx
      scay=(ymax-ymin)/ary
      dxl=1.*scax
      dyl=1.*scay
      if(where.eq.'ru'.or.where.eq.'ur')then
         x0=(xmax+xmin)/2.+0.5*dxl+shift_fac
         y0=ymax-0.8*dyl
         x1=x0
         y1=y0
         do nl=1,10
            print *,' Enter :'
            print *,' col,curve,kind,text for legenda'
            read(i5,'(a)') command
            if(command.eq.'end') return
            read(command,'(a1,1x,a1)') col,curv
            do n=0,8
              if(col.eq.mycl(n)) goto 104
            enddo
  104       continue
            call pgsci(n)
            if(curv.eq.'c')then
               read(command(4:50),'(a2)') fline
               xl(1)=x1
               xl(2)=x1+1.3*dxl
               yl(1)=y1
               yl(2)=y1
               npt=2
               if(fline.eq.'--')then
                  iline=1
                  call pgsls(iline)
                  call pgline(npt,xl,yl)
                  iline=1
                  call pgsls(iline)
                 else if(fline.eq.'- '.or.fline.eq.' -') then
                  iline=2
                  call pgsls(iline)
                  call pgline(npt,xl,yl)
                  iline=1
                  call pgsls(iline)
                 else if(fline.eq.'-.'.or.fline.eq.'.-') then
                  iline=3
                  call pgsls(iline)
                  call pgline(npt,xl,yl)
                  iline=1
                  call pgsls(iline)
                 else if(fline.eq.'..') then
                  iline=4
                  call pgsls(iline)
                  call pgline(npt,xl,yl)
                  iline=1
                  call pgsls(iline)
               endif
               read(command(7:50),'(a)') textr
               call rosetta(textr,text)
               ltext=index(text,'   ')-1
               xtex=x1+1.6*dxl
               ytex=y1-0.1*dyl
               angle=0.
               fjust=0.
               call pgslw(2)
               call pgsch(fchscal)
               call pgptxt(xtex,ytex,angle,fjust,text(1:ltext))
               call pgsch(1.)
               call pgslw(1)
              else if(curv.eq.'p') then
               read(command(4:50),*) kmark
               npoint=1
               xp=x1+0.65*dxl
               yp=y1
c
               call pgpt(1,xp,yp,kmark)
c
               read(command(7:50),'(a)') textr
               call rosetta(textr,text)
               ltext=index(text,'   ')-1
               xtex=x1+1.6*dxl
               ytex=y1-0.1*dyl
               angle=0.
               fjust=0.
               call pgslw(2)
               call pgsch(fchscal)
               call pgptxt(xtex,ytex,angle,fjust,text(1:ltext))
               call pgsch(1.)
               call pgslw(1)
              else if(curv.eq.'s') then
               read(command(4:50),*)kpat
c
               call gr_defpattern(kpat,ifs,angle,sepn,phase)
c
               call pgshs(angle,sepn,phase)
               nvert=4
               iline=1
               call pgsls(iline)
               call pgsfs(ifs)
               dx=1.3*dxl
               dy=0.3*dyl
               xx(1)=x1
               yy(1)=y1-0.5*dy
               xx(2)=xx(1)+dx
               yy(2)=yy(1)
               xx(3)=xx(2)
               yy(3)=yy(2)+dy
               xx(4)=xx(3)-dx
               yy(4)=yy(3)
               xx(5)=xx(1)
               yy(5)=yy(1)
               call pgpoly(nvert,xx,yy)
               nvert=5
               call pgline(nvert,xx,yy)
               read(command(7:50),'(a)') textr
               call rosetta(textr,text)
               ltext=index(text,'   ')-1
               xtex=x1+1.8*dxl
               ytex=y1-0.1*dyl
               angle=0.
               fjust=0.
               call pgslw(2)
               call pgsch(fchscal)
               call pgptxt(xtex,ytex,angle,fjust,text(1:ltext))
               call pgsch(1.)
               call pgslw(1)
            endif
            y1=y1-0.6*dyl
         enddo
c
        else if(where.eq.'lu'.or.where.eq.'ul')then
         x0=xmin+0.8*dxl
         y0=ymax-0.8*dyl
         x1=x0
         y1=y0
         do nl=1,10
            print *,' Enter :'
            print *,' col,curve,kind,text for legenda'
            read(i5,'(a)') command
            if(command.eq.'end') return
            read(command,'(a1,1x,a1)') col,curv
            do n=0,8
               if(col.eq.mycl(n)) goto 103
            enddo
  103       continue
            call pgsci(n)
            if(curv.eq.'c')then
               read(command(4:50),'(a2)') fline
               xl(1)=x1
               xl(2)=x1+1.3*dxl
               yl(1)=y1
               yl(2)=y1
               npt=2
               print *, ' sono qui in lu'
               if(fline.eq.'--')then
                  iline=1
                  call pgsls(iline)
                  call pgline(npt,xl,yl)
                  iline=1
                  call pgsls(iline)
                 else if(fline.eq.'- '.or.fline.eq.' -') then
                  iline=2
                  call pgsls(iline)
                  call pgline(npt,xl,yl)
                  iline=1
                  call pgsls(iline)
                 else if(fline.eq.'-.'.or.fline.eq.'.-') then
                  iline=3
                  call pgsls(iline)
                  call pgline(npt,xl,yl)
                  iline=1
                  call pgsls(iline)
                 else if(fline.eq.'..') then
                  iline=4
                  call pgsls(iline)
                  call pgline(npt,xl,yl)
                  iline=1
                  call pgsls(iline)
               endif
               read(command(7:50),'(a)') textr
               call rosetta(textr,text)
               ltext=index(text,'   ')-1
               xtex=x1+1.6*dxl
               ytex=y1-0.1*dyl
               angle=0.
               fjust=0.
               call pgslw(2)
               call pgsch(fchscal)
               call pgptxt(xtex,ytex,angle,fjust,text(1:ltext))
               call pgsch(1.)
               call pgslw(1)
              else if(curv.eq.'p') then
               read(command(4:50),*) kmark
               npoint=1
               xp=x1+0.65*dxl
               yp=y1
c
               call pgpt(1,xp,yp,kmark)
c
               read(command(7:50),'(a)') textr
               call rosetta(textr,text)
               ltext=index(text,'   ')-1
               xtex=x1+1.6*dxl
               ytex=y1-0.1*dyl
               angle=0.
               fjust=0.
               call pgslw(2)
               call pgsch(fchscal)
               call pgptxt(xtex,ytex,angle,fjust,text(1:ltext))
               call pgsch(1.)
               call pgslw(1)
              else if(curv.eq.'s') then
               read(command(4:50),*)kpat
c
               call gr_defpattern(kpat,ifs,angle,sepn,phase)
c
               call pgshs(angle,sepn,phase)
               nvert=4
               iline=1
               call pgsls(iline)
               call pgsfs(ifs)
               dx=1.3*dxl
               dy=0.3*dyl
               xx(1)=x1
               yy(1)=y1-0.5*dy
               xx(2)=xx(1)+dx
               yy(2)=yy(1)
               xx(3)=xx(2)
               yy(3)=yy(2)+dy
               xx(4)=xx(3)-dx
               yy(4)=yy(3)
               xx(5)=xx(1)
               yy(5)=yy(1)
               call pgpoly(nvert,xx,yy)
               nvert=5
               call pgline(nvert,xx,yy)
               read(command(7:50),'(a)') textr
               call rosetta(textr,text)
               ltext=index(text,'   ')-1
               xtex=x1+1.8*dxl
               ytex=y1-0.1*dyl
               angle=0.
               fjust=0.
               call pgslw(2)
               call pgsch(fchscal)
               call pgptxt(xtex,ytex,angle,fjust,text(1:ltext))
               call pgsch(1.)
               call pgslw(1)
            endif
            y1=y1-0.6*dyl
         enddo
c
        else if(where.eq.'rd'.or.where.eq.'dr')then
         x0=(xmax+xmin)/2.+0.5*dxl+shift_fac
         y0=ymin+0.8*dyl
         x1=x0
         y1=y0
         do nl=1,10
            print *,' Enter :'
            print *,' col,curve,kind,text for legenda'
            read(i5,'(a)') command
            if(command.eq.'end') return
            read(command,'(a1,1x,a1)') col,curv
            do n=0,8
               if(col.eq.mycl(n)) goto 105
            enddo
  105       continue
            call pgsci(n)
            if(curv.eq.'c')then
               read(command(4:50),'(a2)') fline
               xl(1)=x1
               xl(2)=x1+1.3*dxl
               yl(1)=y1
               yl(2)=y1
               npt=2
               if(fline.eq.'--')then
                  iline=1
                  call pgsls(iline)
                  call pgline(npt,xl,yl)
                  iline=1
                  call pgsls(iline)
                 else if(fline.eq.'- '.or.fline.eq.' -') then
                  iline=2
                  call pgsls(iline)
                  call pgline(npt,xl,yl)
                  iline=1
                  call pgsls(iline)
                 else if(fline.eq.'-.'.or.fline.eq.'.-') then
                  iline=3
                  call pgsls(iline)
                  call pgline(npt,xl,yl)
                  iline=1
                  call pgsls(iline)
                 else if(fline.eq.'..') then
                  iline=4
                  call pgsls(iline)
                  call pgline(npt,xl,yl)
                  iline=1
                  call pgsls(iline)
               endif
               read(command(7:50),'(a)') textr
               call rosetta(textr,text)
               ltext=index(text,'   ')-1
               xtex=x1+1.6*dxl
               ytex=y1-0.1*dyl
               angle=0.
               fjust=0.
               call pgslw(2)
               call pgsch(fchscal)
               call pgptxt(xtex,ytex,angle,fjust,text(1:ltext))
               call pgsch(1.)
               call pgslw(1)
              else if(curv.eq.'p') then
               read(command(4:50),*) kmark
               npoint=1
               xp=x1+0.65*dxl
               yp=y1
c
               call pgpt(1,xp,yp,kmark)
c
               read(command(7:50),'(a)') textr
               call rosetta(textr,text)
               ltext=index(text,'   ')-1
               xtex=x1+1.6*dxl
               ytex=y1-0.1*dyl
               angle=0.
               fjust=0.
               call pgslw(2)
               call pgsch(fchscal)
               call pgptxt(xtex,ytex,angle,fjust,text(1:ltext))
               call pgsch(1.)
               call pgslw(1)
              else if(curv.eq.'s') then
               read(command(4:50),*)kpat
c
               call gr_defpattern(kpat,ifs,angle,sepn,phase)
c
               call pgshs(angle,sepn,phase)
               nvert=4
               iline=1
               call pgsls(iline)
               call pgsfs(ifs)
               dx=1.3*dxl
               dy=0.3*dyl
               xx(1)=x1
               yy(1)=y1-0.5*dy
               xx(2)=xx(1)+dx
               yy(2)=yy(1)
               xx(3)=xx(2)
               yy(3)=yy(2)+dy
               xx(4)=xx(3)-dx
               yy(4)=yy(3)
               xx(5)=xx(1)
               yy(5)=yy(1)
               call pgpoly(nvert,xx,yy)
               nvert=5
               call pgline(nvert,xx,yy)
               read(command(7:50),'(a)') textr
               call rosetta(textr,text)
               ltext=index(text,'   ')-1
               xtex=x1+1.8*dxl
               ytex=y1-0.1*dyl
               angle=0.
               fjust=0.
               call pgslw(2)
               call pgsch(fchscal)
               call pgptxt(xtex,ytex,angle,fjust,text(1:ltext))
               call pgsch(1.)
               call pgslw(1)
            endif
            y1=y1+0.6*dyl
         enddo
c
        else if(where.eq.'dl'.or.where.eq.'ld')then
         x0=xmin+0.8*dxl
         y0=ymin+0.8*dyl
         x1=x0
         y1=y0
         do nl=1,10
            print *,' Enter :'
            print *,' col,curve,kind,text for legenda'
            read(i5,'(a)') command
            if(command.eq.'end') return
            read(command,'(a1,1x,a1)') col,curv
            do n=0,8
               if(col.eq.mycl(n)) goto 102
            enddo
  102       continue
            call pgsci(n)
            if(curv.eq.'c')then
               read(command(4:50),'(a2)') fline
               xl(1)=x1
               xl(2)=x1+1.3*dxl
               yl(1)=y1
               yl(2)=y1
               npt=2
               if(fline.eq.'--')then
                  iline=1
                  call pgsls(iline)
                  call pgline(npt,xl,yl)
                  iline=1
                  call pgsls(iline)
                 else if(fline.eq.'- '.or.fline.eq.' -') then
                  iline=2
                  call pgsls(iline)
                  call pgline(npt,xl,yl)
                  iline=1
                  call pgsls(iline)
                 else if(fline.eq.'-.'.or.fline.eq.'.-') then
                  iline=3
                  call pgsls(iline)
                  call pgline(npt,xl,yl)
                  iline=1
                  call pgsls(iline)
                 else if(fline.eq.'..') then
                  iline=4
                  call pgsls(iline)
                  call pgline(npt,xl,yl)
                  iline=1
                  call pgsls(iline)
               endif
               read(command(7:50),'(a)') textr
               call rosetta(textr,text)
               ltext=index(text,'   ')-1
               xtex=x1+1.6*dxl
               ytex=y1-0.1*dyl
               angle=0.
               fjust=0.
               call pgslw(2)
               call pgsch(fchscal)
               call pgptxt(xtex,ytex,angle,fjust,text(1:ltext))
               call pgsch(1.)
               call pgslw(1)
              else if(curv.eq.'p') then
               read(command(4:50),*) kmark
               npoint=1
               xp=x1+0.65*dxl
               yp=y1
c
               call pgpt(1,xp,yp,kmark)
c
               read(command(7:50),'(a)') textr
               call rosetta(textr,text)
               ltext=index(text,'   ')-1
               xtex=x1+1.6*dxl
               ytex=y1-0.1*dyl
               angle=0.
               fjust=0.
               call pgslw(2)
               call pgsch(fchscal)
               call pgptxt(xtex,ytex,angle,fjust,text(1:ltext))
               call pgsch(1.)
               call pgslw(1)
              else if(curv.eq.'s') then
               read(command(4:50),*)kpat
c
               call gr_defpattern(kpat,ifs,angle,sepn,phase)
c
               call pgshs(angle,sepn,phase)
               nvert=4
               iline=1
               call pgsls(iline)
               call pgsfs(ifs)
               dx=1.3*dxl
               dy=0.3*dyl
               xx(1)=x1
               yy(1)=y1-0.5*dy
               xx(2)=xx(1)+dx
               yy(2)=yy(1)
               xx(3)=xx(2)
               yy(3)=yy(2)+dy
               xx(4)=xx(3)-dx
               yy(4)=yy(3)
               xx(5)=xx(1)
               yy(5)=yy(1)
               call pgpoly(nvert,xx,yy)
               nvert=5
               call pgline(nvert,xx,yy)
               read(command(7:50),'(a)') textr
               call rosetta(textr,text)
               ltext=index(text,'   ')-1
               xtex=x1+1.8*dxl
               ytex=y1-0.1*dyl
               angle=0.
               fjust=0.
               call pgslw(2)
               call pgsch(fchscal)
               call pgptxt(xtex,ytex,angle,fjust,text(1:ltext))
               call pgsch(1.)
               call pgslw(1)
            endif
            y1=y1+0.6*dyl
         enddo
      endif
      return
      end
c
c     ------------------------------------------------------------
c
      subroutine gr_palett(type, contra, bright)
c
c      set a "palette" of colors in the range of color indices used by
c      pgimag.
c
      integer type
      real contra, bright
c
      real gl(2), gr(2), gg(2), gb(2)
      real rl(9), rr(9), rg(9), rb(9)
      real hl(5), hr(5), hg(5), hb(5)
      real wl(10), wr(10), wg(10), wb(10)
      real al(20), ar(20), ag(20), ab(20)
c
      data gl /0.0, 1.0/
      data gr /0.0, 1.0/
      data gg /0.0, 1.0/
      data gb /0.0, 1.0/
c
      data rl /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      data rr / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      data rg / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      data rb / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
c
      data hl /0.0, 0.2, 0.4, 0.6, 1.0/
      data hr /0.0, 0.5, 1.0, 1.0, 1.0/
      data hg /0.0, 0.0, 0.5, 1.0, 1.0/
      data hb /0.0, 0.0, 0.0, 0.3, 1.0/
c
      data wl /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      data wr /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      data wg /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      data wb /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
c
      data al /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
     :         0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      data ar /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
     :         0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      data ag /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
     :         0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      data ab /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
     :         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
c
      if (type.eq.1) then
c        -- gray scale
         call pgctab(gl, gr, gg, gb, 2, contra, bright)
      else if (type.eq.2) then
c        -- rainbow
         call pgctab(rl, rr, rg, rb, 9, contra, bright)
      else if (type.eq.3) then
c        -- heat
         call pgctab(hl, hr, hg, hb, 5, contra, bright)
      else if (type.eq.4) then
c        -- weird iraf
         call pgctab(wl, wr, wg, wb, 10, contra, bright)
      else if (type.eq.5) then
c        -- aips
         call pgctab(al, ar, ag, ab, 20, contra, bright)
      end if
      return
      end
c
c     -----------------------------------------------------
c
      subroutine rosetta(string,string2)
c
      character*72  string
      character*200 string2
      character*78  temp,str,str0
      character*2   up,down,clup,cldown
      character*1 backs
c
	backs=char(92)
      up  =backs//'u'
      down=backs//'d'
      clup=backs//'d'
      cldown=backs//'u'
      str0='                                                   '
c
c
      lp2=index(string,'    ')-1
      lp1=1
c  100 continue
      ls=index(string,'    ')
      ls=ls-1
      l1=1
      string2='@'
  150 continue
      ju=index(string(l1:ls),'^')
      jd=index(string(l1:ls),'_')
      if(ju.eq.0.and.jd.eq.0) then
         call dictionary(string(l1:l2),string2)
         goto 999
      endif
      if(jd.eq.0)jd=ls+1
      if(ju.eq.0)ju=ls+1
      if(ju.lt.jd) then
         key=0
         l2=l1+ju-2
         str=str0
         str=string(l1:l2)
         if(ju.ne.1)call dictionary(str,string2)
         l1=l2+3
         j=index(string(l1:ls),'}')
         l2=l1+j-2
         str=str0
         str=string(l1:l2)
         lf=index(string2,'   ')
         string2=string2(1:lf-1)//up   
         call dictionary(str,string2)
         lf=index(string2,'   ')
         string2=string2(1:lf-1)//clup               
         l1=l2+2
         if(l1.le.ls) goto 150
      else
         l2=l1+jd-2
         str=str0
         str=string(l1:l2)
         if(jd.ne.1)call dictionary(str,string2)
         l1=l2+3
         j=index(string(l1:ls),'}')
         l2=l1+j-2
         str=str0
         str=string(l1:l2)
         lf=index(string2,'   ')
         string2=string2(1:lf-1)//down 
         call dictionary(str,string2)
         lf=index(string2,'   ')
         string2=string2(1:lf-1)//cldown    
         l1=l2+2
         if(l1.le.ls) goto 150
      endif
  999 continue
c
      ls=index(string2,'   ')-1
      string2=string2(2:ls)//' '
c
c     from ? to \fs (scrpt)
      lp2=index(string2,'   ')-1
      lp1=1
  100 continue
      lp=index(string2(lp1:lp2),'?')
      if(lp.ne.0) then
         if(lp1.ne.1) lp=lp1+lp-1
         string2=string2(1:lp-1)//backs//'fs'//string2(lp+1:lp2)
         lp1=lp+3
         lp2=index(string2,'   ') - 1
         go to 100
      endif      
c
c     from * to \g (greek)
      lp2=index(string2,'   ')-1
      lp1=1
  110 continue
      lp=index(string2(lp1:lp2),'*')
      if(lp.ne.0) then
         if(lp1.ne.1) lp=lp1+lp-1
         string2=string2(1:lp-1)//backs//'g'//string2(lp+1:lp2)
         lp1=lp+2
         lp2=index(string2,'   ') - 1
         go to 110
      endif      
c
c
c     from ~ to \  (blank)
      lp2=index(string2,'   ')-1
      lp1=1
  120 continue
      lp=index(string2(lp1:lp2),'~')
      if(lp.ne.0) then
         if(lp1.ne.1) lp=lp1+lp-1
         string2=string2(1:lp-1)//' '//string2(lp+1:lp2)
         lp1=lp+1
         lp2=index(string2,'   ') - 1
         go to 120
      endif      
c
c     from # to \  (blank)
      lp2=index(string2,'   ')-1
      lp1=1
  130 continue
      lp=index(string2(lp1:lp2),'#')
      if(lp.ne.0) then
         if(lp1.ne.1) lp=lp1+lp-1
         string2=string2(1:lp-1)//backs//string2(lp+1:lp2)
         lp1=lp+1
         lp2=index(string2,'   ') - 1
         go to 130
      endif
c
      return
      end      
c
c     -------------------
c 
      subroutine dictionary(str1,str2)
      character str1*72,str2*200
      character temp*72,gg*10
      character*10  g(105) 
      character*1 greek
      character*1 backs
      character*3 norm
      character*8 gp(105)
c               1234567890  1234567890  1234567890  1234567890
      data g  /'alpha'    ,'beta'     ,'xi'      ,'delta'    
     1        ,'epsilon'  ,'phi'      ,'gamma'    ,'theta'
     2        ,'iota'     ,'kappa'    ,'lambda'   ,'mu'
     3        ,'nu'       ,'omicron'  ,'pi'       ,'psi'
     4        ,'rho'      ,'sigma'    ,'tau'      ,'upsilon'
     5        ,'omega'    ,'chi'      ,'eta'      ,'zeta'
     6        ,'Delta'    ,'Phi'      ,'Gamma'    ,'Lambda'
     7        ,'Pi'       ,'Theta'    ,'Sigma'    ,'Upsilon'
     8        ,'Omega'    ,'Xi'      ,'Psi'      ,'pm'
     9        ,'mp'       ,'leq'      ,'geq'      ,'propto'
     $        ,'infty'    ,'sim'      ,'partial'  ,'righta'
     $        ,'upa'      ,'lefta'    ,'downa'    ,'sqrt'
     $        ,'int'      ,'Int'      ,'prod'     ,'sum'  
     $        ,'cala'     ,'calb'     ,'calc'     ,'cald'
     $        ,'cale'     ,'calf'     ,'calg'     ,'calh'
     $        ,'cali'     ,'calj'     ,'calk'     ,'call' 
     $        ,'calm'     ,'caln'     ,'calo'     ,'calp' 
     $        ,'calq'     ,'calr'     ,'cals'     ,'calt'   
     $        ,'calu'     ,'calv'     ,'calw'     ,'calx'
     $        ,'caly'     ,'calz'     ,'gzp'
     $        ,'calA'     ,'calB'     ,'calC'     ,'calD' 
     $        ,'calE'     ,'calF'     ,'calG'     ,'calH'
     $        ,'calI'     ,'calJ'     ,'calK'     ,'calL' 
     $        ,'calM'     ,'calN'     ,'calO'     ,'calP' 
     $        ,'calQ'     ,'calR'     ,'calS'     ,'calT'
     $        ,'calU'     ,'calV'     ,'calW'     ,'calX'
     $        ,'calY'     ,'calZ'/
      data gp /'*a'       ,'*b'       ,'*c'       ,'*d'
     1        ,'*e'       ,'*f'       ,'*g'       ,'*h'
     2        ,'*i'       ,'*k'       ,'*l'       ,'*m'
     3        ,'*n'       ,'*o'       ,'*p'       ,'*q' 
     4        ,'*r'       ,'*s'       ,'*t'       ,'*u'
     5        ,'*w'       ,'*x'       ,'*y'       ,'*z'
     6        ,'*D'       ,'*F'       ,'*G'       ,'*L'
     7        ,'*P'       ,'*H'       ,'*S'       ,'*U'
     8        ,'*W'       ,'*C'       ,'*Q'       ,'#(2233)'
     9        ,'#(2234)'  ,'#(2343)'  ,'#(2345)'  ,'#(2245)'
     $        ,'#(2270)'  ,'#(2246)'  ,'#(2265)'  ,'#(2261)'
     $        ,'#(2262)'  ,'#(2263)'  ,'#(2264)'  ,'#(2267)'
     $        ,'#(2268)'  ,'#(2412)'  ,'#(2401)'  ,'#(2402)'
     $        ,'?a'       ,'?b'       ,'?c'       ,'?d'
     $        ,'?e'       ,'?f'       ,'?g'       ,'?h'
     $        ,'?i'       ,'?j'       ,'?k'       ,'?l' 
     $        ,'?m'       ,'?n'       ,'?o'       ,'?p' 
     $        ,'?q'       ,'?r'       ,'?s'       ,'?t'   
     $        ,'?u'       ,'?v'       ,'?w'       ,'?x'
     $        ,'?y'       ,'?z'       ,'?GzP'
     $        ,'?A'       ,'?B'       ,'?C'       ,'?D' 
     $        ,'?E'       ,'?F'       ,'?G'       ,'?H'
     $        ,'?I'       ,'?J'       ,'?K'       ,'?L' 
     $        ,'?M'       ,'?N'       ,'?O'       ,'?P' 
     $        ,'?Q'       ,'?R'       ,'?S'       ,'?T'
     $        ,'?U'       ,'?V'       ,'?W'       ,'?X'
     $        ,'?Y'       ,'?Z' /
c
      data greek /'*'/
      backs =char(92)
      norm=backs//'fr'
      ls=index(str1,'    ')
      ls=ls-1
      l1=1
      l2=ls
  150 continue
      jg=index(str1(l1:l2),backs)
      if(jg.eq.0) then
         lf=index(str2,'   ')
         str2=str2(1:lf-1)//str1(l1:l2)
         goto 999
      endif
      if(jg.eq.1)then
         l1=l1+1
         n=0
 200     continue
         n=n+1
         gg=g(n)
         lc=index(gg,' ')
         lc=lc-1
         j=index(str1(l1:l2),gg(1:lc))
         if(j.ne.1) goto 200
         lf=index(str2,'   ')
c         str2=str2(1:lf-1)//greek//gp(n)//norm
         lgp=index(gp(n),' ')-1
         str2=str2(1:lf-1)//gp(n)(:lgp)//norm
         l1=l1+lc
         if(l1.gt.l2)goto 999
         goto 150
         else
         lp=l1+jg-2
         lf=index(str2,'   ')
         str2=str2(1:lf-1)//str1(l1:lp)
         l1=lp+1
         if(l1.gt.l2)goto 999
         goto 150
      endif
  999 continue
c
      return
      end      
c
c     ----------------------------------------------------
c
      subroutine gr_defpattern(kpat,ifs,angle,sepn,phase)
c
c     definition of filling pattern
      if(kpat.eq.0) then              ! outline
         ifs=2
         angle=0.
         sepn=0.
         phase=0.
        else if(kpat.eq.1) then        ! filled
         ifs=1
         angle=0.
         sepn=0.
         phase=0.
        else if(kpat.eq.2) then        ! hutched 45 degrees
         ifs=3
         angle=45.
         sepn=1.2
         phase=0.
        else if(kpat.eq.3) then        ! hutched 90
         ifs=3
         angle=90.
         sepn=1.2
         phase=0.
        else if(kpat.eq.4)then        ! hutched 135
         ifs=3
         angle=135.
         sepn=1.2
         phase=0.
        elseif(kpat.eq.5) then        ! cross-hutched 45
         ifs=4
         angle=45.
         sepn=1.2
         phase=0.
        else if(kpat.eq.6) then        ! cross-hutched 90
         ifs=4
         angle=90.
         sepn=1.2
         phase=0.
        else if(kpat.eq.7) then        ! cross-hutched 135
         ifs=4
         angle=60.
         sepn=1.2
         phase=0.
      endif   
c
      return
      end
c
c     ---------------------------------------------------------
c
      subroutine gr_text_cm(command)
c
      character*80  command 
      character*78  textr 
      character*200 text
      character*8   dev
c
      common /device/ dev
      common /runit/ i5
      common /limit/ xmin,xmax,ymin,ymax
      common /dims/  pgx ,pgy ,arx ,ary ,pox ,poy,   scx,  scy
      common /frame/ dx,dy,itickx,iticky
c
      scax=(xmax-xmin)/arx
      scay=(ymax-ymin)/ary
c
      hite=1.
      iblank=0
      angle=0.
      read(command(5:50),*,err=200,end=200) hite,angle
  200 continue
c
      print *,' Enter your text'
      read(i5,'(A)') textr
      call rosetta(textr,text)
      ltext=index(text,'   ')-1
      print *,' Enter [x,y] of text in cm from the phys. origin'
      read(i5,*)xcm,ycm
      xxx=xmin+xcm*scax/scx
      yyy=ymin+ycm*scay/scy
c
      fjust=0.
      call pgslw(3)
      call pgsch(hite)
      call pgptxt(xxx,yyy,angle,fjust,text(1:ltext))
      call pgsch(1.)
      call pgslw(1)
c
      return
      end
c
c     ---------------------------------------------------------
c
      subroutine gr_shadow(command)
c
      character*1   curve
      character*2   fline
      character*8   dev
      character*80  filename  ,command
      dimension xv1(500),yv1(500),aaa(20)
      dimension xv2(500),yv2(500)
      character*12 optionx,optiony
c
      common /runit/ i5
      common /device/ dev
      common /logsc/ ilogx,ilogy
      common /frame/ dx,dy,itickx,iticky
      common /opt/ optionx,optiony
c
      kpat=1
      read(command(4:50),'(a2)') fline
      ncol1=1
      ncol2=2
      scale_fac=1.
      shiftx_fac=0.
      shifty_fac=0.
      read(command(7:50),*,end=200,err=200) ncol1,ncol2,scale_fac
     1,shiftx_fac,shifty_fac
  200 continue
c
  300 continue
      print *,'Enter file-name to be plotted '
      read(i5,'(a)') filename
      print *, filename
      if(filename.ne.'sys$input')then
         iread=69
         open(unit=iread,file=filename,status='old',err=300)
        else
         iread=i5
      endif
      ncolmax=max(ncol1,ncol2)
      do  n=1,10000
         read(iread,*,end=150,err=150) (aaa(m),m=1,ncolmax)
         xv1(n)=shiftx_fac+aaa(ncol1)
         if(ilogx.eq.1) xv1(n)=log10(abs(xv1(n)))
         yv1(n)=shifty_fac+scale_fac*aaa(ncol2)
         if(ilogy.eq.1) yv1(n)=log10(abs(yv1(n)))
      enddo
  150 continue
      close(unit=69)
      npoint1=n-1
c
      print *,' Enter (for the second curve):'
      print *,' col,curve,kind,[icol1,icol2,...]'
      read(i5,'(a)')command
c
      read(command(4:50),'(a2)') fline
      ncol1=1
      ncol2=2
      scale_fac=1.
      shiftx_fac=0.
      shifty_fac=0.
      read(command(7:50),*,end=201,err=201) ncol1,ncol2,scale_fac
     1,shiftx_fac,shifty_fac
  201 continue
c
  301 continue
      print *,'Enter file-name to be plotted '
      read(i5,'(a)') filename
      if(filename.ne.'sys$input')then
         iread=69
         open(unit=iread,file=filename,status='old',err=301)
        else
         iread=i5
      endif
      ncolmax=max(ncol1,ncol2)
      do  n=1,10000
         read(iread,*,end=151,err=151) (aaa(m),m=1,ncolmax)
         xv2(n)=shiftx_fac+aaa(ncol1)
         if(ilogx.eq.1) xv2(n)=log10(abs(xv2(n)))
         yv2(n)=shifty_fac+scale_fac*aaa(ncol2)
         if(ilogy.eq.1) yv2(n)=log10(abs(yv2(n)))
      enddo
  151 continue
      close(unit=69)
      npoint2=n-1
c
      print *,' Enter pattern parameter 0...7'
      read (i5,*) kpat
      if(kpat.eq.0) then
         kpat=1
         print *,' warning kpat set to 1'
      endif
c
      npoint=npoint1+npoint2
c
      do n=npoint2,1,-1
         xv1(npoint-n+1)=xv2(n)
         yv1(npoint-n+1)=yv2(n)
      enddo
c
cc      call gr_defpattern(kpat,ifs,angle,sepn,phase,ilpat)
cc
cc    attenzione a cosa e' ilpat che non e mai stato definito
cc
      call gr_defpattern(kpat,ifs,angle,sepn,phase)
      call pgshs(angle,sepn,phase)
      call pgsls(ilpat)
      call pgsfs(ifs)
      call pgpoly(npoint,xv1,yv1)
c
c     redraw the frame after filling
      call pgsci(1)
      call pgsls(1)
      call pgslw(2)
      call pgbox(optionx,dx,itickx,optiony,dy,iticky)
      call pgslw(1)
c
      return
      end
