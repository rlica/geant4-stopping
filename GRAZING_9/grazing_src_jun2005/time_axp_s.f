c      Subroutine delivers time in the format HH:MM
c-------------------------------------------
      subroutine time_s(ttt)
      implicit real*8(a-h,o-z)
      character*8  date
      character*10 time
      character*5  zone 
      integer*4 value(8)
      logical add0
      character*5 ttt
      character*10 ttt1
      character*3 hh
      character*3 mm
      add0 = .false.
      call date_and_time(date,time,zone,value)
      ihh = value(5)
      call convert_s(ihh,hh,ierr)
      imm = value(6)
      if(imm.lt.10)add0 = .true.
      call convert_s(imm,mm,ierr)
      if(add0)then
         ttt1 = hh//':'//'0'//mm
        else
         ttt1 = hh//':'//mm
      endif
      call condens_s(ttt1,ttt,len)
c
      return
      end 






