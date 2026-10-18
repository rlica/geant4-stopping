c      Subroutine delivers time in the format HH:MM
c-------------------------------------------
      subroutine time_s(ttt)
cnan      implicit real*8(a-h,o-z)
cnan      integer*4 value(8)
cnan      logical add0
cnan      character*5 ttt
cnan      character*10 ttt1
cnan      character*3 hh
cnan      character*3 mm
cnan      add0 = .false.
cnan      call date_and_time(,,,value)
cnan      ihh = value(5)
cnan      call convert_s(ihh,hh,ierr)
cnan      imm = value(6)
cnan      if(imm.lt.10)add0 = .true.
cnan      call convert_s(imm,mm,ierr)
cnan      if(add0)then
cnan       ttt1 = hh//':'//'0'//mm
cnan      else
cnan       ttt1 = hh//':'//mm
cnan      endif
cnan      call condens_s(ttt1,ttt,len)
cnan
cnan  ---------------------------------------------------------------
cnan  NB:
cnan  un-comment all the  above lines if you un on alpha and comment
cnan  all the lines below that are for linux
cnan  ---------------------------------------------------------------
cnan
      character*50 result
      character*5 ttt
c
      call CTime(result,time8())
      ttt=result(11:16)
c
      return
      end 






