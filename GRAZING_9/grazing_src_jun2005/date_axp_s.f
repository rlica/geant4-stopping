c      Subroutine delivers date in the format DD-MMM-YYYY
c-------------------------------------------
      subroutine date_s(ddd)
      implicit real*8(a-h,o-z)
      character*8  date
      character*10 time
      character*5  zone 
      integer*4 value(8)
      character*11 ddd
      character*50 ddd1
      character*5 dd
      character*3 mm(12)
      character*3 mmm
      character*5 yyyy
      data mm /
     *  'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct'
     1 ,'Nov','Dec'/
      call date_and_time(date,time,zone,value)
      idd = value(3)
      call convert_s(idd,dd,ierr)
      imm = value(2)
      mmm = mm(imm)
      iyy = value(1)
      call convert_s(iyy,yyyy,ierr)
      ddd1 = dd//'-'//mmm//'-'//yyyy
      call condens_s(ddd1,ddd,len)
c
      return
      end 
