c      Subroutine delivers date in the format DD-MMM-YYYY
c-------------------------------------------
      subroutine date_s(ddd)
      character*50 result
      character*11 ddd
cnan      implicit real*8(a-h,o-z)
cnan      integer*4 value(8)
cnan      character*11 ddd
cnan      character*50 ddd1
cnan      character*5 dd
cnan      character*3 mm(12)
cnan      character*3 mmm
cnan      character*5 yyyy
cnan      data mm /
cnan     *  'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct'
cnan     1 ,'Nov','Dec'/
cnan      call date_and_time(,,,value)
cnan      idd = value(3)
cnan      call convert_s(idd,dd,ierr)
cnan      imm = value(2)
cnan      mmm = mm(imm)
cnan      iyy = value(1)
cnan      call convert_s(iyy,yyyy,ierr)
cnan      ddd1 = dd//'-'//mmm//'-'//yyyy
cnan      call condens_s(ddd1,ddd,len)
cnan
cnan  ---------------------------------------------------------------
cnan  NB:
cnan  un-comment all the  above lines if you un on alpha and comment
cnan  all the lines below that are for linux
cnan  ---------------------------------------------------------------
cnan
      character*50 resultcnan
cnan  ---------------------------------------------------------------
cnan  NB:
cnan  un-comment all the  above lines if you un on alpha and comment
cnan  all the lines below that are for linux
cnan  ---------------------------------------------------------------
cnan
c
      call CTime(result,time8())
      ddd=result(5:7)//'-'//result(9:10)//'-'//result(21:24)
c
      return
      end 
