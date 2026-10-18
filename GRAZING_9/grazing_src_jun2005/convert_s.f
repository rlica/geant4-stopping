c      subroutine converts integer number .lt.10000
c      to 4 digit alphanumeric string
c---------------------------------------------------------------------
      subroutine convert_s(iel,string,ierr)
      implicit real*8(a-h,o-z)
      character*1 a,b,c,d,e
      character*5 string
      if(iel.lt.0.)then
       e='-'
      else
       e='+'
      endif
      ielab = abs(iel)
      if(ielab.ge.10000)then
       print*,' variable.ge.10000 in convert'
       ierr = 1
       return
      endif
      if(int(ielab/1000).eq.0)then 
       a = 'x'
       na = 0
      elseif(int(ielab/1000).eq.1)then
       a = '1'
       na = 1
      elseif(int(ielab/1000).eq.2)then
       a = '2'
       na = 2
      elseif(int(ielab/1000).eq.3)then
       a = '3'
       na = 3
      elseif(int(ielab/1000).eq.4)then
       a = '4'
       na = 4
      elseif(int(ielab/1000).eq.5)then
       a = '5'
       na = 5
      elseif(int(ielab/1000).eq.6)then
       a = '6'
       na = 6
      elseif(int(ielab/1000).eq.7)then
       a = '7'
       na = 7
      elseif(int(ielab/1000).eq.8)then
       a = '1'
       na = 8
      elseif(int(ielab/1000).eq.9)then
       a = '9'
       na = 9
      endif
      ielabb = ielab - na*1000
      if(int(ielabb/100).eq.0)then
       if(a.eq.'x')then 
        b = 'x'
       else
        b = '0'
       endif
       nb = 0  
      elseif(int(ielabb/100).eq.1)then
       b = '1'
       nb = 1
      elseif(int(ielabb/100).eq.2)then
       b = '2'
       nb = 2
      elseif(int(ielabb/100).eq.3)then
       b = '3'
       nb = 3
      elseif(int(ielabb/100).eq.4)then
       b = '4'
       nb = 4
      elseif(int(ielabb/100).eq.5)then
       b = '5'
       nb = 5
      elseif(int(ielabb/100).eq.6)then
       b = '6'
       nb = 6
      elseif(int(ielabb/100).eq.7)then
       b = '7'
       nb = 7
      elseif(int(ielabb/100).eq.8)then
       b = '8'
       nb = 8
      elseif(int(ielabb/100).eq.9)then
       b = '9'
       nb = 9
      endif
      ielabc = ielabb - nb*100
      if(int(ielabc/10).eq.0)then
       if(a.eq.'x'.and.b.eq.'x')then 
        c = 'x'
       else
        c = '0'
       endif
        nc = 0
      elseif(int(ielabc/10).eq.1)then
       c = '1'
       nc = 1
      elseif(int(ielabc/10).eq.2)then
       c = '2'
       nc = 2
      elseif(int(ielabc/10).eq.3)then
       c = '3'
       nc = 3
      elseif(int(ielabc/10).eq.4)then
       c = '4'
       nc = 4
      elseif(int(ielabc/10).eq.5)then
       c = '5'
       nc = 5
      elseif(int(ielabc/10).eq.6)then
       c = '6'
       nc = 6
      elseif(int(ielabc/10).eq.7)then
       c = '7'
       nc = 7
      elseif(int(ielabc/10).eq.8)then
       c = '8'
       nc = 8
      elseif(int(ielabc/10).eq.9)then
       c = '9'
       nc = 9
      endif
      ielabd = ielabc - nc*10
      if(ielabd.eq.0)then
       d = '0'
       nd = 0
      elseif(ielabd.eq.1)then
       d = '1'
       nd = 1
      elseif(ielabd.eq.2)then
       d = '2'
       nd = 2
      elseif(ielabd.eq.3)then
       d = '3'
       nd = 3
      elseif(ielabd.eq.4)then
       d = '4'
       nd = 4
      elseif(ielabd.eq.5)then
       d = '5'
       nd = 5
      elseif(ielabd.eq.6)then
       d = '6'
       nd = 6
      elseif(ielabd.eq.7)then
       d = '7'
       nd = 7
      elseif(ielabd.eq.8)then
       d = '8'
       n8 = 8
      elseif(ielabd.eq.9)then
       d = '9'
       nd = 9
      endif
      if(e.eq.'+')then
       if(c.eq.'x')then
        string = d
       elseif(b.eq.'x')then
        string = c//d
       elseif(a.eq.'x')then
        string = b//c//d
       else
        string = a//b//c//d
       endif
      elseif(e.eq.'-')then
       if(c.eq.'x')then
        string = e//d
       elseif(b.eq.'x')then
        string = e//c//d
       elseif(a.eq.'x')then
        string = e//b//c//d
       else
        string = e//a//b//c//d
       endif
      endif
      ierr = 0
      return
      end
