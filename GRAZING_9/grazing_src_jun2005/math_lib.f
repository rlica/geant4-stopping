c      Library contains following mathematical functions:
c      sqrt_c,gncmplx,ein2,gamma,aloggamma,gammainc 
c      and subroutines:
c      ierfc,gn_s,gn_cs,ibessel_s,fourier_s,prob_s,probd_s,probdtab_s
c      besseldist_s,besseldis_cs,besseldist_cs,zero_s,interpolate_s
c      to be used in program grazing.
c      Uses nag_shr/opt subroutines
c      All in double precision and no common blocks
c----------------------------------------------------------------
c      Function calculates squareroot of complex number with
c      branchcut along (just below) real axis.
c----------------------------------------------------------------
      complex*16 function sqrt_c(z)
      implicit real*8(a-h,o-z)
      complex*16 z, sqrt
      sqrt_c = sqrt(z)
      if(dimag(z).lt.0.d0)sqrt_c = -sqrt_c
      return
      end
c---------------------------------------------------------------
c      Subroutine for multiple integrals of error function
c      Uses S15ADF from nag_shr/opt
c_____________________________________________________________________________

      subroutine ierfc(n,x,func)
      implicit real*8(a-h,o-z)
      ifail = 0
      gh = S15ADF(x,ifail)
      gl = 2.d0/sqrt(3.14159265358979d0)*exp(-x**2)
      if (n.eq.0) goto 2
      nn = 1
1     gll = gh
      gh = (-x*gh + .5d0*gl)/float(nn)
      gl=gll
      if(nn.eq.n) goto 2
      nn=nn+1        
      goto 1
2     func = gh
      end      
c---------------------------------------------------------------------------
c      Function for spherical harmonic Y_lm(pi/2,0)
c--------------------------------------------------------------

      real*8 function ylm(l,m,ifail)
      implicit real*8(a-h,o-z)
      if(abs(m).gt.l)then
       ifail = 2
       print*,' abs(m).gt.l in ylm'
       ylm = 0.d0
       return
      endif
      if(mod(l+m,2).ne.0)then
       ylm = 0.d0
       return
      endif
      iarg1 = l+m
      iarg2 = l-m
      ylm = sqrt((2.*dble(l)+1.)/4.d0/3.141592653589793238d0)
     &*facrat(iarg1)*facrat(iarg2)*(-1)**(iarg1/2)
      return
      end

      real*8 function facrat(n)
      implicit real*8(a-h,o-z)
      facrat = 1.d0
      nn = 0
1     if(nn.ge.n)goto 2
      nn = nn + 2
      facrat = facrat*float(nn-1)/float(nn)
      goto 1
2     facrat = sqrt(facrat)
      return
      end
c----------------------------------------------------------------------
c      function calculates G_n(z) of complex argument in double prec.
c      accuracy decreasing for large n=no.
c      Overflow for real(z)<-7 close to real axis
c      At entry: ierr=0 (standard), ierr=3 (Taylor),
c      ierr=4 (continued fraction), ierr=2 (Assymptotic exp.)
c      At exit ierr = 2 if no convergence in assymptotic expansion
c      ierr = 1 if overflow, ierr = 3 if no convergence in Taylor exp.
c      ierr = 4 if no convergence in continued fraction
c__________________________________________________________________

      complex*16 function gncmplx(xx,yy,no,accur,ierr)
      implicit real*8(a-h,o-z)
      complex*16 z,sn,sum,z2,cw,cde,gh,gl,gll,fab
      complex*16 aaa,bbb,aanm,aan,bbnm,bbn,fn,fnm
      ierr2 = 0
      z=dcmplx(xx,yy)
      abz = abs(z)
      z2=z**2
      nrz = 1    
      if(dreal(z).ge.0.d0)goto 101
      nrz = -1
      cde=cdexp(z2)
      z = - z      
101   if(ierr.eq.0)then
       if(abz.le.3.d0)goto 240
       if(abz.gt.10.d0)goto 222      
       goto 230
      endif
      if(ierr.eq.2)goto 222
      if(ierr.eq.3)goto 240
      if(ierr.eq.4)goto 230
      print*,' ierr.eq.1.or.gt.3 in gn_c'
      ierr = 1
      goto 3
      
c      Taylor expansion
240   cde = cdexp(z2)
      n=0
      sum= -2.d0/1.7724538509055160272d0*z
      sn = sum
100   n=n+1
       ab = abs(sn)
       if(ab.gt.1.d28)then
        ierr2 = 3
        goto 230
       endif
       sn = 2.d0*z2/(2.d0*n+1.d0)*sn
       sum = sum + sn
       if(n.gt.300)then
        ierr2 = 3
        goto 230
       endif
      if(ab.ge.accur.and.ab.ge.abs(cde+sum)*accur)goto 100
      cw = cde + sum
      goto 300

c      Assymptotic expansion
222   sum = 1.d0
      sn = 1.d0
      m = 0
223   m=m+1
       fab = (2.d0*m-1.d0)/(-2.d0*z2)
       if(abs(fab).gt.1.d0)then
        if(ierr2.eq.3)then
         print*,' No convergence in gn_c'
         ierr = 1
         goto 3
         else
         ierr2 = 2
         goto 230
        endif
       endif
       sn=sn*fab
       if(abs(sn).lt.accur*abs(sum))goto 224
       sum = sum + sn
      goto 223
224   cw = sum/1.7724538509055160272d0/z
      goto 300

c      Continued fraction
230   n=0
      aanm = 1.d0
      bbnm = 0.d0
      aan = 0.d+0
      bbn = 1.d+0
      fnm = 0.d+0
231   n=n+1
       aaa = aan*z + aanm*n/2.d0
       bbb = bbn*z + bbnm*n/2.d0
       if(abs(aan)+ abs(bbn).gt.0.1d+30)then
c       correction if realz.lt.0
        if(nrz.eq.-1)z = - z
        ierr2 = 4
        goto 240
       endif
       if(abs(bbb).le.1.d-30)then
        ierr2 = 4
        goto 240
       endif
       fn = aaa/bbb
       if(cdabs(fn-fnm).le.accur*cdabs(z+fn))goto 232
       if(n.ge.100)then
        ierr2 = 4
        goto 240
       endif
       aanm = aan
       bbnm = bbn
       aan = aaa
       bbn = bbb
       fnm = fn
      goto 231
232   cw = 1.d0/1.7724538509055160272d0/(z + fn)      
c      correction if realz.lt.0
300   if(nrz.eq.-1)then
       cw = 2.d0*cde - cw
       z = - z
      endif
      gh = cw
      gl = 2.d0/1.7724538509055160272d0 
      if (no.eq.0) goto 2

c      Recursion for n.ge.1
      nn = 1
1     gll = gh
      gh = (-z*gh + .5d0*gl)/nn
      gl=gll
      if(nn.eq.no) goto 2
      nn=nn+1        
      goto 1
2     gncmplx = gh
      ierr = 0
      if(ierr2.ne.0)ierr = ierr2
      return
3     gncmplx = dcmplx(0.d0,0.d0)
      return
      end 
c----------------------------------------------------------------------
c      subroutine calculates G_n(z) of real argument in double prec.
c      G_0....G_no are stored in gn(1)...gn(no+1)
c      with no_max = 10. Accuracy decreasing for large no.
c      Overflow for real(z)<-7 close to real axis
c      At entry: ierr=0 (standard), ierr=3 (Taylor),
c      ierr=4 (continued fraction), ierr=2 (Assymptotic exp.)
c      At exit ierr = 2 if no convergence in assymptotic expansion
c      ierr = 1 if overflow, ierr = 3 if no convergence in Taylor exp.
c      ierr = 4 if no convergence in continued fraction
c__________________________________________________________________

      subroutine gn_s(z,nomax,gn,accur,ierr)
      implicit real*8(a-h,o-z)
      dimension gn(11)
      if(nomax.lt.0.or.nomax.gt.10)then
       print*,' wrong input in gn_s'
       ierr = 1
       return
      endif
      ierr2 = 0
      abz = abs(z)
      z2=z**2
      nrz = 1    
      if(z.ge.0.d0)goto 101
      nrz = -1
      if(z2.gt.75.d0)then 
       ierr = 1                     
       print*,'overflow in gn_s'
       return
      else
       cde=exp(z2)
      endif
      z = - z      
101   if(ierr.eq.0)then
       if(abz.le.3.d0)goto 240
       if(abz.gt.10.d0)goto 222      
       goto 230
      endif
      if(ierr.eq.2)goto 222
      if(ierr.eq.3)goto 240
      if(ierr.eq.4)goto 230
       print*,' ierr.eq.1.or.gt.3 in input in gn_s'
      ierr = 1
      goto 3
      
c      Taylor expansion
240   cde = exp(z2)
      n=0
      sum= -2.d0/1.7724538509055160272d0*z
      sn = sum
100   n=n+1
       ab = abs(sn)
       if(ab.gt.1.d+28)then
        ierr2 = 3
        goto 230
       endif
       sn = 2.d0*z2/(2.d0*n+1.d0)*sn
       sum = sum + sn
       if(n.gt.300)then
        ierr2 = 3
        goto 230
       endif
      if(ab.ge.accur.and.ab.ge.abs(cde+sum)*accur)goto 100
      cw = cde + sum
      goto 300

c      Assymptotic expansion
222   sum = 1.d0
      sn = 1.d0
      m = 0
223   m=m+1
       fab = (2.d0*m-1.d0)/(-2.d0*z2)
       if(abs(fab).gt.1.d0)then
        if(ierr2.eq.3)then
         print*,' No convergence in gn_s'
         ierr = 1
         goto 3
         else
         ierr2 = 2
         goto 230
        endif
       endif
       sn=sn*fab
       if(abs(sn).lt.accur*abs(sum))goto 224
       sum = sum + sn
      goto 223
224   cw = sum/1.7724538509055160272d0/z
      goto 300

c      Continued fraction
230   n=0
      aanm = 1.d0
      bbnm = 0.d0
      aan = 0.d+0
      bbn = 1.d0
      fnm = 0.d0
231   n=n+1
       aaa = aan*z + aanm*n/2.d0
       bbb = bbn*z + bbnm*n/2.d0
       if(abs(aan)+ abs(bbn).gt.0.1d+30)then
c       correction if realz.lt.0
        if(nrz.eq.-1)z = - z
        ierr2 = 4
        goto 240
       endif
       if(abs(bbb).le.1.d-30)then
        ierr2 = 4
        goto 240
       endif
       fn = aaa/bbb
       if(abs(fn-fnm).le.accur*abs(z+fn))goto 232
       if(n.ge.100)then
        ierr2 = 4
        goto 240
       endif
       aanm = aan
       bbnm = bbn
       aan = aaa
       bbn = bbb
       fnm = fn
      goto 231
232   cw = 1.d0/1.7724538509055160272d0/(z + fn)      
c      correction if realz.lt.0
300   if(nrz.eq.-1)then
       cw = 2.d0*cde - cw
       z = - z
      endif
      gh = cw
      gl = 2.d0/1.7724538509055160272d0 
      gn(1) = gh
      if(nomax.eq.0)goto 4
c
c      Recursion for n.ge.1
      n = 0
1     n=n+1        
      gll = gh
      gh = (-z*gh + .5d0*gl)/n
      gl=gll
2     gn(n+1) = gh
      if(n.lt.nomax)goto 1
4     ierr = 0
      if(ierr2.ne.0)ierr = ierr2
3     return
      end       
c----------------------------------------------------------------------
c      subroutine calculates G_n(z) of complex argument in double prec.
c      G_0....G_no are stored in func(1)...func(no+1)
c      with no_max = 10. Accuracy decreasing for large n.
c      Overflow for real(z)<-7 close to real axis
c      At entry: ierr=0 (standard), ierr=1 (Taylor),
c      ierr=2 (continued fraction), ierr=3 (Assymptotic exp.)
c      At exit ierr = 2 if no convergence in assymptotic expansion
c      ierr = 1 if overflow, ierr = 3 if no convergence in Taylor exp.
c      ierr = 4 if no convergence in continued fraction
c__________________________________________________________________

      subroutine gn_cs(xx,yy,no,func,accur,ierr) 
      implicit real*8(a-h,o-z)
      complex*16 z,sn,sum,z2,cw,cde,gh,gl,gll,func(11),fab
      complex*16 aaa,bbb,aanm,aan,bbnm,bbn,fn,fnm 
      do nn = 1,11
       func(nn) = dcmplx(0.d0,0.d0)
      enddo
      ierr2 = 0
      z=dcmplx(xx,yy)
      abz = abs(z)
      z2=z**2
      nrz = 1    
      if(dreal(z).ge.0.)goto 101
      nrz = -1
      cde=cdexp(z2)
      z = - z      
101   if(ierr.eq.0)then
       if(abz.le.3.)goto 240
       if(abz.gt.10.d0*dble(no+1))goto 222      
       goto 230
      endif
      if(ierr.eq.1)then
       if(abz.ge.7.d0)then
        ierr = 3
        goto 3 
        else
       endif
       goto 240
      endif
      if(ierr.eq.2)then
       if(abz.le.1.d-3)then
        ierr = 4
        goto 3 
       else
       endif
       goto 230
      endif
      if(ierr.eq.3)then
       if(abz.le.1.d-3)then
        ierr = 2
        goto 3 
       else
       endif
       goto 222
      endif
      print*,' ierr.gt.3 in gn_cs'
      ierr = 1
      goto 3
      
c      Taylor expansion
240   cde = cdexp(z2)
      n=0
      sum= -2.d0/1.7724538509055160272d0*z
      sn = sum
100   n=n+1
      ab = abs(sn)
      if(ab.gt.1.d+28)then
       ierr = 3
       goto 3
      endif
      sn = 2.d0*z2/(2.d0*dble(n)+1.d0)*sn
      sum = sum + sn
c      print*,'sn =',sn
      if(n.gt.300)then
       print*,' No convergence in Taylor expansion of gn_cs'
       ierr = 3
       goto 3
      endif
      if(ab.ge.accur.and.ab.ge.abs(cde+sum)*accur)goto 100
      cw = cde + sum
      goto 300

c      Assymptotic expansion
222   sum = 1.d0
      sn = 1.d0
      m = 0
223   m=m+1
      fab = (2.d0*dble(m)-1.d0)/(-2.d0*z2)
      if(abs(fab).gt.1.d0)then
       abserr = abs(sn/z)
       write(*,*)
     &' Accuracy in assymptotic expansion insufficient: abserr =',
     &abserr
       ierr2 = 2
       goto 230
      endif
      sn=sn*fab
      if(abs(sn).lt.accur*abs(sum))goto 224
      sum = sum + sn
      goto 223
224   cw = sum/1.7724538509055160272d0/z
      goto 300

c      Continued fraction
230   n=0
      aanm = 1.d0
      bbnm = 0.d0
      aan = 0.d0
      bbn = 1.d0
      fnm = 0.d0
231   n=n+1
      aaa = aan*z + aanm*dble(n)/2.d0
      bbb = bbn*z + bbnm*n/2.d0
      if(abs(aan)+ abs(bbn).gt.0.1d+30)then
c      correction if realz.lt.0
       if(nrz.eq.-1)z = - z
       ierr2 = 4
       goto 240
      endif
      if(abs(bbb).le.1.d-30)then
       ierr2 = 4
       goto 240
      endif
      fn = aaa/bbb
      if(cdabs(fn-fnm).le.accur*cdabs(z+fn))goto 232
c      print*,'fn=',fn
      if(n.ge.100)then
       ierr2 = 4
       goto 240
      endif
      aanm = aan
      bbnm = bbn
      aan = aaa
      bbn = bbb
      fnm = fn
      goto 231
232   cw = 1.d0/1.7724538509055160272d0/(z + fn)      
c      correction if realz.lt.0
300   if(nrz.eq.-1)then
       cw = 2.d0*cde - cw
       z = - z
      endif
      gh = cw
      gl = 2.d0/1.7724538509055160272d0
      func(1) = gh 
      if (no.eq.0) goto 2

c      Recursion for n.ge.1
      nn = 1
1     gll = gh
      gh = (-z*gh + .5d0*gl)/nn
      gl=gll
      nn=nn+1        
      func(nn) = gh
      if(nn.eq.no+1) goto 2
      goto 1
2     ierr = 0
      if(ierr2.ne.0)ierr = ierr2
3     return
      end 
c--------------------------------------------------------------------
c      Subroutine for modified Bessel function I_n(x)
c      uses subroutines S18AEF,S18AFF from nag_shr/opt
c_______________________________________________________________________
c
      subroutine ibessel_s(aibess,n,x,ifa)
      implicit real*8(a-h,o-z)
      if(n.lt.0)n=-n
      ifail = 0
      ai0 = S18AEF(x,ifail)
      if(ifail.ne.0)goto 4
      ai1 = S18AFF(x,ifail)
      if(ifail.ne.0)goto 4
      if(n.eq.0)aibess = ai0
      if(n.eq.1)aibess = ai1
      if(n.lt.2)goto 1
      if(x**2/n.lt.1.)goto 5
      nn=2
3     ainn = ai0 - 2.d0*(dble(nn)-1.d0)/x*ai1
      if(nn.eq.n)goto 2
      ai0 = ai1
      ai1 = ainn
      nn=nn+1
      goto 3
2     aibess = ainn
1     return
5     fact = 1.d0
      cc=x**2/4.d0
      i=1
      ci = 1.d0
6     ci = ci*cc/dble(n+i)/dble(i)           
      if(ci.lt.1.d-10)goto 7
      fact = fact + ci
      i = i + 1
      goto 6
7     nfac = 1
      do 8, j=1,n
8     nfac = nfac*j
      aibess = (x/2.d0)**n/dble(nfac)*fact 
      return
4     ifa = 1
      print*,'problems with ibessel'
      end                  
      
c---------------------------------------------------------------------
c      Calculates function Ein2(x)=int_0^x{(exp(-t)-1+t)/t**2}dt
c
c----------------------------------------------------------------

      real*8 function ein2(x,ierr)
      implicit real*8(a-h,o-z)
      if(x.gt.15)goto 4
      if(x.lt.-25)goto 3
      ein2 = x/2.d0
      n=1
      term = x/2.d0
1     term = -term*x*dble(n/(n+1)/(n+2))
      ein2 = ein2 + term
      if(abs(term).lt.0.0000001d0)goto 2
      n=n+1
      goto 1
2     return
4     ein2 = log(x)-1.d0+1.d0/x + .5772156649d0
      return
3     ierr = 1
      print*,'Too large argument in ein2'
      ein2 = 0.d0
      return
      end      
c--------------------------------------------------------------------
c      Subroutine calculates Fourier transform 
c      1/(2 pi)\int_-\infty^infty\exp(-i omega y)G(x+iy)dy
c      The function G should be defined through functions 
c      gre(y)=real(G) and gim=imag(G)
c      Uses D01ASF from nag_shr/opt
c---------------------------------------------------------------

      subroutine fourier_s(result,gre,gim,omega,accur,ierr)
      implicit real*8(a-h,o-z)
      external gre,gim
      dimension ierlst(50),iw(1000),erlst(50),rslst(50),w(2000)

      key = 1
      ifail=-1
      acc = accur
      call D01ASF(gre,0.d0,omega,key,acc,res1,abserr,50,lst,
     &  erlst,rslst,ierlst,w,2000,iw,1000,ifail)
      if(ifail.ne.0)then
       print*,'warning cosine fourier abserr =',abserr
       ierr = 1
       result = 0.d0
       return
      endif
      key = 2
      ifail=-1
      acc = accur
      call D01ASF(gim,0.d0,omega,key,acc,res2,abserr,50,lst,
     &  erlst,rslst,ierlst,w,2000,iw,1000,ifail)
      if(ifail.ne.0)then
       print*,'warning sine fourier abserr =',abserr
       ierr = 1
       result = 0.d0
       return
      endif
      result = 1.d0/3.1415926535897932384d0*(res1+res2)
      ierr = 0
      return
      end
c-------------------------------------------------------------------
c      Function calculates gamma function of real argument
c      uses S14AAF from nag_shr/opt
c----------------------------------------------------------------------
      real*8 function gamma(x,ierr)
      implicit real*8(a-h,o-z)
      ierr = -1
      gamma = S14AAF(x,ierr)
      if(ierr.eq.1)then
       print*,' Too large argument in gamma'
       goto 2
       elseif(ierr.eq.2)then
       print*,' Too large negative argument in gamma'
       goto 2
       elseif(ierr.eq.3)then
       print*,' Argument in gamma too close to 0'
       goto 2
       elseif(ierr.eq.4)then
       print*,' Argument in gamma is negative integer'
       goto 2
      endif
      ierr  = 0
      return
2     ierr = 1
      return
      end
c---------------------------------------------------------------
c      Function calculates logarithm of gamma function 
c      of real argument x>0
c      uses S14ABF from nag_shr/opt
c----------------------------------------------------------------------
      real*8 function aloggamma(x,ierr)
      implicit real*8(a-h,o-z)
      if(x.le.0.)goto 2
      aloggamma = S14ABF(x,ierr)
      if(ierr.eq.2)goto 22
      ierr = 0
      return
2     print*,' Error: argument in loggamma .le. 0'
      aloggamma = 0.d0
22    ierr = 1
      return
      end      
c---------------------------------------------------------------------
c      Function calculates incomplete gamma function
c      gammmainc(a,x)=int_x^infty{exp(-t)*t**(a-1)dt}
c      uses D01AMF from nag_shr/opt library
c---------------------------------------------------------------

      real*8 function gammainc(a,x,ierr)
      implicit real*8(a-h,o-z)
      external funna
      common/locgam/aa
      dimension w(800),iw(200)
      aa=a
      if(aa.gt.33.d0)goto 2
      accur0 = 1.d-9 
      accur = accur0
      inf = 1
      ifail = -1
      lw = 800
      liw = 200
      call D01AMF(funna,x,inf,accur,accur,gam,abserr,w
     & ,lw,iw,liw,ifail)
      if(ifail.ne.0)goto 1
      gammainc = gam
      return
2     print*,'too large argument in gammaincompl.'
1     ierr = 1
      gammainc = 0.d0
      return
      end

      real*8 function funna(t)
      implicit real*8(a-h,o-z)
      common/locgam/aa
      if(t.gt.1.d0)goto 3
      funna = exp(-t)*t**(aa-1.d0)
      return
3     tt=-t+(aa-1.d0)*log(t)
      funna = exp(tt)
      return
      end      
c----------------------------------------------------------------
c      Function calculates incomplete gamma function
c      gammmainc(a,x)=int_x^infty{exp(-t)*t**(a-1)dt}
c      divided by gamma(a) i.e. the chi-distribution
c      uses D01AMF from nag_shr/opt library
c---------------------------------------------------------------

      real*8 function chidist(a,x,ierr)
      implicit real*8(a-h,o-z)
      external funnb
      common/locchi/aa,alof
      dimension w(800),iw(200)
      aa=a
      ierr = 0
      alof = aloggamma(a,ierr)
      if(ierr.eq.1)then
       chidist = 0.d0
       return
      endif
      accur0 = 1.d-9 
      accur = accur0
      inf = 1
      ifail = -1
      lw = 800
      liw = 200
      call D01AMF(funnb,x,inf,accur,accur,chidist,abserr,w
     & ,lw,iw,liw,ifail)
      if(ifail.ne.0)goto 1
      ierr = 0
      return
1     ierr = 1
      return
      end

      real*8 function funnb(t)
      implicit real*8(a-h,o-z)
      common/locchi/aa,alof
      tt=-t+(aa-1.d0)*log(t)
      funnb = exp(tt-alof)
      return
      end      
c-------------------------------------------------------------     
c      Subroutine calculates probability with gammadistribution
c      for given average and spread to have the varible
c      larger than a given amount.
c      Uses gammainc,gamma i.e. D01AMF and S14AAF from nag_shr/opt 
c----------------------------------------------------------------
c
      subroutine prob_s(eav,sig2,e0,prob,ierr)
      implicit real*8(a-h,o-z)
      if(e0.le.0.)then
       prob = 1.d+0
       goto 3
      endif
      gamnu = (eav)**2/sig2
      gamalf = eav/sig2
      gamarg = e0*gamalf
      prob = chidist(gamnu,gamarg,ierr)
      if(ierr.eq.1)goto 2
3     ierr = 0
      return
2     print*,'Error in prob_s'
      return
      end        
c---------------------------------------------------------------
c      Subroutine calculates probability density with 
c      gammadistribution for given average and spread 
c      Uses gammainc,gamma i.e. D01AMF and S14AAF from nag_shr/opt 
c      At exit ierr = 2 if gamnu<1 but probd = 0 (actually diverging)
c----------------------------------------------------------------
c
      subroutine probd_s(eav,sig2,e0,probd,ierr)
      implicit real*8(a-h,o-z)
      if(e0.le.0.)then
       probd = 0.d+0
       goto 3
      endif
      gamnu = (eav)**2/sig2
      gamalf = eav/sig2
      gamarg = e0*gamalf
      ierr1 = 0
      alogcoeff = log(gamalf)*gamnu-aloggamma(gamnu,ierr1)
      if(ierr1.eq.1)goto 2
      prolog = (gamnu-1.d0)*log(e0)-gamarg +alogcoeff
      probd = exp(prolog)
3     ierr = 0
      if(gamnu.lt.1.d0)ierr = 2
      return
2     ierr = 1
      return
      end        
c------------------------------------------------------------------------
c      Subroutine calculates table of probability density with 
c      gammadistribution for given energystep, average and spread 
c      Prob(E=ein+i*dele) > accur/eav are stored in gammad(i) (i=1,100)
c      Uses gammainc,gamma i.e. D01AMF and S14AAF from nag_shr/opt 
c      At exit ierr = 2 if gamnu < 1 (probability density diverges)
c              ierr = 3 if dimension too small
c----------------------------------------------------------------
c
      subroutine probdtab_s(eav,sig2,ein,dele,gammad,accur,ierr)
      implicit real*8(a-h,o-z)
      dimension gammad(100)
      if(abs(eav).lt.accur.or.sig2.lt.accur.or.dele.lt.accur.
     &or.accur.le.0.d0)then
       write(*,'(a,/a,4e18.11)')' wrong input in probdtab_s ',
     &' eav,sig2,dele,accur = ',eav,sig2,dele,accur
       goto 2
      endif
      do i=1,100
       gammad(i)=0.d0
      enddo
      acc = accur/eav
      gamnu = (eav)**2/sig2
      if(gamnu.lt.1.d0)iwarn = 1
      gamalf = eav/sig2
      ierr1 = 0
      alogcoeff = log(gamalf)*gamnu-aloggamma(gamnu,ierr1)
      if(ierr1.eq.1)goto 2
      e0 = ein
      n = 0
1     e0 = e0+dele
      n = n + 1
      if(e0.lt.0.d0)then
       gammad(n) = 0.d0
       goto 1
      endif  
      if(n.gt.100)then
       write(*,*)' Number of steps in probdtab_s exceeds dimension'
       write(*,*)'   dele = ',dele,' eav = ',eav,' sig2 = ',sig2
       ierr = 3
       return
      endif
      gamarg = e0*gamalf
      prolog = (gamnu-1.d0)*log(e0)-gamarg +alogcoeff
      if(prolog.gt.72.d0)then
       write(*,*)' probability diverges in probdtab_s'
       goto 2
      elseif(prolog.lt.-72.d0)then
       gammad(n) = 0.d0
      else 
       gammad(n) = exp(prolog)
      endif
      if(e0.lt.eav.or.gammad(n).gt.acc)goto 1
3     ierr = 0
      if(iwarn.eq.1)ierr = 2
      return
2     ierr = 1
      return
      end        
c------------------------------------------------------------------------
c      Subroutine calculates the Bessel distribution
c      within given accuracy limit. 
c      At exit nmin and nmax indicate the limits of n calculated.
c      Result stored in besdis array, with n=0 at besdis((ibesrange+1)/2),
c      Maximum |nmin| and nmax given by ibesrange to be 100
c---------------------------------------------------------------------
      subroutine besseldist_s(x,y,nmin,nmax,besdis,accur,ierr)
      implicit real*8 (a-h,o-z)
      parameter (nrange=100,ibesrange=2*nrange+1,nn0=nrange+1)
      dimension besdis(ibesrange),besdisx(ibesrange)
      acc = accur/10.d0
      if(x+y.gt.70.d0)then
       write(*,*)' x and y out of range in besdist_s'
       ierr = 1
       return
      endif
      deln = .3d0*(-log(accur))*sqrt(x+y)
      if(x+y.lt..1)deln = log(accur)/log(x+accur)
      nstart = x-y+ deln
      if(nstart.ge.nrange)then
       nstart = nrange-1
      endif
      do i=1,ibesrange
       besdisx(i) = 0.d0
       besdis(i) = 0.d0
      enddo
      if(nstart.lt.0.and.x-y.lt.0.d0)then
       nstart = 0
      endif
      xe = x
      ye = y
      ne = nstart
      nmaxe = nrange
      nmine = -nrange
      nistop = nn0
c  calculate first from expansion for n=ne and n=ne+1
      if(ne.gt.nrange-1)ne = nrange-1
3     sk = 1.d0
      do i = 1,ne
       sk = sk*xe/dble(i)
      enddo
      skp = sk*xe/dble(ne+1)
      if(skp.eq.0.d0)then
       ne = ne - .2*ne
       goto 3
      endif
      sum = sk
      sump = skp
      do k = 1,400
       ss = ye/dble(k)*xe
       sk = sk*ss/dble(k+ne)
       skp = skp*ss/dble(k+ne+1)
       sum = sum+sk
       sump = sump+skp 
       if(sk.le.accur*sum.and.skp.le.accur*sump)goto 4
      enddo
      write(*,*)' No convergence in besseldist_s'
c      print*,' sum,sump,sk,skp = ',sum,sump,sk,skp
      ierr = 1
      return
4     eex =exp(-xe-ye)
      besdisx(ne+nn0) = sum*eex
      besdisx(ne+nn0+1) = sump*eex
c  Now backwards recursion to n=-1
      do ni= ne+nn0,nistop,-1
       besdisxx=(ye*besdisx(ni+1)+dble(ni-nn0)*besdisx(ni))/xe
       besdisx(ni-1)=besdisxx
       nn = ni-nn0-1
       if(besdisx(ni+1).lt.acc.and.ni-nn0+1.gt.xe-ye)then
        nmaxe = ni-nn0+1
       endif
      enddo
c if nstart was not large enough:
      if(nmaxe.eq.nrange)then
       if(ne.eq.nrange-1)then
        write(*,*)' range in besdist_s is too small'
        ierr = 2
        goto 6
       endif
       nistop = ne+1+nn0
c  new nstart
       delne =.1d0*(-log(accur))*sqrt(xe+ye)
       if(delne.lt.0.2d0)then
        nmaxe = ne+1
        goto 6
       elseif(delne.lt.2.d0)then 
        delne = 2.001d0
       endif
       ne = dble(ne)+delne
       if(ne.gt.nrange-1)ne=nrange-1
       goto 3
      else
       if(nistop.ne.nn0.
     &and.abs(besdisxx-besdisx(nistop-1)).gt.accur)then 
        write(*,'(2(a,d14.8,/))')
     &' trouble: besdisx = ',besdisx(nistop-1),
     &'         besdisxx = ',besdisxx
        ierr = 1
        return
       endif
      endif
6     if(nmine.eq.-nrange)then
c  store  the positive n values
       nmax = nmaxe
       besmin1 = besdisx(nn0-1)
       bes0 = besdisx(nn0)
       do i=-1,nrange
        besdis(nn0+i)=besdisx(nn0+i)
        besdisx(nn0+i)=0.d0
       enddo
       if(besmin1.lt.accur.and.xe-ye.ge.0.d0)then
        nmin = -1
        ierr = 0
        return
       else
        xe = y
        ye = x
        deln = .3d0*(-log(accur))*sqrt(xe+ye)
        if(xe+ye.lt..1)deln = log(accur)/log(xe)
        ne = xe-ye+ deln
        nmine = -1
        nmaxe = +nrange
        nistop = nn0
        goto 3
       endif
      else
c calculation finished for negative n
       if(abs(besdisx(nn0+1)-besmin1).gt.accur.or.
     &abs(besdisx(nn0)-bes0).gt.accur)then
        write(*,'(4(a,d14.8,/))')
     &' trouble: besdisx_(-1) = ',besdisx(nn0+1),
     &'           besdis_(-1) = ',besmin1,
     &'             besdisx_0 = ',besdisx(nn0),
     &'              besdis_0 = ',bes0
        ierr = 1
        return
       endif
       nmin = -nmaxe
c store the negative n values
       do i=1,nrange
        besdis(nn0-i)=besdisx(nn0+i)
       enddo
       ierr = 0
       return
      endif
      end
c--------------------------------------------------------------------
c      Complex subroutine calculates the Bessel distribution
c      for complex argument for given n = no. 
c      Uses S18DEF from nag_shr/opt 
c      zs,zp and besdis must be declared complex*16 in
c      calling program even if zs and zp are real
c------------------------------------------------------------------
      subroutine besseldis_cs(zp,zs,no,besdis,ierr)
      implicit real*8(a-h,o-z)
      complex*16 z,zs,zp,arg,rat,szs,szp,besdis,cy
      complex*16 fff,ffn,z2,expo
      dimension cy(1)
      character*1 scale
c
      scale = 's'
      nno = abs(no)
      fnu = dble(nno)
      szs = sqrt(zs)
      szp = sqrt(zp)
      z = 2.d+0*szs*szp
      if(abs(z).lt.1.d-10)goto 100
      if(10.d0*abs(z)**2.lt.fnu+1)goto 100
      rat = log(szp/szs)
      arg = -zs-zp+abs(dble(z))
      ifail = 0
      n=1
      call S18DEF(fnu,z,n,scale,cy,nz,ifail)
      if(ifail.ge.1)then
       write(*,*)' Error in besseldis_cs, ifail=',ifail,'in S18DEF'
       ierr = 1
       goto 2
      endif
      expo=rat*no+arg
      if(dble(expo).lt.-75.d0)then
       besdis=0.d0
      elseif(dble(expo).gt.75.d0)then
       write(*,*)' Overflow in besseldis_cs'
       ierr = 1
       goto 2
      else
       besdis=exp(expo)*cy(1)
      endif
      ierr = 0
2     return
100   alogfac = 0.d0
      if(nno.eq.0)goto 101
      do n=1,nno
       alogfac = alogfac + log(dble(n))
      enddo
101   ffn = 1.d0
      fff = ffn
      nnn=0
      z2 = zp*zs
102   nnn=nnn+1
      ffn = ffn*z2/dble(nnn)/dble(nno+nnn)
      fff = fff + ffn
      if(abs(ffn).gt.1.d-14)goto 102
      if(abs(zp).lt.1.d-30)zp=1.d-30
      if(abs(zs).lt.1.d-30)zs=1.d-30
      if(no.ge.0)then
       expo=-zp-zs+log(zp)*fnu-alogfac
       if(dble(expo).lt.-50.d0)then
        besdis=0.d0
       elseif(dble(expo).gt.50.d0)then
        write(*,*)' Overflow in besseldis_cs'
        ierr = 1
        goto 2
       else
        besdis=exp(expo)*fff
       endif
      elseif(no.lt.0)then
       expo=-zp-zs+log(zs)*fnu-alogfac 
       if(dble(expo).lt.-50.d0)then
        besdis=0.d0
       elseif(dble(expo).gt.50.d0)then
        write(*,*)' Overflow in besseldis_cs'
        ierr = 1
        goto 2
       else
        besdis=exp(expo)*fff
       endif
      endif
      ierr = 0
      return
      end
c---------------------------------------------------------------------
c      Complex subroutine calculates the Bessel distribution from nmin
c      to nmax for complex argument. 
c      Uses S18DEF from nag_shr/opt 
c      zs,zp and besdis(81) must be declared complex*16 in
c      calling program even if zs and zp are real
c------------------------------------------------------------------
      subroutine besseldist_cs(zs,zp,nmin,nmax,bdis,ierr)
      implicit real*8(a-h,o-z)
      parameter (nrange=100,ibesrange=2*nrange+1,nn0=nrange+1)
      complex*16 z,zs,zp,fact,rat,szs,szp,bdis(ibesrange),cz(nn0)
      character*1 scale
c
      scale = 's'
      if(nmin.lt.-nrange.or.nmax.gt.nrange.or.nmax.lt.nmin)then
       write(*,*)' nmin or nmax outside range in besseldist_cs'
       ierr = 1
       goto 2            
      endif
      szs = sqrt(zs)
      szp = sqrt(zp)
      z = 2.d0*szs*szp
      rat = szs*szp**(-1)
      fact = exp(-zs-zp+abs(dble(z)))
      ifail = 0
      if(nmin.gt.0)goto 20
      if(nmax.lt.0)goto 21
      fnu = 0.d+0
      nnmax = max(nmax,-nmin)
      n=nnmax+1
      call S18DEF(fnu,z,n,scale,cz,nz,ifail)
      if(ifail.ge.1)then
       write(*,*)' Error in Besseldist_cs'
       ierr = 1
       goto 2
      endif
      do 10 ii = 0,nmax
10    bdis(nn0+ii)=fact*rat**ii*cz(ii+1)
      do 11 ii =nmin,-1
11    bdis(nn0+ii)=fact/rat**(-ii)*cz(-ii+1)           
      ierr = 0
2     return
20    fnu = dble(nmin)
      n = nmax - nmin + 1
      call S18DEF(fnu,z,n,scale,cz,nz,ifail)
      if(ifail.ge.1)then
       write(*,*)' Error in Besseldist_cs'
       ierr = 1
       goto 2
      endif
      do 12 ii = nmin,nmax
12    bdis(nn0+ii)=fact*rat**ii*cz(ii-nmin+1)
      ierr = 0
      return
21    fnu = dble(-nmax)
      n = nmax-nmin + 1
      call S18DEF(fnu,z,n,scale,cz,nz,ifail)
      if(ifail.ge.1)then
       write(*,*)' Error in Besseldist_cs'
       ierr = 1
       goto 2
      endif
      do 13 ii = nmin,nmax
13    bdis(nn0+ii)=fact/rat**(-ii)*cz(-ii+nmax+1)
      ierr = 0
      return
      end        
c___________________________________________________________
c      Subroutine with reverse communication
c      calculates zero of a real function in an interval.
c      The function (func) and the derivative (dfunc)
c      in the point x are used as input.
c      At entry x should be initial guess for the zero,
c      while func and dfunc are the values in this point.
c      xneg and xpos are values of x where func is known 
c      to be negative and positive respectively
c      At exit ierr = 2 if new try is to be made. Then:
c      x contains the new guess for the zero, and func,dfunc
c      should be calculated and inserted. xneg and xpos are changed
c      and contain upper and lower limits for x 
c      When, at exit ierr = 0, 
c      x contains resulting zero to specified accuracy
c      nmax is the maximum number of calls, while n contains the 
c      actual number of calls. Must not be touched but
c      n must be 0 at first entry.
c----------------------------------------------------------------------

      subroutine zero_s(func,dfunc,xneg,xpos,x,n,nmax,accur,ierr)
      implicit real*8(a-h,o-z)
      xin = x
      if(abs(xin-xneg).lt.accur.or.abs(xin-xpos).lt.accur)then
       ierr = 0
       return
      endif
      xnext = x - func/dfunc
      if(xpos.lt.xneg)then
       if(x.gt.xneg.or.x.lt.xpos)then
        ierr = 1
        print*,' Error in zero_s: x outside interval'
        return
       endif 
       if(func.gt.0.d0)xpos=x
       if(func.lt.0.d0)xneg=x       
       if(xnext.gt.xpos.and.xnext.lt.xneg)then
        x=xnext
       else
        x = (xneg + xpos)/2.d0
       endif
      elseif(xpos.gt.xneg)then
       if(x.lt.xneg.or.x.gt.xpos)then
        ierr = 1
        print*,' Error in zero_s: x outside interval'
        return 
       endif
       if(func.gt.0.d0)xpos=x
       if(func.lt.0.d0)xneg=x       
       if(xnext.lt.xpos.and.xnext.gt.xneg)then
        x=xnext
       else
        x = (xneg + xpos)/2.d0
       endif
      else
       ierr = 1
       print*,' Error in zero_s: xneg=xpos'
      endif
      ierr = 2
      n=n+1
      if(n.gt.nmax)then
       write(*,'(a,i5,a/a)')
     &' Error in zero_s: zero not found after',nmax,' iterations',
     &' Check whether func(xneg)<0 and func(xpos)>0'
       ierr = 1
       return
      endif
      if(abs(xin-xnext).lt.accur)ierr = 0
      return
      end
c------------------------------------------------------------
c      Subroutine manipulates tabulated functions of one variable
c       with nso value of increasing or decreasing x values in 
c       sox(1000) and corresponding y values in soy(1000) 
c       to produce an equidistant table with steplength delx
c       in rex(1000), rey(1000) with increasing x.
c------------------------------------------------------------------
      subroutine interpolate_s(sox,soy,nso,rex,rey,delx,ierr)      
      implicit real*8(a-h,o-z)
      dimension x(1000),y(1000),rex(1000),rey(1000),dydx(1000)
     &,sox(1000),soy(1000)
      nre = int(abs(sox(nso)-sox(1))/delx + 1.01)
      if(nre.gt.1000.or.nso.gt.1000)then
       write(*,*)' nre.or.nso.gt.1000 in interpolate_s'
       goto 2
      endif
      do n=1,1000
       rex(n)=0.d0
       rey(n)=0.d0
      enddo
      if(sox(nso).lt.sox(1))then
       do n=1,nso
        x(n)=sox(nso-n+1)
        y(n)=soy(nso-n+1)
       enddo
      else
       do n=1,nso
        x(n)=sox(n)
        y(n)=soy(n)
       enddo
      endif
      ifail = 0
      call E01BEF(nso,x,y,dydx,ifail)
      if(ifail.ne.0)goto 2
      do n=1,nre-1
       rex(n)=x(1)+dble(n-1)*delx
      enddo
      rex(nre)=x(nso)
      ifail=0
      call E01BFF(nso,x,y,dydx,nre,rex,rey,ifail)
      if(ifail.ne.0)goto 2
      ierr = 0
      return
2     ierr = 1
      return
      end
