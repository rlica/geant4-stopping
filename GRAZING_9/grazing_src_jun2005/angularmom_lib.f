c      Library contains functions and subroutines for
c      angularmomentum algebra.
c      They have no common blocks      
c=================================================================
c      Function calculates Clebsch-Gordon coefficient
c      <j1,m1,j2,m2|j,m>. Uses threej function
c      If used several times in same program ierr=1 as input
c      avoids recalculation of factorials
c----------------------------------------------------------------

      real*8 function clebsch(rj1,rm1,rj2,rm2,rj,rm,ierr)
      implicit real*8(a-h,o-z)
      i = abs(rj1-rj2+rm)+.0001d0
      clebsch = (-1)**i*sqrt(2.d0*rj+1.d0)
     &*threej(rj1,rm1,rj2,rm2,rj,-rm,ierr)
      return
      end
c--------------------------------------------------------------------
c      function computes the three-j-symbol
c      input ierr = 1 avoids the calculation of factorials array b(200)
c          from subroutine fac_a
c      sets threej = 0 if any abs(m).gt.j
c-----------------------------------------------------------------------------
      real*8 function threej(rj1,rm1,rj2,rm2,rj3,rm3,ierr)
      implicit real*8(a-h,o-z)
      dimension b(200),iex(200)
      if(rj1.lt.0..or.rj2.lt.0..or.rj3.lt.0.)goto 201
      if(ierr.eq.1)goto 20
      a = 85.d0
      call fac_a(a,b,iex)
      asqrt = sqrt(a)
20    j1=2.001*rj1
      j2=2.001*rj2
      j =2.001*rj3
      m1=2.001*rm1
      m2=2.001*rm2
      m =2.001*rm3
      if(mod(j1-m1,2).ne.0.or.mod(j2-m2,2).ne.0.or.mod(j-m,2).ne.0)
     &goto 201
      if(mod(j1+j2+j,2).ne.0.or.mod(m1+m2+m,2).ne.0)goto 201
      sum=0.0d0
      if((m1+m2+m.ne.0).or.(j.gt.j1+j2).or.(j.lt.iabs(j1-j2))) goto 400
      if((iabs(m1).gt.j1).or.(iabs(m2).gt.j2).or.(iabs(m).gt.j))goto 400
      k2min=0
      if(j-j2+m1.lt.0) k2min=-j+j2-m1
      if(j-j1-m2+k2min.lt.0) k2min=-j+j1+m2
      k2max=j1+j2-j
      if(j2+m2-k2max.lt.0) k2max=j2+m2
      if(j1-m1-k2max.lt.0) k2max=j1-m1
      ja=(j1+m1)/2+1
      jb=ja-m1
      jc=(j2+m2)/2+1
      jd=jc-m2
      je=(j +m )/2+1
      jf=je-m
      jg=(j1+j2-j)/2+1
      jh=ja+jb-jg
      ji=jc+jd-jg
      jj=je+jf+jg-1
      if (jj .gt. 200) goto 600
c      600 stops the computation if the range of factorials is exceeded
      b1=b(ja)*b(jb)*b(jc)*b(jd)*b(je)*b(jf)*b(jg)*b(jh)*b(ji)/b(jj)
      iex1=iex(ja)+iex(jb)+iex(jc)+iex(jd)+iex(je)+iex(jf)+iex(jg)+
     &iex(jh)+iex(ji)-iex(jj)
      ia=k2min/2
      ib=jg-ia+1
      ic=jb-ia+1
      id=jc-ia+1
      ie=ja-jg+ia
      if=jd-jg+ia
      fase=1.0d0
      if(mod(ia,2).eq.0) fase=-fase
      k2 =k2min
c
c     summation over k starts here
100   ia=ia+1
      ib=ib-1
      ic=ic-1
      id=id-1
      ie=ie+1
      if=if+1
      fase=-fase
      b2=b(ia)*b(ib)*b(ic)*b(id)*b(ie)*b(if)
      iex2=iex(ia)+iex(ib)+iex(ic)+iex(id)+iex(ie)+iex(if)
      stor=sqrt(b1/(b2*b2) * 10.0d0**(iex1-iex2-iex2))
      sum = sum + fase*stor
      k2=k2+2
      if(k2.le.k2max) goto 100
      n = abs(rj1-rj2-rm3)+0.01d0
      sum = (-1)**n * sum
400   threej=sum/asqrt
      ierr = 0
      return
c
c     error exit
600   print*,'error - factorials exceeded in threej'
      print*,'j1=',rj1,'j2=',rj2,'j3=',rj3
      threej = 0.d0
      ierr = 1
      return
201   print*,' wrong input in threej'
      threej = 0.d0
      ierr = 1
      return
      end

c     subroutine calculates factorial arrays  b  and  iex
c--------------------------------------------------------------------
      subroutine fac_a(a,b,iex)
      implicit real*8(a-h,o-z)
      dimension b(200),iex(200)
      asqrt = sqrt(a)
      f = 1.0d0
      b(1) = 1.d0
      iex(1) = 0
      do 80 i=1,199
      f = f * float(i)/a
      iex(i+1) = int(log10(f))
      b(i+1) = f * 10.0d0**(-iex(i+1))
80    continue
      return
      end
c--------------------------------------------------------------------------
C   This function calculates the  d-function using the following
C  relations (J1,M1 and M2 should have twice their actual values
C  when the function is called):
C
C  *********************************************************************
C  *                                                                   *
C  *  J1            M1-M2       (J1+M1)!*(J1-M1)!      THET M1+M2      *
C  * D     (THET)=(-)     *SQRT(-----------------)*COS(----)     *     *
C  *  M1 M2                     (J1+M2)!*(J1-M2)!        2             *
C  *                                                                   *
C  *                              -----                                *
C  *                                                                   *
C  *        THET M1-M2  (M1-M2,M1+M2)                                  *
C  *    SIN(----)     *P             (COS(THET))                       *
C  *          2         J1-M1                                          *
C  *                                                                   *
C  *********************************************************************
C
C       This is for M1.GE.M2. for M1.LT.M2 one has:
C
C  *********************************************************************
C  *                                                                   *
C  *             J1             M1-M2  J1                              *
C  *            D     (THET)=(-)     *D    (THET)                      *
C  *             M1 M2                 M2 M1                           *
C  *                                                                   *
C  *********************************************************************
C
C     The Jacobi polynomials P are calculated by recursion according
C       to Abramowitz & Stegun 22.7.1
C
C     Code by B.S. Nilsson, NBI, Copenhagen
C     This version is from 1976-08-23 but modified 1985-01-25.
c     Modified by A.Winther 14.10.1992 (accuracy e-6 up to rj=40)

      real*8 function ddfunc(rj,rm,rmm,theta,ierr)
      implicit real*8(a-h,o-z)
      thet=theta
      x=cos(thet)
      j1=2.0001d0*rj
      m12=2.0001d0*rm
      m22=2.0001d0*rmm
      if(mod(j1-m12,2).ne.0.or.mod(j1-m22,2).ne.0)goto 31
      m1=m12
      m2=m22
      if(iabs(m1).gt.j1.or.iabs(m2).gt.j1) then
        ddfunc=0.d0
        return
       endif
      if(m1+m2.ge.0) goto 110
      m1=-m22
      m2=-m12
110   iph=0
      if(m1.ge.m2) goto 120
      iph=m1
      m1=m2
      m2=iph
      iph=1
120   ia=(m1-m2)/2
      ib=(m1+m2)/2
      l=(j1-m1)/2
      p1=1.d0
      if(l.ne.0) goto 130
      p3=1.d0
      go to 300
130   p2=x*(float(m1)*0.5d0+1.d0)-float(m2)*0.5d0
      if(l.ne.1) goto 140
      p3=p2
      go to 300
140   m12n=m1+2
      do 200 n=2,l
      m12n=m12n+2
      a1=n*(n+m1)*(m12n-2)*2
      a2=-m1*m2*(m12n-1)
      a3=(m12n-2)*(m12n-1)*m12n
      a4=(n+ia-1)*(n-1+ib)*m12n*2
      p3=((a2+a3*x)*p2-a4*p1)/a1
      p1=p2
200   p2=p3
300   i1=(j1+m1)/2
      i2=(j1-m1)/2
      i3=(j1+m2)/2
      i4=(j1-m2)/2
      f1=1.d0
      f2=1.d0
      if(ib.ne.0) f1=cos(0.5d0*thet)**ib
      if(ia.ne.0) f2=sin(0.5d0*thet)**ia
      dux1 = duxfct(i3,i1)
      dux2 = duxfct(i4,i2)
      sq = sqrt(dux1*dux2)
      ddfunc=p3*sq*f1*f2
c      print*,'sq,dux1,dux2 =',sq,dux1,dux2
      if(mod(m1-m2,4).ne.0.and.iph.eq.0) ddfunc=-ddfunc
      return
31    print*,' Error in input in ddfunc'
      ierr = 1
      ddfunc=0.d0
      return
      end

      real*8 function duxfct(m,n)
      implicit real*8 (a-h,o-z)
c      calculates n!/m!*60**(n-m)
      duxfct=1.d0
      numax=m-n
      if(numax.lt.0) then
        ictrl=1
        numax=-numax
        fctin=m
       else if(numax.gt.0) then
        ictrl=0
        fctin=n
       else
        return
       endif
      do 100 nu=1,numax
      fctor=(fctin+nu)/60.d0
100   duxfct=duxfct*fctor
      if(ictrl.eq.0.and.abs(duxfct).gt.1.e-20) duxfct=1.d0/duxfct
      return
      end
c---------------------------------------------------------------------------
c      subroutine for spherical harmonic Y_lm(pi/2,0)
c
c--------------------------------------------------------------

      subroutine ylm_s(ylm,l,m,ifail)
      implicit real*8(a-h,o-z)
      if(abs(m).gt.l)then
       ifail = 1
       print*,' abs(m).gt.l in ylm_s'
       return
      endif
      if(mod(l+m,2).ne.0)then
       ylm = 0.d0
       return
      endif
      iarg1 = l+m
      iarg2 = l-m
      ylm = sqrt((2.*l+1.)/4./3.141592653589793238d+0)
     &*facratt(iarg1)*facratt(iarg2)*(-1)**(iarg1/2)
      return
      end

      real*8 function facratt(n)
      implicit real*8(a-h,o-z)
      facratt = 1.d0
      nn = 0
1     if(nn.ge.n)goto 2
      nn = nn + 2
      facratt = facratt*float(nn-1)/float(nn)
      goto 1
2     facratt = sqrt(facratt)
      return
      end
c----------------------------------------------------------------
c      Function calculates six-j symbol for l_1 = 1/2
c----------------------------------------------------------------

      real*8 function sixj1(rj1,rj2,rj3,rl2,rl3,ierr)
      implicit real*8(a-h,o-z)
      jd1 = 2.0001d0*rj1
      jd2 = 2.0001d0*rj2
      jd3 = 2.0001d0*rj3
      ld2 = 2.0001d0*rl2
      ld3 = 2.0001d0*rl3
      is = (jd1+jd2+jd3)/2
      if(ld2.eq.jd3-1)goto 20
      if(ld2.eq.jd3+1)goto 21
      goto 3
20    if(ld3.eq.jd2+1)goto 22
      if(ld3.eq.jd2-1)goto 23
      goto 3
22    sixj1 = (-1)**is*sqrt((rj1+rj3-rj2)*(rj1+rj2-rj3+1.d0)
     &/((2.d0*rj2+1.d0)*(2.d0*rj2+2.d0)*2.d0*rj3*(2.d0*rj3+1.d0)))
      goto 4
23    sixj1 = (-1)**is*sqrt((rj1+rj3+rj2+1.d0)*(rj2+rj3-rj1)
     &/(2.d0*rj2*(2.d0*rj2+1.d0)*2.d0*rj3*(2.d0*rj3+1.d0)))
      goto 4
21    if(ld3.eq.jd2-1)goto 24
      if(ld3.eq.jd2+1)goto 25
      goto 3
24    sixj1 = (-1)**is*sqrt((rj1+rj3-rj2+1.d0)*(rj1+rj2-rj3)
     &/(2.d0*rj2*(2.d0*rj2+1.d0)*(2.d0*rj3+1.d0)*(2.d0*rj3+2.d0)))
      goto 4
25    sixj1 = (-1)**(is+1)*sqrt((rj1+rj3+rj2+2.d0)*(rj2+rj3-rj1+1.d0)
     &/((2.d0*rj2+1.d0)*(2.d0*rj2+2.d0)*(2.d0*rj3+1.d0)
     &*(2.d0*rj3+2.d0)))
      goto 4
3     print*,' Error in input in sixj1'
      sixj1 = 0.d0
      ierr = 1
      return
4     ierr = 0
      return
      end

