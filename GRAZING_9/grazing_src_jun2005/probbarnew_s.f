c      Subroutine calculates energy widths at capture distance  and
c       probability distribution of barriers
c---------------------------------------------------------------
      subroutine probbar_s(inelast,iout,rcap,alcap,smin,rdur,phid,eta2
     &,ierr)
      implicit real*8(a-h,o-z)
      common/za/za,aa,zaa,aaa
      common/ecm/ecm
      common/probc/adia,capfac,ecmcapt,widcbar,sigcest,icapt
      common/capt/eb,rbar,i5,intp
      common/capr/capred,pow
      common/trans/rturn,rgraz,algraz,homegg,alatr
      common/vdr/vdr,accur
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      common/mml/ssl,qvalcom
      common/iza/iza,iaa,izaa,iaaa
      common/q/qns,qnp,qps,qpp
      common/cap/capn,capp,rr
      common/ba/bna,bnaa,bpa,bpaa
      common/ks/akbarn,akbarp
      common/ebe23/eex(2,3,2)
      common/zp23/zp(2,3,2)
      common/widcin/wi2cin(2,3,2)
      common/widcns/wid2cn,wid2cp
      common/widcc/widc,widcin
      common/bhom/bhomn,bhomp
      common/probbar/delebar(481),barprob(481),ebarstep
      dimension forms(2,2,3,4),formst(2,2,3,4),sfuncst(2,2,3,4)
     &,eavst(2,2,3,4),wid2st(2,2,3,4),eta2(2,3,2,4),eta2a(2,3,2,4)
     &,rdur(2,3)
      dimension gn(11),eacin(2,3,2)
      dimension eben(10),alfn(10),ecomp(10),nselect(10),iselect(10),
     &lamselect(10),mselect(10),wid2sp(10)
c      
      time1 = X05BAF()
c  For the calculation of transfer widths analytically the surface 
c   distance is increasedby a small amount to make a reasonable 
c   collision time. Should be estimated in capture_s (like smin)
      vdrstore = vdr
      vdr = vdr-.5d0
      call constant_s
      ecmstore = ecm
      ecm1 = ecmcapt
      rcapt = rcap
      alcapt = alcap
      am0 = aa*aaa/(aa+aaa)/hc2om      
      nmax = 2
      nmin = 1
      if(inelast.eq.0)nmin = 2
      if(inelast.eq.1.or.inelast.eq.2)nmax = 1
      if(inelast.eq.-1)nmax = 0
c      
      zeta = exp(-.4d0/50.d0*ecm1*(1.d0/aa+1.d0/aaa))
      zeta=1.d0
c   calculate the hbar*omega for transfer          
      cappan = capn/2.d0
      cappap = capp/2.d0
      ierr1 = 1
      dur = dueff(alcapt,rcapt,ierr1)
      hom2 = - cappan/am0*dur
      if(ierr1.eq.1)goto 2
      if(hom2.le.0.)then
       write(*,*)' Error in probbar hom2 = ',hom2,' set = .01'
        hom2 = .01d0
c       ierr = 1
c       goto 2
      endif
      homegan = sqrt(hom2)*adia
      homegap = homegan*sqrt(cappap/cappan)
      homegain = sqrt(hom2/cappan*ainv)
c  calculate phidot and alpha3
      if(icapt.eq.2)then
       phidot = 0.d0
      else
       ierr1 = 1
       phidot2 = 2.d0/am0*(ecm1-unuc(rcapt,ierr1)
     &-za*zaa*esq/rcapt)/rcapt**2
       if(ierr1.eq.1)goto 2
       if(phidot2.lt.0.d0)then
        write(*,*)' Warning: Error in probbar_s phidot^2 = ',phidot2
        phidot2 = 0.d0
       endif
       phidot = sqrt(phidot2)
      endif
c
      alpha3n = phidot/homegan
      alpha3p = phidot/homegap
c  calculate ktrans
      transk = rcapt*phidot/hc2om
      cn = transk/cappan
      cp = transk/cappap
c calculate ba, baa, bb, dd and alphs3
      alf32n = alpha3n**2
      alf32p = alpha3p**2
      bn = bna+bnaa
      bp = bpa+bpaa
      bbn = sqrt(1.d0+ 2.d0*bn*alf32n) 
      bbp = sqrt(1.d0+ 2.d0*bp*alf32p) 
      ddn = 1.d0 + .5d0*bn*ssl
      ddp = 1.d0 + .5d0*bp*ssl
      c3n2 = (1.d0 + 2.d0*bna*bnaa/bn*ssl)/ddn
      c3p2 = (1.d0 + 2.d0*bpa*bpaa/bp*ssl)/ddp
      alphs3n = alpha3n/bbn/ddn
      alphs3p = alpha3p/bbp/ddp
      bhomn = bbn*homegan
      bhomp = bbp*homegap
c  calculate k_n and k_p 
      eeen = (capn+cn**2/4.d0*cappan)*smin
      eeep = (capp+cp**2/4.d0*cappap)*smin
      if(eeen.gt.70..or.eeep.gt.70.)then
       en = 0.d0
       ep = 0.d0
      else
       en = 1.d0/(1.d0+exp(eeen))
       ep = 1.d0/(1.d0+exp(eeep))
      endif
c
      akn = adia*akbarn*bbn/sqrt(ddn)*zeta**2*en
      akp = adia*akbarp*bbp/sqrt(ddp)*zeta**2*ep
c calculate Q_opt
      qoptns = 0.
      qoptps = qoptns + (za-zaa)*esq/rcapt
c calculate Qs
      qsns = -(qns-qoptns)/bhomn
      qsnp = -(qnp+qoptns)/bhomn
      qsps = -(qps-qoptps)/bhomp
      qspp = -(qpp+qoptps)/bhomp
c calculate error functions
      acc = accur/1000.d0
      ejns = qsns + (bna-bnaa)*cn*alphs3n
      ejcns = ejns/capfac
      if(abs(ejcns).gt.8.)then
       efaccns = 0.d0
      else
       efacns = exp(-ejns**2)
       efaccns = exp(-ejcns**2)
      endif
      if(ejcns.lt.-5.)then
       call ierfc(1,ejcns,gcns)
      else
       ierr = 0
       call gn_s(ejcns,1,gn,acc,ierr)
       if(ierr.eq.1)goto 2
       if(ierr.ne.0)write(*,*)
     &' Warning: no convergence for capture in gn_s for nsc'
       gg1cns = gn(2)     
       gcns = gg1cns*efaccns
      endif
      ejnp = qsnp + (bnaa-bna)*cn*alphs3n
      ejcnp = ejnp/capfac
      if(abs(ejcnp).gt.8.)then
       efaccnp = 0.d0
      else
       efaccnp = exp(-ejcnp**2)
      endif
      if(ejcnp.lt.-5.)then
       call ierfc(1,ejcnp,gcnp)
      else
       ierr = 0
       call gn_s(ejcnp,1,gn,acc,ierr)
       if(ierr.eq.1)goto 2
       if(ierr.ne.0)write(*,*)
     &' Warning: no convergence for capture in gn_s for npc'
       gg1cnp = gn(2)     
       gcnp = gg1cnp*efaccnp
      endif
      ejps = qsps + (bpa-bpaa)*cp*alphs3p
      ejcps = ejps/capfac
      if(abs(ejcps).gt.8.)then
       efaccps = 0.d0
      else
       efaccps = exp(-ejcps**2)
      endif
      if(ejcps.lt.-5.)then
       call ierfc(1,ejcps,gcps)
      else
       ierr = 0
       call gn_s(ejcps,1,gn,acc,ierr)
       if(ierr.eq.1)goto 2
       if(ierr.ne.0)write(*,*)
     &' Warning: no convergence for capture in gn_s for psc'
       gg1cps = gn(2)     
       gcps = gg1cps*efaccps
      endif
      ejpp = qspp + (bpaa-bpa)*cp*alphs3p
      ejcpp = ejpp/capfac
      if(abs(ejcpp).gt.8.)then
       efaccpp = 0.d0
      else
       efaccpp = exp(-ejcpp**2)
      endif
      if(ejcpp.lt.-5.)then
       call ierfc(1,ejcpp,gcpp)
      else
       ierr = 0
       call gn_s(ejcpp,1,gn,acc,ierr)
       if(ierr.eq.1)goto 2
       if(ierr.ne.0)write(*,*)
     &' Warning: no convergence for capture in gn_s for ppc'
       gg1cpp = gn(2)     
       gcpp = gg1cpp*efaccpp
      endif
      c0n = exp((bn/2.d0/ddn-cappan/4.d0*(rr+vdr))*cn**2)
      c0p = exp((bp/2.d0/ddp-cappap/4.d0*(rr+vdr))*cp**2)
c   rr is sum of nuclear radii for the densities (not equal to r0)
c calculate F_0's
      f0cns = akn*c0n*gcns
      f0cnp = akn*c0n*gcnp 
      f0cps = akp*c0p*gcps
      f0cpp = akp*c0p*gcpp
c  Spread at distance of closest approach
      wid2cn = homegan**2/(2.d0*pi)*(f0cns+f0cnp)
      wid2cp = homegap**2/(2.d0*pi)*(f0cps+f0cpp)
      widcn = sqrt(wid2cn)
      widcp = sqrt(wid2cp)
      eavcina = 0.d0
      eavcinaa = 0.d0
      wid2cina = 0.d0
      wid2cinaa = 0.d0
      vdr = vdrstore
c  restoring original vdr
      call constant_s
c
c  alternatively from subroutine capture_s:
      phidot = phid
      do n = nmin,nmax
       do lam = 2,3
        do m=-lam,lam,2
         mm = (lam+m)/2+1
         ierr1 = 0
         yl = (-1)**m*ylm(lam,m,ierr1)
         if(ierr1.eq.1)goto 2
         forms(1,n,lam,mm) = -(-1)**lam*zp(1,lam,n)*yl*rdur(1,lam)
         forms(2,n,lam,mm) = -zp(2,lam,n)*yl*rdur(2,lam)
        enddo
       enddo
      enddo
      do n=nmin,nmax
       do i=1,2
        do lam = 2,3
          eacin(i,lam,n)=0.d0
          wi2cin(i,lam,n)=0.d0
        enddo
       enddo
      enddo
      eavcina = 0.d0
      eavcinaa = 0.d0
      do n=nmin,nmax
       do i=1,2
        do lam = 2,3
         do m = -lam,lam,2
          mm = (lam+m)/2+1
          eeff = eex(i,lam,n)-dble(m)*phidot
          if(eeff.eq.0.d0)eeff = 1.d-20
          eavst(i,n,lam,mm)= eeff*
     &(eta2(i,lam,n,mm)-(forms(i,n,lam,mm)/eeff)**2)
          eacin(i,lam,n) = eacin(i,lam,n)+eavst(i,n,lam,mm)
          wid2st(i,n,lam,mm)= eeff**2*eta2(i,lam,n,mm)
          wi2cin(i,lam,n) = wi2cin(i,lam,n)+wid2st(i,n,lam,mm)
         enddo
c         write(*,'(a,i1,a,i1,a,i1,a,4f12.5)')
c     &' eta2(',i,',',lam,',',n,',m) = ',(eta2(i,lam,n,mm),mm=1,lam+1)
        enddo
       enddo
      enddo
      do n = nmin,nmax
       do lam = 2,3
        eavcina = eavcina + eacin(1,lam,n)
        eavcinaa = eavcinaa + eacin(2,lam,n)
        wid2cina = wid2cina + wi2cin(1,lam,n)
        wid2cinaa = wid2cinaa + wi2cin(2,lam,n)
       enddo
      enddo
c
      eavcin = eavcina+eavcinaa
      wid2cin = wid2cina+wid2cinaa
      widcin = sqrt(wid2cin)
      wid2ctr = wid2cn+wid2cp
      widctr = sqrt(wid2ctr)
      wid2c = wid2ctr+wid2cin
      widc = sqrt(wid2c)
c
      iwarn = 0
c
c  Now start selecting states for discrete barrier distribution
      iselectmax = 10
c  this number can be le. 10 (number of do-loops below)
      nma = 24
      esmin1 = .4d0
      esmin = esmin1
      probmin1 = .001d0
      probmin = probmin1
      nselect(10) = 0
c  Select collective (m-)states to be treated explicitly (at most 10!)
9     phidotmin = esmin/3.d0
      do ise = 1,10
       ecomp(ise) = 0.d0
       eben(ise) = 0.d0
       alfn(ise) = 1.d-8      
      enddo
      ise = 10
      wid2cred = wid2c
      do n=nmin,nmax
       do i=1,2
        do lam=2,3
c there are two main cases:
         if(phidot.gt.phidotmin)then
c  m-states included explicitly:
          do m=-lam,lam,2
           mm = (lam+m)/2+1
c    energystep in Poisson:
           alf = eex(i,lam,n)-dble(m)*phidot

c    formfactor:
           fos = forms(i,n,lam,mm)
c    width:
           wi2 = wid2st(i,n,lam,mm)
c    adiabatic energy shift:
           comp = fos**2/alf
c    "probability" in Poisson (mean number of phonons): 

           probb = eta2(i,lam,n,mm) 
c    shift of ground state to make mean energy = 0 (-mean energy in Poisson):
           shiftgs = -probb*alf
           if(alf.gt.esmin.and.probb.lt.2.5d0.
     &and.probb.gt.probmin)then
c    i.e. if energy splitting is large and probability not too large
            if(ise.le.10-iselectmax)then
             iwarn = 1
             if(probmin.ge.probmin1*32.1d0)then
              esmin = esmin + esmin1
              probmin = probmin1
             else
              probmin = probmin*2.d0
             endif
             goto 9
            endif
            nselect(ise) = n
            iselect(ise) = i
            lamselect(ise) = lam
c    mselect ne. 99 indicates that m-states are treated explicitly:
            mselect(ise) = m
c    reduction in total width (inelastic+transfer) due to states treated 
c    explicitly, ie. the width**2 of each of these states:
            wid2cred = wid2cred - wi2
c    spreading due to m-states not applicable in this case:
            wid2sp(ise) = 0.d0
            eben(ise) = alf
            alfn(ise) = probb
            ecomp(ise)= shiftgs
            ise = ise - 1
           endif
          enddo
c  second possibility:
         elseif(phidot.le.phidotmin)then
c    i.e.no need to treat each m-state separately
          probspr = 0.d0
c    Width**2 of states included explicitly:
          wid2 = 0.d0
c    shift of ground state makeing mean energy = 0 (-mean energy in Poisson):
          shiftgs = 0.d0
c    Total probability for lowest state:
          probb = 0.d0
c    Mean excitation energy of lowest state:
          alfm = 0.d0
c    here wid2spr (double-sum) and wid2 (single sum) are calculated 
          do m = -lam,lam,2
           mm = (lam+m)/2+1
           alf1 = eex(i,lam,n)-dble(m)*phidot
           probm1 = eta2(i,lam,n,mm)
           fosfu = probm1*alf1**2
           wid2 = wid2+ fosfu
           probb = probb + probm1
           alfm = alfm+alf1*probm1
           shiftgs = shiftgs - fosfu/alf1
           do m2 = -lam,lam,2
            mm2 = (lam+m2)/2+1
            alf2 = eex(i,lam,n)-dble(m2)*phidot
            probm2 = eta2(i,lam,n,mm2)
            probspr = probspr + probm1*probm2*dble((mm-mm2)**2)            
           enddo
          enddo
          if(probb.ne.0.d0)then
c    Width**2 of included states due to non-vanishing phidot:
           wid2spr = probspr*phidot**2/2.d0/probb**2
           alf = alfm/probb
          endif
          if(alf.gt.esmin.and.probb.lt.4.d0
     &.and.probb.gt.probmin)then
           if(ise.le.10-iselectmax)then
            iwarn = 1
            if(probmin.gt..05d0)then
             esmin = esmin + esmin1
             probmin = probmin1
            else
             probmin = probmin*2.d0
            endif
            goto 9
           endif
           nselect(ise) = n
           iselect(ise) = i
           lamselect(ise) = lam
c  this indicates that m-states are not included explicitly:
           mselect(ise) = 99
c  this is the total width**2 apart from total of that below 
           wid2cred = wid2cred - wid2
c  this is the final width ascribed to each state due to phidot-width
           wid2sp(ise) = wid2spr
           eben(ise) = alf
           alfn(ise) = probb
           ecomp(ise)= shiftgs
           ise = ise - 1
          endif
         endif
        enddo
       enddo
      enddo
      if(iwarn.eq.1)then
       write(*,'(/a)')
     &' Warning: more states need special treatment in probbar_s' 
       write(*,'(a,2f6.3)')' esmin,probmin increased to ',esmin,probmin
      endif
      isemin = ise+1
      if(iout.eq.1.and.nselect(10).ne.0)then
       if(mselect(10).eq.99)then
        write(*,*)' selected states for special treatment in probbar_s:'
        do is = isemin,10
         write(*,'(a,3i2)')' i, lam, n = ', iselect(is),
     &lamselect(is),nselect(is)
         write(*,'(a,3f9.5)')' probability,shift_gs,estep ='
     &,alfn(is),ecomp(is),eben(is)
        enddo
       elseif(mselect(10).ne.99)then
      write(*,*)' selected m-states for special treatment in probbar_s:'
        do is = isemin,10
         write(*,'(a,4i2)')' i, lam, n, m = ', iselect(is),
     &lamselect(is),nselect(is),mselect(is)
         write(*,'(a,3f9.5)')' probability,shift_gs,estep ='
     &,alfn(is),ecomp(is),eben(is)
        enddo
       endif
      endif
c  this sets the range and steps in the probability distribution
      
      range = widc*sqrt(2.d0*log(100.d0/accur))
      ebarstep = range/120.d0
      do i = 1,481
       barprob(i) = 0.d0
       delebar(i) = -range+dble(i-1)*ebarstep
      enddo
c  i.e. mean barrier at i=121  (delebar(121) = 0)
      write(*,'(a)')' Calculating barrier distribution (probbar_s)'
      widcred = sqrt(wid2cred)
c
c  Here starts calculation of barrier distribution
      decmr0 = 0.d0
      probtot = 0.d0
c  For each selected mode (up to 10) one calculates
      iii = 0
      acc100 = accur/100.d0
      acc1000 = accur/1000.d0
c  for factorial:
      fac1 = 1.d0
c  mean energy in 1. mode:
      decmr10=decmr0-ecomp(1)
c  for each number of phonons in this mode:
      do n1 = 0,nma
c  the factorial:
       fac1 = fac1*dble(max(n1,1))
c  the negative of the energy ie. the change in kinetic energy      
       decmr1 = decmr10-eben(1)*dble(n1)
c  multiplied by its probability
       probn1 = alfn(1)**n1/fac1*exp(-alfn(1))
       probnn1 = probn1
       if(n1.gt.alfn(1).and.probnn1.le.acc100)goto 11
c  ie. if finished here all modes are finished
c  else go to second mode:
       fac2 = 1
       decmr20=decmr1-ecomp(2)
       do n2=0,nma
        fac2 = fac2*dble(max(n2,1))
c  here total change in kinetic energy from n1 phonons in 1. and n2 in 2.
        decmr2 = decmr20-eben(2)*dble(n2)
        probn2 = alfn(2)**n2/fac2*exp(-alfn(2))
c  here the probability for a number of phonons n1 in 1. and n2 in 2.
        probnn2 = probnn1*probn2
        if(n2.gt.alfn(2).and.probnn2.le.acc100)goto 12
        fac3 = 1.d0
        decmr30=decmr2-ecomp(3)
        do n3=0,nma
         fac3 = fac3*dble(max(n3,1))
         decmr3 = decmr30-eben(3)*dble(n3)
         probn3 = alfn(3)**n3/fac3*exp(-alfn(3))
         probnn3 = probnn2*probn3
         if(n3.gt.alfn(3).and.probnn3.le.acc100)goto 13
         fac4 = 1.d0
         decmr40=decmr3-ecomp(4)
         do n4=0,nma
          fac4 = fac4*dble(max(n4,1))
          decmr4 = decmr40-eben(4)*dble(n4)
          probn4 = alfn(4)**n4/fac4*exp(-alfn(4))
          probnn4 = probnn3*probn4
          if(n4.gt.alfn(4).and.probnn4.le.acc100)goto 43
          fac5 = 1.d0
          decmr50=decmr4-ecomp(5)
          do n5=0,nma
           fac5 = fac5*dble(max(n5,1))
           decmr5 = decmr50-eben(5)*dble(n5)
           probn5 = alfn(5)**n5/fac5*exp(-alfn(5))
           probnn5 = probnn4*probn5
           if(n5.gt.alfn(5).and.probnn5.le.acc100)goto 53
           fac6 = 1.d0
           decmr60=decmr5-ecomp(6)
           do n6=0,nma
            fac6 = fac6*dble(max(n6,1))
            decmr6 = decmr60-eben(6)*dble(n6)
            probn6 = alfn(6)**n6/fac6*exp(-alfn(6))
            probnn6 = probnn5*probn6
            if(n6.gt.alfn(6).and.probnn6.le.acc100)goto 63
            fac7 = 1.d0
            decmr70=decmr6-ecomp(7)
            do n7=0,nma
             fac7 = fac7*dble(max(n7,1))
             decmr7 = decmr70-eben(7)*dble(n7)
             probn7 = alfn(7)**n7/fac7*exp(-alfn(7))
             probnn7 = probnn6*probn7
             if(n7.gt.alfn(7).and.probnn7.le.acc100)goto 73
             fac8 = 1.d0
             decmr80=decmr7-ecomp(8)
             do n8=0,nma
              fac8 = fac8*dble(max(n8,1))
              decmr8 = decmr80-eben(8)*dble(n8)
              probn8 = alfn(8)**n8/fac8*exp(-alfn(8))
              probnn8 = probnn7*probn8
              if(n8.gt.alfn(8).and.probnn8.le.acc100)goto 83
              fac9 = 1.d0
              decmr90=decmr8-ecomp(9)
              do n9=0,nma
               fac9 = fac9*dble(max(n9,1))
               decmr9 = decmr90-eben(9)*dble(n9)
               probn9 = alfn(9)**n9/fac9*exp(-alfn(9))
               probnn9 = probnn8*probn9
               if(n9.gt.alfn(9).and.probnn9.le.acc100)goto 84
               fac10 = 1.d0
               decmr100=decmr9-ecomp(10)
c  one allways gets to this point at least for n1=n2=...=0
               do n10=0,nma
                fac10 = fac10*dble(max(n10,1))
                decmr = decmr100-eben(10)*dble(n10)
                probn10 = alfn(10)**n10/fac10*exp(-alfn(10))
                probn = probnn9*probn10
c  the contribution from n1 in 1., n2 in 2. ect are added for check
                probtot = probtot + probn
                if(n10.gt.alfn(10).and.probn.le.acc1000)goto 93
c  Calculate contribution to barrier probability if probn.gt.acc100
                if(probn.lt.accur)then
                 prcom = 2.d0*accur
                else
                 prcom = probn
                endif
c   to check total number of steps:
                iii = iii + 1
c  if m-states are not treated explicitly:
                if(mselect(10).eq.99)then
c    then this is the final width ascribed to n1 in 1., n2 in 2. etc.
                 widcredd = sqrt(n1*wid2sp(1)+n2*wid2sp(2)+n3*wid2sp(3)+
     &n4*wid2sp(4)+n5*wid2sp(5)+n6*wid2sp(6)+n7*wid2sp(7)+n8*wid2sp(8)+
     &n9*wid2sp(9)+n10*wid2sp(10)+wid2cred)
                else
c  if m-states are treated explicitly:
                 widcredd = widcred
                endif
                declim = widcredd*sqrt(2.d0*log(prcom/acc100))
                do ii = 1,481
c  this is then the construction of the barier distribution
c   decmr (change in kinetic energy) appears with minus 
c   ie. delebar is actually height of barrier
                 if(delebar(ii).gt.-declim-decmr)then
                  barprob(ii) = barprob(ii) + sqrt(.5d0/pi)/widcredd*
     &exp(-((-decmr-delebar(ii))/widcredd)**2/2.d0)*probn
                  if(delebar(ii).gt.declim-decmr)goto 18
                 endif
                enddo
                write(*,*)' Error: More than 481 steps in probbar_s'
18              continue          
               enddo
               write(*,*)' Error in probbar_s: n10 too large'
               ierr = 1
               goto 2
93             continue
              enddo
              write(*,*)' Error in probbar_s: n9 too large'
              ierr = 1
              goto 2
84            continue
             enddo
             write(*,*)' Error in probbar_s: n8 too large'
             ierr = 1
             goto 2
83           continue
            enddo
            write(*,*)' Error in probbar_s: n7 too large'
            ierr = 1
            goto 2
73          continue
           enddo
           write(*,*)' Error in probbar_s: n6 too large'
           ierr = 1
           goto 2
63         continue
          enddo
          write(*,*)' Error in probbar_s: n5 too large'
          ierr = 1
          goto 2
53        continue
         enddo
         write(*,*)' Error in probbar_s: n4 too large'
         ierr = 1
         goto 2
43       continue  
        enddo
        write(*,*)' Error in probbar_s: n3 too large'
        ierr = 1
        goto 2
13      continue  
       enddo
       write(*,*)' Error in probbar_s: n2 too large'
       ierr = 1
       goto 2
12     continue  
      enddo
      write(*,*)' Error in probbar_s: n1 too large'
      ierr = 1
      goto 2
11    if(abs(probtot-1.d0).gt.20.d0*accur)then
       write(*,'(a,f10.7)')
     &' Warning in probbar_s: total Poisson probability = ',probtot 
       iwarn = 1
      endif
      probetot = 0.d0
      avbar = 0.d0
      do i = 1,481
       probetot = probetot + barprob(i)*ebarstep
       avbar = avbar + delebar(i)*barprob(i)*ebarstep
      enddo
      if(abs(probetot-1.d0).gt.20.d0*accur)then
       write(*,'(a,f10.7)')
     &' Warning in probbar_s: Total barrier probability= ',probetot
       iwarn = 1
      endif
      do i = 1,481
       barprob(i) = barprob(i)*1.d0/probetot
      enddo
      if(abs(avbar).gt.20.d0*accur)then
       write(*,'(a,f10.7)')
     &' Warning in probbar_s: Average barrier position not 0 = ',avbar
       iwarn = 1
      endif
      do i = 1,481
       delebar(i) = delebar(i)-avbar
      enddo
c
      if(iout.ne.0)then
       wt = 2.35482d0*widc
       wi = 2.35482d0*widcin
       write(*,'(a/a,f5.2,a,f5.2,a/a,f6.2,a)')
     &' Widths(FWHM) at capture distance:',
     &'   From inelastic = ',wi,' MeV, Total = ',wt,'MeV',
     &' Mean radial energy loss at capture: ',eavcin,' MeV'
       write(98,'(a/a,f5.2,a,f5.2,a/a,f6.2,a)')
     &' Widths(FWHM) at capture distance:',
     &'   From inelastic = ',wi,' MeV, Total = ',wt,'MeV',
     &' Mean radial energy loss at capture: ',eavcin,' MeV'
      endif
c  calculate estimated capture cross section
      write(*,*)' Total number of steps in probbar_s =',iii
      widcbar = widc
      ierr = 1
      ucap1 = za*zaa*esq/rcapt+unuc(rcapt,ierr)
      if(ierr.eq.1)goto 2
      sq = sqrt(2.d0)
      aarg1 = (ucap1-ecmstore)/sq/widc
      call ierfc(1,aarg1,eerf)
      sigcest = pi*rcapt**2*widcbar/sq/ecmstore*eerf*10.d0
      write(*,'(a,d15.8,a)')
     &' 1. estimate for capture cross section = ',sigcest,' mb'
      ucap2 = ecmcapt- alcapt**2/2.d0/am0/rcapt**2
      aarg2 = (ucap2-ecmstore)/sq/widc
      call ierfc(1,aarg2,eerf)
      sigcest2 = pi*rcapt**2*widcbar/sq/ecmstore*eerf*10.d0
      ierr = 0
      if(iwarn.eq.1)ierr = 2
      time2 = X05BAF()
      time = time2-time1
      write(*,'(a,f7.2,a)')' time for probbar_s = ',time,' sec.'
2     ecm = ecmstore
      return 
      end
c
      real*8 function gfunc(xi,w,ierr)
      implicit real*8(a-h,o-z)
      complex*16 ai,gncmplx
      common/vdr/vdr,accur
      acc = accur/10.d0
      ai=dcmplx(0.,1.)
      x = -w/sqrt(2.d0)
      y = 2.d0*xi/sqrt(2.d0)
      ierr = 0
      n=0
      sp=abs(1.2533141373155d0*gncmplx(x,y,n,acc,ierr))
      if(ierr.eq.1)then
       gfunc = 0.d0
       ierr = 1
       return
      endif
      gfunc=exp(-w**2)*sp**2
      ierr = 0
      return
      end
c
      real*8 function efunc(xi,w,ierr)
      implicit real*8(a-h,o-z)
      complex*16 ai,gncmplx,spp
      common/vdr/vdr,accur
      ai=dcmplx(0.d0,1.d0)
      x = -w/sqrt(2.d0)
      y = 2.d0*xi/sqrt(2.d0)
      ierr = 0
      n=0
      spp = gncmplx(x,y,n,accur,ierr)
      if(ierr.eq.1)then
       efunc = 0.d0
       ierr = 1
       return
      endif
      sp=abs(1.2533141373155d0*spp+ai/2.d0/xi)
      efunc=exp(-w**2)*(2.d0*xi*sp**2-1.d0/(2.d0*xi))
      ierr = 0
      return
      end
c
      real*8 function sfunc(xi,w,ierr)
      implicit real*8(a-h,o-z)
      complex*16 ai,gncmplx
      common/vdr/vdr,accur
      ai=dcmplx(0.d0,1.d0)
      x = -w/sqrt(2.d0)
      y = 2.d0*xi/sqrt(2.d0)
      ierr = 0
      n = 0
      sp=abs(1.2533141373155d0*gncmplx(x,y,n,accur,ierr)+ai/2.d0/xi)
      if(ierr.eq.1)then
       sfunc = 0.d0
       ierr = 1
       return
      endif
      sfunc=exp(-w**2/2.d0)*abs(2.d0*xi)*sp
      ierr = 0
      return
      end
