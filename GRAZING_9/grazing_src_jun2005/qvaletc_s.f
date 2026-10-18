c      Reads and calculates Q-value and BE2,BE3 for schematic model          
c      Uses subroutines mass and be23
c      iopt = 1,2,3,4 indicates the option for the Q-values
c      If at entry iopt = 0 option 3 is chosen in program
c      At entry ierr=istand=1 if parameters are to be changed
c      If at entry inelas = -1 no inelastic channels are considered
c           if inelas = 0 only giant modes are included
c           if inelas = 1(3) only low-lying modes are included
c           if inelas = 2(4)both low- and high-lying modes are included
c         1 or 2 when inelastic parameters for low modes are to be found 
c                  from the be23.dat file or otherwise asked for
c         3 or 4 when inelastic parameters are found in be23.dat
c                  if available, otherwise skipped
c         in both cases ierr = 3 at exit if no data are found
c-------------------------------------------------------------

      subroutine qvaletc_s(iopt,iout,inelas,ierr)
      implicit real*8(a-h,o-z)
c      common:akyuz
      common/za/za,aa,zaa,aaa
      common/iza/iza,iaa,izaa,iaaa
      common/constant/pi,hbarc,c2m,esq,hc2om
      common/parm/ra,raa,aa3,aaa3,r0,v0,ainv
      common/mml/ssl,qvalcom,threshn,threshp
      common/emas/delans,delanp,delaps,delapp,eneut,eprot
     &,delaans,delaanp,delaaps,delaapp,corma,cormaa
      common/q/qns,qnp,qps,qpp
      common/cap/capn,capp,rr
      common/ba/bna,bnaa,bpa,bpaa
      common/ks/akbarn,akbarp
      common/ebe23/eex(2,3,2)
      common/be23/be23(2,3,2)
      common/giantwid/giantwid(3)
      common/zp23/zp(2,3,2)
      common/CD/Cin(2,3,2),Din(2,3,2)
      common/capt/eb,rbar,i5,intp
      common/g_0/g0na,g0naa,g0pa,g0paa
      common/salt/amml,cckn,cckp,cden,cdep,r00
      common/char/elemp,elemt
      dimension efac(3),sfac(3)
      character*5 elemp
      character*5 elemt
      character*5 intp
      character*5 version
      common/vers/version

c
      if(version.eq.'vms  '.or.version.eq.'unix* ')then
       open (unit=98,file='grazing_out.dat',status='old')
      else
      endif
      istand = ierr
c  calculate Q
c  first neutron and proton masses
      iz=0
      ia=1
      ierr1 = 0
      call masse_s(iz,ia,eneut,dd1,ierr1)
      if(ierr1.eq.1)goto 33
      iz=1
      ia=1
      ierr1 = 0
      call masse_s(iz,ia,eprot,dd1,ierr1)
      if(ierr1.eq.1)goto 33
c   then initial masses:
      iz1 = iza
      ia1 = iaa
      iz2 = izaa
      ia2 = iaaa
      imaf = 0
      imaaf = 0
      ierr1 = 2
      call masse_s(iz1,ia1,ea2,d1,ierr1)
      if(ierr1.eq.1)imaf = 3
      ierr1 = 3
      call masse_s(iz1,ia1,ea3,d1,ierr1)
      if(ierr1.eq.1.and.imaf.eq.3)goto 33
      ierr2 = 2
      call masse_s(iz2,ia2,eaa2,d2,ierr2)
      if(ierr2.eq.1)imaaf = 3
      ierr2 = 3
      call masse_s(iz2,ia2,eaa3,d2,ierr2)
      if(ierr2.eq.1.and.imaaf.eq.3)goto 33
c
c  Q-value for compound
      izc = iza+izaa
      iac = iaa+iaaa
      imcf = 0
      ierr6 = 2
      call masse_s(izc,iac,ec2,dc,ierr6)
      if(ierr6.eq.1)imcf = 3
      ierr6 = 3
      call masse_s(izc,iac,ec3,dc,ierr6)
      if(ierr6.eq.1.and.imcf.eq.3)then
       write(*,*)' Could not find mass of compound nucleus'
       qualcom = 0.d0
      elseif(imcf.eq.3.or.imaf.eq.3.or.imaaf.eq.3)then
       qvalcom = ea3 + eaa3 - ec3
       imcf = 3
      else
       qvalcom = ea2 + eaa2 - ec2
      endif
c  evaporation theshhold
      iz=izc
      ia=iac
      ierr1 = imcf
      call masse_s(iz,ia-1,emns,dd1,ierr1)
      if(ierr1.eq.1.and.imcf.eq.3)then
       write(*,*)' Could not find mass of compound nucleus - N'
       threshn = 0.d0
      elseif(ierr1.eq.3.or.imcf.eq.3)then
       threshn = emns + eneut - ec3
      else
       threshn = emns + eneut - ec2
      endif
      ierr1 = imcf
      call masse_s(iz-1,ia-1,emps,dd1,ierr1)
      if(ierr1.eq.1.and.imcf.eq.3)then
       write(*,*)' Could not find mass of compound nucleus - P'
       threshn = 0.d0
      elseif(ierr1.eq.3.or.imcf.eq.3)then
       threshp = emps + eprot - ec3
      else
       threshp = emps + eprot - ec2
      endif
c
c  single neutron transfer Q-values
      ins=1                                              
13    iz3 = iza                                
      ia3 = iaa - ins
      iz4 = izaa 
      ia4 = iaaa + ins
      ierr3 = imaf
      call masse_s(iz3,ia3,e3,d3,ierr3)
      if(ierr3.eq.1)goto 33
      ierr4 = imaaf
      call masse_s(iz4,ia4,e4,d4,ierr4)
      if(ierr4.eq.1)goto 33
      if (ins.eq.-1) goto 14
      ea = ea2
      if(ierr3.eq.3.or.imaf.eq.3)ea = ea3
      eaa = eaa2
      if(ierr4.eq.3.or.imaaf.eq.3)eaa = eaa3
      dans = ea-e3
      bindan = -dans+eneut
      if(iout.eq.1)write(*,'(a,f10.5)')
     &' Binding energy of neutron in projectile = ',bindan 
      daans = eaa-e4
      bindaa1n = daans+eneut
      if(iout.eq.1)write(*,'(a,f10.5)')
     &' Binding energy of neutron in target +1 neutron = ',bindaa1n 
      qqns  = ea+eaa-e3-e4
      write(98,'(a,f7.2,a)')
     &'  Q-value for neutron stripping = ', qqns,' MeV'
      write(*,'(a,f7.2,a)')
     & '  Q-value for neutron stripping = ', qqns,' MeV'
      if(imaf.eq.3.or.imaaf.eq.3.or.ierr3.eq.3.or.ierr4.eq.3)
     &write(98,*)'      Calculated values were used'
      ins = -1
      goto 13 
14    ea = ea2
      if(ierr3.eq.3.or.imaf.eq.3)ea = ea3
      eaa = eaa2
      if(ierr4.eq.3.or.imaaf.eq.3)eaa = eaa3
      danp = ea-e3
      binda1n = danp+eneut
      if(iout.eq.1)write(*,'(a,f10.5)')
     &' Binding energy of neutron in projectile + 1neutron = ',binda1n 
      daanp = eaa-e4
      bindaan = -daanp+eneut
      if(iout.eq.1)write(*,'(a,f10.5)')
     &' Binding energy of neutron in target = ',bindaan 
      qqnp  = ea+eaa-e3-e4
      write(98,'(a,f7.2,a)')
     & '  Q-value for neutron pick-up =   ', qqnp,' MeV' 
      write(*,'(a,f7.2,a)')
     & '  Q-value for neutron pick-up =   ', qqnp,' MeV' 
      if(imaf.eq.3.or.imaaf.eq.3.or.ierr3.eq.3.or.ierr4.eq.3)
     &write(98,*)'      Calculated values were used'
c
c  single proton transfer Q-values
      ips = 1                                  
11    iz3 = iza - ips                          
      ia3 = iaa - ips
      iz4 = izaa + ips
      ia4 = iaaa + ips
      ierr3 = imaf
      call masse_s(iz3,ia3,e3,d3,ierr3)
      if(ierr3.eq.1)goto 33
      ierr4 = imaaf
      call masse_s(iz4,ia4,e4,d4,ierr4)
      if(ierr4.eq.1)goto 33
      if(ips.eq.-1)goto 12
      ea = ea2
      if(ierr3.eq.3.or.imaf.eq.3)ea = ea3
      eaa = eaa2
      if(ierr4.eq.3.or.imaaf.eq.3)eaa = eaa3
      daps = ea-e3
      bindap = -daps+eprot
      if(iout.eq.1)write(*,'(a,f10.5)')
     &' Binding energy of proton in projectile = ',bindap 
      daaps = eaa-e4
      bindaa1p = daaps+eprot
      if(iout.eq.1)write(*,'(a,f10.5)')
     &' Binding energy of proton in target +1 proton = ',bindaa1p 
      qqps  = ea+eaa-e3-e4
      write(98,'(a,f7.2,a)') 
     &'  Q-value for proton stripping =  ', qqps,' MeV'
      write(*,'(a,f7.2,a)') 
     &'  Q-value for proton stripping =  ', qqps,' MeV'
      if(imaf.eq.3.or.imaaf.eq.3.or.ierr3.eq.3.or.ierr4.eq.3)
     &write(98,*)'      Calculated values were used'
      ips = -1
      goto 11
12    ea = ea2
      if(ierr3.eq.3.or.imaf.eq.3)ea = ea3
      eaa = eaa2
      if(ierr4.eq.3.or.imaaf.eq.3)eaa = eaa3
      dapp = ea-e3
      binda1p = dapp+eprot
      if(iout.eq.1)write(*,'(a,f10.5)')
     &' Binding energy of proton in projectile + 1 proton  = ',binda1p 
      daapp = eaa-e4
      bindaap = -daapp+eprot
      if(iout.eq.1)write(*,'(a,f10.5)')
     &' Binding energy of proton in target = ',bindaap 
      qqpp  = ea+eaa-e3-e4
      write(98,'(a,f7.2,a)') 
     &'  Q-value for proton pick-up =    ', qqpp,' MeV'
      write(*,'(a,f7.2,a)') 
     &'  Q-value for proton pick-up =    ', qqpp,' MeV'
      if(imaf.eq.3.or.imaaf.eq.3.or.ierr3.eq.3.or.ierr4.eq.3)
     &write(98,*)'      Calculated values were used'
c
      fermian = -(bindan+binda1n)/2.d0
      fermiaan = -(bindaan+bindaa1n)/2.d0
      fermiap = -(bindap+binda1p)/2.d0
      fermiaap = -(bindaap+bindaa1p)/2.d0
      corma = ea+(fermian+bindan+fermiap+bindap)/2.d0
      cormaa = eaa+(fermiaan+bindaan+fermiaap+bindaap)/2.d0
c
c  two neutron transfer Q-values
      if(iopt.eq.1.or.iopt.eq.3)then
       ins=2                                              
513    iz3 = iza                                
       ia3 = iaa - ins
       iz4 = izaa 
       ia4 = iaaa + ins
       ierr3 = imaf
       call masse_s(iz3,ia3,e3,d3,ierr3)
       if(ierr3.eq.1)goto 33
       ierr4 = imaaf
       call masse_s(iz4,ia4,e4,d4,ierr4)
       if(ierr4.eq.1)goto 33
       if(ins.eq.-2)goto 514
       ea = ea2
       if(ierr3.eq.3.or.imaf.eq.3)ea = ea3
       eaa = eaa2
       if(ierr4.eq.3.or.imaaf.eq.3)eaa = eaa3
       delans = (ea-e3)/2.d0
       delaans = (eaa-e4)/2.d0
       aqqns  = (ea+eaa-e3-e4)/2.d0
       if(iout.eq.1)write(98,'(a,f7.2,a)')
     &' Average Q-value for neutron stripping ',aqqns,' MeV'
       if(iout.eq.1)write(*,'(a,f7.2,a)')
     &' Average Q-value for neutron stripping ',aqqns,' MeV'
       if(ierr1.eq.3.or.ierr2.eq.3.or.ierr3.eq.3.or.ierr4.eq.3)
     &write(98,*)'      Calculated values were used'
       ins = -2
       goto 513 
514    ea = ea2
       if(ierr3.eq.3.or.imaf.eq.3)ea = ea3
       eaa = eaa2
       if(ierr4.eq.3.or.imaaf.eq.3)eaa = eaa3
       delanp = (ea-e3)/2.d0
       delaanp = (eaa-e4)/2.d0
       aqqnp  = (ea+eaa-e3-e4)/2.d0
       if(iout.eq.1)write(98,'(a,f7.2,a)')
     &' Average Q-value for neutron pick-up ', aqqnp,' MeV' 
       if(iout.eq.1)write(*,'(a,f7.2,a)') 
     &' Average Q-value for neutron pick-up ', aqqnp,' MeV' 
       if(ierr1.eq.3.or.ierr2.eq.3.or.ierr3.eq.3.or.ierr4.eq.3)
     &write(98,*)'      Calculated values were used'
c
c  two proton transfer Q-values
       ips = 2                                  
511    iz3 = iza - ips                          
       ia3 = iaa - ips
       iz4 = izaa + ips
       ia4 = iaaa + ips
       ierr3 = imaf
       call masse_s(iz3,ia3,e3,d3,ierr3)
       if(ierr3.eq.1)goto 33
       ierr4 = imaaf
       call masse_s(iz4,ia4,e4,d4,ierr4)
       if(ierr4.eq.1)goto 33
       if(ips.eq.-2)goto 512
       ea = ea2
       if(ierr3.eq.3.or.imaf.eq.3)ea = ea3
       eaa = eaa2
       if(ierr4.eq.3.or.imaaf.eq.3)eaa = eaa3
       delaps = (ea-e3)/2.d0
       delaaps = (eaa-e4)/2.d0
       aqqps  = (ea+eaa-e3-e4)/2.d0
       if(iout.eq.1)write(98,'(a,f7.2,a)') 
     &' Average Q-value for proton stripping ', aqqps,' MeV'
       if(iout.eq.1)write(*,'(a,f7.2,a)') 
     &' Average Q-value for proton stripping ', aqqps,' MeV'
       if(imaf.eq.3.or.imaaf.eq.3.or.ierr3.eq.3.or.ierr4.eq.3)
     &write(98,*)'      Calculated values were used'
       ips = -2
       goto 511
512    ea = ea2
       if(ierr3.eq.3.or.imaf.eq.3)ea = ea3
       eaa = eaa2
       if(ierr4.eq.3.or.imaaf.eq.3)eaa = eaa3
       delapp = (ea-e3)/2.d0
       delaapp = (eaa-e4)/2.d0
       aqqpp  = (ea+eaa-e3-e4)/2.d0
       if(iout.eq.1)write(98,'(a,f7.2,a)') 
     &' Average Q-value for proton pick-up ', aqqpp,' MeV'
       if(iout.eq.1)write(*,'(a,f7.2,a)') 
     &' Average Q-value for proton pick-up ', aqqpp,' MeV'
       if(ierr1.eq.3.or.ierr2.eq.3.or.ierr3.eq.3.or.ierr4.eq.3)
     &write(98,*)'      Calculated values were used'
      endif
c
c  Slope of densities and formfactor-radius
      rra = r00*aa3
      rraa = r00*aaa3
      rr = rra + rraa
      capn = 1.16d0 + .5d0*(1.d0/aa3+1.d0/aaa3)
      capp = 1.65d0 + .1d0*(1.d0/aa3+1.d0/aaa3)
      if(iout.eq.1)then
       write(98,'(a,f6.3,a)')
     &' Average exponential slope of neutron density =',capn,' 1/fm'
       write(98,'(a,f6.3,a)')
     &' Average exponential slope of proton density =',capp,' 1/fm'
c
       write(*,'(a,f10.5,a)')
     &' Fermi level for neutrons in projectile = ',fermian,' MeV'
       write(*,'(a,f10.5,a)')
     &' Fermi level for protons in projectile = ',fermiap,' MeV'
       write(*,'(a,f10.5,a)')
     &' Fermi level for neutrons in target = ',fermiaan,' MeV'
       write(*,'(a,f10.5,a)')
     &' Fermi level for protons in target = ',fermiaap,' MeV'
       write(98,'(a,f10.5,a)')
     &' Fermi level for neutrons in projectile = ',fermian,' MeV'
       write(98,'(a,f10.5,a)')
     &' Fermi level for protons in projectile = ',fermiap,' MeV'
       write(98,'(a,f10.5,a)')
     &' Fermi level for neutrons in target = ',fermiaan,' MeV'
       write(98,'(a,f10.5,a)')
     &' Fermi level for protons in target = ',fermiaap,' MeV'
       write(*,'(a,f10.5,a)')
     &' Corrected mass of projectile = ',corma,' MeV'
       write(*,'(a,f10.5,a)')
     &' Corrected mass of target = ',cormaa,' MeV'
      endif
c
c      Level densities
      epsna = 52.d0 - 22.d0*(aa - 2.d0*za)/aa -8.d0
      epsnaa = 52.d0 - 22.d0*(aaa - 2.d0*zaa)/aaa - 8.d0
      epspa = 52.d0+22.d0*(aa-2.d0*za)/aa-za*esq/(1.2d0*aa3)-8.d0
      epspaa = 52.d0+22.d0*(aaa-2.d0*zaa)/aaa-zaa*esq/(1.2d0*aaa3)-8.d0
      if(iout.eq.1)then
       write(98,'(a,f5.1,a/a,f5.1,a)')
     &' Fermi energies in projectile, for neutrons',epsna,' MeV',
     &'                            and for protons',epspa,' MeV'
       write(98,'(a,f5.1,a/a,f5.1,a)')
     &' Fermi energies in target, for neutrons',epsnaa,' MeV',
     &'                        and for protons',epspaa,' MeV'
      endif
      g0na=1.5d0*(aa-za)/epsna*(15.d0/cden)
      g0naa=1.5d0*(aaa-zaa)/epsnaa*(15.d0/cden)      
      g0pa=1.5d0*za/epspa*(15.d0/cdep)
      g0paa=1.5d0*zaa/epspaa*(15.d0/cdep)
      if(iout.eq.1)then
       write(98,'(a,2f5.2,a)')
     &' Level density for neutrons in a and A',g0na,g0naa,' 1/MeV'
       write(98,'(a,2f5.2,a)')
     &' Level density for protons  in a and A',g0pa,g0paa,' 1/MeV'
      endif
c
c  Moments of inertia/hbar^2 
      airiga=.4d0*aa/hc2om*(1.2d0*aa3)**2
      airigaa=.4d0*aaa/hc2om*(1.2d0*aaa3)**2      
c  Parameters b for <m**2> and S_l
      bna=(aa-za)/aa*airiga/g0na
      bnaa=(aaa-zaa)/aaa*airigaa/g0naa
      bpa=za/aa*airiga/g0pa
      bpaa=zaa/aaa*airigaa/g0paa 
      ssl=amml/(2.d0*(sqrt(bna)+sqrt(bnaa)))              
      if(iout.eq.1)then
       write(98,'(a,4f5.1)')
     &' <m**2> (bn_a,bp_a,bn_A,bp_A) for a and A =',bna,bpa,bnaa,bpaa
       write(98,'(a,f4.3)')' S_l =' ,ssl
      endif
c
c  ksn and ksp parameters
      ckn = cckn*.25d0
      ckp = cckp*.3d0
      akbarn = pi**1.5d0/4.d0*g0na*g0naa*(hc2om/2.d0)**2*(1.6d0/aa)*
     &(1.6d0/aaa)*(rra*rraa/(rra+rraa))**2*ckn
      akbarp = pi**1.5d0/4.d0*g0pa*g0paa*(hc2om/2.d0)**2*(1./aa)*
     &(1./aaa)*(rra*rraa/(rra+rraa))**2*ckp
      if(iout.eq.1)then
       write(98,'(a,2f5.2)')
     &' kbar coefficients for neutrons and protons',akbarn,akbarp
      endif
c
c  corrected Q-value(4 options)
      if(iopt.eq.1)goto 60
      if(iopt.eq.2)goto 61
      if(iopt.eq.3)goto 62
      if(iopt.eq.4)goto 63
60    qns = aqqns
      qnp = aqqnp
      qps = aqqps
      qpp = aqqpp
      goto 550
61    qns = fermian-fermiaan
      qnp = -qns
      qps = fermiap-fermiaap
      qpp = - qps
      delans = fermian+eneut
      delanp = -delans
      delaans = -fermiaan-eneut
      delaanp = -delaans
      delaps = fermiap+eprot
      delapp = -delaps
      delaaps = -fermiaap-eprot
      delaapp = -delaaps
      goto 550
62    qns = (.5d0*(qqns-qqnp)+aqqns)/2.d0
      qnp = (.5d0*(qqnp-qqns)+aqqnp)/2.d0
      qps = (.5d0*(qqps - qqpp)+aqqps)/2.d0
      qpp = (.5d0*(qqpp - qqps)+aqqpp)/2.d0
      delans = ((dans - danp)/2.d0+delans)/2.d0
      delanp = ((danp - dans)/2.d0+delanp)/2.d0
      delaans = ((daans - daanp)/2.d0+delaans)/2.d0
      delaanp = ((daanp - daans)/2.d0+delaanp)/2.d0
      delaps = ((daps - dapp)/2.d0+delaps)/2.d0
      delapp = ((dapp - daps)/2.d0+delapp)/2.d0
      delaaps = ((daaps - daapp)/2.d0+delaaps)/2.d0
      delaapp = ((daapp - daaps)/2.d0+delaapp)/2.d0
      goto 550
63    qns = qqns
      qnp = qqnp
      qps = qqps
      qpp = qqpp
      delans = dans
      delanp = danp
      delaans = daans
      delaanp = daanp
      delaps = daps 
      delapp =  dapp
      delaaps = daaps
      delaapp =  daapp
      corma = ea
      cormaa = eaa
550   continue
      write(98,'(a,4f7.3,a)')
     &' Q-values used: qns,qnp,qps,qpp =',qns,qnp,qps,qpp,' MeV'
      write(*,'(a,4f7.3,a)')
     &' Q-values used: qns,qnp,qps,qpp =',qns,qnp,qps,qpp,' MeV'
c
c   Inelastic parameters
      do i=1,2
       do lam=2,3
        do n=1,2
         eex(i,lam,n) = 0.d0
        enddo
       enddo
      enddo
68    if(inelas.eq.-1)then
       write (*,'(a/a)')' You are neglecting inelastic scattering',
     &'    Capture unrealistic'
       write (98,'(a/a)')' You are neglecting inelastic scattering',
     &'    Capture unrealistic'
       goto 50
      endif
      if(inelas.eq.1.or.inelas.eq.3)then
       write (*,'(a)')' You are neglecting giant resonances'
       write (98,'(a)')' You are neglecting giant resonances'
       goto 51 
      endif
      if(inelas.eq.0)then
       write (*,'(a/a)')' You are neglecting low-lying states',
     &'    Capture unrealistic'
       write (98,'(a/a)')' You are neglecting low-lying states',
     &'    Capture unrealistic'
      endif
c  default parameters for giant resonances  (or read)
      efac(2) = 57.d0
      efac(3) = 94.d0
      sfac(2) = .8d0
      sfac(3) = .4d0
      giantwid(2) = 2.5d0/2.3548
      giantwid(3) = 6.d0/2.3548
      if(istand.eq.1)then
168    if(intp.eq.'inter')write(*,'(6(/a))')
     &' The parameters for the giant resonances are calculated from:',
     &'   Isoscalar quadrupole: energy = 57/A**(1/3),',
     &'     width(FWHM) = 2.5 MeV, strength = 80% of sumrule,',
     &'   Isoscalar octupole: energy = 94/A**(1/3)',
     &'     width(FWHM) = 6.0 MeV, strength = 40% of sumrule,',
     &' You may change this by entering new factors for the energy',
     &'   i.e. Enter: 57, 94 if no change'
       read(i5,*,err=1100)efac(2),efac(3)
       icheck = 0
       goto 1110
1100   write(*,*)' you must provide two numbers'
       icheck = icheck + 1
       if(icheck.gt.5)goto 33
       goto 168
1110   if(intp.eq.'inter')write(10,*)efac(2),efac(3)
268    if(intp.eq.'inter')write(*,'(/a/a)')
     &'  and new fractions of the sum rule',
     &'    i.e. Enter: .8, .4 if no change'
       read(i5,*,err=1101)sfac(2),sfac(3)
       icheck = 0
       goto 1111
1101   write(*,*)' you must provide two numbers'
       icheck = icheck + 1
       if(icheck.gt.5)goto 33
       goto 268
1111   if(intp.eq.'inter')write(10,*)sfac(2),sfac(3)
368    if(intp.eq.'inter')write(*,'(/a/a)')
     &'  and new widths',
     &'    i.e. Enter: 2.5, 6.0 if no change'
       read(i5,*,err=1102)gwidth2,gwidth3
       icheck = 0
       goto 1112
1102   write(*,*)' you must provide two numbers'
       icheck = icheck + 1
       if(icheck.gt.5)goto 33
       goto 368
1112   if(intp.eq.'inter')write(10,*)gwidth2,gwidth3
       giantwid(2) = gwidth2/2.3548
       giantwid(3) = gwidth3/2.3548
      endif
c  calculate collective  modes
51    icount = 0
      do i = 1,2
c  (i=1=projectile 2=target)
       ncount = 1
       if(i.eq.1)then
        zint = za
        aaint3 = aa3
        write (98,'(a,a)') ' Collective states for ',elemp
        write (*,'(a,a)') ' Collective states for ',elemp
       elseif(i.eq.2)then
        zint = zaa
        aaint3 = aaa3
        write (98,'(a,a)') ' Collective states for ',elemt
        write (*,'(a,a)') ' Collective states for ',elemt
       endif
       do lam = 2,3
        if(inelas.eq.1.or.inelas.eq.3)goto 52
c  giant resonances
        slam = dble(3*(2*lam+1)*lam)/8.d0/pi*hc2om*aaint3**(2*lam-5)
     &*(1.20d0)**(2*lam-2)*zint**2
        dirrot = .75d0/pi/dble(lam)*aaint3**5*(1.20d0)**2/hc2om        
        eex(i,lam,2)= efac(lam)/aaint3
        Din(i,lam,2)= dirrot/sfac(lam)
        Cin(i,lam,2)= Din(i,lam,2)*eex(i,lam,2)**2
        zp(i,lam,2)= sqrt(eex(i,lam,2)/2.d0/Cin(i,lam,2))
        be23(i,lam,2)= slam*sfac(lam)/eex(i,lam,2)/100.d0**lam
        gwid = giantwid(lam)*2.3548
        write (98,'(a,f7.3,a,i1,a,e12.3,a,i1,a,f5.2,a)') 
     &' Energy =',eex(i,lam,2),' MeV   BE',lam,' =',be23(i,lam,2),
     &' e^2 b^',lam,' width = ',gwid,' (giant)'
        write (*,'(a,f7.3,a,i1,a,e12.3,a,i1,a,f5.2,a)') 
     &' Energy =',eex(i,lam,2),' MeV   BE',lam,' =',be23(i,lam,2),
     &' e^2 b^',lam,' width = ',gwid,' (giant)'
        if(inelas.eq.0)goto 21
52      if(inelas.ge.3)then
         ierr1 = 1
        else
         ierr1 = 0
        endif
c  Low-lying modes
        if(i.eq.1)then
         iz = iza
         ia = iaa
        elseif(i.eq.2)then
         iz = izaa
         ia = iaaa
        endif
        call be23c_s(iz,ia,lam,eex(i,lam,1),be23(i,lam,1)
     &,Cin(i,lam,1),Din(i,lam,1),zp(i,lam,1),ierr1)
        if(ierr1.eq.2)then
         ncount = ncount+1
         goto 21
        elseif(ierr1.eq.1)then
         goto 33
        endif
        if(ierr1.eq.3)then
         write(*,*)' Interpolated data:'
         write(98,*)' Interpolated data:'
        endif
        write (98,'(a,f7.3,a,i1,a,e12.3,a,i1)') 
     &' Energy =',eex(i,lam,1),' MeV   BE',lam,' =',be23(i,lam,1),
     &' e^2 b^',lam
        write (*,'(a,f7.3,a,i1,a,e12.3,a,i1)') 
     &' Energy =',eex(i,lam,1),' MeV   BE',lam,' =',be23(i,lam,1),
     &' e^2 b^',lam
21      continue
       enddo
       if(ncount.eq.3)then
        write(*,*)' Low-lying states are not included'
        write(98,*)' Low-lying states are not included'
        icount = icount+1
       endif
      enddo
      if(icount.eq.2)then
       write(*,*)' No low-lying data were found'
       write(98,*)' No low-lying data were found'
       ierr = 3
       return
      endif
      if(iout.eq.1)then
       write(98,*)' Parameters for the collective states included:'
       do i=1,2
        do lam= 2,3
         do n=1,2
          if(zp(i,lam,n).gt.1.d-8)then
           if(mod(lam,2).eq.0)
     &write(98,'(a,i1,a,f5.2,a,e10.3,a,e10.3,a,f5.4)')
     &' ',lam,'+: E = ',eex(i,lam,n),' C = ',Cin(i,lam,n),' D = ',
     &Din(i,lam,n),' zerop = ',zp(i,lam,n)
           if(mod(lam,2).eq.1)
     &write(98,'(a,i1,a,f5.2,a,e10.3,a,e10.3,a,f5.4)')
     &' ',lam,'-: E = ',eex(i,lam,n),' C = ',Cin(i,lam,n),' D = ',
     &Din(i,lam,n),' zerop = ',zp(i,lam,n)
          endif
         enddo
        enddo
       enddo
      endif
50    ierr = 0
      return
33    ierr = 1
      return
      end
