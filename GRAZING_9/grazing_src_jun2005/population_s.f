c     Subroutine tabulates (in FOR019) single particle population
c      in projectile (columns 1 and 2) and target (columns 3 and 4)
c      for neutrons (1 and 3 and protons (2 and 4).
c------------------------------------------------------------
      subroutine population_s(ierr)
      implicit real*8(a-h,o-z)
      common/pop/dtann,dtanp,dtaann,dtaanp
      common/ej/ejns,ejnp,ejps,ejpp
      common/bhom/bhomn,bhomp
      write(19,'(a)')
     &' energy above Fermi level, pop_an, pop_ap, pop_An, pop_Ap'
      do 20 n=0,20
       eps = -10.d0+n*.5d0
       arg1 = eps/bhomn-ejns
       arg2 = eps/bhomn-ejnp
       call ierfc(0,arg1,func1)
       call ierfc(0,arg2,func2)
       anan = 1.d0 - dtann*(2.d0-func1)
       anaan = 1.d0 - dtaann*(2.d0-func2) 
       arg1 = eps/bhomp-ejps
       arg2 = eps/bhomp-ejpp
       call ierfc(0,arg1,func1)
       call ierfc(0,arg2,func2)
       anap = 1.d0 - dtanp*(2.d0-func1)
       anaap = 1.d0 - dtaanp*(2.d0-func2) 
       write(19,'(f26.4,4f8.4)')eps,anan,anaan,anap,anaap
20    continue
      do 21 n=0,20
       eps =n*.5d0
       arg1 = eps/bhomn-ejns
       arg2 = eps/bhomn-ejnp
       call ierfc(0,arg1,func1)
       call ierfc(0,arg2,func2)
       anan = dtann*func1
       anaan = dtaann*func2
       arg1 = eps/bhomp-ejps
       arg2 = eps/bhomp-ejpp
       call ierfc(0,arg1,func1)
       call ierfc(0,arg2,func2)
       anap = dtanp*func1
       anaap = dtaanp*func2
       write(19,'(f26.4,4f8.4)')eps,anan,anaan,anap,anaap
21    continue
      ierr = 0
      return 
      end
