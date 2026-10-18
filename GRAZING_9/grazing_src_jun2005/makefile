###############################################################################
#
# the simple makefile for    grazing_9 
#
###############################################################################
#
### some macro-name (alias) definition
#
PROJECTS = grazing_9
#
##############################################################################
#  supported achitecture: LINUX, ALPHA
#  uncomment the appropriate block
##############################################################################
# LINUX (the following lines must be uncommended)
##############################################################################
   NAGLIB = -L/usr/local/lib -lnag18 
   OPTION = -O2 -finit-local-zero -fno-automatic -fno-backslash -fdollar-ok
   COMP   = gfortran
### end LINUX (more at the end) #############################################
#
#############################################################################
# ALPHA (the following lines must me uncommended)
#############################################################################
#   NAGLIB = -L/teosw/lib -lnag18	
#   OPTION =
#   COMP = f77
### end ALPHA (more at the end) #############################################
#
PROGRAM  = grazing_9.o                \
           qvaletc_s.o                \
           capturenew_s.o             \
           population_s.o             \
           probbarnew_s.o             \
           transmnew_s.o              \
           barrierpos_s.o             \
           nzdistri_s.o               \
           convert_s.o                \
           condens_s.o                \
           date_s.o                   \
           time_s.o                   \
           akyuz_lib.o                \
           angularmom_lib.o           \
           fys_lib.o                  \
           interpolate_lib.o          \
           math_lib.o                 \
           leg1new_s.o                \
           treslim_s.o                \
           scatangleavnew_s.o         \
           leg2new_s.o                \
           leg3new_s.o                \
           leg41new_s.o               \
           leg42new_s.o               \
           leg43new_s.o               \
           leg44new_s.o               \
           edista_s.o                 \
           tdistri_s.o                \
           edistapp_s.o               \
           driftnztdis_s.o            \
           mavtdist_s.o               \
           wilsh_s.o                  \
           wilsh_plotnew_s.o          \
           leg4_s.o                   \
           emass_s.o                  \
           enzcor_s.o                 \
           evap_s.o                   \
           nzedisa_s.o                \
           selectnz_s.o               \
           nztdis_s.o                 \
           mavnzdis_s.o               \
           wilsh_nzplotnew_s.o        \
           leg45new_s.o
#
### the main command for the makefile
#
grazing_9:  $(PROGRAM)     ;\
       $(COMP) $(OPTION) $(NAGLIB) $(PROGRAM) $(NAGLIB) -o grazing_9   
#
### the commands to produce all the obj
#
grazing_9.o: grazing_9.f              ;\
       $(COMP) -c $(OPTION) grazing_9.f
qvaletc_s.o:  qvaletc_s.f              ;\
       $(COMP) -c $(OPTION) qvaletc_s.f
capturenew_s.o: capturenew_s.f         ;\
       $(COMP) -c $(OPTION) capturenew_s.f
population_s.o:   population_s.f         ;\
       $(COMP) -c $(OPTION)  population_s.f
probbarnew_s.o:   probbarnew_s.f         ;\
       $(COMP) -c $(OPTION)  probbarnew_s.f
transmnew_s.o:  transmnew_s.f            ;\
       $(COMP) -c $(OPTION)  transmnew_s.f
barrierpos_s.o:  barrierpos_s.f          ;\
       $(COMP) -c $(OPTION)  barrierpos_s.f
nzdistri_s.o:    nzdistri_s.f            ;\
       $(COMP) -c $(OPTION) nzdistri_s.f
convert_s.o:   convert_s.f               ;\
       $(COMP) -c $(OPTION) convert_s.f
condens_s.o:   condens_s.f               ;\
       $(COMP) -c $(OPTION) condens_s.f
akyuz_lib.o:   akyuz_lib.f               ;\
       $(COMP) -c $(OPTION) akyuz_lib.f
angularmom_lib.o:  angularmom_lib.f      ;\
       $(COMP) -c $(OPTION) angularmom_lib.f
fys_lib.o:   fys_lib.f                   ;\
       $(COMP) -c $(OPTION) fys_lib.f
interpolate_lib.o: interpolate_lib.f     ;\
       $(COMP) -c $(OPTION)  interpolate_lib.f
math_lib.o:   math_lib.f                 ;\
       $(COMP) -c $(OPTION)  math_lib.f
leg1new_s.o:  leg1new_s.f                ;\
       $(COMP) -c $(OPTION)  leg1new_s.f
treslim_s.o:  treslim_s.f                ;\
       $(COMP) -c $(OPTION)  treslim_s.f
scatangleavnew_s.o:  scatangleavnew_s.f  ;\
       $(COMP) -c $(OPTION)  scatangleavnew_s.f
leg2new_s.o:   leg2new_s.f               ;\
       $(COMP) -c $(OPTION) leg2new_s.f
leg3new_s.o:  leg3new_s.f                ;\
       $(COMP) -c $(OPTION) leg3new_s.f
leg4_s.o:   leg4_s.f                     ;\
       $(COMP) -c $(OPTION) leg4_s.f
leg41new_s.o:   leg41new_s.f             ;\
       $(COMP) -c $(OPTION) leg41new_s.f
leg42new_s.o:    leg42new_s.f            ;\
       $(COMP) -c $(OPTION) leg42new_s.f
leg43new_s.o:  leg43new_s.f              ;\
       $(COMP) -c $(OPTION) leg43new_s.f
leg44new_s.o:   leg44new_s.f             ;\
       $(COMP) -c $(OPTION) leg44new_s.f
edista_s.o:   edista_s.f                 ;\
       $(COMP) -c $(OPTION) edista_s.f
tdistri_s.o:  tdistri_s.f                ;\
       $(COMP) -c $(OPTION) tdistri_s.f
edistapp_s.o:  edistapp_s.f              ;\
       $(COMP) -c $(OPTION) edistapp_s.f
driftnztdis_s.o:   driftnztdis_s.f       ;\
       $(COMP) -c $(OPTION) driftnztdis_s.f
mavtdist_s.o:  mavtdist_s.f              ;\
       $(COMP) -c $(OPTION) mavtdist_s.f
wilsh_s.o:    wilsh_s.f                  ;\
       $(COMP) -c $(OPTION) wilsh_s.f
wilsh_plotnew_s.o:  wilsh_plotnew_s.f    ;\
       $(COMP) -c $(OPTION) wilsh_plotnew_s.f
emass_s.o:   emass_s.f                   ;\
       $(COMP) -c $(OPTION) emass_s.f
enzcor_s.o:   enzcor_s.f                 ;\
       $(COMP) -c $(OPTION) enzcor_s.f
evap_s.o:   evap_s.f                     ;\
       $(COMP) -c $(OPTION) evap_s.f
nzedisa_s.o:   nzedisa_s.f               ;\
       $(COMP) -c $(OPTION) nzedisa_s.f
selectnz_s.o:    selectnz_s.f            ;\
       $(COMP) -c $(OPTION) selectnz_s.f
nztdis_s.o:   nztdis_s.f                 ;\
       $(COMP) -c $(OPTION) nztdis_s.f
mavnzdis_s.o:     mavnzdis_s.f           ;\
       $(COMP) -c $(OPTION) mavnzdis_s.f
wilsh_nzplotnew_s.o:  wilsh_nzplotnew_s.f  ;\
       $(COMP) -c $(OPTION) wilsh_nzplotnew_s.f
leg45new_s.o:     leg45new_s.f          ;\
       $(COMP) -c $(OPTION) leg45new_s.f
###############################################################################
# LINUX (the following lines must be uncommended) 
###############################################################################
   date_s.o:  date_s.f                      ;\
       $(COMP) -c $(OPTION)  date_s.f
   time_s.o:   time_s.f                     ;\
       $(COMP) -c $(OPTION) time_s.f
###############################################################################
# ALPHA (the following lines must be uncommended)
###############################################################################
#date_s.o:  date_axp_s.f                      ;\
#	$(COMP) -c $(OPTION)  date_axp_s.f -o date_s.o
#time_s.o:   time_s.f                     ;\
#	$(COMP) -c $(OPTION)  time_axp_s.f -o time_s.o
###############################################################################
# end of the makefile
#









