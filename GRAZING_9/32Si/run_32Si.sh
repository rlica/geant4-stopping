#!/bin/bash

export GRAZING_DIR=/mnt/c/Users/licar/cernbox/DAQ/GRAZING_9/grazing_exe
EXE=/mnt/c/Users/licar/cernbox/DAQ/GRAZING_9/grazing_exe/grazing_9_64bit

OUTFILE="cross_sections_18O_30Si_2n_Estar_5p5_6p5.txt"

echo "# 18O + 30Si -> 16O + 32Si" > $OUTFILE
echo "# 2n transfer integrated in E* window 5.5-6.5 MeV" >> $OUTFILE
echo "# Elab(MeV)   total_2n_mb   sigma_Estar_5p5_6p5_mb" >> $OUTFILE

for energy in $(seq 20 0.2 26)
do
    mkdir -p "E_${energy}"
    cd "E_${energy}"

    $EXE << fine > screen_output.txt
OK
8 18 14 30
a
0.00001
no
1
0
$energy
4
0 0
lab
4
yes
no
no
no
1
-2 0
1. -15.
sr
no
no
fine

    TOTAL=$(awk '
    /for deln =  *-2 and delz =  *0/ {flag=1}
    flag && /Total cross section/ {print $5; exit}
    ' fort.98)

    SIG_WINDOW=$(awk '
    /^[ ]*[-+]?[0-9]/ {
        E=$1
        sigma=$2
        gsub("D","E",sigma)
        if (E >= 5.5 && E <= 6.5) {
            sum += sigma
        }
    }
    END {print sum}
    ' fort.95)

    echo "$energy   $TOTAL   $SIG_WINDOW" >> "../$OUTFILE"

    cd ..
done
