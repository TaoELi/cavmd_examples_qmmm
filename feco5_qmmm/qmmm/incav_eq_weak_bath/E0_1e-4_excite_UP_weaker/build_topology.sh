#!/bin/bash

NAT_IMPURITY=11
NAT_SOLVENT=38

#geom_pos=$(cat geom.txt)
geom_pos="$1"

#echo "$geom_pos"

N_LINE=$(echo "$geom_pos" | wc -l)
#echo "This geom file contains $N_LINE lines"

N_SOLVENT=$((($N_LINE - $NAT_IMPURITY) / $NAT_SOLVENT))
#echo "This geom file contains $N_SOLVENT solvent molecules"

MM_TOPOLOGY_IMPURITY="   -1  0   0   0   0
   -2  0   0   0   0 
   -2  0   0   0   0
   -2  0   0   0   0
   -2  0   0   0   0
   -2  0   0   0   0
   -3  0   0   0   0
   -3  0   0   0   0 
   -3  0   0   0   0
   -3  0   0   0   0 
   -3  0   0   0   0
"

MM_TOPOLOGY_SOLVENT=""
for i in $(seq 0 $(($N_SOLVENT-1)))
do
    ntmp=$(($NAT_SOLVENT*i+$NAT_IMPURITY))
    MM_TOPOLOGY_SOLVENT="$MM_TOPOLOGY_SOLVENT    6     $(($ntmp+2))    0    0    0 
    1     $(($ntmp+1))    $(($ntmp+3))    $(($ntmp+4))    $(($ntmp+5))
    6     $(($ntmp+2))    0    0    0
    6     $(($ntmp+2))    0    0    0
    2     $(($ntmp+2))    $(($ntmp+6))    $(($ntmp+7))    $(($ntmp+8))
    6     $(($ntmp+5))    0    0    0
    6     $(($ntmp+5))    0    0    0
    2     $(($ntmp+5))    $(($ntmp+9))   $(($ntmp+10))   $(($ntmp+11))
    6     $(($ntmp+8))    0    0    0
    6     $(($ntmp+8))    0    0    0
    2     $(($ntmp+8))   $(($ntmp+12))   $(($ntmp+13))   $(($ntmp+14))
    6    $(($ntmp+11))    0    0    0
    6    $(($ntmp+11))    0    0    0
    2    $(($ntmp+11))   $(($ntmp+15))   $(($ntmp+16))   $(($ntmp+17))
    6    $(($ntmp+14))    0    0    0
    6    $(($ntmp+14))    0    0    0
    2    $(($ntmp+14))   $(($ntmp+18))   $(($ntmp+19))   $(($ntmp+20))
    6    $(($ntmp+17))    0    0    0
    6    $(($ntmp+17))    0    0    0
    2    $(($ntmp+17))   $(($ntmp+21))   $(($ntmp+22))   $(($ntmp+23))
    6    $(($ntmp+20))    0    0    0
    6    $(($ntmp+20))    0    0    0
    2    $(($ntmp+20))   $(($ntmp+24))   $(($ntmp+25))   $(($ntmp+26))
    6    $(($ntmp+23))    0    0    0
    6    $(($ntmp+23))    0    0    0
    2    $(($ntmp+23))   $(($ntmp+27))   $(($ntmp+28))   $(($ntmp+29))
    6    $(($ntmp+26))    0    0    0
    6    $(($ntmp+26))    0    0    0
    2    $(($ntmp+26))   $(($ntmp+30))   $(($ntmp+31))   $(($ntmp+32))
    6    $(($ntmp+29))    0    0    0
    6    $(($ntmp+29))    0    0    0
    2    $(($ntmp+29))   $(($ntmp+33))   $(($ntmp+34))   $(($ntmp+35))
    6    $(($ntmp+32))    0    0    0
    6    $(($ntmp+32))    0    0    0
    1    $(($ntmp+32))   $(($ntmp+36))   $(($ntmp+37))   $(($ntmp+38))
    6    $(($ntmp+35))    0    0    0
    6    $(($ntmp+35))    0    0    0
    6    $(($ntmp+35))    0    0    0
"
done

MM_TOPOLOGY="$MM_TOPOLOGY_IMPURITY $MM_TOPOLOGY_SOLVENT"

#echo "$MM_TOPOLOGY"

geom_pos=$(paste -d' ' <(echo "$geom_pos") <(echo "$MM_TOPOLOGY"))

# finally append qm_atoms list and force field params to the topology


echo "$geom_pos"
