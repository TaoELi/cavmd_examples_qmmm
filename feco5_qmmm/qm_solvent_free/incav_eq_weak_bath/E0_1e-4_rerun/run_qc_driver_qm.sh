#!/bin/bash

BASIS="def2-msvp"
METHOD="bp86"
NTHREAD=2

# check whether ECHO has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

geom_pos="$1"
$ECHO "The current geometry position is"
$ECHO $geom_pos

cell_length="$2"
$ECHO "The primitive cell length is $cell_length"

do_dipder="$3" # try to do dipder calculation
$ECHO "Dipder calculation: $do_dipder"

IPI_DRIVER_TEMP="$4"
$ECHO "Q-Chem will run in $IPI_DRIVER_TEMP"

QCHEM=$HOME/qchem_dev/qchem/bin/qchem
TMP_DIR=./

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

$ECHO
$ECHO "$EXAMPLE_DIR : starting driver"
$ECHO

for DIR in "$EXAMPLE_DIR/$IPI_DRIVER_TEMP" ; do
    if test ! -d $DIR ; then
        mkdir $DIR
    fi
done
cd $EXAMPLE_DIR/$IPI_DRIVER_TEMP

# create Q-Chem input file for force
cat > ipi.scf.in << EOF
\$molecule
0   1
$geom_pos
\$end

\$rem
     JOBTYPE  force
     METHOD   $METHOD
     BASIS    $BASIS
     MAX_SCF_CYCLES   800
     INPUT_BOHR TRUE
     SYM_IGNORE TRUE
\$end
EOF

if [[ $do_dipder == "dipder_yes" ]];
then
    $ECHO "Append freq calculation..."

cat >> ipi.scf.in << EOF

@@@

\$molecule
READ
\$end

\$rem
     JOBTYPE  freq
     METHOD   $METHOD
     BASIS    $BASIS
     MAX_SCF_CYCLES   800
     INPUT_BOHR TRUE
     SYM_IGNORE TRUE
     SCF_GUESS     read 
\$end
EOF

fi

$ECHO "  running the SCF for IPI driver...\c"
$QCHEM -nt $NTHREAD ipi.scf.in ipi.scf.out
wait
$ECHO "Done SCF calculation..."

f=ipi.scf.out

# after simulation, I need to get energy, forces, and dipole moment from output
if test -f "energy.au"; then
    $ECHO "$geom_pos" > traj_temp
    sed -i '$ d' traj_temp
    nat=$(wc -l traj_temp | awk '{print $1}')
    $ECHO "$nat\n" >> traj_bohr.xyz
    cat traj_temp >> traj_bohr.xyz

    cat energy.au >> energy_traj
    cat dipole.debye >> dipole_traj
    cat egrad.au >> egrad_traj
    wait

    rm energy.au egrad.au dipole.debye
fi
wait

grep -m 1 -i "Total energy in the final basis set" ipi.scf.out | awk '{print $9}' > energy.au
wait
grep -A 1 -i "Dipole Moment (Debye)" ipi.scf.out | tail -n1 | awk '{print $2" " $4" " $6}' > dipole.debye
wait
sed -n '/Gradient of SCF Energy/, /Max gradient/{ /Gradient of SCF Energy/! { /Max gradient/! p } }' ipi.scf.out | awk 'NR % 4 != 1' | grep -i " 1 " | cut -d' ' -f6- | tr -d '\n' > egrad_1.au
wait
sed -n '/Gradient of SCF Energy/, /Max gradient/{ /Gradient of SCF Energy/! { /Max gradient/! p } }' ipi.scf.out | awk 'NR % 4 != 1' | grep -i " 2 " | cut -d' ' -f6- | tr -d '\n' > egrad_2.au
wait
sed -n '/Gradient of SCF Energy/, /Max gradient/{ /Gradient of SCF Energy/! { /Max gradient/! p } }' ipi.scf.out | awk 'NR % 4 != 1' | grep -i " 3 " | cut -d' ' -f6- | tr -d '\n' > egrad_3.au
echo "" >> egrad_1.au
echo "" >> egrad_2.au
echo "" >> egrad_3.au
wait
cat egrad_1.au egrad_2.au egrad_3.au >> egrad.au
wait
rm egrad_*.au

if [[ $do_dipder == "dipder_yes" ]];
then
    if test -f "dipder.au"; then
	cat dipder.au >> dipder_traj
	rm dipder.au
    fi
    wait
    sed -n '/Derivative of Dipole Moment Matrix/, /------/{ /Derivative of Dipole Moment Matrix/! { /------/! p } }' ipi.scf.out | tail -n +2 | cut -d' ' -f6- > dipder.au
    wait
fi

$ECHO "Done information gathering"

cd ..
