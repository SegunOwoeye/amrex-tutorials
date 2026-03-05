#!/bin/bash

EXEC=./main2d.gnu.MPI.ex
INPUT=configs/lax_and_liu_tests/laxliu_2d

OUTDIR=../../../../../../CPP/data/pilot_timing
mkdir -p $OUTDIR

CSV=$OUTDIR/timing_pilot.csv
echo "mode,resolution,amr_level,cpus,time" > $CSV


# Pilot resolutions
RESOLUTIONS=(64 128 256)
CPUS=(1 2 4 8)

# Uniform grid tests
for RES in "${RESOLUTIONS[@]}"
do
    echo "Running uniform N=$RES"

    for C in "${CPUS[@]}" 
    do
        echo "CPUs=$C"

        START=$(date +%s.%N)

        mpirun -np $C $EXEC $INPUT \
            prob.type=3 \
            geometry.is_periodic="0 0" \
            amr.n_cell="$RES $RES" \
            amr.max_level=0 \
            amr.plot_files_output=0 \
            adv.v=0 \
            amr.v=0

        END=$(date +%s.%N)

        TIME=$(echo "$END - $START" | bc)

        echo "uniform,$RES,0,$C,$TIME" >> $CSV
    done
done



# AMR pilot test
BASE_RES=64

for LEVEL in 1 2
do
    EFFECTIVE=$(( BASE_RES * (2**LEVEL) ))

    echo "Running AMR base=$BASE_RES level=$LEVEL effective=$EFFECTIVE"

    for C in "${CPUS[@]}"
    do

        START=$(date +%s.%N)

        mpirun -np $C $EXEC $INPUT \
            prob.type=3 \
            geometry.is_periodic="0 0" \
            amr.n_cell="$BASE_RES $BASE_RES" \
            amr.max_level=$LEVEL \
            amr.regrid_int=2 \
            amr.n_error_buf=2 \
            adv.max_rhograd_lev=2 \
            amr.plot_files_output=0 \
            adv.v=0 \
            amr.v=0

        END=$(date +%s.%N)

        TIME=$(echo "$END - $START" | bc)

        echo "amr,$EFFECTIVE,$LEVEL,$C,$TIME" >> $CSV
    done
done


echo "Pilot run finished"
echo "Results in $CSV"

